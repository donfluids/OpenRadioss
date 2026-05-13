module time_integration
   use kinds,         only : wp
   use constants,     only : NVAR
   use mpi_f08
   use mesh_types,    only : t_mesh
   use fields,        only : t_state
   use bc_types,      only : t_bc_data
   use eos_ideal_gas, only : max_wave_speed
   use flux_assembly, only : residual_begin, residual_pure_interior, &
                             residual_partition, residual_boundary
   use gradients,     only : compute_primitives, compute_gradients
   use limiters,      only : compute_venkat_limiter
   use halo_exchange, only : halo_pack_and_start, halo_wait, &
                             halo_pack_and_start_gp, halo_wait_gp
   use mpi_runtime,   only : t_mpi_ctx
   implicit none
   private

   public :: compute_dt, rk3_step

contains

   function compute_dt(mesh, s, cfl, ctx) result(dt)
      type(t_mesh),    intent(in) :: mesh
      type(t_state),   intent(in) :: s
      real(wp),        intent(in) :: cfl
      type(t_mpi_ctx), intent(in) :: ctx
      real(wp) :: dt

      integer :: f, c
      real(wp) :: lam_face, area, lam_o, lam_n, cell_dt, dt_local, dt_global
      real(wp), allocatable :: lam_sum(:)
      real(wp) :: nrml(3)

      allocate(lam_sum(mesh%nc_total), source=0.0_wp)

      do f = 1, mesh%nf
         nrml = mesh%face_normal(:, f)
         area = mesh%face_area(f)
         c = mesh%face_owner(f)
         lam_o = max_wave_speed(s%U(:, c), nrml)
         if (mesh%face_neighbor(f) > 0) then
            lam_n = max_wave_speed(s%U(:, mesh%face_neighbor(f)), nrml)
            lam_face = 0.5_wp * (lam_o + lam_n) * area
         else
            lam_face = lam_o * area
         end if
         lam_sum(c) = lam_sum(c) + lam_face
         if (mesh%face_neighbor(f) > 0) then
            lam_sum(mesh%face_neighbor(f)) = lam_sum(mesh%face_neighbor(f)) + lam_face
         end if
      end do

      dt_local = huge(1.0_wp)
      do c = 1, mesh%nc_internal
         if (lam_sum(c) > 0.0_wp) then
            cell_dt = cfl * mesh%cell_volume(c) / lam_sum(c)
            if (cell_dt < dt_local) dt_local = cell_dt
         end if
      end do
      deallocate(lam_sum)

      if (ctx%nproc > 1) then
         call MPI_Allreduce(dt_local, dt_global, 1, MPI_DOUBLE_PRECISION, &
                            MPI_MIN, ctx%comm)
         dt = dt_global
      else
         dt = dt_local
      end if
   end function compute_dt

   ! SSP-RK3 (Shu-Osher) full step.
   ! Per stage: halo U → primitives → (grads, limiters) → halo GP overlapped
   ! with pure-interior residual → partition + boundary residuals → RK update.
   subroutine rk3_step(mesh, bc_dat, s, dt, ctx, muscl, viscous, K_venkat)
      type(t_mesh),    intent(inout) :: mesh
      type(t_bc_data), intent(in)    :: bc_dat(:)
      type(t_state),   intent(inout) :: s
      real(wp),        intent(in)    :: dt
      type(t_mpi_ctx), intent(in)    :: ctx
      logical,         intent(in)    :: muscl, viscous
      real(wp),        intent(in)    :: K_venkat

      integer :: c
      real(wp) :: inv_V

      ! Stash U^n
      s%U0 = s%U

      call one_stage(mesh, bc_dat, s, ctx, muscl, viscous, K_venkat)
      do c = 1, mesh%nc_internal
         inv_V = 1.0_wp / mesh%cell_volume(c)
         s%U(:, c) = s%U0(:, c) - dt * inv_V * s%R(:, c)
      end do

      call one_stage(mesh, bc_dat, s, ctx, muscl, viscous, K_venkat)
      do c = 1, mesh%nc_internal
         inv_V = 1.0_wp / mesh%cell_volume(c)
         s%U(:, c) = 0.75_wp * s%U0(:, c) &
                   + 0.25_wp * ( s%U(:, c) - dt * inv_V * s%R(:, c) )
      end do

      call one_stage(mesh, bc_dat, s, ctx, muscl, viscous, K_venkat)
      do c = 1, mesh%nc_internal
         inv_V = 1.0_wp / mesh%cell_volume(c)
         s%U(:, c) = (1.0_wp/3.0_wp) * s%U0(:, c) &
                   + (2.0_wp/3.0_wp) * ( s%U(:, c) - dt * inv_V * s%R(:, c) )
      end do
   end subroutine rk3_step

   subroutine one_stage(mesh, bc_dat, s, ctx, muscl, viscous, K_venkat)
      type(t_mesh),    intent(inout) :: mesh
      type(t_bc_data), intent(in)    :: bc_dat(:)
      type(t_state),   intent(inout) :: s
      type(t_mpi_ctx), intent(in)    :: ctx
      logical,         intent(in)    :: muscl, viscous
      real(wp),        intent(in)    :: K_venkat
      logical :: need_grads

      need_grads = muscl .or. viscous

      ! 1) Exchange U so ghosts have current state.
      call halo_pack_and_start(mesh, s%U, ctx)
      call halo_wait(mesh, s%U)

      ! 2) Primitives from cons (all cells incl ghosts).
      call compute_primitives(mesh, s)

      ! 3) Per-cell gradients and (if MUSCL) limiters on local cells.
      if (need_grads) then
         call compute_gradients(mesh, bc_dat, s)
         if (muscl) then
            call compute_venkat_limiter(mesh, s, K_venkat)
         else
            s%psi = 1.0_wp
         end if
         ! 4) Exchange grads + limiters; overlap with pure-interior residual.
         call halo_pack_and_start_gp(mesh, s, ctx)
         call residual_begin(s)
         call residual_pure_interior(mesh, s, muscl, viscous)
         call halo_wait_gp(mesh, s)
         call residual_partition(mesh, s, muscl, viscous)
         call residual_boundary(mesh, bc_dat, s, muscl, viscous)
      else
         call residual_begin(s)
         call residual_pure_interior(mesh, s, .false., .false.)
         call residual_partition(mesh, s, .false., .false.)
         call residual_boundary(mesh, bc_dat, s, .false., .false.)
      end if
   end subroutine one_stage

end module time_integration
