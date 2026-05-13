module time_integration
   use kinds,         only : wp
   use constants,     only : NVAR
   use mpi_f08
   use mesh_types,    only : t_mesh
   use fields,        only : t_state
   use bc_types,      only : t_bc_data
   use eos_ideal_gas, only : max_wave_speed
   use flux_assembly, only : compute_residual
   use mpi_runtime,   only : t_mpi_ctx
   implicit none
   private

   public :: compute_dt, rk3_step

contains

   ! Global stable time step (CFL). Local min over each rank's cells, then
   ! MPI_Allreduce(MIN) to get the global minimum.
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

   ! SSP-RK3 (Shu-Osher) full step. Each stage exchanges halos before residual.
   subroutine rk3_step(mesh, bc_dat, s, dt, ctx)
      type(t_mesh),    intent(inout) :: mesh
      type(t_bc_data), intent(in)    :: bc_dat(:)
      type(t_state),   intent(inout) :: s
      real(wp),        intent(in)    :: dt
      type(t_mpi_ctx), intent(in)    :: ctx

      integer :: c
      real(wp) :: inv_V

      ! Stash U^n
      s%U0 = s%U

      call compute_residual(mesh, bc_dat, s, ctx)
      do c = 1, mesh%nc_internal
         inv_V = 1.0_wp / mesh%cell_volume(c)
         s%U(:, c) = s%U0(:, c) - dt * inv_V * s%R(:, c)
      end do

      call compute_residual(mesh, bc_dat, s, ctx)
      do c = 1, mesh%nc_internal
         inv_V = 1.0_wp / mesh%cell_volume(c)
         s%U(:, c) = 0.75_wp * s%U0(:, c) &
                   + 0.25_wp * ( s%U(:, c) - dt * inv_V * s%R(:, c) )
      end do

      call compute_residual(mesh, bc_dat, s, ctx)
      do c = 1, mesh%nc_internal
         inv_V = 1.0_wp / mesh%cell_volume(c)
         s%U(:, c) = (1.0_wp/3.0_wp) * s%U0(:, c) &
                   + (2.0_wp/3.0_wp) * ( s%U(:, c) - dt * inv_V * s%R(:, c) )
      end do
   end subroutine rk3_step

end module time_integration
