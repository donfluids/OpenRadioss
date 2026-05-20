module time_integration
   use kinds,         only : wp
   use constants,     only : NVAR
   use mpi_f08
   use mesh_types,    only : t_mesh
   use fields,        only : t_state
   use bc_types,      only : t_bc_data
   use eos_ideal_gas, only : max_wave_speed, prim_from_cons, temperature, mu_sutherland
   use gas_properties, only : gas, SGS_NONE
   use sgs_model,      only : sgs_nu_t, sgs_filter_width
   use constants,      only : NPRIM
   use flux_assembly, only : residual_begin, residual_pure_interior, &
                             residual_partition, residual_boundary, &
                             residual_fragments
   use gradients,     only : compute_primitives, compute_gradients
   use limiters,      only : compute_venkat_limiter
   use halo_exchange, only : halo_pack_and_start, halo_wait, &
                             halo_pack_and_start_gp, halo_wait_gp
   use mpi_runtime,   only : t_mpi_ctx
   implicit none
   private

   public :: compute_dt, rk3_step

contains

   ! Combined acoustic + (optionally) viscous CFL.
   !   dt_acoustic(c) = CFL * V_c / Σ_faces 0.5(|u·n|+c)_avg A
   !   dt_visc(c)     = CFL * V_c² / Σ_faces κ A²
   !     with κ = max(4μ/(3ρ),  γ μ/(ρ Pr)) — the larger of viscous and thermal
   !     diffusivities.
   function compute_dt(mesh, s, cfl, ctx, viscous) result(dt)
      type(t_mesh),    intent(in) :: mesh
      type(t_state),   intent(in) :: s
      real(wp),        intent(in) :: cfl
      type(t_mpi_ctx), intent(in) :: ctx
      logical,         intent(in) :: viscous
      real(wp) :: dt

      integer :: f, c, c_n
      real(wp) :: lam_face, area, lam_o, lam_n, cell_dt, dt_local, dt_global
      real(wp), allocatable :: lam_sum(:), visc_sum(:)
      real(wp) :: nrml(3)
      real(wp) :: rho, u, v, w, p, T, mu, kappa, kappa_o, kappa_n
      real(wp) :: dt_visc_cell

      allocate(lam_sum(mesh%nc_total),  source=0.0_wp)
      allocate(visc_sum(mesh%nc_total), source=0.0_wp)

      do f = 1, mesh%nf
         nrml = mesh%face_normal(:, f)
         area = mesh%face_area_eff(f)
         c    = mesh%face_owner(f)
         c_n  = mesh%face_neighbor(f)
         lam_o = max_wave_speed(s%U(:, c), nrml)
         if (c_n > 0) then
            lam_n    = max_wave_speed(s%U(:, c_n), nrml)
            lam_face = 0.5_wp * (lam_o + lam_n) * area
         else
            lam_face = lam_o * area
         end if
         lam_sum(c) = lam_sum(c) + lam_face
         if (c_n > 0) lam_sum(c_n) = lam_sum(c_n) + lam_face

         if (viscous) then
            call prim_from_cons(s%U(:, c), rho, u, v, w, p)
            T = temperature(rho, p)
            mu = mu_sutherland(T)
            ! Include SGS eddy viscosity using the gradient from the previous
            ! step (s%gradW, valid after the first complete RK step; zero on
            ! the first step which is the most conservative direction anyway).
            if (gas%sgs_kind /= SGS_NONE) then
               mu = mu + rho * sgs_nu_t(s%gradW(:, :, c), sgs_filter_width(mesh%cell_volume(c)))
            end if
            kappa_o = max(4.0_wp*mu/(3.0_wp*rho), 1.4_wp*mu/(rho*gas%Pr))
            if (c_n > 0) then
               call prim_from_cons(s%U(:, c_n), rho, u, v, w, p)
               T = temperature(rho, p)
               mu = mu_sutherland(T)
               if (gas%sgs_kind /= SGS_NONE) then
                  mu = mu + rho * sgs_nu_t(s%gradW(:, :, c_n), sgs_filter_width(mesh%cell_volume(c_n)))
               end if
               kappa_n = max(4.0_wp*mu/(3.0_wp*rho), 1.4_wp*mu/(rho*gas%Pr))
               kappa = 0.5_wp * (kappa_o + kappa_n)
            else
               kappa = kappa_o
            end if
            visc_sum(c) = visc_sum(c) + kappa * area * area
            if (c_n > 0) visc_sum(c_n) = visc_sum(c_n) + kappa * area * area
         end if
      end do

      dt_local = huge(1.0_wp)
      do c = 1, mesh%nc_internal
         if (lam_sum(c) > 0.0_wp) then
            cell_dt = cfl * mesh%cell_vol_eff(c) / lam_sum(c)
            if (cell_dt < dt_local) dt_local = cell_dt
         end if
         if (viscous .and. visc_sum(c) > 0.0_wp) then
            dt_visc_cell = cfl * mesh%cell_vol_eff(c)**2 / visc_sum(c)
            if (dt_visc_cell < dt_local) dt_local = dt_visc_cell
         end if
      end do
      deallocate(lam_sum, visc_sum)

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
         inv_V = 1.0_wp / mesh%cell_vol_eff(c)
         s%U(:, c) = s%U0(:, c) - dt * inv_V * s%R(:, c)
      end do

      call one_stage(mesh, bc_dat, s, ctx, muscl, viscous, K_venkat)
      do c = 1, mesh%nc_internal
         inv_V = 1.0_wp / mesh%cell_vol_eff(c)
         s%U(:, c) = 0.75_wp * s%U0(:, c) &
                   + 0.25_wp * ( s%U(:, c) - dt * inv_V * s%R(:, c) )
      end do

      call one_stage(mesh, bc_dat, s, ctx, muscl, viscous, K_venkat)
      do c = 1, mesh%nc_internal
         inv_V = 1.0_wp / mesh%cell_vol_eff(c)
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
         call residual_fragments(mesh, s)
      else
         call residual_begin(s)
         call residual_pure_interior(mesh, s, .false., .false.)
         call residual_partition(mesh, s, .false., .false.)
         call residual_boundary(mesh, bc_dat, s, .false., .false.)
         call residual_fragments(mesh, s)
      end if
   end subroutine one_stage

end module time_integration
