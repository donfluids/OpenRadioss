module time_integration
   use kinds,         only : wp
   use constants,     only : NVAR
   use mesh_types,    only : t_mesh
   use fields,        only : t_state
   use bc_types,      only : t_bc_data
   use eos_ideal_gas, only : max_wave_speed
   use flux_assembly, only : compute_residual
   implicit none
   private

   public :: compute_dt, rk3_step

contains

   ! Global stable time step from CFL condition:
   !   dt = CFL * min_c V_c / sum_{faces of c} (|u·n| + c) A_f
   function compute_dt(mesh, s, cfl) result(dt)
      type(t_mesh),  intent(in) :: mesh
      type(t_state), intent(in) :: s
      real(wp),      intent(in) :: cfl
      real(wp) :: dt

      integer :: f, c
      real(wp) :: lam_face, n(3), area
      real(wp), allocatable :: lam_sum(:)
      real(wp) :: cell_dt, lam_o, lam_n

      allocate(lam_sum(mesh%nc_total), source=0.0_wp)

      do f = 1, mesh%nf
         n    = mesh%face_normal(:, f)
         area = mesh%face_area(f)
         c = mesh%face_owner(f)
         lam_o = max_wave_speed(s%U(:, c), n)
         if (mesh%face_neighbor(f) > 0) then
            lam_n = max_wave_speed(s%U(:, mesh%face_neighbor(f)), n)
            lam_face = 0.5_wp * (lam_o + lam_n) * area
         else
            lam_face = lam_o * area
         end if
         lam_sum(c) = lam_sum(c) + lam_face
         if (mesh%face_neighbor(f) > 0) then
            lam_sum(mesh%face_neighbor(f)) = lam_sum(mesh%face_neighbor(f)) + lam_face
         end if
      end do

      dt = huge(1.0_wp)
      do c = 1, mesh%nc_internal
         if (lam_sum(c) > 0.0_wp) then
            cell_dt = cfl * mesh%cell_volume(c) / lam_sum(c)
            if (cell_dt < dt) dt = cell_dt
         end if
      end do

      deallocate(lam_sum)
   end function compute_dt

   ! One full step of SSP-RK3 (Shu-Osher form). U is updated in place.
   subroutine rk3_step(mesh, bc_dat, s, dt)
      type(t_mesh),    intent(in)    :: mesh
      type(t_bc_data), intent(in)    :: bc_dat(:)
      type(t_state),   intent(inout) :: s
      real(wp),        intent(in)    :: dt

      integer :: c
      real(wp) :: inv_V

      ! Stash U^n
      s%U0 = s%U

      ! Stage 1: U^(1) = U^n - dt/V R(U^n)
      call compute_residual(mesh, bc_dat, s)
      do c = 1, mesh%nc_internal
         inv_V = 1.0_wp / mesh%cell_volume(c)
         s%U(:, c) = s%U0(:, c) - dt * inv_V * s%R(:, c)
      end do

      ! Stage 2: U^(2) = 3/4 U^n + 1/4 (U^(1) - dt/V R(U^(1)))
      call compute_residual(mesh, bc_dat, s)
      do c = 1, mesh%nc_internal
         inv_V = 1.0_wp / mesh%cell_volume(c)
         s%U(:, c) = 0.75_wp * s%U0(:, c) &
                   + 0.25_wp * ( s%U(:, c) - dt * inv_V * s%R(:, c) )
      end do

      ! Stage 3: U^{n+1} = 1/3 U^n + 2/3 (U^(2) - dt/V R(U^(2)))
      call compute_residual(mesh, bc_dat, s)
      do c = 1, mesh%nc_internal
         inv_V = 1.0_wp / mesh%cell_volume(c)
         s%U(:, c) = (1.0_wp/3.0_wp) * s%U0(:, c) &
                   + (2.0_wp/3.0_wp) * ( s%U(:, c) - dt * inv_V * s%R(:, c) )
      end do
   end subroutine rk3_step

end module time_integration
