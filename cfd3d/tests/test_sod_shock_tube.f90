! Integration acceptance test: 3D Sod shock tube vs. exact Riemann solution.
program test_sod_shock_tube
   use, intrinsic :: iso_fortran_env, only : error_unit
   use kinds,            only : wp
   use constants,        only : GAMMA, NVAR
   use mesh_types,       only : t_mesh
   use mesh_module,      only : load_mesh
   use fields,           only : t_state, alloc_state, free_state
   use bc_types,         only : t_bc_data
   use eos_ideal_gas,    only : prim_from_cons
   use time_integration, only : compute_dt, rk3_step
   use solver_control,   only : t_run_params, read_namelist
   use solver_driver,    only : assign_patch_bcs, set_initial_condition
   implicit none

   character(len=512) :: nml_path
   integer :: nargs
   type(t_run_params) :: p
   type(t_mesh)       :: mesh
   type(t_state)      :: s
   type(t_bc_data), allocatable :: bc_dat(:)
   real(wp) :: t, dt, x, x_d
   real(wp) :: rho_n, u_n, v_n, w_n, p_n
   real(wp) :: rho_e, u_e, p_e
   real(wp) :: err_rho, err_u, err_p
   real(wp) :: tot_rho, tot_u, tot_p
   real(wp) :: rel_rho, rel_u, rel_p
   integer :: c, step
   integer :: fails

   nargs = command_argument_count()
   if (nargs < 1) then
      write(error_unit,'(A)') 'usage: test_sod_shock_tube <input.nml>'
      stop 1
   end if
   call get_command_argument(1, nml_path)

   call read_namelist(trim(nml_path), p)
   call load_mesh(trim(p%mesh_file), mesh)
   call assign_patch_bcs(p, mesh, bc_dat)
   call alloc_state(s, mesh)
   call set_initial_condition(p, mesh, s)

   t = 0.0_wp
   step = 0
   do while (t < p%t_end .and. step < p%max_steps)
      dt = compute_dt(mesh, s, p%cfl)
      if (t + dt > p%t_end) dt = p%t_end - t
      call rk3_step(mesh, bc_dat, s, dt)
      t = t + dt
      step = step + 1
   end do

   ! Compare to exact Sod at t = p%t_end along the diaphragm axis.
   x_d = p%diaphragm_pos
   err_rho = 0.0_wp; err_u = 0.0_wp; err_p = 0.0_wp
   tot_rho = 0.0_wp; tot_u = 0.0_wp; tot_p = 0.0_wp

   do c = 1, mesh%nc_internal
      x = mesh%cell_centroid(p%diaphragm_axis, c)
      call prim_from_cons(s%U(:, c), rho_n, u_n, v_n, w_n, p_n)
      call exact_sod(p%rho_L, p%u_L, p%p_L, p%rho_R, p%u_R, p%p_R, &
                     x - x_d, t, rho_e, u_e, p_e)
      err_rho = err_rho + abs(rho_n - rho_e) * mesh%cell_volume(c)
      err_u   = err_u   + abs(u_n   - u_e)   * mesh%cell_volume(c)
      err_p   = err_p   + abs(p_n   - p_e)   * mesh%cell_volume(c)
      tot_rho = tot_rho + abs(rho_e)         * mesh%cell_volume(c)
      tot_u   = tot_u   + abs(u_e)           * mesh%cell_volume(c)
      tot_p   = tot_p   + abs(p_e)           * mesh%cell_volume(c)
   end do

   rel_rho = err_rho / max(tot_rho, 1.0e-30_wp)
   rel_p   = err_p   / max(tot_p,   1.0e-30_wp)
   ! u reference can be zero almost everywhere; use absolute L1 against domain
   rel_u   = err_u / max(tot_u, 1.0e-30_wp)

   write(*,'(A,1PE12.4)') 'L1_rel(rho) = ', rel_rho
   write(*,'(A,1PE12.4)') 'L1_rel(p)   = ', rel_p
   write(*,'(A,1PE12.4)') 'L1_abs(u)/L1(u_exact) = ', rel_u

   fails = 0
   if (rel_rho > 0.06_wp) then
      write(error_unit,'(A,1PE12.4)') 'rel L1 rho exceeds 6% : ', rel_rho
      fails = fails + 1
   end if
   if (rel_p > 0.06_wp) then
      write(error_unit,'(A,1PE12.4)') 'rel L1 p exceeds 6% : ', rel_p
      fails = fails + 1
   end if

   call free_state(s)
   if (allocated(bc_dat)) deallocate(bc_dat)

   if (fails == 0) then
      write(*,'(A)') 'test_sod_shock_tube: PASS'
   else
      write(*,'(A,I0,A)') 'test_sod_shock_tube: FAIL (', fails, ')'
      stop 1
   end if

contains

   ! Exact solution to the 1D Riemann problem (Sod-type) for an ideal gas.
   ! Returns (rho, u, p) at offset xi = x - x_d at time t (>0).
   subroutine exact_sod(rhoL, uL, pL, rhoR, uR, pR, xi, t, rho, u, pres)
      real(wp), intent(in)  :: rhoL, uL, pL, rhoR, uR, pR, xi, t
      real(wp), intent(out) :: rho, u, pres

      real(wp) :: cL, cR, p_star, u_star
      real(wp) :: SL, SR, S_HL, S_TL, S_HR, S_TR, S_C
      real(wp) :: rho_starL, rho_starR
      real(wp) :: s
      real(wp) :: c
      real(wp), parameter :: gm1 = GAMMA - 1.0_wp
      real(wp), parameter :: gp1 = GAMMA + 1.0_wp

      cL = sqrt(GAMMA * pL / rhoL)
      cR = sqrt(GAMMA * pR / rhoR)

      call solve_p_star(rhoL, uL, pL, cL, rhoR, uR, pR, cR, p_star, u_star)

      if (t <= 0.0_wp) then
         if (xi < 0.0_wp) then
            rho = rhoL; u = uL; pres = pL
         else
            rho = rhoR; u = uR; pres = pR
         end if
         return
      end if

      s = xi / t

      ! Left wave
      if (p_star > pL) then
         ! Left shock
         SL = uL - cL * sqrt( 0.5_wp*gp1/GAMMA * p_star/pL + 0.5_wp*gm1/GAMMA )
         if (s < SL) then
            rho = rhoL; u = uL; pres = pL; return
         end if
         rho_starL = rhoL * ( p_star/pL + gm1/gp1 ) / ( (gm1/gp1) * p_star/pL + 1.0_wp )
      else
         ! Left rarefaction
         S_HL = uL - cL
         S_TL = u_star - cL * (p_star/pL) ** (0.5_wp*gm1/GAMMA)
         if (s < S_HL) then
            rho = rhoL; u = uL; pres = pL; return
         end if
         if (s < S_TL) then
            ! Inside fan
            u   = (2.0_wp/gp1) * ( cL + 0.5_wp*gm1*uL + s )
            c   = (2.0_wp/gp1) * ( cL + 0.5_wp*gm1*(uL - s) )
            rho = rhoL * (c/cL) ** (2.0_wp/gm1)
            pres = pL  * (c/cL) ** (2.0_wp*GAMMA/gm1)
            return
         end if
         rho_starL = rhoL * (p_star/pL) ** (1.0_wp/GAMMA)
      end if

      S_C = u_star

      ! Right wave
      if (p_star > pR) then
         SR = uR + cR * sqrt( 0.5_wp*gp1/GAMMA * p_star/pR + 0.5_wp*gm1/GAMMA )
         if (s > SR) then
            rho = rhoR; u = uR; pres = pR; return
         end if
         rho_starR = rhoR * ( p_star/pR + gm1/gp1 ) / ( (gm1/gp1) * p_star/pR + 1.0_wp )
      else
         S_HR = uR + cR
         S_TR = u_star + cR * (p_star/pR) ** (0.5_wp*gm1/GAMMA)
         if (s > S_HR) then
            rho = rhoR; u = uR; pres = pR; return
         end if
         if (s > S_TR) then
            u   = (2.0_wp/gp1) * ( -cR + 0.5_wp*gm1*uR + s )
            c   = (2.0_wp/gp1) * (  cR - 0.5_wp*gm1*(uR - s) )
            rho = rhoR * (c/cR) ** (2.0_wp/gm1)
            pres = pR  * (c/cR) ** (2.0_wp*GAMMA/gm1)
            return
         end if
         rho_starR = rhoR * (p_star/pR) ** (1.0_wp/GAMMA)
      end if

      ! Between waves: star region
      if (s < S_C) then
         rho = rho_starL; u = u_star; pres = p_star
      else
         rho = rho_starR; u = u_star; pres = p_star
      end if
   end subroutine exact_sod

   subroutine solve_p_star(rhoL, uL, pL, cL, rhoR, uR, pR, cR, p_star, u_star)
      real(wp), intent(in)  :: rhoL, uL, pL, cL, rhoR, uR, pR, cR
      real(wp), intent(out) :: p_star, u_star

      real(wp) :: p_old, fL, fR, fLp, fRp, delta
      integer  :: iter
      integer, parameter :: MAX_ITER = 100
      real(wp), parameter :: TOL = 1.0e-12_wp

      ! Initial guess: PVRS
      p_star = max(1.0e-10_wp, 0.5_wp*(pL + pR) - 0.5_wp*(uR-uL)*0.5_wp*(rhoL+rhoR)*0.5_wp*(cL+cR))

      do iter = 1, MAX_ITER
         call f_wave(p_star, pL, rhoL, cL, fL, fLp)
         call f_wave(p_star, pR, rhoR, cR, fR, fRp)
         delta = (fL + fR + (uR - uL)) / (fLp + fRp)
         p_old = p_star
         p_star = max(1.0e-10_wp, p_star - delta)
         if (abs(p_star - p_old) / max(p_star, 1.0e-10_wp) < TOL) exit
      end do

      u_star = 0.5_wp * (uL + uR) + 0.5_wp * (fR - fL)
   end subroutine solve_p_star

   subroutine f_wave(p, pK, rhoK, cK, f, fp)
      real(wp), intent(in)  :: p, pK, rhoK, cK
      real(wp), intent(out) :: f, fp
      real(wp) :: A, B, rat, gm1, gp1
      gm1 = GAMMA - 1.0_wp
      gp1 = GAMMA + 1.0_wp
      if (p > pK) then
         ! Shock
         A = 2.0_wp / (gp1 * rhoK)
         B = pK * gm1 / gp1
         f  = (p - pK) * sqrt(A / (p + B))
         fp = sqrt(A / (p + B)) * (1.0_wp - 0.5_wp*(p - pK)/(p + B))
      else
         ! Rarefaction
         rat = p / pK
         f  = (2.0_wp * cK / gm1) * ( rat ** (0.5_wp*gm1/GAMMA) - 1.0_wp )
         fp = (1.0_wp / (rhoK*cK)) * rat ** (-0.5_wp*gp1/GAMMA)
      end if
   end subroutine f_wave

end program test_sod_shock_tube
