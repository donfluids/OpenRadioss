program test_riemann_hllc
   use, intrinsic :: iso_fortran_env, only : error_unit
   use kinds,         only : wp
   use constants,     only : NVAR, IRHO, IRHOU, IRHOE
   use eos_ideal_gas, only : cons_from_prim
   use riemann_hllc,  only : hllc_flux
   implicit none

   integer :: fails
   real(wp) :: UL(NVAR), UR(NVAR), F(NVAR), n(3)
   real(wp) :: rho, u, v, w, p
   real(wp), parameter :: TOL = 1.0e-12_wp

   fails = 0

   ! Test 1: Identical states → flux equals the physical flux through n.
   call cons_from_prim(1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, UL)
   UR = UL
   n  = [1.0_wp, 0.0_wp, 0.0_wp]
   F  = hllc_flux(UL, UR, n)
   ! rho*u = 0, rho*u^2 + p = 1, rho*u*v = 0, rho*u*w = 0, rho*u*H = 0
   if (abs(F(IRHO))  > TOL) then; fails = fails + 1; write(error_unit,*) 'identical: F(rho) /= 0', F(IRHO); end if
   if (abs(F(IRHOU) - 1.0_wp) > TOL) then; fails = fails + 1; write(error_unit,*) 'identical: F(rhou) /= p', F(IRHOU); end if
   if (abs(F(IRHOE)) > TOL) then; fails = fails + 1; write(error_unit,*) 'identical: F(rhoE) /= 0', F(IRHOE); end if

   ! Test 2: Rotational invariance — same problem rotated should give rotated fluxes.
   ! Setup: Sod problem in +x, then in +y; the rho flux magnitude must match.
   block
      real(wp) :: UL1(NVAR), UR1(NVAR), F1(NVAR)
      real(wp) :: UL2(NVAR), UR2(NVAR), F2(NVAR)
      call cons_from_prim(1.0_wp,    0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, UL1)
      call cons_from_prim(0.125_wp,  0.0_wp, 0.0_wp, 0.0_wp, 0.1_wp, UR1)
      F1 = hllc_flux(UL1, UR1, [1.0_wp, 0.0_wp, 0.0_wp])

      ! Same problem, but velocities rotated x→y; normal also rotated
      call cons_from_prim(1.0_wp,    0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, UL2)
      call cons_from_prim(0.125_wp,  0.0_wp, 0.0_wp, 0.0_wp, 0.1_wp, UR2)
      F2 = hllc_flux(UL2, UR2, [0.0_wp, 1.0_wp, 0.0_wp])
      if (abs(F1(IRHO) - F2(IRHO)) > 1.0e-12_wp) then
         fails = fails + 1
         write(error_unit,*) 'rot invariance F(rho)', F1(IRHO), F2(IRHO)
      end if
      if (abs(F1(IRHOE) - F2(IRHOE)) > 1.0e-12_wp) then
         fails = fails + 1
         write(error_unit,*) 'rot invariance F(rhoE)', F1(IRHOE), F2(IRHOE)
      end if
   end block

   ! Test 3: Wall-reflection — UR is mirror of UL across normal → mass flux = 0
   block
      real(wp) :: UR_ref(NVAR)
      n = [1.0_wp, 0.0_wp, 0.0_wp]
      call cons_from_prim(1.0_wp, 0.5_wp, 0.2_wp, -0.1_wp, 1.0_wp, UL)
      call cons_from_prim(1.0_wp,-0.5_wp, 0.2_wp, -0.1_wp, 1.0_wp, UR_ref)
      F = hllc_flux(UL, UR_ref, n)
      ! Mass flux should be approximately 0 by symmetry across a slip wall.
      if (abs(F(IRHO)) > 1.0e-10_wp) then
         fails = fails + 1
         write(error_unit,*) 'wall-mirror mass flux', F(IRHO)
      end if
      ! Energy flux should also vanish.
      if (abs(F(IRHOE)) > 1.0e-10_wp) then
         fails = fails + 1
         write(error_unit,*) 'wall-mirror energy flux', F(IRHOE)
      end if
   end block

   ! Suppress unused
   if (.false.) then
      rho = 0.0_wp; u = 0.0_wp; v = 0.0_wp; w = 0.0_wp; p = 0.0_wp
   end if

   if (fails == 0) then
      write(*,'(A)') 'test_riemann_hllc: PASS'
   else
      write(*,'(A,I0,A)') 'test_riemann_hllc: FAIL (', fails, ' issues)'
      stop 1
   end if
end program test_riemann_hllc
