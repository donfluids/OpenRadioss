! Unit test for the JWL EOS module.
!
! Checks the published reference state for TNT (Dobratz-Crawford 1985):
!   At V = 1 (ρ = ρ₀_TNT), e = E₀/ρ₀, the JWL pressure should be ≈ 8.4 GPa.
! Also verifies positivity of (∂p/∂ρ)_e + Grüneisen term (sound speed²)
! at the reference state.
program test_eos_jwl
   use, intrinsic :: iso_fortran_env, only : error_unit
   use kinds,     only : wp
   use eos_jwl,   only : jwl_tnt, jwl_pressure, jwl_sound_speed_squared, &
                         jwl_state_at_reference
   implicit none

   integer :: fails
   real(wp) :: rho, e, p, c2
   real(wp), parameter :: P_REF_TNT     = 8.4e9_wp     ! Pa, expected at V=1
   real(wp), parameter :: P_REF_TOL_REL = 0.05_wp      ! 5 % tolerance on the
                                                       ! literature value

   fails = 0

   ! TNT reference state: ρ = ρ₀, e = E₀/ρ₀.
   call jwl_state_at_reference(jwl_tnt, rho, e, p)
   write(*,'(A,F8.1,A,1PE12.4,A,1PE12.4)') &
      'TNT V=1: rho=', rho, '  e=', e, '  p=', p

   if (abs(p - P_REF_TNT) / P_REF_TNT > P_REF_TOL_REL) then
      write(error_unit,'(A,1PE12.4,A,1PE12.4)') &
         'JWL p at reference too far from 8.4 GPa: got ', p, ' expected ~', P_REF_TNT
      fails = fails + 1
   end if

   c2 = jwl_sound_speed_squared(jwl_tnt, rho, e)
   write(*,'(A,1PE12.4,A,F8.1,A)') 'TNT V=1: c²=', c2, '  c=', sqrt(c2), ' m/s'
   if (c2 <= 0.0_wp) then
      write(error_unit,'(A,1PE12.4)') 'JWL c² non-positive at reference: ', c2
      fails = fails + 1
   end if

   ! At very low density (ρ → 0), p should approach ω·ρ·e → 0.
   p = jwl_pressure(jwl_tnt, 1.0e-3_wp * jwl_tnt%rho0, e)
   write(*,'(A,1PE12.4)') 'JWL p at V → ∞ (low ρ): ', p
   if (p < 0.0_wp) then
      write(error_unit,'(A)') 'JWL p went negative at low density (unphysical)'
      fails = fails + 1
   end if

   if (fails == 0) then
      write(*,'(A)') 'test_eos_jwl: PASS'
   else
      write(*,'(A,I0,A)') 'test_eos_jwl: FAIL (', fails, ')'
      stop 1
   end if
end program test_eos_jwl
