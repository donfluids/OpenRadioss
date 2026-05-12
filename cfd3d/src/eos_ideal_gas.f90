module eos_ideal_gas
   use kinds,     only : wp
   use constants, only : GAMMA, GM1, NVAR, IRHO, IRHOU, IRHOV, IRHOW, IRHOE, &
                         TINY_RHO, TINY_P
   implicit none
   private

   public :: prim_from_cons, cons_from_prim, sound_speed, pressure_from_cons, &
             max_wave_speed

contains

   pure subroutine prim_from_cons(Q, rho, u, v, w, p)
      real(wp), intent(in)  :: Q(NVAR)
      real(wp), intent(out) :: rho, u, v, w, p
      real(wp) :: ke, inv_rho
      rho = max(Q(IRHO), TINY_RHO)
      inv_rho = 1.0_wp / rho
      u = Q(IRHOU) * inv_rho
      v = Q(IRHOV) * inv_rho
      w = Q(IRHOW) * inv_rho
      ke = 0.5_wp * rho * (u*u + v*v + w*w)
      p = max(GM1 * (Q(IRHOE) - ke), TINY_P)
   end subroutine prim_from_cons

   pure subroutine cons_from_prim(rho, u, v, w, p, Q)
      real(wp), intent(in)  :: rho, u, v, w, p
      real(wp), intent(out) :: Q(NVAR)
      real(wp) :: e_int
      e_int = p / GM1
      Q(IRHO)  = rho
      Q(IRHOU) = rho * u
      Q(IRHOV) = rho * v
      Q(IRHOW) = rho * w
      Q(IRHOE) = e_int + 0.5_wp * rho * (u*u + v*v + w*w)
   end subroutine cons_from_prim

   pure function sound_speed(rho, p) result(c)
      real(wp), intent(in) :: rho, p
      real(wp) :: c
      c = sqrt(GAMMA * max(p, TINY_P) / max(rho, TINY_RHO))
   end function sound_speed

   pure function pressure_from_cons(Q) result(p)
      real(wp), intent(in) :: Q(NVAR)
      real(wp) :: p, rho, u, v, w
      call prim_from_cons(Q, rho, u, v, w, p)
   end function pressure_from_cons

   ! Maximum local wave speed |u·n| + c for state Q projected on normal n.
   pure function max_wave_speed(Q, n) result(s)
      real(wp), intent(in) :: Q(NVAR), n(3)
      real(wp) :: s, rho, u, v, w, p, c, vn
      call prim_from_cons(Q, rho, u, v, w, p)
      c = sound_speed(rho, p)
      vn = u*n(1) + v*n(2) + w*n(3)
      s = abs(vn) + c
   end function max_wave_speed

end module eos_ideal_gas
