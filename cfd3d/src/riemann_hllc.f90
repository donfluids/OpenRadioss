module riemann_hllc
   use kinds,         only : wp
   use constants,     only : NVAR, IRHO, IRHOU, IRHOV, IRHOW, IRHOE, &
                             GAMMA, GM1, TINY_RHO, TINY_P
   use eos_ideal_gas, only : prim_from_cons
   implicit none
   private

   public :: hllc_flux

contains

   ! HLLC Riemann flux through a face with unit normal n, given left/right
   ! conservative states. Reference: Toro, "Riemann Solvers and Numerical
   ! Methods for Fluid Dynamics", 3rd ed., Sec. 10.6.
   pure function hllc_flux(QL, QR, n) result(F)
      real(wp), intent(in) :: QL(NVAR), QR(NVAR), n(3)
      real(wp) :: F(NVAR)

      real(wp) :: rhoL, uL, vL, wL, pL, EL, cL, vnL
      real(wp) :: rhoR, uR, vR, wR, pR, ER, cR, vnR
      real(wp) :: rho_avg, c_avg, p_pvrs, p_star
      real(wp) :: SL, SR, SM
      real(wp) :: rhoLs, rhoRs, factL, factR
      real(wp) :: vnL_factor, vnR_factor
      real(wp) :: QsL(NVAR), QsR(NVAR), FL(NVAR), FR(NVAR)

      call prim_from_cons(QL, rhoL, uL, vL, wL, pL)
      call prim_from_cons(QR, rhoR, uR, vR, wR, pR)

      EL = QL(IRHOE) / max(rhoL, TINY_RHO)
      ER = QR(IRHOE) / max(rhoR, TINY_RHO)

      cL = sqrt(GAMMA * max(pL, TINY_P) / max(rhoL, TINY_RHO))
      cR = sqrt(GAMMA * max(pR, TINY_P) / max(rhoR, TINY_RHO))

      vnL = uL*n(1) + vL*n(2) + wL*n(3)
      vnR = uR*n(1) + vR*n(2) + wR*n(3)

      ! Pressure estimate via PVRS (Toro 10.59)
      rho_avg = 0.5_wp * (rhoL + rhoR)
      c_avg   = 0.5_wp * (cL + cR)
      p_pvrs  = 0.5_wp * (pL + pR) - 0.5_wp * (vnR - vnL) * rho_avg * c_avg
      p_star  = max(0.0_wp, p_pvrs)

      ! Wave speed estimates (Toro 10.59-10.60, Einfeldt-style)
      if (p_star <= pL) then
         vnL_factor = 1.0_wp
      else
         vnL_factor = sqrt(1.0_wp + 0.5_wp*(GAMMA+1.0_wp)/GAMMA * (p_star/pL - 1.0_wp))
      end if
      if (p_star <= pR) then
         vnR_factor = 1.0_wp
      else
         vnR_factor = sqrt(1.0_wp + 0.5_wp*(GAMMA+1.0_wp)/GAMMA * (p_star/pR - 1.0_wp))
      end if
      SL = vnL - cL * vnL_factor
      SR = vnR + cR * vnR_factor

      ! Contact wave speed
      SM = ( pR - pL + rhoL*vnL*(SL - vnL) - rhoR*vnR*(SR - vnR) ) &
           / ( rhoL*(SL - vnL) - rhoR*(SR - vnR) )

      call phys_flux(QL, n, FL)
      call phys_flux(QR, n, FR)

      if (SL >= 0.0_wp) then
         F = FL
      else if (SR <= 0.0_wp) then
         F = FR
      else if (SM >= 0.0_wp) then
         factL = rhoL * (SL - vnL) / (SL - SM)
         rhoLs = factL
         QsL(IRHO)  = rhoLs
         QsL(IRHOU) = rhoLs * (uL + (SM - vnL)*n(1))
         QsL(IRHOV) = rhoLs * (vL + (SM - vnL)*n(2))
         QsL(IRHOW) = rhoLs * (wL + (SM - vnL)*n(3))
         QsL(IRHOE) = rhoLs * ( EL + (SM - vnL) * (SM + pL / (rhoL * (SL - vnL))) )
         F = FL + SL * (QsL - QL)
      else
         factR = rhoR * (SR - vnR) / (SR - SM)
         rhoRs = factR
         QsR(IRHO)  = rhoRs
         QsR(IRHOU) = rhoRs * (uR + (SM - vnR)*n(1))
         QsR(IRHOV) = rhoRs * (vR + (SM - vnR)*n(2))
         QsR(IRHOW) = rhoRs * (wR + (SM - vnR)*n(3))
         QsR(IRHOE) = rhoRs * ( ER + (SM - vnR) * (SM + pR / (rhoR * (SR - vnR))) )
         F = FR + SR * (QsR - QR)
      end if
   end function hllc_flux

   pure subroutine phys_flux(Q, n, F)
      real(wp), intent(in)  :: Q(NVAR), n(3)
      real(wp), intent(out) :: F(NVAR)
      real(wp) :: rho, u, v, w, p, vn, H, E
      call prim_from_cons(Q, rho, u, v, w, p)
      vn = u*n(1) + v*n(2) + w*n(3)
      E  = Q(IRHOE) / max(rho, TINY_RHO)
      H  = E + p / max(rho, TINY_RHO)
      F(IRHO)  = rho * vn
      F(IRHOU) = rho * u * vn + p * n(1)
      F(IRHOV) = rho * v * vn + p * n(2)
      F(IRHOW) = rho * w * vn + p * n(3)
      F(IRHOE) = rho * H * vn
   end subroutine phys_flux

end module riemann_hllc
