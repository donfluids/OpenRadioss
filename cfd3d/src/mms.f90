! Method of Manufactured Solutions for spatial-order verification.
!
! Manufactured solution (steady, satisfies slip wall BCs u·n = 0 at all box
! boundaries of [0,1]^3):
!     ρ_a(x)   = ρ0 = 1
!     u_a(x)   = sin(πx)
!     v_a(x)   = 0
!     w_a(x)   = 0
!     p_a(x)   = p0 = 10
!     T_a       = p0 / (ρ0 R)   constant
!
! Body-force source terms (derived from steady NS with this W_a):
!     S_ρ      = π cos(πx)
!     S_(ρu)   = π sin(2πx)       + (4π²μ/3) sin(πx)
!     S_(ρv)   = 0
!     S_(ρw)   = 0
!     S_(ρE)   = π cos(πx) (γ p₀/(γ-1) + ρ₀ sin²(πx))    (inviscid energy advection)
!              + π cos(πx) ρ₀ sin²(πx) / 2               (kinetic-energy gradient)
!              - (4π²μ/3) cos(2πx)                       (viscous heating divergence)
!
! At T = const, ∇T = 0 so heat flux is zero.
module mms
   use kinds,          only : wp
   use constants,      only : NVAR, NPRIM, GAMMA, GM1, IRHO, IRHOU, IRHOV, IRHOW, IRHOE, &
                              IP_RHO, IP_U, IP_V, IP_W, IP_P, PI
   use gas_properties, only : gas
   implicit none
   private

   public :: mms_primitives, mms_source
   public :: MMS_RHO0, MMS_P0

   real(wp), parameter :: MMS_RHO0 = 1.0_wp
   real(wp), parameter :: MMS_P0   = 10.0_wp

contains

   pure subroutine mms_primitives(x, y, z, W)
      real(wp), intent(in)  :: x, y, z
      real(wp), intent(out) :: W(NPRIM)
      W(IP_RHO) = MMS_RHO0
      W(IP_U)   = sin(PI * x)
      W(IP_V)   = 0.0_wp
      W(IP_W)   = 0.0_wp
      W(IP_P)   = MMS_P0
      ! Suppress unused
      if (.false.) W(IP_RHO) = W(IP_RHO) + y + z
   end subroutine mms_primitives

   ! S(x,y,z): per-conservative-variable body force (units: conservative / time).
   pure subroutine mms_source(x, y, z, S)
      real(wp), intent(in)  :: x, y, z
      real(wp), intent(out) :: S(NVAR)
      real(wp) :: mu, c2px, c2px2, s2px, spx, cpx

      spx  = sin(PI * x)
      cpx  = cos(PI * x)
      s2px = sin(2.0_wp * PI * x)
      c2px = cos(2.0_wp * PI * x)

      ! Constant-μ approximation: viscosity at T = p0/(ρ0*R) is the same
      ! everywhere for this manufactured solution.
      mu = mu_at_state()

      S(IRHO)  = PI * cpx
      S(IRHOU) = PI * s2px + (4.0_wp * PI * PI * mu / 3.0_wp) * spx
      S(IRHOV) = 0.0_wp
      S(IRHOW) = 0.0_wp

      ! Energy source — derived from ∇·(ρu(E + p/ρ)) - ∇·(τ·u - q):
      !   ∇·(ρu H) = π cos(πx) * (γ p₀/(γ-1) + 1.5 ρ₀ sin²(πx))
      !   ∇·(τ·u)  = (4π²μ/3) cos(2πx)
      !   q = 0 (T constant)
      ! (We split S_E into the convective + kinetic-energy-gradient + viscous parts.)
      c2px2 = 0.0_wp     ! unused placeholder removed below
      S(IRHOE) = PI * cpx * ( GAMMA * MMS_P0 / GM1 + 1.5_wp * MMS_RHO0 * spx * spx ) &
               - (4.0_wp * PI * PI * mu / 3.0_wp) * c2px

      ! Suppress unused
      if (.false.) S(1) = S(1) + y + z + c2px2
   end subroutine mms_source

   pure function mu_at_state() result(mu)
      real(wp) :: mu, T
      if (gas%sutherland) then
         T = MMS_P0 / (MMS_RHO0 * gas%R_gas)
         mu = gas%mu_ref * (T / gas%T_ref) ** 1.5_wp &
              * (gas%T_ref + gas%S_S) / (T + gas%S_S)
      else
         mu = gas%mu_const
      end if
   end function mu_at_state

end module mms
