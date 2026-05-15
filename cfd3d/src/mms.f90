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
!
! M4 SGS extension: when gas%sgs_kind != SGS_NONE, mms_source adds an
! additional contribution -∂(τ_xx_SGS)/∂x to S_(ρu) and -∂(u·τ_xx_SGS)/∂x to
! S_(ρE). Both are computed via 4th-order central finite differences of the
! analytical τ_xx_SGS along x (the only non-zero stress component for this
! 1-D manufactured velocity field). FD step h=1e-5 keeps FD error ~ h⁴ ≈
! 1e-20, well below any mesh-discretization error of interest.
module mms
   use kinds,          only : wp
   use constants,      only : NVAR, NPRIM, GAMMA, GM1, IRHO, IRHOU, IRHOV, IRHOW, IRHOE, &
                              IP_RHO, IP_U, IP_V, IP_W, IP_P, PI
   use gas_properties, only : gas, SGS_NONE, SGS_SMAG, SGS_WALE
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

   ! S(x,y,z; Delta): per-conservative-variable body force.
   ! Delta is the per-cell filter width (V_cell^(1/3)); only used when the
   ! SGS model is active and is otherwise ignored.
   pure subroutine mms_source(x, y, z, Delta, S)
      real(wp), intent(in)  :: x, y, z, Delta
      real(wp), intent(out) :: S(NVAR)
      real(wp) :: mu, c2px, s2px, spx, cpx
      real(wp) :: dtau_dx, du_tau_dx

      spx  = sin(PI * x)
      cpx  = cos(PI * x)
      s2px = sin(2.0_wp * PI * x)
      c2px = cos(2.0_wp * PI * x)

      mu = mu_at_state()

      S(IRHO)  = PI * cpx
      S(IRHOU) = PI * s2px + (4.0_wp * PI * PI * mu / 3.0_wp) * spx
      S(IRHOV) = 0.0_wp
      S(IRHOW) = 0.0_wp
      S(IRHOE) = PI * cpx * ( GAMMA * MMS_P0 / GM1 + 1.5_wp * MMS_RHO0 * spx * spx ) &
               - (4.0_wp * PI * PI * mu / 3.0_wp) * c2px

      ! SGS contribution (added when sgs_kind != NONE).
      if (gas%sgs_kind /= SGS_NONE) then
         call sgs_source_at(x, Delta, dtau_dx, du_tau_dx)
         S(IRHOU) = S(IRHOU) - dtau_dx
         S(IRHOE) = S(IRHOE) - du_tau_dx
      end if

      ! Suppress unused
      if (.false.) S(1) = S(1) + y + z
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

   ! ν_t evaluated analytically from the manufactured solution at point x,
   ! using the same model formula the discretization sees, with filter width
   ! Delta. For our 1-D manufactured velocity (sin(πx), 0, 0) only ∂u/∂x is
   ! non-zero, so |S| and (S^d:S^d), (S:S) all simplify to closed form.
   pure function nu_t_at(x, Delta) result(nu_t)
      real(wp), intent(in) :: x, Delta
      real(wp) :: nu_t
      real(wp) :: cpx_abs, SmagS, SijSij, SdSd, num, den
      cpx_abs = abs(cos(PI * x))
      SmagS = sqrt(2.0_wp) * PI * cpx_abs    ! √(2 S_ij S_ij)
      SijSij = (PI * cpx_abs)**2             ! S_ij S_ij (only S_11^2)
      ! WALE invariants:
      !   trace(g²) = π² cos²
      !   Sd_11 = 2π² cos² / 3, Sd_22 = Sd_33 = -π² cos² / 3, off-diagonals = 0
      !   S^d : S^d = 2 π⁴ cos⁴ / 3
      SdSd = (2.0_wp / 3.0_wp) * (PI * cpx_abs)**4
      select case (gas%sgs_kind)
      case (SGS_SMAG)
         nu_t = (gas%C_s * Delta)**2 * SmagS
      case (SGS_WALE)
         num = SdSd ** 1.5_wp
         den = (SijSij ** 2.5_wp) + (SdSd ** 1.25_wp)
         if (den <= 1.0e-30_wp) then
            nu_t = 0.0_wp
         else
            nu_t = (gas%C_w * Delta)**2 * num / den
         end if
      case default
         nu_t = 0.0_wp
      end select
   end function nu_t_at

   pure function tau_xx_sgs_at(x, Delta) result(tau_xx)
      real(wp), intent(in) :: x, Delta
      real(wp) :: tau_xx, nu_t, dudx, div_u
      nu_t  = nu_t_at(x, Delta)
      dudx  = PI * cos(PI * x)
      div_u = dudx
      ! Boussinesq: τ_SGS = ρ ν_t (2 S_ij - (2/3)(∇·u) δ_ij). For our 1-D
      ! velocity, S_11 = ∂u/∂x and ∇·u = ∂u/∂x, so τ_xx = ρ ν_t (4/3) ∂u/∂x.
      tau_xx = MMS_RHO0 * nu_t * (4.0_wp / 3.0_wp) * dudx
   end function tau_xx_sgs_at

   pure function u_tau_xx_sgs_at(x, Delta) result(uTau)
      real(wp), intent(in) :: x, Delta
      real(wp) :: uTau
      ! u · τ_xx_SGS — used for the energy source ∂(u·τ_xx)/∂x.
      uTau = sin(PI * x) * tau_xx_sgs_at(x, Delta)
   end function u_tau_xx_sgs_at

   ! 4th-order central FD of ∂τ_xx_SGS/∂x and ∂(u·τ_xx_SGS)/∂x at point x.
   pure subroutine sgs_source_at(x, Delta, dtau_dx, du_tau_dx)
      real(wp), intent(in)  :: x, Delta
      real(wp), intent(out) :: dtau_dx, du_tau_dx
      real(wp), parameter :: h = 1.0e-5_wp
      real(wp) :: t_p2, t_p1, t_m1, t_m2
      real(wp) :: ut_p2, ut_p1, ut_m1, ut_m2

      t_p2 = tau_xx_sgs_at(x + 2.0_wp*h, Delta)
      t_p1 = tau_xx_sgs_at(x +        h, Delta)
      t_m1 = tau_xx_sgs_at(x -        h, Delta)
      t_m2 = tau_xx_sgs_at(x - 2.0_wp*h, Delta)
      dtau_dx = (-t_p2 + 8.0_wp*t_p1 - 8.0_wp*t_m1 + t_m2) / (12.0_wp * h)

      ut_p2 = u_tau_xx_sgs_at(x + 2.0_wp*h, Delta)
      ut_p1 = u_tau_xx_sgs_at(x +        h, Delta)
      ut_m1 = u_tau_xx_sgs_at(x -        h, Delta)
      ut_m2 = u_tau_xx_sgs_at(x - 2.0_wp*h, Delta)
      du_tau_dx = (-ut_p2 + 8.0_wp*ut_p1 - 8.0_wp*ut_m1 + ut_m2) / (12.0_wp * h)
   end subroutine sgs_source_at

end module mms
