! Jones-Wilkins-Lee (JWL) equation of state for detonation products.
!
! The standard JWL form for the pressure of explosive products is
!
!     p(ρ, e) = A·(1 − ω/(R₁·V))·exp(−R₁·V)
!             + B·(1 − ω/(R₂·V))·exp(−R₂·V)
!             + ω·ρ·e
!
! where V = ρ₀/ρ is the relative volume, ρ₀ is the un-reacted explosive
! density, and e is the specific internal energy (J/kg). Constants for TNT
! from the Dobratz-Crawford LLNL Explosives Handbook (1985):
!
!     A   = 371.2 GPa
!     B   = 3.231 GPa
!     R₁  = 4.15
!     R₂  = 0.95
!     ω   = 0.30
!     ρ₀  = 1630 kg/m³
!     E₀  = 7.0 GJ/m³  (energy per unit reference volume)
!
! At the reference state (V = 1, ρ = ρ₀), the JWL pressure of the
! detonation products is the immediate post-detonation value used as an
! initial condition for the surrounding fluid simulation. For TNT this
! evaluates to roughly 8.4 GPa.
!
! M5.2 uses this module only to *initialize* the explosive-products
! region; evolution is single-γ ideal gas. The functions below are
! general (any ρ, any e ≥ 0) so they can be reused for full multi-gas
! evolution in a later sub-phase.
module eos_jwl
   use kinds, only : wp
   implicit none
   private

   public :: jwl_tnt
   public :: jwl_pressure, jwl_sound_speed_squared
   public :: jwl_state_at_reference

   type :: t_jwl_params
      real(wp) :: A      ! Pa
      real(wp) :: B      ! Pa
      real(wp) :: R1     ! dimensionless
      real(wp) :: R2     ! dimensionless
      real(wp) :: omega  ! dimensionless (Grüneisen)
      real(wp) :: rho0   ! kg/m³  un-reacted explosive density
      real(wp) :: E0     ! J/m³   energy per unit reference volume
   end type t_jwl_params

   ! Dobratz-Crawford 1985 TNT parameters.
   type(t_jwl_params), parameter :: jwl_tnt = t_jwl_params( &
      A     = 371.2e9_wp, &
      B     = 3.231e9_wp, &
      R1    = 4.15_wp,    &
      R2    = 0.95_wp,    &
      omega = 0.30_wp,    &
      rho0  = 1630.0_wp,  &
      E0    = 7.0e9_wp )

   real(wp), parameter :: V_MIN = 1.0e-3_wp   ! relative-volume floor

contains

   pure function jwl_pressure(jp, rho, e) result(p)
      type(t_jwl_params), intent(in) :: jp
      real(wp),           intent(in) :: rho, e
      real(wp) :: p, V

      V = max(jp%rho0 / max(rho, 1.0e-30_wp), V_MIN)
      p = jp%A * (1.0_wp - jp%omega / (jp%R1 * V)) * exp(-jp%R1 * V) &
        + jp%B * (1.0_wp - jp%omega / (jp%R2 * V)) * exp(-jp%R2 * V) &
        + jp%omega * rho * e
   end function jwl_pressure

   ! ∂p/∂ρ|_e at constant specific energy. For the full thermodynamic
   ! sound speed (constant entropy), include the (∂p/∂e)_ρ · (p/ρ²) term
   ! from the Grüneisen relation.
   !
   !   c² = (∂p/∂ρ)_e  +  (∂p/∂e)_ρ · p / ρ²
   !      = ∂_ρ[A·f₁(V)·exp(−R₁V) + B·f₂(V)·exp(−R₂V)]  +  ω·e  +  ω·p/ρ
   !
   ! with V = ρ₀/ρ, ∂V/∂ρ = −ρ₀/ρ², ∂_ρ[g(V)] = g'(V)·∂V/∂ρ.
   pure function jwl_sound_speed_squared(jp, rho, e) result(c2)
      type(t_jwl_params), intent(in) :: jp
      real(wp),           intent(in) :: rho, e
      real(wp) :: c2, V, p, dpdrho
      real(wp) :: f1, df1dV, term1, f2, df2dV, term2

      V = max(jp%rho0 / max(rho, 1.0e-30_wp), V_MIN)

      f1     = 1.0_wp - jp%omega / (jp%R1 * V)
      df1dV  = jp%omega / (jp%R1 * V * V)
      term1  = jp%A * exp(-jp%R1 * V) * (df1dV - jp%R1 * f1)

      f2     = 1.0_wp - jp%omega / (jp%R2 * V)
      df2dV  = jp%omega / (jp%R2 * V * V)
      term2  = jp%B * exp(-jp%R2 * V) * (df2dV - jp%R2 * f2)

      ! d/drho of cold-curve = (term1 + term2) · dV/drho.
      dpdrho = -(term1 + term2) * jp%rho0 / (max(rho, 1.0e-30_wp)**2)

      p  = jwl_pressure(jp, rho, e)
      c2 = dpdrho + jp%omega * e + jp%omega * p / max(rho, 1.0e-30_wp)
      c2 = max(c2, 1.0e-10_wp)
   end function jwl_sound_speed_squared

   ! Reference initial state of detonation products at the un-expanded
   ! charge volume: ρ = ρ₀, e = E₀/ρ₀ (specific internal energy from the
   ! volumetric reference E₀). Returns the JWL pressure at that state.
   subroutine jwl_state_at_reference(jp, rho_out, e_out, p_out)
      type(t_jwl_params), intent(in)  :: jp
      real(wp),           intent(out) :: rho_out, e_out, p_out
      rho_out = jp%rho0
      e_out   = jp%E0 / jp%rho0
      p_out   = jwl_pressure(jp, rho_out, e_out)
   end subroutine jwl_state_at_reference

end module eos_jwl
