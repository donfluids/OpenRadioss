! Gas thermophysical properties — set once at startup from the namelist
! and used by EOS and viscous-flux modules. Module-level state keeps the
! per-face hot path branch-free.
module gas_properties
   use kinds,     only : wp
   use constants, only : GAMMA, GM1
   implicit none
   private

   public :: gas, init_gas
   public :: SGS_NONE, SGS_SMAG, SGS_WALE

   integer, parameter :: SGS_NONE = 0
   integer, parameter :: SGS_SMAG = 1
   integer, parameter :: SGS_WALE = 2

   type :: t_gas
      real(wp) :: R_gas = 1.0_wp          ! specific gas constant (non-dim default)
      real(wp) :: cv    = 1.0_wp / GM1
      real(wp) :: cp    = GAMMA / GM1

      ! Viscosity: Sutherland's law μ(T) = μ_ref (T/T_ref)^1.5 (T_ref + S) / (T + S).
      ! If sutherland=.false., μ ≡ mu_const.
      logical  :: sutherland = .true.
      real(wp) :: mu_const   = 0.0_wp     ! when sutherland=.false.
      real(wp) :: mu_ref     = 1.716e-5_wp
      real(wp) :: T_ref      = 273.15_wp
      real(wp) :: S_S        = 110.4_wp

      real(wp) :: Pr         = 0.72_wp    ! Prandtl number (constant)

      ! Sub-grid scale (M4 algebraic LES) configuration.
      integer  :: sgs_kind   = SGS_NONE
      real(wp) :: C_s        = 0.17_wp    ! Smagorinsky constant
      real(wp) :: C_w        = 0.5_wp     ! WALE constant
      real(wp) :: Pr_t       = 0.9_wp     ! turbulent Prandtl
   end type t_gas

   ! Module-level singleton. Mutable, but written only by init_gas at startup.
   type(t_gas) :: gas

contains

   subroutine init_gas(R_gas, sutherland, mu_const, mu_ref, T_ref, S_S, Pr, &
                       sgs_kind, C_s, C_w, Pr_t)
      real(wp), intent(in) :: R_gas, mu_const, mu_ref, T_ref, S_S, Pr
      logical,  intent(in) :: sutherland
      integer,  intent(in), optional :: sgs_kind
      real(wp), intent(in), optional :: C_s, C_w, Pr_t
      gas%R_gas      = R_gas
      gas%cv         = R_gas / GM1
      gas%cp         = R_gas * GAMMA / GM1
      gas%sutherland = sutherland
      gas%mu_const   = mu_const
      gas%mu_ref     = mu_ref
      gas%T_ref      = T_ref
      gas%S_S        = S_S
      gas%Pr         = Pr
      if (present(sgs_kind)) gas%sgs_kind = sgs_kind
      if (present(C_s     )) gas%C_s      = C_s
      if (present(C_w     )) gas%C_w      = C_w
      if (present(Pr_t    )) gas%Pr_t     = Pr_t
   end subroutine init_gas

end module gas_properties
