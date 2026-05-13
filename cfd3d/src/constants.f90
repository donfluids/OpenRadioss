module constants
   use kinds, only : wp
   implicit none
   private

   public :: PI, GAMMA, GM1, GP1, ONE_OVER_GM1
   public :: SMALL, TINY_RHO, TINY_P
   public :: BC_INTERIOR, BC_SLIP_WALL, BC_SYMMETRY, &
             BC_SUPERSONIC_INLET, BC_SUPERSONIC_OUTLET, BC_FARFIELD, &
             BC_NO_SLIP_WALL, BC_DIRICHLET
   public :: NVAR, IRHO, IRHOU, IRHOV, IRHOW, IRHOE
   public :: NPRIM, IP_RHO, IP_U, IP_V, IP_W, IP_P

   real(wp), parameter :: PI = 3.141592653589793238462643383279502884_wp

   ! Ratio of specific heats (ideal gas, air). Override via solver_control.
   real(wp), parameter :: GAMMA = 1.4_wp
   real(wp), parameter :: GM1   = GAMMA - 1.0_wp
   real(wp), parameter :: GP1   = GAMMA + 1.0_wp
   real(wp), parameter :: ONE_OVER_GM1 = 1.0_wp / GM1

   real(wp), parameter :: SMALL    = 1.0e-30_wp
   real(wp), parameter :: TINY_RHO = 1.0e-12_wp
   real(wp), parameter :: TINY_P   = 1.0e-12_wp

   ! Boundary condition tags. Stored on patches; dispatched in bc_apply.
   integer, parameter :: BC_INTERIOR          = 0
   integer, parameter :: BC_SLIP_WALL         = 1
   integer, parameter :: BC_SYMMETRY          = 2
   integer, parameter :: BC_SUPERSONIC_INLET  = 3
   integer, parameter :: BC_SUPERSONIC_OUTLET = 4
   integer, parameter :: BC_FARFIELD          = 5
   integer, parameter :: BC_NO_SLIP_WALL      = 6   ! viscous wall: u=v=w=0, dT/dn=0 (adiabatic)
   integer, parameter :: BC_DIRICHLET         = 7   ! all primitives prescribed (MMS, IO tests)

   ! Conservative state index layout: U = [rho, rho*u, rho*v, rho*w, rho*E]
   integer, parameter :: NVAR  = 5
   integer, parameter :: IRHO  = 1
   integer, parameter :: IRHOU = 2
   integer, parameter :: IRHOV = 3
   integer, parameter :: IRHOW = 4
   integer, parameter :: IRHOE = 5

   ! Primitive state index layout: W = [rho, u, v, w, p]
   integer, parameter :: NPRIM  = 5
   integer, parameter :: IP_RHO = 1
   integer, parameter :: IP_U   = 2
   integer, parameter :: IP_V   = 3
   integer, parameter :: IP_W   = 4
   integer, parameter :: IP_P   = 5
end module constants
