module bc_types
   use kinds,     only : wp
   use constants, only : NVAR
   implicit none
   private

   public :: t_bc_data, set_default_freestream

   ! Per-patch BC parameters. Far-field / supersonic-inlet primitives stored here.
   type :: t_bc_data
      real(wp) :: rho = 1.0_wp
      real(wp) :: u   = 0.0_wp
      real(wp) :: v   = 0.0_wp
      real(wp) :: w   = 0.0_wp
      real(wp) :: p   = 1.0_wp
   end type t_bc_data

contains

   pure subroutine set_default_freestream(b)
      type(t_bc_data), intent(out) :: b
      b%rho = 1.0_wp
      b%u   = 0.0_wp
      b%v   = 0.0_wp
      b%w   = 0.0_wp
      b%p   = 1.0_wp
   end subroutine set_default_freestream

end module bc_types
