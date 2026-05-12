module fields
   use kinds,     only : wp
   use constants, only : NVAR
   use mesh_types, only : t_mesh
   implicit none
   private

   public :: t_state, alloc_state, free_state

   type :: t_state
      integer :: nc = 0
      real(wp), allocatable :: U(:,:)        ! (NVAR, nc) conservative
      real(wp), allocatable :: R(:,:)        ! (NVAR, nc) residual
      real(wp), allocatable :: U0(:,:)       ! (NVAR, nc) RK stage buffer
   end type t_state

contains

   subroutine alloc_state(s, mesh)
      type(t_state), intent(out) :: s
      type(t_mesh),  intent(in)  :: mesh
      s%nc = mesh%nc_total
      allocate(s%U (NVAR, s%nc), source=0.0_wp)
      allocate(s%R (NVAR, s%nc), source=0.0_wp)
      allocate(s%U0(NVAR, s%nc), source=0.0_wp)
   end subroutine alloc_state

   subroutine free_state(s)
      type(t_state), intent(inout) :: s
      if (allocated(s%U )) deallocate(s%U )
      if (allocated(s%R )) deallocate(s%R )
      if (allocated(s%U0)) deallocate(s%U0)
      s%nc = 0
   end subroutine free_state

end module fields
