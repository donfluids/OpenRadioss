module fields
   use kinds,      only : wp
   use constants,  only : NVAR, NPRIM
   use mesh_types, only : t_mesh
   implicit none
   private

   public :: t_state, alloc_state, free_state

   type :: t_state
      integer :: nc = 0

      ! Conservative state and residual.
      real(wp), allocatable :: U (:,:)        ! (NVAR, nc) conservative
      real(wp), allocatable :: R (:,:)        ! (NVAR, nc) residual
      real(wp), allocatable :: U0(:,:)        ! (NVAR, nc) RK stage stash

      ! M3 additions: primitives, primitive gradients, limiters (Venkatakrishnan).
      real(wp), allocatable :: W    (:,:)     ! (NPRIM, nc) primitives [rho, u, v, w, p]
      real(wp), allocatable :: gradW(:,:,:)   ! (3, NPRIM, nc) ∇W per cell
      real(wp), allocatable :: psi  (:,:)     ! (NPRIM, nc) per-variable limiter ψ ∈ [0,1]
   end type t_state

contains

   subroutine alloc_state(s, mesh)
      type(t_state), intent(out) :: s
      type(t_mesh),  intent(in)  :: mesh
      s%nc = mesh%nc_total
      allocate(s%U    (NVAR,        s%nc), source=0.0_wp)
      allocate(s%R    (NVAR,        s%nc), source=0.0_wp)
      allocate(s%U0   (NVAR,        s%nc), source=0.0_wp)
      allocate(s%W    (NPRIM,       s%nc), source=0.0_wp)
      allocate(s%gradW(3,    NPRIM, s%nc), source=0.0_wp)
      allocate(s%psi  (NPRIM,       s%nc), source=1.0_wp)
   end subroutine alloc_state

   subroutine free_state(s)
      type(t_state), intent(inout) :: s
      if (allocated(s%U    )) deallocate(s%U    )
      if (allocated(s%R    )) deallocate(s%R    )
      if (allocated(s%U0   )) deallocate(s%U0   )
      if (allocated(s%W    )) deallocate(s%W    )
      if (allocated(s%gradW)) deallocate(s%gradW)
      if (allocated(s%psi  )) deallocate(s%psi  )
      s%nc = 0
   end subroutine free_state

end module fields
