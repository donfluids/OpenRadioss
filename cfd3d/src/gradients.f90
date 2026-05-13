! Cell-centered least-squares gradient computation for primitive variables.
! Weighted by 1/|dr|^2. Includes face-neighbor cells (interior, ghost) AND
! BC ghost states at physical-boundary faces (half-step stencil).
module gradients
   use kinds,         only : wp
   use constants,     only : NVAR, NPRIM, IP_RHO, IP_U, IP_V, IP_W, IP_P, &
                             BC_SLIP_WALL, BC_SYMMETRY, BC_SUPERSONIC_OUTLET, &
                             BC_FARFIELD
   use mesh_types,    only : t_mesh
   use fields,        only : t_state
   use bc_types,      only : t_bc_data
   use bc_apply,      only : ghost_state
   use eos_ideal_gas, only : prim_from_cons
   implicit none
   private

   public :: compute_primitives, compute_gradients

contains

   ! Skip boundary types that don't carry meaningful gradient information
   ! (slip walls reflect the normal velocity, producing spurious LSQ
   ! gradients at edge cells; outflow boundaries simply copy state).
   pure function bc_contributes_to_gradient(bc_type) result(yes)
      integer, intent(in) :: bc_type
      logical :: yes
      select case (bc_type)
      case (BC_SLIP_WALL, BC_SYMMETRY, BC_SUPERSONIC_OUTLET, BC_FARFIELD)
         yes = .false.
      case default
         yes = .true.
      end select
   end function bc_contributes_to_gradient

   ! Fill s%W (primitives) from s%U (conservative) for ALL cells incl. ghosts.
   subroutine compute_primitives(mesh, s)
      type(t_mesh),  intent(in)    :: mesh
      type(t_state), intent(inout) :: s
      integer  :: c
      real(wp) :: rho, u, v, w, p
      do c = 1, mesh%nc_total
         call prim_from_cons(s%U(:, c), rho, u, v, w, p)
         s%W(IP_RHO, c) = rho
         s%W(IP_U,   c) = u
         s%W(IP_V,   c) = v
         s%W(IP_W,   c) = w
         s%W(IP_P,   c) = p
      end do
   end subroutine compute_primitives

   subroutine compute_gradients(mesh, bc_dat, s)
      type(t_mesh),    intent(in)    :: mesh
      type(t_bc_data), intent(in)    :: bc_dat(:)
      type(t_state),   intent(inout) :: s

      integer  :: c, f, c_o, c_n, ip, vv
      real(wp), allocatable :: A_cells(:,:,:)
      real(wp), allocatable :: b_cells(:,:,:)
      real(wp) :: A(3,3), Ainv(3,3), det
      real(wp) :: dr(3), w2, dW(NPRIM), Wbc(NPRIM)

      allocate(A_cells(3,    3, mesh%nc_internal), source=0.0_wp)
      allocate(b_cells(3, NPRIM, mesh%nc_internal), source=0.0_wp)

      do f = 1, mesh%nf
         c_o = mesh%face_owner(f)
         c_n = mesh%face_neighbor(f)
         if (c_n > 0) then
            dr = mesh%cell_centroid(:, c_n) - mesh%cell_centroid(:, c_o)
            w2 = 1.0_wp / max(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3), 1.0e-30_wp)
            dW = s%W(:, c_n) - s%W(:, c_o)
            if (c_o <= mesh%nc_internal) &
               call accumulate(A_cells(:,:,c_o), b_cells(:,:,c_o), dr, w2, dW)
            if (c_n <= mesh%nc_internal) &
               call accumulate(A_cells(:,:,c_n), b_cells(:,:,c_n), -dr, w2, -dW)
         else
            ! Physical boundary face. Only include in the LSQ stencil for BC
            ! types that carry meaningful information (Dirichlet, no-slip,
            ! supersonic inlet). Slip walls / symmetry / outflow are excluded
            ! to avoid the reflection-induced spurious normal gradient.
            ip = mesh%face_patch(f)
            if (.not. bc_contributes_to_gradient(mesh%patches(ip)%bc_type)) cycle
            block
               real(wp) :: UR(NVAR), QL(NVAR), rho_b, u_b, v_b, w_b, p_b
               QL = s%U(:, c_o)
               call ghost_state(mesh%patches(ip)%bc_type, bc_dat(ip), &
                                QL, mesh%face_normal(:, f), UR)
               call prim_from_cons(UR, rho_b, u_b, v_b, w_b, p_b)
               Wbc(IP_RHO) = rho_b
               Wbc(IP_U)   = u_b
               Wbc(IP_V)   = v_b
               Wbc(IP_W)   = w_b
               Wbc(IP_P)   = p_b
            end block
            dr = mesh%face_centroid(:, f) - mesh%cell_centroid(:, c_o)
            w2 = 1.0_wp / max(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3), 1.0e-30_wp)
            dW = 0.5_wp * (Wbc - s%W(:, c_o))
            if (c_o <= mesh%nc_internal) &
               call accumulate(A_cells(:,:,c_o), b_cells(:,:,c_o), dr, w2, dW)
         end if
      end do

      do c = 1, mesh%nc_internal
         A = A_cells(:,:,c)
         call invert_3x3(A, Ainv, det)
         if (abs(det) < 1.0e-30_wp) then
            s%gradW(:,:,c) = 0.0_wp
         else
            do vv = 1, NPRIM
               s%gradW(:, vv, c) = matmul(Ainv, b_cells(:, vv, c))
            end do
         end if
      end do

      deallocate(A_cells, b_cells)
   end subroutine compute_gradients

   subroutine accumulate(A, b, dr, w2, dW)
      real(wp), intent(inout) :: A(3,3), b(3, NPRIM)
      real(wp), intent(in)    :: dr(3), w2, dW(NPRIM)
      integer  :: i, j, v
      do i = 1, 3
         do j = 1, 3
            A(i,j) = A(i,j) + w2 * dr(i) * dr(j)
         end do
      end do
      do v = 1, NPRIM
         do i = 1, 3
            b(i, v) = b(i, v) + w2 * dr(i) * dW(v)
         end do
      end do
   end subroutine accumulate

   pure subroutine invert_3x3(A, Ainv, det)
      real(wp), intent(in)  :: A(3,3)
      real(wp), intent(out) :: Ainv(3,3), det
      det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
          - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
          + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
      if (abs(det) < 1.0e-30_wp) then
         Ainv = 0.0_wp
         return
      end if
      Ainv(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2)) / det
      Ainv(1,2) = (A(1,3)*A(3,2) - A(1,2)*A(3,3)) / det
      Ainv(1,3) = (A(1,2)*A(2,3) - A(1,3)*A(2,2)) / det
      Ainv(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3)) / det
      Ainv(2,2) = (A(1,1)*A(3,3) - A(1,3)*A(3,1)) / det
      Ainv(2,3) = (A(1,3)*A(2,1) - A(1,1)*A(2,3)) / det
      Ainv(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1)) / det
      Ainv(3,2) = (A(1,2)*A(3,1) - A(1,1)*A(3,2)) / det
      Ainv(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1)) / det
   end subroutine invert_3x3

end module gradients
