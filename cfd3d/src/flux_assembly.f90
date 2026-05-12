module flux_assembly
   use kinds,         only : wp
   use constants,     only : NVAR
   use mesh_types,    only : t_mesh, t_patch
   use bc_types,      only : t_bc_data
   use bc_apply,      only : ghost_state
   use riemann_hllc,  only : hllc_flux
   use fields,        only : t_state
   implicit none
   private

   public :: compute_residual

contains

   subroutine compute_residual(mesh, bc_dat, s)
      type(t_mesh),    intent(in)    :: mesh
      type(t_bc_data), intent(in)    :: bc_dat(:)         ! (np)
      type(t_state),   intent(inout) :: s

      integer  :: ifc, c_o, c_n, ip
      real(wp) :: Fflx(NVAR), QL(NVAR), QR(NVAR), nrml(3), area

      s%R = 0.0_wp

      ! Interior faces
      do ifc = 1, mesh%nf_interior
         c_o  = mesh%face_owner(ifc)
         c_n  = mesh%face_neighbor(ifc)
         nrml = mesh%face_normal(:, ifc)
         area = mesh%face_area(ifc)
         QL   = s%U(:, c_o)
         QR   = s%U(:, c_n)
         Fflx = hllc_flux(QL, QR, nrml)
         s%R(:, c_o) = s%R(:, c_o) + Fflx * area
         s%R(:, c_n) = s%R(:, c_n) - Fflx * area
      end do

      ! Boundary faces — ghost state, no scatter to neighbor
      do ifc = mesh%nf_interior + 1, mesh%nf
         c_o  = mesh%face_owner(ifc)
         ip   = mesh%face_patch(ifc)
         nrml = mesh%face_normal(:, ifc)
         area = mesh%face_area(ifc)
         QL   = s%U(:, c_o)
         call ghost_state(mesh%patches(ip)%bc_type, bc_dat(ip), QL, nrml, QR)
         Fflx = hllc_flux(QL, QR, nrml)
         s%R(:, c_o) = s%R(:, c_o) + Fflx * area
      end do
   end subroutine compute_residual

end module flux_assembly
