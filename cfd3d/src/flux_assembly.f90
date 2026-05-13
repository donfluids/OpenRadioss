module flux_assembly
   use kinds,         only : wp
   use constants,     only : NVAR
   use mesh_types,    only : t_mesh, t_patch
   use bc_types,      only : t_bc_data
   use bc_apply,      only : ghost_state
   use riemann_hllc,  only : hllc_flux
   use fields,        only : t_state
   use mpi_runtime,   only : t_mpi_ctx
   use halo_exchange, only : halo_pack_and_start, halo_wait
   implicit none
   private

   public :: compute_residual            ! MPI-aware orchestrator
   public :: residual_begin
   public :: residual_pure_interior
   public :: residual_partition
   public :: residual_boundary

contains

   ! MPI-aware residual assembly with comm/compute overlap.
   subroutine compute_residual(mesh, bc_dat, s, ctx)
      type(t_mesh),    intent(inout) :: mesh
      type(t_bc_data), intent(in)    :: bc_dat(:)
      type(t_state),   intent(inout) :: s
      type(t_mpi_ctx), intent(in)    :: ctx

      call residual_begin(s)
      call halo_pack_and_start(mesh, s%U, ctx)
      call residual_pure_interior(mesh, s)
      call halo_wait(mesh, s%U)
      call residual_partition(mesh, s)
      call residual_boundary(mesh, bc_dat, s)
   end subroutine compute_residual

   subroutine residual_begin(s)
      type(t_state), intent(inout) :: s
      s%R = 0.0_wp
   end subroutine residual_begin

   ! Faces with both cells local; safe to compute while halos are in flight.
   subroutine residual_pure_interior(mesh, s)
      type(t_mesh),  intent(in)    :: mesh
      type(t_state), intent(inout) :: s
      integer  :: ifc, c_o, c_n
      real(wp) :: Fflx(NVAR), nrml(3), area
      do ifc = 1, mesh%nf_pure_interior
         c_o  = mesh%face_owner(ifc)
         c_n  = mesh%face_neighbor(ifc)
         nrml = mesh%face_normal(:, ifc)
         area = mesh%face_area(ifc)
         Fflx = hllc_flux(s%U(:, c_o), s%U(:, c_n), nrml)
         s%R(:, c_o) = s%R(:, c_o) + Fflx * area
         s%R(:, c_n) = s%R(:, c_n) - Fflx * area
      end do
   end subroutine residual_pure_interior

   ! Partition faces: one local, one ghost. Requires halo data (call after halo_wait).
   subroutine residual_partition(mesh, s)
      type(t_mesh),  intent(in)    :: mesh
      type(t_state), intent(inout) :: s
      integer  :: ifc, c_o, c_n
      real(wp) :: Fflx(NVAR), nrml(3), area
      do ifc = mesh%nf_pure_interior + 1, mesh%nf_interior
         c_o  = mesh%face_owner(ifc)
         c_n  = mesh%face_neighbor(ifc)
         nrml = mesh%face_normal(:, ifc)
         area = mesh%face_area(ifc)
         Fflx = hllc_flux(s%U(:, c_o), s%U(:, c_n), nrml)
         ! Scatter to both — local cell uses it for its residual; ghost-cell residual
         ! is harmlessly written but won't be used in the RK update (loop runs only
         ! over [1..nc_internal]).
         s%R(:, c_o) = s%R(:, c_o) + Fflx * area
         s%R(:, c_n) = s%R(:, c_n) - Fflx * area
      end do
   end subroutine residual_partition

   ! Physical boundary faces: ghost state from BC, then HLLC.
   subroutine residual_boundary(mesh, bc_dat, s)
      type(t_mesh),    intent(in)    :: mesh
      type(t_bc_data), intent(in)    :: bc_dat(:)
      type(t_state),   intent(inout) :: s
      integer  :: ifc, c_o, ip
      real(wp) :: Fflx(NVAR), QL(NVAR), QR(NVAR), nrml(3), area
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
   end subroutine residual_boundary

end module flux_assembly
