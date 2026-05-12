module mesh_module
   use kinds,         only : wp
   use mesh_types,    only : t_mesh
   use mesh_io_gmsh,  only : gmsh_read, t_bface_raw
   use mesh_topology, only : build_topology
   use mesh_metrics,  only : build_metrics
   implicit none
   private

   public :: load_mesh

contains

   subroutine load_mesh(filename, mesh)
      character(len=*), intent(in)  :: filename
      type(t_mesh),     intent(out) :: mesh

      type(t_bface_raw), allocatable :: bfaces(:)
      integer :: nb

      call gmsh_read(filename, mesh, bfaces, nb)
      call build_topology(mesh, bfaces, nb)
      call build_metrics(mesh)

      if (allocated(bfaces)) deallocate(bfaces)
   end subroutine load_mesh

end module mesh_module
