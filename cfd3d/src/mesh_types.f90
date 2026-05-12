module mesh_types
   use kinds, only : wp, i4
   implicit none
   private

   public :: t_patch, t_mesh
   public :: CELL_TET, CELL_HEX, CELL_PRISM, CELL_PYRAMID
   public :: GMSH_TYPE_TRI, GMSH_TYPE_QUAD
   public :: GMSH_TYPE_TET, GMSH_TYPE_HEX, GMSH_TYPE_PRISM, GMSH_TYPE_PYRAMID
   public :: n_vtx_per_cell, n_faces_per_cell, max_vpf

   ! Internal cell-type tags (we re-tag Gmsh element ids onto these)
   integer, parameter :: CELL_TET     = 1
   integer, parameter :: CELL_HEX     = 2
   integer, parameter :: CELL_PRISM   = 3
   integer, parameter :: CELL_PYRAMID = 4

   ! Gmsh element type ids (msh2 spec)
   integer, parameter :: GMSH_TYPE_TRI     = 2
   integer, parameter :: GMSH_TYPE_QUAD    = 3
   integer, parameter :: GMSH_TYPE_TET     = 4
   integer, parameter :: GMSH_TYPE_HEX     = 5
   integer, parameter :: GMSH_TYPE_PRISM   = 6
   integer, parameter :: GMSH_TYPE_PYRAMID = 7

   ! Max vertices per face across supported cell types (hex/prism/pyramid quad = 4)
   integer, parameter :: max_vpf = 4

   type :: t_patch
      character(len=64) :: name = ''
      integer           :: bc_type = 0       ! one of constants::BC_*
      integer           :: gmsh_tag = -1     ! physical group tag for matching
      integer           :: face_start = 0    ! first boundary face index
      integer           :: face_count = 0    ! number of boundary faces
   end type t_patch

   type :: t_mesh
      ! Vertices
      integer  :: nv = 0
      real(wp), allocatable :: xv(:,:)                 ! (3, nv)

      ! Cells (SoA, CSR vertex list)
      integer  :: nc_internal = 0                     ! cells owned locally (M1: all)
      integer  :: nc_total    = 0                     ! nc_internal + nc_ghost (M1: equal)
      integer,  allocatable :: cell_type(:)            ! (nc_total) CELL_*
      integer,  allocatable :: cell_vtx_ptr(:)         ! (nc_total+1) CSR ptr
      integer,  allocatable :: cell_vtx(:)             ! flat CSR data
      real(wp), allocatable :: cell_volume(:)          ! (nc_total)
      real(wp), allocatable :: cell_centroid(:,:)      ! (3, nc_total)
      integer,  allocatable :: cell_perm(:)            ! (nc_total) identity in M1; RCM in M2

      ! Faces — interior first, then boundary grouped by patch
      integer  :: nf = 0
      integer  :: nf_interior = 0
      integer  :: nf_boundary = 0
      integer,  allocatable :: face_owner(:)           ! (nf) local cell index
      integer,  allocatable :: face_owner_lf(:)        ! (nf) local face index on owner
      integer,  allocatable :: face_neighbor(:)        ! (nf) 0 for boundary in M1
      real(wp), allocatable :: face_area(:)            ! (nf)
      real(wp), allocatable :: face_normal(:,:)        ! (3, nf) unit, owner->neighbor
      real(wp), allocatable :: face_centroid(:,:)      ! (3, nf)
      integer,  allocatable :: face_patch(:)           ! (nf) patch index or 0 for interior

      ! Patches
      integer  :: np = 0
      type(t_patch), allocatable :: patches(:)         ! (np)
   end type t_mesh

contains

   pure function n_vtx_per_cell(ct) result(n)
      integer, intent(in) :: ct
      integer :: n
      select case (ct)
      case (CELL_TET)
         n = 4
      case (CELL_HEX)
         n = 8
      case (CELL_PRISM)
         n = 6
      case (CELL_PYRAMID)
         n = 5
      case default
         n = 0
      end select
   end function n_vtx_per_cell

   pure function n_faces_per_cell(ct) result(n)
      integer, intent(in) :: ct
      integer :: n
      select case (ct)
      case (CELL_TET)
         n = 4
      case (CELL_HEX)
         n = 6
      case (CELL_PRISM)
         n = 5
      case (CELL_PYRAMID)
         n = 5
      case default
         n = 0
      end select
   end function n_faces_per_cell

end module mesh_types
