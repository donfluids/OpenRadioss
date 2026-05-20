module mesh_types
   use kinds, only : wp, i4
   use mpi_f08, only : MPI_Request
   implicit none
   private

   public :: t_patch, t_mesh, t_halo
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

   ! Halo descriptor for inter-rank cell-state exchange.
   ! Maintained per-mesh; built once after partition, reused per RK stage.
   type :: t_halo
      integer :: n_neighbors = 0
      integer, allocatable :: neighbor_ranks(:)        ! (n_neighbors)

      ! CSR layout — send side: for each neighbor, list of LOCAL cell indices to pack.
      integer, allocatable :: send_offset(:)           ! (n_neighbors+1)
      integer, allocatable :: send_cells (:)           ! (send_offset(n_neighbors+1)-1)

      ! CSR layout — recv side: for each neighbor, list of GHOST cell indices to unpack into.
      integer, allocatable :: recv_offset(:)           ! (n_neighbors+1)
      integer, allocatable :: recv_ghosts(:)           ! (recv_offset(n_neighbors+1)-1)

      ! Channel U: NVAR floats per cell (conservative state).
      real(wp), allocatable :: send_buf(:)             ! (NVAR * total_send)
      real(wp), allocatable :: recv_buf(:)             ! (NVAR * total_recv)
      type(MPI_Request), allocatable :: send_req(:)
      type(MPI_Request), allocatable :: recv_req(:)
      logical :: persistent_inited = .false.

      ! Channel GP: gradients + limiters (3*NPRIM + NPRIM = 20 floats per cell).
      real(wp), allocatable :: send_buf_gp(:)
      real(wp), allocatable :: recv_buf_gp(:)
      type(MPI_Request), allocatable :: send_req_gp(:)
      type(MPI_Request), allocatable :: recv_req_gp(:)
      logical :: persistent_inited_gp = .false.
   end type t_halo

   type :: t_mesh
      ! Vertices
      integer  :: nv = 0
      real(wp), allocatable :: xv(:,:)                 ! (3, nv)

      ! Cells (SoA, CSR vertex list). Ordering:
      !   [1..nc_internal] are locally-owned cells
      !   [nc_internal+1..nc_total] are halo/ghost cells mirroring neighbor-rank state
      integer  :: nc_internal = 0
      integer  :: nc_total    = 0
      integer,  allocatable :: cell_type(:)            ! (nc_total) CELL_*
      integer,  allocatable :: cell_vtx_ptr(:)         ! (nc_total+1) CSR ptr
      integer,  allocatable :: cell_vtx(:)             ! flat CSR data
      real(wp), allocatable :: cell_volume(:)          ! (nc_total)
      real(wp), allocatable :: cell_centroid(:,:)      ! (3, nc_total)
      integer,  allocatable :: cell_perm(:)            ! (nc_total) identity in M1; RCM in M2

      ! Faces, ordered:
      !   [1..nf_pure_interior]                       — both cells local
      !   [nf_pure_interior+1..nf_interior]           — partition faces (1 local + 1 ghost)
      !   [nf_interior+1..nf]                         — physical boundary faces, grouped by patch
      integer  :: nf = 0
      integer  :: nf_pure_interior = 0
      integer  :: nf_interior      = 0                 ! = nf_pure_interior + n_partition_faces
      integer  :: nf_boundary      = 0
      integer,  allocatable :: face_owner(:)           ! (nf) local cell index
      integer,  allocatable :: face_owner_lf(:)        ! (nf) local face index on owner
      integer,  allocatable :: face_neighbor(:)        ! (nf) 0 for boundary
      real(wp), allocatable :: face_area(:)            ! (nf)
      real(wp), allocatable :: face_normal(:,:)        ! (3, nf) unit, owner->neighbor
      real(wp), allocatable :: face_centroid(:,:)      ! (3, nf)
      integer,  allocatable :: face_patch(:)           ! (nf) patch index or 0 for interior

      ! Patches
      integer  :: np = 0
      type(t_patch), allocatable :: patches(:)         ! (np)

      ! Halo exchange (only meaningful when nproc > 1)
      type(t_halo) :: halo

      ! Global cell IDs and owning rank — for each cell in [1..nc_total].
      ! Local cells [1..nc_internal] have owner_rank = my_rank.
      ! Ghost cells [nc_internal+1..nc_total] have owner_rank = neighbor rank.
      integer, allocatable :: cell_global_id (:)       ! (nc_total)
      integer, allocatable :: cell_owner_rank(:)       ! (nc_total)

      ! Cut-cell tables (Phase 1b). Populated by cut_cell::build_cut_cell_tables
      ! after build_metrics. For meshes with no obstacle cubes these are
      ! initialised to copies of cell_volume / face_area, and n_frags = 0.
      ! For meshes with obstacle cubes, per-cell V_eff < V_full where
      ! the cube clips the cell, per-face A_eff < A_full where the cube
      ! covers the face, and fragments carry the obstacle-surface
      ! geometry (outward normal, area, centroid) inside each cut cell.
      real(wp), allocatable :: cell_vol_eff (:)        ! (nc_total)
      real(wp), allocatable :: face_area_eff(:)        ! (nf)
      integer  :: n_frags = 0
      integer,  allocatable :: cell_frag_offset(:)     ! (nc_total+1) CSR
      real(wp), allocatable :: frag_normal  (:,:)      ! (3, n_frags) outward from cell
      real(wp), allocatable :: frag_area    (:)        ! (n_frags)
      real(wp), allocatable :: frag_centroid(:,:)      ! (3, n_frags)

      ! Cut-cell merging (Phase 1d). Sliver cells (cell_vol_eff/cell_volume
      ! below a threshold) are linked to a larger non-sliver face-neighbor
      ! ("host"). The group shares a single conserved state advanced with
      ! the combined group volume, which restores a sane CFL time step.
      !   cell_merge_root(c) = representative cell of c's group
      !                        (= c for un-merged cells and for hosts)
      !   cell_merge_vol (c) = total V_eff of c's group (= V_eff for
      !                        un-merged cells)
      integer,  allocatable :: cell_merge_root(:)      ! (nc_total)
      real(wp), allocatable :: cell_merge_vol (:)      ! (nc_total)
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
