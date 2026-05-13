module mesh_topology
   use kinds,        only : wp
   use mesh_types,   only : t_mesh, t_patch, max_vpf, &
                            CELL_TET, CELL_HEX, CELL_PRISM, CELL_PYRAMID, &
                            n_faces_per_cell
   use mesh_io_gmsh, only : t_bface_raw
   implicit none
   private

   public :: build_topology

   type :: t_face_cand
      integer :: key(max_vpf) = 0   ! sorted vertex ids; trailing 0 if tri
      integer :: nvf  = 0
      integer :: cell = 0
      integer :: lface = 0
   end type t_face_cand

   ! Local face vertex templates (1-based local cell-vertex indices)
   integer, parameter :: TET_FACE_NV(4) = [3, 3, 3, 3]
   integer, parameter :: TET_FACE_LV(4,4) = reshape( [ &
      1, 2, 3, 0,  &
      1, 2, 4, 0,  &
      2, 3, 4, 0,  &
      1, 3, 4, 0 ], shape=[4,4])  ! columns are faces; rows are local vertex slots

   integer, parameter :: HEX_FACE_NV(6) = [4, 4, 4, 4, 4, 4]
   integer, parameter :: HEX_FACE_LV(4,6) = reshape( [ &
      1, 2, 3, 4,  &
      5, 6, 7, 8,  &
      1, 2, 6, 5,  &
      2, 3, 7, 6,  &
      3, 4, 8, 7,  &
      4, 1, 5, 8 ], shape=[4,6])

   integer, parameter :: PRISM_FACE_NV(5) = [3, 3, 4, 4, 4]
   integer, parameter :: PRISM_FACE_LV(4,5) = reshape( [ &
      1, 2, 3, 0,  &
      4, 5, 6, 0,  &
      1, 2, 5, 4,  &
      2, 3, 6, 5,  &
      3, 1, 4, 6 ], shape=[4,5])

   integer, parameter :: PYR_FACE_NV(5) = [4, 3, 3, 3, 3]
   integer, parameter :: PYR_FACE_LV(4,5) = reshape( [ &
      1, 2, 3, 4,  &
      1, 2, 5, 0,  &
      2, 3, 5, 0,  &
      3, 4, 5, 0,  &
      4, 1, 5, 0 ], shape=[4,5])

contains

   ! Get local vertex indices and count for face lf of a cell of type ct.
   pure subroutine cell_face_local(ct, lf, nvf, lv)
      integer, intent(in)  :: ct, lf
      integer, intent(out) :: nvf
      integer, intent(out) :: lv(max_vpf)
      lv = 0
      select case (ct)
      case (CELL_TET)
         nvf = TET_FACE_NV(lf)
         lv(1:nvf) = TET_FACE_LV(1:nvf, lf)
      case (CELL_HEX)
         nvf = HEX_FACE_NV(lf)
         lv(1:nvf) = HEX_FACE_LV(1:nvf, lf)
      case (CELL_PRISM)
         nvf = PRISM_FACE_NV(lf)
         lv(1:nvf) = PRISM_FACE_LV(1:nvf, lf)
      case (CELL_PYRAMID)
         nvf = PYR_FACE_NV(lf)
         lv(1:nvf) = PYR_FACE_LV(1:nvf, lf)
      case default
         nvf = 0
      end select
   end subroutine cell_face_local

   subroutine build_topology(mesh, bfaces, nb)
      type(t_mesh),      intent(inout) :: mesh
      type(t_bface_raw), intent(in)    :: bfaces(:)
      integer,           intent(in)    :: nb

      integer :: c, ct, nf_c, lf, k, ncand
      integer :: nvf
      integer :: lv(max_vpf), gv(max_vpf)
      integer :: vptr
      type(t_face_cand), allocatable :: cand(:)
      type(t_face_cand), allocatable :: bcand(:)
      integer :: nf_interior, nf_boundary, nf_total
      integer :: i, j

      ! Phase 1: enumerate candidate faces
      ncand = 0
      do c = 1, mesh%nc_total
         nf_c = n_faces_per_cell(mesh%cell_type(c))
         ncand = ncand + nf_c
      end do
      allocate(cand(ncand))

      k = 0
      do c = 1, mesh%nc_total
         ct = mesh%cell_type(c)
         nf_c = n_faces_per_cell(ct)
         vptr = mesh%cell_vtx_ptr(c) - 1   ! base
         do lf = 1, nf_c
            call cell_face_local(ct, lf, nvf, lv)
            do i = 1, nvf
               gv(i) = mesh%cell_vtx(vptr + lv(i))
            end do
            do i = nvf+1, max_vpf
               gv(i) = 0
            end do
            call sort_key(gv, nvf)
            k = k + 1
            cand(k)%key   = gv
            cand(k)%nvf   = nvf
            cand(k)%cell  = c
            cand(k)%lface = lf
         end do
      end do

      ! Phase 2: sort candidates lexicographically by key
      call qsort_cands(cand, 1, ncand)

      ! Phase 3: walk through, pair consecutive equal keys → interior; lone → boundary
      ! Pre-scan to count
      nf_interior = 0
      nf_boundary = 0
      i = 1
      do while (i <= ncand)
         if (i < ncand) then
            if (cmp_key(cand(i)%key, cand(i+1)%key) == 0) then
               nf_interior = nf_interior + 1
               i = i + 2
               cycle
            end if
         end if
         nf_boundary = nf_boundary + 1
         i = i + 1
      end do

      nf_total = nf_interior + nf_boundary
      mesh%nf               = nf_total
      mesh%nf_interior      = nf_interior
      mesh%nf_pure_interior = nf_interior   ! serial / pre-partition: no partition faces yet
      mesh%nf_boundary      = nf_boundary

      allocate(mesh%face_owner   (nf_total))
      allocate(mesh%face_owner_lf(nf_total))
      allocate(mesh%face_neighbor(nf_total))
      allocate(mesh%face_patch   (nf_total))
      allocate(mesh%face_area    (nf_total))
      allocate(mesh%face_normal  (3, nf_total))
      allocate(mesh%face_centroid(3, nf_total))
      mesh%face_owner    = 0
      mesh%face_owner_lf = 0
      mesh%face_neighbor = 0
      mesh%face_patch    = 0
      mesh%face_area     = 0.0_wp
      mesh%face_normal   = 0.0_wp
      mesh%face_centroid = 0.0_wp

      ! Phase 4a: emit interior faces first (write to indices 1..nf_interior),
      ! then boundary faces (nf_interior+1..nf_total). Within boundary, we'll
      ! re-order again after matching to patches.
      allocate(bcand(nf_boundary))
      block
         integer :: ii, jj, owner_c, neigh_c, owner_lf
         ii = 0
         jj = 0
         i = 1
         do while (i <= ncand)
            if (i < ncand) then
               if (cmp_key(cand(i)%key, cand(i+1)%key) == 0) then
                  ii = ii + 1
                  if (cand(i)%cell <= cand(i+1)%cell) then
                     owner_c  = cand(i)%cell
                     owner_lf = cand(i)%lface
                     neigh_c  = cand(i+1)%cell
                  else
                     owner_c  = cand(i+1)%cell
                     owner_lf = cand(i+1)%lface
                     neigh_c  = cand(i)%cell
                  end if
                  mesh%face_owner   (ii) = owner_c
                  mesh%face_owner_lf(ii) = owner_lf
                  mesh%face_neighbor(ii) = neigh_c
                  i = i + 2
                  cycle
               end if
            end if
            jj = jj + 1
            bcand(jj)%key   = cand(i)%key
            bcand(jj)%nvf   = cand(i)%nvf
            bcand(jj)%cell  = cand(i)%cell
            bcand(jj)%lface = cand(i)%lface
            i = i + 1
         end do
      end block

      ! Phase 4b: match boundary candidates to Gmsh boundary elements by sorted key
      block
         type(t_face_cand), allocatable :: bf_keys(:)
         integer :: nb_local, m, kk, lo, hi, mid
         integer :: pkey(max_vpf)
         integer :: face_phys_tag(nf_boundary)

         face_phys_tag = 0
         nb_local = nb
         allocate(bf_keys(nb_local))
         do m = 1, nb_local
            do kk = 1, bfaces(m)%nvf
               pkey(kk) = bfaces(m)%vtx(kk)
            end do
            do kk = bfaces(m)%nvf+1, max_vpf
               pkey(kk) = 0
            end do
            call sort_key(pkey, bfaces(m)%nvf)
            bf_keys(m)%key   = pkey
            bf_keys(m)%nvf   = bfaces(m)%nvf
            bf_keys(m)%cell  = bfaces(m)%gmsh_tag    ! reuse 'cell' as tag carrier
            bf_keys(m)%lface = m
         end do
         call qsort_cands(bf_keys, 1, nb_local)

         do j = 1, nf_boundary
            ! binary search in bf_keys for matching key
            lo = 1; hi = nb_local
            do while (lo <= hi)
               mid = (lo + hi) / 2
               select case (cmp_key(bf_keys(mid)%key, bcand(j)%key))
               case (-1)
                  lo = mid + 1
               case (1)
                  hi = mid - 1
               case (0)
                  face_phys_tag(j) = bf_keys(mid)%cell   ! gmsh_tag
                  exit
               end select
            end do
         end do

         ! Reorder boundary candidates by patch (gmsh_tag) order in mesh%patches.
         ! For any boundary face with unknown tag (face_phys_tag == 0), tag = -1 (orphan).
         call reorder_and_emit_boundary(mesh, bcand, face_phys_tag, nf_boundary, nf_interior)
      end block

   end subroutine build_topology

   subroutine reorder_and_emit_boundary(mesh, bcand, face_phys_tag, nfb, nfi)
      type(t_mesh),      intent(inout) :: mesh
      type(t_face_cand), intent(in)    :: bcand(:)
      integer,           intent(in)    :: face_phys_tag(:)
      integer,           intent(in)    :: nfb, nfi

      integer :: ip, p, np, j, idx

      np = mesh%np
      ! Reset patches' face_count
      do p = 1, np
         mesh%patches(p)%face_count = 0
      end do

      ! First pass: count per patch
      do j = 1, nfb
         ip = 0
         do p = 1, np
            if (face_phys_tag(j) == mesh%patches(p)%gmsh_tag) then
               ip = p
               exit
            end if
         end do
         if (ip > 0) mesh%patches(ip)%face_count = mesh%patches(ip)%face_count + 1
      end do

      ! Assign face_start for each patch
      idx = nfi
      do p = 1, np
         mesh%patches(p)%face_start = idx + 1
         idx = idx + mesh%patches(p)%face_count
      end do

      ! Second pass: emit
      ! Use running counters per patch; orphans get unique slots after all patches.
      block
         integer, allocatable :: pcnt(:)
         integer :: orphan_start, orphan_cursor
         allocate(pcnt(np))
         pcnt = 0
         orphan_start = nfi
         do p = 1, np
            orphan_start = orphan_start + mesh%patches(p)%face_count
         end do
         orphan_cursor = orphan_start
         do j = 1, nfb
            ip = 0
            do p = 1, np
               if (face_phys_tag(j) == mesh%patches(p)%gmsh_tag) then
                  ip = p
                  exit
               end if
            end do
            if (ip > 0) then
               pcnt(ip) = pcnt(ip) + 1
               idx = mesh%patches(ip)%face_start + pcnt(ip) - 1
            else
               orphan_cursor = orphan_cursor + 1
               idx = orphan_cursor
            end if
            mesh%face_owner   (idx) = bcand(j)%cell
            mesh%face_owner_lf(idx) = bcand(j)%lface
            mesh%face_neighbor(idx) = 0
            mesh%face_patch   (idx) = ip
         end do
         deallocate(pcnt)
      end block

      ! Orphans (boundary candidates with no matching gmsh bface) are expected
      ! during MPI partition (ghost-only faces). They get unique slots but
      ! patch=0; the partition classifier in partition.f90 drops them.

      ! Tag interior faces with patch=0
      do j = 1, nfi
         mesh%face_patch(j) = 0
      end do
   end subroutine reorder_and_emit_boundary

   pure subroutine sort_key(key, n)
      integer, intent(inout) :: key(:)
      integer, intent(in)    :: n
      integer :: i, j, t
      ! Insertion sort over first n entries; n <= 4 so this is optimal.
      do i = 2, n
         t = key(i)
         j = i - 1
         do while (j >= 1)
            if (key(j) <= t) exit
            key(j+1) = key(j)
            j = j - 1
         end do
         key(j+1) = t
      end do
   end subroutine sort_key

   pure function cmp_key(a, b) result(s)
      integer, intent(in) :: a(:), b(:)
      integer :: s, i
      s = 0
      do i = 1, size(a)
         if (a(i) < b(i)) then
            s = -1
            return
         else if (a(i) > b(i)) then
            s = 1
            return
         end if
      end do
   end function cmp_key

   recursive subroutine qsort_cands(a, lo, hi)
      type(t_face_cand), intent(inout) :: a(:)
      integer,           intent(in)    :: lo, hi
      integer :: i, j
      integer :: pivot(max_vpf)
      type(t_face_cand) :: tmp
      if (lo >= hi) return
      pivot = a((lo+hi)/2)%key
      i = lo
      j = hi
      do
         do while (cmp_key(a(i)%key, pivot) < 0)
            i = i + 1
         end do
         do while (cmp_key(a(j)%key, pivot) > 0)
            j = j - 1
         end do
         if (i <= j) then
            tmp = a(i); a(i) = a(j); a(j) = tmp
            i = i + 1
            j = j - 1
         end if
         if (i > j) exit
      end do
      if (lo < j) call qsort_cands(a, lo, j)
      if (i < hi) call qsort_cands(a, i, hi)
   end subroutine qsort_cands

end module mesh_topology
