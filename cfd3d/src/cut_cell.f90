! Cut-cell table construction. After build_metrics, this populates
! cell_vol_eff, face_area_eff, and the per-cell obstacle-fragment CSR
! by intersecting each cell / face with the list of axis-aligned cube
! obstacles supplied via the namelist.
!
! For meshes with no obstacle cubes, the tables are still allocated
! and initialised to the un-cut volume / face area / zero fragments,
! so the solver can always read mesh%cell_vol_eff / face_area_eff and
! iterate mesh%frag_* without a special-case branch.
!
! Phase 1b: structures only. The flux-loop integration (use V_eff in
! the cell update, add Σ-fragments slip-wall flux) lands in Phase 1c.
!
! Cells are assumed axis-aligned hexes for cut-cell purposes; the cell
! extent is taken as the AABB of its 8 vertices. Tet / prism / pyramid
! cells will simply see no overlap (or wrong overlap) and need a
! follow-up generalisation if mixed-element cut-cell is ever needed.
module cut_cell
   use kinds,         only : wp
   use mesh_types,    only : t_mesh
   use cut_cell_geom, only : MAX_FRAGS_PER_CELL, cube_cell_intersection
   implicit none
   private

   public :: build_cut_cell_tables
   public :: get_axis_aligned_extent, classify_face_axis_side, axis_side_idx

contains

   subroutine build_cut_cell_tables(mesh, n_cubes, cube_lo, cube_hi)
      type(t_mesh), intent(inout) :: mesh
      integer,      intent(in)    :: n_cubes
      real(wp),     intent(in)    :: cube_lo(:,:)   ! (3, n_cubes)
      real(wp),     intent(in)    :: cube_hi(:,:)   ! (3, n_cubes)

      integer  :: c, f, k, j, ax, sgn, idx
      real(wp) :: cell_lo_c(3), cell_hi_c(3)
      real(wp) :: V_full, V_eff_c, A_full(6), A_eff_c(6)
      integer  :: nfrag_c
      real(wp) :: fn_c  (3, MAX_FRAGS_PER_CELL)
      real(wp) :: fa_c  (   MAX_FRAGS_PER_CELL)
      real(wp) :: fcen_c(3, MAX_FRAGS_PER_CELL)
      real(wp), allocatable :: A_eff_cell(:,:)         ! (6, nc_total)
      integer,  allocatable :: cell_nfrag(:)
      integer,  allocatable :: pos(:)
      integer  :: n_frags_total

      ! 1) Allocate / initialise the per-cell and per-face tables to the
      ! un-cut quantities. Subsequent passes only subtract from these.
      if (allocated(mesh%cell_vol_eff))  deallocate(mesh%cell_vol_eff)
      if (allocated(mesh%face_area_eff)) deallocate(mesh%face_area_eff)
      allocate(mesh%cell_vol_eff (mesh%nc_total))
      allocate(mesh%face_area_eff(mesh%nf))
      mesh%cell_vol_eff  = mesh%cell_volume
      mesh%face_area_eff = mesh%face_area

      allocate(cell_nfrag(mesh%nc_total), source=0)
      allocate(A_eff_cell(6, mesh%nc_total))

      ! 2) Per-cell pass: extent → for each cube, accumulate V_overlap +
      ! per-axis-side A_overlap, count fragments.
      do c = 1, mesh%nc_total
         call get_axis_aligned_extent(mesh, c, cell_lo_c, cell_hi_c)
         A_eff_cell(1, c) = (cell_hi_c(2)-cell_lo_c(2)) * (cell_hi_c(3)-cell_lo_c(3))
         A_eff_cell(2, c) = A_eff_cell(1, c)
         A_eff_cell(3, c) = (cell_hi_c(1)-cell_lo_c(1)) * (cell_hi_c(3)-cell_lo_c(3))
         A_eff_cell(4, c) = A_eff_cell(3, c)
         A_eff_cell(5, c) = (cell_hi_c(1)-cell_lo_c(1)) * (cell_hi_c(2)-cell_lo_c(2))
         A_eff_cell(6, c) = A_eff_cell(5, c)

         do k = 1, n_cubes
            call cube_cell_intersection(cell_lo_c, cell_hi_c, cube_lo(:,k), cube_hi(:,k), &
                 V_full, V_eff_c, A_full, A_eff_c, nfrag_c, fn_c, fa_c, fcen_c)
            if (V_full - V_eff_c > 0.0_wp) then
               mesh%cell_vol_eff(c) = mesh%cell_vol_eff(c) - (V_full - V_eff_c)
               do j = 1, 6
                  A_eff_cell(j, c) = A_eff_cell(j, c) - (A_full(j) - A_eff_c(j))
               end do
               cell_nfrag(c) = cell_nfrag(c) + nfrag_c
            end if
         end do
      end do

      ! 3) Face pass: copy the owner's per-axis-side A_eff to face_area_eff.
      ! Each face is shared between (owner, neighbor) for interior faces but
      ! the obstacle clipping is the same from either side, so reading off
      ! the owner only is correct (no double-counting).
      do f = 1, mesh%nf
         c = mesh%face_owner(f)
         call classify_face_axis_side(mesh%face_normal(:,f), ax, sgn)
         idx = axis_side_idx(ax, sgn)
         mesh%face_area_eff(f) = A_eff_cell(idx, c)
      end do

      ! 4) Build the fragment CSR offsets, allocate fragment arrays.
      if (allocated(mesh%cell_frag_offset)) deallocate(mesh%cell_frag_offset)
      allocate(mesh%cell_frag_offset(mesh%nc_total + 1))
      mesh%cell_frag_offset(1) = 1
      do c = 1, mesh%nc_total
         mesh%cell_frag_offset(c+1) = mesh%cell_frag_offset(c) + cell_nfrag(c)
      end do
      n_frags_total = mesh%cell_frag_offset(mesh%nc_total + 1) - 1
      mesh%n_frags = n_frags_total

      if (allocated(mesh%frag_normal))   deallocate(mesh%frag_normal)
      if (allocated(mesh%frag_area))     deallocate(mesh%frag_area)
      if (allocated(mesh%frag_centroid)) deallocate(mesh%frag_centroid)
      allocate(mesh%frag_normal  (3, max(n_frags_total, 1)))
      allocate(mesh%frag_area    (   max(n_frags_total, 1)))
      allocate(mesh%frag_centroid(3, max(n_frags_total, 1)))

      ! 5) Fill the fragment CSR by re-iterating (cube_cell_intersection
      ! is cheap; caching the results from pass 2 would save flops at the
      ! cost of a 7×n_cubes×nc temp buffer — premature for n_cubes ~ 1).
      allocate(pos(mesh%nc_total))
      do c = 1, mesh%nc_total
         pos(c) = mesh%cell_frag_offset(c)
      end do
      do c = 1, mesh%nc_total
         call get_axis_aligned_extent(mesh, c, cell_lo_c, cell_hi_c)
         do k = 1, n_cubes
            call cube_cell_intersection(cell_lo_c, cell_hi_c, cube_lo(:,k), cube_hi(:,k), &
                 V_full, V_eff_c, A_full, A_eff_c, nfrag_c, fn_c, fa_c, fcen_c)
            do j = 1, nfrag_c
               mesh%frag_normal  (:, pos(c)) = fn_c(:, j)
               mesh%frag_area    (   pos(c)) = fa_c(j)
               mesh%frag_centroid(:, pos(c)) = fcen_c(:, j)
               pos(c) = pos(c) + 1
            end do
         end do
      end do

      deallocate(pos, cell_nfrag, A_eff_cell)
   end subroutine build_cut_cell_tables

   ! Axis-aligned bounding box of a cell's vertex set.
   pure subroutine get_axis_aligned_extent(mesh, c, cell_lo, cell_hi)
      type(t_mesh), intent(in)  :: mesh
      integer,      intent(in)  :: c
      real(wp),     intent(out) :: cell_lo(3), cell_hi(3)
      integer  :: pt, ed, k, vid
      real(wp) :: x(3)
      pt = mesh%cell_vtx_ptr(c)
      ed = mesh%cell_vtx_ptr(c+1) - 1
      cell_lo = mesh%xv(:, mesh%cell_vtx(pt))
      cell_hi = cell_lo
      do k = pt + 1, ed
         vid = mesh%cell_vtx(k)
         x = mesh%xv(:, vid)
         cell_lo(1) = min(cell_lo(1), x(1))
         cell_lo(2) = min(cell_lo(2), x(2))
         cell_lo(3) = min(cell_lo(3), x(3))
         cell_hi(1) = max(cell_hi(1), x(1))
         cell_hi(2) = max(cell_hi(2), x(2))
         cell_hi(3) = max(cell_hi(3), x(3))
      end do
   end subroutine get_axis_aligned_extent

   ! Find the dominant axis (1..3) and its sign (±1) for a unit normal
   ! that is approximately axis-aligned (true for axis-aligned hex cell
   ! faces; falls back to the largest |component| for any other case).
   pure subroutine classify_face_axis_side(normal, ax, sgn)
      real(wp), intent(in)  :: normal(3)
      integer,  intent(out) :: ax, sgn
      real(wp) :: amax
      integer  :: k
      amax = -1.0_wp
      ax = 1
      do k = 1, 3
         if (abs(normal(k)) > amax) then
            amax = abs(normal(k))
            ax = k
         end if
      end do
      if (normal(ax) >= 0.0_wp) then
         sgn = +1
      else
         sgn = -1
      end if
   end subroutine classify_face_axis_side

   ! Convert (axis, sign) to the per-cell A_eff index used by
   ! cube_cell_intersection. Ordering matches: 1=-x, 2=+x, 3=-y, 4=+y,
   ! 5=-z, 6=+z.
   pure integer function axis_side_idx(ax, sgn) result(idx)
      integer, intent(in) :: ax, sgn
      if (sgn > 0) then
         idx = 2 * ax        ! 2, 4, 6
      else
         idx = 2 * ax - 1    ! 1, 3, 5
      end if
   end function axis_side_idx

end module cut_cell
