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
   use kinds,          only : wp
   use constants,      only : NVAR
   use mesh_types,     only : t_mesh
   use fields,         only : t_state
   use cut_cell_geom,  only : MAX_FRAGS_PER_CELL, cube_cell_intersection
   use cut_cell_shapes, only : sampled_cell_deficit
   implicit none
   private

   public :: build_cut_cell_tables
   public :: get_axis_aligned_extent, classify_face_axis_side, axis_side_idx
   public :: build_merge_groups, merge_fold_residual, merge_sync_state

   ! A cell with cut volume fraction below this is a "sliver" and gets
   ! merged into a larger neighbour for time-step stability.
   real(wp), parameter, public :: MERGE_ALPHA = 0.5_wp

   ! A cell with cut volume fraction below this is "dead" (fully inside an
   ! obstacle): flux-isolated, residual identically zero, kept frozen.
   real(wp), parameter, public :: V_DEAD_FRAC = 1.0e-8_wp

   ! Curved-obstacle subsampling resolution (nsub^3 per cell).
   integer, parameter :: NSUB = 16

contains

   subroutine build_cut_cell_tables(mesh, n_cubes, cube_lo, cube_hi, &
                                    n_shapes, shape_kind, shape_p)
      type(t_mesh), intent(inout) :: mesh
      integer,      intent(in)    :: n_cubes
      real(wp),     intent(in)    :: cube_lo(:,:)   ! (3, n_cubes)
      real(wp),     intent(in)    :: cube_hi(:,:)   ! (3, n_cubes)
      integer,  intent(in), optional :: n_shapes
      integer,  intent(in), optional :: shape_kind(:)    ! (n_shapes)
      real(wp), intent(in), optional :: shape_p(:,:)     ! (6, n_shapes)

      integer  :: c, f, k, j, ax, sgn, idx, nsh
      real(wp) :: cell_lo_c(3), cell_hi_c(3)
      real(wp) :: V_full, V_eff_c, A_full(6), A_eff_c(6)
      integer  :: nfrag_c
      real(wp) :: fn_c  (3, MAX_FRAGS_PER_CELL)
      real(wp) :: fa_c  (   MAX_FRAGS_PER_CELL)
      real(wp) :: fcen_c(3, MAX_FRAGS_PER_CELL)
      real(wp) :: V_def, A_def(6), s_nrm(3), s_area, s_cen(3)
      logical  :: s_has
      real(wp), allocatable :: A_eff_cell(:,:)         ! (6, nc_total)
      integer,  allocatable :: cell_nfrag(:)
      integer,  allocatable :: pos(:)
      integer  :: n_frags_total

      nsh = 0
      if (present(n_shapes)) nsh = n_shapes

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

         ! Curved shapes: subsampled deficit + single closure fragment.
         do k = 1, nsh
            call sampled_cell_deficit(cell_lo_c, cell_hi_c, shape_kind(k), &
                 shape_p(:,k), NSUB, V_def, A_def)
            if (V_def > 0.0_wp) then
               mesh%cell_vol_eff(c) = mesh%cell_vol_eff(c) - V_def
               do j = 1, 6
                  A_eff_cell(j, c) = A_eff_cell(j, c) - A_def(j)
               end do
               call closure_fragment(A_def, cell_lo_c, cell_hi_c, s_has, s_nrm, s_area, s_cen)
               if (s_has) cell_nfrag(c) = cell_nfrag(c) + 1
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
         do k = 1, nsh
            call sampled_cell_deficit(cell_lo_c, cell_hi_c, shape_kind(k), &
                 shape_p(:,k), NSUB, V_def, A_def)
            if (V_def > 0.0_wp) then
               call closure_fragment(A_def, cell_lo_c, cell_hi_c, s_has, s_nrm, s_area, s_cen)
               if (s_has) then
                  mesh%frag_normal  (:, pos(c)) = s_nrm
                  mesh%frag_area    (   pos(c)) = s_area
                  mesh%frag_centroid(:, pos(c)) = s_cen
                  pos(c) = pos(c) + 1
               end if
            end if
         end do
      end do

      deallocate(pos, cell_nfrag, A_eff_cell)

      ! 6) Merge groups: default each cell to its own group, then link
      ! slivers to hosts. (No-op default when there are no cut cells.)
      if (allocated(mesh%cell_merge_root)) deallocate(mesh%cell_merge_root)
      if (allocated(mesh%cell_merge_vol))  deallocate(mesh%cell_merge_vol)
      allocate(mesh%cell_merge_root(mesh%nc_total))
      allocate(mesh%cell_merge_vol (mesh%nc_total))
      do c = 1, mesh%nc_total
         mesh%cell_merge_root(c) = c
         mesh%cell_merge_vol (c) = mesh%cell_vol_eff(c)
      end do
      if (n_cubes > 0 .or. nsh > 0) call build_merge_groups(mesh, MERGE_ALPHA)
   end subroutine build_cut_cell_tables

   ! Closure fragment for a cell from its per-axis-side area deficit
   ! A_def (1:-x 2:+x 3:-y 4:+y 5:-z 6:+z). The obstacle-surface area
   ! vector that closes the cell is Σ_faces deficit·n_face_outward; for a
   ! slip wall at constant cell pressure this single planar fragment
   ! carries the exact net pressure force on the (possibly curved)
   ! obstacle surface inside the cell. Centroid is the cell centre
   ! (unused at first order on cut cells).
   pure subroutine closure_fragment(A_def, cell_lo, cell_hi, has_frag, nrm, area, cen)
      real(wp), intent(in)  :: A_def(6), cell_lo(3), cell_hi(3)
      logical,  intent(out) :: has_frag
      real(wp), intent(out) :: nrm(3), area, cen(3)
      real(wp) :: vec(3), mag
      vec(1) = A_def(2) - A_def(1)
      vec(2) = A_def(4) - A_def(3)
      vec(3) = A_def(6) - A_def(5)
      mag = sqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
      if (mag > 1.0e-14_wp) then
         has_frag = .true.
         nrm  = vec / mag
         area = mag
         cen  = 0.5_wp * (cell_lo + cell_hi)
      else
         has_frag = .false.
         nrm = 0.0_wp; area = 0.0_wp; cen = 0.0_wp
      end if
   end subroutine closure_fragment

   ! Link each sliver cell (cut volume fraction < alpha) to the largest
   ! non-sliver LOCAL face-neighbour, then recompute per-group volumes.
   ! Single-level star merging — hosts are themselves non-slivers, so no
   ! merge chains form. A sliver with no eligible host stays un-merged
   ! (it keeps its small dt; a warning is printed once).
   subroutine build_merge_groups(mesh, alpha)
      type(t_mesh), intent(inout) :: mesh
      real(wp),     intent(in)    :: alpha

      integer  :: c, f, o, n, root, n_unmerged
      real(wp), allocatable :: best_host_vol(:)
      logical,  allocatable :: is_sliver(:)

      allocate(is_sliver(mesh%nc_total), source=.false.)
      do c = 1, mesh%nc_total
         if (mesh%cell_volume(c) > 0.0_wp) then
            is_sliver(c) = (mesh%cell_vol_eff(c) / mesh%cell_volume(c)) < alpha
         end if
      end do

      ! Single pass over interior faces: for each sliver, track the
      ! best (largest V_eff) non-sliver local neighbour as its host.
      allocate(best_host_vol(mesh%nc_total), source=-1.0_wp)
      do f = 1, mesh%nf_interior
         o = mesh%face_owner(f)
         n = mesh%face_neighbor(f)
         if (n <= 0) cycle
         call consider(o, n)
         call consider(n, o)
      end do

      ! Recompute group volumes: accumulate each cell's V_eff onto its
      ! root, then propagate the group total back to every member.
      ! Count only LIVE slivers (not dead cells) that found no host —
      ! those are the ones that keep a reduced dt.
      n_unmerged = 0
      do c = 1, mesh%nc_internal
         if (is_sliver(c) .and. mesh%cell_merge_root(c) == c .and. &
             mesh%cell_vol_eff(c) > V_DEAD_FRAC * mesh%cell_volume(c)) &
            n_unmerged = n_unmerged + 1
      end do
      mesh%cell_merge_vol = 0.0_wp
      do c = 1, mesh%nc_total
         root = mesh%cell_merge_root(c)
         mesh%cell_merge_vol(root) = mesh%cell_merge_vol(root) + mesh%cell_vol_eff(c)
      end do
      do c = 1, mesh%nc_total
         mesh%cell_merge_vol(c) = mesh%cell_merge_vol(mesh%cell_merge_root(c))
      end do

      ! Dead cells (fully inside an obstacle, V_eff ~ 0) are flux-isolated
      ! — all their faces have A_eff = 0 and they carry no fragments, so
      ! their residual is identically zero. They cannot find a host (all
      ! neighbours are dead or cut), so give them a safe nonzero
      ! denominator (the full cell volume); with R = 0 they simply stay
      ! frozen at their initial state and contribute nothing.
      do c = 1, mesh%nc_total
         if (mesh%cell_volume(c) > 0.0_wp) then
            if (mesh%cell_vol_eff(c) <= V_DEAD_FRAC * mesh%cell_volume(c)) &
               mesh%cell_merge_vol(c) = mesh%cell_volume(c)
         end if
      end do

      if (n_unmerged > 0) &
         write(*,'(A,I0,A)') 'cut_cell: warning — ', n_unmerged, &
            ' sliver cell(s) had no eligible host; they keep a reduced dt'

      deallocate(best_host_vol, is_sliver)

   contains

      ! If cell s is a sliver and cell h is a non-sliver local cell with
      ! a larger V_eff than the current best, adopt h as s's host.
      subroutine consider(s, h)
         integer, intent(in) :: s, h
         if (.not. is_sliver(s)) return
         if (h > mesh%nc_internal) return     ! host must be locally owned
         if (is_sliver(h)) return             ! host must be a non-sliver
         if (mesh%cell_vol_eff(h) > best_host_vol(s)) then
            best_host_vol(s) = mesh%cell_vol_eff(h)
            mesh%cell_merge_root(s) = h
         end if
      end subroutine consider

   end subroutine build_merge_groups

   ! Fold each sliver's residual into its host (in place) so the host
   ! carries the whole group's residual and slivers carry zero. The RK
   ! update then advances every group member with R(root)/V_group.
   subroutine merge_fold_residual(mesh, s)
      type(t_mesh),  intent(in)    :: mesh
      type(t_state), intent(inout) :: s
      integer :: c, root
      do c = 1, mesh%nc_internal
         root = mesh%cell_merge_root(c)
         if (root /= c) then
            s%R(:, root) = s%R(:, root) + s%R(:, c)
            s%R(:, c)    = 0.0_wp
         end if
      end do
   end subroutine merge_fold_residual

   ! Volume-weighted average each group's conserved state and assign it
   ! to every member. Conserves Σ U·V_eff within each group. Call once
   ! after the initial condition so merged members start coherent.
   subroutine merge_sync_state(mesh, s)
      type(t_mesh),  intent(in)    :: mesh
      type(t_state), intent(inout) :: s
      integer :: c, root
      real(wp), allocatable :: Usum(:,:)
      logical :: any_merge

      any_merge = .false.
      do c = 1, mesh%nc_internal
         if (mesh%cell_merge_root(c) /= c) then
            any_merge = .true.
            exit
         end if
      end do
      if (.not. any_merge) return

      allocate(Usum(NVAR, mesh%nc_total), source=0.0_wp)
      do c = 1, mesh%nc_internal
         root = mesh%cell_merge_root(c)
         Usum(:, root) = Usum(:, root) + s%U(:, c) * mesh%cell_vol_eff(c)
      end do
      do c = 1, mesh%nc_internal
         root = mesh%cell_merge_root(c)
         if (mesh%cell_merge_vol(root) > 0.0_wp) &
            s%U(:, c) = Usum(:, root) / mesh%cell_merge_vol(root)
      end do
      deallocate(Usum)
   end subroutine merge_sync_state

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
