! Analytic geometry for axis-aligned cube obstacles cutting an
! axis-aligned hex cell. Building block for cut-cell handling of
! arbitrary-shaped obstacle surfaces (Phase 1: cubes only).
!
! Given an air cell extent [cell_lo, cell_hi] and a cube obstacle
! extent [cube_lo, cube_hi], compute:
!   V_full / V_eff   — original / clipped cell volumes
!   A_full / A_eff   — original / clipped face areas, ordered
!                      1: -x, 2: +x, 3: -y, 4: +y, 5: -z, 6: +z
!                      (face area still covered by cube is removed)
!   n_frags           — number of obstacle-surface fragments inside the
!                      cell, 0..6
!   frag_normal(:, k) — outward unit normal FROM THE CELL (i.e. pointing
!                      into the cube), one of ±e_x, ±e_y, ±e_z
!   frag_area(k)      — area of the k-th rectangular fragment
!   frag_centroid(:,k)— geometric centroid of the fragment
!
! All formulas are exact analytic clipping — no numerical integration.
! Single-cube; cells affected by multiple non-overlapping cubes can be
! handled by superposition (sum V_overlap, sum A overlaps, concatenate
! fragment lists) in the caller.
module cut_cell_geom
   use kinds, only : wp
   implicit none
   private

   public :: MAX_FRAGS_PER_CELL, cube_cell_intersection

   integer, parameter :: MAX_FRAGS_PER_CELL = 6

contains

   pure subroutine cube_cell_intersection( &
        cell_lo, cell_hi, cube_lo, cube_hi, &
        V_full, V_eff, A_full, A_eff, &
        n_frags, frag_normal, frag_area, frag_centroid)
      real(wp), intent(in)  :: cell_lo(3), cell_hi(3)
      real(wp), intent(in)  :: cube_lo(3), cube_hi(3)
      real(wp), intent(out) :: V_full, V_eff
      real(wp), intent(out) :: A_full(6), A_eff(6)
      integer,  intent(out) :: n_frags
      real(wp), intent(out) :: frag_normal  (3, MAX_FRAGS_PER_CELL)
      real(wp), intent(out) :: frag_area    (   MAX_FRAGS_PER_CELL)
      real(wp), intent(out) :: frag_centroid(3, MAX_FRAGS_PER_CELL)

      real(wp) :: d_cell(3), inter_lo(3), inter_hi(3), d_inter(3)
      integer  :: ax, ax_o1, ax_o2
      real(wp) :: area_frag, cen(3), nrm(3)

      d_cell = cell_hi - cell_lo
      V_full = d_cell(1) * d_cell(2) * d_cell(3)
      A_full(1) = d_cell(2) * d_cell(3)
      A_full(2) = d_cell(2) * d_cell(3)
      A_full(3) = d_cell(1) * d_cell(3)
      A_full(4) = d_cell(1) * d_cell(3)
      A_full(5) = d_cell(1) * d_cell(2)
      A_full(6) = d_cell(1) * d_cell(2)

      V_eff   = V_full
      A_eff   = A_full
      n_frags = 0
      frag_normal   = 0.0_wp
      frag_area     = 0.0_wp
      frag_centroid = 0.0_wp

      inter_lo(1) = max(cell_lo(1), cube_lo(1))
      inter_lo(2) = max(cell_lo(2), cube_lo(2))
      inter_lo(3) = max(cell_lo(3), cube_lo(3))
      inter_hi(1) = min(cell_hi(1), cube_hi(1))
      inter_hi(2) = min(cell_hi(2), cube_hi(2))
      inter_hi(3) = min(cell_hi(3), cube_hi(3))
      d_inter(1) = max(inter_hi(1) - inter_lo(1), 0.0_wp)
      d_inter(2) = max(inter_hi(2) - inter_lo(2), 0.0_wp)
      d_inter(3) = max(inter_hi(3) - inter_lo(3), 0.0_wp)

      ! Disjoint? Done — return defaults (no overlap).
      if (d_inter(1) <= 0.0_wp .or. d_inter(2) <= 0.0_wp .or. &
          d_inter(3) <= 0.0_wp) return

      ! Clip the cell volume.
      V_eff = V_full - d_inter(1) * d_inter(2) * d_inter(3)

      ! Clip the six cell face areas: on a cell face whose plane lies
      ! inside the cube's perpendicular-axis extent, the 2D intersection
      ! with the cube has area = d_inter(other_axes_product).
      if (cube_lo(1) <= cell_lo(1) .and. cell_lo(1) <= cube_hi(1)) &
         A_eff(1) = A_full(1) - d_inter(2) * d_inter(3)
      if (cube_lo(1) <= cell_hi(1) .and. cell_hi(1) <= cube_hi(1)) &
         A_eff(2) = A_full(2) - d_inter(2) * d_inter(3)
      if (cube_lo(2) <= cell_lo(2) .and. cell_lo(2) <= cube_hi(2)) &
         A_eff(3) = A_full(3) - d_inter(1) * d_inter(3)
      if (cube_lo(2) <= cell_hi(2) .and. cell_hi(2) <= cube_hi(2)) &
         A_eff(4) = A_full(4) - d_inter(1) * d_inter(3)
      if (cube_lo(3) <= cell_lo(3) .and. cell_lo(3) <= cube_hi(3)) &
         A_eff(5) = A_full(5) - d_inter(1) * d_inter(2)
      if (cube_lo(3) <= cell_hi(3) .and. cell_hi(3) <= cube_hi(3)) &
         A_eff(6) = A_full(6) - d_inter(1) * d_inter(2)

      ! Obstacle-surface fragments: a cube face contributes a fragment to
      ! this cell iff the cube face's plane passes strictly through the
      ! cell interior (cell_lo(ax) < cube_face_coord < cell_hi(ax)).
      ! Outward normal direction (from cell into cube interior):
      !   cube's "lo" face on axis ax → cell sits on -ax side → +e_ax
      !   cube's "hi" face on axis ax → cell sits on +ax side → -e_ax
      do ax = 1, 3
         ax_o1 = mod(ax,     3) + 1
         ax_o2 = mod(ax + 1, 3) + 1
         area_frag = d_inter(ax_o1) * d_inter(ax_o2)
         if (area_frag <= 0.0_wp) cycle

         if (cell_lo(ax) < cube_lo(ax) .and. cube_lo(ax) < cell_hi(ax)) then
            n_frags = n_frags + 1
            nrm = 0.0_wp; nrm(ax) = +1.0_wp
            cen(ax)    = cube_lo(ax)
            cen(ax_o1) = 0.5_wp * (inter_lo(ax_o1) + inter_hi(ax_o1))
            cen(ax_o2) = 0.5_wp * (inter_lo(ax_o2) + inter_hi(ax_o2))
            frag_normal  (:, n_frags) = nrm
            frag_area    (   n_frags) = area_frag
            frag_centroid(:, n_frags) = cen
         end if
         if (cell_lo(ax) < cube_hi(ax) .and. cube_hi(ax) < cell_hi(ax)) then
            n_frags = n_frags + 1
            nrm = 0.0_wp; nrm(ax) = -1.0_wp
            cen(ax)    = cube_hi(ax)
            cen(ax_o1) = 0.5_wp * (inter_lo(ax_o1) + inter_hi(ax_o1))
            cen(ax_o2) = 0.5_wp * (inter_lo(ax_o2) + inter_hi(ax_o2))
            frag_normal  (:, n_frags) = nrm
            frag_area    (   n_frags) = area_frag
            frag_centroid(:, n_frags) = cen
         end if
      end do
   end subroutine cube_cell_intersection

end module cut_cell_geom
