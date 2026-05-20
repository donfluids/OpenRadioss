! Implicit (curved) obstacle shapes for cut-cell handling: spheres and
! axis-aligned cylinders. Unlike the analytic cube path, these use
! subsampling to estimate the per-cell cut volume and per-face cut area,
! then the caller derives a single closure fragment from the area
! deficit. The closure fragment is exact for the NET pressure force on
! the curved surface at constant cell pressure (p·∮n dA = p·net), so the
! scheme stays conservative regardless of subsampling accuracy — only
! the absolute cut volume / area carry the O(1/nsub) sampling error.
module cut_cell_shapes
   use kinds, only : wp
   implicit none
   private

   public :: SHAPE_SPHERE, SHAPE_CYL_X, SHAPE_CYL_Y, SHAPE_CYL_Z
   public :: point_in_shape, shape_aabb, sampled_cell_deficit

   integer, parameter :: SHAPE_SPHERE = 1
   integer, parameter :: SHAPE_CYL_X  = 2
   integer, parameter :: SHAPE_CYL_Y  = 3
   integer, parameter :: SHAPE_CYL_Z  = 4

   ! Parameter conventions (p(1:6) per shape):
   !   SPHERE : p(1:3) = centre,         p(4) = radius
   !   CYL_X  : p(2:3) = (cy,cz) axis,   p(4) = radius, p(5)=x_lo, p(6)=x_hi
   !   CYL_Y  : p(1),p(3) = (cx,cz) axis,p(4) = radius, p(5)=y_lo, p(6)=y_hi
   !   CYL_Z  : p(1:2) = (cx,cy) axis,   p(4) = radius, p(5)=z_lo, p(6)=z_hi

contains

   pure logical function point_in_shape(x, kind, p) result(inside)
      real(wp), intent(in) :: x(3)
      integer,  intent(in) :: kind
      real(wp), intent(in) :: p(:)
      real(wp) :: dx, dy, dz, r2
      inside = .false.
      select case (kind)
      case (SHAPE_SPHERE)
         dx = x(1) - p(1); dy = x(2) - p(2); dz = x(3) - p(3)
         inside = (dx*dx + dy*dy + dz*dz) <= p(4)*p(4)
      case (SHAPE_CYL_X)
         dy = x(2) - p(2); dz = x(3) - p(3)
         r2 = dy*dy + dz*dz
         inside = (r2 <= p(4)*p(4)) .and. (x(1) >= p(5)) .and. (x(1) <= p(6))
      case (SHAPE_CYL_Y)
         dx = x(1) - p(1); dz = x(3) - p(3)
         r2 = dx*dx + dz*dz
         inside = (r2 <= p(4)*p(4)) .and. (x(2) >= p(5)) .and. (x(2) <= p(6))
      case (SHAPE_CYL_Z)
         dx = x(1) - p(1); dy = x(2) - p(2)
         r2 = dx*dx + dy*dy
         inside = (r2 <= p(4)*p(4)) .and. (x(3) >= p(5)) .and. (x(3) <= p(6))
      end select
   end function point_in_shape

   ! Axis-aligned bounding box of a shape, for cheap cell rejection.
   pure subroutine shape_aabb(kind, p, lo, hi)
      integer,  intent(in)  :: kind
      real(wp), intent(in)  :: p(:)
      real(wp), intent(out) :: lo(3), hi(3)
      select case (kind)
      case (SHAPE_SPHERE)
         lo = p(1:3) - p(4)
         hi = p(1:3) + p(4)
      case (SHAPE_CYL_X)
         lo = [p(5),      p(2)-p(4), p(3)-p(4)]
         hi = [p(6),      p(2)+p(4), p(3)+p(4)]
      case (SHAPE_CYL_Y)
         lo = [p(1)-p(4), p(5),      p(3)-p(4)]
         hi = [p(1)+p(4), p(6),      p(3)+p(4)]
      case (SHAPE_CYL_Z)
         lo = [p(1)-p(4), p(2)-p(4), p(5)]
         hi = [p(1)+p(4), p(2)+p(4), p(6)]
      case default
         lo = 0.0_wp; hi = 0.0_wp
      end select
   end subroutine shape_aabb

   ! Subsample a cell to estimate the cut volume and the per-face cut
   ! area for a single implicit shape. A_def is ordered 1:-x 2:+x 3:-y
   ! 4:+y 5:-z 6:+z (matching cut_cell_geom). Returns zeros when the
   ! cell's AABB does not meet the shape's AABB.
   pure subroutine sampled_cell_deficit(cell_lo, cell_hi, kind, p, nsub, V_def, A_def)
      real(wp), intent(in)  :: cell_lo(3), cell_hi(3)
      integer,  intent(in)  :: kind, nsub
      real(wp), intent(in)  :: p(:)
      real(wp), intent(out) :: V_def, A_def(6)

      real(wp) :: d(3), V_full, slo(3), shi(3)
      real(wp) :: x(3)
      integer  :: i, j, k, cnt
      real(wp) :: inv

      V_def = 0.0_wp
      A_def = 0.0_wp

      ! Cheap reject: cell AABB vs shape AABB.
      call shape_aabb(kind, p, slo, shi)
      if (cell_hi(1) < slo(1) .or. cell_lo(1) > shi(1) .or. &
          cell_hi(2) < slo(2) .or. cell_lo(2) > shi(2) .or. &
          cell_hi(3) < slo(3) .or. cell_lo(3) > shi(3)) return

      d = cell_hi - cell_lo
      V_full = d(1) * d(2) * d(3)
      inv = 1.0_wp / real(nsub, wp)

      ! Volume fraction inside the shape.
      cnt = 0
      do k = 1, nsub
         x(3) = cell_lo(3) + (real(k, wp) - 0.5_wp) * inv * d(3)
         do j = 1, nsub
            x(2) = cell_lo(2) + (real(j, wp) - 0.5_wp) * inv * d(2)
            do i = 1, nsub
               x(1) = cell_lo(1) + (real(i, wp) - 0.5_wp) * inv * d(1)
               if (point_in_shape(x, kind, p)) cnt = cnt + 1
            end do
         end do
      end do
      V_def = (real(cnt, wp) / real(nsub, wp)**3) * V_full

      ! Per-face area fractions. Each face is sampled nsub^2 in-plane.
      A_def(1) = face_frac(1, cell_lo(1)) * d(2) * d(3)   ! -x
      A_def(2) = face_frac(1, cell_hi(1)) * d(2) * d(3)   ! +x
      A_def(3) = face_frac(2, cell_lo(2)) * d(1) * d(3)   ! -y
      A_def(4) = face_frac(2, cell_hi(2)) * d(1) * d(3)   ! +y
      A_def(5) = face_frac(3, cell_lo(3)) * d(1) * d(2)   ! -z
      A_def(6) = face_frac(3, cell_hi(3)) * d(1) * d(2)   ! +z

   contains

      ! Fraction of a face (perpendicular to axis `ax`, at coordinate
      ! `coord`) whose subsample points lie inside the shape.
      pure real(wp) function face_frac(ax, coord) result(frac)
         integer,  intent(in) :: ax
         real(wp), intent(in) :: coord
         real(wp) :: xx(3)
         integer  :: a, b, m, n2
         integer  :: u_ax, v_ax
         frac = 0.0_wp
         ! the two in-plane axes
         u_ax = mod(ax,     3) + 1
         v_ax = mod(ax + 1, 3) + 1
         m = 0
         do a = 1, nsub
            do b = 1, nsub
               xx(ax)   = coord
               xx(u_ax) = cell_lo(u_ax) + (real(a, wp) - 0.5_wp) * inv * d(u_ax)
               xx(v_ax) = cell_lo(v_ax) + (real(b, wp) - 0.5_wp) * inv * d(v_ax)
               if (point_in_shape(xx, kind, p)) m = m + 1
            end do
         end do
         n2 = nsub * nsub
         frac = real(m, wp) / real(n2, wp)
      end function face_frac

   end subroutine sampled_cell_deficit

end module cut_cell_shapes
