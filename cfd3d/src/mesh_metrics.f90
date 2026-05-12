module mesh_metrics
   use kinds,      only : wp
   use mesh_types, only : t_mesh, max_vpf, n_vtx_per_cell, &
                          CELL_TET, CELL_HEX, CELL_PRISM, CELL_PYRAMID
   implicit none
   private

   public :: build_metrics

   ! Local-face templates (duplicated from mesh_topology so we don't need to
   ! expose its private data). Kept lockstep by tests.
   integer, parameter :: TET_NV(4) = [3, 3, 3, 3]
   integer, parameter :: TET_LV(4,4) = reshape( [ &
      1, 2, 3, 0,  &
      1, 2, 4, 0,  &
      2, 3, 4, 0,  &
      1, 3, 4, 0 ], shape=[4,4])
   integer, parameter :: HEX_NV(6) = [4,4,4,4,4,4]
   integer, parameter :: HEX_LV(4,6) = reshape( [ &
      1, 2, 3, 4,  &
      5, 6, 7, 8,  &
      1, 2, 6, 5,  &
      2, 3, 7, 6,  &
      3, 4, 8, 7,  &
      4, 1, 5, 8 ], shape=[4,6])
   integer, parameter :: PRI_NV(5) = [3,3,4,4,4]
   integer, parameter :: PRI_LV(4,5) = reshape( [ &
      1, 2, 3, 0,  &
      4, 5, 6, 0,  &
      1, 2, 5, 4,  &
      2, 3, 6, 5,  &
      3, 1, 4, 6 ], shape=[4,5])
   integer, parameter :: PYR_NV(5) = [4,3,3,3,3]
   integer, parameter :: PYR_LV(4,5) = reshape( [ &
      1, 2, 3, 4,  &
      1, 2, 5, 0,  &
      2, 3, 5, 0,  &
      3, 4, 5, 0,  &
      4, 1, 5, 0 ], shape=[4,5])

contains

   subroutine build_metrics(mesh)
      type(t_mesh), intent(inout) :: mesh

      integer :: c, f, nc, nf, i
      real(wp) :: tmp

      nc = mesh%nc_total
      nf = mesh%nf

      allocate(mesh%cell_volume(nc), source=0.0_wp)
      allocate(mesh%cell_centroid(3, nc), source=0.0_wp)
      allocate(mesh%cell_perm(nc))
      do i = 1, nc
         mesh%cell_perm(i) = i
      end do

      call compute_cell_centroids(mesh)
      call compute_face_metrics(mesh)
      call orient_face_normals(mesh)

      ! Cell volumes via divergence theorem using oriented faces.
      ! For each face f: contribution = (1/3) * (centroid · normal) * area
      ! Owner gets +; neighbor gets - (because normal points away from owner).
      do f = 1, nf
         tmp = ( mesh%face_centroid(1,f) * mesh%face_normal(1,f) &
               + mesh%face_centroid(2,f) * mesh%face_normal(2,f) &
               + mesh%face_centroid(3,f) * mesh%face_normal(3,f) ) &
               * mesh%face_area(f) / 3.0_wp
         c = mesh%face_owner(f)
         mesh%cell_volume(c) = mesh%cell_volume(c) + tmp
         if (mesh%face_neighbor(f) > 0) then
            c = mesh%face_neighbor(f)
            mesh%cell_volume(c) = mesh%cell_volume(c) - tmp
         end if
      end do
   end subroutine build_metrics

   subroutine compute_cell_centroids(mesh)
      type(t_mesh), intent(inout) :: mesh
      integer :: c, k, vptr, nv
      real(wp) :: s(3)
      do c = 1, mesh%nc_total
         nv = n_vtx_per_cell(mesh%cell_type(c))
         vptr = mesh%cell_vtx_ptr(c) - 1
         s = 0.0_wp
         do k = 1, nv
            s = s + mesh%xv(:, mesh%cell_vtx(vptr + k))
         end do
         mesh%cell_centroid(:, c) = s / real(nv, wp)
      end do
   end subroutine compute_cell_centroids

   subroutine compute_face_metrics(mesh)
      type(t_mesh), intent(inout) :: mesh
      integer :: f, c, lf, nvf, k, vptr
      integer :: lv(max_vpf)
      real(wp) :: x(3, max_vpf), cen(3), nrml(3), area

      do f = 1, mesh%nf
         c  = mesh%face_owner(f)
         lf = mesh%face_owner_lf(f)
         call get_local_face(mesh%cell_type(c), lf, nvf, lv)
         vptr = mesh%cell_vtx_ptr(c) - 1
         do k = 1, nvf
            x(:, k) = mesh%xv(:, mesh%cell_vtx(vptr + lv(k)))
         end do
         cen = 0.0_wp
         do k = 1, nvf
            cen = cen + x(:, k)
         end do
         cen = cen / real(nvf, wp)
         mesh%face_centroid(:, f) = cen
         if (nvf == 3) then
            call tri_normal(x(:,1), x(:,2), x(:,3), nrml, area)
         else
            call quad_normal(x(:,1), x(:,2), x(:,3), x(:,4), nrml, area)
         end if
         mesh%face_normal(:, f) = nrml
         mesh%face_area(f) = area
      end do
   end subroutine compute_face_metrics

   subroutine get_local_face(ct, lf, nvf, lv)
      integer, intent(in)  :: ct, lf
      integer, intent(out) :: nvf
      integer, intent(out) :: lv(max_vpf)
      lv = 0
      select case (ct)
      case (CELL_TET)
         nvf = TET_NV(lf)
         lv(1:nvf) = TET_LV(1:nvf, lf)
      case (CELL_HEX)
         nvf = HEX_NV(lf)
         lv(1:nvf) = HEX_LV(1:nvf, lf)
      case (CELL_PRISM)
         nvf = PRI_NV(lf)
         lv(1:nvf) = PRI_LV(1:nvf, lf)
      case (CELL_PYRAMID)
         nvf = PYR_NV(lf)
         lv(1:nvf) = PYR_LV(1:nvf, lf)
      case default
         nvf = 0
      end select
   end subroutine get_local_face

   pure subroutine tri_normal(a, b, c, n, area)
      real(wp), intent(in)  :: a(3), b(3), c(3)
      real(wp), intent(out) :: n(3), area
      real(wp) :: e1(3), e2(3), cr(3), nrm
      e1 = b - a
      e2 = c - a
      cr(1) = e1(2)*e2(3) - e1(3)*e2(2)
      cr(2) = e1(3)*e2(1) - e1(1)*e2(3)
      cr(3) = e1(1)*e2(2) - e1(2)*e2(1)
      nrm = sqrt(cr(1)**2 + cr(2)**2 + cr(3)**2)
      area = 0.5_wp * nrm
      if (nrm > 0.0_wp) then
         n = cr / nrm
      else
         n = 0.0_wp
      end if
   end subroutine tri_normal

   pure subroutine quad_normal(a, b, c, d, n, area)
      real(wp), intent(in)  :: a(3), b(3), c(3), d(3)
      real(wp), intent(out) :: n(3), area
      real(wp) :: d1(3), d2(3), cr(3), nrm
      d1 = c - a
      d2 = d - b
      cr(1) = d1(2)*d2(3) - d1(3)*d2(2)
      cr(2) = d1(3)*d2(1) - d1(1)*d2(3)
      cr(3) = d1(1)*d2(2) - d1(2)*d2(1)
      nrm = sqrt(cr(1)**2 + cr(2)**2 + cr(3)**2)
      area = 0.5_wp * nrm
      if (nrm > 0.0_wp) then
         n = cr / nrm
      else
         n = 0.0_wp
      end if
   end subroutine quad_normal

   subroutine orient_face_normals(mesh)
      type(t_mesh), intent(inout) :: mesh
      integer :: f, c
      real(wp) :: d(3), dot
      do f = 1, mesh%nf
         c = mesh%face_owner(f)
         d = mesh%face_centroid(:, f) - mesh%cell_centroid(:, c)
         dot = d(1)*mesh%face_normal(1,f) + d(2)*mesh%face_normal(2,f) + d(3)*mesh%face_normal(3,f)
         if (dot < 0.0_wp) then
            mesh%face_normal(:, f) = -mesh%face_normal(:, f)
         end if
      end do
   end subroutine orient_face_normals

end module mesh_metrics
