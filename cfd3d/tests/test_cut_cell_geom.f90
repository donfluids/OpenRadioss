! Unit tests for cut_cell_geom: axis-aligned cube ∩ axis-aligned hex
! cell. Tests every topological case (no overlap, single-face cut,
! edge cut, corner cut, cube-inside-cell, cell-inside-cube).
program test_cut_cell_geom
   use kinds, only : wp
   use cut_cell_geom, only : MAX_FRAGS_PER_CELL, cube_cell_intersection
   implicit none

   integer :: fails

   fails = 0

   call test_no_overlap                  (fails)
   call test_single_face_cut             (fails)
   call test_edge_cut                    (fails)
   call test_corner_cut                  (fails)
   call test_cube_inside_cell            (fails)
   call test_cell_inside_cube            (fails)
   call test_cube_face_flush_with_cell   (fails)

   if (fails == 0) then
      write(*,'(A)') 'test_cut_cell_geom: PASS'
   else
      write(*,'(A,I0,A)') 'test_cut_cell_geom: FAIL (', fails, ')'
      stop 1
   end if

contains

   subroutine check(name, ok, fails)
      character(len=*), intent(in)    :: name
      logical,          intent(in)    :: ok
      integer,          intent(inout) :: fails
      if (.not. ok) then
         write(*,'(A,A)') 'FAIL: ', name
         fails = fails + 1
      end if
   end subroutine check

   pure logical function close_to(a, b, tol)
      real(wp), intent(in) :: a, b, tol
      close_to = abs(a - b) <= tol
   end function close_to

   ! ----------------------------------------------------------------
   subroutine test_no_overlap(fails)
      integer, intent(inout) :: fails
      real(wp) :: cell_lo(3), cell_hi(3), cube_lo(3), cube_hi(3)
      real(wp) :: V_full, V_eff, A_full(6), A_eff(6)
      integer  :: n_frags
      real(wp) :: fn(3, MAX_FRAGS_PER_CELL), fa(MAX_FRAGS_PER_CELL)
      real(wp) :: fc(3, MAX_FRAGS_PER_CELL)
      real(wp), parameter :: TOL = 1.0e-12_wp

      cell_lo = [0.0_wp, 0.0_wp, 0.0_wp]
      cell_hi = [1.0_wp, 1.0_wp, 1.0_wp]
      cube_lo = [2.0_wp, 2.0_wp, 2.0_wp]
      cube_hi = [3.0_wp, 3.0_wp, 3.0_wp]
      call cube_cell_intersection(cell_lo, cell_hi, cube_lo, cube_hi, &
           V_full, V_eff, A_full, A_eff, n_frags, fn, fa, fc)
      call check('no_overlap: V_eff = V_full',    close_to(V_eff, V_full, TOL),  fails)
      call check('no_overlap: A_eff(1) unchanged',close_to(A_eff(1), 1.0_wp,TOL),fails)
      call check('no_overlap: A_eff(6) unchanged',close_to(A_eff(6), 1.0_wp,TOL),fails)
      call check('no_overlap: n_frags = 0',       n_frags == 0,                  fails)
   end subroutine test_no_overlap

   ! ----------------------------------------------------------------
   ! Cell [0,1]^3; cube cuts the +x face only (cube starts at x=0.5,
   ! extends beyond the cell on all sides except -x). One fragment,
   ! normal = +x_hat (cell sits on -x side of the cube).
   subroutine test_single_face_cut(fails)
      integer, intent(inout) :: fails
      real(wp) :: cell_lo(3), cell_hi(3), cube_lo(3), cube_hi(3)
      real(wp) :: V_full, V_eff, A_full(6), A_eff(6)
      integer  :: n_frags
      real(wp) :: fn(3, MAX_FRAGS_PER_CELL), fa(MAX_FRAGS_PER_CELL)
      real(wp) :: fc(3, MAX_FRAGS_PER_CELL)
      real(wp), parameter :: TOL = 1.0e-12_wp

      cell_lo = [0.0_wp, 0.0_wp, 0.0_wp]
      cell_hi = [1.0_wp, 1.0_wp, 1.0_wp]
      cube_lo = [0.5_wp, -1.0_wp, -1.0_wp]
      cube_hi = [2.0_wp,  2.0_wp,  2.0_wp]
      call cube_cell_intersection(cell_lo, cell_hi, cube_lo, cube_hi, &
           V_full, V_eff, A_full, A_eff, n_frags, fn, fa, fc)
      call check('single_face: V_eff = 0.5',        close_to(V_eff,   0.5_wp,TOL), fails)
      call check('single_face: -x face unchanged',  close_to(A_eff(1),1.0_wp,TOL), fails)
      call check('single_face: +x face fully covered', close_to(A_eff(2),0.0_wp,TOL), fails)
      call check('single_face: -y face cut in half',close_to(A_eff(3),0.5_wp,TOL), fails)
      call check('single_face: +y face cut in half',close_to(A_eff(4),0.5_wp,TOL), fails)
      call check('single_face: -z face cut in half',close_to(A_eff(5),0.5_wp,TOL), fails)
      call check('single_face: +z face cut in half',close_to(A_eff(6),0.5_wp,TOL), fails)
      call check('single_face: n_frags = 1',        n_frags == 1, fails)
      if (n_frags >= 1) then
         call check('single_face: normal = +x', &
            close_to(fn(1,1), 1.0_wp, TOL) .and. &
            close_to(fn(2,1), 0.0_wp, TOL) .and. &
            close_to(fn(3,1), 0.0_wp, TOL), fails)
         call check('single_face: area = 1',  close_to(fa(1), 1.0_wp, TOL), fails)
         call check('single_face: centroid_x = 0.5', close_to(fc(1,1), 0.5_wp, TOL), fails)
         call check('single_face: centroid_y = 0.5', close_to(fc(2,1), 0.5_wp, TOL), fails)
         call check('single_face: centroid_z = 0.5', close_to(fc(3,1), 0.5_wp, TOL), fails)
      end if
   end subroutine test_single_face_cut

   ! ----------------------------------------------------------------
   ! Cell [0,1]^3; cube cuts +x and +y faces (cube's lo corner at
   ! (0.5,0.5,-1), hi far outside in +x/+y, +z full). Two fragments:
   ! one with normal +x, one with normal +y. Quarter-cube removed.
   subroutine test_edge_cut(fails)
      integer, intent(inout) :: fails
      real(wp) :: cell_lo(3), cell_hi(3), cube_lo(3), cube_hi(3)
      real(wp) :: V_full, V_eff, A_full(6), A_eff(6)
      integer  :: n_frags
      real(wp) :: fn(3, MAX_FRAGS_PER_CELL), fa(MAX_FRAGS_PER_CELL)
      real(wp) :: fc(3, MAX_FRAGS_PER_CELL)
      real(wp), parameter :: TOL = 1.0e-12_wp

      cell_lo = [0.0_wp, 0.0_wp, 0.0_wp]
      cell_hi = [1.0_wp, 1.0_wp, 1.0_wp]
      cube_lo = [0.5_wp, 0.5_wp, -1.0_wp]
      cube_hi = [2.0_wp, 2.0_wp,  2.0_wp]
      call cube_cell_intersection(cell_lo, cell_hi, cube_lo, cube_hi, &
           V_full, V_eff, A_full, A_eff, n_frags, fn, fa, fc)
      call check('edge_cut: V_eff = 0.75',        close_to(V_eff, 0.75_wp, TOL),  fails)
      call check('edge_cut: A_eff(+x) = 0.5',     close_to(A_eff(2), 0.5_wp, TOL),fails)
      call check('edge_cut: A_eff(+y) = 0.5',     close_to(A_eff(4), 0.5_wp, TOL),fails)
      call check('edge_cut: A_eff(-z) = 0.75',    close_to(A_eff(5), 0.75_wp,TOL),fails)
      call check('edge_cut: A_eff(+z) = 0.75',    close_to(A_eff(6), 0.75_wp,TOL),fails)
      call check('edge_cut: n_frags = 2',         n_frags == 2,                   fails)
      ! sum of fragment areas
      if (n_frags == 2) then
         call check('edge_cut: total frag area = 1', &
            close_to(fa(1) + fa(2), 1.0_wp, TOL), fails)
         ! both frags should have area 0.5
         call check('edge_cut: each frag area = 0.5', &
            close_to(fa(1), 0.5_wp, TOL) .and. close_to(fa(2), 0.5_wp, TOL), fails)
      end if
   end subroutine test_edge_cut

   ! ----------------------------------------------------------------
   ! Cell [0,1]^3; cube hi corner at (0.5,0.5,0.5), lo far outside.
   ! Three cube faces (the +x, +y, +z faces at coords 0.5, 0.5, 0.5)
   ! all pass through the cell interior. Three fragments. The
   ! intersection volume = 0.5^3 = 0.125.
   subroutine test_corner_cut(fails)
      integer, intent(inout) :: fails
      real(wp) :: cell_lo(3), cell_hi(3), cube_lo(3), cube_hi(3)
      real(wp) :: V_full, V_eff, A_full(6), A_eff(6)
      integer  :: n_frags
      real(wp) :: fn(3, MAX_FRAGS_PER_CELL), fa(MAX_FRAGS_PER_CELL)
      real(wp) :: fc(3, MAX_FRAGS_PER_CELL)
      real(wp), parameter :: TOL = 1.0e-12_wp
      integer  :: k
      real(wp) :: total_area

      cell_lo = [0.0_wp, 0.0_wp, 0.0_wp]
      cell_hi = [1.0_wp, 1.0_wp, 1.0_wp]
      cube_lo = [-1.0_wp, -1.0_wp, -1.0_wp]
      cube_hi = [ 0.5_wp,  0.5_wp,  0.5_wp]
      call cube_cell_intersection(cell_lo, cell_hi, cube_lo, cube_hi, &
           V_full, V_eff, A_full, A_eff, n_frags, fn, fa, fc)
      call check('corner_cut: V_eff = 0.875',     close_to(V_eff, 0.875_wp, TOL), fails)
      call check('corner_cut: n_frags = 3',       n_frags == 3, fails)
      total_area = 0.0_wp
      do k = 1, n_frags
         total_area = total_area + fa(k)
         ! Each fragment is 0.5 x 0.5 = 0.25 m^2
         call check('corner_cut: each frag area = 0.25', close_to(fa(k), 0.25_wp, TOL), fails)
      end do
      call check('corner_cut: total frag area = 0.75', close_to(total_area, 0.75_wp, TOL), fails)
      ! All normals point in -x, -y, -z (cube hi faces; cell on +axis side)
      if (n_frags == 3) then
         call check('corner_cut: each normal is -e_ax', &
            ((close_to(fn(1,1),-1.0_wp,TOL) .and. close_to(fn(2,1),0.0_wp,TOL) .and. close_to(fn(3,1),0.0_wp,TOL)) .or. &
             (close_to(fn(2,1),-1.0_wp,TOL) .and. close_to(fn(1,1),0.0_wp,TOL) .and. close_to(fn(3,1),0.0_wp,TOL)) .or. &
             (close_to(fn(3,1),-1.0_wp,TOL) .and. close_to(fn(1,1),0.0_wp,TOL) .and. close_to(fn(2,1),0.0_wp,TOL))), fails)
      end if
   end subroutine test_corner_cut

   ! ----------------------------------------------------------------
   ! Small cube entirely inside the cell: 6 fragments (all cube
   ! faces), V_overlap = cube volume.
   subroutine test_cube_inside_cell(fails)
      integer, intent(inout) :: fails
      real(wp) :: cell_lo(3), cell_hi(3), cube_lo(3), cube_hi(3)
      real(wp) :: V_full, V_eff, A_full(6), A_eff(6)
      integer  :: n_frags, k
      real(wp) :: fn(3, MAX_FRAGS_PER_CELL), fa(MAX_FRAGS_PER_CELL)
      real(wp) :: fc(3, MAX_FRAGS_PER_CELL)
      real(wp) :: total_area
      real(wp), parameter :: TOL = 1.0e-12_wp

      cell_lo = [0.0_wp, 0.0_wp, 0.0_wp]
      cell_hi = [1.0_wp, 1.0_wp, 1.0_wp]
      cube_lo = [0.4_wp, 0.4_wp, 0.4_wp]
      cube_hi = [0.6_wp, 0.6_wp, 0.6_wp]
      call cube_cell_intersection(cell_lo, cell_hi, cube_lo, cube_hi, &
           V_full, V_eff, A_full, A_eff, n_frags, fn, fa, fc)
      call check('cube_in_cell: V_eff = 1 - 0.008', &
         close_to(V_eff, 1.0_wp - 0.008_wp, TOL), fails)
      call check('cube_in_cell: all 6 faces unchanged', &
         all([(close_to(A_eff(k),A_full(k),TOL), k=1,6)]), fails)
      call check('cube_in_cell: n_frags = 6', n_frags == 6, fails)
      total_area = 0.0_wp
      do k = 1, n_frags
         total_area = total_area + fa(k)
         call check('cube_in_cell: each frag area = 0.04', close_to(fa(k), 0.04_wp, TOL), fails)
      end do
      call check('cube_in_cell: total frag area = cube surface = 0.24', &
         close_to(total_area, 0.24_wp, TOL), fails)
   end subroutine test_cube_inside_cell

   ! ----------------------------------------------------------------
   ! Cube engulfs the cell. Degenerate but valid: V_eff = 0, all face
   ! areas zero, no fragments (cube faces don't pass through cell
   ! interior since cell is fully enclosed).
   subroutine test_cell_inside_cube(fails)
      integer, intent(inout) :: fails
      real(wp) :: cell_lo(3), cell_hi(3), cube_lo(3), cube_hi(3)
      real(wp) :: V_full, V_eff, A_full(6), A_eff(6)
      integer  :: n_frags, k
      real(wp) :: fn(3, MAX_FRAGS_PER_CELL), fa(MAX_FRAGS_PER_CELL)
      real(wp) :: fc(3, MAX_FRAGS_PER_CELL)
      real(wp), parameter :: TOL = 1.0e-12_wp

      cell_lo = [0.0_wp, 0.0_wp, 0.0_wp]
      cell_hi = [1.0_wp, 1.0_wp, 1.0_wp]
      cube_lo = [-1.0_wp, -1.0_wp, -1.0_wp]
      cube_hi = [ 2.0_wp,  2.0_wp,  2.0_wp]
      call cube_cell_intersection(cell_lo, cell_hi, cube_lo, cube_hi, &
           V_full, V_eff, A_full, A_eff, n_frags, fn, fa, fc)
      call check('cell_in_cube: V_eff = 0', close_to(V_eff, 0.0_wp, TOL), fails)
      call check('cell_in_cube: all face areas = 0', &
         all([(close_to(A_eff(k),0.0_wp,TOL), k=1,6)]), fails)
      call check('cell_in_cube: n_frags = 0', n_frags == 0, fails)
   end subroutine test_cell_inside_cube

   ! ----------------------------------------------------------------
   ! Cube face EXACTLY at a cell face (grid-aligned obstacle). The
   ! cube and the cell share a 2D plane but the 3D intersection volume
   ! is zero. By design, cut-cell reports "no overlap" in this case:
   ! the flush cube face is handled by the existing obstacle-boundary-
   ! quad mechanism (cell-center cutout emits this face as a boundary
   ! quad tagged obstacle_faces, which the slip-wall BC dispatcher
   ! already handles). Cut-cell only kicks in for non-grid-aligned
   ! cubes whose faces cross strictly inside a cell — avoiding double-
   ! counting of the wall load.
   subroutine test_cube_face_flush_with_cell(fails)
      integer, intent(inout) :: fails
      real(wp) :: cell_lo(3), cell_hi(3), cube_lo(3), cube_hi(3)
      real(wp) :: V_full, V_eff, A_full(6), A_eff(6)
      integer  :: n_frags
      real(wp) :: fn(3, MAX_FRAGS_PER_CELL), fa(MAX_FRAGS_PER_CELL)
      real(wp) :: fc(3, MAX_FRAGS_PER_CELL)
      real(wp), parameter :: TOL = 1.0e-12_wp
      integer  :: k

      cell_lo = [0.0_wp, 0.0_wp, 0.0_wp]
      cell_hi = [1.0_wp, 1.0_wp, 1.0_wp]
      cube_lo = [1.0_wp, 0.0_wp, 0.0_wp]
      cube_hi = [2.0_wp, 1.0_wp, 1.0_wp]
      call cube_cell_intersection(cell_lo, cell_hi, cube_lo, cube_hi, &
           V_full, V_eff, A_full, A_eff, n_frags, fn, fa, fc)
      call check('flush: V_eff = V_full',    close_to(V_eff, V_full, TOL), fails)
      call check('flush: A_eff = A_full (handled by existing BC, not cut-cell)', &
         all([(close_to(A_eff(k), A_full(k), TOL), k=1,6)]), fails)
      call check('flush: n_frags = 0',       n_frags == 0, fails)
   end subroutine test_cube_face_flush_with_cell

end program test_cut_cell_geom
