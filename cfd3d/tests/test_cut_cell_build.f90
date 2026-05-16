! Phase 1b unit test for cut_cell::build_cut_cell_tables.
!
! Constructs a tiny synthetic t_mesh consisting of a single axis-aligned
! unit-hex cell [0,1]^3 with its six boundary faces, then drives the
! cut-cell table builder against a list of axis-aligned cube obstacles
! and verifies the resulting cell_vol_eff, face_area_eff, and fragment
! CSR. No MPI, no Gmsh I/O — exercises the data-structure plumbing in
! isolation.
program test_cut_cell_build
   use kinds,         only : wp
   use mesh_types,    only : t_mesh, CELL_HEX
   use cut_cell_geom, only : MAX_FRAGS_PER_CELL
   use cut_cell,      only : build_cut_cell_tables, classify_face_axis_side, axis_side_idx
   implicit none

   integer :: fails

   fails = 0

   call test_no_obstacles            (fails)
   call test_single_cube_corner_cut  (fails)
   call test_grid_aligned_cube       (fails)

   if (fails == 0) then
      write(*,'(A)') 'test_cut_cell_build: PASS'
   else
      write(*,'(A,I0,A)') 'test_cut_cell_build: FAIL (', fails, ')'
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

   ! Build a synthetic unit-hex mesh: one [0,1]^3 cell, six axis-aligned
   ! boundary faces, no neighbor. All fields populated to the minimum
   ! needed by build_cut_cell_tables.
   subroutine build_unit_hex_mesh(mesh)
      type(t_mesh), intent(out) :: mesh
      integer :: k

      mesh%nv = 8
      allocate(mesh%xv(3, 8))
      mesh%xv(:,1) = [0.0_wp, 0.0_wp, 0.0_wp]
      mesh%xv(:,2) = [1.0_wp, 0.0_wp, 0.0_wp]
      mesh%xv(:,3) = [1.0_wp, 1.0_wp, 0.0_wp]
      mesh%xv(:,4) = [0.0_wp, 1.0_wp, 0.0_wp]
      mesh%xv(:,5) = [0.0_wp, 0.0_wp, 1.0_wp]
      mesh%xv(:,6) = [1.0_wp, 0.0_wp, 1.0_wp]
      mesh%xv(:,7) = [1.0_wp, 1.0_wp, 1.0_wp]
      mesh%xv(:,8) = [0.0_wp, 1.0_wp, 1.0_wp]

      mesh%nc_internal = 1
      mesh%nc_total    = 1
      allocate(mesh%cell_type(1));        mesh%cell_type(1) = CELL_HEX
      allocate(mesh%cell_vtx_ptr(2));     mesh%cell_vtx_ptr = [1, 9]
      allocate(mesh%cell_vtx(8));         mesh%cell_vtx = [(k, k=1,8)]
      allocate(mesh%cell_volume(1));      mesh%cell_volume(1)      = 1.0_wp
      allocate(mesh%cell_centroid(3, 1)); mesh%cell_centroid(:, 1) = 0.5_wp

      ! Six boundary faces, all owner = 1, no neighbor.
      mesh%nf               = 6
      mesh%nf_pure_interior = 0
      mesh%nf_interior      = 0
      mesh%nf_boundary      = 6
      allocate(mesh%face_owner   (6))
      allocate(mesh%face_owner_lf(6))
      allocate(mesh%face_neighbor(6))
      allocate(mesh%face_area    (6))
      allocate(mesh%face_normal  (3, 6))
      allocate(mesh%face_centroid(3, 6))
      allocate(mesh%face_patch   (6))
      mesh%face_owner    = 1
      mesh%face_owner_lf = [(k, k=1,6)]
      mesh%face_neighbor = 0
      mesh%face_area     = 1.0_wp
      mesh%face_patch    = 1
      mesh%face_normal(:,1) = [-1.0_wp, 0.0_wp, 0.0_wp]   ! -x
      mesh%face_normal(:,2) = [+1.0_wp, 0.0_wp, 0.0_wp]   ! +x
      mesh%face_normal(:,3) = [0.0_wp, -1.0_wp, 0.0_wp]   ! -y
      mesh%face_normal(:,4) = [0.0_wp, +1.0_wp, 0.0_wp]   ! +y
      mesh%face_normal(:,5) = [0.0_wp, 0.0_wp, -1.0_wp]   ! -z
      mesh%face_normal(:,6) = [0.0_wp, 0.0_wp, +1.0_wp]   ! +z
      mesh%face_centroid(:,1) = [0.0_wp, 0.5_wp, 0.5_wp]
      mesh%face_centroid(:,2) = [1.0_wp, 0.5_wp, 0.5_wp]
      mesh%face_centroid(:,3) = [0.5_wp, 0.0_wp, 0.5_wp]
      mesh%face_centroid(:,4) = [0.5_wp, 1.0_wp, 0.5_wp]
      mesh%face_centroid(:,5) = [0.5_wp, 0.5_wp, 0.0_wp]
      mesh%face_centroid(:,6) = [0.5_wp, 0.5_wp, 1.0_wp]
   end subroutine build_unit_hex_mesh

   subroutine free_mesh_local(mesh)
      type(t_mesh), intent(inout) :: mesh
      if (allocated(mesh%xv))               deallocate(mesh%xv)
      if (allocated(mesh%cell_type))        deallocate(mesh%cell_type)
      if (allocated(mesh%cell_vtx_ptr))     deallocate(mesh%cell_vtx_ptr)
      if (allocated(mesh%cell_vtx))         deallocate(mesh%cell_vtx)
      if (allocated(mesh%cell_volume))      deallocate(mesh%cell_volume)
      if (allocated(mesh%cell_centroid))    deallocate(mesh%cell_centroid)
      if (allocated(mesh%face_owner))       deallocate(mesh%face_owner)
      if (allocated(mesh%face_owner_lf))    deallocate(mesh%face_owner_lf)
      if (allocated(mesh%face_neighbor))    deallocate(mesh%face_neighbor)
      if (allocated(mesh%face_area))        deallocate(mesh%face_area)
      if (allocated(mesh%face_normal))      deallocate(mesh%face_normal)
      if (allocated(mesh%face_centroid))    deallocate(mesh%face_centroid)
      if (allocated(mesh%face_patch))       deallocate(mesh%face_patch)
      if (allocated(mesh%cell_vol_eff))     deallocate(mesh%cell_vol_eff)
      if (allocated(mesh%face_area_eff))    deallocate(mesh%face_area_eff)
      if (allocated(mesh%cell_frag_offset)) deallocate(mesh%cell_frag_offset)
      if (allocated(mesh%frag_normal))      deallocate(mesh%frag_normal)
      if (allocated(mesh%frag_area))        deallocate(mesh%frag_area)
      if (allocated(mesh%frag_centroid))    deallocate(mesh%frag_centroid)
   end subroutine free_mesh_local

   ! ----------------------------------------------------------------
   ! n_cubes = 0 → tables defaulted to un-cut values, no fragments.
   subroutine test_no_obstacles(fails)
      integer, intent(inout) :: fails
      type(t_mesh) :: mesh
      real(wp) :: cube_lo(3,1), cube_hi(3,1)
      real(wp), parameter :: TOL = 1.0e-12_wp
      integer :: k

      call build_unit_hex_mesh(mesh)
      cube_lo = 0.0_wp;  cube_hi = 0.0_wp
      call build_cut_cell_tables(mesh, 0, cube_lo, cube_hi)

      call check('no_obs: cell_vol_eff = cell_volume', &
         close_to(mesh%cell_vol_eff(1), 1.0_wp, TOL), fails)
      call check('no_obs: all face_area_eff = 1', &
         all([(close_to(mesh%face_area_eff(k), 1.0_wp, TOL), k=1,6)]), fails)
      call check('no_obs: n_frags = 0', mesh%n_frags == 0, fails)
      call check('no_obs: cell_frag_offset = [1,1]', &
         mesh%cell_frag_offset(1) == 1 .and. mesh%cell_frag_offset(2) == 1, fails)

      call free_mesh_local(mesh)
   end subroutine test_no_obstacles

   ! ----------------------------------------------------------------
   ! One cube cutting the +x/+y/+z corner of the unit cell. Removes
   ! 0.5^3 = 0.125 volume; 3 fragments (each 0.25 m^2). Three faces
   ! see their area halved by an in-cube quarter, the other three
   ! see a 1/4 quadrant covered (A_eff = 0.75).
   subroutine test_single_cube_corner_cut(fails)
      integer, intent(inout) :: fails
      type(t_mesh) :: mesh
      real(wp) :: cube_lo(3,1), cube_hi(3,1)
      real(wp), parameter :: TOL = 1.0e-12_wp
      integer  :: k
      real(wp) :: total_frag_area

      call build_unit_hex_mesh(mesh)
      cube_lo(:,1) = [0.5_wp, 0.5_wp, 0.5_wp]
      cube_hi(:,1) = [2.0_wp, 2.0_wp, 2.0_wp]
      call build_cut_cell_tables(mesh, 1, cube_lo, cube_hi)

      call check('corner: cell_vol_eff = 0.875', &
         close_to(mesh%cell_vol_eff(1), 0.875_wp, TOL), fails)
      ! The -x face at x=0 is OUTSIDE the cube's x-extent [0.5, 2]
      ! (face-coverage test cube_lo(1) <= 0 <= cube_hi(1) is false),
      ! so A_eff(-x) = A_full = 1. Same for -y, -z. The +x face at
      ! x=1 IS inside the cube's x-extent, with a 0.5x0.5 quadrant
      ! covered, so A_eff(+x) = 1 - 0.25 = 0.75. Same for +y, +z.
      call check('corner: -x face A_eff = 1.0',  close_to(mesh%face_area_eff(1), 1.0_wp,  TOL), fails)
      call check('corner: +x face A_eff = 0.75', close_to(mesh%face_area_eff(2), 0.75_wp, TOL), fails)
      call check('corner: -y face A_eff = 1.0',  close_to(mesh%face_area_eff(3), 1.0_wp,  TOL), fails)
      call check('corner: +y face A_eff = 0.75', close_to(mesh%face_area_eff(4), 0.75_wp, TOL), fails)
      call check('corner: -z face A_eff = 1.0',  close_to(mesh%face_area_eff(5), 1.0_wp,  TOL), fails)
      call check('corner: +z face A_eff = 0.75', close_to(mesh%face_area_eff(6), 0.75_wp, TOL), fails)

      call check('corner: n_frags = 3', mesh%n_frags == 3, fails)
      call check('corner: CSR [1,4]', &
         mesh%cell_frag_offset(1) == 1 .and. mesh%cell_frag_offset(2) == 4, fails)
      total_frag_area = 0.0_wp
      do k = 1, 3
         total_frag_area = total_frag_area + mesh%frag_area(k)
      end do
      call check('corner: total frag area = 0.75', close_to(total_frag_area, 0.75_wp, TOL), fails)

      call free_mesh_local(mesh)
   end subroutine test_single_cube_corner_cut

   ! ----------------------------------------------------------------
   ! Cube flush with cell +x face (grid-aligned, no 3D overlap).
   ! Tables should match the no-obstacle case — the existing
   ! boundary-quad mechanism handles the flush face.
   subroutine test_grid_aligned_cube(fails)
      integer, intent(inout) :: fails
      type(t_mesh) :: mesh
      real(wp) :: cube_lo(3,1), cube_hi(3,1)
      real(wp), parameter :: TOL = 1.0e-12_wp
      integer  :: k

      call build_unit_hex_mesh(mesh)
      cube_lo(:,1) = [1.0_wp, 0.0_wp, 0.0_wp]
      cube_hi(:,1) = [2.0_wp, 1.0_wp, 1.0_wp]
      call build_cut_cell_tables(mesh, 1, cube_lo, cube_hi)

      call check('flush: cell_vol_eff = 1.0',   close_to(mesh%cell_vol_eff(1), 1.0_wp, TOL), fails)
      call check('flush: all A_eff = 1.0',      &
         all([(close_to(mesh%face_area_eff(k), 1.0_wp, TOL), k=1,6)]), fails)
      call check('flush: n_frags = 0',          mesh%n_frags == 0, fails)

      call free_mesh_local(mesh)
   end subroutine test_grid_aligned_cube

end program test_cut_cell_build
