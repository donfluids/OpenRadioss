program test_mesh_topology
   use, intrinsic :: iso_fortran_env, only : error_unit
   use kinds,        only : wp
   use mesh_types,   only : t_mesh
   use mesh_module,  only : load_mesh
   implicit none

   type(t_mesh) :: mesh
   character(len=512) :: meshfile
   integer :: nargs, f, c, fails
   real(wp) :: vol_sum, dot
   real(wp), parameter :: TOL = 1.0e-10_wp

   nargs = command_argument_count()
   if (nargs < 1) then
      write(error_unit,'(A)') 'usage: test_mesh_topology <mesh.msh>'
      stop 1
   end if
   call get_command_argument(1, meshfile)

   call load_mesh(trim(meshfile), mesh)

   fails = 0

   ! For an Nx=2, Ny=1, Nz=1 hex mesh:
   !   2 cells, 1 interior face, 10 boundary faces (= 11 total).
   if (mesh%nc_internal /= 2) then
      write(error_unit,*) 'expected 2 cells, got', mesh%nc_internal; fails = fails + 1
   end if
   if (mesh%nf_interior /= 1) then
      write(error_unit,*) 'expected 1 interior face, got', mesh%nf_interior; fails = fails + 1
   end if
   if (mesh%nf_boundary /= 10) then
      write(error_unit,*) 'expected 10 boundary faces, got', mesh%nf_boundary; fails = fails + 1
   end if
   if (mesh%nf /= 11) then
      write(error_unit,*) 'expected 11 total faces, got', mesh%nf; fails = fails + 1
   end if

   ! Interior face: owner < neighbor
   do f = 1, mesh%nf_interior
      if (mesh%face_owner(f) >= mesh%face_neighbor(f)) then
         write(error_unit,*) 'interior face', f, ': owner >= neighbor', &
            mesh%face_owner(f), mesh%face_neighbor(f)
         fails = fails + 1
      end if
   end do

   ! Boundary face: neighbor == 0, patch > 0
   do f = mesh%nf_interior + 1, mesh%nf
      if (mesh%face_neighbor(f) /= 0) then
         write(error_unit,*) 'boundary face', f, 'neighbor /= 0', mesh%face_neighbor(f)
         fails = fails + 1
      end if
      if (mesh%face_patch(f) <= 0) then
         write(error_unit,*) 'boundary face', f, 'patch <= 0', mesh%face_patch(f)
         fails = fails + 1
      end if
   end do

   ! Normal orientation: for every face, (face_centroid - owner_centroid)·normal > 0
   do f = 1, mesh%nf
      c = mesh%face_owner(f)
      dot = (mesh%face_centroid(1,f) - mesh%cell_centroid(1,c)) * mesh%face_normal(1,f) &
          + (mesh%face_centroid(2,f) - mesh%cell_centroid(2,c)) * mesh%face_normal(2,f) &
          + (mesh%face_centroid(3,f) - mesh%cell_centroid(3,c)) * mesh%face_normal(3,f)
      if (dot <= 0.0_wp) then
         write(error_unit,*) 'face', f, ' normal not pointing owner-outward, dot=', dot
         fails = fails + 1
      end if
   end do

   ! Volume sum equals domain volume
   vol_sum = sum(mesh%cell_volume(1:mesh%nc_internal))
   ! For Nx=2, Ny=1, Nz=1, Lx=2, Ly=1, Lz=1, total volume = 2.
   if (abs(vol_sum - 2.0_wp) > TOL) then
      write(error_unit,'(A,1PE15.6)') 'volume sum mismatch (expected 2.0):', vol_sum
      fails = fails + 1
   end if

   ! All patches have face_count > 0 (6 patches × 1 face each for this mesh)
   if (mesh%np /= 6) then
      write(error_unit,*) 'expected 6 patches, got', mesh%np
      fails = fails + 1
   end if

   if (fails == 0) then
      write(*,'(A)') 'test_mesh_topology: PASS'
   else
      write(*,'(A,I0,A)') 'test_mesh_topology: FAIL (', fails, ' issues)'
      stop 1
   end if
end program test_mesh_topology
