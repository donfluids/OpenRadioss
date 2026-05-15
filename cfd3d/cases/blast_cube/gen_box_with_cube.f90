! Standalone Gmsh msh2 ASCII generator: axis-aligned hex grid over
! [0,Lx]x[0,Ly]x[0,Lz] with an axis-aligned cube cut out.
!
! Cells whose centers lie inside the cube are omitted. Faces between an
! air cell and the cube interior are tagged "cube_faces" (slip-wall).
! The six outer-domain faces are tagged "farfield". Cube edges should
! align with cell edges — the caller picks Nx,Ny,Nz so that hx/dx etc.
! are integers; off-grid cubes still work but produce a staircase.
!
! Usage:
!   gen_box_with_cube Nx Ny Nz Lx Ly Lz cx cy cz hx hy hz out.msh
program gen_box_with_cube
   use, intrinsic :: iso_fortran_env, only : real64
   implicit none

   integer, parameter :: wp = real64

   integer :: Nx, Ny, Nz
   real(wp) :: Lx, Ly, Lz
   real(wp) :: cx, cy, cz, hx, hy, hz
   character(len=512) :: outfile
   character(len=64)  :: arg

   integer, parameter :: TAG_FARFIELD = 1, TAG_CUBE = 2

   integer :: u, ios
   integer :: i, j, k, eid
   integer :: nvx, nvy, nvz, nv
   integer :: nc_air, nbf, ne_total
   integer :: q(4)
   real(wp) :: xcen, ycen, zcen, dx, dy, dz
   logical, allocatable :: is_air(:,:,:)
   logical :: L_air, R_air
   integer :: face_tag

   if (command_argument_count() < 13) then
      write(*,'(A)') 'usage: gen_box_with_cube Nx Ny Nz Lx Ly Lz cx cy cz hx hy hz out.msh'
      stop 1
   end if
   call get_command_argument(1,  arg); read(arg,*) Nx
   call get_command_argument(2,  arg); read(arg,*) Ny
   call get_command_argument(3,  arg); read(arg,*) Nz
   call get_command_argument(4,  arg); read(arg,*) Lx
   call get_command_argument(5,  arg); read(arg,*) Ly
   call get_command_argument(6,  arg); read(arg,*) Lz
   call get_command_argument(7,  arg); read(arg,*) cx
   call get_command_argument(8,  arg); read(arg,*) cy
   call get_command_argument(9,  arg); read(arg,*) cz
   call get_command_argument(10, arg); read(arg,*) hx
   call get_command_argument(11, arg); read(arg,*) hy
   call get_command_argument(12, arg); read(arg,*) hz
   call get_command_argument(13, outfile)

   nvx = Nx + 1
   nvy = Ny + 1
   nvz = Nz + 1
   nv  = nvx * nvy * nvz
   dx  = Lx / real(Nx, wp)
   dy  = Ly / real(Ny, wp)
   dz  = Lz / real(Nz, wp)

   ! Per-cell air mask
   allocate(is_air(Nx, Ny, Nz))
   nc_air = 0
   do k = 1, Nz
      zcen = (real(k, wp) - 0.5_wp) * dz
      do j = 1, Ny
         ycen = (real(j, wp) - 0.5_wp) * dy
         do i = 1, Nx
            xcen = (real(i, wp) - 0.5_wp) * dx
            if (xcen >= cx .and. xcen <= cx + hx .and. &
                ycen >= cy .and. ycen <= cy + hy .and. &
                zcen >= cz .and. zcen <= cz + hz) then
               is_air(i,j,k) = .false.
            else
               is_air(i,j,k) = .true.
               nc_air = nc_air + 1
            end if
         end do
      end do
   end do

   ! Count boundary faces: a face is a boundary iff exactly one side is air.
   ! Cells outside [1..Nx]x[1..Ny]x[1..Nz] are treated as non-air (out of domain);
   ! such faces are outer-boundary -> farfield. Air-vs-void within domain -> cube.
   nbf = 0
   ! x-faces at node plane i_node in 1..nvx
   do k = 1, Nz
      do j = 1, Ny
         do i = 1, nvx
            L_air = cell_is_air(i-1, j, k)
            R_air = cell_is_air(i,   j, k)
            if (L_air .neqv. R_air) nbf = nbf + 1
         end do
      end do
   end do
   ! y-faces
   do k = 1, Nz
      do j = 1, nvy
         do i = 1, Nx
            L_air = cell_is_air(i, j-1, k)
            R_air = cell_is_air(i, j,   k)
            if (L_air .neqv. R_air) nbf = nbf + 1
         end do
      end do
   end do
   ! z-faces
   do k = 1, nvz
      do j = 1, Ny
         do i = 1, Nx
            L_air = cell_is_air(i, j, k-1)
            R_air = cell_is_air(i, j, k  )
            if (L_air .neqv. R_air) nbf = nbf + 1
         end do
      end do
   end do

   ne_total = nc_air + nbf

   open(newunit=u, file=trim(outfile), status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      write(*,'(A,A)') 'gen_box_with_cube: cannot open ', trim(outfile)
      stop 1
   end if

   write(u,'(A)') '$MeshFormat'
   write(u,'(A)') '2.2 0 8'
   write(u,'(A)') '$EndMeshFormat'

   write(u,'(A)') '$PhysicalNames'
   write(u,'(I0)') 2
   write(u,'(I0,1X,I0,1X,A)') 2, TAG_FARFIELD, '"farfield"'
   write(u,'(I0,1X,I0,1X,A)') 2, TAG_CUBE,     '"cube_faces"'
   write(u,'(A)') '$EndPhysicalNames'

   write(u,'(A)') '$Nodes'
   write(u,'(I0)') nv
   do k = 1, nvz
      do j = 1, nvy
         do i = 1, nvx
            write(u,'(I0,3(1X,1PE22.14))') node_id(i,j,k,nvx,nvy), &
               real(i-1, wp)*dx, real(j-1, wp)*dy, real(k-1, wp)*dz
         end do
      end do
   end do
   write(u,'(A)') '$EndNodes'

   write(u,'(A)') '$Elements'
   write(u,'(I0)') ne_total
   eid = 0

   ! Hexahedra for air cells
   do k = 1, Nz
      do j = 1, Ny
         do i = 1, Nx
            if (is_air(i,j,k)) then
               eid = eid + 1
               write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,8(1X,I0))') &
                  eid, 5, 2, 0, 0, &
                  node_id(i,   j,   k,   nvx, nvy), &
                  node_id(i+1, j,   k,   nvx, nvy), &
                  node_id(i+1, j+1, k,   nvx, nvy), &
                  node_id(i,   j+1, k,   nvx, nvy), &
                  node_id(i,   j,   k+1, nvx, nvy), &
                  node_id(i+1, j,   k+1, nvx, nvy), &
                  node_id(i+1, j+1, k+1, nvx, nvy), &
                  node_id(i,   j+1, k+1, nvx, nvy)
            end if
         end do
      end do
   end do

   ! Boundary quads: x-faces
   do k = 1, Nz
      do j = 1, Ny
         do i = 1, nvx
            L_air = cell_is_air(i-1, j, k)
            R_air = cell_is_air(i,   j, k)
            if (L_air .neqv. R_air) then
               face_tag = pick_tag(i == 1 .or. i == nvx)
               q(1) = node_id(i, j,   k,   nvx, nvy)
               q(2) = node_id(i, j+1, k,   nvx, nvy)
               q(3) = node_id(i, j+1, k+1, nvx, nvy)
               q(4) = node_id(i, j,   k+1, nvx, nvy)
               eid = eid + 1
               write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,4(1X,I0))') &
                  eid, 3, 2, face_tag, face_tag, q(1), q(2), q(3), q(4)
            end if
         end do
      end do
   end do
   ! y-faces
   do k = 1, Nz
      do j = 1, nvy
         do i = 1, Nx
            L_air = cell_is_air(i, j-1, k)
            R_air = cell_is_air(i, j,   k)
            if (L_air .neqv. R_air) then
               face_tag = pick_tag(j == 1 .or. j == nvy)
               q(1) = node_id(i,   j, k,   nvx, nvy)
               q(2) = node_id(i+1, j, k,   nvx, nvy)
               q(3) = node_id(i+1, j, k+1, nvx, nvy)
               q(4) = node_id(i,   j, k+1, nvx, nvy)
               eid = eid + 1
               write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,4(1X,I0))') &
                  eid, 3, 2, face_tag, face_tag, q(1), q(2), q(3), q(4)
            end if
         end do
      end do
   end do
   ! z-faces
   do k = 1, nvz
      do j = 1, Ny
         do i = 1, Nx
            L_air = cell_is_air(i, j, k-1)
            R_air = cell_is_air(i, j, k  )
            if (L_air .neqv. R_air) then
               face_tag = pick_tag(k == 1 .or. k == nvz)
               q(1) = node_id(i,   j,   k, nvx, nvy)
               q(2) = node_id(i+1, j,   k, nvx, nvy)
               q(3) = node_id(i+1, j+1, k, nvx, nvy)
               q(4) = node_id(i,   j+1, k, nvx, nvy)
               eid = eid + 1
               write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,4(1X,I0))') &
                  eid, 3, 2, face_tag, face_tag, q(1), q(2), q(3), q(4)
            end if
         end do
      end do
   end do

   write(u,'(A)') '$EndElements'
   close(u)

   write(*,'(A,I0,A,I0,A,I0,A)') 'gen_box_with_cube: wrote ', nc_air, &
        ' air cells, ', nbf, ' boundary faces (', ne_total, ' total elements)'

contains

   pure function node_id(i, j, k, nvx, nvy) result(id)
      integer, intent(in) :: i, j, k, nvx, nvy
      integer :: id
      id = (k - 1) * nvx * nvy + (j - 1) * nvx + i
   end function node_id

   pure logical function cell_is_air(ii, jj, kk) result(yes)
      integer, intent(in) :: ii, jj, kk
      if (ii < 1 .or. ii > Nx .or. jj < 1 .or. jj > Ny .or. kk < 1 .or. kk > Nz) then
         yes = .false.
      else
         yes = is_air(ii, jj, kk)
      end if
   end function cell_is_air

   pure integer function pick_tag(is_outer) result(t)
      logical, intent(in) :: is_outer
      if (is_outer) then
         t = TAG_FARFIELD
      else
         t = TAG_CUBE
      end if
   end function pick_tag

end program gen_box_with_cube
