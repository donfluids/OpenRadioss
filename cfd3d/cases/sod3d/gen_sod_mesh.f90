! Standalone Gmsh msh2 ASCII generator for an axis-aligned hex grid.
! Usage: gen_sod_mesh <Nx> <Ny> <Nz> <Lx> <Ly> <Lz> <out.msh>
program gen_sod_mesh
   use, intrinsic :: iso_fortran_env, only : real64
   implicit none

   integer, parameter :: wp = real64

   integer :: Nx, Ny, Nz
   real(wp) :: Lx, Ly, Lz
   character(len=512) :: outfile
   character(len=64) :: arg

   integer :: u, ios
   integer :: i, j, k, eid
   integer :: nvx, nvy, nvz, nv
   integer :: nc, nbf
   integer :: v(8), q(4)
   real(wp) :: x, y, z, dx, dy, dz

   integer, parameter :: TAG_XMIN = 1, TAG_XMAX = 2, &
                         TAG_YMIN = 3, TAG_YMAX = 4, &
                         TAG_ZMIN = 5, TAG_ZMAX = 6

   if (command_argument_count() < 7) then
      write(*,'(A)') 'usage: gen_sod_mesh <Nx> <Ny> <Nz> <Lx> <Ly> <Lz> <out.msh>'
      stop 1
   end if
   call get_command_argument(1, arg); read(arg,*) Nx
   call get_command_argument(2, arg); read(arg,*) Ny
   call get_command_argument(3, arg); read(arg,*) Nz
   call get_command_argument(4, arg); read(arg,*) Lx
   call get_command_argument(5, arg); read(arg,*) Ly
   call get_command_argument(6, arg); read(arg,*) Lz
   call get_command_argument(7, outfile)

   nvx = Nx + 1
   nvy = Ny + 1
   nvz = Nz + 1
   nv  = nvx * nvy * nvz
   nc  = Nx * Ny * Nz
   nbf = 2*(Ny*Nz) + 2*(Nx*Nz) + 2*(Nx*Ny)
   dx = Lx / real(Nx, wp)
   dy = Ly / real(Ny, wp)
   dz = Lz / real(Nz, wp)

   open(newunit=u, file=trim(outfile), status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      write(*,'(A,A)') 'gen_sod_mesh: cannot open ', trim(outfile)
      stop 1
   end if

   write(u,'(A)') '$MeshFormat'
   write(u,'(A)') '2.2 0 8'
   write(u,'(A)') '$EndMeshFormat'

   write(u,'(A)') '$PhysicalNames'
   write(u,'(I0)') 6
   write(u,'(I0,1X,I0,1X,A)') 2, TAG_XMIN, '"x_min"'
   write(u,'(I0,1X,I0,1X,A)') 2, TAG_XMAX, '"x_max"'
   write(u,'(I0,1X,I0,1X,A)') 2, TAG_YMIN, '"y_min"'
   write(u,'(I0,1X,I0,1X,A)') 2, TAG_YMAX, '"y_max"'
   write(u,'(I0,1X,I0,1X,A)') 2, TAG_ZMIN, '"z_min"'
   write(u,'(I0,1X,I0,1X,A)') 2, TAG_ZMAX, '"z_max"'
   write(u,'(A)') '$EndPhysicalNames'

   write(u,'(A)') '$Nodes'
   write(u,'(I0)') nv
   do k = 1, nvz
      z = real(k-1, wp) * dz
      do j = 1, nvy
         y = real(j-1, wp) * dy
         do i = 1, nvx
            x = real(i-1, wp) * dx
            write(u,'(I0,3(1X,1PE22.14))') node_id(i,j,k,nvx,nvy), x, y, z
         end do
      end do
   end do
   write(u,'(A)') '$EndNodes'

   write(u,'(A)') '$Elements'
   write(u,'(I0)') nc + nbf
   eid = 0

   ! Hexahedra
   do k = 1, Nz
      do j = 1, Ny
         do i = 1, Nx
            eid = eid + 1
            v(1) = node_id(i,   j,   k,   nvx, nvy)
            v(2) = node_id(i+1, j,   k,   nvx, nvy)
            v(3) = node_id(i+1, j+1, k,   nvx, nvy)
            v(4) = node_id(i,   j+1, k,   nvx, nvy)
            v(5) = node_id(i,   j,   k+1, nvx, nvy)
            v(6) = node_id(i+1, j,   k+1, nvx, nvy)
            v(7) = node_id(i+1, j+1, k+1, nvx, nvy)
            v(8) = node_id(i,   j+1, k+1, nvx, nvy)
            ! type 5 = hex; 2 tags (phys, elem)
            write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,8(1X,I0))') &
               eid, 5, 2, 0, 0, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8)
         end do
      end do
   end do

   ! Boundary quads — emit per-face
   ! x_min: i=1
   do k = 1, Nz
      do j = 1, Ny
         eid = eid + 1
         q(1) = node_id(1, j,   k,   nvx, nvy)
         q(2) = node_id(1, j+1, k,   nvx, nvy)
         q(3) = node_id(1, j+1, k+1, nvx, nvy)
         q(4) = node_id(1, j,   k+1, nvx, nvy)
         write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,4(1X,I0))') &
            eid, 3, 2, TAG_XMIN, TAG_XMIN, q(1), q(2), q(3), q(4)
      end do
   end do
   ! x_max: i=Nx+1
   do k = 1, Nz
      do j = 1, Ny
         eid = eid + 1
         q(1) = node_id(nvx, j,   k,   nvx, nvy)
         q(2) = node_id(nvx, j+1, k,   nvx, nvy)
         q(3) = node_id(nvx, j+1, k+1, nvx, nvy)
         q(4) = node_id(nvx, j,   k+1, nvx, nvy)
         write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,4(1X,I0))') &
            eid, 3, 2, TAG_XMAX, TAG_XMAX, q(1), q(2), q(3), q(4)
      end do
   end do
   ! y_min: j=1
   do k = 1, Nz
      do i = 1, Nx
         eid = eid + 1
         q(1) = node_id(i,   1, k,   nvx, nvy)
         q(2) = node_id(i+1, 1, k,   nvx, nvy)
         q(3) = node_id(i+1, 1, k+1, nvx, nvy)
         q(4) = node_id(i,   1, k+1, nvx, nvy)
         write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,4(1X,I0))') &
            eid, 3, 2, TAG_YMIN, TAG_YMIN, q(1), q(2), q(3), q(4)
      end do
   end do
   ! y_max: j=Ny+1
   do k = 1, Nz
      do i = 1, Nx
         eid = eid + 1
         q(1) = node_id(i,   nvy, k,   nvx, nvy)
         q(2) = node_id(i+1, nvy, k,   nvx, nvy)
         q(3) = node_id(i+1, nvy, k+1, nvx, nvy)
         q(4) = node_id(i,   nvy, k+1, nvx, nvy)
         write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,4(1X,I0))') &
            eid, 3, 2, TAG_YMAX, TAG_YMAX, q(1), q(2), q(3), q(4)
      end do
   end do
   ! z_min: k=1
   do j = 1, Ny
      do i = 1, Nx
         eid = eid + 1
         q(1) = node_id(i,   j,   1, nvx, nvy)
         q(2) = node_id(i+1, j,   1, nvx, nvy)
         q(3) = node_id(i+1, j+1, 1, nvx, nvy)
         q(4) = node_id(i,   j+1, 1, nvx, nvy)
         write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,4(1X,I0))') &
            eid, 3, 2, TAG_ZMIN, TAG_ZMIN, q(1), q(2), q(3), q(4)
      end do
   end do
   ! z_max: k=Nz+1
   do j = 1, Ny
      do i = 1, Nx
         eid = eid + 1
         q(1) = node_id(i,   j,   nvz, nvx, nvy)
         q(2) = node_id(i+1, j,   nvz, nvx, nvy)
         q(3) = node_id(i+1, j+1, nvz, nvx, nvy)
         q(4) = node_id(i,   j+1, nvz, nvx, nvy)
         write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,4(1X,I0))') &
            eid, 3, 2, TAG_ZMAX, TAG_ZMAX, q(1), q(2), q(3), q(4)
      end do
   end do

   write(u,'(A)') '$EndElements'
   close(u)

contains

   pure function node_id(i, j, k, nvx, nvy) result(id)
      integer, intent(in) :: i, j, k, nvx, nvy
      integer :: id
      id = (k - 1) * nvx * nvy + (j - 1) * nvx + i
   end function node_id

end program gen_sod_mesh
