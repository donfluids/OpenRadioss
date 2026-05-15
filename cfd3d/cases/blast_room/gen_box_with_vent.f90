! Standalone Gmsh msh2 ASCII generator: axis-aligned hex grid over
! [0,Lx]x[0,Ly]x[0,Lz] with a rectangular vent on the x_max face.
!
! All outer faces are tagged "walls" except quads on x=Lx whose center
! lies inside the vent rectangle [vy0, vy0+vdy] x [vz0, vz0+vdz] —
! those are tagged "vent_outside". This matches an internal-blast case
! with one opening on the +x face.
!
! Usage:
!   gen_box_with_vent Nx Ny Nz Lx Ly Lz vy0 vz0 vdy vdz out.msh
program gen_box_with_vent
   use, intrinsic :: iso_fortran_env, only : real64
   implicit none

   integer, parameter :: wp = real64

   integer :: Nx, Ny, Nz
   real(wp) :: Lx, Ly, Lz
   real(wp) :: vy0, vz0, vdy, vdz
   character(len=512) :: outfile
   character(len=64)  :: arg

   integer, parameter :: TAG_WALLS = 1, TAG_VENT = 2

   integer :: u, ios
   integer :: i, j, k, eid
   integer :: nvx, nvy, nvz, nv
   integer :: nc, nbf, n_vent
   integer :: v(8), q(4)
   real(wp) :: x, y, z, dx, dy, dz
   real(wp) :: yc, zc
   integer :: tag

   if (command_argument_count() < 11) then
      write(*,'(A)') 'usage: gen_box_with_vent Nx Ny Nz Lx Ly Lz vy0 vz0 vdy vdz out.msh'
      stop 1
   end if
   call get_command_argument(1,  arg); read(arg,*) Nx
   call get_command_argument(2,  arg); read(arg,*) Ny
   call get_command_argument(3,  arg); read(arg,*) Nz
   call get_command_argument(4,  arg); read(arg,*) Lx
   call get_command_argument(5,  arg); read(arg,*) Ly
   call get_command_argument(6,  arg); read(arg,*) Lz
   call get_command_argument(7,  arg); read(arg,*) vy0
   call get_command_argument(8,  arg); read(arg,*) vz0
   call get_command_argument(9,  arg); read(arg,*) vdy
   call get_command_argument(10, arg); read(arg,*) vdz
   call get_command_argument(11, outfile)

   nvx = Nx + 1
   nvy = Ny + 1
   nvz = Nz + 1
   nv  = nvx * nvy * nvz
   nc  = Nx * Ny * Nz
   nbf = 2*(Ny*Nz) + 2*(Nx*Nz) + 2*(Nx*Ny)
   dx  = Lx / real(Nx, wp)
   dy  = Ly / real(Ny, wp)
   dz  = Lz / real(Nz, wp)

   open(newunit=u, file=trim(outfile), status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      write(*,'(A,A)') 'gen_box_with_vent: cannot open ', trim(outfile)
      stop 1
   end if

   write(u,'(A)') '$MeshFormat'
   write(u,'(A)') '2.2 0 8'
   write(u,'(A)') '$EndMeshFormat'

   write(u,'(A)') '$PhysicalNames'
   write(u,'(I0)') 2
   write(u,'(I0,1X,I0,1X,A)') 2, TAG_WALLS, '"walls"'
   write(u,'(I0,1X,I0,1X,A)') 2, TAG_VENT,  '"vent_outside"'
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
            write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,8(1X,I0))') &
               eid, 5, 2, 0, 0, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8)
         end do
      end do
   end do

   ! Boundary quads — same six-face traversal as gen_sod_mesh, but the
   ! x_max face is split per-quad into walls vs. vent_outside.

   ! x_min — entirely walls
   do k = 1, Nz
      do j = 1, Ny
         eid = eid + 1
         q(1) = node_id(1, j,   k,   nvx, nvy)
         q(2) = node_id(1, j+1, k,   nvx, nvy)
         q(3) = node_id(1, j+1, k+1, nvx, nvy)
         q(4) = node_id(1, j,   k+1, nvx, nvy)
         call emit_quad(u, eid, TAG_WALLS, q)
      end do
   end do
   ! x_max — split into walls / vent_outside per-quad
   n_vent = 0
   do k = 1, Nz
      zc = (real(k, wp) - 0.5_wp) * dz
      do j = 1, Ny
         yc = (real(j, wp) - 0.5_wp) * dy
         if (yc >= vy0 .and. yc <= vy0 + vdy .and. &
             zc >= vz0 .and. zc <= vz0 + vdz) then
            tag = TAG_VENT
            n_vent = n_vent + 1
         else
            tag = TAG_WALLS
         end if
         eid = eid + 1
         q(1) = node_id(nvx, j,   k,   nvx, nvy)
         q(2) = node_id(nvx, j+1, k,   nvx, nvy)
         q(3) = node_id(nvx, j+1, k+1, nvx, nvy)
         q(4) = node_id(nvx, j,   k+1, nvx, nvy)
         call emit_quad(u, eid, tag, q)
      end do
   end do
   ! y_min / y_max / z_min / z_max — entirely walls
   do k = 1, Nz
      do i = 1, Nx
         eid = eid + 1
         q(1) = node_id(i,   1, k,   nvx, nvy)
         q(2) = node_id(i+1, 1, k,   nvx, nvy)
         q(3) = node_id(i+1, 1, k+1, nvx, nvy)
         q(4) = node_id(i,   1, k+1, nvx, nvy)
         call emit_quad(u, eid, TAG_WALLS, q)
      end do
   end do
   do k = 1, Nz
      do i = 1, Nx
         eid = eid + 1
         q(1) = node_id(i,   nvy, k,   nvx, nvy)
         q(2) = node_id(i+1, nvy, k,   nvx, nvy)
         q(3) = node_id(i+1, nvy, k+1, nvx, nvy)
         q(4) = node_id(i,   nvy, k+1, nvx, nvy)
         call emit_quad(u, eid, TAG_WALLS, q)
      end do
   end do
   do j = 1, Ny
      do i = 1, Nx
         eid = eid + 1
         q(1) = node_id(i,   j,   1, nvx, nvy)
         q(2) = node_id(i+1, j,   1, nvx, nvy)
         q(3) = node_id(i+1, j+1, 1, nvx, nvy)
         q(4) = node_id(i,   j+1, 1, nvx, nvy)
         call emit_quad(u, eid, TAG_WALLS, q)
      end do
   end do
   do j = 1, Ny
      do i = 1, Nx
         eid = eid + 1
         q(1) = node_id(i,   j,   nvz, nvx, nvy)
         q(2) = node_id(i+1, j,   nvz, nvx, nvy)
         q(3) = node_id(i+1, j+1, nvz, nvx, nvy)
         q(4) = node_id(i,   j+1, nvz, nvx, nvy)
         call emit_quad(u, eid, TAG_WALLS, q)
      end do
   end do

   write(u,'(A)') '$EndElements'
   close(u)

   write(*,'(A,I0,A,I0,A,I0,A)') 'gen_box_with_vent: wrote ', nc, &
      ' hex cells, ', nbf, ' boundary faces (', n_vent, ' vent quads on x_max)'

contains

   pure function node_id(i, j, k, nvx, nvy) result(id)
      integer, intent(in) :: i, j, k, nvx, nvy
      integer :: id
      id = (k - 1) * nvx * nvy + (j - 1) * nvx + i
   end function node_id

   subroutine emit_quad(u, eid, tag, q)
      integer, intent(in) :: u, eid, tag, q(4)
      write(u,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,4(1X,I0))') &
         eid, 3, 2, tag, tag, q(1), q(2), q(3), q(4)
   end subroutine emit_quad

end program gen_box_with_vent
