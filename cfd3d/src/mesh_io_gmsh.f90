module mesh_io_gmsh
   use kinds,       only : wp
   use mesh_types,  only : t_mesh, t_patch, max_vpf, &
                           CELL_TET, CELL_HEX, CELL_PRISM, CELL_PYRAMID, &
                           GMSH_TYPE_TRI, GMSH_TYPE_QUAD, &
                           GMSH_TYPE_TET, GMSH_TYPE_HEX, GMSH_TYPE_PRISM, GMSH_TYPE_PYRAMID, &
                           n_vtx_per_cell
   implicit none
   private

   public :: t_bface_raw, gmsh_read

   type :: t_bface_raw
      integer :: nvf = 0                       ! 3 (tri) or 4 (quad)
      integer :: vtx(max_vpf) = 0              ! 1-based vertex ids (local)
      integer :: gmsh_tag = 0                  ! physical group tag
   end type t_bface_raw

contains

   subroutine gmsh_read(filename, mesh, bfaces, nb)
      character(len=*),  intent(in)    :: filename
      type(t_mesh),      intent(inout) :: mesh
      type(t_bface_raw), allocatable, intent(out) :: bfaces(:)
      integer,           intent(out)   :: nb

      integer :: u, ios
      character(len=256) :: line

      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(A,A)') 'gmsh_read: cannot open ', trim(filename)
         error stop 1
      end if

      nb = 0
      mesh%nv = 0
      mesh%nc_internal = 0
      mesh%nc_total    = 0

      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         select case (trim(line))
         case ('$MeshFormat')
            call skip_to_end(u, '$EndMeshFormat')
         case ('$PhysicalNames')
            call read_physical_names(u, mesh)
         case ('$Nodes')
            call read_nodes(u, mesh)
         case ('$Elements')
            call read_elements(u, mesh, bfaces, nb)
         case default
            ! Unknown / future section — skip silently
            cycle
         end select
      end do

      close(u)

      if (mesh%nv == 0) then
         write(*,'(A)') 'gmsh_read: no nodes found'
         error stop 1
      end if
      if (mesh%nc_internal == 0) then
         write(*,'(A)') 'gmsh_read: no 3D cells found'
         error stop 1
      end if
   end subroutine gmsh_read

   subroutine skip_to_end(u, endtag)
      integer,          intent(in) :: u
      character(len=*), intent(in) :: endtag
      character(len=256) :: line
      integer :: ios
      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) return
         if (trim(adjustl(line)) == endtag) return
      end do
   end subroutine skip_to_end

   subroutine read_physical_names(u, mesh)
      integer,       intent(in)    :: u
      type(t_mesh),  intent(inout) :: mesh
      integer :: n, i, dim, tag, ios, q1, q2
      character(len=256) :: line, nm

      read(u, *) n
      if (allocated(mesh%patches)) deallocate(mesh%patches)
      allocate(mesh%patches(n))
      mesh%np = n
      do i = 1, n
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         ! Format: dim tag "name"
         read(line, *) dim, tag
         q1 = index(line, '"')
         q2 = index(line(q1+1:), '"')
         if (q1 > 0 .and. q2 > 0) then
            nm = line(q1+1 : q1+q2-1)
         else
            nm = ''
         end if
         mesh%patches(i)%name     = trim(nm)
         mesh%patches(i)%gmsh_tag = tag
      end do
      call skip_to_end(u, '$EndPhysicalNames')
   end subroutine read_physical_names

   subroutine read_nodes(u, mesh)
      integer,       intent(in)    :: u
      type(t_mesh),  intent(inout) :: mesh
      integer :: n, i, id
      real(wp) :: x, y, z

      read(u, *) n
      mesh%nv = n
      if (allocated(mesh%xv)) deallocate(mesh%xv)
      allocate(mesh%xv(3, n))
      do i = 1, n
         read(u, *) id, x, y, z
         ! msh2 node ids are usually 1..n consecutive; trust id when contiguous
         mesh%xv(1, id) = x
         mesh%xv(2, id) = y
         mesh%xv(3, id) = z
      end do
      call skip_to_end(u, '$EndNodes')
   end subroutine read_nodes

   subroutine read_elements(u, mesh, bfaces, nb)
      integer,           intent(in)    :: u
      type(t_mesh),      intent(inout) :: mesh
      type(t_bface_raw), allocatable, intent(inout) :: bfaces(:)
      integer,           intent(inout) :: nb

      integer :: nelem, i, id, etype, ntags, t, k
      integer :: phys_tag, n_cells, n_bf, csr_size
      integer :: tags(16)
      integer :: nodes(32)
      integer, allocatable :: cell_type(:), cell_vtx_ptr(:), cell_vtx(:)
      type(t_bface_raw), allocatable :: tmp_bf(:)
      integer :: cap_bf
      integer :: nv_cell, ct, npos
      integer, allocatable :: line_buf(:)

      read(u, *) nelem

      ! First pass would be cleaner but we read once with growable arrays.
      n_cells = 0
      n_bf    = 0
      csr_size = 0
      cap_bf = 16
      allocate(tmp_bf(cap_bf))
      ! Worst-case bounds (we'll trim at end):
      allocate(cell_type(nelem))
      allocate(cell_vtx_ptr(nelem+1))
      allocate(cell_vtx(8*nelem))  ! 8 = max vpc (hex); resize if exceeded
      cell_vtx_ptr(1) = 1

      do i = 1, nelem
         read(u, *) id, etype, ntags
         ! Re-read with tags + nodes
         ! Easier: read full line and parse
         backspace(u)
         call read_int_line(u, line_buf)
         id    = line_buf(1)
         etype = line_buf(2)
         ntags = line_buf(3)
         do t = 1, ntags
            tags(t) = line_buf(3 + t)
         end do
         phys_tag = tags(1)
         npos = 3 + ntags

         select case (etype)
         case (GMSH_TYPE_TET)
            ct = CELL_TET
            nv_cell = 4
         case (GMSH_TYPE_HEX)
            ct = CELL_HEX
            nv_cell = 8
         case (GMSH_TYPE_PRISM)
            ct = CELL_PRISM
            nv_cell = 6
         case (GMSH_TYPE_PYRAMID)
            ct = CELL_PYRAMID
            nv_cell = 5
         case (GMSH_TYPE_TRI)
            ct = 0
            nv_cell = 3
         case (GMSH_TYPE_QUAD)
            ct = 0
            nv_cell = 4
         case default
            cycle
         end select

         do k = 1, nv_cell
            nodes(k) = line_buf(npos + k)
         end do

         if (ct > 0) then
            ! 3D cell
            n_cells = n_cells + 1
            cell_type(n_cells) = ct
            ! Ensure CSR capacity
            if (csr_size + nv_cell > size(cell_vtx)) then
               call grow_int(cell_vtx, max(size(cell_vtx)*2, csr_size + nv_cell))
            end if
            do k = 1, nv_cell
               cell_vtx(csr_size + k) = nodes(k)
            end do
            csr_size = csr_size + nv_cell
            cell_vtx_ptr(n_cells + 1) = csr_size + 1
         else
            ! 2D boundary element
            n_bf = n_bf + 1
            if (n_bf > cap_bf) then
               call grow_bf(tmp_bf, cap_bf * 2)
               cap_bf = size(tmp_bf)
            end if
            tmp_bf(n_bf)%nvf = nv_cell
            tmp_bf(n_bf)%vtx = 0
            do k = 1, nv_cell
               tmp_bf(n_bf)%vtx(k) = nodes(k)
            end do
            tmp_bf(n_bf)%gmsh_tag = phys_tag
         end if
      end do

      call skip_to_end(u, '$EndElements')

      ! Commit to mesh
      mesh%nc_internal = n_cells
      mesh%nc_total    = n_cells
      allocate(mesh%cell_type(n_cells))
      allocate(mesh%cell_vtx_ptr(n_cells + 1))
      allocate(mesh%cell_vtx(csr_size))
      mesh%cell_type   = cell_type(1:n_cells)
      mesh%cell_vtx_ptr(1:n_cells+1) = cell_vtx_ptr(1:n_cells+1)
      mesh%cell_vtx    = cell_vtx(1:csr_size)

      nb = n_bf
      if (allocated(bfaces)) deallocate(bfaces)
      allocate(bfaces(n_bf))
      do i = 1, n_bf
         bfaces(i) = tmp_bf(i)
      end do
   end subroutine read_elements

   ! Read a whitespace-separated line of integers into a growable buffer
   subroutine read_int_line(u, buf)
      integer,              intent(in)    :: u
      integer, allocatable, intent(inout) :: buf(:)
      character(len=4096) :: line
      integer :: ios, n, k, pos, length
      logical :: in_tok
      character :: ch

      read(u, '(A)', iostat=ios) line
      if (ios /= 0) then
         if (allocated(buf)) deallocate(buf)
         allocate(buf(0))
         return
      end if

      ! Count tokens
      n = 0
      in_tok = .false.
      length = len_trim(line)
      do pos = 1, length
         ch = line(pos:pos)
         if (ch == ' ' .or. ch == char(9)) then
            in_tok = .false.
         else
            if (.not. in_tok) then
               n = n + 1
               in_tok = .true.
            end if
         end if
      end do

      if (allocated(buf)) deallocate(buf)
      allocate(buf(n))
      if (n > 0) read(line, *) (buf(k), k = 1, n)
   end subroutine read_int_line

   subroutine grow_int(arr, new_size)
      integer, allocatable, intent(inout) :: arr(:)
      integer,              intent(in)    :: new_size
      integer, allocatable :: tmp(:)
      integer :: n_old
      n_old = size(arr)
      allocate(tmp(new_size))
      tmp(1:n_old) = arr
      call move_alloc(tmp, arr)
   end subroutine grow_int

   subroutine grow_bf(arr, new_size)
      type(t_bface_raw), allocatable, intent(inout) :: arr(:)
      integer,                        intent(in)    :: new_size
      type(t_bface_raw), allocatable :: tmp(:)
      integer :: n_old
      n_old = size(arr)
      allocate(tmp(new_size))
      tmp(1:n_old) = arr
      call move_alloc(tmp, arr)
   end subroutine grow_bf

end module mesh_io_gmsh
