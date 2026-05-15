! Standalone Gmsh msh2 ASCII generator for an arbitrary-shaped room.
!
! The room is defined as the intersection of:
!   - a polygonal footprint in (x, y), specified by N vertices in order,
!     tested via ray-casting (point-in-polygon)
!   - a vertical extrusion range [ROOM_Z z0 z1]
! plus a list of axis-aligned cube obstacles inside the room.
!
! Cells are classified by center-inclusion:
!     AIR           — inside polygon AND inside [z0,z1] AND not in any obstacle
!     VOID_OBSTACLE — inside any obstacle
!     VOID_OUTSIDE  — outside polygon or outside [z0,z1]
!
! Hex elements are emitted only for AIR cells. Boundary quads are
! emitted for every face with exactly one AIR neighbour; the tag is:
!     "walls"          — the other side is VOID_OUTSIDE (interior room
!                        wall) OR the outer bbox face is in a wall area
!     "vent_outside"   — the outer bbox face falls inside a vent rectangle
!     "obstacle_faces" — the other side is VOID_OBSTACLE
!
! Spec file format (whitespace-separated tokens, '#' starts a comment):
!     GRID    Nx Ny Nz
!     BBOX    Lx Ly Lz                 ! origin at (0,0,0)
!     ROOM_Z  z0 z1
!     POLYGON N
!     VERT    x1 y1
!     VERT    x2 y2
!     ...     (N VERT lines)
!     VENT    face u0 v0 u1 v1         ! face: 0=x_min 1=x_max 2=y_min
!                                      !       3=y_max 4=z_min 5=z_max
!                                      ! (u,v) are the in-plane coords
!                                      ! on that face, in (x,y,z) order
!                                      ! (multiple VENT lines allowed)
!     CUBE    x0 y0 z0 x1 y1 z1        ! optional, multiple allowed
!     END                              ! optional
!
! Usage:
!   gen_room_mesh <spec.txt> <out.msh>
program gen_room_mesh
   use, intrinsic :: iso_fortran_env, only : real64, error_unit, iostat_end
   implicit none

   integer, parameter :: wp = real64

   integer, parameter :: TAG_WALLS    = 1
   integer, parameter :: TAG_VENT     = 2
   integer, parameter :: TAG_OBSTACLE = 3

   integer, parameter :: KIND_AIR           = 1
   integer, parameter :: KIND_VOID_OUTSIDE  = 2
   integer, parameter :: KIND_VOID_OBSTACLE = 3
   integer, parameter :: KIND_OFF_DOMAIN    = 4

   integer, parameter :: MAX_POLY      = 64
   integer, parameter :: MAX_VENTS     = 8
   integer, parameter :: MAX_OBSTACLES = 64

   ! Spec
   integer  :: Nx, Ny, Nz
   real(wp) :: Lx, Ly, Lz
   real(wp) :: room_z0, room_z1
   integer  :: poly_n
   real(wp) :: poly_x(MAX_POLY), poly_y(MAX_POLY)
   integer  :: n_vents
   integer  :: vent_face(MAX_VENTS)
   real(wp) :: vent_u0(MAX_VENTS), vent_v0(MAX_VENTS)
   real(wp) :: vent_u1(MAX_VENTS), vent_v1(MAX_VENTS)
   integer  :: n_obs
   real(wp) :: obs_x0(MAX_OBSTACLES), obs_y0(MAX_OBSTACLES), obs_z0(MAX_OBSTACLES)
   real(wp) :: obs_x1(MAX_OBSTACLES), obs_y1(MAX_OBSTACLES), obs_z1(MAX_OBSTACLES)

   character(len=512) :: spec_file, out_file
   integer, allocatable :: cell_kind(:,:,:)
   real(wp) :: dx, dy, dz
   integer  :: nvx, nvy, nvz, nv
   integer  :: u_out, ios
   integer  :: i, j, k, eid
   integer  :: nc_air
   integer  :: nbf_walls, nbf_vent, nbf_obs, nbf_total

   if (command_argument_count() < 2) then
      write(error_unit,'(A)') 'usage: gen_room_mesh <spec.txt> <out.msh>'
      stop 1
   end if
   call get_command_argument(1, spec_file)
   call get_command_argument(2, out_file)

   Nx = 0;  Ny = 0;  Nz = 0
   Lx = 0.0_wp;  Ly = 0.0_wp;  Lz = 0.0_wp
   room_z0 = 0.0_wp;  room_z1 = 0.0_wp
   poly_n  = 0;  n_vents = 0;  n_obs = 0

   call parse_spec(spec_file)

   if (Nx <= 0 .or. Ny <= 0 .or. Nz <= 0) then
      write(error_unit,'(A)') 'spec: GRID Nx Ny Nz must be positive'; stop 1
   end if
   if (Lx <= 0.0_wp .or. Ly <= 0.0_wp .or. Lz <= 0.0_wp) then
      write(error_unit,'(A)') 'spec: BBOX Lx Ly Lz must be positive'; stop 1
   end if
   if (poly_n < 3) then
      write(error_unit,'(A)') 'spec: POLYGON must have at least 3 vertices'; stop 1
   end if
   if (room_z1 <= room_z0) then
      write(error_unit,'(A)') 'spec: ROOM_Z must have z1 > z0'; stop 1
   end if

   nvx = Nx + 1;  nvy = Ny + 1;  nvz = Nz + 1
   nv  = nvx*nvy*nvz
   dx = Lx / real(Nx, wp)
   dy = Ly / real(Ny, wp)
   dz = Lz / real(Nz, wp)

   allocate(cell_kind(Nx, Ny, Nz))
   call classify_cells()
   nc_air = count(cell_kind == KIND_AIR)

   nbf_walls = 0;  nbf_vent = 0;  nbf_obs = 0
   call visit_boundary_faces(do_emit=.false., write_unit=0)
   nbf_total = nbf_walls + nbf_vent + nbf_obs

   open(newunit=u_out, file=trim(out_file), status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      write(error_unit,'(A,A)') 'cannot open ', trim(out_file); stop 1
   end if

   write(u_out,'(A)') '$MeshFormat'
   write(u_out,'(A)') '2.2 0 8'
   write(u_out,'(A)') '$EndMeshFormat'

   write(u_out,'(A)') '$PhysicalNames'
   write(u_out,'(I0)') 3
   write(u_out,'(I0,1X,I0,1X,A)') 2, TAG_WALLS,    '"walls"'
   write(u_out,'(I0,1X,I0,1X,A)') 2, TAG_VENT,     '"vent_outside"'
   write(u_out,'(I0,1X,I0,1X,A)') 2, TAG_OBSTACLE, '"obstacle_faces"'
   write(u_out,'(A)') '$EndPhysicalNames'

   write(u_out,'(A)') '$Nodes'
   write(u_out,'(I0)') nv
   do k = 1, nvz
      do j = 1, nvy
         do i = 1, nvx
            write(u_out,'(I0,3(1X,1PE22.14))') &
               node_id(i, j, k), &
               real(i-1, wp)*dx, real(j-1, wp)*dy, real(k-1, wp)*dz
         end do
      end do
   end do
   write(u_out,'(A)') '$EndNodes'

   write(u_out,'(A)') '$Elements'
   write(u_out,'(I0)') nc_air + nbf_total
   eid = 0
   call emit_hexes()

   ! Reset counters (so visit_boundary_faces emit pass is clean) and emit
   nbf_walls = 0;  nbf_vent = 0;  nbf_obs = 0
   call visit_boundary_faces(do_emit=.true., write_unit=u_out)
   write(u_out,'(A)') '$EndElements'
   close(u_out)

   write(*,'(A,I0,A)') 'gen_room_mesh: wrote ', nc_air, ' air cells'
   write(*,'(A,I0,A,I0,A,I0,A,I0)') &
      '   walls=', nbf_walls, ' vent=', nbf_vent, &
      ' obstacle=', nbf_obs, ' total_bnd=', nbf_total

contains

   subroutine parse_spec(fname)
      character(len=*), intent(in) :: fname
      integer :: u_in, io, iv, dummy_n
      character(len=512) :: line
      character(len=32)  :: kw

      open(newunit=u_in, file=trim(fname), status='old', action='read', iostat=io)
      if (io /= 0) then
         write(error_unit,'(A,A)') 'cannot open spec file ', trim(fname); stop 1
      end if

      dummy_n = 0
      do
         read(u_in, '(A)', iostat=io) line
         if (io == iostat_end) exit
         if (io /= 0) then
            write(error_unit,'(A)') 'spec: read error'; stop 1
         end if
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle

         read(line, *) kw
         select case (trim(kw))
         case ('GRID')
            read(line, *) kw, Nx, Ny, Nz
         case ('BBOX')
            read(line, *) kw, Lx, Ly, Lz
         case ('ROOM_Z')
            read(line, *) kw, room_z0, room_z1
         case ('POLYGON')
            read(line, *) kw, dummy_n
            if (dummy_n > MAX_POLY) then
               write(error_unit,'(A,I0)') 'spec: polygon vertices exceed MAX_POLY=', MAX_POLY
               stop 1
            end if
            poly_n = 0
            do iv = 1, dummy_n
               read(u_in, '(A)', iostat=io) line
               if (io /= 0) then
                  write(error_unit,'(A)') 'spec: POLYGON truncated'; stop 1
               end if
               line = adjustl(line)
               if (line(1:5) == 'VERT ' .or. line(1:5) == 'VERT'//char(9)) then
                  read(line, *) kw, poly_x(iv), poly_y(iv)
               else
                  read(line, *) poly_x(iv), poly_y(iv)
               end if
               poly_n = iv
            end do
         case ('VENT')
            n_vents = n_vents + 1
            if (n_vents > MAX_VENTS) then
               write(error_unit,'(A)') 'spec: too many vents'; stop 1
            end if
            read(line, *) kw, vent_face(n_vents), &
                              vent_u0(n_vents), vent_v0(n_vents), &
                              vent_u1(n_vents), vent_v1(n_vents)
         case ('CUBE')
            n_obs = n_obs + 1
            if (n_obs > MAX_OBSTACLES) then
               write(error_unit,'(A)') 'spec: too many obstacles'; stop 1
            end if
            read(line, *) kw, obs_x0(n_obs), obs_y0(n_obs), obs_z0(n_obs), &
                              obs_x1(n_obs), obs_y1(n_obs), obs_z1(n_obs)
         case ('END')
            exit
         case default
            write(error_unit,'(A,A)') 'spec: unknown keyword ', trim(kw)
            stop 1
         end select
      end do
      close(u_in)
   end subroutine parse_spec

   subroutine classify_cells()
      integer  :: ic, jc, kc
      real(wp) :: xc, yc, zc
      do kc = 1, Nz
         zc = (real(kc, wp) - 0.5_wp) * dz
         do jc = 1, Ny
            yc = (real(jc, wp) - 0.5_wp) * dy
            do ic = 1, Nx
               xc = (real(ic, wp) - 0.5_wp) * dx
               if (in_any_obstacle(xc, yc, zc)) then
                  cell_kind(ic, jc, kc) = KIND_VOID_OBSTACLE
               else if (zc >= room_z0 .and. zc <= room_z1 .and. &
                        point_in_polygon(xc, yc)) then
                  cell_kind(ic, jc, kc) = KIND_AIR
               else
                  cell_kind(ic, jc, kc) = KIND_VOID_OUTSIDE
               end if
            end do
         end do
      end do
   end subroutine classify_cells

   pure logical function point_in_polygon(px, py) result(inside)
      real(wp), intent(in) :: px, py
      integer :: ii, jj
      inside = .false.
      jj = poly_n
      do ii = 1, poly_n
         if ((poly_y(ii) > py) .neqv. (poly_y(jj) > py)) then
            if (px < (poly_x(jj) - poly_x(ii)) * (py - poly_y(ii)) / &
                     (poly_y(jj) - poly_y(ii)) + poly_x(ii)) then
               inside = .not. inside
            end if
         end if
         jj = ii
      end do
   end function point_in_polygon

   pure logical function in_any_obstacle(xc, yc, zc) result(yes)
      real(wp), intent(in) :: xc, yc, zc
      integer :: ii
      yes = .false.
      do ii = 1, n_obs
         if (xc >= obs_x0(ii) .and. xc <= obs_x1(ii) .and. &
             yc >= obs_y0(ii) .and. yc <= obs_y1(ii) .and. &
             zc >= obs_z0(ii) .and. zc <= obs_z1(ii)) then
            yes = .true.
            return
         end if
      end do
   end function in_any_obstacle

   pure integer function node_id(ii, jj, kk) result(id)
      integer, intent(in) :: ii, jj, kk
      id = (kk - 1) * nvx * nvy + (jj - 1) * nvx + ii
   end function node_id

   pure integer function cell_kind_safe(ii, jj, kk) result(t)
      integer, intent(in) :: ii, jj, kk
      if (ii < 1 .or. ii > Nx .or. jj < 1 .or. jj > Ny .or. kk < 1 .or. kk > Nz) then
         t = KIND_OFF_DOMAIN
      else
         t = cell_kind(ii, jj, kk)
      end if
   end function cell_kind_safe

   pure integer function pick_outer_face_tag(face_id, uc, vc) result(tag)
      integer,  intent(in) :: face_id
      real(wp), intent(in) :: uc, vc
      integer :: v
      do v = 1, n_vents
         if (vent_face(v) == face_id) then
            if (uc >= vent_u0(v) .and. uc <= vent_u1(v) .and. &
                vc >= vent_v0(v) .and. vc <= vent_v1(v)) then
               tag = TAG_VENT
               return
            end if
         end if
      end do
      tag = TAG_WALLS
   end function pick_outer_face_tag

   subroutine emit_hexes()
      integer :: ic, jc, kc
      do kc = 1, Nz
         do jc = 1, Ny
            do ic = 1, Nx
               if (cell_kind(ic, jc, kc) == KIND_AIR) then
                  eid = eid + 1
                  write(u_out,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,8(1X,I0))') &
                     eid, 5, 2, 0, 0, &
                     node_id(ic,   jc,   kc  ), &
                     node_id(ic+1, jc,   kc  ), &
                     node_id(ic+1, jc+1, kc  ), &
                     node_id(ic,   jc+1, kc  ), &
                     node_id(ic,   jc,   kc+1), &
                     node_id(ic+1, jc,   kc+1), &
                     node_id(ic+1, jc+1, kc+1), &
                     node_id(ic,   jc+1, kc+1)
               end if
            end do
         end do
      end do
   end subroutine emit_hexes

   ! Visit every face on the structured candidate grid. For each face
   ! with exactly one air neighbour: increment the per-tag count, and
   ! (if do_emit) also write the quad to write_unit.
   subroutine visit_boundary_faces(do_emit, write_unit)
      logical, intent(in) :: do_emit
      integer, intent(in) :: write_unit
      integer :: ic, jc, kc, tag
      integer :: q(4)

      ! x-faces at node plane i ∈ 1..nvx
      do kc = 1, Nz
         do jc = 1, Ny
            do ic = 1, nvx
               call classify_x_face(ic, jc, kc, tag)
               if (tag == 0) cycle
               call bump_count(tag)
               if (do_emit) then
                  q(1) = node_id(ic, jc,   kc  )
                  q(2) = node_id(ic, jc+1, kc  )
                  q(3) = node_id(ic, jc+1, kc+1)
                  q(4) = node_id(ic, jc,   kc+1)
                  eid = eid + 1
                  write(write_unit,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,4(1X,I0))') &
                     eid, 3, 2, tag, tag, q(1), q(2), q(3), q(4)
               end if
            end do
         end do
      end do
      ! y-faces
      do kc = 1, Nz
         do jc = 1, nvy
            do ic = 1, Nx
               call classify_y_face(ic, jc, kc, tag)
               if (tag == 0) cycle
               call bump_count(tag)
               if (do_emit) then
                  q(1) = node_id(ic,   jc, kc  )
                  q(2) = node_id(ic+1, jc, kc  )
                  q(3) = node_id(ic+1, jc, kc+1)
                  q(4) = node_id(ic,   jc, kc+1)
                  eid = eid + 1
                  write(write_unit,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,4(1X,I0))') &
                     eid, 3, 2, tag, tag, q(1), q(2), q(3), q(4)
               end if
            end do
         end do
      end do
      ! z-faces
      do kc = 1, nvz
         do jc = 1, Ny
            do ic = 1, Nx
               call classify_z_face(ic, jc, kc, tag)
               if (tag == 0) cycle
               call bump_count(tag)
               if (do_emit) then
                  q(1) = node_id(ic,   jc,   kc)
                  q(2) = node_id(ic+1, jc,   kc)
                  q(3) = node_id(ic+1, jc+1, kc)
                  q(4) = node_id(ic,   jc+1, kc)
                  eid = eid + 1
                  write(write_unit,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,4(1X,I0))') &
                     eid, 3, 2, tag, tag, q(1), q(2), q(3), q(4)
               end if
            end do
         end do
      end do
   end subroutine visit_boundary_faces

   subroutine classify_x_face(i_node, j_idx, k_idx, tag)
      integer, intent(in)  :: i_node, j_idx, k_idx
      integer, intent(out) :: tag
      integer  :: kL, kR, void_kind
      real(wp) :: uc, vc
      kL = cell_kind_safe(i_node - 1, j_idx, k_idx)
      kR = cell_kind_safe(i_node,     j_idx, k_idx)
      if ((kL == KIND_AIR) .eqv. (kR == KIND_AIR)) then
         tag = 0; return
      end if
      if (kL == KIND_AIR) then
         void_kind = kR
      else
         void_kind = kL
      end if
      select case (void_kind)
      case (KIND_OFF_DOMAIN)
         uc = (real(j_idx, wp) - 0.5_wp) * dy
         vc = (real(k_idx, wp) - 0.5_wp) * dz
         if (i_node == 1) then
            tag = pick_outer_face_tag(0, uc, vc)
         else
            tag = pick_outer_face_tag(1, uc, vc)
         end if
      case (KIND_VOID_OUTSIDE)
         tag = TAG_WALLS
      case (KIND_VOID_OBSTACLE)
         tag = TAG_OBSTACLE
      case default
         tag = 0
      end select
   end subroutine classify_x_face

   subroutine classify_y_face(i_idx, j_node, k_idx, tag)
      integer, intent(in)  :: i_idx, j_node, k_idx
      integer, intent(out) :: tag
      integer  :: kL, kR, void_kind
      real(wp) :: uc, vc
      kL = cell_kind_safe(i_idx, j_node - 1, k_idx)
      kR = cell_kind_safe(i_idx, j_node,     k_idx)
      if ((kL == KIND_AIR) .eqv. (kR == KIND_AIR)) then
         tag = 0; return
      end if
      if (kL == KIND_AIR) then
         void_kind = kR
      else
         void_kind = kL
      end if
      select case (void_kind)
      case (KIND_OFF_DOMAIN)
         uc = (real(i_idx, wp) - 0.5_wp) * dx
         vc = (real(k_idx, wp) - 0.5_wp) * dz
         if (j_node == 1) then
            tag = pick_outer_face_tag(2, uc, vc)
         else
            tag = pick_outer_face_tag(3, uc, vc)
         end if
      case (KIND_VOID_OUTSIDE)
         tag = TAG_WALLS
      case (KIND_VOID_OBSTACLE)
         tag = TAG_OBSTACLE
      case default
         tag = 0
      end select
   end subroutine classify_y_face

   subroutine classify_z_face(i_idx, j_idx, k_node, tag)
      integer, intent(in)  :: i_idx, j_idx, k_node
      integer, intent(out) :: tag
      integer  :: kL, kR, void_kind
      real(wp) :: uc, vc
      kL = cell_kind_safe(i_idx, j_idx, k_node - 1)
      kR = cell_kind_safe(i_idx, j_idx, k_node    )
      if ((kL == KIND_AIR) .eqv. (kR == KIND_AIR)) then
         tag = 0; return
      end if
      if (kL == KIND_AIR) then
         void_kind = kR
      else
         void_kind = kL
      end if
      select case (void_kind)
      case (KIND_OFF_DOMAIN)
         uc = (real(i_idx, wp) - 0.5_wp) * dx
         vc = (real(j_idx, wp) - 0.5_wp) * dy
         if (k_node == 1) then
            tag = pick_outer_face_tag(4, uc, vc)
         else
            tag = pick_outer_face_tag(5, uc, vc)
         end if
      case (KIND_VOID_OUTSIDE)
         tag = TAG_WALLS
      case (KIND_VOID_OBSTACLE)
         tag = TAG_OBSTACLE
      case default
         tag = 0
      end select
   end subroutine classify_z_face

   subroutine bump_count(tag)
      integer, intent(in) :: tag
      select case (tag)
      case (TAG_WALLS);    nbf_walls = nbf_walls + 1
      case (TAG_VENT);     nbf_vent  = nbf_vent  + 1
      case (TAG_OBSTACLE); nbf_obs   = nbf_obs   + 1
      end select
   end subroutine bump_count

end program gen_room_mesh
