! M2 partitioning + halo construction.
!
! Architecture for M2 first iteration: every rank reads the full Gmsh mesh
! (cheap relative to the solve) and runs build_topology on it once for cell
! adjacency. Rank 0 calls METIS_PartMeshDual to produce a per-cell partition
! vector; this is MPI_Bcast to all ranks. Each rank then locally extracts its
! slice (local cells + one ghost layer) into a new t_mesh, rebuilds local
! topology + metrics, reorders pure-interior / partition faces, and builds
! the halo exchange descriptor via point-to-point neighbor discovery.
!
! True distributed-input ParMETIS is deferred to M2.5 once mesh sizes
! exceed per-rank memory.
module partition
   use, intrinsic :: iso_c_binding
   use kinds,         only : wp
   use constants,     only : NVAR, NPRIM
   use mpi_f08
   use mpi_runtime,   only : t_mpi_ctx, mpi_abort_msg
   use mesh_types,    only : t_mesh, t_patch, t_halo
   use mesh_io_gmsh,  only : gmsh_read, t_bface_raw
   use mesh_topology, only : build_topology
   use mesh_metrics,  only : build_metrics
   use cut_cell,      only : build_cut_cell_tables
   use metis_binding, only : c_idx, METIS_PartMeshDual_C, METIS_OK
   implicit none
   private

   public :: partition_and_load

contains

   subroutine partition_and_load(filename, ctx, mesh)
      character(len=*), intent(in)  :: filename
      type(t_mpi_ctx),  intent(in)  :: ctx
      type(t_mesh),     intent(out) :: mesh

      type(t_mesh) :: gmesh
      type(t_bface_raw), allocatable :: bfaces(:)
      integer :: nb
      integer, allocatable :: epart(:)         ! (nc_global) rank assignment per cell

      ! 1) All ranks read the full mesh (M2 simplification).
      call gmsh_read(filename, gmesh, bfaces, nb)

      ! 2) All ranks build global topology (for cell-to-cell adjacency).
      call build_topology(gmesh, bfaces, nb)
      ! Metrics not needed for partition; we'll build metrics on the LOCAL mesh.

      ! 3) Rank 0 partitions; broadcast.
      allocate(epart(gmesh%nc_internal))
      if (ctx%is_root) then
         if (ctx%nproc == 1) then
            epart = 0
         else
            call compute_metis_partition(gmesh, ctx%nproc, epart)
         end if
      end if
      call MPI_Bcast(epart, gmesh%nc_internal, MPI_INTEGER, 0, ctx%comm)

      ! 4) Each rank extracts its slice into `mesh`.
      call extract_local(gmesh, bfaces, nb, epart, ctx, mesh)

      ! 5) Build local topology + metrics.
      block
         type(t_bface_raw), allocatable :: lbfaces(:)
         integer :: nl
         call select_local_bfaces(gmesh, bfaces, nb, epart, ctx%rank, lbfaces, nl)
         call build_topology(mesh, lbfaces, nl)
         call build_metrics(mesh)
         if (allocated(lbfaces)) deallocate(lbfaces)
      end block

      ! 6) Reclassify interior faces into pure-interior vs partition,
      !    reorder so pure-interior come first.
      call classify_partition_faces(mesh)

      ! 7) Build halo descriptor (CSR send/recv lists + persistent MPI requests).
      call build_halo_desc(mesh, ctx)

      ! 8) Default (identity) cut-cell tables so every mesh — including
      !    those loaded directly by tests that don't go through run_case —
      !    can be read by the solver. run_case overrides with the namelist
      !    obstacle list right after this returns.
      block
         real(wp) :: no_cube_lo(3, 1), no_cube_hi(3, 1)
         no_cube_lo = 0.0_wp
         no_cube_hi = 0.0_wp
         call build_cut_cell_tables(mesh, 0, no_cube_lo, no_cube_hi)
      end block

      ! Cleanup global mesh
      deallocate(epart)
      call free_mesh(gmesh)
      if (allocated(bfaces)) deallocate(bfaces)
   end subroutine partition_and_load

   ! --- METIS partition ------------------------------------------------------

   subroutine compute_metis_partition(gmesh, nparts, epart)
      type(t_mesh), intent(in)  :: gmesh
      integer,      intent(in)  :: nparts
      integer,      intent(out) :: epart(:)        ! (nc_global) 0-based rank assignment

      integer(c_idx) :: ne, nn, ncommon, np_metis, objval
      integer(c_idx), allocatable, target :: eptr(:), eind(:), opts(:), ep(:), npart(:)
      integer(c_int) :: ierr
      integer :: i, csr_size

      ne = int(gmesh%nc_internal, c_idx)
      nn = int(gmesh%nv,          c_idx)
      ncommon = 3_c_idx   ! 3 shared nodes ⇒ neighbors (works for tet/hex/prism/pyramid)
      np_metis = int(nparts, c_idx)

      csr_size = gmesh%cell_vtx_ptr(gmesh%nc_internal + 1) - 1
      allocate(eptr(gmesh%nc_internal + 1))
      allocate(eind(csr_size))
      ! METIS expects 0-based indices.
      eptr = int(gmesh%cell_vtx_ptr - 1, c_idx)
      eind = int(gmesh%cell_vtx     - 1, c_idx)

      allocate(opts(40))   ! METIS_NOPTIONS = 40
      opts = -1_c_idx       ! all options to default

      allocate(ep   (gmesh%nc_internal))
      allocate(npart(gmesh%nv))

      ierr = METIS_PartMeshDual_C(ne, nn, eptr, eind, c_null_ptr, c_null_ptr, &
                                   ncommon, np_metis, c_null_ptr, opts, &
                                   objval, ep, npart)
      if (ierr /= METIS_OK) then
         write(*,'(A,I0)') 'METIS_PartMeshDual returned error ', ierr
         error stop 1
      end if

      do i = 1, gmesh%nc_internal
         epart(i) = int(ep(i))
      end do

      deallocate(eptr, eind, opts, ep, npart)
   end subroutine compute_metis_partition

   ! --- Local mesh extraction ------------------------------------------------

   subroutine extract_local(gmesh, bfaces, nb, epart, ctx, mesh)
      type(t_mesh),      intent(in)    :: gmesh
      type(t_bface_raw), intent(in)    :: bfaces(:)
      integer,           intent(in)    :: nb
      integer,           intent(in)    :: epart(:)
      type(t_mpi_ctx),   intent(in)    :: ctx
      type(t_mesh),      intent(out)   :: mesh

      integer :: c, f, o, nbr, idx, nv_pc, k, vptr_g, vptr_l, csr_size
      integer :: nc_internal_loc, nc_total_loc
      integer, allocatable :: glob_to_loc(:)         ! (nc_global) → local cell idx or 0
      integer, allocatable :: ghost_set(:)           ! list of global cell IDs
      integer :: n_ghost

      ! 1) Count local cells.
      nc_internal_loc = 0
      do c = 1, gmesh%nc_internal
         if (epart(c) == ctx%rank) nc_internal_loc = nc_internal_loc + 1
      end do

      ! 2) Identify ghost cells via global mesh's interior face adjacency.
      allocate(glob_to_loc(gmesh%nc_internal))
      glob_to_loc = 0
      idx = 0
      do c = 1, gmesh%nc_internal
         if (epart(c) == ctx%rank) then
            idx = idx + 1
            glob_to_loc(c) = idx
         end if
      end do
      ! glob_to_loc now contains local IDs for local cells; 0 elsewhere.

      ! Collect ghost candidates (cells on other ranks that touch a local cell).
      allocate(ghost_set(0))     ! growable
      n_ghost = 0
      do f = 1, gmesh%nf_interior
         o   = gmesh%face_owner   (f)
         nbr = gmesh%face_neighbor(f)
         if (epart(o) == ctx%rank .and. epart(nbr) /= ctx%rank) then
            call insert_unique(ghost_set, n_ghost, nbr)
         else if (epart(nbr) == ctx%rank .and. epart(o) /= ctx%rank) then
            call insert_unique(ghost_set, n_ghost, o)
         end if
      end do

      ! 3) Assign local IDs to ghost cells [nc_internal+1 .. nc_total].
      do k = 1, n_ghost
         glob_to_loc(ghost_set(k)) = nc_internal_loc + k
      end do
      nc_total_loc = nc_internal_loc + n_ghost

      ! 4) Allocate local mesh.
      mesh%nv = gmesh%nv
      allocate(mesh%xv(3, gmesh%nv), source=gmesh%xv)   ! replicate global vertex coords

      mesh%nc_internal = nc_internal_loc
      mesh%nc_total    = nc_total_loc
      allocate(mesh%cell_type      (nc_total_loc))
      allocate(mesh%cell_vtx_ptr   (nc_total_loc + 1))
      allocate(mesh%cell_global_id (nc_total_loc))
      allocate(mesh%cell_owner_rank(nc_total_loc))

      ! Compute CSR ptr sizes (sum of nv_per_cell across local + ghost cells)
      csr_size = 0
      do c = 1, gmesh%nc_internal
         if (glob_to_loc(c) > 0) then
            csr_size = csr_size + (gmesh%cell_vtx_ptr(c+1) - gmesh%cell_vtx_ptr(c))
         end if
      end do
      allocate(mesh%cell_vtx(csr_size))

      ! Fill local cells first (preserve global order among locals)
      idx = 0
      vptr_l = 1
      mesh%cell_vtx_ptr(1) = 1
      do c = 1, gmesh%nc_internal
         if (epart(c) == ctx%rank) then
            idx = idx + 1
            mesh%cell_type(idx)       = gmesh%cell_type(c)
            mesh%cell_global_id(idx)  = c
            mesh%cell_owner_rank(idx) = ctx%rank
            nv_pc = gmesh%cell_vtx_ptr(c+1) - gmesh%cell_vtx_ptr(c)
            vptr_g = gmesh%cell_vtx_ptr(c)
            do k = 1, nv_pc
               mesh%cell_vtx(vptr_l + k - 1) = gmesh%cell_vtx(vptr_g + k - 1)
            end do
            vptr_l = vptr_l + nv_pc
            mesh%cell_vtx_ptr(idx + 1) = vptr_l
         end if
      end do
      ! Then ghosts
      do k = 1, n_ghost
         c = ghost_set(k)
         idx = idx + 1
         mesh%cell_type(idx)       = gmesh%cell_type(c)
         mesh%cell_global_id(idx)  = c
         mesh%cell_owner_rank(idx) = epart(c)
         nv_pc = gmesh%cell_vtx_ptr(c+1) - gmesh%cell_vtx_ptr(c)
         vptr_g = gmesh%cell_vtx_ptr(c)
         do o = 1, nv_pc
            mesh%cell_vtx(vptr_l + o - 1) = gmesh%cell_vtx(vptr_g + o - 1)
         end do
         vptr_l = vptr_l + nv_pc
         mesh%cell_vtx_ptr(idx + 1) = vptr_l
      end do

      ! 5) Replicate patches (names + gmsh_tag); face counts will be set by topology.
      mesh%np = gmesh%np
      if (gmesh%np > 0) then
         allocate(mesh%patches(gmesh%np))
         do k = 1, gmesh%np
            mesh%patches(k)%name     = gmesh%patches(k)%name
            mesh%patches(k)%gmsh_tag = gmesh%patches(k)%gmsh_tag
            mesh%patches(k)%bc_type  = 0
            mesh%patches(k)%face_start = 0
            mesh%patches(k)%face_count = 0
         end do
      end if

      ! Suppress unused-warning: bfaces / nb consumed by select_local_bfaces upstream.
      if (.false.) then
         k = nb + size(bfaces)
      end if

      deallocate(ghost_set, glob_to_loc)
   end subroutine extract_local

   subroutine insert_unique(arr, n, val)
      integer, allocatable, intent(inout) :: arr(:)
      integer,              intent(inout) :: n
      integer,              intent(in)    :: val
      integer :: i, cap
      integer, allocatable :: tmp(:)
      do i = 1, n
         if (arr(i) == val) return
      end do
      cap = size(arr)
      if (n + 1 > cap) then
         allocate(tmp(max(8, 2*cap)))
         if (n > 0) tmp(1:n) = arr(1:n)
         call move_alloc(tmp, arr)
      end if
      n = n + 1
      arr(n) = val
   end subroutine insert_unique

   ! Filter global boundary-face list to those incident on a local cell.
   subroutine select_local_bfaces(gmesh, bfaces, nb, epart, my_rank, &
                                  lbfaces, nl)
      type(t_mesh),      intent(in)  :: gmesh
      type(t_bface_raw), intent(in)  :: bfaces(:)
      integer,           intent(in)  :: nb, my_rank
      integer,           intent(in)  :: epart(:)
      type(t_bface_raw), allocatable, intent(out) :: lbfaces(:)
      integer,           intent(out) :: nl

      ! Even simpler: pass all global bfaces. Build_topology matches by vertex set;
      ! bfaces whose vertices aren't on any local boundary candidate are silently
      ! unmatched (and don't show up as orphans because we filter via topology).
      integer :: i
      nl = nb
      allocate(lbfaces(nl))
      do i = 1, nb
         lbfaces(i) = bfaces(i)
      end do
      ! Suppress unused-warnings
      if (.false.) then
         i = gmesh%nv + epart(1) + my_rank
      end if
   end subroutine select_local_bfaces

   ! --- Face classification + reordering ------------------------------------

   ! Classify faces after build_topology runs on local+ghost cells, and
   ! reorder so:
   !   [1 .. nf_pure_interior]   — both cells local
   !   [.. nf_interior]          — one local, one ghost (partition faces)
   !   [.. nf]                   — physical BC faces (owner local), grouped by patch
   ! Faces that don't involve any local cell (ghost-ghost interior, ghost-owned
   ! boundary) are dropped — they're computed by the rank that owns them.
   subroutine classify_partition_faces(mesh)
      type(t_mesh), intent(inout) :: mesh
      integer :: f, nf_old, npure, npart, nb_keep
      integer :: nc_int
      integer, allocatable :: keep(:)
      integer, allocatable :: kind(:)     ! 1=pure 2=partition 3=boundary 0=drop
      integer :: nf_new
      integer, allocatable :: order(:)
      integer, allocatable :: tmp_owner(:), tmp_owner_lf(:), tmp_neigh(:), tmp_patch(:)
      real(wp), allocatable :: tmp_area(:), tmp_normal(:,:), tmp_centroid(:,:)

      nf_old = mesh%nf
      nc_int = mesh%nc_internal

      allocate(kind(nf_old)); kind = 0

      npure = 0
      npart = 0
      do f = 1, mesh%nf_interior
         block
            integer :: o, n
            o = mesh%face_owner(f)
            n = mesh%face_neighbor(f)
            if (o <= nc_int .and. n <= nc_int) then
               kind(f) = 1; npure = npure + 1
            else if (o <= nc_int .or. n <= nc_int) then
               kind(f) = 2; npart = npart + 1
            else
               kind(f) = 0     ! ghost-ghost: drop
            end if
         end block
      end do
      nb_keep = 0
      do f = mesh%nf_interior + 1, nf_old
         if (mesh%face_owner(f) <= nc_int .and. mesh%face_patch(f) > 0) then
            kind(f) = 3; nb_keep = nb_keep + 1
         else
            kind(f) = 0
         end if
      end do
      nf_new = npure + npart + nb_keep

      ! Build new order: pure → partition → boundary (by patch order).
      allocate(order(nf_new))
      block
         integer :: cur_pure, cur_part, cur_bnd
         cur_pure = 0
         cur_part = npure
         cur_bnd  = npure + npart
         ! For pure/partition: emit in input order.
         do f = 1, mesh%nf_interior
            if (kind(f) == 1) then
               cur_pure = cur_pure + 1
               order(cur_pure) = f
            else if (kind(f) == 2) then
               cur_part = cur_part + 1
               order(cur_part) = f
            end if
         end do
         ! Boundary: emit grouped by patch (existing layout already groups them).
         do f = mesh%nf_interior + 1, nf_old
            if (kind(f) == 3) then
               cur_bnd = cur_bnd + 1
               order(cur_bnd) = f
            end if
         end do
      end block

      ! Apply.
      allocate(tmp_owner(nf_new), tmp_owner_lf(nf_new), tmp_neigh(nf_new), tmp_patch(nf_new))
      allocate(tmp_area(nf_new), tmp_normal(3, nf_new), tmp_centroid(3, nf_new))
      do f = 1, nf_new
         tmp_owner    (f) = mesh%face_owner    (order(f))
         tmp_owner_lf (f) = mesh%face_owner_lf (order(f))
         tmp_neigh    (f) = mesh%face_neighbor (order(f))
         tmp_patch    (f) = mesh%face_patch    (order(f))
         tmp_area     (f) = mesh%face_area     (order(f))
         tmp_normal (:,f) = mesh%face_normal (:, order(f))
         tmp_centroid(:,f)= mesh%face_centroid(:,order(f))
      end do
      deallocate(mesh%face_owner, mesh%face_owner_lf, mesh%face_neighbor, mesh%face_patch)
      deallocate(mesh%face_area, mesh%face_normal, mesh%face_centroid)
      call move_alloc(tmp_owner,    mesh%face_owner)
      call move_alloc(tmp_owner_lf, mesh%face_owner_lf)
      call move_alloc(tmp_neigh,    mesh%face_neighbor)
      call move_alloc(tmp_patch,    mesh%face_patch)
      call move_alloc(tmp_area,     mesh%face_area)
      call move_alloc(tmp_normal,   mesh%face_normal)
      call move_alloc(tmp_centroid, mesh%face_centroid)

      mesh%nf               = nf_new
      mesh%nf_pure_interior = npure
      mesh%nf_interior      = npure + npart
      mesh%nf_boundary      = nb_keep

      ! Recompute patch face_start / face_count for the kept boundary faces.
      block
         integer :: p
         do p = 1, mesh%np
            mesh%patches(p)%face_count = 0
         end do
         do f = mesh%nf_interior + 1, mesh%nf
            block
               integer :: p2
               p2 = mesh%face_patch(f)
               if (p2 > 0) mesh%patches(p2)%face_count = mesh%patches(p2)%face_count + 1
            end block
         end do
         block
            integer :: cur
            cur = mesh%nf_interior + 1
            do p = 1, mesh%np
               mesh%patches(p)%face_start = cur
               cur = cur + mesh%patches(p)%face_count
            end do
         end block
      end block

      deallocate(kind, order)
      if (allocated(keep)) deallocate(keep)
   end subroutine classify_partition_faces

   ! --- Halo descriptor ------------------------------------------------------

   subroutine build_halo_desc(mesh, ctx)
      type(t_mesh),    intent(inout) :: mesh
      type(t_mpi_ctx), intent(in)    :: ctx

      integer :: g, r, k, np, nproc
      integer, allocatable :: send_count(:), recv_count(:)
      integer, allocatable :: recv_global_ids(:)   ! flat: cells we want from each neighbor
      integer, allocatable :: recv_offset_per_rank(:)
      integer, allocatable :: send_global_ids(:)   ! flat: cells they want from us
      integer, allocatable :: neighbor_idx(:)      ! map rank → idx in neighbor list

      nproc = ctx%nproc

      ! 1) Build per-rank recv_count from ghost cells' owners.
      allocate(recv_count(0:nproc-1), source=0)
      do g = mesh%nc_internal + 1, mesh%nc_total
         r = mesh%cell_owner_rank(g)
         recv_count(r) = recv_count(r) + 1
      end do
      ! self-recv is zero by construction

      ! 2) MPI_Alltoall to learn send_count (how many cells each rank wants from us).
      allocate(send_count(0:nproc-1), source=0)
      call MPI_Alltoall(recv_count, 1, MPI_INTEGER, &
                        send_count, 1, MPI_INTEGER, ctx%comm)

      ! 3) Determine neighbor list (ranks with non-zero exchange in either direction).
      np = 0
      do r = 0, nproc - 1
         if (r == ctx%rank) cycle
         if (recv_count(r) > 0 .or. send_count(r) > 0) np = np + 1
      end do
      mesh%halo%n_neighbors = np
      allocate(mesh%halo%neighbor_ranks(np))
      allocate(neighbor_idx(0:nproc-1)); neighbor_idx = -1
      k = 0
      do r = 0, nproc - 1
         if (r == ctx%rank) cycle
         if (recv_count(r) > 0 .or. send_count(r) > 0) then
            k = k + 1
            mesh%halo%neighbor_ranks(k) = r
            neighbor_idx(r) = k
         end if
      end do

      ! 4) Build recv CSR.
      allocate(mesh%halo%recv_offset(np + 1))
      mesh%halo%recv_offset(1) = 1
      do k = 1, np
         r = mesh%halo%neighbor_ranks(k)
         mesh%halo%recv_offset(k+1) = mesh%halo%recv_offset(k) + recv_count(r)
      end do
      allocate(mesh%halo%recv_ghosts(mesh%halo%recv_offset(np+1) - 1))

      ! Per-rank temp cursors to fill recv_ghosts in order
      allocate(recv_offset_per_rank(0:nproc-1), source=0)
      do k = 1, np
         r = mesh%halo%neighbor_ranks(k)
         recv_offset_per_rank(r) = mesh%halo%recv_offset(k) - 1   ! cursor
      end do

      ! Fill recv_ghosts AND build the parallel global_ids buffer to send to owner ranks.
      allocate(recv_global_ids(mesh%halo%recv_offset(np+1) - 1))
      do g = mesh%nc_internal + 1, mesh%nc_total
         r = mesh%cell_owner_rank(g)
         if (neighbor_idx(r) <= 0) cycle    ! shouldn't happen
         recv_offset_per_rank(r) = recv_offset_per_rank(r) + 1
         mesh%halo%recv_ghosts(recv_offset_per_rank(r)) = g
         recv_global_ids   (recv_offset_per_rank(r)) = mesh%cell_global_id(g)
      end do

      ! 5) Build send CSR. We need to know which LOCAL cells to send to each rank.
      !    Use MPI_Alltoallv: send recv_global_ids → other ranks receive as "send list".
      allocate(mesh%halo%send_offset(np + 1))
      mesh%halo%send_offset(1) = 1
      do k = 1, np
         r = mesh%halo%neighbor_ranks(k)
         mesh%halo%send_offset(k+1) = mesh%halo%send_offset(k) + send_count(r)
      end do
      allocate(send_global_ids(mesh%halo%send_offset(np+1) - 1))

      ! Build sendcounts/sdispls / recvcounts/rdispls in a full-nproc layout for Alltoallv.
      block
         integer, allocatable :: sc(:), sd(:), rc(:), rd(:)
         integer :: r2
         allocate(sc(nproc), sd(nproc), rc(nproc), rd(nproc))
         sc = 0; sd = 0; rc = 0; rd = 0
         ! In this Alltoallv, this rank SENDS its want-list (recv_global_ids) to owners.
         ! That means from this rank's PoV: send sizes = recv_count (we send to owner the
         ! cells we want), receive sizes = send_count (owners tell us cells they will send).
         do r2 = 0, nproc - 1
            sc(r2+1) = recv_count(r2)
            rc(r2+1) = send_count(r2)
         end do
         sd(1) = 0
         rd(1) = 0
         do r2 = 1, nproc - 1
            sd(r2+1) = sd(r2) + sc(r2)
            rd(r2+1) = rd(r2) + rc(r2)
         end do

         ! Build a contiguous buffer of recv_global_ids in nproc-rank order.
         block
            integer, allocatable :: full_send_buf(:), full_recv_buf(:)
            integer :: nsend, nrecv, ii
            nsend = sum(sc)
            nrecv = sum(rc)
            allocate(full_send_buf(max(nsend,1)))
            allocate(full_recv_buf(max(nrecv,1)))
            ! Pack recv_global_ids in rank order
            block
               integer :: cursor, neighbor_k
               cursor = 0
               do r2 = 0, nproc - 1
                  neighbor_k = neighbor_idx(r2)
                  if (neighbor_k > 0) then
                     ii = mesh%halo%recv_offset(neighbor_k) - 1
                     do k = 1, recv_count(r2)
                        cursor = cursor + 1
                        full_send_buf(cursor) = recv_global_ids(ii + k)
                     end do
                  end if
               end do
            end block

            call MPI_Alltoallv(full_send_buf, sc, sd, MPI_INTEGER, &
                               full_recv_buf, rc, rd, MPI_INTEGER, ctx%comm)

            ! Unpack full_recv_buf back into send_global_ids in neighbor-list order.
            block
               integer :: cursor, neighbor_k
               cursor = 0
               do r2 = 0, nproc - 1
                  neighbor_k = neighbor_idx(r2)
                  if (neighbor_k > 0) then
                     ii = mesh%halo%send_offset(neighbor_k) - 1
                     do k = 1, send_count(r2)
                        cursor = cursor + 1
                        send_global_ids(ii + k) = full_recv_buf(cursor)
                     end do
                  end if
               end do
            end block

            deallocate(full_send_buf, full_recv_buf)
         end block
         deallocate(sc, sd, rc, rd)
      end block

      ! 6) Translate send_global_ids to local cell indices.
      allocate(mesh%halo%send_cells(mesh%halo%send_offset(np+1) - 1))
      block
         integer, allocatable :: glob_to_local(:)
         integer :: nc_glob, c, lid
         nc_glob = maxval(mesh%cell_global_id)
         allocate(glob_to_local(nc_glob))
         glob_to_local = 0
         do c = 1, mesh%nc_internal
            glob_to_local(mesh%cell_global_id(c)) = c
         end do
         do k = 1, size(send_global_ids)
            lid = glob_to_local(send_global_ids(k))
            if (lid == 0) then
               call mpi_abort_msg(ctx, 'halo: requested global cell not local', 99)
            end if
            mesh%halo%send_cells(k) = lid
         end do
         deallocate(glob_to_local)
      end block

      ! 7) Allocate packed buffers. Channel U: NVAR per cell. Channel GP:
      !    gradients (3*NPRIM) + limiters (NPRIM) = 4*NPRIM per cell.
      block
         integer, parameter :: GP_SIZE = 4 * NPRIM
         integer :: ntotal_send, ntotal_recv
         ntotal_send = mesh%halo%send_offset(np+1) - 1
         ntotal_recv = mesh%halo%recv_offset(np+1) - 1
         allocate(mesh%halo%send_buf   (NVAR    * ntotal_send))
         allocate(mesh%halo%recv_buf   (NVAR    * ntotal_recv))
         allocate(mesh%halo%send_buf_gp(GP_SIZE * ntotal_send))
         allocate(mesh%halo%recv_buf_gp(GP_SIZE * ntotal_recv))
      end block

      ! Persistent request init happens in halo_exchange:halo_init_persistent
      ! (needs the matching MPI_Send_init / MPI_Recv_init wired to the bufs above).
      mesh%halo%persistent_inited = .false.

      deallocate(recv_count, send_count, recv_global_ids, send_global_ids)
      deallocate(recv_offset_per_rank, neighbor_idx)
   end subroutine build_halo_desc

   subroutine free_mesh(m)
      type(t_mesh), intent(inout) :: m
      if (allocated(m%xv))            deallocate(m%xv)
      if (allocated(m%cell_type))     deallocate(m%cell_type)
      if (allocated(m%cell_vtx_ptr))  deallocate(m%cell_vtx_ptr)
      if (allocated(m%cell_vtx))      deallocate(m%cell_vtx)
      if (allocated(m%cell_volume))   deallocate(m%cell_volume)
      if (allocated(m%cell_centroid)) deallocate(m%cell_centroid)
      if (allocated(m%cell_perm))     deallocate(m%cell_perm)
      if (allocated(m%face_owner))    deallocate(m%face_owner)
      if (allocated(m%face_owner_lf)) deallocate(m%face_owner_lf)
      if (allocated(m%face_neighbor)) deallocate(m%face_neighbor)
      if (allocated(m%face_area))     deallocate(m%face_area)
      if (allocated(m%face_normal))   deallocate(m%face_normal)
      if (allocated(m%face_centroid)) deallocate(m%face_centroid)
      if (allocated(m%face_patch))    deallocate(m%face_patch)
      if (allocated(m%patches))       deallocate(m%patches)
      if (allocated(m%cell_vol_eff))     deallocate(m%cell_vol_eff)
      if (allocated(m%face_area_eff))    deallocate(m%face_area_eff)
      if (allocated(m%cell_frag_offset)) deallocate(m%cell_frag_offset)
      if (allocated(m%frag_normal))      deallocate(m%frag_normal)
      if (allocated(m%frag_area))        deallocate(m%frag_area)
      if (allocated(m%frag_centroid))    deallocate(m%frag_centroid)
   end subroutine free_mesh

end module partition
