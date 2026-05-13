! Persistent-request halo exchange of cell state across rank boundaries.
!
! Usage:
!   call halo_init_persistent(mesh, ctx)             ! once after partition
!   ...
!   call halo_pack_and_start(mesh, U, ctx)           ! kick off Isend/Irecv each stage
!   call compute_residual_pure_interior(...)         ! overlap window
!   call halo_wait(mesh)
!   call compute_residual_partition(...)
!   call compute_residual_boundary(...)
!   ...
!   call halo_free_persistent(mesh)                  ! at shutdown
module halo_exchange
   use kinds,       only : wp
   use constants,   only : NVAR
   use mpi_f08
   use mpi_runtime, only : t_mpi_ctx
   use mesh_types,  only : t_mesh
   implicit none
   private

   public :: halo_init_persistent, halo_free_persistent
   public :: halo_pack_and_start, halo_wait

   integer, parameter :: HALO_TAG = 7777

contains

   subroutine halo_init_persistent(mesh, ctx)
      type(t_mesh),    intent(inout) :: mesh
      type(t_mpi_ctx), intent(in)    :: ctx

      integer :: k, np, r, ns, nr, send_off, recv_off

      np = mesh%halo%n_neighbors
      if (np == 0) then
         mesh%halo%persistent_inited = .true.
         return
      end if

      allocate(mesh%halo%send_req(np))
      allocate(mesh%halo%recv_req(np))

      do k = 1, np
         r = mesh%halo%neighbor_ranks(k)
         ns = mesh%halo%send_offset(k+1) - mesh%halo%send_offset(k)
         nr = mesh%halo%recv_offset(k+1) - mesh%halo%recv_offset(k)
         send_off = (mesh%halo%send_offset(k) - 1) * NVAR + 1
         recv_off = (mesh%halo%recv_offset(k) - 1) * NVAR + 1

         if (ns > 0) then
            call MPI_Send_init(mesh%halo%send_buf(send_off : send_off + ns*NVAR - 1), &
                               ns*NVAR, MPI_DOUBLE_PRECISION, r, HALO_TAG, &
                               ctx%comm, mesh%halo%send_req(k))
         end if
         if (nr > 0) then
            call MPI_Recv_init(mesh%halo%recv_buf(recv_off : recv_off + nr*NVAR - 1), &
                               nr*NVAR, MPI_DOUBLE_PRECISION, r, HALO_TAG, &
                               ctx%comm, mesh%halo%recv_req(k))
         end if
      end do

      mesh%halo%persistent_inited = .true.
   end subroutine halo_init_persistent

   subroutine halo_free_persistent(mesh)
      type(t_mesh), intent(inout) :: mesh
      integer :: k
      if (.not. mesh%halo%persistent_inited) return
      if (allocated(mesh%halo%send_req)) then
         do k = 1, size(mesh%halo%send_req)
            call MPI_Request_free(mesh%halo%send_req(k))
         end do
         deallocate(mesh%halo%send_req)
      end if
      if (allocated(mesh%halo%recv_req)) then
         do k = 1, size(mesh%halo%recv_req)
            call MPI_Request_free(mesh%halo%recv_req(k))
         end do
         deallocate(mesh%halo%recv_req)
      end if
      mesh%halo%persistent_inited = .false.
   end subroutine halo_free_persistent

   ! Pack the local cell state into the send buffer (in send-CSR order) and
   ! start all persistent send/recv. Receives unpack into ghost cells later
   ! (in halo_wait).
   subroutine halo_pack_and_start(mesh, U, ctx)
      type(t_mesh),    intent(inout) :: mesh
      real(wp),        intent(in)    :: U(:, :)         ! (NVAR, nc_total)
      type(t_mpi_ctx), intent(in)    :: ctx

      integer :: k, np, idx, cell_idx, var
      np = mesh%halo%n_neighbors
      if (np == 0) return

      ! Pack send_cells → send_buf
      do k = 1, size(mesh%halo%send_cells)
         cell_idx = mesh%halo%send_cells(k)
         do var = 1, NVAR
            mesh%halo%send_buf((k - 1) * NVAR + var) = U(var, cell_idx)
         end do
      end do

      ! Start receives first (good practice for matching)
      do k = 1, np
         if (mesh%halo%recv_offset(k+1) > mesh%halo%recv_offset(k)) then
            call MPI_Start(mesh%halo%recv_req(k))
         end if
      end do
      do k = 1, np
         if (mesh%halo%send_offset(k+1) > mesh%halo%send_offset(k)) then
            call MPI_Start(mesh%halo%send_req(k))
         end if
      end do

      ! Suppress unused
      if (.false.) idx = ctx%rank
   end subroutine halo_pack_and_start

   ! Wait for all pending sends/recvs and unpack received data into ghost cells.
   subroutine halo_wait(mesh, U)
      type(t_mesh), intent(inout) :: mesh
      real(wp),     intent(inout) :: U(:, :)             ! (NVAR, nc_total)
      integer :: k, np, ghost_idx, var
      np = mesh%halo%n_neighbors
      if (np == 0) return

      do k = 1, np
         if (mesh%halo%recv_offset(k+1) > mesh%halo%recv_offset(k)) then
            call MPI_Wait(mesh%halo%recv_req(k), MPI_STATUS_IGNORE)
         end if
      end do
      do k = 1, np
         if (mesh%halo%send_offset(k+1) > mesh%halo%send_offset(k)) then
            call MPI_Wait(mesh%halo%send_req(k), MPI_STATUS_IGNORE)
         end if
      end do

      ! Unpack recv_buf → ghost cells
      do k = 1, size(mesh%halo%recv_ghosts)
         ghost_idx = mesh%halo%recv_ghosts(k)
         do var = 1, NVAR
            U(var, ghost_idx) = mesh%halo%recv_buf((k - 1) * NVAR + var)
         end do
      end do
   end subroutine halo_wait

end module halo_exchange
