! Persistent-request halo exchange across rank boundaries.
!
! Two channels:
!   Channel U  — NVAR floats per cell (conservative state).
!   Channel GP — 4*NPRIM floats per cell: 3*NPRIM gradients of primitives,
!                followed by NPRIM Venkat limiters.
!
! Both channels share the send_cells / recv_ghosts CSR pattern built by
! `partition.build_halo_desc`. Each channel has its own packed buffers and
! persistent MPI request handles.
module halo_exchange
   use kinds,       only : wp
   use constants,   only : NVAR, NPRIM
   use mpi_f08
   use mpi_runtime, only : t_mpi_ctx
   use mesh_types,  only : t_mesh
   use fields,      only : t_state
   implicit none
   private

   public :: halo_init_persistent, halo_free_persistent
   public :: halo_pack_and_start, halo_wait
   public :: halo_init_persistent_gp, halo_free_persistent_gp
   public :: halo_pack_and_start_gp, halo_wait_gp

   integer, parameter :: HALO_TAG_U  = 7777
   integer, parameter :: HALO_TAG_GP = 7778
   integer, parameter :: GP_SIZE     = 4 * NPRIM     ! 3 grad + 1 ψ per primitive

contains

   ! --- Channel U (conservative state) ---------------------------------------

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
         if (ns > 0) call MPI_Send_init( &
            mesh%halo%send_buf(send_off : send_off + ns*NVAR - 1), &
            ns*NVAR, MPI_DOUBLE_PRECISION, r, HALO_TAG_U, ctx%comm, mesh%halo%send_req(k))
         if (nr > 0) call MPI_Recv_init( &
            mesh%halo%recv_buf(recv_off : recv_off + nr*NVAR - 1), &
            nr*NVAR, MPI_DOUBLE_PRECISION, r, HALO_TAG_U, ctx%comm, mesh%halo%recv_req(k))
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

   subroutine halo_pack_and_start(mesh, U, ctx)
      type(t_mesh),    intent(inout) :: mesh
      real(wp),        intent(in)    :: U(:, :)         ! (NVAR, nc_total)
      type(t_mpi_ctx), intent(in)    :: ctx
      integer :: k, np, cell_idx, var
      np = mesh%halo%n_neighbors
      if (np == 0) return
      do k = 1, size(mesh%halo%send_cells)
         cell_idx = mesh%halo%send_cells(k)
         do var = 1, NVAR
            mesh%halo%send_buf((k - 1) * NVAR + var) = U(var, cell_idx)
         end do
      end do
      do k = 1, np
         if (mesh%halo%recv_offset(k+1) > mesh%halo%recv_offset(k)) &
            call MPI_Start(mesh%halo%recv_req(k))
      end do
      do k = 1, np
         if (mesh%halo%send_offset(k+1) > mesh%halo%send_offset(k)) &
            call MPI_Start(mesh%halo%send_req(k))
      end do
      if (.false.) k = ctx%rank
   end subroutine halo_pack_and_start

   subroutine halo_wait(mesh, U)
      type(t_mesh), intent(inout) :: mesh
      real(wp),     intent(inout) :: U(:, :)
      integer :: k, np, ghost_idx, var
      np = mesh%halo%n_neighbors
      if (np == 0) return
      do k = 1, np
         if (mesh%halo%recv_offset(k+1) > mesh%halo%recv_offset(k)) &
            call MPI_Wait(mesh%halo%recv_req(k), MPI_STATUS_IGNORE)
      end do
      do k = 1, np
         if (mesh%halo%send_offset(k+1) > mesh%halo%send_offset(k)) &
            call MPI_Wait(mesh%halo%send_req(k), MPI_STATUS_IGNORE)
      end do
      do k = 1, size(mesh%halo%recv_ghosts)
         ghost_idx = mesh%halo%recv_ghosts(k)
         do var = 1, NVAR
            U(var, ghost_idx) = mesh%halo%recv_buf((k - 1) * NVAR + var)
         end do
      end do
   end subroutine halo_wait

   ! --- Channel GP (gradients + limiters) ------------------------------------

   subroutine halo_init_persistent_gp(mesh, ctx)
      type(t_mesh),    intent(inout) :: mesh
      type(t_mpi_ctx), intent(in)    :: ctx
      integer :: k, np, r, ns, nr, send_off, recv_off

      np = mesh%halo%n_neighbors
      if (np == 0) then
         mesh%halo%persistent_inited_gp = .true.
         return
      end if
      allocate(mesh%halo%send_req_gp(np))
      allocate(mesh%halo%recv_req_gp(np))
      do k = 1, np
         r = mesh%halo%neighbor_ranks(k)
         ns = mesh%halo%send_offset(k+1) - mesh%halo%send_offset(k)
         nr = mesh%halo%recv_offset(k+1) - mesh%halo%recv_offset(k)
         send_off = (mesh%halo%send_offset(k) - 1) * GP_SIZE + 1
         recv_off = (mesh%halo%recv_offset(k) - 1) * GP_SIZE + 1
         if (ns > 0) call MPI_Send_init( &
            mesh%halo%send_buf_gp(send_off : send_off + ns*GP_SIZE - 1), &
            ns*GP_SIZE, MPI_DOUBLE_PRECISION, r, HALO_TAG_GP, ctx%comm, mesh%halo%send_req_gp(k))
         if (nr > 0) call MPI_Recv_init( &
            mesh%halo%recv_buf_gp(recv_off : recv_off + nr*GP_SIZE - 1), &
            nr*GP_SIZE, MPI_DOUBLE_PRECISION, r, HALO_TAG_GP, ctx%comm, mesh%halo%recv_req_gp(k))
      end do
      mesh%halo%persistent_inited_gp = .true.
   end subroutine halo_init_persistent_gp

   subroutine halo_free_persistent_gp(mesh)
      type(t_mesh), intent(inout) :: mesh
      integer :: k
      if (.not. mesh%halo%persistent_inited_gp) return
      if (allocated(mesh%halo%send_req_gp)) then
         do k = 1, size(mesh%halo%send_req_gp)
            call MPI_Request_free(mesh%halo%send_req_gp(k))
         end do
         deallocate(mesh%halo%send_req_gp)
      end if
      if (allocated(mesh%halo%recv_req_gp)) then
         do k = 1, size(mesh%halo%recv_req_gp)
            call MPI_Request_free(mesh%halo%recv_req_gp(k))
         end do
         deallocate(mesh%halo%recv_req_gp)
      end if
      mesh%halo%persistent_inited_gp = .false.
   end subroutine halo_free_persistent_gp

   ! Pack layout for one cell: [∇W_rho(3), ∇W_u(3), ∇W_v(3), ∇W_w(3), ∇W_p(3),
   !                            ψ_rho, ψ_u, ψ_v, ψ_w, ψ_p]
   subroutine halo_pack_and_start_gp(mesh, s, ctx)
      type(t_mesh),    intent(inout) :: mesh
      type(t_state),   intent(in)    :: s
      type(t_mpi_ctx), intent(in)    :: ctx
      integer :: k, np, cell_idx, v, off
      np = mesh%halo%n_neighbors
      if (np == 0) return
      do k = 1, size(mesh%halo%send_cells)
         cell_idx = mesh%halo%send_cells(k)
         off = (k - 1) * GP_SIZE
         do v = 1, NPRIM
            mesh%halo%send_buf_gp(off + (v-1)*3 + 1) = s%gradW(1, v, cell_idx)
            mesh%halo%send_buf_gp(off + (v-1)*3 + 2) = s%gradW(2, v, cell_idx)
            mesh%halo%send_buf_gp(off + (v-1)*3 + 3) = s%gradW(3, v, cell_idx)
         end do
         do v = 1, NPRIM
            mesh%halo%send_buf_gp(off + 3*NPRIM + v) = s%psi(v, cell_idx)
         end do
      end do
      do k = 1, np
         if (mesh%halo%recv_offset(k+1) > mesh%halo%recv_offset(k)) &
            call MPI_Start(mesh%halo%recv_req_gp(k))
      end do
      do k = 1, np
         if (mesh%halo%send_offset(k+1) > mesh%halo%send_offset(k)) &
            call MPI_Start(mesh%halo%send_req_gp(k))
      end do
      if (.false.) k = ctx%rank
   end subroutine halo_pack_and_start_gp

   subroutine halo_wait_gp(mesh, s)
      type(t_mesh),  intent(inout) :: mesh
      type(t_state), intent(inout) :: s
      integer :: k, np, ghost_idx, v, off
      np = mesh%halo%n_neighbors
      if (np == 0) return
      do k = 1, np
         if (mesh%halo%recv_offset(k+1) > mesh%halo%recv_offset(k)) &
            call MPI_Wait(mesh%halo%recv_req_gp(k), MPI_STATUS_IGNORE)
      end do
      do k = 1, np
         if (mesh%halo%send_offset(k+1) > mesh%halo%send_offset(k)) &
            call MPI_Wait(mesh%halo%send_req_gp(k), MPI_STATUS_IGNORE)
      end do
      do k = 1, size(mesh%halo%recv_ghosts)
         ghost_idx = mesh%halo%recv_ghosts(k)
         off = (k - 1) * GP_SIZE
         do v = 1, NPRIM
            s%gradW(1, v, ghost_idx) = mesh%halo%recv_buf_gp(off + (v-1)*3 + 1)
            s%gradW(2, v, ghost_idx) = mesh%halo%recv_buf_gp(off + (v-1)*3 + 2)
            s%gradW(3, v, ghost_idx) = mesh%halo%recv_buf_gp(off + (v-1)*3 + 3)
         end do
         do v = 1, NPRIM
            s%psi(v, ghost_idx) = mesh%halo%recv_buf_gp(off + 3*NPRIM + v)
         end do
      end do
   end subroutine halo_wait_gp

end module halo_exchange
