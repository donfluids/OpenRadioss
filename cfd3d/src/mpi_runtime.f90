module mpi_runtime
   use mpi_f08
   implicit none
   private

   public :: t_mpi_ctx, mpi_init_ctx, mpi_finalize_ctx, mpi_abort_msg

   type :: t_mpi_ctx
      type(MPI_Comm) :: comm = MPI_COMM_WORLD
      integer :: rank   = 0
      integer :: nproc  = 1
      logical :: is_root = .true.
   end type t_mpi_ctx

contains

   subroutine mpi_init_ctx(ctx)
      type(t_mpi_ctx), intent(out) :: ctx
      logical :: initialized
      call MPI_Initialized(initialized)
      if (.not. initialized) call MPI_Init()
      ctx%comm = MPI_COMM_WORLD
      call MPI_Comm_rank(ctx%comm, ctx%rank)
      call MPI_Comm_size(ctx%comm, ctx%nproc)
      ctx%is_root = (ctx%rank == 0)
   end subroutine mpi_init_ctx

   subroutine mpi_finalize_ctx()
      logical :: finalized
      call MPI_Finalized(finalized)
      if (.not. finalized) call MPI_Finalize()
   end subroutine mpi_finalize_ctx

   subroutine mpi_abort_msg(ctx, msg, code)
      type(t_mpi_ctx),  intent(in) :: ctx
      character(len=*), intent(in) :: msg
      integer,          intent(in) :: code
      write(*,'(A,I0,A,A)') 'rank ', ctx%rank, ' abort: ', trim(msg)
      call MPI_Abort(ctx%comm, code)
   end subroutine mpi_abort_msg

end module mpi_runtime
