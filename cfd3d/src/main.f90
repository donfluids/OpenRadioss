program cfd3d
   use solver_control, only : t_run_params, read_namelist
   use solver_driver,  only : run_case
   use mpi_runtime,    only : t_mpi_ctx, mpi_init_ctx, mpi_finalize_ctx
   implicit none

   type(t_run_params) :: p
   type(t_mpi_ctx)    :: ctx
   character(len=512) :: arg
   integer :: nargs

   call mpi_init_ctx(ctx)

   nargs = command_argument_count()
   if (nargs < 1) then
      if (ctx%is_root) write(*,'(A)') 'usage: cfd3d_solver <input.nml>'
      call mpi_finalize_ctx()
      stop 1
   end if
   call get_command_argument(1, arg)

   call read_namelist(trim(arg), p)
   call run_case(p, ctx)

   call mpi_finalize_ctx()
end program cfd3d
