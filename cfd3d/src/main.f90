program cfd3d
   use solver_control, only : t_run_params, read_namelist
   use solver_driver,  only : run_case
   implicit none

   type(t_run_params) :: p
   character(len=512) :: arg
   integer :: nargs

   nargs = command_argument_count()
   if (nargs < 1) then
      write(*,'(A)') 'usage: cfd3d_solver <input.nml>'
      stop 1
   end if
   call get_command_argument(1, arg)

   call read_namelist(trim(arg), p)
   call run_case(p)
end program cfd3d
