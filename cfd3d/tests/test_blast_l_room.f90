! Arbitrary-shaped-room capability smoke test (M5.4 follow-up).
!
! A 5 g TNT-equivalent point energy is detonated at the far end of the
! top arm of an L-shaped 1 m^3 room. The bottom arm contains a 0.15 m
! rigid cube obstacle. A 0.2 m x 0.2 m vent on the x_max face (in the
! bottom arm) bleeds the chamber to ambient.
!
! Geometry capability check: probes record p(t) in the same arm as
! the charge, at the L's inside corner, and at the far end of the
! opposite arm — confirming the shock both reflects off walls and
! diffracts around the bend.
!
! Pass criteria (smoke-test scope; full validation is a follow-up):
!   1. top_far peak overpressure > 30 kPa — direct shock + back-wall
!      reflection in the charge arm.
!   2. corner  peak overpressure > 30 kPa — shock arrives at the L
!      inside corner (still strong here, having only diffracted around
!      a quarter-turn).
!   3. bot_far peak overpressure >  5 kPa — shock turns the corner
!      and reaches the far end of the opposite arm.
program test_blast_l_room
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mpi_f08
   use kinds,            only : wp
   use constants,        only : IP_P
   use mesh_types,       only : t_mesh
   use partition,        only : partition_and_load
   use halo_exchange,    only : halo_init_persistent, halo_free_persistent, &
                                halo_init_persistent_gp, halo_free_persistent_gp
   use fields,           only : t_state, alloc_state, free_state
   use bc_types,         only : t_bc_data
   use time_integration, only : compute_dt, rk3_step
   use solver_control,   only : t_run_params, read_namelist, sgs_string_to_int
   use solver_driver,    only : assign_patch_bcs, set_initial_condition
   use mpi_runtime,      only : t_mpi_ctx, mpi_init_ctx, mpi_finalize_ctx
   use gas_properties,   only : init_gas
   use probes,           only : t_probe_set, probes_load, probes_locate, &
                                probes_sample, probes_free
   implicit none

   character(len=512) :: nml_path
   integer :: nargs
   type(t_run_params) :: p
   type(t_mesh)       :: mesh
   type(t_state)      :: s
   type(t_bc_data), allocatable :: bc_dat(:)
   type(t_mpi_ctx)    :: ctx
   type(t_probe_set)  :: probes
   real(wp) :: t, dt
   integer  :: step, i, k, fails

   real(wp), parameter :: TOPFAR_MIN_DP_KPA  = 30.0_wp
   real(wp), parameter :: CORNER_MIN_DP_KPA  = 30.0_wp
   real(wp), parameter :: BOTFAR_MIN_DP_KPA  = 5.0_wp

   real(wp) :: p_peak_local, p_peak_global, dp_peak_kpa
   real(wp) :: dp_top_kpa, dp_corner_kpa, dp_bot_kpa
   real(wp) :: p_amb

   call mpi_init_ctx(ctx)

   nargs = command_argument_count()
   if (nargs < 1) then
      if (ctx%is_root) write(error_unit,'(A)') 'usage: test_blast_l_room <input.nml>'
      call mpi_finalize_ctx()
      stop 1
   end if
   call get_command_argument(1, nml_path)

   call read_namelist(trim(nml_path), p)
   call init_gas(p%R_gas, p%sutherland, p%mu_const, p%mu_ref, p%T_ref, p%S_S, p%Pr, &
                 sgs_kind=sgs_string_to_int(p%sgs_model), &
                 C_s=p%C_s, C_w=p%C_w, Pr_t=p%Pr_t)
   call partition_and_load(trim(p%mesh_file), ctx, mesh)
   call assign_patch_bcs(p, mesh, bc_dat)
   call alloc_state(s, mesh)
   call set_initial_condition(p, mesh, s, ctx)
   call halo_init_persistent(mesh, ctx)
   call halo_init_persistent_gp(mesh, ctx)
   call probes_load(trim(p%probes_file), trim(p%case_name), trim(p%output_dir), ctx, probes)
   call probes_locate(mesh, ctx, probes)

   if (ctx%is_root) write(*,'(A,F8.4,A,F8.4,A,1PE12.4)') &
      'blast_l_room: tnt_mass=', p%tnt_mass, ' kg  blast_radius=', p%blast_radius, &
      ' m  E_total=', p%tnt_mass * p%tnt_specific_E

   t = 0.0_wp
   step = 0
   call probes_sample(mesh, s, t, ctx, probes)
   do while (t < p%t_end .and. step < p%max_steps)
      dt = compute_dt(mesh, s, p%cfl, ctx, p%viscous_enabled)
      if (t + dt > p%t_end) dt = p%t_end - t
      call rk3_step(mesh, bc_dat, s, dt, ctx, p%muscl_enabled, p%viscous_enabled, p%venkat_K)
      t = t + dt
      step = step + 1
      call probes_sample(mesh, s, t, ctx, probes)
      if (ctx%is_root .and. mod(step, 200) == 0) &
         write(*,'(A,I7,A,1PE12.5,A,1PE12.5)') 'step ', step, '  t=', t, '  dt=', dt
   end do

   p_amb         = p%ambient_p
   dp_top_kpa    = -huge(1.0_wp)
   dp_corner_kpa = -huge(1.0_wp)
   dp_bot_kpa    = -huge(1.0_wp)

   do i = 1, probes%nprobes
      p_peak_local  = -huge(1.0_wp)
      p_peak_global = -huge(1.0_wp)
      if (probes%pr(i)%owner_rank == ctx%rank) then
         do k = 1, probes%nsamples
            if (probes%w_hist(IP_P, k, i) > p_peak_local) &
               p_peak_local = probes%w_hist(IP_P, k, i)
         end do
      end if
      if (ctx%nproc > 1) then
         call MPI_Allreduce(p_peak_local, p_peak_global, 1, &
                            MPI_DOUBLE_PRECISION, MPI_MAX, ctx%comm)
      else
         p_peak_global = p_peak_local
      end if
      dp_peak_kpa = (p_peak_global - p_amb) * 1.0e-3_wp

      if (ctx%is_root) write(*,'(A,A,A,F8.2,A)') &
         '  probe ', trim(probes%pr(i)%name), &
         '  Δp_peak = ', dp_peak_kpa, ' kPa'

      select case (trim(probes%pr(i)%name))
      case ('top_far'); dp_top_kpa    = dp_peak_kpa
      case ('corner');  dp_corner_kpa = dp_peak_kpa
      case ('bot_far'); dp_bot_kpa    = dp_peak_kpa
      end select
   end do

   fails = 0
   if (dp_top_kpa < TOPFAR_MIN_DP_KPA) then
      if (ctx%is_root) write(error_unit,'(A,F8.2,A,F8.2,A)') &
         'FAIL: top_far Δp_peak = ', dp_top_kpa, &
         ' kPa < threshold ', TOPFAR_MIN_DP_KPA, ' kPa'
      fails = fails + 1
   end if
   if (dp_corner_kpa < CORNER_MIN_DP_KPA) then
      if (ctx%is_root) write(error_unit,'(A,F8.2,A,F8.2,A)') &
         'FAIL: corner Δp_peak = ', dp_corner_kpa, &
         ' kPa < threshold ', CORNER_MIN_DP_KPA, ' kPa'
      fails = fails + 1
   end if
   if (dp_bot_kpa < BOTFAR_MIN_DP_KPA) then
      if (ctx%is_root) write(error_unit,'(A,F8.2,A,F8.2,A)') &
         'FAIL: bot_far Δp_peak = ', dp_bot_kpa, &
         ' kPa < threshold ', BOTFAR_MIN_DP_KPA, ' kPa (shock did not turn the L corner)'
      fails = fails + 1
   end if

   call probes_free(probes)
   call halo_free_persistent_gp(mesh)
   call halo_free_persistent(mesh)
   call free_state(s)
   if (allocated(bc_dat)) deallocate(bc_dat)

   if (ctx%is_root) then
      if (fails == 0) then
         write(*,'(A)') 'test_blast_l_room: PASS'
      else
         write(*,'(A,I0,A)') 'test_blast_l_room: FAIL (', fails, ')'
      end if
   end if
   call mpi_finalize_ctx()
   if (fails /= 0) stop 1
end program test_blast_l_room
