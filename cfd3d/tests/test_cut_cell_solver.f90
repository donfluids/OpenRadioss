! Phase 1e end-to-end cut-cell test: conservation in a closed box with
! an embedded non-grid-aligned cube obstacle.
!
! All six domain faces and the cube surface are slip walls, so the box
! is closed and adiabatic. Total mass (Σ ρ V_eff) and total energy
! (Σ ρE V_eff) must therefore be conserved to machine precision as the
! blast waves bounce off the obstacle and the walls. Any leak through
! the cut faces or the obstacle fragments would show up as drift.
!
! This exercises the full cut-cell path: build_cut_cell_tables with a
! namelist cube (cut cells + dead interior cells + fragments), merge
! groups for sliver stability, merge_sync_state, and the V_eff / A_eff /
! fragment-flux time integration.
program test_cut_cell_solver
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mpi_f08
   use kinds,            only : wp
   use constants,        only : IRHO, IRHOE
   use mesh_types,       only : t_mesh
   use partition,        only : partition_and_load
   use cut_cell,         only : build_cut_cell_tables, merge_sync_state, V_DEAD_FRAC
   use halo_exchange,    only : halo_init_persistent, halo_free_persistent, &
                                halo_init_persistent_gp, halo_free_persistent_gp
   use fields,           only : t_state, alloc_state, free_state
   use bc_types,         only : t_bc_data
   use time_integration, only : compute_dt, rk3_step
   use solver_control,   only : t_run_params, read_namelist, sgs_string_to_int
   use solver_driver,    only : assign_patch_bcs, set_initial_condition
   use mpi_runtime,      only : t_mpi_ctx, mpi_init_ctx, mpi_finalize_ctx
   use gas_properties,   only : init_gas
   use eos_ideal_gas,    only : prim_from_cons
   implicit none

   character(len=512) :: nml_path
   type(t_run_params) :: p
   type(t_mesh)       :: mesh
   type(t_state)      :: s
   type(t_bc_data), allocatable :: bc_dat(:)
   type(t_mpi_ctx)    :: ctx
   real(wp) :: t, dt
   integer  :: step, c, fails, nobs
   integer  :: n_cut, n_dead, n_full
   real(wp) :: frac
   real(wp) :: mass0, ener0, mass1, ener1, dmass, dener
   real(wp) :: rho, u, v, w, pr, pmin

   real(wp), parameter :: TOL_CONS = 1.0e-10_wp

   call mpi_init_ctx(ctx)
   if (command_argument_count() < 1) then
      if (ctx%is_root) write(error_unit,'(A)') 'usage: test_cut_cell_solver <input.nml>'
      call mpi_finalize_ctx(); stop 1
   end if
   call get_command_argument(1, nml_path)

   call read_namelist(trim(nml_path), p)
   call init_gas(p%R_gas, p%sutherland, p%mu_const, p%mu_ref, p%T_ref, p%S_S, p%Pr, &
                 sgs_kind=sgs_string_to_int(p%sgs_model), &
                 C_s=p%C_s, C_w=p%C_w, Pr_t=p%Pr_t)
   call partition_and_load(trim(p%mesh_file), ctx, mesh)
   nobs = max(p%n_obstacles, 1)
   call build_cut_cell_tables(mesh, p%n_obstacles, &
        p%obstacle_cube_lo(:, 1:nobs), p%obstacle_cube_hi(:, 1:nobs))
   call assign_patch_bcs(p, mesh, bc_dat)
   call alloc_state(s, mesh)
   call set_initial_condition(p, mesh, s, ctx)
   call merge_sync_state(mesh, s)
   call halo_init_persistent(mesh, ctx)
   call halo_init_persistent_gp(mesh, ctx)

   ! Cut-cell census (serial test — just rank 0's local cells).
   n_cut = 0; n_dead = 0; n_full = 0
   do c = 1, mesh%nc_internal
      frac = mesh%cell_vol_eff(c) / mesh%cell_volume(c)
      if (frac <= V_DEAD_FRAC) then
         n_dead = n_dead + 1
      else if (frac < 1.0_wp - 1.0e-12_wp) then
         n_cut = n_cut + 1
      else
         n_full = n_full + 1
      end if
   end do
   if (ctx%is_root) write(*,'(A,I0,A,I0,A,I0,A,I0)') &
      'cut_cell census: full=', n_full, ' cut=', n_cut, &
      ' dead=', n_dead, ' fragments=', mesh%n_frags

   call totals(mesh, s, mass0, ener0)
   if (ctx%is_root) write(*,'(A,1PE20.12,A,1PE20.12)') &
      'initial  mass=', mass0, '  energy=', ener0

   t = 0.0_wp; step = 0
   do while (t < p%t_end .and. step < p%max_steps)
      dt = compute_dt(mesh, s, p%cfl, ctx, p%viscous_enabled)
      if (t + dt > p%t_end) dt = p%t_end - t
      call rk3_step(mesh, bc_dat, s, dt, ctx, p%muscl_enabled, p%viscous_enabled, p%venkat_K)
      t = t + dt; step = step + 1
      if (ctx%is_root .and. mod(step, 100) == 0) &
         write(*,'(A,I7,A,1PE12.5)') 'step ', step, '  t=', t
   end do

   call totals(mesh, s, mass1, ener1)
   if (ctx%is_root) write(*,'(A,1PE20.12,A,1PE20.12)') &
      'final    mass=', mass1, '  energy=', ener1

   ! Minimum pressure over live cells (stability / positivity check).
   pmin = huge(1.0_wp)
   do c = 1, mesh%nc_internal
      if (mesh%cell_vol_eff(c) <= V_DEAD_FRAC * mesh%cell_volume(c)) cycle
      call prim_from_cons(s%U(:, c), rho, u, v, w, pr)
      if (pr < pmin) pmin = pr
   end do

   dmass = abs(mass1 - mass0) / mass0
   dener = abs(ener1 - ener0) / ener0
   if (ctx%is_root) write(*,'(A,1PE10.3,A,1PE10.3,A,1PE12.5)') &
      'rel drift: mass=', dmass, '  energy=', dener, '  pmin=', pmin

   fails = 0
   if (n_cut == 0) then
      if (ctx%is_root) write(error_unit,'(A)') 'FAIL: no cut cells produced (cube not engaged)'
      fails = fails + 1
   end if
   if (mesh%n_frags == 0) then
      if (ctx%is_root) write(error_unit,'(A)') 'FAIL: no obstacle fragments produced'
      fails = fails + 1
   end if
   if (dmass > TOL_CONS) then
      if (ctx%is_root) write(error_unit,'(A,1PE10.3)') 'FAIL: mass not conserved, rel drift=', dmass
      fails = fails + 1
   end if
   if (dener > TOL_CONS) then
      if (ctx%is_root) write(error_unit,'(A,1PE10.3)') 'FAIL: energy not conserved, rel drift=', dener
      fails = fails + 1
   end if
   if (.not. (pmin > 0.0_wp)) then
      if (ctx%is_root) write(error_unit,'(A,1PE12.5)') 'FAIL: non-positive pressure pmin=', pmin
      fails = fails + 1
   end if

   call halo_free_persistent_gp(mesh)
   call halo_free_persistent(mesh)
   call free_state(s)
   if (allocated(bc_dat)) deallocate(bc_dat)

   if (ctx%is_root) then
      if (fails == 0) then
         write(*,'(A)') 'test_cut_cell_solver: PASS'
      else
         write(*,'(A,I0,A)') 'test_cut_cell_solver: FAIL (', fails, ')'
      end if
   end if
   call mpi_finalize_ctx()
   if (fails /= 0) stop 1

contains

   subroutine totals(mesh, s, mass, ener)
      type(t_mesh),  intent(in)  :: mesh
      type(t_state), intent(in)  :: s
      real(wp),      intent(out) :: mass, ener
      integer  :: cc
      real(wp) :: ml, el, mg, eg
      ml = 0.0_wp; el = 0.0_wp
      do cc = 1, mesh%nc_internal
         ml = ml + s%U(IRHO,  cc) * mesh%cell_vol_eff(cc)
         el = el + s%U(IRHOE, cc) * mesh%cell_vol_eff(cc)
      end do
      if (ctx%nproc > 1) then
         call MPI_Allreduce(ml, mg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ctx%comm)
         call MPI_Allreduce(el, eg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ctx%comm)
         mass = mg; ener = eg
      else
         mass = ml; ener = el
      end if
   end subroutine totals

end program test_cut_cell_solver
