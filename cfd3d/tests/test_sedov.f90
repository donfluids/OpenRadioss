! Sedov-Taylor self-similar blast verification.
!
! Deposit total energy E into a small sphere at blast_center. The analytical
! self-similar solution gives the shock-front radius at time t as
!     R(t) = ξ₀ · (E / ρ_∞)^{1/5} · t^{2/5}
! with ξ₀ ≈ 1.033 in 3-D for γ = 1.4 (Sedov 1959 / Landau-Lifshitz §99).
!
! We run the case to t_end and compare the radius at which the cell-centred
! pressure is maximum to the analytical value. Acceptance: within 10 %.
program test_sedov
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mpi_f08
   use kinds,            only : wp
   use mesh_types,       only : t_mesh
   use partition,        only : partition_and_load
   use halo_exchange,    only : halo_init_persistent, halo_free_persistent, &
                                halo_init_persistent_gp, halo_free_persistent_gp
   use fields,           only : t_state, alloc_state, free_state
   use bc_types,         only : t_bc_data
   use eos_ideal_gas,    only : prim_from_cons
   use time_integration, only : compute_dt, rk3_step
   use solver_control,   only : t_run_params, read_namelist, sgs_string_to_int
   use solver_driver,    only : assign_patch_bcs, set_initial_condition
   use mpi_runtime,      only : t_mpi_ctx, mpi_init_ctx, mpi_finalize_ctx
   use gas_properties,   only : init_gas
   implicit none

   character(len=512) :: nml_path
   integer :: nargs
   type(t_run_params) :: p
   type(t_mesh)       :: mesh
   type(t_state)      :: s
   type(t_bc_data), allocatable :: bc_dat(:)
   type(t_mpi_ctx)    :: ctx
   real(wp) :: t, dt
   real(wp) :: rho, ux, uy, uz, pres
   real(wp) :: dx, dy, dz, r
   real(wp) :: r_at_peak_local, r_at_peak_global, R_analytic, err
   integer  :: c, step, fails
   ! MPI MAXLOC pair (pressure, rank-encoded r metadata) — we cheat by sending
   ! r as the value to be maxed-on-rank-of-max-pressure: see below.
   type :: dpair
      real(wp) :: x
      integer  :: rk
   end type dpair
   real(wp) :: p_peak_local, p_peak_global
   type(dpair) :: pl, pg

   call mpi_init_ctx(ctx)

   nargs = command_argument_count()
   if (nargs < 1) then
      if (ctx%is_root) write(error_unit,'(A)') 'usage: test_sedov <input.nml>'
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

   t = 0.0_wp
   step = 0
   do while (t < p%t_end .and. step < p%max_steps)
      dt = compute_dt(mesh, s, p%cfl, ctx, p%viscous_enabled)
      if (t + dt > p%t_end) dt = p%t_end - t
      call rk3_step(mesh, bc_dat, s, dt, ctx, p%muscl_enabled, p%viscous_enabled, p%venkat_K)
      t = t + dt
      step = step + 1
   end do

   ! Find the cell with maximum cell-centred pressure on this rank; record
   ! its radius from the blast centre. Then MPI MAXLOC to find the
   ! globally-maximum-pressure rank, and that rank broadcasts its r.
   p_peak_local    = -huge(1.0_wp)
   r_at_peak_local = 0.0_wp
   do c = 1, mesh%nc_internal
      call prim_from_cons(s%U(:, c), rho, ux, uy, uz, pres)
      if (pres > p_peak_local) then
         p_peak_local = pres
         dx = mesh%cell_centroid(1, c) - p%blast_center(1)
         dy = mesh%cell_centroid(2, c) - p%blast_center(2)
         dz = mesh%cell_centroid(3, c) - p%blast_center(3)
         r_at_peak_local = sqrt(dx*dx + dy*dy + dz*dz)
      end if
   end do

   pl%x  = p_peak_local
   pl%rk = ctx%rank
   if (ctx%nproc > 1) then
      call MPI_Allreduce(pl, pg, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, ctx%comm)
   else
      pg = pl
   end if
   p_peak_global = pg%x
   ! Broadcast the r_at_peak from the winning rank.
   if (ctx%nproc > 1) then
      if (ctx%rank == pg%rk) r_at_peak_global = r_at_peak_local
      call MPI_Bcast(r_at_peak_global, 1, MPI_DOUBLE_PRECISION, pg%rk, ctx%comm)
   else
      r_at_peak_global = r_at_peak_local
   end if

   ! Analytical Sedov shock radius (3-D, γ=1.4).
   R_analytic = 1.033_wp * (p%blast_energy / p%ambient_rho)**(0.2_wp) * t**(0.4_wp)
   err = abs(r_at_peak_global - R_analytic) / max(R_analytic, 1.0e-30_wp)

   if (ctx%is_root) then
      write(*,'(A,I0,A,1PE12.4)') 'sedov: steps=', step, '  t=', t
      write(*,'(A,1PE12.4)') 'peak pressure (global) = ', p_peak_global
      write(*,'(A,F8.4,A,F8.4)') 'R_num = ', r_at_peak_global, '   R_analytic = ', R_analytic
      write(*,'(A,F6.2,A)')      '|R_num - R_analytic| / R_analytic = ', err*100.0_wp, ' %'
   end if

   fails = 0
   if (err > 0.10_wp) then
      if (ctx%is_root) write(error_unit,'(A,F6.2)') 'sedov R-shock error exceeds 10 % : ', err*100.0_wp
      fails = fails + 1
   end if

   call halo_free_persistent_gp(mesh)
   call halo_free_persistent(mesh)
   call free_state(s)
   if (allocated(bc_dat)) deallocate(bc_dat)

   if (ctx%is_root) then
      if (fails == 0) then
         write(*,'(A)') 'test_sedov: PASS'
      else
         write(*,'(A,I0,A)') 'test_sedov: FAIL (', fails, ')'
      end if
   end if
   call mpi_finalize_ctx()
   if (fails /= 0) stop 1
end program test_sedov
