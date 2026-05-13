! MMS spatial-order consistency test.
!
! At t=0 we initialise U = U_a (manufactured solution) at every cell centroid.
! Then we evaluate the discrete residual R = Σ_faces F·A (inviscid HLLC + viscous)
! and subtract the analytical body source V·S(centroid). What remains is the
! truncation error of the spatial discretization on the manufactured solution.
!
! The L2 norm of (R - V·S)/V should decrease at order p with mesh refinement
! for a p-th-order spatial scheme. We run on two mesh sizes (N=8 and N=16
! cubes), compute the L2 norm on each, and report
!     p = log2( L2(N) / L2(2N) ).
program test_mms_order
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mpi_f08
   use kinds,            only : wp
   use constants,        only : NVAR, NPRIM, IP_RHO, IP_U, IP_V, IP_W, IP_P
   use mesh_types,       only : t_mesh
   use partition,        only : partition_and_load
   use halo_exchange,    only : halo_init_persistent, halo_free_persistent, &
                                halo_init_persistent_gp, halo_free_persistent_gp, &
                                halo_pack_and_start, halo_wait, &
                                halo_pack_and_start_gp, halo_wait_gp
   use fields,           only : t_state, alloc_state, free_state
   use bc_types,         only : t_bc_data
   use eos_ideal_gas,    only : cons_from_prim
   use solver_control,   only : t_run_params, read_namelist
   use solver_driver,    only : assign_patch_bcs
   use mpi_runtime,      only : t_mpi_ctx, mpi_init_ctx, mpi_finalize_ctx
   use gas_properties,   only : init_gas
   use gradients,        only : compute_primitives, compute_gradients
   use limiters,         only : compute_venkat_limiter
   use flux_assembly,    only : residual_begin, residual_pure_interior, &
                                residual_partition, residual_boundary
   use mms,              only : mms_primitives, mms_source
   implicit none

   integer :: nargs
   character(len=512) :: nml_path
   real(wp) :: L2_coarse, L2_fine, order
   type(t_mpi_ctx) :: ctx
   integer :: fails

   call mpi_init_ctx(ctx)

   nargs = command_argument_count()
   if (nargs < 1) then
      if (ctx%is_root) write(error_unit,'(A)') 'usage: test_mms_order <input.nml>'
      call mpi_finalize_ctx()
      stop 1
   end if
   call get_command_argument(1, nml_path)

   call run_one(nml_path, 'coarse', ctx, L2_coarse)
   call run_one(nml_path, 'fine',   ctx, L2_fine)

   order = log(L2_coarse / max(L2_fine, 1.0e-30_wp)) / log(2.0_wp)
   if (ctx%is_root) then
      write(*,'(A,1PE12.4,A,1PE12.4)') 'MMS L2: coarse=', L2_coarse, '  fine=', L2_fine
      write(*,'(A,F6.3)')              'observed spatial order p = ', order
   end if

   fails = 0
   ! Expect p ≈ 2 for MUSCL on smooth flow. Accept p > 1.5 (the order may be
   ! reduced by the boundary stencil's first-order behaviour at slip walls
   ! when dW=0 is enforced; the bulk dominates the L2 norm).
   if (order < 1.5_wp) then
      if (ctx%is_root) write(error_unit,'(A,F6.3)') 'MMS order below 1.5: ', order
      fails = fails + 1
   end if
   if (L2_fine > 0.5_wp * L2_coarse) then
      if (ctx%is_root) write(error_unit,'(A)') 'MMS L2 did not decrease enough with refinement'
      fails = fails + 1
   end if

   if (ctx%is_root) then
      if (fails == 0) then
         write(*,'(A)') 'test_mms_order: PASS'
      else
         write(*,'(A,I0,A)') 'test_mms_order: FAIL (', fails, ')'
      end if
   end if
   call mpi_finalize_ctx()
   if (fails /= 0) stop 1

contains

   subroutine run_one(base_nml, label, ctx, L2_global)
      character(len=*),  intent(in)    :: base_nml, label
      type(t_mpi_ctx),   intent(in)    :: ctx
      real(wp),          intent(out)   :: L2_global

      character(len=512) :: this_nml
      type(t_run_params) :: p
      type(t_mesh)       :: mesh
      type(t_state)      :: s
      type(t_bc_data), allocatable :: bc_dat(:)
      integer :: c
      real(wp) :: W(NPRIM), U_a(NVAR), Svec(NVAR)
      real(wp) :: err_local, err_global, vol_local, vol_global
      real(wp) :: rho, ux, uy, uz, pres

      this_nml = trim(base_nml) // '.' // trim(label)
      call read_namelist(trim(this_nml), p)
      call init_gas(p%R_gas, p%sutherland, p%mu_const, p%mu_ref, p%T_ref, p%S_S, p%Pr)
      call partition_and_load(trim(p%mesh_file), ctx, mesh)
      call assign_patch_bcs(p, mesh, bc_dat)
      call alloc_state(s, mesh)
      call halo_init_persistent(mesh, ctx)
      call halo_init_persistent_gp(mesh, ctx)

      ! Initialise U from analytical W_a at every cell centroid.
      do c = 1, mesh%nc_total
         call mms_primitives(mesh%cell_centroid(1, c), &
                              mesh%cell_centroid(2, c), &
                              mesh%cell_centroid(3, c), W)
         call cons_from_prim(W(IP_RHO), W(IP_U), W(IP_V), W(IP_W), W(IP_P), U_a)
         s%U(:, c) = U_a
      end do

      ! One residual evaluation, with MUSCL on and viscous on.
      call halo_pack_and_start(mesh, s%U, ctx)
      call halo_wait(mesh, s%U)
      call compute_primitives(mesh, s)
      call compute_gradients(mesh, bc_dat, s)
      call compute_venkat_limiter(mesh, s, p%venkat_K)
      call halo_pack_and_start_gp(mesh, s, ctx)
      call residual_begin(s)
      call residual_pure_interior(mesh, s, p%muscl_enabled, p%viscous_enabled)
      call halo_wait_gp(mesh, s)
      call residual_partition(mesh, s, p%muscl_enabled, p%viscous_enabled)
      call residual_boundary(mesh, bc_dat, s, p%muscl_enabled, p%viscous_enabled)

      ! Compute L2 norm of ||R - V·S||_2 / Σ V on local cells.
      err_local = 0.0_wp
      vol_local = 0.0_wp
      do c = 1, mesh%nc_internal
         call mms_source(mesh%cell_centroid(1, c), &
                          mesh%cell_centroid(2, c), &
                          mesh%cell_centroid(3, c), Svec)
         err_local = err_local + sum((s%R(:, c) - mesh%cell_volume(c) * Svec)**2)
         vol_local = vol_local + mesh%cell_volume(c)
      end do

      if (ctx%nproc > 1) then
         call MPI_Allreduce(err_local, err_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ctx%comm)
         call MPI_Allreduce(vol_local, vol_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ctx%comm)
      else
         err_global = err_local
         vol_global = vol_local
      end if
      L2_global = sqrt(err_global / max(vol_global, 1.0e-30_wp))

      if (ctx%is_root) then
         write(*,'(A,A,A,I0,A,I0,A,1PE12.4)') &
            '[', trim(label), '] nc=', mesh%nc_internal, &
            ' nf=', mesh%nf, '  L2(R-VS)/sqrt(V) = ', L2_global
      end if

      call halo_free_persistent_gp(mesh)
      call halo_free_persistent(mesh)
      call free_state(s)
      if (allocated(bc_dat)) deallocate(bc_dat)

      ! Suppress unused warnings
      if (.false.) then
         rho = 0.0_wp; ux = 0.0_wp; uy = 0.0_wp; uz = 0.0_wp; pres = 0.0_wp
      end if
   end subroutine run_one

end program test_mms_order
