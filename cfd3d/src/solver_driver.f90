module solver_driver
   use kinds,            only : wp
   use constants,        only : NVAR
   use mesh_types,       only : t_mesh, t_patch
   use partition,        only : partition_and_load
   use fields,           only : t_state, alloc_state, free_state
   use bc_types,         only : t_bc_data
   use eos_ideal_gas,    only : cons_from_prim
   use time_integration, only : compute_dt, rk3_step
   use io_vtk_legacy,    only : write_vtk
   use solver_control,   only : t_run_params, bc_string_to_int
   use mpi_runtime,      only : t_mpi_ctx
   use halo_exchange,    only : halo_init_persistent, halo_free_persistent, &
                                halo_init_persistent_gp, halo_free_persistent_gp
   use gas_properties,   only : init_gas
   implicit none
   private

   public :: run_case, assign_patch_bcs, set_initial_condition

contains

   subroutine run_case(p, ctx)
      type(t_run_params), intent(in) :: p
      type(t_mpi_ctx),    intent(in) :: ctx

      type(t_mesh)    :: mesh
      type(t_state)   :: s
      type(t_bc_data), allocatable :: bc_dat(:)

      real(wp) :: t, dt, next_out
      integer  :: step, out_idx
      character(len=512) :: outfile

      call init_gas(p%R_gas, p%sutherland, p%mu_const, p%mu_ref, p%T_ref, p%S_S, p%Pr)

      call partition_and_load(trim(p%mesh_file), ctx, mesh)
      call assign_patch_bcs(p, mesh, bc_dat)
      call alloc_state(s, mesh)
      call set_initial_condition(p, mesh, s)
      call halo_init_persistent(mesh, ctx)
      call halo_init_persistent_gp(mesh, ctx)

      write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
         'rank=', ctx%rank, &
         '  nc_loc=', mesh%nc_internal, &
         '  nc_ghost=', mesh%nc_total - mesh%nc_internal, &
         '  nf_pure=', mesh%nf_pure_interior, &
         '  nf_part=', mesh%nf_interior - mesh%nf_pure_interior, &
         '  nf_bnd=', mesh%nf_boundary, &
         '  np_halo_neighbors=', mesh%halo%n_neighbors, &
         '  np=', mesh%np

      out_idx = 0
      call write_one(p, mesh, s, 0.0_wp, ctx, out_idx, outfile)

      t = 0.0_wp
      next_out = p%output_interval
      step = 0
      do while (t < p%t_end .and. step < p%max_steps)
         dt = compute_dt(mesh, s, p%cfl, ctx)
         if (t + dt > p%t_end) dt = p%t_end - t
         call rk3_step(mesh, bc_dat, s, dt, ctx, p%muscl_enabled, p%viscous_enabled, p%venkat_K)
         t = t + dt
         step = step + 1
         if (ctx%is_root .and. mod(step, 50) == 0) then
            write(*,'(A,I7,A,1PE12.5,A,1PE12.5)') &
               'step ', step, '  t=', t, '  dt=', dt
         end if
         if (t >= next_out - 1.0e-15_wp) then
            call write_one(p, mesh, s, t, ctx, out_idx, outfile)
            next_out = next_out + p%output_interval
         end if
      end do

      ! Final dump
      call write_one(p, mesh, s, t, ctx, out_idx, outfile)

      call halo_free_persistent_gp(mesh)
      call halo_free_persistent(mesh)
      call free_state(s)
      if (allocated(bc_dat)) deallocate(bc_dat)
   end subroutine run_case

   subroutine assign_patch_bcs(p, mesh, bc_dat)
      type(t_run_params),           intent(in)    :: p
      type(t_mesh),                 intent(inout) :: mesh
      type(t_bc_data), allocatable, intent(out)   :: bc_dat(:)

      integer :: ip, ii
      allocate(bc_dat(mesh%np))
      do ip = 1, mesh%np
         bc_dat(ip)%rho = 1.0_wp
         bc_dat(ip)%u   = 0.0_wp
         bc_dat(ip)%v   = 0.0_wp
         bc_dat(ip)%w   = 0.0_wp
         bc_dat(ip)%p   = 1.0_wp
         mesh%patches(ip)%bc_type = 0
      end do

      do ip = 1, mesh%np
         do ii = 1, p%patch_count
            if (trim(mesh%patches(ip)%name) == trim(p%patch_name(ii))) then
               mesh%patches(ip)%bc_type = bc_string_to_int(p%patch_bc(ii))
               bc_dat(ip)%rho = p%patch_rho(ii)
               bc_dat(ip)%u   = p%patch_u(ii)
               bc_dat(ip)%v   = p%patch_v(ii)
               bc_dat(ip)%w   = p%patch_w(ii)
               bc_dat(ip)%p   = p%patch_p(ii)
               exit
            end if
         end do
         if (mesh%patches(ip)%bc_type == 0) then
            mesh%patches(ip)%bc_type = bc_string_to_int('slip_wall')
         end if
      end do
   end subroutine assign_patch_bcs

   subroutine set_initial_condition(p, mesh, s)
      type(t_run_params), intent(in)    :: p
      type(t_mesh),       intent(in)    :: mesh
      type(t_state),      intent(inout) :: s

      integer :: c, ax
      real(wp) :: pos, U(NVAR)

      select case (trim(p%init_type))
      case ('uniform')
         call cons_from_prim(p%rho_L, p%u_L, p%v_L, p%w_L, p%p_L, U)
         do c = 1, mesh%nc_total
            s%U(:, c) = U
         end do
      case ('riemann')
         ax = p%diaphragm_axis
         if (ax < 1 .or. ax > 3) then
            write(*,'(A)') 'set_initial_condition: invalid diaphragm_axis'
            error stop 1
         end if
         do c = 1, mesh%nc_total
            pos = mesh%cell_centroid(ax, c)
            if (pos < p%diaphragm_pos) then
               call cons_from_prim(p%rho_L, p%u_L, p%v_L, p%w_L, p%p_L, s%U(:, c))
            else
               call cons_from_prim(p%rho_R, p%u_R, p%v_R, p%w_R, p%p_R, s%U(:, c))
            end if
         end do
      case default
         write(*,'(A,A)') 'set_initial_condition: unknown init_type ', trim(p%init_type)
         error stop 1
      end select
   end subroutine set_initial_condition

   subroutine write_one(p, mesh, s, t, ctx, idx, outfile)
      type(t_run_params),  intent(in)    :: p
      type(t_mesh),        intent(in)    :: mesh
      type(t_state),       intent(in)    :: s
      real(wp),            intent(in)    :: t
      type(t_mpi_ctx),     intent(in)    :: ctx
      integer,             intent(inout) :: idx
      character(len=512),  intent(out)   :: outfile
      character(len=8)  :: idxs
      character(len=8)  :: ranks
      write(idxs,'(I8.8)') idx
      write(ranks,'(I0)') ctx%rank
      outfile = trim(p%output_dir) // '/' // trim(p%case_name) // &
                '_' // trim(idxs) // '_r' // trim(ranks) // '.vtk'
      call write_vtk(trim(outfile), mesh, s, t)
      idx = idx + 1
   end subroutine write_one

end module solver_driver
