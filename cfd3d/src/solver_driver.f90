module solver_driver
   use kinds,            only : wp
   use constants,        only : NVAR
   use mesh_types,       only : t_mesh, t_patch
   use partition,        only : partition_and_load
   use cut_cell,         only : build_cut_cell_tables, merge_sync_state
   use fields,           only : t_state, alloc_state, free_state
   use bc_types,         only : t_bc_data
   use eos_ideal_gas,    only : cons_from_prim
   use time_integration, only : compute_dt, rk3_step
   use io_vtk_legacy,    only : write_vtk
   use solver_control,   only : t_run_params, bc_string_to_int, sgs_string_to_int
   use mpi_runtime,      only : t_mpi_ctx
   use halo_exchange,    only : halo_init_persistent, halo_free_persistent, &
                                halo_init_persistent_gp, halo_free_persistent_gp
   use gas_properties,   only : init_gas
   use probes,           only : t_probe_set, probes_load, probes_locate, &
                                probes_sample, probes_flush, probes_free
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
      type(t_probe_set) :: probes

      real(wp) :: t, dt, next_out
      integer  :: step, out_idx
      character(len=512) :: outfile

      call init_gas(p%R_gas, p%sutherland, p%mu_const, p%mu_ref, p%T_ref, p%S_S, p%Pr, &
                    sgs_kind=sgs_string_to_int(p%sgs_model), &
                    C_s=p%C_s, C_w=p%C_w, Pr_t=p%Pr_t)

      call partition_and_load(trim(p%mesh_file), ctx, mesh)
      call build_cut_cell_tables(mesh, p%n_obstacles, &
           p%obstacle_cube_lo(:, 1:max(p%n_obstacles,1)), &
           p%obstacle_cube_hi(:, 1:max(p%n_obstacles,1)))
      call assign_patch_bcs(p, mesh, bc_dat)
      call alloc_state(s, mesh)
      call set_initial_condition(p, mesh, s, ctx)
      call merge_sync_state(mesh, s)
      call halo_init_persistent(mesh, ctx)
      call halo_init_persistent_gp(mesh, ctx)
      call probes_load(trim(p%probes_file), trim(p%case_name), trim(p%output_dir), ctx, probes)
      call probes_locate(mesh, ctx, probes)

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
      call probes_sample(mesh, s, t, ctx, probes)
      do while (t < p%t_end .and. step < p%max_steps)
         dt = compute_dt(mesh, s, p%cfl, ctx, p%viscous_enabled)
         if (t + dt > p%t_end) dt = p%t_end - t
         call rk3_step(mesh, bc_dat, s, dt, ctx, p%muscl_enabled, p%viscous_enabled, p%venkat_K)
         t = t + dt
         step = step + 1
         call probes_sample(mesh, s, t, ctx, probes)
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
      call probes_flush(probes, ctx)

      call probes_free(probes)
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

   subroutine set_initial_condition(p, mesh, s, ctx)
      use mpi_f08
      use constants,     only : GM1
      type(t_run_params), intent(in)    :: p
      type(t_mesh),       intent(in)    :: mesh
      type(t_state),      intent(inout) :: s
      type(t_mpi_ctx),    intent(in)    :: ctx

      integer :: c, ax
      real(wp) :: pos, U(NVAR), dx, dy, dz, r2, rblast2

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
      case ('blast_bubble')
         ! Hot region: rho=blast_rho, p=blast_p inside sphere; ambient outside.
         rblast2 = p%blast_radius * p%blast_radius
         do c = 1, mesh%nc_total
            dx = mesh%cell_centroid(1, c) - p%blast_center(1)
            dy = mesh%cell_centroid(2, c) - p%blast_center(2)
            dz = mesh%cell_centroid(3, c) - p%blast_center(3)
            r2 = dx*dx + dy*dy + dz*dz
            if (r2 <= rblast2) then
               call cons_from_prim(p%blast_rho, 0.0_wp, 0.0_wp, 0.0_wp, p%blast_p, s%U(:, c))
            else
               call cons_from_prim(p%ambient_rho, 0.0_wp, 0.0_wp, 0.0_wp, p%ambient_p, s%U(:, c))
            end if
         end do
      case ('sedov')
         ! Sedov-Taylor: deposit total energy E into the sphere of radius
         ! blast_radius around blast_center. p_hot = (γ-1)·E / V_region with
         ! V_region globally summed over MPI ranks for consistency.
         call init_sedov_state(p, mesh, s, ctx)
      case ('tnt_charge')
         ! TNT free-field charge. Map physical charge mass to deposited
         ! energy E = m_TNT · 4.184 MJ/kg, then use the Sedov machinery
         ! with E auto-computed. The deposit uses ambient density at
         ! blast_radius (the "fluid-energy equivalent" / balloon-analogy
         ! approach; full JWL evolution lives in a follow-up sub-phase).
         block
            type(t_run_params) :: pl
            pl = p
            pl%blast_energy = p%tnt_mass * p%tnt_specific_E
            call init_sedov_state(pl, mesh, s, ctx)
         end block
      case default
         write(*,'(A,A)') 'set_initial_condition: unknown init_type ', trim(p%init_type)
         error stop 1
      end select
   end subroutine set_initial_condition

   subroutine init_sedov_state(p, mesh, s, ctx)
      use mpi_f08
      use constants, only : GM1
      type(t_run_params), intent(in)    :: p
      type(t_mesh),       intent(in)    :: mesh
      type(t_state),      intent(inout) :: s
      type(t_mpi_ctx),    intent(in)    :: ctx

      integer  :: c
      real(wp) :: dx, dy, dz, r2, rblast2, V_local, V_global, p_hot

      rblast2 = p%blast_radius * p%blast_radius

      ! First pass: count hot-region volume on local cells.
      V_local = 0.0_wp
      do c = 1, mesh%nc_internal
         dx = mesh%cell_centroid(1, c) - p%blast_center(1)
         dy = mesh%cell_centroid(2, c) - p%blast_center(2)
         dz = mesh%cell_centroid(3, c) - p%blast_center(3)
         r2 = dx*dx + dy*dy + dz*dz
         if (r2 <= rblast2) V_local = V_local + mesh%cell_volume(c)
      end do
      if (ctx%nproc > 1) then
         call MPI_Allreduce(V_local, V_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ctx%comm)
      else
         V_global = V_local
      end if
      if (V_global <= 0.0_wp) then
         if (ctx%is_root) write(*,'(A)') &
            'init_sedov_state: no cells inside blast_radius; refine the mesh or increase blast_radius'
         error stop 1
      end if
      p_hot = GM1 * p%blast_energy / V_global

      ! Second pass: write state. ghosts get the same treatment so partition
      ! faces see consistent ICs before the first halo exchange.
      do c = 1, mesh%nc_total
         dx = mesh%cell_centroid(1, c) - p%blast_center(1)
         dy = mesh%cell_centroid(2, c) - p%blast_center(2)
         dz = mesh%cell_centroid(3, c) - p%blast_center(3)
         r2 = dx*dx + dy*dy + dz*dz
         if (r2 <= rblast2) then
            call cons_from_prim(p%ambient_rho, 0.0_wp, 0.0_wp, 0.0_wp, p_hot, s%U(:, c))
         else
            call cons_from_prim(p%ambient_rho, 0.0_wp, 0.0_wp, 0.0_wp, p%ambient_p, s%U(:, c))
         end if
      end do
   end subroutine init_sedov_state

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
