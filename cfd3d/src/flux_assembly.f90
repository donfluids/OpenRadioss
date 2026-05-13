module flux_assembly
   use kinds,          only : wp
   use constants,      only : NVAR, NPRIM, IP_RHO, IP_U, IP_V, IP_W, IP_P, &
                              BC_SLIP_WALL, BC_SYMMETRY, BC_SUPERSONIC_OUTLET, &
                              BC_SUPERSONIC_INLET, BC_FARFIELD, &
                              BC_NO_SLIP_WALL, BC_DIRICHLET
   use mesh_types,     only : t_mesh, t_patch
   use bc_types,       only : t_bc_data
   use bc_apply,       only : ghost_state
   use eos_ideal_gas,  only : cons_from_prim, prim_from_cons
   use riemann_hllc,   only : hllc_flux
   use viscous_fluxes, only : viscous_flux_face
   use sgs_model,      only : sgs_filter_width
   use fields,         only : t_state
   use mpi_runtime,    only : t_mpi_ctx
   use halo_exchange,  only : halo_pack_and_start, halo_wait
   implicit none
   private

   public :: compute_residual
   public :: residual_begin
   public :: residual_pure_interior, residual_partition, residual_boundary

contains

   ! Slip and symmetry walls are inviscid by construction — no shear can be
   ! transmitted by the wall. Returning .false. for these prevents the LSQ
   ! reflection trick from manufacturing a spurious viscous stress.
   pure function bc_has_viscous_flux(bc_type) result(yes)
      integer, intent(in) :: bc_type
      logical :: yes
      select case (bc_type)
      case (BC_SLIP_WALL, BC_SYMMETRY, BC_SUPERSONIC_OUTLET, BC_FARFIELD)
         yes = .false.
      case (BC_NO_SLIP_WALL, BC_DIRICHLET, BC_SUPERSONIC_INLET)
         yes = .true.
      case default
         yes = .true.
      end select
   end function bc_has_viscous_flux

   ! MPI-aware orchestrator. (Halos for U and grads must already be
   ! exchanged by the time-integration driver; this routine only handles
   ! the state-halo overlap with pure-interior compute.)
   subroutine compute_residual(mesh, bc_dat, s, ctx, muscl, viscous)
      type(t_mesh),    intent(inout) :: mesh
      type(t_bc_data), intent(in)    :: bc_dat(:)
      type(t_state),   intent(inout) :: s
      type(t_mpi_ctx), intent(in)    :: ctx
      logical,         intent(in)    :: muscl, viscous

      call residual_begin(s)
      call halo_pack_and_start(mesh, s%U, ctx)
      call residual_pure_interior(mesh, s, muscl, viscous)
      call halo_wait(mesh, s%U)
      call residual_partition(mesh, s, muscl, viscous)
      call residual_boundary(mesh, bc_dat, s, muscl, viscous)
   end subroutine compute_residual

   subroutine residual_begin(s)
      type(t_state), intent(inout) :: s
      s%R = 0.0_wp
   end subroutine residual_begin

   subroutine residual_pure_interior(mesh, s, muscl, viscous)
      type(t_mesh),  intent(in)    :: mesh
      type(t_state), intent(inout) :: s
      logical,       intent(in)    :: muscl, viscous
      integer  :: ifc
      do ifc = 1, mesh%nf_pure_interior
         call interior_face_residual(mesh, s, ifc, muscl, viscous)
      end do
   end subroutine residual_pure_interior

   subroutine residual_partition(mesh, s, muscl, viscous)
      type(t_mesh),  intent(in)    :: mesh
      type(t_state), intent(inout) :: s
      logical,       intent(in)    :: muscl, viscous
      integer  :: ifc
      do ifc = mesh%nf_pure_interior + 1, mesh%nf_interior
         call interior_face_residual(mesh, s, ifc, muscl, viscous)
      end do
   end subroutine residual_partition

   subroutine interior_face_residual(mesh, s, ifc, muscl, viscous)
      type(t_mesh),  intent(in)    :: mesh
      type(t_state), intent(inout) :: s
      integer,       intent(in)    :: ifc
      logical,       intent(in)    :: muscl, viscous
      integer  :: c_o, c_n
      real(wp) :: nrml(3), area, drO(3), drN(3)
      real(wp) :: W_L(NPRIM), W_R(NPRIM), Q_L(NVAR), Q_R(NVAR), Fflx(NVAR)
      real(wp) :: Wf(NPRIM), gWf(3, NPRIM), Fvis(NVAR)
      integer  :: v, j

      c_o  = mesh%face_owner(ifc)
      c_n  = mesh%face_neighbor(ifc)
      nrml = mesh%face_normal(:, ifc)
      area = mesh%face_area(ifc)

      drO = mesh%face_centroid(:, ifc) - mesh%cell_centroid(:, c_o)
      drN = mesh%face_centroid(:, ifc) - mesh%cell_centroid(:, c_n)

      if (muscl) then
         do v = 1, NPRIM
            W_L(v) = s%W(v, c_o) + s%psi(v, c_o) * &
                     (s%gradW(1, v, c_o) * drO(1) + s%gradW(2, v, c_o) * drO(2) + &
                      s%gradW(3, v, c_o) * drO(3))
            W_R(v) = s%W(v, c_n) + s%psi(v, c_n) * &
                     (s%gradW(1, v, c_n) * drN(1) + s%gradW(2, v, c_n) * drN(2) + &
                      s%gradW(3, v, c_n) * drN(3))
         end do
      else
         W_L = s%W(:, c_o)
         W_R = s%W(:, c_n)
      end if

      call cons_from_prim(W_L(IP_RHO), W_L(IP_U), W_L(IP_V), W_L(IP_W), W_L(IP_P), Q_L)
      call cons_from_prim(W_R(IP_RHO), W_R(IP_U), W_R(IP_V), W_R(IP_W), W_R(IP_P), Q_R)
      Fflx = hllc_flux(Q_L, Q_R, nrml)

      if (viscous) then
         Wf = 0.5_wp * (W_L + W_R)
         do v = 1, NPRIM
            do j = 1, 3
               gWf(j, v) = 0.5_wp * (s%gradW(j, v, c_o) + s%gradW(j, v, c_n))
            end do
         end do
         block
            real(wp) :: Delta_face
            Delta_face = 0.5_wp * ( sgs_filter_width(mesh%cell_volume(c_o)) &
                                  + sgs_filter_width(mesh%cell_volume(c_n)) )
            call viscous_flux_face(Wf, gWf, nrml, Delta_face, Fvis)
         end block
         Fflx = Fflx - Fvis
      end if

      s%R(:, c_o) = s%R(:, c_o) + Fflx * area
      s%R(:, c_n) = s%R(:, c_n) - Fflx * area
   end subroutine interior_face_residual

   subroutine residual_boundary(mesh, bc_dat, s, muscl, viscous)
      type(t_mesh),    intent(in)    :: mesh
      type(t_bc_data), intent(in)    :: bc_dat(:)
      type(t_state),   intent(inout) :: s
      logical,         intent(in)    :: muscl, viscous

      integer  :: ifc, c_o, ip, v
      real(wp) :: nrml(3), area, drO(3)
      real(wp) :: W_L(NPRIM), Q_L(NVAR), Q_R(NVAR), Fflx(NVAR)
      real(wp) :: W_R(NPRIM), Wf(NPRIM), gWf(3, NPRIM), Fvis(NVAR)
      real(wp) :: rho_b, u_b, v_b, w_b, p_b

      do ifc = mesh%nf_interior + 1, mesh%nf
         c_o  = mesh%face_owner(ifc)
         ip   = mesh%face_patch(ifc)
         nrml = mesh%face_normal(:, ifc)
         area = mesh%face_area(ifc)

         drO = mesh%face_centroid(:, ifc) - mesh%cell_centroid(:, c_o)

         if (muscl) then
            do v = 1, NPRIM
               W_L(v) = s%W(v, c_o) + s%psi(v, c_o) * &
                        (s%gradW(1, v, c_o) * drO(1) + s%gradW(2, v, c_o) * drO(2) + &
                         s%gradW(3, v, c_o) * drO(3))
            end do
         else
            W_L = s%W(:, c_o)
         end if
         call cons_from_prim(W_L(IP_RHO), W_L(IP_U), W_L(IP_V), W_L(IP_W), W_L(IP_P), Q_L)
         call ghost_state(mesh%patches(ip)%bc_type, bc_dat(ip), Q_L, nrml, Q_R)
         Fflx = hllc_flux(Q_L, Q_R, nrml)

         if (viscous .and. bc_has_viscous_flux(mesh%patches(ip)%bc_type)) then
            call prim_from_cons(Q_R, rho_b, u_b, v_b, w_b, p_b)
            W_R(IP_RHO) = rho_b; W_R(IP_U) = u_b; W_R(IP_V) = v_b
            W_R(IP_W)   = w_b;   W_R(IP_P) = p_b
            Wf = 0.5_wp * (W_L + W_R)
            ! Over-relaxed gradient correction at the boundary face: replace
            ! the along-dr component of the cell gradient with the local
            ! one-sided difference (W_face - W_cell)/|dr|. This is critical
            ! at no-slip walls where the cell-averaged LSQ gradient would
            ! over-state the wall stress in transient.
            block
               real(wp) :: ndr(3), dr_mag, local_deriv, dot_gn
               integer  :: v
               ndr = drO   ! face_centroid - owner_centroid
               dr_mag = sqrt(ndr(1)**2 + ndr(2)**2 + ndr(3)**2)
               if (dr_mag > 0.0_wp) ndr = ndr / dr_mag
               do v = 1, NPRIM
                  gWf(:, v) = s%gradW(:, v, c_o)
                  if (dr_mag > 0.0_wp) then
                     local_deriv = (Wf(v) - s%W(v, c_o)) / dr_mag
                     dot_gn = gWf(1, v)*ndr(1) + gWf(2, v)*ndr(2) + gWf(3, v)*ndr(3)
                     gWf(:, v) = gWf(:, v) + (local_deriv - dot_gn) * ndr
                  end if
               end do
            end block
            call viscous_flux_face(Wf, gWf, nrml, sgs_filter_width(mesh%cell_volume(c_o)), Fvis)
            Fflx = Fflx - Fvis
         end if

         s%R(:, c_o) = s%R(:, c_o) + Fflx * area
      end do
   end subroutine residual_boundary

end module flux_assembly
