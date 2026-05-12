module bc_apply
   use kinds,         only : wp
   use constants,     only : NVAR, IRHO, IRHOU, IRHOV, IRHOW, IRHOE, &
                             BC_SLIP_WALL, BC_SYMMETRY, &
                             BC_SUPERSONIC_INLET, BC_SUPERSONIC_OUTLET, BC_FARFIELD
   use eos_ideal_gas, only : prim_from_cons, cons_from_prim, sound_speed
   use bc_types,      only : t_bc_data
   implicit none
   private

   public :: ghost_state

contains

   ! Construct the ghost (right) state for a boundary face given the interior
   ! (left) state and the unit outward normal. The face normal stored in
   ! t_mesh points from the owner cell outward for boundary faces.
   pure subroutine ghost_state(bc_type, bc_dat, UL, n, UR)
      integer,         intent(in)  :: bc_type
      type(t_bc_data), intent(in)  :: bc_dat
      real(wp),        intent(in)  :: UL(NVAR), n(3)
      real(wp),        intent(out) :: UR(NVAR)

      real(wp) :: rho, u, v, w, p, vn
      real(wp) :: ur_rho, ur_u, ur_v, ur_w, ur_p, c

      select case (bc_type)
      case (BC_SLIP_WALL, BC_SYMMETRY)
         call prim_from_cons(UL, rho, u, v, w, p)
         vn = u*n(1) + v*n(2) + w*n(3)
         ur_rho = rho
         ur_u   = u - 2.0_wp * vn * n(1)
         ur_v   = v - 2.0_wp * vn * n(2)
         ur_w   = w - 2.0_wp * vn * n(3)
         ur_p   = p
         call cons_from_prim(ur_rho, ur_u, ur_v, ur_w, ur_p, UR)

      case (BC_SUPERSONIC_INLET)
         call cons_from_prim(bc_dat%rho, bc_dat%u, bc_dat%v, bc_dat%w, bc_dat%p, UR)

      case (BC_SUPERSONIC_OUTLET)
         UR = UL   ! zero-gradient extrapolation

      case (BC_FARFIELD)
         ! Classify via normal Mach number; supersonic-in -> Dirichlet,
         ! supersonic-out -> extrapolation, subsonic -> simple zero-gradient
         ! (refined in M3 with characteristic BCs).
         call prim_from_cons(UL, rho, u, v, w, p)
         c = sound_speed(rho, p)
         vn = u*n(1) + v*n(2) + w*n(3)
         if (vn <= -c) then
            ! Supersonic inflow
            call cons_from_prim(bc_dat%rho, bc_dat%u, bc_dat%v, bc_dat%w, bc_dat%p, UR)
         else if (vn >= c) then
            ! Supersonic outflow
            UR = UL
         else
            ! Subsonic — coarse zero-gradient
            UR = UL
         end if

      case default
         UR = UL
      end select
   end subroutine ghost_state

end module bc_apply
