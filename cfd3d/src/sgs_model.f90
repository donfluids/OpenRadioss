! M4 algebraic sub-grid scale (SGS) eddy-viscosity models.
!
!   Smagorinsky:   ν_t = (C_s Δ)² |S|,           |S| = √(2 S_ij S_ij)
!   WALE:          ν_t = (C_w Δ)² (S^d:S^d)^{3/2} / ((S:S)^{5/2} + (S^d:S^d)^{5/4})
!
! where S_ij = ½(∂u_i/∂x_j + ∂u_j/∂x_i) is the strain-rate tensor and S^d the
! traceless symmetric part of (∇u)·(∇u). WALE self-damps near walls (∝ y³)
! without an explicit Van Driest factor — Nicoud & Ducros 1999.
!
! Inputs are the 9 velocity-gradient components packed in `gW(3, NPRIM)` rows
! IP_U/IP_V/IP_W (consistent with `t_state%gradW` / face-averaged gradient
! passed by `viscous_flux_face`). Filter width Δ is per-face (V_cell^(1/3)).
module sgs_model
   use kinds,          only : wp
   use constants,      only : NPRIM, IP_U, IP_V, IP_W
   use gas_properties, only : gas, SGS_NONE, SGS_SMAG, SGS_WALE
   implicit none
   private

   public :: sgs_nu_t, sgs_filter_width

contains

   pure function sgs_filter_width(V_cell) result(d)
      real(wp), intent(in) :: V_cell
      real(wp) :: d
      d = max(V_cell, 1.0e-30_wp) ** (1.0_wp / 3.0_wp)
   end function sgs_filter_width

   ! Eddy viscosity from velocity-gradient tensor. Branches on gas%sgs_kind.
   pure function sgs_nu_t(gW, Delta) result(nu_t)
      real(wp), intent(in) :: gW(3, NPRIM)
      real(wp), intent(in) :: Delta
      real(wp) :: nu_t
      real(wp) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
      real(wp) :: S11, S22, S33, S12, S13, S23
      real(wp) :: SijSij, SmagS

      if (gas%sgs_kind == SGS_NONE) then
         nu_t = 0.0_wp
         return
      end if

      dudx = gW(1, IP_U); dudy = gW(2, IP_U); dudz = gW(3, IP_U)
      dvdx = gW(1, IP_V); dvdy = gW(2, IP_V); dvdz = gW(3, IP_V)
      dwdx = gW(1, IP_W); dwdy = gW(2, IP_W); dwdz = gW(3, IP_W)

      S11 = dudx
      S22 = dvdy
      S33 = dwdz
      S12 = 0.5_wp * (dudy + dvdx)
      S13 = 0.5_wp * (dudz + dwdx)
      S23 = 0.5_wp * (dvdz + dwdy)
      SijSij = S11*S11 + S22*S22 + S33*S33 + 2.0_wp * (S12*S12 + S13*S13 + S23*S23)
      SmagS  = sqrt(2.0_wp * SijSij)

      select case (gas%sgs_kind)
      case (SGS_SMAG)
         nu_t = (gas%C_s * Delta)**2 * SmagS
      case (SGS_WALE)
         nu_t = wale_nu_t(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, &
                          SijSij, Delta)
      case default
         nu_t = 0.0_wp
      end select
   end function sgs_nu_t

   pure function wale_nu_t(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, &
                            SijSij, Delta) result(nu_t)
      real(wp), intent(in) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
      real(wp), intent(in) :: SijSij, Delta
      real(wp) :: nu_t
      real(wp) :: g2_11, g2_22, g2_33, g2_12, g2_13, g2_23, g2_21, g2_31, g2_32
      real(wp) :: Sd11, Sd22, Sd33, Sd12, Sd13, Sd23, trace_g2
      real(wp) :: SdSd, num, den

      ! g²_ij = g_ik g_kj where g_ij = ∂u_i/∂x_j
      g2_11 = dudx*dudx + dudy*dvdx + dudz*dwdx
      g2_22 = dvdx*dudy + dvdy*dvdy + dvdz*dwdy
      g2_33 = dwdx*dudz + dwdy*dvdz + dwdz*dwdz
      g2_12 = dudx*dudy + dudy*dvdy + dudz*dwdy
      g2_21 = dvdx*dudx + dvdy*dvdx + dvdz*dwdx
      g2_13 = dudx*dudz + dudy*dvdz + dudz*dwdz
      g2_31 = dwdx*dudx + dwdy*dvdx + dwdz*dwdx
      g2_23 = dvdx*dudz + dvdy*dvdz + dvdz*dwdz
      g2_32 = dwdx*dudy + dwdy*dvdy + dwdz*dwdy

      trace_g2 = g2_11 + g2_22 + g2_33
      ! Symmetric part: Sd_ij = 0.5(g²_ij + g²_ji) - (1/3) δ_ij trace(g²)
      Sd11 = g2_11 - trace_g2 / 3.0_wp
      Sd22 = g2_22 - trace_g2 / 3.0_wp
      Sd33 = g2_33 - trace_g2 / 3.0_wp
      Sd12 = 0.5_wp * (g2_12 + g2_21)
      Sd13 = 0.5_wp * (g2_13 + g2_31)
      Sd23 = 0.5_wp * (g2_23 + g2_32)

      SdSd = Sd11*Sd11 + Sd22*Sd22 + Sd33*Sd33 &
           + 2.0_wp * (Sd12*Sd12 + Sd13*Sd13 + Sd23*Sd23)

      num = SdSd ** 1.5_wp
      den = (SijSij ** 2.5_wp) + (SdSd ** 1.25_wp)
      if (den <= 1.0e-30_wp) then
         nu_t = 0.0_wp
      else
         nu_t = (gas%C_w * Delta)**2 * num / den
      end if
   end function wale_nu_t

end module sgs_model
