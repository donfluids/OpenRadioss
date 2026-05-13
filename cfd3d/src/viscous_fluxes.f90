! Viscous flux contribution for compressible Navier-Stokes.
!
!   τ_ij = μ ( ∂u_i/∂x_j + ∂u_j/∂x_i ) - (2/3) μ (∇·u) δ_ij
!   q_i  = -k ∂T/∂x_i,  k = μ cp / Pr
!
! Per-face fluxes through the unit normal n:
!   F_mass = 0
!   F_mom_i = τ_ij n_j
!   F_E    = u_i τ_ij n_j - q_j n_j
!
! Face-centered gradient uses the arithmetic mean of the two cells' gradients
! ("face average"). For first M3 iteration we skip the directional correction
! along (x_n - x_o); good enough for the order test on smooth meshes and the
! Couette case.
module viscous_fluxes
   use kinds,          only : wp
   use constants,      only : NVAR, NPRIM, IP_RHO, IP_U, IP_V, IP_W, IP_P, &
                              IRHOU, IRHOV, IRHOW, IRHOE
   use gas_properties, only : gas
   use eos_ideal_gas,  only : temperature, mu_sutherland
   implicit none
   private

   public :: viscous_flux_face

contains

   ! Viscous flux through a face given the two cell-side primitive states and
   ! their primitive gradients (face-averaged). dT_face is the face-centered
   ! temperature gradient (derived from p and rho gradients via the ideal gas).
   pure subroutine viscous_flux_face(W_face, gW_face, n, F)
      real(wp), intent(in)  :: W_face(NPRIM)
      real(wp), intent(in)  :: gW_face(3, NPRIM)         ! (3, NPRIM)
      real(wp), intent(in)  :: n(3)
      real(wp), intent(out) :: F(NVAR)

      real(wp) :: rho, u, v, w, p, T, mu, k_th
      real(wp) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
      real(wp) :: drho_dx, drho_dy, drho_dz, dp_dx, dp_dy, dp_dz
      real(wp) :: dT_dx, dT_dy, dT_dz
      real(wp) :: div_u, txx, tyy, tzz, txy, txz, tyz
      real(wp) :: qx, qy, qz
      real(wp) :: tau_n_x, tau_n_y, tau_n_z
      real(wp) :: inv_rho_R

      rho = W_face(IP_RHO)
      u   = W_face(IP_U)
      v   = W_face(IP_V)
      w   = W_face(IP_W)
      p   = W_face(IP_P)
      T   = temperature(rho, p)
      mu  = mu_sutherland(T)
      k_th = mu * gas%cp / gas%Pr

      dudx = gW_face(1, IP_U); dudy = gW_face(2, IP_U); dudz = gW_face(3, IP_U)
      dvdx = gW_face(1, IP_V); dvdy = gW_face(2, IP_V); dvdz = gW_face(3, IP_V)
      dwdx = gW_face(1, IP_W); dwdy = gW_face(2, IP_W); dwdz = gW_face(3, IP_W)
      drho_dx = gW_face(1, IP_RHO); drho_dy = gW_face(2, IP_RHO); drho_dz = gW_face(3, IP_RHO)
      dp_dx   = gW_face(1, IP_P);   dp_dy   = gW_face(2, IP_P);   dp_dz   = gW_face(3, IP_P)

      ! ∇T = ∂(p/(ρR))/∂x = (dp - p dρ/ρ) / (ρ R)
      inv_rho_R = 1.0_wp / (max(rho, 1.0e-30_wp) * gas%R_gas)
      dT_dx = (dp_dx - p * drho_dx / max(rho, 1.0e-30_wp)) * inv_rho_R
      dT_dy = (dp_dy - p * drho_dy / max(rho, 1.0e-30_wp)) * inv_rho_R
      dT_dz = (dp_dz - p * drho_dz / max(rho, 1.0e-30_wp)) * inv_rho_R

      div_u = dudx + dvdy + dwdz

      txx = mu * (2.0_wp * dudx - (2.0_wp/3.0_wp) * div_u)
      tyy = mu * (2.0_wp * dvdy - (2.0_wp/3.0_wp) * div_u)
      tzz = mu * (2.0_wp * dwdz - (2.0_wp/3.0_wp) * div_u)
      txy = mu * (dudy + dvdx)
      txz = mu * (dudz + dwdx)
      tyz = mu * (dvdz + dwdy)

      qx = -k_th * dT_dx
      qy = -k_th * dT_dy
      qz = -k_th * dT_dz

      tau_n_x = txx * n(1) + txy * n(2) + txz * n(3)
      tau_n_y = txy * n(1) + tyy * n(2) + tyz * n(3)
      tau_n_z = txz * n(1) + tyz * n(2) + tzz * n(3)

      F(1)     = 0.0_wp
      F(IRHOU) = tau_n_x
      F(IRHOV) = tau_n_y
      F(IRHOW) = tau_n_z
      F(IRHOE) = u * tau_n_x + v * tau_n_y + w * tau_n_z &
               - (qx * n(1) + qy * n(2) + qz * n(3))
   end subroutine viscous_flux_face

end module viscous_fluxes
