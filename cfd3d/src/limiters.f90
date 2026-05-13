! Venkatakrishnan slope limiter — smooth, differentiable, robust for
! compressible flow with shocks. ψ ∈ [0, 1] per cell per primitive variable.
!
! Reference: Venkatakrishnan, AIAA-93-0880, "On the accuracy of limiters and
! convergence to steady state solutions".
module limiters
   use kinds,      only : wp
   use constants,  only : NPRIM
   use mesh_types, only : t_mesh
   use fields,     only : t_state
   implicit none
   private

   public :: compute_venkat_limiter

contains

   subroutine compute_venkat_limiter(mesh, s, K)
      type(t_mesh),  intent(in)    :: mesh
      type(t_state), intent(inout) :: s
      real(wp),      intent(in)    :: K       ! Venkat tuning constant (typ. 1-10)

      integer  :: c, f, c_o, c_n, v
      real(wp), allocatable :: W_min(:,:), W_max(:,:)
      real(wp) :: dW_face, Wcell, eps2, h
      real(wp) :: delta_pos, delta_neg, psi_face

      allocate(W_min(NPRIM, mesh%nc_internal))
      allocate(W_max(NPRIM, mesh%nc_internal))
      ! Seed min/max with the cell's own state.
      do c = 1, mesh%nc_internal
         W_min(:, c) = s%W(:, c)
         W_max(:, c) = s%W(:, c)
      end do

      ! Pass 1: sweep faces to compute per-cell stencil min/max over all neighbors.
      do f = 1, mesh%nf_interior          ! pure-interior + partition
         c_o = mesh%face_owner(f)
         c_n = mesh%face_neighbor(f)
         if (c_o <= mesh%nc_internal) then
            do v = 1, NPRIM
               W_min(v, c_o) = min(W_min(v, c_o), s%W(v, c_n))
               W_max(v, c_o) = max(W_max(v, c_o), s%W(v, c_n))
            end do
         end if
         if (c_n <= mesh%nc_internal) then
            do v = 1, NPRIM
               W_min(v, c_n) = min(W_min(v, c_n), s%W(v, c_o))
               W_max(v, c_n) = max(W_max(v, c_n), s%W(v, c_o))
            end do
         end if
      end do
      ! Note: physical-boundary face neighbor states would normally also feed
      ! the stencil, but doing so robustly requires evaluating BC ghost states
      ! here. We use the cell-only contribution for boundary faces, which gives
      ! a slightly tighter stencil — safer for shocks.

      ! Pass 2: per-cell, find limiting ψ over its faces.
      ! Initialise ψ to 1 (no limiting).
      do c = 1, mesh%nc_internal
         s%psi(:, c) = 1.0_wp
      end do

      do f = 1, mesh%nf
         c_o = mesh%face_owner(f)
         c_n = mesh%face_neighbor(f)
         if (c_o <= mesh%nc_internal) then
            h = max((mesh%cell_volume(c_o))**(1.0_wp/3.0_wp), 1.0e-30_wp)
            eps2 = (K * h)**3
            do v = 1, NPRIM
               Wcell  = s%W(v, c_o)
               dW_face = s%gradW(1, v, c_o) * (mesh%face_centroid(1, f) - mesh%cell_centroid(1, c_o)) &
                       + s%gradW(2, v, c_o) * (mesh%face_centroid(2, f) - mesh%cell_centroid(2, c_o)) &
                       + s%gradW(3, v, c_o) * (mesh%face_centroid(3, f) - mesh%cell_centroid(3, c_o))
               if (dW_face > 0.0_wp) then
                  delta_pos = W_max(v, c_o) - Wcell
                  psi_face = venkat_psi(delta_pos, dW_face, eps2)
               else if (dW_face < 0.0_wp) then
                  delta_neg = W_min(v, c_o) - Wcell
                  psi_face = venkat_psi(delta_neg, dW_face, eps2)
               else
                  psi_face = 1.0_wp
               end if
               s%psi(v, c_o) = min(s%psi(v, c_o), psi_face)
            end do
         end if
         if (c_n > 0 .and. c_n <= mesh%nc_internal) then
            h = max((mesh%cell_volume(c_n))**(1.0_wp/3.0_wp), 1.0e-30_wp)
            eps2 = (K * h)**3
            do v = 1, NPRIM
               Wcell  = s%W(v, c_n)
               dW_face = s%gradW(1, v, c_n) * (mesh%face_centroid(1, f) - mesh%cell_centroid(1, c_n)) &
                       + s%gradW(2, v, c_n) * (mesh%face_centroid(2, f) - mesh%cell_centroid(2, c_n)) &
                       + s%gradW(3, v, c_n) * (mesh%face_centroid(3, f) - mesh%cell_centroid(3, c_n))
               if (dW_face > 0.0_wp) then
                  delta_pos = W_max(v, c_n) - Wcell
                  psi_face = venkat_psi(delta_pos, dW_face, eps2)
               else if (dW_face < 0.0_wp) then
                  delta_neg = W_min(v, c_n) - Wcell
                  psi_face = venkat_psi(delta_neg, dW_face, eps2)
               else
                  psi_face = 1.0_wp
               end if
               s%psi(v, c_n) = min(s%psi(v, c_n), psi_face)
            end do
         end if
      end do

      deallocate(W_min, W_max)
   end subroutine compute_venkat_limiter

   ! Venkat limiter formula (Venkatakrishnan 1993, eq. 9):
   !   ψ = (Δ_pos² + 2 Δ_neg Δ_pos + ε²) / (Δ_pos² + 2 Δ_neg² + Δ_pos Δ_neg + ε²)
   ! where Δ_pos = (max/min - cell_value) [signed], Δ_neg = ∇W·dr (signed).
   pure function venkat_psi(delta_p, delta_m, eps2) result(psi)
      real(wp), intent(in) :: delta_p, delta_m, eps2
      real(wp) :: psi, num, den
      num = (delta_p*delta_p + eps2) * delta_m + 2.0_wp * delta_m*delta_m * delta_p
      den =  delta_m * (delta_p*delta_p + 2.0_wp*delta_m*delta_m + delta_p*delta_m + eps2)
      if (abs(den) < 1.0e-30_wp) then
         psi = 1.0_wp
      else
         psi = num / den
      end if
      psi = max(0.0_wp, min(1.0_wp, psi))
   end function venkat_psi

end module limiters
