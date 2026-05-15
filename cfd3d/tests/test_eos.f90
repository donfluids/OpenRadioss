program test_eos
   use, intrinsic :: iso_fortran_env, only : real64, error_unit
   use kinds,         only : wp
   use constants,     only : NVAR
   use eos_ideal_gas, only : prim_from_cons, cons_from_prim, sound_speed, max_wave_speed
   implicit none

   integer :: i, fails
   real(wp) :: rho, u, v, w, p
   real(wp) :: rho2, u2, v2, w2, p2
   real(wp) :: Q(NVAR), n(3), c, s
   real(wp), parameter :: TOL = 1.0e-12_wp

   fails = 0

   ! Round-trip: prim -> cons -> prim
   do i = 1, 20
      rho = 0.1_wp + 0.5_wp * real(i, wp)
      u   = -2.0_wp + 0.3_wp * real(i, wp)
      v   = 0.5_wp * real(i, wp) - 4.0_wp
      w   = real(mod(i, 7), wp) * 0.25_wp
      p   = 0.2_wp + 0.7_wp * real(i, wp)
      call cons_from_prim(rho, u, v, w, p, Q)
      call prim_from_cons(Q, rho2, u2, v2, w2, p2)
      if (abs(rho-rho2) > TOL * abs(rho) .or. &
          abs(u-u2)     > TOL * (1.0_wp+abs(u)) .or. &
          abs(v-v2)     > TOL * (1.0_wp+abs(v)) .or. &
          abs(w-w2)     > TOL * (1.0_wp+abs(w)) .or. &
          abs(p-p2)     > TOL * abs(p) ) then
         write(error_unit,'(A,I0,5(1X,1PE12.4))') 'roundtrip fail i=', i, &
            rho-rho2, u-u2, v-v2, w-w2, p-p2
         fails = fails + 1
      end if
   end do

   ! Sound speed positivity & monotone in p
   c = sound_speed(1.0_wp, 1.0_wp)
   if (c <= 0.0_wp) then
      write(error_unit,'(A)') 'sound_speed not positive'
      fails = fails + 1
   end if
   if (sound_speed(1.0_wp, 2.0_wp) <= sound_speed(1.0_wp, 1.0_wp)) then
      write(error_unit,'(A)') 'sound_speed not monotone in p'
      fails = fails + 1
   end if

   ! max_wave_speed >= 0, equals c when u=0
   call cons_from_prim(1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, Q)
   n = [1.0_wp, 0.0_wp, 0.0_wp]
   s = max_wave_speed(Q, n)
   if (abs(s - sound_speed(1.0_wp, 1.0_wp)) > TOL) then
      write(error_unit,'(A,2(1X,1PE12.4))') 'max_wave_speed wrong: ', s, sound_speed(1.0_wp, 1.0_wp)
      fails = fails + 1
   end if

   if (fails == 0) then
      write(*,'(A)') 'test_eos: PASS'
   else
      write(*,'(A,I0,A)') 'test_eos: FAIL (', fails, ' issues)'
      stop 1
   end if
end program test_eos
