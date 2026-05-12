module kinds
   use, intrinsic :: iso_fortran_env, only : real64, int32, int64
   implicit none
   private

   public :: wp, i4, i8

   integer, parameter :: wp = real64
   integer, parameter :: i4 = int32
   integer, parameter :: i8 = int64
end module kinds
