! ISO C binding interface to METIS_PartMeshDual.
! METIS idx_t defaults to 32-bit int and real_t defaults to float.
! See /usr/include/metis.h on Ubuntu (libmetis-dev 5.x).
module metis_binding
   use, intrinsic :: iso_c_binding
   implicit none
   private

   public :: c_idx, c_real, METIS_PartMeshDual_C, METIS_OK

   ! idx_t is 32-bit int by default in the Ubuntu libmetis package.
   integer, parameter :: c_idx  = c_int32_t
   integer, parameter :: c_real = c_float

   integer(c_int), parameter :: METIS_OK = 1

   interface
      function METIS_PartMeshDual_C(ne, nn, eptr, eind, vwgt, vsize, &
                                    ncommon, nparts, tpwgts, options, &
                                    objval, epart, npart) &
                                    bind(C, name="METIS_PartMeshDual") result(ierr)
         import :: c_int, c_idx, c_real, c_ptr
         integer(c_idx), intent(inout) :: ne, nn
         integer(c_idx), intent(inout) :: eptr(*), eind(*)
         type(c_ptr),    value         :: vwgt, vsize
         integer(c_idx), intent(inout) :: ncommon, nparts
         type(c_ptr),    value         :: tpwgts
         integer(c_idx), intent(inout) :: options(*)
         integer(c_idx), intent(out)   :: objval
         integer(c_idx), intent(out)   :: epart(*), npart(*)
         integer(c_int) :: ierr
      end function METIS_PartMeshDual_C
   end interface

end module metis_binding
