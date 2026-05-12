module io_vtk_legacy
   use kinds,         only : wp
   use constants,     only : NVAR, IRHO, IRHOU, IRHOV, IRHOW
   use mesh_types,    only : t_mesh, &
                             CELL_TET, CELL_HEX, CELL_PRISM, CELL_PYRAMID, &
                             n_vtx_per_cell
   use eos_ideal_gas, only : prim_from_cons
   use fields,        only : t_state
   implicit none
   private

   public :: write_vtk

contains

   subroutine write_vtk(filename, mesh, s, time)
      character(len=*), intent(in) :: filename
      type(t_mesh),     intent(in) :: mesh
      type(t_state),    intent(in) :: s
      real(wp),         intent(in) :: time

      integer :: u, ios, i, c, k, vptr, nv
      integer :: total_int, vtk_type
      real(wp) :: rho, ux, uy, uz, p

      open(newunit=u, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) then
         write(*,'(A,A)') 'write_vtk: cannot open ', trim(filename)
         return
      end if

      write(u,'(A)') '# vtk DataFile Version 3.0'
      write(u,'(A,1PE15.6)') 'cfd3d output, t = ', time
      write(u,'(A)') 'ASCII'
      write(u,'(A)') 'DATASET UNSTRUCTURED_GRID'

      ! Points
      write(u,'(A,I0,A)') 'POINTS ', mesh%nv, ' double'
      do i = 1, mesh%nv
         write(u,'(3(1X,1PE22.14))') mesh%xv(1,i), mesh%xv(2,i), mesh%xv(3,i)
      end do

      ! Cells
      total_int = 0
      do c = 1, mesh%nc_internal
         total_int = total_int + 1 + n_vtx_per_cell(mesh%cell_type(c))
      end do
      write(u,'(A,I0,1X,I0)') 'CELLS ', mesh%nc_internal, total_int
      do c = 1, mesh%nc_internal
         nv = n_vtx_per_cell(mesh%cell_type(c))
         vptr = mesh%cell_vtx_ptr(c) - 1
         write(u,'(I0)', advance='no') nv
         do k = 1, nv
            write(u,'(1X,I0)', advance='no') mesh%cell_vtx(vptr + k) - 1
         end do
         write(u,*)
      end do

      write(u,'(A,I0)') 'CELL_TYPES ', mesh%nc_internal
      do c = 1, mesh%nc_internal
         select case (mesh%cell_type(c))
         case (CELL_TET)
            vtk_type = 10
         case (CELL_HEX)
            vtk_type = 12
         case (CELL_PRISM)
            vtk_type = 13
         case (CELL_PYRAMID)
            vtk_type = 14
         case default
            vtk_type = 0
         end select
         write(u,'(I0)') vtk_type
      end do

      ! Cell data
      write(u,'(A,I0)') 'CELL_DATA ', mesh%nc_internal

      ! rho
      write(u,'(A)') 'SCALARS rho double 1'
      write(u,'(A)') 'LOOKUP_TABLE default'
      do c = 1, mesh%nc_internal
         call prim_from_cons(s%U(:, c), rho, ux, uy, uz, p)
         write(u,'(1PE22.14)') rho
      end do

      ! p
      write(u,'(A)') 'SCALARS pressure double 1'
      write(u,'(A)') 'LOOKUP_TABLE default'
      do c = 1, mesh%nc_internal
         call prim_from_cons(s%U(:, c), rho, ux, uy, uz, p)
         write(u,'(1PE22.14)') p
      end do

      ! velocity
      write(u,'(A)') 'VECTORS velocity double'
      do c = 1, mesh%nc_internal
         call prim_from_cons(s%U(:, c), rho, ux, uy, uz, p)
         write(u,'(3(1X,1PE22.14))') ux, uy, uz
      end do

      close(u)
   end subroutine write_vtk

end module io_vtk_legacy
