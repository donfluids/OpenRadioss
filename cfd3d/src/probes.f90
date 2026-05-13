! Point-gauge ("probe") output for blast diagnostics.
!
! Reads a CSV file of probe locations:
!     name, x, y, z
!     wall_centre, 0.5, 0.5, 0.0
!     ...
! Each rank computes the nearest local cell centroid to each probe; the rank
! whose distance is smallest globally is declared the probe's owner. Owners
! record the per-step primitive state at the probe's nearest cell into a
! per-probe history buffer. At shutdown each owner flushes its probes to a
! per-probe CSV file (`<case>_<probe_name>.csv`).
module probes
   use kinds,         only : wp
   use constants,     only : NPRIM, IP_RHO, IP_U, IP_V, IP_W, IP_P
   use mpi_f08
   use mpi_runtime,   only : t_mpi_ctx
   use mesh_types,    only : t_mesh
   use fields,        only : t_state
   use eos_ideal_gas, only : prim_from_cons
   implicit none
   private

   public :: t_probe_set
   public :: probes_load, probes_locate, probes_sample, probes_flush, probes_free

   integer, parameter :: PROBE_NAME_LEN = 32
   integer, parameter :: INIT_CAPACITY  = 1024

   type :: t_probe
      character(len=PROBE_NAME_LEN) :: name = ''
      real(wp) :: x(3)        = 0.0_wp
      integer  :: owner_rank  = -1
      integer  :: cell        = 0
   end type t_probe

   type :: t_probe_set
      integer :: nprobes   = 0
      integer :: nsamples  = 0
      integer :: ncapacity = 0
      type(t_probe), allocatable :: pr(:)
      real(wp), allocatable :: t_hist(:)         ! (ncapacity)
      real(wp), allocatable :: w_hist(:, :, :)   ! (NPRIM, ncapacity, nprobes)
      character(len=256) :: case_name = 'case'
      character(len=256) :: out_dir   = '.'
   end type t_probe_set

contains

   ! Parse the probes CSV — comments start with '#'. If the path is empty,
   ! the probe set stays empty.
   subroutine probes_load(filename, case_name, out_dir, ctx, ps)
      character(len=*),   intent(in)  :: filename, case_name, out_dir
      type(t_mpi_ctx),    intent(in)  :: ctx
      type(t_probe_set),  intent(out) :: ps

      integer :: u, ios, n, i, comma1, comma2, comma3
      character(len=256) :: line, namepart
      type(t_probe), allocatable :: tmp(:)
      real(wp) :: x, y, z

      ps%case_name = case_name
      ps%out_dir   = out_dir
      ps%nprobes   = 0

      if (len_trim(filename) == 0) return

      open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         if (ctx%is_root) &
            write(*,'(A,A)') 'probes_load: no probes file, skipping: ', trim(filename)
         return
      end if

      n = 0
      allocate(tmp(64))
      do
         read(u, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         comma1 = index(line, ',')
         if (comma1 <= 1) cycle
         comma2 = index(line(comma1+1:), ',')
         if (comma2 <= 0) cycle
         comma2 = comma1 + comma2
         comma3 = index(line(comma2+1:), ',')
         if (comma3 <= 0) cycle
         comma3 = comma2 + comma3
         namepart = adjustl(line(:comma1-1))
         read(line(comma1+1 : comma2-1), *, iostat=ios) x
         if (ios /= 0) cycle
         read(line(comma2+1 : comma3-1), *, iostat=ios) y
         if (ios /= 0) cycle
         read(line(comma3+1:), *, iostat=ios) z
         if (ios /= 0) cycle
         n = n + 1
         if (n > size(tmp)) call grow(tmp, 2 * size(tmp))
         tmp(n)%name = trim(namepart)
         tmp(n)%x = [x, y, z]
      end do
      close(u)

      ps%nprobes = n
      if (n > 0) then
         allocate(ps%pr(n))
         ps%pr = tmp(1:n)
         ps%ncapacity = INIT_CAPACITY
         allocate(ps%t_hist(ps%ncapacity))
         allocate(ps%w_hist(NPRIM, ps%ncapacity, n))
         ps%t_hist = 0.0_wp
         ps%w_hist = 0.0_wp
      end if
      deallocate(tmp)

      if (ctx%is_root) &
         write(*,'(A,I0,A,A)') 'probes_load: loaded ', n, ' probes from ', trim(filename)
   end subroutine probes_load

   subroutine grow(arr, new_size)
      type(t_probe), allocatable, intent(inout) :: arr(:)
      integer,                    intent(in)    :: new_size
      type(t_probe), allocatable :: tmp(:)
      integer :: n_old
      n_old = size(arr)
      allocate(tmp(new_size))
      tmp(1:n_old) = arr
      call move_alloc(tmp, arr)
   end subroutine grow

   ! For each probe, find the closest local cell centroid; MPI_MINLOC on the
   ! (distance, rank) pair determines the global owner.
   subroutine probes_locate(mesh, ctx, ps)
      type(t_mesh),       intent(in)    :: mesh
      type(t_mpi_ctx),    intent(in)    :: ctx
      type(t_probe_set),  intent(inout) :: ps

      integer  :: i, c, c_best
      real(wp) :: d2, d2_best, dx, dy, dz
      type :: minloc_t
         real(wp) :: d
         integer  :: rank
      end type minloc_t
      type(minloc_t) :: local, global

      do i = 1, ps%nprobes
         c_best = 0
         d2_best = huge(1.0_wp)
         do c = 1, mesh%nc_internal
            dx = mesh%cell_centroid(1, c) - ps%pr(i)%x(1)
            dy = mesh%cell_centroid(2, c) - ps%pr(i)%x(2)
            dz = mesh%cell_centroid(3, c) - ps%pr(i)%x(3)
            d2 = dx*dx + dy*dy + dz*dz
            if (d2 < d2_best) then
               d2_best = d2
               c_best  = c
            end if
         end do
         local%d    = sqrt(max(d2_best, 0.0_wp))
         local%rank = ctx%rank
         if (ctx%nproc > 1) then
            call MPI_Allreduce(local, global, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, ctx%comm)
         else
            global = local
         end if
         ps%pr(i)%owner_rank = global%rank
         if (global%rank == ctx%rank) then
            ps%pr(i)%cell = c_best
         else
            ps%pr(i)%cell = 0
         end if
      end do
   end subroutine probes_locate

   ! Each owner rank appends (t, W_at_cell) to its history buffer.
   subroutine probes_sample(mesh, s, t, ctx, ps)
      type(t_mesh),       intent(in)    :: mesh
      type(t_state),      intent(in)    :: s
      real(wp),           intent(in)    :: t
      type(t_mpi_ctx),    intent(in)    :: ctx
      type(t_probe_set),  intent(inout) :: ps

      integer  :: i, c
      real(wp) :: rho, u, v, w, p

      if (ps%nprobes == 0) return
      if (ps%nsamples + 1 > ps%ncapacity) call grow_history(ps)

      ps%nsamples = ps%nsamples + 1
      ps%t_hist(ps%nsamples) = t
      do i = 1, ps%nprobes
         if (ps%pr(i)%owner_rank == ctx%rank .and. ps%pr(i)%cell > 0) then
            c = ps%pr(i)%cell
            call prim_from_cons(s%U(:, c), rho, u, v, w, p)
            ps%w_hist(IP_RHO, ps%nsamples, i) = rho
            ps%w_hist(IP_U,   ps%nsamples, i) = u
            ps%w_hist(IP_V,   ps%nsamples, i) = v
            ps%w_hist(IP_W,   ps%nsamples, i) = w
            ps%w_hist(IP_P,   ps%nsamples, i) = p
         end if
      end do
      ! Suppress unused-warning
      if (.false.) i = mesh%nc_total
   end subroutine probes_sample

   subroutine grow_history(ps)
      type(t_probe_set), intent(inout) :: ps
      real(wp), allocatable :: tmp_t(:), tmp_w(:,:,:)
      integer :: new_cap
      new_cap = max(2 * ps%ncapacity, INIT_CAPACITY)
      allocate(tmp_t(new_cap))
      tmp_t(1:ps%nsamples) = ps%t_hist(1:ps%nsamples)
      call move_alloc(tmp_t, ps%t_hist)
      allocate(tmp_w(NPRIM, new_cap, ps%nprobes))
      tmp_w(:, 1:ps%nsamples, :) = ps%w_hist(:, 1:ps%nsamples, :)
      call move_alloc(tmp_w, ps%w_hist)
      ps%ncapacity = new_cap
   end subroutine grow_history

   ! Each owner rank writes its probes to per-probe CSV files. Files are
   ! named `<case_name>_probe_<name>.csv` and live in out_dir.
   subroutine probes_flush(ps, ctx)
      type(t_probe_set), intent(in) :: ps
      type(t_mpi_ctx),   intent(in) :: ctx
      integer :: i, u, ios, k
      character(len=512) :: fname

      do i = 1, ps%nprobes
         if (ps%pr(i)%owner_rank /= ctx%rank) cycle
         fname = trim(ps%out_dir) // '/' // trim(ps%case_name) // &
                 '_probe_' // trim(ps%pr(i)%name) // '.csv'
         open(newunit=u, file=trim(fname), status='replace', action='write', iostat=ios)
         if (ios /= 0) then
            write(*,'(A,A)') 'probes_flush: cannot open ', trim(fname)
            cycle
         end if
         write(u,'(A)') '# time, rho, u, v, w, p'
         do k = 1, ps%nsamples
            write(u,'(1PE16.8,5(",",1PE16.8))') &
               ps%t_hist(k), ps%w_hist(IP_RHO,k,i), ps%w_hist(IP_U,k,i), &
               ps%w_hist(IP_V,k,i),   ps%w_hist(IP_W,k,i), ps%w_hist(IP_P,k,i)
         end do
         close(u)
      end do
   end subroutine probes_flush

   subroutine probes_free(ps)
      type(t_probe_set), intent(inout) :: ps
      if (allocated(ps%pr))     deallocate(ps%pr)
      if (allocated(ps%t_hist)) deallocate(ps%t_hist)
      if (allocated(ps%w_hist)) deallocate(ps%w_hist)
      ps%nprobes   = 0
      ps%nsamples  = 0
      ps%ncapacity = 0
   end subroutine probes_free

end module probes
