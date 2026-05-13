module solver_control
   use kinds,     only : wp
   use constants, only : BC_SLIP_WALL, BC_SYMMETRY, &
                         BC_SUPERSONIC_INLET, BC_SUPERSONIC_OUTLET, BC_FARFIELD, &
                         BC_NO_SLIP_WALL, BC_DIRICHLET
   implicit none
   private

   public :: t_run_params, read_namelist, bc_string_to_int

   integer, parameter, public :: MAX_PATCHES = 32

   type :: t_run_params
      character(len=256) :: mesh_file = ''
      character(len=64)  :: case_name = 'cfd3d'
      character(len=256) :: output_dir = '.'
      real(wp) :: t_end = 1.0_wp
      real(wp) :: cfl   = 0.5_wp
      integer  :: max_steps = 100000
      real(wp) :: output_interval = 0.1_wp

      ! Initial condition
      character(len=32) :: init_type = 'uniform'  ! 'uniform' | 'riemann' | 'mms'
      integer  :: diaphragm_axis = 1
      real(wp) :: diaphragm_pos  = 0.0_wp
      real(wp) :: rho_L = 1.0_wp, u_L = 0.0_wp, v_L = 0.0_wp, w_L = 0.0_wp, p_L = 1.0_wp
      real(wp) :: rho_R = 1.0_wp, u_R = 0.0_wp, v_R = 0.0_wp, w_R = 0.0_wp, p_R = 1.0_wp

      ! Scheme toggles (M3)
      logical  :: muscl_enabled    = .true.    ! 2nd-order MUSCL with Venkat limiter
      real(wp) :: venkat_K         = 5.0_wp    ! Venkat tuning constant
      logical  :: viscous_enabled  = .false.   ! compressible NS vs Euler
      logical  :: mms_enabled      = .false.   ! manufactured solution source + ref

      ! Gas thermophysical properties (M3)
      real(wp) :: R_gas       = 1.0_wp         ! specific gas constant (non-dim default)
      logical  :: sutherland  = .false.        ! false = constant μ
      real(wp) :: mu_const    = 0.0_wp
      real(wp) :: mu_ref      = 1.716e-5_wp    ! Sutherland μ at T_ref
      real(wp) :: T_ref       = 273.15_wp
      real(wp) :: S_S         = 110.4_wp
      real(wp) :: Pr          = 0.72_wp

      ! Patch BC mapping
      integer :: patch_count = 0
      character(len=64) :: patch_name(MAX_PATCHES) = ''
      character(len=32) :: patch_bc  (MAX_PATCHES) = ''
      real(wp) :: patch_rho(MAX_PATCHES) = 1.0_wp
      real(wp) :: patch_u  (MAX_PATCHES) = 0.0_wp
      real(wp) :: patch_v  (MAX_PATCHES) = 0.0_wp
      real(wp) :: patch_w  (MAX_PATCHES) = 0.0_wp
      real(wp) :: patch_p  (MAX_PATCHES) = 1.0_wp
   end type t_run_params

contains

   subroutine read_namelist(filename, p)
      character(len=*),   intent(in)  :: filename
      type(t_run_params), intent(out) :: p

      ! Locals mirror the public type for namelist binding
      character(len=256) :: mesh_file
      character(len=64)  :: case_name
      character(len=256) :: output_dir
      real(wp) :: t_end, cfl, output_interval
      integer  :: max_steps
      character(len=32) :: init_type
      integer  :: diaphragm_axis
      real(wp) :: diaphragm_pos
      real(wp) :: rho_L, u_L, v_L, w_L, p_L
      real(wp) :: rho_R, u_R, v_R, w_R, p_R
      logical  :: muscl_enabled, viscous_enabled, mms_enabled, sutherland
      real(wp) :: venkat_K
      real(wp) :: R_gas, mu_const, mu_ref, T_ref, S_S, Pr
      integer  :: patch_count
      character(len=64) :: patch_name(MAX_PATCHES)
      character(len=32) :: patch_bc  (MAX_PATCHES)
      real(wp) :: patch_rho(MAX_PATCHES)
      real(wp) :: patch_u  (MAX_PATCHES)
      real(wp) :: patch_v  (MAX_PATCHES)
      real(wp) :: patch_w  (MAX_PATCHES)
      real(wp) :: patch_p  (MAX_PATCHES)

      namelist /cfd3d/ mesh_file, case_name, output_dir, t_end, cfl, &
         output_interval, max_steps, init_type, diaphragm_axis, diaphragm_pos, &
         rho_L, u_L, v_L, w_L, p_L, rho_R, u_R, v_R, w_R, p_R, &
         muscl_enabled, venkat_K, viscous_enabled, mms_enabled, &
         R_gas, sutherland, mu_const, mu_ref, T_ref, S_S, Pr, &
         patch_count, patch_name, patch_bc, &
         patch_rho, patch_u, patch_v, patch_w, patch_p

      integer :: u, ios

      ! Seed defaults from a default-constructed t_run_params
      mesh_file       = p%mesh_file
      case_name       = p%case_name
      output_dir      = p%output_dir
      t_end           = p%t_end
      cfl             = p%cfl
      output_interval = p%output_interval
      max_steps       = p%max_steps
      init_type       = p%init_type
      diaphragm_axis  = p%diaphragm_axis
      diaphragm_pos   = p%diaphragm_pos
      rho_L = p%rho_L; u_L = p%u_L; v_L = p%v_L; w_L = p%w_L; p_L = p%p_L
      rho_R = p%rho_R; u_R = p%u_R; v_R = p%v_R; w_R = p%w_R; p_R = p%p_R
      muscl_enabled   = p%muscl_enabled
      venkat_K        = p%venkat_K
      viscous_enabled = p%viscous_enabled
      mms_enabled     = p%mms_enabled
      R_gas      = p%R_gas
      sutherland = p%sutherland
      mu_const   = p%mu_const
      mu_ref     = p%mu_ref
      T_ref      = p%T_ref
      S_S        = p%S_S
      Pr         = p%Pr
      patch_count = p%patch_count
      patch_name  = p%patch_name
      patch_bc    = p%patch_bc
      patch_rho   = p%patch_rho
      patch_u     = p%patch_u
      patch_v     = p%patch_v
      patch_w     = p%patch_w
      patch_p     = p%patch_p

      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(A,A)') 'read_namelist: cannot open ', trim(filename)
         error stop 1
      end if
      read(u, nml=cfd3d, iostat=ios)
      close(u)
      if (ios /= 0) then
         write(*,'(A,A,I0)') 'read_namelist: failed to parse ', trim(filename), ios
         error stop 1
      end if

      p%mesh_file       = mesh_file
      p%case_name       = case_name
      p%output_dir      = output_dir
      p%t_end           = t_end
      p%cfl             = cfl
      p%output_interval = output_interval
      p%max_steps       = max_steps
      p%init_type       = init_type
      p%diaphragm_axis  = diaphragm_axis
      p%diaphragm_pos   = diaphragm_pos
      p%rho_L = rho_L; p%u_L = u_L; p%v_L = v_L; p%w_L = w_L; p%p_L = p_L
      p%rho_R = rho_R; p%u_R = u_R; p%v_R = v_R; p%w_R = w_R; p%p_R = p_R
      p%muscl_enabled   = muscl_enabled
      p%venkat_K        = venkat_K
      p%viscous_enabled = viscous_enabled
      p%mms_enabled     = mms_enabled
      p%R_gas      = R_gas
      p%sutherland = sutherland
      p%mu_const   = mu_const
      p%mu_ref     = mu_ref
      p%T_ref      = T_ref
      p%S_S        = S_S
      p%Pr         = Pr
      p%patch_count = patch_count
      p%patch_name  = patch_name
      p%patch_bc    = patch_bc
      p%patch_rho   = patch_rho
      p%patch_u     = patch_u
      p%patch_v     = patch_v
      p%patch_w     = patch_w
      p%patch_p     = patch_p
   end subroutine read_namelist

   pure function bc_string_to_int(s) result(it)
      character(len=*), intent(in) :: s
      integer :: it
      select case (trim(adjustl(s)))
      case ('slip_wall', 'wall')
         it = BC_SLIP_WALL
      case ('symmetry')
         it = BC_SYMMETRY
      case ('supersonic_inlet', 'inlet')
         it = BC_SUPERSONIC_INLET
      case ('supersonic_outlet', 'outlet')
         it = BC_SUPERSONIC_OUTLET
      case ('farfield')
         it = BC_FARFIELD
      case ('no_slip_wall', 'no_slip')
         it = BC_NO_SLIP_WALL
      case ('dirichlet')
         it = BC_DIRICHLET
      case default
         it = BC_SLIP_WALL
      end select
   end function bc_string_to_int

end module solver_control
