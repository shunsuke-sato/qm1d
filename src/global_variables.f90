!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module global_variables

! parameter
  integer,parameter :: dp = kind(0d0), zp = kind((1d0,1d0))
  complex(zp),parameter :: zI = (0d0,1d0)
  real(dp),parameter :: pi=3.14159265359d0

! mesh
  integer :: Nx
  real(dp) :: length_x,dx
  real(dp),allocatable :: xn(:),xyn(:,:)

  real(dp),parameter :: gN0 = 0d0
  real(dp),parameter :: gN1 = 4.d0/5.d0
  real(dp),parameter :: gN2 = -1.d0/5.d0
  real(dp),parameter :: gN3 = 4.d0/105.d0
  real(dp),parameter :: gN4 = -1.d0/280.d0
  real(dp),parameter :: cN0 = -205.d0/72.d0
  real(dp),parameter :: cN1 = 8.d0/5.d0
  real(dp),parameter :: cN2 = -1.d0/5.d0
  real(dp),parameter :: cN3 = 8.d0/315.d0
  real(dp),parameter :: cN4 = -1.d0/560.d0


! GS: wave-function, density, potential and so on
  integer :: Ncg
  real(dp),allocatable :: wfn(:,:), rho(:), v_ext(:), v_int(:,:), v_all(:,:)
  real(dp),allocatable :: jx(:)

! TD
  character(4) :: RT_mode
  integer :: Nt_iter
  real(dp) :: T_calc,dt,kick_mom
  complex(zp),allocatable :: zwfn(:,:)
  real(dp),allocatable :: dipole_t(:),norm_t(:)
  real(dp) :: field_max,field_duration,field_omega
  real(dp) :: field_max_eV_per_AA,field_duration_fs,field_omega_eV
  real(dp),allocatable :: field_t(:)
! TD-Lanczos
  integer,parameter :: NLanczos = 6
  complex(zp),allocatable :: zwfn_Lanczos(:,:,:)

! collision
  real(dp) :: Ekin0, x0_col, delta_x_col

! H2
  real(dp),parameter :: sigma0 = sqrt(2d0)
  real(dp),parameter :: r_dist = 2d0

! single-particle problem
  real(dp),allocatable :: wfn_sp(:)

! exact exchange-correlation potential
!! TDHF type
  integer,parameter :: Nc_exc = 10
  real(dp),allocatable :: alpha_xc(:,:), v_xc_tdhf(:)
  complex(dp),allocatable :: zwfn_st(:,:)

end module global_variables
