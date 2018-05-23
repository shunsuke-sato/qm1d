module global_variables

! parameter
  integer,parameter :: dp = kind(0d0), zp = kind((1d0,1d0))
  complex(zp),parameter :: zI = (0d0,1d0)
  real(dp),parameter :: pi=3.14159265359d0

! mesh
  integer :: Nx
  real(dp) :: length_x,dx
  real(dp),allocatable :: xn(:)

! TD/GS: wave-function, density, potential and so on
  integer :: Nt_im
  real(dp) :: dt_im
  real(dp),allocatable :: rho(:), v_ext(:), v_int(:)
  complex(zp),allocatable :: zwfn(:,:), zw_mat(:,:,:), zSmat(:,:), zAmat(:,:)

  character(4) :: RT_mode
  integer :: Nt_iter
  real(dp) :: T_calc,dt,kick_mom
  real(dp),allocatable :: dipole_t(:),norm_t(:)
  real(dp) :: field_max,field_duration,field_omega
  real(dp) :: field_max_eV_per_AA,field_duration_fs,field_omega_eV
  real(dp),allocatable :: field_t(:)

! collision
  real(dp) :: Ekin0, x0_col, delta_x_col
  real(dp) :: V0_delta ! V(x) = V0_delta*delta(x-x')

! temporary
  complex(zp),allocatable :: ztmp_wfn(:,:),ztmp_hwfn(:,:),ztmp_swfn(:,:),ztmp_wfn_b(:,:)

! single-particle problem
  real(dp),allocatable :: wfn_sp(:)

  real(dp),parameter :: V0 = -0.73485720d0, sigma_V = 1.56906200d0

end module global_variables
!=======10========20========30========40========50========60========70========80========90=======100
program main
  use global_variables
  implicit none

  write(*,'(A)')'Start qm1d'
  call input
  call mesh

  call preparation_GS

  call single_particle_GS
!  call IMA_TIME_PROP

  call preparation_RT_collision
  call RT_prop

!  call write_results

end program main
!=======10========20========30========40========50========60========70========80========90=======100
subroutine input
  use global_variables
  implicit none

! == input parameter == !                                                                                                                       
  write(*,'(A)')'===== Input parameter ============================================================='
  write(*,'(A)')

  Nx = 2000
  length_x = 400d0
! GS
  Nt_im = 30000
  dt_im = 0.001d0
! TD
  T_calc = 400d0
  dt = 0.01d0
  Nt_iter = aint(T_calc/dt)+1


! collision
  Ekin0 = 0.60d0
  x0_col = -100d0
  delta_x_col = 15d0

! collision-option
!  V0_delta = 0d0

  RT_mode='kick' ! kick or gs
  kick_mom = 1e-3

  field_max_eV_per_AA = 1d0
  field_duration_fs = 10d0
  field_omega_eV = 1.55d0

  field_max = field_max_eV_per_AA * (0.529d0/27.2d0)
  field_duration = field_duration_fs/0.02418d0
  field_omega = field_omega_eV/27.2d0

  write(*,'(A,2x,I4)')'Nx =',Nx
  write(*,'(A,2x,e26.16e3)')'length_x =',length_x
  write(*,'(A,2x,e26.16e3)')'T_calc =',T_calc
  write(*,'(A,2x,e26.16e3)')'dt =',dt
  write(*,'(A,2x,I10)')'Nt_iter =',Nt_iter


! temporary array
  allocate(ztmp_wfn(0:Nx,2),ztmp_hwfn(0:Nx,2),ztmp_swfn(0:Nx,2),ztmp_wfn_b(-4:Nx+4,2))

  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== Complete Input parameter ===================================================' 
  return
end subroutine input
!=======10========20========30========40========50========60========70========80========90=======100
subroutine mesh
  use global_variables
  implicit none
  integer :: ix

  write(*,'(A)')'===== Making mesh ================================================================'
  write(*,'(A)')
  write(*,'(A)')

  allocate(xn(0:Nx))
  dx = length_x/dble(Nx)
  
  do ix = 0,Nx
     xn(ix) = dx*dble(ix) - 0.5d0*length_x
  end do

  write(*,'(A)')'===== Complete Making mesh ========================================================'
  return
end subroutine mesh
!=======10========20========30========40========50========60========70========80========90=======100
subroutine preparation_GS
  use global_variables
  implicit none
  integer :: ix,iy
  real(dp) :: tmp,r


  write(*,'(A)')'===== preparatin_GS =============================================================='
  write(*,'(A)')
  write(*,'(A)')

  allocate(v_ext(0:Nx), v_int(0:Nx))
  allocate(rho(0:Nx))

  allocate(zwfn(0:Nx,2), zw_mat(0:Nx,2,2), zSmat(2,2), zAmat(2,2))

  write(*,'(A)')'=== preparing initial wave-function ===='
! wfn(0,:) == wfn(Nx,:) == 0 
  zwfn(:,:) = 0d0

  do ix = 1,Nx-1
    zwfn(ix,1) = exp(-(xn(ix)-0.6d0)**2)
    zwfn(ix,2) = exp(-(xn(ix)+0.4d0)**2)
  end do


! Normalize
  zSmat(1,1) = sum(abs(zwfn(:,1))**2)*dx
  zSmat(2,2) = sum(abs(zwfn(:,2))**2)*dx
  zwfn(:,2) = zwfn(:,2)*sqrt(real(zSmat(1,1)/zSmat(2,2)))
  zSmat(2,2) = zSmat(1,1)
  zSmat(1,2) = sum(conjg(zwfn(:,1))*zwfn(:,2))*dx
  zSmat(2,1) = conjg(zSmat(1,2))
  tmp = zSmat(1,1)*zSmat(2,2) + abs(zSmat(1,2))**2
  zwfn(:,:)=zwfn(:,:)/sqrt(sqrt(tmp))


  write(*,'(A)')'=== preparing external potential ===='
  do ix=0,Nx
!     v_ext(ix) = 0.5d0*xn(ix)**2
!     v_ext(ix) = -2d0/sqrt(1d0+xn(ix)**2)
     v_ext(ix) = -1d0/sqrt(2d0+xn(ix)**2)
!     v_ext(ix) = V0*exp(-0.5d0*(xn(ix)/sigma_V)**2)
  end do

  write(*,'(A)')'=== preparing interaction potential ===='
  do ix=0,Nx
     r = dx*dble(ix)
!     v_int(ix) = 1d0 !1d0/sqrt(1d0+r**2)
     v_int(ix) = 1d0/sqrt(2d0+r**2)
!     v_int(ix) = exp(-r**2)
!     v_int(ix) = -V0*exp(-0.5d0*(r/sigma_V)**2)
  end do

  call wfn_rho

  write(*,'(A)')'=== preparing total potential ===='
  
  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== Complete preparatin_GS ====================================================='

  return
end subroutine preparation_GS
!=======10========20========30========40========50========60========70========80========90=======100
subroutine wfn_rho
  use global_variables
  implicit none
  integer :: ix

  zSmat(1,1) = sum(abs(zwfn(:,1))**2)*dx
  zSmat(2,2) = sum(abs(zwfn(:,2))**2)*dx
  zSmat(1,2) = sum(conjg(zwfn(:,1))*zwfn(:,2))*dx
  zSmat(2,1) = conjg(zSmat(1,2))

  rho(:)=0d0
  do ix = 0,Nx
     rho(ix) = 0.5d0*(abs(zwfn(ix,1))**2*zSmat(2,2) + abs(zwfn(ix,2))**2*zSmat(1,1))
  end do

  rho(:) = rho(:) + real( conjg(zwfn(:,2))*zwfn(:,1)*zSmat(1,2) )

  return
end subroutine wfn_rho
!=======10========20========30========40========50========60========70========80========90=======100
subroutine calc_zw_mat
  use global_variables
  implicit none
  integer :: ix,iy,ixy,i,j
  complex(zp) :: ztmp

! (1,1)
  do ix = 0,Nx
    ztmp = 0d0
    do iy = 0,Nx
      ixy = abs(ix-iy)
      ztmp = ztmp + v_int(ixy)*abs(zwfn(iy,1))**2
    end do
    zw_mat(ix,1,1) = ztmp*dx
  end do

! (1,2)
  do ix = 0,Nx
    ztmp = 0d0
    do iy = 0,Nx
      ixy = abs(ix-iy)
      ztmp = ztmp + v_int(ixy)*conjg(zwfn(iy,1))*zwfn(iy,2)
    end do
    zw_mat(ix,1,2) = ztmp*dx
  end do

! (2,2)
  do ix = 0,Nx
    ztmp = 0d0
    do iy = 0,Nx
      ixy = abs(ix-iy)
      ztmp = ztmp + v_int(ixy)*abs(zwfn(iy,2))**2
    end do
    zw_mat(ix,2,2) = ztmp*dx
  end do

  zw_mat(:,2,1) = conjg(zw_mat(:,1,2))

!delta-type additional force
  do i = 1,2
  do j = 1,2
    zw_mat(:,i,j) = zw_mat(:,i,j) !+ V0_delta*conjg(zwfn(:,i))*zwfn(:,j)
  end do
  end do

  return
end subroutine Calc_zw_mat
!=======10========20========30========40========50========60========70========80========90=======100
subroutine zhpsi_pre
  use global_variables
  implicit none
  real(dp) :: sigma
  complex(zp) :: zw_ijkl(2,2,2,2)
  integer :: i,j,k,l
  complex(zp) :: zt1,zt2

! Normalize
  zSmat(1,1) = sum(abs(zwfn(:,1))**2)*dx
  zSmat(2,2) = sum(abs(zwfn(:,2))**2)*dx
  zSmat(1,2) = sum(conjg(zwfn(:,1))*zwfn(:,2))*dx
  zSmat(2,1) = conjg(zSmat(1,2))
  sigma = zSmat(1,1)*zSmat(2,2)-abs(zSmat(1,2))**2

  do i = 1,2
  do j = 1,2
  do k = 1,2
  do l = 1,2
!    if (i == 1 .and. j == 1 .and. k == 1 .and. l == 1)cycle
!    if (i == 2 .and. j == 2 .and. k == 2 .and. l == 2)cycle
    zw_ijkl(i,j,k,l) = sum(conjg(zwfn(:,i))*zwfn(:,k)*zw_mat(:,j,l))*dx
  end do
  end do
  end do
  end do



  zt1 = zSmat(1,1)*zSmat(2,2)*(zw_ijkl(1,2,2,1) + zw_ijkl(1,2,1,2))
  zt2 = zSmat(2,1)*zSmat(1,1)*zw_ijkl(1,1,1,2) + zSmat(1,2)*zSmat(1,1)*zw_ijkl(2,2,1,2)
  zAmat(1,1) = (zt1 - zt2)/((zSmat(1,1)+zSmat(2,2))*sigma)
  zAmat(2,2) = zAmat(1,1)


  zt1 = (zSmat(1,1)*zSmat(2,2)+zSmat(2,2)**2-abs(zSmat(1,2))**2)*zw_ijkl(1,1,1,2) &
    +zSmat(1,2)**2*zw_ijkl(2,2,1,2)
  zt2 = zSmat(2,2)*zSmat(1,2)*(zw_ijkl(1,2,2,1) + zw_ijkl(1,2,1,2))
  zAmat(1,2) = (zt1 - zt2)/((zSmat(1,1)+zSmat(2,2))*sigma)

  zt1 = (zSmat(1,1)*zSmat(2,2)+zSmat(1,1)**2-abs(zSmat(1,2))**2)*zw_ijkl(2,2,1,2) &
    +zSmat(2,1)**2*zw_ijkl(1,1,1,2)
  zt2 = zSmat(1,1)*zSmat(2,1)*(zw_ijkl(1,2,2,1) + zw_ijkl(1,2,1,2))
  zAmat(2,1) = (zt1 - zt2)/((zSmat(1,1)+zSmat(2,2))*sigma)


  return
end subroutine zhpsi_pre
!=======10========20========30========40========50========60========70========80========90=======100
subroutine zhpsi(ft)
  use global_variables
  implicit none
! finite difference
  real(dp),parameter :: cN0=-205d0/72d0,cN1=8d0/5d0
  real(dp),parameter :: cN2=-1d0/5d0,cN3=8d0/315d0
  real(dp),parameter :: cN4=-1d0/560d0    
  real(dp) :: ft
  integer :: ix,iy
  real(dp) :: c0,c1,c2,c3,c4
  real(dp) :: sigma
  complex(zp) :: zSmat_i(2,2)
  integer :: n
! nine-points formula  
  c0=-0.5d0*cN0/(dx**2)
  c1=-0.5d0*cN1/(dx**2)
  c2=-0.5d0*cN2/(dx**2)
  c3=-0.5d0*cN3/(dx**2)
  c4=-0.5d0*cN4/(dx**2)
  sigma = zSmat(1,1)*zSmat(2,2)-abs(zSmat(1,2))**2
  zSmat_i(1,1) = zSmat(2,2)/sigma
  zSmat_i(1,2) = -zSmat(1,2)/sigma
  zSmat_i(2,1) = -zSmat(2,1)/sigma
  zSmat_i(2,2) = zSmat(1,1)/sigma

  ztmp_wfn_b(:,:)=0d0
  ztmp_wfn_b(0:Nx,:) = ztmp_wfn(0:Nx,:)
  do n = 1,2
    do ix=0,Nx
      ztmp_hwfn(ix,n)=c0*ztmp_wfn_b(ix,n) &
        + c1*(ztmp_wfn_b(ix+1,n) + ztmp_wfn_b(ix-1,n)) &
        + c2*(ztmp_wfn_b(ix+2,n) + ztmp_wfn_b(ix-2,n)) &
        + c3*(ztmp_wfn_b(ix+3,n) + ztmp_wfn_b(ix-3,n)) &
        + c4*(ztmp_wfn_b(ix+4,n) + ztmp_wfn_b(ix-4,n)) 
    end do
  ztmp_hwfn(:,n) = ztmp_hwfn(:,n) + (v_ext(:)+ft*xn(:))*ztmp_wfn(:,n)
  end do


  do ix = 0,Nx
    ztmp_swfn(ix,1) = (zAmat(1,1)-zw_mat(ix,1,1))*ztmp_wfn(ix,2) &
      & + (zAmat(1,2)-zw_mat(ix,1,2))*ztmp_wfn(ix,1)  
    ztmp_swfn(ix,2) = (zAmat(2,1)-zw_mat(ix,2,1))*ztmp_wfn(ix,2) &
      & + (zAmat(2,2)-zw_mat(ix,2,2))*ztmp_wfn(ix,1)  
  end do

  ztmp_wfn(:,1) = zSmat_i(1,1)*ztmp_swfn(:,1) + zSmat_i(1,2)*ztmp_swfn(:,2)
  ztmp_wfn(:,2) = zSmat_i(2,1)*ztmp_swfn(:,1) + zSmat_i(2,2)*ztmp_swfn(:,2)

  ztmp_hwfn(:,1) = ztmp_hwfn(:,1) - ztmp_wfn(:,2)
  ztmp_hwfn(:,2) = ztmp_hwfn(:,2) - ztmp_wfn(:,1)

  return
end subroutine zhpsi
!=======10========20========30========40========50========60========70========80========90=======100
subroutine Ima_time_prop
  use global_variables
  implicit none
  integer :: it,ix
  real(dp) :: esp,esp_res,diff_rho,Etot,ft,tmp
  real(dp),allocatable :: rho_old(:)
  
  write(*,'(A)')'===== Ground state calculation ===================================================='
  write(*,'(A)')
  write(*,'(A)')
  allocate(rho_old(0:Nx))
  ft = 0d0
  call wfn_rho
  rho_old(:) = rho(:)
  call calc_zw_mat
  call zhpsi_pre

  open(90,file = 'GS_data.out')
  do it = 1,Nt_im
    write(*,'(A,2x,I7,2x,A,I7)')'it =',it,'/',Nt_im

    ztmp_wfn(:,:) = zwfn(:,:)
    call zhpsi(ft)
    zwfn(:,:) = zwfn(:,:) -dt_im*ztmp_hwfn(:,:)

! Normalize
    zSmat(1,1) = sum(abs(zwfn(:,1))**2)*dx
    zSmat(2,2) = sum(abs(zwfn(:,2))**2)*dx
    zwfn(:,2) = zwfn(:,2)*sqrt(real(zSmat(1,1)/zSmat(2,2)))
    zSmat(2,2) = zSmat(1,1)
    zSmat(1,2) = sum(conjg(zwfn(:,1))*zwfn(:,2))*dx
    zSmat(2,1) = conjg(zSmat(1,2))
    tmp = zSmat(1,1)*zSmat(2,2) + abs(zSmat(1,2))**2
    zwfn(:,:)=zwfn(:,:)/sqrt(sqrt(tmp))
    call calc_zw_mat
    call zhpsi_pre

    call wfn_rho
    diff_rho = sum((rho(:)-rho_old(:))**2)*dx
    rho_old = rho
    call Total_Energy(Etot,ft)
    write(90,'(I6,2x,3e28.16e3)')it,Etot,diff_rho
    write(*,'(I6,2x,3e28.16e3)')it,Etot,diff_rho

  end do
  close(90)


  open(90,file = 'GS_wfn.out')
  write(90,'(A)')'# x, wfn(x,y)'
  do ix = 0,Nx
     write(90,'(100e26.16e3)')xn(ix),abs(zwfn(ix,1))**2,abs(zwfn(ix,2))**2,v_ext(ix),rho(ix)
  end do
  close(90)

  call pair_distribution_GS

  write(*,'(A,2x,e26.16e3)')'Total Energy (a.u.) =',Etot
  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== End Ground state calculation ================================================'  
  
  return
end subroutine IMA_TIME_PROP
!=======10========20========30========40========50========60========70========80========90=======100
subroutine Total_Energy(Etot,ft)
  use global_variables
  implicit none
  real(dp) :: Etot,ft
  complex(zp) :: zHmat(2,2)

! finite difference
  real(dp),parameter :: cN0=-205d0/72d0,cN1=8d0/5d0
  real(dp),parameter :: cN2=-1d0/5d0,cN3=8d0/315d0
  real(dp),parameter :: cN4=-1d0/560d0    
  integer :: ix,iy,n
  real(dp) :: c0,c1,c2,c3,c4
  real(dp) :: sigma
  complex(zp) :: zSmat_i(2,2)
! nine-points formula  
  c0=-0.5d0*cN0/(dx**2)
  c1=-0.5d0*cN1/(dx**2)
  c2=-0.5d0*cN2/(dx**2)
  c3=-0.5d0*cN3/(dx**2)
  c4=-0.5d0*cN4/(dx**2)
  sigma = zSmat(1,1)*zSmat(2,2)-abs(zSmat(1,2))**2
  zSmat_i(1,1) = zSmat(2,2)/sigma
  zSmat_i(1,2) = -zSmat(1,2)/sigma
  zSmat_i(2,1) = -zSmat(2,1)/sigma
  zSmat_i(2,2) = zSmat(1,1)/sigma

  ztmp_wfn_b(:,:)=0d0
  ztmp_wfn_b(0:Nx,:) = zwfn(0:Nx,:)
  do n = 1,2
    do ix=0,Nx
      ztmp_hwfn(ix,n)=c0*ztmp_wfn_b(ix,n) &
        + c1*(ztmp_wfn_b(ix+1,n) + ztmp_wfn_b(ix-1,n)) &
        + c2*(ztmp_wfn_b(ix+2,n) + ztmp_wfn_b(ix-2,n)) &
        + c3*(ztmp_wfn_b(ix+3,n) + ztmp_wfn_b(ix-3,n)) &
        + c4*(ztmp_wfn_b(ix+4,n) + ztmp_wfn_b(ix-4,n)) 
    end do
  ztmp_hwfn(:,n) = ztmp_hwfn(:,n) + (v_ext(:)+ft*xn(:))*zwfn(:,n)
  end do


  zHmat(1,1) = sum(conjg(zwfn(:,1))*ztmp_hwfn(:,1))*dx
  zHmat(2,2) = sum(conjg(zwfn(:,2))*ztmp_hwfn(:,2))*dx
  zHmat(1,2) = sum(conjg(zwfn(:,1))*ztmp_hwfn(:,2))*dx
  zHmat(2,1) = conjg(zHmat(1,2))

  Etot = zHmat(1,1)*zSmat(2,2) + zHmat(2,2)*zSmat(1,1) &
     & + zHmat(1,2)*zSmat(2,1) + zHmat(2,1)*zSmat(1,2) &
     & + sum(conjg(zwfn(:,1))*zwfn(:,1)*zw_mat(:,2,2))*dx &
     & + sum(conjg(zwfn(:,1))*zwfn(:,2)*zw_mat(:,2,1))*dx 

  return
end subroutine Total_Energy
!=======10========20========30========40========50========60========70========80========90=======100
!=======10========20========30========40========50========60========70========80========90=======100
subroutine pair_distribution_GS
  use global_variables
  implicit none 
  real(dp),parameter :: length_limit_pair = 20d0
  integer :: Nlength_limit_pair,N0
  real(dp),allocatable :: Pair_dist(:,:)
  integer :: ix,iy,ir
  complex(zp),allocatable :: zwfn_2d(:,:)

  Nlength_limit_pair = aint(length_limit_pair /dx)+1
  allocate(Pair_dist(0:Nx,0:4))
  allocate(zwfn_2d(0:Nx,0:Nx))
  call wfn_rho
  rho = rho * 2d0

  do ix = 1,Nx
  do iy = 1,Nx
    zwfn_2d(ix,iy) = 1d0/sqrt(2d0)*( &
      &zwfn(ix,1)*zwfn(iy,2)+zwfn(ix,2)*zwfn(iy,1))
  end do
  end do


! (x+y)/2 = 0.0 a.u.
  N0 = Nx/2  
  write(*,*)'Pair 0 =',xn(N0)
  do ir = 0,Nx
    ix = N0 + ir
    iy = N0 - ir
    Pair_dist(ix-iy,0) = 2d0*abs(zwfn_2d(ix,iy))**2/(rho(ix)*rho(iy))
    if(dble(ix-iy)*dx > length_limit_pair)exit
  end do
  write(*,*)ix-iy

! (x+y)/2 = 1.0 a.u.
  N0 = Nx/2 + aint(1d0/dx)  
  write(*,*)'Pair 1 =',xn(N0)
  do ir = 0,Nx
    ix = N0 + ir
    iy = N0 - ir
    Pair_dist(ix-iy,1) = 2d0*abs(zwfn_2d(ix,iy))**2/(rho(ix)*rho(iy))
    if(dble(ix-iy)*dx > length_limit_pair)exit
  end do

! (x+y)/2 = 2.0 a.u.
  N0 = Nx/2 + aint(2d0/dx)  
  write(*,*)'Pair 1 =',xn(N0)
  do ir = 0,Nx
    ix = N0 + ir
    iy = N0 - ir
    Pair_dist(ix-iy,2) = 2d0*abs(zwfn_2d(ix,iy))**2/(rho(ix)*rho(iy))
    if(dble(ix-iy)*dx > length_limit_pair)exit
  end do

! (x+y)/2 = 3.0 a.u.
  N0 = Nx/2 + aint(3d0/dx)  
  write(*,*)'Pair 1 =',xn(N0)
  do ir = 0,Nx
    ix = N0 + ir
    iy = N0 - ir
    Pair_dist(ix-iy,3) = 2d0*abs(zwfn_2d(ix,iy))**2/(rho(ix)*rho(iy))
    if(dble(ix-iy)*dx > length_limit_pair)exit
  end do

! (x+y)/2 = 4.0 a.u.
  N0 = Nx/2 + aint(4d0/dx)  
  write(*,*)'Pair 1 =',xn(N0)
  do ir = 0,Nx
    ix = N0 + ir
    iy = N0 - ir
    Pair_dist(ix-iy,4) = 2d0*abs(zwfn_2d(ix,iy))**2/(rho(ix)*rho(iy))
    if(dble(ix-iy)*dx > length_limit_pair)exit
  end do

  open(90,file='Pair_dist_rel.out')
  write(90,'(A)')'# x-y,  P(R=0 au), P(R=1 au), P(R=2 au), P(R=3 au), P(R=4 au)'
  do ir = 0,Nx,2
    write(90,'(100e26.16e3)')dble(ir)*dx,Pair_dist(ir,:)
    if(dble(ir)*dx > length_limit_pair)exit
  end do
  close(90)
  write(*,*)ir
  return
end subroutine pair_distribution_GS
!=======10========20========30========40========50========60========70========80========90=======100
subroutine single_particle_GS
  use global_variables
  implicit none
  integer,parameter :: Ncg = 2000
  real(8),parameter :: epsilon_res = 1d-8
! CG
  real(dp),allocatable :: xvec(:),pvec(:),rvec(:)
  real(dp),allocatable :: hxvec(:),gvec(:),hpvec(:)

  real(dp) :: xx,pp,xp,xhx,php,xhp,esp,esp_res,gg,gg0
  real(dp) :: ss,lambda,alpha,beta,aa,bb,cc,res

  integer :: ix, icg
  real(dp) :: rnd


  allocate (wfn_sp(0:Nx))
  allocate( xvec(0:Nx),pvec(0:Nx),rvec(0:Nx) )
  allocate( hxvec(0:Nx),gvec(0:Nx),hpvec(0:Nx) )


  do ix = 0,Nx
    call random_number(rnd)
    wfn_sp(ix) = rnd
  end do

  ss = sum(wfn_sp(:)**2)*dx
  wfn_sp(:) = wfn_sp(:)/sqrt(ss)

  xvec(:) = wfn_sp(:)
  call hpsi_sp(xvec,hxvec,v_ext)
  xx=sum(xvec(:)**2)*dx
  xhx=sum(xvec(:)*hxvec(:))*dx
  lambda=xhx/xx

  do icg = 1, Ncg
    gvec(:)=2d0*(hxvec(:)-lambda*xvec(:))/xx
     
    gg0=sum(gvec(:)**2)*dx
    select case(icg)
    case(1)
      pvec(:)=-gvec(:)
    case default
      beta=gg0/gg
      pvec(:)=-gvec(:)+beta*pvec(:)
    end select
    gg=gg0
    
    call hpsi_sp(pvec,hpvec,v_ext)
    
    pp=sum(pvec(:)**2)*dx
    php=sum(pvec(:)*hpvec(:))*dx
    xp=sum(xvec(:)*pvec(:))*dx
    xhp=sum(hxvec(:)*pvec(:))*dx
    
    aa=php*xp-xhp*pp
    bb=php*xx-xhx*pp
    cc=xhp*xx-xhx*xp
    ss=bb**2-4d0*aa*cc
    if(ss > 0d0)then
      alpha=(-bb+sqrt(ss))/(2d0*aa)
    else
      exit
    end if
    
    xvec(:)=xvec(:)+alpha*pvec(:)
    
    call hpsi_sp(xvec,hxvec,v_ext)
    xx=sum(xvec(:)**2)*dx
    xhx=sum(xvec(:)*hxvec(:))*dx
    lambda=xhx/xx
    res = sum((hxvec(:)-lambda*xvec(:))**2)*dx
    write(*,'(I7,2x,2e26.16e3)')icg-1,lambda,res
    if(res < epsilon_res)exit
  end do


  xx=sum(xvec(:)**2)*dx
  wfn_sp(:) = xvec(:)/sqrt(xx)

  return
end subroutine single_particle_GS
!=======10========20========30========40========50========60========70========80========90=======100
subroutine hpsi_sp(u,hu,v)
  use global_variables
  implicit none
  real(dp) :: u(0:Nx),hu(0:Nx),tu(0-4:Nx+4),v(0:Nx)
! finite difference
  real(dp),parameter :: cN0=-205d0/72d0,cN1=8d0/5d0
  real(dp),parameter :: cN2=-1d0/5d0,cN3=8d0/315d0
  real(dp),parameter :: cN4=-1d0/560d0    
  real(dp) :: c0,c1,c2,c3,c4
  integer :: ix

  c0=-0.5d0*cN0/(dx**2)
  c1=-0.5d0*cN1/(dx**2)
  c2=-0.5d0*cN2/(dx**2)
  c3=-0.5d0*cN3/(dx**2)
  c4=-0.5d0*cN4/(dx**2)

  tu(0:Nx) = u(0:Nx)
  do ix = 0,Nx
    hu(ix) = c0*tu(ix) &
      + c1*(tu(ix+1)+tu(ix-1)) &
      + c2*(tu(ix+2)+tu(ix-2)) &
      + c3*(tu(ix+3)+tu(ix-3)) &
      + c4*(tu(ix+4)+tu(ix-4)) 
  end do

  hu(:) = hu(:) + v(:)*u(:)

  return
end subroutine hpsi_sp
!=======10========20========30========40========50========60========70========80========90=======100
subroutine preparation_RT_collision
  use global_variables
  implicit none
  integer :: ix
  real(dp) :: k_mom, tmp

  k_mom = sqrt(2d0*Ekin0)
  write(*,*)'k_mom =',k_mom
  zwfn(0:Nx,1) = wfn_sp(0:Nx)

  do ix = 0,Nx
    zwfn(ix,2) = exp(zI*k_mom*xn(ix))*exp(-0.5d0*((xn(ix)-x0_col)/delta_x_col)**2)
  end do

! Normalize
  zSmat(1,1) = sum(abs(zwfn(:,1))**2)*dx
  zSmat(2,2) = sum(abs(zwfn(:,2))**2)*dx
  zwfn(:,2) = zwfn(:,2)*sqrt(real(zSmat(1,1)/zSmat(2,2)))
  zSmat(2,2) = zSmat(1,1)
  zSmat(1,2) = sum(conjg(zwfn(:,1))*zwfn(:,2))*dx
  zSmat(2,1) = conjg(zSmat(1,2))
  tmp = zSmat(1,1)*zSmat(2,2) + abs(zSmat(1,2))**2
  zwfn(:,:)=zwfn(:,:)/sqrt(sqrt(tmp))

  call wfn_rho

  open(21,file='initial_wf_rho.out')
  do ix = 0,Nx
    write(21,'(999e26.16e3)')xn(ix),real(zwfn(ix,1)),aimag(zwfn(ix,1)) &
      ,real(zwfn(ix,2)),aimag(zwfn(ix,2)),rho(ix)
  end do
  close(21)

!  stop 'end prep wf'

  return
end subroutine preparation_RT_collision
!=======10========20========30========40========50========60========70========80========90=======100
subroutine RT_prop
  use global_variables
  implicit none
  integer :: it,iexp, ix
  real(dp) :: ft, Etot,P_impact_excitation(0:2)
  complex(zp) :: zcoef
  character(50) :: citer
! work (pred-corr)  
  complex(zp),allocatable :: zwfn_pre_cor(:,:)

  allocate( zwfn_pre_cor(0:Nx,2))

  write(*,'(A)')'===== Real time propagation ======================================================'
  write(*,'(A)')
  write(*,'(A)')

  ft = 0d0
  call calc_zw_mat
  call zhpsi_pre
  call update_v_ext

  open(21,file='consv_chk.out')
  call wfn_rho
  call Total_Energy(Etot,ft)
  write(21,'(A)')'# Time (a.u.), Norm = 2, Etot (a.u.)'
  write(21,'(99e26.16e3)')dt*dble(it),sum(rho)*dx,Etot
  write(*,'(A)')'# Time (a.u.), Norm = 2, Etot (a.u.)'
  write(*,'(99e26.16e3)')dt*dble(it),sum(rho)*dx,Etot

  it = 0
  write(citer,'(I7.7)')it
  open(22,file='wfn_rho_'//trim(citer)//'.out')
  do ix = 0,Nx
    write(22,'(999e26.16e3)')xn(ix),real(zwfn(ix,1)),aimag(zwfn(ix,1)) &
      ,real(zwfn(ix,2)),aimag(zwfn(ix,2)),rho(ix)
  end do
  close(22)

  open(23,file='impact_excitation.out')
  call impact_ionization(P_impact_excitation); write(23,'(999e26.16e3)')it*dt,P_impact_excitation
  
  do it = 1,Nt_iter
    if(mod(it,1) == 0)  then
      write(*,'(A,2x,I7,2x,f7.2,2x,A)')'it =',it,dble(it)/dble(Nt_iter)*100.0,'%'
    end if


! == predictor-ccorector
    zwfn_pre_cor(:,:) = zwfn(:,:)
! prediction t + dt/2
    zcoef = 1d0
    ztmp_wfn(:,:) = zwfn(:,:)
    do iexp = 1,4
      zcoef = zcoef * (-zI*dt*0.5d0)/dble(iexp)
      call zhpsi(ft)
      zwfn(:,:) = zwfn(:,:) +zcoef*ztmp_hwfn(:,:)
      ztmp_wfn(:,:) = ztmp_hwfn(:,:)
    end do
    
    call calc_zw_mat
    call zhpsi_pre
    call update_v_ext

    zwfn(:,:) = zwfn_pre_cor(:,:) 
! correction t + dt
    zcoef = 1d0
    ztmp_wfn(:,:) = zwfn(:,:)
    do iexp = 1,4
      zcoef = zcoef * (-zI*dt)/dble(iexp)
      call zhpsi(ft)
      zwfn(:,:) = zwfn(:,:) +zcoef*ztmp_hwfn(:,:)
      ztmp_wfn(:,:) = ztmp_hwfn(:,:)
    end do
    
    call calc_zw_mat
    call zhpsi_pre
    call update_v_ext
! == predictor-ccorector

     if(mod(it,50) == 0)then
       call impact_ionization(P_impact_excitation); write(23,'(999e26.16e3)')it*dt,P_impact_excitation
     end if

    if (mod(it,200) == 0)then
      call wfn_rho
      call Total_Energy(Etot,ft)
      write(21,'(99e26.16e3)')dt*dble(it),sum(rho)*dx,Etot
      write(*,'(A)')'conserv chk'
      write(*,'(99e26.16e3)')dt*dble(it),sum(rho)*dx,Etot

      write(citer,'(I7.7)')it
      open(22,file='wfn_rho_'//trim(citer)//'.out')
      do ix = 0,Nx
        write(22,'(999e26.16e3)')xn(ix),real(zwfn(ix,1)),aimag(zwfn(ix,1)) &
          ,real(zwfn(ix,2)),aimag(zwfn(ix,2)),rho(ix)
      end do
      close(22)
    end if

  end do

  close(21)
  close(23)
  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== complete Real time propagation ============================================'

  return
end subroutine RT_prop
!=======10========20========30========40========50========60========70========80========90=======100
subroutine update_v_ext
  use global_variables
  implicit none
  real(dp),parameter :: V0_pair =0d0
  real(dp) :: rho_pair(0:Nx)
  integer :: ix

  do ix = 0,Nx
    rho_pair(ix) = (abs(zwfn(ix,1))**2)* (abs(zwfn(ix,2))**2)
  end do

  do ix=0,Nx
     v_ext(ix) = V0*exp(-0.5d0*(xn(ix)/sigma_V)**2)
!     v_ext(ix) = -1d0/sqrt(1d0+xn(ix)**2) + V0_pair * rho_pair(ix)
  end do
  
  return
end subroutine update_v_ext
!=======10========20========30========40========50========60========70========80========90=======100
subroutine impact_ionization(P_impact_excitation)
  use global_variables
  implicit none
  real(dp) :: P_impact_excitation(0:2)
  complex(zp) :: zv_t0(0:Nx),zv_t(0:Nx),zc_t(0:2)
  integer :: ix,iy

  do ix = 0,Nx
    zv_t(ix) = (sum(wfn_sp(:)*zwfn(:,1))*dx*zwfn(ix,2) + sum(wfn_sp(:)*zwfn(:,2))*dx*zwfn(ix,1))/sqrt(2d0)
  end do

  zc_t(0) = sum(wfn_sp(:)*zv_t(:))*dx
  
  zv_t0(:) = sqrt(2d0)*(zv_t(:)-zc_t(0)*wfn_sp(:))

  zc_t(1) = sqrt(sum(abs(zv_t0(:))**2)*dx)

  P_impact_excitation(0) = abs(zc_t(0))**2
  P_impact_excitation(1) = abs(zc_t(1))**2
  P_impact_excitation(2) =   1d0 - ( P_impact_excitation(1) + P_impact_excitation(0) )

  return
end subroutine impact_ionization
!=======10========20========30========40========50========60========70========80========90=======100
!=======10========20========30========40========50========60========70========80========90=======100
!=======10========20========30========40========50========60========70========80========90=======100
!=======10========20========30========40========50========60========70========80========90=======100
