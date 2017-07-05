!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine dt_evolve(ft)
  use global_variables
  use hpsi
  implicit none
  integer :: iexp
  real(dp) :: ft
  complex(zp) :: zs
  complex(zp) :: ztmp_wfn(0:Nx,0:Nx),ztmp_hwfn(0:Nx,0:Nx)

  ztmp_wfn(:,:) = zwfn(:,:)
  zs = 1d0
  do iexp = 1,4
     zs = zs*(-zI*dt)/dble(iexp)
     call zhpsi(ztmp_wfn,ztmp_hwfn,ft)
     zwfn(:,:) = zwfn(:,:) + zs*ztmp_hwfn(:,:)
     ztmp_wfn(:,:) = ztmp_hwfn(:,:)
  end do
  return
end subroutine dt_evolve
!=======10========20========30========40========50========60========70========80========90=======100
subroutine dt_evolve_Lanczos(ft)
  use global_variables
  use hpsi
  implicit none
  integer :: iexp,j
  real(dp) :: ft
  complex(zp) :: zs
  complex(zp) :: Hamiltonian_L(NLanczos,NLanczos), Exp_Ham_L(NLanczos,NLanczos), Umat_L(NLanczos,NLanczos)
  complex(zp) :: zvec(NLanczos),zvec_t(NLanczos)

  real(dp) :: ss
!LAPACK
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info
  integer :: ix,iy
  complex(zp) :: ztmp_wfn(0:Nx,0:Nx),ztmp_hwfn(0:Nx,0:Nx)

  lwork=6*NLanczos
  allocate(work_lp(lwork),rwork(3*NLanczos-2),w(NLanczos))
  ft = 0d0

  ss = 0d0
!$omp parallel do private(ix,iy) reduction(+:ss)
  do iy = 0,Nx
     do ix = 0,Nx
        ss = ss + abs(zwfn(ix,iy))**2
     end do
  end do
  ss = dx**2/sqrt(ss)

!$omp parallel do private(ix,iy)
  do iy = 0,Nx
     do ix = 0,Nx
        zwfn(ix,iy) = zwfn(ix,iy) * ss
        zwfn_Lanczos(ix,iy,1) = zwfn(ix,iy)
     end do
  end do

  Hamiltonian_L = 0d0

  do j = 1,NLanczos
!$omp parallel do private(ix,iy)
    do iy = 0,Nx
      do ix = 0,Nx
        ztmp_wfn(ix,iy) = zwfn_Lanczos(ix,iy,j)
      end do
    end do

    call zhpsi(ztmp_wfn,ztmp_hwfn,ft)
    ss = 0d0
!$omp parallel do private(ix,iy) reduction(+:ss)
    do iy = 0,Nx
       do ix = 0,Nx
          ss = ss + conjg(ztmp_wfn(ix,iy))*ztmp_hwfn(ix,iy)
       end do
    end do

    Hamiltonian_L(j,j) = ss*dx**2

    if(j == NLanczos)exit
    if(j == 1) then
!$omp parallel do private(ix,iy) 
      do iy = 0,Nx
        do ix = 0,Nx
          ztmp_wfn(ix,iy) = ztmp_hwfn(ix,iy) - Hamiltonian_L(j,j)*zwfn_Lanczos(ix,iy,j)
        end do
      end do
    else
!$omp parallel do private(ix,iy) 
      do iy = 0,Nx
        do ix = 0,Nx
          ztmp_wfn(ix,iy) = ztmp_hwfn(ix,iy) &
               - Hamiltonian_L(j,j)*zwfn_Lanczos(ix,iy,j)  &
               & -Hamiltonian_L(j,j-1)*zwfn_Lanczos(ix,iy,j-1)
        end do
      end do
    end if

    ss = 0d0
!$omp parallel do private(ix,iy) reduction(+:ss)
    do iy = 0,Nx
       do ix = 0,Nx
          ss = ss + abs(ztmp_wfn(ix,iy))**2
       end do
    end do
    ss = sqrt(ss*dx**2)
    Hamiltonian_L(j,j+1) = ss
    Hamiltonian_L(j+1,j) = ss

!$omp parallel do private(ix,iy) 
    do iy = 0,Nx
       do ix = 0,Nx
          zwfn_Lanczos(ix,iy,j+1) = ztmp_wfn(ix,iy)/ss
       end do
    end do

  end do

  Call zheev('V', 'U', NLanczos, Hamiltonian_L, NLanczos, w, work_lp, lwork, rwork, info)

  Exp_Ham_L= 0d0
  do j = 1,NLanczos
    Exp_Ham_L(j,j) = exp(-zI*dt*w(j))
  end do

  zvec = 0d0; zvec(1) = 1d0
  do j = 1,NLanczos
    zvec_t(j) = sum(conjg(Hamiltonian_L(:,j))*zvec(:))
  end do
  do j = 1,NLanczos
    zvec(j) = sum(exp_Ham_L(j,:)*zvec_t(:))
  end do
  do j = 1,NLanczos
    zvec_t(j) = sum(Hamiltonian_L(j,:)*zvec(:))
  end do
  write(*,*)'norm =',sum(abs(zvec_t)**2)

!$omp parallel do private(ix,iy) 
  do ix = 0,Nx
     do iy = 0,Nx
        zwfn(ix,iy) = zvec_t(1)*zwfn_Lanczos(ix,iy,1)
     end do
  end do

  do j = 1,NLanczos
!$omp parallel do private(ix,iy) 
    do iy = 0,Nx
      do ix = 0,Nx
           zwfn(ix,iy) = zwfn(ix,iy) + zvec_t(j)*zwfn_Lanczos(ix,iy,j)
        end do
     end do
  end do

  return
end subroutine dt_evolve_Lanczos
