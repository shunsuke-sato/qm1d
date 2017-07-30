!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine RT_prop_Born_exp
  use global_variables
  use timer
  use density_matrix
  implicit none
  real(dp) :: ft,dipole,norm
  integer :: it,ix
  character(50) :: citer
  real(dp) :: time_ini,time_end,time_elapse
  real(dp) :: occ_nat(Nx+1)
  complex(zp) :: zwfn_0th(0:Nx,0:Nx)
  complex(zp) :: zwfn_1st(0:Nx,0:Nx)

  call elapse_time(time_ini)

  write(*,'(A)')'===== Real time propagation (Born exp.) ======================================================'
  write(*,'(A)')
  write(*,'(A)')


  it = 0
  zwfn_0th = zwfn
  zwfn_1st = 0d0
  zwfn = zwfn_0th + zwfn_1st

  call zwfn_rho
  call zwfn_jx

  write(citer,'(I7.7)')it
  open(22,file='rho_'//trim(citer)//'.out')
  do ix = 0,Nx
    write(22,'(999e26.16e3)')xn(ix),rho(ix),jx(ix)
  end do
  close(22)

  open(45,file='occ_nat.out')
  call calc_occ_natural_orbital(occ_nat)
  write(45,"(999e26.16e3)")it*dt,occ_nat(1:100)

  write(citer,'(I7.7)')it

  do it = 1, Nt_iter

     write(*,"(A,2x,I0)")"iter =",it

     ft = 0d0

     zwfn_1st = zwfn_1st -0.5d0*zI*dt*v_all*zwfn_0th
     call dt_evolve_Lanczos_Born_exp(zwfn_1st,ft)
     call dt_evolve_Lanczos_Born_exp(zwfn_0th,ft)
     zwfn_1st = zwfn_1st -0.5d0*zI*dt*v_all*zwfn_0th

     if (mod(it,1000) == 0)then
       zwfn = zwfn_0th + zwfn_1st        
       call calc_occ_natural_orbital(occ_nat)
       write(45,"(999e26.16e3)")it*dt,occ_nat(1:100)
     end if

     if (mod(it,100) == 0 .or. it == 1)then
!     if (1 == 0)then
       zwfn = zwfn_0th + zwfn_1st        
       call zwfn_rho
       call zwfn_jx

       write(citer,'(I7.7)')it
       open(22,file='rho_'//trim(citer)//'.out')
       do ix = 0,Nx
         write(22,'(999e26.16e3)')xn(ix),rho(ix),jx(ix)
       end do
       close(22)

     end if
  end do

  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== complete Real time propagation ============================================'

  close(45)

  call elapse_time(time_end)
  time_elapse = time_end - time_ini
  write(*,"(A,2x,e16.6e3,A)")'Elapse time for RT_prop',time_elapse,"[sec]"
  

end subroutine RT_prop_Born_exp
!----------------------------------------------------------------------------------
!=======10========20========30========40========50========60========70========80========90=======100
subroutine dt_evolve_Lanczos_Born_exp(zwfn_t)
  use global_variables
  use hpsi
  implicit none
  complex(zp),intent(inout) :: zwfn_t(0:Nx,0:Nx)
  real(dp) :: ft
  integer :: iexp,j
  complex(zp) :: zs
  complex(zp) :: Hamiltonian_L(NLanczos,NLanczos), Exp_Ham_L(NLanczos,NLanczos), Umat_L(NLanczos,NLanczos)
  complex(zp) :: zvec(NLanczos),zvec_t(NLanczos)

  real(dp) :: ss
!LAPACK
  integer :: lwork
  complex(zp),allocatable :: work_lp(:)
  real(dp),allocatable :: rwork(:),w(:)
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
        ss = ss + abs(zwfn_t(ix,iy))**2
     end do
  end do
  ss = dx**2*ss
  ss = 1d0/sqrt(ss)

!$omp parallel do private(ix,iy)
  do iy = 0,Nx
     do ix = 0,Nx
        zwfn_t(ix,iy) = zwfn_t(ix,iy) * ss
        zwfn_Lanczos(ix,iy,1) = zwfn_t(ix,iy)
     end do
  end do

  Hamiltonian_L = 0d0

  do j = 1,NLanczos

    call zhpsi0(zwfn_Lanczos(:,:,j),ztmp_hwfn,ft)
    ss = 0d0
!$omp parallel do private(ix,iy) reduction(+:ss)
    do iy = 0,Nx
      do ix = 0,Nx
!          ss = ss + conjg(ztmp_wfn(ix,iy))*ztmp_hwfn(ix,iy)
        ss = ss + conjg(zwfn_Lanczos(ix,iy,j))*ztmp_hwfn(ix,iy)
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

    ss = 1d0/ss
!$omp parallel do private(ix,iy) 
    do iy = 0,Nx
       do ix = 0,Nx
          zwfn_Lanczos(ix,iy,j+1) = ztmp_wfn(ix,iy)*ss
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


  call Lanczos_sum


contains
  subroutine Lanczos_sum
    implicit none

    select case(NLanczos)
    case(3)
!$omp parallel do private(ix,iy) 
      do iy = 0,Nx
        do ix = 0,Nx
          zwfn_t(ix,iy) = zvec_t(1)*zwfn_Lanczos(ix,iy,1) &
                       +zvec_t(2)*zwfn_Lanczos(ix,iy,2) &
                       +zvec_t(3)*zwfn_Lanczos(ix,iy,3)
        end do
      end do
    case(4)
!$omp parallel do private(ix,iy) 
      do iy = 0,Nx
        do ix = 0,Nx
          zwfn_t(ix,iy) = zvec_t(1)*zwfn_Lanczos(ix,iy,1) &
                       +zvec_t(2)*zwfn_Lanczos(ix,iy,2) &
                       +zvec_t(3)*zwfn_Lanczos(ix,iy,3) &
                       +zvec_t(4)*zwfn_Lanczos(ix,iy,4)
        end do
      end do
    case(5)
!$omp parallel do private(ix,iy) 
      do iy = 0,Nx
        do ix = 0,Nx
          zwfn_t(ix,iy) = zvec_t(1)*zwfn_Lanczos(ix,iy,1) &
                       +zvec_t(2)*zwfn_Lanczos(ix,iy,2) &
                       +zvec_t(3)*zwfn_Lanczos(ix,iy,3) &
                       +zvec_t(4)*zwfn_Lanczos(ix,iy,4) &
                       +zvec_t(5)*zwfn_Lanczos(ix,iy,5)
        end do
      end do
    case(6)
!$omp parallel do private(ix,iy) 
      do iy = 0,Nx
        do ix = 0,Nx
          zwfn_t(ix,iy) = zvec_t(1)*zwfn_Lanczos(ix,iy,1) &
                       +zvec_t(2)*zwfn_Lanczos(ix,iy,2) &
                       +zvec_t(3)*zwfn_Lanczos(ix,iy,3) &
                       +zvec_t(4)*zwfn_Lanczos(ix,iy,4) &
                       +zvec_t(5)*zwfn_Lanczos(ix,iy,5) &
                       +zvec_t(6)*zwfn_Lanczos(ix,iy,6)
        end do
      end do
    case default
!$omp parallel do private(ix,iy) 
      do iy = 0,Nx
        do ix = 0,Nx
          zwfn_t(ix,iy) = zvec_t(1)*zwfn_Lanczos(ix,iy,1)
        end do
      end do

      do j = 2,NLanczos
!$omp parallel do private(ix,iy) 
        do iy = 0,Nx
          do ix = 0,Nx
            zwfn_t(ix,iy) = zwfn_t(ix,iy) + zvec_t(j)*zwfn_Lanczos(ix,iy,j)
          end do
        end do
      end do
    end select
  end subroutine Lanczos_sum

end subroutine dt_evolve_Lanczos_Born_exp
