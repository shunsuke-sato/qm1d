!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine preparation_RT_collision
  use global_variables
  implicit none
  integer :: ix,iy
  real(dp) :: k_mom, tmp
  complex(zp) :: ztmp1,ztmp2
  complex(zp) :: zwfn_prj(0:Nx)

  k_mom = sqrt(2d0*Ekin0)
  write(*,*)'k_mom =',k_mom

  allocate(zwfn(0:Nx,0:Nx),  zwfn_Lanczos(0:Nx,0:Nx,NLanczos))      
  allocate(field_t(0:Nt_iter+1))
  allocate(dipole_t(0:Nt_iter),norm_t(0:Nt_iter))

  do ix = 0,Nx
    zwfn_prj(ix) = exp(zI*k_mom*xn(ix))*exp(-0.5d0*((xn(ix)-x0_col)/delta_x_col)**2)
  end do
  ztmp1 = sum(zwfn_prj(:)*wfn_sp(:))*dx
  zwfn_prj = zwfn_prj - ztmp1 * wfn_sp

  tmp = sum(abs(zwfn_prj(:))**2)*dx
  zwfn_prj = zwfn_prj /sqrt(tmp)


  do ix = 0,Nx
  do iy = 0,Nx
    ztmp1 = zwfn_prj(ix)*wfn_sp(iy)
    ztmp2 = zwfn_prj(iy)*wfn_sp(ix)
!    ztmp1 = exp(zI*k_mom*xn(ix))*exp(-0.5d0*((xn(ix)-x0_col)/delta_x_col)**2)*wfn_sp(iy)
!    ztmp2 = exp(zI*k_mom*xn(iy))*exp(-0.5d0*((xn(iy)-x0_col)/delta_x_col)**2)*wfn_sp(ix)
    zwfn(ix,iy) = ztmp1 + ztmp2
  end do
  end do

! Normalize
  tmp = sum(abs(zwfn(:,:))**2)*dx**2
  zwfn(:,:) = zwfn(:,:)/sqrt(tmp)

  call zwfn_rho
  call zwfn_jx

  open(21,file='initial_rho.out')
  do ix = 0,Nx
    write(21,'(999e26.16e3)')xn(ix),rho(ix),jx(ix)
  end do
  close(21)

  return
end subroutine preparation_RT_collision
