!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine wfn_rho
  use global_variables
  implicit none
  integer :: ix,iy
  real(dp) :: tmp

  rho(:)=0d0

!  do ix = 0,Nx
!    tmp = 0d0
!    do iy = 0,Nx
!      tmp = tmp + wfn(ix,iy)**2
!    end do
!    rho(ix) = tmp
!  end do

  do ix = 0,Nx
    rho(ix) = sum(wfn(:,ix)**2)
  end do
  rho = rho*2d0*dx

  return
end subroutine wfn_rho
!=======10========20========30========40========50========60========70========80========90=======100
subroutine zwfn_rho
  use global_variables
  implicit none
  integer :: ix,iy
  real(dp) :: tmp

  rho(:)=0d0

!  do ix = 0,Nx
!    tmp = 0d0
!    do iy = 0,Nx
!      tmp = tmp + wfn(ix,iy)**2
!    end do
!    rho(ix) = tmp
!  end do

  do ix = 0,Nx
    rho(ix) = sum(abs(zwfn(:,ix))**2)
  end do
  rho = rho*2d0*dx

  return
end subroutine zwfn_rho
