!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine zwfn_jx
  use global_variables
  implicit none
  integer :: ix,iy
  real(dp) :: tmp
  complex(zp) :: ztmp

  jx(:)=0d0

!  do ix = 0,Nx
!    tmp = 0d0
!    do iy = 0,Nx
!      tmp = tmp + wfn(ix,iy)**2
!    end do
!    rho(ix) = tmp
!  end do

  ztmp_wfn_b(:,:)=0d0
  ztmp_wfn_b(0:Nx,0:Nx) = zwfn(0:Nx,0:Nx)

  do ix = 0,Nx
    ztmp = 0d0
    do iy = 0,Nx
      ztmp = ztmp + conjg(zwfn(iy,ix)) * ( &
        + gN1*(ztmp_wfn_b(iy,ix+1) - ztmp_wfn_b(iy,ix-1)) &
        + gN2*(ztmp_wfn_b(iy,ix+2) - ztmp_wfn_b(iy,ix-2)) &
        + gN3*(ztmp_wfn_b(iy,ix+3) - ztmp_wfn_b(iy,ix-3)) &
        + gN4*(ztmp_wfn_b(iy,ix+4) - ztmp_wfn_b(iy,ix-4)) ) 
    end do
    jx(ix) = real(2d0 * ztmp /zI) !*dx/dx
  end do

  return
end subroutine zwfn_jx
