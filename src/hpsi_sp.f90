!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine hpsi_sp(u,hu,v)
  use global_variables
  implicit none
  real(dp) :: u(0:Nx),hu(0:Nx),tu(0-4:Nx+4),v(0:Nx)
! finite difference
  real(dp) :: c0,c1,c2,c3,c4
  integer :: ix

  c0=-0.5d0*cN0/(dx**2)
  c1=-0.5d0*cN1/(dx**2)
  c2=-0.5d0*cN2/(dx**2)
  c3=-0.5d0*cN3/(dx**2)
  c4=-0.5d0*cN4/(dx**2)

  tu = 0d0
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
