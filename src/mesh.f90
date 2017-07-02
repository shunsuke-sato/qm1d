!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine mesh
  use global_variables
  use hpsi
  implicit none
  integer :: ix, iy 

  write(*,'(A)')'===== Making mesh ================================================================'
  write(*,'(A)')
  write(*,'(A)')

  allocate(xn(0:Nx),xyn(0:Nx,0:Nx))
  dx = length_x/dble(Nx)
  
  do ix = 0,Nx
     xn(ix) = dx*dble(ix) - 0.5d0*length_x
  end do
  
  do ix = 0,Nx
  do iy = 0,Nx
     xyn(ix,iy)= xn(ix)+xn(iy)
  end do
  end do

  call initialize_hpsi

  write(*,'(A)')'===== Complete Making mesh ========================================================'
  return
end subroutine mesh
