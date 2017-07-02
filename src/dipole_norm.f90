!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine dipole_norm(dipole,norm)
  use global_variables
  implicit none
  real(dp) :: dipole,norm

  dipole = sum(xyn(:,:)*abs(zwfn(:,:))**2)*dx**2
  norm = sum(abs(zwfn(:,:))**2)*dx**2
  
  return
end subroutine dipole_norm
