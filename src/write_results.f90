!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine write_results
  use global_variables
  implicit none
  integer :: it
  real(dp) :: tt

  open(20,file="dipole_norm_t.out")
  write(20,"(A)")"time (a.u.), dipole (a.u.), norm (a.u.)"
  do it = 0,Nt_iter
     tt = dt*dble(it)
     write(20,"(2x,3(e26.16e3,2x))")tt, dipole_t(it), norm_t(it)
  end do
  close(20)
  return
end subroutine write_results
