!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine RT_prop
  use global_variables
  implicit none
  real(dp) :: ft,dipole,norm
  integer :: it,ix
  character(50) :: citer

  write(*,'(A)')'===== Real time propagation ======================================================'
  write(*,'(A)')
  write(*,'(A)')

  call dipole_norm(dipole,norm)
  dipole_t(0) = dipole
  norm_t(0) = norm


  it = 0
  call zwfn_rho
  call zwfn_jx

  write(citer,'(I7.7)')it

  do it = 1, Nt_iter
!     if(mod(it, 1000) == 0)write(*,"(A,2x,I0)")"iter =",it
     write(*,"(A,2x,I0)")"iter =",it
!     ft = 0.5d0*(field_t(it)+field_t(it-1))
     ft = 0d0
!     call dt_evolve(ft)
     call dt_evolve_Lanczos(ft)
     call dipole_norm(dipole,norm)
     dipole_t(it) = dipole
     norm_t(it) = norm

     call zwfn_rho
     call zwfn_jx
       

     if (mod(it,10) == 0 .or. it == 1)then
!     if (1 == 0)then

       write(citer,'(I7.7)')it
       open(22,file='rho_'//trim(citer)//'.out')
       do ix = 0,Nx
         write(22,'(999e26.16e3)')xn(ix),rho(ix),jx(ix)
       end do
       close(22)


       call zwfn_rho
       call zwfn_jx

     end if
  end do

  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== complete Real time propagation ============================================'

  
  return
end subroutine RT_prop
