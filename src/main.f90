!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
program main
use global_variables
  implicit none

  write(*,'(A)')'Start qm1d'
  call input
  call mesh

  call preparation_GS
!  call GS_CG

  call single_particle_GS
!  call IMA_TIME_PROP
  write(*,*)'End single_particle_GS'

  call preparation_RT_collision
!  stop
!  call preparation_RT
  call RT_prop

  call write_results

end program main

