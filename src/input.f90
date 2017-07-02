!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine input
  use global_variables
  implicit none

! == input parameter == !                                                                                                                       
  write(*,'(A)')'===== Input parameter ============================================================='
  write(*,'(A)')

  Nx = 2000
  length_x = 400d0
! GS
  Ncg = 800
! TD
  T_calc = 400d0
  dt = 0.01d0 ! 0.01d0
  Nt_iter = aint(T_calc/dt)+1

  RT_mode='kick' ! kick or gs
  kick_mom = 1d-3

! collision
  Ekin0 = 0.6d0
  x0_col = -100d0
  delta_x_col = 15d0

  field_max_eV_per_AA = 1d0
  field_duration_fs = 10d0
  field_omega_eV = 1.55d0

  field_max = field_max_eV_per_AA * (0.529d0/27.2d0)
  field_duration = field_duration_fs/0.02418d0
  field_omega = field_omega_eV/27.2d0

  write(*,'(A,2x,I4)')'Nx =',Nx
  write(*,'(A,2x,e26.16e3)')'length_x =',length_x
  write(*,'(A,2x,e26.16e3)')'T_calc =',T_calc
  write(*,'(A,2x,e26.16e3)')'dt =',dt
  write(*,'(A,2x,I10)')'Nt_iter =',Nt_iter


! temporary array
  allocate(tmp_wfn(0:Nx,0:Nx),tmp_hwfn(0:Nx,0:Nx),tmp_wfn_b(-4:Nx+4,-4:Nx+4))
  allocate(ztmp_wfn(0:Nx,0:Nx),ztmp_hwfn(0:Nx,0:Nx),ztmp_wfn_b(-4:Nx+4,-4:Nx+4))

  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== Complete Input parameter ===================================================' 
  return
end subroutine input
