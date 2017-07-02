!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine preparation_RT
  use global_variables
  implicit none
  integer :: ix,iy,it
  real(dp) :: tt

  write(*,'(A)')'===== preparatin_RT =============================================================='
  write(*,'(A)')
  write(*,'(A)')
  allocate(zwfn(0:Nx,0:Nx))
  allocate(field_t(0:Nt_iter+1))
  allocate(dipole_t(0:Nt_iter),norm_t(0:Nt_iter))
  select case(RT_mode)
     case('kick')
        do iy=0,Nx
        do ix=0,Nx
           zwfn(ix,iy)=wfn(ix,iy)*exp(-zI*(xn(ix)+xn(iy))*kick_mom)
        end do
        end do
     case('gs')
        zwfn(:,:)=wfn(:,:)
     case default 
        stop 'Error RT_mode'
     end select


  field_max = field_max_eV_per_AA * (0.529d0/27.2d0)
  field_duration = field_duration_fs/0.02418d0
  field_omega = field_omega_eV/27.2d0

     do it = 0,Nt_iter
        tt = dt*dble(it)
        if(tt < field_duration)then
           field_t(it) = field_max*sin(field_omega*(tt-0.5d0*field_duration)) &
                &*cos(pi*(tt-0.5d0*field_duration)/field_duration)**2
        else
           field_t(it) = 0d0
        end if
     end do
     if(RT_mode == 'kick')field_t(:) = 0d0

  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== complete preparatin_RT ====================================================='

        
  return
end subroutine preparation_RT
