!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine preparation_GS
  use global_variables
  implicit none
  integer :: ix,iy
  real(dp) :: tmp,r

  write(*,'(A)')'===== preparatin_GS =============================================================='
  write(*,'(A)')
  write(*,'(A)')

  allocate(wfn(0:Nx,0:Nx), v_ext(0:Nx), v_int(0:Nx,0:Nx), v_all(0:Nx,0:Nx))
  allocate(rho(0:Nx), jx(0:Nx))

  write(*,'(A)')'=== preparing initial wave-function ===='
! wfn(0,:) == wfn(Nx,:) == 0 
  wfn(:,:) = 0d0

  do iy = 1,Nx-1  
  do ix = 1,Nx-1
     call random_number(tmp)
     tmp = tmp-0.5d0
     wfn(ix,iy) = tmp
  end do
  end do

  tmp_wfn(:,:) = wfn(:,:)

! Symmetrize
  do iy = 0,Nx
  do ix = 0,Nx
     wfn(ix,iy)=0.5d0*(tmp_wfn(ix,iy)+tmp_wfn(iy,ix))
  end do
  end do

! Normalize
  tmp = sum(wfn(:,:)**2)*dx**2
  wfn(:,:)=wfn(:,:)/sqrt(tmp)


  write(*,'(A)')'=== preparing external potential ===='
  do ix=0,Nx
!     v_ext(ix) = 0.5d0*xn(ix)**2
!     v_ext(ix) = -2d0/sqrt(1d0+xn(ix)**2)
     v_ext(ix) = -1d0/sqrt(1d0+xn(ix)**2)
  end do
!  v_KS_tdehf_old = v_ext
!  v_KS_tdehf = v_ext

  write(*,'(A)')'=== preparing interaction potential ===='
  do iy=0,Nx
  do ix=0,Nx
     r = xn(ix) - xn(iy)
     v_int(ix,iy) = 0d0 !1d0/sqrt(1d0+r**2) ! debut
!     v_int(ix,iy) = 1d0/sqrt(1d0+r**2)
  end do
  end do

  write(*,'(A)')'=== preparing total potential ===='
  do iy=0,Nx
  do ix=0,Nx
     v_all(ix,iy) = v_ext(ix) + v_ext(iy) + v_int(ix,iy)
  end do
  end do
  
  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== Complete preparatin_GS ====================================================='

  return
end subroutine preparation_GS
