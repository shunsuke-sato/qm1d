!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module hpsi
  use global_variables
  implicit none

  private
! temporary
  real(dp),allocatable :: tmp_wfn_b(:,:)
  complex(zp),allocatable :: ztmp_wfn_b(:,:)

  public :: initialize_hpsi &
           ,rhpsi &
           ,zhpsi
  
contains
!-----------------------------------------------------------------------------------------
  subroutine initialize_hpsi
    implicit none

    allocate(tmp_wfn_b(-4:Nx+4,-4:Nx+4))
    allocate(ztmp_wfn_b(-4:Nx+4,-4:Nx+4))
    tmp_wfn_b = 0d0
    ztmp_wfn_b = 0d0

  end subroutine initialize_hpsi
!-----------------------------------------------------------------------------------------
  subroutine rhpsi(tmp_wfn,tmp_hwfn)
    implicit none
    real(dp),intent(in) :: tmp_wfn(0:Nx,0:Nx)
    real(dp),intent(out) :: tmp_hwfn(0:Nx,0:Nx)
! finite difference
    integer :: ix,iy
    real(dp) :: c0,c1,c2,c3,c4
! nine-points formula  
    c0=-0.5d0*cN0/(dx**2)
    c1=-0.5d0*cN1/(dx**2)
    c2=-0.5d0*cN2/(dx**2)
    c3=-0.5d0*cN3/(dx**2)
    c4=-0.5d0*cN4/(dx**2)

    tmp_wfn_b(0:Nx,0:Nx) = tmp_wfn(0:Nx,0:Nx)
    do iy=0,Nx
      do ix=0,Nx
        tmp_hwfn(ix,iy)=2d0*c0*tmp_wfn_b(ix,iy) &
          + c1*(tmp_wfn_b(ix+1,iy) + tmp_wfn_b(ix-1,iy) + tmp_wfn_b(ix,iy+1) + tmp_wfn_b(ix,iy-1)) &
          + c2*(tmp_wfn_b(ix+2,iy) + tmp_wfn_b(ix-2,iy) + tmp_wfn_b(ix,iy+2) + tmp_wfn_b(ix,iy-2)) &
          + c3*(tmp_wfn_b(ix+3,iy) + tmp_wfn_b(ix-3,iy) + tmp_wfn_b(ix,iy+3) + tmp_wfn_b(ix,iy-3)) &
          + c4*(tmp_wfn_b(ix+4,iy) + tmp_wfn_b(ix-4,iy) + tmp_wfn_b(ix,iy+4) + tmp_wfn_b(ix,iy-4)) 
      end do
    end do
    
    tmp_hwfn(:,:) = tmp_hwfn(:,:) + v_all(:,:)*tmp_wfn(:,:)
    return
  end subroutine rhpsi
  !-----------------------------------------------------------------------------------------
  subroutine zhpsi(ztmp_wfn,ztmp_hwfn,ft)
    implicit none
    complex(zp),intent(in) :: ztmp_wfn(0:Nx,0:Nx)
    complex(zp),intent(out) :: ztmp_hwfn(0:Nx,0:Nx)
    real(dp),intent(in) :: ft
    integer :: ix,iy
    real(dp) :: c0,c1,c2,c3,c4

! nine-points formula  
    c0=-0.5d0*cN0/(dx**2)
    c1=-0.5d0*cN1/(dx**2)
    c2=-0.5d0*cN2/(dx**2)
    c3=-0.5d0*cN3/(dx**2)
    c4=-0.5d0*cN4/(dx**2)

    ztmp_wfn_b(0:Nx,0:Nx) = ztmp_wfn(0:Nx,0:Nx)
    do iy=0,Nx
      do ix=0,Nx
        ztmp_hwfn(ix,iy)=2d0*c0*ztmp_wfn_b(ix,iy) &
          + c1*(ztmp_wfn_b(ix+1,iy) + ztmp_wfn_b(ix-1,iy) + ztmp_wfn_b(ix,iy+1) + ztmp_wfn_b(ix,iy-1)) &
          + c2*(ztmp_wfn_b(ix+2,iy) + ztmp_wfn_b(ix-2,iy) + ztmp_wfn_b(ix,iy+2) + ztmp_wfn_b(ix,iy-2)) &
          + c3*(ztmp_wfn_b(ix+3,iy) + ztmp_wfn_b(ix-3,iy) + ztmp_wfn_b(ix,iy+3) + ztmp_wfn_b(ix,iy-3)) &
          + c4*(ztmp_wfn_b(ix+4,iy) + ztmp_wfn_b(ix-4,iy) + ztmp_wfn_b(ix,iy+4) + ztmp_wfn_b(ix,iy-4)) 
      end do
    end do
    
    ztmp_hwfn(:,:) = ztmp_hwfn(:,:) + (v_all(:,:)+ft*xyn(:,:))*ztmp_wfn(:,:)
    
    return
  end subroutine zhpsi
!-----------------------------------------------------------------------------------------

end module hpsi
