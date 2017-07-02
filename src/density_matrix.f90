!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module density_matrix
  use global_variables
  implicit none

  private
! temporary
  complex(zp),allocatable :: zdensity_matrix(:,:)

  public :: initialize_density_matrix &
           ,calc_occ_natural_orbital

  
contains
!-----------------------------------------------------------------------------------------
  subroutine initialize_density_matrix
    implicit none

    allocate(zdensity_matrix(0:Nx,0:Nx))

  end subroutine initialize_density_matrix
!-----------------------------------------------------------------------------------------
  subroutine calc_occ_natural_orbital(occ_nat)
    implicit none
    real(8),intent(out) :: occ_nat(Nx+1)
    integer :: ix,iy,iz
    complex(zp) :: zs
!LAPACK
    integer :: lwork
    integer :: Nmat
    complex(8),allocatable :: work_lp(:)
    real(8),allocatable :: rwork(:),w(:)
    integer :: info

    lwork = 60*(Nx+1)
    Nmat = Nx+1
    allocate(work_lp(lwork),rwork(3*Nmat-2),w(Nmat))

!$omp parallel do private(ix,iy,iz,zs)
    do ix = 0,Nx
      do iy = 0,Nx
        zs = 0d0
        do iz = 0,Nx
          zs = zs + zwfn(ix,iz)*conjg(zwfn(iy,iz))
        end do
        zs = zs * dx**2
        zdensity_matrix(ix,iy) = zs
      end do
    end do

    call zheev('N', 'U', Nx+1, zdensity_matrix(0:Nx,0:Nx), Nx+1, w, work_lp, lwork, rwork, info)

    do ix = 1,Nx+1
      occ_nat(ix) = w(Nx+2-ix)
    end do

  end subroutine calc_occ_natural_orbital
!-----------------------------------------------------------------------------------------
end module density_matrix
