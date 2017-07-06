!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine single_particle_GS
  use global_variables
  use timer
  implicit none
  integer,parameter :: Ncg_t = 500
! CG
  real(dp),allocatable :: xvec(:),pvec(:),rvec(:)
  real(dp),allocatable :: hxvec(:),gvec(:),hpvec(:)

  real(dp) :: xx,pp,xp,xhx,php,xhp,esp,esp_res,gg,gg0
  real(dp) :: ss,lambda,alpha,beta,aa,bb,cc

  integer :: ix, icg
  real(dp) :: rnd

  real(8) :: time_ini,time_end,time_elapse

  call elapse_time(time_ini)

  allocate (wfn_sp(0:Nx))
  allocate( xvec(0:Nx),pvec(0:Nx),rvec(0:Nx) )
  allocate( hxvec(0:Nx),gvec(0:Nx),hpvec(0:Nx) )


  do ix = 0,Nx
    call random_number(rnd)
    wfn_sp(ix) = rnd
  end do

  ss = sum(wfn_sp(:)**2)*dx
  wfn_sp(:) = wfn_sp(:)/sqrt(ss)

  xvec(:) = wfn_sp(:)
  call hpsi_sp(xvec,hxvec,v_ext)
  xx=sum(xvec(:)**2)*dx
  xhx=sum(xvec(:)*hxvec(:))*dx
  lambda=xhx/xx

  do icg = 1, Ncg_t
    gvec(:)=2d0*(hxvec(:)-lambda*xvec(:))/xx
     
    gg0=sum(gvec(:)**2)*dx
    select case(icg)
    case(1)
      pvec(:)=-gvec(:)
    case default
      beta=gg0/gg
      pvec(:)=-gvec(:)+beta*pvec(:)
    end select
    gg=gg0
    
    call hpsi_sp(pvec,hpvec,v_ext)
    
    pp=sum(pvec(:)**2)*dx
    php=sum(pvec(:)*hpvec(:))*dx
    xp=sum(xvec(:)*pvec(:))*dx
    xhp=sum(hxvec(:)*pvec(:))*dx
    
    aa=php*xp-xhp*pp
    bb=php*xx-xhx*pp
    cc=xhp*xx-xhx*xp
    ss=bb**2-4d0*aa*cc
    if(ss > 0d0)then
      alpha=(-bb+sqrt(ss))/(2d0*aa)
    else
      exit
    end if
    
    xvec(:)=xvec(:)+alpha*pvec(:)
    
    call hpsi_sp(xvec,hxvec,v_ext)
    xx=sum(xvec(:)**2)*dx
    xhx=sum(xvec(:)*hxvec(:))*dx
    lambda=xhx/xx

    write(*,'(I7,2x,3e26.16e3)')icg-1,lambda,sum((hxvec(:)-lambda*xvec(:))**2)*dx
  end do

  xx=sum(xvec(:)**2)*dx
  wfn_sp(:) = xvec(:)/sqrt(xx)

  call elapse_time(time_end)
  time_elapse = time_end - time_ini
  write(*,"(A,2x,e16.6e3,A)")'Elapse time for single_particle_GS',time_elapse,"[sec]"

end subroutine single_particle_GS
