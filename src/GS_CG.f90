!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine GS_CG
  use global_variables
  use hpsi
  implicit none
  real(dp),allocatable :: xvec(:,:),pvec(:,:),rvec(:,:)
  real(dp),allocatable :: hxvec(:,:),gvec(:,:),hpvec(:,:)

  real(dp) :: xx,pp,xp,xhx,php,xhp,esp,esp_res,gg,gg0
  real(dp) :: ss,lambda,alpha,beta,aa,bb,cc
  integer :: iorb,iorb_t,ix,iy,iter_cg
  real(dp) :: esp_iter(Ncg),esp_res_iter(Ncg)

  allocate( xvec(0:Nx,0:Nx),pvec(0:Nx,0:Nx),rvec(0:Nx,0:Nx) )
  allocate( hxvec(0:Nx,0:Nx),gvec(0:Nx,0:Nx),hpvec(0:Nx,0:Nx) )

  write(*,'(A)')'===== Ground state calculation ===================================================='
  write(*,'(A)')
  write(*,'(A)')

  xvec(:,:)=wfn(:,:)

  call rhpsi(xvec,hxvec)
  
  xx=sum(xvec(:,:)**2)*dx**2
  xhx=sum(xvec(:,:)*hxvec(:,:))*dx**2
  lambda=xhx/xx
  do iter_cg=1,Ncg
     gvec(:,:)=2d0*(hxvec(:,:)-lambda*xvec(:,:))/xx
     
     gg0=sum(gvec(:,:)**2)*dx**2
     select case(iter_cg)
     case(1)
        pvec(:,:)=-gvec(:,:)
     case default
        beta=gg0/gg
        pvec(:,:)=-gvec(:,:)+beta*pvec(:,:)
     end select
     gg=gg0

     call rhpsi(pvec,hpvec)
        
     pp=sum(pvec(:,:)**2)*dx**2
     php=sum(pvec(:,:)*hpvec(:,:))*dx**2
     xp=sum(xvec(:,:)*pvec(:,:))*dx**2
     xhp=sum(hxvec(:,:)*pvec(:,:))*dx**2

     aa=php*xp-xhp*pp
     bb=php*xx-xhx*pp
     cc=xhp*xx-xhx*xp
     ss=bb**2-4d0*aa*cc
     if(ss > 0d0)then
        alpha=(-bb+sqrt(ss))/(2d0*aa)
     else
        exit
     end if
     
     xvec(:,:)=xvec(:,:)+alpha*pvec(:,:)

     call rhpsi(xvec,hxvec)
     xx=sum(xvec(:,:)**2)*dx**2
     xhx=sum(xvec(:,:)*hxvec(:,:))*dx**2
     lambda=xhx/xx
     esp_iter(iter_cg)=lambda
     esp_res_iter(iter_cg)=sum((hxvec(:,:)-lambda*xvec(:,:))**2)*dx**2
     
  end do

  xvec(:,:)=xvec(:,:)/sqrt(xx)
  call rhpsi(xvec,hxvec)
  esp=sum(xvec(:,:)*hxvec(:,:))*dx**2
  esp_res=sum((hxvec(:,:)-esp*xvec(:,:))**2)*dx**2
  wfn(:,:)=xvec(:,:)
  if(wfn(Nx/2,Nx/2) < 0d0) wfn(:,:) = -wfn(:,:)


  write(*,'(A)')'esp,     esp_res'
  write(*,'(e16.6e3,3x,e16.6e3)')esp,esp_res
  write(*,'(A)')'esp + ion-ion, r_dist'
  write(*,'(e16.6e3,3x,e16.6e3)')esp+1d0/sqrt(sigma0**2+r_dist**2),r_dist

  ss = 0d0
  do iy =0,Nx
  do ix =0,Nx
     ss = ss + abs(wfn(ix,iy)-wfn(iy,ix))
  end do
  end do
  ss = ss *dx**2
  write(*,'(A,2x,e26.16e3)')'Symmerty error =',ss


  open(90,file = 'GS_data.out')
  write(90,'(A)')'# iter,  esp(iter),  esp_res(iter)'
  do iter_cg =1,Ncg
     write(90,'(I6,2x,e26.16e3,2x,e26.16e3)')iter_cg,esp_iter(iter_cg),esp_res_iter(iter_cg)
  end do
  close(90)

  open(90,file = 'GS_wfn.out')
  write(90,'(A)')'# x, y, wfn(x,y)'
  do ix =1,Nx
  do iy =1,Nx
     write(90,'(100e26.16e3)')xn(ix),xn(iy),wfn(ix,iy)
  end do
  end do
  close(90)


  open(90,file = 'GS_wfn_matrix.out')
  do iy =0,Nx
     write(90,'(9999e26.16e3)')wfn(:,iy)
  end do
  close(90)

  call wfn_rho
  open(90,file = 'GS_rho.out')
  write(90,'(A)')'# x,  rho(x)'
  do ix =0,Nx
     write(90,'(100e26.16e3)')xn(ix),rho(ix)
  end do
  close(90)
  write(*,*)'rho',sum(rho(:))*dx


  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== End Ground state calculation ================================================'  

  return
end subroutine GS_CG
