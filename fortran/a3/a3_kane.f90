program kadai3
    implicit none

    double precision,parameter::h=10.0d0
    double precision,parameter::a=1.2d0
    double precision,parameter::b=(a+1.0d0)/(a-1.0d0)
    integer,parameter::xm=30
    integer,parameter::ym=30
    double precision,parameter::r=2.871d2       !ガス定数
    double precision,parameter::gamma=1.4d0     !比熱比
    double precision,parameter::dt=1.0d-6
    double precision,parameter::ii=1.0d0
    double precision,dimension(0:xm,0:ym)::x
    double precision,dimension(0:xm,0:ym)::y
    double precision,dimension(0:xm,0:ym)::gx,gy,ex,ey,xg,yg,xe,ye,jj,uu,vv
    double precision,dimension(0:xm,0:ym,1:4)::qq,ee,ff
    double precision,dimension(0:xm,0:ym)::rho=1.2d0,p=0.0d0,u=0.0d0,v=0.0d0,&
    e=0.0d0,eb,t,c
    integer i,j,n,k
    double precision dx,dy,m,alpha,jm

    do j=0,ym
       do i=0,xm
          dx=(3.0d0*h)/dble(xm)
          x(i,j)=dble(i)*dx
          jm=dble(j)/(dble(ym))
          y(i,j)=h*(b**jm-1.0d0)/(b-1.0d0)
       end do
    end do

    do j=0,ym
       do i=0,xm
          t(i,j)=2.930d2
          c(i,j)=dsqrt(gamma*r*t(i,j)) !音速
       end do
    end do
    m=2.9d0                         !マッハ数

    do j=0,ym
       u(0,j)=m*c(0,j)
       u(1,j)=u(0,j)
    end do

    do j=0,ym
       do i=2,xm
          alpha=(x(i,j)-x(i-2,j))/(x(i-2,j)-x(i-1,j))
          u(i,j)=u(i-2,j)+alpha*(u(i-2,j)-u(i-1,j))
       end do
    end do


    do j=0,ym
       do i=0,xm
    p(i,j)=rho(i,j)*r*t(i,j)
    e(i,j)=rho(i,j)*((r*t(i,j)/(gamma-1.0d0))+((u(i,j)**2.0)+v(i,j)**2.0)/2.0d0)
       end do
    end do

    !座標変換
    do j=1,ym-1
       do i=0,xm
          ye(i,j)=((y(i,j+1)-y(i,j-1))/2.0d0)
          xe(i,j)=((x(i,j+1)-x(i,j-1))/2.0d0)
       end do
    end do
    do j=0,ym
       do i=1,xm-1
          yg(i,j)=((y(i+1,j)-y(i-1,j))/2.0d0)
          xg(i,j)=((x(i+1,j)-x(i-1,j))/2.0d0)
       end do
    end do
    do i=0,xm
       ye(i,0)=(-3.0d0*y(i,0)+4.0d0*y(i,1)-y(i,2))/(2.0d0*ii)
       ye(i,ym)=(3.0d0*y(i,ym)-4.0d0*y(i,ym-1)+y(i,ym-2))/(2.0d0*ii)
       xe(i,0)=(-3.0d0*x(i,0)+4.0d0*x(i,1)-x(i,2))/(2.0d0**ii)
       xe(i,ym)=(3.0d0*x(i,ym)-4.0d0*x(i,ym-1)+x(i,ym-2))/(2.0d0*ii)
    end do
    do j=0,ym
       yg(0,j)=(-3.0d0*y(0,j)+4.0d0*y(1,j)-y(2,j))/(2.0d0*ii)
       yg(xm,j)=(3.0d0*y(xm,j)-4.0d0*y(xm-1,j)+y(xm-2,j))/(2.0d0*ii)
       xg(0,j)=(-3.0d0*x(0,j)+4.0d0*x(1,j)-x(2,j))/(2.0d0*ii)
       xg(xm,j)=(3.0d0*x(xm,j)-4.0d0*x(xm-1,j)+x(xm-2,j))/(2.0d0*ii)
    end do

    do j=0,ym
       do i=0,xm
          jj(i,j)=1.0d0/(xg(i,j)*ye(i,j)-yg(i,j)*xe(i,j))
          gx(i,j)=ye(i,j)*jj(i,j)
          gy(i,j)=-xe(i,j)*jj(i,j)
          ex(i,j)=-yg(i,j)*jj(i,j)
          ey(i,j)=xg(i,j)*jj(i,j)
       end do
    end do


    do n=0,1000           !計算回数
    do j=0,ym
       do i=0,xm
          uu(i,j)=gx(i,j)*u(i,j)+gy(i,j)*v(i,j)
          vv(i,j)=ex(i,j)*u(i,j)+ey(i,j)*v(i,j)
          qq(i,j,1)=rho(i,j)/jj(i,j)
          qq(i,j,2)=(rho(i,j)*u(i,j))/jj(i,j)
          qq(i,j,3)=(rho(i,j)*v(i,j))/jj(i,j)
          qq(i,j,4)=e(i,j)/jj(i,j)
          ee(i,j,1)=(rho(i,j)*uu(i,j))/jj(i,j)
          ee(i,j,2)=(rho(i,j)*u(i,j)*uu(i,j)+gx(i,j)*p(i,j))/jj(i,j)
          ee(i,j,3)=(rho(i,j)*v(i,j)*uu(i,j)+gy(i,j)*p(i,j))/jj(i,j)
          ee(i,j,4)=((e(i,j)+p(i,j))*uu(i,j))/jj(i,j)
          ff(i,j,1)=(rho(i,j)*vv(i,j))/jj(i,j)
          ff(i,j,2)=(rho(i,j)*u(i,j)*vv(i,j)+ex(i,j)*p(i,j))/jj(i,j)
          ff(i,j,3)=(rho(i,j)*v(i,j)*vv(i,j)+ey(i,j)*p(i,j))/jj(i,j)
          ff(i,j,4)=((e(i,j)+p(i,j))*vv(i,j))/jj(i,j)
       end do
    end do

    do j=1,ym-1
       do i=1,xm-1
          do k=1,4
             qq(i,j,k)=qq(i,j,k)+dt*(((ee(i-1,j,k)-ee(i+1,j,k))/2.0d0)+&
                       ((ff(i,j-1,k)-ff(i,j+1,k))/2.0d0))
          end do
       end do
    end do


    do j=1,ym-1
       do i=1,xm-1
          rho(i,j)=jj(i,j)*qq(i,j,1)
          u(i,j)=jj(i,j)*qq(i,j,2)/rho(i,j)
          v(i,j)=jj(i,j)*qq(i,j,3)/rho(i,j)
          e(i,j)=jj(i,j)*qq(i,j,4)
       end do
    end do

    do i=1,xm-1
       alpha=(y(i,ym)-y(i,ym-2))/(y(i,ym-1)-y(i,ym-2))
       rho(i,ym)=rho(i,ym-2)+alpha*(rho(i,ym-1)-rho(i,ym-2))
       u(i,ym)=u(i,ym-2)+alpha*(u(i,ym-1)-u(i,ym-2))
       v(i,ym)=v(i,ym-2)+alpha*(v(i,ym-1)-v(i,ym-2))
       e(i,ym)=e(i,ym-2)+alpha*(e(i,ym-1)-e(i,ym-2))
       alpha=(y(i,0)-y(i,2))/(y(i,1)-y(i,2))
       rho(i,0)=rho(i,2)+alpha*(rho(i,1)-rho(i,2))
       u(i,0)=u(i,2)+alpha*(u(i,1)-u(i,2))
       v(i,0)=v(i,2)+alpha*(v(i,1)-v(i,2))
       e(i,0)=e(i,2)+alpha*(e(i,1)-e(i,2))
    end do
    do j=0,ym
       alpha=(x(xm,j)-x(xm-2,j))/(x(xm-1,j)-x(xm-2,j))
       rho(xm,j)=rho(xm-2,j)+alpha*(rho(xm-1,j)-rho(xm-2,j))
       u(xm,j)=u(xm-2,j)+alpha*(u(xm-1,j)-u(xm-2,j))
       v(xm,j)=v(xm-2,j)+alpha*(v(xm-1,j)-v(xm-2,j))
       e(xm,j)=e(xm-2,j)+alpha*(e(xm-1,j)-e(xm-2,j))
    end do

    do i=0,xm
       do j=0,ym
          eb(i,j)=e(i,j)/rho(i,j)-((u(i,j)**2.0)/2.0d0)
          t(i,j)=eb(i,j)*(gamma-1.0d0)/r
          p(i,j)=(gamma-1.0d0)*rho(i,j)*eb(i,j)
       end do
    end do


    end do

    open(17,file='kadai3.dat',status='replace')
       do j=0,ym
          do i=0,xm
             write(17,'(f13.6,$)') x(i,j)
             write(17,'(f13.6,$)') y(i,j)
             write(17,'(f13.6,$)') u(i,j)
             write(17,'(f13.6)') v(i,j)
          end do
       end do
    close(17)
    end program kadai3
