program assyuku4
implicit none
!ŠiqğŒ
double precision,parameter :: H = 10.0d0
integer,parameter :: IM = 30,JM = 30 !Šiq”
double precision,parameter :: dx = 3.0*H/30.0 !x•ûŒü‚ÌŠiq•
double precision,parameter :: b = 1.2d0 !Œö”äb
double precision,parameter :: bt = 30.0d0/180.0d0*dacos(-1.0d0)
double precision :: a,ix

!—¬‚êê‰ŠúğŒ
double precision,dimension(0:IM,0:JM) :: R = 287.1d0,g = 1.4d0,dt = 1.0d-6 !ƒ}ƒbƒn”,‰·“x,ƒKƒX’è”,”ä”M”ä,–§“x,ŠÔ‚İ

!‚»‚Ì‘¼‚ÌğŒ
integer :: i,j,n
double precision,dimension(0:IM,0:JM) :: x,y,u,v,p,e,T,rho
double precision,dimension(0:IM,0:JM) :: M1,M2,T1,T2,rho1,rho2,p1,p2,e1,e2,u1,u2,v1,v2,vv1,vv2,vn1,vn2,vt,theta

!Šiq¶¬
do i = 0,IM
  do j = 0,JM
    x(i,j) = dx*real(i) !x•ûŒüi“™ŠÔŠuŠiqj
  end do
end do

do j = 0,JM
  do i = 0,IM
    a = (b+1)/(b-1)
    ix = dble(j)/JM
    y(i,j) = H*(a**ix-1)/(a-1) !y•ûŒüi•s“™ŠÔŠuŠiqj
  end do
end do

!”ğŒ
do j = 0,JM
  do i = 0,IM
    M1 = 2.9d0
    rho1 = 1.2d0
    T1 = 293.0d0
    u1(i,j) = M1(i,j)*dsqrt(g(i,j)*R(i,j)*T1(i,j))
    v1(i,j) = 0.0
    vv1(i,j) = u1(i,j)
    vn1(i,j) = u1(i,j)*dsin(bt)
    vt(i,j) = vn1(i,j)/dtan(bt)
    p1(i,j) = rho1(i,j)*R(i,j)*T1(i,j)
    e1(i,j) = rho1(i,j)*((R(i,j)*T1(i,j)/g(i,j)-1.0d0))+((u1(i,j)**2.0+v1(i,j)**2.0)/2.0d0)

    vn2 = vn1*(((g-1.0d0)*(M1**2)*(dsin(bt))**2)+2.0d0)/((g+1.0d0)*(M1**2)*(dsin(bt))**2)
    vv2 = dsqrt(vn2**2+vt**2)
    theta = bt-datan(vn2/vt)
    u2 = vv2*dcos(theta)
    v2 = vv2*dsin(-theta)
    T2 = T1*((2.0d0*g*(M1**2)*(dsin(bt))**2)-(g-1.0d0))*(((g-1.0d0)*(M1**2)*(dsin(bt))**2)+2.0d0)/&
           &(((g+1.0d0)**2)*(M1**2)*((dsin(bt))**2))
    rho2 = rho1*vn1/vn2
    p2 = p1*((2.0d0*g*(M1**2)*(dsin(bt))**2)-(g-1.0d0))/(g+1.0d0)
    e2 = rho2*((R*T2/g-1.0d0))+((u2**2.0+v2**2.0)/2.0d0)
  end do
end do

do j = 0,JM
  do i = 0,IM
    if ((y(0,JM-2)-dtan(bt)*x(i,j))>y(i,j)) then
      rho(i,j) = rho1(i,j) !ÕŒ‚”g‘O
      u(i,j) = u1(i,j)
      v(i,j) = v1(i,j)
      e(i,j) = e1(i,j)
      p(i,j) = p1(i,j)
    else
      rho(i,j) = rho2(i,j) !ÕŒ‚”gŒã
      u(i,j) = u2(i,j)
      v(i,j) = v2(i,j)
      e(i,j) = e2(i,j)
      p(i,j) = p2(i,j)
    end if
    print*,p,e
  end do
end do

!ƒOƒ‰ƒtì¬
open (1,file='assyuku4.dat',status='replace')
do j = 0,JM
  do i = 0,IM
    write (1,*) x(i,j),y(i,j),u(i,j),v(i,j)
  end do
end do
close (1)

end program assyuku4
