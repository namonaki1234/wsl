program kadai4
implicit none

double precision,parameter :: h = 10.0d0
double precision,parameter :: a = 1.2d0
double precision,parameter :: b = (a+1.0d0)/(a-1.0d0)
integer,parameter :: xm = 30
integer,parameter :: ym = 30
double precision,parameter :: r = 2.871d2 !ガス定数
double precision,parameter :: gamma = 1.4d0 !比熱比
double precision,parameter :: dt = 1.0d-6
double precision,parameter :: pi = dacos(0.0d0)*2.0d0
double precision,parameter :: beta = 3.0d1*pi/1.80d2
double precision,dimension(0:xm) :: x
double precision,dimension(0:ym) :: y
double precision,dimension(0:xm,0:ym) :: rho = 1.2d0,p = 0.0d0,u = 0.0d0,v = 0.0d0, &
  e = 0.0d0,eb,t,c,vn1,vn2,vt,v2
integer,dimension(0:xm,0:ym) :: flag
integer i,j,n
double precision dx,dy,m,alpha,ex,fy,jm,tanb,sinb,p1,rho1,cosb,theta

!for flag=1:after wave  flag=0:before wave
do i = 0,xm
  dx = (3.0d0*h)/dble(xm)
  x(i) = dble(i)*dx
end do
do j = 0,ym
  jm = dble(j)/(dble(ym))
  y(j) = h*(b**jm-1.0d0)/(b-1.0d0)
end do

tanb = dtan(beta)

do j = 0,ym
  do i = 0,xm
    if ((-tanb*x(i)+y(ym-2))<=y(j)) then
      flag(i,j) = 1
    else
      flag(i,j) = 0
    end if
  end do
end do

do j = 0,ym
  do i = 0,xm
    t(i,j) = 2.930d2
    c(i,j) = dsqrt(gamma*r*t(i,j)) !音速
  end do
end do
m = 2.9d0 !マッハ数

do j = 0,ym
  u(0,j) = m*c(0,j)
  u(1,j) = u(0,j)
end do

do j = 0,ym
  do i = 2,xm
    alpha = (x(i)-x(i-2))/(x(i-2)-x(i-1))
    u(i,j) = u(i-2,j)+alpha*(u(i-2,j)-u(i-1,j))
  end do
end do

do j = 0,ym
  do i = 0,xm
    if (flag(i,j)==0) then
      p(i,j) = rho(i,j)*r*t(i,j)
      e(i,j) = rho(i,j)*((r*t(i,j)/(gamma-1.0d0))+((u(i,j)**2.0)+v(i,j)**2.0)/2.0d0)
    end if
  end do
end do

sinb = dsin(beta)
cosb = dcos(beta)
p1 = 1.2d0*r*293.0d0
rho1 = 1.2d0
do j = 0,ym
  do i = 0,xm
    if (flag(i,j)==1) then
      t(i,j) = 293.0d0*(2.0d0*gamma*(m**2.0)*(sinb**2.0)-(gamma-1.0d0))* &
               ((gamma-1.0d0)*(m**2.0)*(sinb**2.0)+2.0d0)/(((gamma+1.0d0)**2.0)*(m**2.0)* &
                                                           (sinb**2.0))
      p(i,j) = p1*(2.0d0*gamma*(m**2.0)*(sinb**2.0)-(gamma-1.0d0))/(gamma+1.0d0)
      rho(i,j) = rho1*((gamma+1.0d0)*(m**2.0)*(sinb**2.0))/((gamma-1.0d0)* &
                                                            (m**2.0)*(sinb**2.0)+2.0d0)
      vn1(i,j) = u(i,j)*sinb
      vt(i,j) = vn1(i,j)/tanb
      vn2(i,j) = vn1(i,j)*((gamma-1.0d0)*(m**2.0)*(sinb**2.0)+2.0d0)/((gamma+1.0d0)* &
                                                                      (m**2.0)*(sinb**2.0))
      v2(i,j) = dsqrt((vn2(i,j)**2.0)+(vt(i,j)**2.0))
      theta = beta-datan((vn2(i,j)/vt(i,j)))
      u(i,j) = v2(i,j)*dcos(-theta)
      v(i,j) = v2(i,j)*dsin(-theta)
      e(i,j) = rho(i,j)*((r*t(i,j)/(gamma-1.0d0))+(u(i,j)**2.0))
    end if
  end do
end do

open (17,file='kadai4.dat',status='replace')
do j = 0,ym
  do i = 0,xm
    write (17,'(f13.6,$)') x(i)
    write (17,'(f13.6,$)') y(j)
    write (17,'(f13.6,$)') u(i,j)
    write (17,'(f13.6)') v(i,j)
  end do
end do
close (17)

end program kadai4
