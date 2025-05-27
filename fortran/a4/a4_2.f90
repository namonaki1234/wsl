program a4
implicit none

double precision,parameter :: H = 10.0d0
integer,parameter :: IM = 30,JM = 30
double precision :: dx,r_y,a,beta

double precision,parameter :: M0 = 2.9d0,T0 = 293.0d0,R0 = 287.1d0
double precision,parameter :: gam0 = 1.4d0,rho0 = 1.2d0,dt = 1.0d-6

integer :: i,j
double precision,dimension(0:IM,0:JM) :: x,y
double precision,dimension(0:IM,0:JM) :: rho,u,v,p,energy,T

double precision,dimension(0:IM,0:JM) :: x_xi,y_xi,x_eta,y_eta,Jac,Jac_inv
double precision,dimension(0:IM,0:JM) :: xi_x,xi_y,eta_x,eta_y

double precision :: u_init,v_init,p_init,energy_init,c_sound_init,slope
double precision,parameter::pi=3.14159265358979d0 !π

dx = 3.0d0*H/dble(IM)
r_y = 1.1d0

! === 格子生成 ===
do i = 0,IM
  do j = 0,JM
    x(i,j) = dx*dble(i)
  end do
end do

do j = 0,JM
  a = H*(r_y-1.0d0)/(r_y**dble(JM)-1.0d0)
  do i = 0,IM
    y(i,j) = a*(r_y**dble(j)-1.0d0)/(r_y-1.0d0)
  end do
end do

!初期条件
beta=pi/6.0d0 !衝撃波角β
slope=((gam0+1.0d0)*(M0*dsin(beta))**2.0)/((gam0-1.0d0)*(M0*dsin(beta))**2.0+2.0d0) !係数
c_sound_init = dsqrt(gam0*R0*T0)
u_init = M0*c_sound_init
v_init = 0.0d0
p_init = rho0*R0*T0
energy_init = rho0*(R0*T0/(gam0-1.0d0)+0.5d0*(u_init**2+v_init**2))

do j = 0,JM
  do i = 0,IM
    rho(i,j) = rho0
    u(i,j) = u_init
    v(i,j) = v_init
    p(i,j) = p_init
    energy(i,j) = energy_init
    T(i,j) = T0
  end do
end do

do i = 0,IM
  do j = 1,JM-1
    x_eta(i,j) = (x(i,j+1)-x(i,j-1))/2.0d0
    y_eta(i,j) = (y(i,j+1)-y(i,j-1))/2.0d0
  end do
end do
do j = 0,JM
  do i = 1,IM-1
    x_xi(i,j) = (x(i+1,j)-x(i-1,j))/2.0d0
    y_xi(i,j) = (y(i+1,j)-y(i-1,j))/2.0d0
  end do
end do
do j = 0,JM
  x_xi(0,j) = (-3.0d0*x(0,j)+4.0d0*x(1,j)-x(2,j))/2.0d0
  x_xi(IM,j) = -(-3.0d0*x(IM,j)+4.0d0*x(IM-1,j)-x(IM-2,j))/2.0d0
  y_xi(0,j) = (-3.0d0*y(0,j)+4.0d0*y(1,j)-y(2,j))/2.0d0
  y_xi(IM,j) = -(-3.0d0*y(IM,j)+4.0d0*y(IM-1,j)-y(IM-2,j))/2.0d0
end do
do i = 0,IM
  x_eta(i,0) = (-3.0d0*x(i,0)+4.0d0*x(i,1)-x(i,2))/2.0d0
  x_eta(i,JM) = -(-3.0d0*x(i,JM)+4.0d0*x(i,JM-1)-x(i,JM-2))/2.0d0
  y_eta(i,0) = (-3.0d0*y(i,0)+4.0d0*y(i,1)-y(i,2))/2.0d0
  y_eta(i,JM) = -(-3.0d0*y(i,JM)+4.0d0*y(i,JM-1)-y(i,JM-2))/2.0d0
end do

do i = 0,IM
  do j = 0,JM
    Jac_inv(i,j) = x_xi(i,j)*y_eta(i,j)-y_xi(i,j)*x_eta(i,j)
    if (abs(Jac_inv(i,j))<1.0d-12) then
      print*,"エラー：ヤコビアンがゼロ i=",i," j=",j," → jac_inv=",Jac_inv(i,j)
      stop
    end if
    Jac(i,j) = 1.0d0/(x_xi(i,j)*y_eta(i,j)-y_xi(i,j)*x_eta(i,j))
    xi_x(i,j) = Jac(i,j)*y_eta(i,j)
    xi_y(i,j) = -Jac(i,j)*x_eta(i,j)
    eta_x(i,j) = -Jac(i,j)*y_xi(i,j)
    eta_y(i,j) = Jac(i,j)*x_xi(i,j)
  end do
end do

! === データ出力（MicroAVS用 DAT）===
open (10,file='a3.dat',status='replace')
do j = 0,JM
  do i = 0,IM
    write (10,*) x(i,j),y(i,j),u(i,j),v(i,j)
  end do
end do
close (10)

! === MicroAVSのFLDヘッダ出力 ===
open (11,file='a4.fld',status='replace')
write (11,'(A)') '# AVS field file'
write (11,'(A)') 'ndim = 2'
write (11,'(A,I5)') 'dim1 =',IM+1
write (11,'(A,I5)') 'dim2 =',JM+1
write (11,'(A)') 'nspace = 2'
write (11,'(A)') 'veclen = 2'
write (11,'(A)') 'data = double'
write (11,'(A)') 'field = irregular'
write (11,'(A)') 'label = u v'
write (11,'(A)') 'variable 1 file=a4.dat filetype=ascii skip=0 offset=2 stride=4'
write (11,'(A)') 'variable 2 file=a4.dat filetype=ascii skip=0 offset=3 stride=4'
write (11,'(A)') 'coord 1 file=a4.dat filetype=ascii skip=0 offset=0 stride=4'
write (11,'(A)') 'coord 2 file=a4.dat filetype=ascii skip=0 offset=1 stride=4'
close (11)

print*,"→ MicroAVS用の .dat および .fld を出力しました。"

end program a4
