program a3
implicit none

double precision,parameter :: H = 10.0d0
integer,parameter :: IM = 30,JM = 30
double precision :: dx,r_y,a

double precision,parameter :: M0 = 2.9d0,T0 = 293.0d0,R0 = 287.1d0
double precision,parameter :: gam0 = 1.4d0,rho0 = 1.2d0,dt = 1.0d-6

integer :: i,j,vec_index,n_step
double precision,dimension(0:IM,0:JM) :: x,y
double precision,dimension(0:IM,0:JM) :: rho,u,v,p,energy,T,E_bar

double precision,dimension(0:IM,0:JM,4) :: Q,E,F

double precision,dimension(0:IM,0:JM) :: x_xi,y_xi,x_eta,y_eta,Jac,Jac_inv
double precision,dimension(0:IM,0:JM) :: xi_x,xi_y,eta_x,eta_y
double precision,dimension(0:IM,0:JM) :: U_cvc,V_cvc

double precision :: u_init,v_init,p_init,energy_init,c_sound_init

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

do n_step = 1,10000
  call calc_QEF(IM,JM,rho,u,v,energy,p,Jac,xi_x,xi_y,eta_x,eta_y,U_cvc,V_cvc,Q,E,F)

  do vec_index = 1,4
    do j = 1,JM-1
      do i = 1,IM-1
        Q(i,j,vec_index) = Q(i,j,vec_index)-dt*( &
                           (E(i+1,j,vec_index)-E(i-1,j,vec_index))/(2.0d0)+ &
                           (F(i,j+1,vec_index)-F(i,j-1,vec_index))/(2.0d0))
      end do
    end do
  end do

  do j = 1,JM-1
    do i = 1,IM-1
      rho(i,j) = Q(i,j,1)*Jac(i,j)
      u(i,j) = (Q(i,j,2)*Jac(i,j))/rho(i,j)
      v(i,j) = (Q(i,j,3)*Jac(i,j))/rho(i,j)
      energy(i,j) = Q(i,j,4)*Jac(i,j)
    end do
  end do

!   一次外挿で流入を固定
  do vec_index = 1,4
    do i = 1,IM-1
      Q(i,JM,vec_index) = 2.0d0*Q(i,JM-1,vec_index)-Q(i,JM-2,vec_index)
    end do
  end do
  do vec_index = 1,4
    do i = 1,IM-1
      Q(i,0,vec_index) = 2.0d0*Q(i,1,vec_index)-Q(i,2,vec_index)
    end do
  end do
  do vec_index = 1,4
    do j = 0,JM
      Q(IM,j,vec_index) = 2.0d0*Q(IM-1,j,vec_index)-Q(IM-2,j,vec_index)
    end do
  end do

  do i = 1,IM-1
    rho(i,JM) = Q(i,JM,1)*Jac(i,JM)
    u(i,JM) = (Q(i,JM,2)*Jac(i,JM))/rho(i,JM)
    v(i,JM) = 0.0d0
    energy(i,JM) = (Q(i,JM,4)*Jac(i,JM))
    E_bar(i,JM) = energy(i,JM)/rho(i,JM)-0.5d0*(u(i,JM)**2+v(i,JM)**2)
    p(i,JM) = (gam0-1.0d0)*rho(i,JM)*E_bar(i,JM)
    energy(i,JM) = rho(i,JM)*(E_bar(i,JM)+0.5d0*(u(i,JM)**2+v(i,JM)**2))
  end do
  do i = 1,IM-1
    rho(i,0) = Q(i,0,1)*Jac(i,0)
    u(i,0) = (Q(i,0,2)*Jac(i,0))/rho(i,0)
    v(i,0) = 0.0d0
    energy(i,0) = (Q(i,0,4)*Jac(i,0))
    E_bar(i,0) = energy(i,0)/rho(i,0)-0.5d0*(u(i,0)**2+v(i,0)**2)
    p(i,0) = (gam0-1.0d0)*rho(i,0)*E_bar(i,0)
    energy(i,0) = rho(i,0)*(E_bar(i,0)+0.5d0*(u(i,0)**2+v(i,0)**2))
  end do
  do j = 0,JM
    rho(IM,j) = Q(IM,j,1)*Jac(IM,j)
    u(IM,j) = (Q(IM,j,2)*Jac(IM,j))/rho(IM,j)
    v(IM,j) = (Q(IM,j,3)*Jac(IM,j))/rho(IM,j)
    energy(IM,j) = Q(IM,j,4)*Jac(IM,j)
  end do

!   状態方程式の再計算
  do j = 0,JM
    do i = 0,IM
      E_bar(i,j) = energy(i,j)/rho(i,j)-(u(i,j)**2+v(i,j)**2)/2.0d0
      p(i,j) = (gam0-1.0d0)*rho(i,j)*E_bar(i,j)
      T(i,j) = p(i,j)/(rho(i,j)*R0)
    end do
  end do

  if (mod(n_step,100)==0) then
    print*,'n_step = ',n_step,', time = ',dble(n_step)*dt
  end if
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
open (11,file='a3.fld',status='replace')
write (11,'(A)') '# AVS field file'
write (11,'(A)') 'ndim = 2'
write (11,'(A,I5)') 'dim1 =',IM+1
write (11,'(A,I5)') 'dim2 =',JM+1
write (11,'(A)') 'nspace = 2'
write (11,'(A)') 'veclen = 2'
write (11,'(A)') 'data = double'
write (11,'(A)') 'field = irregular'
write (11,'(A)') 'label = u v'
write (11,'(A)') 'variable 1 file=a3.dat filetype=ascii skip=0 offset=2 stride=4'
write (11,'(A)') 'variable 2 file=a3.dat filetype=ascii skip=0 offset=3 stride=4'
write (11,'(A)') 'coord 1 file=a3.dat filetype=ascii skip=0 offset=0 stride=4'
write (11,'(A)') 'coord 2 file=a3.dat filetype=ascii skip=0 offset=1 stride=4'
close (11)

print*,"→ MicroAVS用の .dat および .fld を出力しました。"

contains

subroutine calc_QEF(IM,JM,rho,u,v,energy,p, &
                    Jac,xi_x,xi_y,eta_x,eta_y, &
                    U_cvc,V_cvc,Q,E,F)
  implicit none
  integer,intent(in) :: IM,JM
  double precision,intent(in) :: rho(0:IM,0:JM),u(0:IM,0:JM),v(0:IM,0:JM)
  double precision,intent(in) :: energy(0:IM,0:JM),p(0:IM,0:JM)
  double precision,intent(in) :: Jac(0:IM,0:JM)
  double precision,intent(in) :: xi_x(0:IM,0:JM),xi_y(0:IM,0:JM)
  double precision,intent(in) :: eta_x(0:IM,0:JM),eta_y(0:IM,0:JM)

  double precision,intent(out) :: U_cvc(0:IM,0:JM),V_cvc(0:IM,0:JM)
  double precision,intent(out) :: Q(0:IM,0:JM,4)
  double precision,intent(out) :: E(0:IM,0:JM,4)
  double precision,intent(out) :: F(0:IM,0:JM,4)
  integer :: i,j

  do j = 0,JM
    do i = 0,IM
      U_cvc(i,j) = xi_x(i,j)*u(i,j)+xi_y(i,j)*v(i,j)
      V_cvc(i,j) = eta_x(i,j)*u(i,j)+eta_y(i,j)*v(i,j)

      Q(i,j,1) = rho(i,j)/Jac(i,j)
      Q(i,j,2) = (rho(i,j)*u(i,j))/Jac(i,j)
      Q(i,j,3) = (rho(i,j)*v(i,j))/Jac(i,j)
      Q(i,j,4) = energy(i,j)/Jac(i,j)

      E(i,j,1) = (rho(i,j)*U_cvc(i,j))/Jac(i,j)
      E(i,j,2) = (rho(i,j)*u(i,j)*U_cvc(i,j)+xi_x(i,j)*p(i,j))/Jac(i,j)
      E(i,j,3) = (rho(i,j)*v(i,j)*U_cvc(i,j)+xi_y(i,j)*p(i,j))/Jac(i,j)
      E(i,j,4) = (energy(i,j)+p(i,j))*U_cvc(i,j)/Jac(i,j)

      F(i,j,1) = (rho(i,j)*V_cvc(i,j))/Jac(i,j)
      F(i,j,2) = (rho(i,j)*u(i,j)*V_cvc(i,j)+eta_x(i,j)*p(i,j))/Jac(i,j)
      F(i,j,3) = (rho(i,j)*v(i,j)*V_cvc(i,j)+eta_y(i,j)*p(i,j))/Jac(i,j)
      F(i,j,4) = (energy(i,j)+p(i,j))*V_cvc(i,j)/Jac(i,j)
    end do
  end do
end subroutine calc_QEF

end program a3
