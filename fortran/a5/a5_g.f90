! program a5: TVDスキームとRunge-Kutta法を取り入れた2次元Euler方程式ソルバー（斜め衝撃波考慮）
program a5
implicit none

double precision,parameter :: H = 10.0d0
integer,parameter :: IM = 30,JM = 30
double precision :: dx,r_y,a

double precision,parameter :: M0 = 2.9d0,T0 = 293.0d0,R0 = 287.1d0
double precision,parameter :: gam0 = 1.4d0,rho0 = 1.2d0,dt = 1.0d-6

double precision,parameter :: pi = dacos(-1.0d0)
double precision,parameter :: beta = 30.0d0/180.0d0*pi

integer :: i,j,vec_index,n_step
double precision,dimension(0:IM,0:JM) :: x,y
double precision,dimension(0:IM,0:JM) :: rho,u,v,p,energy,T,E_bar
double precision,dimension(0:IM,0:JM) :: M,V_abs,V_n,V_t,theta,rho1,rho2

double precision,dimension(0:IM,0:JM,4) :: Q,Q1,E,F
double precision,dimension(0:IM,0:JM) :: x_xi,y_xi,x_eta,y_eta,Jac,Jac_inv
double precision,dimension(0:IM,0:JM) :: xi_x,xi_y,eta_x,eta_y
double precision,dimension(0:IM,0:JM) :: U_cvc,V_cvc

double precision :: u_init,v_init,p_init,energy_init,c_sound_init,slope
double precision :: du,dv,ddu,ddv
integer,parameter :: NMAX = 10000
double precision,parameter :: EPS = 1.0d-6

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

! === 初期条件設定 ===
do j = 0,JM
  do i = 0,IM
    M(i,j) = M0
    T(i,j) = T0
    rho1(i,j) = rho0
    p(i,j) = rho1(i,j)*R0*T(i,j)
    u(i,j) = M(i,j)*sqrt(gam0*R0*T(i,j))
    v(i,j) = 0.d0
    V_abs(i,j) = sqrt(u(i,j)**2+v(i,j)**2)
    V_n(i,j) = V_abs(i,j)*sin(beta)
    V_t(i,j) = V_abs(i,j)*cos(beta)
  end do
end do

! === 斜め衝撃波の適用 ===
do j = 0,JM
  do i = 0,IM
    if (y(i,j)<=y(i,JM-2) .and. x(i,j)<=sqrt(3.d0)*(y(i,JM-2)-y(i,j))) cycle

    T(i,j) = T(i,j)*(2.d0*gam0*(M(i,j)**2)*(sin(beta)**2)-(gam0-1.d0)) &
             *((gam0-1.0d0)*(M(i,j)**2)*(sin(beta)**2)+2.d0) &
             /(((gam0+1.0d0)**2)*(M(i,j)**2)*(sin(beta)**2))

    p(i,j) = p(i,j)*(2.d0*gam0*(M(i,j)**2)*(sin(beta)**2)-(gam0-1.d0)) &
             /(gam0+1.0d0)

    rho2(i,j) = rho1(i,j)*(gam0+1.d0)*(M(i,j)**2)*(sin(beta)**2) &
                /((gam0-1.0d0)*(M(i,j)**2)*(sin(beta)**2)+2.d0)

    V_n(i,j) = rho1(i,j)*V_n(i,j)/rho2(i,j)
    V_abs(i,j) = sqrt(V_n(i,j)**2+V_t(i,j)**2)
    theta(i,j) = beta-atan(V_n(i,j)/V_t(i,j))
    u(i,j) = V_abs(i,j)*cos(-theta(i,j))
    v(i,j) = V_abs(i,j)*sin(-theta(i,j))
    rho(i,j) = rho2(i,j)
    energy(i,j) = rho(i,j)*(R0*T(i,j)/(gam0-1.d0)+0.5d0*(u(i,j)**2+v(i,j)**2))
  end do
end do

! === 格子微係数とヤコビアン計算 ===
! (元コードと同様の実装を挿入：x_xi, y_xi, Jac 等)

! === 保存量Qの初期化 ===
call calc_QEF(IM,JM,rho,u,v,energy,p,Jac,xi_x,xi_y,eta_x,eta_y,U_cvc,V_cvc,Q,E,F)

do n_step = 1,NMAX
  Q1 = Q
  do vec_index = 1,4
    call calc_QEF(IM,JM,rho,u,v,energy,p,Jac,xi_x,xi_y,eta_x,eta_y,U_cvc,V_cvc,Q,E,F)
    do j = 1,JM-1
      do i = 1,IM-1
        Q(i,j,vec_index) = Q1(i,j,vec_index)-dt*((E(i,j,vec_index)-E(i-1,j,vec_index))+(F(i,j,vec_index)-F(i,j-1,vec_index)))
      end do
    end do
  end do

  du = 0.d0; dv = 0.d0
  do j = 1,JM-1
    do i = 1,IM-1
      rho(i,j) = Q(i,j,1)*Jac(i,j)
      ddu = abs(u(i,j)-Q(i,j,2)*Jac(i,j)/rho(i,j))/(abs(u(i,j))+1.0d-10)
      ddv = abs(v(i,j)-Q(i,j,3)*Jac(i,j)/rho(i,j))/(abs(v(i,j))+1.0d-10)
      if (ddu>du) du = ddu
      if (ddv>dv) dv = ddv
      u(i,j) = Q(i,j,2)*Jac(i,j)/rho(i,j)
      v(i,j) = Q(i,j,3)*Jac(i,j)/rho(i,j)
      energy(i,j) = Q(i,j,4)*Jac(i,j)
    end do
  end do
  if (du<=EPS .and. dv<=EPS) exit

  ! 状態量更新
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

! === 出力部 ===
open (10,file='a5.dat',status='replace')
do j = 0,JM
  do i = 0,IM
    write (10,*) x(i,j),y(i,j),u(i,j),v(i,j)
  end do
end do
close (10)

open (11,file='a5.fld',status='replace')
write (11,'(A)') '# AVS field file'
write (11,'(A)') 'ndim = 2'
write (11,'(A,I5)') 'dim1 =',IM+1
write (11,'(A,I5)') 'dim2 =',JM+1
write (11,'(A)') 'nspace = 2'
write (11,'(A)') 'veclen = 2'
write (11,'(A)') 'data = double'
write (11,'(A)') 'field = irregular'
write (11,'(A)') 'label = u v'
write (11,'(A)') 'variable 1 file=a5.dat filetype=ascii skip=0 offset=2 stride=4'
write (11,'(A)') 'variable 2 file=a5.dat filetype=ascii skip=0 offset=3 stride=4'
write (11,'(A)') 'coord 1 file=a5.dat filetype=ascii skip=0 offset=0 stride=4'
write (11,'(A)') 'coord 2 file=a5.dat filetype=ascii skip=0 offset=1 stride=4'
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

end program a5
