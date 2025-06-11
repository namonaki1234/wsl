program a6

implicit none
integer :: i,j,k,k_runge_kutta,n
double precision :: a
double precision :: M =2.9d0 !マッハ数
! double precision :: rho1,rho2,u1,u2,v1,v2,p1,p2,T1,T2,M,V_abs1,V_abs2,V_n1,V_n2,V_t,theta,energy1,energy2
!********************格子空間・計算条件・初期条件********************
double precision,parameter :: H = 10.0d0
double precision :: du,dv,COE
! double precision :: slope
integer,parameter :: IM = 30,JM = 30
double precision,parameter :: r_y = 1.1d0,dx = 3.0d0*H/dble(IM)
double precision,parameter :: dt = 1.0d-8
double precision,parameter :: g = 1.4d0
double precision,parameter :: R = 287.1d0
double precision,parameter :: pi = dacos(-1.0d0)
double precision,parameter :: beta = pi/6.0d0

!********************配列定義********************
!-------衝撃波-------
! double precision,dimension(0:IM,0:JM) :: E_bar = 0.0d0
double precision,dimension(0:IM,0:JM) :: rho = 0.0d0 !密度
double precision,dimension(0:IM,0:JM) :: u = 0.0d0 !x方向速度
double precision,dimension(0:IM,0:JM) :: v = 0.0d0 !y方向速度
double precision,dimension(0:IM,0:JM) :: p = 0.0d0 !圧力
double precision,dimension(0:IM,0:JM) :: T = 0.0d0 !温度
double precision,dimension(0:IM,0:JM) :: energy = 0.0d0 !エネルギー
double precision,dimension(0:IM,0:JM,4) :: Q,Q1,E,F
double precision :: us !速度uの値保持用

!-------ヤコビアンJa-------
double precision,dimension(0:IM,0:JM) :: x = 0.0d0,y = 0.0d0,x_xi = 0.0d0,x_eta = 0.0d0,y_xi = 0.0d0,y_eta = 0.0d0
double precision,dimension(0:IM,0:JM) :: xi_x,xi_y,eta_x,eta_y,Jac,Jac_inv
integer,parameter :: NMAX = 100000
double precision,parameter :: EPS = 1.0d-6
double precision,dimension(1:4,0:IM,0:JM) :: Rv = 0.0d0,Sv = 0.0d0 !粘性項ベクトルR,S
double precision :: E_s !Eの値保持用

!********************不等間隔格子の生成・ファイル出力********************
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

!********************一般座標へ変換・初期値設定********************
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

!********************ヤコビアンJaの計算********************
do i = 0,IM
  do j = 0,JM
    Jac_inv(i,j) = x_xi(i,j)*y_eta(i,j)-y_xi(i,j)*x_eta(i,j)
    if (abs(Jac_inv(i,j))<1.0d-12) then
      print*,"エラー：ヤコビアンの分母がゼロ i=",i," j=",j," → jac_inv=",Jac_inv(i,j)
      stop
    end if
    Jac(i,j) = 1.0d0/(x_xi(i,j)*y_eta(i,j)-y_xi(i,j)*x_eta(i,j))
    xi_x(i,j) = Jac(i,j)*y_eta(i,j)
    xi_y(i,j) = -Jac(i,j)*x_eta(i,j)
    eta_x(i,j) = -Jac(i,j)*y_xi(i,j)
    eta_y(i,j) = Jac(i,j)*x_xi(i,j)
  end do
end do

!初期条件***************************************************
COE = ((g+1.0d0)*(M*dsin(beta))**2.0)/((g-1.0d0)*(M*dsin(beta))**2.0+2.0d0) !係数
do i = 0,IM
  do j = 0,JM
    rho(i,j) = 1.2d0
    T(i,j) = 293.0d0
    u(i,j) = M*dsqrt(g*R*T(i,j))
    p(i,j) = rho(i,j)*R*T(i,j)
    energy(i,j) = rho(i,j)*(R*T(i,j)/(g-1.0d0)+(u(i,j)**2.0+v(i,j)**2.0)/2.0d0)

    if (y(i,j)>=y(0,JM-2)-x(i,j)/dsqrt(3.0d0)) then !衝撃波後の物理量
      us = u(i,j)
      rho(i,j) = COE*rho(i,j)
      u(i,j) = u(i,j)*dsin(beta)*dsqrt((1.0d0/COE)**2.0+(1.0d0/dtan(beta))**2d0)*dcos(-beta+datan(dtan(beta)/COE))
      v(i,j) = us*dsin(beta)*dsqrt((1.0d0/COE)**2.0+(1.0d0/dtan(beta))**2.0)*dsin(-beta+datan(dtan(beta)/COE))
      p(i,j) = (2.0d0*g*(M*dsin(beta))**2.0-(g-1.0d0))/(g+1.0d0)*p(i,j)
     T(i,j) = (2.0d0*g*(M*dsin(beta))**2.0d0-(g-1.0d0))*((g-1.0d0)*(M*dsin(beta))**2.0+2.0d0)/(((g+1.0d0)*M*dsin(beta))**2.0)*T(i,j)
      energy(i,j) = rho(i,j)*(R*T(i,j)/(g-1.0d0)+(u(i,j)**2.0+v(i,j)**2.0)/2.0d0)
    end if

    Q(i,j,1) = rho(i,j)
    Q(i,j,2) = rho(i,j)*u(i,j)
    Q(i,j,3) = rho(i,j)*v(i,j)
    Q(i,j,4) = energy(i,j)
    E(i,j,1) = rho(i,j)*u(i,j)
    E(i,j,2) = p(i,j)+rho(i,j)*u(i,j)**2.0
    E(i,j,3) = rho(i,j)*u(i,j)*v(i,j)
    E(i,j,4) = (energy(i,j)+p(i,j))*u(i,j)
    F(i,j,1) = rho(i,j)*v(i,j)
    F(i,j,2) = rho(i,j)*u(i,j)*v(i,j)
    F(i,j,3) = p(i,j)+rho(i,j)*v(i,j)**2.0
    F(i,j,4) = (energy(i,j)+p(i,j))*v(i,j)

    do k = 1,4 !各流束を一般座標系に変換
      Q(i,j,k) = Q(i,j,k)/Jac(i,j)
      Q1(i,j,k) = Q(i,j,k)
      E_s = E(i,j,k)
      E(i,j,k) = (xi_x(i,j)*E(i,j,k)+xi_y(i,j)*F(i,j,k))/Jac(i,j)
      F(i,j,k) = (eta_x(i,j)*E_s+eta_y(i,j)*F(i,j,k))/Jac(i,j)
    end do
  end do
end do
!********************ルンゲクッタ法／メイン計算********************
do n = 1,NMAX
  do k = 1,4
    do j = 0,JM
      do i = 0,IM
        Q1(i,j,k) = Q(i,j,k)
      end do
    end do
  end do
  do k_runge_kutta = 1,4
    call tvd_xi
    call tvd_eta
    call viscosity_calc
    call runge_kutta(1,IM-1,1,JM-1)
    call boundary_condition
  end do
  if (mod(n,1000)==0) then
    print*,"n du dv",n,du,dv
  end if
  ! if (du<=EPS .and. dv<=EPS) exit
end do
call Save_data

contains
!********************サブルーチン5　境界条件********************
subroutine boundary_condition
  do j = 0,JM
    do i = 0,IM
      !境界条件の更新
      if (i==IM) then !流出→ノイマン境界条件（境界法線方向差分0）
        rho(IM,j) = rho(IM-1,j)
        u(IM,j) = u(IM-1,j)
        v(IM,j) = v(IM-1,j)
        p(IM,j) = p(IM-1,j)
        T(IM,j) = p(IM,j)/(rho(IM,j)*R)
        energy(IM,j) = rho(IM,j)*(R*T(IM,j)/(g-1.0d0)+(u(IM,j)**2.0+v(IM,j)**2.0)/2.0d0)
      end if
      if (j==JM) then !上端→ノイマン境界条件（境界法線方向差分0）
        rho(i,JM) = rho(i,JM-1)
        u(i,JM) = u(i,JM-1)
        v(i,JM) = v(i,JM-1)
        p(i,JM) = p(i,JM-1)
        T(i,JM) = p(i,JM)/(rho(i,JM)*R)
        energy(i,JM) = rho(i,JM)*(R*T(i,JM)/(g-1.0d0)+(u(i,JM)**2.0+v(i,JM)**2.0)/2.0d0)
      end if
      if (j==0) then !下端（壁）→粘性を考慮するためu=0，他はノイマン境界条件（境界法線方向差分0）
        rho(i,0) = rho(i,1)
        u(i,0) = 0.0d0
        v(i,0) = 0.0d0
        p(i,0) = p(i,1)
        T(i,0) = p(i,0)/(rho(i,0)*R)
        energy(i,0) = rho(i,0)*(R*T(i,0)/(g-1.0d0)+(u(i,0)**2.0+v(i,0)**2.0)/2.0d0)
      end if
      !流束ベクトルの更新
      Q(i,j,1) = rho(i,j)
      Q(i,j,2) = rho(i,j)*u(i,j)
      Q(i,j,3) = rho(i,j)*v(i,j)
      Q(i,j,4) = energy(i,j)
      E(i,j,1) = rho(i,j)*u(i,j)
      E(i,j,2) = p(i,j)+rho(i,j)*u(i,j)**2.0
      E(i,j,3) = rho(i,j)*u(i,j)*v(i,j)
      E(i,j,4) = (energy(i,j)+p(i,j))*u(i,j)
      F(i,j,1) = rho(i,j)*v(i,j)
      F(i,j,2) = rho(i,j)*u(i,j)*v(i,j)
      F(i,j,3) = p(i,j)+rho(i,j)*v(i,j)**2.0
      F(i,j,4) = (energy(i,j)+p(i,j))*v(i,j)

      do k = 1,4 !各流束を一般座標系に変換
        Q(i,j,k) = Q(i,j,k)/Jac(i,j)
        E_s = E(i,j,k)
        E(i,j,k) = (xi_x(i,j)*E(i,j,k)+xi_y(i,j)*F(i,j,k))/Jac(i,j)
        F(i,j,k) = (eta_x(i,j)*E_s+eta_y(i,j)*F(i,j,k))/Jac(i,j)
      end do
    end do

  end do
end subroutine boundary_condition

!********************ファンクション　FPSI制限関数の定義********************
double precision function FPSI(Z,DELTA)
  double precision DELTA,Z
  if (abs(Z)>=DELTA) then
    FPSI = abs(Z)
  else
    FPSI = 0.5d0*(Z**2+DELTA**2)/DELTA
  end if
  return
end function FPSI

!********************サブルーチン1　TVDスキーム_ξ方向********************
!*****Eを求める*****
subroutine tvd_xi
  double precision RWL,RWR,AKX,AKY,AJACM,UM,VM,DUM,CM
  double precision PHI,BETA,AKXT,AKYT,THIT
  double precision,dimension(1:4,1:4,0:IM) :: RR = 0.0d0,RI = 0.0d0
  double precision,dimension(1:4,0:IM) :: EIGM,GG,PHIM,EH
  double precision,dimension(1:4) :: D
  double precision,dimension(1:4,-1:IM) :: ALPHA
  double precision S,DELTA,GAMMA,UCONT

  do j = 1,JM-1
    do i = 0,IM-1
      RWL = dsqrt(Q(i,j,1)*Jac(i,j))/(dsqrt(Q(i+1,j,1)*Jac(i+1,j))+dsqrt(Q(i,j,1)*Jac(i,j)))
      RWR = dsqrt(Q(i+1,j,1)*Jac(i+1,j))/(dsqrt(Q(i+1,j,1)*Jac(i+1,j))+dsqrt(Q(i,j,1)*Jac(i,j)))
      AKX = 0.5d0*(xi_x(i,j)+xi_x(i+1,j))
      AKY = 0.5d0*(xi_y(i,j)+xi_y(i+1,j))
      AJACM = 0.5d0*(Jac(i,j)+Jac(i+1,j))
      UM = RWL*u(i,j)+RWR*u(i+1,j)
      VM = RWL*v(i,j)+RWR*v(i+1,j)
      DUM = RWL*Q(i,j,4)/Q(i,j,1)+RWR*Q(i+1,j,4)/Q(i+1,j,1)
      CM = dsqrt(g*(g-1.0d0)*abs(DUM-0.5d0*(UM**2+VM**2)))
      !-------NOMENCLATURE-------
      PHI = 0.5d0*(g-1.0d0)*(UM**2+VM**2)
      BETA = 1.0d0/(2.0d0*CM**2)
      AKXT = AKX/dsqrt(AKX**2+AKY**2)
      AKYT = AKY/dsqrt(AKX**2+AKY**2)
      THIT = AKXT*UM+AKYT*VM
      !-------RIGHT ENGEN-VECTORS-------　P.88 R_k
      RR(1,1,i) = 1.0d0
      RR(1,2,i) = 0.0d0
      RR(1,3,i) = 1.0d0
      RR(1,4,i) = 1.0d0
      RR(2,1,i) = UM
      RR(2,2,i) = AKYT
      RR(2,3,i) = UM+AKXT*CM
      RR(2,4,i) = UM-AKXT*CM
      RR(3,1,i) = VM
      RR(3,2,i) = -AKXT
      RR(3,3,i) = VM+AKYT*CM
      RR(3,4,i) = VM-AKYT*CM
      RR(4,1,i) = PHI/(g-1.0d0)
      RR(4,2,i) = AKYT*UM-AKXT*VM
      RR(4,3,i) = (PHI+CM**2)/(g-1.0d0)+CM*THIT
      RR(4,4,i) = (PHI+CM**2)/(g-1.0d0)-CM*THIT
      !-------INVERS OR RIGHT EIGEN-VECTORS-------　 P.88 R_k^-1
      RI(1,1,i) = 1.0d0-PHI/CM**2
      RI(1,2,i) = (g-1.0d0)*UM/CM**2
      RI(1,3,i) = (g-1.0d0)*VM/CM**2
      RI(1,4,i) = -(g-1.0d0)/CM**2
      RI(2,1,i) = -AKYT*UM+AKXT*VM
      RI(2,2,i) = AKYT
      RI(2,3,i) = -AKXT
      RI(2,4,i) = 0.0d0
      RI(3,1,i) = BETA*(PHI-CM*THIT)
      RI(3,2,i) = BETA*(AKXT*CM-(g-1.0d0)*UM)
      RI(3,3,i) = BETA*(AKYT*CM-(g-1.0d0)*VM)
      RI(3,4,i) = BETA*(g-1.0d0)
      RI(4,1,i) = BETA*(PHI+CM*THIT)
      RI(4,2,i) = -BETA*(AKXT*CM+(g-1.0d0)*UM)
      RI(4,3,i) = -BETA*(AKYT*CM+(g-1.0d0)*VM)
      RI(4,4,i) = BETA*(g-1.0d0)
      !-------ENGEN-VALUES AT INTERMEDIATE-------　 P.87  固有値
      EIGM(1,i) = AKX*UM+AKY*VM
      EIGM(2,i) = EIGM(1,i)
      EIGM(3,i) = EIGM(1,i)+CM*dsqrt(AKX**2+AKY**2)
      EIGM(4,i) = EIGM(1,i)-CM*dsqrt(AKX**2+AKY**2)
      !-------ALPHA-------
      do k = 1,4
        D(k) = (Q(i+1,j,k)*Jac(i+1,j)-Q(i,j,k)*Jac(i,j))/AJACM
      end do
      do k = 1,4
        ALPHA(k,i) = RI(k,1,i)*D(1)+RI(k,2,i)*D(2)+RI(k,3,i)*D(3)+RI(k,4,i)*D(4)
      end do
    end do
    do k = 1,4
      ALPHA(k,-1) = ALPHA(k,0)
      ALPHA(k,IM) = ALPHA(k,IM-1)
    end do

    do k = 1,4
      do i = 0,IM
        !-------MINMOD LIMITER-------
        S = DSIGN(1.0d0,ALPHA(k,i))
        GG(k,i) = S*DMAX1(0.0d0,DMIN1(abs(ALPHA(k,i)),S*ALPHA(k,i-1)))
      end do
    end do
    do k = 1,4
      do i = 0,IM-1
        !-------GAMMA-------
        DELTA = 1.0d-10
        !DELTA=AMAX1(0.0,EIGM(K,J)-EIG(K,J),EIG(K,J+1)-EIGM(K,J))
        if (ALPHA(k,i)/=0.0d0) then
          GAMMA = 0.5d0*FPSI(EIGM(k,i),DELTA)*(GG(k,i+1)-GG(k,i))/ALPHA(k,i)
        else
          GAMMA = 0.0d0
        end if
        !-------PHI-------
        PHIM(k,i) = 0.5d0*FPSI(EIGM(k,i),DELTA)*(GG(k,i+1)+GG(k,i))-FPSI(EIGM(k,i)+GAMMA,DELTA)*ALPHA(k,i)
      end do
    end do
    !-------CONVECTION COMPORNENTS-------
    do i = 0,IM
      UCONT = xi_x(i,j)*u(i,j)+xi_y(i,j)*v(i,j)
      EH(1,i) = Q(i,j,1)*UCONT
      EH(2,i) = Q(i,j,2)*UCONT+p(i,j)*xi_x(i,j)/Jac(i,j)
      EH(3,i) = Q(i,j,3)*UCONT+p(i,j)*xi_y(i,j)/Jac(i,j)
      EH(4,i) = Q(i,j,4)*UCONT+p(i,j)*UCONT/Jac(i,j)
    end do
    !-------XI-DIRECTION CONVECTION FLUX-------
    do k = 1,4
      do i = 0,IM-1
        E(i,j,k) = 0.5d0*(EH(k,i)+EH(k,i+1)+RR(k,1,i)*PHIM(1,i)+RR(k,2,i)*PHIM(2,i) &
                          +RR(k,3,i)*PHIM(3,i)+RR(k,4,i)*PHIM(4,i))
      end do
    end do
  end do
end subroutine tvd_xi

!********************サブルーチン2　TVDスキーム_η方向********************
!*****Fを求める*****
subroutine tvd_eta
  double precision RWL,RWR,AKX,AKY,AJACM,UM,VM,DUM,CM
  double precision PHI,BETA,AKXT,AKYT,THIT
  double precision,dimension(1:4,1:4,0:JM) :: RR,RJ
  double precision,dimension(1:4,0:JM) :: EIGM,GG,PHIM,FH
  !   double precision,dimension(1:4,0:JM) :: EIG
  double precision,dimension(1:4) :: D
  double precision,dimension(1:4,-1:JM) :: ALPHA
  double precision S,DELTA,GAMMA,VCONT

  do i = 1,IM-1
    do j = 0,JM-1
      RWL = dsqrt(Q(i,j,1)*Jac(i,j))/(dsqrt(Q(i,j+1,1)*Jac(i,j+1))+dsqrt(Q(i,j,1)*Jac(i,j)))
      RWR = dsqrt(Q(i,j+1,1)*Jac(i,j+1))/(dsqrt(Q(i,j+1,1)*Jac(i,j+1))+dsqrt(Q(i,j,1)*Jac(i,j)))
      AKX = 0.5d0*(eta_x(i,j)+eta_x(i,j+1))
      AKY = 0.5d0*(eta_y(i,j)+eta_y(i,j+1))
      AJACM = 0.5d0*(Jac(i,j)+Jac(i,j+1))
      UM = RWL*u(i,j)+RWR*u(i,j+1)
      VM = RWL*v(i,j)+RWR*v(i,j+1)
      DUM = RWL*Q(i,j,4)/Q(i,j,1)+RWR*Q(i,j+1,4)/Q(i,j+1,1)
      CM = dsqrt(g*(g-1.0d0)*abs(DUM-0.5d0*(UM**2+VM**2)))
      !NOMENCLATURE
      PHI = 0.5d0*(g-1.0d0)*(UM**2+VM**2)
      BETA = 1.0d0/(2.0d0*CM**2)
      AKXT = AKX/dsqrt(AKX**2+AKY**2)
      AKYT = AKY/dsqrt(AKX**2+AKY**2)
      THIT = AKXT*UM+AKYT*VM
      !RIGHT ENGEN-VECTORS
      RR(1,1,j) = 1.0d0
      RR(1,2,j) = 0.0d0
      RR(1,3,j) = 1.0d0
      RR(1,4,j) = 1.0d0
      RR(2,1,j) = UM
      RR(2,2,j) = AKYT
      RR(2,3,j) = UM+AKXT*CM
      RR(2,4,j) = UM-AKXT*CM
      RR(3,1,j) = VM
      RR(3,2,j) = -AKXT
      RR(3,3,j) = VM+AKYT*CM
      RR(3,4,j) = VM-AKYT*CM
      RR(4,1,j) = PHI/(g-1.0d0)
      RR(4,2,j) = AKYT*UM-AKXT*VM
      RR(4,3,j) = (PHI+CM**2)/(g-1.0d0)+CM*THIT
      RR(4,4,j) = (PHI+CM**2)/(g-1.0d0)-CM*THIT
      !INVERS OR RIGHT EIGEN-VECTORS
      RJ(1,1,j) = 1.0d0-PHI/CM**2
      RJ(1,2,j) = (g-1.0d0)*UM/CM**2
      RJ(1,3,j) = (g-1.0d0)*VM/CM**2
      RJ(1,4,j) = -(g-1.0d0)/CM**2
      RJ(2,1,j) = -AKYT*UM+AKXT*VM
      RJ(2,2,j) = AKYT
      RJ(2,3,j) = -AKXT
      RJ(2,4,j) = 0.0d0
      RJ(3,1,j) = BETA*(PHI-CM*THIT)
      RJ(3,2,j) = BETA*(AKXT*CM-(g-1.0d0)*UM)
      RJ(3,3,j) = BETA*(AKYT*CM-(g-1.0d0)*VM)
      RJ(3,4,j) = BETA*(g-1.0d0)
      RJ(4,1,j) = BETA*(PHI+CM*THIT)
      RJ(4,2,j) = -BETA*(AKXT*CM+(g-1.0d0)*UM)
      RJ(4,3,j) = -BETA*(AKYT*CM+(g-1.0d0)*VM)
      RJ(4,4,j) = BETA*(g-1.0d0)
      !ENGEN-VALUES AT INTERMEDIATE
      EIGM(1,j) = AKX*UM+AKY*VM
      EIGM(2,j) = EIGM(1,j)
      EIGM(3,j) = EIGM(1,j)+CM*dsqrt(AKX**2+AKY**2)
      EIGM(4,j) = EIGM(1,j)-CM*dsqrt(AKX**2+AKY**2)
      !ALPHA
      do k = 1,4
        D(k) = (Q(i,j+1,k)*Jac(i,j+1)-Q(i,j,k)*Jac(i,j))/AJACM
      end do
      do k = 1,4
        ALPHA(k,j) = RJ(k,1,j)*D(1)+RJ(k,2,j)*D(2)+RJ(k,3,j)*D(3)+RJ(k,4,j)*D(4)
      end do
    end do
    do k = 1,4
      ALPHA(k,-1) = ALPHA(k,0)
      ALPHA(k,JM) = ALPHA(k,JM-1)
    end do

    do k = 1,4
      do j = 0,JM
        !MINMOD LIMITER
        S = DSIGN(1.0d0,ALPHA(k,j))
        GG(k,j) = S*DMAX1(0.0d0,DMIN1(abs(ALPHA(k,j)),S*ALPHA(k,j-1)))
      end do
    end do
    do k = 1,4
      do j = 0,JM-1
        !GAMMA
        DELTA = 1.0d-10
        !DELTA=AMAX1(0.0,EIGM(K,J)-EIG(K,J),EIG(K,J+1)-EIGM(K,J))
        if (ALPHA(k,j)/=0.0) then
          GAMMA = 0.5d0*FPSI(EIGM(k,j),DELTA)*(GG(k,j+1)-GG(k,j))/ALPHA(k,j)
        else
          GAMMA = 0.0d0
        end if
        !PHI
        PHIM(k,j) = 0.5d0*FPSI(EIGM(k,j),DELTA)*(GG(k,j+1)+GG(k,j))-FPSI(EIGM(k,j)+GAMMA,DELTA)*ALPHA(k,j)
      end do
    end do
    !CONVECTION COMPORNENTS
    do j = 0,JM
      VCONT = eta_x(i,j)*u(i,j)+eta_y(i,j)*v(i,j)
      FH(1,j) = Q(i,j,1)*VCONT
      FH(2,j) = Q(i,j,2)*VCONT+p(i,j)*eta_x(i,j)/Jac(i,j)
      FH(3,j) = Q(i,j,3)*VCONT+p(i,j)*eta_y(i,j)/Jac(i,j)
      FH(4,j) = Q(i,j,4)*VCONT+p(i,j)*VCONT/Jac(i,j)
    end do
    !XI-DIRECTION CONVECTION FLUX
    do k = 1,4
      do j = 0,JM-1
        F(i,j,k) = 0.5d0*(FH(k,j)+FH(k,j+1)+RR(k,1,j)*PHIM(1,j)+RR(k,2,j)*PHIM(2,j) &
                          +RR(k,3,j)*PHIM(3,j)+RR(k,4,j)*PHIM(4,j))
      end do
    end do
  end do
end subroutine tvd_eta

!********************サブルーチン3　粘性項の計算*****************
subroutine viscosity_calc
  double precision :: txx,tyy,txy
  double precision :: mu0,T0,Sl
  double precision,dimension(0:IM,0:JM) :: u_xi,v_xi,T_xi,u_eta,v_eta,T_eta
  double precision,dimension(0:IM,0:JM) :: mu
  double precision :: Rv_s !Rvの値保持用
  double precision,parameter :: Pr = 0.72d0 !定数

  mu0 = 1.0d-3 !粘性係数の基準値
  T0 = 293.0d0 !基準温度
  Sl = 110.4d0 !サザーランド定数

  do i = 0,IM
    do j = 0,JM
      mu(i,j) = ((T(i,j)/T0)**1.5)*(T0+Sl)/(T(i,j)+Sl)*mu0 !粘性係数（サザーランドの式）
    end do
  end do
  do j = 0,JM
    do i = 0,IM
      if (i>0 .and. i<IM) then
        u_xi(i,j) = (u(i+1,j)-u(i-1,j))/2.0d0
        v_xi(i,j) = (v(i+1,j)-v(i-1,j))/2.0d0
        T_xi(i,j) = (T(i+1,j)-T(i-1,j))/2.0d0
      end if
      if (j>0 .and. j<JM) then
        u_eta(i,j) = (u(i,j+1)-u(i,j-1))/2.0d0
        v_eta(i,j) = (v(i,j+1)-v(i,j-1))/2.0d0
        T_eta(i,j) = (T(i,j+1)-T(i,j-1))/2.0d0
      end if
      if (i==0) then
        u_xi(i,j) = (-3.0d0*u(i,j)+4.0d0*u(i+1,j)-u(i+2,j))/2.0d0
        v_xi(i,j) = (-3.0d0*v(i,j)+4.0d0*v(i+1,j)-v(i+2,j))/2.0d0
        T_xi(i,j) = (-3.0d0*T(i,j)+4.0d0*T(i+1,j)-T(i+2,j))/2.0d0
      end if
      if (i==IM) then
        u_xi(i,j) = (3.0d0*u(i,j)-4.0d0*u(i-1,j)+u(i-2,j))/2.0d0
        v_xi(i,j) = (3.0d0*v(i,j)-4.0d0*v(i-1,j)+v(i-2,j))/2.0d0
        T_xi(i,j) = (3.0d0*T(i,j)-4.0d0*T(i-1,j)+T(i-2,j))/2.0d0
      end if
      if (j==0) then
        u_eta(i,j) = (-3.0d0*u(i,j)+4.0d0*u(i,j+1)-u(i,j+2))/2.0d0
        v_eta(i,j) = (-3.0d0*v(i,j)+4.0d0*v(i,j+1)-v(i,j+2))/2.0d0
        T_eta(i,j) = (-3.0d0*T(i,j)+4.0d0*T(i,j+1)-T(i,j+2))/2.0d0
      end if
      if (j==JM) then
        u_eta(i,j) = (3.0d0*u(i,j)-4.0d0*u(i,j-1)+u(i,j-2))/2.0d0
        v_eta(i,j) = (3.0d0*v(i,j)-4.0d0*v(i,j-1)+v(i,j-2))/2.0d0
        T_eta(i,j) = (3.0d0*T(i,j)-4.0d0*T(i,j-1)+T(i,j-2))/2.0d0
      end if
    end do
  end do
  do i = 0,IM !粘性項ベクトル
    do j = 0,JM
      txx = 2.0d0*mu(i,j)*(2.0d0*(u_xi(i,j)*xi_x(i,j)+u_eta(i,j)*eta_x(i,j))-(v_xi(i,j)*xi_y(i,j)+v_eta(i,j)*eta_y(i,j)))/3.0d0
      txy = mu(i,j)*((u_xi(i,j)*xi_y(i,j)+u_eta(i,j)*eta_y(i,j))+(v_xi(i,j)*xi_x(i,j)+v_eta(i,j)*eta_x(i,j)))
      tyy = 2.0d0*mu(i,j)*(2.0d0*(v_xi(i,j)*xi_y(i,j)+v_eta(i,j)*eta_y(i,j))-(u_xi(i,j)*xi_x(i,j)+u_eta(i,j)*eta_x(i,j)))/3.0d0
      Rv(2,i,j) = txx
      Rv(3,i,j) = txy
      Rv(4,i,j) = txx*u(i,j)+txy*v(i,j)+mu(i,j)*g*R*(T_xi(i,j)*xi_x(i,j)+T_eta(i,j)*eta_x(i,j))/(Pr*(g-1.0d0))
      Sv(2,i,j) = txy
      Sv(3,i,j) = tyy
      Sv(4,i,j) = txy*u(i,j)+tyy*v(i,j)+mu(i,j)*g*R*(T_xi(i,j)*xi_y(i,j)+T_eta(i,j)*eta_y(i,j))/(Pr*(g-1.0d0))
      do k = 1,4 !各流束を一般座標系に変換
        Rv_s = Rv(k,i,j)
        Rv(k,i,j) = (xi_x(i,j)*Rv(k,i,j)+xi_y(i,j)*Sv(k,i,j))/Jac(i,j)
        Sv(k,i,j) = (eta_x(i,j)*Rv_s+eta_y(i,j)*Sv(k,i,j))/Jac(i,j)
      end do
    end do
  end do
end subroutine viscosity_calc

!********************サブルーチン4　ルンゲクッタ********************
subroutine runge_kutta(is,ie,js,je)
  integer is,ie,js,je
  ! double precision us,ddu,vs,ddv
  ! du = 0.d0
  ! dv = 0.d0
  do j = js,je
    do i = is,ie
      do k = 1,4
        Q(i,j,k) = Q1(i,j,k)-1.0d0/(5.0d0-dble(k_runge_kutta))*dt*((E(i,j,k)-E(i-1,j,k)) &
                                                     +(F(i,j,k)-F(i,j-1,k)-(Rv(k,i+1,j)-Rv(k,i-1,j)+Sv(k,i,j+1)-Sv(k,i,j-1))/2.0d0))
        Q(i,j,k)=Q1(i,j,k)
      end do
      ! us = u(i,j)
      ! vs = v(i,j)
      ! rho(i,j) = Q(i,j,1)*Jac(i,j)
      ! u(i,j) = Q(i,j,2)*Jac(i,j)/rho(i,j)
      ! v(i,j) = Q(i,j,3)*Jac(i,j)/rho(i,j)
      ! U_cvc(i,j) = xi_x(i,j)*u(i,j)+xi_y(i,j)*v(i,j)
      ! V_cvc(i,j) = eta_x(i,j)*u(i,j)+eta_y(i,j)*v(i,j)
      ! energy(i,j) = Q(i,j,4)*Jac(i,j)
      ! T(i,j) = (g-1.0d0)*(energy(i,j)/rho(i,j)-(u(i,j)**2+v(i,j)**2)/2.0d0)/R
      ! p(i,j) = rho(i,j)*R*T(i,j)

      ! ddu = abs(u(i,j)-us)/us
      ! ddv = abs(v(i,j)-vs)/vs

      ! if (ddu>du) du = ddu
      ! if (ddv>dv) dv = ddv

    end do
  end do
  !保存量ベクトルQ^から諸量を求める
  do i = 1,IM-1
    do j = 1,JM-1
      rho(i,j) = Q(i,j,1)*Jac(i,j)
      u(i,j) = Q(i,j,2)*Jac(i,j)/rho(i,j)
      v(i,j) = Q(i,j,3)*Jac(i,j)/rho(i,j)
      energy(i,j) = Q(i,j,4)*Jac(i,j)
      p(i,j) = (g-1.0d0)*(energy(i,j)-(u(i,j)**2.0+v(i,j)**2.0)*rho(i,j)/2.0d0)
      T(i,j) = p(i,j)/(rho(i,j)*R)
    end do
  end do
end subroutine runge_kutta

!サブルーチン6　出力
subroutine Save_data
  open (10,file='a6.dat',status='replace')
  do j = 0,JM
    do i = 0,IM
      write (10,'(f13.6,$)') x(i,j)
      write (10,'(f13.6,$)') y(i,j)
      write (10,'(f13.6,$)') u(i,j)
      write (10,'(f13.6,$)') v(i,j)
      write (10,'(f13.3,$)') p(i,j)
      write (10,'(f13.6)') T(i,j)
    end do
  end do
  close (10)

  open (11,file='a6.fld',status='replace')
  write (11,'(a)') '# AVS field file'
  write (11,'(a)') 'ndim=2'
  write (11,'(a,i0)') 'dim1=',IM+1
  write (11,'(a,i0)') 'dim2=',JM+1
  write (11,'(a)') 'nspace=2'
  write (11,'(a)') 'veclen=4'
  write (11,'(a)') 'data=float'
  write (11,'(a)') 'field=irregular'
  write (11,'(a)') 'label=u,v,p,T'
  write (11,'(a)') 'variable 1 file=a6.dat filetype=ascii skip=0 offset=2 stride=6'
  write (11,'(a)') 'variable 2 file=a6.dat filetype=ascii skip=0 offset=3 stride=6'
  write (11,'(a)') 'variable 3 file=a6.dat filetype=ascii skip=0 offset=4 stride=6'
  write (11,'(a)') 'variable 4 file=a6.dat filetype=ascii skip=0 offset=5 stride=6'
  write (11,'(a)') 'coord 1 file=a6.dat filetype=ascii skip=0 offset=0 stride=6'
  write (11,'(a)') 'coord 2 file=a6.dat filetype=ascii skip=0 offset=1 stride=6'
  close (11)
  print*,"→ MicroAVS用の .dat および .fld を出力しました。"
end subroutine Save_data

end program a6
