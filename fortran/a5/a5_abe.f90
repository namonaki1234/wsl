program comp5

!変数の定義*************************************************
implicit none

double precision,parameter::pi=dacos(0.0d0)*2 !π

double precision,parameter::H=10.0d0 !基準長さ（高さ）
integer,parameter::IM=30 !x方向領域数
integer,parameter::JM=30 !y方向領域数
double precision,parameter::dt=1.0d-6 !時刻刻み幅

double precision,dimension(0:IM,0:JM)::x=0.0d0 !格子点(i,j)におけるx座標
double precision,dimension(0:IM,0:JM)::y=0.0d0 !格子点(i,j)におけるy座標
double precision,parameter::cr=1.1d0 !等比数列公比
double precision::ft !等比数列初項

double precision,dimension(0:IM,0:JM)::xxi=0.0d0,xeta=0.0d0,yxi=0.0d0,yeta=0.0d0 !ξ→xi，η→eta
double precision,dimension(0:IM,0:JM)::JJ !ヤコビアン
double precision,dimension(0:IM,0:JM)::xix=0.0d0,xiy=0.0d0,etax=0.0d0,etay=0.0d0

integer i,j,n,k,nn
integer,parameter::nmax=5000
double precision,parameter::eps=1.0d-6
double precision dq,ddq
double precision,dimension(0:IM,0:JM)::qp=1.0d0

double precision,parameter::R=287.1d0 !ガス定数
double precision,parameter::g=1.4d0 !比熱比

double precision::M=2.9d0 !衝撃波前マッハ数
double precision::b !衝撃波角β
double precision::AA !係数
double precision::un !速度uの値保持用

double precision,dimension(0:IM,0:JM)::rho=0.0d0 !密度
double precision,dimension(0:IM,0:JM)::u=0.0d0 !x方向速度
double precision,dimension(0:IM,0:JM)::v=0.0d0 !y方向速度
double precision,dimension(0:IM,0:JM)::p=0.0d0 !圧力
double precision,dimension(0:IM,0:JM)::T=0.0d0 !温度
double precision,dimension(0:IM,0:JM)::ene=0.0d0 !エネルギー

double precision,dimension(1:4,0:IM,0:JM)::Q=0.0d0 !保存量ベクトル
double precision,dimension(1:4,0:IM,0:JM)::E=0.0d0 !x方向流束ベクトル
double precision::Eh
double precision,dimension(1:4,0:IM,0:JM)::F=0.0d0 !y方向流束ベクトル
double precision,dimension(1:4,0:IM,0:JM,0:4)::QQ=0.0d0 !保存量ベクトル（ルンゲクッタ）
double precision,dimension(1:4,0:IM,0:JM)::EE=0.0d0 !x方向流束ベクトル（TVDスキーム）E~
double precision,dimension(1:4,0:IM,0:JM)::FF=0.0d0 !y方向流束ベクトル（TVDスキーム）F~

double precision kx,ky,JJm,um,vm,dum,cm,kxt,kyt,phi,thetat,beta,S,delta,gamma
double precision,dimension(1:4,1:4,0:IM)::RR,RI !行列
double precision,dimension(1:4,0:IM)::a !固有値
double precision,dimension(1:4,0:IM)::gg !制限関数
double precision,dimension(1:4,-1:IM)::alpha !α
double precision,dimension(1:4)::D !α計算用
double precision,dimension(1:4,0:IM)::phim !Φ


!物理面座標の設定*******************************************
ft=H*(1.0d0-cr)/(1.0d0-cr**real(JM)) !等比数列初項
do i=0,IM
 do j=0,JM
  x(i,j)=dble(i)*3.0d0*H/dble(IM)
  y(i,j)=ft*(cr**real(j)-1.0d0)/(cr-1.0d0)
 end do
end do


!座標変換***************************************************
do i=0,IM
 do j=0,JM
  if (i>0.and.i<IM) then
   xxi(i,j)=(x(i+1,j)-x(i-1,j))/2.0d0
   yxi(i,j)=(y(i+1,j)-y(i-1,j))/2.0d0
  end if
  if(j>0.and.j<JM) then
   xeta(i,j)=(x(i,j+1)-x(i,j-1))/2.0d0
   yeta(i,j)=(y(i,j+1)-y(i,j-1))/2.0d0
  end if
  if(i==0) then
   xxi(i,j)=(-3.0d0*x(i,j)+4.0d0*x(i+1,j)-x(i+2,j))/2.0d0
   yxi(i,j)=(-3.0d0*y(i,j)+4.0d0*y(i+1,j)-y(i+2,j))/2.0d0
  end if
  if(i==IM) then
   xxi(i,j)=(3.0d0*x(i,j)-4.0d0*x(i-1,j)+x(i-2,j))/2.0d0
   yxi(i,j)=(3.0d0*y(i,j)-4.0d0*y(i-1,j)+y(i-2,j))/2.0d0
  end if
  if(j==0) then
   xeta(i,j)=(-3.0d0*x(i,j)+4.0d0*x(i,j+1)-x(i,j+2))/2.0d0
   yeta(i,j)=(-3.0d0*y(i,j)+4.0d0*y(i,j+1)-y(i,j+2))/2.0d0
  end if
  if(j==JM) then
   xeta(i,j)=(3.0d0*x(i,j)-4.0d0*x(i,j-1)+x(i,j-2))/2.0d0
   yeta(i,j)=(3.0d0*y(i,j)-4.0d0*y(i,j-1)+y(i,j-2))/2.0d0
  end if

  JJ(i,j)=1.0d0/(xxi(i,j)*yeta(i,j)-xeta(i,j)*yxi(i,j))
  xix(i,j)=JJ(i,j)*yeta(i,j)
  xiy(i,j)=-JJ(i,j)*xeta(i,j)
  etax(i,j)=-JJ(i,j)*yxi(i,j)
  etay(i,j)=JJ(i,j)*xxi(i,j)

 end do
end do


!初期条件***************************************************
b=pi/6.0d0 !衝撃波角β
AA=((g+1.0d0)*(M*dsin(b))**2.0)/((g-1.0d0)*(M*dsin(b))**2.0+2.0d0) !係数
do i=0,IM
 do j=0,JM
  rho(i,j)=1.2d0
  T(i,j)=293.0d0
  u(i,j)=M*dsqrt(g*R*T(i,j))
  p(i,j)=rho(i,j)*R*T(i,j)
  ene(i,j)=rho(i,j)*(R*T(i,j)/(g-1.0d0)+(u(i,j)**2.0+v(i,j)**2.0)/2.0d0)

  if(y(i,j)>=y(0,JM-2)-x(i,j)/dsqrt(3.0d0)) then !衝撃波後の物理量
   un=u(i,j)
   rho(i,j)=AA*rho(i,j)
   u(i,j)=u(i,j)*dsin(b)*dsqrt((1.0d0/AA)**2.0+(1.0d0/dtan(b))**2d0)*dcos(-b+datan(dtan(b)/AA))
   v(i,j)=un*dsin(b)*dsqrt((1.0d0/AA)**2.0+(1.0d0/dtan(b))**2.0)*dsin(-b+datan(dtan(b)/AA))
   p(i,j)=(2.0d0*g*(M*dsin(b))**2.0-(g-1.0d0))/(g+1.0d0)*p(i,j)
   T(i,j)=(2.0d0*g*(M*dsin(b))**2.0d0-(g-1.0d0))*((g-1.0d0)*(M*dsin(b))**2.0+2.0d0)/(((g+1.0d0)*M*dsin(b))**2.0)*T(i,j)
   ene(i,j)=rho(i,j)*(R*T(i,j)/(g-1.0d0)+(u(i,j)**2.0+v(i,j)**2.0)/2.0d0)
  end if

  Q(1,i,j)=rho(i,j)
  Q(2,i,j)=rho(i,j)*u(i,j)
  Q(3,i,j)=rho(i,j)*v(i,j)
  Q(4,i,j)=ene(i,j)
  E(1,i,j)=rho(i,j)*u(i,j)
  E(2,i,j)=p(i,j)+rho(i,j)*u(i,j)**2.0
  E(3,i,j)=rho(i,j)*u(i,j)*v(i,j)
  E(4,i,j)=(ene(i,j)+p(i,j))*u(i,j)
  F(1,i,j)=rho(i,j)*v(i,j)
  F(2,i,j)=rho(i,j)*u(i,j)*v(i,j)
  F(3,i,j)=p(i,j)+rho(i,j)*v(i,j)**2.0
  F(4,i,j)=(ene(i,j)+p(i,j))*v(i,j)

  do k=1,4 !各流束を一般座標系に変換
   Q(k,i,j)=Q(k,i,j)/JJ(i,j)
   QQ(k,i,j,0)=Q(k,i,j)
   Eh=E(k,i,j)
   E(k,i,j)=(xix(i,j)*E(k,i,j)+xiy(i,j)*F(k,i,j))/JJ(i,j)
   F(k,i,j)=(etax(i,j)*Eh+etay(i,j)*F(k,i,j))/JJ(i,j)
  end do

 end do
end do


!計算*******************************************************
n=0
dq=EPS+1.0d0
do while(dq>EPS.and.n<nmax)
 n=n+1
 dq=0.0d0
 do nn=1,4
  !ξ方向TVDスキーム----------------------------------------
  do j=1,JM-1
   do i=0,IM-1
    kx=0.5d0*(xix(i,j)+xix(i+1,j)) !中間点のk_x=ξ_x
    ky=0.5d0*(xiy(i,j)+xiy(i+1,j)) !中間点のk_y=ξ_y
    JJm=0.5d0*(JJ(i,j)+JJ(i+1,j)) !中間点のJJ
    um=(dsqrt(rho(i,j))*u(i,j)+dsqrt(rho(i+1,j))*u(i+1,j))/(dsqrt(rho(i,j))+dsqrt(rho(i+1,j))) !中間点のu
    vm=(dsqrt(rho(i,j))*v(i,j)+dsqrt(rho(i+1,j))*v(i+1,j))/(dsqrt(rho(i,j))+dsqrt(rho(i+1,j))) !中間点のv
    dum=(dsqrt(rho(i,j))*Q(4,i,j)/Q(1,i,j)+dsqrt(rho(i+1,j))*Q(4,i+1,j)/Q(1,i+1,j))&
        /(dsqrt(rho(i,j))+dsqrt(rho(i+1,j))) !中間点のe/ρ
    cm=sqrt(g*(g-1.0d0)*dabs(dum-0.5d0*(um**2.0+vm**2.0))) !中間点のc（音速）

    kxt=kx/dsqrt(kx**2.0+ky**2.0) !k_x~
    kyt=ky/dsqrt(kx**2.0+ky**2.0) !k_y~
    phi=0.5d0*(g-1.0d0)*(um**2.0+vm**2.0) !φ
    thetat=kxt*um+kyt*vm !θ~
    beta=1.0d0/(2.0d0*cm**2.0)

    !行列R（中間点）
    RR(1,1,i)=1.0d0
    RR(1,2,i)=0.0d0
    RR(1,3,i)=1.0d0
    RR(1,4,i)=1.0d0
    RR(2,1,i)=um
    RR(2,2,i)=kyt
    RR(2,3,i)=um+kxt*cm
    RR(2,4,i)=um-kxt*cm
    RR(3,1,i)=vm
    RR(3,2,i)=-kxt
    RR(3,3,i)=vm+kyt*cm
    RR(3,4,i)=vm-kyt*cm
    RR(4,1,i)=phi/(g-1.0d0)
    RR(4,2,i)=kyt*um-kxt*vm
    RR(4,3,i)=(phi+cm**2.0)/(g-1.0d0)+cm*thetat
    RR(4,4,i)=(phi+cm**2.0)/(g-1.0d0)-cm*thetat

    !Rの逆行列R^-1（中間点）
    RI(1,1,i)=1.0d0-phi/cm**2.0
    RI(1,2,i)=(g-1.0d0)*um/cm**2.0
    RI(1,3,i)=(g-1.0d0)*vm/cm**2.0
    RI(1,4,i)=-(g-1.0d0)/cm**2.0
    RI(2,1,i)=-kyt*um+kxt*vm
    RI(2,2,i)=kyt
    RI(2,3,i)=-kxt
    RI(2,4,i)=0.0d0
    RI(3,1,i)=beta*(phi-cm*thetat)
    RI(3,2,i)=beta*(kxt*cm-(g-1.0d0)*um)
    RI(3,3,i)=beta*(kyt*cm-(g-1.0d0)*vm)
    RI(3,4,i)=beta*(g-1.0d0)
    RI(4,1,i)=beta*(phi+cm*thetat)
    RI(4,2,i)=-beta*(kxt*cm+(g-1.0d0)*um)
    RI(4,3,i)=-beta*(kyt*cm+(g-1.0d0)*vm)
    RI(4,4,i)=beta*(g-1.0d0)

    !固有値a
    a(1,i)=kx*um+ky*vm
    a(2,i)=a(1,i)
    a(3,i)=a(1,i)+cm*dsqrt(kx**2.0+ky**2.0)
    a(4,i)=a(1,i)-cm*dsqrt(kx**2.0+ky**2.0)

    !α
    do k=1,4
     D(k)=(Q(k,i+1,j)*JJ(i+1,j)-Q(k,i,j)*JJ(i,j))/JJm
    end do
    do k=1,4
     alpha(k,i)=RI(k,1,i)*D(1)+RI(k,2,i)*D(2)+RI(k,3,i)*D(3)+RI(k,4,i)*D(4)
    end do
   end do
   do k=1,4
    alpha(k,-1)=alpha(k,0)
    alpha(k,IM)=alpha(k,IM-1)
   end do

  !制限関数g
   do k=1,4
    do i=0,IM
     S=dsign(1.0d0,alpha(k,i))
     gg(k,i)=S*dmax1(0.0d0,dmin1(dabs(alpha(k,i)),S*alpha(k,i-1)))
    end do
   end do

   do k=1,4
    do i=0,IM-1
     !γ（比熱比とは別）
     delta=1.0d-10
     if(alpha(k,i)/=0.0d0)then
      gamma=0.5d0*FPSI(a(k,i),delta)*(gg(k,i+1)-gg(k,i))/alpha(k,i)
     else
      gamma=0.0d0
     end if
     !Φ（φとは別）
     phim(k,i)=0.5d0*FPSI(a(k,i),delta)*(gg(k,i)+gg(k,i+1))-FPSI(a(k,i)+gamma,delta)*alpha(k,i)
    end do
   end do

   !数値流束E~
   do k=1,4
    do i=0,IM-1
      EE(k,i,j)=0.5d0*(E(k,i,j)+E(k,i+1,j)&
                +RR(k,1,i)*phim(1,i)+RR(k,2,i)*phim(2,i)&
                +RR(k,3,i)*phim(3,i)+RR(k,4,i)*phim(4,i))
    end do
   end do
  end do

  !η方向TVDスキーム-----------------------------------------
  do i=1,IM-1
   do j=0,JM-1
    kx=0.5d0*(etax(i,j)+etax(i,j+1)) !中間点のk_x=η_x
    ky=0.5d0*(etay(i,j)+etay(i,j+1)) !中間点のk_y=η_y
    JJm=0.5d0*(JJ(i,j)+JJ(i,j+1)) !中間点のJJ
    um=(dsqrt(rho(i,j))*u(i,j)+dsqrt(rho(i,j+1))*u(i,j+1))/(dsqrt(rho(i,j))+dsqrt(rho(i,j+1))) !中間点のu
    vm=(dsqrt(rho(i,j))*v(i,j)+dsqrt(rho(i,j+1))*v(i,j+1))/(dsqrt(rho(i,j))+dsqrt(rho(i,j+1))) !中間点のv
    dum=(dsqrt(rho(i,j))*Q(4,i,j)/Q(1,i,j)+dsqrt(rho(i,j+1))*Q(4,i,j+1)/Q(1,i,j+1))&
        /(dsqrt(rho(i,j))+dsqrt(rho(i,j+1))) !中間点のe/ρ
    cm=sqrt(g*(g-1.0d0)*dabs(dum-0.5d0*(um**2.0+vm**2.0))) !中間点のc（音速）

    kxt=kx/dsqrt(kx**2.0+ky**2.0) !k_x~
    kyt=ky/dsqrt(kx**2.0+ky**2.0) !k_y~
    phi=0.5d0*(g-1.0d0)*(um**2.0+vm**2.0) !φ
    thetat=kxt*um+kyt*vm !θ~
    beta=1.0d0/(2.0d0*cm**2.0)

    !行列R（中間点）
    RR(1,1,j)=1.0d0
    RR(1,2,j)=0.0d0
    RR(1,3,j)=1.0d0
    RR(1,4,j)=1.0d0
    RR(2,1,j)=um
    RR(2,2,j)=kyt
    RR(2,3,j)=um+kxt*cm
    RR(2,4,j)=um-kxt*cm
    RR(3,1,j)=vm
    RR(3,2,j)=-kxt
    RR(3,3,j)=vm+kyt*cm
    RR(3,4,j)=vm-kyt*cm
    RR(4,1,j)=phi/(g-1.0d0)
    RR(4,2,j)=kyt*um-kxt*vm
    RR(4,3,j)=(phi+cm**2.0)/(g-1.0d0)+cm*thetat
    RR(4,4,j)=(phi+cm**2.0)/(g-1.0d0)-cm*thetat

    !Rの逆行列R^-1（中間点）
    RI(1,1,j)=1.0d0-phi/cm**2.0
    RI(1,2,j)=(g-1.0d0)*um/cm**2.0
    RI(1,3,j)=(g-1.0d0)*vm/cm**2.0
    RI(1,4,j)=-(g-1.0d0)/cm**2.0
    RI(2,1,j)=-kyt*um+kxt*vm
    RI(2,2,j)=kyt
    RI(2,3,j)=-kxt
    RI(2,4,j)=0.0d0
    RI(3,1,j)=beta*(phi-cm*thetat)
    RI(3,2,j)=beta*(kxt*cm-(g-1.0d0)*um)
    RI(3,3,j)=beta*(kyt*cm-(g-1.0d0)*vm)
    RI(3,4,j)=beta*(g-1.0d0)
    RI(4,1,j)=beta*(phi+cm*thetat)
    RI(4,2,j)=-beta*(kxt*cm+(g-1.0d0)*um)
    RI(4,3,j)=-beta*(kyt*cm+(g-1.0d0)*vm)
    RI(4,4,j)=beta*(g-1.0d0)

    !固有値a
    a(1,j)=kx*um+ky*vm
    a(2,j)=a(1,j)
    a(3,j)=a(1,j)+cm*dsqrt(kx**2.0+ky**2.0)
    a(4,j)=a(1,j)-cm*dsqrt(kx**2.0+ky**2.0)

    !α
    do k=1,4
     D(k)=(Q(k,i,j+1)*JJ(i,j+1)-Q(k,i,j)*JJ(i,j))/JJm
    end do
    do k=1,4
     alpha(k,j)=RI(k,1,j)*D(1)+RI(k,2,j)*D(2)+RI(k,3,j)*D(3)+RI(k,4,j)*D(4)
    end do
   end do
   do k=1,4
    alpha(k,-1)=alpha(k,0)
    alpha(k,JM)=alpha(k,JM-1)
   end do

   !制限関数g
   do k=1,4
    do j=0,JM
     S=dsign(1.0d0,alpha(k,j))
     gg(k,j)=S*dmax1(0.0d0,dmin1(dabs(alpha(k,j)),S*alpha(k,j-1)))
    end do
   end do

   do k=1,4
    do j=0,JM-1
     !γ（比熱比とは別）
     delta=1.0d-10
     if(alpha(k,j)/=0.0d0)then
      gamma=0.5d0*FPSI(a(k,j),delta)*(gg(k,j+1)-gg(k,j))/alpha(k,j)
     else
      gamma=0.0d0
     end if

     !Φ（φとは別）
     phim(k,j)=0.5d0*FPSI(a(k,j),delta)*(gg(k,j)+gg(k,j+1))-FPSI(a(k,j)+gamma,delta)*alpha(k,j)
    end do
   end do

   !数値流束F~
   do k=1,4
    do j=0,JM-1
      FF(k,i,j)=0.5d0*(F(k,i,j)+F(k,i,j+1)&
                +RR(k,1,j)*phim(1,j)+RR(k,2,j)*phim(2,j)&
                +RR(k,3,j)*phim(3,j)+RR(k,4,j)*phim(4,j))
    end do
   end do
  end do

  !Runge-Kutta法--------------------------------------------
  do j=1,JM-1
   do i=1,IM-1
    do k=1,4
     QQ(k,i,j,nn)=QQ(k,i,j,0)-dt*(EE(k,i,j)-EE(k,i-1,j)+FF(k,i,j)-FF(k,i,j-1))/dble(5-nn)
     Q(k,i,j)=QQ(k,i,j,nn)
    end do
   end do
  end do

  !諸量の計算-----------------------------------------------
  !保存量ベクトルQ^から諸量を求める
  do i=1,IM-1
   do j=1,JM-1
    rho(i,j)=Q(1,i,j)*JJ(i,j)
    u(i,j)=Q(2,i,j)*JJ(i,j)/rho(i,j)
    v(i,j)=Q(3,i,j)*JJ(i,j)/rho(i,j)
    ene(i,j)=Q(4,i,j)*JJ(i,j)
    p(i,j)=(g-1.0d0)*(ene(i,j)-(u(i,j)**2.0+v(i,j)**2.0)*rho(i,j)/2.0d0)
    T(i,j)=p(i,j)/(rho(i,j)*R)
   end do
  end do
  do j=0,JM
   do i=0,IM
    !境界条件の更新
    if(i==IM)then !流出
     rho(IM,j)=2.0d0*rho(IM-1,j)-rho(IM-2,j)
     u(IM,j)=2.0d0*u(IM-1,j)-u(IM-2,j)
     v(IM,j)=2.0d0*v(IM-1,j)-v(IM-2,j)
     p(IM,j)=2.0d0*p(IM-1,j)-p(IM-2,j)
     T(IM,j)=p(IM,j)/(rho(IM,j)*R)
     ene(IM,j)=rho(IM,j)*(R*T(IM,j)/(g-1.0d0)+(u(IM,j)**2.0+v(IM,j)**2.0)/2.0d0)
    end if
    if(j==JM)then !上端
     rho(i,JM)=rho(i,JM-1)+(rho(i,JM-1)-rho(i,JM-2))*cr
     u(i,JM)=u(i,JM-1)+(u(i,JM-1)-u(i,JM-2))*cr
     v(i,JM)=v(i,JM-1)+(v(i,JM-1)-v(i,JM-2))*cr
     p(i,JM)=p(i,JM-1)+(p(i,JM-1)-p(i,JM-2))*cr
     T(i,JM)=p(i,JM)/(rho(i,JM)*R)
     ene(i,JM)=rho(i,JM)*(R*T(i,JM)/(g-1.0d0)+(u(i,JM)**2.0+v(i,JM)**2.0)/2.0d0)
    end if
    if(j==0)then !下端（壁）
     rho(i,0)=rho(i,1)+(rho(i,1)-rho(i,2))/cr
     u(i,0)=u(i,1)+(u(i,1)-u(i,2))/cr
     v(i,0)=0.0d0
     p(i,0)=p(i,1)+(p(i,1)-p(i,2))/cr
     T(i,0)=p(i,0)/(rho(i,0)*R)
     ene(i,0)=rho(i,0)*(R*T(i,0)/(g-1.0d0)+(u(i,0)**2.0+v(i,0)**2.0)/2.0d0)
    end if
    !流束ベクトルの更新
    Q(1,i,j)=rho(i,j)
    Q(2,i,j)=rho(i,j)*u(i,j)
    Q(3,i,j)=rho(i,j)*v(i,j)
    Q(4,i,j)=ene(i,j)
    E(1,i,j)=rho(i,j)*u(i,j)
    E(2,i,j)=p(i,j)+rho(i,j)*u(i,j)**2.0
    E(3,i,j)=rho(i,j)*u(i,j)*v(i,j)
    E(4,i,j)=(ene(i,j)+p(i,j))*u(i,j)
    F(1,i,j)=rho(i,j)*v(i,j)
    F(2,i,j)=rho(i,j)*u(i,j)*v(i,j)
    F(3,i,j)=p(i,j)+rho(i,j)*v(i,j)**2.0
    F(4,i,j)=(ene(i,j)+p(i,j))*v(i,j)

    do k=1,4 !各流束を一般座標系に変換
     Q(k,i,j)=Q(k,i,j)/JJ(i,j)
     Eh=E(k,i,j)
     E(k,i,j)=(xix(i,j)*E(k,i,j)+xiy(i,j)*F(k,i,j))/JJ(i,j)
     F(k,i,j)=(etax(i,j)*Eh+etay(i,j)*F(k,i,j))/JJ(i,j)
    end do

   end do
  end do

 end do

 !次のステップの初期値を設定/最大残差の計算
 do i=0,IM
  do j=0,JM
   do k=1,4
    QQ(k,i,j,0)=Q(k,i,j)
   end do
   ddq=abs((Q(1,j,k)+Q(2,j,k)+Q(3,j,k)+Q(4,j,k)-qp(i,j))/qp(i,j))
   qp(i,j)=Q(1,j,k)+Q(2,j,k)+Q(3,j,k)+Q(4,j,k)
   if(ddq>dq) dq=ddq
  end do
 end do

 print*,'count =',n,'dq =',dq

end do


!.datファイル書き出し***************************************
open(1,file='comp5.dat',status='replace')
do j=0,JM
 do i=0,IM
  write(1,*) x(i,j),y(j,j),u(i,j),v(i,j),p(i,j),T(i,j)
 end do
end do
close(1)


!.fldファイル作成*******************************************
open(2,file='comp5.fld',status='replace')
write(2,'(A)') '# AVS field file'
write(2,'(A)') 'ndim = 2'
write(2,'(A)',advance='no') 'dim1 ='
write(2,'(I5)') IM+1
write(2,'(A)',advance='no') 'dim2 ='
write(2,'(I5)') JM+1
write(2,'(A)') 'nspace = 2'
write(2,'(A)') 'veclen = 4'
write(2,'(A)') 'data = double'
write(2,'(A)') 'field = irregular'
write(2,'(A)') 'label = u v p T'
write(2,'(A)') 'variable 1 file=comp5.dat filetype=ascii skip=0 offset=2 stride=6'
write(2,'(A)') 'variable 2 file=comp5.dat filetype=ascii skip=0 offset=3 stride=6'
write(2,'(A)') 'variable 3 file=comp5.dat filetype=ascii skip=0 offset=4 stride=6'
write(2,'(A)') 'variable 4 file=comp5.dat filetype=ascii skip=0 offset=5 stride=6'
write(2,'(A)') 'coord 1 file=comp5.dat filetype=ascii skip=0 offset=0 stride=6'
write(2,'(A)') 'coord 2 file=comp5.dat filetype=ascii skip=0 offset=1 stride=6'
close(2)


contains

!FPSI関数の定義*********************************************
double precision function FPSI(z,delta)
double precision z,delta
if(abs(z)>=delta)then
 FPSI=dabs(z)
else
 FPSI=0.5d0*(z**2.0+delta**2.0)/delta
end if
end function


end program comp5