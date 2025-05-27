program comp3

!変数の宣言
implicit none
double precision,parameter::H=10.d0 !基準長さ（高さ）
integer,parameter::IM=30 !x方向領域数
integer,parameter::JM=30 !y方向領域数
double precision,parameter::dt=1.0d-6

double precision,dimension(0:IM,0:JM)::x=0 !格子点(i,j)におけるx座標
double precision,dimension(0:IM,0:JM)::y=0 !格子点(i,j)におけるy座標
double precision,parameter::cr=1.1 !等比数列公比
double precision::a !等比数列初項

double precision,dimension(0:IM,0:JM)::xxi=0,xeta=0,yxi=0,yeta=0 !ξ→xi，η→eta
double precision,dimension(0:IM,0:JM)::JJ !ヤコビアン
double precision,dimension(0:IM,0:JM)::xix=0,xiy=0,etax=0,etay=0

integer i,j,n,k
integer,parameter::nmax=1000
double precision,parameter::EPS=1.0d-5
double precision qp,dq,ddq

double precision,parameter::R=287.1 !ガス定数
double precision,parameter::g=1.4 !比熱比

double precision::M=2.9 !マッハ数
double precision::T=293.0 !温度

double precision,dimension(0:IM,0:JM)::rho=0 !密度
double precision,dimension(0:IM,0:JM)::u=0 !x方向速度
double precision,dimension(0:IM,0:JM)::v=0 !y方向速度
double precision,dimension(0:IM,0:JM)::p=0 !圧力
double precision,dimension(0:IM,0:JM)::ene=0 !エネルギー

double precision,dimension(0:IM,0:JM,1:4)::Q=0 !保存量ベクトル
double precision,dimension(0:IM,0:JM,1:4)::E=0 !x方向流束ベクトル
double precision,dimension(0:IM,0:JM,1:4)::F=0 !y方向流束ベクトル


!物理面座標の設定
a=H*(1-cr)/(1-cr**JM) !等比数列初項
do i=0,IM
 do j=0,JM
  x(i,j)=dble(i)*3.d0*H/dble(IM)
  y(i,j)=a*(cr**j-1)/(cr-1)
 end do
end do


!座標変換
do i=0,IM
 do j=0,JM
  if (i>0.and.i<IM) then
   xxi(i,j)=(x(i+1,j)-x(i-1,j))/2.d0
   yxi(i,j)=(y(i+1,j)-y(i-1,j))/2.d0
  end if
  if(j>0.and.j<JM) then
   xeta(i,j)=(x(i,j+1)-x(i,j-1))/2.d0
   yeta(i,j)=(y(i,j+1)-y(i,j-1))/2.d0
  end if
  if(i==0) then
   xxi(i,j)=(-3.d0*x(i,j)+4.d0*x(i+1,j)-x(i+2,j))/2.d0
   yxi(i,j)=(-3.d0*y(i,j)+4.d0*y(i+1,j)-y(i+2,j))/2.d0
  end if
  if(i==IM) then
   xxi(i,j)=(3.d0*x(i,j)-4.d0*x(i-1,j)+x(i-2,j))/2.d0
   yxi(i,j)=(3.d0*y(i,j)-4.d0*y(i-1,j)+y(i-2,j))/2.d0
  end if
  if(j==0) then
   xeta(i,j)=(-3.d0*x(i,j)+4.d0*x(i,j+1)-x(i,j+2))/2.d0
   yeta(i,j)=(-3.d0*y(i,j)+4.d0*y(i,j+1)-y(i,j+2))/2.d0
  end if
  if(j==JM) then
   xeta(i,j)=(3.d0*x(i,j)-4.d0*x(i,j-1)+x(i,j-2))/2.d0
   yeta(i,j)=(3.d0*y(i,j)-4.d0*y(i,j-1)+y(i,j-2))/2.d0
  end if

  JJ(i,j)=1.d0/(xxi(i,j)*yeta(i,j)-xeta(i,j)*yxi(i,j))
  xix(i,j)=JJ(i,j)*yeta(i,j)
  xiy(i,j)=-JJ(i,j)*xeta(i,j)
  etax(i,j)=-JJ(i,j)*yxi(i,j)
  etay(i,j)=JJ(i,j)*xxi(i,j)

 end do
end do


!初期条件
do i=0,IM
 do j=0,JM
  rho(i,j)=1.2
  u(i,j)=M*sqrt(g*R*T)
  p(i,j)=rho(i,j)*R*T
  ene(i,j)=rho(i,j)*(R*T/(g-1.d0)+(u(i,j)**2.d0+v(i,j)**2.d0)/2.d0)

  Q(i,j,1)=rho(i,j)
  Q(i,j,2)=rho(i,j)*u(i,j)
  Q(i,j,3)=rho(i,j)*v(i,j)
  Q(i,j,4)=ene(i,j)
  E(i,j,1)=rho(i,j)*u(i,j)
  E(i,j,2)=p(i,j)+rho(i,j)*u(i,j)**2.d0
  E(i,j,3)=rho(i,j)*u(i,j)*v(i,j)
  E(i,j,4)=(ene(i,j)+p(i,j))*u(i,j)
  F(i,j,1)=rho(i,j)*v(i,j)
  F(i,j,2)=rho(i,j)*u(i,j)*v(i,j)
  F(i,j,3)=p(i,j)+rho(i,j)*v(i,j)**2.d0
  F(i,j,4)=(ene(i,j)+p(i,j))*v(i,j)

  do k=1,4 !各流束を一般座標系に変換
   Q(i,j,k)=Q(i,j,k)/JJ(i,j)
   E(i,j,k)=(xix(i,j)*E(i,j,k)+xiy(i,j)*F(i,j,k))/JJ(i,j)
   F(i,j,k)=(etax(i,j)*E(i,j,k)+etay(i,j)*F(i,j,k))/JJ(i,j)
  end do

 end do
end do


!計算
n=0
dq=EPS+1.d0
do while(dq>EPS.and.n<nmax)
 n=n+1
 dq=0.d0
 !境界条件の更新
 do i=1,IM !上下端
  do k=1,4
   Q(i,0,k)=2*Q(i,1,k)-Q(i,2,k)
   Q(i,JM,k)=2*Q(i,JM-1,k)-Q(i,JM-2,k)
   E(i,0,k)=2*E(i,1,k)-E(i,2,k)
   E(i,JM,k)=2*E(i,JM-1,k)-E(i,JM-2,k)
   F(i,0,k)=2*F(i,1,k)-F(i,2,k)
   F(i,JM,k)=2*F(i,JM-1,k)-F(i,JM-2,k)
  end do
 end do
 do j=0,JM !流出
  do k=1,4
   Q(IM,j,k)=2*Q(IM-1,j,k)-Q(IM-2,j,k)
   E(IM,j,k)=2*E(IM-1,j,k)-E(IM-2,j,k)
   F(IM,j,k)=2*F(IM-1,j,k)-F(IM-2,j,k)
  end do
 end do
 !オイラー方程式による計算
 do i=0,IM
  do j=0,JM
   qp=0
   if(i>0.and.i<IM.and.j>1.and.j<JM) then
    do k=1,4
     qp=qp+Q(i,j,k)
     Q(i,j,k)=Q(i,j,k)+dt*((E(i+1,j,k)-E(i-1,j,k))/(2.d0)+(F(i,j+1,k)-F(i,j-1,k))/(2.d0))
    end do
    ddq=abs((Q(i,j,1)+Q(i,j,2)+Q(i,j,3)+Q(i,j,4)-qp)/qp)
    if(ddq>dq) dq=ddq
   end if
   !保存量ベクトルQから各値を求める
   rho(i,j)=JJ(i,j)*Q(i,j,1)
   u(i,j)=JJ(i,j)*Q(i,j,2)/rho(i,j)
   v(i,j)=JJ(i,j)*Q(i,j,3)/rho(i,j)
   ene(i,j)=JJ(i,j)*Q(i,j,4)
   p(i,j)=ene(i,j)-(u(i,j)**2.d0+v(i,j)**2.d0)*rho(i,j)/2.d0
   !流束ベクトルの更新
    E(i,j,1)=rho(i,j)*u(i,j)
    E(i,j,2)=p(i,j)+rho(i,j)*u(i,j)**2.d0
    E(i,j,3)=rho(i,j)*u(i,j)*v(i,j)
    E(i,j,4)=(ene(i,j)+p(i,j))*u(i,j)
    F(i,j,1)=rho(i,j)*v(i,j)
    F(i,j,2)=rho(i,j)*u(i,j)*v(i,j)
    F(i,j,3)=p(i,j)+rho(i,j)*v(i,j)**2.d0
    F(i,j,4)=(ene(i,j)+p(i,j))*v(i,j)
    do k=1,4
     E(i,j,k)=(xix(i,j)*E(i,j,k)+xiy(i,j)*F(i,j,k))/JJ(i,j)
     F(i,j,k)=(etax(i,j)*E(i,j,k)+etay(i,j)*F(i,j,k))/JJ(i,j)
    end do
  end do
 end do
 print*,'Count',n,'dq = ',dq
end do


!.datファイル書き出し
open(1,file='comp3.dat',status='replace')
do j=0,JM
 do i=0,IM
  write(1,*) x(i,j),y(j,j),rho(i,j),u(i,j),v(i,j),p(i,j),ene(i,j)
 end do
end do
close(1)


!.fldファイル作成
open(2,file='comp3.fld',status='replace')
write(2,'(A)') '# AVS field file'
write(2,'(A)') 'ndim = 2'
write(2,'(A)',advance='no') 'dim1 ='
write(2,'(I5)') IM+1
write(2,'(A)',advance='no') 'dim2 ='
write(2,'(I5)') JM+1
write(2,'(A)') 'nspace = 2'
write(2,'(A)') 'veclen = 5'
write(2,'(A)') 'data = double'
write(2,'(A)') 'field = irregular'
write(2,'(A)') 'label = rho u v p ene'
write(2,'(A)') 'variable 1 file=comp3.dat filetype=ascii skip=0 offset=2 stride=7'
write(2,'(A)') 'variable 2 file=comp3.dat filetype=ascii skip=0 offset=3 stride=7'
write(2,'(A)') 'variable 3 file=comp3.dat filetype=ascii skip=0 offset=4 stride=7'
write(2,'(A)') 'variable 4 file=comp3.dat filetype=ascii skip=0 offset=5 stride=7'
write(2,'(A)') 'variable 5 file=comp3.dat filetype=ascii skip=0 offset=6 stride=7'
write(2,'(A)') 'coord 1 file=comp3.dat filetype=ascii skip=0 offset=0 stride=7'
write(2,'(A)') 'coord 2 file=comp3.dat filetype=ascii skip=0 offset=1 stride=7'
close(2)


end program comp3
