program a3
implicit none
!格子条件
double precision,parameter::H=10.0d0
integer,parameter::IM=30,JM=30 !格子数
double precision,parameter::dx=3.0*H/30.0 !x方向の格子幅
double precision,parameter::r_y=1.1d0 !yの公比
double precision::a

!流れ場初期条件
double precision,dimension(0:IM,0:JM)::M=2.9d0,T=293.0d0,R=287.1d0,g=1.4d0,rho=1.2d0,dt=1.0d-6 !マッハ数,温度,ガス定数,比熱比,密度,時間刻み

!その他の条件
integer::i,j,n
double precision,dimension(0:IM,0:JM)::c,x,y,u,v,p,e,Q1,Q2,Q3,Q4,E1,E2,E3,E4,F1,F2,F3,F4,E_bar&
,xh,yh,xe,ye,Jn,xix,xiy,ex,ey,uu,vv

  ! === 格子生成 ===
do i = 0, IM
  do j = 0, JM
    x(i,j) = dx * dble(i)
  end do
end do

do j = 0, JM
  a = H * (r_y - 1.0d0) / (r_y**dble(JM) - 1.0d0)
  do i = 0, IM
    y(i,j) = a * (r_y**dble(j) - 1.0d0) / (r_y - 1.0d0)
  end do
end do

!流入条件設定
do i=0,IM
 do j=0,JM
 c(i,j)=dsqrt(g(i,j)*R(i,j)*T(i,j)) !音速
 u(0,j)=M(0,j)*c(0,j) !u=Mc
 u(0,j)=u(1,j)
 end do
end do

 !初期値以外の速度計算
 do i=2,IM
  do j=0,JM
   u(i,j)=M(i,j)*c(i,j)
   v(i,j)=0.0
  end do
 end do

 !圧力とエネルギーの計算
 do j=0,JM
  do i=0,IM
   p(i,j)=rho(i,j)*R(i,j)*T(i,j) !圧力p
   e(i,j)=rho(i,j)*((R(i,j)*T(i,j)/g(i,j)-1.0))+((u(i,j)**2+v(i,j)**2.0)/2.0d0) !全エネルギーe
  end do
 end do


!座標変換(境界以外)
do i=0,IM
 do j=1,JM-1
  xh(i,j)=(x(i,j+1)-x(i,j-1))/2.0d0
  yh(i,j)=(y(i,j+1)-y(i,j-1))/2.0d0
 end do
end do

do j=0,JM
 do i=1,IM-1
  xe(i,j)=(x(i+1,j)-x(i-1,j))/2.0d0
  ye(i,j)=(y(i+1,j)-y(i-1,j))/2.0d0
  end do
end do

!座標変換（境界）
do j=0,JM
 xe(0,j)=(-3.0d0*x(0,j)+4.0d0*x(1,j)-x(2,j))/2.0d0
 xe(IM,j)=(3.0d0*x(IM,j)-4.0d0*x(IM-1,j)+x(IM-2,j))/2.0d0
 ye(0,j)=(-3.0d0*y(0,j)+4.0d0*y(1,j)-y(2,j))/2.0d0
 ye(IM,j)=(3.0d0*y(IM,j)-4.0d0*y(IM-1,j)+y(IM-2,j))/2.0d0
end do

do i=0,IM
 xh(i,0)=(-3.0d0*x(i,0)+4.0d0*x(i,1)-x(i,2))/2.0d0
 xh(i,JM)=(3.0d0*x(i,JM)-4.0d0*x(i,JM-1)+x(i,JM-2))/2.0d0
 yh(i,0)=(-3.0d0*y(i,0)+4.0d0*y(i,1)-y(i,2))/2.0d0
 yh(i,JM)=-(-3.0d0*y(i,JM)+4.0d0*y(i,JM-1)-y(i,JM-2))/2.0d0
end do

!ヤコビアン計算
do i=0,IM
 do j=0,JM
Jn(i,j)=1.0d0/(xe(i,j)*yh(i,j)-ye(i,j)*xh(i,j))
xix(i,j)=Jn(i,j)*yh(i,j)
xiy(i,j)=Jn(i,j)*xh(i,j)
ex(i,j)=-Jn(i,j)*ye(i,j)
ey(i,j)=Jn(i,j)*xe(i,j)
 end do
end do

!メイン計算
do n=1,10000
 do i=0,IM
  do j=0,JM
   uu(i,j)=xix(i,j)*u(i,j)+xiy(i,j)*v(i,j) !反変速度
   vv(i,j)=ex(i,j)*u(i,j)+ey(i,j)*v(i,j)
   Q1(i,j)=rho(i,j)/Jn(i,j) !保存変数ベクトル
   Q2(i,j)=(rho(i,j)*u(i,j))/Jn(i,j)
   Q3(i,j)=(rho(i,j)*v(i,j))/Jn(i,j)
   Q4(i,j)=e(i,j)/Jn(i,j)
   E1(i,j)=(rho(i,j)*uu(i,j))/Jn(i,j) !x方向流束ベクトル
   E2(i,j)=(rho(i,j)*u(i,j)*uu(i,j)+xix(i,j)*p(i,j))/Jn(i,j)
   E3(i,j)=(rho(i,j)*v(i,j)*uu(i,j)+xiy(i,j)*p(i,j))/Jn(i,j)
   E4(i,j)=(e(i,j)+p(i,j))*uu(i,j)/Jn(i,j)
   F1(i,j)=rho(i,j)*vv(i,j)/Jn(i,j) !y方向流束ベクトル
   F2(i,j)=(rho(i,j)*u(i,j)*vv(i,j)+ex(i,j)*p(i,j))/Jn(i,j)
   F3(i,j)=(rho(i,j)*v(i,j)*vv(i,j)+ey(i,j)*p(i,j))/Jn(i,j)
   F4(i,j)=(e(i,j)+p(i,j))*vv(i,j)/Jn(i,j)
  end do
 end do

!境界条件以外における式(7.16)
 do i=1,IM-1
  do j=1,JM-1
   Q1(i,j)=Q1(i,j)-dt(i,j)*((E1(i+1,j)-E1(i-1,j))/(2.0d0)+(F1(i,j+1)-F1(i,j-1))/2.0d0)
   Q2(i,j)=Q2(i,j)-dt(i,j)*((E2(i+1,j)-E2(i-1,j))/(2.0d0)+(F2(i,j+1)-F2(i,j-1))/2.0d0)
   Q3(i,j)=Q3(i,j)-dt(i,j)*((E3(i+1,j)-E3(i-1,j))/(2.0d0)+(F3(i,j+1)-F3(i,j-1))/2.0d0)
   Q4(i,j)=Q4(i,j)-dt(i,j)*((E4(i+1,j)-E4(i-1,j))/(2.0d0)+(F4(i,j+1)-F4(i,j-1))/2.0d0)
  end do
 end do

 do j=1,JM-1
  do i=1,IM-1
   rho(i,j)=Q1(i,j)*Jn(i,j)
   u(i,j)=Q2(i,j)*Jn(i,j)/rho(i,j)
   v(i,j)=Q3(i,j)*Jn(i,j)/rho(i,j)
   e(i,j)=Q4(i,j)*Jn(i,j)
  end do
 end do

!境界の計算
do i=1,IM-1
!上壁
  rho(i,JM)=rho(i,JM-1)+(y(i,JM)-y(i,JM-1))/(y(i,JM-1)-y(i,JM-2))*(rho(i,JM-1)-rho(i,JM-2))
  u(i,JM)=u(i,JM-1)+(y(i,JM)-y(i,JM-1))/(y(i,JM-1)-y(i,JM-2))*(u(i,JM-1)-u(i,JM-2))
  v(i,JM)=v(i,JM-1)+(y(i,JM)-y(i,JM-1))/(y(i,JM-1)-y(i,JM-2))*(v(i,JM-1)-v(i,JM-2))
  e(i,JM)=e(i,JM-1)+(y(i,JM)-y(i,JM-1))/(y(i,JM-1)-y(i,JM-2))*(e(i,JM-1)-e(i,JM-2))
  p(i,JM)=p(i,JM-1)+(y(i,JM)-y(i,JM-1))/(y(i,JM-1)-y(i,JM-2))*(p(i,JM-1)-p(i,JM-2))
!下壁
  rho(i,0)=rho(i,1)-(y(i,1)-y(i,0))/(y(i,2)-y(i,1))*(rho(i,2)-rho(i,1))
  u(i,0)=u(i,1)-(y(i,1)-y(i,0))/(y(i,2)-y(i,1))*(u(i,2)-u(i,1))
  v(i,0)=v(i,1)-(y(i,1)-y(i,0))/(y(i,2)-y(i,1))*(v(i,2)-v(i,1))
  e(i,0)=e(i,1)-(y(i,1)-y(i,0))/(y(i,2)-y(i,1))*(e(i,2)-e(i,1))
  p(i,0)=p(i,1)-(y(i,1)-y(i,0))/(y(i,2)-y(i,1))*(p(i,2)-p(i,1))
end do

do j=0,JM
!流出
  rho(IM,j)=rho(IM-1,j)+(x(IM,j)-x(IM-1,j))/(x(IM-1,j)-x(IM-2,j))*(rho(IM-1,j)-rho(IM-2,j))
  u(IM,j)=u(IM-1,j)+(x(IM,j)-x(IM-1,j))/(x(IM-1,j)-x(IM-2,j))*(u(IM-1,j)-u(IM-2,j))
  v(IM,j)=v(IM-1,j)+(x(IM,j)-x(IM-1,j))/(x(IM-1,j)-x(IM-2,j))*(v(IM-1,j)-v(IM-2,j))
  e(IM,j)=e(IM-1,j)+(x(IM,j)-x(IM-1,j))/(x(IM-1,j)-x(IM-2,j))*(e(IM-1,j)-e(IM-2,j))
  p(IM,j)=p(IM-1,j)+(x(IM,j)-x(IM-1,j))/(x(IM-1,j)-x(IM-2,j))*(p(IM-1,j)-p(IM-2,j))
end do

!更新
 do j=0,JM
  do i=0,IM
   T(i,j)=E_bar(i,j)*(g(i,j)-1.0d0)/r(i,j)
   E_bar(i,j)=e(i,j)/rho(i,j)-(u(i,j)**2.0+v(i,j)**2.0)/2.0d0 !式(5.7)
   p(i,j)=(g(i,j)-1.0d0)*rho(i,j)*E_bar(i,j) !式(5.10)
  end do
 end do
 !print*,n
end do

! === データ出力（MicroAVS用 DAT）===
open(10,file='a3.dat',status='replace')
do j=0,JM
 do i=0,IM
  write(10,*) x(i,j),y(i,j),u(i,j),v(i,j)
 end do
end do
close(10)


  ! === MicroAVSのFLDヘッダ出力 ===
open (11, file='a3.fld', status='replace')
write (11, '(A)') '# AVS field file'
write (11, '(A)') 'ndim = 2'
write (11, '(A,I5)') 'dim1 =', IM + 1
write (11, '(A,I5)') 'dim2 =', JM + 1
write (11, '(A)') 'nspace = 2'
write (11, '(A)') 'veclen = 2'
write (11, '(A)') 'data = double'
write (11, '(A)') 'field = irregular'
write (11, '(A)') 'label = u v'
write (11, '(A)') 'variable 1 file=a3.dat filetype=ascii skip=0 offset=2 stride=4'
write (11, '(A)') 'variable 2 file=a3.dat filetype=ascii skip=0 offset=3 stride=4'
write (11, '(A)') 'coord 1 file=a3.dat filetype=ascii skip=0 offset=0 stride=4'
write (11, '(A)') 'coord 2 file=a3.dat filetype=ascii skip=0 offset=1 stride=4'
close (11)

print *, "→ MicroAVS用の .dat および .fld を出力しました。"


end program a3
