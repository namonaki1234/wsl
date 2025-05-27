program comp4

!変数の定義
implicit none

double precision,parameter::pi=3.14159265358979d0 !π

double precision,parameter::H=10.0d0 !基準長さ（高さ）
integer,parameter::IM=30 !x方向領域数
integer,parameter::JM=30 !y方向領域数

double precision,dimension(0:IM,0:JM)::x=0.0d0 !格子点(i,j)におけるx座標
double precision,dimension(0:IM,0:JM)::y=0.0d0 !格子点(i,j)におけるy座標
double precision,parameter::cr=1.1d0 !等比数列公比
double precision::a !等比数列初項

double precision,dimension(0:IM,0:JM)::xxi=0.0d0,xeta=0.0d0,yxi=0.0d0,yeta=0.0d0 !ξ→xi，η→eta
double precision,dimension(0:IM,0:JM)::JJ !ヤコビアン
double precision,dimension(0:IM,0:JM)::xix=0.0d0,xiy=0.0d0,etax=0.0d0,etay=0.0d0

integer i,j

double precision,parameter::R=287.1d0 !ガス定数
double precision,parameter::g=1.4d0 !比熱比

double precision::M=2.9d0 !衝撃波前マッハ数
double precision::b !衝撃波角β
double precision::alpha !係数
double precision::un !速度uの値保持用

double precision,dimension(0:IM,0:JM)::rho=0.0d0 !密度
double precision,dimension(0:IM,0:JM)::u=0.0d0 !x方向速度
double precision,dimension(0:IM,0:JM)::v=0.0d0 !y方向速度
double precision,dimension(0:IM,0:JM)::p=0.0d0 !圧力
double precision,dimension(0:IM,0:JM)::T=0.0d0 !温度
double precision,dimension(0:IM,0:JM)::ene=0.0d0 !エネルギー

!物理面座標の設定
a=H*(1.0d0-cr)/(1.0d0-cr**JM) !等比数列初項
do i=0,IM
 do j=0,JM
  x(i,j)=dble(i)*3.0d0*H/dble(IM)
  y(i,j)=a*(cr**real(j-1))/(cr-1.0d0)
 end do
end do


!座標変換
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


!初期条件
b=pi/6.0d0 !衝撃波角β
alpha=((g+1.0d0)*(M*dsin(b))**2.0)/((g-1.0d0)*(M*dsin(b))**2.0+2.0d0) !係数
do i=0,IM
 do j=0,JM
  rho(i,j)=1.2d0
  T(i,j)=293.0d0
  u(i,j)=M*dsqrt(g*R*T(i,j))
  p(i,j)=rho(i,j)*R*T(i,j)
  ene(i,j)=rho(i,j)*(R*T(i,j)/(g-10.d0)+(u(i,j)**2.0+v(i,j)**2.0)/2.0d0)

  if(y(i,j)>=y(0,JM-2)-x(i,j)/dsqrt(3.0d0)) then !衝撃波後の物理量
   un=u(i,j)
   rho(i,j)=alpha*rho(i,j)
   u(i,j)=u(i,j)*dsin(b)*dsqrt((1.0d0/alpha)**2.0+(1.0d0/dtan(b))**2d0)*dcos(-b+datan(dtan(b)/alpha))
   v(i,j)=un*dsin(b)*dsqrt((1.0d0/alpha)**2.0+(1.0d0/dtan(b))**2.0)*dsin(-b+datan(dtan(b)/alpha))
   p(i,j)=(2.0d0*g*(M*dsin(b))**2.0-(g-1.0d0))/(g+1.0d0)*p(i,j)
   T(i,j)=(2.0d0*g*(M*dsin(b))**2.0d0-(g-1.0d0))*((g-1.0d0)*(M*dsin(b))**2.0+2.0d0)/(((g+1.0d0)*M*dsin(b))**2.0)*T(i,j)
   ene(i,j)=rho(i,j)*(R*T(i,j)/(g-1.0d0)+(u(i,j)**2.0+v(i,j)**2.0)/2.0d0)
  end if

 end do
end do


!.datファイル書き出し
open(1,file='comp4.dat',status='replace')
do j=0,JM
 do i=0,IM
  write(1,*) x(i,j),y(j,j),rho(i,j),u(i,j),v(i,j),p(i,j),ene(i,j)
 end do
end do
close(1)


!.fldファイル作成
open(2,file='comp4.fld',status='replace')
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
write(2,'(A)') 'variable 1 file=comp4.dat filetype=ascii skip=0 offset=2 stride=7'
write(2,'(A)') 'variable 2 file=comp4.dat filetype=ascii skip=0 offset=3 stride=7'
write(2,'(A)') 'variable 3 file=comp4.dat filetype=ascii skip=0 offset=4 stride=7'
write(2,'(A)') 'variable 4 file=comp4.dat filetype=ascii skip=0 offset=5 stride=7'
write(2,'(A)') 'variable 5 file=comp4.dat filetype=ascii skip=0 offset=6 stride=7'
write(2,'(A)') 'coord 1 file=comp4.dat filetype=ascii skip=0 offset=0 stride=7'
write(2,'(A)') 'coord 2 file=comp4.dat filetype=ascii skip=0 offset=1 stride=7'
close(2)


end program comp4