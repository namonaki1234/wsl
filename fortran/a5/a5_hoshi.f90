program kadai05

implicit none

integer::i,j,k,kk,mm,n
!********************格子空間・計算条件・初期条件********************
double precision,parameter		::H=10.0d0		!格子数
double precision,parameter		::XM=3.0d0*H		!長さ_x方向
double precision,parameter		::YM=H			!長さ_y方向
double precision			::du,dv,ddu,ddv,jy
integer,parameter			::IM=30			!終点_x方向
integer,parameter			::JM=30			!終点_y方向
double precision,parameter		::a=1.2d0
double precision,parameter		::b=(a+1.0d0)/(a-1.0d0)
double precision,parameter		::dx=XM/dble(IM)
double precision,parameter		::dy=YM/dble(JM)
double precision,parameter		::dxi=1.0d0,deta=1.0d0
double precision,parameter		::dt=1.0d-6
double precision,parameter		::gamma_h=1.4d0
double precision,parameter		::R=287.1d0
double precision,parameter		::pi=dacos(-1.0d0)
double precision,parameter		::beta_rad=pi/6.0d0

!********************配列定義********************
!-------衝撃波-------
double precision,dimension(0:IM,0:JM)	::T=0.0d0
double precision,dimension(0:IM,0:JM)	::M=2.9d0
double precision,dimension(0:IM,0:JM)	::rho=0.0d0
double precision,dimension(0:IM,0:JM)	::rho2=0.0d0
double precision,dimension(0:IM,0:JM)	::u=0.0d0
double precision,dimension(0:IM,0:JM)	::v=0.0d0
double precision,dimension(0:IM,0:JM)	::p=0.0d0
double precision,dimension(0:IM,0:JM)	::ee=0.0d0
double precision,dimension(0:IM,0:JM)	::bj=0.0d0
double precision,dimension(0:IM,0:JM,4)	::q=0.0d0
double precision,dimension(0:IM,0:JM,4)	::q1=0.0d0
double precision,dimension(0:IM,0:JM,4)	::e=0.0d0
double precision,dimension(0:IM,0:JM,4)	::e1=0.0d0
double precision,dimension(0:IM,0:JM,4)	::f=0.0d0
double precision,dimension(0:IM,0:JM,4)	::f1=0.0d0
double precision,dimension(0:IM,0:JM)	::u2=0.0d0
double precision,dimension(0:IM,0:JM)	::v2=0.0d0
double precision,dimension(0:IM,0:JM)	::T2=0.0d0
double precision,dimension(0:IM,0:JM)	::p2=0.0d0
double precision,dimension(0:IM,0:JM)	::Vt=0.0d0
double precision,dimension(0:IM,0:JM)	::Vn1=0.0d0
double precision,dimension(0:IM,0:JM)	::Vn2=0.0d0
double precision,dimension(0:IM,0:JM)	::lv2=0.0d0
double precision,dimension(0:IM,0:JM)	::theta=0.0d0
!-------ヤコビアンJa-------
double precision,dimension(0:IM,0:JM)	::x=0.0d0
double precision,dimension(0:IM,0:JM)	::y=0.0d0
double precision,dimension(0:IM,0:JM)	::xxi=0.0d0
double precision,dimension(0:IM,0:JM)	::xeta=0.0d0
double precision,dimension(0:IM,0:JM)	::yxi=0.0d0
double precision,dimension(0:IM,0:JM)	::yeta=0.0d0
double precision,dimension(0:IM,0:JM)	::xix=0.0d0
double precision,dimension(0:IM,0:JM)	::xiy=0.0d0
double precision,dimension(0:IM,0:JM)	::etax=0.0d0
double precision,dimension(0:IM,0:JM)	::etay=0.0d0
double precision,dimension(0:IM,0:JM)	::bu=0.0d0
double precision,dimension(0:IM,0:JM)	::bv=0.0d0
double precision,parameter::EPS=1.0d-6
integer,parameter::NMAX=200000

!********************不等間隔格子の生成・ファイル出力********************
do i=0,IM
	do j=0,JM
	x(i,j)=XM/dble(IM)*dble(i)
	end do
end do

do j=0,JM
	do i=0,IM
	jy=dble(j)/dble(JM)
	y(i,j)=YM*(b**jy-1.0d0)/(b-1.0d0)
	end do
end do

open(2,file='kadai05x.dat', status='replace')
do i=0,IM
write(2,*) x(i,0)
end do
close(2)

open(3,file='kadai05y.dat', status='replace')
do j=0,JM
write(3,*) y(0,j)
end do
close(3)


!********************一般座標へ変換・初期値設定********************
do j=0,JM
	do i=1,IM-1
	xxi(i,j)=(x(i+1,j)-x(i-1,j))/2.0d0
	yxi(i,j)=(y(i+1,j)-y(i-1,j))/2.0d0
	end do
end do

do i=0,IM
	do j=1,JM-1
	xeta(i,j)=(x(i,j+1)-x(i,j-1))/2.0d0
	yeta(i,j)=(y(i,j+1)-y(i,j-1))/2.0d0
	end do
end do

do i=0,IM
xeta(i,0)=(3.0d0*(x(i,1)-x(i,0))-(x(i,2)-x(i,1)))/(2.0d0*deta)
xeta(i,JM)=-(3.0d0*(x(i,JM-1)-x(i,JM))-(x(i,JM-2)-x(i,JM-1)))/(2.0d0*deta)
yeta(i,0)=(3.0d0*(y(i,1)-y(i,0))-(y(i,2)-y(i,1)))/(2.0d0*deta)
yeta(i,JM)=-(3.0d0*(y(i,JM-1)-y(i,JM))-(y(i,JM-2)-y(i,JM-1)))/(2.0d0*deta)
end do

do j=0,JM
xxi(0,j)=(3.0d0*(x(1,j)-x(0,j))-(x(2,j)-x(1,j)))/(2.0d0*dxi)
xxi(IM,j)=-(3.0d0*(x(IM-1,j)-x(IM,j))-(x(IM-2,j)-x(IM-1,j)))/(2.0d0*dxi)
yxi(0,j)=(3.0d0*(y(1,j)-y(0,j))-(y(2,j)-y(1,j)))/(2.0d0*dxi)
yxi(IM,j)=-(3.0d0*(y(IM-1,j)-y(IM,j))-(y(IM-2,j)-y(IM-1,j)))/(2.0d0*dxi)
end do


!********************ヤコビアンJaの計算********************
do i=0,IM
	do j=0,JM
	bj(i,j)=1.0d0/(xxi(i,j)*yeta(i,j)-xeta(i,j)*yxi(i,j))
	xix(i,j)=bj(i,j)*yeta(i,j)
	xiy(i,j)=-bj(i,j)*xeta(i,j)
	etax(i,j)=-bj(i,j)*yxi(i,j)
	etay(i,j)=bj(i,j)*xxi(i,j)
	end do
end do


!********************初期条件設定********************
do j=0,JM
	do i=0,IM
	T(i,j)=293.0d0
	rho(i,j)=1.2d0
	p(i,j)=rho(i,j)*R*T(i,j)
	u(i,j)=M(i,j)*dsqrt(gamma_h*R*T(i,j))
	v(i,j)=0.0d0
	ee(i,j)=rho(i,j)*(R*T(i,j)/(gamma_h-1.0d0)+(u(i,j)**2+v(i,j)**2)/2.0d0)
	end do
end do

!********************衝撃波条件（衝撃波通過後の変化）********************
do j=0,JM
	do i=0,IM
	if (y(i,j)>=1.0d0/(-dsqrt(3.0d0))*x(i,j)+y(0,JM-2)) then
	T2(i,j)=(2.0d0*gamma_h*M(i,j)**2*dsin(beta_rad)**2-(gamma_h-1.0d0))*&
	((gamma_h-1.0d0)*M(i,j)**2*dsin(beta_rad)**2+2.0d0)/((gamma_h+1.0d0)**2&
	*M(i,j)**2*dsin(beta_rad)**2)*T(i,j)
	p2(i,j)=(2.0d0*gamma_h*M(i,j)**2*dsin(beta_rad)**2-(gamma_h-1.0d0))/&
	(gamma_h+1.0d0)*p(i,j)
	rho2(i,j)=((gamma_h+1.0d0)*M(i,j)**2*dsin(beta_rad)**2)/&
	((gamma_h-1.0d0)*M(i,j)**2*dsin(beta_rad)**2+2.0d0)*rho(i,j)
	Vn1(i,j)=u(i,j)*dsin(beta_rad)
	Vn2(i,j)=((gamma_h-1.0d0)*M(i,j)**2*dsin(beta_rad)**2+2.0d0)/&
	((gamma_h+1.0d0)*M(i,j)**2*dsin(beta_rad)**2)*Vn1(i,j)
	Vt(i,j)=Vn1(i,j)/dtan(beta_rad)
	lv2(i,j)=dsqrt(Vn2(i,j)**2+Vt(i,j)**2)
	theta(i,j)=beta_rad-datan(Vn2(i,j)/Vt(i,j))
	u2(i,j)=lv2(i,j)*dcos(theta(i,j))
	v2(i,j)=lv2(i,j)*dsin(-theta(i,j))
	T(i,j)=T2(i,j)
	rho(i,j)=rho2(i,j)
	p(i,j)=p2(i,j)
	u(i,j)=u2(i,j)
	v(i,j)=v2(i,j)
	ee(i,j)=rho(i,j)*(R*T(i,j)/(gamma_h-1.0d0)+(u(i,j)**2+v(i,j)**2)/2.0d0)
	
	end if

	bu(i,j)=xix(i,j)*u(i,j)+xiy(i,j)*v(i,j)
	bv(i,j)=etax(i,j)*u(i,j)+etay(i,j)*v(i,j)
	end do
end do

!********************ベクトルQ　初期条件********************
do j=0,JM
	do i=0,IM
	q(i,j,1)=rho(i,j)/bj(i,j)
	q(i,j,2)=rho(i,j)*u(i,j)/bj(i,j)
	q(i,j,3)=rho(i,j)*v(i,j)/bj(i,j)
	q(i,j,4)=ee(i,j)/bj(i,j)
	end do
end do

!********************ルンゲクッタ法／メイン計算********************
do n=1,NMAX
	do k=1,4
		do j=0,JM
			do i=0,IM
			q1(i,j,k)=q(i,j,k)
			end do
		end do
	end do

	do kk=1,4
	call tvdxi
	call tvdeta
	call runge(1,IM-1,1,JM-1)
	call boundary
	end do

	print*,n,du,dv

	if (du<=EPS.and.dv<=EPS) exit
end do

call Output

contains

!********************サブルーチン①　境界条件********************
subroutine boundary
do j=0,JM
rho(IM,j)=rho(IM-1,j)+(rho(IM-1,j)-rho(IM-2,j))/dx*dx
u(IM,j)=u(IM-1,j)+(u(IM-1,j)-u(IM-2,j))/dx*dx
v(IM,j)=v(IM-1,j)+(v(IM-1,j)-v(IM-2,j))/dx*dx
ee(IM,j)=ee(IM-1,j)+(ee(IM-1,j)-ee(IM-2,j))/dx*dx
T(IM,j)=T(IM-1,j)+(T(IM-1,j)-T(IM-2,j))/dx*dx
p(IM,j)=p(IM-1,j)+(p(IM-1,j)-p(IM-2,j))/dx*dx
q(IM,j,1)=rho(IM,j)/bj(IM,j)
q(IM,j,2)=rho(IM,j)*u(IM,j)/bj(IM,j)
q(IM,j,3)=rho(IM,j)*v(IM,j)/bj(IM,j)
q(IM,j,4)=ee(IM,j)/bj(IM,j)
end do

do i=1,IM-1
rho(i,JM)=rho(i,JM-1)
u(i,0)=u(i,1)
u(i,JM)=u(i,JM-1)
v(i,0)=0.0d0
v(i,JM)=v(i,JM-1)
ee(i,0)=ee(i,1)
ee(i,JM)=ee(i,JM-1)
T(i,0)=T(i,1)
T(i,JM)=T(i,JM-1)
p(i,0)=p(i,1)
p(i,JM)=p(i,JM-1)
q(i,0,1)=rho(i,0)/bj(i,0)
q(i,0,2)=rho(i,0)*u(i,0)/bj(i,0)
q(i,0,3)=rho(i,0)*v(i,0)/bj(i,0)
q(i,0,4)=ee(i,0)/bj(i,0)
q(i,JM,1)=rho(i,JM)/bj(i,JM)
q(i,JM,2)=rho(i,JM)*u(i,JM)/bj(i,JM)
q(i,JM,3)=rho(i,JM)*v(i,JM)/bj(i,JM)
q(i,JM,4)=ee(i,JM)/bj(i,JM)
end do
end subroutine boundary


!********************ファンクション　FPSI関数の定義********************
double precision FUNCTION FPSI(Z,DELTA)
double precision DELTA,Z
IF(ABS(Z).GE.DELTA)THEN
FPSI=ABS(Z)
      ELSE
        FPSI=0.5d0*(Z**2+DELTA**2)/DELTA
      END IF
      return
      END function FPSI

!********************サブルーチン②　TVDスキーム_ξ方向********************
!*****EEが求まる*****
subroutine tvdxi
double precision RWL,RWR,AKX,AKY,AJACM,UM,VM,DUM,CM
double precision PHI,BETA,AKXT,AKYT,THIT
double precision,dimension(1:4,1:4,0:IM)	::RR=0.0d0
double precision,dimension(1:4,1:4,0:IM)	::RI=0.0d0
double precision,dimension(1:4,0:IM)		::EIGM,EIG,GG,PHIM,EH
double precision,dimension(1:4)::D
double precision,dimension(1:4,-1:IM)::ALPHA
double precision S,DELTA,GAMMA,UCONT

DO j=1,JM-1
     DO i=0,IM-1
        RWL=DSQRT(q(i,j,1)*bj(i,j))/(DSQRT(q(i+1,j,1)*bj(i+1,j))+DSQRT(q(i,j,1)*bj(i,j)))
        RWR=DSQRT(q(i+1,j,1)*bj(i+1,j))/(DSQRT(q(i+1,j,1)*bj(i+1,j))+DSQRT(q(i,j,1)*bj(i,j)))
        AKX=0.5d0*(xix(i,j)+xix(i+1,j))
        AKY=0.5d0*(xiy(i,j)+xiy(i+1,j))
        AJACM=0.5d0*(bj(i,j)+bj(i+1,j))
        UM=RWL*u(i,j)+RWR*u(i+1,j)
        VM=RWL*v(i,j)+RWR*v(i+1,j)
        DUM=RWL*q(i,j,4)/q(i,j,1)+RWR*q(i+1,j,4)/q(i+1,j,1)
        CM=DSQRT(gamma_h*(gamma_h-1.0d0)*ABS(DUM-0.5d0*(UM**2+VM**2)))
!-------NOMENCLATURE-------
        PHI=0.5d0*(gamma_h-1.0d0)*(UM**2+VM**2)
        BETA=1.0d0/(2.0d0*CM**2)
        AKXT=AKX/DSQRT(AKX**2+AKY**2)
        AKYT=AKY/DSQRT(AKX**2+AKY**2)
        THIT=AKXT*UM+AKYT*VM
!-------RIGHT ENGEN-VECTORS-------
        RR(1,1,i)=1.0d0
        RR(1,2,i)=0.0d0
        RR(1,3,i)=1.0d0
        RR(1,4,i)=1.0d0
        RR(2,1,i)=UM
        RR(2,2,i)=AKYT
        RR(2,3,i)=UM+AKXT*CM
        RR(2,4,i)=UM-AKXT*CM
        RR(3,1,i)=VM
        RR(3,2,i)=-AKXT
        RR(3,3,i)=VM+AKYT*CM
        RR(3,4,i)=VM-AKYT*CM
        RR(4,1,i)=PHI/(gamma_h-1.0d0)
        RR(4,2,i)=AKYT*UM-AKXT*VM
        RR(4,3,i)=(PHI+CM**2)/(gamma_h-1.0d0)+CM*THIT
        RR(4,4,i)=(PHI+CM**2)/(gamma_h-1.0d0)-CM*THIT
!-------INVERS OR RIGHT EIGEN-VECTORS-------
        RI(1,1,i)=1.0d0-PHI/CM**2
        RI(1,2,i)=(gamma_h-1.0d0)*UM/CM**2
        RI(1,3,i)=(gamma_h-1.0d0)*VM/CM**2
        RI(1,4,i)=-(gamma_h-1.0d0)/CM**2
        RI(2,1,i)=-AKYT*UM+AKXT*VM
        RI(2,2,i)=AKYT
        RI(2,3,i)=-AKXT
        RI(2,4,i)=0.0d0
        RI(3,1,i)=BETA*(PHI-CM*THIT)
        RI(3,2,i)=BETA*(AKXT*CM-(gamma_h-1.0d0)*UM)
        RI(3,3,i)=BETA*(AKYT*CM-(gamma_h-1.0d0)*VM)
        RI(3,4,i)=BETA*(gamma_h-1.0d0)
        RI(4,1,i)=BETA*(PHI+CM*THIT)
        RI(4,2,i)=-BETA*(AKXT*CM+(gamma_h-1.0d0)*UM)
        RI(4,3,i)=-BETA*(AKYT*CM+(gamma_h-1.0d0)*VM)
        RI(4,4,i)=BETA*(gamma_h-1.0d0)
!-------ENGEN-VALUES AT INTERMEDIATE-------
        EIGM(1,i)=AKX*UM+AKY*VM
        EIGM(2,i)=EIGM(1,i)
        EIGM(3,i)=EIGM(1,i)+CM*DSQRT(AKX**2+AKY**2)
        EIGM(4,i)=EIGM(1,i)-CM*DSQRT(AKX**2+AKY**2)
!-------ALPHA-------
do k=1,4
        D(k)=(q(i+1,j,k)*bj(i+1,j)-q(i,j,k)*bj(i,j))/AJACM
end do
do k=1,4
        ALPHA(k,i)=RI(k,1,i)*D(1)+RI(k,2,i)*D(2)+RI(k,3,i)*D(3)+RI(k,4,i)*D(4)
end do
end do
do k=1,4
        ALPHA(k,-1)=ALPHA(k,0)
        ALPHA(k,IM)=ALPHA(k,IM-1)
end do

do k=1,4
do i=0,IM
!-------MINMOD LIMITER-------
        S=DSIGN(1.0d0,ALPHA(k,i))
        GG(k,i)=S*DMAX1(0.0d0,DMIN1(ABS(ALPHA(k,i)),S*ALPHA(k,i-1)))
end do
end do
do k=1,4
do i=0,IM-1
!-------GAMMA-------
        DELTA=1.0d-10
        IF(ALPHA(k,i).NE.0.0d0)THEN
        GAMMA=0.5d0*FPSI(EIGM(k,i),DELTA)*(GG(k,i+1)-GG(k,i))/ALPHA(k,i)
        ELSE
        GAMMA=0.0d0
        END IF
!-------PHI-------
        PHIM(k,i)=0.5d0*FPSI(EIGM(k,i),DELTA)*(GG(k,i+1)+GG(k,i))-FPSI(EIGM(k,i)+GAMMA,DELTA)*ALPHA(k,i)
end do
end do
!-------CONVECTION COMPORNENTS-------
do i=0,IM
        UCONT=xix(i,j)*u(i,j)+xiy(i,j)*v(i,j)
        EH(1,i)=q(i,j,1)*UCONT
        EH(2,i)=q(i,j,2)*UCONT+p(i,j)*xix(i,j)/bj(i,j)
        EH(3,i)=q(i,j,3)*UCONT+p(i,j)*xiy(i,j)/bj(i,j)
        EH(4,i)=q(i,j,4)*UCONT+p(i,j)*UCONT/bj(i,j)
end do
!-------XI-DIRECTION CONVECTION FLUX-------
        do k=1,4
        do i=0,IM-1
        e(i,j,k)=0.5d0*(EH(k,i)+EH(k,i+1)+RR(k,1,i)*PHIM(1,i)+RR(k,2,i)*PHIM(2,i)&
+RR(k,3,i)*PHIM(3,i)+RR(k,4,i)*PHIM(4,i))
end do
end do
end do
end subroutine tvdxi


!********************サブルーチン③　TVDスキーム_η方向********************
!*****FFが求まる*****
subroutine tvdeta
double precision RWL,RWR,AKX,AKY,AJACM,UM,VM,DUM,CM
double precision PHI,BETA,AKXT,AKYT,THIT
double precision,dimension(1:4,1:4,0:JM)::RR,RJ
double precision,dimension(1:4,0:JM)::EIGM,EIG,GG,PHIM,FH
double precision,dimension(1:4)::D
double precision,dimension(1:4,-1:JM)::ALPHA
double precision S,DELTA,GAMMA,VCONT

 DO i=1,IM-1
 DO j=0,JM-1
      
        RWL=DSQRT(q(i,j,1)*bj(i,j))/(DSQRT(q(i,j+1,1)*bj(i,j+1))+DSQRT(q(i,j,1)*bj(i,j)))
        RWR=DSQRT(q(i,j+1,1)*bj(i,j+1))/(DSQRT(q(i,j+1,1)*bj(i,j+1))+DSQRT(q(i,j,1)*bj(i,j)))
        AKX=0.5d0*(etax(i,j)+etax(i,j+1))
        AKY=0.5d0*(etay(i,j)+etay(i,j+1))
        AJACM=0.5d0*(bj(i,j)+bj(i,j+1))
        UM=RWL*u(i,j)+RWR*u(i,j+1)
        VM=RWL*v(i,j)+RWR*v(i,j+1)
        DUM=RWL*q(i,j,4)/q(i,j,1)+RWR*q(i,j+1,4)/q(i,j+1,1)
        CM=DSQRT(gamma_h*(gamma_h-1.0d0)*ABS(DUM-0.5d0*(UM**2+VM**2)))
!NOMENCLATURE
        PHI=0.5d0*(gamma_h-1.0d0)*(UM**2+VM**2)
        BETA=1.0d0/(2.0d0*CM**2)
        AKXT=AKX/DSQRT(AKX**2+AKY**2)
        AKYT=AKY/DSQRT(AKX**2+AKY**2)
        THIT=AKXT*UM+AKYT*VM
!RIGHT ENGEN-VECTORS
        RR(1,1,j)=1.0d0
        RR(1,2,j)=0.0d0
        RR(1,3,j)=1.0d0
        RR(1,4,j)=1.0d0
        RR(2,1,j)=UM
        RR(2,2,j)=AKYT
        RR(2,3,j)=UM+AKXT*CM
        RR(2,4,j)=UM-AKXT*CM
        RR(3,1,j)=VM
        RR(3,2,j)=-AKXT
        RR(3,3,j)=VM+AKYT*CM
        RR(3,4,j)=VM-AKYT*CM
        RR(4,1,j)=PHI/(gamma_h-1.0d0)
        RR(4,2,j)=AKYT*UM-AKXT*VM
        RR(4,3,j)=(PHI+CM**2)/(gamma_h-1.0d0)+CM*THIT
        RR(4,4,j)=(PHI+CM**2)/(gamma_h-1.0d0)-CM*THIT
!INVERS OR RIGHT EIGEN-VECTORS
        RJ(1,1,j)=1.0d0-PHI/CM**2
        RJ(1,2,j)=(gamma_h-1.0d0)*UM/CM**2
        RJ(1,3,j)=(gamma_h-1.0d0)*VM/CM**2
        RJ(1,4,j)=-(gamma_h-1.0d0)/CM**2
        RJ(2,1,j)=-AKYT*UM+AKXT*VM
        RJ(2,2,j)=AKYT
        RJ(2,3,j)=-AKXT
        RJ(2,4,j)=0.0d0
        RJ(3,1,j)=BETA*(PHI-CM*THIT)
        RJ(3,2,j)=BETA*(AKXT*CM-(gamma_h-1.0d0)*UM)
        RJ(3,3,j)=BETA*(AKYT*CM-(gamma_h-1.0d0)*VM)
        RJ(3,4,j)=BETA*(gamma_h-1.0d0)
        RJ(4,1,j)=BETA*(PHI+CM*THIT)
        RJ(4,2,j)=-BETA*(AKXT*CM+(gamma_h-1.0d0)*UM)
        RJ(4,3,j)=-BETA*(AKYT*CM+(gamma_h-1.0d0)*VM)
        RJ(4,4,j)=BETA*(gamma_h-1.0d0)
!ENGEN-VALUES AT INTERMEDIATE
        EIGM(1,j)=AKX*UM+AKY*VM
        EIGM(2,j)=EIGM(1,j)
        EIGM(3,j)=EIGM(1,j)+CM*DSQRT(AKX**2+AKY**2)
        EIGM(4,j)=EIGM(1,j)-CM*DSQRT(AKX**2+AKY**2)
!ALPHA
do k=1,4
        D(k)=(q(i,j+1,k)*bj(i,j+1)-q(i,j,k)*bj(i,j))/AJACM
end do
do k=1,4
        ALPHA(k,j)=RJ(k,1,j)*D(1)+RJ(k,2,j)*D(2)+RJ(k,3,j)*D(3)+RJ(k,4,j)*D(4)
end do
end do
do k=1,4
        ALPHA(k,-1)=ALPHA(k,0)
        ALPHA(k,JM)=ALPHA(k,JM-1)
end do

do k=1,4
do j=0,JM
!MINMOD LIMITER
        S=DSIGN(1.0d0,ALPHA(k,j))
        GG(k,j)=S*DMAX1(0.0d0,DMIN1(ABS(ALPHA(k,j)),S*ALPHA(k,j-1)))
end do
end do
do k=1,4
do j=0,JM-1
!GAMMA
        DELTA=1.0d-10
        !DELTA=AMAX1(0.0,EIGM(K,J)-EIG(K,J),EIG(K,J+1)-EIGM(K,J))
        IF(ALPHA(k,j).NE.0.0)THEN
        GAMMA=0.5d0*FPSI(EIGM(k,j),DELTA)*(GG(k,j+1)-GG(k,j))/ALPHA(k,j)
        ELSE
        GAMMA=0.0d0
        END IF
!PHI
        PHIM(k,j)=0.5d0*FPSI(EIGM(k,j),DELTA)*(GG(k,j+1)+GG(k,j))-FPSI(EIGM(k,j)+GAMMA,DELTA)*ALPHA(k,j)
end do
end do
!CONVECTION COMPORNENTS
do j=0,JM
        VCONT=etax(i,j)*u(i,j)+etay(i,j)*v(i,j)
        FH(1,j)=q(i,j,1)*VCONT
        FH(2,j)=q(i,j,2)*VCONT+p(i,j)*etax(i,j)/bj(i,j)
        FH(3,j)=q(i,j,3)*VCONT+p(i,j)*etay(i,j)/bj(i,j)
        FH(4,j)=q(i,j,4)*VCONT+p(i,j)*VCONT/bj(i,j)
end do
!XI-DIRECTION CONVECTION FLUX
        do k=1,4
        do j=0,JM-1
        f(i,j,k)=0.5d0*(FH(k,j)+FH(k,j+1)+RR(k,1,j)*PHIM(1,j)+RR(k,2,j)*PHIM(2,j)&
+RR(k,3,j)*PHIM(3,j)+RR(k,4,j)*PHIM(4,j))
end do
end do
end do
end subroutine tvdeta


!********************サブルーチン④　ルンゲクッタ********************
subroutine runge(is,ie,js,je)
integer ik,is,ie,js,je
double precision uu,ddu,vv,ddv
du=0.d0
dv=0.d0
do j=js,je
	do i=is,ie
		do k=1,4
		q(i,j,k)=q1(i,j,k)-1.0d0/(5.0d0-dble(kk))*dt*((e(i,j,k)&
		-e(i-1,j,k))/dxi+(f(i,j,k)-f(i,j-1,k))/deta)
		end do
	uu=u(i,j)
	vv=v(i,j)
	rho(i,j)=q(i,j,1)*bj(i,j)
	u(i,j)=q(i,j,2)*bj(i,j)/rho(i,j)
	v(i,j)=q(i,j,3)*bj(i,j)/rho(i,j)
	bu(i,j)=xix(i,j)*u(i,j)+xiy(i,j)*v(i,j)
	bv(i,j)=etax(i,j)*u(i,j)+etay(i,j)*v(i,j)
	ee(i,j)=q(i,j,4)*bj(i,j)
	T(i,j)=(gamma_h-1.0d0)*(ee(i,j)/rho(i,j)-(u(i,j)**2+v(i,j)**2)/2.0d0)/R
	p(i,j)=rho(i,j)*R*T(i,j)

	ddu=abs(u(i,j)-uu)/uu

	if (ddu>du) then
	du=ddu
	
	end if
	
	ddv=abs(v(i,j)-vv)/vv
	if (ddv>dv) then
	dv=ddv
	end if
	end do
end do
end subroutine runge


!出力
subroutine Output
open(1,file='kadai05.dat', status='replace')
do j=0,JM
do i=0,IM
write (1,*) u(i,j),v(i,j),p(i,j),T(i,j)
end do
end do
close(1)
end subroutine Output



end program kadai05
