program assyuku8_1
implicit none

!----------定義----------
double precision,parameter::H=0.2d0
integer,parameter::IM=50,JM=50
double precision,parameter::EPS=1.0d-6,R=287.1d0,g=1.4d0,Pr=0.9d0,Ret=500.0d0,T0=293.0d0
double precision,parameter::myu0=1.0d-3,SS=110.0d0
double precision,parameter::Cmyu=0.09d0,sigmak=1.0d0,sigmaepsi=1.3d0
double precision::du,dv,ddu,ddv,slope
integer::i,j,n,kk,rk
double precision,parameter::M=2.9d0
double precision,dimension(0:IM,0:JM)::x,y,u,v,p,t
double precision,dimension(0:IM,0:JM)::xe,ye,xh,yh,xix,xiy,ex,ey,Jac
double precision,dimension(1:6,0:IM,0:JM)::Q,EE,FF,SSS,RRR,PPP
double precision,dimension(1:6,0:IM,0:JM,0:5)::QQ
double precision,dimension(0:IM,0:JM)::kwall,epsiwall,ldt=1.0d-8,enebar,c,un,vn
double precision,dimension(0:IM,0:JM)::rho,ene,k,epsi,myu,nyu,myut


!------メイン計算------
call mesh
call conditions
call coordinate
n=0
do n=0,100000
	do j=0,JM
		do i=0,IM
			do kk=1,6
				QQ(kk,i,j,0)=Q(kk,i,j)
			end do
		end do
	end do
	do rk=1,4						!Runge-Kuttaの段階数
		call tvd_xi
		call tvd_eta
		call parameter
		call viscosity_xi
		call viscosity_eta
		call production
		call wall_law
		call Runge_Kutta
		call calculation
		call local_time_step_method
	end do
	print*,n
	if(du<EPS.and.dv<EPS)exit
end do
call output_yplus
call output

!------------------------------------------------------------
contains

!------格子生成------
subroutine mesh
double precision::dx,ix,iy
double precision,parameter::a=1.2d0,b=(a+1.0d0)/(a-1.0d0)
do j=0,JM
	do i=0,IM
		if(i<=5)then
		dx=50.0d0*H/50.0d0
		x(i,j)=dble(i)*dx
		else
		ix=dble(i)/IM
		x(i,j)=x(5,j)+45.0d0*H/(b-1.0d0)*(b**((dble(i-5))/((dble(IM-5))))-1.0d0)
		end if
		iy=dble(j)/dble(JM)
		y(i,j)=H*(b**iy-1.0d0)/(b-1.0d0)
	end do
end do
end subroutine mesh

!------初期条件------
subroutine conditions
do j=0,JM
	do i=0,IM
	T(i,j)=293.15d0
	rho(i,j)=1.2d0
	c(i,j)=dsqrt(g*R*T(i,j))
	u(i,j)=M*c(i,j)
	v(i,j)=0.0d0
	p(i,j)=rho(i,j)*R*T(i,j)
	ene(i,j)=p(i,j)/(g-1.0d0)+rho(i,j)*((u(i,j)**2.0+v(i,j)**2.0))/2.0d0
	myu(i,j)=1.82d-5
	k(i,j)=1.5d0*(1.0d-2*u(i,j))**2.0
	nyu(i,j)=myu(i,j)/rho(i,j)
	epsi(i,j)=(k(i,j)**2.0)/(nyu(i,j)*Ret)
	end do
end do

call wall_law

do i=5,IM
	u(i,0)=0.0d0
	v(i,0)=0.0d0
	ene(i,0)=p(i,0)/(g-1.0d0)+rho(i,0)*((u(i,0)**2.0+v(i,0)**2.0))/2.0d0
	nyu(i,0)=myu(i,0)/rho(i,0)
	k(i,0)=1.5d0*(1.0d-2*u(i,0))**2.0
	epsi(i,0)=(k(i,0)**2.0)/(nyu(i,0)*Ret)
end do

do j=0,JM
	do i=0,IM
	myut(i,j)=Cmyu*rho(i,j)*(k(i,j)**2.0)/epsi(i,j)
	end do
end do

end subroutine conditions

!------座標変換------
subroutine coordinate
do j=0,JM
	do i=1,IM-1
		xe(i,j)=(x(i+1,j)-x(i-1,j))/2.0d0
		ye(i,j)=(y(i+1,j)-y(i-1,j))/2.0d0
	end do
end do

do j=1,JM-1
	do i=0,IM
		xh(i,j)=(x(i,j+1)-x(i,j-1))/2.0d0
		yh(i,j)=(y(i,j+1)-y(i,j-1))/2.0d0
	end do
end do

do j=0,JM
	xe(0,j)=(-3.0*x(0,j)+4.0*x(1,j)-x(2,j))/2.0d0
	xe(IM,j)=(3.0*x(IM,j)-4.0*x(IM-1,j)+x(IM-2,j))/2.0d0
	ye(0,j)=(-3.0*y(0,j)+4.0*y(1,j)-y(2,j))/2.0d0
	ye(IM,j)=(3.0*y(IM,j)-4.0*y(IM-1,j)+y(IM-2,j))/2.0d0
end do

do i=0,IM
	xh(i,0)=(-3.0*x(i,0)+4.0*x(i,1)-x(i,2))/2.0d0
	xh(i,JM)=(3.0*x(i,JM)-4.0*x(i,JM-1)+x(i,JM-2))/2.0d0
	yh(i,0)=(-3.0*y(i,0)+4.0*y(i,1)-y(i,2))/2.0d0
	yh(i,JM)=(3.0*y(i,JM)-4.0*y(i,JM-1)+y(i,JM-2))/2.0d0
end do

do j=0,JM
	do i=0,IM
		Jac(i,j)=1.0d0/(xe(i,j)*yh(i,j)-ye(i,j)*xh(i,j))
		xix(i,j)=yh(i,j)*Jac(i,j)
		xiy(i,j)=-xh(i,j)*Jac(i,j)
		ex(i,j)=-ye(i,j)*Jac(i,j)
		ey(i,j)=xe(i,j)*Jac(i,j)
	end do
end do

do j=0,JM
	do i=0,IM
		Q(1,i,j)=rho(i,j)/Jac(i,j)
		Q(2,i,j)=(rho(i,j)*u(i,j))/Jac(i,j)
		Q(3,i,j)=(rho(i,j)*v(i,j))/Jac(i,j)
		Q(4,i,j)=ene(i,j)/Jac(i,j)
		Q(5,i,j)=rho(i,j)*k(i,j)/Jac(i,j)
		Q(6,i,j)=rho(i,j)*epsi(i,j)/Jac(i,j)
	end do
end do
end subroutine coordinate

!------ξ方向TVD------
subroutine tvd_xi
double precision:: RWL,RWR,AKL,AKX,AKY,AJACM,UM,VM,DUM,CM,PHI,BETA,AKXT,AKYT,THIT,AKM,S
double precision,dimension(1:4,1:4,0:IM)::RR,RI
double precision,dimension(1:6,0:IM)::EIGM,EIG,GG,E
double precision,dimension(1:6)::D
double precision,dimension(1:6,-1:IM)::ALPHA
double precision,dimension(1:6,0:IM)::PHIM
double precision::DELTA,GAMMA,UCONT

do j=1,JM-1
	do i=0,IM-1
		RWL=dSQRT(Q(1,I,J)*Jac(I,J))/(dSQRT(Q(1,I+1,J)*Jac(I+1,J))+dSQRT(Q(1,I,J)*Jac(I,J)))
		RWR=dSQRT(Q(1,I+1,J)*Jac(I+1,J))/(dSQRT(Q(1,I+1,J)*Jac(I+1,J))+dSQRT(Q(1,I,J)*Jac(I,J)))
		AKX=0.5d0*(xix(I,J)+xix(I+1,J))!中間点のkx
		AKY=0.5d0*(xiy(I,J)+xiy(I+1,J))!中間点のky
		AJACM=0.5d0*(Jac(I,J)+Jac(I+1,J))!中間点のJac

		UM=RWL*U(I,J)+RWR*U(I+1,J)!中間点のu
		VM=RWL*V(I,J)+RWR*V(I+1,J)!中間点のv
		DUM=RWL*Q(4,I,J)/Q(1,I,J)+RWR*Q(4,I+1,J)/Q(1,I+1,J)!e/ρ
		CM=dSQRT(G*(G-1.0d0)*ABS(DUM-0.5d0*(UM**2.0+VM**2.0)))!中間点の音速c
		!NOMENCLATURE
		PHI=0.5d0*(G-1.0d0)*(UM**2.0+VM**2.0)!Φ
		BETA=1.0d0/(2.0d0*CM**2.0)
		AKXT=AKX/dSQRT(AKX**2.0+AKY**2.0)
		AKYT=AKY/dSQRT(AKX**2.0+AKY**2.0)
		THIT=AKXT*UM+AKYT*VM

		!RIGHT ENGEN-VECTORS(右固有ベクトル）
		RR(1,1,I)=1.0d0
		RR(1,2,I)=0.0d0
		RR(1,3,I)=1.0d0
		RR(1,4,I)=1.0d0
		RR(2,1,I)=UM
		RR(2,2,I)=AKYT
		RR(2,3,I)=UM+AKXT*CM
		RR(2,4,I)=UM-AKXT*CM
		RR(3,1,I)=VM
		RR(3,2,I)=-AKXT
		RR(3,3,I)=VM+AKYT*CM
		RR(3,4,I)=VM-AKYT*CM
		RR(4,1,I)=PHI/(G-1.0d0)
		RR(4,2,I)=AKYT*UM-AKXT*VM
		RR(4,3,I)=(PHI+CM**2.0)/(G-1.0d0)+CM*THIT
		RR(4,4,I)=(PHI+CM**2.0)/(G-1.0d0)-CM*THIT

		!INVERS OR RIGHT EIGEN-VECTORS(右固有ベクトルの逆行列）
		RI(1,1,I)=1.0d0-PHI/CM**2.0
		RI(1,2,I)=(G-1.0d0)*UM/CM**2.0
		RI(1,3,I)=(G-1.0d0)*VM/CM**2.0
		RI(1,4,I)=-(G-1.0d0)/CM**2.0
		RI(2,1,I)=-AKYT*UM+AKXT*VM
		RI(2,2,I)=AKYT
		RI(2,3,I)=-AKXT
		RI(2,4,I)=0.0d0
		RI(3,1,I)=BETA*(PHI-CM*THIT)
		RI(3,2,I)=BETA*(AKXT*CM-(G-1.0d0)*UM)
		RI(3,3,I)=BETA*(AKYT*CM-(G-1.0d0)*VM)
		RI(3,4,I)=BETA*(G-1.0d0)
		RI(4,1,I)=BETA*(PHI+CM*THIT)
		RI(4,2,I)=-BETA*(AKXT*CM+(G-1.0d0)*UM)
		RI(4,3,I)=-BETA*(AKYT*CM+(G-1.0d0)*VM)
		RI(4,4,I)=BETA*(G-1.0d0)

		!ENGEN-VALUES AT INTERMEDATE(固有値）
		EIGM(1,I)=AKX*UM+AKY*VM
		EIGM(2,I)=EIGM(1,I)
		EIGM(3,I)=EIGM(1,I)+CM*dSQRT(AKX**2.0+AKY**2.0)
		EIGM(4,I)=EIGM(1,I)-CM*dSQRT(AKX**2.0+AKY**2.0)
		EIGM(5,I)=EIGM(1,I)
		EIGM(6,I)=EIGM(1,I)

		!ALPHAの導出
		do kk=1,6
			D(kk)=(Q(kk,I+1,J)*Jac(I+1,J)-Q(kk,I,J)*Jac(I,J))/AJACM
		end do

		do kk=1,6
			if(Kk.LE.4)then
				ALPHA(kk,I)=RI(kk,1,I)*D(1)+RI(kk,2,I)*D(2)+RI(kk,3,I)*D(3)+RI(kk,4,I)*D(4)
			else
				ALPHA(kk,I)=D(kk)
			end if
		end do
	end do

	do kk=1,6
		ALPHA(kk,-1)=ALPHA(kk,0)
		ALPHA(kk,IM)=ALPHA(kk,IM-1)
	end do

	do i=0,IM
		DUM=dSQRT(G*R*T(I,J)*(xix(I,J)**2+xiy(I,J)**2))
		EIG(1,I)=xix(I,J)*U(I,J)+xiy(I,J)*V(I,J)
		EIG(2,I)=EIG(1,I)
		EIG(3,I)=EIG(1,I)+DUM
		EIG(4,I)=EIG(1,I)-DUM
		EIG(5,I)=EIG(1,I)
		EIG(6,I)=EIG(1,I)
	end do

	!CAL LIMITER FUNCTION G(制御関数)
	do kk=1,6
		do I=0,IM
			!MINMOD LIMITER
			!if(ALPHA(Kk,I)>=0.0)b=1.0
			!if(ALPHA(Kk,I)<0.0)b=-1.0
			s=dSIGN(1.0d0,ALPHA(kk,I))
			GG(kk,I)=s*dMAX1(0.0,dMIN1(dabs(ALPHA(kk,I)),s*ALPHA(kk,I-1)))
		end do
	end do

	do kk=1,6
		do I=0,IM-1
		!GAMMA
			DELTA=1.0d-10
			!DELTA=dMAX1(0.0,EIGM(kk,I)-EIG(kk,I),EIG(kk,I+1)-EIGM(kk,I))
			IF(ALPHA(kk,I).NE.0.0d0)THEN
				GAMMA=0.5d0*FPSI(EIGM(kk,I),DELTA)*(GG(kk,I+1)-GG(kk,I))/ALPHA(kk,I)
			ELSE
				GAMMA=0.0d0
			END IF
			!PHI(Φ)
			PHIM(kk,I)=0.5d0*FPSI(EIGM(kk,I),DELTA)*(GG(kk,I+1)+GG(kk,I))-FPSI(EIGM(kk,I)+GAMMA,DELTA)*ALPHA(kk,I)
		end do
	end do
	!CONVECTION COMPORNENTS(計算空間上のE計算)
	do I=0,IM
		UCONT=xix(I,J)*U(I,J)+xiy(I,J)*V(I,J)!(反変速度)
		E(1,I)=Q(1,I,J)*UCONT
		E(2,I)=Q(2,I,J)*UCONT+P(I,J)*xix(I,J)/Jac(i,j)
		E(3,I)=Q(3,I,J)*UCONT+P(I,J)*xiy(I,J)/Jac(i,j)
		E(4,I)=Q(4,I,J)*UCONT+P(I,J)*UCONT/Jac(i,j)
		E(5,I)=Q(5,I,J)*UCONT
		E(6,I)=Q(6,I,J)*UCONT
	end do
	!XI-DRECTION CONVECTION FLUX
	do kk=1,6
		do I=0,IM-1
			if(Kk.LE.4.0d0)then
				EE(kk,I,J)=0.5d0*(E(kk,I)+E(kk,I+1)+RR(kk,1,I)*PHIM(1,I)&
							&+RR(kk,2,I)*PHIM(2,I)+RR(kk,3,I)*PHIM(3,I)+RR(kk,4,I)*PHIM(4,I))
			else
				EE(kk,I,J)=0.5d0*(E(kk,I)+E(kk,I+1)+PHIM(kk,I))
			end if
		end do
	end do
end do
return
end subroutine tvd_xi

!------η方向TVD------
subroutine tvd_eta
double precision:: RWL,RWR,AKL,AKX,AKY,AJACM,UM,VM,DUM,CM,PHI,BETA,AKXT,AKYT,THIT,AKM,S
double precision,dimension(1:4,1:4,0:JM)::RR,RI
double precision,dimension(1:6,0:JM)::FIGM,FIG,GG,F
double precision,dimension(1:6)::D
double precision,dimension(1:6,-1:JM)::ALPHA
double precision,dimension(1:6,0:JM)::PHIM
double precision::DELTA,GAMMA,UCONT

do I=1,IM-1
	do J=0,JM-1
		RWL=dSQRT(Q(1,I,J)*Jac(I,J))/(dSQRT(Q(1,I,J+1)*Jac(I,J+1))+dSQRT(Q(1,I,J)*Jac(I,J)))
		RWR=dSQRT(Q(1,I,J+1)*Jac(I,J+1))/(dSQRT(Q(1,I,J+1)*Jac(I,J+1))+dSQRT(Q(1,I,J)*Jac(I,J)))
		AKX=0.5d0*(ex(I,J)+ex(I,J+1))
		AKY=0.5d0*(ey(I,J)+ey(I,J+1))
		AJACM=0.5d0*(Jac(I,J)+Jac(I,J+1))
		UM=RWL*U(I,J)+RWR*U(I,J+1)
		VM=RWL*V(I,J)+RWR*V(I,J+1)
		DUM=RWL*Q(4,I,J)/Q(1,I,J)+RWR*Q(4,I,J+1)/Q(1,I,J+1)
		CM=dSQRT(G*(G-1.0d0)*dabs(DUM-0.5d0*(UM**2.0+VM**2.0)))

		!NOMENCLATURE
		PHI=0.5d0*(G-1.0d0)*(UM**2.0+VM**2.0)
		BETA=1.0d0/(2.0d0*CM**2.0)
		AKXT=AKX/dSQRT(AKX**2.0+AKY**2.0)
		AKYT=AKY/dSQRT(AKX**2.0+AKY**2.0)
		THIT=AKXT*UM+AKYT*VM

		!RIGHT ENGEN-VECTORS
		RR(1,1,J)=1.0d0
		RR(1,2,J)=0.0d0
		RR(1,3,J)=1.0d0
		RR(1,4,J)=1.0d0
		RR(2,1,J)=UM
		RR(2,2,J)=AKYT
		RR(2,3,J)=UM+AKXT*CM
		RR(2,4,J)=UM-AKXT*CM
		RR(3,1,J)=VM
		RR(3,2,J)=-AKXT
		RR(3,3,J)=VM+AKYT*CM
		RR(3,4,J)=VM-AKYT*CM
		RR(4,1,J)=PHI/(G-1.0d0)
		RR(4,2,J)=AKYT*UM-AKXT*VM
		RR(4,3,J)=(PHI+CM**2.0)/(G-1.0d0)+CM*THIT
		RR(4,4,J)=(PHI+CM**2.0)/(G-1.0d0)-CM*THIT

		!INVERS OR RIGHT EIGEN-VECTORS
		RI(1,1,J)=1.0d0-PHI/CM**2.0
		RI(1,2,J)=(G-1.0d0)*UM/CM**2.0
		RI(1,3,J)=(G-1.0d0)*VM/CM**2.0
		RI(1,4,J)=-(G-1.0d0)/CM**2.0
		RI(2,1,J)=-AKYT*UM+AKXT*VM
		RI(2,2,J)=AKYT
		RI(2,3,J)=-AKXT
		RI(2,4,J)=0.0d0
		RI(3,1,J)=BETA*(PHI-CM*THIT)
		RI(3,2,J)=BETA*(AKXT*CM-(G-1.0d0)*UM)
		RI(3,3,J)=BETA*(AKYT*CM-(G-1.0d0)*VM)
		RI(3,4,J)=BETA*(G-1.0d0)
		RI(4,1,J)=BETA*(PHI+CM*THIT)
		RI(4,2,J)=-BETA*(AKXT*CM+(G-1.0d0)*UM)
		RI(4,3,J)=-BETA*(AKYT*CM+(G-1.0d0)*VM)
		RI(4,4,J)=BETA*(G-1.0d0)

		!ENGEN-VALUES AT INTERMEDIATE
		FIGM(1,J)=AKX*UM+AKY*VM
		FIGM(2,J)=FIGM(1,J)
		FIGM(3,J)=FIGM(1,J)+CM*dSQRT(AKX**2.0+AKY**2.0)
		FIGM(4,J)=FIGM(1,J)-CM*dSQRT(AKX**2.0+AKY**2.0)
		FIGM(5,J)=FIGM(1,J)
		FIGM(6,J)=FIGM(1,J)

		!ALPHA
		do kk=1,6
			D(kk)=(Q(kk,I,J+1)*Jac(I,J+1)-Q(kk,I,J)*Jac(I,J))/AJACM
		end do

		do kk=1,6
			IF(Kk.LE.4) THEN
				ALPHA(kk,J)=RI(kk,1,J)*D(1)+RI(kk,2,J)*D(2)+RI(kk,3,J)*D(3)+RI(kk,4,J)*D(4)
			else
				ALPHA(kk,J)=D(kk)
			end if
		end do
	end do

	do kk=1,6
		ALPHA(kk,-1)=ALPHA(kk,0)
		ALPHA(kk,JM)=ALPHA(kk,JM-1)
	end do

	do J=0,JM
		DUM=dSQRT(G*R*T(I,J)*(ex(I,J)**2.0+ey(I,J)**2.0))
		FIG(1,J)=ex(I,J)*U(I,J)+ey(I,J)*V(I,J)
		FIG(2,J)=FIG(1,J)
		FIG(3,J)=FIG(1,J)+DUM
		FIG(4,J)=FIG(1,J)-DUM
		FIG(5,I)=FIG(1,I)
		FIG(6,I)=FIG(1,I)
	end do
	!CAL LIMITER FUNCTION G
	do kk=1,6
		do J=0,JM
	!MINMOD LIMITER 制限関数
			!if(ALPHA(Kk,I)>=0.0)b=1.0
			!if(ALPHA(Kk,I)<0.0)b=-1.0
			s=dSIGN(1.0d0,ALPHA(kk,J))
			GG(kk,J)=s*dMAX1(0.0d0,dMIN1(dabs(ALPHA(kk,J)),s*ALPHA(kk,J-1)))
		end do
	end do

	do kk=1,6
		do J=0,JM-1
			!GAMMA
			DELTA=1.0d-10
			!DELTA=dMAX1(0.0,FIGM(kk,I)-FIG(kk,I),FIG(kk,I+1)-FIGM(kk,I))
			if(ALPHA(kk,J).NE.0.0d0)then
				GAMMA=0.5d0*FPSI(FIGM(kk,J),DELTA)*(GG(kk,J+1)-GG(kk,J))/ALPHA(kk,J)
			else
				GAMMA=0.0d0
			end if
			!PHI
			PHIM(kk,J)=0.5d0*FPSI(FIGM(kk,J),DELTA)*(GG(kk,J+1)+GG(kk,J))-FPSI(FIGM(kk,J)+GAMMA,DELTA)*ALPHA(kk,J)
		end do
	end do

	!CONVECTION COMPORNENTS
	do J=0,JM
		UCONT=ex(I,J)*U(I,J)+ey(I,J)*V(I,J)
		F(1,J)=Q(1,I,J)*UCONT
		F(2,J)=Q(2,I,J)*UCONT+P(I,J)/Jac(i,j)*ex(I,J)
		F(3,J)=Q(3,I,J)*UCONT+P(I,J)/Jac(i,j)*ey(I,J)
		F(4,J)=Q(4,I,J)*UCONT+P(I,J)/Jac(i,j)*UCONT
		F(5,J)=Q(5,I,J)*UCONT
		F(6,J)=Q(6,I,J)*UCONT
	end do

	!XI-DRECTION CONVECTION FLUX
	do Kk=1,6
		do J=0,JM-1
			if(kk<=4)then
				FF(KK,I,J)=0.5d0*(F(Kk,J)+F(Kk,J+1)+RR(Kk,1,J)*PHIM(1,J)+RR(Kk,2,J)*PHIM(2,J)&
				&+RR(Kk,3,J)*PHIM(3,J)+RR(Kk,4,J)*PHIM(4,J))
			else
				FF(KK,I,J)=0.5d0*(F(Kk,J)+F(Kk,J+1)+PHIM(Kk,J))
			end if
		end do
	end do
end do
return
end subroutine tvd_eta

!*****FPSI関数の定義*****
double precision function FPSI(Z,DELTA)
double precision Z,DELTA
	if(dabs(Z).GE.DELTA)then
		FPSI=dabs(Z)
	else
		FPSI=0.5d0*(Z**2.0+DELTA**2.0)/DELTA
	end if
return
end function

!------サザーランド------
double precision function Sutherland(T)
double precision T
Sutherland=myu0*((T/T0)**1.5)*(T0+SS)/(T+SS)
end function

!------諸量計算------
subroutine parameter
do j=0,JM
	do i=0,IM
		myu(i,j)=Sutherland(T(i,j))
		nyu(i,j)=myu(i,j)/rho(i,j)
		myut(i,j)=Cmyu*rho(i,j)*(k(i,j)**2.0)/epsi(i,j)
	end do
end do
end subroutine parameter

!------粘性項（ξ方向）------
subroutine viscosity_xi
double precision,dimension(0:IM-1,0:JM-1)::xix_r,xiy_r,ex_r,ey_r,Jac_r,u_r,v_r,T_r,k_r,epsi_r,rho_r
double precision,dimension(0:IM-1,0:JM-1)::ue_r,ve_r,Te_r,ke_r,epsie_r,uh_r,vh_r,Th_r,kh_r,epsih_r,myut_r
double precision::myu,Sxx,Syy,Sxy,Syx,txx,tyy,txy,tyx,R4,R5,R6,S4,S5,S6

do j=1,JM-1
	do i=0,IM-1
		xix_r(i,j)=(xix(i,j)+xix(i+1,j))/2.0d0
		xiy_r(i,j)=(xiy(i,j)+xiy(i+1,j))/2.0d0
		ex_r(i,j)=(ex(i,j)+ex(i+1,j))/2.0d0
		ey_r(i,j)=(ey(i,j)+ey(i+1,j))/2.0d0
		Jac_r(i,j)=(Jac(i,j)+Jac(i+1,j))/2.0d0
		u_r(i,j)=(u(i,j)+u(i+1,j))/2.0d0
		v_r(i,j)=(v(i,j)+v(i+1,j))/2.0d0
		T_r(i,j)=(t(i,j)+t(i+1,j))/2.0d0
		k_r(i,j)=(k(i,j)+k(i+1,j))/2.0d0
		epsi_r(i,j)=(epsi(i,j)+epsi(i+1,j))/2.0d0
		rho_r(i,j)=(Q(1,i,j)*Jac(i,j)+Q(1,i+1,j)*Jac(i+1,j))/(Jac(i,j)+Jac(i+1,j))*Jac_r(i,j)

		ue_r(i,j)=u(i+1,j)-u(i,j)
		ve_r(i,j)=v(i+1,j)-v(i,j)
		Te_r(i,j)=t(i+1,j)-t(i,j)
		ke_r(i,j)=k(i+1,j)-k(i,j)
		epsie_r(i,j)=epsi(i+1,j)-epsi(i,j)

		uh_r(i,j)=(((u(i,j+1)+u(i+1,j+1))/2.0d0)-((u(i,j-1)+u(i+1,j-1))/2.0d0))/2.0d0
		vh_r(i,j)=(((v(i,j+1)+v(i+1,j+1))/2.0d0)-((v(i,j-1)+v(i+1,j-1))/2.0d0))/2.0d0
		Th_r(i,j)=(((t(i,j+1)+t(i+1,j+1))/2.0d0)-((t(i,j-1)+t(i+1,j-1))/2.0d0))/2.0d0
		kh_r(i,j)=(((k(i,j+1)+k(i+1,j+1))/2.0d0)-((k(i,j-1)+k(i+1,j-1))/2.0d0))/2.0d0
		epsih_r(i,j)=(((epsi(i,j+1)+epsi(i+1,j+1))/2.0d0)-((epsi(i,j-1)+epsi(i+1,j-1))/2.0d0))/2.0d0

		myut_r(i,j)=Cmyu*rho_r(i,j)*(k_r(i,j)**2.0)/epsi_r(i,j)
	end do
end do

do j=1,JM-1
	do i=0,IM-1
		Sxx=2.0d0/3.0d0*(2.0d0*(ue_r(i,j)*xix_r(i,j)+uh_r(i,j)*ex_r(i,j))-(ve_r(i,j)*xiy_r(i,j)+vh_r(i,j)*ey_r(i,j)))
		Syy=2.0d0/3.0d0*(2.0d0*(ve_r(i,j)*xiy_r(i,j)+vh_r(i,j)*ey_r(i,j))-(ue_r(i,j)*xix_r(i,j)+uh_r(i,j)*ex_r(i,j)))
		Sxy=(ue_r(i,j)*xiy_r(i,j)+uh_r(i,j)*ey_r(i,j))+(ve_r(i,j)*xix_r(i,j)+vh_r(i,j)*ex_r(i,j))
		Syx=Sxy
		txx=myut_r(i,j)*Sxx-2.0d0/3.0d0*rho_r(i,j)*k_r(i,j)
		tyy=myut_r(i,j)*Syy-2.0d0/3.0d0*rho_r(i,j)*k_r(i,j)
		txy=myut_r(i,j)*Sxy
		tyx=txy

		R4=txx*u_r(i,j)+txy*v_r(i,j)+(myu/(Pr*(g-1.0d0)))*g*R*(Te_r(i,j)*xix_r(i,j)+Th_r(i,j)*ex_r(i,j))
		R5=myut_r(i,j)*(ke_r(i,j)*xix_r(i,j)+kh_r(i,j)*ex_r(i,j))/sigmak
		R6=myut_r(i,j)*(epsie_r(i,j)*xix_r(i,j)+epsih_r(i,j)*ex_r(i,j))/sigmaepsi
		S4=tyx*u_r(i,j)+tyy*v_r(i,j)+(myu/(Pr*(g-1.0d0)))*g*R*(Te_r(i,j)*xiy_r(i,j)+Th_r(i,j)*ey_r(i,j))
		S5=myut_r(i,j)*(ke_r(i,j)*xiy_r(i,j)+kh_r(i,j)*ey_r(i,j))/sigmak
		S6=myut_r(i,j)*(epsie_r(i,j)*xiy_r(i,j)+epsih_r(i,j)*ey_r(i,j))/sigmaepsi

		RRR(1,i,j)=0.0d0
		RRR(2,i,j)=(xix_r(i,j)*txx+xiy_r(i,j)*txy)/Jac_r(i,j)
		RRR(3,i,j)=(xix_r(i,j)*tyx+xiy_r(i,j)*tyy)/Jac_r(i,j)
		RRR(4,i,j)=(xix_r(i,j)*R4+xiy_r(i,j)*S4)/Jac_r(i,j)
		RRR(5,i,j)=(xix_r(i,j)*R5+xiy_r(i,j)*S5)/Jac_r(i,j)
		RRR(6,i,j)=(xix_r(i,j)*R6+xiy_r(i,j)*S6)/Jac_r(i,j)
	end do
end do
end subroutine viscosity_xi

!----------η方向粘性項----------
subroutine viscosity_eta
double precision,dimension(0:IM-1,0:JM-1)::xix_s,xiy_s,ex_s,ey_s,Jac_s,u_s,v_s,T_s,ue_s,ve_s,Te_s,uh_s,vh_s,Th_s
double precision,dimension(0:IM-1,0:JM-1)::k_s,epsi_s,rho_s,myut_s,ke_s,epsie_s,kh_s,epsih_s
double precision::R4,R5,R6,S4,S5,S6,Sxx,Syy,Sxy,Syx,txx,tyy,txy,tyx

do i=1,IM-1
	do j=0,JM-1
		xix_s(i,j)=(xix(i,j)+xix(i,j+1))/2.0d0
		xiy_s(i,j)=(xiy(i,j)+xiy(i,j+1))/2.0d0
		ex_s(i,j)=(ex(i,j)+ex(i,j+1))/2.0d0
		ey_s(i,j)=(ey(i,j)+ey(i,j+1))/2.0d0
		Jac_s(i,j)=(Jac(i,j)+Jac(i,j+1))/2.0d0
		u_s(i,j)=(u(i,j)+u(i,j+1))/2.0d0
		v_s(i,j)=(v(i,j)+v(i,j+1))/2.0d0
		T_s(i,j)=(T(i,j)+T(i,j+1))/2.0d0
		k_s(i,j)=(k(i,j)+k(i,j+1))/2.0d0
		epsi_s(i,j)=(epsi(i,j)+epsi(i,j+1))/2.0d0
		rho_s(i,j)=(Q(1,i,j)*Jac(i,j)+Q(1,i,j+1)*Jac(i,j+1))/(Jac(i,j)+Jac(i,j+1))*Jac_s(i,j)

		uh_s(i,j)=u(i,j+1)-u(i,j)
		vh_s(i,j)=v(i,j+1)-v(i,j)
		Th_s(i,j)=T(i,j+1)-T(i,j)
		kh_s(i,j)=k(i,j+1)-k(i,j)
		epsih_s(i,j)=epsi(i,j+1)-epsi(i,j)

		ue_s(i,j)=(((u(i+1,j)+u(i+1,j+1))/2.0d0)-((u(i-1,j)+u(i-1,j+1))/2.0d0))/2.0d0
		ve_s(i,j)=(((v(i+1,j)+v(i+1,j+1))/2.0d0)-((v(i-1,j)+v(i-1,j+1))/2.0d0))/2.0d0
		Te_s(i,j)=(((T(i+1,j)+T(i+1,j+1))/2.0d0)-((T(i-1,j)+T(i-1,j+1))/2.0d0))/2.0d0
		ke_s(i,j)=(((k(i+1,j)+k(i+1,j+1))/2.0d0)-((k(i-1,j)+k(i-1,j+1))/2.0d0))/2.0d0
		epsie_s(i,j)=(((epsi(i+1,j)+epsi(i+1,j+1))/2.0d0)-((epsi(i-1,j)+epsi(i-1,j+1))/2.0d0))/2.0d0

		myut_s(i,j)=Cmyu*rho_s(i,j)*(k_s(i,j)**2.0)/epsi_s(i,j)
	end do
end do

do i=1,IM-1
	do j=0,JM-1
		Sxx=2.0d0/3.0d0*(2.0d0*(ue_s(i,j)*xix_s(i,j)+uh_s(i,j)*ex_s(i,j))-(ve_s(i,j)*xiy_s(i,j)+vh_s(i,j)*ey_s(i,j)))
		Syy=2.0d0/3.0d0*(2.0d0*(ve_s(i,j)*xiy_s(i,j)+vh_s(i,j)*ey_s(i,j))-(ue_s(i,j)*xix_s(i,j)+uh_s(i,j)*ex_s(i,j)))
		Sxy=(ue_s(i,j)*xiy_s(i,j)+uh_s(i,j)*ey_s(i,j))+(ve_s(i,j)*xix_s(i,j)+vh_s(i,j)*ex_s(i,j))
		Syx=Sxy
		txx=myut_s(i,j)*Sxx-2.0d0/3.0d0*rho_s(i,j)*k_s(i,j)
		tyy=myut_s(i,j)*Syy-2.0d0/3.0d0*rho_s(i,j)*k_s(i,j)
		txy=myut_s(i,j)*Sxy
		tyx=txy
		
		R4=txx*u_s(i,j)+txy*v_s(i,j)+(myut_s(i,j)/(Pr*(g-1.0d0)))*g*R*(Te_s(i,j)*xix_s(i,j)+Th_s(i,j)*ex_s(i,j))
		R5=myut_s(i,j)*(ke_s(i,j)*xix_s(i,j)+kh_s(i,j)*ex_s(i,j))/sigmak
		R6=myut_s(i,j)*(epsie_s(i,j)*xix_s(i,j)+epsih_s(i,j)*ex_s(i,j))/sigmaepsi
		S4=tyx*u_s(i,j)+tyy*v_s(i,j)+(myut_s(i,j)/(Pr*(g-1.0d0)))*g*R*(Te_s(i,j)*xiy_s(i,j)+Th_s(i,j)*ey_s(i,j))
		S5=myut_s(i,j)*(ke_s(i,j)*xiy_s(i,j)+kh_s(i,j)*ey_s(i,j))/sigmak
		S6=myut_s(i,j)*(epsie_s(i,j)*xiy_s(i,j)+epsih_s(i,j)*ey_s(i,j))/sigmaepsi

		SSS(1,i,j)=0.0d0
		SSS(2,i,j)=(ex_s(i,j)*txx+ey_s(i,j)*txy)/Jac_s(i,j)
		SSS(3,i,j)=(ex_s(i,j)*tyx+ey_s(i,j)*tyy)/Jac_s(i,j)
		SSS(4,i,j)=(ex_s(i,j)*R4+ey_s(i,j)*S4)/Jac_s(i,j)
		SSS(5,i,j)=(ex_s(i,j)*R5+ey_s(i,j)*S5)/Jac_s(i,j)
		SSS(6,i,j)=(ex_s(i,j)*R6+ey_s(i,j)*S6)/Jac_s(i,j)
	end do
end do

end subroutine viscosity_eta

!----------生産項----------
subroutine production
double precision,dimension(0:IM,0:JM)::ue,ve,uh,vh
double precision::Sxx,Syy,Sxy,Syx,txx,tyy,txy,tyx,PRO
double precision,parameter::Cepsi1=1.44d0,Cepsi2=1.92d0
do j=1,JM-1
	do i=1,IM-1
		ue(i,j)=(u(i+1,j)-u(i-1,j))/2.0d0
		ve(i,j)=(v(i+1,j)-v(i-1,j))/2.0d0
		uh(i,j)=(u(i,j+1)-u(i,j-1))/2.0d0
		vh(i,j)=(v(i,j+1)-v(i,j-1))/2.0d0

		Sxx=2.0d0/3.0d0*(2.0d0*(ue(i,j)*xix(i,j)+uh(i,j)*ex(i,j))-(ve(i,j)*xiy(i,j)+vh(i,j)*ey(i,j)))
		Syy=2.0d0/3.0d0*(2.0d0*(ve(i,j)*xiy(i,j)+vh(i,j)*ey(i,j))-(ue(i,j)*xix(i,j)+uh(i,j)*ex(i,j)))
		Sxy=(ue(i,j)*xiy(i,j)+uh(i,j)*ey(i,j))+(ve(i,j)*xix(i,j)+vh(i,j)*ex(i,j))
		Syx=Sxy

		txx=myut(i,j)*Sxx-2.0d0/3.0d0*rho(i,j)*k(i,j)
		tyy=myut(i,j)*Syy-2.0d0/3.0d0*rho(i,j)*k(i,j)
		txy=myut(i,j)*Sxy
		tyx=txy

		PRO=txx*(ue(i,j)*xix(i,j)+uh(i,j)*ex(i,j))&
			&+txy*((ve(i,j)*xix(i,j)+vh(i,j)*ex(i,j))+(ue(i,j)*xiy(i,j)+uh(i,j)*ey(i,j)))&
			&+tyy*(ve(i,j)*xiy(i,j)+vh(i,j)*ey(i,j))

		PPP(1,i,j)=0.0d0
		PPP(2,i,j)=0.0d0
		PPP(3,i,j)=0.0d0
		PPP(4,i,j)=0.0d0
		PPP(5,i,j)=(PRO-rho(i,j)*epsi(i,j))/Jac(i,j)
		PPP(6,i,j)=(Cepsi1*epsi(i,j)/k(i,j)*PRO-Cepsi2*rho(i,j)*(epsi(i,j)**2.0)/k(i,j))/Jac(i,j)
	end do
end do
end subroutine production

!----------壁法則----------
subroutine wall_law
double precision::utau,ypp,yp,kp,epsip
double precision,parameter::kap=0.42d0,EEE=9.7d0

do i=5,IM
	yp=y(i,1)
	kp=k(i,1)
	epsip=epsi(i,1)

	utau=Cmyu*(kp**2.0)/(kap*yp*epsip)
	ypp=yp*utau/nyu(i,1)
	if(ypp>11.63d0)then
		utau=kap*u(i,1)/dlog(ypp*EEE)
	else
		utau=dsqrt(nyu(i,1)*u(i,1)/yp)
	end if

	k(i,1)=(utau**2.0)/dsqrt(Cmyu)
	epsi(i,1)=(utau**3.0)/(kap*yp)
	kwall(i,1)=k(i,1)
	epsiwall(i,1)=epsi(i,1)
end do
end subroutine wall_law

!------Runge-Kutta------
subroutine Runge_Kutta
do j=1,JM-1
	do i=1,IM-1
		do kk=1,6
			QQ(kk,i,j,rk)=QQ(kk,i,j,0)-(1.0d0/dble(5-rk))*ldt(i,j)*(-EE(kk,i-1,j)+EE(kk,i,j)&
							&-FF(kk,i,j-1)+FF(kk,i,j)+((RRR(kk,i-1,j)-RRR(kk,i,j)))+&
							&((SSS(kk,i,j-1)-SSS(kk,i,j)))-PPP(kk,i,j))
			Q(kk,i,j)=QQ(kk,i,j,rk)
		end do
	end do
end do
end subroutine Runge_Kutta

!------諸量の計算------
subroutine calculation
do j=1,JM-1
	do i=1,IM-1
		un(i,j)=u(i,j)
		vn(i,j)=v(i,j)
		rho(i,j)=Jac(i,j)*Q(1,i,j)
		u(i,j)=Jac(i,j)*Q(2,i,j)/rho(i,j)
		v(i,j)=Jac(i,j)*Q(3,i,j)/rho(i,j)
		k(i,j)=Q(5,i,j)/Q(1,i,j)
		epsi(i,j)=Q(6,i,j)/Q(1,i,j)
		ene(i,j)=Jac(i,j)*Q(4,i,j)
		enebar(i,j)=ene(i,j)/rho(i,j)-((u(i,j)**2.0+v(i,j)**2.0)/2.0d0)
		T(i,j)=enebar(i,j)*(g-1.0d0)/R
		p(i,j)=(g-1.0d0)*rho(i,j)*enebar(i,j)
	end do
end do

do j=1,JM-1 !流出
	slope=(x(IM,j)-x(IM-2,j))/(x(IM-2,j)-x(IM-1,j))
	rho(IM,j)=rho(IM-2,j)+slope*(rho(IM-2,j)-rho(IM-1,j))
	u(IM,j)=u(IM-2,j)+slope*(u(IM-2,j)-u(IM-1,j))
	v(IM,j)=v(IM-2,j)+slope*(v(IM-2,j)-v(IM-1,j))
	p(IM,j)=p(IM-2,j)+slope*(p(IM-2,j)-p(IM-1,j))
	ene(IM,j)=p(IM,j)/(g-1.0d0)+rho(IM,j)*((u(IM,j)**2.0)+(v(IM,j)**2.0))/2.0d0
	T(IM,j)=T(IM-2,j)+slope*(T(IM-2,j)-T(IM-1,j))
	k(IM,j)=k(IM-2,j)+slope*(k(IM-2,j)-k(IM-1,j))
	epsi(IM,j)=epsi(IM-2,j)+slope*(epsi(IM-2,j)-epsi(IM-1,j))
	
	Q(1,IM,j)=rho(IM,j)/Jac(IM,j)
	Q(2,IM,j)=rho(IM,j)*u(IM,j)/Jac(IM,j)
	Q(3,IM,j)=rho(IM,j)*v(IM,j)/Jac(IM,j)
	Q(4,IM,j)=ene(IM,j)/Jac(IM,j)
	Q(5,IM,j)=rho(IM,j)*k(IM,j)/Jac(IM,j)
	Q(6,IM,j)=rho(IM,j)*epsi(IM,j)/Jac(IM,j)
end do

do i=0,IM !上壁
	rho(i,JM)=rho(i,JM-1)
	u(i,JM)=u(i,JM-1)
	v(i,JM)=v(i,JM-1)
	p(i,JM)=p(i,JM-1)
	ene(i,JM)=ene(i,JM-1)
	T(i,JM)=T(i,JM-1)
	k(i,JM)=k(i,JM-1)
	epsi(i,JM)=epsi(i,JM-1)

	Q(1,i,JM)=rho(i,JM)/Jac(i,JM)
	Q(2,i,JM)=rho(i,JM)*u(i,JM)/Jac(i,JM)
	Q(3,i,JM)=rho(i,JM)*v(i,JM)/Jac(i,JM)
	Q(4,i,JM)=ene(i,JM)/Jac(i,JM)
	Q(5,i,JM)=rho(i,JM)*k(i,JM)/Jac(i,JM)
	Q(6,i,JM)=rho(i,JM)*epsi(i,JM)/Jac(i,JM)
end do

do i=0,IM !下壁
	rho(i,0)=rho(i,1)
	v(i,0)=0.0d0
	p(i,0)=p(i,1)
	T(i,0)=T(i,1)
	ene(i,0)=p(i,0)/(g-1.0d0)+rho(i,0)*((u(i,0)**2.0+v(i,0)**2.0))/2.0d0

	if(i<=5)then
		u(i,0)=u(i,1)
		k(i,0)=k(i,1)
		epsi(i,0)=epsi(i,1)
	else
		u(i,0)=0.0d0
		k(i,1)=kwall(i,1)
		k(i,0)=k(i,1)
		epsi(i,1)=epsiwall(i,1)
		epsi(i,0)=epsi(i,1)
	end if	
	
	Q(1,i,0)=rho(i,0)/Jac(i,0)
	Q(2,i,0)=rho(i,0)*u(i,0)/Jac(i,0)
	Q(3,i,0)=rho(i,0)*v(i,0)/Jac(i,0)
	Q(4,i,0)=ene(i,0)/Jac(i,0)
	Q(5,i,0)=rho(i,0)*k(i,0)/Jac(i,0)
	Q(6,i,0)=rho(i,0)*epsi(i,0)/Jac(i,0)
end do

du=0.0d0
dv=0.0d0
do j=1,JM
	do i=0,IM
		ddu=0.0d0
		ddv=0.0d0
		ddu=(abs(u(i,j)-un(i,j)))/abs(un(i,j))
		ddv=(abs(v(i,j)-vn(i,j)))/abs(vn(i,j))
		if(du<ddu)then
			du=ddu
		end if
		if(dv<ddv)then
			dv=ddv
		end if
	end do
end do

end subroutine calculation

!------局所時間刻み法------
subroutine local_time_step_method
double precision,parameter::CN=0.2
double precision,dimension(0:IM,0:JM)::uu,vv
double precision::lcxi,lce,ldxi,lde,lepsixi,lepsie,Lxi,Le

do j=0,JM
	do i=0,IM
		c(i,j)=dsqrt(g*r*T(i,j))
		uu(i,j)=xix(i,j)*u(i,j)+xiy(i,j)*v(i,j) !反変速度
		vv(i,j)=ex(i,j)*u(i,j)+ey(i,j)*v(i,j)

		!対流項
		lcxi=dabs(uu(i,j))+c(i,j)*dsqrt(xix(i,j)**2.0+xiy(i,j)**2.0)
		lce=dabs(vv(i,j))+c(i,j)*dsqrt(ex(i,j)**2.0+ey(i,j)**2.0)

		!拡散項
		ldxi=2.0d0*(myu(i,j)+myut(i,j))/rho(i,j)*(xix(i,j)**2.0+xiy(i,j)**2.0)
		lde=2.0d0*(myu(i,j)+myut(i,j))/rho(i,j)*(ex(i,j)**2.0+ey(i,j)**2.0)

		!流れの散逸項
		lepsixi=epsi(i,j)/k(i,j)
		lepsie=lepsixi

		Lxi=lcxi+ldxi+lepsixi
		Le=lce+lde+lepsie

		ldt(i,j)=CN*dmin1(1.0d0/Lxi,1.0d0/Le)
	end do
end do
end subroutine local_time_step_method


!----------y＋とu+計算・出力----------
subroutine output_yplus
double precision::utau,tauw
double precision,dimension(0:IM,0:JM)::yplus,uplus

do j=0,JM
	do i=5,IM
		tauw=myu(i,0)*(u(i,1)-u(i,0))/(y(i,1)-y(i,0))
		utau=dsqrt(tauw/rho(i,0))
		uplus(i,j)=u(i,j)/utau
		yplus(i,j)=y(i,j)*utau/nyu(i,j)
	end do
end do
open(10,file='assyuku8_1_yu.dat',status='replace')
do j=0,JM
	write(10,*) yplus(45,j),uplus(45,j)
end do
close(10)
end subroutine output_yplus

!----------計算結果出力----------
subroutine output
open(1,file='assyuku8_1.dat',status='replace')
do j=0,JM
	do i=0,IM
		write(1,*) x(i,j),y(i,j),u(i,j),v(i,j),p(i,j),T(i,j)
	end do
end do
close(1)

open(11,file='assyuku8_1.fld',status='replace')
write(11,'(a)')'# AVS field file'
write(11,'(a)')'ndIM=2'
write(11,'(a,i0)')'dIM1=',IM+1
write(11,'(a,i0)')'dIM2=',JM+1
write(11,'(a)')'nspace=2'
write(11,'(a)')'veclen=4'
write(11,'(a)')'data=float'
write(11,'(a)')'field=irregular'
write(11,'(a)')'label=u,v,p,t'
write(11,'(a)')'coord 1 file=assyuku8_1.dat filetype=ascii skip=0 offset=0 stride=6'
write(11,'(a)')'coord 2 file=assyuku8_1.dat filetype=ascii skip=0 offset=1 stride=6'
write(11,'(a)')'variable 1 file=assyuku8_1.dat filetype=ascii skip=0 offset=2 stride=6'
write(11,'(a)')'variable 2 file=assyuku8_1.dat filetype=ascii skip=0 offset=3 stride=6'
write(11,'(a)')'variable 3 file=assyuku8_1.dat filetype=ascii skip=0 offset=4 stride=6'
write(11,'(a)')'variable 4 file=assyuku8_1.dat filetype=ascii skip=0 offset=5 stride=6'
close(11)
end subroutine output

!--------------------

end program assyuku8_1