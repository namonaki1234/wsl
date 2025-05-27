program kadai5
implicit none

double precision,parameter::h=10.0d0
double precision,parameter::a=1.2d0
double precision,parameter::b=(a+1.0d0)/(a-1.0d0)
integer,parameter::im=30
integer,parameter::jm=30
double precision,parameter::r=2.871d2       !ガス定数
double precision,parameter::g=1.4d0     !比熱比
double precision,parameter::dt=1.0d-6
double precision,parameter::pi=dacos(0.0d0)*2.0d0
double precision,parameter::betas=3.0d1*pi/1.80d2
double precision,parameter::ii=1.0d0
double precision,parameter::eps=1.0d-4
double precision,dimension(0:im,0:jm)::x
double precision,dimension(0:im,0:jm)::y
double precision,dimension(0:im,0:jm)::rho=1.2d0,p=0.0d0,u=0.0d0,v=0.0d0,&
energy=0.0d0,eb,t,c,vn1,vn2,vt,v2,un,vn
double precision,dimension(0:im,0:jm)::xix,xiy,etax,etay,xxi,yxi,xeta,yeta,&
jj,uu,vv
double precision,dimension(1:4,0:im,0:jm)::q,ee,ff
double precision,dimension(1:4,0:im,0:jm,0:5)::qq
double precision,dimension(4,4,0:im)::RR,RI
double precision,dimension(4,0:im)::EIGM,FIGM,GG,PHIM,E,F,EIG,FIG
double precision,dimension(4,-1:im)::ALPHA
double precision,dimension(4)::D
integer,dimension(0:im,0:jm)::flag
integer i,j,n,mm,k
double precision dx,dy,m,jmm,tanb,sinb,p1,rho1,cosb,theta,time,slope,veps,dd,&
AKX,AKY,AJACM,UM,VM,DUM,CM,BETA,AKXT,THIT,RWL,RWR,sij,GAMMA,UCONT,AKYT,DELTA,&
PHI


!flag=1:after wave  flag=0:before wave


!main 
print*,"*****START CALCULATION*****"
time=0.0d0
call mesh
call initialsetting
call metrix

do n=0,1000000  !計算回数
do j=0,jm
   do i=0,im
      do k=1,4
QQ(k,i,j,0)=Q(k,i,j)
      end do
   end do
end do

do mm=1,4
call xi_TVD
call eta_TVD
call runge_kutta
call calculate
end do
print*,time,dd
time=time+dt
if(dd<eps)exit
end do
print*,"*****END CALCULATION*****"



open(17,file='kadai5.dat',status='replace')
   do j=0,jm
      do i=0,im
         write(17,'(f13.6,$)') x(i,j)
         write(17,'(f13.6,$)') y(i,j)
         write(17,'(f13.6,$)') u(i,j)
         write(17,'(f13.6,$)') v(i,j)
         write(17,'(f13.3,$)') p(i,j)
         write(17,'(f13.6)') t(i,j)
      end do
   end do
close(17)

open(18,file='kadai5.fld',status='replace')
write(18,'(a)')'# AVS field file'
write(18,'(a)')'ndim=2'
write(18,'(a,i0)')'dim1=',im+1
write(18,'(a,i0)')'dim2=',jm+1
write(18,'(a)')'nspace=2'
write(18,'(a)')'veclen=4'
write(18,'(a)')'data=float'
write(18,'(a)')'field=irregular'
write(18,'(a)')'label=u,v,p,t'
write(18,'(a)')'variable 1 file=kadai5.dat filetype=ascii skip=0 offset=2 stride=6'
write(18,'(a)')'variable 2 file=kadai5.dat filetype=ascii skip=0 offset=3 stride=6'
write(18,'(a)')'variable 3 file=kadai5.dat filetype=ascii skip=0 offset=4 stride=6'
write(18,'(a)')'variable 4 file=kadai5.dat filetype=ascii skip=0 offset=5 stride=6'
write(18,'(a)')'coord 1 file=kadai5.dat filetype=ascii skip=0 offset=0 stride=6'
write(18,'(a)')'coord 2 file=kadai5.dat filetype=ascii skip=0 offset=1 stride=6'
close(18)


!スキーム内変数まとめ
!AKYT=k~_y
!CM=c
!AKXT=k~_x
!PHI=φ
!G=γ
!THIT=θ~
!BETA=β
!EIGM=α_j+1/2,k::特性変数
!ALPHA=α
!GG=g_j,k::制限関数
!EE=E~_j+1/2,k::ξ方向数値流束
!PHIM=φ_j+1/2,k
!XIX=ξ_x
!XIY=ξ_y


!ここまで

contains

subroutine mesh !不等間隔格子生成

do j=0,jm
   do i=0,im
      dx=(3.0d0*h)/dble(im)
      x(i,j)=dble(i)*dx
      jmm=dble(j)/(dble(jm))
      y(i,j)=h*(b**jmm-1.0d0)/(b-1.0d0)
   end do
end do

tanb=dtan(betas)

do j=0,jm
   do i=0,im
      if((-tanb*x(i,j)+y(i,jm-2))<=y(i,j))then
         flag(i,j)=1
      else
         flag(i,j)=0
      end if
   end do
end do

end subroutine mesh

subroutine initialsetting !初期設定
do j=0,jm
   do i=0,im
      t(i,j)=2.930d2
      c(i,j)=dsqrt(g*r*t(i,j)) !音速
   end do
end do
m=2.9d0                         !マッハ数

do j=0,jm
   u(0,j)=m*c(0,j)
   u(1,j)=u(0,j)
end do

do j=0,jm
   do i=2,im
      slope=(x(i,j)-x(i-2,j))/(x(i-2,j)-x(i-1,j))
      u(i,j)=u(i-2,j)+slope*(u(i-2,j)-u(i-1,j))
   end do
end do


do j=0,jm
   do i=0,im
      if(flag(i,j)==0)then
p(i,j)=rho(i,j)*r*t(i,j)
energy(i,j)=rho(i,j)*((r*t(i,j)/(g-1.0d0))+((u(i,j)**2.0)+v(i,j)**2.0)/2.0d0)
      end if
   end do
end do

sinb=dsin(betas)
p1=1.2d0*r*293.0d0
rho1=1.2d0
do j=0,jm
   do i=0,im
      if(flag(i,j)==1)then
t(i,j)=293.0d0*(2.0d0*g*(m**2.0)*(sinb**2.0)-(g-1.0d0))*&
   ((g-1.0d0)*(m**2.0)*(sinb**2.0)+2.0d0)/(((g+1.0d0)**2.0)*(m**2.0)*&
   (sinb**2.0))
p(i,j)=p1*(2.0d0*g*(m**2.0)*(sinb**2.0)-(g-1.0d0))/(g+1.0d0)
rho(i,j)=rho1*((g+1.0d0)*(m**2.0)*(sinb**2.0))/((g-1.0d0)*&
   (m**2.0)*(sinb**2.0)+2.0d0)
vn1(i,j)=u(i,j)*sinb
vt(i,j)=vn1(i,j)/tanb
vn2(i,j)=vn1(i,j)*((g-1.0d0)*(m**2.0)*(sinb**2.0)+2.0d0)/((g+1.0d0)*&
   (m**2.0)*(sinb**2.0))
v2(i,j)=dsqrt((vn2(i,j)**2.0)+(vt(i,j)**2.0))
theta=betas-datan((vn2(i,j)/vt(i,j)))
u(i,j)=v2(i,j)*dcos(-theta)
v(i,j)=v2(i,j)*dsin(-theta)
energy(i,j)=rho(i,j)*((r*t(i,j)/(g-1.0d0))+((u(i,j)**2.0)+v(i,j)**2.0)/2.0d0)
      end if
   end do
end do
end subroutine initialsetting

subroutine metrix !座標変換

do j=1,jm-1
   do i=0,im
      yeta(i,j)=((y(i,j+1)-y(i,j-1))/2.0d0)
      xeta(i,j)=((x(i,j+1)-x(i,j-1))/2.0d0)
   end do
end do
do j=0,jm
   do i=1,im-1
      yxi(i,j)=((y(i+1,j)-y(i-1,j))/2.0d0)
      xxi(i,j)=((x(i+1,j)-x(i-1,j))/2.0d0)
   end do
end do
do i=0,im
   yeta(i,0)=(-3.0d0*y(i,0)+4.0d0*y(i,1)-y(i,2))/(2.0d0*ii)
   yeta(i,jm)=(3.0d0*y(i,jm)-4.0d0*y(i,jm-1)+y(i,jm-2))/(2.0d0*ii)
   xeta(i,0)=(-3.0d0*x(i,0)+4.0d0*x(i,1)-x(i,2))/(2.0d0**ii)
   xeta(i,jm)=(3.0d0*x(i,jm)-4.0d0*x(i,jm-1)+x(i,jm-2))/(2.0d0*ii)
end do
do j=0,jm
   yxi(0,j)=(-3.0d0*y(0,j)+4.0d0*y(1,j)-y(2,j))/(2.0d0*ii)
   yxi(im,j)=(3.0d0*y(im,j)-4.0d0*y(im-1,j)+y(im-2,j))/(2.0d0*ii)
   xxi(0,j)=(-3.0d0*x(0,j)+4.0d0*x(1,j)-x(2,j))/(2.0d0*ii)
   xxi(im,j)=(3.0d0*x(im,j)-4.0d0*x(im-1,j)+x(im-2,j))/(2.0d0*ii)
end do

do j=0,jm
   do i=0,im
      jj(i,j)=1.0d0/(xxi(i,j)*yeta(i,j)-yxi(i,j)*xeta(i,j))
      xix(i,j)=yeta(i,j)*jj(i,j)
      xiy(i,j)=-xeta(i,j)*jj(i,j)
      etax(i,j)=-yxi(i,j)*jj(i,j)
      etay(i,j)=xxi(i,j)*jj(i,j)
   end do
end do

do j=0,jm
   do i=0,im
      q(1,i,j)=rho(i,j)/jj(i,j)
      q(2,i,j)=(rho(i,j)*u(i,j))/jj(i,j)
      q(3,i,j)=(rho(i,j)*v(i,j))/jj(i,j)
      q(4,i,j)=energy(i,j)/jj(i,j)
   end do
end do

end subroutine metrix



subroutine xi_TVD
!*****ξ方向*****
      DO  J=1,JM-1
       DO  I=0,IM-1

        RWL=dSQRT(Q(1,I,J)*jj(I,J))&
            /(dSQRT(Q(1,I+1,J)*jj(I+1,J))+dSQRT(Q(1,I,J)*jj(I,J)))
        RWR=dSQRT(Q(1,I+1,J)*jj(I+1,J))&
            /(dSQRT(Q(1,I+1,J)*jj(I+1,J))+dSQRT(Q(1,I,J)*jj(I,J)))
        AKX=0.5d0*(XIX(I,J)+XIX(I+1,J))      !k_x
        AKY=0.5d0*(XIY(I,J)+XIY(I+1,J))      !k_y
        AJACM=0.5d0*(jj(I,J)+jj(I+1,J))

        UM=RWL*U(I,J)+RWR*U(I+1,J)
        VM=RWL*V(I,J)+RWR*V(I+1,J)
        DUM=RWL*Q(4,I,J)/Q(1,I,J)+RWR*Q(4,I+1,J)/Q(1,I+1,J)
        CM=dSQRT(G*(G-1.0d0)*ABS(DUM-0.5d0*(UM**2.0+VM**2.0)))

!       NOMENCLATURE
        PHI=0.5d0*(G-1.0d0)*(UM**2.0+VM**2.0)
        BETA=1.0d0/(2.0d0*CM**2.0)
        AKXT=AKX/dSQRT(AKX**2.0+AKY**2.0)
        AKYT=AKY/dSQRT(AKX**2.0+AKY**2.0)
        THIT=AKXT*UM+AKYT*VM
!
!       RIGHT ENGEN-VECTORS     P.88 R_k
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
        RR(4,3,I)=(PHI+CM**2)/(G-1.0d0)+CM*THIT
        RR(4,4,I)=(PHI+CM**2)/(G-1.0d0)-CM*THIT
!
!       INVERS OR RIGHT EIGEN-VECTORS    P.88 R_k^-1
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
!
!       ENGEN-VALUES AT INTERMEDIATE                      P.87  固有値
        EIGM(1,I)=AKX*UM+AKY*VM
        EIGM(2,I)=EIGM(1,I)
        EIGM(3,I)=EIGM(1,I)+CM*dSQRT(AKX**2.0+AKY**2.0)
        EIGM(4,I)=EIGM(1,I)-CM*dSQRT(AKX**2.0+AKY**2.0)
!
!       ALPHA
        DO K=1,4
        D(K)=(Q(K,I+1,J)*jj(I+1,J)-Q(K,I,J)*jj(I,J))/AJACM
end do
!
        DO K=1,4
        ALPHA(K,I)=RI(K,1,I)*D(1)+RI(K,2,I)*D(2)+RI(K,3,I)*D(3)+RI(K,4,I)*D(4)
end do
end do
!
        DO K=1,4
        ALPHA(K,-1)=ALPHA(K,0)
        ALPHA(K,IM)=ALPHA(K,IM-1)
end do
!
        DO I=0,IM
        DUM=dSQRT(G*R*T(I,J)*(XIX(I,J)**2+XIY(I,J)**2))
        EIG(1,I)=XIX(I,J)*U(I,J)+XIY(I,J)*V(I,J)
        EIG(2,I)=EIG(1,I)
        EIG(3,I)=EIG(1,I)+DUM
        EIG(4,I)=EIG(1,I)-DUM
end do
!
!       CAL LIMITER FUNCTION G
        DO K=1,4
        DO I=0,IM
!
!       MINMOD LIMITER         制限関数
        sij=dSIGN(1.0d0,ALPHA(K,I))
        GG(K,I)=sij*dMAX1(0.0,dMIN1(dabs(ALPHA(K,I)),sij*ALPHA(K,I-1)))
end do
end do

        DO K=1,4
        DO I=0,IM-1

!       GAMMA
        DELTA=1.0d-10
!        DELTA=dMAX1(0.0,EIGM(K,I)-EIG(K,I),EIG(K,I+1)-EIGM(K,I))
        IF(ALPHA(K,I).NE.0.0)THEN
        GAMMA=0.5d0*FPSI(EIGM(K,I),DELTA)*(GG(K,I+1)-GG(K,I))/ALPHA(K,I)
        ELSE
        GAMMA=0.0d0
        ENDIF

!       PHI
        PHIM(K,I)=0.5d0*FPSI(EIGM(K,I),DELTA)*(GG(K,I+1)+GG(K,I))&
                   -FPSI(EIGM(K,I)+GAMMA,DELTA)*ALPHA(K,I)
end do
end do

!       CONVECTION COMPORNENTS    P.83 E^
        DO I=0,IM
        UCONT=XIX(I,J)*U(I,J)+XIY(I,J)*V(I,J)
        E(1,I)=Q(1,I,J)*UCONT
        E(2,I)=Q(2,I,J)*UCONT+P(I,J)*XIX(I,J)/jj(i,j)
        E(3,I)=Q(3,I,J)*UCONT+P(I,J)*XIY(I,J)/jj(i,j)
        E(4,I)=Q(4,I,J)*UCONT+P(I,J)*UCONT/jj(i,j)
end do

!       XI-DIRECTION CONVECTION FLUX　　　　　E~_j+1/2,k
        DO K=1,4
        DO I=0,IM-1
        EE(K,I,J)=0.5d0*(E(K,I)+E(K,I+1)&
                  +RR(K,1,I)*PHIM(1,I)+RR(K,2,I)*PHIM(2,I)&
                  +RR(K,3,I)*PHIM(3,I)+RR(K,4,I)*PHIM(4,I))
end do
end do

end do

      RETURN
end subroutine xi_TVD

!*****FPSI関数の定義*****           Ψ(z)の計算
double precision FUNCTION FPSI(Z,DELTA)
   double precision Z,DELTA
      IF(dabs(Z).GE.DELTA)THEN       !.ge. >=
        FPSI=dabs(Z)
      ELSE
        FPSI=0.5d0*(Z**2.0+DELTA**2.0)/DELTA
      END IF
      RETURN
      END function

subroutine eta_TVD
!*****η方向*****
      DO I=1,IM-1
       do J=0,JM-1

        RWL=dSQRT(Q(1,I,J)*jj(I,J))&
            /(dSQRT(Q(1,I,J+1)*jj(I,J+1))+dSQRT(Q(1,I,J)*jj(I,J)))
        RWR=dSQRT(Q(1,I,J+1)*jj(I,J+1))&
            /(dSQRT(Q(1,I,J+1)*jj(I,J+1))+dSQRT(Q(1,I,J)*jj(I,J)))
        AKX=0.5d0*(etaX(I,J)+etaX(I,J+1))      !k_x
        AKY=0.5d0*(etaY(I,J)+etaY(I,J+1))      !k_y
        AJACM=0.5d0*(jj(I,J)+jj(I,J+1))

        UM=RWL*U(I,J)+RWR*U(I,J+1)
        VM=RWL*V(I,J)+RWR*V(I,J+1)
        DUM=RWL*Q(4,I,J)/Q(1,I,J)+RWR*Q(4,I,J+1)/Q(1,I,J+1)
        CM=dSQRT(G*(G-1.0d0)*dabs(DUM-0.5d0*(UM**2.0+VM**2.0)))

!       NOMENCLATURE
        PHI=0.5d0*(G-1.0d0)*(UM**2.0+VM**2.0)
        BETA=1.0d0/(2.0d0*CM**2.0)
        AKXT=AKX/dSQRT(AKX**2.0+AKY**2.0)
        AKYT=AKY/dSQRT(AKX**2.0+AKY**2.0)
        THIT=AKXT*UM+AKYT*VM
!
!       RIGHT ENGEN-VECTORS     P.88 R_k
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
!
!       INVERS OR RIGHT EIGEN-VECTORS    P.88 R_k^-1
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
!
!       ENGEN-VALUES AT INTERMEDIATE                      P.87  固有値
        FIGM(1,J)=AKX*UM+AKY*VM
        FIGM(2,J)=FIGM(1,J)
        FIGM(3,J)=FIGM(1,J)+CM*dSQRT(AKX**2.0+AKY**2.0)
        FIGM(4,J)=FIGM(1,J)-CM*dSQRT(AKX**2.0+AKY**2.0)
!
!       ALPHA
        DO K=1,4
        D(K)=(Q(K,I,J+1)*jj(I,J+1)-Q(K,I,J)*jj(I,J))/AJACM
end do
!
        DO K=1,4
        ALPHA(K,J)=RI(K,1,J)*D(1)+RI(K,2,J)*D(2)+RI(K,3,J)*D(3)+RI(K,4,J)*D(4)
end do
end do
!
        DO K=1,4
        ALPHA(K,-1)=ALPHA(K,0)
        ALPHA(K,JM)=ALPHA(K,JM-1)
end do
!
        DO J=0,JM
        DUM=dSQRT(G*R*T(I,J)*(etaX(I,J)**2.0+etaY(I,J)**2.0))
        FIG(1,J)=etaX(I,J)*U(I,J)+etaY(I,J)*V(I,J)
        FIG(2,J)=FIG(1,J)
        FIG(3,J)=FIG(1,J)+DUM
        FIG(4,J)=FIG(1,J)-DUM
end do
!
!       CAL LIMITER FUNCTION G
        DO K=1,4
        DO J=0,JM
!
!       MINMOD LIMITER         制限関数
        sij=dSIGN(1.0d0,ALPHA(K,J))
        GG(K,J)=sij*dMAX1(0.0d0,dMIN1(dabs(ALPHA(K,J)),sij*ALPHA(K,J-1)))
end do
end do

        DO K=1,4
        DO J=0,JM-1

!       GAMMA
        DELTA=1.0d-10
!        DELTA=dMAX1(0.0,EIGM(K,I)-EIG(K,I),EIG(K,I+1)-EIGM(K,I))
        IF(ALPHA(K,J).NE.0.0)THEN
        GAMMA=0.5d0*FPSI(FIGM(K,J),DELTA)*(GG(K,J+1)-GG(K,J))/ALPHA(K,J)
        ELSE
        GAMMA=0.0d0
        ENDIF

!       PHI
        PHIM(K,J)=0.5d0*FPSI(FIGM(K,J),DELTA)*(GG(K,J+1)+GG(K,J))&
                   -FPSI(FIGM(K,J)+GAMMA,DELTA)*ALPHA(K,J)
end do
end do


!       CONVECTION COMPORNENTS    P.83 F^
        DO J=0,JM
        UCONT=etaX(I,J)*U(I,J)+etaY(I,J)*V(I,J)
        F(1,J)=Q(1,I,J)*UCONT
        F(2,J)=Q(2,I,J)*UCONT+P(I,J)*etaX(I,J)/jj(i,j)
        F(3,J)=Q(3,I,J)*UCONT+P(I,J)*etaY(I,J)/jj(i,j)
        F(4,J)=Q(4,I,J)*UCONT+P(I,J)*UCONT/jj(i,j)
end do

!       XI-DIRECTION CONVECTION FLUX　　　　　F~_j+1/2,k
        DO K=1,4
        DO J=0,JM-1
        FF(K,I,J)=0.5d0*(F(K,J)+F(K,J+1)&
                  +RR(K,1,J)*PHIM(1,J)+RR(K,2,J)*PHIM(2,J)&
                  +RR(K,3,J)*PHIM(3,J)+RR(K,4,J)*PHIM(4,J))
end do
end do

end do

      RETURN
end subroutine eta_TVD

subroutine runge_kutta

do j=1,jm-1
   do i=1,im-1
      do k=1,4
QQ(k,i,j,mm)=QQ(k,i,j,0)-(1.0d0/dble(5-mm))*dt*(-EE(k,i-1,j)+EE(k,i,j)-&
            FF(k,i,j-1)+FF(k,i,j))
q(k,i,j)=qq(k,i,j,mm)
      end do
   end do
end do


end subroutine runge_kutta

subroutine calculate    !諸量の計算
un=u
vn=v

do j=1,jm-1
   do i=1,im-1
      rho(i,j)=jj(i,j)*q(1,i,j)
      u(i,j)=jj(i,j)*q(2,i,j)/rho(i,j)
      v(i,j)=jj(i,j)*q(3,i,j)/rho(i,j)
      energy(i,j)=jj(i,j)*q(4,i,j)
      eb(i,j)=energy(i,j)/rho(i,j)-(((u(i,j)**2.0)+v(i,j)**2.0)/2.0d0)
      t(i,j)=eb(i,j)*(g-1.0d0)/r
      p(i,j)=(g-1.0d0)*rho(i,j)*eb(i,j)
   end do
end do

do i=0,im
   slope=(y(i,jm)-y(i,jm-2))/(y(i,jm-2)-y(i,jm-1))
   rho(i,jm)=rho(i,jm-2)+slope*(rho(i,jm-2)-rho(i,jm-1))
   u(i,jm)=u(i,jm-2)+slope*(u(i,jm-2)-u(i,jm-1))
   v(i,jm)=v(i,jm-2)+slope*(v(i,jm-2)-v(i,jm-1))
   p(i,jm)=p(i,jm-2)+slope*(p(i,jm-2)-p(i,jm-1))
   t(i,jm)=t(i,jm-2)+slope*(t(i,jm-2)-t(i,jm-1))
energy(i,jm)=p(i,jm)/(g-1.0d0)+rho(i,jm)*((u(i,jm)**2.0)+(v(i,jm)**2.0))/2.0d0

   q(1,i,jm)=rho(i,jm)/jj(i,jm)
   q(2,i,jm)=rho(i,jm)*u(i,jm)/jj(i,jm)
   q(3,i,jm)=rho(i,jm)*v(i,jm)/jj(i,jm)
   q(4,i,jm)=energy(i,jm)/jj(i,jm)

    slope=(y(i,0)-y(i,2))/(y(i,2)-y(i,1))
    rho(i,0)=rho(i,2)+slope*(rho(i,2)-rho(i,1))
    p(i,0)=p(i,2)+slope*(p(i,2)-p(i,1))
    u(i,0)=u(i,2)+slope*(u(i,2)-u(i,1))
    v(i,0)=0.0d0
    t(i,0)=t(i,2)+slope*(t(i,2)-t(i,1))
   energy(i,0)=p(i,0)/(g-1.0d0)+rho(i,0)*((u(i,0)**2.0+v(i,0)**2.0))/2.0d0

   q(1,i,0)=rho(i,0)/jj(i,0)
   q(2,i,0)=rho(i,0)*u(i,0)/jj(i,0)
   q(3,i,0)=rho(i,0)*v(i,0)/jj(i,0)
   q(4,i,0)=energy(i,0)/jj(i,0)
end do

do j=1,jm-1
 slope=(x(im,j)-x(im-2,j))/(x(im-2,j)-x(im-1,j))
   rho(im,j)=rho(im-2,j)+slope*(rho(im-2,j)-rho(im-1,j))
   p(im,j)=p(im-2,j)+slope*(p(im-2,j)-p(im-1,j))
   u(im,j)=u(im-2,j)+slope*(u(im-2,j)-u(im-1,j))
   v(im,j)=v(im-2,j)+slope*(v(im-2,j)-v(im-1,j))
   t(im,j)=t(im-2,j)+slope*(t(im-2,j)-t(im-1,j))
energy(im,j)=p(im,j)/(g-1.0d0)+rho(im,j)*((u(im,j)**2.0)+(v(im,j)**2.0))/2.0d0

   q(1,im,j)=rho(im,j)/jj(im,j)
   q(2,im,j)=rho(im,j)*u(im,j)/jj(im,j)
   q(3,im,j)=rho(im,j)*v(im,j)/jj(im,j)
   q(4,im,j)=energy(im,j)/jj(im,j)
end do

dd=0.0d0
do j=1,jm
   do i=0,im
      veps=(abs(u(i,j)-un(i,j)))/abs(un(i,j))
      if(dd<veps)then
         dd=veps
      end if
      veps=(abs(v(i,j)-vn(i,j)))/abs(vn(i,j))
      if(dd<veps)then
         dd=veps
      end if
   end do
end do

end subroutine calculate



end program kadai5
