program assyuku5_1
implicit none

!------初期条件------
double precision,parameter :: H = 10.0d0
double precision,parameter :: dx = 3.0d0*H/30.0d0
double precision,parameter :: b = 1.2d0 !公比b
double precision,parameter :: bt = 30.0d0/180.0d0*dacos(-1.0d0)
double precision :: a,ix
integer,parameter :: IM = 30,JM = 30
double precision,parameter :: R = 287.1d0,g = 1.4d0,dt = 1.0d-6 !マッハ数,温度,ガス定数,比熱比,密度,時間刻み
integer :: i,j,n,k,rk
double precision :: M1,rho1,rho2,T1,T2,p1,p2,u1,u2,v1,v2,ene1,ene2,vt,theta,vv1,vv2,vn1,vn2
double precision,dimension(0:IM,0:JM) :: xh,yh,xe,ye
double precision,dimension(0:IM,0:JM) :: xix,xiy,ex,ey,Jac
double precision,dimension(0:IM,0:JM) :: M,M2
double precision,dimension(0:IM,0:JM) :: p = 0.0d0
double precision,dimension(0:IM,0:JM) :: ene = 0.0d0,enebar
double precision,dimension(0:IM,0:JM) :: T
double precision,dimension(0:IM,0:JM) :: u = 0.0d0,un
double precision,dimension(0:IM,0:JM) :: v = 0.0d0,vn
double precision,dimension(0:IM,0:JM) :: rho
double precision,dimension(0:IM,0:JM) :: x,y
double precision,dimension(1:4,0:IM,0:JM,0:4) :: QQ
double precision,dimension(1:4,0:IM,0:JM) :: Q,EE,FF
double precision,parameter :: EPS = 1.0d-6
double precision :: du,dv,ddu,ddv,slope

!------メイン計算------
call mesh
call coordinate
call conditions
do n = 1,100000
  do j = 0,JM
    do i = 0,IM
      do k = 1,4
        QQ(k,i,j,0) = Q(k,i,j)
      end do
    end do
  end do

  do rk = 1,4 !Runge-Kuttaの段階数
    call tvd_xi
    call tvd_eta
    call Runge_Kutta
    call calquation
    call boundary
  end do

  print*,n
  if (dv<EPS) exit
end do

!グラフ作成
open (1,file='assyuku5_1.dat',status='replace')
do j = 0,JM
  do i = 0,IM
    write (1,*) x(i,j),y(i,j),u(i,j),v(i,j),T(i,j),p(i,j)
  end do
end do
close (1)

!------------------------------------------------------------
contains

!------格子生成------
subroutine mesh
  do i = 0,IM
    do j = 0,JM
      x(i,j) = dx*real(i)
    end do
  end do

  do j = 0,JM
    do i = 0,IM
      a = (b+1.0d0)/(b-1.0d0)
      ix = dble(j)/JM
      y(i,j) = H*(a**ix-1.0d0)/(a-1.0d0)
    end do
  end do
end subroutine mesh

!------座標変換(境界以外)------
subroutine coordinate
  do i = 0,IM
    do j = 1,JM-1
      xh(i,j) = (x(i,j+1)-x(i,j-1))/2.0d0
      yh(i,j) = (y(i,j+1)-y(i,j-1))/2.0d0
    end do
  end do

  do j = 0,JM
    do i = 1,IM-1
      xe(i,j) = (x(i+1,j)-x(i-1,j))/2.0d0
      ye(i,j) = (y(i+1,j)-y(i-1,j))/2.0d0
    end do
  end do

!------座標変換（境界）------
  do j = 0,JM
    xe(0,j) = (-3.0d0*x(0,j)+4.0d0*x(1,j)-x(2,j))/2.0d0
    xe(IM,j) = (3.0d0*x(IM,j)-4.0d0*x(IM-1,j)+x(IM-2,j))/2.0d0
    ye(0,j) = (-3.0d0*y(0,j)+4.0d0*y(1,j)-y(2,j))/2.0d0
    ye(IM,j) = (3.0d0*y(IM,j)-4.0d0*y(IM-1,j)+y(IM-2,j))/2.0d0
  end do

  do i = 0,IM
    xh(i,0) = (-3.0d0*x(i,0)+4.0d0*x(i,1)-x(i,2))/2.0d0
    xh(i,JM) = (3.0d0*x(i,JM)-4.0d0*x(i,JM-1)+x(i,JM-2))/2.0d0
    yh(i,0) = (-3.0d0*y(i,0)+4.0d0*y(i,1)-y(i,2))/2.0d0
    yh(i,JM) = -(-3.0d0*y(i,JM)+4.0d0*y(i,JM-1)-y(i,JM-2))/2.0d0
  end do

!------ヤコビアン計算------
  do j = 0,JM
    do i = 0,IM
      Jac(i,j) = 1.0d0/(xe(i,j)*yh(i,j)-ye(i,j)*xh(i,j))
      xix(i,j) = Jac(i,j)*yh(i,j)
      xiy(i,j) = Jac(i,j)*xh(i,j)
      ex(i,j) = -Jac(i,j)*ye(i,j)
      ey(i,j) = Jac(i,j)*xe(i,j)
    end do
  end do
end subroutine coordinate

!------諸条件------
subroutine conditions
  do j = 0,JM
    do i = 0,IM
      M1 = 2.9d0
      rho1 = 1.2d0
      T1 = 293.0d0
      u1 = M1*dsqrt(g*R*T1)
      v1 = 0.0d0
      vv1 = u1
      vn1 = u1*dsin(bt)
      vt = vn1/dtan(bt)
      p1 = rho1*R*T1
      ene1 = rho1*((R*T1/(g-1.0d0))+((u1**2.0+v1**2.0)/2.0d0))

      vn2 = vn1*(((g-1.0d0)*(M1**2)*(dsin(bt))**2)+2.0d0)/((g+1.0d0)*(M1**2)*(dsin(bt))**2)
      vv2 = dsqrt(vn2**2.0+vt**2.0)
      theta = bt-datan(vn2/vt)
      u2 = vv2*dcos(theta)
      v2 = vv2*dsin(-theta)
      T2 = T1*((2.0d0*g*(M1**2)*(dsin(bt))**2)-(g-1.0d0))*(((g-1.0d0)*(M1**2)*(dsin(bt))**2)+2.0d0)/&
             &(((g+1.0d0)**2)*(M1**2)*((dsin(bt))**2))
      rho2 = rho1*vn1/vn2
      p2 = p1*((2.0d0*g*(M1**2)*(dsin(bt))**2)-(g-1.0d0))/(g+1.0d0)
      ene2 = rho2*((R*T2/(g-1.0d0))+((u2**2.0+v2**2.0)/2.0d0))
    end do
  end do

  do i = 0,IM
    do j = 0,JM
      if ((y(0,JM-2)-dtan(bt)*x(i,j))>y(i,j)) then
        rho(i,j) = rho1 !衝撃波前
        u(i,j) = u1
        v(i,j) = v1
        ene(i,j) = ene1
        p(i,j) = p1
        T(i,j) = T1
      else
        rho(i,j) = rho2 !衝撃波後
        u(i,j) = u2
        v(i,j) = v2
        ene(i,j) = ene2
        p(i,j) = p2
        T(i,j) = T2
      end if
      Q(1,i,j) = rho(i,j)/Jac(i,j) !保存変数ベクトル
      Q(2,i,j) = (rho(i,j)*u(i,j))/Jac(i,j)
      Q(3,i,j) = (rho(i,j)*v(i,j))/Jac(i,j)
      Q(4,i,j) = ene(i,j)/Jac(i,j)
    end do
  end do
end subroutine conditions

!------tvdスキーム（ξ方向）------
subroutine tvd_xi
  double precision :: RWL,RWR,AKX,AKY,AJACM,UM,VM,CM,PHI,BETA,AKXT,AKYT,THIT
  double precision,dimension(1:4,1:4,0:IM) :: RR,RI
  double precision,dimension(1:4,0:IM) :: EIGM,EIG,GG,PHIM,EH
  double precision,dimension(0:IM,0:JM) :: DUM
  double precision,dimension(1:4,0:IM) :: E
  double precision,dimension(1:4) :: D
  double precision,dimension(1:4,-1:IM) :: ALPHA
  double precision :: S,DELTA,GAMMA,UCONT

  do j = 1,JM-1
    do i = 0,IM-1
      RWL = dSQRT(Q(1,i,j)*Jac(i,j))/(dSQRT(Q(1,i+1,j)*Jac(i+1,j))+dSQRT(Q(1,i,j)*Jac(i,j)))
      RWR = dSQRT(Q(1,i+1,j)*Jac(i+1,j))/(dSQRT(Q(1,i+1,j)*Jac(i+1,j))+dSQRT(Q(1,i,j)*Jac(i,j)))
      AKX = 0.5d0*(xix(i,j)+xix(i+1,j)) !中間点のkx
      AKY = 0.5d0*(xiy(i,j)+xiy(i+1,j)) !中間点のky
      AJACM = 0.5d0*(Jac(i,j)+Jac(i+1,j)) !中間点のJac

      UM = RWL*u(i,j)+RWR*u(i+1,j) !中間点のu
      VM = RWL*v(i,j)+RWR*v(i+1,j) !中間点のv
      DUM(i,j) = RWL*Q(4,i,j)/Q(1,i,j)+RWR*Q(4,i+1,j)/Q(1,i+1,j) !e/ρ
      CM = dSQRT(g*(g-1.0d0)*dabs(DUM(i,j)-0.5d0*(UM**2.0+VM**2.0))) !中間点の音速c
      !NOMENCLATURE
      PHI = 0.5d0*(g-1.0d0)*(UM**2.0+VM**2.0) !Φ
      BETA = 1.0d0/(2.0d0*CM**2.0)
      AKXT = AKX/dSQRT(AKX**2.0+AKY**2.0)
      AKYT = AKY/dSQRT(AKX**2.0+AKY**2.0)
      THIT = AKXT*UM+AKYT*VM

      !RIGHT ENGEN-VECTORS(右固有ベクトル）
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
      RR(4,3,i) = (PHI+CM**2.0)/(g-1.0d0)+CM*THIT
      RR(4,4,i) = (PHI+CM**2.0)/(g-1.0d0)-CM*THIT
      !print *,RR(4,1,i),RR(4,2,i),RR(4,3,i),RR(4,4,i)

      !INVERS OR RIGHT EIGEN-VECTORS(右固有ベクトルの逆行列）
      RI(1,1,i) = 1.0d0-PHI/CM**2.0
      RI(1,2,i) = (g-1.0d0)*UM/CM**2.0
      RI(1,3,i) = (g-1.0d0)*VM/CM**2.0
      RI(1,4,i) = -(g-1.0d0)/CM**2.0
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

      !ENGEN-VALUES AT INTERMEDIATE(固有値）
      EIGM(1,i) = AKX*UM+AKY*VM
      EIGM(2,i) = EIGM(1,i)
      EIGM(3,i) = EIGM(1,i)+CM*dSQRT(AKX**2.0+AKY**2.0)
      EIGM(4,i) = EIGM(1,i)-CM*dSQRT(AKX**2.0+AKY**2.0)

      !ALPHAの導出
      do k = 1,4
        D(k) = (Q(k,i+1,j)*Jac(i+1,j)-Q(k,i,j)*Jac(i,j))/AJACM
      end do

      do k = 1,4
        ALPHA(k,i) = RI(k,1,i)*D(1)+RI(k,2,i)*D(2)+RI(k,3,i)*D(3)+RI(k,4,i)*D(4)
      end do
    end do

    do k = 1,4
      ALPHA(k,-1) = ALPHA(k,0)
      ALPHA(k,IM) = ALPHA(k,IM-1)
    end do

    do i = 0,IM
      DUM = dSQRT(g*R*T(i,j)*(xix(i,j)**2.0+xiy(i,j)**2.0))
      EIG(1,i) = xix(i,j)*u(i,j)+xiy(i,j)*v(i,j)
      EIG(2,i) = EIG(1,i)
      EIG(3,i) = EIG(1,i)+DUM(i,j)
      EIG(4,i) = EIG(1,i)-DUM(i,j)
    end do

    !CAL LIMITER FUNCTION G(制御関数)
    do k = 1,4
      do i = 0,IM
        !MINMOD LIMITER
        S = dSIGN(1.0d0,ALPHA(k,i))
        GG(k,i) = S*AMAX1(0.0d0,AMIN1(abs(ALPHA(k,i)),S*ALPHA(k,i-1)))
      end do
    end do

    do k = 1,4
      do i = 0,IM-1
        !GAMMA
        DELTA = 1.0d-10
        !DELTA=AMAX1(0.0,EIGM(k,i)-EIG(k,i),EIG(k,i+1)-EIGM(k,i))
        if (ALPHA(k,i)/=0.0) then
          GAMMA = 0.5d0*FPSI(EIGM(k,i),DELTA)*(GG(k,i+1)-GG(k,i))/ALPHA(k,i)
        else
          GAMMA = 0.0d0
        end if
        !PHI(Φ)
        PHIM(k,i) = 0.5d0*FPSI(EIGM(k,i),DELTA)*(GG(k,i+1)+GG(k,i))-FPSI(EIGM(k,i)+GAMMA,DELTA)*ALPHA(k,i)
      end do
    end do

    !CONVECTION COMPORNENTS(計算空間上のE計算)
    do i = 0,IM
      UCONT = xix(i,j)*u(i,j)+xiy(i,j)*v(i,j) !(反変速度)
      E(1,i) = Q(1,i,j)*UCONT
      E(2,i) = Q(2,i,j)*UCONT+p(i,j)*xix(i,j)/Jac(i,j)
      E(3,i) = Q(3,i,j)*UCONT+p(i,j)*xiy(i,j)/Jac(i,j)
      E(4,i) = Q(4,i,j)*UCONT+p(i,j)*UCONT/Jac(i,j)
    end do

    !XI-DIRECTION CONVECTION FLUX
    do k = 1,4
      do i = 0,IM-1
        EE(k,i,j) = 0.5d0*(E(k,i)+E(k,i+1)+RR(k,1,i)*PHIM(1,i)+RR(k,2,i)*PHIM(2,i)&
                                                &+RR(k,3,i)*PHIM(3,i)+RR(k,4,i)*PHIM(4,i))
      end do
    end do
  end do
  return
end subroutine tvd_xi

!------tvdスキーム（η方向）------
subroutine tvd_eta
  double precision :: RWL,RWR,AKX,AKY,AJACM,UM,VM,CM,PHI,BETA,AKXT,AKYT,THIT
  double precision,dimension(1:4,1:4,0:JM) :: RR,RI
  double precision,dimension(0:IM,0:JM) :: DUM
  double precision,dimension(1:4,0:IM) :: FIGM,FIG,GG,PHIM
  double precision,dimension(1:4,0:JM) :: F
  double precision,dimension(1:4) :: D
  double precision,dimension(1:4,-1:IM) :: ALPHA
  double precision :: S,DELTA,GAMMA,UCONT

  do i = 1,IM-1
    do j = 0,JM-1
      RWL = dSQRT(Q(1,i,j)*Jac(i,j))/(dSQRT(Q(1,i,j+1)*Jac(i,j+1))+dSQRT(Q(1,i,j)*Jac(i,j)))
      RWR = dSQRT(Q(1,i,j+1)*Jac(i,j+1))/(dSQRT(Q(1,i,j+1)*Jac(i,j+1))+dSQRT(Q(1,i,j)*Jac(i,j)))
      AKX = 0.5d0*(ex(i,j)+ex(i,j+1))
      AKY = 0.5d0*(ey(i,j)+ey(i,j+1))
      AJACM = 0.5d0*(Jac(i,j)+Jac(i,j+1))

      UM = RWL*u(i,j)+RWR*u(i,j+1)
      VM = RWL*v(i,j)+RWR*v(i,j+1)
      DUM(i,j) = RWL*Q(4,i,j)/Q(1,i,j)+RWR*Q(4,i,j+1)/Q(1,i,j+1)
      CM = dSQRT(g*(g-1.0d0)*dabs(DUM(i,j)-0.5d0*(UM**2.0+VM**2.0)))

      !NOMENCLATURE
      PHI = 0.5d0*(g-1.0d0)*(UM**2.0+VM**2.0)
      BETA = 1.0d0/(2.0d0*CM**2.0)
      AKXT = AKX/dSQRT(AKX**2.0+AKY**2.0)
      AKYT = AKY/dSQRT(AKX**2.0+AKY**2.0)
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
      RR(4,3,j) = (PHI+CM**2.0)/(g-1.0d0)+CM*THIT
      RR(4,4,j) = (PHI+CM**2.0)/(g-1.0d0)-CM*THIT

      !INVERS OR RIGHT EIGEN-VECTORS
      RI(1,1,j) = 1.0d0-PHI/CM**2.0
      RI(1,2,j) = (g-1.0d0)*UM/CM**2.0
      RI(1,3,j) = (g-1.0d0)*VM/CM**2.0
      RI(1,4,j) = -(g-1.0d0)/CM**2.0
      RI(2,1,j) = -AKYT*UM+AKXT*VM
      RI(2,2,j) = AKYT
      RI(2,3,j) = -AKXT
      RI(2,4,j) = 0.0d0
      RI(3,1,j) = BETA*(PHI-CM*THIT)
      RI(3,2,j) = BETA*(AKXT*CM-(g-1.0d0)*UM)
      RI(3,3,j) = BETA*(AKYT*CM-(g-1.0d0)*VM)
      RI(3,4,j) = BETA*(g-1.0d0)
      RI(4,1,j) = BETA*(PHI+CM*THIT)
      RI(4,2,j) = -BETA*(AKXT*CM+(g-1.0d0)*UM)
      RI(4,3,j) = -BETA*(AKYT*CM+(g-1.0d0)*VM)
      RI(4,4,j) = BETA*(g-1.0d0)

      !ENGEN-VALUES AT INTERMEDIATE
      FIGM(1,j) = AKX*UM+AKY*VM
      FIGM(2,j) = FIGM(1,j)
      FIGM(3,j) = FIGM(1,j)+CM*dSQRT(AKX**2.0+AKY**2.0)
      FIGM(4,j) = FIGM(1,j)-CM*dSQRT(AKX**2.0+AKY**2.0)

      !ALPHA
      do k = 1,4
        D(k) = (Q(k,i,j+1)*Jac(i,j+1)-Q(k,i,j)*Jac(i,j))/AJACM
      end do

      do k = 1,4
        ALPHA(k,j) = RI(k,1,j)*D(1)+RI(k,2,j)*D(2)+RI(k,3,j)*D(3)+RI(k,4,j)*D(4)
      end do
    end do

    do k = 1,4
      ALPHA(k,-1) = ALPHA(k,0)
      ALPHA(k,JM) = ALPHA(k,JM-1)
    end do

    do j = 0,JM
      DUM(i,j) = dSQRT(g*R*T(i,j)*(ex(i,j)**2+ey(i,j)**2))
      FIG(1,j) = ex(i,j)*u(i,j)+ey(i,j)*v(i,j)
      FIG(2,j) = FIG(1,j)
      FIG(3,j) = FIG(1,j)+DUM(i,j)
      FIG(4,j) = FIG(1,j)-DUM(i,j)
    end do

    !CAL LIMITER FUNCTION G
    do k = 1,4
      do j = 0,JM
        !MINMOD LIMITER
        S = dSIGN(1.0d0,ALPHA(k,j))
        GG(k,j) = S*MAX1(0.0d0,AMIN1(abs(ALPHA(k,j)),S*ALPHA(k,j-1)))
      end do
    end do

    do k = 1,4
      do j = 0,JM-1
        !GAMMA
        DELTA = 1.0d-10
        !DELTA=AMAX1(0.0d0,FIGM(k,j)-FIG(k,j),FIG(k,j+1)-FIGM(k,j))
        if (ALPHA(k,j)/=0.0) then
          GAMMA = 0.5d0*FPSI(FIGM(k,i),DELTA)*(GG(k,j+1)-GG(k,j))/ALPHA(k,j)
        else
          GAMMA = 0.0d0
        end if
        !PHI
        PHIM(k,j) = 0.5d0*FPSI(FIGM(k,j),DELTA)*(GG(k,j+1)+GG(k,j))-FPSI(FIGM(k,j)+GAMMA,DELTA)*ALPHA(k,j)
      end do
    end do

    !CONVECTION COMPORNENTS
    do j = 0,JM
      UCONT = ex(i,j)*u(i,j)+ey(i,j)*v(i,j)
      F(1,j) = Q(1,i,j)*UCONT
      F(2,j) = Q(2,i,j)*UCONT+p(i,j)*ex(i,j)/Jac(i,j)
      F(3,j) = Q(3,i,j)*UCONT+p(i,j)*ey(i,j)/Jac(i,j)
      F(4,j) = Q(4,i,j)*UCONT+p(i,j)*UCONT/Jac(i,j)
    end do

    !XI-DIRECTION CONVECTION FLUX
    do k = 1,4
      do j = 0,JM-1
        FF(k,i,j) = 0.5d0*(F(k,j)+F(k,j+1)+RR(k,1,j)*PHIM(1,j)&
                  &+RR(k,2,j)*PHIM(2,j)+RR(k,3,j)*PHIM(3,j)+RR(k,4,j)*PHIM(4,j))
      end do
    end do
  end do
  return
end subroutine tvd_eta

!*****FPSI関数の定義*****
function FPSI(z,DELTA) !Ψ(z)
  double precision :: DELTA,z,FPSI
  if (abs(z)>=DELTA) then
    FPSI = abs(z)
  else
    FPSI = 0.5d0*(z**2.0+DELTA**2.0)/DELTA
  end if
  return
end function FPSI

!------Runge-Kutta------
subroutine Runge_Kutta
  do j = 1,JM-1
    do i = 1,IM-1
      do k = 1,4
        QQ(k,i,j,rk) = QQ(k,i,j,0)-(1.0d0/dble(5-rk))*dt*((EE(k,i,j)-EE(k,i-1,j))/1.0d0+(FF(k,i,j)-FF(k,i,j-1))/1.0d0)
        Q(k,i,j) = QQ(k,i,j,rk)
      end do
    end do
  end do
end subroutine Runge_Kutta

!------諸量の計算------
subroutine calquation
  do j = 1,JM-1
    do i = 1,IM-1
      un = u(i,j)
      vn = v(i,j)
      rho(i,j) = Q(1,i,j)*Jac(i,j)
      u(i,j) = Q(2,i,j)*Jac(i,j)/rho(i,j)
      v(i,j) = Q(3,i,j)*Jac(i,j)/rho(i,j)
      ene(i,j) = Q(4,i,j)*Jac(i,j)
      enebar(i,j) = ene(i,j)/rho(i,j)-(u(i,j)**2.0+v(i,j)**2.0)/2.0d0
      T(i,j) = enebar(i,j)*(g-1.0d0)/R
      p(i,j) = rho(i,j)*enebar(i,j)*(g-1.0d0)
      ddu = dabs((u(i,j)-un(i,j))/un(i,j))
      if (ddu>du) then
        du = ddu
      end if
      ddv = dabs((v(i,j)-vn(i,j))/vn(i,j))
      if (ddv>dv) then
        dv = ddv
      end if
    end do
  end do
end subroutine calquation

!------境界条件------
subroutine boundary

  do i = 0,IM !上下壁
    slope = (y(i,0)-y(i,2))/(y(i,2)-y(i,1))
    rho(i,0) = rho(i,2)+slope*(rho(i,2)-rho(i,1))
    u(i,0) = u(i,2)+slope*(u(i,2)-u(i,1))
    v(i,0) = 0.0d0
    p(i,0) = p(i,2)+slope*(p(i,2)-p(i,1))
    ene(i,0) = p(i,0)/(g-1.0d0)+rho(i,0)*((u(i,0)**2.0+v(i,0)**2.0))/2.0d0
    T(i,0) = T(i,2)+slope*(T(i,2)-T(i,1))

    slope = (y(i,JM)-y(i,JM-2))/(y(i,JM-2)-y(i,JM-1))
    rho(i,JM) = rho(i,JM-2)+slope*(rho(i,JM-2)-rho(i,JM-1))
    u(i,JM) = u(i,JM-2)+slope*(u(i,JM-2)-u(i,JM-1))
    v(i,JM) = v(i,JM-2)+slope*(v(i,JM-2)-v(i,JM-1))
    p(i,JM) = p(i,JM-2)+slope*(p(i,JM-2)-p(i,JM-1))
    ene(i,JM) = p(i,jm)/(g-1.0d0)+rho(i,jm)*((u(i,jm)**2.0)+(v(i,jm)**2.0))/2.0d0
    T(i,JM) = T(i,JM-2)+slope*(T(i,JM-2)-T(i,JM-1))
  end do

  do j = 1,JM-1 !流出
    slope = (x(IM,j)-x(IM-2,j))/(x(IM-2,j)-x(IM-1,j))
    rho(IM,j) = rho(IM-2,j)+slope*(rho(IM-2,j)-rho(IM-1,j))
    u(IM,j) = u(IM-2,j)+slope*(u(IM-2,j)-u(IM-1,j))
    v(IM,j) = v(IM-2,j)+slope*(v(IM-2,j)-v(IM-1,j))
    p(IM,j) = p(IM-2,j)+slope*(p(IM-2,j)-p(IM-1,j))
    ene(IM,j) = p(im,j)/(g-1.0d0)+rho(im,j)*((u(im,j)**2.0)+(v(im,j)**2.0))/2.0d0
    T(IM,j) = T(IM-2,j)+slope*(T(IM-2,j)-T(IM-1,j))
  end do

  do j = 0,JM
    do i = 0,IM
      if (i==IM .or. j==0 .or. j==JM) then
        Q(1,i,j) = rho(i,j)/Jac(i,j)
        Q(2,i,j) = rho(i,j)*u(i,j)/Jac(i,j)
        Q(3,i,j) = rho(i,j)*v(i,j)/Jac(i,j)
        Q(4,i,j) = ene(i,j)/Jac(i,j)
      end if
    end do
  end do
end subroutine boundary
!------------
end program assyuku5_1
