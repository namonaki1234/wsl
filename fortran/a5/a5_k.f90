program a5
implicit none

double precision,parameter :: H = 10.0d0
double precision,parameter :: r_y = 1.1d0
integer,parameter :: IM = 30,JM = 30
double precision,parameter :: dx = 3.0d0*H/dble(IM)
double precision :: a,time
double precision,parameter :: r = 287.1d0 !ガス定数
double precision,parameter :: gam = 1.4d0 !比熱比
double precision,parameter :: dt = 1.0d-6
double precision,parameter :: pi = dacos(-1.0d0)
double precision,parameter :: betas = 3.0d1*pi/1.80d2
double precision,parameter :: EPS = 1.0d-4
double precision,dimension(0:IM,0:JM) :: x,y
double precision,dimension(0:IM,0:JM) :: energy = 0.0d0
double precision,dimension(0:IM,0:JM) :: rho1,rho2,u,v,p,T,M,V_abs,V_n,V_t,theta,us,vs,E_bar
double precision,dimension(0:IM,0:JM) :: xi_x,xi_y,eta_x,eta_y,x_xi,y_xi,x_eta,y_eta,Jac,Jac_inv,U_cvc,V_cvc
double precision,dimension(1:4,0:IM,0:JM) :: Q,EE,FF
double precision,dimension(1:4,0:IM,0:JM,0:5) :: QQ
double precision,dimension(4,4,0:IM) :: RR,RI
double precision,dimension(4,0:IM) :: EIGM,FIGM,GG,PHIM,E,F,FIG
! double precision,dimension(4,0:IM) :: EIG
double precision,dimension(4,-1:IM) :: ALPHA
double precision,dimension(4) :: D
integer :: i,j,n,k_runge_kutta,k
double precision :: slope,dd,duv, &
  AKX,AKY,AJACM,UM,VM,DUM,CM,BETA,AKXT,THIT,RWL,RWR,S,GAMMA,UCONT,AKYT,DELTA,PHI

print*,"*****START CALCULATION*****"
time = 0.0d0
call generate_grid
call init_condition
call coordinate_transformation

do n = 0,1000000 !計算回数
  do j = 0,JM
    do i = 0,IM
      do k = 1,4
        QQ(k,i,j,0) = Q(k,i,j)
      end do
    end do
  end do

  do k_runge_kutta = 1,4
    call xi_TVD
    call eta_TVD
    call runge_kutta
    call update_Q
  end do
  if (mod(n,1000)==0) then
    print*,"n=",n
  end if
  time = time+dt
  if (dd<EPS) exit
end do
print*,"*****END CALCULATION*****"

open (10,file='a5.dat',status='replace')
do j = 0,JM
  do i = 0,IM
    write (10,'(f13.6,$)') x(i,j)
    write (10,'(f13.6,$)') y(i,j)
    write (10,'(f13.6,$)') u(i,j)
    write (10,'(f13.6,$)') v(i,j)
    write (10,'(f13.3,$)') p(i,j)
    write (10,'(f13.6)') t(i,j)
  end do
end do
close (10)

open (11,file='a5.fld',status='replace')
write (11,'(a)') '# AVS field file'
write (11,'(a)') 'ndim=2'
write (11,'(a,i0)') 'dim1=',IM+1
write (11,'(a,i0)') 'dim2=',JM+1
write (11,'(a)') 'nspace=2'
write (11,'(a)') 'veclen=4'
write (11,'(a)') 'data=float'
write (11,'(a)') 'field=irregular'
write (11,'(a)') 'label=u,v,p,t'
write (11,'(a)') 'variable 1 file=a5.dat filetype=ascii skip=0 offset=2 stride=6'
write (11,'(a)') 'variable 2 file=a5.dat filetype=ascii skip=0 offset=3 stride=6'
write (11,'(a)') 'variable 3 file=a5.dat filetype=ascii skip=0 offset=4 stride=6'
write (11,'(a)') 'variable 4 file=a5.dat filetype=ascii skip=0 offset=5 stride=6'
write (11,'(a)') 'coord 1 file=a5.dat filetype=ascii skip=0 offset=0 stride=6'
write (11,'(a)') 'coord 2 file=a5.dat filetype=ascii skip=0 offset=1 stride=6'
close (11)

contains

subroutine generate_grid !不等間隔格子生成

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

end subroutine generate_grid

subroutine init_condition !初期設定
  do j = 0,JM
    do i = 0,IM
      M(i,j) = 2.9d0
      T(i,j) = 293.0d0
      rho1(i,j) = 1.2d0
      p(i,j) = rho1(i,j)*R*T(i,j)
      u(i,j) = M(i,j)*sqrt(gam*R*T(i,j))
      v(i,j) = 0.d0
      V_abs(i,j) = sqrt((u(i,j)**2)+(v(i,j)**2))
      V_n(i,j) = V_abs(i,j)*sin(beta)
      V_t(i,j) = V_abs(i,j)*cos(beta)
      energy(i, j) = rho1(i, j)*((r*t(i, j)/(gam - 1.0d0)) + ((u(i, j)**2.0) + v(i, j)**2.0)/2.0d0)
    end do
  end do

  !********************衝撃波条件（衝撃波通過後の変化）********************
  do j = 0,JM
    do i = 0,IM
      if (y(i,j)<=y(i,JM-2) .and. x(i,j)<=sqrt(3.d0)*(y(i,JM-2)-y(i,j))) cycle

      T(i,j) = T(i,j)*(2.d0*gam*(M(i,j)**2)*(sin(beta)**2)-(gam-1.d0)) &
               *((gam-1.0d0)*(M(i,j)**2)*(sin(beta)**2)+2.d0) &
               /(((gam+1.0d0)**2)*(M(i,j)**2)*(sin(beta)**2))

      p(i,j) = p(i,j)*(2.d0*gam*(M(i,j)**2)*(sin(beta)**2)-(gam-1.d0)) &
               /(gam+1.0d0)

      rho2(i,j) = rho1(i,j)*(gam+1.d0)*(M(i,j)**2)*(sin(beta)**2) &
                  /((gam-1.0d0)*(M(i,j)**2)*(sin(beta)**2)+2.d0)

      V_n(i,j) = rho1(i,j)*V_n(i,j)/rho2(i,j)

      V_abs(i,j) = sqrt((V_n(i,j)**2)+(V_t(i,j)**2))

      theta(i,j) = beta-atan(V_n(i,j)/V_t(i,j))
      u(i,j) = V_abs(i,j)*cos(-theta(i,j))
      v(i,j) = V_abs(i,j)*sin(-theta(i,j))

      U_cvc(i,j) = xi_x(i,j)*u(i,j)+xi_y(i,j)*v(i,j)
      V_cvc(i,j) = eta_x(i,j)*u(i,j)+eta_y(i,j)*v(i,j)
      energy(i, j) = rho1(i, j)*((r*t(i, j)/(gam - 1.0d0)) + ((u(i, j)**2.0) + v(i, j)**2.0)/2.0d0)
    end do
  end do
end subroutine init_condition

subroutine coordinate_transformation !座標変換
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

  do j = 0,JM
    do i = 0,IM
      Q(1,i,j) = rho1(i,j)/Jac(i,j)
      Q(2,i,j) = (rho1(i,j)*u(i,j))/Jac(i,j)
      Q(3,i,j) = (rho1(i,j)*v(i,j))/Jac(i,j)
      Q(4,i,j) = energy(i,j)/Jac(i,j)
    end do
  end do

end subroutine coordinate_transformation

subroutine xi_TVD
!*****ξ方向*****
  do J = 1,JM-1
    do I = 0,IM-1

      RWL = dSQRT(Q(1,I,J)*Jac(I,J)) &
            /(dSQRT(Q(1,I+1,J)*Jac(I+1,J))+dSQRT(Q(1,I,J)*Jac(I,J)))
      RWR = dSQRT(Q(1,I+1,J)*Jac(I+1,J)) &
            /(dSQRT(Q(1,I+1,J)*Jac(I+1,J))+dSQRT(Q(1,I,J)*Jac(I,J)))
      AKX = 0.5d0*(xi_x(I,J)+xi_x(I+1,J)) !k_x
      AKY = 0.5d0*(xi_y(I,J)+xi_y(I+1,J)) !k_y
      AJACM = 0.5d0*(Jac(I,J)+Jac(I+1,J))

      UM = RWL*U(I,J)+RWR*U(I+1,J)
      VM = RWL*V(I,J)+RWR*V(I+1,J)
      DUM = RWL*Q(4,I,J)/Q(1,I,J)+RWR*Q(4,I+1,J)/Q(1,I+1,J)
      CM = dSQRT(gam*(gam-1.0d0)*abs(DUM-0.5d0*(UM**2.0+VM**2.0)))

!       NOMENCLATURE
      PHI = 0.5d0*(gam-1.0d0)*(UM**2.0+VM**2.0)
      BETA = 1.0d0/(2.0d0*CM**2.0)
      AKXT = AKX/dSQRT(AKX**2.0+AKY**2.0)
      AKYT = AKY/dSQRT(AKX**2.0+AKY**2.0)
      THIT = AKXT*UM+AKYT*VM
!
!       RIGHT ENGEN-VECTORS     P.88 R_k
      RR(1,1,I) = 1.0d0
      RR(1,2,I) = 0.0d0
      RR(1,3,I) = 1.0d0
      RR(1,4,I) = 1.0d0
      RR(2,1,I) = UM
      RR(2,2,I) = AKYT
      RR(2,3,I) = UM+AKXT*CM
      RR(2,4,I) = UM-AKXT*CM
      RR(3,1,I) = VM
      RR(3,2,I) = -AKXT
      RR(3,3,I) = VM+AKYT*CM
      RR(3,4,I) = VM-AKYT*CM
      RR(4,1,I) = PHI/(gam-1.0d0)
      RR(4,2,I) = AKYT*UM-AKXT*VM
      RR(4,3,I) = (PHI+CM**2)/(gam-1.0d0)+CM*THIT
      RR(4,4,I) = (PHI+CM**2)/(gam-1.0d0)-CM*THIT
!
!       INVERS OR RIGHT EIGEN-VECTORS    P.88 R_k^-1
      RI(1,1,I) = 1.0d0-PHI/CM**2.0
      RI(1,2,I) = (gam-1.0d0)*UM/CM**2.0
      RI(1,3,I) = (gam-1.0d0)*VM/CM**2.0
      RI(1,4,I) = -(gam-1.0d0)/CM**2.0
      RI(2,1,I) = -AKYT*UM+AKXT*VM
      RI(2,2,I) = AKYT
      RI(2,3,I) = -AKXT
      RI(2,4,I) = 0.0d0
      RI(3,1,I) = BETA*(PHI-CM*THIT)
      RI(3,2,I) = BETA*(AKXT*CM-(gam-1.0d0)*UM)
      RI(3,3,I) = BETA*(AKYT*CM-(gam-1.0d0)*VM)
      RI(3,4,I) = BETA*(gam-1.0d0)
      RI(4,1,I) = BETA*(PHI+CM*THIT)
      RI(4,2,I) = -BETA*(AKXT*CM+(gam-1.0d0)*UM)
      RI(4,3,I) = -BETA*(AKYT*CM+(gam-1.0d0)*VM)
      RI(4,4,I) = BETA*(gam-1.0d0)
!
!       ENGEN-VALUES AT INTERMEDIATE     P.87  固有値
      EIGM(1,I) = AKX*UM+AKY*VM
      EIGM(2,I) = EIGM(1,I)
      EIGM(3,I) = EIGM(1,I)+CM*dSQRT(AKX**2.0+AKY**2.0)
      EIGM(4,I) = EIGM(1,I)-CM*dSQRT(AKX**2.0+AKY**2.0)
!
!       ALPHA
      do K = 1,4
        D(K) = (Q(K,I+1,J)*Jac(I+1,J)-Q(K,I,J)*Jac(I,J))/AJACM
      end do
!
      do K = 1,4
        ALPHA(K,I) = RI(K,1,I)*D(1)+RI(K,2,I)*D(2)+RI(K,3,I)*D(3)+RI(K,4,I)*D(4)
      end do
    end do
!
    do K = 1,4
      ALPHA(K,-1) = ALPHA(K,0)
      ALPHA(K,IM) = ALPHA(K,IM-1)
    end do
!
    do I = 0,IM
      DUM = dSQRT(gam*R*T(I,J)*(xi_x(I,J)**2+xi_y(I,J)**2))
      ! EIG(1,I) = xi_x(I,J)*U(I,EIG(1,I) = xi_x(I,J)*U(I,J)+xi_y(I,J)*J)
      ! EIG(2,I) = EIG(1,I)
      ! ! EIG(
      ! EIG(3,I) = EIG(1,I)+DUM
      ! ! EIG(EIG(1,I)
      ! EIG(4,I) = EIG(1,I)-DUM
      ! EIG(EIG(1,I)
    end do
!
!       CAL LIMITER FUNCTION gam
    do K = 1,4
      do I = 0,IM
!
!       MINMOD LIMITER         制限関数
        S = dSIGN(1.0d0,ALPHA(K,I))
        GG(K,I) = S*dMAX1(0.0,dMIN1(dabs(ALPHA(K,I)),S*ALPHA(K,I-1)))
      end do
    end do

    do K = 1,4
      do I = 0,IM-1

!       GAMMA
        DELTA = 1.0d-10
!        DELTA=dMAX1(0.0,EIGM(K,I)-EIG(K,I),EIG(K,I+1)-EIGM(K,I))
        if (ALPHA(K,I)/=0.0) then
          GAMMA = 0.5d0*FPSI(EIGM(K,I),DELTA)*(GG(K,I+1)-GG(K,I))/ALPHA(K,I)
        else
          GAMMA = 0.0d0
        end if

!       PHI
        PHIM(K,I) = 0.5d0*FPSI(EIGM(K,I),DELTA)*(GG(K,I+1)+GG(K,I)) &
                    -FPSI(EIGM(K,I)+GAMMA,DELTA)*ALPHA(K,I)
      end do
    end do

!       CONVECTION COMPORNENTS    P.83 E^
    do I = 0,IM
      UCONT = xi_x(I,J)*U(I,J)+xi_y(I,J)*V(I,J)
      E(1,I) = Q(1,I,J)*UCONT
      E(2,I) = Q(2,I,J)*UCONT+P(I,J)*xi_x(I,J)/Jac(i,j)
      E(3,I) = Q(3,I,J)*UCONT+P(I,J)*xi_y(I,J)/Jac(i,j)
      E(4,I) = Q(4,I,J)*UCONT+P(I,J)*UCONT/Jac(i,j)
    end do

!       XI-DIRECTION CONVECTION FLUX　　　　　E~_j+1/2,k
    do K = 1,4
      do I = 0,IM-1
        EE(K,I,J) = 0.5d0*(E(K,I)+E(K,I+1) &
                           +RR(K,1,I)*PHIM(1,I)+RR(K,2,I)*PHIM(2,I) &
                           +RR(K,3,I)*PHIM(3,I)+RR(K,4,I)*PHIM(4,I))
      end do
    end do
  end do
  return
end subroutine xi_TVD

!*****FPSI関数の定義*****           Ψ(z)の計算
double precision function FPSI(Z,DELTA)
  double precision Z,DELTA
  if (dabs(Z)>=DELTA) then
    FPSI = dabs(Z)
  else
    FPSI = 0.5d0*(Z**2.0+DELTA**2.0)/DELTA
  end if
  return
end function

subroutine eta_TVD
!*****η方向*****
  do I = 1,IM-1
    do J = 0,JM-1
      RWL = dSQRT(Q(1,I,J)*Jac(I,J)) &
            /(dSQRT(Q(1,I,J+1)*Jac(I,J+1))+dSQRT(Q(1,I,J)*Jac(I,J)))
      RWR = dSQRT(Q(1,I,J+1)*Jac(I,J+1)) &
            /(dSQRT(Q(1,I,J+1)*Jac(I,J+1))+dSQRT(Q(1,I,J)*Jac(I,J)))
      AKX = 0.5d0*(eta_x(I,J)+eta_x(I,J+1)) !k_x
      AKY = 0.5d0*(eta_y(I,J)+eta_y(I,J+1)) !k_y
      AJACM = 0.5d0*(Jac(I,J)+Jac(I,J+1))
      UM = RWL*U(I,J)+RWR*U(I,J+1)
      VM = RWL*V(I,J)+RWR*V(I,J+1)
      DUM = RWL*Q(4,I,J)/Q(1,I,J)+RWR*Q(4,I,J+1)/Q(1,I,J+1)
      CM = dSQRT(gam*(gam-1.0d0)*dabs(DUM-0.5d0*(UM**2.0+VM**2.0)))

!       NOMENCLATURE
      PHI = 0.5d0*(gam-1.0d0)*(UM**2.0+VM**2.0)
      BETA = 1.0d0/(2.0d0*CM**2.0)
      AKXT = AKX/dSQRT(AKX**2.0+AKY**2.0)
      AKYT = AKY/dSQRT(AKX**2.0+AKY**2.0)
      THIT = AKXT*UM+AKYT*VM
!
!       RIGHT ENGEN-VECTORS     P.88 R_k
      RR(1,1,J) = 1.0d0
      RR(1,2,J) = 0.0d0
      RR(1,3,J) = 1.0d0
      RR(1,4,J) = 1.0d0
      RR(2,1,J) = UM
      RR(2,2,J) = AKYT
      RR(2,3,J) = UM+AKXT*CM
      RR(2,4,J) = UM-AKXT*CM
      RR(3,1,J) = VM
      RR(3,2,J) = -AKXT
      RR(3,3,J) = VM+AKYT*CM
      RR(3,4,J) = VM-AKYT*CM
      RR(4,1,J) = PHI/(gam-1.0d0)
      RR(4,2,J) = AKYT*UM-AKXT*VM
      RR(4,3,J) = (PHI+CM**2.0)/(gam-1.0d0)+CM*THIT
      RR(4,4,J) = (PHI+CM**2.0)/(gam-1.0d0)-CM*THIT
!
!       INVERS OR RIGHT EIGEN-VECTORS    P.88 R_k^-1
      RI(1,1,J) = 1.0d0-PHI/CM**2.0
      RI(1,2,J) = (gam-1.0d0)*UM/CM**2.0
      RI(1,3,J) = (gam-1.0d0)*VM/CM**2.0
      RI(1,4,J) = -(gam-1.0d0)/CM**2.0
      RI(2,1,J) = -AKYT*UM+AKXT*VM
      RI(2,2,J) = AKYT
      RI(2,3,J) = -AKXT
      RI(2,4,J) = 0.0d0
      RI(3,1,J) = BETA*(PHI-CM*THIT)
      RI(3,2,J) = BETA*(AKXT*CM-(gam-1.0d0)*UM)
      RI(3,3,J) = BETA*(AKYT*CM-(gam-1.0d0)*VM)
      RI(3,4,J) = BETA*(gam-1.0d0)
      RI(4,1,J) = BETA*(PHI+CM*THIT)
      RI(4,2,J) = -BETA*(AKXT*CM+(gam-1.0d0)*UM)
      RI(4,3,J) = -BETA*(AKYT*CM+(gam-1.0d0)*VM)
      RI(4,4,J) = BETA*(gam-1.0d0)
!
!       ENGEN-VALUES AT INTERMEDIATE   P.87  固有値
      FIGM(1,J) = AKX*UM+AKY*VM
      FIGM(2,J) = FIGM(1,J)
      FIGM(3,J) = FIGM(1,J)+CM*dSQRT(AKX**2.0+AKY**2.0)
      FIGM(4,J) = FIGM(1,J)-CM*dSQRT(AKX**2.0+AKY**2.0)
!
!       ALPHA
      do K = 1,4
        D(K) = (Q(K,I,J+1)*Jac(I,J+1)-Q(K,I,J)*Jac(I,J))/AJACM
      end do
!
      do K = 1,4
        ALPHA(K,J) = RI(K,1,J)*D(1)+RI(K,2,J)*D(2)+RI(K,3,J)*D(3)+RI(K,4,J)*D(4)
      end do
    end do
!
    do K = 1,4
      ALPHA(K,-1) = ALPHA(K,0)
      ALPHA(K,JM) = ALPHA(K,JM-1)
    end do
!
    do J = 0,JM
      DUM = dSQRT(gam*R*T(I,J)*(eta_x(I,J)**2.0+eta_y(I,J)**2.0))
      FIG(1,J) = eta_x(I,J)*U(I,J)+eta_y(I,J)*V(I,J)
      FIG(2,J) = FIG(1,J)
      FIG(3,J) = FIG(1,J)+DUM
      FIG(4,J) = FIG(1,J)-DUM
    end do
!
!       CAL LIMITER FUNCTION gam
    do K = 1,4
      do J = 0,JM
!
!       MINMOD LIMITER  制限関数
        S = dSIGN(1.0d0,ALPHA(K,J))
        GG(K,J) = S*dMAX1(0.0d0,dMIN1(dabs(ALPHA(K,J)),S*ALPHA(K,J-1)))
      end do
    end do

    do K = 1,4
      do J = 0,JM-1

!       GAMMA
        DELTA = 1.0d-10
        !   DELTA = dMAX1(0.0,EIGM(K,I)-EIG(K,I),EIG(K,I+1)-EIGM(K,I))
        if (ALPHA(K,J)/=0.0) then
          GAMMA = 0.5d0*FPSI(FIGM(K,J),DELTA)*(GG(K,J+1)-GG(K,J))/ALPHA(K,J)
        else
          GAMMA = 0.0d0
        end if

!       PHI
        PHIM(K,J) = 0.5d0*FPSI(FIGM(K,J),DELTA)*(GG(K,J+1)+GG(K,J)) &
                    -FPSI(FIGM(K,J)+GAMMA,DELTA)*ALPHA(K,J)
      end do
    end do

!       CONVECTION COMPORNENTS    P.83 F^
    do J = 0,JM
      UCONT = eta_x(I,J)*U(I,J)+eta_y(I,J)*V(I,J)
      F(1,J) = Q(1,I,J)*UCONT
      F(2,J) = Q(2,I,J)*UCONT+P(I,J)*eta_x(I,J)/Jac(i,j)
      F(3,J) = Q(3,I,J)*UCONT+P(I,J)*eta_y(I,J)/Jac(i,j)
      F(4,J) = Q(4,I,J)*UCONT+P(I,J)*UCONT/Jac(i,j)
    end do

!       XI-DIRECTION CONVECTION FLUX　　　　　F~_j+1/2,k
    do K = 1,4
      do J = 0,JM-1
        FF(K,I,J) = 0.5d0*(F(K,J)+F(K,J+1) &
                           +RR(K,1,J)*PHIM(1,J)+RR(K,2,J)*PHIM(2,J) &
                           +RR(K,3,J)*PHIM(3,J)+RR(K,4,J)*PHIM(4,J))
      end do
    end do
  end do

  return
end subroutine eta_TVD

subroutine runge_kutta
  do j = 1,JM-1
    do i = 1,IM-1
      do k = 1,4
        QQ(k,i,j,k_runge_kutta) = QQ(k,i,j,0)-(1.0d0/dble(5-k_runge_kutta))*dt*(-EE(k,i-1,j)+EE(k,i,j)-FF(k,i,j-1)+FF(k,i,j))
        Q(k,i,j) = QQ(k,i,j,k_runge_kutta)
      end do
    end do
  end do

end subroutine runge_kutta

subroutine update_Q !ベクトルQの更新
  us = u
  vs = v

  do j = 1,JM-1
    do i = 1,IM-1
      rho1(i,j) = Jac(i,j)*Q(1,i,j)
      u(i,j) = Jac(i,j)*Q(2,i,j)/rho1(i,j)
      v(i,j) = Jac(i,j)*Q(3,i,j)/rho1(i,j)
      energy(i,j) = Jac(i,j)*Q(4,i,j)
      E_bar(i, j) = energy(i, j)/rho1(i, j) - (((u(i, j)**2.0) + v(i, j)**2.0)/2.0d0)
      t(i, j) = E_bar(i, j)*(gam - 1.0d0)/r
      p(i, j) = (gam - 1.0d0)*rho1(i, j)*E_bar(i, j)
    end do
  end do

  do i = 0,IM
    slope = (y(i,JM)-y(i,JM-2))/(y(i,JM-2)-y(i,JM-1))
    rho1(i,JM) = rho1(i,JM-2)+slope*(rho1(i,JM-2)-rho1(i,JM-1))
    u(i,JM) = u(i,JM-2)+slope*(u(i,JM-2)-u(i,JM-1))
    v(i,JM) = v(i,JM-2)+slope*(v(i,JM-2)-v(i,JM-1))
    p(i,JM) = p(i,JM-2)+slope*(p(i,JM-2)-p(i,JM-1))
    t(i,JM) = t(i,JM-2)+slope*(t(i,JM-2)-t(i,JM-1))
    energy(i,JM) = p(i,JM)/(gam-1.0d0)+rho1(i,JM)*((u(i,JM)**2.0)+(v(i,JM)**2.0))/2.0d0

    Q(1,i,JM) = rho1(i,JM)/Jac(i,JM)
    Q(2,i,JM) = rho1(i,JM)*u(i,JM)/Jac(i,JM)
    Q(3,i,JM) = rho1(i,JM)*v(i,JM)/Jac(i,JM)
    Q(4,i,JM) = energy(i,JM)/Jac(i,JM)

    slope = (y(i,0)-y(i,2))/(y(i,2)-y(i,1))
    rho1(i,0) = rho1(i,2)+slope*(rho1(i,2)-rho1(i,1))
    p(i,0) = p(i,2)+slope*(p(i,2)-p(i,1))
    u(i,0) = u(i,2)+slope*(u(i,2)-u(i,1))
    v(i,0) = 0.0d0
    t(i,0) = t(i,2)+slope*(t(i,2)-t(i,1))
    energy(i,0) = p(i,0)/(gam-1.0d0)+rho1(i,0)*((u(i,0)**2.0+v(i,0)**2.0))/2.0d0

    Q(1,i,0) = rho1(i,0)/Jac(i,0)
    Q(2,i,0) = rho1(i,0)*u(i,0)/Jac(i,0)
    Q(3,i,0) = rho1(i,0)*v(i,0)/Jac(i,0)
    Q(4,i,0) = energy(i,0)/Jac(i,0)
  end do

  do j = 1,JM-1
    slope = (x(IM,j)-x(IM-2,j))/(x(IM-2,j)-x(IM-1,j))
    rho1(IM,j) = rho1(IM-2,j)+slope*(rho1(IM-2,j)-rho1(IM-1,j))
    p(IM,j) = p(IM-2,j)+slope*(p(IM-2,j)-p(IM-1,j))
    u(IM,j) = u(IM-2,j)+slope*(u(IM-2,j)-u(IM-1,j))
    v(IM,j) = v(IM-2,j)+slope*(v(IM-2,j)-v(IM-1,j))
    t(IM,j) = t(IM-2,j)+slope*(t(IM-2,j)-t(IM-1,j))
    energy(IM,j) = p(IM,j)/(gam-1.0d0)+rho1(IM,j)*((u(IM,j)**2.0)+(v(IM,j)**2.0))/2.0d0

    Q(1,IM,j) = rho1(IM,j)/Jac(IM,j)
    Q(2,IM,j) = rho1(IM,j)*u(IM,j)/Jac(IM,j)
    Q(3,IM,j) = rho1(IM,j)*v(IM,j)/Jac(IM,j)
    Q(4,IM,j) = energy(IM,j)/Jac(IM,j)
  end do

  dd = 0.0d0
  do j = 1,JM
    do i = 0,IM
      duv = (abs(u(i,j)-us(i,j)))/abs(us(i,j))
      if (dd<duv) then
        dd = duv
      end if
      duv = (abs(v(i,j)-vs(i,j)))/abs(vs(i,j))
      if (dd<duv) then
        dd = duv
      end if
    end do
  end do

end subroutine update_Q

end program a5
