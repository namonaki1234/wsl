program assyuku7
implicit none

!------��������------
double precision,parameter :: H = 0.01d0
integer,parameter :: IM = 30,JM = 30
double precision,parameter :: dx = 3.0d0*H/30.0d0
double precision,parameter :: b = 1.2d0,EPS = 1.0d-4
double precision,parameter :: bt = 30.0d0/180.0d0*dacos(-1.0d0)
double precision :: a,ix,slope
double precision,parameter :: R = 287.1d0,g = 1.4d0,dt = 1.0d-8,T0 = 293.0d0,myu0 = 1.0d-3,SS = 110.0d0,Pr = 0.72d0
integer :: i,j,n,k,rk
double precision,dimension(0:IM,0:JM) :: x,y
double precision,dimension(0:IM,0:JM) :: rho = 1.2d0,p = 0.0d0,u = 0.0d0,v = 0.0d0,ene = 0.0d0,enebar,T,c
double precision :: M1,rho1 = 1.2d0,rho2,T1,T2,p1,p2,u1,u2,v1,v2,vv1,vv2,vt,vn1,vn2,ene1,ene2,theta,un,vn,du,ddu,dv,ddv,myu
double precision :: AJACM,AKX,AKXT,AKY,AKYT,CM,DELTA,DUM,GAMMA,PHI,RWL,RWR,THIT,UCONT,UM,VM,S,BETA,R4,S4,txx,tyy,txy,tyx
double precision,dimension(0:IM,0:JM) :: xh,yh,xe,ye
double precision,dimension(0:IM,0:JM) :: xix,xiy,ex,ey,Jac
double precision,dimension(0:IM,0:JM) :: uu,vv,ue,uh,ve,vh,Te,Th
double precision,dimension(1:4,0:IM,0:JM) :: Q,EE,FF,SSS,RRR
double precision,dimension(1:4,0:IM,0:JM,0:5) :: QQ
double precision,dimension(4,4,0:IM) :: RR,RI
double precision,dimension(4,0:IM) :: EIGM,FIGM,GG,PHIM,E,F,EIG,FIG
double precision,dimension(4) :: D
double precision,dimension(4,-1:IM) :: ALPHA

!------���C���v�Z------
call mesh
call conditions
call coordinate
n = 0
du = EPS+1.0d0
dv = EPS+1.0d0
do while (((du>EPS) .and. (dv>EPS)) .and. (n<100000))
  print*,n
  n = n+1
  du = 0
  dv = 0
  do j = 0,JM
    do i = 0,IM
      do k = 1,4
        QQ(k,i,j,0) = Q(k,i,j)
      end do
    end do
  end do
  do rk = 1,4 !Runge-Kutta�̒i�K��
    call tvd_xi
    call tvd_eta
    call Runge_Kutta
    call viscosity_xi
    call viscosity_eta
    call calculation
    call boundary
  end do
end do

!�O���t�쐬
open (1,file='assyuku7.dat',status='replace')
do j = 0,JM
  do i = 0,IM
    write (1,*) x(i,j),y(i,j),u(i,j),v(i,j),T(i,j),p(i,j)
  end do
end do
close (1)

!------------------------------------------------------------
contains

!------�i�q����------
subroutine mesh
  do j = 0,JM
    do i = 0,IM
      x(i,j) = dble(i)*dx !���Ԋu�i�q
      a = (b+1.0d0)/(b-1.0d0)
      ix = dble(j)/JM
      y(i,j) = H*(a**ix-1.0d0)/(a-1.0d0) !�s���Ԋu�i�q
    end do
  end do
end subroutine mesh

!------������------
subroutine conditions
  do j = 0,JM
    do i = 0,IM
      M1 = 2.9d0
      T1 = 293.0d0
      p1 = rho1*R*T1
      u1 = M1*dsqrt(g*R*T1)
      v1 = 0.0d0
      vv1 = u1
      vt = vv1*dcos(bt)
      vn1 = u1*dsin(bt)
      ene1 = rho1*((R*T1/(g-1.0d0))+((u1**2.0d0)+(v1**2.0d0))/2.0d0)

      vn2 = vn1*((g-1.0d0)*(M1*dsin(bt))**2+2.0d0)/((g+1.0d0)*((M1*dsin(bt))**2.0))
      vv2 = dsqrt(vn2**2.0+vt**2.0)
      theta = bt-datan(vn2/vt)
      u2 = vv2*dcos(theta)
      v2 = vv2*dsin(-theta)
      T2 = T1*((2.0d0*g*(M1**2)*(dsin(bt))**2)-(g-1.0d0))*(((g-1.0d0)*(M1**2)*(dsin(bt))**2.0)+2.0d0)/&
             &(((g+1.0d0)**2)*(M1**2)*((dsin(bt))**2))
      rho2 = rho1*vn1/vn2
      p2 = p1*(2.0d0*g*((M1*dsin(bt))**2.0)-(g-1.0d0))/(g+1.0d0)
      ene2 = rho2*((R*T2/(g-1.0d0))+(u2**2.0+v2**2.0)/2.0d0)
    end do
  end do

  do j = 0,JM
    do i = 0,IM
      if ((y(0,JM-2)-dtan(bt)*x(i,j))>y(i,j)) then
        rho(i,j) = rho1 !�Ռ��g�O
        u(i,j) = u1
        v(i,j) = v1
        ene(i,j) = ene1
        T(i,j) = T1
        p(i,j) = p1
      else
        rho(i,j) = rho2 !�Ռ��g��
        u(i,j) = u2
        v(i,j) = v2
        ene(i,j) = ene2
        T(i,j) = T2
        p(i,j) = p2
      end if
      c(i,j) = dsqrt(g*R*T(i,j))
    end do
  end do
end subroutine conditions

!------���W�ϊ�------
subroutine coordinate
  do j = 0,JM
    do i = 1,IM-1
      xe(i,j) = (x(i+1,j)-x(i-1,j))/2.0d0
      ye(i,j) = (y(i+1,j)-y(i-1,j))/2.0d0
    end do
  end do

  do j = 1,JM-1
    do i = 0,IM
      xh(i,j) = (x(i,j+1)-x(i,j-1))/2.0d0
      yh(i,j) = (y(i,j+1)-y(i,j-1))/2.0d0
    end do
  end do

  do j = 0,JM
    xe(0,j) = (-3.0*x(0,j)+4.0*x(1,j)-x(2,j))/2.0d0
    xe(IM,j) = (3.0*x(IM,j)-4.0*x(IM-1,j)+x(IM-2,j))/2.0d0
    ye(0,j) = (-3.0*y(0,j)+4.0*y(1,j)-y(2,j))/2.0d0
    ye(IM,j) = (3.0*y(IM,j)-4.0*y(IM-1,j)+y(IM-2,j))/2.0d0
  end do

  do i = 0,IM
    xh(i,0) = (-3.0*x(i,0)+4.0*x(i,1)-x(i,2))/2.0d0
    xh(i,JM) = (3.0*x(i,JM)-4.0*x(i,JM-1)+x(i,JM-2))/2.0d0
    yh(i,0) = (-3.0*y(i,0)+4.0*y(i,1)-y(i,2))/2.0d0
    yh(i,JM) = (3.0*y(i,JM)-4.0*y(i,JM-1)+y(i,JM-2))/2.0d0
  end do

  do j = 0,JM
    do i = 0,IM
      Jac(i,j) = 1.0/(xe(i,j)*yh(i,j)-ye(i,j)*xh(i,j))
      xix(i,j) = yh(i,j)*Jac(i,j)
      xiy(i,j) = -xh(i,j)*Jac(i,j)
      ex(i,j) = -ye(i,j)*Jac(i,j)
      ey(i,j) = xe(i,j)*Jac(i,j)
    end do
  end do

  do j = 0,JM
    do i = 0,IM
      Q(1,i,j) = rho(i,j)/Jac(i,j) !�ۑ��ϐ��x�N�g��
      Q(2,i,j) = (rho(i,j)*u(i,j))/Jac(i,j)
      Q(3,i,j) = (rho(i,j)*v(i,j))/Jac(i,j)
      Q(4,i,j) = ene(i,j)/Jac(i,j)
    end do
  end do
end subroutine coordinate

!----------�̕���TVD----------
subroutine tvd_xi
  do j = 1,JM-1
    do i = 0,IM-1
      RWL = dSQRT(Q(1,I,J)*Jac(I,J))/(dSQRT(Q(1,I+1,J)*Jac(I+1,J))+dSQRT(Q(1,I,J)*Jac(I,J)))
      RWR = dSQRT(Q(1,I+1,J)*Jac(I+1,J))/(dSQRT(Q(1,I+1,J)*Jac(I+1,J))+dSQRT(Q(1,I,J)*Jac(I,J)))
      AKX = 0.5d0*(xix(I,J)+xix(I+1,J)) !���ԓ_��kx
      AKY = 0.5d0*(xiy(I,J)+xiy(I+1,J)) !���ԓ_��ky
      AJACM = 0.5d0*(Jac(I,J)+Jac(I+1,J)) !���ԓ_��Jac

      UM = RWL*U(I,J)+RWR*U(I+1,J) !���ԓ_��u
      VM = RWL*V(I,J)+RWR*V(I+1,J) !���ԓ_��v
      DUM = RWL*Q(4,I,J)/Q(1,I,J)+RWR*Q(4,I+1,J)/Q(1,I+1,J) !e/��
      CM = dSQRT(G*(G-1.0d0)*abs(DUM-0.5d0*(UM**2.0+VM**2.0))) !���ԓ_�̉���c
      !NOMENCLATURE
      PHI = 0.5d0*(G-1.0d0)*(UM**2.0+VM**2.0) !��
      BETA = 1.0d0/(2.0d0*CM**2.0)
      AKXT = AKX/dSQRT(AKX**2.0+AKY**2.0)
      AKYT = AKY/dSQRT(AKX**2.0+AKY**2.0)
      THIT = AKXT*UM+AKYT*VM

      !RIGHT ENGEN-VECTORS(�E�ŗL�x�N�g���j
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
      RR(4,1,I) = PHI/(G-1.0d0)
      RR(4,2,I) = AKYT*UM-AKXT*VM
      RR(4,3,I) = (PHI+CM**2)/(G-1.0d0)+CM*THIT
      RR(4,4,I) = (PHI+CM**2)/(G-1.0d0)-CM*THIT
      !print *,RR(4,1,i),RR(4,2,i),RR(4,3,i),RR(4,4,i)

      !INVERS OR RIGHT EIGEN-VECTORS(�E�ŗL�x�N�g���̋t�s��j
      RI(1,1,I) = 1.0d0-PHI/CM**2.0
      RI(1,2,I) = (G-1.0d0)*UM/CM**2.0
      RI(1,3,I) = (G-1.0d0)*VM/CM**2.0
      RI(1,4,I) = -(G-1.0d0)/CM**2.0
      RI(2,1,I) = -AKYT*UM+AKXT*VM
      RI(2,2,I) = AKYT
      RI(2,3,I) = -AKXT
      RI(2,4,I) = 0.0d0
      RI(3,1,I) = BETA*(PHI-CM*THIT)
      RI(3,2,I) = BETA*(AKXT*CM-(G-1.0d0)*UM)
      RI(3,3,I) = BETA*(AKYT*CM-(G-1.0d0)*VM)
      RI(3,4,I) = BETA*(G-1.0d0)
      RI(4,1,I) = BETA*(PHI+CM*THIT)
      RI(4,2,I) = -BETA*(AKXT*CM+(G-1.0d0)*UM)
      RI(4,3,I) = -BETA*(AKYT*CM+(G-1.0d0)*VM)
      RI(4,4,I) = BETA*(G-1.0d0)

      !ENGEN-VALUES AT INTERMEDIATE(�ŗL�l�j
      EIGM(1,I) = AKX*UM+AKY*VM
      EIGM(2,I) = EIGM(1,I)
      EIGM(3,I) = EIGM(1,I)+CM*dSQRT(AKX**2.0+AKY**2.0)
      EIGM(4,I) = EIGM(1,I)-CM*dSQRT(AKX**2.0+AKY**2.0)

      !ALPHA�̓��o
      do K = 1,4
        D(K) = (Q(K,I+1,J)*Jac(I+1,J)-Q(K,I,J)*Jac(I,J))/AJACM
      end do

      do K = 1,4
        ALPHA(K,I) = RI(K,1,I)*D(1)+RI(K,2,I)*D(2)+RI(K,3,I)*D(3)+RI(K,4,I)*D(4)
      end do
    end do

    do k = 1,4
      ALPHA(K,-1) = ALPHA(K,0)
      ALPHA(K,IM) = ALPHA(K,IM-1)
    end do

    do i = 0,IM
      DUM = dSQRT(G*R*T(I,J)*(xix(I,J)**2+xiy(I,J)**2))
      EIG(1,I) = xix(I,J)*U(I,J)+xiy(I,J)*V(I,J)
      EIG(2,I) = EIG(1,I)
      EIG(3,I) = EIG(1,I)+DUM
      EIG(4,I) = EIG(1,I)-DUM
    end do

    !CAL LIMITER FUNCTION G(����֐�)
    do K = 1,4
      do I = 0,IM
        !MINMOD LIMITER
        s = dSIGN(1.0d0,ALPHA(K,I))
        GG(K,I) = s*dMAX1(0.0,dMIN1(dabs(ALPHA(K,I)),s*ALPHA(K,I-1)))
      end do
    end do
    do K = 1,4
      do I = 0,IM-1
        !GAMMA
        DELTA = 1.0d-10
        !DELTA=dMAX1(0.0,EIGM(K,I)-EIG(K,I),EIG(K,I+1)-EIGM(K,I))
        if (ALPHA(K,I)/=0.0) then
          GAMMA = 0.5d0*FPSI(EIGM(K,I),DELTA)*(GG(K,I+1)-GG(K,I))/ALPHA(K,I)
        else
          GAMMA = 0.0d0
        end if
        !PHI(��)
        PHIM(K,I) = 0.5d0*FPSI(EIGM(K,I),DELTA)*(GG(K,I+1)+GG(K,I))-FPSI(EIGM(K,I)+GAMMA,DELTA)*ALPHA(K,I)
      end do
    end do
    !CONVECTION COMPORNENTS(�v�Z��ԏ��E�v�Z)
    do I = 0,IM
      UCONT = xix(I,J)*U(I,J)+xiy(I,J)*V(I,J) !(���ϑ��x)
      E(1,I) = Q(1,I,J)*UCONT
      E(2,I) = Q(2,I,J)*UCONT+P(I,J)*xix(I,J)/Jac(i,j)
      E(3,I) = Q(3,I,J)*UCONT+P(I,J)*xiy(I,J)/Jac(i,j)
      E(4,I) = Q(4,I,J)*UCONT+P(I,J)*UCONT/Jac(i,j)
    end do
    !XI-DIRECTION CONVECTION FLUX
    do K = 1,4
      do I = 0,IM-1
        EE(K,I,J) = 0.5d0*(E(K,I)+E(K,I+1)+RR(K,1,I)*PHIM(1,I)+RR(K,2,I)*PHIM(2,I)&
                                &+RR(K,3,I)*PHIM(3,I)+RR(K,4,I)*PHIM(4,I))
      end do
    end do
  end do
  return
end subroutine tvd_xi

!*****FPSI�֐��̒�`*****
double precision function FPSI(Z,DELTA)
  double precision Z,DELTA
  if (dabs(Z)>=DELTA) then
    FPSI = dabs(Z)
  else
    FPSI = 0.5d0*(Z**2.0+DELTA**2.0)/DELTA
  end if
  return
end function

!----------�ŕ���TVD----------
subroutine tvd_eta
  do I = 1,IM-1
    do J = 0,JM-1
      RWL = dSQRT(Q(1,I,J)*Jac(I,J))/(dSQRT(Q(1,I,J+1)*Jac(I,J+1))+dSQRT(Q(1,I,J)*Jac(I,J)))
      RWR = dSQRT(Q(1,I,J+1)*Jac(I,J+1))/(dSQRT(Q(1,I,J+1)*Jac(I,J+1))+dSQRT(Q(1,I,J)*Jac(I,J)))
      AKX = 0.5d0*(ex(I,J)+ex(I,J+1))
      AKY = 0.5d0*(ey(I,J)+ey(I,J+1))
      AJACM = 0.5d0*(Jac(I,J)+Jac(I,J+1))
      UM = RWL*U(I,J)+RWR*U(I,J+1)
      VM = RWL*V(I,J)+RWR*V(I,J+1)
      DUM = RWL*Q(4,I,J)/Q(1,I,J)+RWR*Q(4,I,J+1)/Q(1,I,J+1)
      CM = dSQRT(G*(G-1.0d0)*dabs(DUM-0.5d0*(UM**2.0+VM**2.0)))

      !NOMENCLATURE
      PHI = 0.5d0*(G-1.0d0)*(UM**2.0+VM**2.0)
      BETA = 1.0d0/(2.0d0*CM**2.0)
      AKXT = AKX/dSQRT(AKX**2.0+AKY**2.0)
      AKYT = AKY/dSQRT(AKX**2.0+AKY**2.0)
      THIT = AKXT*UM+AKYT*VM

      !RIGHT ENGEN-VECTORS
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
      RR(4,1,J) = PHI/(G-1.0d0)
      RR(4,2,J) = AKYT*UM-AKXT*VM
      RR(4,3,J) = (PHI+CM**2.0)/(G-1.0d0)+CM*THIT
      RR(4,4,J) = (PHI+CM**2.0)/(G-1.0d0)-CM*THIT

      !INVERS OR RIGHT EIGEN-VECTORS
      RI(1,1,J) = 1.0d0-PHI/CM**2.0
      RI(1,2,J) = (G-1.0d0)*UM/CM**2.0
      RI(1,3,J) = (G-1.0d0)*VM/CM**2.0
      RI(1,4,J) = -(G-1.0d0)/CM**2.0
      RI(2,1,J) = -AKYT*UM+AKXT*VM
      RI(2,2,J) = AKYT
      RI(2,3,J) = -AKXT
      RI(2,4,J) = 0.0d0
      RI(3,1,J) = BETA*(PHI-CM*THIT)
      RI(3,2,J) = BETA*(AKXT*CM-(G-1.0d0)*UM)
      RI(3,3,J) = BETA*(AKYT*CM-(G-1.0d0)*VM)
      RI(3,4,J) = BETA*(G-1.0d0)
      RI(4,1,J) = BETA*(PHI+CM*THIT)
      RI(4,2,J) = -BETA*(AKXT*CM+(G-1.0d0)*UM)
      RI(4,3,J) = -BETA*(AKYT*CM+(G-1.0d0)*VM)
      RI(4,4,J) = BETA*(G-1.0d0)

      !ENGEN-VALUES AT INTERMEDIATE
      FIGM(1,J) = AKX*UM+AKY*VM
      FIGM(2,J) = FIGM(1,J)
      FIGM(3,J) = FIGM(1,J)+CM*dSQRT(AKX**2.0+AKY**2.0)
      FIGM(4,J) = FIGM(1,J)-CM*dSQRT(AKX**2.0+AKY**2.0)

      !ALPHA
      do K = 1,4
        D(K) = (Q(K,I,J+1)*Jac(I,J+1)-Q(K,I,J)*Jac(I,J))/AJACM
      end do

      do K = 1,4
        ALPHA(K,J) = RI(K,1,J)*D(1)+RI(K,2,J)*D(2)+RI(K,3,J)*D(3)+RI(K,4,J)*D(4)
      end do
    end do

    do K = 1,4
      ALPHA(K,-1) = ALPHA(K,0)
      ALPHA(K,JM) = ALPHA(K,JM-1)
    end do

    do J = 0,JM
      DUM = dSQRT(G*R*T(I,J)*(ex(I,J)**2.0+ey(I,J)**2.0))
      FIG(1,J) = ex(I,J)*U(I,J)+ey(I,J)*V(I,J)
      FIG(2,J) = FIG(1,J)
      FIG(3,J) = FIG(1,J)+DUM
      FIG(4,J) = FIG(1,J)-DUM
    end do
    !CAL LIMITER FUNCTION G
    do K = 1,4
      do J = 0,JM
        !MINMOD LIMITER �����֐�
        s = dSIGN(1.0d0,ALPHA(K,J))
        GG(K,J) = s*dMAX1(0.0d0,dMIN1(dabs(ALPHA(K,J)),s*ALPHA(K,J-1)))
      end do
    end do

    do K = 1,4
      do J = 0,JM-1
        !GAMMA
        DELTA = 1.0d-10
        !DELTA=dMAX1(0.0,EIGM(K,I)-EIG(K,I),EIG(K,I+1)-EIGM(K,I))
        if (ALPHA(K,J)/=0.0) then
          GAMMA = 0.5d0*FPSI(FIGM(K,J),DELTA)*(GG(K,J+1)-GG(K,J))/ALPHA(K,J)
        else
          GAMMA = 0.0d0
        end if
        !PHI
        PHIM(K,J) = 0.5d0*FPSI(FIGM(K,J),DELTA)*(GG(K,J+1)+GG(K,J))-FPSI(FIGM(K,J)+GAMMA,DELTA)*ALPHA(K,J)
      end do
    end do

    !CONVECTION COMPORNENTS
    do J = 0,JM
      UCONT = ex(I,J)*U(I,J)+ey(I,J)*V(I,J)
      F(1,J) = Q(1,I,J)*UCONT
      F(2,J) = Q(2,I,J)*UCONT+P(I,J)*ex(I,J)/Jac(i,j)
      F(3,J) = Q(3,I,J)*UCONT+P(I,J)*ey(I,J)/Jac(i,j)
      F(4,J) = Q(4,I,J)*UCONT+P(I,J)*UCONT/Jac(i,j)
    end do

    !XI-DIRECTION CONVECTION FLUX
    do K = 1,4
      do J = 0,JM-1
        FF(K,I,J) = 0.5d0*(F(K,J)+F(K,J+1)+RR(K,1,J)*PHIM(1,J)+RR(K,2,J)*PHIM(2,J)&
                                &+RR(K,3,J)*PHIM(3,J)+RR(K,4,J)*PHIM(4,J))
      end do
    end do
  end do
  return
end subroutine tvd_eta

!------�T�U�[�����h------
double precision function Sutherland(T)
  double precision T
  Sutherland = myu0*((T/T0)**1.5)*(T0+SS)/(T+SS)
end function

!------�S�����i�̕����j------
subroutine viscosity_xi
  double precision,dimension(0:IM-1,0:JM-1) :: xix_r,xiy_r,ex_r,ey_r,Jac_r,u_r,v_r,T_r,ue_r,ve_r,Te_r,uh_r,vh_r,Th_r
  do j = 1,JM-1
    do i = 0,IM-1
      xix_r(i,j) = (xix(i,j)+xix(i+1,j))/2.0d0
      xiy_r(i,j) = (xiy(i,j)+xiy(i+1,j))/2.0d0
      ex_r(i,j) = (ex(i,j)+ex(i+1,j))/2.0d0
      ey_r(i,j) = (ey(i,j)+ey(i+1,j))/2.0d0
      Jac_r(i,j) = (Jac(i,j)+Jac(i+1,j))/2.0d0
      u_r(i,j) = (u(i,j)+u(i+1,j))/2.0d0
      v_r(i,j) = (v(i,j)+v(i+1,j))/2.0d0
      T_r(i,j) = (t(i,j)+t(i+1,j))/2.0d0

      ue_r(i,j) = u(i+1,j)-u(i,j)
      ve_r(i,j) = v(i+1,j)-v(i,j)
      Te_r(i,j) = t(i+1,j)-t(i,j)

      uh_r(i,j) = (((u(i,j+1)+u(i+1,j+1))/2.0d0)-((u(i,j-1)+u(i+1,j-1))/2.0d0))/2.0d0
      vh_r(i,j) = (((v(i,j+1)+v(i+1,j+1))/2.0d0)-((v(i,j-1)+v(i+1,j-1))/2.0d0))/2.0d0
      Th_r(i,j) = (((t(i,j+1)+t(i+1,j+1))/2.0d0)-((t(i,j-1)+t(i+1,j-1))/2.0d0))/2.0d0
      !print *,Jac_r(i,j)
    end do
  end do
  do j = 1,JM-1
    do i = 0,IM-1
      myu = Sutherland(T_r(i,j))
      txx = 2.0d0/3.0d0*myu*(2.0d0*(ue_r(i,j)*xix_r(i,j)+uh_r(i,j)*ex_r(i,j))-(ve_r(i,j)*xiy_r(i,j)+vh_r(i,j)*ey_r(i,j)))
      tyy = 2.0d0/3.0d0*myu*(2.0d0*(ve_r(i,j)*xiy_r(i,j)+vh_r(i,j)*ey_r(i,j))-(ue_r(i,j)*xix_r(i,j)+uh_r(i,j)*ex_r(i,j)))
      txy = myu*((ue_r(i,j)*xiy_r(i,j)+uh_r(i,j)*ey_r(i,j))+(ve_r(i,j)*xix_r(i,j)+vh_r(i,j)*ex_r(i,j)))
      tyx = txy
      R4 = txx*u_r(i,j)+txy*v_r(i,j)+(myu/(Pr*(g-1.0d0)))*g*R*(Te_r(i,j)*xix_r(i,j)+Th_r(i,j)*ex_r(i,j))
      S4 = tyx*u_r(i,j)+tyy*v_r(i,j)+(myu/(Pr*(g-1.0d0)))*g*R*(Te_r(i,j)*xiy_r(i,j)+Th_r(i,j)*ey_r(i,j))
      RRR(1,i,j) = 0.0d0
      RRR(2,i,j) = (xix_r(i,j)*txx+xiy_r(i,j)*txy)/Jac_r(i,j)
      RRR(3,i,j) = (xix_r(i,j)*tyx+xiy_r(i,j)*tyy)/Jac_r(i,j)
      RRR(4,i,j) = (xix_r(i,j)*r4+xiy_r(i,j)*s4)/Jac_r(i,j)
    end do
  end do
end subroutine viscosity_xi

!------�S�����i�ŕ����j------
subroutine viscosity_eta
  double precision,dimension(0:IM-1,0:JM-1) :: xix_s,xiy_s,ex_s,ey_s,Jac_s,u_s,v_s,T_s,ue_s,ve_s,Te_s,uh_s,vh_s,Th_s
  do i = 1,IM-1
    do j = 0,JM-1
      xix_s(i,j) = (xix(i,j)+xix(i,j+1))/2.0d0
      xiy_s(i,j) = (xiy(i,j)+xiy(i,j+1))/2.0d0
      ex_s(i,j) = (ex(i,j)+ex(i,j+1))/2.0d0
      ey_s(i,j) = (ey(i,j)+ey(i,j+1))/2.0d0
      Jac_s(i,j) = (Jac(i,j)+Jac(i,j+1))/2.0d0
      u_s(i,j) = (u(i,j)+u(i,j+1))/2.0d0
      v_s(i,j) = (v(i,j)+v(i,j+1))/2.0d0
      T_s(i,j) = (t(i,j)+t(i,j+1))/2.0d0

      uh_s(i,j) = u(i,j+1)-u(i,j)
      vh_s(i,j) = v(i,j+1)-v(i,j)
      Th_s(i,j) = t(i,j+1)-t(i,j)

      ue_s(i,j) = (((u(i+1,j)+u(i+1,j+1))/2.0d0)-((u(i-1,j)+u(i-1,j+1))/2.0d0))/2.0d0
      ve_s(i,j) = (((v(i+1,j)+v(i+1,j+1))/2.0d0)-((v(i-1,j)+v(i-1,j+1))/2.0d0))/2.0d0
      Te_s(i,j) = (((t(i+1,j)+t(i+1,j+1))/2.0d0)-((t(i-1,j)+t(i-1,j+1))/2.0d0))/2.0d0
    end do
  end do

  do i = 1,IM-1
    do j = 0,JM-1
      myu = Sutherland(T_s(i,j))
      txx = 2.0d0/3.0d0*myu*(2.0d0*(ue_s(i,j)*xix_s(i,j)+uh_s(i,j)*ex_s(i,j))-(ve_s(i,j)*xiy_s(i,j)+vh_s(i,j)*ey_s(i,j)))
      tyy = 2.0d0/3.0d0*myu*(2.0d0*(ve_s(i,j)*xiy_s(i,j)+vh_s(i,j)*ey_s(i,j))-(ue_s(i,j)*xix_s(i,j)+uh_s(i,j)*ex_s(i,j)))
      txy = myu*((ue_s(i,j)*xiy_s(i,j)+uh_s(i,j)*ey_s(i,j))+(ve_s(i,j)*xix_s(i,j)+vh_s(i,j)*ex_s(i,j)))
      tyx = txy
      R4 = txx*u_s(i,j)+txy*v_s(i,j)+(myu/(Pr*(g-1.0d0)))*g*R*(Te_s(i,j)*xix_s(i,j)+Th_s(i,j)*ex_s(i,j))
      S4 = tyx*u_s(i,j)+tyy*v_s(i,j)+(myu/(Pr*(g-1.0d0)))*g*R*(Te_s(i,j)*xiy_s(i,j)+Th_s(i,j)*ey_s(i,j))
      SSS(1,i,j) = 0.0d0
      SSS(2,i,j) = (ex_s(i,j)*txx+ey_s(i,j)*txy)/Jac_s(i,j)
      SSS(3,i,j) = (ex_s(i,j)*tyx+ey_s(i,j)*tyy)/Jac_s(i,j)
      SSS(4,i,j) = (ex_s(i,j)*r4+ey_s(i,j)*s4)/Jac_s(i,j)
    end do
  end do

end subroutine viscosity_eta

!------Runge-Kutta------
subroutine Runge_Kutta
  do j = 1,JM-1
    do i = 1,IM-1
      do k = 1,4
        QQ(k,i,j,rk) = QQ(k,i,j,0)-(1.0d0/dble(5-rk))*dt*(-EE(k,i-1,j)+EE(k,i,j)-&
                                        &FF(k,i,j-1)+FF(k,i,j)+(RRR(k,i-1,j)-RRR(k,i,j))+&
                                        &(SSS(k,i,j-1)-SSS(k,i,j)))
        Q(k,i,j) = QQ(k,i,j,rk)
      end do
    end do
  end do
end subroutine Runge_Kutta

!------���ʂ̌v�Z------
subroutine calculation
  do j = 1,JM-1
    do i = 1,IM-1
      un = u(i,j)
      vn = v(i,j)
      rho(i,j) = Jac(i,j)*Q(1,i,j)
      u(i,j) = Jac(i,j)*Q(2,i,j)/rho(i,j)
      v(i,j) = Jac(i,j)*Q(3,i,j)/rho(i,j)
      ene(i,j) = Jac(i,j)*Q(4,i,j)
      enebar(i,j) = ene(i,j)/rho(i,j)-((u(i,j)**2.0+v(i,j)**2.0)/2.0d0)
      T(i,j) = enebar(i,j)*(g-1.0d0)/R
      p(i,j) = (g-1.0d0)*rho(i,j)*enebar(i,j)
      ddu = dabs((u(i,j)-un)/un)
      ddv = dabs((v(i,j)-vn)/vn)
      if (ddu>du) then
        du = ddu
      end if
      if (ddv>dv) then
        dv = ddv
      end if
    end do
  end do
end subroutine calculation

!------���E����------
subroutine boundary
  do i = 0,IM !��ǁi�m�C�}�������j
    rho(i,JM) = rho(i,JM-1)
    u(i,JM) = u(i,JM-1)
    v(i,JM) = v(i,JM-1)
    p(i,JM) = p(i,JM-1)
    T(i,JM) = T(i,JM-1)
    ene(i,JM) = p(i,JM)/(g-1.0d0)+rho(i,JM)*((u(i,JM)**2.0)+(v(i,JM)**2.0))/2.0d0
    Q(1,i,JM) = rho(i,JM)/Jac(i,JM)
    Q(2,i,JM) = rho(i,JM)*u(i,JM)/Jac(i,JM)
    Q(3,i,JM) = rho(i,JM)*v(i,JM)/Jac(i,JM)
    Q(4,i,JM) = ene(i,JM)/Jac(i,JM)

    !����(�S�����l�����邽�߁Cu=0���̓m�C�}������)
    rho(i,0) = rho(i,1)
    u(i,0) = 0.0d0
    v(i,0) = 0.0d0
    p(i,0) = p(i,1)
    T(i,0) = T(i,1)
    ene(i,0) = p(i,0)/(g-1.0d0)+rho(i,0)*((u(i,0)**2.0+v(i,0)**2.0))/2.0d0
    Q(1,i,0) = rho(i,0)/Jac(i,0)
    Q(2,i,0) = rho(i,0)*u(i,0)/Jac(i,0)
    Q(3,i,0) = rho(i,0)*v(i,0)/Jac(i,0)
    Q(4,i,0) = ene(i,0)/Jac(i,0)
  end do

  do j = 1,JM-1 !���o(�m�C�}������)
    slope = (x(IM,j)-x(IM-2,j))/(x(IM-2,j)-x(IM-1,j))
    rho(IM,j) = rho(IM-2,j)+slope*(rho(IM-2,j)-rho(IM-1,j))
    u(IM,j) = u(IM-2,j)+slope*(u(IM-2,j)-u(IM-1,j))
    v(IM,j) = v(IM-2,j)+slope*(v(IM-2,j)-v(IM-1,j))
    p(IM,j) = p(IM-2,j)+slope*(p(IM-2,j)-p(IM-1,j))
    ene(IM,j) = p(im,j)/(g-1.0d0)+rho(im,j)*((u(im,j)**2.0)+(v(im,j)**2.0))/2.0d0
    T(IM,j) = T(IM-2,j)+slope*(T(IM-2,j)-T(IM-1,j))

    Q(1,IM,j) = rho(IM,j)/Jac(IM,j)
    Q(2,IM,j) = rho(IM,j)*u(IM,j)/Jac(IM,j)
    Q(3,IM,j) = rho(IM,j)*v(IM,j)/Jac(IM,j)
    Q(4,IM,j) = ene(IM,j)/Jac(IM,j)
  end do
end subroutine boundary

!--------------------
end program assyuku7
