program comp8
!�ϐ��̒�`*************************************************
implicit none
double precision,parameter::pi=dacos(0.0d0)*2 !��

double precision,parameter::H=0.2d0 !������i�����j
integer,parameter::IM=50 !x�����̈搔
integer,parameter::JM=50 !y�����̈搔

double precision,dimension(0:IM,0:JM)::x=0.0d0 !�i�q�_(i,j)�ɂ�����x���W
double precision,dimension(0:IM,0:JM)::y=0.0d0 !�i�q�_(i,j)�ɂ�����y���W
double precision,parameter::cr=1.05d0 !���䐔�����
double precision::ftx,fty !���䐔�񏉍�

double precision,dimension(0:IM,0:JM)::xxi=0.0d0,xeta=0.0d0,yxi=0.0d0,yeta=0.0d0 !�́�xi�C�Ł�eta
double precision,dimension(0:IM,0:JM)::JJ !���R�r�A��
double precision,dimension(0:IM,0:JM)::xix=0.0d0,xiy=0.0d0,etax=0.0d0,etay=0.0d0 !���W�ϊ��p�Δ���

integer i,j,n,l,nn
integer,parameter::nmax=100000 !�ő�v�Z��
double precision,parameter::epsln=1.0d-4 !��������l
double precision dq,ddq
double precision,dimension(0:IM,0:JM)::qp=1.0d0

double precision,parameter::R=287.1d0 !�K�X�萔
double precision,parameter::g=1.4d0 !��M��

double precision::M=2.9d0 !�����}�b�n��
double precision::un !���xu�̒l�ێ��p

double precision,dimension(0:IM,0:JM)::rho=0.0d0 !���x
double precision,dimension(0:IM,0:JM)::u=0.0d0 !x�������x
double precision,dimension(0:IM,0:JM)::v=0.0d0 !y�������x
double precision,dimension(0:IM,0:JM)::p=0.0d0 !����
double precision,dimension(0:IM,0:JM)::T=0.0d0 !���x
double precision,dimension(0:IM,0:JM)::ene=0.0d0 !�G�l���M�[
double precision,dimension(0:Im,0:JM)::k=0.0d0 !�����G�l���M�[
double precision,dimension(0:IM,0:JM)::ep=0.0d0 !�U�헦

double precision,dimension(0:IM,0:JM)::mu !�S���W��
double precision,dimension(0:IM,0:JM)::nu !���S���W��
double precision,parameter::Ret=500.d0 !�������C�m���Y��

double precision,dimension(1:6,0:IM,0:JM)::Q=0.0d0 !�ۑ��ʃx�N�g��
double precision,dimension(1:6,0:IM,0:JM)::E=0.0d0 !x���������x�N�g��
double precision::Eh !E�̒l�ێ��p
double precision,dimension(1:6,0:IM,0:JM)::F=0.0d0 !y���������x�N�g��
double precision,dimension(1:6,0:IM,0:JM,0:4)::QQ=0.0d0 !�ۑ��ʃx�N�g���i�����Q�N�b�^�j
double precision,dimension(1:6,0:IM,0:JM)::EE=0.0d0 !x���������x�N�g���iTVD�X�L�[���jE~
double precision,dimension(1:6,0:IM,0:JM)::FF=0.0d0 !y���������x�N�g���iTVD�X�L�[���jF~

double precision::lCxi,lCeta,lDxi,lDeta,lep !��C�́C��C�ŁC��D�́C��D�ŁC�Ƀ�
double precision,parameter::CN=0.2d0 !�N�[������
double precision,dimension(0:IM,0:JM)::dt !�������ݕ�

double precision kx,ky,JJm,um,vm,dum,cm,kxt,kyt,phi,thetat,beta,S,delta,gamma
double precision,dimension(1:4,1:4,0:IM)::RR,RI !�s��
double precision,dimension(1:6,0:IM)::a !�ŗL�l
double precision,dimension(1:6,0:IM)::gg !�����֐�
double precision,dimension(1:6,-1:IM)::alpha !��
double precision,dimension(1:6)::D !���v�Z�p
double precision,dimension(1:6,0:IM)::phim !��

double precision,dimension(0:IM,0:JM)::muS !�S���W���i�X�^�K�[�h�i�q�j
double precision,parameter::T0=293.15d0 !���x�i�萔�j
double precision,parameter::mu0=1.0d-3 !�S���W���i�萔�j
double precision,parameter::Sl=110.0d0 !�T�U�[�����h�萔
double precision,parameter::Cmu=0.09d0 !�������萔
double precision,dimension(0:IM,0:JM)::mutS !�Q�S���W���i�X�^�K�[�h�i�q�j
double precision,dimension(0:IM,0:JM)::xixS,xiyS,etaxS,etayS,JJS,rhoS,uS,vS,TS,kS,epS,&
                                       uxiS,vxiS,TxiS,kxiS,epxiS,uetaS,vetaS,TetaS,ketaS,epetaS !��ij�p�Δ����i�X�^�K�[�h�i�q�j
double precision::SxxS,SxyS,SyyS !Sij�i�X�^�K�[�h�i�q�j
double precision::txxS,txyS,tyyS !��xx,��xy,��yy�i�X�^�K�[�h�i�q�j
double precision,dimension(1:6,0:IM,0:JM)::RvS=0.0d0,SvS=0.0d0,RvSS=0.0d0,SvSS=0.0d0 !�S�����x�N�g��R,S�i�X�^�K�[�h�i�q�j
double precision,parameter::Prt=0.9d0 !�萔
double precision,parameter::sigk=1.0d0,sigep=1.3d0 !�萔

double precision,dimension(0:IM,0:JM)::mut !�Q�S���W��
double precision,dimension(0:IM,0:JM)::uxi,vxi,ueta,veta
double precision::Sxx,Sxy,Syy !Sij
double precision::txx,txy,tyy !��xx,��xy,��yy
double precision,dimension(1:6,0:IM,0:JM)::PP=0.0d0 !�������x�N�g��P
double precision::PRO
double precision,parameter::Cep1=1.44d0,Cep2=1.92d0 !�萔

double precision::kp,epp,yp,ypp !�ǂ�����_�ڂ̒l
double precision::ut !���C���x
double precision::kap=0.42d0,Ec=9.7d0 !�萔

double precision::utau !u��
double precision,dimension(0:IM,0:JM)::uplus,yplus !u^+,y^+


!�����ʍ��W�̐ݒ�*******************************************
ftx=45.0d0*H*(cr-1.0d0)/(cr**40.0d0-1.0d0) !���䐔�񏉍��ix�����j
fty=H*(cr-1.0d0)/(cr**real(JM)-1.0d0) !���䐔�񏉍��iy�����j
do i=0,IM
 do j=0,JM
  if(i<=10) then
   x(i,j)=dble(i)*50.0d0*H/dble(IM)/2.0d0
  else
   x(i,j)=x(10,j)+ftx*(cr**real(i-10)-1.0d0)/(cr-1.0d0)
  end if
  y(i,j)=fty*(cr**real(j)-1.0d0)/(cr-1.0d0)
 end do
end do


!���W�ϊ�***************************************************
do j=0,IM
 do i=0,JM
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


!��������***************************************************
do j=0,IM
 do i=0,JM
  rho(i,j)=1.2d0
  T(i,j)=293.15d0
  if(i<10.or.j/=0) then
   u(i,j)=M*dsqrt(g*R*T(i,j))
  end if
  p(i,j)=rho(i,j)*R*T(i,j)
  ene(i,j)=rho(i,j)*(R*T(i,j)/(g-1.0d0)+(u(i,j)**2.0+v(i,j)**2.0)/2.0d0)
  mu(i,j)=1.82d-5
  nu(i,j)=mu(i,j)/rho(i,j)
  k(i,j)=1.5d0*(0.01d0*u(i,j))**2.0
  ep(i,j)=k(i,j)**2.0/(nu(i,j)*Ret)

  Q(1,i,j)=rho(i,j)
  Q(2,i,j)=rho(i,j)*u(i,j)
  Q(3,i,j)=rho(i,j)*v(i,j)
  Q(4,i,j)=ene(i,j)
  Q(5,i,j)=rho(i,j)*k(i,j)
  Q(6,i,j)=rho(i,j)*ep(i,j)
  E(1,i,j)=rho(i,j)*u(i,j)
  E(2,i,j)=p(i,j)+rho(i,j)*u(i,j)**2.0
  E(3,i,j)=rho(i,j)*u(i,j)*v(i,j)
  E(4,i,j)=(ene(i,j)+p(i,j))*u(i,j)
  E(5,i,j)=rho(i,j)*k(i,j)*u(i,j)
  E(6,i,j)=rho(i,j)*ep(i,j)*u(i,j)
  F(1,i,j)=rho(i,j)*v(i,j)
  F(2,i,j)=rho(i,j)*u(i,j)*v(i,j)
  F(3,i,j)=p(i,j)+rho(i,j)*v(i,j)**2.0
  F(4,i,j)=(ene(i,j)+p(i,j))*v(i,j)
  F(5,i,j)=rho(i,j)*k(i,j)*v(i,j)
  F(6,i,j)=rho(i,j)*ep(i,j)*v(i,j)

  do l=1,6 !�e��������ʍ��W�n�ɕϊ�
   Q(l,i,j)=Q(l,i,j)/JJ(i,j)
   QQ(l,i,j,0)=Q(l,i,j)
   Eh=E(l,i,j)
   E(l,i,j)=(xix(i,j)*E(l,i,j)+xiy(i,j)*F(l,i,j))/JJ(i,j)
   F(l,i,j)=(etax(i,j)*Eh+etay(i,j)*F(l,i,j))/JJ(i,j)
  end do

 end do
end do


!�v�Z*******************************************************
n=0
dq=epsln+1.0d0
do while(dq>epsln.and.n<nmax)
 n=n+1
 dq=0.0d0
 !�Ǐ����ԍ��ݖ@--------------------------------------------
 do j=0,JM
  do i=0,IM
   mu(i,j)=((T(i,j)/T0)**1.5)*(T0+Sl)/(T(i,j)+Sl)*mu0
   mut(i,j)=(Cmu*rho(i,j)*k(i,j)**2.0)/ep(i,j)
   !�Η���
   lCxi=dabs(E(1,i,j)*JJ(i,j)/rho(i,j))+dsqrt(g*R*T(i,j)*(xix(i,j)**2.0+xiy(i,j)**2.0))
   lCeta=dabs(F(1,i,j)*JJ(i,j)/rho(i,j))+dsqrt(g*R*T(i,j)*(etax(i,j)**2.0+etay(i,j)**2.0))
   !�g�U��
   lDxi=2.0d0*(mu(i,j)+mut(i,j))*(xix(i,j)**2.0+xiy(i,j)**2.0)/rho(i,j)
   lDeta=2.0d0*(mu(i,j)+mut(i,j))*(etax(i,j)**2.0+etay(i,j)**2.0)/rho(i,j)
   !����̎U�퍀
   lep=ep(i,j)/k(i,j)
   !���ԍ���
   dt(i,j)=CN*dmin1(1.0d0/(lCxi+lDxi+lep),1.0d0/(lCeta+lDeta+lep))
  end do
 end do

 do nn=1,4
  !�̕���TVD�X�L�[��----------------------------------------
  do j=1,JM-1
   do i=0,IM-1
    kx=0.5d0*(xix(i,j)+xix(i+1,j)) !���ԓ_��k_x=��_x
    ky=0.5d0*(xiy(i,j)+xiy(i+1,j)) !���ԓ_��k_y=��_y
    JJm=0.5d0*(JJ(i,j)+JJ(i+1,j)) !���ԓ_��JJ
    um=(dsqrt(rho(i,j))*u(i,j)+dsqrt(rho(i+1,j))*u(i+1,j))/(dsqrt(rho(i,j))+dsqrt(rho(i+1,j))) !���ԓ_��u
    vm=(dsqrt(rho(i,j))*v(i,j)+dsqrt(rho(i+1,j))*v(i+1,j))/(dsqrt(rho(i,j))+dsqrt(rho(i+1,j))) !���ԓ_��v
    dum=(dsqrt(rho(i,j))*Q(4,i,j)/Q(1,i,j)+dsqrt(rho(i+1,j))*Q(4,i+1,j)/Q(1,i+1,j))&
        /(dsqrt(rho(i,j))+dsqrt(rho(i+1,j))) !���ԓ_��e/��
    cm=sqrt(g*(g-1.0d0)*dabs(dum-0.5d0*(um**2.0+vm**2.0))) !���ԓ_��c�i�����j

    kxt=kx/dsqrt(kx**2.0+ky**2.0) !k_x~
    kyt=ky/dsqrt(kx**2.0+ky**2.0) !k_y~
    phi=0.5d0*(g-1.0d0)*(um**2.0+vm**2.0) !��
    thetat=kxt*um+kyt*vm !��~
    beta=1.0d0/(2.0d0*cm**2.0)

    !�s��R�i���ԓ_�j
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

    !R�̋t�s��R^-1�i���ԓ_�j
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

    !�ŗL�la
    a(1,i)=kx*um+ky*vm
    a(2,i)=a(1,i)
    a(3,i)=a(1,i)+cm*dsqrt(kx**2.0+ky**2.0)
    a(4,i)=a(1,i)-cm*dsqrt(kx**2.0+ky**2.0)
    a(5,i)=a(1,i)
    a(6,i)=a(1,i)

    !��
    do l=1,6
     D(l)=(Q(l,i+1,j)*JJ(i+1,j)-Q(l,i,j)*JJ(i,j))/JJm
    end do
    do l=1,6
     if(l<=4) then
      alpha(l,i)=RI(l,1,i)*D(1)+RI(l,2,i)*D(2)+RI(l,3,i)*D(3)+RI(l,4,i)*D(4)
     else
      alpha(l,i)=D(l)
     end if
    end do
   end do
   do l=1,6
    alpha(l,-1)=alpha(l,0)
    alpha(l,IM)=alpha(l,IM-1)
   end do

  !�����֐�g
   do l=1,6
    do i=0,IM
     S=dsign(1.0d0,alpha(l,i))
     gg(l,i)=S*dmax1(0.0d0,dmin1(dabs(alpha(l,i)),S*alpha(l,i-1)))
    end do
   end do

   do l=1,6
    do i=0,IM-1
     !���i��M��Ƃ͕ʁj
     delta=1.0d-10
     if(alpha(l,i)/=0.0d0)then
      gamma=0.5d0*FPSI(a(l,i),delta)*(gg(l,i+1)-gg(l,i))/alpha(l,i)
     else
      gamma=0.0d0
     end if
     !���i�ӂƂ͕ʁj
     phim(l,i)=0.5d0*FPSI(a(l,i),delta)*(gg(l,i)+gg(l,i+1))-FPSI(a(l,i)+gamma,delta)*alpha(l,i)
    end do
   end do

   !���l����E~
   do l=1,6
    do i=0,IM-1
     if(l<=4) then
      EE(l,i,j)=0.5d0*(E(l,i,j)+E(l,i+1,j)&
                +RR(l,1,i)*phim(1,i)+RR(l,2,i)*phim(2,i)&
                +RR(l,3,i)*phim(3,i)+RR(l,4,i)*phim(4,i))
     else
      EE(l,i,j)=0.5d0*(E(l,i,j)+E(l,i+1,j)+phim(l,i))
     end if
    end do
   end do
  end do

  !�ŕ���TVD�X�L�[��-----------------------------------------
  do i=1,IM-1
   do j=0,JM-1
    kx=0.5d0*(etax(i,j)+etax(i,j+1)) !���ԓ_��k_x=��_x
    ky=0.5d0*(etay(i,j)+etay(i,j+1)) !���ԓ_��k_y=��_y
    JJm=0.5d0*(JJ(i,j)+JJ(i,j+1)) !���ԓ_��JJ
    um=(dsqrt(rho(i,j))*u(i,j)+dsqrt(rho(i,j+1))*u(i,j+1))/(dsqrt(rho(i,j))+dsqrt(rho(i,j+1))) !���ԓ_��u
    vm=(dsqrt(rho(i,j))*v(i,j)+dsqrt(rho(i,j+1))*v(i,j+1))/(dsqrt(rho(i,j))+dsqrt(rho(i,j+1))) !���ԓ_��v
    dum=(dsqrt(rho(i,j))*Q(4,i,j)/Q(1,i,j)+dsqrt(rho(i,j+1))*Q(4,i,j+1)/Q(1,i,j+1))&
        /(dsqrt(rho(i,j))+dsqrt(rho(i,j+1))) !���ԓ_��e/��
    cm=sqrt(g*(g-1.0d0)*dabs(dum-0.5d0*(um**2.0+vm**2.0))) !���ԓ_��c�i�����j

    kxt=kx/dsqrt(kx**2.0+ky**2.0) !k_x~
    kyt=ky/dsqrt(kx**2.0+ky**2.0) !k_y~
    phi=0.5d0*(g-1.0d0)*(um**2.0+vm**2.0) !��
    thetat=kxt*um+kyt*vm !��~
    beta=1.0d0/(2.0d0*cm**2.0)

    !�s��R�i���ԓ_�j
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

    !R�̋t�s��R^-1�i���ԓ_�j
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

    !�ŗL�la
    a(1,j)=kx*um+ky*vm
    a(2,j)=a(1,j)
    a(3,j)=a(1,j)+cm*dsqrt(kx**2.0+ky**2.0)
    a(4,j)=a(1,j)-cm*dsqrt(kx**2.0+ky**2.0)
    a(5,j)=a(1,j)
    a(6,j)=a(1,j)

    !��
    do l=1,6
     D(l)=(Q(l,i,j+1)*JJ(i,j+1)-Q(l,i,j)*JJ(i,j))/JJm
    end do
    do l=1,6
     if(l<=4) then
      alpha(l,j)=RI(l,1,j)*D(1)+RI(l,2,j)*D(2)+RI(l,3,j)*D(3)+RI(l,4,j)*D(4)
     else
      alpha(l,j)=D(l)
     end if
    end do
   end do
   do l=1,6
    alpha(l,-1)=alpha(l,0)
    alpha(l,JM)=alpha(l,JM-1)
   end do

   !�����֐�g
   do l=1,6
    do j=0,JM
     S=dsign(1.0d0,alpha(l,j))
     gg(l,j)=S*dmax1(0.0d0,dmin1(dabs(alpha(l,j)),S*alpha(l,j-1)))
    end do
   end do

   do l=1,6
    do j=0,JM-1
     !���i��M��Ƃ͕ʁj
     delta=1.0d-10
     if(alpha(l,j)/=0.0d0)then
      gamma=0.5d0*FPSI(a(l,j),delta)*(gg(l,j+1)-gg(l,j))/alpha(l,j)
     else
      gamma=0.0d0
     end if
     !���i�ӂƂ͕ʁj
     phim(l,j)=0.5d0*FPSI(a(l,j),delta)*(gg(l,j)+gg(l,j+1))-FPSI(a(l,j)+gamma,delta)*alpha(l,j)
    end do
   end do

   !���l����F~
   do l=1,6
    do j=0,JM-1
     if(l<=4) then
      FF(l,i,j)=0.5d0*(F(l,i,j)+F(l,i,j+1)&
                +RR(l,1,j)*phim(1,j)+RR(l,2,j)*phim(2,j)&
                +RR(l,3,j)*phim(3,j)+RR(l,4,j)*phim(4,j))
     else
      FF(l,i,j)=0.5d0*(F(l,i,j)+F(l,i,j+1)+phim(l,j))
     end if
    end do
   end do
  end do

  !�S�����v�Z-----------------------------------------------
  !�̕����i�X�^�K�[�h�i�q��R^�����߂�j
  do j=1,JM-1 !S�t���ϐ��̓X�^�K�[�h�i�q�_��̒l�i�Y����i+1/2,j�j
   do i=0,IM-1
    xixS(i,j)=(xix(i,j)+xix(i+1,j))/2.0d0
    xiyS(i,j)=(xiy(i,j)+xiy(i+1,j))/2.0d0
    etaxS(i,j)=(etax(i,j)+etax(i+1,j))/2.0d0
    etayS(i,j)=(etay(i,j)+etay(i+1,j))/2.0d0
    JJS(i,j)=(JJ(i,j)+JJ(i+1,j))/2.0d0
    rhoS(i,j)=(rho(i,j)+rho(i+1,j))/2.0d0
    uS(i,j)=(u(i,j)+u(i+1,j))/2.0d0
    vS(i,j)=(v(i,j)+v(i+1,j))/2.0d0
    TS(i,j)=(T(i,j)+T(i+1,j))/2.0d0
    kS(i,j)=(k(i,j)+k(i+1,j))/2.0d0
    epS(i,j)=(ep(i,j)+ep(i+1,j))/2.0d0
    mutS(i,j)=(Cmu*rhoS(i,j)*kS(i,j)**2.0)/epS(i,j)
    uxiS(i,j)=u(i+1,j)-u(i,j)
    vxiS(i,j)=v(i+1,j)-v(i,j)
    TxiS(i,j)=T(i+1,j)-T(i,j)
    kxiS(i,j)=k(i+1,j)-k(i,j)
    epxiS(i,j)=ep(i+1,j)-ep(i,j)
    uetaS(i,j)=(u(i,j+1)+u(i+1,j+1)-u(i,j-1)-u(i+1,j-1))/4.0d0
    vetaS(i,j)=(v(i,j+1)+v(i+1,j+1)-v(i,j-1)-v(i+1,j-1))/4.0d0
    TetaS(i,j)=(T(i,j+1)+T(i+1,j+1)-T(i,j-1)-T(i+1,j-1))/4.0d0
    ketaS(i,j)=(k(i,j+1)+k(i+1,j+1)-k(i,j-1)-k(i+1,j-1))/4.0d0
    epetaS(i,j)=(ep(i,j+1)+ep(i+1,j+1)-ep(i,j-1)-ep(i+1,j-1))/4.0d0
    muS(i,j)=((TS(i,j)/T0)**1.5)*(T0+Sl)/(TS(i,j)+Sl)*mu0 !�S���W���i�T�U�[�����h�̎��j
    !�S�����x�N�g��
    SxxS=2.0d0*(2.0d0*(uxiS(i,j)*xixS(i,j)+uetaS(i,j)*etaxS(i,j))-(vxiS(i,j)*xiyS(i,j)+vetaS(i,j)*etayS(i,j)))/3.0d0
    SxyS=(uxiS(i,j)*xiyS(i,j)+uetaS(i,j)*etayS(i,j))+(vxiS(i,j)*xixS(i,j)+vetaS(i,j)*etaxS(i,j))
    SyyS=2.0d0*(2.0d0*(vxiS(i,j)*xiyS(i,j)+vetaS(i,j)*etayS(i,j))-(uxiS(i,j)*xixS(i,j)+uetaS(i,j)*etaxS(i,j)))/3.0d0
    txxS=mutS(i,j)*SxxS-2.0d0*rhoS(i,j)*kS(i,j)/3.0d0
    txyS=mutS(i,j)*SxyS
    tyyS=mutS(i,j)*SyyS-2.0d0*rhoS(i,j)*kS(i,j)/3.0d0
    RvS(2,i,j)=txxS
    RvS(3,i,j)=txyS
    RvS(4,i,j)=txxS*uS(i,j)+txyS*vS(i,j)+mutS(i,j)*g*R*(TxiS(i,j)*xixS(i,j)+TetaS(i,j)*etaxS(i,j))/(Prt*(g-1.0d0))
    RvS(5,i,j)=mutS(i,j)*(kxiS(i,j)*xixS(i,j)+ketaS(i,j)*etaxS(i,j))/sigk
    RvS(6,i,j)=mutS(i,j)*(epxiS(i,j)*xixS(i,j)+epetaS(i,j)*etaxS(i,j))/sigep
    SvS(2,i,j)=txyS
    SvS(3,i,j)=tyyS
    SvS(4,i,j)=txyS*uS(i,j)+tyyS*vS(i,j)+mutS(i,j)*g*R*(TxiS(i,j)*xiyS(i,j)+TetaS(i,j)*etayS(i,j))/(Prt*(g-1.0d0))
    SvS(5,i,j)=mutS(i,j)*(kxiS(i,j)*xiyS(i,j)+ketaS(i,j)*etayS(i,j))/sigk
    SvS(6,i,j)=mutS(i,j)*(epxiS(i,j)*xiyS(i,j)+epetaS(i,j)*etayS(i,j))/sigep
    do l=1,6 !��������ʍ��W�n�ɕϊ�
     RvSS(l,i,j)=(xixS(i,j)*RvS(l,i,j)+xiyS(i,j)*SvS(l,i,j))/JJS(i,j)
    end do
   end do
  end do

  !�ŕ����i�X�^�K�[�h�i�q��S^�����߂�j
  do j=0,JM-1 !S�t���ϐ��̓X�^�K�[�h�i�q�_��̒l�i�Y����i,j+1/2�j
   do i=1,IM-1
    xixS(i,j)=(xix(i,j)+xix(i,j+1))/2.0d0
    xiyS(i,j)=(xiy(i,j)+xiy(i,j+1))/2.0d0
    etaxS(i,j)=(etax(i,j)+etax(i,j+1))/2.0d0
    etayS(i,j)=(etay(i,j)+etay(i,j+1))/2.0d0
    JJS(i,j)=(JJ(i,j)+JJ(i,j+1))/2.0d0
    rhoS(i,j)=(rho(i,j)+rho(i,j+1))/2.0d0
    uS(i,j)=(u(i,j)+u(i,j+1))/2.0d0
    vS(i,j)=(v(i,j)+v(i,j+1))/2.0d0
    TS(i,j)=(T(i,j)+T(i,j+1))/2.0d0
    kS(i,j)=(k(i,j)+k(i,j+1))/2.0d0
    epS(i,j)=(ep(i,j)+ep(i,j+1))/2.0d0
    mutS(i,j)=(Cmu*rhoS(i,j)*kS(i,j)**2.0)/epS(i,j)
    uxiS(i,j)=(u(i+1,j)+u(i+1,j+1)-u(i-1,j)-u(i-1,j+1))/4.0d0
    vxiS(i,j)=(v(i+1,j)+v(i+1,j+1)-v(i-1,j)-v(i-1,j+1))/4.0d0
    TxiS(i,j)=(T(i+1,j)+T(i+1,j+1)-T(i-1,j)-T(i-1,j+1))/4.0d0
    kxiS(i,j)=(k(i+1,j)+k(i+1,j+1)-k(i-1,j)-k(i-1,j+1))/4.0d0
    epxiS(i,j)=(ep(i+1,j)+ep(i+1,j+1)-ep(i-1,j)-ep(i-1,j+1))/4.0d0
    uetaS(i,j)=u(i,j+1)-u(i,j)
    vetaS(i,j)=v(i,j+1)-v(i,j)
    TetaS(i,j)=T(i,j+1)-T(i,j)
    ketaS(i,j)=k(i,j+1)-k(i,j)
    epetaS(i,j)=ep(i,j+1)-ep(i,j)
    muS(i,j)=((TS(i,j)/T0)**1.5)*(T0+Sl)/(TS(i,j)+Sl)*mu0 !�S���W���i�T�U�[�����h�̎��j
    !�S�����x�N�g��
    SxxS=2.0d0*(2.0d0*(uxiS(i,j)*xixS(i,j)+uetaS(i,j)*etaxS(i,j))-(vxiS(i,j)*xiyS(i,j)+vetaS(i,j)*etayS(i,j)))/3.0d0
    SxyS=(uxiS(i,j)*xiyS(i,j)+uetaS(i,j)*etayS(i,j))+(vxiS(i,j)*xixS(i,j)+vetaS(i,j)*etaxS(i,j))
    SyyS=2.0d0*(2.0d0*(vxiS(i,j)*xiyS(i,j)+vetaS(i,j)*etayS(i,j))-(uxiS(i,j)*xixS(i,j)+uetaS(i,j)*etaxS(i,j)))/3.0d0
    txxS=mutS(i,j)*SxxS-2.0d0*rhoS(i,j)*kS(i,j)/3.0d0
    txyS=mutS(i,j)*SxyS
    tyyS=mutS(i,j)*SyyS-2.0d0*rhoS(i,j)*kS(i,j)/3.0d0
    RvS(2,i,j)=txxS
    RvS(3,i,j)=txyS
    RvS(4,i,j)=txxS*uS(i,j)+txyS*vS(i,j)+mutS(i,j)*g*R*(TxiS(i,j)*xixS(i,j)+TetaS(i,j)*etaxS(i,j))/(Prt*(g-1.0d0))
    RvS(5,i,j)=mutS(i,j)*(kxiS(i,j)*xixS(i,j)+ketaS(i,j)*etaxS(i,j))/sigk
    RvS(6,i,j)=mutS(i,j)*(epxiS(i,j)*xixS(i,j)+epetaS(i,j)*etaxS(i,j))/sigep
    SvS(2,i,j)=txyS
    SvS(3,i,j)=tyyS
    SvS(4,i,j)=txyS*uS(i,j)+tyyS*vS(i,j)+mutS(i,j)*g*R*(TxiS(i,j)*xiyS(i,j)+TetaS(i,j)*etayS(i,j))/(Prt*(g-1.0d0))
    SvS(5,i,j)=mutS(i,j)*(kxiS(i,j)*xiyS(i,j)+ketaS(i,j)*etayS(i,j))/sigk
    SvS(6,i,j)=mutS(i,j)*(epxiS(i,j)*xiyS(i,j)+epetaS(i,j)*etayS(i,j))/sigep
    do l=1,6 !�e��������ʍ��W�n�ɕϊ�
     SvSS(l,i,j)=(etaxS(i,j)*RvS(l,i,j)+etayS(i,j)*SvS(l,i,j))/JJS(i,j)
    end do
   end do
  end do

  !�������v�Z-----------------------------------------------
  do j=1,JM-1
   do i=1,IM-1
    !mut(i,j)=(Cmu*rho(i,j)*k(i,j)**2.0)/ep(i,j)
    uxi(i,j)=(u(i+1,j)-u(i-1,j))/2.0d0
    vxi(i,j)=(v(i+1,j)-v(i-1,j))/2.0d0
    ueta(i,j)=(u(i,j+1)-u(i,j-1))/2.0d0
    veta(i,j)=(v(i,j+1)-v(i,j-1))/2.0d0
    !�������x�N�g��
    Sxx=2.0d0*(2.0d0*(uxi(i,j)*xix(i,j)+ueta(i,j)*etax(i,j))-(vxi(i,j)*xiy(i,j)+veta(i,j)*etay(i,j)))/3.0d0
    Sxy=(uxi(i,j)*xiy(i,j)+ueta(i,j)*etay(i,j))+(vxi(i,j)*xix(i,j)+veta(i,j)*etax(i,j))
    Syy=2.0d0*(2.0d0*(vxi(i,j)*xiy(i,j)+veta(i,j)*etay(i,j))-(uxi(i,j)*xix(i,j)+ueta(i,j)*etax(i,j)))/3.0d0
    txx=mut(i,j)*Sxx-2.0d0*rho(i,j)*k(i,j)/3.0d0
    txy=mut(i,j)*Sxy
    tyy=mut(i,j)*Syy-2.0d0*rho(i,j)*k(i,j)/3.0d0
    PRO=txx*(uxi(i,j)*xix(i,j)+ueta(i,j)*etax(i,j))&
        +txy*(vxi(i,j)*xix(i,j)+veta(i,j)*etax(i,j)+uxi(i,j)*xiy(i,j)+ueta(i,j)*etay(i,j))&
        +tyy*(vxi(i,j)*xiy(i,j)+veta(i,j)*etay(i,j))
    PP(5,i,j)=(PRO-rho(i,j)*ep(i,j))/JJ(i,j)
    PP(6,i,j)=(Cep1*ep(i,j)*PRO/k(i,j)-(Cep2*rho(i,j)*ep(i,j)**2.0)/k(i,j))/JJ(i,j)
   end do
  end do

  !�ǖ@��---------------------------------------------------
  do i=10,IM
   nu(i,1)=mu(i,1)/rho(i,1)
   kp=k(i,1)
   epp=ep(i,1)
   yp=y(i,1)
   ut=(Cmu*kp**2.0)/(kap*yp*epp)
   ypp=yp*ut/nu(i,1)
   if(ypp>11.63d0) then
    ut=kap*u(i,1)/dlog(ypp*Ec)
   else
    ut=dsqrt(nu(i,1)*u(i,1)/yp)
   end if
   k(i,1)=ut**2.0/dsqrt(Cmu)
   ep(i,1)=ut**3.0/(kap*yp)
  end do

  !Runge-Kutta�@--------------------------------------------
  do j=1,JM-1
   do i=1,IM-1
    do l=1,6
     QQ(l,i,j,nn)=QQ(l,i,j,0)-dt(i,j)*(EE(l,i,j)-EE(l,i-1,j)+FF(l,i,j)-FF(l,i,j-1)&
                  -(RvSS(l,i,j)-RvSS(l,i-1,j)+SvSS(l,i,j)-SvSS(l,i,j-1)+PP(l,i,j)))/dble(5-nn)
     Q(l,i,j)=QQ(l,i,j,nn)
    end do
   end do
  end do

  !���ʂ̌v�Z-----------------------------------------------
  !�ۑ��ʃx�N�g��Q^���珔�ʂ����߂�
  do i=1,IM-1
   do j=1,JM-1
    rho(i,j)=Q(1,i,j)*JJ(i,j)
    u(i,j)=Q(2,i,j)*JJ(i,j)/rho(i,j)
    v(i,j)=Q(3,i,j)*JJ(i,j)/rho(i,j)
    ene(i,j)=Q(4,i,j)*JJ(i,j)
    p(i,j)=(g-1.0d0)*(ene(i,j)-(u(i,j)**2.0+v(i,j)**2.0)*rho(i,j)/2.0d0)
    T(i,j)=p(i,j)/(rho(i,j)*R)
    if(i<10.or.j/=1) then
     k(i,j)=Q(5,i,j)*JJ(i,j)/rho(i,j)
     ep(i,j)=Q(6,i,j)*JJ(i,j)/rho(i,j)
    end if
   end do
  end do
  do j=0,JM
   do i=0,IM
    !���E�����̍X�V
    if(i==IM)then !���o���ꎟ�O�}
     rho(IM,j)=rho(IM-1,j)+(rho(IM-1,j)-rho(IM-2,j))*(x(IM,j)-x(IM-1,j))/(x(IM-1,j)-x(IM-2,j))
     u(IM,j)=u(IM-1,j)+(u(IM-1,j)-u(IM-2,j))*(x(IM,j)-x(IM-1,j))/(x(IM-1,j)-x(IM-2,j))
     v(IM,j)=v(IM-1,j)+(v(IM-1,j)-v(IM-2,j))*(x(IM,j)-x(IM-1,j))/(x(IM-1,j)-x(IM-2,j))
     p(IM,j)=p(IM-1,j)+(p(IM-1,j)-p(IM-2,j))*(x(IM,j)-x(IM-1,j))/(x(IM-1,j)-x(IM-2,j))
     T(IM,j)=p(IM,j)/(rho(IM,j)*R)
     ene(IM,j)=rho(IM,j)*(R*T(IM,j)/(g-1.0d0)+(u(IM,j)**2.0+v(IM,j)**2.0)/2.0d0)
     k(IM,j)=k(IM-1,j)+(k(IM-1,j)-k(IM-2,j))*(x(IM,j)-x(IM-1,j))/(x(IM-1,j)-x(IM-2,j))
     ep(IM,j)=ep(IM-1,j)+(ep(IM-1,j)-ep(IM-2,j))*(x(IM,j)-x(IM-1,j))/(x(IM-1,j)-x(IM-2,j))
    end if
    if(j==JM)then !��[���ꎟ�O�}
     rho(i,JM)=rho(i,JM-1) !+(rho(i,JM-1)-rho(i,JM-2))*(y(i,JM)-y(i,JM-1))/(y(i,JM-1)-y(i,JM-2))
     u(i,JM)=u(i,JM-1) !+(u(i,JM-1)-u(i,JM-2))*(y(i,JM)-y(i,JM-1))/(y(i,JM-1)-y(i,JM-2))
     v(i,JM)=v(i,JM-1) !+(v(i,JM-1)-v(i,JM-2))*(y(i,JM)-y(i,JM-1))/(y(i,JM-1)-y(i,JM-2))
     p(i,JM)=p(i,JM-1) !+(p(i,JM-1)-p(i,JM-2))*(y(i,JM)-y(i,JM-1))/(y(i,JM-1)-y(i,JM-2))
     T(i,JM)=p(i,JM)/(rho(i,JM)*R)
     ene(i,JM)=rho(i,JM)*(R*T(i,JM)/(g-1.0d0)+(u(i,JM)**2.0+v(i,JM)**2.0)/2.0d0)
     k(i,JM)=k(i,JM-1) !+(k(i,JM-1)-k(i,JM-2))*(y(i,JM)-y(i,JM-1))/(y(i,JM-1)-y(i,JM-2))
     ep(i,JM)=ep(i,JM-1) !+(ep(i,JM-1)-ep(i,JM-2))*(y(i,JM)-y(i,JM-1))/(y(i,JM-1)-y(i,JM-2))
    end if
    if(j==0)then !���[���ǂŊ���Ȃ��C���͈ꎟ�O�}
     rho(i,0)=rho(i,1) !+(rho(2,j)-rho(1,j))*(y(i,1)-y(i,0))/(y(i,2)-y(i,1))
     if(i<10) then
      u(i,0)=u(i,1) !+(u(2,j)-u(1,j))*(y(i,1)-y(i,0))/(y(i,2)-y(i,1))
     else
      u(i,0)=0.0d0
     end if
     v(i,0)=0.0d0
     p(i,0)=p(i,1) !+(p(2,j)-p(1,j))*(y(i,1)-y(i,0))/(y(i,2)-y(i,1))
     T(i,0)=p(i,0)/(rho(i,0)*R)
     ene(i,0)=rho(i,0)*(R*T(i,0)/(g-1.0d0)+(u(i,0)**2.0+v(i,0)**2.0)/2.0d0)
     k(i,0)=k(i,1) !+(k(2,j)-k(1,j))*(y(i,1)-y(i,0))/(y(i,2)-y(i,1))
     ep(i,0)=ep(i,1) !+(ep(2,j)-ep(1,j))*(y(i,1)-y(i,0))/(y(i,2)-y(i,1))
    end if
    !�����x�N�g���̍X�V
    Q(1,i,j)=rho(i,j)
    Q(2,i,j)=rho(i,j)*u(i,j)
    Q(3,i,j)=rho(i,j)*v(i,j)
    Q(4,i,j)=ene(i,j)
    Q(5,i,j)=rho(i,j)*k(i,j)
    Q(6,i,j)=rho(i,j)*ep(i,j)
    E(1,i,j)=rho(i,j)*u(i,j)
    E(2,i,j)=p(i,j)+rho(i,j)*u(i,j)**2.0
    E(3,i,j)=rho(i,j)*u(i,j)*v(i,j)
    E(4,i,j)=(ene(i,j)+p(i,j))*u(i,j)
    E(5,i,j)=rho(i,j)*k(i,j)*u(i,j)
    E(6,i,j)=rho(i,j)*ep(i,j)*u(i,j)
    F(1,i,j)=rho(i,j)*v(i,j)
    F(2,i,j)=rho(i,j)*u(i,j)*v(i,j)
    F(3,i,j)=p(i,j)+rho(i,j)*v(i,j)**2.0
    F(4,i,j)=(ene(i,j)+p(i,j))*v(i,j)
    F(5,i,j)=rho(i,j)*k(i,j)*v(i,j)
    F(6,i,j)=rho(i,j)*ep(i,j)*v(i,j)

    do l=1,6 !�e��������ʍ��W�n�ɕϊ�
     Q(l,i,j)=Q(l,i,j)/JJ(i,j)
     Eh=E(l,i,j)
     E(l,i,j)=(xix(i,j)*E(l,i,j)+xiy(i,j)*F(l,i,j))/JJ(i,j)
     F(l,i,j)=(etax(i,j)*Eh+etay(i,j)*F(l,i,j))/JJ(i,j)
    end do

   end do
  end do

 end do

 !���̃X�e�b�v�̏����l��ݒ�/�ő�c���̌v�Z
 do i=0,IM
  do j=0,JM
   do l=1,6
    QQ(l,i,j,0)=Q(l,i,j)
   end do
   ddq=abs((Q(1,i,j)+Q(2,i,j)+Q(3,i,j)+Q(4,i,j)+Q(5,i,j)+Q(6,i,j)-qp(i,j))/qp(i,j))
   qp(i,j)=Q(1,i,j)+Q(2,i,j)+Q(3,i,j)+Q(4,i,j)+Q(5,i,j)+Q(6,i,j)
   if(ddq>dq) dq=ddq
  end do
 end do
 print*,'count =',n,'dq =',dq

end do


!y+,u+�̏����o��
do i=10,IM
 utau=dsqrt((mu(i,0)*(u(i,1)-u(i,0))/(y(i,1)-y(i,0)))/rho(i,0))
 do j=0,JM
  nu(i,j)=mu(i,j)/rho(i,j)
  uplus(i,j)=u(i,j)/utau
  yplus(i,j)=utau*y(i,j)/nu(i,j)
 end do
end do
open(1,file='comp8.csv',status='replace')
do j=0,JM
 write(1,*) yplus(45,j),uplus(45,j)
end do
close(1)


!.dat�t�@�C�������o��***************************************
open(2,file='comp8.dat',status='replace')
do j=0,JM
 do i=0,IM
  write(2,*) x(i,j),y(j,j),u(i,j),v(i,j),p(i,j),T(i,j) !'(F6.3,F6.3)'
 end do
end do
close(2)


!.fld�t�@�C���쐬*******************************************
open(3,file='comp8.fld',status='replace')
write(3,'(A)') '# AVS field file'
write(3,'(A)') 'ndim = 2'
write(3,'(A)',advance='no') 'dim1 ='
write(3,'(I5)') IM+1
write(3,'(A)',advance='no') 'dim2 ='
write(3,'(I5)') JM+1
write(3,'(A)') 'nspace = 2'
write(3,'(A)') 'veclen = 4'
write(3,'(A)') 'data = double'
write(3,'(A)') 'field = irregular'
write(3,'(A)') 'label = u v p T'
write(3,'(A)') 'variable 1 file=comp8.dat filetype=ascii skip=0 offset=2 stride=6'
write(3,'(A)') 'variable 2 file=comp8.dat filetype=ascii skip=0 offset=3 stride=6'
write(3,'(A)') 'variable 3 file=comp8.dat filetype=ascii skip=0 offset=4 stride=6'
write(3,'(A)') 'variable 4 file=comp8.dat filetype=ascii skip=0 offset=5 stride=6'
write(3,'(A)') 'coord 1 file=comp8.dat filetype=ascii skip=0 offset=0 stride=6'
write(3,'(A)') 'coord 2 file=comp8.dat filetype=ascii skip=0 offset=1 stride=6'
close(3)


contains
!�֐��̒�`
!FPSI�֐�---------------------------------------------------
double precision function FPSI(z,delta)
double precision z,delta
if(abs(z)>=delta)then
 FPSI=dabs(z)
else
 FPSI=0.5d0*(z**2.0+delta**2.0)/delta
end if
end function


end program comp8