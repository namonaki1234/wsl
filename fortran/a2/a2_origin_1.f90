!圧縮性課題2
program cka21
  implicit none

  !格子条件
  double precision,parameter::H=10.0d0		!代表長さ
  double precision,parameter::XM=3.0d0*H		!x方向長さ
  double precision,parameter::YM=H		!y方向長さ
  integer,parameter::IM=90			!x方向終点
  integer,parameter::JM=30			!y方向終点
  integer,parameter::IL=0				!x方向終点
  integer,parameter::JL=0				!y方向終点
  double precision x				!x座標
  double precision y				!y座標
  double precision,parameter::dx=XM/IM		!x方向格子刻み
  double precision,parameter::dy=YM/JM		!y方向格子刻み

  !流れ場初期条件
  double precision,parameter::R=287.1d0		!ガス定数
  double precision,parameter::G=1.4d0		!比熱比
  double precision,parameter::dt=1.0E-6		!時間刻み

  !配列宣言
  double precision u(IL-1:IM+1,JL-1:JM+1)		!x方向速度
  double precision v(IL-1:IM+1,JL-1:JM+1)		!y方向速度
  double precision p(IL-1:IM+1,JL-1:JM+1)		!圧力
  double precision T(IL-1:IM+1,JL-1:JM+1)		!絶対温度
  double precision rho(IL-1:IM+1,JL-1:JM+1)	!密度
  double precision M(IL-1:IM+1,JL-1:JM+1)		!マッハ数
  double precision en(IL-1:IM+1,JL-1:JM+1)	!全エネルギー
  double precision ENE(IL-1:IM+1,JL-1:JM+1)	!単位質量あたりの内部エネルギー
  double precision Q(IL-1:IM+1,JL-1:JM+1,1:4)	!保存量ベクトル
  double precision Q1(IL-1:IM+1,JL-1:JM+1,1:4)	!次の時刻の保存量ベクトル
  double precision E(IL-1:IM+1,JL-1:JM+1,1:4)	!x方向流束ベクトル
  double precision E1(IL-1:IM+1,JL-1:JM+1,1:4)	!次の時刻のx方向流束ベクトル
  double precision F(IL-1:IM+1,JL-1:JM+1,1:4)	!y方向流束ベクトル
  double precision F1(IL-1:IM+1,JL-1:JM+1,1:4)	!次の時刻のy方向流束ベクトル

  !計算条件
  integer::n=0					!カウンタ
  integer,parameter::NMAX=10			!繰り返し数
  integer i,j,k

  !初期条件計算
  do j=JL-1,JM+1
      do i=IL-1,IM+1
      M(i,j)=2.9
      T(i,j)=293.0
      RHO(i,j)=1.2
      p(i,j)=RHO(i,j)*R*T(i,j)
      u(i,j)=M(i,j)*sqrt(G*R*T(i,j))
      v(i,j)=0.0d0
      ENE(i,j)=(R*T(i,j))/(G-1)
      en(i,j)=RHO(i,j)*(ENE(i,j)+u(i,j)**2.0/2.0)

      Q(i,j,1)=RHO(i,j)
      Q(i,j,2)=RHO(i,j)*u(i,j)
      Q(i,j,3)=RHO(i,j)*v(i,j)
      Q(i,j,4)=en(i,j)

      E(i,j,1)=RHO(i,j)*u(i,j)
      E(i,j,2)=p(i,j)+RHO(i,j)*u(i,j)**2.0
      E(i,j,3)=RHO(i,j)*u(i,j)*v(i,j)
      E(i,j,4)=(en(i,j)+p(i,j))*u(i,j)

      F(i,j,1)=RHO(i,j)*v(i,j)
      F(i,j,2)=RHO(i,j)*v(i,j)*u(i,j)
      F(i,j,3)=p(i,j)+RHO(i,j)*v(i,j)**2.0
      F(i,j,4)=(en(i,j)+p(i,j))*v(i,j)
      end do
    end do

  !call outputX
  !call outputY

  !メイン計算
  do while(n<NMAX)
  n=n+1
  call Boundary
  call CalQ1(IL+1,IM-1,JL+1,JM-1)
  call CalQ2(IL+1,IM-1,JL+1,JM-1)
  print*,'繰返し数',n
  end do

  call output

  contains

  subroutine Boundary
  implicit none
  integer::i,j

  !上壁・下壁条件(1次外挿)
  do i=IL,IM
  u(i,JM+1)=(u(i,JM)-u(i,JM-1))/dy*dy+u(i,JM)
  u(i,JL-1)=(u(i,JL+1)-u(i,JL))/dy*dy+u(i,JL)
  v(i,JM+1)=(v(i,JM)-v(i,JM-1))/dy*dy+v(i,JM)
  v(i,JL-1)=(v(i,JL+1)-v(i,JL))/dy*dy+v(i,JL)
  rho(i,JM+1)=(rho(i,JM)-rho(i,JM-1))/dy*dy+rho(i,JM)
  rho(i,JL-1)=(rho(i,JL+1)-rho(i,JL))/dy*dy+rho(i,JL)
  en(i,JM+1)=(en(i,JM)-en(i,JM-1))/dy*dy+en(i,JM)
  en(i,JL-1)=(en(i,JL+1)-en(i,JL))/dy*dy+en(i,JL)
  end do

  !流入条件
  do j=JL,JM
  M(IL,j)=2.9d0
  T(IL,j)=293.0d0
  rho(IL,j)=1.2d0
  p(IL,j)=rho(IL,j)*R*T(IL,j)
  u(IL,j)=M(IL,j)*sqrt(G*R*T(IL,j))
  v(IL,j)=0.0d0
  ENE(IL,j)=(R*T(IL,j))/(G-1.0d0)
  en(IL,j)=rho(IL,j)*(ENE(IL,j)+u(IL,j)**2.0d0/2.0d0)

  !流出条件(1次外挿)
  u(IM+1,j)=(u(IM,j)-u(IM-1,j))/dx*dx+u(IM,j)
  v(IM+1,j)=(v(IM,j)-v(IM-1,j))/dx*dx+v(IM,j)
  rho(IM+1,j)=(rho(IM,j)-rho(IM-1,j))/dx*dx+rho(IM,j)
  en(IM+1,j)=(en(IM,j)-en(IM-1,j))/dx*dx+en(IM,j)
  end do
  end subroutine Boundary

  !x,y方向流束ベクトルE,Fの計算
  subroutine CalQ1(is,ie,js,je)
  implicit none
  integer is,ie,js,je
  do j=js,je
    do i=is,ie
    E1(i,j,1)=rho(i,j)*u(i,j)
    E1(i,j,2)=p(i,j)+rho(i,j)*(u(i,j)**2.0d0)
    E1(i,j,3)=rho(i,j)*u(i,j)*v(i,j)
    E1(i,j,4)=(en(i,j)+p(i,j))*u(i,j)

    F1(i,j,1)=rho(i,j)*v(i,j)
    F1(i,j,2)=rho(i,j)*u(i,j)*v(i,j)
    F1(i,j,3)=p(i,j)+rho(i,j)*(v(i,j)**2.0d0)
    F1(i,j,4)=(en(i,j)+p(i,j))*v(i,j)
      do k=1,4
      E(i,j,k)=E1(i,j,k)
      F(i,j,k)=F1(i,j,k)
      end do
    end do
  end do
  end subroutine CalQ1

  !保存量Qの計算
  subroutine CalQ2(is,ie,js,je)
  implicit none
  integer is,ie,js,je
  do j=js,je
    do i=is,ie
      do k=1,4
      Q1(i,j,k)=Q(i,j,k)-dt*((E(i+1,j,k)-E(i-1,j,k))/(2.0d0*dx)+(F(i,j+1,k)-F(i,j-1,k))/(2.0d0*dy))
      end do
    rho(i,j)=Q1(i,j,1)
    u(i,j)=Q1(i,j,2)/rho(i,j)
    v(i,j)=Q1(i,j,3)/rho(i,j)
    en(i,j)=Q1(i,j,4)
    ENE(i,j)=en(i,j)/rho(i,j)-((u(i,j)**2.0d0)/2.0d0)
    T(i,j)=ENE(i,j)*(G-1)/R
    p(i,j)=rho(i,j)*R*T(i,j)
    M(i,j)=u(i,j)/sqrt(G*R*T(i,j))
      do k=1,4
      Q(i,j,k)=Q1(i,j,k)
      end do
    end do
  end do
  end subroutine CalQ2

  subroutine output
  open(17,file='a.dat',status='replace')
  do j=JL,JM
  do i=IL,IM
  write(17,'(f13.6,$)')u(i,j)
  write(17,'(f13.6)')v(i,j)
  end do
  end do
  close(17)
  end subroutine output

  subroutine outputX
  open(18,file='cka21x.dat',status='replace')
  do i=IL,IM
  x=(XM/dble(IM))*dble(i)
  write(18,*)x
  end do
  close(18)
  end subroutine outputX

  subroutine outputY
  open(19,file='cka21y.dat',status='replace')
  do j=JL,JM
  y=(YM/dble(JM))*dble(j)
  write(19,*)y
  end do
  close(19)
  end subroutine outputY

  end program cka21

