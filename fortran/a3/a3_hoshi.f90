program kadai3
implicit none

!---------------格子・計算空間条件・刻み---------------
double precision,parameter :: H = 10.0d0 !格子数
double precision,parameter :: XM = 3.0d0*H !長さ_x方向
double precision,parameter :: YM = H !長さ_y方向

integer,parameter :: IM = 30 !終点_x方向
integer,parameter :: JM = 30 !終点_y方向
integer,parameter :: IL = 0 !始点_x方向
integer,parameter :: JL = 0 !始点_y方向

double precision,dimension(IL-1:IM+1,JL-1:JM+1) :: x = 0.0d0 !不等間隔_x
double precision,dimension(IL-1:IM+1,JL-1:JM+1) :: y = 0.0d0 !不等間隔_y
double precision jy !格子生成_べき乗
double precision,parameter :: a = 1.2d0 !格子生成_定数
double precision,parameter :: b = (a+1.0d0)/(a-1.0d0)

double precision,parameter :: R = 287.1d0 !ガス定数
double precision,parameter :: GAMMA = 1.4d0 !比熱比

double precision,parameter :: dt = 0.000001d0 !時間刻み

double precision,parameter :: dx = XM/dble(IM) !x方向微小長さ
double precision dy !y方向微小長さ
double precision,parameter :: dxi = 1.0d0 !ξ方向微小長さ
double precision,parameter :: deta = 1.0d0 !η方向微小長さ

integer n !計算回数
integer,parameter :: NMAX = 10 !計算回数_上限

!---------------配列定義---------------
double precision,dimension(IL:IM,JL:JM) :: u = 0.0d0 !速度_x方向
double precision,dimension(IL:IM,JL:JM) :: v = 0.0d0 !速度_y方向
double precision,dimension(IL:IM,JL:JM) :: p = 0.0d0 !圧力
double precision,dimension(IL:IM,JL:JM) :: T = 0.0d0 !絶対温度[k]
double precision,dimension(IL:IM,JL:JM) :: RHO = 0.0d0 !密度
double precision,dimension(IL:IM,JL:JM) :: M = 0.0d0 !マッハ数
double precision,dimension(IL:IM,JL:JM) :: e = 0.0d0 !全エネ
double precision,dimension(IL:IM,JL:JM) :: EB = 0.0d0 !エネ／単位質量

double precision,dimension(IL:IM,JL:JM) :: UL = 0.0d0 !ξ方向反変速度
double precision,dimension(IL:IM,JL:JM) :: VL = 0.0d0 !η方向反変速度

double precision,dimension(IL:IM,JL:JM,1:4) :: QH = 0.0d0 !保存量ベク
double precision,dimension(IL:IM,JL:JM,1:4) :: QH1 = 0.0d0 !QH_n-1
double precision,dimension(IL:IM,JL:JM,1:4) :: EH = 0.0d0 !x方向流速ベクトル
double precision,dimension(IL:IM,JL:JM,1:4) :: EH1 = 0.0d0 !EH_n-1
double precision,dimension(IL:IM,JL:JM,1:4) :: FH = 0.0d0 !y方向流速ベク
double precision,dimension(IL:IM,JL:JM,1:4) :: FH1 = 0.0d0 !FH_n-1

double precision,dimension(IL:IM,JL:JM) :: Ja = 0.0d0 !ヤコビアン
double precision,dimension(IL:IM,JL:JM) :: Yeta = 0.0d0 !dy/dη
double precision,dimension(IL:IM,JL:JM) :: Xeta = 0.0d0 !dx/dη
double precision,dimension(IL:IM,JL:JM) :: Yxi = 0.0d0 !dy/dξ
double precision,dimension(IL:IM,JL:JM) :: Xxi = 0.0d0 !dx/dξ

double precision,dimension(IL:IM,JL:JM) :: XIx = 0.0d0 !dξ/dx
double precision,dimension(IL:IM,JL:JM) :: XIy = 0.0d0 !dξ/dy
double precision,dimension(IL:IM,JL:JM) :: ETAx = 0.0d0 !dη/dx
double precision,dimension(IL:IM,JL:JM) :: ETAy = 0.0d0 !dη/dy

!---------------初期条件設定---------------
integer :: i,j !x,y方向座標
integer :: k !行列成分
do i = IL,IM
  do j = JL,JM
    M(i,j) = 2.9d0
    T(i,j) = 293.d0
    RHO(i,j) = 1.2d0
    p(i,j) = RHO(i,j)*R*T(i,j)
    u(i,j) = M(i,j)*sqrt(GAMMA*R*T(i,j))
    v(i,j) = 0.0d0
    EB(i,j) = (R*T(i,j))/(GAMMA-1)
    e(i,j) = RHO(i,j)*(EB(i,j)+(u(i,j)**2+v(i,j)**2)*0.5d0)
  end do
end do

!---------------x,y座標設定（不等間隔格子生成）---------------
do i = JL,JM
  do j = IL,IM
    x(i,j) = XM/dble(IM)*dble(i)
  end do
end do
do j = JL,JM
  do i = IL,IM
    jy = dble(j)/(dble(JM+1)-1.d0)
    y(i,j) = YM/(b-1.d0)*(b**jy-1.d0)
  end do
end do

!---------------x,y座標書き出し（不等間隔格子）---------------
open (12,file='kadai3_xy.dat',status='replace')

do j = JL,JM
  do i = IL,IM
    write (12,'(f13.2,$)') x(i,j)
    write (12,'(f13.2)') y(i,j)
  end do
end do

close (12)

!---------------一般座標へ変換・初期値設定---------------

!********①_領域内部→中心差分********
do i = IL+1,IM-1
  do j = JL+1,JM-1
    Yeta(i,j) = 0.5d0*(y(i,j+1)-y(i,j-1))
    Xeta(i,j) = 0.5d0*(x(i,j+1)-x(i,j-1))
    Yxi(i,j) = 0.5d0*(y(i+1,j)-y(i-1,j))
    Xxi(i,j) = 0.5d0*(x(i+1,j)-x(i-1,j))
  end do
end do

!********②_領域境界上_上下→公式・中心差分********
do i = IL+1,IM-1
  Yeta(i,JL) = (3.d0*(y(i,JL+1)-y(i,JL))-(y(i,JL+2)-y(i,JL+1)))/(2.d0*deta)
  Xeta(i,JL) = (3.d0*(x(i,JL+1)-x(i,JL))-(x(i,JL+2)-x(i,JL+1)))/(2.d0*deta)
  Yxi(i,JL) = 0.5d0*(y(i+1,JL)-y(i-1,JL))
  Xxi(i,JL) = 0.5d0*(x(i+1,JL)-x(i-1,JL))

  Yeta(i,JM) = -(3.d0*(y(i,JM-1)-y(i,JM))-(y(i,JM-2)-y(i,JM-1)))/(2.d0*deta)
  Xeta(i,JM) = -(3.d0*(x(i,JM-1)-x(i,JM))-(x(i,JM-2)-x(i,JM-1)))/(2.d0*deta)
  Yxi(i,JM) = 0.5d0*(y(i+1,JM)-y(i-1,JM))
  Xxi(i,JM) = 0.5d0*(x(i+1,JM)-x(i-1,JM))
end do

!********③_領域境界上_左右→公式・中心差分********
do j = JL+1,JM-1
  Yeta(IM,j) = 0.5d0*(y(i,j+1)-y(i,j-1))
  Xeta(IM,j) = 0.5d0*(x(i,j+1)-x(i,j-1))
  Yxi(IM,j) = -(3.d0*(y(IM-1,j)-y(IM,j))-(y(IM-2,j)-y(IM-1,j)))/(2.d0*dxi)
  Xxi(IM,j) = -(3.d0*(x(IM-1,j)-x(IM,j))-(x(IM-2,j)-x(IM-1,j)))/(2.d0*dxi)

  Yeta(IL,j) = 0.5d0*(y(i,j+1)-y(i,j-1))
  Xeta(IL,j) = 0.5d0*(x(i,j+1)-x(i,j-1))
  Yxi(IL,j) = (3.d0*(y(IL+1,j)-y(IL,j))-(y(IL+2,j)-y(IL+1,j)))/(2.d0*dxi)
  Xxi(IL,j) = (3.d0*(x(IL+1,j)-x(IL,j))-(x(IL+2,j)-x(IL+1,j)))/(2.d0*dxi)
end do

!********④_頂点×4→公式********
Yeta(IM,JM) = -(3.d0*(y(IM,JM-1)-y(IM,JM))-(y(IM,JM-2)-y(IM,JM-1)))/(2.d0*deta)
Xeta(IM,JM) = -(3.d0*(x(IM,JM-1)-x(IM,JM))-(x(IM,JM-2)-x(IM,JM-1)))/(2.d0*deta)
Yxi(IM,JM) = -(3.d0*(y(IM-1,JM)-y(IM,JM))-(y(IM-2,JM)-y(IM-1,JM)))/(2.d0*dxi)
Xxi(IM,JM) = -(3.d0*(x(IM-1,JM)-x(IM,JM))-(x(IM-2,JM)-x(IM-1,JM)))/(2.d0*dxi)

Yeta(IM,JL) = (3.0d0*(y(IM,JL+1)-y(IM,JL))-(y(IM,JL+2)-y(IM,JL+1)))/(2.0d0*deta)
Xeta(IM,JL) = (3.0d0*(x(IM,JL+1)-x(IM,JL))-(x(IM,JL+2)-x(IM,JL+1)))/(2.0d0*deta)
Yxi(IM,JL) = -(3.0d0*(y(IM-1,JL)-y(IM,JL))-(y(IM-2,JL)-y(IM-1,JL)))/(2.0d0*dxi)
Xxi(IM,JL) = -(3.0d0*(x(IM-1,JL)-x(IM,JL))-(x(IM-2,JL)-x(IM-1,JL)))/(2.0d0*dxi)

Yeta(IL,JL) = (3.0d0*(y(IL,JL+1)-y(IL,JL))-(y(IL,JL+2)-y(IL,JL+1)))/(2.0d0*deta)
Xeta(IL,JL) = (3.0d0*(x(IL,JL+1)-x(IL,JL))-(x(IL,JL+2)-x(IL,JL+1)))/(2.0d0*deta)
Yxi(IL,JL) = (3.0d0*(y(IL+1,JL)-y(IL,JL))-(y(IL+2,JL)-y(IL+1,JL)))/(2.0d0*dxi)
Xxi(IL,JL) = (3.0d0*(x(IL+1,JL)-x(IL,JL))-(x(IL+2,JL)-x(IL+1,JL)))/(2.0d0*dxi)

Yeta(IL,JM) = -(3.0d0*(y(IL,JM-1)-y(IL,JM))-(y(IL,JM-2)-y(IL,JM-1)))/(2.0d0*deta)
Xeta(IL,JM) = -(3.0d0*(x(IL,JM-1)-x(IL,JM))-(x(IL,JM-2)-x(IL,JM-1)))/(2.0d0*deta)
Yxi(IL,JM) = (3.0d0*(y(IL+1,JM)-y(IL,JM))-(y(IL+2,JM)-y(IL+1,JM)))/(2.0d0*dxi)
Xxi(IL,JM) = (3.0d0*(x(IL+1,JM)-x(IL,JM))-(x(IL+2,JM)-x(IL+1,JM)))/(2.0d0*dxi)

!********上記を用いて計算（測度から・ベクトル成分）********
do j = JL,JM
  do i = IL,IM
    Ja(i,j) = 1.0d0/(Xxi(i,j)*Yeta(i,j)-Yxi(i,j)*Xeta(i,j))
    XIx(i,j) = Ja(i,j)*Yeta(i,j)
    XIy(i,j) = -Ja(i,j)*Xeta(i,j)
    ETAx(i,j) = -Ja(i,j)*Yxi(i,j)
    ETAy(i,j) = Ja(i,j)*Xxi(i,j)

    UL(i,j) = XIx(i,j)*u(i,j)+XIy(i,j)*v(i,j)
    VL(i,j) = ETAx(i,j)*u(i,j)+ETAy(i,j)*v(i,j)

    QH(i,j,1) = RHO(i,j)/Ja(i,j)
    QH(i,j,2) = RHO(i,j)*u(i,j)/Ja(i,j)
    QH(i,j,3) = RHO(i,j)*v(i,j)/Ja(i,j)
    QH(i,j,4) = e(i,j)/Ja(i,j)

    EH(i,j,1) = RHO(i,j)*UL(i,j)/Ja(i,j)
    EH(i,j,2) = (p(i,j)*XIx(i,j)+RHO(i,j)*u(i,j)*UL(i,j))/Ja(i,j)
    EH(i,j,3) = (RHO(i,j)*v(i,j)*UL(i,j)+XIy(i,j)*p(i,j))/Ja(i,j)
    EH(i,j,4) = ((e(i,j)+p(i,j))*UL(i,j))/Ja(i,j)

    FH(i,j,1) = RHO(i,j)*VL(i,j)/Ja(i,j)
    FH(i,j,2) = (RHO(i,j)*u(i,j)*VL(i,j)+ETAx(i,j)*p(i,j))/Ja(i,j)
    FH(i,j,3) = (ETAy(i,j)*p(i,j)+RHO(i,j)*v(i,j)*VL(i,j))/Ja(i,j)
    FH(i,j,4) = ((e(i,j)+p(i,j))*VL(i,j))/Ja(i,j)
  end do
end do

!---------------メイン計算---------------
do while (n<=NMAX)

  call Boundary
  QH1 = QH
  EH1 = EH
  FH1 = FH
  call Euler(IL+1,IM-1,JL+1,JM-1)
  print '("時刻 n",I6)',n
  n = n+1
end do
call Output

contains

!---------------サブルーチン①_境界条件---------------
subroutine Boundary
  integer i,j

  !******** [A] 流入部（左端）→固定 ********
  do j = JL,JM
    M(IL,j) = 2.9d0
    T(IL,j) = 293.0d0
    RHO(IL,j) = 1.2d0
    p(IL,j) = RHO(0,j)*R*T(0,j)
    u(IL,j) = M(0,j)*sqrt(GAMMA*R*T(0,j))
    v(IL,j) = 0.0
    EB(IL,j) = (R*T(0,j))/(GAMMA-1)
    e(IL,j) = RHO(0,j)*(EB(0,j)+u(0,j)**2/2.0d0)
  end do

  !******** [B] 境界上（上_下）→1次外挿 ********
  do i = IL+1,IM-1
    v(i,JM) = v(i,JM-1)+dy*(v(i,JM-1)-v(i,JM-2))/dy
    u(i,JM) = u(i,JM-1)+dy*(u(i,JM-1)-u(i,JM-2))/dy
    RHO(i,JM) = RHO(i,JM-1)+dy*(RHO(i,JM-1)-RHO(i,JM-2))/dy
    e(i,JM) = e(i,JM-1)+dy*(e(i,JM-1)-e(i,JM-2))/dy

    v(i,JL) = v(i,JL+1)-dy*(v(i,JL+2)-v(i,JL+1))/dy
    u(i,JL) = u(i,JL+1)-dy*(u(i,JL+2)-u(i,JL+1))/dy
    RHO(i,JL) = RHO(i,JL+1)-dy*(RHO(i,JL+2)-RHO(i,JL+1))/dy
    e(i,JL) = e(i,JL+1)-dy*(e(i,JL+2)-e(i,JL+1))/dy
  end do

  !******** [C] 流出部（右端）→1次外挿 ********
  do j = JL,JM
    RHO(IM,j) = RHO(IM-1,j)-dx*(RHO(IM-2,j)-RHO(IM-1,j))/dx
    u(IM,j) = u(IM-1,j)-dx*(u(IM-2,j)-u(IM-1,j))/dx
    v(IM,j) = v(IM-1,j)-dx*(v(IM-2,j)-v(IM-1,j))/dx
    p(IM,j) = p(IM-1,j)-dx*(p(IM-2,j)-p(IM-1,j))/dx
    e(IM,j) = e(IM-1,j)-dx*(e(IM-2,j)-e(IM-1,j))/dx
  end do

end subroutine Boundary

!---------------サブルーチン②_オイラー方程式---------------
subroutine Euler(is,ie,js,je)
  integer is,ie,js,je
  do j = js,je
    do i = is,ie
      do k = 1,4
        QH(i,j,k) = QH1(i,j,k)-dt*((EH1(i+1,j,k)-EH1(i-1,j,k))/(2.0d0*dxi)+(FH1(i,j+1,k)-FH1(i,j-1,k))/(2.0d0*deta))
      end do
    end do
  end do

  do j = js,je
    do i = is,ie
      RHO(i,j) = QH(i,j,1)*Ja(i,j)
      u(i,j) = QH(i,j,2)*Ja(i,j)/RHO(i,j)
      v(i,j) = QH(i,j,3)*Ja(i,j)/RHO(i,j)

      UL(i,j) = EH(i,j,1)*Ja(i,j)/RHO(i,j)
      VL(i,j) = FH(i,j,1)*Ja(i,j)/RHO(i,j)

      e(i,j) = QH(i,j,4)*Ja(i,j)
      EB(i,j) = e(i,j)/RHO(i,j)-(u(i,j)**2+v(i,j)**2)*0.5d0
      T(i,j) = (GAMMA-1)*EB(i,j)/R
      p(i,j) = RHO(i,j)*R*T(i,j)
      M(i,j) = u(i,j)/sqrt(GAMMA*R*T(i,j))

      EH(i,j,1) = RHO(i,j)*UL(i,j)/Ja(i,j)
      EH(i,j,2) = (p(i,j)*XIx(i,j)+RHO(i,j)*u(i,j)*UL(i,j))/Ja(i,j)
      EH(i,j,3) = (RHO(i,j)*v(i,j)*UL(i,j)+XIy(i,j)*p(i,j))/Ja(i,j)
      EH(i,j,4) = ((e(i,j)+p(i,j))*UL(i,j))/Ja(i,j)

      FH(i,j,1) = RHO(i,j)*VL(i,j)/Ja(i,j)
      FH(i,j,2) = (RHO(i,j)*u(i,j)*VL(i,j)+ETAx(i,j)*p(i,j))/Ja(i,j)
      FH(i,j,3) = (ETAy(i,j)*p(i,j)+RHO(i,j)*v(i,j)*VL(i,j))/Ja(i,j)
      FH(i,j,4) = ((e(i,j)+p(i,j))*VL(i,j))/Ja(i,j)

    end do
  end do
end subroutine Euler

!---------------サブルーチン③_ファイル書き出し---------------
subroutine Output
  open (11,file='kadai3.dat',status='replace')

  do j = JL,JM
    do i = IL,IM
      write (11,'(f13.6,$)') u(i,j)
      write (11,'(f13.6)') v(i,j)
    end do
  end do

  close (11)
end subroutine Output

end program kadai3
