program kadai4
implicit none

!--------------格子空間・計算空間条件--------------
double precision,parameter :: H = 10.0d0 !格子数
double precision,parameter :: XM = 3.0d0*H !長さ_x方向
double precision,parameter :: YM = H !長さ_y方向

integer,parameter :: IM = 30 !終点_x方向
integer,parameter :: JM = 30 !終点_y方向
integer,parameter :: IL = 0 !始点_x方向
integer,parameter :: JL = 0 !始点_y方向

double precision,dimension(IL-1:IM+1,JL-1:JM+1) :: x = 0.0d0 !等間隔座標_x
double precision,dimension(IL-1:IM+1,JL-1:JM+1) :: y = 0.0d0 !不等間隔座標_y
double precision jy !bのべき乗数
double precision,parameter :: a = 1.2d0 !格子間隔_始点近傍
double precision,parameter :: b = (a+1.0d0)/(a-1.0d0) !格子間隔_始点近傍
integer :: i,j !座標_x方向,y方向

double precision,parameter :: dx = XM/dble(IM) !微小長さ_x方向
double precision :: dy !微小長さ_y方向

!--------------流れ場の初期条件--------------
double precision,parameter :: R = 287.1d0 !ガス定数
double precision,parameter :: GAMMA = 1.4d0 !比熱比
double precision,parameter :: pi = dacos(-1.0d0) !円周率π
double precision,parameter :: beta = 30.0d0/180.0d0*pi !衝撃波角

!--------------配列定義--------------
double precision,dimension(IL:IM,JL:JM) :: u = 0.0d0 !速度_x方向
double precision,dimension(IL:IM,JL:JM) :: v = 0.0d0 !速度_y方向
double precision,dimension(IL:IM,JL:JM) :: p = 0.0d0 !圧力
double precision,dimension(IL:IM,JL:JM) :: T = 0.0d0 !絶対温度
double precision,dimension(IL:IM,JL:JM) :: RHO = 0.0d0 !密度
double precision,dimension(IL:IM,JL:JM) :: RHO2 = 0.0d0 !密度2
double precision,dimension(IL:IM,JL:JM) :: M = 0.0d0 !マッハ数
double precision,dimension(IL:IM,JL:JM) :: VELO = 0.0d0 !速度ベク
double precision,dimension(IL:IM,JL:JM) :: VELOn = 0.0d0 !速度ベク_垂直
double precision,dimension(IL:IM,JL:JM) :: VELOt = 0.0d0 !速度ベク_平行
double precision,dimension(IL:IM,JL:JM) :: theta = 0.0d0 !偏角θ
double precision,dimension(IL:IM,JL:JM) :: c = 0.0d0 !音速c
double precision,dimension(IL:IM,JL:JM) :: c_s = 0.0d0 !臨界音速c*

!--------------不等間隔格子の生成・ファイル出力--------------
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

open (12,file='kadai4_xy.dat',status='replace')

do j = JL,JM
  do i = IL,IM
    write (12,'(f13.2,$)') x(i,j)
    write (12,'(f13.2)') y(i,j)
  end do
end do
close (12)

!--------------初期条件設定--------------
do j = JL,JM
  do i = IL,IM
    M(i,j) = 2.9d0
    T(i,j) = 293.d0
    RHO(i,j) = 1.2d0
    p(i,j) = RHO(i,j)*R*T(i,j)
    u(i,j) = M(i,j)*sqrt(GAMMA*R*T(i,j))
    v(i,j) = 0.d0
    VELO(i,j) = sqrt((u(i,j)**2)+(v(i,j)**2))
    VELOn(i,j) = VELO(i,j)*sin(beta)
    VELOt(i,j) = VELOn(i,j)/tan(beta)
    c(i,j) = VELO(i,j)/M(i,j)
  end do
end do

!--------------衝撃波条件（衝撃波通過後の変化）--------------
do j = JL,JM
  do i = IL,IM
    if (y(i,j)<y(i,JM-2) .and. x(i,j)<=sqrt(3.d0)*(y(i,JM-2)-y(i,j))) then

    else

      T(i,j) = T(i,j)*(2.d0*GAMMA*(M(i,j)**2)*(sin(beta)**2)-(GAMMA-1.d0)) &
               *((GAMMA-1.0d0)*(M(i,j)**2)*(sin(beta)**2)+2.d0) &
               /(((GAMMA+1.0d0)**2)*(M(i,j)**2)*(sin(beta)**2))

      p(i,j) = p(i,j)*(2.d0*GAMMA*(M(i,j)**2)*(sin(beta)**2)-(GAMMA-1.d0)) &
               /(GAMMA+1.0d0)

      RHO2(i,j) = RHO(i,j)*(GAMMA+1.d0)*(M(i,j)**2)*(sin(beta)**2) &
                  /((GAMMA-1.0d0)*(M(i,j)**2)*(sin(beta)**2)+2.d0)

      VELOn(i,j) = RHO(i,j)*VELOn(i,j)/RHO2(i,j)

      VELO(i,j) = sqrt((VELOn(i,j)**2)+(VELOt(i,j)**2))

      theta(i,j) = beta-atan(VELOn(i,j)/VELOt(i,j))
      u(i,j) = VELO(i,j)*cos(-theta(i,j))
      v(i,j) = VELO(i,j)*sin(-theta(i,j))

    end if
  end do
end do

!--------------ファイル書き出し--------------
open (11,file='kadai4.dat',status='replace')

do j = JL,JM
  do i = IL,IM
    write (11,'(f13.6,$)') u(i,j)
    write (11,'(f13.6)') v(i,j)
  end do
end do

close (11)

end program kadai4
