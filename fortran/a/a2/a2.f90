program a2
  implicit none

  ! ==========================
  ! パラメータ設定
  !格子条件
  double precision, parameter :: H = 10.0d0
  integer, parameter :: IM = 30, JM = 30
  double precision, parameter :: dx = 3.0d0*H/dble(IM), dy = H/dble(JM)

  ! 流れ場の初期条件
  double precision, dimension(0:IM, 0:JM) :: M = 2.9d0, T = 293.0d0, R = 287.1d0, k = 1.4d0, rho = 1.2d0, dt = 1.0d-6 !マッハ数,温度,ガス定数,比熱比,密度,時間刻み

  ! 変数宣言
  integer :: i, j, n
  double precision, dimension(0:IM, 0:JM) :: a, x, y, u, v, p, energy, E_bar
  double precision, dimension(0:IM, 0:JM, 4) :: Q, E, F

  ! ==========================
  ! 等間隔格子作成
  do j = 0, JM
    do i = 0, IM
      x(i, j) = dx*dble(i)
      y(i, j) = dy*dble(j)
    end do
  end do

  ! ==========================
  ! 流入の初期条件
  do j = 0, JM
    do i = 0, IM
      a(i, j) = sqrt(k(i, j)*R(i, j)*T(i, j)) !音速
      u(0, j) = M(0, j)*a(0, j) !u=Ma
      u(1, j) = u(0, j) !一次外挿
      v(i, j) = 0.0d0
    end do
  end do

  !流入付近以外の速度計算
  do j = 0, JM
    do i = 2, IM
      u(i, j) = M(i, j)*a(i, j)
    end do
  end do

  !圧力とエネルギーの計算
  do j = 0, JM
    do i = 0, IM
      p(i, j) = rho(i, j)*R(i, j)*T(i, j) !圧力p,式(5.10)にE_barを代入し、変形した式
      energy(i, j) = rho(i, j)*((R(i, j)*T(i, j)/(k(i, j) - 1.0d0)) + (u(i, j)**2 + v(i, j)**2)/2.0d0) !全エネルギーe 式(5.7)の2次元ver
    end do
  end do

  ! ==========================
  ! メインループ、式(5.26)の計算
  do n = 1, 10000
    call calc_QEF(Q, E, F, rho, u, v, energy, p, IM, JM)

    !境界条件を除く内部における式(5.25)
    do j = 1, JM - 1
      do i = 1, IM - 1
        Q(i, j, 1) = Q(i, j, 1) - dt(i, j) * ((E(i+1, j, 1) - E(i-1, j, 1)) &
        / (2.0d0 * dx) + (F(i, j+1, 1) - F(i, j-1, 1)) / (2.0d0 * dy))
        Q(i, j, 2) = Q(i, j, 2) - dt(i, j) * ((E(i+1, j, 2) - E(i-1, j, 2)) &
        / (2.0d0 * dx) + (F(i, j+1, 2) - F(i, j-1, 2)) / (2.0d0 * dy))
        Q(i, j, 3) = Q(i, j, 3) - dt(i, j) * ((E(i+1, j, 3) - E(i-1, j, 3)) &
        / (2.0d0 * dx) + (F(i, j+1, 3) - F(i, j-1, 3)) / (2.0d0 * dy))
        Q(i, j, 4) = Q(i, j, 4) - dt(i, j) * ((E(i+1, j, 4) - E(i-1, j, 4)) &
        / (2.0d0 * dx) + (F(i, j+1, 4) - F(i, j-1, 4)) / (2.0d0 * dy))
      end do
    end do

    ! ==========================
    !境界条件の設定（一次外挿）
    !上壁
    do i = 1, IM
      Q(i, JM,1) = Q(i, JM - 2,1) + (y(i, JM) - y(i, JM - 2))/(y(i, JM - 2) - y(i, JM - 1))*(Q(i, JM - 2,1) - Q(i, JM - 1,1))
      Q(i, JM,2) = Q(i, JM - 2,2) + (y(i, JM) - y(i, JM - 2))/(y(i, JM - 2) - y(i, JM - 1))*(Q(i, JM - 2,2) - Q(i, JM - 1,2))
      Q(i, JM,3) = Q(i, JM - 2,3) + (y(i, JM) - y(i, JM - 2))/(y(i, JM - 2) - y(i, JM - 1))*(Q(i, JM - 2,3) - Q(i, JM - 1,3))
      Q(i, JM,4) = Q(i, JM - 2,4) + (y(i, JM) - y(i, JM - 2))/(y(i, JM - 2) - y(i, JM - 1))*(Q(i, JM - 2,4) - Q(i, JM - 1,4))
    end do

    !下壁
    do i = 1, IM
      Q(i, 0,1) = Q(i, 2,1) + (y(i, 0) - y(i, 2))/(y(i, 2) - y(i, 1))*(Q(i, 2,1) - Q(i, 1,1))
      Q(i, 0,2) = Q(i, 2,2) + (y(i, 0) - y(i, 2))/(y(i, 2) - y(i, 1))*(Q(i, 2,2) - Q(i, 1,2))
      Q(i, 0,3) = Q(i, 2,3) + (y(i, 0) - y(i, 2))/(y(i, 2) - y(i, 1))*(Q(i, 2,3) - Q(i, 1,3))
      Q(i, 0,4) = Q(i, 2,4) + (y(i, 0) - y(i, 2))/(y(i, 2) - y(i, 1))*(Q(i, 2,4) - Q(i, 1,4))
    end do

    !流出境界
    do j = i, JM - 1
      Q(IM, j,1) = Q(IM - 2, j,1) + (x(IM, j) - x(IM - 2, j))/(x(IM - 2, j) - x(IM - 1, j))*(Q(IM - 2, j,1) - Q(IM - 1, j,1))
      Q(IM, j,2) = Q(IM - 2, j,2) + (x(IM, j) - x(IM - 2, j))/(x(IM - 2, j) - x(IM - 1, j))*(Q(IM - 2, j,2) - Q(IM - 1, j,2))
      Q(IM, j,3) = Q(IM - 2, j,3) + (x(IM, j) - x(IM - 2, j))/(x(IM - 2, j) - x(IM - 1, j))*(Q(IM - 2, j,3) - Q(IM - 1, j,3))
      Q(IM, j,4) = Q(IM - 2, j,4) + (x(IM, j) - x(IM - 2, j))/(x(IM - 2, j) - x(IM - 1, j))*(Q(IM - 2, j,4) - Q(IM - 1, j,4))
    end do

    ! 時間発展に伴う変数の再計算(更新)
    do j = 0, JM
      do i = 0, IM
        rho(i, j) = Q(i, j, 1)
        u(i, j) = Q(i, j, 2)/rho(i, j)
        v(i, j) = Q(i, j, 3)/rho(i, j)
        energy(i, j) = Q(i, j, 4)
        E_bar(i, j) = energy(i, j)/rho(i, j) - (u(i, j)**2 + v(i, j)**2)/2.0d0
        T(i, j) = E_bar(i, j)*(k(i, j) - 1.0d0)/R(i, j)
        p(i, j) = (k(i, j) - 1.0d0)*rho(i, j)*E_bar(i, j)
      end do
    end do
  end do

  ! === データ出力（MicroAVS用 DAT）===
  open (10, file='a2.dat', status='replace')
  do j = 0, JM
    do i = 0, IM
      write (10, *) x(i, j), y(i, j), u(i, j), v(i, j)
    end do
  end do
  close (10)

  ! === MicroAVSのFLDヘッダ出力 ===
  open (11, file='a2.fld', status='replace')
  write (11, '(A)') '# AVS field file'
  write (11, '(A)') 'ndim = 2'
  write (11, '(A,I5)') 'dim1 =', IM + 1
  write (11, '(A,I5)') 'dim2 =', JM + 1
  write (11, '(A)') 'nspace = 2'
  write (11, '(A)') 'veclen = 2'
  write (11, '(A)') 'data = double'
  write (11, '(A)') 'field = irregular'
  write (11, '(A)') 'label = u v'
  write (11, '(A)') 'variable 1 file=a2.dat filetype=ascii skip=0 offset=2 stride=4'
  write (11, '(A)') 'variable 2 file=a2.dat filetype=ascii skip=0 offset=3 stride=4'
  write (11, '(A)') 'coord 1 file=a2.dat filetype=ascii skip=0 offset=0 stride=4'
  write (11, '(A)') 'coord 2 file=a2.dat filetype=ascii skip=0 offset=1 stride=4'
  close (11)

  print *, "→ MicroAVS用の .dat および .fld を出力しました。"


contains

  subroutine calc_QEF(Q, E, F, rho, u, v, energy, p, IM, JM)
    implicit none
    integer, intent(in) :: IM, JM
    double precision, intent(in) :: rho(0:IM,0:JM), u(0:IM,0:JM), v(0:IM,0:JM)
    double precision, intent(in) :: energy(0:IM,0:JM), p(0:IM,0:JM)
    double precision, intent(out) :: Q(0:IM,0:JM,4), E(0:IM,0:JM,4), F(0:IM,0:JM,4)
    integer :: i, j

    do j = 0, JM
      do i = 0, IM
        !保存変数ベクトル
        Q(i,j,1) = rho(i,j)
        Q(i,j,2) = rho(i,j)*u(i,j)
        Q(i,j,3) = rho(i,j)*v(i,j)
        Q(i,j,4) = energy(i,j)
        !x方向流束ベクトル
        E(i,j,1) = Q(i,j,2)
        E(i,j,2) = p(i,j) + rho(i,j)*u(i,j)**2
        E(i,j,3) = rho(i,j)*u(i,j)*v(i,j)
        E(i,j,4) = (energy(i,j) + p(i,j))*u(i,j)
        !y方向流束ベクトル
        F(i,j,1) = Q(i,j,3)
        F(i,j,2) = rho(i,j)*v(i,j)*u(i,j)
        F(i,j,3) = p(i,j) + rho(i,j)*v(i,j)**2
        F(i,j,4) = (energy(i,j) + p(i,j))*v(i,j)
      end do
    end do
  end subroutine calc_QEF

end program a2
