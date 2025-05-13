program a2
  implicit none

  !格子条件
  double precision, parameter::H = 10.0d0
  integer, parameter::IM = 30, JM = 30 !格子数
  double precision, parameter::dx = 3.0d0*H/dble(IM), dy = H/dble(JM) !x,y方向の格子幅

  !流れ場初期条件
  double precision, dimension(0:IM, 0:JM)::M = 2.9d0, T = 293.0d0, R = 287.1d0, k = 1.4d0, rho = 1.2d0, dt = 1.0d-6 !マッハ数,温度,ガス定数,比熱比,密度,時間刻み

  !その他の条件
  integer::i, j, n
  double precision, dimension(0:IM, 0:JM)::a, x, y, u, v, p, e, Q1, Q2, Q3, Q4, E1, E2, E3, E4, F1, F2, F3, F4, E_bar

  !格子作成
  do j = 0, JM
    do i = 0, IM
      x(i, j) = dx*dble(i)
      y(i, j) = dy*dble(j)
    end do
  end do

  !初期条件
  do j = 0, JM
    do i = 0, IM
      a(i, j) = sqrt(k(i, j)*R(i, j)*T(i, j)) !音速

      !流入条件
      u(0, j) = M(0, j)*a(0, j) !u=Ma
      u(1, j) = u(0, j) !一次外挿
      v(i, j) = 0.0
    end do
  end do

  !初期値以外の速度計算
  do j = 0, JM
    do i = 2, IM
      u(i, j) = M(i, j)*a(i, j) 
    end do
  end do

  !圧力とエネルギーの計算
  do j = 0, JM
    do i = 0, IM
      p(i, j) = rho(i, j)*R(i, j)*T(i, j) !圧力p 式(5.10)にE_barを代入し、変形した式
      e(i, j) = rho(i, j)*((R(i, j)*T(i, j)/k(i, j) - 1.0)) + ((u(i, j)**2 + v(i, j)**2)/2.0) !全エネルギーe 式(5.7)の2次元ver
    end do
  end do

  !式(5.26)の要素計算
  do n = 1, 10000
    do j = 0, JM
      do i = 0, IM
        !保存変数ベクトル
        Q1(i, j) = rho(i, j) 
        Q2(i, j) = rho(i, j)*u(i, j)
        Q3(i, j) = rho(i, j)*v(i, j)
        Q4(i, j) = e(i, j)
        !x方向流束ベクトル
        E1(i, j) = rho(i, j)*u(i, j) 
        E2(i, j) = p(i, j) + rho(i, j)*(u(i, j)**2.0)
        E3(i, j) = rho(i, j)*u(i, j)*v(i, j)
        E4(i, j) = (e(i, j) + p(i, j))*u(i, j)
        !y方向流束ベクトル
        F1(i, j) = rho(i, j)*v(i, j) 
        F2(i, j) = rho(i, j)*v(i, j)*u(i, j)
        F3(i, j) = p(i, j) + rho(i, j)*(v(i, j)**2.0)
        F4(i, j) = (e(i, j) + p(i, j))*v(i, j)
      end do
    end do

    !境界条件を除く内部における式(5.25)
    do j = 1, JM - 1
      do i = 1, IM - 1
        Q1(i, j) = Q1(i, j) - dt(i, j)*(((E1(i + 1, j) - E1(i - 1, j))/(2.0*dx) + (F1(i, j + 1) - F1(i, j - 1))/(2.0*dy)))
        Q2(i, j) = Q2(i, j) - dt(i, j)*(((E2(i + 1, j) - E2(i - 1, j))/(2.0*dx) + (F2(i, j + 1) - F2(i, j - 1))/(2.0*dy)))
        Q3(i, j) = Q3(i, j) - dt(i, j)*(((E3(i + 1, j) - E3(i - 1, j))/(2.0*dx) + (F3(i, j + 1) - F3(i, j - 1))/(2.0*dy)))
        Q4(i, j) = Q4(i, j) - dt(i, j)*(((E4(i + 1, j) - E4(i - 1, j))/(2.0*dx) + (F4(i, j + 1) - F4(i, j - 1))/(2.0*dy)))
      end do
    end do

    !境界の計算（一次外挿）
    !上壁
    do i = 1, IM
      Q1(i, JM) = Q1(i, JM - 2) + (y(i, JM) - y(i, JM - 2))/(y(i, JM - 2) - y(i, JM - 1))*(Q1(i, JM - 2) - Q1(i, JM - 1))
      Q2(i, JM) = Q2(i, JM - 2) + (y(i, JM) - y(i, JM - 2))/(y(i, JM - 2) - y(i, JM - 1))*(Q2(i, JM - 2) - Q2(i, JM - 1))
      Q3(i, JM) = Q3(i, JM - 2) + (y(i, JM) - y(i, JM - 2))/(y(i, JM - 2) - y(i, JM - 1))*(Q3(i, JM - 2) - Q3(i, JM - 1))
      Q4(i, JM) = Q4(i, JM - 2) + (y(i, JM) - y(i, JM - 2))/(y(i, JM - 2) - y(i, JM - 1))*(Q4(i, JM - 2) - Q4(i, JM - 1))
    end do

    !下壁
    do i = 1, IM
      Q1(i, 0) = Q1(i, 2) + (y(i, 0) - y(i, 2))/(y(i, 2) - y(i, 1))*(Q1(i, 2) - Q1(i, 1))
      Q2(i, 0) = Q2(i, 2) + (y(i, 0) - y(i, 2))/(y(i, 2) - y(i, 1))*(Q2(i, 2) - Q2(i, 1))
      Q3(i, 0) = Q3(i, 2) + (y(i, 0) - y(i, 2))/(y(i, 2) - y(i, 1))*(Q3(i, 2) - Q3(i, 1))
      Q4(i, 0) = Q4(i, 2) + (y(i, 0) - y(i, 2))/(y(i, 2) - y(i, 1))*(Q4(i, 2) - Q4(i, 1))
    end do

    !流出
    do j = i, JM - 1
      Q1(IM, j) = Q1(IM - 2, j) + (x(IM, j) - x(IM - 2, j))/(x(IM - 2, j) - x(IM - 1, j))*(Q1(IM - 2, j) - Q1(IM - 1, j))
      Q2(IM, j) = Q2(IM - 2, j) + (x(IM, j) - x(IM - 2, j))/(x(IM - 2, j) - x(IM - 1, j))*(Q2(IM - 2, j) - Q2(IM - 1, j))
      Q3(IM, j) = Q3(IM - 2, j) + (x(IM, j) - x(IM - 2, j))/(x(IM - 2, j) - x(IM - 1, j))*(Q3(IM - 2, j) - Q3(IM - 1, j))
      Q4(IM, j) = Q4(IM - 2, j) + (x(IM, j) - x(IM - 2, j))/(x(IM - 2, j) - x(IM - 1, j))*(Q4(IM - 2, j) - Q4(IM - 1, j))
    end do

    !更新
    do j = 0, JM
      do i = 0, IM
        rho(i, j) = Q1(i, j)
        u(i, j) = Q2(i, j)/rho(i, j)
        v(i, j) = Q3(i, j)/rho(i, j)
        e(i, j) = Q4(i, j)
        T(i, j) = E_bar(i, j)*(k(i, j) - 1.0)/r(i, j)
        E_bar(i, j) = e(i, j)/rho(i, j) - (u(i, j)**2.0 + v(i, j)**2.0)/2.0 !式(5.7)
        p(i, j) = (k(i, j) - 1.0)*rho(i, j)*E_bar(i, j) !式(5.10)
      end do
    end do
    print *, n
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

end program a2
