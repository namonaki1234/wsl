program generate_grid
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, parameter :: NX = 30, NY = 30
    real(dp), parameter :: H = 1.0_dp       ! 高さ
    real(dp), parameter :: L = 1.0_dp       ! 幅
    real(dp), parameter :: r = 0.9_dp       ! Y方向の等比数列の公比（下寄せ）
    real(dp) :: dx, x(NX), y(NY)
    real(dp) :: sum_r, base
    integer :: i, j
    character(len=20) :: grid_file

    grid_file = 'grid_data.dat'
    open(10, file=grid_file, status='replace')

    ! X方向：等間隔格子生成
    dx = L / real(NX - 1, dp)
    do i = 1, NX
      x(i) = (i - 1) * dx
    end do

    ! Y方向：等比級数和を計算して基準長を求める
    sum_r = 0.0_dp
    do j = 1, NY
      sum_r = sum_r + r**real(j - 1, dp)
    end do
    base = H / sum_r

    y(1) = 0.0_dp
    do j = 2, NY
      y(j) = y(j-1) + base * r**real(j - 2, dp)
    end do

    ! 格子点の出力
    do j = 1, NY
      do i = 1, NX
        write(10,'(F12.6,1x,F12.6)') x(i), y(j)
      end do
      write(10,*)
    end do

    close(10)
    print *, "格子データを書き出しました → ", grid_file

  end program generate_grid
