program plt_to_g
  implicit none

  character(len=256) :: line
  integer :: ni, nj, nk, nblocks
  double precision :: val_x, val_y, val_z
  integer :: i, j, k
  
  ! --- 1. .plt ファイルの読み込み（ヘッダー情報とグリッドサイズ） ---
  open(unit=10, file='airfoil_grid.plt', status='old', form='formatted')

  ! ヘッダー行を読み飛ばす
  read(10, '(a)') line  ! TITLE
  read(10, '(a)') line  ! VARIABLES
  read(10, '(a)') line  ! ZONE

  ! ZONE行からサイズを読み取る
  read(line, *)
  read(line(index(line, 'I=')+2:), *) ni
  read(line(index(line, 'J=')+2:), *) nj
  read(line(index(line, 'K=')+2:), *) nk

  print *, 'Detected grid size from .plt file: I=', ni, ', J=', nj, ', K=', nk

  ! --- 2. Plot3D形式（.g）への変換と書き込み ---
  ! Plot3Dグリッドファイルは非フォーマット（unformatted）でなければならない
  open(unit=7, file='airfoil_grid.g', status='replace', form='unformatted')

  ! Plot3Dヘッダー情報を書き込む
  ! ブロック数
  nblocks = 1
  write(7) nblocks
  
  ! 各ブロックのni, nj, nkを書き込む
  write(7) ni, nj, nk
  
  ! 座標データを読み込み、Plot3D形式で書き込む
  ! Plot3DはまずX座標をすべて、次にY座標、最後にZ座標を書き込む
  
  ! X座標を書き込む
  do k = 1, nk
    do j = 1, nj
      do i = 1, ni
        read(10, *) val_x, val_y, val_z ! この行は、x,y,zのセットとして読み込まれる
        write(7) val_x
      end do
    end do
  end do
  
  rewind(10)
  ! ヘッダー行を再読み込み
  read(10, '(a)') line
  read(10, '(a)') line
  read(10, '(a)') line

  ! Y座標を書き込む
  do k = 1, nk
    do j = 1, nj
      do i = 1, ni
        read(10, *) val_x, val_y, val_z
        write(7) val_y
      end do
    end do
  end do
  
  rewind(10)
  ! ヘッダー行を再読み込み
  read(10, '(a)') line
  read(10, '(a)') line
  read(10, '(a)') line

  ! Z座標を書き込む
  do k = 1, nk
    do j = 1, nj
      do i = 1, ni
        read(10, *) val_x, val_y, val_z
        write(7) val_z
      end do
    end do
  end do

  close(10)
  close(7)

  print *, 'Conversion successful. Output written to airfoil_grid.g'

end program plt_to_g
