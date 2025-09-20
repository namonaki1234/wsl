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
  ! この行は、ZONE I=..., J=..., K=... の形式を想定しています
  ! 内部ファイル機能と文字列操作を組み合わせて、より堅牢な読み取りを実現
  read(line, *)
  read(line(index(line, 'I=')+2:), *) ni
  read(line(index(line, 'J=')+2:), *) nj
  read(line(index(line, 'K=')+2:), *) nk

  print *, 'Detected grid size from .plt file: I=', ni, ', J=', nj, ', K=', nk

  ! --- 2. .g ファイルへの変換と書き込み ---
  open(unit=7, file='airfoil_grid.g', status='replace', form='unformatted')

  ! .g ファイル形式に合わせてヘッダー情報を書き込む
  nblocks = 1
  write(7) nblocks
  write(7) ni, nj, nk

  ! 座標データを読み込み、すぐに.gファイルに書き込む
  do k = 1, nk
    do j = 1, nj
      do i = 1, ni
        read(10, *) val_x, val_y, val_z
        write(7) val_x, val_y, val_z
      end do
    end do
  end do

  close(10)
  close(7)

  print *, 'Conversion successful. Output written to airfoil_grid.g'

end program plt_to_g
