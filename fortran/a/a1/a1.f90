program a1
  implicit none

  integer :: i, j
  integer, parameter :: IM = 30, JM = 30
  double precision, parameter :: H = 20.0d0, dx = 3.0d0 * H / IM
  double precision, parameter :: r = 1.1d0
  double precision, dimension(0:IM,0:JM) :: x, y, z
  double precision :: a

  z = 0.0d0  ! 初期化

  ! === 格子生成 ===
  do i = 0, IM
    do j = 0, JM
      x(i,j) = dx * dble(i)
    end do
  end do

  do j = 0, JM
    a = H * (r - 1.0d0) / (r**dble(JM) - 1.0d0)
    do i = 0, IM
      y(i,j) = a * (r**dble(j) - 1.0d0) / (r - 1.0d0)
    end do
  end do

  ! === データ出力（MicroAVS用 DAT）===
  open(10, file='a1.dat', status='replace')
  do j = 0, JM
    do i = 0, IM
      write(10,'(3F15.8)') x(i,j), y(i,j), z(i,j)
    end do
  end do
  close(10)

  ! === MicroAVSのFLDヘッダ出力 ===
  open(11, file='a1.fld', status='replace')
  write(11,'(A)') '# AVS field file'
  write(11,'(A)') 'ndim = 2'
  write(11,'(A,I5)') 'dim1 =', IM+1
  write(11,'(A,I5)') 'dim2 =', JM+1
  write(11,'(A)') 'nspace = 2'
  write(11,'(A)') 'veclen = 1'
  write(11,'(A)') 'data = double'
  write(11,'(A)') 'field = irregular'
  write(11,'(A)') 'label = z'
  write(11,'(A)') 'variable 1 file=a1.dat filetype=ascii skip=0 offset=2 stride=3'
  write(11,'(A)') 'coord 1 file=a1.dat filetype=ascii skip=0 offset=0 stride=3'
  write(11,'(A)') 'coord 2 file=a1.dat filetype=ascii skip=0 offset=1 stride=3'
  close(11)

  print *, "→ MicroAVS用の .dat および .fld を出力しました。"

end program a1
