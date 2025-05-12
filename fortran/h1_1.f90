program assyuku1_avs
    use iso_fortran_env
    implicit none

    integer, parameter :: IM = 30, JM = 30
    real(real64), parameter :: H = 20.0, dx = 3.0 * H / IM
    real(real64), parameter :: r = 1.2

    real(real64), dimension(0:IM,0:JM) :: x, y, z = 0.0
    real(real64) :: a
    integer :: i, j

    ! === 格子生成 ===
    do i = 0, IM
      do j = 0, JM
        x(i,j) = dx * real(i, kind=real64)
      end do
    end do

    do j = 0, JM
      a = H * (r - 1.0) / (r**real(JM, real64) - 1.0)
      do i = 0, IM
        y(i,j) = a * (r**real(j, real64) - 1.0) / (r - 1.0)
      end do
    end do

    ! === データ出力（MicroAVS用 DAT）===
    open(10, file='assyuku1.dat', status='replace')
    do j = 0, JM
      do i = 0, IM
        write(10,'(3F15.8)') x(i,j), y(i,j), z(i,j)
      end do
    end do
    close(10)

    ! === MicroAVSのFLDヘッダ出力 ===
    open(11, file='assyuku1.fld', status='replace')
    write(11,'(A)') '# AVS field file'
    write(11,'(A)') 'ndim = 2'
    write(11,'(A,I5)') 'dim1 =', IM+1
    write(11,'(A,I5)') 'dim2 =', JM+1
    write(11,'(A)') 'nspace = 2'
    write(11,'(A)') 'veclen = 1'
    write(11,'(A)') 'data = double'
    write(11,'(A)') 'field = uniform'
    write(11,'(A)') 'label = z'
    write(11,'(A)') 'variable 1 file=assyuku1.dat filetype=ascii skip=0 offset=2 stride=3'
    close(11)

    print *, "→ MicroAVS用の .dat および .fld を出力しました。"

  end program assyuku1_avs
