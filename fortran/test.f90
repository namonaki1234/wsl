! hw2_5
program mac_method
  implicit none
  double precision, parameter :: dh = 0.01
  double precision, parameter :: dt = 0.005
  double precision, parameter :: H = 1.0
  integer, parameter :: IM = int(11*H/dh)
  integer, parameter :: JM = int(2*H/dh)
  integer, parameter :: I1 = int(H/dh)
  integer, parameter :: J1 = int(H/dh)
  double precision, parameter :: UMAX = 1.0
  double precision, parameter :: Re = 100.0
  double precision, parameter :: EPS = 1.0e-6
  double precision, parameter :: EPSP = 1.0e-3
  integer, parameter :: NMAX = 100000
  integer, parameter :: MMAX = 100000

  double precision, dimension(0:IM, 0:JM) :: u = 0.0, v = 0.0, p = 0.0
  double precision, dimension(0:IM, 0:JM) :: us = 0.0, vs = 0.0, ps = 0.0
  integer :: i, j, n, m
  double precision :: dp, du, dv, ddp, ddu, ddv
  double precision :: u_x, u_y, v_x, v_y, p_x, p_y
  double precision :: u_u_x, v_u_y, u_v_x, v_v_y, lap_u, lap_v

  n = 0
  du = EPS + 1.0
  dv = EPS + 1.0

  do while (((du .gt. EPS) .or. (dv .gt. EPS)) .and. n .lt. NMAX)
    m = 0
    dp = EPSP + 1.0

    do while (dp .gt. EPSP .and. m .lt. MMAX)
      call bc
      m = m + 1
      dp = 0.0
      ps = p
      do i = 1, I1
        do j = J1 + 1, JM - 1
          call cal_p
        end do
      end do
      do i = I1 + 1, IM - 1
        do j = 1, JM - 1
          call cal_p
        end do
      end do
    end do

    call bc
    n = n + 1
    du = 0.0
    dv = 0.0
    us = u
    vs = v

    do i = 1, I1
      do j = J1 + 1, JM - 1
        call cal_uv
      end do
    end do
    do i = I1 + 1, IM - 1
      do j = 1, JM - 1
        call cal_uv
      end do
    end do

    print *, n, dp, du, dv
  end do

  open (10, file='mac_method.dat', status='replace')
  do j = 0, JM
    do i = 0, IM
      write (10, *) p(i, j), u(i, j), v(i, j)
    end do
  end do
  close (10)

  open (11, file='mac_method.fld', status='replace')
  write (11, '(A)') '# AVS field file'
  write (11, '(A)') 'ndim = 2'
  write (11, '(A,I5)') 'dim1 =', IM + 1
  write (11, '(A,I5)') 'dim2 =', JM + 1
  write (11, '(A)') 'nspace = 2'
  write (11, '(A)') 'veclen = 3'
  write (11, '(A)') 'data = double'
  write (11, '(A)') 'field = uniform'
  write (11, '(A)') 'label = p u v'
  write (11, '(A)') 'variable 1 file=mac_method.dat filetype=ascii skip=0 offset=0 stride=3'
  write (11, '(A)') 'variable 2 file=mac_method.dat filetype=ascii skip=0 offset=1 stride=3'
  write (11, '(A)') 'variable 3 file=mac_method.dat filetype=ascii skip=0 offset=2 stride=3'
  close (11)

  open (12, file='mac_method.mgf', status='replace')
  write (12, '(A)') '# Micro AVS Geom:1.00'
  write (12, '(A)') 'polyline'
  write (12, '(A)') 'pline'
  write (12, '(A)') 'vertex'
  write (12, '(A)') '5'
  write (12, *) 0, 0, 0
  write (12, *) 0, H/dh - 1, 0
  write (12, *) H/dh - 1, H/dh - 1, 0
  write (12, *) H/dh - 1, 0, 0
  write (12, *) 0, 0, 0
  close (12)

contains

  subroutine bc
    do j = J1, JM
      u(0, j) = -4.0*UMAX*(j - J1)*(j - JM)/(J1 - JM)**2
      v(0, j) = 0.0
      p(0, j) = p(1, j)
    end do
    do j = 0, JM
      u(IM, j) = u(IM - 1, j)
      v(IM, j) = v(IM - 1, j)
      p(IM, j) = p(IM - 1, j)
    end do
    do i = 0, IM
      u(i, JM) = 0.0
      v(i, JM) = 0.0
      p(i, JM) = p(i, JM - 1)
    end do
    do i = I1, IM
      u(i, 0) = 0.0
      v(i, 0) = 0.0
      p(i, 0) = p(i, 1)
    end do
    do i = 0, I1
      u(i, J1) = 0.0
      v(i, J1) = 0.0
      p(i, J1) = p(i, J1 + 1)
    end do
    do j = 0, J1
      u(I1, j) = 0.0
      v(I1, j) = 0.0
      p(I1, j) = p(I1 + 1, j)
    end do
    p(I1, J1) = (p(I1 + 1, J1) + p(I1, J1 + 1))/2.0
    p(I1, 0) = (p(I1, 0) + p(I1, 1))/2.0
  end subroutine bc

  subroutine cal_p
    u_x = (us(i + 1, j) - us(i - 1, j))/(2*dh)
    u_y = (us(i, j + 1) - us(i, j - 1))/(2*dh)
    v_x = (vs(i + 1, j) - vs(i - 1, j))/(2*dh)
    v_y = (vs(i, j + 1) - vs(i, j - 1))/(2*dh)
    ps(i, j) = p(i, j)
    p(i, j) = (ps(i + 1, j) + p(i - 1, j) + ps(i, j + 1) + p(i, j - 1))/4.0 + &
              (dh**2)/4.0*((u_x**2 + v_y**2 + 2*v_x*u_y) - (u_x + v_y)/dt)

    if (abs(ps(i, j)) > 1.0e-12) then
      ddp = abs((p(i, j) - ps(i, j)) / max(abs(ps(i, j)), 1.0e-12))
    else
    ddp = 0.0
    endif
    if (ddp .gt. dp) dp = ddp
  end subroutine cal_p

  subroutine cal_uv
    u_x = (us(i + 1, j) - us(i - 1, j))/(2*dh)
    u_y = (us(i, j + 1) - us(i, j - 1))/(2*dh)
    v_x = (vs(i + 1, j) - vs(i - 1, j))/(2*dh)
    v_y = (vs(i, j + 1) - vs(i, j - 1))/(2*dh)
    p_x = (p(i + 1, j) - p(i - 1, j))/(2*dh)
    p_y = (p(i, j + 1) - p(i, j - 1))/(2*dh)
    lap_u = (us(i + 1, j) + us(i - 1, j) + us(i, j + 1) + us(i, j - 1) - 4.0*us(i, j))/dh**2
    lap_v = (vs(i + 1, j) + vs(i - 1, j) + vs(i, j + 1) + vs(i, j - 1) - 4.0*vs(i, j))/dh**2
    u_u_x = us(i, j)*u_x - (abs(us(i, j))/2)*(us(i + 1, j) - 2*us(i, j) + us(i - 1, j))/dh
    v_u_y = vs(i, j)*u_y - (abs(vs(i, j))/2)*(us(i, j + 1) - 2*us(i, j) + us(i, j - 1))/dh
    u_v_x = us(i, j)*v_x - (abs(us(i, j))/2)*(vs(i + 1, j) - 2*vs(i, j) + vs(i - 1, j))/dh
    v_v_y = vs(i, j)*v_y - (abs(vs(i, j))/2)*(vs(i, j + 1) - 2*vs(i, j) + vs(i, j - 1))/dh

    us(i, j) = u(i, j)  ! 修正: 単一要素への代入
    if (abs(us(i, j)) > 1.0e-12) then
      ddu = abs((-(u_u_x + v_u_y + p_x)*dt + dt/Re*lap_u)/us(i, j))
    else
      ddu = 0.0
    endif
    u(i, j) = -(u_u_x + v_u_y + p_x)*dt + dt/Re*lap_u + us(i, j)
    if (ddu .gt. du) du = ddu

    vs(i, j) = v(i, j)
    if (abs(vs(i, j)) > 1.0e-12) then
      ddv = abs((-(u_v_x + v_v_y + p_y)*dt + dt/Re*lap_v)/vs(i, j))
    else
      ddv = 0.0
    endif
    v(i, j) = -(u_v_x + v_v_y + p_y)*dt + dt/Re*lap_v + vs(i, j)
    if (ddv .gt. dv) dv = ddv
  end subroutine cal_uv

end program mac_method
