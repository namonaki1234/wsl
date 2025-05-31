program a4
   implicit none

   double precision, parameter :: H = 10.0d0
   integer, parameter :: IM = 30, JM = 30
   double precision :: dx, r_y, a
   double precision, parameter :: pi = dacos(-1.0d0)
   double precision, parameter :: beta = 30.0d0/180.0d0*pi
   double precision, parameter :: R = 287.1d0
   double precision, parameter :: gam = 1.4d0

   integer :: i, j
   double precision, dimension(0:IM, 0:JM) :: x, y
   double precision, dimension(0:IM, 0:JM) :: rho1, rho2, u, v, p, T, M, V_abs, V_n, V_t, theta

   dx = 3.0d0*H/dble(IM)
   r_y = 1.1d0

   do i = 0, IM
      do j = 0, JM
         x(i, j) = dx*dble(i)
      end do
   end do

   do j = 0, JM
      a = H*(r_y - 1.0d0)/(r_y**dble(JM) - 1.0d0)
      do i = 0, IM
         y(i, j) = a*(r_y**dble(j) - 1.0d0)/(r_y - 1.0d0)
      end do
   end do

   do j = 0, JM
      do i = 0, IM
         M(i, j) = 2.9d0
         T(i, j) = 293.0d0
         rho1(i, j) = 1.2d0
         p(i, j) = rho1(i, j)*R*T(i, j)
         u(i, j) = M(i, j)*sqrt(gam*R*T(i, j))
         v(i, j) = 0.d0
         V_abs(i, j) = sqrt((u(i, j)**2) + (v(i, j)**2))
         V_n(i, j) = V_abs(i, j)*sin(beta)
         V_t(i, j) = V_abs(i, j)*cos(beta)
      end do
   end do

   do j = 0, JM
      do i = 0, IM
         if (y(i, j) <= y(i, JM - 2) .and. x(i, j) <= sqrt(3.d0)*(y(i, JM - 2) - y(i, j))) cycle

         T(i, j) = T(i, j)*(2.d0*gam*(M(i, j)**2)*(sin(beta)**2) - (gam - 1.d0)) &
                   *((gam - 1.0d0)*(M(i, j)**2)*(sin(beta)**2) + 2.d0) &
                   /(((gam + 1.0d0)**2)*(M(i, j)**2)*(sin(beta)**2))

         p(i, j) = p(i, j)*(2.d0*gam*(M(i, j)**2)*(sin(beta)**2) - (gam - 1.d0)) &
                   /(gam + 1.0d0)

         rho2(i, j) = rho1(i, j)*(gam + 1.d0)*(M(i, j)**2)*(sin(beta)**2) &
                      /((gam - 1.0d0)*(M(i, j)**2)*(sin(beta)**2) + 2.d0)

         V_n(i, j) = rho1(i, j)*V_n(i, j)/rho2(i, j)

         V_abs(i, j) = sqrt((V_n(i, j)**2) + (V_t(i, j)**2))

         theta(i, j) = beta - atan(V_n(i, j)/V_t(i, j))
         u(i, j) = V_abs(i, j)*cos(-theta(i, j))
         v(i, j) = V_abs(i, j)*sin(-theta(i, j))
      end do
   end do

! === データ出力（MicroAVS用 DAT）===
   open (10, file='a4.dat', status='replace')
   do j = 0, JM
      do i = 0, IM
         write (10, *) x(i, j), y(i, j), u(i, j), v(i, j)
      end do
   end do
   close (10)

   ! === MicroAVSのFLDヘッダ出力 ===
   open (11, file='a4.fld', status='replace')
   write (11, '(A)') '# AVS field file'
   write (11, '(A)') 'ndim = 2'
   write (11, '(A,I5)') 'dim1 =', IM + 1
   write (11, '(A,I5)') 'dim2 =', JM + 1
   write (11, '(A)') 'nspace = 2'
   write (11, '(A)') 'veclen = 2'
   write (11, '(A)') 'data = double'
   write (11, '(A)') 'field = irregular'
   write (11, '(A)') 'label = u v'
   write (11, '(A)') 'variable 1 file=a4.dat filetype=ascii skip=0 offset=2 stride=4'
   write (11, '(A)') 'variable 2 file=a4.dat filetype=ascii skip=0 offset=3 stride=4'
   write (11, '(A)') 'coord 1 file=a4.dat filetype=ascii skip=0 offset=0 stride=4'
   write (11, '(A)') 'coord 2 file=a4.dat filetype=ascii skip=0 offset=1 stride=4'
   close (11)

   print *, "→ MicroAVS用の .dat および .fld を出力しました。"
end program a4
