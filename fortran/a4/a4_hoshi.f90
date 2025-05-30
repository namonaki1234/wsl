program kadai4
implicit none

!--------------�i�q��ԁE�v�Z��ԏ���--------------
double precision,parameter :: H = 10.0d0 !�i�q��
double precision,parameter :: XM = 3.0d0*H !����_x����
double precision,parameter :: YM = H !����_y����

integer,parameter :: IM = 30 !�I�__x����
integer,parameter :: JM = 30 !�I�__y����
integer,parameter :: IL = 0 !�n�__x����
integer,parameter :: JL = 0 !�n�__y����

double precision,dimension(IL-1:IM+1,JL-1:JM+1) :: x = 0.0d0 !���Ԋu���W_x
double precision,dimension(IL-1:IM+1,JL-1:JM+1) :: y = 0.0d0 !�s���Ԋu���W_y
double precision jy !b�ׂ̂��搔
double precision,parameter :: a = 1.2d0 !�i�q�Ԋu_�n�_�ߖT
double precision,parameter :: b = (a+1.0d0)/(a-1.0d0) !�i�q�Ԋu_�n�_�ߖT
integer :: i,j !���W_x����,y����

double precision,parameter :: dx = XM/dble(IM) !��������_x����
double precision :: dy !��������_y����

!--------------�����̏�������--------------
double precision,parameter :: R = 287.1d0 !�K�X�萔
double precision,parameter :: GAMMA = 1.4d0 !��M��
double precision,parameter :: pi = dacos(-1.0d0) !�~������
double precision,parameter :: beta = 30.0d0/180.0d0*pi !�Ռ��g�p

!--------------�z���`--------------
double precision,dimension(IL:IM,JL:JM) :: u = 0.0d0 !���x_x����
double precision,dimension(IL:IM,JL:JM) :: v = 0.0d0 !���x_y����
double precision,dimension(IL:IM,JL:JM) :: p = 0.0d0 !����
double precision,dimension(IL:IM,JL:JM) :: T = 0.0d0 !��Ή��x
double precision,dimension(IL:IM,JL:JM) :: RHO = 0.0d0 !���x
double precision,dimension(IL:IM,JL:JM) :: RHO2 = 0.0d0 !���x2
double precision,dimension(IL:IM,JL:JM) :: M = 0.0d0 !�}�b�n��
double precision,dimension(IL:IM,JL:JM) :: VELO = 0.0d0 !���x�x�N
double precision,dimension(IL:IM,JL:JM) :: VELOn = 0.0d0 !���x�x�N_����
double precision,dimension(IL:IM,JL:JM) :: VELOt = 0.0d0 !���x�x�N_���s
double precision,dimension(IL:IM,JL:JM) :: theta = 0.0d0 !�Ίp��
double precision,dimension(IL:IM,JL:JM) :: c = 0.0d0 !����c

!--------------�s���Ԋu�i�q�̐����E�t�@�C���o��--------------
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

!--------------���������ݒ�--------------
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

!--------------�Ռ��g�����i�Ռ��g�ʉߌ�̕ω��j--------------
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

!--------------�t�@�C�������o��--------------
open (11,file='kadai4.dat',status='replace')

do j = JL,JM
  do i = IL,IM
    write (11,'(f13.6,$)') u(i,j)
    write (11,'(f13.6)') v(i,j)
  end do
end do

close (11)

end program kadai4
