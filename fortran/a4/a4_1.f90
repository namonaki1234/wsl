program a4

!�ϐ��̒�`
implicit none

double precision,parameter::pi=3.14159265358979d0 !��

!�i�q����
double precision,parameter::H=10.0d0 !�����
integer,parameter::IM=30,JM=30 !�i�q��
double precision,parameter::dx=3.0*H/30.0 !x�����̊i�q��
double precision,parameter::r_y=1.1d0 !y�̌���
double precision::a !����
double precision,dimension(0:IM,0:JM) :: x,y

double precision,dimension(0:IM,0:JM)::xxi=0.0d0,xeta=0.0d0,yxi=0.0d0,yeta=0.0d0 !�́�xi�C�Ł�eta
double precision,dimension(0:IM,0:JM)::JJ !���R�r�A��
double precision,dimension(0:IM,0:JM)::xix=0.0d0,xiy=0.0d0,etax=0.0d0,etay=0.0d0

integer i,j

double precision,parameter::R=287.1d0 !�K�X�萔
double precision,parameter::g=1.4d0 !��M��

double precision::M=2.9d0 !�Ռ��g�O�}�b�n��
double precision::b !�Ռ��g�p��
double precision::alpha !�W��
double precision::un !���xu�̒l�ێ��p

double precision,dimension(0:IM,0:JM)::rho=0.0d0 !���x
double precision,dimension(0:IM,0:JM)::u=0.0d0 !x�������x
double precision,dimension(0:IM,0:JM)::v=0.0d0 !y�������x
double precision,dimension(0:IM,0:JM)::p=0.0d0 !����
double precision,dimension(0:IM,0:JM)::T=0.0d0 !���x
double precision,dimension(0:IM,0:JM)::ene=0.0d0 !�G�l���M�[

  ! === �i�q���� ===
do i = 0, IM
    do j = 0, JM
      x(i,j) = dx * dble(i)
    end do
  end do

  do j = 0, JM
    a = H * (r_y - 1.0d0) / (r_y**dble(JM) - 1.0d0)
    do i = 0, IM
      y(i,j) = a * (r_y**dble(j) - 1.0d0) / (r_y - 1.0d0)
    end do
  end do


!���W�ϊ�
do i=0,IM
 do j=0,JM
  if (i>0.and.i<IM) then
   xxi(i,j)=(x(i+1,j)-x(i-1,j))/2.0d0
   yxi(i,j)=(y(i+1,j)-y(i-1,j))/2.0d0
  end if
  if(j>0.and.j<JM) then
   xeta(i,j)=(x(i,j+1)-x(i,j-1))/2.0d0
   yeta(i,j)=(y(i,j+1)-y(i,j-1))/2.0d0
  end if
  if(i==0) then
   xxi(i,j)=(-3.0d0*x(i,j)+4.0d0*x(i+1,j)-x(i+2,j))/2.0d0
   yxi(i,j)=(-3.0d0*y(i,j)+4.0d0*y(i+1,j)-y(i+2,j))/2.0d0
  end if
  if(i==IM) then
   xxi(i,j)=(3.0d0*x(i,j)-4.0d0*x(i-1,j)+x(i-2,j))/2.0d0
   yxi(i,j)=(3.0d0*y(i,j)-4.0d0*y(i-1,j)+y(i-2,j))/2.0d0
  end if
  if(j==0) then
   xeta(i,j)=(-3.0d0*x(i,j)+4.0d0*x(i,j+1)-x(i,j+2))/2.0d0
   yeta(i,j)=(-3.0d0*y(i,j)+4.0d0*y(i,j+1)-y(i,j+2))/2.0d0
  end if
  if(j==JM) then
   xeta(i,j)=(3.0d0*x(i,j)-4.0d0*x(i,j-1)+x(i,j-2))/2.0d0
   yeta(i,j)=(3.0d0*y(i,j)-4.0d0*y(i,j-1)+y(i,j-2))/2.0d0
  end if

  JJ(i,j)=1.0d0/(xxi(i,j)*yeta(i,j)-xeta(i,j)*yxi(i,j))
  xix(i,j)=JJ(i,j)*yeta(i,j)
  xiy(i,j)=-JJ(i,j)*xeta(i,j)
  etax(i,j)=-JJ(i,j)*yxi(i,j)
  etay(i,j)=JJ(i,j)*xxi(i,j)

 end do
end do


!��������
b=pi/6.0d0 !�Ռ��g�p��
alpha=((g+1.0d0)*(M*dsin(b))**2.0)/((g-1.0d0)*(M*dsin(b))**2.0+2.0d0) !�W��
do i=0,IM
 do j=0,JM
  rho(i,j)=1.2d0
  T(i,j)=293.0d0
  u(i,j)=M*dsqrt(g*R*T(i,j))
  p(i,j)=rho(i,j)*R*T(i,j)
  ene(i,j)=rho(i,j)*(R*T(i,j)/(g-10.d0)+(u(i,j)**2.0+v(i,j)**2.0)/2.0d0)

  if(y(i,j)>=y(0,JM-2)-x(i,j)/dsqrt(3.0d0)) then !�Ռ��g��̕�����
   un=u(i,j)
   rho(i,j)=alpha*rho(i,j)
   u(i,j)=u(i,j)*dsin(b)*dsqrt((1.0d0/alpha)**2.0+(1.0d0/dtan(b))**2d0)*dcos(-b+datan(dtan(b)/alpha))
   v(i,j)=un*dsin(b)*dsqrt((1.0d0/alpha)**2.0+(1.0d0/dtan(b))**2.0)*dsin(-b+datan(dtan(b)/alpha))
   p(i,j)=(2.0d0*g*(M*dsin(b))**2.0-(g-1.0d0))/(g+1.0d0)*p(i,j)
   T(i,j)=(2.0d0*g*(M*dsin(b))**2.0d0-(g-1.0d0))*((g-1.0d0)*(M*dsin(b))**2.0+2.0d0)/(((g+1.0d0)*M*dsin(b))**2.0)*T(i,j)
   ene(i,j)=rho(i,j)*(R*T(i,j)/(g-1.0d0)+(u(i,j)**2.0+v(i,j)**2.0)/2.0d0)
  end if

 end do
end do


! === �f�[�^�o�́iMicroAVS�p DAT�j===
open(10,file='a4.dat',status='replace')
do j=0,JM
 do i=0,IM
  write(1,*) x(i,j),y(i,j),u(i,j),v(i,j)
 end do
end do
close(10)


  ! === MicroAVS��FLD�w�b�_�o�� ===
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

print *, "�� MicroAVS�p�� .dat ����� .fld ���o�͂��܂����B"





end program a4