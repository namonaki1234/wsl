! ファイル名kadai3.3
program kadai3
implicit none
!---------------�i�q�E�v�Z��ԏ����E����---------------
double precision,parameter	::H=10.0d0	!�i�q��
double precision,parameter	::XM=3.0d0*H	!����_x����
double precision,parameter	::YM=H		!����_y����

integer,parameter		::IM=30		!�I�__x����
integer,parameter		::JM=30		!�I�__y����
integer,parameter		::IL=0		!�n�__x����
integer,parameter		::JL=0		!�n�__y����

double precision ,dimension(IL-1:IM+1,JL-1:JM+1)::x=0.0d0	!�s���Ԋu_x
double precision ,dimension(IL-1:IM+1,JL-1:JM+1)::y=0.0d0	!�s���Ԋu_y
double precision jy				!�i�q����_�ׂ���
double precision,parameter	::a=1.2d0	!�i�q����_�萔
double precision,parameter	::b=(a+1.0d0)/(a-1.0d0)

double precision,parameter	::R=287.1d0	!�K�X�萔
double precision,parameter	::GAMMA=1.4d0	!��M��

double precision,parameter	::dt=0.000001d0	!���ԍ���

double precision,parameter	::dx=XM/dble(IM)	!x������������
double precision dy				!y������������
double precision,parameter	::dxi=1.0d0	!�̕�����������
double precision,parameter	::deta=1.0d0	!�ŕ�����������

integer n					!�v�Z��
integer,parameter		::NMAX=10	!�v�Z��_���

!---------------�z���`---------------
double precision ,dimension(IL:IM,JL:JM)::u=0.0d0	!���x_x����
double precision ,dimension(IL:IM,JL:JM)::v=0.0d0	!���x_y����
double precision ,dimension(IL:IM,JL:JM)::p=0.0d0	!����
double precision ,dimension(IL:IM,JL:JM)::T=0.0d0	!��Ή��x[k]
double precision ,dimension(IL:IM,JL:JM)::RHO=0.0d0		!���x
double precision ,dimension(IL:IM,JL:JM)::M=0.0d0	!�}�b�n��
double precision ,dimension(IL:IM,JL:JM)::e=0.0d0	!�S�G�l
double precision ,dimension(IL:IM,JL:JM)::EB=0.0d0	!�G�l�^�P�ʎ���

double precision ,dimension(IL:IM,JL:JM)::UL=0.0d0	!�̕������ϑ��x
double precision ,dimension(IL:IM,JL:JM)::VL=0.0d0	!�ŕ������ϑ��x

double precision ,dimension(IL:IM,JL:JM,1:4)::QH=0.0d0	!�ۑ��ʃx�N
double precision ,dimension(IL:IM,JL:JM,1:4)::QH1=0.0d0	!QH_n-1
double precision ,dimension(IL:IM,JL:JM,1:4)::EH=0.0d0	!x���������x�N�g��
double precision ,dimension(IL:IM,JL:JM,1:4)::EH1=0.0d0	!EH_n-1
double precision ,dimension(IL:IM,JL:JM,1:4)::FH=0.0d0	!y���������x�N
double precision ,dimension(IL:IM,JL:JM,1:4)::FH1=0.0d0	!FH_n-1

double precision ,dimension(IL:IM,JL:JM)::Ja=0.0d0	!���R�r�A��
double precision ,dimension(IL:IM,JL:JM)::Yeta=0.0d0	!dy/d��
double precision ,dimension(IL:IM,JL:JM)::Xeta=0.0d0	!dx/d��
double precision ,dimension(IL:IM,JL:JM)::Yxi=0.0d0	!dy/d��
double precision ,dimension(IL:IM,JL:JM)::Xxi=0.0d0	!dx/d��

double precision ,dimension(IL:IM,JL:JM)::XIx=0.0d0	!d��/dx
double precision ,dimension(IL:IM,JL:JM)::XIy=0.0d0	!d��/dy
double precision ,dimension(IL:IM,JL:JM)::ETAx=0.0d0	!d��/dx
double precision ,dimension(IL:IM,JL:JM)::ETAy=0.0d0	!d��/dy

!---------------���������ݒ�---------------
integer				::i,j		!x,y�������W
integer				::k		!�s�񐬕�
do i=IL,IM
	do j=JL,JM
		M(i,j)=2.9d0
		T(i,j)=293.d0
		RHO(i,j)=1.2d0
		p(i,j)=RHO(i,j)*R*T(i,j)
		u(i,j)=M(i,j)*sqrt(GAMMA*R*T(i,j))
		v(i,j)=0.0d0
		EB(i,j)=(R*T(i,j))/(GAMMA-1)
		e(i,j)=RHO(i,j)*(EB(i,j)+(u(i,j)**2+v(i,j)**2)*0.5d0)
		end do
	end do

!---------------x,y���W�ݒ�i�s���Ԋu�i�q�����j---------------
do i=JL,JM
		do j=IL,IM
		x(i,j)=XM/dble(IM)*dble(i)
		end do
	end do
	do j=JL,JM
		do i=IL,IM
		jy=dble(j)/(dble(JM+1)-1.d0)
		y(i,j)=YM/(b-1.d0)*(b**jy-1.d0)
		end do
	end do

!---------------x,y���W�����o���i�s���Ԋu�i�q�j---------------
open(12,file='kadai3_xy.dat',status='replace')

	do j=JL,JM
		do i=IL,IM
	write(12,'(f13.2,$)')x(i,j)
	write(12,'(f13.2)')y(i,j)
		end do
	end do

close(12)

!---------------��ʍ��W�֕ϊ��E�����l�ݒ�---------------

!********�@_�̈���������S����********
do i=IL+1,IM-1
		do j=JL+1,JM-1
		Yeta(i,j)=0.5d0*(y(i,j+1)-y(i,j-1))
		Xeta(i,j)=0.5d0*(x(i,j+1)-x(i,j-1))
		Yxi(i,j)=0.5d0*(y(i+1,j)-y(i-1,j))
		Xxi(i,j)=0.5d0*(x(i+1,j)-x(i-1,j))
		end do
	end do

!********�A_�̈拫�E��_�㉺�������E���S����********
	do i=IL+1,IM-1
	Yeta(i,JL)=(3.d0*(y(i,JL+1)-y(i,JL))-(y(i,JL+2)-y(i,JL+1)))/(2.d0*deta)
	Xeta(i,JL)=(3.d0*(x(i,JL+1)-x(i,JL))-(x(i,JL+2)-x(i,JL+1)))/(2.d0*deta)
	Yxi(i,JL)=0.5d0*(y(i+1,JL)-y(i-1,JL))
	Xxi(i,JL)=0.5d0*(x(i+1,JL)-x(i-1,JL))

	Yeta(i,JM)=-(3.d0*(y(i,JM-1)-y(i,JM))-(y(i,JM-2)-y(i,JM-1)))/(2.d0*deta)
	Xeta(i,JM)=-(3.d0*(x(i,JM-1)-x(i,JM))-(x(i,JM-2)-x(i,JM-1)))/(2.d0*deta)
	Yxi(i,JM)=0.5d0*(y(i+1,JM)-y(i-1,JM))
	Xxi(i,JM)=0.5d0*(x(i+1,JM)-x(i-1,JM))
	end do

!********�B_�̈拫�E��_���E�������E���S����********
	do j=JL+1,JM-1
	Yeta(IM,j)=0.5d0*(y(i,j+1)-y(i,j-1))
	Xeta(IM,j)=0.5d0*(x(i,j+1)-x(i,j-1))
	Yxi(IM,j)=-(3.d0*(y(IM-1,j)-y(IM,j))-(y(IM-2,j)-y(IM-1,j)))/(2.d0*dxi)
	Xxi(IM,j)=-(3.d0*(x(IM-1,j)-x(IM,j))-(x(IM-2,j)-x(IM-1,j)))/(2.d0*dxi)

	Yeta(IL,j)=0.5d0*(y(i,j+1)-y(i,j-1))
	Xeta(IL,j)=0.5d0*(x(i,j+1)-x(i,j-1))
	Yxi(IL,j)=(3.d0*(y(IL+1,j)-y(IL,j))-(y(IL+2,j)-y(IL+1,j)))/(2.d0*dxi)
	Xxi(IL,j)=(3.d0*(x(IL+1,j)-x(IL,j))-(x(IL+2,j)-x(IL+1,j)))/(2.d0*dxi)
	end do

!********�C_���_�~4������********
	Yeta(IM,JM)=-(3.d0*(y(IM,JM-1)-y(IM,JM))-(y(IM,JM-2)-y(IM,JM-1)))/(2.d0*deta)
	Xeta(IM,JM)=-(3.d0*(x(IM,JM-1)-x(IM,JM))-(x(IM,JM-2)-x(IM,JM-1)))/(2.d0*deta)
	Yxi(IM,JM)=-(3.d0*(y(IM-1,JM)-y(IM,JM))-(y(IM-2,JM)-y(IM-1,JM)))/(2.d0*dxi)
	Xxi(IM,JM)=-(3.d0*(x(IM-1,JM)-x(IM,JM))-(x(IM-2,JM)-x(IM-1,JM)))/(2.d0*dxi)

	Yeta(IM,JL)=(3.0d0*(y(IM,JL+1)-y(IM,JL))-(y(IM,JL+2)-y(IM,JL+1)))/(2.0d0*deta)
	Xeta(IM,JL)=(3.0d0*(x(IM,JL+1)-x(IM,JL))-(x(IM,JL+2)-x(IM,JL+1)))/(2.0d0*deta)
	Yxi(IM,JL)=-(3.0d0*(y(IM-1,JL)-y(IM,JL))-(y(IM-2,JL)-y(IM-1,JL)))/(2.0d0*dxi)
	Xxi(IM,JL)=-(3.0d0*(x(IM-1,JL)-x(IM,JL))-(x(IM-2,JL)-x(IM-1,JL)))/(2.0d0*dxi)

	Yeta(IL,JL)=(3.0d0*(y(IL,JL+1)-y(IL,JL))-(y(IL,JL+2)-y(IL,JL+1)))/(2.0d0*deta)
	Xeta(IL,JL)=(3.0d0*(x(IL,JL+1)-x(IL,JL))-(x(IL,JL+2)-x(IL,JL+1)))/(2.0d0*deta)
	Yxi(IL,JL)=(3.0d0*(y(IL+1,JL)-y(IL,JL))-(y(IL+2,JL)-y(IL+1,JL)))/(2.0d0*dxi)
	Xxi(IL,JL)=(3.0d0*(x(IL+1,JL)-x(IL,JL))-(x(IL+2,JL)-x(IL+1,JL)))/(2.0d0*dxi)

	Yeta(IL,JM)=-(3.0d0*(y(IL,JM-1)-y(IL,JM))-(y(IL,JM-2)-y(IL,JM-1)))/(2.0d0*deta)
	Xeta(IL,JM)=-(3.0d0*(x(IL,JM-1)-x(IL,JM))-(x(IL,JM-2)-x(IL,JM-1)))/(2.0d0*deta)
	Yxi(IL,JM)=(3.0d0*(y(IL+1,JM)-y(IL,JM))-(y(IL+2,JM)-y(IL+1,JM)))/(2.0d0*dxi)
	Xxi(IL,JM)=(3.0d0*(x(IL+1,JM)-x(IL,JM))-(x(IL+2,JM)-x(IL+1,JM)))/(2.0d0*dxi)

!********��L��p���Čv�Z�i���x����E�x�N�g�������j********
	do j=JL,JM
		do i=IL,IM
		Ja(i,j)=1.0d0/(Xxi(i,j)*Yeta(i,j)-Yxi(i,j)*Xeta(i,j))
		XIx(i,j)=Ja(i,j)*Yeta(i,j)
		XIy(i,j)=-Ja(i,j)*Xeta(i,j)
		ETAx(i,j)=-Ja(i,j)*Yxi(i,j)
		ETAy(i,j)=Ja(i,j)*Xxi(i,j)

		UL(i,j)=XIx(i,j)*u(i,j)+XIy(i,j)*v(i,j)
		VL(i,j)=ETAx(i,j)*u(i,j)+ETAy(i,j)*v(i,j)

		QH(i,j,1)=RHO(i,j)/Ja(i,j)
		QH(i,j,2)=RHO(i,j)*u(i,j)/Ja(i,j)
		QH(i,j,3)=RHO(i,j)*v(i,j)/Ja(i,j)
		QH(i,j,4)=e(i,j)/Ja(i,j)

		EH(i,j,1)=RHO(i,j)*UL(i,j)/Ja(i,j)
		EH(i,j,2)=(p(i,j)*XIx(i,j)+RHO(i,j)*u(i,j)*UL(i,j))/Ja(i,j)
		EH(i,j,3)=(RHO(i,j)*v(i,j)*UL(i,j)+XIy(i,j)*p(i,j))/Ja(i,j)
		EH(i,j,4)=((e(i,j)+p(i,j))*UL(i,j))/Ja(i,j)

		FH(i,j,1)=RHO(i,j)*VL(i,j)/Ja(i,j)
		FH(i,j,2)=(RHO(i,j)*u(i,j)*VL(i,j)+ETAx(i,j)*p(i,j))/Ja(i,j)
		FH(i,j,3)=(ETAy(i,j)*p(i,j)+RHO(i,j)*v(i,j)*VL(i,j))/Ja(i,j)
		FH(i,j,4)=((e(i,j)+p(i,j))*VL(i,j))/Ja(i,j)
		end do
	end do

!---------------���C���v�Z---------------
do while (n<=NMAX)

call Boundary
QH1=QH
EH1=EH
FH1=FH
call Euler(IL+1,IM-1,JL+1,JM-1)
print '("���� n",I6)',n
n=n+1
end do
call Output

contains

!---------------�T�u���[�`���@_���E����---------------
subroutine Boundary
integer i,j

!******** [A] �������i���[�j���Œ� ********
do j=JL,JM
	M(IL,j)=2.9d0
	T(IL,j)=293.0d0
	RHO(IL,j)=1.2d0
	p(IL,j)=RHO(0,j)*R*T(0,j)
	u(IL,j)=M(0,j)*sqrt(GAMMA*R*T(0,j))
	v(IL,j)=0.0
	EB(IL,j)=(R*T(0,j))/(GAMMA-1)
	e(IL,j)=RHO(0,j)*(EB(0,j)+u(0,j)**2/2.0d0)
end do

!******** [B] ���E��i��_���j��1���O�} ********
do i=IL+1,IM-1
	v(i,JM)=v(i,JM-1)+dy*(v(i,JM-1)-v(i,JM-2))/dy
	u(i,JM)=u(i,JM-1)+dy*(u(i,JM-1)-u(i,JM-2))/dy
	RHO(i,JM)=RHO(i,JM-1)+dy*(RHO(i,JM-1)-RHO(i,JM-2))/dy
	e(i,JM)=e(i,JM-1)+dy*(e(i,JM-1)-e(i,JM-2))/dy

	v(i,JL)=v(i,JL+1)-dy*(v(i,JL+2)-v(i,JL+1))/dy
	u(i,JL)=u(i,JL+1)-dy*(u(i,JL+2)-u(i,JL+1))/dy
	RHO(i,JL)=RHO(i,JL+1)-dy*(RHO(i,JL+2)-RHO(i,JL+1))/dy
	e(i,JL)=e(i,JL+1)-dy*(e(i,JL+2)-e(i,JL+1))/dy
end do

!******** [C] ���o���i�E�[�j��1���O�} ********
do j=JL,JM
        RHO(IM,j)=RHO(IM-1,j)-dx*(RHO(IM-2,j)-RHO(IM-1,j))/dx
        u(IM,j)=u(IM-1,j)-dx*(u(IM-2,j)-u(IM-1,j))/dx
        v(IM,j)=v(IM-1,j)-dx*(v(IM-2,j)-v(IM-1,j))/dx
        p(IM,j)=p(IM-1,j)-dx*(p(IM-2,j)-p(IM-1,j))/dx
	e(IM,j)=e(IM-1,j)-dx*(e(IM-2,j)-e(IM-1,j))/dx
end do

end subroutine Boundary

!---------------�T�u���[�`���A_�I�C���[������---------------
subroutine Euler(is,ie,js,je)
integer	is,ie,js,je
do j=js,je
		do i=is,ie
			do k=1,4
			QH(i,j,k)=QH1(i,j,k)-dt*((EH1(i+1,j,k)-EH1(i-1,j,k))/(2.0d0*dxi)+(FH1(i,j+1,k)-FH1(i,j-1,k))/(2.0d0*deta))
			end do
		end do
	end do

	do j=js,je
		do i=is,ie
		RHO(i,j)=QH(i,j,1)*Ja(i,j)
		u(i,j)=QH(i,j,2)*Ja(i,j)/RHO(i,j)
		v(i,j)=QH(i,j,3)*Ja(i,j)/RHO(i,j)

		UL(i,j)=EH(i,j,1)*Ja(i,j)/RHO(i,j)
		VL(i,j)=FH(i,j,1)*Ja(i,j)/RHO(i,j)

		e(i,j)=QH(i,j,4)*Ja(i,j)
		EB(i,j)=e(i,j)/RHO(i,j)-(u(i,j)**2+v(i,j)**2)*0.5d0
		T(i,j)=(GAMMA-1)*EB(i,j)/R
		p(i,j)=RHO(i,j)*R*T(i,j)
		M(i,j)=u(i,j)/sqrt(GAMMA*R*T(i,j))

		EH(i,j,1)=RHO(i,j)*UL(i,j)/Ja(i,j)
		EH(i,j,2)=(p(i,j)*XIx(i,j)+RHO(i,j)*u(i,j)*UL(i,j))/Ja(i,j)
		EH(i,j,3)=(RHO(i,j)*v(i,j)*UL(i,j)+XIy(i,j)*p(i,j))/Ja(i,j)
		EH(i,j,4)=((e(i,j)+p(i,j))*UL(i,j))/Ja(i,j)

		FH(i,j,1)=RHO(i,j)*VL(i,j)/Ja(i,j)
		FH(i,j,2)=(RHO(i,j)*u(i,j)*VL(i,j)+ETAx(i,j)*p(i,j))/Ja(i,j)
		FH(i,j,3)=(ETAy(i,j)*p(i,j)+RHO(i,j)*v(i,j)*VL(i,j))/Ja(i,j)
		FH(i,j,4)=((e(i,j)+p(i,j))*VL(i,j))/Ja(i,j)

		end do
	end do
end subroutine Euler

!---------------�T�u���[�`���B_�t�@�C�������o��---------------
subroutine Output
open(11,file='kadai3.dat',status='replace')

	do j=JL,JM
		do i=IL,IM
		write(11,'(f13.6,$)')u(i,j)
		write(11,'(f13.6)')v(i,j)
		end do
	end do

close(11)
end subroutine Output

end program kadai3
