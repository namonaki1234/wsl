program assyuku1
implicit none

double precision,parameter::H=10.0d0
integer,parameter::IM=30,JM=30 !�i�q��
double precision,parameter::dx=3.0*H/30.0 !x�����̊i�q��
double precision,parameter::b=1.2d0 !����b
double precision::a,ix
integer::i,j
real,dimension(0:IM,0:JM)::x,y,z=0.0

do i=0,IM
 do j=0,JM
 x(i,j)=dx*real(i) !x�����i���Ԋu�i�q�j
 !print *,x(i,j)
end do
end do

do j=0,JM
 do i=0,IM
 a=(b+1)/(b-1)
 ix=dble(j)/JM
 y(i,j)=H*(a**ix-1)/(a-1) !y�����i�s���Ԋu�i�q�j
end do
end do

!�O���t�쐬
open(1,file='assyuku1.dat',status='replace')
do j=0,JM
 do i=0,IM
		write(1,*) x(i,j),y(i,j),z(i,j)
 end do
end do
close(1)

end program assyuku1
