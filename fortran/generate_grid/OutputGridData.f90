!***************************************************
!*****           OutputGridData.f90            *****
!***************************************************
!格子の座標を出力するプログラム

program grid
implicit none

!変数宣言
integer i,j,k,m,e,outlet						!変数
integer nb,imax,jmax,kmax						!マルチブロック数,各座標の最大数
integer,allocatable::ni(:),nj(:),nk(:)														!各座標の最大値
double precision,allocatable::x(:,:,:,:),y(:,:,:,:),z(:,:,:,:)		!xyz座標
character(len=50) fn								!読み込むファイル名

print*,'***************************************************'
print*,'*****           OutputGridData.f90            *****'
print*,'***************************************************' 


!*****格子データの読み込み*****
open (unit=10, form='unformatted', file='grid.g')
read(10) nb
allocate(ni(nb),nj(nb),nk(nb))
read(10) (ni(m),nj(m),nk(m),m=1,nb)
imax=maxval(ni)
jmax=maxval(nj)
kmax=maxval(nk)
allocate(x(imax,jmax,kmax,nb),y(imax,jmax,kmax,nb),z(imax,jmax,kmax,nb))
do m=1,nb
 read(10) &
  ((( x(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)), &
  ((( y(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)), &
  ((( z(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m))

!print *,((( x(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)), &
!  ((( y(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)), &
!  ((( z(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m))

end do
close(10)


!*****格子データの出力*****
open (unit = 11, form='formatted',file='codelength.csv')
	do i=1,ni(3)
		write(11,*) i,',',sqrt((x(i,1,47,3)-x(1,1,47,1))**2.0d0+(y(i,1,47,3)-y(1,1,47,1))**2.0d0+(z(i,1,47,3)-z(1,1,47,1))**2.0d0)
	end do
close(11)

open (unit = 12, form='formatted',file='codelengths.csv')
	do i=1,ni(30)
  do k=1,nk(30)
		write(12,*) i,',',k,',',sqrt((x(i,30,k,30)-x(1,1,93,1))**2.0d0+(y(i,30,k,30)-y(1,1,93,1))**2.0d0+ &
                             (z(i,30,k,30)-z(1,1,93,1))**2.0d0)
	end do
	end do
close(12)


!*****特定の格子データの出力*****
!open (unit=10,form='formatted',file='log.txt')
!	write(10,*) nb
!	do m=29,29
!  	do i = 1,ni(m)
!    	do j = 1,nj(m)
!      	do k = 1,nk(m)
!					if (64<=i .and. i<=66 .and. 1<=j .and. j<=2 .and. 199<=k .and. k<=201) then
!        		write(10,'(i4)',advance='no')m
!	        	write(10,'(i4)',advance='no')i
!  	    	  write(10,'(i4)',advance='no')j
!    		    write(10,'(i4)',advance='no')k
!    	  	  write(10,'(f30.20)',advance='no')x(i,j,k,m)
!  	      	write(10,'(f30.20)',advance='no')y(i,j,k,m)
!	    	    write(10,'(f30.20)')z(i,j,k,m)
!					end if
!	      end do
!  	  end do
!	  end do
!	end do
!close(10)


end program grid