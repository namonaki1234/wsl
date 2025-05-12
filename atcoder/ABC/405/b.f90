program test
  implicit none

  integer :: n,m,i,count,j,challenge_divide
  integer,allocatable,dimension(:):: a
  read (5, *) n,m
  allocate(a(N))
  read (5, *) (a(i), i=1,n)

count=0
 
do j=1,1000
    challenge_divide=0
 
    do i=1, n
        if (mod(a(i),2)==0) then
            challenge_divide=challenge_divide+1
        endif
    enddo
 
    if (challenge_divide==n) then
        count=count+1
        do i=1,n
            a(i)=a(i)/2
        enddo
    else
        exit
    endif
enddo
 
write(6,*) count
  ! if (mod(a*b, 2) .eq. 0) then
  !   write (6, *) "Even"
  ! else
  !   write (6, *) "Odd"
  ! end if

end program test
