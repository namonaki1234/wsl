program test
  implicit none

  integer :: R, X
  ! read (5, *) R, X
  R = 2000
  X = 1
  if (X == 1) then
    if ((1600 <= R) .and. (2999 >= R)) then
      write (6, *) "Yes"
    else
      write (6, *) "No"
    end if
  end if
  if (X == 2) then
    if ((1200 <= R) .and. (2399 >= R)) then
      write (6, *) "Yes"
    else
      write (6, *) "No"
    end if
  end if
  ! if (mod(a*b, 2) .eq. 0) then
  !   write (6, *) "Even"
  ! else
  !   write (6, *) "Odd"
  ! end if

end program test
