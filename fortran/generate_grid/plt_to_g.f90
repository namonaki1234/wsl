program plt_to_g
  implicit none

  character(len=256) :: line
  integer :: ni, nj, nk, nblocks
  integer :: i, j, k
  
  ! Use allocatable arrays to store all data in memory
  double precision, allocatable :: x_coords(:,:,:), y_coords(:,:,:), z_coords(:,:,:)
  double precision :: val_x, val_y, val_z

  ! --- 1. Read .plt file into memory arrays ---
  open(unit=10, file='airfoil_grid.plt', status='old', form='formatted')

  ! Skip header lines
  read(10, '(a)') line  ! TITLE
  read(10, '(a)') line  ! VARIABLES
  read(10, '(a)') line  ! ZONE

  ! Read dimensions from the ZONE line
  read(line, *)
  read(line(index(line, 'I=')+2:), *) ni
  read(line(index(line, 'J=')+2:), *) nj
  read(line(index(line, 'K=')+2:), *) nk
  
  print *, 'Detected grid size from .plt file: I=', ni, ', J=', nj, ', K=', nk

  ! Allocate memory for coordinates
  allocate(x_coords(ni, nj, nk), y_coords(ni, nj, nk), z_coords(ni, nj, nk))

  ! Read all coordinate data into the arrays
  do k = 1, nk
    do j = 1, nj
      do i = 1, ni
        read(10, *) val_x, val_y, val_z
        x_coords(i, j, k) = val_x
        y_coords(i, j, k) = val_y
        z_coords(i, j, k) = val_z
      end do
    end do
  end do

  close(10)

  ! --- 2. Write Plot3D file (.g) from memory ---
  open(unit=7, file='airfoil_grid.g', status='replace', form='unformatted')

  ! Plot3D header (write block count)
  nblocks = 1
  write(7) nblocks

  ! Write ni, nj, nk for each block
  write(7) ni, nj, nk
  
  ! Write coordinate data in Plot3D's specific XYZ order
  ! X coordinates
  write(7) x_coords
  ! Y coordinates
  write(7) y_coords
  ! Z coordinates
  write(7) z_coords

  close(7)

  print *, 'Conversion successful. Output written to airfoil_grid.g'

end program plt_to_g
