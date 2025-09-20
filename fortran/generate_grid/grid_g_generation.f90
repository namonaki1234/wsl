program grid_generator
  implicit none

  ! Parameters based on the provided drawings
  double precision, parameter :: CHORD = 105.0d0  ! 翼の弦長 [mm]
  double precision, parameter :: WING_THICKNESS = 0.18d0 ! NACA0018翼厚比
  double precision, parameter :: WING_SPAN = 250.0d0 ! 翼のスパン長 [mm]
  double precision, parameter :: PI = acos(-1.0d0)

  ! Grid generation parameters
  integer, parameter :: NI_UPSTREAM = 51   ! 上流ブロックのi方向格子点数
  integer, parameter :: NJ_UPSTREAM = 51   ! 上流ブロックのj方向格子点数
  integer, parameter :: NK_UPSTREAM = 51   ! 上流ブロックのk方向格子点数
  
  integer, parameter :: NI_AIRFOIL = 101   ! 翼周りブロックのi方向格子点数
  integer, parameter :: NJ_AIRFOIL = 51    ! 翼周りブロックのj方向格子点数
  integer, parameter :: NK_AIRFOIL = 51    ! 翼周りブロックのk方向格子点数
  
  integer, parameter :: NI_DOWNSTREAM = 101 ! 下流ブロックのi方向格子点数
  integer, parameter :: NJ_DOWNSTREAM = 51  ! 下流ブロックのj方向格子点数
  integer, parameter :: NK_DOWNSTREAM = 51  ! 下流ブロックのk方向格子点数

  ! Multiblock data
  integer, parameter :: NBLOCKS = 3
  integer :: ni(NBLOCKS), nj(NBLOCKS), nk(NBLOCKS)
  
  ! Variables for grid coordinates
  double precision, allocatable :: x_up(:,:,:), y_up(:,:,:), z_up(:,:,:)
  double precision, allocatable :: x_af(:,:,:), y_af(:,:,:), z_af(:,:,:)
  double precision, allocatable :: x_down(:,:,:), y_down(:,:,:), z_down(:,:,:)
  
  double precision, allocatable :: wing_x(:), wing_y(:)
  integer :: i, j, k

  ! Set block dimensions
  ni(1) = NI_UPSTREAM; nj(1) = NJ_UPSTREAM; nk(1) = NK_UPSTREAM
  ni(2) = NI_AIRFOIL; nj(2) = NJ_AIRFOIL; nk(2) = NK_AIRFOIL
  ni(3) = NI_DOWNSTREAM; nj(3) = NJ_DOWNSTREAM; nk(3) = NK_DOWNSTREAM

  ! Allocate arrays for each block
  allocate(x_up(ni(1), nj(1), nk(1)), y_up(ni(1), nj(1), nk(1)), z_up(ni(1), nj(1), nk(1)))
  allocate(x_af(ni(2), nj(2), nk(2)), y_af(ni(2), nj(2), nk(2)), z_af(ni(2), nj(2), nk(2)))
  allocate(x_down(ni(3), nj(3), nk(3)), y_down(ni(3), nj(3), nk(3)), z_down(ni(3), nj(3), nk(3)))

  allocate(wing_x(ni(2)), wing_y(ni(2)))

  ! --- Generate grid for each block ---

  ! Block 1: Upstream H-grid
  call generate_h_grid_upstream(CHORD, ni(1), nj(1), nk(1), x_up, y_up, z_up)

  ! Block 2: Airfoil C-O grid
  call generate_naca0018(CHORD, WING_THICKNESS, wing_x, wing_y)
  call generate_co_grid(wing_x, wing_y, ni(2), nj(2), nk(2), x_af, y_af, z_af)

  ! Block 3: Downstream H-grid
  call generate_h_grid_downstream(CHORD, ni(3), nj(3), nk(3), x_down, y_down, z_down)

  ! --- Write Plot3D multiblock file (.g) ---
  open(unit=7, file='airfoil_grid.g', status='replace', form='unformatted')

  ! Header: Number of blocks
  write(7) NBLOCKS

  ! Header: Dimensions of each block
  write(7) ni, nj, nk
  
  ! Data: Write coordinates for each block in a multiblock-specific order
  write(7) (((x_up(i,j,k), i=1,ni(1)), j=1,nj(1)), k=1,nk(1)), &
           (((y_up(i,j,k), i=1,ni(1)), j=1,nj(1)), k=1,nk(1)), &
           (((z_up(i,j,k), i=1,ni(1)), j=1,nj(1)), k=1,nk(1))
           
  write(7) (((x_af(i,j,k), i=1,ni(2)), j=1,nj(2)), k=1,nk(2)), &
           (((y_af(i,j,k), i=1,ni(2)), j=1,nj(2)), k=1,nk(2)), &
           (((z_af(i,j,k), i=1,ni(2)), j=1,nj(2)), k=1,nk(2))
           
  write(7) (((x_down(i,j,k), i=1,ni(3)), j=1,nj(3)), k=1,nk(3)), &
           (((y_down(i,j,k), i=1,ni(3)), j=1,nj(3)), k=1,nk(3)), &
           (((z_down(i,j,k), i=1,ni(3)), j=1,nj(3)), k=1,nk(3))
  
  close(7)

  print *, 'Multiblock grid generation complete. Output written to airfoil_grid.g'

contains

  subroutine generate_naca0018(chord, t, x, y)
    double precision, intent(in) :: chord, t
    double precision, intent(out) :: x(:), y(:)
    integer :: n, i
    double precision :: xt, yt, dx

    n = size(x)
    dx = chord / (dble(n) - 1.0d0)

    do i = 1, n
      xt = dx * (dble(i) - 1.0d0)
      yt = 5.0d0 * t * chord * (0.2969d0*dsqrt(xt/chord) - 0.1260d0*(xt/chord) - &
                                0.3516d0*(xt/chord)**2 + 0.2843d0*(xt/chord)**3 - &
                                0.1015d0*(xt/chord)**4)
      x(i) = xt
      y(i) = yt
    end do

  end subroutine generate_naca0018

  subroutine generate_h_grid_upstream(chord, ni, nj, nk, x_h, y_h, z_h)
    double precision, intent(in) :: chord
    integer, intent(in) :: ni, nj, nk
    double precision, intent(inout) :: x_h(:,:,:), y_h(:,:,:), z_h(:,:,:)
    integer :: i, j, k
    double precision :: dx, dy, dz

    dx = (-1.0d0 * chord) / (dble(ni) - 1.0d0)
    dy = CHORD / (dble(nj) - 1.0d0)
    dz = WING_SPAN / (dble(nk) - 1.0d0)

    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          x_h(i, j, k) = (-1.0d0 * chord) + dx * (dble(i) - 1.0d0)
          y_h(i, j, k) = (-0.5d0 * CHORD) + dy * (dble(j) - 1.0d0)
          z_h(i, j, k) = dz * (dble(k) - 1.0d0)
        end do
      end do
    end do
  end subroutine generate_h_grid_upstream

  subroutine generate_h_grid_downstream(chord, ni, nj, nk, x_h, y_h, z_h)
    double precision, intent(in) :: chord
    integer, intent(in) :: ni, nj, nk
    double precision, intent(inout) :: x_h(:,:,:), y_h(:,:,:), z_h(:,:,:)
    integer :: i, j, k
    double precision :: dx, dy, dz

    dx = (4.0d0 * chord) / (dble(ni) - 1.0d0)
    dy = CHORD / (dble(nj) - 1.0d0)
    dz = WING_SPAN / (dble(nk) - 1.0d0)

    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          x_h(i, j, k) = chord + dx * (dble(i) - 1.0d0)
          y_h(i, j, k) = (-0.5d0 * CHORD) + dy * (dble(j) - 1.0d0)
          z_h(i, j, k) = dz * (dble(k) - 1.0d0)
        end do
      end do
    end do
  end subroutine generate_h_grid_downstream

  subroutine generate_co_grid(wing_x, wing_y, ni, nj, nk, x_co, y_co, z_co)
    double precision, intent(in) :: wing_x(:), wing_y(:)
    integer, intent(in) :: ni, nj, nk
    double precision, intent(inout) :: x_co(:,:,:), y_co(:,:,:), z_co(:,:,:)
    double precision :: r_outer, r_inner, ratio
    integer :: i, j, k
    double precision :: dy_dist, dz

    r_outer = 1.5d0 * CHORD
    ratio = 1.05d0
    dy_dist = 0.5d0 * CHORD

    dz = WING_SPAN / (dble(nk) - 1.0d0)

    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          r_inner = dsqrt(wing_x(i)**2 + wing_y(i)**2)
          x_co(i,j,k) = wing_x(i) + (r_outer - r_inner) * (dble(j-1)/dble(nj-1))
          y_co(i,j,k) = wing_y(i) + (dy_dist - wing_y(i)) * (dble(j-1)/dble(nj-1))
          z_co(i,j,k) = dz * (dble(k) - 1.0d0)
        end do
      end do
    end do
  end subroutine generate_co_grid

end program grid_generator
