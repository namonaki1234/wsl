program grid_generator
  implicit none

  ! Parameters based on the provided drawings
  double precision, parameter :: CHORD = 105.0d0  ! 翼の弦長 [mm]
  double precision, parameter :: WING_THICKNESS = 0.18d0 ! NACA0018翼厚比
  double precision, parameter :: WING_SPAN = 250.0d0 ! 翼のスパン長 [mm]
  double precision, parameter :: FLAP_WIDTH = 105.0d0 ! 翼のフラップ幅 [mm]
  double precision, parameter :: MEASUREMENT_SECTION_LENGTH = 370.0d0 ! 測定部の長さ [mm]
  double precision, parameter :: INLET_DIA = 260.0d0 ! 入口直径 [mm]
  double precision, parameter :: TEST_SECTION_SIDE = 110.0d0 ! 測定部一辺 [mm]
  double precision, parameter :: PI = acos(-1.0d0)

  ! Grid generation parameters (tuned for demonstration)
  integer, parameter :: NI_O = 101   ! O-grid: 翼周り (i方向)
  integer, parameter :: NJ_O = 51    ! O-grid: 翼から外側へ (j方向)
  integer, parameter :: NI_H_UPSTREAM = 51 ! H-grid: 翼上流 (i方向)
  integer, parameter :: NI_H_DOWNSTREAM = 101 ! H-grid: 翼下流 (i方向)
  integer, parameter :: NJ_H = 51    ! H-grid: 翼から外側へ (j方向)
  integer, parameter :: NK = 51      ! スパン方向 (k方向)

  ! Variables
  double precision, allocatable :: x_grid(:,:,:)
  double precision, allocatable :: y_grid(:,:,:)
  double precision, allocatable :: z_grid(:,:,:)
  double precision, allocatable :: wing_x(:), wing_y(:)
  
  integer :: i, j, k

  ! Allocate arrays
  allocate(x_grid(NI_H_UPSTREAM + NI_O + NI_H_DOWNSTREAM - 2, NJ_H, NK))
  allocate(y_grid(NI_H_UPSTREAM + NI_O + NI_H_DOWNSTREAM - 2, NJ_H, NK))
  allocate(z_grid(NI_H_UPSTREAM + NI_O + NI_H_DOWNSTREAM - 2, NJ_H, NK))
  allocate(wing_x(NI_O), wing_y(NI_O))

  ! Subroutine to define NACA0018 coordinates
  call generate_naca0018(CHORD, WING_THICKNESS, wing_x, wing_y)

  ! Subroutine for O-grid generation around the airfoil
  call generate_o_grid(wing_x, wing_y, NI_O, NJ_O, NK, &
                        x_grid(NI_H_UPSTREAM : NI_H_UPSTREAM + NI_O - 1, 1:NJ_O, :), &
                        y_grid(NI_H_UPSTREAM : NI_H_UPSTREAM + NI_O - 1, 1:NJ_O, :))

  ! Subroutine for H-grid generation (upstream & downstream)
  call generate_h_grid(CHORD, NI_H_UPSTREAM, NI_H_DOWNSTREAM, NJ_H, NK, &
                       x_grid(1:NI_H_UPSTREAM, 1:NJ_H, :), &
                       y_grid(1:NI_H_UPSTREAM, 1:NJ_H, :), &
                       x_grid(NI_H_UPSTREAM + NI_O - 1 : NI_H_UPSTREAM + NI_O + NI_H_DOWNSTREAM - 2, 1:NJ_H, :), &
                       y_grid(NI_H_UPSTREAM + NI_O - 1 : NI_H_UPSTREAM + NI_O + NI_H_DOWNSTREAM - 2, 1:NJ_H, :))

  ! Subroutine to set Z-coordinates (spanwise direction)
  call set_z_coordinates(NK, WING_SPAN, z_grid)

  ! Output grid to a file
  open(unit=10, file='airfoil_grid.plt', status='replace')
  write(10, *) 'TITLE = "Airfoil C-O Grid"'
  write(10, *) 'VARIABLES = "X", "Y", "Z"'
  write(10, *) 'ZONE I=', NI_H_UPSTREAM + NI_O + NI_H_DOWNSTREAM - 2, &
               ', J=', NJ_H, ', K=', NK, ' ZONETYPE=Ordered'
  do k = 1, NK
    do j = 1, NJ_H
      do i = 1, (NI_H_UPSTREAM + NI_O + NI_H_DOWNSTREAM - 2)
        write(10, '(3f15.6)') x_grid(i, j, k), y_grid(i, j, k), z_grid(i, j, k)
      end do
    end do
  end do
  close(10)

  print *, 'Grid generation complete. Output written to airfoil_grid.plt'

contains

  ! Subroutine to generate NACA0018 coordinates
  subroutine generate_naca0018(chord, t, x, y)
    double precision, intent(in) :: chord, t
    double precision, intent(out) :: x(:), y(:)
    integer :: n, i
    double precision :: xt, yt, dx

    n = size(x)
    dx = chord / (dble(n) - 1.0d0)

    do i = 1, n
      ! 線形補間で翼弦方向の座標を生成
      xt = dx * (dble(i) - 1.0d0)
      
      ! NACA 4桁シリーズ翼型の厚さ分布
      yt = 5.0d0 * t * chord * (0.2969d0*dsqrt(xt/chord) - 0.1260d0*(xt/chord) - &
                                0.3516d0*(xt/chord)**2 + 0.2843d0*(xt/chord)**3 - &
                                0.1015d0*(xt/chord)**4)

      x(i) = xt
      y(i) = yt
    end do

  end subroutine generate_naca0018

  ! Subroutine for O-grid generation around the airfoil
  subroutine generate_o_grid(wing_x, wing_y, ni, nj, nk, x_o, y_o)
    double precision, intent(in) :: wing_x(:), wing_y(:)
    integer, intent(in) :: ni, nj, nk
    double precision, intent(inout) :: x_o(:,:,:), y_o(:,:,:)
    double precision :: r_outer, r_inner, r_log, ratio
    integer :: i, j, k

    ! Outer boundary (circular shape)
    r_outer = 1.5d0 * CHORD ! 適当な外側境界の半径
    ratio = 1.05d0           ! 指数関数的補間の比率

    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          ! Normalize j
          r_inner = dsqrt(wing_x(i)**2 + wing_y(i)**2)
          r_log = dlog(1.0d0 + dble(j-1)/dble(nj-1) * (ratio-1.0d0)) / dlog(ratio)
          
          ! Linear interpolation for now; need to map to O-grid coordinates
          x_o(i, j, k) = wing_x(i) + (r_outer - r_inner) * r_log
          y_o(i, j, k) = wing_y(i) + (r_outer - r_inner) * r_log
        end do
      end do
    end do
  end subroutine generate_o_grid

  ! Subroutine for H-grid generation
  subroutine generate_h_grid(chord, ni_up, ni_down, nj, nk, x_up, y_up, x_down, y_down)
    double precision, intent(in) :: chord
    integer, intent(in) :: ni_up, ni_down, nj, nk
    double precision, intent(inout) :: x_up(:,:,:), y_up(:,:,:), x_down(:,:,:), y_down(:,:,:)
    integer :: i, j, k

    ! Upstream H-grid
    ! (Linear interpolation)
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni_up
          x_up(i, j, k) = -3.0d0 * chord + (3.0d0 * chord) * dble(i-1) / dble(ni_up-1)
          y_up(i, j, k) = -TEST_SECTION_SIDE/2.0d0 + TEST_SECTION_SIDE * dble(j-1) / dble(nj-1)
        end do
      end do
    end do

    ! Downstream H-grid
    ! (Linear interpolation)
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni_down
          x_down(i, j, k) = chord + (3.0d0 * chord) * dble(i-1) / dble(ni_down-1)
          y_down(i, j, k) = -TEST_SECTION_SIDE/2.0d0 + TEST_SECTION_SIDE * dble(j-1) / dble(nj-1)
        end do
      end do
    end do

  end subroutine generate_h_grid

  ! Subroutine to set Z-coordinates
  subroutine set_z_coordinates(nk, span, z)
    integer, intent(in) :: nk
    double precision, intent(in) :: span
    double precision, intent(inout) :: z(:,:,:)
    double precision :: dz
    integer :: k

    dz = span / (dble(nk)-1.0d0)

    do k = 1, nk
      z(:,:,k) = dz * (dble(k)-1.0d0)
    end do
  end subroutine set_z_coordinates

end program grid_generator