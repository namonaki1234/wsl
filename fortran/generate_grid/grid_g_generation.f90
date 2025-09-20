program grid_airfoil_naca0018
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)

  ! --- 翼形状パラメータ ---
  real(dp), parameter :: CHORD = 105.0_dp       ! 弦長 [mm]
  real(dp), parameter :: THICKNESS = 0.18_dp    ! 翼厚比 NACA0018
  real(dp), parameter :: SPAN = 250.0_dp        ! スパン長 [mm]
  real(dp), parameter :: PI = acos(-1.0_dp)

  ! --- 格子分割数 ---
  integer, parameter :: NI_UP = 51,  NJ_UP = 51,  NK_UP = 51
  integer, parameter :: NI_AF = 101, NJ_AF = 51, NK_AF = 51
  integer, parameter :: NI_DN = 101, NJ_DN = 51, NK_DN = 51
  integer, parameter :: NBLOCKS = 3

  ! --- 配列 ---
  real(dp), allocatable :: x_up(:,:,:), y_up(:,:,:), z_up(:,:,:)
  real(dp), allocatable :: x_af(:,:,:), y_af(:,:,:), z_af(:,:,:)
  real(dp), allocatable :: x_dn(:,:,:), y_dn(:,:,:), z_dn(:,:,:)
  real(dp), allocatable :: wing_x(:), wing_y(:)
  integer :: ni(NBLOCKS), nj(NBLOCKS), nk(NBLOCKS)
  integer :: i,j,k

  ! --- ブロック寸法 ---
  ni(1)=NI_UP; nj(1)=NJ_UP; nk(1)=NK_UP
  ni(2)=NI_AF; nj(2)=NJ_AF; nk(2)=NK_AF
  ni(3)=NI_DN; nj(3)=NJ_DN; nk(3)=NK_DN

  ! --- メモリ確保 ---
  allocate(x_up(NI_UP,NJ_UP,NK_UP), y_up(NI_UP,NJ_UP,NK_UP), z_up(NI_UP,NJ_UP,NK_UP))
  allocate(x_af(NI_AF,NJ_AF,NK_AF), y_af(NI_AF,NJ_AF,NK_AF), z_af(NI_AF,NJ_AF,NK_AF))
  allocate(x_dn(NI_DN,NJ_DN,NK_DN), y_dn(NI_DN,NJ_DN,NK_DN), z_dn(NI_DN,NJ_DN,NK_DN))
  allocate(wing_x(NI_AF), wing_y(NI_AF))

  ! --- 格子生成 ---
  call generate_h_up(CHORD, NI_UP,NJ_UP,NK_UP, x_up,y_up,z_up)
  call generate_naca0018(CHORD, THICKNESS, wing_x, wing_y)
  call generate_co_grid(wing_x, wing_y, NI_AF,NJ_AF,NK_AF, x_af,y_af,z_af)
  call generate_h_dn(CHORD, NI_DN,NJ_DN,NK_DN, x_dn,y_dn,z_dn)

  ! --- 出力 (Plot3D) ---
  open(10,file="airfoil_grid.g",status="replace",form="unformatted")
  write(10) NBLOCKS
  write(10) (ni(i),nj(i),nk(i), i=1,NBLOCKS)

  write(10) (((real(x_up(i,j,k),4), i=1,NI_UP), j=1,NJ_UP), k=1,NK_UP)
  write(10) (((real(y_up(i,j,k),4), i=1,NI_UP), j=1,NJ_UP), k=1,NK_UP)
  write(10) (((real(z_up(i,j,k),4), i=1,NI_UP), j=1,NJ_UP), k=1,NK_UP)

  write(10) (((real(x_af(i,j,k),4), i=1,NI_AF), j=1,NJ_AF), k=1,NK_AF)
  write(10) (((real(y_af(i,j,k),4), i=1,NI_AF), j=1,NJ_AF), k=1,NK_AF)
  write(10) (((real(z_af(i,j,k),4), i=1,NI_AF), j=1,NJ_AF), k=1,NK_AF)

  write(10) (((real(x_dn(i,j,k),4), i=1,NI_DN), j=1,NJ_DN), k=1,NK_DN)
  write(10) (((real(y_dn(i,j,k),4), i=1,NI_DN), j=1,NJ_DN), k=1,NK_DN)
  write(10) (((real(z_dn(i,j,k),4), i=1,NI_DN), j=1,NJ_DN), k=1,NK_DN)
  close(10)

  print *, "Grid written to airfoil_grid.g"

contains

  !--- NACA0018 翼型生成 ---
  subroutine generate_naca0018(chord, thick, x, y)
    real(dp), intent(in) :: chord, thick
    real(dp), intent(out) :: x(:), y(:)
    integer :: n, i, half
    real(dp) :: xt, yt, dx

    n = size(x)
    half = n/2        ! 上面と下面に分割

    dx = chord / (half - 1)

    ! 上面 (前縁→後縁)
    do i = 1, half
        xt = dx * (i - 1)
        yt = 5*thick*chord * ( 0.2969*sqrt(xt/chord) - 0.1260*(xt/chord) - &
                            0.3516*(xt/chord)**2 + 0.2843*(xt/chord)**3 - &
                            0.1015*(xt/chord)**4 )
        x(i) = xt
        y(i) =  yt
    end do

    ! 下面 (後縁→前縁)
    do i = 1, half
        xt = dx * (half - i)
        yt = 5*thick*chord * ( 0.2969*sqrt(xt/chord) - 0.1260*(xt/chord) - &
                            0.3516*(xt/chord)**2 + 0.2843*(xt/chord)**3 - &
                            0.1015*(xt/chord)**4 )
        x(half+i) = xt
        y(half+i) = -yt
    end do
  end subroutine generate_naca0018

  !--- 上流 H 格子 ---
  subroutine generate_h_up(chord,ni,nj,nk,xh,yh,zh)
    real(dp), intent(in) :: chord
    integer,intent(in)::ni,nj,nk
    real(dp),intent(out)::xh(ni,nj,nk),yh(ni,nj,nk),zh(ni,nj,nk)
    real(dp)::dx,dy,dz
    integer::i,j,k
    dx = -chord/(ni-1)
    dy = chord/(nj-1)
    dz = SPAN/(nk-1)
    do k=1,nk; do j=1,nj; do i=1,ni
      xh(i,j,k)=-chord + dx*(i-1)
      yh(i,j,k)=-0.5*chord + dy*(j-1)
      zh(i,j,k)=dz*(k-1)
    end do; end do; end do
  end subroutine

  !--- 下流 H 格子 ---
  subroutine generate_h_dn(chord,ni,nj,nk,xh,yh,zh)
    real(dp), intent(in) :: chord
    integer,intent(in)::ni,nj,nk
    real(dp),intent(out)::xh(ni,nj,nk),yh(ni,nj,nk),zh(ni,nj,nk)
    real(dp)::dx,dy,dz
    integer::i,j,k
    dx = 4*chord/(ni-1)
    dy = chord/(nj-1)
    dz = SPAN/(nk-1)
    do k=1,nk; do j=1,nj; do i=1,ni
      xh(i,j,k)=chord + dx*(i-1)
      yh(i,j,k)=-0.5*chord + dy*(j-1)
      zh(i,j,k)=dz*(k-1)
    end do; end do; end do
  end subroutine

  !--- 翼周り C-O 格子 ---
  subroutine generate_co_grid(wx, wy, ni, nj, nk, xc, yc, zc)
    real(dp), intent(in) :: wx(:), wy(:)
    integer, intent(in) :: ni, nj, nk
    real(dp), intent(out) :: xc(ni,nj,nk), yc(ni,nj,nk), zc(ni,nj,nk)
    real(dp) :: r_out, theta, eta, r_in, dz
    integer :: i, j, k

    r_out = 5.0_dp * CHORD
    dz = SPAN / real(nk-1,dp)

    do k = 1, nk
        do j = 1, nj
        eta = real(j-1,dp)/real(nj-1,dp)   ! 0→1
        do i = 1, ni
            theta = atan2(wy(i), wx(i))
            r_in  = sqrt(wx(i)**2 + wy(i)**2)
            ! 内側: 翼表面、外側: 半径 r_out の円周
            xc(i,j,k) = (1.0-eta)*wx(i) + eta*r_out*cos(theta)
            yc(i,j,k) = (1.0-eta)*wy(i) + eta*r_out*sin(theta)
            zc(i,j,k) = dz * real(k-1,dp)
        end do
        end do
    end do
  end subroutine

end program grid_airfoil_naca0018
