!*********************************************************************
!*** mps_icing_simulation                                          ***
!*** ver. 2014.09.05                                               ***
!*** written by Yuki Koji                                          ***
!*** E-MPS base program is written by Takahashi Ryohei             ***
!*** Grid-based Method program is written by Hayashi Ryosuke       ***
!*********************************************************************
!*********************************************************************

!*********************************************************************
!***********  粒子法と格子法を無理に合わせたので変数ごちゃごちゃ******
!*********************************************************************

program main
  ! module_setting ***************************************************
  use mod_mps
  use mod_mps_dim
  use mod_icing_dim
  use Package_NACA
  use Package_FileIO
  use Package_Droplet
  use Package_Grid
  use Package_Flow
!  use mt95

  implicit none
  ! variable_definition **********************************************
  integer :: i,j,k,l,m,n
  integer :: nn,idim,ifluid
  integer :: nn_grid,neighbor
  integer :: output_limit
  integer :: phase,phase_min,phase_max
  real    :: time_sim,time_max
  double precision :: time_sim_dble
  real    :: cfl,c_vel
  real    :: dis
  real    :: grav,delta,sigma
  integer :: ite_max
  integer :: nump,nump_out
  real    :: dt,dt_fix,dtout
  real    :: rlambda,rlambda_ref
  real    :: dis_rat
  integer :: grid_num,wall_num,initial_num
  integer :: phase_ref
  real    :: dist_min
  real    :: coll_rat
  integer :: ww_num
  integer :: wall_n
  real    :: num_ref = 0.0
  real    :: beta    = 0.8
  real    :: sum_g
  integer :: sum_n
  real    :: col_wall
  real    :: j_temp


  ! variable_definition(with_initial_setting) ************************
  integer :: ite     = 0
  integer :: fcount  = 0
  integer :: nn_wall = 50

  ! array_definition *************************************************
  real   , allocatable :: x(:,:)
  real   , allocatable :: v(:,:)
  real   , allocatable :: t(:)
  real   , allocatable :: an(:)
  real   , allocatable :: n0(:,:)
  real   , allocatable :: surf(:,:)
  real   , allocatable :: visc(:,:)
  real   , allocatable :: accel(:,:)
  real   , allocatable :: press(:,:)
  real   , allocatable :: p_out(:)
  real   , allocatable :: out(:,:)
  real   , allocatable :: p(:)
  real   , allocatable :: rho(:)
  real   , allocatable :: nu(:)
  integer, allocatable :: itypep(:)
  integer, allocatable :: cal_sw(:,:)
  integer, allocatable :: grid_p(:)
  integer, allocatable :: neigh(:,:)
  integer, allocatable :: g_neigh(:,:)
  real   , allocatable :: den(:)
  real   , allocatable :: vis(:)
  integer, allocatable :: type(:)
  integer, allocatable :: grid_n(:)
  real   , allocatable :: non_cal(:)
  integer, allocatable :: f_type(:)
  real   , allocatable :: ker_c(:)
  real   , allocatable :: grid(:,:,:)
  real   , allocatable :: sol_vel(:)
  real   , allocatable :: ker_r(:)
  real   , allocatable :: ww_wall(:,:)
  real   , allocatable :: grid_dist(:,:)
  integer, allocatable :: near_wall(:,:)
  real   , allocatable :: wall_anc(:,:,:)
  real   , allocatable :: wall_nor(:,:)
  integer, allocatable :: g_wall(:,:)
  integer, allocatable :: solid_sur(:)
  real, allocatable    :: tb(:)
  real, allocatable    :: tdd(:)
  real, allocatable    :: Q1(:)
  real, allocatable    :: Q2(:)
  real, allocatable    :: Qa(:)
  real, allocatable    :: Qk(:)
  real, allocatable    :: Qc(:)
  real, allocatable    :: Qd(:)
  real, allocatable    :: Qr(:)
  real, allocatable    :: Mimw(:)
  real   , allocatable :: ganma(:)
  integer   , allocatable :: swi_koeki(:)
  real, allocatable    :: Qhc(:)      !追加：熱伝達考慮
  real, allocatable    :: hc(:)
  real   , allocatable :: uf(:,:,:)
  real   , allocatable :: vf(:,:,:)
  real   , allocatable :: wf(:,:,:)
  real   , allocatable :: xf(:,:,:)
  real   , allocatable :: yf(:,:,:)
  real   , allocatable :: zf(:,:,:)
  real   , allocatable :: rhof(:,:,:)
  real   , allocatable :: muf(:,:,:)
  real   , allocatable :: drag(:,:)

  real   , allocatable :: grid_x(:,:,:)    !追加
  real   , allocatable :: grid_y(:,:,:)    !追加
  ! 粒子投入用パラメータ *******************************************
  double precision     :: random1
  integer              :: random2
  real                 :: randomx
  real                 :: randomy
  real                 :: randomz
  real                 :: dom_xmax
  real                 :: dom_xmin
  real                 :: dom_ymax
  real                 :: dom_ymin
  real                 :: dom_zmax
  real                 :: dom_zmin
  real                 :: vxmax
  real                 :: vymax
  real                 :: vzmax
  real                 :: vmax
  real                 :: interval
  real                 :: dtin
  real                 :: nextin
  integer              :: incount
  integer              :: inptnum
  integer              :: temp_nump
  real                 :: vx
  real                 :: vy
  real                 :: vz
  integer              :: data_out
   real                 :: lf
  ! 着氷計算用パラメータ *******************************************
!  real                 :: mvd
!  real,parameter       :: lwc    = 1.2
!  real, parameter      :: pi     = acos(0.0)*2.0.
!  real                 :: AOA
  ! 上記4つは格子法の変数宣言と競合するので、格子法を使う場合はコメントアウト
  real                 :: in_vel
  real                 :: AOAdeg
  integer              :: ini_nump
  integer, allocatable :: icing_col(:)
  real, allocatable    :: ice_an(:)
  integer              :: changep
  real                 :: t_wall
  real                 :: ca
  real                 :: ice_dist_min
  real                 :: ice_an_min

  real                 :: ganma_ini
  ! 着氷量計算パラメータ *******************************************
  INTEGER :: WATER_NUM
  INTEGER :: ICE_NUM
  REAL :: INITIAL_VOL
  REAL :: ICE_VOL
  REAL :: WATER_VOL
  REAL :: VOL_RATE

  ! 時間出力用パラメータ *******************************************
  integer              :: date_time(8)
  character*10         :: date_time_cha(3)

  REAL, PARAMETER      :: vmax_lim     = 0.05
  !REAL, PARAMETER      :: interval_lim = 0.015
  REAL, PARAMETER      :: interval_lim = 0.015
  ! 並列化用 *******************************************************
  real   , allocatable :: x_temp(:,:)
  real   , allocatable :: v_temp(:,:)
  real   , allocatable :: p_temp(:)
  real   , allocatable :: t_temp(:)
  real   , allocatable :: ganma_temp(:)

  ! 格子法用 *******************************************************
  integer, parameter :: mini     =   2         ! 液滴初期ブロック
  integer, parameter :: CalNum   =   3         ! 1 : 衝突　2 : 可視化　3 : 両方
  integer, parameter :: nDrpView = 500         ! 可視化液滴数
  real   , parameter :: SplLim   = 40.0e-6     ! スプラッシュ判定
  ! 共有変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  logical, pointer :: fOver(:, :, :)         ! 重合格子部のフラグ
  integer :: DrpNum                          ! 液滴ファイル番号
  real    :: xmin, ymin, ymax, zmin, zmax    ! 液滴投入範囲
  real    :: uini, vini, wini                ! 液滴初期速度
  real    :: xp, yp, zp                      ! 液滴位置 (物理空間)
  real    :: up, vp, wp                      ! 液滴速度 (物理空間)
  real    :: xip, etp, zep                   ! 液滴位置 (計算空間)
  real    :: uxp, vep, wzp                   ! 液滴速度 (計算空間)
  real    :: dp                              ! 液滴直径
  real    :: fx, fy, fz                      ! 液滴に働く力
  integer :: Nlost                           ! 液滴ロスト数
  logical :: fSearch                         ! 液滴探索成功フラグ
  logical :: fSplash                         ! 液滴スプラッシュフラグ
  logical :: fBounce                         ! 液滴バウンドフラグ
  integer, ALLOCATABLE :: ip(:)
  integer, ALLOCATABLE :: jp(:)
  integer, ALLOCATABLE :: kp(:)
  character :: para_name * (40)

  ! 作成した翼に応じて、適宜翼位置の変更を行う
  real   , parameter :: xpos = 0.0
  real   , parameter :: ypos = 0.0

  ! カップリング用
  integer,parameter :: droplet_nn = 100000
  integer, allocatable :: droplet_num(:)
  integer, allocatable :: droplet_swi(:)
  integer, allocatable :: droplet_mm(:)
  integer, allocatable :: cal_type0(:)
  integer :: mps_swi
  integer :: change_incount
  real, allocatable :: droplet_xaccel(:)
  real, allocatable :: droplet_yaccel(:)
  real, allocatable :: droplet_zaccel(:)
  real, allocatable :: droplet_x(:)
  real, allocatable :: droplet_y(:)
  real, allocatable :: droplet_z(:)
  real, allocatable :: droplet_u(:)
  real, allocatable :: droplet_v(:)
  real, allocatable :: droplet_w(:)
  real :: dp_temp_x
  real :: dp_temp_y
  real :: dp_temp_z
  real :: dx
  real :: dy
  real :: dz
  real :: dp_dt
  integer :: in
  integer, allocatable :: dp_grid_p(:)
  integer, allocatable :: dp_neigh(:,:)
  integer, allocatable :: dp_g_neigh(:,:)
  integer :: swi_lost
  integer :: in_retry
  integer, parameter :: retry_lim = 20
  integer :: swi_cal_mps
  integer :: cal_dp_num

  real, allocatable :: tau(:,:)

  integer,allocatable :: flag_temp(:)

  integer,allocatable :: swi_hc(:)
  integer,allocatable :: swi_solid(:)

  real,dimension(0:130,0:1) :: wall
  real,allocatable	:: wall_dist(:)

  !流れ場リメッシュ用
  real :: time_remesh
  integer :: count_remesh = 0

  real :: finish_time

  !対流熱伝達用
  real	:: fflat
  real,allocatable :: sur_temp(:)

  real  :: real_mvd

  integer :: int_time = 0

  !ICMデータ用
  integer :: num_surf
  real,allocatable :: ICMsurf(:,:)

  integer :: step


  ! 高速化/計算安定用 **********************************************
  ! この値以下の粒子数密度で粒子を取り除く
  real       :: limit_low_an

  ! 計算スイッチ
  ! 初期粒子数密度計算 0:計算, 1:手打ち *基本的に0
  integer,parameter       :: swi_ini_n = 0

  ! 初回の粒子テキスト投入 0:なし, 1:あり
  integer,parameter       :: swi_intxt = 0

  ! 流入のランダム化
  ! 0:指定位置(単一液滴) 1:ランダム(複数液滴) 2:平板着氷 3:液滴滴下
  integer,parameter       :: swi_random = 1

  ! 時間のスライス(凍結)   0:なし, 1:あり
  integer,parameter       :: swi_time_slice = 1

  ! ParaOutput時の時間出力   0:なし, 1:あり
  integer,parameter       :: swi_date_and_time = 1

  ! 着氷計算 0:なし 1:あり 2:0℃でスイッチ
  integer,parameter       :: swi_icing = 2

  ! パラメータ選択 0:単一液滴着氷(NACA) 1:長時間着氷  2:平板着氷 3:液滴滴下
  !                4:格子法カップリング
  integer,parameter       :: swi_parameter = 4

  ! リメッシュ 0:なし 1:あり(0秒から計算) 2:あり(途中から計算)
  integer                 :: swi_remesh

  !対流熱伝達の考慮 0:なし 1:あり
  integer,parameter	  :: swi_heattransfer = 1
  
   ! *******************************************************************
  ! **********************  計算開始  *********************************
  ! *******************************************************************

  open(1,file = './data/STEP.txt',status = 'old')
   read(1,*) step
  close(1)
  if(step .eq. 1) then
   swi_remesh = 1
  else
   swi_remesh = 2
  end if

  ! こいつだけ1からなのに注意
  call date_and_time(date_time_cha(1),date_time_cha(2),date_time_cha(3), &
  &                  date_time)
  write(*,*)
  write(*,*)"**************  CALUCATION START  *********************"
  write(*,'(A,I4,A,I2,A,I2)')" DATE = ",date_time(1),"/",date_time(2),&
  &                          "/",date_time(3)
  write(*,'(A,I2,A,I2,A,I2)')" TIME = ",date_time(5),":",date_time(6),&
  &                          ":",date_time(7)
  write(*,*)



  ! ****************************************************************
  ! *****************  イニシャルデータの読み込み  *****************
  ! ****************************************************************

  ini_nump = 0
  nump = 0
  nextin = 0.0
  incount = 0
  swi_cal_mps = 0

  ! output_setting *************************************************
  open(900,file='./log/MPS/mps_log.txt')
  open(800,file='./data/MPS/random_table.txt')


  write(*,*)"INPUT INITIAL DATA START"

  ! read_initial_file **********************************************
  call input_file                       &
  &    (nn,idim,ifluid,neighbor,type,   &
  &     den,itypep,x,v,an,grav,vis,     &
  &     time_max,ite_max,dt_fix,cfl,    &
  &     dtout,dis,ker_c,                &
  &     delta,sigma,n0,                 &
  &     ite,time_sim,nump,rlambda,      &
  &     grid_n,grid_num,grid,near_wall, &
  &     non_cal,dis_rat,c_vel,          &
  &     f_type,num_ref,p,               &
  &     nn_grid,output_limit,t,         &
  &     dist_min,coll_rat,surf,         &
  &     wall_n,grid_dist,wall_anc,      &
  &     wall_nor)

  ! read_file_for_weight_funtion_calculation ***********************
  call input_ww_wall(ww_num,ww_wall,dis)

  ! read_wall_position *********************************************
  call input_wall(wall)

  ! array_allocation ***********************************************
  call allocation                          &
  &    (idim,nn,neighbor,grid_num,nn_grid, &
  &     visc,accel,press,cal_sw,           &
  &     p_out,grid_p,out,neigh,solid_sur,  &
  &     g_neigh,rho,nu,sol_vel,ker_r)

  allocate(ice_an(0:nn-1))
  allocate(icing_col(0:nn-1))
  allocate(x_temp(0:nn-1,0:idim-1))
  allocate(v_temp(0:nn-1,0:idim-1))
  allocate(p_temp(0:nn-1))
  allocate(t_temp(0:nn-1))
  allocate(tb(0:nn-1))
  allocate(tdd(0:nn-1))
  allocate(Q2(0:nn-1))
  allocate(Q1(0:nn-1))
  allocate(Qa(0:nn-1))
  allocate(Qk(0:nn-1))
  allocate(Qc(0:nn-1))
  allocate(Qd(0:nn-1))
  allocate(Qr(0:nn-1))
  allocate(Mimw(0:nn-1))
  allocate(ganma_temp(0:nn-1))
  allocate(ganma(0:nn-1))
  allocate(swi_koeki(0:nn-1))
  allocate(ip(0:droplet_nn))
  allocate(jp(0:droplet_nn))
  allocate(kp(0:droplet_nn))
  allocate(droplet_num(0:nn-1))
  allocate(droplet_swi(0:droplet_nn))
  allocate(droplet_mm(0:droplet_nn))
  allocate(droplet_xaccel(0:droplet_nn))
  allocate(droplet_yaccel(0:droplet_nn))
  allocate(droplet_zaccel(0:droplet_nn))
  allocate(droplet_x(0:droplet_nn))
  allocate(droplet_y(0:droplet_nn))
  allocate(droplet_z(0:droplet_nn))
  allocate(droplet_u(0:droplet_nn))
  allocate(droplet_v(0:droplet_nn))
  allocate(droplet_w(0:droplet_nn))
  allocate(cal_type0(0:nn-1))
  allocate(dp_g_neigh(1:grid_num,0:nn_grid))
  allocate(dp_grid_p(0:nn-1))
  allocate(dp_neigh(0:nn-1,0:neighbor))
  allocate(Qhc(0:nn-1))

  allocate(flag_temp(0:nn-1))
  flag_temp = 1
  
  allocate(swi_hc(0:nn-1))
  swi_hc = 0
  allocate(swi_solid(0:nn-1))
  swi_solid = -1
  allocate(wall_dist(0:nn-1))
  wall_dist = -1.0

  droplet_num(:)  = 0
  dp_g_neigh(:,:) = 0
  neigh(:,:)      = 0
  Qhc(:)=0


  ! 0:計算なし 1:MPS 2:軌道計算(MPS領域外) 3:軌道計算(MPS領域内)
  droplet_swi(:) = 0

  ! 0:計算なし 1:MPS領域外軌道計算用type(0) 2:MPS領域内
  cal_type0(:) = 0

  ! 液滴初期ブロック
  droplet_mm(:) = mini

  mps_swi = 0

  ! ******************************************************************
  ! ******************   格子法用 初期条件   *************************
  ! ******************************************************************

  ! 検証実験ケース選択 ++++++++++++++++++++++++++++++++++++++++++++++++
  write(*, '(a)'  ) "<< Exp. Case Selection >>"
  call SelectExpCase
  ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  write(*, '(/,a)') "<< Initial Condition Set >>"
  call InitialSetting
   write(*, '(/,a)') "<< CHANGE AOA >>"
  call changeAOAandPOS(xpos,ypos,dis)

  do m =ms,me
    if(m .eq. 1)then
      para_name = "./data/MPS/grid_main"
    else if(m .eq. 2)then
      para_name = "./data/MPS/grid_sub"
    end if

    call Output_Grid3D_para(para_name,  &
    &    Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &    Flw(m)%x, Flw(m)%y, Flw(m)%z)

    call Output_Q3D_para(para_name,  &
    &    Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &    Flw(m)%rho, Flw(m)%u, Flw(m)%v,  Flw(m)%w, Flw(m)%p, &
    &    Flw(m)%t, Flw(m)%mu, Flw(m)%kin, Flw(m)%eps)

    call Output_Function3D_para(para_name,  &
    &    Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &    Flw(m)%rho, Flw(m)%u, Flw(m)%v,  Flw(m)%w, Flw(m)%p, &
    &    Flw(m)%t, Flw(m)%mu, Flw(m)%kin, Flw(m)%eps)





    write(*,*)'---------------------------'
    write(*,*)'m  = ',m
    write(*,*)'is = ',Flw(m)%is
    write(*,*)'ie = ',Flw(m)%ie
    write(*,*)'js = ',Flw(m)%js
    write(*,*)'je = ',Flw(m)%je
    write(*,*)'ks = ',Flw(m)%ks
    write(*,*)'ke = ',Flw(m)%ke
    write(*,*)'---------------------------'

  end do

!   write(*, '(/,a)') "<< Cal Tau >>"

!   m=2
!
!allocate(grid_x(Flw(m)%is:Flw(m)%ie,Flw(m)%js:Flw(m)%je,Flw(m)%ks:Flw(m)%ke))
!allocate(grid_y(Flw(m)%is:Flw(m)%ie,Flw(m)%js:Flw(m)%je,Flw(m)%ks:Flw(m)%ke))
!
!grid_x(:,:,:) = Flw(m)%x(:,:,:)*lRef
!grid_y(:,:,:) = Flw(m)%y(:,:,:)*lRef
!
!    call CalTau( &
!    &    Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
!    &    Flw(m)%rho, tau, idim, grid_x,grid_y, &
!    &           Flw(m)%u,Flw(m)%v)
!
!    allocate(hc(Flw(m)%is:Flw(m)%ie))
!    open(500,file='./HeatTransferCoefficient.txt')
!     do i = Flw(m)%is,Flw(m)%ie
!      read(500,*) hc(i)
!     end do
!    close(500)
!
!    allocate(sur_temp(Flw(m)%is:Flw(m)%ie))
!    open(501,file = 'SurfaceTemperature.txt',status = 'old')
!    do i = Flw(m)%is,Flw(m)%ie
!      read(501,*) sur_temp(i)
!     end do
!    close(501)

  write(*, '(/,a)') "<< Input Grid Flow Data >>"
  open(1,file = './data/GridFlow/GridFlowData.dat',status = 'old')
   read(1,*) num_surf
   allocate(ICMsurf(0:num_surf-1,0:1))
   allocate(tau(0:num_surf-1,0:1))
   allocate(hc(0:num_surf-1))
   allocate(sur_temp(0:num_surf-1))
   do i = 0,num_surf-1
    read(1,*) ICMsurf(i,0),ICMsurf(i,1),tau(i,0),tau(i,1),hc(i),sur_temp(i)
   end do
  close(1)


  ! kernel_radius_setting ******************************************
  do i=0,3
    ker_r(i) = dis*ker_c(i)
  end do

  phase_min = 4
  phase_max = 4

  ! number_of_wall_particle_counting *******************************
  wall_num = 0
  if(swi_remesh .eq. 2)then
  else
    do i=0,nump-1
      if((itypep(i) .eq. type(3)).or. &
      &  (itypep(i) .eq. type(4)).or. &
      &  (itypep(i) .eq. type(5)))then
        wall_num = i
        exit
      end if
    end do
  end if


  ! *****************************************************************
  ! ********************* パラメ―タの強制変更  *********************
  ! *****************************************************************
  if (idim.eq.2)then
    limit_low_an = 0.2
    ice_dist_min = 0.98
    ice_an_min   = 2.2
  else if(idim.eq.3)then
    limit_low_an = 0.5
    ice_dist_min = 0.98
    ice_an_min   = 4.0
  end if

  dt = 1.0e-6           !1.0e-6（元）
  dt_fix = 1.0e-6

  if(swi_parameter.eq.1)then
    output_limit = 999
    DTOUT  = 0.2
    mvd    = dis*10.0
    in_vel = 50.0
    AOAdeg    = 0.0
    AOA = AOAdeg * PI / 180
    t(:)   = 0.0
    t_wall = -10.0
    data_out = 2000
    ca       = 0.0
    Vx = in_vel
    Vy = 0.0
    Vz = 0.0
  else if(swi_parameter.eq.0)then
    output_limit = 999
    DTOUT  = 1.0E-5
    mvd    = dis*10.0
    in_vel = 50.0
    AOAdeg    = 4.0
    AOA = AOAdeg * PI / 180
    t(:)   = 1.0+273.15
    t_wall = -1.0+273.15
    data_out = 100
    ca       = 0.0
    lf=3.334*1.0e5*917.0*(1.0/2.0)*pi*(((0.5*(dis))**2.0))
    Vx = in_vel
    Vy = 0.0
    Vz = 0.0
  else if(swi_parameter.eq.2)then
    DTOUT = 5.0E-6
    OUTPUT_LIMIT = 500
    mvd    = dis*10.0
    in_vel = 50.0
!    AOAdeg    = 10.0
    open(90,file='./AOA.txt')
      read(90,*)ice_an_min
      read(90,*)AOAdeg
      read(90,*)CA
    close(90)
    AOA = AOAdeg * PI / 180
    t(:)   = 1.0
    t_wall = -1.0
    data_out = 100
    Vx = in_vel*cos(AOA)
    Vy = -in_vel*sin(AOA)
    Vz = 0.0
  else if(swi_parameter.eq.3)then
    DTOUT = 5.0E-3
    OUTPUT_LIMIT = 999
    mvd    = dis*10.0
    in_vel = 0.0
    AOAdeg    = 0.0
    AOA = AOAdeg * PI / 180
    t(:)   = 5.0
    t_wall = -10.0
    data_out = 500
    ca       = 120.0
    Vx = in_vel*cos(AOA)
    Vy = -in_vel*sin(AOA)
    Vz = 0.0
  else if(swi_parameter.eq.4)then

    real_mvd = 20.0e-6 !1000.0e-6 !20.0e-6

    ganma_ini = 1.0 - (4.2*1.0e3)*(273.15-(-11.0+273.15))/(3.344*1.0e5)!0.925!(-6度)     !0.937(-5℃のとき)
    ganma(:) = ganma_ini
    !過冷却度．液滴の温度(大気温度)によってここの値は変わる．この算出式は河野さんの修論参照
    swi_koeki(:) = 0

    output_limit = 99999
    DTOUT = 0.1 !1.0E-4 !0.2
    mvd    = dis*10.0
    in_vel = 50.93 !58.09
    AOAdeg    = 0.0 !4.0
    AOA = AOAdeg * PI / 180
    t(:)   = 273.15
    t_wall = 273.15 -11.0 !-8.0    !温度変更するときはここも変更
    data_out = 20000
    ca       = 0.0
    lf=3.334*1.0e5*917.0*(1.0/2.0)*pi*(((0.5*(dis))**2.0))
    Vx = 0.0
    Vy = 0.0
    Vz = 0.0

    lwc = 1.19e-3 !2.1E-3

  else
    write(*,*)'switch error (random)'
    stop
  end if

  if(swi_remesh .ne. 0)then
    time_remesh = 30.0    ![s]    リメッシュ（リスタートファイルを作成）する間隔
  end if

  if(swi_remesh .eq. 1) then
   finish_time = time_remesh
  end if

  write(*,*)'dimension    =',idim
  write(*,*)'ice_dist_min =',ice_dist_min
  write(*,*)'ice_an_min   =',ice_an_min


  time_max = 500.0
  ite_max  = 1000000000

  ! *****************************************************************
  ! ****************  パラメータの強制変更ここまで  *****************
  ! *****************************************************************
  write(*,*)'---------------------------------------------------------'
  write(*,*)'                     PARAMETER                           '
  write(*,*)'dis=',dis
  write(*,*)'mvd=',mvd
  write(*,*)'---------------------------------------------------------'

  ! *****************************************************************
  ! ********************* 初期粒子数密度の計算  *********************
  ! *****************************************************************

  if(swi_ini_n .eq. 0)then
    write(*,*)"START Calculation Initial Particle Number Density"

    ! 2次元計算
    call box_positioning         &
    &    (0.0,20.0*dis,0.0,20.0*dis,0.0,0.0, &
    &     1.0*dis,0.9,type(4),ini_nump,    &
    &     nn,idim,ifluid,x,v,type,itypep, &
    &     dis,grid,grid_num,wall_num,     &
    &     nn_grid,g_neigh,grid_n,         &
    &     f_type,phase_min,phase_max)


    ! 3次元計算
!    call box_positioning               &
!    &    (0.0,0.0,0.0,                  &
!    &    20.0*dis,20.0*dis,20.0*dis,   &
!    &    0.0,0.0,0.0,                  &
!    &    1.0*dis,type(4),ini_nump,    &
!    &    nn,idim,ifluid,x,v,type,itypep, &
!    &    dis,grid,grid_num,wall_num,   &
!    &    nn_grid,g_neigh,grid_n,       &
!    &    f_type,phase_min,phase_max    &
!    &    )

    ! initial_particle_number_setting *********************************
    initial_num = ini_nump

    ! divide_particles_into_grid **************************************
    call grid_division                  &
    &    (nn,nn_grid,idim,ifluid,ini_nump,  &
    &     grid_num,grid_n,grid,         &
    &     g_neigh,x,itypep,type,grid_p)

    ! wall_grid_setting ************************************************
    call wall_grid_setting               &
    &    (idim,grid_num,grid,wall_anc,   &
    &     wall_nor,g_wall,ker_r,nn_wall, &
    &     wall_n)

    ! initial_setting_for_calculation **********************************
    if((rlambda .eq. 0.0) .and. &
    &  (num_ref .eq. 0.0)) then
      ! search_neighboring_particles
      !$omp parallel do default(shared) private(i) schedule(dynamic,64)
      do i=wall_num,ini_nump-1
        if((itypep(i) .ge. phase_min).and. &
        &  (itypep(i) .le. phase_max)) then
          call set_neigh                   &
          &    (i,                         &
          &     nn,idim,ifluid,neighbor,x, &
          &     type,neigh,itypep,         &
          &     ker_r,nn_grid,             &
          &     grid_num,grid_p,g_neigh)
        end if
      end do
      !$omp end parallel do

      ! count_standart_particle_number_density(num_ref)
      do i=wall_num,ini_nump-1
        if((itypep(i) .ge. phase_min).and. &
        &  (itypep(i) .le. phase_max)) then
          call cal_par_num               &
          &    (i,2,ker_r,               &
          &     nn,idim,ifluid,neighbor, &
          &     an,neigh,x,itypep,grid,  &
          &     grid_num,ww_wall,ww_num, &
          &     grid_p,wall_n,wall_anc)
          num_ref = max(an(i),num_ref)
        end if
      end do

      ! count_standart_particle_number_density(n0)
      do i=wall_num,ini_nump-1
        if((itypep(i) .ge. phase_min).and. &
        &  (itypep(i) .le. phase_max)) then
          call cal_n_wall                    &
          &    (i,0,ker_r,                   &
          &     nn,idim,ifluid,neighbor,     &
          &     grid_num,ww_num,grid,        &
          &     ww_wall,wall_n,grid_p,       &
          &     an,neigh,x,itypep,wall_anc)
          n0(0,0) = max(an(i),n0(0,0))
        end if
      end do

      ! count_standart_particle_number_density(n0)
      do i=wall_num,ini_nump-1
        if((itypep(i) .ge. phase_min).and. &
        &  (itypep(i) .le. phase_max)) then
          call cal_n_wall                    &
          &    (i,2,ker_r,                   &
          &     nn,idim,ifluid,neighbor,     &
          &     grid_num,ww_num,grid,        &
          &     ww_wall,wall_n,grid_p,       &
          &     an,neigh,x,itypep,wall_anc)
          n0(0,2) = max(an(i),n0(0,2))
        end if
      end do

      ! count_standart_particle_number_density(n0)
      do i=wall_num,ini_nump-1
        if((itypep(i) .ge. phase_min).and. &
        &  (itypep(i) .le. phase_max)) then
          call cal_n_wall                    &
          &    (i,1,ker_r,                   &
          &     nn,idim,ifluid,neighbor,     &
          &     grid_num,ww_num,grid,        &
          &     ww_wall,wall_n,grid_p,       &
          &     an,neigh,x,itypep,wall_anc)
          n0(0,1) = max(an(i),n0(0,1))
        end if
      end do

      ! data_copy
      do i=0,nn-1
        n0(i,0:3) = n0(0,0:3)
      end do

      ! count_standart_particle_number_density(rlambda)
      rlambda = -1.0e10
      do i=wall_num,ini_nump-1
        call cal_lambda                     &
        &    (i,1,ker_r,rlambda_ref,        &
        &     nn,idim,ifluid,neighbor,      &
        &     grid_num,ww_num,grid,ww_wall, &
        &     wall_n,grid_p,neigh,x,n0,wall_anc)
        rlambda = max(rlambda,rlambda_ref)
      end do
    end if

    ! rho, nu_setting **************************************************
!   !$omp parallel do default(shared) private(i) schedule(dynamic,64)
!   do i=0,nump-1
!     call set_rhonu                   &
!     &    (i,nn,idim,ifluid,neighbor,   &
!     &     itypep,den,vis,rho,nu,nump)
!   end do
!   !$omp end parallel do
! [call set_rhonu] is manually inlined - begin
    do i=0,nump-1
      nu(i)  = vis(itypep(i))
      rho(i) = den(itypep(i))
    end do
! [call set_rhonu] is manually inlined - end

!    call output_file                   &
!    &    (idim,nn,ini_nump,ite,time_sim,   &
!    &     type,ifluid,itypep,          &
!    &     den,vis,grav,delta,sigma,    &
!    &     time_max,ite_max,            &
!    &     dt,dtout,cfl,                &
!    &     dis,ker_c,n0,                &
!    &     x,non_cal,num_ref,rlambda,   &
!    &     grid,grid_n,grid_num,        &
!    &     dis_rat,c_vel,               &
!    &     f_type,output_limit,         &
!    &     nn_grid,neighbor,t,wall_nor, &
!    &     dist_min,coll_rat,near_wall, &
!    &     wall_n,grid_dist,wall_anc)
!
!    ! output :vtk_file *************************************************
!    call output_para_bin            &
!    &    (nn,idim,ifluid,fcount,t,  &
!    &     x,v,itypep,type,an,p,surf,cal_type0)

    do i=wall_num,ini_nump-1
      itypep(i) = type(0)
      an(i)     = 0.0
    end do

    call calculation_fi_flat(idim,dis,fflat)    !|F_i_flat|の計算
    write(*,*) 'f_flat :',fflat

    write(*,*)"END Calculation Initial Particle Number Density"
    write(*,*) n0(0,0),n0(0,1),n0(0,2),n0(0,3),rlambda,num_ref
    write(*,*)'---------------------------------------------------------'

  else if(swi_ini_n .eq.1)then

    ! ****************************************************************
    ! ********************手打ち用************************************
    ! ****************************************************************

    write(*,*)"INPUT Initial Particle Number Density"


    ! divide_particles_into_grid **************************************
    call grid_division                  &
    &    (nn,nn_grid,idim,ifluid,nump,  &
    &     grid_num,grid_n,grid,         &
    &     g_neigh,x,itypep,type,grid_p)

    ! wall_grid_setting ************************************************
    call wall_grid_setting               &
    &    (idim,grid_num,grid,wall_anc,   &
    &     wall_nor,g_wall,ker_r,nn_wall, &
    &     wall_n)


    if(idim .eq. 2)then
      n0(0,0) = 9.963232
      n0(0,1) = 9.963232
      n0(0,2) = 27.98962
      n0(0,3) = 0.0
      rlambda = 1.8682274E-08
      num_ref = 29.00000
    else if(idim .eq. 3)then
      n0(0,0) = 21.71523
      n0(0,1) = 21.71523
      n0(0,2) = 65.82861
      n0(0,3) = 0.0
      rlambda = 2.1100842E-08
      num_ref = 142.0000
    end if

    ! data_copy
    do i=0,nn-1
      n0(i,0:3) = n0(0,0:3)
    end do

    write(*,*) n0(0,0:3),rlambda
    write(*,*)'---------------------------------------------------------'

    ! rho, nu_setting **************************************************
!   !$omp parallel do default(shared) private(i) schedule(dynamic,64)
!   do i=0,nump-1
!     call set_rhonu                   &
!     &    (i,nn,idim,ifluid,neighbor,   &
!     &     itypep,den,vis,rho,nu,nump)
!   end do
!   !$omp end parallel do
    do i=0,nump-1
      nu(i)  = vis(itypep(i))
      rho(i) = den(itypep(i))
    end do

  end if

  call dom_search                  &
  &    (nn,nn_grid,idim,ifluid,nump,  &
  &     grid_num,grid_n,grid,         &
  &     g_neigh,x,itypep,type,grid_p, &
  &     dom_xmax,dom_xmin,dom_ymax,dom_ymin,dom_zmax,dom_zmin)

  if(idim .eq. 2)then
    dom_zmin = - dis
    dom_zmax = + dis
  end if

  ! ********************************************************************
  ! *************************  初期流入  *******************************
  ! ********************************************************************

  if((swi_remesh .eq. 0) .or. (swi_remesh .eq. 1))then
    write(*,*)'---------------------------------------------------------'
    write(*,*)"INITIAL FLOW SETTING START"
    
    in_retry = 0
    do
      if(swi_random .eq. 1)then
        ! divide_particles_into_grid **************************************
        
        RANDOMX = -0.1 * chord * lRef !-0.053
        
        ! 乱数テーブルの使用
        ! 乱数による投入位置決定で余計なことをしている感がある
        ! 下のコメントアウトしてるので十分
        
        !call random_number(random1)
        !read(800,*)random1
        !random1 = random1 * 10e8
        !random2 = nint(random1)
        !random2 = mod(random2,100000)
        !randomy = dom_ymin + mvd*0.5 + ((dom_ymax - dom_ymin) - mvd*0.5) * random2 * 0.00001
        
        
        !call random_number(random1)
        !read(800,*)random1
        !random1 = random1 * 10e8
        !random2 = nint(random1)
        !random2 = mod(random2,100000)
        !randomz = dom_zmin + mvd*0.5 + ((dom_zmax - dom_zmin) - mvd*0.5) * random2 * 0.00001
        
        read(800,*)random1
        randomy = dom_ymin + mvd*0.5 + ((dom_ymax - dom_ymin) - mvd*0.5) * real(random1)
      !  read(800,*)random1
      !  randomz = dom_zmin + mvd*0.5 + ((dom_zmax - dom_zmin) - mvd*0.5) * real(random1)
        
      else if(swi_random .eq. 0)then
        
    !    RANDOMX = 0.027
    !    RANDOMY = 0.031
    !    RANDOMZ = 0.0025
        
    !    RANDOMX = 0.023
    !    RANDOMY = 0.029
    !    RANDOMZ = 0.0025
        
    !    RANDOMX = 0.005
    !    RANDOMY = 0.0255
    !    RANDOMZ = 0.0025
        
        RANDOMX = 0.01
        RANDOMY = 0.037
        RANDOMZ = 0.0025
        
    !    RANDOMX = 0.040
    !    RANDOMY = 0.052
    !    RANDOMZ = 0.01
        
      else if(swi_random .eq. 2)then
        
        RANDOMX = 0.0015
        RANDOMY = 0.0015
        RANDOMZ = 0.005
        
    !    RANDOMX = 0.005
    !    RANDOMY = 0.005
    !    RANDOMZ = 0.005
    !    Vx = 0.0
    !    vy =0.0
    !    Vz = 0.0
    !    grav = 0.0
    !!    DTOUT = 5.0E-6
    !    DTOUT = 500.0E-6
        
      else if(swi_random .eq. 3)then
        
        RANDOMX = 0.005
        RANDOMY = 0.003
        RANDOMZ = 0.005
      
      end if
  !    randomy = 0.037 + MVD*20.0 - MVD * real(incount+1)*4.0
      !randomy = 0.037
      
      incount = incount + 1
      droplet_swi(incount) = 2
      in = incount
      
      if(idim .eq. 2)then
        droplet_x(in) = randomx
        droplet_y(in) = randomy
        droplet_z(in) = 0.0
      else if(idim .eq. 3)then
        droplet_x(in) = randomx
        droplet_y(in) = randomy
        droplet_z(in) = randomz
      end if
      
      droplet_u(in) = Vx
      droplet_v(in) = Vy
      droplet_w(in) = Vz
      
      swi_lost = 0
        
      ! 無次元化
      dp = MVD / lRef
      xp = droplet_x(in) / lRef
      yp = droplet_y(in) / lRef
      zp = droplet_z(in) / lRef
      up = droplet_u(in) / aRef
      vp = droplet_v(in) / aRef
      wp = droplet_w(in) / aRef
      
      call ini_DropletTimeStep(ip(in),jp(in),kp(in),droplet_mm(in))
      
      ! 有次元化
      droplet_x(in) = xp * lRef
      droplet_y(in) = yp * lRef
      droplet_z(in) = zp * lRef
      droplet_u(in) = up * aRef
      droplet_v(in) = vp * aRef
      droplet_w(in) = wp * aRef
     
      if(swi_lost .eq. 1)then
        incount = incount - 1
      else if(swi_lost .eq. 0)then
        exit
      end if
      
      in_retry = in_retry + 1
      if(in_retry .ge. retry_lim)then
        write(*,*)"error inflow position"
        stop
      end if
    end do
    
    temp_nump = nump
    
    if(swi_intxt .eq. 0)then
      
      ! 2次元計算
         call circle_positioning         &
         &    (randomx,randomy,mvd,Vx,Vy, &
         &     0.25*dis,2.0,type(0),nump,      &
         &     nn,idim,ifluid,x,v,type,itypep, &
         &     dis,grid,grid_num,wall_num,     &
         &     nn_grid,g_neigh,grid_n,         &
         &     f_type,phase_min,phase_max,inptnum, &
         &     incount,droplet_num,cal_type0)
         
         call write_circle_position(   &
        &    nn,idim,dis,mvd,x,nump,temp_nump, &
        &    randomx,randomy,randomz           &
        &    )
      
      ! 3次元計算
  !      call sphere_positioning             &
  !      &    (randomx,randomy,randomz,mvd,     &
  !      &     Vx,Vy,Vz, &
  !      &     0.25*dis,2.0,type(4),nump,  &
  !      &     nn,idim,ifluid,x,v,type,itypep, &
  !      &     dis,grid,grid_num,wall_num,   &
  !      &     nn_grid,g_neigh,grid_n,       &
  !      &     f_type,phase_min,phase_max,inptnum)
  !      call write_sphere_position(   &
  !      &    nn,idim,dis,mvd,x,nump,temp_nump, &
  !      &    randomx,randomy,randomz           &
  !      &    )
      
    else if(swi_intxt .eq. 1)then
      
      ! 2次元
      call circle_positioning_txt         &
      &    (randomx,randomy,mvd,Vx,Vy, &
      &     0.25*dis,2.0,type(0),nump,      &
      &     nn,idim,ifluid,x,v,type,itypep, &
      &     dis,grid,grid_num,wall_num,     &
      &     nn_grid,g_neigh,grid_n,         &
      &     f_type,phase_min,phase_max,inptnum, &
      &     incount,droplet_num,cal_type0,swi_solid)
      
      ! 3次元
  !      call sphere_positioning_txt             &
  !      &    (randomx,randomy,randomz,mvd,     &
  !      &     Vx,Vy,Vz, &
  !      &     0.25*dis,2.0,type(4),nump,  &
  !      &     nn,idim,ifluid,x,v,type,itypep, &
  !      &     dis,grid,grid_num,wall_num,   &
  !      &     nn_grid,g_neigh,grid_n,       &
  !      &     f_type,phase_min,phase_max,inptnum)
      
    end if
    
    if(swi_random .eq. 1)then
    
      write(*,*)'---------------------------------------------------------'
      WRITE(*,'(A,I4)')" inflow droplet : number = ",incount
      WRITE(*,*)"x = ",randomx
      WRITE(*,*)"y = ",randomy
      WRITE(*,*)"z = ",randomz
      write(*,*)
          WRITE(*,*)"Y_MAX = ",dom_ymax
      WRITE(*,*)"Y_MIN = ",dom_ymin
      WRITE(*,*)"LWC = ",lwc
     WRITE(*,*)"z = ",((ABS(dom_ymax - dom_ymin) - MVD) * Vx) * 1000
      write(*,*)
      
      if(idim .eq. 2)then
        DTIN    = DEN(4) / LWC * PI * (DIS * 0.5)**2 * INPTNUM / &
        &         ((ABS(dom_ymax - dom_ymin) - MVD) * in_vel)
        !DTIN    = DEN(4) / LWC * PI * (DIS * 0.5)**3 *4.0 / 3.0 * 560 / &
        !&         ((ABS(dom_ymax - dom_ymin) - MVD) * MVD * in_vel)
      else if(idim .eq. 3)then
        dtin    = den(4) / lwc * pi *4.0/3.0 * (dis * 0.5)**3 * inptnum / &
        &         ((abs(dom_ymax - dom_ymin) - mvd) * (abs(dom_zmax - dom_zmin) - mvd) * &
        &         in_vel) * 1000
      end if
      
      nextin = nextin + dtin
      
      WRITE(*,*)'INPTNUM  = ', INPTNUM
      WRITE(*,*)'DTIN     = ', DTIN
      WRITE(*,*)'NEXTIN   = ', NEXTIN
      write(*,*)'---------------------------------------------------------'
      
    end if
  
  
    ! *******************************************************************
    ! *********************  初期流入終わり  ****************************
    ! *******************************************************************
    
    ! initial_particle_number_setting *********************************
    initial_num = nump
    
    ! output_initial_file **********************************************
    call output_file                   &
    &    (idim,nn,nump,ite,time_sim,   &
    &     type,ifluid,itypep,          &
    &     den,vis,grav,delta,sigma,    &
    &     time_max,ite_max,            &
    &     dt,dtout,cfl,                &
    &     dis,ker_c,n0,                &
    &     x,non_cal,num_ref,rlambda,   &
    &     grid,grid_n,grid_num,        &
    &     dis_rat,c_vel,               &
    &     f_type,output_limit,         &
    &     nn_grid,neighbor,t,wall_nor, &
    &     dist_min,coll_rat,near_wall, &
    &     wall_n,grid_dist,wall_anc)
  
    ! output :vtk_file *************************************************
!    call output_para_bin            &
!    &    (nn,idim,ifluid,fcount,t,  &
!    &     x,v,itypep,type,an,p,surf,cal_type0,swi_koeki,ganma)
    call output_para_bin            &
    &    (nn,idim,ifluid,fcount,t,  &
    &     x,v,itypep,type,an,p,surf,cal_type0,swi_koeki,ganma, &
    &     swi_solid,swi_hc)
  
    ! count_up
    fcount = fcount+1
  
  else if(swi_remesh .eq. 2)then
    write(*,*) '---------------------------------------------------------'
    write(*,*) 'START READING RESTART FILE'

    call read_restart_file         &                    !リスタートするためのファイルの読み込み
&           (nump,droplet_nn,time_sim_dble,    &
&           nn,idim,ifluid,x,v,type,itypep, &
&           dis,grid,grid_num,wall_num,     &
&           nn_grid,g_neigh,grid_n,         &
&           f_type,phase_min,phase_max,inptnum, &
&           incount,droplet_num,cal_type0,  &
&           t,ganma,swi_koeki,time_sim,nextin, &
&           fcount,count_remesh, &
&           droplet_swi,droplet_x,droplet_y,droplet_z, &
&           droplet_u,droplet_v,droplet_w, &
&           flag_temp)

  count_remesh = count_remesh + 1
  swi_cal_mps = 1
!  droplet_swi(incount) = 3

  finish_time = time_sim + time_remesh

  ! particle_search **********************************************
  call grid_division                  &
  &    (nn,nn_grid,idim,ifluid,nump,  &
  &     grid_num,grid_n,grid,         &
  &     g_neigh,x,itypep,type,grid_p)


  !$omp parallel do default(shared) private(j) schedule(dynamic,64)
  do j=0,nump-1
    if((itypep(j) .ge. phase_min).and. &
    &  (itypep(j) .le. phase_max))then
      cal_sw(j,0:1) = 1
    else
      cal_sw(j,0:1) = 0
    end if
 
    if(itypep(j) .eq. type(3))then
      if(solid_sur(j) .eq. 0)then
        cal_sw(j,0) = 0
      end if
    end if
  end do
  !$omp end parallel do

  
  !$omp parallel do default(shared) private(i) schedule(dynamic,64)
  do i=wall_num,nump-1
    if(cal_sw(i,1) .eq. 1)then
      call set_neigh                   &
      &    (i,                         &
      &     nn,idim,ifluid,neighbor,x, &
      &     type,neigh,itypep,         &
      &     ker_r,nn_grid,             &
      &     grid_num,grid_p,g_neigh)
    end if
  end do
  !$omp end parallel do

    write(*,*) 'FINISH READY TO RESTART'
    write(*,*) 'Restart time  =',time_sim
    write(*,*) '---------------------------------------------------------'
    write(*,*)

  else
    write(*,*) 'Switch Error swi_remesh'
    write(*,*) 'Swi_remesh   =',swi_remesh
    stop
  end if


  time_sim_dble = dble(time_sim)

  write(*,*)"INITIAL FLOW SETTING COMPLETED"
  write(*,*)
  
  phase_min = 4
  phase_max = 4
  phase = type(4)

  ! *******************************************************************
  ! *********************  メインループ開始  **************************
  ! *******************************************************************
  write(*,*)"**********************************************************"
  write(*,*)"****************  MAIN LOOP START  ***********************"
  write(*,*)"**********************************************************"

  ! main_loop_start **************************************************
  1000 ite= ite+1

  ! calculation_dt ***************************************************

!if(time_sim .gt. real(int_time)*0.001) then
! write(*,*) time_sim_dble
! int_time = int_time + 1
!end if

  call cal_dt                          &
  &    (nn,idim,ifluid,neighbor,       &
  &     neigh,nump,itypep,type,cfl,x,  &
  &     dt_fix,dt,time_sim_dble,c_vel,dis,v,vmax)

  time_sim = real(time_sim_dble)

  interval = interval + dt

!  if ((swi_time_slice .eq. 1) .and. &
!  &   (swi_cal_mps .eq. 1))then
!    if((vmax .lt. vmax_lim) .or. &
!    &  (interval .ge. interval_lim))then
!      ! 格子法とのカップリングなので、時間のスライスではなく
!      ! MPS法パートの時間の凍結

!      swi_cal_mps = 0
!      write(*,*)"mps computation freeze"
!      write(*,*)"time = ",time_sim

!      ! 残っている液体粒子を固体粒子に
!      do i = 0, nump-1
!        if(itypep(i) .eq. type(4))then
!          itypep(i) = type(2)
!        end if
!      end do

!    end if
!  end if


  if((time_sim .ge. nextin) .and. &
  &  (swi_random .eq. 1))then

    in_retry = 0
    do
      RANDOMX = -0.1 * chord * lRef !-0.053

      ! 乱数テーブルの使用
      ! 乱数による投入位置決定で余計なことをしている感がある
      ! 下のコメントアウトしてるので十分

      !call random_number(random1)
      read(800,*)random1
      !random1 = random1 * 10e8
      !random2 = nint(random1)
      !random2 = mod(random2,100000)
      !randomy = dom_ymin + mvd*0.5 + ((dom_ymax - dom_ymin) - mvd*0.5) * random2 * 0.00001


      !call random_number(random1)
      !read(800,*)random1
      !random1 = random1 * 10e8
      !random2 = nint(random1)
      !random2 = mod(random2,100000)
      !randomz = dom_zmin + mvd*0.5 + ((dom_zmax - dom_zmin) - mvd*0.5) * random2 * 0.00001

      read(800,*)random1
      randomy = dom_ymin + mvd*0.5 + ((dom_ymax - dom_ymin) - mvd*0.5) * real(random1)
  !    read(800,*)random1
  !    randomz = dom_zmin + mvd*0.5 + ((dom_zmax - dom_zmin) - mvd*0.5) * real(random1)

      incount = incount + 1
      droplet_swi(incount) = 2
      in = incount

      if(idim .eq. 2)then
        droplet_x(in) = randomx
        droplet_y(in) = randomy
        droplet_z(in) = 0.0
      else if(idim .eq. 3)then
        droplet_x(in) = randomx
        droplet_y(in) = randomy
        droplet_z(in) = randomz
      end if

      droplet_u(in) = Vx
      droplet_v(in) = Vy
      droplet_w(in) = Vz

      swi_lost = 0

      ! 無次元化
      dp = MVD / lRef
      xp = droplet_x(in) / lRef
      yp = droplet_y(in) / lRef
      zp = droplet_z(in) / lRef
      up = droplet_u(in) / aRef
      vp = droplet_v(in) / aRef
      wp = droplet_w(in) / aRef

      call ini_DropletTimeStep(ip(in),jp(in),kp(in),droplet_mm(in))

      ! 有次元化
      droplet_x(in) = xp * lRef
      droplet_y(in) = yp * lRef
      droplet_z(in) = zp * lRef
      droplet_u(in) = up * aRef
      droplet_v(in) = vp * aRef
      droplet_w(in) = wp * aRef

      if(swi_lost .eq. 1)then
        incount = incount - 1
      else if(swi_lost .eq. 0)then
        exit
      end if

      in_retry = in_retry + 1
      if(in_retry .ge. retry_lim)then
        write(*,*)"error inflow position"
        stop
      end if
    end do
      ! 2次元
      call circle_positioning_txt         &
      &    (randomx,randomy,mvd,Vx,Vy, &
      &     0.25*dis,2.0,type(0),nump,      &
      &     nn,idim,ifluid,x,v,type,itypep, &
      &     dis,grid,grid_num,wall_num,     &
      &     nn_grid,g_neigh,grid_n,         &
      &     f_type,phase_min,phase_max,inptnum, &
      &    incount,droplet_num,cal_type0,swi_solid)

    ! 3次元
!      call sphere_positioning_txt             &
!      &    (randomx,randomy,randomz,mvd,     &
!      &     Vx,Vy,Vz, &
!      &     0.25*dis,2.0,type(4),nump,  &
!      &     nn,idim,ifluid,x,v,type,itypep, &
!      &     dis,grid,grid_num,wall_num,   &
!      &     nn_grid,g_neigh,grid_n,       &
!      &     f_type,phase_min,phase_max,inptnum)



    if(idim .eq. 2)then
      DTIN    = DEN(4) / LWC * PI * (DIS * 0.5)**2 * INPTNUM / &
      &         ((ABS(dom_ymax - dom_ymin) - MVD) * in_vel)
      !DTIN    = DEN(4) / LWC * PI * (DIS * 0.5)**3 *4.0 / 3.0 * 560 / &
      !&         ((ABS(dom_ymax - dom_ymin) - MVD) * MVD * in_vel)
    else if(idim .eq. 3)then
      dtin    = den(4) / lwc * pi *4/3 * (dis * 0.5)**3 * inptnum / &
      &         ((abs(dom_ymax - dom_ymin) - mvd) * (abs(dom_zmax - dom_zmin) - mvd) * &
      &         in_vel) * 1000
    end if

    nextin = nextin + dtin

    write(*,*)'---------------------------------------------------------'
    WRITE(*,'(A,I4)')"inflow doroplet : number = ",incount
    WRITE(*,*)"x = ",randomx
    WRITE(*,*)"y = ",randomy
    WRITE(*,*)"z = ",randomz
    WRITE(*,*)"x = ",dom_ymax
    WRITE(*,*)"y = ",dom_ymin
    WRITE(*,*)"z = ",MVD
    WRITE(*,*)"z = ",lwc
   WRITE(*,*)"z = ",((ABS(dom_ymax - dom_ymin) - MVD) * in_vel) * 1000
    write(*,*)
    write(*,*)'inptnum = ', inptnum
    write(*,*)'dtin = ', dtin
    write(*,*)'nextin = ', nextin
    write(*,*)'---------------------------------------------------------'


    do i=0,nump-1
      if((itypep(i) .ge. phase_min).and. &
      &  (itypep(i) .le. phase_max))then
        cal_sw(i,0:1) = 1
      else
        cal_sw(i,0:1) = 0
      end if

      if(itypep(i) .eq. type(3))then
        if(solid_sur(i) .eq. 0)then
          cal_sw(i,0) = 0
        end if
      end if
    end do

  end if

!!ここ絶対いらない
!  do
!  !if(fcount .ge. 100*incount)then
!  if(10 .ge. incount)then
!
!  in_retry = 0
!
!  do
!
!    randomy = 0.037 + MVD*20.0 - MVD * real(incount+1)*4.0
!
!    incount = incount + 1
!    droplet_swi(incount) = 2
!    in = incount
!
!    if(idim .eq. 2)then
!      droplet_x(in) = randomx
!      droplet_y(in) = randomy
!      droplet_z(in) = 0.0
!    else if(idim .eq. 3)then
!      droplet_x(in) = randomx
!      droplet_y(in) = randomy
!      droplet_z(in) = randomz
!    end if
!
!    droplet_u(in) = Vx
!    droplet_v(in) = Vy
!    droplet_w(in) = Vz
!
!    swi_lost = 0
!
!    ! 無次元化
!    dp = MVD / lRef
!    xp = droplet_x(in) / lRef
!    yp = droplet_y(in) / lRef
!    zp = droplet_z(in) / lRef
!    up = droplet_u(in) / aRef
!    vp = droplet_v(in) / aRef
!    wp = droplet_w(in) / aRef
!
!    call ini_DropletTimeStep(ip(in),jp(in),kp(in))
!
!    ! 有次元化
!    droplet_x(in) = xp * lRef
!    droplet_y(in) = yp * lRef
!    droplet_z(in) = zp * lRef
!    droplet_u(in) = up * aRef
!    droplet_v(in) = vp * aRef
!    droplet_w(in) = wp * aRef
!
!    if(swi_lost .eq. 1)then
!      incount = incount - 1
!    else if(swi_lost .eq. 0)then
!      exit
!    end if
!    in_retry = in_retry + 1
!    if(in_retry .ge. retry_lim)then
!      write(*,*)"error inflow position"
!      stop
!    end if
!
!  end do
!
!  temp_nump = nump
!
!  if(swi_intxt .eq. 0)then
!
!    ! 2次元計算
!       call circle_positioning         &
!       &    (randomx,randomy,mvd,Vx,Vy, &
!       &     0.25*dis,2.0,type(0),nump,      &
!       &     nn,idim,ifluid,x,v,type,itypep, &
!       &     dis,grid,grid_num,wall_num,     &
!       &     nn_grid,g_neigh,grid_n,         &
!       &     f_type,phase_min,phase_max,inptnum, &
!       &     incount,droplet_num,cal_type0)
!
!       call write_circle_position(   &
!      &    nn,idim,dis,mvd,x,nump,temp_nump, &
!      &    randomx,randomy,randomz           &
!      &    )
!
!    ! 3次元計算
!!      call sphere_positioning             &
!!      &    (randomx,randomy,randomz,mvd,     &
!!      &     Vx,Vy,Vz, &
!!      &     0.25*dis,2.0,type(4),nump,  &
!!      &     nn,idim,ifluid,x,v,type,itypep, &
!!      &     dis,grid,grid_num,wall_num,   &
!!      &     nn_grid,g_neigh,grid_n,       &
!!      &     f_type,phase_min,phase_max,inptnum)
!!      call write_sphere_position(   &
!!      &    nn,idim,dis,mvd,x,nump,temp_nump, &
!!      &    randomx,randomy,randomz           &
!!      &    )
!
!  else if(swi_intxt .eq. 1)then
!
!    ! 2次元
!    call circle_positioning_txt         &
!    &    (randomx,randomy,mvd,Vx,Vy, &
!    &     0.25*dis,2.0,type(0),nump,      &
!    &     nn,idim,ifluid,x,v,type,itypep, &
!    &     dis,grid,grid_num,wall_num,     &
!    &     nn_grid,g_neigh,grid_n,         &
!    &     f_type,phase_min,phase_max,inptnum, &
!    &     incount,droplet_num,cal_type0)
!
!    ! 3次元
!!      call sphere_positioning_txt             &
!!      &    (randomx,randomy,randomz,mvd,     &
!!      &     Vx,Vy,Vz, &
!!      &     0.25*dis,2.0,type(4),nump,  &
!!      &     nn,idim,ifluid,x,v,type,itypep, &
!!      &     dis,grid,grid_num,wall_num,   &
!!      &     nn_grid,g_neigh,grid_n,       &
!!      &     f_type,phase_min,phase_max,inptnum)
!
!  end if
!
!  else
!    exit
!  end if
!  end do

  ! **********************  格子法粒子属性チェック  ******************

  do i = 1,incount
    if(droplet_swi(i) .eq. 2)then
      if((dom_XMIN .LE. droplet_x(i)) .AND. (droplet_x(i) .LE. dom_XMAX) .AND. &
      &  (dom_YMIN .LE. droplet_y(i)) .AND. (droplet_y(i) .LE. dom_YMAX) .and. &
      &  (dom_ZMIN .LE. droplet_z(i)) .AND. (droplet_z(i) .LE. dom_ZMAX))then
        do j = 0,nump-1
          if(droplet_num(j) .eq. i)then
            cal_type0(j)   = 2
            droplet_swi(i) = 3
          end if
        end do
      end if
    else if(droplet_swi(i) .eq. 3)then
      if((dom_XMIN .LE. droplet_x(i)) .AND. (droplet_x(i) .LE. dom_XMAX) .AND. &
      &  (dom_YMIN .LE. droplet_y(i)) .AND. (droplet_y(i) .LE. dom_YMAX) .and. &
      &  (dom_ZMIN .LE. droplet_z(i)) .AND. (droplet_z(i) .LE. dom_ZMAX))then
      else
        do j = 0,nump-1
          if(droplet_num(j) .eq. i)then
            cal_type0(j)   = 1
            droplet_swi(i) = 2
          end if
        end do
      end if
    end if
  end do

  ! ******************************************************************
  ! *******************  格子法 計算パート開始  **********************
  ! ******************************************************************
  do i = 1,incount
    if(droplet_swi(i) .ge. 2)then

      swi_lost = 0

      dp_temp_x = droplet_x(i)
      dp_temp_y = droplet_y(i)
      dp_temp_z = droplet_z(i)

      ! 無次元化
      dp = MVD / lRef !real_mvd / lRef !MVD / lRef
      dp_dt = dt / (lRef/aRef)

      xp = droplet_x(i) / lRef
      yp = droplet_y(i) / lRef
      zp = droplet_z(i) / lRef
      up = droplet_u(i) / aRef
      vp = droplet_v(i) / aRef
      wp = droplet_w(i) / aRef

      call DropletTimeStep(ip(i),jp(i),kp(i),droplet_mm(in))

      ! 有次元化
      droplet_x(i) = xp * lRef
      droplet_y(i) = yp * lRef
      droplet_z(i) = zp * lRef
      droplet_u(i) = up * aRef
      droplet_v(i) = vp * aRef
      droplet_w(i) = wp * aRef

      if(swi_lost .eq. 1) then

        droplet_swi(i) = 0
  !$omp parallel default(shared) private(j)
  !$omp do schedule(dynamic,64)
        do j = 0,nump-1
          if((itypep(j) .eq. type(0)) .and. &
          &  (cal_type0(j) .ge. 1)    .and. &
          &  (droplet_num(j) .eq. i))then
            cal_type0(j) = 0
            droplet_num(j) = 0

            if(idim .eq. 2)then
              v(j,0) = 0.0
              v(j,1) = 0.0
            else if(idim .eq. 3)then
              v(j,0) = 0.0
              v(j,1) = 0.0
              v(j,2) = 0.0
            end if
          end if
        end do
  !$omp end do
  !$omp end parallel

      else if(swi_lost .eq. 0) then

        dx = droplet_x(i) - dp_temp_x
        dy = droplet_y(i) - dp_temp_y
        dz = droplet_z(i) - dp_temp_z
  !$omp parallel default(shared) private(j)
  !$omp do schedule(dynamic,64)
        do j = 0,nump-1
          if((itypep(j) .eq. type(0)) .and. &
          &  (cal_type0(j) .ge. 1)    .and. &
          &  (droplet_num(j) .eq. i))then
            if(idim .eq. 2)then
              x(j,0) = x(j,0) + dx
              x(j,1) = x(j,1) + dy
            else if(idim .eq. 3)then
              x(j,0) = x(j,0) + dx
              x(j,1) = x(j,1) + dy
              x(j,2) = x(j,2) + dz
            end if
          end if
        end do
  !$omp end do
  !$omp end parallel
      end if
    end if
  end do

  ! ******************************************************************
  ! *******************  格子法 計算パート終了  **********************
  ! ******************************************************************

  if(swi_cal_mps .eq. 1)then

  ! ******************************************************************
  ! *********************  MPS 計算パート開始  ***********************
  ! ******************************************************************

  ! calculation_surface_tension_term *********************************
!  call surface_tension_pot        &
!  &    (nn,idim,ifluid,neighbor,  &
!  &     neigh,nump,itypep,type,x, &
!  &     dis,den,surf,sigma,       &
!  &     wall_num,ker_r,an,n0,     &
!  &     grid_num,ww_num,grid,     &
!  &     ww_wall,wall_n,grid_p,    &
!  &     grid_dist,wall_anc,       &
!  &     near_wall,wall_nor,       &
!  &     nn_wall,g_wall,CA)

  !$omp parallel default(shared) private(i) 
  !$omp do schedule(dynamic,64)
  do j=0,idim-1
   do i=0,nump-1
    surf(i,j) = 0.0
   end do
  end do
  !$omp end do
  !$omp do schedule(dynamic,64)
  do i=0,nump-1
   swi_hc(i) = 0
  end do  
  !$omp end do
  !$omp end parallel

  !$omp parallel default(shared) private(i) 
  !$omp do schedule(dynamic,64)
  do i=wall_num,nump-1
   call surface_tension_pot2          &
   &    (i,nn,idim,ifluid,neighbor,  &
   &     neigh,nump,itypep,type,x,   &
   &     dis,den,surf,sigma,         &
   &     wall_num,ker_r,an,n0,       &
   &     grid_num,ww_num,grid,       &
   &     ww_wall,wall_n,grid_p,      &
   &     grid_dist,wall_anc,         &
   &     near_wall,wall_nor,         &
   &     nn_wall,g_wall,CA,          &
   &     swi_solid,swi_hc,fflat)
   if(flag_temp(i) .eq. 0) swi_hc(i) = 0
  end do
  !$omp end do
  !$omp end parallel

  !$omp parallel default(shared) private(i) 
  !$omp do schedule(dynamic,64)
  ! data copy for viscosity
  do j=0,idim-1
   do i=0,nump-1
    x_temp(i,j) = x(i,j)
    v_temp(i,j) = v(i,j)
   end do
  end do
  !$omp end do
  !$omp end parallel

 !$omp parallel default(shared) private(i)
 !$omp do schedule(dynamic,64)
  do i=wall_num,nump-1
    if(((itypep(i) .ge. phase_min).and. &
    &   (itypep(i) .le. phase_max)).or. &
    &   (itypep(i) .eq. type(2)))then
      call cal_n_2                     &
      &    (i,1,ker_r,               &
      &     nn,idim,ifluid,neighbor, &
      &     grid_num,ww_num,         &
      &     ww_wall,wall_n,grid_p,   &
      &     an,neigh,x,itypep,       &
      &     wall_anc,dis,            &
      &     nn_wall,g_wall,wall_nor, &
      &     swi_hc,swi_solid)
     ! call cal_n                     &
     ! &    (i,1,ker_r,               &
     ! &     nn,idim,ifluid,neighbor, &
     ! &     grid_num,ww_num,         &
     ! &     ww_wall,wall_n,grid_p,   &
     ! &     an,neigh,x,itypep,       &
     ! &     wall_anc,dis,            &
     ! &     nn_wall,g_wall,wall_nor)
    end if
  end do
  !$omp end do
  !$omp end parallel

  ! calculation_external_force_term **********************************

  !$omp parallel default(shared) private(i)
  !$omp do schedule(dynamic,64)
  do i=wall_num,nump-1
    if(cal_sw(i,1) .eq. 1)then
!    if(itypep(i) .ge. type(3))then
    if(swi_solid(i) .eq. 0) then
       call external_force_term       &
&          (i,nn,idim,ifluid, &
&           accel,grav,nump,x,  &
&           num_surf,tau,ICMsurf, &
&           dis,den,phase, &
&           swi_hc)
    end if
    end if
  end do
  !$omp end do

  !$omp do schedule(dynamic,64)
  do i=0,nump-1
!    if(itypep(i) .ge. type(3))then
    if(swi_solid(i) .eq. 0) then
      call viscous_term                       &
      &    (i,nn,idim,ifluid,neighbor,wall_num, &
      &     nump,visc,itypep,type,cal_sw,     &
      &     neigh,x_temp,v_temp,nu,ker_r,rlambda,n0,    &
      &     grid_num,ww_num,grid,ww_wall,     &
      &     wall_n,grid_p,grid_dist,dis,      &
      &     near_wall,nn_wall,g_wall,wall_anc,&
      &     wall_nor,solid_sur)

!      call cal_dragforce2d( &
!      & i,nn,idim,dis,rho,x,v,&
!      & Flw(m)%x,Flw(m)%x,rhof,muf,Flw(m)%x,Flw(m)%x,drag,is,ie,js,je,ks,ke)
    end if

    if(itypep(i) .eq. 2 .and. swi_solid(i) .eq. 1) then
     v(i,0) = 0.0
     v(i,1) = 0.0
    end if

!  end do
!  !$omp end do
!
!  !$omp do schedule(dynamic,64)
!  do i=0,nump-1
    call prediction                         &
    &    (i,nn,idim,ifluid,neighbor,wall_num, &
    &     itypep,nump,type,v,x,dt,visc,     &
    &     accel,cal_sw,surf,time_sim)

  end do
  !$omp end do
  !$omp end parallel

  ! particle_search **************************************************

  call grid_division                  &
  &    (nn,nn_grid,idim,ifluid,nump,  &
  &     grid_num,grid_n,grid,         &
  &     g_neigh,x,itypep,type,grid_p)

  ! data copy for collision
  !$omp parallel default(shared) private(i)
  !$omp do schedule(dynamic,64)
  do j=0,idim-1
   do i=0,nump-1
    x_temp(i,j) = x(i,j)
    v_temp(i,j) = v(i,j)
   end do
  end do
  !$omp end do
  !$omp end parallel

  ! neighboring_particle_search **************************************

  !$omp parallel default(shared) private(i)
  !$omp do schedule(dynamic,64)
  do i=wall_num,nump-1
    call set_neigh                   &
    &    (i,                         &
    &     nn,idim,ifluid,neighbor,x, &
    &     type,neigh,itypep,         &
    &     ker_r,nn_grid,             &
    &     grid_num,grid_p,g_neigh)
  end do
  !$omp end do

  !$omp do schedule(dynamic,64)
  do i=0,nump-1
    ! delete *********************************************************
!    if(itypep(i) .eq. phase)then
    if(swi_solid(i) .eq. 0) then
      cal_sw(i,0) = 1
    else
      cal_sw(i,0) = 0
    end if
    ! delete *********************************************************

    if(itypep(i) .le. phase)then
      cal_sw(i,1) = 1
    else
      cal_sw(i,1) = 0
    end if
  end do
  !$omp end do

  ! particle_collision *********************************************
  !$omp do schedule(dynamic,64)
  do i=0,nump-1
    call collision                         &
    &    (i,nn,idim,ifluid,neighbor,nump,    &
    &     itypep,type,cal_sw,              &
    &     den,neigh,x,v,dist_min,coll_rat, &
    &     dt,phase,dis,wall_num,x_temp,v_temp)
  end do
  !$omp end do

  ! wall_collision ***************************************************
  !$omp do schedule(dynamic,64)
  do i=wall_num,nump-1
    call wall_collision                         &
    &    (i,idim,ifluid,nn,dis,dist_min,x,v,type, &
    &     nump,wall_num,wall_n,wall_anc,itypep, &
    &     nn_wall,g_wall,wall_nor,grid_p,grid_num)
  end do
  !$omp end do
  !$omp end parallel

  ! particle_search **************************************************
  call grid_division                  &
  &    (nn,nn_grid,idim,ifluid,nump,  &
  &     grid_num,grid_n,grid,         &
  &     g_neigh,x,itypep,type,grid_p)

  ! data copy for pressure
  !$omp parallel default(shared) private(i)
  !$omp do schedule(dynamic,64)
  do j=0,idim-1
   do i=0,nump-1
    x_temp(i,j) = x(i,j)
    v_temp(i,j) = v(i,j)
   end do
  end do
  !$omp end do
  !$omp end parallel

  ! neighboring_particle_search ***************************************

  !$omp parallel default(shared) private(i)
  !$omp do schedule(dynamic,64)
  do i=wall_num,nump-1
    call set_neigh                   &
    &    (i,                         &
    &     nn,idim,ifluid,neighbor,x, &
    &     type,neigh,itypep,         &
    &     ker_r,nn_grid,             &
    &     grid_num,grid_p,g_neigh)
  end do
  !$omp end do

  !$omp do schedule(dynamic,64)
  do i=wall_num,nump-1
    if(((itypep(i) .ge. phase_min).and. &
    &   (itypep(i) .le. phase_max)).or. &
    &   (itypep(i) .eq. type(2)))then
      call cal_n_div                       &
      &    (i,1,ker_r,cal_sw,              &
      &     nn,idim,ifluid,neighbor,phase, &
      &     an,neigh,x,itypep,type,ww_num, &
      &     grid_num,grid,ww_wall,wall_n,  &
      &     grid_p,grid_dist,wall_anc,dis, &
      &     near_wall,nn_wall,g_wall,wall_nor)
    end if

    ! rho, nu_calculation ******************************************
    call cal_rhonu                               &
    &    (i,nn,idim,ifluid,neighbor,phase,num_ref, &
    &     den,vis,rho,nu,nump,n0,an,wall_num)

    call cal_pressure                            &
    &    (i,nn,idim,ifluid,neighbor,nump,          &
    &     wall_num,cal_sw,itypep,type,rho,p,     &
    &     neigh,x,den,ker_r,n0,press,c_vel,phase,&
    &     grid_num,ww_num,grid,ww_wall,wall_n,   &
    &     grid_p,grid_dist,near_wall,wall_anc,   &
    &     dis,nn_wall,g_wall,wall_nor)
  end do
  !$omp end do

  !$omp do schedule(dynamic,64)
  do i=wall_num,nump-1
!    if(itypep(i) .eq. phase)then
    if(swi_solid(i) .eq. 0) then
      cal_sw(i,0:1) = 1
    else
      cal_sw(i,0:1) = 0
    end if

    if(itypep(i) .eq. type(3))then
      if(solid_sur(i) .eq. 0)then
        cal_sw(i,0) = 0
      end if
    end if

  ! calculation_pressure_term ************************************

    call pressure_term                            &
    &    (i,nn,idim,ifluid,neighbor,nump,           &
    &     wall_num,cal_sw,itypep,type,rho,p,      &
    &     neigh,x_temp,den,ker_r,n0,press,c_vel,phase, &
    &     grid_num,ww_num,grid,ww_wall,wall_n,    &
    &     grid_p,grid_dist,near_wall,wall_anc,    &
    &     dis,nn_wall,g_wall,wall_nor,swi_solid)

  ! correction_term **********************************************
    call correction                         &
    &    (i,nn,idim,ifluid,neighbor,wall_num, &
    &     nump,type,v,x,dt,cal_sw,press,itypep)
  end do
  !$omp end do
  !$omp end parallel

  ! particle_search **********************************************

  call grid_division                  &
  &    (nn,nn_grid,idim,ifluid,nump,  &
  &     grid_num,grid_n,grid,         &
  &     g_neigh,x,itypep,type,grid_p)

  !$omp parallel do default(shared) private(i) schedule(dynamic,64)
  do i=wall_num,nump-1
    if(cal_sw(i,1) .eq. 1)then
      call set_neigh                   &
      &    (i,                         &
      &     nn,idim,ifluid,neighbor,x, &
      &     type,neigh,itypep,         &
      &     ker_r,nn_grid,             &
      &     grid_num,grid_p,g_neigh)
    end if
  end do
  !$omp end parallel do

  if(swi_icing .eq.1)then
    do
      phase=type(4)
      changep = 0

      !$omp parallel do default(shared) private(i) schedule(dynamic,64)
      do i=wall_num,nump-1
        if(itypep(i) .eq. type(4))then
          icing_col(i) = 0

          ! 着氷衝突判定
!          call icing_judge_collision             &
!          &    (i,nn,idim,ifluid,neighbor,nump,    &
!          &     itypep,type,cal_sw,              &
!          &     den,neigh,x,v,dist_min,coll_rat, &
!          &     dt,phase,dis,wall_num,wall_n,wall_anc, &
!          &     grid_p,nn_wall,grid_num,g_wall,wall_nor, &
!          &     ice_dist_min,icing_col)

          ! 着氷密度判定
          if((icing_col(i) .eq. 1))then
            call cal_n_ice                      &
            &   (i,1,ker_r,ice_an,nump,        &
            &    nn,idim,ifluid,neighbor, &
            &    type,grid_num,ww_num,         &
            &    ww_wall,wall_n,grid_p,   &
            &    neigh,x,itypep,       &
            &    wall_anc,dis,            &
            &    nn_wall,g_wall,wall_nor)
          end if

          ! 着氷判定
          call icing                      &
          &    (i,nn,idim,ifluid,nump,x,v, &
          &     itypep,type,dis,ice_an,icing_col, &
          &     ice_dist_min,ice_an_min,changep)
        end if
      end do
      !$omp end parallel do

      if(changep .eq. 0)then
        exit
      end if
    end do

  else if(swi_icing .eq. 2)then

   !$omp parallel default(shared) private(i)
   !$omp do schedule(dynamic,64)
   do i=0,nump-1
    t_temp(i) =t(i)
    ganma_temp(i) =ganma(i)
    Qhc(i) = 0.0
   end do
  !$omp end do
  !$omp end parallel

      ! 温度場計算

    !$omp parallel default(shared) private(i)
    !$omp do schedule(dynamic,64)
    do i=wall_num,nump-1
      if(((itypep(i) .ge. phase_min).and. &
      &   (itypep(i) .le. phase_max)).or. &
      &   (itypep(i) .eq. type(2)))then
        call cal_n                     &
        &    (i,1,ker_r,               &
        &     nn,idim,ifluid,neighbor, &
        &     grid_num,ww_num,         &
        &     ww_wall,wall_n,grid_p,   &
        &     an,neigh,x,itypep,       &
        &     wall_anc,dis,            &
        &     nn_wall,g_wall,wall_nor)
      end if
    end do
    !$omp end do
    !$omp end parallel

  if(swi_heattransfer .eq. 1) then
    !$omp parallel do default(shared) private(i) schedule(dynamic,64)
    do i=wall_num,nump-1
!      if (an(i) .le. 10.02) then
      if(swi_hc(i) .eq. 1) then
        if(((itypep(i) .ge. type(2)).and.  &
        &  (itypep(i) .le. type(4)))) then
          call HeatTransfer(&
          & i,nn,idim,nump,num_surf,ICMsurf,&
          & x,t_temp,sur_temp,dis,hc,Qhc)
        end if
      end if
!      end if
    end do

    !$omp end parallel do
  end if

  !write(*,*) Qhc(10)

      !$omp parallel do default(shared) private(i) schedule(dynamic,64)
      !温度計算
      do i=0,nump-1
        if(((itypep(i) .ge. type(2)).and. &
        &  (itypep(i) .le. type(4)))) then
          if(flag_temp(i) .ne. 1) cycle
          call temperature                       &
          &    (i,nn,idim,ifluid,neighbor,wall_num, &
          &     nump,visc,itypep,type,cal_sw,     &
          &     neigh,x,v,nu,ker_r,rlambda,n0,    &
          &     grid_num,ww_num,grid,ww_wall,     &
          &     wall_n,grid_p,grid_dist,dis,      &
          &     near_wall,nn_wall,g_wall,wall_anc,&
          &     wall_nor,solid_sur,den,t,dt,      &
          &     t_wall,t_temp,ganma,ganma_temp,swi_koeki,Qhc, &
          &     swi_solid,flag_temp,wall,wall_dist)   !Qhc追加！！！
        end if
      end do
      !$omp end parallel do

  end if

  !$omp parallel default(shared) private(i)
  !$omp do schedule(dynamic,64)
  do i=wall_num,nump-1
!    if((itypep(i) .ge. phase_min).and. &
!    &  (itypep(i) .le. phase_max)) then
    if(swi_solid(i) .eq. 0) then
      call cal_n                     &
      &    (i,1,ker_r,               &
      &     nn,idim,ifluid,neighbor, &
      &     grid_num,ww_num,         &
      &     ww_wall,wall_n,grid_p,   &
      &     an,neigh,x,itypep,       &
      &     wall_anc,dis,            &
      &     nn_wall,g_wall,wall_nor)
    end if

    ! 低粒子数密度の粒子を除去
    if(an(i) .le. limit_low_an)then
      itypep(i) = type(0)
    end if
  end do
  !$omp end do

  !$omp do schedule(dynamic, 64)
  do i=0,nump-1
!    if((itypep(i) .ge. phase_min).and. &
!    &  (itypep(i) .le. phase_max))then
    if(swi_solid(i) .eq. 0) then

      cal_sw(i,0:1) = 1
    else
      cal_sw(i,0:1) = 0
    end if

    if(itypep(i) .eq. type(3))then
      if(solid_sur(i) .eq. 0)then
        cal_sw(i,0) = 0
      end if
    end if
  end do
 !$omp end do

! !$omp do schedule(dynamic,64)
! do i=0,nump-1
!   call set_rhonu                   &
!   &    (i,nn,idim,ifluid,neighbor,   &
!   &     itypep,den,vis,rho,nu,nump)
! end do
! !$omp end do
! !$omp end parallel

  !$omp do schedule(dynamic,64)
  do i=0,nump-1
    nu(i)  = vis(itypep(i))
    rho(i) = den(itypep(i))
  end do
  !$omp end do
  !$omp end parallel

  ! ******************************************************************
  ! *********************  MPS 計算パート終了  ***********************
  ! ******************************************************************
  end if 

  ! ******************************************************************
  ! *********************    格子法 to MPS     ***********************
  ! ******************************************************************

  call grid_division_type0            &
  &    (nn,nn_grid,idim,ifluid,nump,  &
  &     grid_num,grid_n,grid,         &
  &     dp_g_neigh,x,itypep,type,dp_grid_p,cal_type0)

  !$omp parallel default(shared) private(i)
  !$omp do schedule(dynamic,64)
  do i=0,nump-1
    call set_neigh_type0                   &
    &    (i,                         &
    &     nn,idim,ifluid,neighbor,x, &
    &     type,dp_neigh,itypep,         &
    &     ker_r,nn_grid,             &
    &     grid_num,dp_grid_p,dp_g_neigh,cal_type0)
  end do
  !$omp end do
  !$omp end parallel

  do i = 0,nump-1
    if((itypep(i) .eq. type(0)) .and. &
    &  (cal_type0(i) .eq. 2))then
      call switch_grid                       &
      &    (i,nn,idim,ifluid,neighbor,nump,  &
      &     itypep,type,cal_sw,              &
      &     den,dp_neigh,x,v,dist_min,coll_rat, &
      &     dt,phase,dis,wall_num,wall_n,wall_anc, &
      &     dp_grid_p,nn_wall,grid_num,g_wall,wall_nor, &
      &     ker_r,mps_swi,droplet_num,change_incount)
    end if

    if(mps_swi .eq. 1)then
      droplet_swi(change_incount) = 1
      write(*,*)'swi_mps'
      write(*,*)"time = ",time_sim
      swi_cal_mps = 1
      interval = 0.0

      !$omp parallel default(shared) private(j)
      !$omp do schedule(dynamic,64)
      do j = 0, nump-1
        if((itypep(j) .eq. type(0)) .and. &
        &  (cal_type0(j) .eq. 2)    .and. &
        &  (droplet_num(j) .eq. change_incount))then
          itypep(j) = type(4)
          swi_solid(j) = 0
          an(j) = 0.0
          p(j) = 0.0
          cal_type0(j) = 0
          droplet_num(j) = 0
          ganma(j) = ganma_ini

          if(idim .eq. 2)then
            v(j,0) = droplet_u(change_incount)
            v(j,1) = droplet_v(change_incount)
          else if(idim .eq. 3)then
            v(j,0) = droplet_u(change_incount)
            v(j,1) = droplet_v(change_incount)
            v(j,2) = droplet_w(change_incount)
          end if
        end if
      end do
      !$omp end do
      !$omp end parallel

      ! particle_search **************************************************
      call grid_division            &
      &    (nn,nn_grid,idim,ifluid,nump,  &
      &     grid_num,grid_n,grid,         &
      &     g_neigh,x,itypep,type,grid_p)

      !$omp parallel default(shared) private(j)
      !$omp do schedule(dynamic,64)
      do j=0,nump-1
        if((itypep(j) .ge. phase_min).and. &
        &  (itypep(j) .le. phase_max))then
          cal_sw(j,0:1) = 1
        else
          cal_sw(j,0:1) = 0
        end if

        if(itypep(j) .eq. type(3))then
          if(solid_sur(j) .eq. 0)then
            cal_sw(j,0) = 0
          end if
        end if
      end do
      !$omp end do
      !$omp end parallel

      !$omp parallel default(shared) private(j)
      !$omp do schedule(dynamic,64)
      do j=wall_num,nump-1
        call set_neigh                   &
        &    (j,                         &
        &     nn,idim,ifluid,neighbor,x, &
        &     type,neigh,itypep,         &
        &     ker_r,nn_grid,             &
        &     grid_num,grid_p,g_neigh)
      end do
      !$omp end do
      !$omp end parallel

!     !$omp parallel default(shared) private(j)
!     !$omp do schedule(dynamic,64)
!     do j=0,nump-1
!       call set_rhonu                   &
!       &    (j,nn,idim,ifluid,neighbor,   &
!       &     itypep,den,vis,rho,nu,nump)
!     end do
!     !$omp end do
!     !$omp end parallel

      do j=0,nump-1
        nu(j)  = vis(itypep(j))
        rho(j) = den(itypep(j))
      end do

      mps_swi = 0

    end if
  end do
  ! ******************************************************************
  ! *********************    格子法 to MPS     ***********************
  ! ******************************************************************

  ! ******************************************************************
  ! ************************  繰り返し計算  **************************
  ! ******************************************************************

  ! ******************************************************************
  ! ********************** データのアウトプット **********************
  ! ******************************************************************
  if(mod(ite,data_out) .eq. 0) then

    ! count : nump_out
    nump_out = 0
    do i=wall_num,nump-1
     if(itypep(i) .ge. phase_min)nump_out = nump_out+1
    end do

    ! output_calculation_state(file)
    write(900,'(a,i15  )') 'iteration       = ',ite
    write(900,'(a,e15.4)') 'dt              = ',dt
    write(900,'(a,f15.8)') 'time            = ',time_sim
    write(900,'(a,i15  )') 'num_out         = ',nump_out
    write(900,'(a,i15  )') 'wall_num        = ',wall_num
    write(900,'(a,i15  )') 'nump            = ',nump
    write(900,'(a,i15  )') 'neigh_max       = ',  &
    &                    maxval(neigh(0:nn-1,0))
    write(900,'(a,i15  )') 'g_neigh_max     = ',  &
    &                    maxval(g_neigh(1:grid_num,0))
    write(900,'(a,i15  )') 'fcount          = ',fcount-1
    write(900,'(a,e15.4)') 'Vmax            = ',vmax
    write(900,'(a,e15.5)') 'Interval        = ',interval
    write(900,'(a,e15.4)') 'Tmax            = ',  &
    &                    maxval(t(0:nump-1))
    write(900,'(a,e15.4)') 'Tmin            = ',  &
    &                    minval(t(0:nump-1))
    write(900,'(a      )') '---------------------------------'
  end if

  ! output_paraview_file *********************************************
  if(time_sim .ge. dtout*real(fcount))then
    if(fcount .le. output_limit)then
    ! output :vtk_file
!      call output_para_bin            &
!      &    (nn,idim,ifluid,fcount,t,  &
!      &     x,v,itypep,type,an,p,surf,cal_type0,swi_koeki,ganma)
      call output_para_bin            &
      &    (nn,idim,ifluid,fcount,t,  &
      &     x,v,itypep,type,an,p,surf,cal_type0,swi_koeki,ganma, &
      &     swi_solid,swi_hc)
      if (swi_date_and_time .eq.1)then
        write(*,*)'---------------------------------------------------------'
        write(*,*)'                 ParaViewData Output                     '
        write(*,*)'file number = ',fcount
        write(*,*)'sim time    = ',time_sim
        call date_and_time(date_time_cha(1),date_time_cha(2),date_time_cha(3), &
        &                  date_time)
        write(*,'(A,I4,A,I2,A,I2)')" DATE = ",date_time(1),"/",date_time(2),&
        &                          "/",date_time(3)
        write(*,'(A,I2,A,I2,A,I2)')" TIME = ",date_time(5),":",date_time(6),&
        &                          ":",date_time(7)
        write(*,*)'---------------------------------------------------------'

  ! output_initial_file **********************************************
!  call output_file                   &
!  &    (idim,nn,nump,ite,time_sim,   &
!  &     type,ifluid,itypep,          &
!  &     den,vis,grav,delta,sigma,    &
!  &     time_max,ite_max,            &
!  &     dt,dtout,cfl,                &
!  &     dis,ker_c,n0,                &
!  &     x,non_cal,num_ref,rlambda,   &
!  &     grid,grid_n,grid_num,        &
!  &     dis_rat,c_vel,               &
!  &     f_type,output_limit,         &
!  &     nn_grid,neighbor,t,wall_nor, &
!  &     dist_min,coll_rat,near_wall, &
!  &     wall_n,grid_dist,wall_anc)

      end if

      ! count_up
      fcount = fcount+1

      if ((swi_remesh .eq. 1) .or. (swi_remesh .eq. 2))then
        if (time_sim .ge. (real(count_remesh + 1) *time_remesh) &
          & .or. time_sim .ge. finish_time)then
          call write_restart_file(   &
          &    nn,idim,dis,mvd,nump,wall_num,droplet_nn,time_sim_dble, &
          &    time_sim,nextin,incount,fcount,count_remesh,           &
          &    x,v,t,itypep,ganma,swi_koeki,cal_type0,droplet_num,    &
          &    droplet_swi,droplet_x,droplet_y,droplet_z,             &
          &    droplet_u,droplet_v,droplet_w)
  
          write(*,*)'**********************************************'
          write(*,*)'********** Time reach to remesh **************'
          write(*,*)'**********************************************'
          write(*,*)'Simulation time  =',time_sim
          write(*,*)'Remesh count     =',count_remesh
          write(*,*)'Please recalculate flow field'
          call date_and_time(date_time_cha(1),date_time_cha(2),date_time_cha(3), &
          &                  date_time)
          write(*,'(A,I4,A,I2,A,I2)')" DATE = ",date_time(1),"/",date_time(2),&
          &                          "/",date_time(3)
          write(*,'(A,I2,A,I2,A,I2)')" TIME = ",date_time(5),":",date_time(6),&
          &                          ":",date_time(7)
          write(*,*)'Calculation stop'
          stop
        end if
      end if

      if(fcount .gt. output_limit) then
        write(*,*)'**********************************************'
        write(*,*)'************ data output overflow ************'
        write(*,*)'**********************************************'
        call date_and_time(date_time_cha(1),date_time_cha(2),date_time_cha(3), &
        &                  date_time)

        if(swi_parameter.eq.2)then
          ICE_NUM = 0
          DO I=WALL_NUM,NUMP-1
            IF(ITYPEP(I) .EQ. TYPE(2))THEN
              ICE_NUM = ICE_NUM + 1
            END IF
          END DO

          WATER_NUM = INITIAL_NUM-ICE_NUM
          INITIAL_VOL = (DIS**3)/6*REAL(INITIAL_NUM)
          ICE_VOL   = (DIS**3)/6*REAL(ICE_NUM)
          WATER_VOL = (DIS**3)/6*REAL(WATER_NUM)

          VOL_RATE  = WATER_VOL/INITIAL_VOL

          WRITE(*,*)
          WRITE(*,*)"******** CALCULATION RESULT ***********"
          WRITE(*,*)"ICE_DIST_MIN            = ",ICE_DIST_MIN
          WRITE(*,*)"ICE_AN_MIN              = ",ICE_AN_MIN
          WRITE(*,*)"ANGLE OF ATTACK         = ",AOAdeg
          WRITE(*,*)"CONTACT ANGLE           = ",CA
          WRITE(*,*)"INITIAL PARTICLE NUMBER = ",INITIAL_NUM
          WRITE(*,*)"ICING PARTICLE NUMBER   = ",ICE_NUM
          WRITE(*,*)"WATER PARTICLE NUMBER   = ",WATER_NUM
          WRITE(*,*)"VOL_RATE                = ",VOL_RATE
          WRITE(*,*)"***************************************"

          open(90,file='./VOL_RATE.txt',position='append')
            write(90,*)ICE_AN_MIN,CA,AOAdeg,VOL_RATE
          close(90)

          AOAdeg = AOAdeg+5.0

          if(AOAdeg .eq. 95.0)then
          !  CA = CA + 15.0
            ice_an_min = ice_an_min + 0.5
            AOAdeg = 5.0
          end if

          open(90,file='./AOA.txt')
            write(90,*)ice_an_min
            write(90,*)AOAdeg
            write(90,*)CA
          close(90)
        end if


        ! 計算終了後形状の出力
!        call write_circle_position(   &
!        &    nn,idim,dis,mvd,x,nump,temp_nump, &
!        &    randomx,randomy,randomz           &
!        &    )

!        call write_sphere_position(   &
!        &    nn,idim,dis,mvd,x,nump,temp_nump, &
!        &    randomx,randomy,randomz           &
!        &    )


        write(*,*)
        write(*,*)"**************  CALUCATION END  *********************"
        write(*,'(A,I4,A,I2,A,I2)')" DATE = ",date_time(1),"/",date_time(2),&
        &                          "/",date_time(3)
        write(*,'(A,I2,A,I2,A,I2)')" TIME = ",date_time(5),":",date_time(6),&
        &                          ":",date_time(7)
        write(*,*)
       stop
      end if
    end if
  end if

  if((time_sim .le. time_max) .and.   &
  &  (ite      .le. ite_max)) goto 1000
  close(80)

  write(*,'(a)')'**********************************************'
  write(*,'(a)')'**************** main loop end ***************'
  write(*,'(a)')'**********************************************'
  call date_and_time(date_time_cha(1),date_time_cha(2),date_time_cha(3), &
  &                  date_time)
  write(*,*)
  write(*,*)"**************  CALUCATION END  *********************"
  write(*,'(A,I4,A,I2,A,I2)')"DATE = ",date_time(1),"/",date_time(2),&
  &                          "/",date_time(3)
  write(*,'(A,I2,A,I2,A,I2)')"TIME = ",date_time(5),":",date_time(6),&
  &                          ":",date_time(7)
  write(*,*)
  close(800)
  close(900)

  contains

  ! ******************************************************************
  ! ** モジュールに持って行きたいけど変数ぐちゃぐちゃ 超めんどい *****
  ! ******************************************************************

  !*******************************************************************
  !********               検証実験ケース選択                  ********
  !*******************************************************************
  subroutine SelectExpCase
   ! 変数宣言 ********************************************************
   implicit none
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   character :: fname * 20
   ! 処理開始 ********************************************************
   ! 計算条件ファイル入力
   call Input_CalSetting( './data/GridFlow/' // trim(ND_CalSetFile) // strtxt )
   ! ディレクトリ設定
   if( IceStep == 0 ) then
     GrdInDir    = './result/GridFlow/grid/clean/'
     OSGDir      = './result/GridFlow/overset/clean/'
     FlwIniDir   = './result/GridFlow/flow/initial/clean/'
     FlwCalInDir = './result/GridFlow/flow/cal/clean/'
     DrpImpDir   = './result/GridFlow/droplet/impingement/clean/'
     DrpTraDir   = './result/GridFlow/droplet/trajectory/clean/'
    else
     GrdInDir    = './result/GridFlow/grid/icing/'
     OSGDir      = './result/GridFlow/overset/icing/'
     FlwIniDir   = './result/GridFlow/flow/initial//icing/'
     FlwCalInDir = './result/GridFlow/flow/cal/clean/'
     DrpImpDir   = './result/GridFlow/droplet/impingement/icing/'
     DrpTraDir   = './result/GridFlow/droplet/trajectory/icing/'
   endif
!   if( IceStep == 0 ) then
!     GrdInDir    = './data/grid-based/grid/clean/'
!     OSGDir      = './data/grid-based/overset/clean/'
!     FlwIniDir   = './data/grid-based/flow/initial/clean/'
!     FlwCalInDir = './data/grid-based/flow/cal/clean/'
!     DrpImpDir   = './data/grid-based/droplet/impingement/clean/'
!     DrpTraDir   = './data/grid-based/droplet/trajectory/clean/'
!    else
!     GrdInDir    = './data/grid-based/grid/icing/'
!     OSGDir      = './data/grid-based/overset/icing/'
!     FlwIniDir   = './data/grid-based/flow/initial//icing/'
!     FlwCalInDir = './data/grid-based/flow/cal/clean/'
!     DrpImpDir   = './data/grid-based/droplet/impingement/icing/'
!     DrpTraDir   = './data/grid-based/droplet/trajectory/icing/'
!   endif
   write(*, '(a,i2)') '* Ice step      = ', IceStep
   write(*, '(a,i2)') '* Ice step max. = ', IceStepMax
   write(*, '(a,e16.8e3)') '* Ts    = ', TsExp * aRef**2
   write(*, '(a,e16.8e3)') '* Ps    = ', PsExp * (rhoRef * aRef**2)
   write(*, '(a,e16.8e3)') '* V     = ', VelExp * aRef
   write(*, '(a,e16.8e3)') '* LWC   = ', LWC * RhoRef
   write(*, '(a,e16.8e3)') '* MVD   = ', MVD * LRef
   write(*, '(a,e16.8e3)') '* Rho   = ', Rhod * RhoRef
   write(*, '(a,e16.8e3)') '* Chord = ', Chord * LRef
   write(*, '(a,e16.8e3)') '* AOA   = ', AOA * 180.0 / pi
   ! 処理終了 ********************************************************
   return
  end subroutine SelectExpCase

  !***********************************************
  !********      格子パラメータ変更       ********
  !***********************************************
  subroutine changeAOAandPOS(xpos,ypos,dis)
   ! 変数宣言 ********************************************************************************************
   implicit none

   real   , intent(in) :: xpos
   real   , intent(in) :: ypos
   real   , intent(in) :: dis
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer   :: i, j, k, m
   real      :: theta

   type(Type_Flow), pointer :: Flw2(:)
   ! 処理開始 ********************************************************************************************
   ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   theta = - AOA
   allocate( Flw2(ms:me))

   ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! 一時変数の作成
   do m = ms, me
    allocate( Flw2(m)%x   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw2(m)%y   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw2(m)%z   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw2(m)%u   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw2(m)%v   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw2(m)%w   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke))
   enddo

   do m = ms,me
   do k = Flw(m)%ks,Flw(m)%ke
   do j = Flw(m)%js,Flw(m)%je
   do i = Flw(m)%is,Flw(m)%ie
    Flw2(m)%x(i,j,k) = (Flw(m)%x(i,j,k)*lRef)*cos(theta) - (Flw(m)%y(i,j,k)*lRef) * sin(theta)
    Flw2(m)%y(i,j,k) = (Flw(m)%x(i,j,k)*lRef)*sin(theta) + (Flw(m)%y(i,j,k)*lRef) * cos(theta)
    Flw2(m)%u(i,j,k) = (Flw(m)%u(i,j,k)*aRef)*cos(theta) - (Flw(m)%v(i,j,k)*aRef) * sin(theta)
    Flw2(m)%v(i,j,k) = (Flw(m)%u(i,j,k)*aRef)*sin(theta) + (Flw(m)%v(i,j,k)*aRef) * cos(theta)
   end do
   end do
   end do
   end do

   do m = ms,me
   do k = Flw(m)%ks,Flw(m)%ke
   do j = Flw(m)%js,Flw(m)%je
   do i = Flw(m)%is,Flw(m)%ie
    Flw(m)%x(i,j,k) = (Flw2(m)%x(i,j,k) + xpos)/lRef
    Flw(m)%y(i,j,k) = (Flw2(m)%y(i,j,k) + ypos)/lRef
    Flw(m)%u(i,j,k) = (Flw2(m)%u(i,j,k))/aRef
    Flw(m)%v(i,j,k) = (Flw2(m)%v(i,j,k))/aRef
   end do
   end do
   end do
   end do

   ! 応急処置
   ! 翼格子が薄いので伸ばす
   ! 格子法が2次元計算の場合のみ可能

   do m = ms,me
   do k = Flw(m)%ks,Flw(m)%ke
   do j = Flw(m)%js,Flw(m)%je
   do i = Flw(m)%is,Flw(m)%ie
     if(idim .eq. 2)then
       Flw(m)%z(i,j,k) = (-0.5*dis + (real(k)/(Flw(m)%ke-Flw(m)%ks))*dis) /lRef
     !  可視化用
     !  Flw(m)%z(i,j,k) = (-5.0*dis + (real(k)/(Flw(m)%ke-Flw(m)%ks))*dis) /lRef
     else if(idim .eq. 3)then
       Flw(m)%z(i,j,k) = (0.00 + (real(k)/(Flw(m)%ke-Flw(m)%ks))*30.0*dis) /lRef
     end if

   end do
   end do
   end do
   end do

   m = 2
   do k = Flw(m)%ks,Flw(m)%ke
   do j = Flw(m)%js,Flw(m)%je
   do i = Flw(m)%is,Flw(m)%ie
    Flw(m)%x(i,j,k) = (Flw2(m)%x(i,j,k) + xpos)/lRef
    Flw(m)%y(i,j,k) = (Flw2(m)%y(i,j,k) + ypos)/lRef
    Flw(m)%u(i,j,k) = (Flw2(m)%u(i,j,k))/aRef
    Flw(m)%v(i,j,k) = (Flw2(m)%v(i,j,k))/aRef
   end do
   end do
   end do
   deallocate( Flw2)

   ! 処理終了 ********************************************************************************************
   return
  end subroutine changeAOAandPOS

  !***********************************************
  !********        Paraview出力           ********
  !***********************************************
  ! 値が無次元化されているので、有次元に戻してあげる
  ! aRef   : 音速                      速度、温度などの無次元化に用いる
  ! rhoRef : 密度の最大値              密度などの無次元化に用いる
  ! lRef   : 代表長さ(翼ではコード長)  長さなどの無次元化に用いる

  subroutine Output_Grid3D_para( &
  &            strname, is, ie, js, je, ks, ke, x, y, z )
   ! 変数宣言 ********************************************************************************************
   implicit none
   ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   character, intent(in) :: strname*(*)
   integer  , intent(in) :: is, ie, js, je, ks, ke
   real     , intent(in) :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke), z(is:ie, js:je, ks:ke)
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer :: i, j, k
   character :: strg * 4 = '.g'
   character :: strq * 4 = '.q'
   character :: strfun * 4 = '.fun'
   ! 処理開始 ********************************************************************************************
   open(1, file = trim(strname) // strg , form = 'unformatted', status = 'replace')

    write(1)ie-is+1, je-js+1, ke-ks+1

    write(1) (((x(i,j,k)*lRef, i = is, ie), j = js, je), k = ks, ke), &
    &        (((y(i,j,k)*lRef, i = is, ie), j = js, je), k = ks, ke), &
    &        (((z(i,j,k)*lRef, i = is, ie), j = js, je), k = ks, ke)


   close(1)
   ! 処理終了 ********************************************************************************************
   return
  end subroutine Output_Grid3D_para
  subroutine Output_Q3D_para( &
  &          strname, is, ie, js, je, ks, ke, rho, u, v, w, &
  &          p, t, mu, kin, eps )
   ! 変数宣言 ********************************************************************************************
   implicit none
   ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   character, intent(in) :: strname*(*)
   integer  , intent(in) :: is, ie, js, je, ks, ke
   real     , intent(in) :: rho(is:ie, js:je, ks:ke), u(is:ie, js:je, ks:ke), v(is:ie, js:je, ks:ke)
   real     , intent(in) :: w(is:ie, js:je, ks:ke), p(is:ie, js:je, ks:ke), t(is:ie, js:je, ks:ke)
   real     , intent(in) :: mu(is:ie, js:je, ks:ke), kin(is:ie, js:je, ks:ke), eps(is:ie, js:je, ks:ke)
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer :: i, j, k
   character :: strg * 4 = '.g'
   character :: strq * 4 = '.q'
   character :: strfun * 4 = '.fun'
   integer :: fun_num = 9        !変数の個数

   integer   :: count
   real      :: mu_ave
   real      :: out_rho (is:ie, js:je, ks:ke)
   real      :: out_urho(is:ie, js:je, ks:ke)
   real      :: out_vrho(is:ie, js:je, ks:ke)
   real      :: out_wrho(is:ie, js:je, ks:ke)
   real      :: out_e   (is:ie, js:je, ks:ke)

   ! 処理開始 ********************************************************************************************

   count = 0
   mu_ave = 0.0

   do k = ks,ke
   do j = js,je
   do i = is,ie
    count = count + 1
    mu_ave = mu_ave + mu(i,j,k)
    out_rho(i,j,k)  = rho(i,j,k)*rhoRef
    out_urho(i,j,k) = u(i,j,k)*aRef*rho(i,j,k)*rhoRef
    out_vrho(i,j,k) = v(i,j,k)*aRef*rho(i,j,k)*rhoRef
    out_wrho(i,j,k) = w(i,j,k)*aRef*rho(i,j,k)*rhoRef
    out_e(i,j,k)    = ((t(i,j,k)*aref)**2)*Rg / (gamma-1)

   end do
   end do
   end do

   mu_ave = mu_ave / count

   open(1, file = trim(strname) // strq , form = 'unformatted', status = 'replace')

    ! 格子数
    write(1)ie-is+1, je-js+1, ke-ks+1
    ! マッハ数、AOA、レイノルズ数、時間
    write(1)VelExp * aRef / aRef, AOA * 180.0 / pi, &
    &       VelExp*Rhod*Chord / mu_ave, 0.0

    write(1) (((out_rho(i,j,k),   i = is, ie), j = js, je), k = ks, ke), &
    &        (((out_urho(i,j,k),  i = is, ie), j = js, je), k = ks, ke), &
    &        (((out_vrho(i,j,k),  i = is, ie), j = js, je), k = ks, ke), &
    &        (((out_wrho(i,j,k),  i = is, ie), j = js, je), k = ks, ke), &
    &        (((out_e(i,j,k),     i = is, ie), j = js, je), k = ks, ke)

   close(1)
   ! 処理終了 ********************************************************************************************
   return
  end subroutine Output_Q3D_para
  subroutine Output_Function3D_para( &
  &          strname, is, ie, js, je, ks, ke, rho, u, v, w, &
  &          p, t, mu, kin, eps )
   ! 変数宣言 ********************************************************************************************
   implicit none
   ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   character, intent(in) :: strname*(*)
   integer  , intent(in) :: is, ie, js, je, ks, ke
   real     , intent(in) :: rho(is:ie, js:je, ks:ke), u(is:ie, js:je, ks:ke), v(is:ie, js:je, ks:ke)
   real     , intent(in) :: w(is:ie, js:je, ks:ke), p(is:ie, js:je, ks:ke), t(is:ie, js:je, ks:ke)
   real     , intent(in) :: mu(is:ie, js:je, ks:ke), kin(is:ie, js:je, ks:ke), eps(is:ie, js:je, ks:ke)
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer :: i, j, k
   character :: strg * 4 = '.g'
   character :: strq * 4 = '.q'
   character :: strfun * 4 = '.fun'
   integer :: fun_num = 9        !変数の個数

   ! 処理開始 ********************************************************************************************
   open(1, file = trim(strname) // strfun , form = 'unformatted', status = 'replace')

    write(1)ie-is+1, je-js+1, ke-ks+1, fun_num

    write(1) (((rho(i,j,k)*rhoRef,          i = is, ie), j = js, je), k = ks, ke), &
    &        (((u(i,j,k)*aRef,              i = is, ie), j = js, je), k = ks, ke), &
    &        (((v(i,j,k)*aRef,              i = is, ie), j = js, je), k = ks, ke), &
    &        (((w(i,j,k)*aRef,              i = is, ie), j = js, je), k = ks, ke), &
    &        (((p(i,j,k)*rhoRef*aRef**2,    i = is, ie), j = js, je), k = ks, ke), &
    &        (((t(i,j,k)*aRef**2,           i = is, ie), j = js, je), k = ks, ke), &
    &        (((mu(i,j,k)*rhoRef*aRef*lRef, i = is, ie), j = js, je), k = ks, ke), &
    &        (((kin(i,j,k)*aref**2,         i = is, ie), j = js, je), k = ks, ke), &
    &        (((eps(i,j,k)*aref**3/lRef,    i = is, ie), j = js, je), k = ks, ke)


   close(1)
   ! 処理終了 ********************************************************************************************
   return
  end subroutine Output_Function3D_para


  subroutine DropletTimeStep(ip,jp,kp,mp)
  ! 変数宣言 ********************************************************************************************
  implicit none
  ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  integer, intent(inout) :: ip, jp, kp,mp
  integer, parameter :: nRetry = 100         ! 液滴探索最大再試行回数
  integer, parameter :: mRetry = 10          ! 液滴補間最大再試行回数
  integer, parameter :: SR     = 2           ! 液滴探索範囲
  real   , parameter :: MGN    = 1.0e-3          ! 補間経緯数許容誤差
  ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  integer :: mp
  real    :: xp0, yp0, zp0, up0, vp0, wp0
  integer :: n, m
  integer :: ip0, jp0, kp0
  real(8) :: rr
  real    :: mr
  logical :: fMove, fImpi, fExit
  ! 処理開始 ********************************************************************************************
  ! 液滴初期位置及び速度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  mp   = mini
  mr   = 1.0
  fSplash = .false.
  fBounce = .false.

  ! 液滴位置探索 ---------------------------------------------------------------------------------------
  ip0 = ip; jp0 = jp; kp0 = kp
  do n = 1, nRetry
   xp0 = xp; yp0 = yp; zp0 = zp
   up0 = up; vp0 = vp; wp0 = wp
   ! 液滴探索
   do m = 1, mRetry
    call SearchDroplet( &
    &      mp, &
    &      max(Flw(mp)%is, ip0 - SR), min(Flw(mp)%ie - 1, ip0 + SR), &
    &      max(Flw(mp)%js, jp0 - SR), min(Flw(mp)%je - 1, jp0 + SR), &
    &      max(Flw(mp)%ks, kp0 - SR), min(Flw(mp)%ke - 1, kp0 + SR), &
    &      MGN * real(m - 1), ip, jp, kp )
    if(fSearch) exit
   enddo
   xp0 = xp; yp0 = yp; zp0 = zp
   up0 = up; vp0 = vp; wp0 = wp
   ! 時間進行
   call TimeEuler3D( &
   &      dp_dt, fx, fy, fz, xp0, yp0, zp0, up0, vp0, wp0, &
   &      xp, yp, zp, up, vp, wp )

   if(fSearch) then
     exit
    else
     call BoundaryBladeSurface( ip0, jp0, kp0, Flw(mp)%js, mp, fImpi, up, vp, wp, mr )
     if(fImpi) return
   endif
  enddo
  if(.not. fSearch) then
    write(*, '(a)') '!!!!! Error : Droplet search !!!!!'
    write(*, '(a, i6)')  '* Iteration number = ', nCount
    write(*, '(a, i11)') '* Droplet number   = ', incount
    write(*, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip0, ',', jp0, ',', kp0, ')'
    Nlost = Nlost + 1
    swi_lost = 1
    return                ! ロスト．次の液滴へ
  endif
  ! 境界条件 -------------------------------------------------------------------------------------------
  fMove = .false.; fImpi = .false.; fExit = .false.
  select case(mp)
   ! Main Grid
   case(1)
    ! 流出境界
    if(.not. fMove) call BoundaryOutlet( ip, jp, mp, fExit )
    if(fExit) then
    end if
    ! Sub Grid 境界
    if(.not. fMove) call BoundaryMGtoSG( ip, jp, kp, mp, fMove )
    ! C 型格子ブランチ・カット
    if(.not. fMove) call BoundaryCtypeBranchCut( ip, jp, mp, fMove )
    ! 翼表面
    if(.not. fMove) call BoundaryBladeSurface( ip, jp, kp, Flw(mp)%js, mp, fImpi, up, vp, wp, mr )
    if(fImpi) return
    ! 周期境界
    if(.not. fMove) call BoundaryPeriodic( ip, jp, kp, mp, fMove )
   ! Sub Grid
   case(2)
    ! Main Grid 境界
    if(.not. fMove) call BoundarySGtoMG( ip, jp, kp, mp, fMove )
    ! 翼表面
    if(.not. fMove) call BoundaryBladeSurface( ip, jp, kp, Flw(mp)%js, mp, fImpi, up, vp, wp, mr )
    if(fImpi) return
    ! 周期境界
    if(.not. fMove) call BoundaryPeriodic( ip, jp, kp, mp, fMove )
   case default
    write(*, '(a)') '!!!!! Error : Block index !!!!!'
    stop
  end select

  ! 格子移動後の液滴位置探索 ---------------------------------------------------------------------------
  if(fMove) then
    do n = 1, nRetry
     ! 液滴位置探索
     do m = 1, mRetry
      call SearchDroplet( &
      &      mp, &
      &      Flw(mp)%is, Flw(mp)%ie - 1, &
      &      Flw(mp)%js, Flw(mp)%je - 1, &
      &      Flw(mp)%ks, Flw(mp)%ke - 1, &
      &      MGN * real(m - 1), ip, jp, kp )
      if(fSearch) exit
     enddo
     if(fSearch) exit
     ! 時間進行 (探索失敗時)
     xp0 = xp; yp0 = yp; zp0 = zp
     up0 = up; vp0 = vp; wp0 = wp
     call TimeEuler3D( &
     &      0.1 * dp_dt, fx, fy, fz, xp0, yp0, zp0, up0, vp0, wp0, &
     &      xp, yp, zp, up, vp, wp )
     ! 格子を突き抜けた場合の処理
     call BoundaryBladeSurface( ip, jp, kp, Flw(mp)%js, mp, fImpi, up, vp, wp, mr )
     if(fImpi) return
    enddo
    if(.not. fSearch) then
      write(*, '(a)') '!!!!! Error : Droplet search !!!!!'
      write(*, '(a, i6)')  '* Iteration number = ', nCount
     ! write(*, '(a, i11)') '* Droplet number   = ', nDrop
      write(*, '(a, i11)') '* Droplet number   = ', incount
      write(*, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip0, ',', jp0, ',', kp0, ')'
      Nlost = Nlost + 1
      swi_lost = 1
      return                ! ロスト．次の液滴へ
    endif
  endif
! ! 液滴軌道データログ ---------------------------------------------------------------------------------
! Drp(DrpNum)%step       = nCount
! Drp(DrpNum)%m (nCount) = mp
! Drp(DrpNum)%x (nCount) = xp
! Drp(DrpNum)%y (nCount) = yp
! Drp(DrpNum)%z (nCount) = zp
! Drp(DrpNum)%u (nCount) = up
! Drp(DrpNum)%v (nCount) = vp
! Drp(DrpNum)%w (nCount) = wp
! Drp(DrpNum)%mr(nCount) = mr


  ! 処理終了 ********************************************************************************************
  return
  end subroutine DropletTimeStep
  !*******************************************************************************************************
  !******** 液滴探索                    ********
  !*******************************************************************************************************
  subroutine SearchDroplet( &
  &            mp, isp, iep, jsp, jep, ksp, kep, MGN, ip, jp, kp )
   ! 変数宣言 ********************************************************************************************
   implicit none
   ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer, intent(in)  :: mp           ! 補間点探索ブロック
   integer, intent(in)  :: isp, iep, jsp, jep, ksp, kep     ! 補間点探索範囲
   real   , intent(in)  :: MGN            ! 補間係数許容誤差
   integer, intent(out) :: ip, jp, kp         ! 補間点格子番号
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer :: i1, i2, i3, i4, i5, i6, i7, i8, &
   &          j1, j2, j3, j4, j5, j6, j7, j8, &
   &          k1, k2, k3, k4, k5, k6, k7, k8        ! 補間点
   real    :: alp, bet, gam           ! 補間係数
   integer :: fInter              ! 補間パターン
   real    :: term1, term2, term3, term4, term5, term6, term7, term8  ! 補間係数
   real    :: rhof, uf, vf, wf, muf                                       ! 周囲流体の物理量
   real    :: fdx, fdy, fdz           ! 液滴に働く抗力
   real    :: fgx, fgy, fgz           ! 液滴に働く重力
   real    :: gx, gy, gz              ! 重力係数
   real    :: xixp, xiyp, xizp, etxp, etyp, etzp, zexp, zeyp, zezp
   ! 処理開始 ********************************************************************************************
   ! 探索開始 (物理空間) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   fInter = 0
   do kp = ksp, kep
   do jp = jsp, jep
   do ip = isp, iep
    ! ワイド・サーチ -------------------------------------------------------------------------------------
    i1 = ip    ; j1 = jp    ; k1 = kp    ; i2 = ip + 1; j2 = jp    ; k2 = kp
    i3 = ip    ; j3 = jp + 1; k3 = kp    ; i4 = ip + 1; j4 = jp + 1; k4 = kp
    i5 = ip    ; j5 = jp    ; k5 = kp + 1; i6 = ip + 1; j6 = jp    ; k6 = kp + 1
    i7 = ip    ; j7 = jp + 1; k7 = kp + 1; i8 = ip + 1; j8 = jp + 1; k8 = kp + 1
    call WideSearch3D8Point( &
    &      xp, yp, zp, &
    &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
    &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
    &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
    &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
    &      Flw(mp)%x(i5,j5,k5), Flw(mp)%y(i5,j5,k5), Flw(mp)%z(i5,j5,k5), &
    &      Flw(mp)%x(i6,j6,k6), Flw(mp)%y(i6,j6,k6), Flw(mp)%z(i6,j6,k6), &
    &      Flw(mp)%x(i7,j7,k7), Flw(mp)%y(i7,j7,k7), Flw(mp)%z(i7,j7,k7), &
    &      Flw(mp)%x(i8,j8,k8), Flw(mp)%y(i8,j8,k8), Flw(mp)%z(i8,j8,k8), &
    &      MGN, fSearch )
    if(.not. fSearch) then
     cycle
    endif
    ! 三重線形補間 ---------------------------------------------------------------------------------------
    alp = 0.5; bet = 0.5; gam = 0.5
    call Interpolation3D8Point( &
    &      xp, yp, zp, &
    &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
    &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
    &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
    &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
    &      Flw(mp)%x(i5,j5,k5), Flw(mp)%y(i5,j5,k5), Flw(mp)%z(i5,j5,k5), &
    &      Flw(mp)%x(i6,j6,k6), Flw(mp)%y(i6,j6,k6), Flw(mp)%z(i6,j6,k6), &
    &      Flw(mp)%x(i7,j7,k7), Flw(mp)%y(i7,j7,k7), Flw(mp)%z(i7,j7,k7), &
    &      Flw(mp)%x(i8,j8,k8), Flw(mp)%y(i8,j8,k8), Flw(mp)%z(i8,j8,k8), &
    &    MGN, alp, bet, gam, fSearch )
    if(fSearch) then
      fInter = 1
      exit
    endif
    ! 三次元線形補間 -------------------------------------------------------------------------------------
    ! Pattern-1
    i1 = ip    ; j1 = jp    ; k1 = kp    ; i2 = ip    ; j2 = jp    ; k2 = kp + 1
    i3 = ip + 1; j3 = jp    ; k3 = kp + 1; i4 = ip    ; j4 = jp + 1; k4 = kp + 1
    alp = 0.5; bet = 0.5; gam = 0.5
    call Interpolation3D4Point( &
    &      xp, yp, zp, &
    &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
    &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
    &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
    &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
    &  MGN, alp, bet, gam, fSearch )
    if(fSearch) then
      fInter = 2
      exit
    endif
    ! Pattern-2
    i1 = ip    ; j1 = jp    ; k1 = kp    ; i2 = ip    ; j2 = jp + 1; k2 = kp
    i3 = ip + 1; j3 = jp + 1; k3 = kp    ; i4 = ip    ; j4 = jp + 1; k4 = kp + 1
    alp = 0.5; bet = 0.5; gam = 0.5
    call Interpolation3D4Point( &
    &      xp, yp, zp, &
    &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
    &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
    &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
    &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
    &  MGN, alp, bet, gam, fSearch )
    if(fSearch) then
      fInter = 3
      exit
    endif
    ! Pattern-3
    i1 = ip    ; j1 = jp    ; k1 = kp    ; i2 = ip + 1; j2 = jp + 1; k2 = kp
    i3 = ip + 1; j3 = jp    ; k3 = kp + 1; i4 = ip    ; j4 = jp + 1; k4 = kp + 1
    alp = 0.5; bet = 0.5; gam = 0.5
    call Interpolation3D4Point( &
    &      xp, yp, zp, &
    &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
    &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
    &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
    &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
    &  MGN, alp, bet, gam, fSearch )
    if(fSearch) then
      fInter = 4
      exit
    endif
    ! Pattern-4
    i1 = ip    ; j1 = jp    ; k1 = kp    ; i2 = ip + 1; j2 = jp    ; k2 = kp
    i3 = ip + 1; j3 = jp + 1; k3 = kp    ; i4 = ip + 1; j4 = jp    ; k4 = kp + 1
    alp = 0.5; bet = 0.5; gam = 0.5
    call Interpolation3D4Point( &
    &      xp, yp, zp, &
    &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
    &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
    &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
    &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
    &  MGN, alp, bet, gam, fSearch )
    if(fSearch) then
      fInter = 5
      exit
    endif
    ! Pattern-5
    i1 = ip + 1; j1 = jp + 1; k1 = kp    ; i2 = ip + 1; j2 = jp    ; k2 = kp + 1
    i3 = ip    ; j3 = jp + 1; k3 = kp + 1; i4 = ip + 1; j4 = jp + 1; k4 = kp + 1
    alp = 0.5; bet = 0.5; gam = 0.5
    call Interpolation3D4Point( &
    &      xp, yp, zp, &
    &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
    &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
    &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
    &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
    &  MGN, alp, bet, gam, fSearch )
    if(fSearch) then
      fInter = 6
      exit
    endif
   enddo
   if(fInter /= 0) exit
   enddo
   if(fInter /= 0) exit
   enddo
   ! 補間係数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   select case(fInter)
    case(0)
     return
    case(1)
     term1 = (1.0 - alp) * (1.0 - bet) * (1.0 - gam)
     term2 =        alp  * (1.0 - bet) * (1.0 - gam)
     term3 = (1.0 - alp) *        bet  * (1.0 - gam)
     term4 =        alp  *        bet  * (1.0 - gam)
     term5 = (1.0 - alp) * (1.0 - bet) *        gam
     term6 =        alp  * (1.0 - bet) *        gam
     term7 = (1.0 - alp) *        bet  *        gam
     term8 =        alp  *        bet  *        gam
    case(2:6)
     term1 = 1.0
     term2 = alp
     term3 = bet
     term4 = gam
    case default
     write(*, '(a)') '!!!!! Error : fInter number !!!!!'
     stop
   end select
   ! 周囲流体の情報を補間 (計算空間) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   select case(fInter)
    case(1)
     xixp = term1 * Flw(mp)%xix(i1,j1,k1) + term2 * Flw(mp)%xix(i2,j2,k2) &
     &    + term3 * Flw(mp)%xix(i3,j3,k3) + term4 * Flw(mp)%xix(i4,j4,k4) &
     &    + term5 * Flw(mp)%xix(i5,j5,k5) + term6 * Flw(mp)%xix(i6,j6,k6) &
     &    + term7 * Flw(mp)%xix(i7,j7,k7) + term8 * Flw(mp)%xix(i8,j8,k8)
     xiyp = term1 * Flw(mp)%xiy(i1,j1,k1) + term2 * Flw(mp)%xiy(i2,j2,k2) &
     &    + term3 * Flw(mp)%xiy(i3,j3,k3) + term4 * Flw(mp)%xiy(i4,j4,k4) &
     &    + term5 * Flw(mp)%xiy(i5,j5,k5) + term6 * Flw(mp)%xiy(i6,j6,k6) &
     &    + term7 * Flw(mp)%xiy(i7,j7,k7) + term8 * Flw(mp)%xiy(i8,j8,k8)
     xizp = term1 * Flw(mp)%xiz(i1,j1,k1) + term2 * Flw(mp)%xiz(i2,j2,k2) &
     &    + term3 * Flw(mp)%xiz(i3,j3,k3) + term4 * Flw(mp)%xiz(i4,j4,k4) &
     &    + term5 * Flw(mp)%xiz(i5,j5,k5) + term6 * Flw(mp)%xiz(i6,j6,k6) &
     &    + term7 * Flw(mp)%xiz(i7,j7,k7) + term8 * Flw(mp)%xiz(i8,j8,k8)
     etxp = term1 * Flw(mp)%etx(i1,j1,k1) + term2 * Flw(mp)%etx(i2,j2,k2) &
     &    + term3 * Flw(mp)%etx(i3,j3,k3) + term4 * Flw(mp)%etx(i4,j4,k4) &
     &    + term5 * Flw(mp)%etx(i5,j5,k5) + term6 * Flw(mp)%etx(i6,j6,k6) &
     &    + term7 * Flw(mp)%etx(i7,j7,k7) + term8 * Flw(mp)%etx(i8,j8,k8)
     etyp = term1 * Flw(mp)%ety(i1,j1,k1) + term2 * Flw(mp)%ety(i2,j2,k2) &
     &    + term3 * Flw(mp)%ety(i3,j3,k3) + term4 * Flw(mp)%ety(i4,j4,k4) &
     &    + term5 * Flw(mp)%ety(i5,j5,k5) + term6 * Flw(mp)%ety(i6,j6,k6) &
     &    + term7 * Flw(mp)%ety(i7,j7,k7) + term8 * Flw(mp)%ety(i8,j8,k8)
     etzp = term1 * Flw(mp)%etz(i1,j1,k1) + term2 * Flw(mp)%etz(i2,j2,k2) &
     &    + term3 * Flw(mp)%etz(i3,j3,k3) + term4 * Flw(mp)%etz(i4,j4,k4) &
     &    + term5 * Flw(mp)%etz(i5,j5,k5) + term6 * Flw(mp)%etz(i6,j6,k6) &
     &    + term7 * Flw(mp)%etz(i7,j7,k7) + term8 * Flw(mp)%etz(i8,j8,k8)
     zexp = term1 * Flw(mp)%zex(i1,j1,k1) + term2 * Flw(mp)%zex(i2,j2,k2) &
     &    + term3 * Flw(mp)%zex(i3,j3,k3) + term4 * Flw(mp)%zex(i4,j4,k4) &
     &    + term5 * Flw(mp)%zex(i5,j5,k5) + term6 * Flw(mp)%zex(i6,j6,k6) &
     &    + term7 * Flw(mp)%zex(i7,j7,k7) + term8 * Flw(mp)%zex(i8,j8,k8)
     zeyp = term1 * Flw(mp)%zey(i1,j1,k1) + term2 * Flw(mp)%zey(i2,j2,k2) &
     &    + term3 * Flw(mp)%zey(i3,j3,k3) + term4 * Flw(mp)%zey(i4,j4,k4) &
     &    + term5 * Flw(mp)%zey(i5,j5,k5) + term6 * Flw(mp)%zey(i6,j6,k6) &
     &    + term7 * Flw(mp)%zey(i7,j7,k7) + term8 * Flw(mp)%zey(i8,j8,k8)
     zezp = term1 * Flw(mp)%zez(i1,j1,k1) + term2 * Flw(mp)%zez(i2,j2,k2) &
     &    + term3 * Flw(mp)%zez(i3,j3,k3) + term4 * Flw(mp)%zez(i4,j4,k4) &
     &    + term5 * Flw(mp)%zez(i5,j5,k5) + term6 * Flw(mp)%zez(i6,j6,k6) &
     &    + term7 * Flw(mp)%zez(i7,j7,k7) + term8 * Flw(mp)%zez(i8,j8,k8)
     rhof = term1 * Flw(mp)%rho(i1,j1,k1) + term2 * Flw(mp)%rho(i2,j2,k2) &
     &    + term3 * Flw(mp)%rho(i3,j3,k3) + term4 * Flw(mp)%rho(i4,j4,k4) &
     &    + term5 * Flw(mp)%rho(i5,j5,k5) + term6 * Flw(mp)%rho(i6,j6,k6) &
     &    + term7 * Flw(mp)%rho(i7,j7,k7) + term8 * Flw(mp)%rho(i8,j8,k8)
     uf   = term1 * Flw(mp)%u  (i1,j1,k1) + term2 * Flw(mp)%u  (i2,j2,k2) &
     &    + term3 * Flw(mp)%u  (i3,j3,k3) + term4 * Flw(mp)%u  (i4,j4,k4) &
     &    + term5 * Flw(mp)%u  (i5,j5,k5) + term6 * Flw(mp)%u  (i6,j6,k6) &
     &    + term7 * Flw(mp)%u  (i7,j7,k7) + term8 * Flw(mp)%u  (i8,j8,k8)
     vf   = term1 * Flw(mp)%v  (i1,j1,k1) + term2 * Flw(mp)%v  (i2,j2,k2) &
     &    + term3 * Flw(mp)%v  (i3,j3,k3) + term4 * Flw(mp)%v  (i4,j4,k4) &
     &    + term5 * Flw(mp)%v  (i5,j5,k5) + term6 * Flw(mp)%v  (i6,j6,k6) &
     &    + term7 * Flw(mp)%v  (i7,j7,k7) + term8 * Flw(mp)%v  (i8,j8,k8)
     wf   = term1 * Flw(mp)%w  (i1,j1,k1) + term2 * Flw(mp)%w  (i2,j2,k2) &
     &    + term3 * Flw(mp)%w  (i3,j3,k3) + term4 * Flw(mp)%w  (i4,j4,k4) &
     &    + term5 * Flw(mp)%w  (i5,j5,k5) + term6 * Flw(mp)%w  (i6,j6,k6) &
     &    + term7 * Flw(mp)%w  (i7,j7,k7) + term8 * Flw(mp)%w  (i8,j8,k8)
     muf  = term1 * Flw(mp)%mu (i1,j1,k1) + term2 * Flw(mp)%mu (i2,j2,k2) &
     &    + term3 * Flw(mp)%mu (i3,j3,k3) + term4 * Flw(mp)%mu (i4,j4,k4) &
     &    + term5 * Flw(mp)%mu (i5,j5,k5) + term6 * Flw(mp)%mu (i6,j6,k6) &
     &    + term7 * Flw(mp)%mu (i7,j7,k7) + term8 * Flw(mp)%mu (i8,j8,k8)
    case(2:6)
     xixp = term1 * Flw(mp)%xix(i1,j1,k1) + term2 * (Flw(mp)%xix(i2,j2,k2) - Flw(mp)%xix(i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%xix(i3,j3,k3) - Flw(mp)%xix(i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%xix(i4,j4,k4) - Flw(mp)%xix(i1,j1,k1))
     xiyp = term1 * Flw(mp)%xiy(i1,j1,k1) + term2 * (Flw(mp)%xiy(i2,j2,k2) - Flw(mp)%xiy(i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%xiy(i3,j3,k3) - Flw(mp)%xiy(i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%xiy(i4,j4,k4) - Flw(mp)%xiy(i1,j1,k1))
     xizp = term1 * Flw(mp)%xiz(i1,j1,k1) + term2 * (Flw(mp)%xiz(i2,j2,k2) - Flw(mp)%xiz(i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%xiz(i3,j3,k3) - Flw(mp)%xiz(i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%xiz(i4,j4,k4) - Flw(mp)%xiz(i1,j1,k1))
     etxp = term1 * Flw(mp)%etx(i1,j1,k1) + term2 * (Flw(mp)%etx(i2,j2,k2) - Flw(mp)%etx(i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%etx(i3,j3,k3) - Flw(mp)%etx(i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%etx(i4,j4,k4) - Flw(mp)%etx(i1,j1,k1))
     etyp = term1 * Flw(mp)%ety(i1,j1,k1) + term2 * (Flw(mp)%ety(i2,j2,k2) - Flw(mp)%ety(i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%ety(i3,j3,k3) - Flw(mp)%ety(i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%ety(i4,j4,k4) - Flw(mp)%ety(i1,j1,k1))
     etzp = term1 * Flw(mp)%etz(i1,j1,k1) + term2 * (Flw(mp)%etz(i2,j2,k2) - Flw(mp)%etz(i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%etz(i3,j3,k3) - Flw(mp)%etz(i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%etz(i4,j4,k4) - Flw(mp)%etz(i1,j1,k1))
     zexp = term1 * Flw(mp)%zex(i1,j1,k1) + term2 * (Flw(mp)%zex(i2,j2,k2) - Flw(mp)%zex(i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%zex(i3,j3,k3) - Flw(mp)%zex(i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%zex(i4,j4,k4) - Flw(mp)%zex(i1,j1,k1))
     zeyp = term1 * Flw(mp)%zey(i1,j1,k1) + term2 * (Flw(mp)%zey(i2,j2,k2) - Flw(mp)%zey(i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%zey(i3,j3,k3) - Flw(mp)%zey(i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%zey(i4,j4,k4) - Flw(mp)%zey(i1,j1,k1))
     zezp = term1 * Flw(mp)%zez(i1,j1,k1) + term2 * (Flw(mp)%zez(i2,j2,k2) - Flw(mp)%zez(i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%zez(i3,j3,k3) - Flw(mp)%zez(i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%zez(i4,j4,k4) - Flw(mp)%zez(i1,j1,k1))
     rhof = term1 * Flw(mp)%rho(i1,j1,k1) + term2 * (Flw(mp)%rho(i2,j2,k2) - Flw(mp)%rho(i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%rho(i3,j3,k3) - Flw(mp)%rho(i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%rho(i4,j4,k4) - Flw(mp)%rho(i1,j1,k1))
     uf   = term1 * Flw(mp)%u  (i1,j1,k1) + term2 * (Flw(mp)%u  (i2,j2,k2) - Flw(mp)%u  (i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%u  (i3,j3,k3) - Flw(mp)%u  (i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%u  (i4,j4,k4) - Flw(mp)%u  (i1,j1,k1))
     vf   = term1 * Flw(mp)%v  (i1,j1,k1) + term2 * (Flw(mp)%v  (i2,j2,k2) - Flw(mp)%v  (i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%v  (i3,j3,k3) - Flw(mp)%v  (i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%v  (i4,j4,k4) - Flw(mp)%v  (i1,j1,k1))
     wf   = term1 * Flw(mp)%w  (i1,j1,k1) + term2 * (Flw(mp)%w  (i2,j2,k2) - Flw(mp)%w  (i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%w  (i3,j3,k3) - Flw(mp)%w  (i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%w  (i4,j4,k4) - Flw(mp)%w  (i1,j1,k1))
     muf  = term1 * Flw(mp)%mu (i1,j1,k1) + term2 * (Flw(mp)%mu (i2,j2,k2) - Flw(mp)%mu (i1,j1,k1)) &
     &                                    + term3 * (Flw(mp)%mu (i3,j3,k3) - Flw(mp)%mu (i1,j1,k1)) &
     &                                    + term4 * (Flw(mp)%mu (i4,j4,k4) - Flw(mp)%mu (i1,j1,k1))
    case default
     write(*, '(a)') '!!!!! Error : fInter number !!!!!'
     stop
   end select
   ! 液滴速度 (計算空間) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   call Velocity3D( &
   &      xixp, xiyp, xizp, etxp, etyp, etzp, zexp, zeyp, zezp, up, vp, wp, &
   &      uxp, vep, wzp )
   ! 液滴に働く抗力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   call Drag3D( &
   &      dp, Rhod, up, vp, wp, Rhof, uf, vf, wf, muf, &
   &      fdx, fdy, fdz )
   ! 液滴に働く重力と浮力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   gx = 0.0
   gy = -9.806650
   gz = 0.0
   gx = gx / (aRef**2 / lRef)
   gy = gy / (aRef**2 / lRef)
   gz = gz / (aRef**2 / lRef)
   call GravityBuoyancy3D( &
   &      Rhod, Rhof, gx, gy, gz, &
   &      fgx, fgy, fgz )
   ! 合力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   fx = fdx + fgx; fy = fdy + fgy; fz = fdz + fgz
   ! 処理終了 ********************************************************************************************
   return
  end subroutine SearchDroplet


  subroutine ini_DropletTimeStep(ip,jp,kp,mp)
    ! 変数宣言 ********************************************************************************************
    implicit none
    integer, intent(inout) :: ip, jp, kp
    integer, intent(in)    :: mp
    integer, parameter :: nRetry = 100         ! 液滴探索最大再試行回数
    integer, parameter :: mRetry = 10          ! 液滴補間最大再試行回数
    integer, parameter :: SR     = 2         ! 液滴探索範囲
    real   , parameter :: MGN    = 1.0e-3          ! 補間経緯数許容誤差
    ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    integer :: mp
    real    :: xp0, yp0, zp0, up0, vp0, wp0
    integer :: n, m
    integer :: ip0, jp0, kp0
    real(8) :: rr
    real    :: mr
    logical :: fMove, fImpi, fExit
    ! 処理開始 ********************************************************************************************
    ! 液滴初期位置及び速度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    mp   = mini
    mr   = 1.0
    fSplash = .false.
    fBounce = .false.
   do n = 1, nRetry
    ! 初期位置探索
    do m = 1, mRetry
     call SearchDroplet( &
     &      mp, &
     &      Flw(mp)%is, Flw(mp)%ie - 1, &
     &      Flw(mp)%js, Flw(mp)%je - 1, &
     &      Flw(mp)%ks, Flw(mp)%ke - 1, &
     &      MGN * real(m - 1), ip, jp, kp )
     if(fSearch) exit
    enddo
    if(fSearch) exit
   enddo
   if(.not. fSearch) then
     write(*, '(a)') '!!!!! Error : Droplet initial position !!!!!'
   !  write(*, '(a, i11)') '* Droplet number   = ', nDrop
     write(*, '(a, i11)') '* Droplet number   = ', incount
     write(*, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
     swi_lost = 1
     return
   endif


  ! 処理終了 ********************************************************************************************
  return
  end subroutine ini_DropletTimeStep



!subroutine CalTau(is,ie,js,je,ks,ke,rho,tau,idim,xx,yy,uu,vv)
!!変数宣言
!implicit none
!
!!引数変数
!
!integer,intent(in) :: is,ie,js,je,ks,ke,idim
!real   ,intent(in) :: rho(is:ie,js:je,ks:ke)
!real   ,intent(in) :: xx(is:ie,js:je,ks:ke)    !格子法の格子点のx座標
!real   ,intent(in) :: yy(is:ie,js:je,ks:ke)    !格子法の格子点のy座標
!real   ,intent(in) :: uu(is:ie,js:je,ks:ke)    !格子法の速度u
!real   ,intent(in) :: vv(is:ie,js:je,ks:ke)    !格子法の速度v
!real   ,intent(out), allocatable :: tau(:)
!
!!局所変数
!real    :: utau(is:ie),grid_1(0:idim-1,0:idim-1)
!real    :: vector_u(0:idim-1),inner(0:idim-1)
!integer :: i,j,m
!real    :: grid_sca,u_sca
!integer :: jj1
!real, parameter :: zero = 1.0e-10
!
!!array_allocation
!allocate(tau(is:ie))
!
!open(700,file='./utau.txt')
!read(700,*) utau
!close(700)
!
!
!tau(:)=0.0
!do i=is,ie
!  tau(i) = (utau(i)*aRef)**2 * rho(i,1,1)*rhoRef
!
!  do j = js+1,je
!    if((abs(uu(i,j,ks)) .gt. zero) .and. (abs(vv(i,j,ks)) .gt. zero)) then
!      jj1 = j-1    !壁面上の格子点の座標を選定．jj1+1が第一格子点のj座標
!      exit
!    end if
!  end do
!!write(*,*) jj1   !確認用
!
!!  tau(i) = 10.0
!  if(i .eq. is) then
!    tau(i) = -1.0*tau(i)
!  else if((i .ne. is) .and. (i .ne. ie)) then
!    grid_1(0,0) = xx(i-1,jj1,ks) - xx(i,jj1,ks)     !ここから下は二次元のみ対応．三次元計算の時には変更の必要あり
!    grid_1(1,0) = yy(i-1,jj1,ks) - yy(i,jj1,ks)
!    grid_1(0,1) = xx(i+1,jj1,ks) - xx(i,jj1,ks)
!    grid_1(1,1) = yy(i+1,jj1,ks) - yy(i,jj1,ks)
!    do m = 0,idim-1
!      grid_sca = (grid_1(0,m)**2+grid_1(1,m)**2)
!      grid_1(0,m) = grid_1(0,m)/grid_sca
!      grid_1(1,m) = grid_1(1,m)/grid_sca
!    end do
!    vector_u(0) = uu(i,jj1+1,ks)
!    vector_u(1) = vv(i,jj1+1,ks)
!    u_sca = vector_u(0)**2+vector_u(1)**2
!    vector_u(0) = vector_u(0)/u_sca
!    vector_u(1) = vector_u(1)/u_sca
!    do m=0,idim-1
!      inner(m) = grid_1(0,m)*vector_u(0)+grid_1(1,m)*vector_u(1)
!    end do
!    if(inner(0) .gt. inner(1)) then   !なんかとりあえず向きは合ってそう…？
!      tau(i) = -1.0*tau(i)
!    end if
!  end if
!
!!if(i .lt. ie/2)then
!!tau(i) = -10.0   !追加
!!else if(i .gt. ie/2)then
!!tau(i) = 10.0
!!end if
!end do
!
!!check
!write(*,*)
!
!open(600,file='./tau.txt',status = 'replace')
!write(600,'(f15.8)') tau
!close(600)
!
!return
!end subroutine CalTau

subroutine input_wall(wall)
 real, intent(out) :: wall(0:130,0:1)
 character(200)	:: dummy
 integer	:: i,j
 real,dimension(0:1) :: temp_x
 real,parameter	:: zero = 1.0e-10

 open(1,file = './data/MPSini/wall.vtk',status = 'old')
  do i = 0,4
   read(1,*) dummy
  end do
  i = 0
  do
   read(1,*) temp_x(0),temp_x(1)
   if(i .eq. 0) then
    wall(i,0) = temp_x(0)
    wall(i,1) = temp_x(1)
    i = i + 1
   else
    do j = 0,i-1
     if((abs(wall(j,0)-temp_x(0)) .lt. zero) .and. (abs(wall(j,1)-temp_x(1)) .lt. zero)) exit
    end do
    if(j .eq. i) then
     wall(i,0) = temp_x(0)
     wall(i,1) = temp_x(1)
     i = i + 1
    end if
   end if
   if(i .gt. 130) exit
  end do
 close(1)

end subroutine input_wall

subroutine HeatTransfer(&
      &i,nn,idim,nump,num_surf,ICMsurf, &
      &x,t_temp,sur_temp,dis,hc,Qhc)
implicit none
  integer, intent(in) :: i,nn,idim,nump
  integer, intent(in) :: num_surf
  real   , intent(in) :: ICMsurf(0:num_surf-1,0:1)
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  real   , intent(in) :: t_temp(0:nn-1)
  real   , intent(in) :: sur_temp(0:num_surf-1)
  real   , intent(in) :: dis
  real   , intent(in) :: hc(0:num_surf-1)
  real   , intent(inout) :: Qhc(0:nn-1)
!  integer, intent(in) :: swi_hc(0:nn-1)

  real   , parameter :: pi = acos(-1.0)
  integer :: i1
  real		:: min
  integer	:: imin
  real		:: dista

! if(swi_hc(i) .eq. 1) then
  min = 1.0e10
  do i1 = 0,num_surf-1
   dista = (x(i,0)-ICMsurf(i1,0))**2+(x(i,1)-ICMsurf(i1,1))**2
   if(dista .lt. min) then
    min = dista
    imin = i1
   end if
  end do
  Qhc(i) = - hc(imin) * (t_temp(i) - sur_temp(imin)) * dis / (pi*(dis*0.5)**2)
!
! else
!  Qhc(i) = 0.0
!
! end if

end subroutine HeatTransfer



  ! *****************************************************************************************************
  ! *****************************  ここから先は全くいじってない、分からん  *****************************
  ! *****************************************************************************************************

  !*******************************************************************************************************
  !******** 初期設定                    ********
  !*******************************************************************************************************
  subroutine InitialSetting
   ! 変数宣言 ********************************************************************************************
   implicit none
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer   :: i, j, k, m
   integer   :: jp
   character :: fname * 30
   ! 処理開始 ********************************************************************************************
   ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   allocate( Flw(ms:me), OSG(ms:me), Drp(1:nDrpFile), Ice(ms:me) )
   ! ブロック名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   call Set_BlockName
   ! 格子解像度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do m = ms, me
     call Input_Resolution3D( &
     &      trim(GrdInDir) // trim(BlkName(m)) // trim(RslFile), strtxt, &
     &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke )
   enddo
   ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do m = ms, me
    allocate( Flw(m)%x   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%y   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%z   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%rho (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%u   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%v   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%w   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%p   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%t   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%mu  (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%kin (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%eps (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%xix (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%xiy (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%xiz (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%etx (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%ety (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%etz (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%zex (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%zey (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%zez (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%jac (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%qh  (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le) )
   enddo
   ! 格子座標 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do m = ms, me
    call Input_Grid3D( &
    &      trim(FlwIniDir) // trim(BlkName(m)) // trim(ND_GrdFile), strbin, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
   enddo
   ! C 型格子分割点 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   m = 1
   call Input_CtypeGridPoint( &
   &      trim(GrdInDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
   &      Flw(m)%i1, Flw(m)%i2, Flw(m)%i3 )
   ! メトリックス ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do m = ms, me
    call Input_Metrics3D( &
    &      trim(FlwIniDir) // trim(BlkName(m)) // trim(ND_MetFile), strbin, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      Flw(m)%jac, &
    &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
    &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
    &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez )
   enddo
   ! 流束関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do m = ms, me
    call Input_Flux3D( &
    &      trim(FlwCalInDir) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      ls, le, Flw(m)%qh )
  !  call Input_Flux3D( &
  !  &      trim(FlwCalInDir) // trim(fname) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
  !  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  !  &      ls, le, Flw(m)%qh )
   enddo
   ! 物理量 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do m = ms, me
    call SetPhysics3DKEM( &
    &      Rg, gamma, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
    &      Flw(m)%qh, Flw(m)%jac, &
    &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps )
    call ViscosityCoefficient3D( &
    &      muSth, TsSth, s1, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      Flw(m)%t, &
    &      Flw(m)%mu )
    Flw(m)%w(:,:,:) = 0.0
   enddo
   ! 重合格子補間係数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do m = ms, me
    ! 補間探索範囲 ---------------------------------------------------------------------------------------
    call Input_Resolution3D( &
    &      trim(OSGDir) // trim(BlkName(m)) // trim(OverAreaFile), strtxt, &
    &      OSG(m)%is, OSG(m)%ie, OSG(m)%js, OSG(m)%je, OSG(m)%ks, OSG(m)%ke )
    ! メモリ確保 -----------------------------------------------------------------------------------------
    allocate( OSG(m)%ip   (OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
    &         OSG(m)%jp   (OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
    &         OSG(m)%kp   (OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
    &         OSG(m)%fOver(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
    &         OSG(m)%term1(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
    &         OSG(m)%term2(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
    &         OSG(m)%term3(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
    &         OSG(m)%term4(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
    &         OSG(m)%term5(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
    &         OSG(m)%term6(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
    &         OSG(m)%term7(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
    &         OSG(m)%term8(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke) )
    ! 補間係数 -------------------------------------------------------------------------------------------
    call Input_OversetCoe3D( &
    &      trim(OSGDir) // trim(BlkName(m)) // trim(OverCoeFile), strbin, &
    &      OSG(m)%is, OSG(m)%ie, OSG(m)%js, OSG(m)%je, OSG(m)%ks, OSG(m)%ke, &
    &      OSG(m)%fOver, OSG(m)%ip, OSG(m)%jp, OSG(m)%kp, &
    &      OSG(m)%term1, OSG(m)%term2, OSG(m)%term3, OSG(m)%term4, &
    &      OSG(m)%term5, OSG(m)%term6, OSG(m)%term7, OSG(m)%term8 )
   enddo
   ! 重合格子部のフラグ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   m = 1
   allocate( fOver(Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke) )
   fOver(:,:,:) = .false.
   do k = OSG(m)%ks, OSG(m)%ke
   do j = OSG(m)%js, OSG(m)%je
   do i = OSG(m)%is, OSG(m)%ie
    if( OSG(m)%fOver(i,j,k) == 0 ) cycle
    fOver(i,j,k) = .true.
   enddo
   enddo
   enddo
   ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do m = ms, me
    !deallocate( Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps, Flw(m)%jac, Flw(m)%qh )
    deallocate( Flw(m)%jac, Flw(m)%qh )
    deallocate( OSG(m)%fOver, OSG(m)%ip   , OSG(m)%jp   , OSG(m)%kp   , &
    &           OSG(m)%term1, OSG(m)%term2, OSG(m)%term3, OSG(m)%term4, &
    &           OSG(m)%term5, OSG(m)%term6, OSG(m)%term7, OSG(m)%term8 )
   enddo
   ! 液滴軌道 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do m = 1, nDrpFile
    allocate( Drp(m)%m(1:nCalMax), &
    &         Drp(m)%x(1:nCalMax), Drp(m)%y(1:nCalMax), Drp(m)%z(1:nCalMax), &
    &         Drp(m)%u(1:nCalMax), Drp(m)%v(1:nCalMax), Drp(m)%w(1:nCalMax), &
    &         Drp(m)%mr(1:nCalMax)  )
    Drp(m)%m(:) = 1
    Drp(m)%u(:) = 0.0; Drp(m)%v(:) = 0.0; Drp(m)%w(:) = 0.0
    Drp(m)%mr(:) = 1.0
   enddo
   ! 着氷計算領域 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do m = ms, me
    call Input_Resolution1D( &
    &      trim(GrdInDir) // trim(BlkName(m)) // trim(IceRslFile), strtxt, &
    &      Ice(m)%is, Ice(m)%ie )
   enddo
   ! 衝突セル ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do m = ms, me
    ! メモリ確保
    allocate( Flw(m)%nimp(Ice(m)%is: Ice(m)%ie, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%uimp(Ice(m)%is: Ice(m)%ie, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%vimp(Ice(m)%is: Ice(m)%ie, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%wimp(Ice(m)%is: Ice(m)%ie, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%simp(Ice(m)%is: Ice(m)%ie, Flw(m)%ks: Flw(m)%ke), &
    &         Flw(m)%CE  (Ice(m)%is: Ice(m)%ie, Flw(m)%ks: Flw(m)%ke) )
    Flw(m)%nimp(:,:) = 0
    Flw(m)%uimp(:,:) = 0.0; Flw(m)%vimp(:,:) = 0.0; Flw(m)%wimp(:,:) = 0.0; Flw(m)%simp(:,:) = 0
   enddo
   ! 液滴計算途中解 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if( nDrop > nDrpOutputCount .and. CalNum /= 2 ) then
     write(fname, '(a,i8.8,a)') 'Drop', nDrop, '_'
     do m = ms, me
      call Input_Impingement3D( &
      &      trim(DrpImpDir) // trim(BlkName(m)) // trim(ND_DrpImpiDataFile), strbin, &
      &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, &
      &      Flw(m)%nimp, Flw(m)%uimp, Flw(m)%vimp, Flw(m)%wimp )
  !    call Input_Impingement3D( &
  !    &      trim(DrpImpDir) // trim(fname) // trim(BlkName(m)) // trim(ND_DrpImpiDataFile), strbin, &
  !    &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, &
  !    &      Flw(m)%nimp, Flw(m)%uimp, Flw(m)%vimp, Flw(m)%wimp )
      nStart = nDrop
     enddo
    else
     nStart = 1
   endif
   ! 処理終了 ********************************************************************************************
   return
  end subroutine InitialSetting


  !*******************************************************************************************************
  !******** メモリ解放                    ********
  !*******************************************************************************************************
  subroutine Deallocating
   ! 変数宣言 ********************************************************************************************
   implicit none
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer :: m
   ! 処理開始 ********************************************************************************************
   do m = ms, me
    deallocate( Flw(m)%x, Flw(m)%y, Flw(m)%z, &
    &           Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%mu, &
    &           Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
    &           Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
    &           Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
    &           Flw(m)%nimp, Flw(m)%uimp, Flw(m)%vimp, Flw(m)%wimp, Flw(m)%CE )
    deallocate( Drp(m)%m, Drp(m)%x, Drp(m)%y, Drp(m)%z, Drp(m)%u, Drp(m)%v, Drp(m)%w, Drp(m)%mr )
   enddo
   ! 処理終了 ********************************************************************************************
   return
  end subroutine Deallocating
  !*******************************************************************************************************
  !******** Main Grid から Sub Grid への移動              ********
  !*******************************************************************************************************
  subroutine BoundaryMGtoSG( ip, jp, kp, mp, fMove )
   ! 変数宣言 ********************************************************************************************
   implicit none
   ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer, intent(in)    :: ip, jp, kp
   integer, intent(inout) :: mp
   logical, intent(out)   :: fMove
   ! 処理開始 ********************************************************************************************
   ! 例外処理 --------------------------------------------------------------------------------------------
   if( ip == Flw(mp)%ie .or. jp == Flw(mp)%je .or. kp == Flw(mp)%ke ) then
    fMove = .false.
    return
   endif
   ! ブロック移動判定 ------------------------------------------------------------------------------------
   if( fOver(ip  , jp  , kp) .and. fOver(ip+1, jp  , kp) .and. &
   &   fOver(ip  , jp+1, kp) .and. fOver(ip+1, jp+1, kp) ) then
     fMove = .true.
     mp    = 2
    else
     fMove = .false.
     mp    = 1
   endif
   ! 処理終了 ********************************************************************************************
   return
  end subroutine BoundaryMGtoSG
  !*******************************************************************************************************
  !******** Sub Grid から Main Grid への移動              ********
  !*******************************************************************************************************
  subroutine BoundarySGtoMG( ip, jp, kp, mp, fMove )
   ! 変数宣言 ********************************************************************************************
   implicit none
   ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer, intent(in)    :: ip, jp, kp
   integer, intent(inout) :: mp
   logical, intent(out)   :: fMove
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   logical :: isOut, ieOut, jeOut
   ! 処理開始 ********************************************************************************************
   isOut = ip <  Flw(mp)%is + 1 .and. uxp < 0.0
   ieOut = ip >= Flw(mp)%ie - 1 .and. uxp > 0.0
   jeOut = jp >= Flw(mp)%je - 1 .and. vep > 0.0
   if( isOut .or. ieOut .or. jeOut ) then
     fMove = .true.
     mp    = 2 !1
    else
     fMove = .false.
     mp    = 2
   endif
   ! 処理終了 ********************************************************************************************
   return
  end subroutine BoundarySGtoMG
  !*******************************************************************************************************
  !******** 翼表面の衝突判定                  ********
  !*******************************************************************************************************
  subroutine BoundaryBladeSurface( ip, jp, kp, j0, mp, fImpi, up, vp, wp, mr )
   ! 変数宣言 ********************************************************************************************
   implicit none
   ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer, intent(in)    :: ip, jp, kp
   integer, intent(in)    :: j0
   integer, intent(in)    :: mp
   logical, intent(out)   :: fImpi
   real   , intent(inout) :: up, vp, wp, mr
   ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer, parameter :: jMin   = 20
   integer, parameter :: nRetry = 100
   real   , parameter :: MGN    = 0.01
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   real    :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
   real    :: alp, bet, gam
   integer :: n
   integer :: i1, i2, i3, i4, j1, j2, j3, j4, k1, k2, k3, k4
   logical :: fDivi
   real    :: wgh1, wgh2, wgh3, wgh4, wgh5, wgh6, wgh7, wgh8
   real    :: r1, r2, r3, r4
   real    :: u1, v1, w1, d1, m1, n1, u2, v2, w2, d2, m2, n2
   real    :: SLDLim
   ! 処理開始 ********************************************************************************************
   ! 例外処理 --------------------------------------------------------------------------------------------
   if( ip == Flw(mp)%ie .or. jp == Flw(mp)%je .or. kp == Flw(mp)%ke ) then
     fImpi = .false.
     return
   endif
   if( (mp == 1) .and. (ip < Flw(mp)%i1 .or. ip >= Flw(mp)%i3) ) then
     fImpi = .false.
     return
   endif
   if( jp > jMin ) then
     fImpi = .false.
     return
   endif
   ! 衝突判定 --------------------------------------------------------------------------------------------
   ! 周囲点のインデックス
   i1 = ip; k1 = kp; i2 = ip; k2 = kp + 1; i3 = ip + 1; k3 = kp
   ! 周囲 3点の座標
   x1 = Flw(mp)%x(i1,j0,k1); y1 = Flw(mp)%y(i1,j0,k1); z1 = Flw(mp)%z(i1,j0,k1)
   x2 = Flw(mp)%x(i2,j0,k2); y2 = Flw(mp)%y(i2,j0,k2); z2 = Flw(mp)%z(i2,j0,k2)
   x3 = Flw(mp)%x(i3,j0,k3); y3 = Flw(mp)%y(i3,j0,k3); z3 = Flw(mp)%z(i3,j0,k3)
   ! 衝突前の物理量
   u1 = up; v1 = vp; w1 = wp; d1 = dp; m1 = mr
   n1 = nDrpSpl + nDrpBou
   SLDLim = SplLim / lRef
   ! 衝突判定
  ! call ImpingementJudge3D( &
  ! &      xp, yp, zp, x1, y1, z1, x2, y2, z2, x3, y3, z3, MVD, &
  ! &      fImpi )
  ! m2 = 0.0
   call SplashLEWICE3D( &
   &      xp, yp, zp, u1, v1, w1, d1, m1, x1, y1, z1, x2, y2, z2, x3, y3, z3, &
   &      LWC, Rhod, Sigd, mud, SLDLim, &
   &      fImpi, u2, v2, w2, d2, m2, nDrpSpl, nDrpBou )
   ! 衝突後の物理量
   up = u2; vp = v2; wp = w2; dp = d2; mr = m2
   n2 = nDrpSpl + nDrpBou
   ! スプラッシュ・バウンドのログ
   if(n1 < n2 ) then
     if(d1 > d2 )then
  !     write(*, '(a)') '----- Splash -----'
  !     write(*, '(a, i11)') '* Droplet number   = ', nDrop
  !     write(*, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
  !     write(9, '(a)') '----- Splash -----'
  !     write(9, '(a, i11)') '* Droplet number   = ', nDrop
  !     write(9, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
       fSplash = .true.
      else
  !     write(*, '(a)') '----- Bounce -----'
  !     write(*, '(a, i11)') '* Droplet number   = ', nDrop
  !     write(*, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
  !     write(9, '(a)') '----- Bounce -----'
  !     write(9, '(a, i11)') '* Droplet number   = ', nDrop
  !     write(9, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
       fBounce = .true.
     endif
    else
     if(fImpi) then
  !     write(*, '(a)') '----- Impingement -----'
  !     write(*, '(a, i11)') '* Droplet number   = ', nDrop
  !     write(*, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
  !     write(9, '(a)') '----- Impingement -----'
  !     write(9, '(a, i11)') '* Droplet number   = ', nDrop
  !     write(9, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
     endif
   endif
   if(.not. fImpi) return
   ! 衝突量の分配 ----------------------------------------------------------------------------------------
   ! 周囲 4点と液滴の距離
   i1 = ip; k1 = kp; i2 = ip; k2 = kp + 1; i3 = ip + 1; k3 = kp; i4 = ip + 1; k4 = kp + 1
   x1 = Flw(mp)%x(i1,j0,k1); y1 = Flw(mp)%y(i1,j0,k1); z1 = Flw(mp)%z(i1,j0,k1)
   x2 = Flw(mp)%x(i2,j0,k2); y2 = Flw(mp)%y(i2,j0,k2); z2 = Flw(mp)%z(i2,j0,k2)
   x3 = Flw(mp)%x(i3,j0,k3); y3 = Flw(mp)%y(i3,j0,k3); z3 = Flw(mp)%z(i3,j0,k3)
   x4 = Flw(mp)%x(i4,j0,k4); y4 = Flw(mp)%y(i4,j0,k4); z4 = Flw(mp)%z(i4,j0,k4)
   r1 = sqrt( (x1 - xp)**2 + (y1 - yp)**2 + (z1 - zp)**2 )
   r2 = sqrt( (x2 - xp)**2 + (y2 - yp)**2 + (z2 - zp)**2 )
   r3 = sqrt( (x3 - xp)**2 + (y3 - yp)**2 + (z3 - zp)**2 )
   r4 = sqrt( (x4 - xp)**2 + (y4 - yp)**2 + (z4 - zp)**2 )
   ! 重み関数
   if( r1 == 0.0 ) then
     wgh1 = 1.0; wgh2 = 0.0; wgh3 = 0.0; wgh4 = 0.0
    else if( r2 == 0.0 ) then
     wgh1 = 0.0; wgh2 = 1.0; wgh3 = 0.0; wgh4 = 0.0
    else if( r3 == 0.0 ) then
     wgh1 = 0.0; wgh2 = 0.0; wgh3 = 1.0; wgh4 = 0.0
    else if( r4 == 0.0 ) then
     wgh1 = 0.0; wgh2 = 0.0; wgh3 = 0.0; wgh4 = 1.0
    else
     alp  = (r2 * r3 * r4 + r1 * r3 * r4 + r1 * r2 * r4 + r1 * r2 * r3) / (r1 * r2 * r3 * r4)
     wgh1 = 1.0 / (r1 * alp)
     wgh2 = 1.0 / (r2 * alp)
     wgh3 = 1.0 / (r3 * alp)
     wgh4 = 1.0 / (r4 * alp)
   endif
   ! 衝突量
   Flw(mp)%nimp(i1,k1) = Flw(mp)%nimp(i1,k1) + (m1 - m2) * wgh1
   Flw(mp)%uimp(i1,k1) = Flw(mp)%uimp(i1,k1) + u1        * wgh1
   Flw(mp)%vimp(i1,k1) = Flw(mp)%vimp(i1,k1) + v1        * wgh1
   Flw(mp)%wimp(i1,k1) = Flw(mp)%wimp(i1,k1) + w1        * wgh1
   Flw(mp)%nimp(i2,k2) = Flw(mp)%nimp(i2,k2) + (m1 - m2) * wgh2
   Flw(mp)%uimp(i2,k2) = Flw(mp)%uimp(i2,k2) + u1        * wgh2
   Flw(mp)%vimp(i2,k2) = Flw(mp)%vimp(i2,k2) + v1        * wgh2
   Flw(mp)%wimp(i2,k2) = Flw(mp)%wimp(i2,k2) + w1        * wgh2
   Flw(mp)%nimp(i3,k3) = Flw(mp)%nimp(i3,k3) + (m1 - m2) * wgh3
   Flw(mp)%uimp(i3,k3) = Flw(mp)%uimp(i3,k3) + u1        * wgh3
   Flw(mp)%vimp(i3,k3) = Flw(mp)%vimp(i3,k3) + v1        * wgh3
   Flw(mp)%wimp(i3,k3) = Flw(mp)%wimp(i3,k3) + w1        * wgh3
   Flw(mp)%nimp(i4,k4) = Flw(mp)%nimp(i4,k4) + (m1 - m2) * wgh4
   Flw(mp)%uimp(i4,k4) = Flw(mp)%uimp(i4,k4) + u1        * wgh4
   Flw(mp)%vimp(i4,k4) = Flw(mp)%vimp(i4,k4) + v1        * wgh4
   Flw(mp)%wimp(i4,k4) = Flw(mp)%wimp(i4,k4) + w1        * wgh4
   ! スプラッシュの場合は計算続行
   if(n1 < n2 ) fImpi = .false.
   ! 処理終了 ********************************************************************************************
   return
  end subroutine BoundaryBladeSurface
  !*******************************************************************************************************
  !******** C 型格子ブランチ・カット                ********
  !*******************************************************************************************************
  subroutine BoundaryCtypeBranchCut( ip, jp, mp, fMove )
   ! 変数宣言 ********************************************************************************************
   implicit none
   ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer, intent(in)  :: ip, jp
   integer, intent(in)  :: mp
   logical, intent(out) :: fMove
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   logical :: jsOut
   ! 処理開始 ********************************************************************************************
   ! 例外処理 --------------------------------------------------------------------------------------------
   if( ip >= Flw(mp)%i1 .and. ip < Flw(mp)%i3 ) return
   ! ブロック移動判定 ------------------------------------------------------------------------------------
   jsOut = jp < Flw(mp)%js + 1 .and. vep < 0.0
   if(jsOut) then
     xp    = xp + up / abs(vep) * 2.0
     yp    = yp + vp / abs(vep) * 2.0
     zp    = zp + wp / abs(vep) * 2.0
     fMove = .true.
    else
     fMove = .false.
   endif
   ! 処理終了 ********************************************************************************************
   return
  end subroutine BoundaryCtypeBranchCut
  !*******************************************************************************************************
  !******** 周期境界                    ********
  !*******************************************************************************************************
  subroutine BoundaryPeriodic( ip, jp, kp, mp, fMove )
   ! 変数宣言 ********************************************************************************************
   implicit none
   ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer, intent(in)  :: ip, jp, kp
   integer, intent(in)  :: mp
   logical, intent(out) :: fMove
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   real    :: dz
   logical :: ksOut, keOut
   ! 処理開始 ********************************************************************************************
   ksOut = kp <  Flw(mp)%ks + 1 .and. wzp < 0.0
   keOut = kp >= Flw(mp)%ke - 1 .and. wzp > 0.0
   dz    = maxval( Flw(mp)%z(ip,jp,:Flw(mp)%ke-1) )  - minval( Flw(mp)%z(ip,jp,:) )
   if(ksOut) then
     zp = zp + dz
     fMove = .true.
     return
   endif
   if(keOut) then
     zp = zp - dz
     fMove = .true.
     return
   endif
   fMove = .false.
   ! 処理終了 ********************************************************************************************
   return
  end subroutine BoundaryPeriodic
  !*******************************************************************************************************
  !******** 流出境界                    ********
  !*******************************************************************************************************
  subroutine BoundaryOutlet( ip, jp, mp, fExit )
   ! 変数宣言 ********************************************************************************************
   implicit none
   ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer, intent(in)  :: ip, jp
   integer, intent(in)  :: mp
   logical, intent(out) :: fExit
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer :: j0
   integer :: ibs, ibe
   logical :: TEOut, isOut, ieOut, jeOut
   ! 処理開始 ********************************************************************************************
   ! ξ方向 ----------------------------------------------------------------------------------------------
   j0 = Flw(1)%js
   ! T.E. 通過
   ibs = Flw(1)%i1 + 1; ibe = Flw(1)%i3 - 1
   TEOut = xp > maxval( Flw(1)%x(ibs:ibe, j0, :) )
   ! 計算領域
   isOut = ip <  Flw(mp)%is + 1 .and. uxp < 0.0
   ieOut = ip >= Flw(mp)%ie - 1 .and. uxp > 0.0
   ! 流出判定
   if( TEOut .or. isOut .or. ieOut ) then
     fExit = .true.
     return
    else
     fExit = .false.
   endif
   ! η方向 ----------------------------------------------------------------------------------------------
   jeOut = jp >= Flw(mp)%je - 1 .and. vep > 0.0
   ! 流出判定
   if( jeOut ) then
     fExit = .true.
     return
    else
     fExit = .false.
   endif
   ! 処理終了 ********************************************************************************************
   return
  end subroutine BoundaryOutlet





end program main
