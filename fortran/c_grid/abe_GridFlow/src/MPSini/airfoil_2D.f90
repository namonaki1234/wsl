!*********************************************************************
!****** mps calculation dam initialfile maker                  ***
!****** ver. 2014.08.12                                            ***
!*********************************************************************
program main
  ! module_setting ***************************************************
  use mod_is
  use mod_is_dim
  use mod_is_airfoil_dim
  implicit none
  ! array_definition *************************************************
  integer, allocatable :: type(:)
  real   , allocatable :: den(:)
  real   , allocatable :: vis(:)
  real   , allocatable :: t(:)
  integer, allocatable :: f_type(:)
  real   , allocatable :: grid(:,:,:)
  real   , allocatable :: grid_dist(:,:)
  integer, allocatable :: near_wall(:,:)
  real   , allocatable :: x(:,:)
  integer, allocatable :: out_type(:)
  real   , allocatable :: dom_size(:,:)
  integer, allocatable :: grid_n(:)
  integer, allocatable :: itypep(:)
  real   , allocatable :: non_cal(:)
  real   , allocatable :: ker_c(:)
  real   , allocatable :: wall_anc(:,:,:)
  real   , allocatable :: wall_out(:,:)
  real   , allocatable :: wall_nor(:,:)

  ! variable_difinition **********************************************
  integer :: i
  integer :: nump
  real    :: delta
  real    :: kernel_max
  integer :: n
  integer :: grid_num
  real    :: interval
  integer :: switch
  integer :: wall_n

  ! variable_definition(with_initial_setting) ************************
  integer, parameter :: idim         = 2

  integer, parameter :: nn           = 1000000 !600000
  integer, parameter :: ifluid       = 6
  integer, parameter :: nn_grid      = 3000
  integer, parameter :: neighbor     = 1000
  integer, parameter :: nwb          = 3
  integer, parameter :: ite          = 0
  real   , parameter :: time_sim     = 0.0
  integer, parameter :: output_limit = 999
  real   , parameter :: pi           = acos(0.0)*2.0
  real   , parameter :: grav         = 9.8065
  real   , parameter :: dis          = 100.0*1.0e-6 !1.35*1.0e-5
  real   , parameter :: dis_rat      = 0.5
  real   , parameter :: dis_min      = 0.9
  real   , parameter :: cfl          = 0.2
  real   , parameter :: dt           = 1.0e-6
  real   , parameter :: dtout        = 200.0e-6
  integer, parameter :: ite_max      = 10000000
  real   , parameter :: time_max     = 10.0
  real   , parameter :: c_vel        = 15.0
  real   , parameter :: sigma        = 7.274e-2
  real   , parameter :: dist_min     = 0.7
  real   , parameter :: coll_rat     = 0.2
  real   , parameter :: wall_temp    = 80.0
  real   , parameter :: fluid_temp   = 0.0
  real   , parameter :: margin       = 3.0
  character(len = 40), parameter :: stl_wall = '../data/NACA0012.dat'
!  character(len = 40), parameter :: stl_part = '../data/dam_par.stl'
  
  ! 翼設定用のパラメータ
  
  integer :: data_n
  real,parameter  :: wscale = 0.267 !0.53 !0.267      ! 翼の倍率変更、1だと翼弦長1m
  
  ! いままでの翼の基準位置
  ! w/ aoa   w/o aoa
  ! x=0.03   x=0.03
  ! y=0.04   y=0.35
  ! 翼の基準位置の変更
!  real,parameter  :: xpos   = 0.03
  real,parameter  :: xpos   = 0.0
  real,parameter  :: ypos   = 0.0
  real  :: zmin
  real  :: zmax
  real,parameter  :: aoa     = 0.0 !4.0  
  real, allocatable     :: xdata(:)
  real, allocatable     :: ydata(:)
  real, allocatable     :: cent(:)
  
  ! 時間用パラメータ
  integer              :: date_time(8)
  character*10         :: date_time_cha(3)
  
  
  ! output_setting ***************************************************
  open(900,file='./log/MPSini/ini_log.txt')
  
  ! こいつだけ1からなのに注意
  call date_and_time(date_time_cha(1),date_time_cha(2),date_time_cha(3), &
  &                  date_time)
  write(*,*)
  write(*,'(A,I4,A,I2,A,I2)')"DATE = ",date_time(1),"/",date_time(2),&
  &                          "/",date_time(3)
  write(*,'(A,I2,A,I2,A,I2)')"TIME = ",date_time(5),":",date_time(6),&
  &                          ":",date_time(7)
  write(*,*)
  
  
  ! allocation *******************************************************
  call allocation                   &
  &    (nn,idim,nump,ifluid,        &
  &     type,den,vis,f_type,itypep, &
  &     dom_size,grid_n,non_cal,    &
  &     ker_c,x,out_type,t,wall_n)
  
  
  ! parameter_setting ************************************************
  ! initialize
  do n=0,ifluid
   type(n) = n
   den(n)  = 998.2
   vis(n)  = 1.004e-4
  end do

  ! liguid_phase
  den(4)  = 998.2
  vis(4)  = 1.004e-4

  ! gas_phase
  den(5)  = 1.205
  vis(5)  = 15.12e-6

  ! solid_phase
  den(3)  = 2.0*den(4)
  vis(3)  = 0.0

  ! calculation_parameter
  nump = 0
  ker_c(0:3) = real(nwb)+0.1
  delta      = 1.0/(0.5*ker_c(3))
  interval = dis*dis_rat
  
  ! 計算領域の設定
  ! しばらく手打ちで！
  ! 最小x座標
  
  ! データ数をとりあえず手打ちで！
  data_n = 132
  
  ! 最小X座標
  DOM_SIZE(0,0) = -0.03 !-0.06 !0.01
  ! 最大X座標
!  DOM_SIZE(0,1) = 0.2
  DOM_SIZE(0,1) = 0.09 !0.06!0.105
  ! 最小Y座標
  DOM_SIZE(1,0) = -0.025 !-0.05!0.0 !0.02
  ! 最大Y座標
!  DOM_SIZE(1,1) = 0.07
  DOM_SIZE(1,1) = 0.025 !0.05!0.06 !0.05
  
!  ! 最小Z座標
!  DOM_SIZE(2,0) = 0.0
!  ! 最大Z座標
!!  DOM_SIZE(2,1) = DIS * 50.0
!  DOM_SIZE(2,1) = DIS * 30.0
  
  ZMIN   = 0.0
  ZMAX   = 0.0
  
  allocate(cent(0:idim-1))
  allocate(xdata(1:data_n))
  allocate(ydata(1:data_n))
 
  xdata(1:data_n) = 0.0
  ydata(1:data_n) = 0.0
  cent(0:idim-1) = 0.0
  
      call change_airfoil_param &
   &          (data_n, xdata, ydata, &
   &           wscale, xpos, ypos,  &
   &           cent,idim, aoa)
   
   
     call make_airfoil_stl &
   &          (data_n, xdata, ydata, &
   &           zmin, zmax, &
   &           cent, idim)
  
!  ! read_stl_file(particle_file) *************************************
!  call read_stl_file                                    &
!  &    (idim,wall_n,wall_nor,wall_anc,stl_part,wall_out)
!
!  ! particle_positioning *********************************************
!  call particle_positioning         &
!  &    (idim,nn,nump,x,dis,         &
!  &     type(4),interval,wall_temp, &
!  &     itypep,dis_min,t,           &
!  &     wall_anc,wall_nor,          &
!  &     wall_n)
  
  ! 2次元計算用 ********************************************************
    wall_n = 131
    ALLOCATE(WALL_ANC(0:IDIM-1,0:WALL_N-1,0:IDIM-1))
    ALLOCATE(WALL_NOR(0:IDIM-1,0:WALL_N-1))
    ALLOCATE(WALL_OUT(0:IDIM-1,0:WALL_N-1))
    
    WALL_ANC(0:IDIM-1,0:WALL_N-1,0:IDIM-1) = 0.0
    WALL_OUT(0:IDIM-1,0:WALL_N-1) = 0.0
    WALL_NOR(0:IDIM-1,0:WALL_N-1) = 0.0
    
    call make_airfoil_wall &
    &    (data_n, xdata, ydata, aoa,idim, &
    &     wall_n,wall_anc,wall_out,wall_nor,dis)
  
  ! 3次元計算用*********************************************************
!    ! read_stl_file(wall_file) *****************************************
!    call read_stl_file                                    &
!    &    (idim,wall_n,wall_nor,wall_anc,stl_wall,wall_out)
!    
!    ! wall_normal_vector ***********************************************
!    call wall_setting                               &
!    &    (idim,wall_n,wall_anc,wall_out,wall_nor,dis)
  
  ! ここまで************************************************************
  
  
  ! grid_and_particle_interference ***********************************
  call grid_and_particle                 &
  &    (idim,nn,ifluid,wall_n,wall_anc,  &
  &     nump,dis,type,x,itypep,wall_out, &
  &     wall_nor,dis_min)

  ! flag_setting *****************************************************
  call flag_setting                      &
  &    (nn,nump,ifluid,itypep,type,f_type)
  
  ! grid_position_calculation ****************************************
  ! 2次元計算用 ******************************************************
    call grid_setting_2d         &
    &    (nn,idim,nump,dis,x,    &
    &     grid_n,grid_num,ker_c, &
    &     grid,dom_size)
  ! 3次元計算用 ******************************************************
!    call grid_setting_3d         &
!    &    (nn,idim,nump,dis,x,    &
!    &     grid_n,grid_num,ker_c, &
!    &     grid,dom_size)
  ! ここまで********************************
  
  ! grid_information_calculation *************************************
  call grid_information                             &
  &    (idim,wall_n,wall_anc,wall_out,wall_nor,dis, &
  &     grid_dist,near_wall,grid,grid_num,ker_c)

  ! paraview_output **************************************************
  call output_para                                  &
  &    (nn,idim,ifluid,nump,x,t,out_type,itypep,type)
  
  call output_para_wall              &
  &    (idim,grid_num,wall_n,wall_anc)
!  call output_para_grid                    &                  !コンパイルするときにtracebackを入れると回るが入れないと回らない．なんで？？
!  &    (idim,grid_num,wall_n,grid,grid_dist,near_wall)        !上記，コア数を1にすると解決する．また，デフォルトだと壁付近のグリッドしか出力しない

  ! initial_file_output **********************************************
  call output_initial_file         &
  &    (idim,nn,nump,ite,time_sim, &
  &     type,ifluid,itypep,        &
  &     den,vis,grav,delta,sigma,  &
  &     time_max,ite_max,          &
  &     dt,dtout,cfl,near_wall,    &
  &     dis,ker_c,                 &
  &     x,dom_size,non_cal,        &
  &     grid,grid_n,grid_num,      &
  &     dis_rat,c_vel,wall_nor,    &
  &     f_type,output_limit,       &
  &     nn_grid,neighbor,t,        &
  &     dist_min,coll_rat,         &
  &     wall_n,wall_anc,grid_dist)

  write(*,'(a    )')'************** calculation end **************'
  write(*,'(a,i21)')'total_particle_number : ',nump
  write(*,'(a,i21)')'total_grid_number     : ',grid_num
  write(*,'(a,i21)')'total_wall_number     : ',wall_n
  write(*,'(a,a34)')'wall     : ',trim(adjustl(stl_wall))
!  write(*,'(a,a34)')'particle : ',trim(adjustl(stl_part))
  write(*,'(a    )')' '
  
  ! 時間の出力
  ! こいつだけ1からなのに注意
  call date_and_time(date_time_cha(1),date_time_cha(2),date_time_cha(3), &
  &                  date_time)
  write(*,*)
  write(*,'(A,I4,A,I2,A,I2)')"DATE = ",date_time(1),"/",date_time(2),&
  &                          "/",date_time(3)
  write(*,'(A,I2,A,I2,A,I2)')"TIME = ",date_time(5),":",date_time(6),&
  &                          ":",date_time(7)
  write(*,*)
  

  stop
end program main
