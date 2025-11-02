!*********************************************************************
!****** module_for_MPS_calculation(2D)                             ***
!****** ver. 2014.09.16                                            ***
!*********************************************************************
module mod_mps_dim
  implicit none
  ! subroutine_definition ********************************************
  public :: grid_division
  public :: arrange_particle
  public :: surface_tension_term

  ! function_definition ***********************************************
  public :: dist_2D
  public :: nan_2D
contains

!*********************************************************************
!*** correction_weight_function_for_2D                             ***
!*********************************************************************
subroutine correction_weight_function_for_viscosity      &
&          (i,drr,dist,nn,idim,x,grid_num,wall_n,grid_p, &
&           dis,ww_w,nn_wall,g_wall,wall_anc,wall_nor)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i,nn
  integer, intent(in) :: idim
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  integer, intent(in) :: grid_num
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: dis
  integer, intent(in) :: nn_wall
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  real   , intent(inout) :: ww_w
  real   , intent(inout) :: dist
  real   , intent(inout) :: drr(0:idim-1)

  ! subroutine_variable **********************************************
  integer :: j,k,l,m,n
  real    :: kappa_all
  real    :: ww_res
  real    :: interval
  real    :: rr,rr0
  real    :: drr_all
  real    :: drr_c
  real, parameter   :: pi = acos(-1.0)
  real    :: x_sub(0:idim-1,0:3)
  real    :: x0(0:idim-1)
  real    :: w_sub(0:3)
  real    :: kappa(0:idim-1)
  real    :: drr0(0:idim-1)

  ! position_setting  ************************************************
  interval = dis

  x_sub(0,0) = x(i,0)-interval
  x_sub(1,0) = x(i,1)

  x_sub(0,1) = x(i,0)+interval
  x_sub(1,1) = x(i,1)

  x_sub(0,2) = x(i,0)
  x_sub(1,2) = x(i,1)-interval

  x_sub(0,3) = x(i,0)
  x_sub(1,3) = x(i,1)+interval

  ! calculation_ww_each_point ****************************************
  do k=0,3
    w_sub(k) = 10000.0
    do n=1,g_wall(grid_p(i),0)
      l = g_wall(grid_p(i),n)
      do m=0,idim-1
        x0(m) = x_sub(m,k)
      end do

      call calc_dist                             &
      &    (x0,l,                                &
      &     drr0,rr,idim,wall_n,wall_anc,wall_nor)

      w_sub(k) = min(rr,w_sub(k))
    end do
  end do

  rr = 10000.0
  do n=1,g_wall(grid_p(i),0)
    l = g_wall(grid_p(i),n)
    do m=0,idim-1
      x0(m) = x(i,m)
    end do

    call calc_dist                &
    &    (x0,l,                   &
    &     drr0,rr0,idim,wall_n,wall_anc,wall_nor)

    rr = min(rr,rr0)
  end do

  ! correction_ww *****************************************************
  ww_res = ww_w

  kappa(0) = (w_sub(0)-2.0*rr+w_sub(1))/interval
  kappa(1) = (w_sub(2)-2.0*rr+w_sub(3))/interval

  kappa_all = (kappa(0)+kappa(1))/2.0

  ww_w = (2.0*acos(kappa_all)/pi)*ww_res

  return
end subroutine correction_weight_function_for_viscosity

!*********************************************************************
!*** correction_weight_function_for_2D                             ***
!*********************************************************************
subroutine correction_weight_function_for_pressure       &
&          (i,drr,dist,nn,idim,x,grid_num,wall_n,grid_p, &
&           dis,ww_w,nn_wall,g_wall,wall_anc,wall_nor)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i,nn
  integer, intent(in) :: idim
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  integer, intent(in) :: grid_num
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: dis
  integer, intent(in) :: nn_wall
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  real   , intent(inout) :: ww_w
  real   , intent(inout) :: dist
  real   , intent(inout) :: drr(0:idim-1)

  ! subroutine_variable **********************************************
  integer :: j,k,l,m,n
  real    :: ww_res
  real    :: interval
  real    :: rr,rr0
  real    :: drr_all
  real    :: drr_c(0:idim-1)
  real    :: x_sub(0:idim-1,0:3)
  real    :: x0(0:idim-1)
  real    :: w_sub(0:3)
  real    :: drr0(0:idim-1)
  real, parameter   :: pi = acos(-1.0)

  ! position_setting  ************************************************
  interval = dis

  x_sub(0,0) = x(i,0)-interval
  x_sub(1,0) = x(i,1)

  x_sub(0,1) = x(i,0)+interval
  x_sub(1,1) = x(i,1)

  x_sub(0,2) = x(i,0)
  x_sub(1,2) = x(i,1)-interval

  x_sub(0,3) = x(i,0)
  x_sub(1,3) = x(i,1)+interval

  ! calculation_ww_each_point ****************************************
  do k=0,3
    w_sub(k) = 10000.0
    do n=1,g_wall(grid_p(i),0)
      l = g_wall(grid_p(i),n)
      do m=0,idim-1
        x0(m) = x_sub(m,k)
      end do

      call calc_dist                             &
      &    (x0,l,                                &
      &     drr0,rr,idim,wall_n,wall_anc,wall_nor)

      w_sub(k) = min(rr,w_sub(k))
    end do
  end do

  ! correction_drr ****************************************************
  drr0(0) = (w_sub(1)-w_sub(0))
  drr0(1) = (w_sub(3)-w_sub(2))

  drr_all = sqrt(drr0(0)**2+drr0(1)**2)

  if(drr_all .ne. 0.0)then
    do m=0,idim-1
      drr0(m) = drr0(m)/drr_all
      drr(m) = drr0(m)
    end do
  else
    do m=0,idim-1
      drr(m) = 0.0
    end do
  end if

  return
end subroutine correction_weight_function_for_pressure

!*********************************************************************
!*** correction_weight_function_for_2D                             ***
!*********************************************************************
subroutine correction_weight_function                &
&          (i,drr,dist,nn,idim,x,grid_num,wall_n,grid_p, &
&           dis,ww_w,nn_wall,g_wall,wall_anc,wall_nor)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i,nn
  integer, intent(in) :: idim
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  integer, intent(in) :: grid_num
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: dis
  integer, intent(in) :: nn_wall
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  real   , intent(inout) :: ww_w
  real   , intent(inout) :: dist
  real   , intent(inout) :: drr(0:idim-1)

  ! subroutine_variable **********************************************
  integer :: j,k,l,m,n
  real    :: kappa_all
  real    :: ww_res
  real    :: interval
  real    :: rr,rr0
  real    :: drr_all
  real    :: drr_c
  real    :: x_sub(0:idim-1,0:3)
  real    :: w_sub(0:3)
  real    :: kappa(0:idim-1)
  real    :: drr0(0:idim-1)
  real    :: x0(0:idim-1)
  real, parameter   :: pi = acos(-1.0)

  ! position_setting  ************************************************
  interval = dis

  x_sub(0,0) = x(i,0)-interval
  x_sub(1,0) = x(i,1)

  x_sub(0,1) = x(i,0)+interval
  x_sub(1,1) = x(i,1)

  x_sub(0,2) = x(i,0)
  x_sub(1,2) = x(i,1)-interval

  x_sub(0,3) = x(i,0)
  x_sub(1,3) = x(i,1)+interval

  ! calculation_ww_each_point ****************************************
  do k=0,3
    w_sub(k) = 10000.0
    do n=1,g_wall(grid_p(i),0)
      l = g_wall(grid_p(i),n)
      do m=0,idim-1
        x0(m) = x_sub(m,k)
      end do

      call calc_dist                             &
      &    (x0,l,             &
      &     drr0,rr,idim,wall_n,wall_anc,wall_nor)

      w_sub(k) = min(rr,w_sub(k))
    end do
  end do

  rr = 10000.0
  do n=1,g_wall(grid_p(i),0)
    l = g_wall(grid_p(i),n)
    do m=0,idim-1
      x0(m) = x(i,m)
    end do

    call calc_dist                              &
    &    (x0,l,                      &
    &     drr0,rr0,idim,wall_n,wall_anc,wall_nor)

    rr = min(rr,rr0)
  end do

  ! correction_ww *****************************************************
  ww_res = ww_w

  kappa(0) = (w_sub(1)-2.0*rr+w_sub(0))/interval
  kappa(1) = (w_sub(3)-2.0*rr+w_sub(2))/interval

  kappa_all = (kappa(0)+kappa(1))/(2.0*sqrt(2.0))

  kappa_all = min(kappa_all,1.0)
  kappa_all = max(kappa_all,-1.0)

  ww_w = (2.0*acos(kappa_all)/pi)*ww_res

  return
end subroutine correction_weight_function

!*********************************************************************
!*** calc_dist                                                     ***
!*********************************************************************
subroutine calc_dist                             &
&          (x0,w_n,drr,rr,idim,wall_n,wall_anc,wall_nor)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in)    :: idim
  real   , intent(in)    :: x0(0:idim-1)
  integer, intent(in)    :: w_n
  integer, intent(in)    :: wall_n
  real   , intent(in)    :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in)    :: wall_nor(0:idim-1,0:wall_n-1)
  real   , intent(inout) :: drr(0:idim-1)
  real   , intent(inout) :: rr

  ! subroutine_variable **********************************************
  integer :: l
  real    :: a,b,c
  real    :: div
  real    :: x1,x2,x3
  real    :: y1,y2,y3
  real    :: x_sub ,y_sub
  real    :: x_sub1,y_sub1
  real    :: x_sub2,y_sub2
  real    :: r_ref1,r_ref2

 ! position_setting *************************************************
  x1 = wall_anc(0,w_n,0)
  x2 = wall_anc(0,w_n,1)

  y1 = wall_anc(1,w_n,0)
  y2 = wall_anc(1,w_n,1)

  ! variable_setting *************************************************
  a = wall_nor(0,w_n)
  b = wall_nor(1,w_n)
  c = -(a*x1+b*y1)

  ! perpendicular_point **********************************************
  rr     = abs(a*x0(0)+b*x0(1)+c)/sqrt(a**2+b**2)
  x_sub1 = x0(0)+rr*a/sqrt(a**2+b**2)
  y_sub1 = x0(1)+rr*b/sqrt(a**2+b**2)
  x_sub2 = x0(0)-rr*a/sqrt(a**2+b**2)
  y_sub2 = x0(1)-rr*b/sqrt(a**2+b**2)

  r_ref1 = abs(a*x_sub1+b*y_sub1+c)/sqrt(a**2+b**2)
  r_ref2 = abs(a*x_sub2+b*y_sub2+c)/sqrt(a**2+b**2)

  if(r_ref1 .le. r_ref2)then
    x_sub = x_sub1
    y_sub = y_sub1
  else
    x_sub = x_sub2
    y_sub = y_sub2
  end if

  ! distance_calc ****************************************************
  drr(0) = x0(0)-x_sub
  drr(1) = x0(1)-y_sub

  return
end subroutine calc_dist

!*********************************************************************
!*** ww_wall_input_for_2D                                          ***
!*********************************************************************
subroutine input_ww_wall(ww_num,ww_wall,dis)
  implicit none
  ! mainroutine_variable *********************************************
  real                , intent(in)  :: dis
  real   , allocatable, intent(out) :: ww_wall(:,:)
  integer             , intent(out) :: ww_num

  ! subroutine_variable **********************************************
  integer :: i

  ! file_name_setting ************************************************
  open(13,file='./data/MPS/ww_function_2d.dat')
  read(13,*)ww_num

  ! array_allocation *************************************************
  allocate(ww_wall(0:7,0:ww_num))

  ! input_data *******************************************************
  do i=0,ww_num
    read(13,*)ww_wall(0:7,i)
  end do

  ! data_correction **************************************************
  ww_wall(0,0:ww_num) = ww_wall(0,0:ww_num)*dis

  return
end subroutine input_ww_wall

!!*********************************************************************
!!*** calc_dist                                                     ***
!!*********************************************************************
!subroutine calc_dist                             &
!&          (p_n,w_n,drr,rr,idim,nn,x,wall_n,wall_anc)
!  implicit none
!  ! mainroutine_variable *********************************************
!  integer, intent(in)    :: idim,nn
!  integer, intent(in)    :: p_n,w_n
!  real   , intent(in)    :: x(0:nn-1,0:idim-1)
!  integer, intent(in)    :: wall_n
!  real   , intent(in)    :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
!  real   , intent(inout) :: drr(0:idim-1)
!  real   , intent(inout) :: rr
!
!  ! subroutine_variable **********************************************
!  integer :: l
!  real    :: a,b,c
!  real    :: div
!  real   , allocatable :: x_sub(:)
!
!  ! allocation *******************************************************
!  allocate(x_sub(0:idim-1))
!
!  ! variable_setting *************************************************
!  a = -(wall_anc(1,w_n,1)-wall_anc(1,w_n,0))
!  b = (wall_anc(0,w_n,1)-wall_anc(0,w_n,0))
!  c = -(wall_anc(0,w_n,1)-wall_anc(0,w_n,0))*wall_anc(1,w_n,0) &
!  &   +(wall_anc(1,w_n,1)-wall_anc(1,w_n,0))*wall_anc(0,w_n,0)
!
!  ! perpendicular_point **********************************************
!  div = (a*x(p_n,0)+b*x(p_n,1)+c)/(a**2+b**2)
!  x_sub(0) = x(p_n,0)-div*a
!  x_sub(1) = x(p_n,1)-div*b
!
!  ! distance_calc ****************************************************
!  do l=0,idim-1
!    drr(l) = x(p_n,l)-x_sub(l)
!  end do
!
!  rr = dist_2D(x_sub(0),x_sub(1),x(p_n,0),x(p_n,1))
!
!  ! dealocation ******************************************************
!  deallocate(x_sub)
!
!  return
!end subroutine calc_dist

!*********************************************************************
!*** wall_collision                                                ***
!*********************************************************************
subroutine wall_collision                         &
&          (i,idim,ifluid,nn,dis,dist_min,x,v,type, &
&           nump,wall_num,wall_n,wall_anc,itypep, &
&           nn_wall,g_wall,wall_nor,grid_p,grid_num)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i,idim
  integer, intent(in) :: ifluid
  integer, intent(in) :: nn
  real   , intent(in) :: dis
  real   , intent(in) :: dist_min
  integer, intent(in) :: nump
  integer, intent(in) :: wall_num
  integer, intent(in) :: wall_n
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  integer, intent(in) :: itypep(0:nn-1)
  integer, intent(in) :: type(0:ifluid)
  integer, intent(in) :: grid_p(0:nn-1)
  integer, intent(in) :: nn_wall
  integer, intent(in) :: grid_num
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  real   , intent(inout) :: x(0:nn-1,0:idim-1)
  real   , intent(inout) :: v(0:nn-1,0:idim-1)

  ! subroutine_variable ***********************************************
  integer :: j,l,k
  real    :: distance
  real    :: vel_a
  real    :: vel_n
  real    :: cos_t
  real    :: theta
  real    :: n_a,h_a
  real    :: x1,x2
  real    :: y1,y2
  real    :: nd
  real    :: v_ref0,v_ref1
  real    :: nx(0:idim-1)
  real    :: hx(0:idim-1)
  real    :: v_sub(0:idim-1)
  real, parameter   :: zero = 1.0e-10
  real, parameter   :: pi   = acos(-1.0)
  real, parameter   :: col_rat = 0.0

  ! particle-wall_collision *******************************************

    if((type(3) .le. itypep(i)) .and. &
    &  (itypep(i) .le. type(5)))then
      do k=1,g_wall(grid_p(i),0)
        j = g_wall(grid_p(i),k)
        ! wall_normal_vector
        x1 = wall_anc(0,j,0)
        x2 = wall_anc(0,j,1)

        y1 = wall_anc(1,j,0)
        y2 = wall_anc(1,j,1)

        nx(0) = -wall_nor(0,j)
        nx(1) = -wall_nor(1,j)
        nd    = -(nx(0)*x1+nx(1)*y1)

        n_a = 1.0/max(zero,sqrt(nx(0)**2+nx(1)**2))

        do l=0,idim-1
          nx(l) = nx(l)*n_a
        end do

        ! particle-wall_distance
        distance = abs(nx(0)*x(i,0)+nx(1)*x(i,1)+nd) &
        &          /sqrt(nx(0)**2+nx(1)**2)

        ! collision_judgment
        if(distance .ge. dist_min*dis)then
          ! wall-particle_distance_is_greater_than_partcle_radius
          exit
        else
          ! wall-particle_distance_is_less_than_particle_radius
          vel_a = sqrt(v(i,0)**2+v(i,1)**2)
          vel_n = v(i,0)*nx(0)+v(i,1)*nx(1)

          cos_t = -vel_n/max(zero,vel_a)
          cos_t = min(1.0,max(0.0,cos_t))
          theta = 0.5*pi-acos(cos_t)

          ! particle_velocity_vector_check
          if(vel_n .gt. 0.0)then
            exit
          end if

          ! collision_vector_component
          v_sub(0) = vel_a*cos(theta)
          v_sub(1) = vel_a*sin(theta)

          ! collision_vector_component(wall_coodinate)
          hx(0) = v(i,0)-vel_n*nx(0)
          hx(1) = v(i,1)-vel_n*nx(1)

          h_a = 1.0/max(zero,sqrt(hx(0)**2+hx(1)**2))

          do l=0,idim-1
            hx(l) = hx(l)*h_a
          end do

          v_ref0 = v(i,0)
          v_ref1 = v(i,1)

          ! rebound_velocity_component
!          v(i,0) = v_sub(0)*hx(0)+v_sub(1)*nx(0)
!          v(i,1) = v_sub(0)*hx(1)+v_sub(1)*nx(1)
          v(i,0) = v_sub(0)*hx(0)+v_sub(1)*nx(0)*col_rat
          v(i,1) = v_sub(0)*hx(1)+v_sub(1)*nx(1)*col_rat
!          v(i,0) = v_sub(0)*hx(0)
!          v(i,1) = v_sub(0)*hx(1)

          ! particle_is_inside_the_wall
          if(distance .lt. dist_min*dis)then
            do l=0,idim-1
              x(i,l) = x(i,l)+(dist_min*dis-distance)*nx(l)
            end do
          end if
        end if
      end do
    end if

  return
end subroutine wall_collision

!*********************************************************************
!*** dam_positioning                                               ***
!*********************************************************************
subroutine dam_positioning               &
&          (poix0,poiy0,poix1,poiy1,      &
&           inter_l,p_type,nump,         &
&           nn,idim,ifluid,x,type,itypep, &
&           dis,grid,grid_num,wall_num,   &
&           nn_grid,g_neigh,grid_n)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: nn,idim,ifluid
  integer, intent(in) :: p_type
  real   , intent(in) :: dis
  real   , intent(in) :: poix0,poix1
  real   , intent(in) :: poiy0,poiy1
  real   , intent(in) :: inter_l
  integer, intent(in) :: type(0:ifluid)
  real   , intent(inout) :: x(0:nn-1,0:idim-1)
  integer, intent(inout) :: itypep(0:nn-1)
  integer, intent(inout) :: nump
  integer, intent(in)    :: wall_num
  integer, intent(in)    :: grid_num
  real,    intent(in)    :: grid(1:grid_num,0:idim-1,0:1)
  integer, intent(in)    :: nn_grid
  integer, intent(inout) :: g_neigh(1:grid_num,0:nn_grid)
  integer, intent(in)    :: grid_n(0:idim-1)

  ! subroutine_variable **********************************************
  integer              :: i,j,n,k,l,g
  real   , parameter   :: pi = acos(0.0)*2.0
  real                 :: rr,rr_min
  integer              :: num(0:1)
  real                 :: x_pos(0:1)
  real :: theta

  ! search_point_number **********************************************
  num(0) = int(abs(poix0-poix1)/inter_l)
  num(1) = int(abs(poiy0-poiy1)/inter_l)

  ! particle_positioning *********************************************
  do i=0,num(0)
    do j=0,num(1)
      ! serch position
      x_pos(0) = min(poix0,poix1)+real(i)*inter_l
      x_pos(1) = min(poiy0,poiy1)+real(j)*inter_l

      ! particle_positioning *****************************************
      do n=wall_num,nn-1
        if(itypep(n) .eq. type(0))then
          do k=0,idim-1
            x(n,k) = x_pos(k)
          end do
          itypep(n) = p_type
          nump = max(nump,n)+1

          ! g_neigh_setting ****************************************
          call g_neigh_setting                        &
          &    (n,g,nn_grid,idim,grid_num,grid_n,g_neigh)
          exit
        end if
      end do
    end do
  end do

  return
end subroutine dam_positioning

!*********************************************************************
!*** line_positioning                                              ***
!*********************************************************************
subroutine line_positioning               &
&          (poix0,poiy0,poix1,poiy1,      &
&           inter_l,p_type,nump,         &
&           nn,idim,ifluid,x,type,itypep, &
&           dis,grid,grid_num,wall_num,   &
&           nn_grid,g_neigh,grid_n)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: nn,idim,ifluid
  integer, intent(in) :: p_type
  real   , intent(in) :: dis
  real   , intent(in) :: poix0,poix1
  real   , intent(in) :: poiy0,poiy1
  real   , intent(in) :: inter_l
  integer, intent(in) :: type(0:ifluid)
  real   , intent(inout) :: x(0:nn-1,0:idim-1)
  integer, intent(inout) :: itypep(0:nn-1)
  integer, intent(inout) :: nump
  integer, intent(in)    :: wall_num
  integer, intent(in)    :: grid_num
  real,    intent(in)    :: grid(1:grid_num,0:idim-1,0:1)
  integer, intent(in)    :: nn_grid
  integer, intent(inout) :: g_neigh(1:grid_num,0:nn_grid)
  integer, intent(in)    :: grid_n(0:idim-1)

  ! subroutine_variable **********************************************
  integer              :: i,j,n,k,l,g
  real   , parameter   :: pi = acos(0.0)*2.0
  real                 :: rr,rr_min
  integer              :: num(0:1)
  real                 :: x_pos(0:1)
  real :: theta

  ! search_point_number **********************************************
  num(0) = int(abs(poix0-poix1)/inter_l)
  num(1) = int(abs(poiy0-poiy1)/inter_l)

  ! particle_positioning *********************************************
  do i=0,num(0)
    do j=0,num(1)
      ! serch position
      x_pos(0) = min(poix0,poix1)+real(i)*inter_l
      x_pos(1) = min(poiy0,poiy1)+real(j)*inter_l

      do g=1,grid_num
        if(x_pos(0) .ge. grid(g,0,0))then
          if(x_pos(0) .le. grid(g,0,1))then
            if(x_pos(1) .ge. grid(g,1,0))then
              if(x_pos(1) .le. grid(g,1,1))then
                exit
              end if
            end if
          end if
        end if
      end do

      ! particle distance calculation ********************************
      rr_min = 1.0e10

      do l=1,g_neigh(g,0)
        n = g_neigh(g,l)
        if(itypep(n) .ne. type(0))then
          rr     = dist_2D(x(n,0),x(n,1),x_pos(0),x_pos(1))
          rr_min = min(rr,rr_min)
        end if
      end do

      ! particle_positioning *****************************************
      if(rr_min > dis)then
        do n=wall_num,nn-1
          if(itypep(n) .eq. type(0))then
            do k=0,idim-1
              x(n,k) = x_pos(k)
            end do
            itypep(n) = p_type
            nump = max(nump,n)+1

            ! g_neigh_setting ****************************************
            call g_neigh_setting                        &
            &    (n,g,nn_grid,idim,grid_num,grid_n,g_neigh)
            exit
          end if
        end do
      end if
    end do
  end do

  return
end subroutine line_positioning

!*********************************************************************
!*** circle_positioning                                            ***
!*********************************************************************
subroutine circle_positioning_solid         &
&          (cent0,cent1,radius,v0,v1,       &
&           inter_l,inter_t,p_type,nump,    &
&           nn,idim,ifluid,x,v,type,itypep, &
&           dis,grid,grid_num,wall_num,     &
&           nn_grid,g_neigh,grid_n,         &
&           f_type,phase_min,phase_max)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: nn,idim,ifluid
  integer, intent(in) :: p_type
  real   , intent(in) :: dis
  real   , intent(in) :: cent0
  real   , intent(in) :: cent1
  real   , intent(in) :: v0
  real   , intent(in) :: v1
  real   , intent(in) :: radius
  real   , intent(in) :: inter_l
  real   , intent(in) :: inter_t
  integer, intent(in) :: type(0:ifluid)
  real   , intent(inout) :: x(0:nn-1,0:idim-1)
  real   , intent(inout) :: v(0:nn-1,0:idim-1)
  integer, intent(inout) :: itypep(0:nn-1)
  integer, intent(inout) :: nump
  integer, intent(in)    :: wall_num
  integer, intent(in)    :: grid_num
  real,    intent(in)    :: grid(1:grid_num,0:idim-1,0:1)
  integer, intent(in)    :: nn_grid
  integer, intent(inout) :: g_neigh(1:grid_num,0:nn_grid)
  integer, intent(inout) :: f_type(1:3)
  integer, intent(inout) :: phase_min
  integer, intent(inout) :: phase_max
  integer, intent(in)    :: grid_n(0:idim-1)

  ! subroutine_variable **********************************************
  integer              :: i,j,n,k,l,g
  real   , parameter   :: pi = acos(0.0)*2.0
  real                 :: rr,rr_min
  integer              :: num(0:1)
  real                 :: x_pos(0:1)
  real :: theta

  ! particle_removing ************************************************
  do i=0,nump-1
    if(sqrt((x(i,0)-cent0)**2+(x(i,1)-cent1)**2) .lt. 0.5*radius+1.0*dis)then
      itypep(i) = type(0)
      do k=0,idim-1
        x(i,k) = 0.0
      end do
    end if
  end do

  ! search_point_number **********************************************
  num(0) = int(radius*0.5/inter_l)
  num(1) = int(360.0/inter_t)

  ! particle_positioning *********************************************
  do i=num(0)-1,0,-1
    do j=0,num(1)-1
      ! serch position
      theta    = (real(j)*inter_t+90.0)*pi/180.0
      x_pos(0) = cent0+real(i)*inter_l*cos(theta)
      x_pos(1) = cent1+real(i)*inter_l*sin(theta)

      do g=1,grid_num
        if(x_pos(0) .ge. grid(g,0,0))then
          if(x_pos(0) .le. grid(g,0,1))then
            if(x_pos(1) .ge. grid(g,1,0))then
              if(x_pos(1) .le. grid(g,1,1))then
                exit
              end if
            end if
          end if
        end if
      end do

      ! particle distance calculation ********************************
      rr_min = 1.0e10

      do l=1,g_neigh(g,0)
        n      = g_neigh(g,l)
        if(itypep(n) .ne. type(0))then
          rr     = dist_2D(x(n,0),x(n,1),x_pos(0),x_pos(1))
          rr_min = min(rr,rr_min)
        end if
      end do

      ! particle_positioning *****************************************
      if(rr_min > dis)then
        do n=wall_num,nn-1
          if(itypep(n) .eq. type(0))then
            do k=0,idim-1
              x(n,k) = x_pos(k)
            end do

            v(n,0) = v0
            v(n,1) = v1
            itypep(n) = p_type
            nump = max(nump,n)+1

            if(p_type .eq. type(3))then
              f_type(3) = type(3)
            else if(p_type .eq. type(4))then
              f_type(1) = type(4)
            else if(p_type .eq. type(5))then
              f_type(2) = type(5)
            end if

            phase_min = min(phase_min,p_type)
            phase_max = max(phase_max,p_type)

            ! g_neigh_setting ****************************************
            call g_neigh_setting                        &
            &    (n,g,nn_grid,idim,grid_num,grid_n,g_neigh)
            exit
          end if
        end do
      end if
    end do
  end do

  return
end subroutine circle_positioning_solid

!*********************************************************************
!*** NaN_checker                                                   ***
!*********************************************************************
integer function nan_2D(a)
  implicit none
  ! mainroutine_variable *********************************************
  real, intent(in) :: a

  ! NaN_checker ******************************************************
  if(a*0.0 /= 0.0)then
    nan_2D = 1
  else
    nan_2D = 0
  end if

  return
end function nan_2D

!*********************************************************************
!*** calculation_surface_tension_term                              ***
!*********************************************************************
subroutine surface_tension_term                &
&          (nn,idim,ifluid,neighbor,wall_num,  &
&           dis,nump,itypep,neigh,beta,type,   &
&           n0,an,ker_r,x,sigma,delta,den,surf)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: nn,idim,ifluid,neighbor
  integer, intent(in) :: wall_num
  real   , intent(in) :: dis
  integer, intent(in) :: nump
  integer, intent(in) :: itypep(0:nn-1)
  integer, intent(in) :: type(0:ifluid)
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  real   , intent(in) :: beta
  real   , intent(in) :: sigma,delta
  real   , intent(in) :: n0(0:nn-1,0:3)
  real   , intent(in) :: an(0:nn-1)
  real   , intent(in) :: ker_r(0:3)
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  real   , intent(in) :: den(0:ifluid)
  real   , intent(inout) :: surf(0:nn-1,0:idim-1)

  ! subroutine_variable **********************************************
  real    :: pi
  integer :: i,j,l,m
  real    :: pdst
  real    :: theta
  real    :: dum
  integer :: nsmooth
  real    :: dist(0:idim-1)
  real    :: diff(0:idim-1)
  real    :: xd(0:idim-1,1:6)
  real   , allocatable :: akk(:)
  real   , allocatable :: pdst1(:)
  real   , allocatable :: anst(:,:)
  real   , allocatable :: dn(:,:)

  ! array_allocation *************************************************
  allocate(akk(0:nn-1))
  allocate(pdst1(0:nn-1))
  allocate(anst(0:nn-1,0:idim-1))
  allocate(dn(0:nn-1,0:idim-1))

  ! initial_setting **************************************************
  pi  = acos(-1.0)
  nsmooth = 0

!  ! judging_surface_particle *****************************************
!  do i=wall_num,nump-1
!    pdst1(i) = 0.0
!    if((itypep(i) .eq. type(4)).and.  &
!    &  (an(i) .lt. beta*n0(i,1)))then
!      do l=1,neigh(i,0)
!        j = neigh(i,l)
!        if((itypep(j) .eq. type(2)).or.  &
!        &  (itypep(j) .eq. type(4)))then
!          pdst1(i) = pdst1(i)+1.0
!        end if
!      end do
!    end if
!  end do

  ! judging_surface_tension(second) **********************************
  do i=wall_num,nump-1
    akk(i)= 0.0
    if((itypep(i) .eq. type(4)).and. &
    &  (an(i) .lt. beta*n0(i,1)))then
      pdst = 0.0
      do l=1,neigh(i,0)
        j = neigh(i,l)
        if((itypep(j) .eq. type(2)).or.  &
        &  (itypep(j) .eq. type(4)))then
!          if(pdst1(j) .ge. pdst1(i)) then
            pdst = pdst+1.0
!          end if
        end if
      end do

      ! curvature
      theta  = 0.5*pdst/n0(i,3)*pi
      akk(i) = 2.0*cos(theta)/ker_r(3)
    end if
  end do

  ! calculation_surface_tension_direction ****************************
  do i=wall_num,nump-1
    ! array_initialize
    anst(i,0:idim-1) = 0.0

    if(itypep(i) .eq. type(4)) then
      if(an(i) .lt. beta*n0(i,1)) then
        ! position_setting
        xd(0,1) = x(i,0)+dis
        xd(0,2) = x(i,0)
        xd(0,3) = x(i,0)-dis
        xd(0,4) = x(i,0)
        xd(1,1) = x(i,1)
        xd(1,2) = x(i,1)+dis
        xd(1,3) = x(i,1)
        xd(1,4) = x(i,1)-dis

        do m=1,4
          pdst1(m) = 0.0
          do l=1,neigh(i,0)
            j = neigh(i,l)
            if(itypep(j).eq.type(4)) then
              dum = sqrt((x(j,0)-xd(0,m))**2 &
              &         +(x(j,1)-xd(1,m))**2)
              if(dum .le. ker_r(3)) then
                pdst1(m) = pdst1(m)+1.0
              end if
            end if
          end do
        end do

        ! normal_vector_component
        do m=0,idim-1
          diff(m) = 0.5*(pdst1(m+1)-pdst1(m+3))/dis
        end do

        ! magnitude_of_vector
        dum = sqrt(diff(0)**2 &
        &         +diff(1)**2)

        ! unit_normal_vector
        if(dum .ne. 0.0) then
          do m=0,idim-1
            anst(i,m) = diff(m)/dum
          end do
        end if
      end if
    end if
  end do

  ! vector_smoothing *************************************************
  do while(nsmooth < 10)
    nsmooth = nsmooth+1
    do i=wall_num,nump-1
      do m=0,idim-1
        dn(i,m) = 0.0
      end do
      if(itypep(i) .eq. type(4)) then
        if(an(i) .le. beta*n0(i,1)) then
          do m=0,idim-1
            dist(m) = 0.0
          end do

          do l=1,neigh(i,0)
            j = neigh(i,l)
            if(itypep(j) .eq. type(4)) then
              do m=0,idim-1
                dist(m) = dist(m)+anst(j,m)
              end do
            end if
          end do

          dum = sqrt(dist(0)**2+dist(1)**2)

          if(dum .ne. 0.0) then
            do m=0,idim-1
              dn(i,m) = dist(m)/dum
            end do
          end if
        end if
      end if
    end do

    do i= wall_num,nump-1
      if(itypep(i) .eq. type(4)) then
        if(an(i) .lt. beta*n0(i,1)) then
          do m=0,idim-1
            anst(i,m) = dn(i,m)
          end do
        end if
      end if
    end do
  end do

  ! calculation_surface_tension_term *********************************
  do i= wall_num,nump-1
    if((itypep(i) .eq. type(4)).and. &
    &  (an(i) .lt. beta*n0(i,1))) then
      do m=0,idim-1
        surf(i,m) = sigma*akk(i)*anst(i,m)*delta/den(itypep(i))
      end do
    else
      do m=0,idim-1
        surf(i,m) = 0.0
      end do
    end if
  end do

  ! array_deallocation ***********************************************
  deallocate(akk)
  deallocate(pdst1)
  deallocate(anst)
  deallocate(dn)

  return
end subroutine surface_tension_term

!*********************************************************************
!*** calculation_distance_for_2D                                   ***
!*********************************************************************
real function dist_2D(x1,y1,x2,y2)
  implicit none
  ! mainroutine_variable *********************************************
  real :: x1,x2
  real :: y1,y2

  ! calculation_distance *********************************************
  dist_2D = (x2-x1)**2+(y2-y1)**2
  dist_2D = sqrt(dist_2D)
end function dist_2D

!*********************************************************************
!*** arrange_particle                                              ***
!*********************************************************************
subroutine arrange_particle            &
&          (type_p,                    &
&           cord,cord_max,dis_rat,     &
&           nn,idim,ifluid,dis,        &
&           nump,type,itypep,x,v,p,an)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in)    :: nn,idim,ifluid
  integer, intent(in)    :: cord_max
  integer, intent(in)    :: type(0:ifluid)
  integer, intent(in)    :: type_p
  real   , intent(in)    :: dis
  real   , intent(in)    :: dis_rat
  real   , intent(in)    :: cord(0:idim-1,0:1,0:cord_max)
  integer, intent(inout) :: nump
  real   , intent(inout) :: an(0:nn-1)
  real   , intent(inout) :: p(0:nn-1)
  integer, intent(inout) :: itypep(0:nn-1)
  real   , intent(inout) :: x(0:nn-1,0:idim-1)
  real   , intent(inout) :: v(0:nn-1,0:idim-1)

  ! subroutine_variable ***********************************************
  integer              :: i,j,n,k
  real                 :: rr,rr_min
  integer              :: count(0:idim-1)
  real                 :: point(0:idim-1)
  real                 :: dist(0:idim-1)

  ! count_searching_point **********************************************
  do i=0,idim-1
    count(i) = int(abs(cord(i,1,1)-cord(i,0,1))/(dis_rat*dis))
  end do

  ! arrange_particle ***************************************************
  do j=0,count(1)
    do i=0,count(0)
      ! initial_setting
      rr_min = 1.0e10
      ! position_setting
      point(0) = cord(0,0,1)+real(i)*dis*dis_rat
      point(1) = cord(1,0,1)+real(j)*dis*dis_rat

      ! interparticle_distance
      do n=0,nump-1
        if((itypep(n) .ne. type(0)))then
          rr = dist_2D(x(n,0),x(n,1),point(0),point(1))
          if(rr .lt. rr_min)rr_min = rr
        end if
      end do

      ! judging_if_particle_can_be_inlet
      if(rr_min .ge. dis*0.99)then
        do n=0,nn-1
          if((itypep(n) .eq. type(0)))then
            if((point(1)  .ge. cord(1,0,1)).and. &
            &  (point(1)  .le. cord(1,1,1)).and. &
            &  (point(0)  .ge. cord(0,0,1)).and. &
            &  (point(0)  .le. cord(0,1,1)))then
              itypep(n) = type_p
              do k=0,idim-1
                x(n,k) = point(k)
                v(n,k) = 0.0
              end do
              an(n) = 0.0
              p(n)  = 0.0
            end if

            if((n .ge. nump).and.(itypep(n).ne.type(0)))nump = nump+1
            if(n .eq. nn-1)then
              write(700,*)'initial_particle : n = nn'
              write(700,*)n,nn
              stop
            end if
            exit
          end if
        end do
      end if
    end do
  end do

  return
end subroutine arrange_particle

!*********************************************************************
!*** grid_particle_setting                                         ***
!*********************************************************************
subroutine g_neigh_setting                        &
&          (i,j,nn_grid,idim,grid_num,grid_n,g_neigh)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in)    :: i,j
  integer, intent(in)    :: nn_grid
  integer, intent(in)    :: idim
  integer, intent(in)    :: grid_num
  integer, intent(in)    :: grid_n(0:idim-1)
  integer, intent(inout) :: g_neigh(1:grid_num,0:nn_grid)

  ! subroutine_variable **********************************************
  integer k

  ! neighboring_grids ************************************************
  ! left
  if(mod(j,grid_n(0)) .ne. 1)then
    k = j-1
    g_neigh(k,0) = g_neigh(k,0)+1
    g_neigh(k,g_neigh(k,0)) = i

    if(j .le. grid_num-grid_n(0))then
      k = j+grid_n(0)-1
      g_neigh(k,0) = g_neigh(k,0)+1
      g_neigh(k,g_neigh(k,0)) = i
    end if

    if(j .gt. grid_n(0))then
      k = j-grid_n(0)-1
      g_neigh(k,0) = g_neigh(k,0)+1
      g_neigh(k,g_neigh(k,0)) = i
    end if
  end if

  ! mid
  if(j .gt. grid_n(0))then
    k = j-grid_n(0)
    g_neigh(k,0) = g_neigh(k,0)+1
    g_neigh(k,g_neigh(k,0)) = i
  end if

  k = j
  g_neigh(k,0) = g_neigh(k,0)+1
  g_neigh(k,g_neigh(k,0)) = i

  if(j .le. grid_num-grid_n(0))then
    k = j+grid_n(0)
    g_neigh(k,0) = g_neigh(k,0)+1
    g_neigh(k,g_neigh(k,0)) = i
  end if

  ! right
  if(mod(j,grid_n(0)) .ne. 0)then
    k = j+1
    g_neigh(k,0) = g_neigh(k,0)+1
    g_neigh(k,g_neigh(k,0)) = i

    if(j .gt. grid_n(0))then
      k = j-grid_n(0)+1
      g_neigh(k,0) = g_neigh(k,0)+1
      g_neigh(k,g_neigh(k,0)) = i
    end if

    if(j .le. grid_num-grid_n(0))then
      k = j+grid_n(0)+1
      g_neigh(k,0) = g_neigh(k,0)+1
      g_neigh(k,g_neigh(k,0)) = i
    end if
  end if

  return
end subroutine g_neigh_setting

!*********************************************************************
!*** wall_grid_settin                                              ***
!*********************************************************************
subroutine wall_grid_setting               &
&          (idim,grid_num,grid,wall_anc,   &
&           wall_nor,g_wall,ker_r,nn_wall, &
&           wall_n)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: idim
  integer, intent(in) :: nn_wall
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_num
  real,    intent(in) :: grid(1:grid_num,0:idim-1,0:1)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  real    ,intent(in) :: ker_r(0:3)
  integer, allocatable, intent(out) :: g_wall(:,:)

  ! subroutine_variable **********************************************
  integer :: i,j,k
  real    :: a,b,c
  real    :: rr
  real    :: x_sub ,y_sub ,z_sub
  real    :: x_sub1,y_sub1,z_sub1
  real    :: x_sub2,y_sub2,z_sub2
  real    :: r_ref1,r_ref2
  real    :: vec1,vec2,vec0
  real    :: d0(0:idim-1,0:2)
  real    :: posi(0:idim-1,0:3)

  ! array_allocation *************************************************
  allocate(g_wall(1:grid_num,0:nn_wall))

  ! wall_grid_setting ************************************************
  do i=1,grid_num
    ! initialize_g_wall
    g_wall(i,0:nn_wall) = 0

    ! grid_corner_position *******************************************
    posi(0,0) = grid(i,0,0)
    posi(1,0) = grid(i,1,0)

    posi(0,1) = grid(i,0,1)
    posi(1,1) = grid(i,1,0)

    posi(0,2) = grid(i,0,1)
    posi(1,2) = grid(i,1,1)

    posi(0,3) = grid(i,0,0)
    posi(1,3) = grid(i,1,1)

    ! wall_search ****************************************************
    do j=0,wall_n-1
      a = wall_nor(0,j)
      b = wall_nor(1,j)
      c = -(a*wall_anc(0,j,0)+b*wall_anc(1,j,0))

      do k=0,3
        rr = abs(a*posi(0,k)+b*posi(1,k)+c)/sqrt(a**2+b**2)

        if(rr .lt. maxval(ker_r(0:3)))then
          x_sub1 = posi(0,k)+rr*a/sqrt(a**2+b**2)
          y_sub1 = posi(1,k)+rr*b/sqrt(a**2+b**2)
          x_sub2 = posi(0,k)-rr*a/sqrt(a**2+b**2)
          y_sub2 = posi(1,k)-rr*b/sqrt(a**2+b**2)

          r_ref1 = abs(a*x_sub1+b*y_sub1+c)/sqrt(a**2+b**2)
          r_ref2 = abs(a*x_sub2+b*y_sub2+c)/sqrt(a**2+b**2)

          if(r_ref1 .le. r_ref2)then
            x_sub = x_sub1
            y_sub = y_sub1
          else
            x_sub = x_sub2
            y_sub = y_sub2
          end if

          ! vectors_between_anchor1_to_anchor2 ***********************
          d0(0,0) = wall_anc(0,j,0)-wall_anc(0,j,1)
          d0(1,0) = wall_anc(1,j,0)-wall_anc(1,j,1)

          d0(0,1) = x_sub-wall_anc(0,j,0)
          d0(1,1) = y_sub-wall_anc(1,j,0)

          d0(0,2) = x_sub-wall_anc(0,j,1)
          d0(1,2) = y_sub-wall_anc(1,j,1)

          ! magnitude_of_vector **************************************
          vec0 = sqrt(d0(0,0)**2+d0(1,0)**2)
          vec1 = sqrt(d0(0,1)**2+d0(1,1)**2)
          vec2 = sqrt(d0(0,2)**2+d0(1,2)**2)

          ! judge_if_the_foot_on_the_segment *************************
          if(vec1 .le. vec0)then
            if(vec2 .le. vec0)then
              g_wall(i,0) = g_wall(i,0)+1

              if(g_wall(i,0) .ge. nn_wall)then
                write(*,*)'Error in wall_grid_setting'
                write(*,*)'neighboring wall is greater than nn_wall'
              end if

              g_wall(i,g_wall(i,0)) = j
              exit
            end if
          end if
        end if
      end do
    end do
  end do

  return
end subroutine wall_grid_setting


!*********************************************************************
!*** devide_particle_into_grids                                    ***
!*********************************************************************
subroutine grid_division                  &
&          (nn,nn_grid,idim,ifluid,nump,  &
&           grid_num,grid_n,grid,         &
&           g_neigh,x,itypep,type,grid_p)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in)    :: nn,nn_grid
  integer, intent(in)    :: idim,ifluid
  integer, intent(in)    :: nump
  integer, intent(in)    :: grid_num
  integer, intent(in)    :: grid_n(0:idim-1)
  real,    intent(in)    :: grid(1:grid_num,0:idim-1,0:1)
  real,    intent(in)    :: x(0:nn-1,0:idim-1)
  integer, intent(inout) :: itypep(0:nn-1)
  integer, intent(in)    :: type(0:ifluid)
  integer, intent(inout) :: g_neigh(1:grid_num,0:nn_grid)
  integer, intent(inout) :: grid_p(0:nn-1)

  ! subroutine_variable **********************************************
  integer :: ip,jg
  integer :: around

  ! この変数達をallocatableにするとバグる、、、なんで？？
  integer :: shift(9)
  integer :: jj(0:nn-1)
  integer :: switch2(0:nn-1)

  logical :: mask(0:nump-1)

  ! local constant variables *****************************************

  shift(1) = -grid_n(0)-1
  shift(2) = -grid_n(0)
  shift(3) = -grid_n(0)+1
  shift(4) =           -1
  shift(5) =            0
  shift(6) =           +1
  shift(7) = +grid_n(0)-1
  shift(8) = +grid_n(0)
  shift(9) = +grid_n(0)+1

  ! array_initialize *************************************************

  switch2(:) = 1

  do jg=1,grid_num
    g_neigh(jg,0) = 0
  end do

  do ip=0,nump-1
    mask(ip) = .false.
  end do
  do ip=0,nump-1
    if((itypep(ip) .eq. type(2)).or.  &
    &  (itypep(ip) .ge. type(3)))then
      mask(ip) = .true.
    end if
  end do

  ! grid_serch *******************************************************

  ! neighboring_grid_search
  do around=1,9
!   !$omp parallel do default(shared) private(ip,jg) schedule(dynamic,64)
!   !$omp do
    do ip=0,nump-1
      jg = grid_p(ip)+shift(around)
      if((mask(ip)).and. &
      &  (jg .ge. 1).and. &
      &  (jg .le. grid_num))then
        if((x(ip,0) .ge. grid(jg,0,0)).and. &
        &  (x(ip,0) .le. grid(jg,0,1)).and. &
        &  (x(ip,1) .ge. grid(jg,1,0)).and. &
        &  (x(ip,1) .le. grid(jg,1,1)))then
          grid_p(ip) = jg
          jj(ip)     = jg
          switch2(ip) = 0
        endif
      endif
    end do
!   !$omp end do
!   !$omp end parallel do
    do ip=0,nump-1
      if(switch2(ip) .eq. 0)then
        mask(ip) = .false.
      end if
    end do
  end do

  ! if_the_particle_in_not_inside_the_neighboring_grids
! !$omp parallel do default(shared) private(ip,jg) schedule(dynamic,64)
! !$omp do
  do ip=0,nump-1
    if(mask(ip))then
      do jg=1,grid_num
        if((x(ip,0) .ge. grid(jg,0,0)).and. &
        &  (x(ip,0) .le. grid(jg,0,1)).and. &
        &  (x(ip,1) .ge. grid(jg,1,0)).and. &
        &  (x(ip,1) .le. grid(jg,1,1)))then
          grid_p(ip) = jg
          jj(ip)     = jg
          switch2(ip) = 0
          exit
        end if
      end do
    end if
  end do
! !$omp end do

! !$omp do
  do ip=0,nump-1
    if(switch2(ip) .eq. 0)then
      mask(ip) = .false.
    end if
  end do
! !$omp end do

  ! if_the_particle_is_not_inside_any_grids
! !$omp do
  do ip=0,nump-1
    if(mask(ip))then
      itypep(ip) = type(0)
      switch2(ip) = 1
    end if
  end do
! !$omp end do
! !$omp end parallel do

  do ip=0,nump-1
    if(switch2(ip) .eq. 0)then
      ! g_neigh_setting
      call g_neigh_setting                        &
      &    (ip,jj(ip),nn_grid,idim,grid_num,grid_n,g_neigh)
    end if
  end do


  ! checking_nn_grid ************************************************
  do jg=1,grid_num
   if(g_neigh(jg,0) .gt. nn_grid)then
    write(700,*)'g_neigh is greater than nn_grid',jg,g_neigh(jg,0)
    stop
   end if
  end do

  return
end subroutine grid_division

!*********************************************************************
!*** cal_ww_point_for_2D                                           ***
!*********************************************************************
subroutine cal_ww_point     &
&          (x_s0,x_s1,dist, &
&           x0,x1,y0,y1,    &
&           w11,w21,w12,w22)
  implicit none
  ! mainroutine_variable *********************************************
  real   , intent(in)    :: x_s0,x_s1
  real   , intent(in)    :: x0,x1
  real   , intent(in)    :: y0,y1
  real   , intent(in)    :: w11,w21,w12,w22
  real   , intent(inout) :: dist

  ! subroutine_variable  *********************************************
  real    :: alph,beta
  integer :: n
  real    :: x2,x3
  real    :: y2,y3

  ! interpolation_cal ************************************************
  alph = (x_s0-x0)/(x1-x0)
  beta = (x_s1-y0)/(y1-y0)

   dist = (1.0-beta)*((1.0-alph)*w11+alph*w21) &
   &     +beta*((1.0-alph)*w12+alph*w22)

  return
end subroutine cal_ww_point
end module mod_mps_dim
