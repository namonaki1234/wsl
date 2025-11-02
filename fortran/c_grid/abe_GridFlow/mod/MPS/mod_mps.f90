!*********************************************************************
!****** module_for_mps_calculation                                 ***
!****** ver. 2014.09.22                                            ***
!*********************************************************************
module mod_mps
use mod_mps_dim
!use prvwutil
  implicit none
  ! subroutine_definition ********************************************
  public :: surface_tension_pot
  public :: file_name
  public :: collision
  public :: solid_phase
  public :: solid_correction
  public :: validation
  public :: correction
  public :: prediction
  public :: viscous_term
  public :: cal_par_num
  public :: cal_n_div
  public :: cal_n_wall
  public :: cal_n
  public :: cal_lambda
  public :: pressure_term
  public :: external_force_term
  public :: cal_rhonu
  public :: set_rhonu
  public :: cal_dt
  public :: output_para_bin
!  public :: output_para_kmgt
  public :: output_term
  public :: output_para
  public :: output_an
  public :: set_neigh
  public :: allocation
  public :: input_file
  public :: output_file
  public :: cal_temperature_wall_out

  ! function_definition **********************************************
  public :: nan
contains

!*********************************************************************
!*** calculation_surface_tension(potential_model)                  ***
!*********************************************************************

subroutine surface_tension_pot        &
&          (nn,idim,ifluid,neighbor,  &
&           neigh,nump,itypep,type,x, &
&           dis,den,surf,sigma,       &
&           wall_num,ker_r,an,n0,     &
&           grid_num,ww_num,grid,     &
&           ww_wall,wall_n,grid_p,    &
&           grid_dist,wall_anc,       &
&           near_wall,wall_nor,       &
&           nn_wall,g_wall,ca)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: nn,idim,ifluid,neighbor
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  integer, intent(in) :: nump
  integer, intent(in) :: type(0:ifluid)
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  real   , intent(in) :: dis
  integer, intent(in) :: itypep(0:nn-1)
  real   , intent(in) :: den(0:ifluid)
  real   , intent(in) :: sigma
  integer, intent(in) :: wall_num
  real   , intent(in) :: an(0:nn-1)
  real   , intent(in) :: n0(0:nn-1,0:3)
  integer, intent(in) :: grid_num
  integer, intent(in) :: ww_num
  real   , intent(in) :: grid(1:grid_num,0:idim-1,0:1)
  real   , intent(in) :: ww_wall(0:7,0:ww_num)
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: grid_dist(1:grid_num,0:4*(idim-1)-1)
  integer, intent(in) :: near_wall(1:grid_num,0:4*(idim-1)-1)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  integer, intent(in) :: nn_wall
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  real   , intent(in) :: ker_r(0:3)
  real   , intent(inout) :: surf(0:nn-1,0:idim-1)
  real   , intent(in) :: ca

  ! subroutine_variable **********************************************
  integer :: i,j,l,m,k
  real    :: rr,rr2,rr0
  real    :: mass
  real    :: l_p
  real    :: s_p
  real    :: l_c
  real    :: pi
  real    :: theta_sf
  real    :: r_kernel
  real    :: drr0(0:idim-1)
  real    :: drr(0:idim-1)
  real    :: x0(0:idim-1)
  real    :: r1,r2,r3,rn
  integer :: l_ref
  real    :: ww

  ! array_initialize *************************************************
  surf(0:nn-1,0:idim-1) = 0.0

  ! initial_setting **************************************************
  l_p  = 2.0*sigma*(dis**2)
  !pi   = acos(-1.0)
  pi   = 3.14159265358979323
  r_kernel = maxval(ker_r(0:3))
  mass = 4.0/3.0*pi*((0.5*dis)**3)*den(itypep(4))

  ! contact_angle ****************************************************
  theta_sf = ca
  theta_sf = theta_sf*2.0*pi/360.0

  ! calculation_surface_tension **************************************
  do i=wall_num,nump-1
    if(itypep(i) .eq. type(4))then
      ! initialize
      s_p = 0.0
      l_c = 0.0

      do l=1,neigh(i,0)
        j = neigh(i,l)
!        if((itypep(j) .eq. type(2)) .or. &
!        &  (itypep(j) .eq. type(4)))then
        if((itypep(j) .ge. type(2)) .or. &
        &  (itypep(j) .le. type(4)))then
          ! interparticle_distance
          rr2 = 0.0
          do m=0,idim-1
           drr(m) = x(j,m)-x(i,m)
           rr2   = rr2+drr(m)**2
          end do
          rr = sqrt(rr2)

          ! potential_function(Kondo_model)
          s_p = s_p-1.0/3.0*(rr-1.5*dis+0.5*ker_r(1))*((rr-ker_r(1))**2)
        end if
      end do

      if(s_p .ne. 0.0)then
        ! potential_coefficient
        l_c = l_p/s_p
      end if

      ! surface_tension_cauclation_for_wall_part *********************
      if(g_wall(grid_p(i),0) .ne. 0)then
        ! calculation_rr
        rr = r_kernel

        do k=1,g_wall(grid_p(i),0)
          l = g_wall(grid_p(i),k)
          x0(0:idim-1) = x(i,0:idim-1)

          call calc_dist                     &
          &    (x0,l,drr0,rr0,               &
          &     idim,wall_n,wall_anc,wall_nor)

          if(rr0 .lt. rr)then
            rr = rr0
            do m=0,idim-1
              drr(m) = drr0(m)
            end do
          end if
        end do

        call calc_ww(1,rr,rn,ww_num,ww_wall,dis)

        ! an_correction
        call correction_weight_function                  &
        &    (i,drr,rr0,nn,idim,x,grid_num,wall_n,grid_p, &
        &     dis,rn,nn_wall,g_wall,wall_anc,wall_nor)

        ! calculation_r**1
        call calc_ww(4,rr,r1,ww_num,ww_wall,dis)

        ! an_correction
        call correction_weight_function                  &
        &    (i,drr,rr0,nn,idim,x,grid_num,wall_n,grid_p, &
        &     dis,r1,nn_wall,g_wall,wall_anc,wall_nor)

        ! calculation_r**1
        call calc_ww(5,rr,r2,ww_num,ww_wall,dis)

        ! an_correction
        call correction_weight_function                  &
        &    (i,drr,rr0,nn,idim,x,grid_num,wall_n,grid_p, &
        &     dis,r2,nn_wall,g_wall,wall_anc,wall_nor)

        ! calculation_r**1
        call calc_ww(6,rr,r3,ww_num,ww_wall,dis)

        ! an_correction
        call correction_weight_function                  &
        &    (i,drr,rr0,nn,idim,x,grid_num,wall_n,grid_p, &
        &     dis,r3,nn_wall,g_wall,wall_anc,wall_nor)

        ! potential_function
        s_p = s_p-1.0/3.0*(r3*(dis**3)                                &
        &                 +r2*(dis**2)*(-1.5*ker_r(1)-1.5*dis)        &
        &                 +r1*(dis**1)*(3.0*dis*ker_r(1))             &
        &                 +rn*(-1.5*dis*(ker_r(1)**2)+0.5*(ker_r(1)**3)))

        if(s_p .ne. 0.0)then
          ! potential_coefficient
          l_c = l_p/s_p
          ! calculation_surface_tension_term(wall_part)
          do m=0,idim-1
            surf(i,m) = surf(i,m)-l_c*0.5*(1.0+cos(theta_sf))     &
            &           /mass*(r2*(dis**2)-(dis+ker_r(1))*r1*dis  &
            &                  +dis*ker_r(1)*rn)*drr(m)/rr
          end do
        end if
      end if

      if(s_p .ne. 0.0)then
      ! surface_tension_calculation_for_particle_part ****************
        do l=1,neigh(i,0)
          j = neigh(i,l)
          if((itypep(j) .eq. type(2)) .or. &
          &  (itypep(j) .eq. type(4)))then
            ! interparticle_distance
            rr2 = 0.0
            do m=0,idim-1
              drr(m) = x(j,m)-x(i,m)
              rr2    = rr2+drr(m)**2
            end do
            rr = sqrt(rr2)

  !          if(itypep(j) .eq. type(4))then
  !            do m=0,idim-1
  !              surf(i,m) = surf(i,m)+l_c/mass*(rr-dis)*(rr-ker_r(1))*drr(m)/rr
  !            end do
  !          else if(itypep(j) .eq. type(2))then
  !            do m=0,idim-1
  !              surf(i,m) = surf(i,m)+l_c*0.5*(1.0+cos(theta_sf)) &
  !              &           /mass*(rr-dis)*(rr-ker_r(1))*drr(m)/rr
  !            end do
  !         end if

           if((itypep(j) .eq. type(2)) .or. &
          &  (itypep(j) .eq. type(4)))then
              do m=0,idim-1
                surf(i,m) = surf(i,m)+l_c/mass*(rr-dis)*(rr-ker_r(1))*drr(m)/rr
              end do
           end if
          end if
        end do
      end if
    end if
  end do

  return
end subroutine surface_tension_pot


subroutine surface_tension_pot2          &
&          (i,nn,idim,ifluid,neighbor,  &
&           neigh,nump,itypep,type,x,   &
&           dis,den,surf,sigma,         &
&           wall_num,ker_r,an,n0,       &
&           grid_num,ww_num,grid,       &
&           ww_wall,wall_n,grid_p,      &
&           grid_dist,wall_anc,         &
&           near_wall,wall_nor,         &
&           nn_wall,g_wall,ca,          &
&           swi_solid,swi_hc,fflat)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i
  integer, intent(in) :: nn,idim,ifluid,neighbor
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  integer, intent(in) :: nump
  integer, intent(in) :: type(0:ifluid)
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  real   , intent(in) :: dis
  integer, intent(in) :: itypep(0:nn-1)
  real   , intent(in) :: den(0:ifluid)
  real   , intent(in) :: sigma
  integer, intent(in) :: wall_num
  real   , intent(in) :: an(0:nn-1)
  real   , intent(in) :: n0(0:nn-1,0:3)
  integer, intent(in) :: grid_num
  integer, intent(in) :: ww_num
  real   , intent(in) :: grid(1:grid_num,0:idim-1,0:1)
  real   , intent(in) :: ww_wall(0:7,0:ww_num)
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: grid_dist(1:grid_num,0:4*(idim-1)-1)
  integer, intent(in) :: near_wall(1:grid_num,0:4*(idim-1)-1)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  integer, intent(in) :: nn_wall
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  real   , intent(in) :: ker_r(0:3)
  real   , intent(inout) :: surf(0:nn-1,0:idim-1)
  real   , intent(in) :: ca
  integer, intent(in) :: swi_solid(0:nn-1)
  integer, intent(inout) :: swi_hc(0:nn-1)
  real   , intent(in) :: fflat

  ! subroutine_variable **********************************************
  integer :: j,l,m,k
  real    :: rr,rr2,rr0
  real    :: mass
  real    :: l_p
  real    :: s_p
  real    :: l_c
  real    :: pi
  real    :: theta_sf
  real    :: r_kernel
  real    :: drr0(0:idim-1)
  real    :: drr(0:idim-1)
  real    :: x0(0:idim-1)
  real    :: r1,r2,r3,rn
  integer :: l_ref
  real    :: ww
  real    :: fi(0:idim-1)
  real    :: ff

  ! array_initialize *************************************************
  fi(0:idim-1) = 0.0

  ! initial_setting **************************************************
  l_p  = 2.0*sigma*(dis**2)
  !pi   = acos(-1.0)
  pi   = 3.14159265358979323
  r_kernel = maxval(ker_r(0:3))
  mass = 4.0/3.0*pi*((0.5*dis)**3)*den(itypep(4))

  ! contact_angle ****************************************************
  theta_sf = ca
  theta_sf = theta_sf*2.0*pi/360.0

  ! calculation_surface_tension **************************************
!    if(itypep(i) .eq. type(4))then
    if(itypep(i) .ne. type(0)) then
      ! initialize
      s_p = 0.0
      l_c = 0.0

      do l=1,neigh(i,0)
        j = neigh(i,l)
!        if((itypep(j) .eq. type(2)) .or. &
!        &  (itypep(j) .eq. type(4)))then
        if((itypep(j) .ge. type(2)) .or. &
        &  (itypep(j) .le. type(4)))then
          ! interparticle_distance
          rr2 = 0.0
          do m=0,idim-1
           drr(m) = x(j,m)-x(i,m)
           rr2   = rr2+drr(m)**2
          end do
          rr = sqrt(rr2)

          ! potential_function(Kondo_model)
          s_p = s_p-1.0/3.0*(rr-1.5*dis+0.5*ker_r(1))*((rr-ker_r(1))**2)
        end if
      end do

      if(s_p .ne. 0.0)then
        ! potential_coefficient
        l_c = l_p/s_p
      end if

      ! surface_tension_cauclation_for_wall_part *********************
      if(g_wall(grid_p(i),0) .ne. 0)then
        ! calculation_rr
        rr = r_kernel

        do k=1,g_wall(grid_p(i),0)
          l = g_wall(grid_p(i),k)
          x0(0:idim-1) = x(i,0:idim-1)

          call calc_dist                     &
          &    (x0,l,drr0,rr0,               &
          &     idim,wall_n,wall_anc,wall_nor)

          if(rr0 .lt. rr)then
            rr = rr0
            do m=0,idim-1
              drr(m) = drr0(m)
            end do
          end if
        end do

        call calc_ww(1,rr,rn,ww_num,ww_wall,dis)

        ! an_correction
        call correction_weight_function                  &
        &    (i,drr,rr0,nn,idim,x,grid_num,wall_n,grid_p, &
        &     dis,rn,nn_wall,g_wall,wall_anc,wall_nor)

        ! calculation_r**1
        call calc_ww(4,rr,r1,ww_num,ww_wall,dis)

        ! an_correction
        call correction_weight_function                  &
        &    (i,drr,rr0,nn,idim,x,grid_num,wall_n,grid_p, &
        &     dis,r1,nn_wall,g_wall,wall_anc,wall_nor)

        ! calculation_r**1
        call calc_ww(5,rr,r2,ww_num,ww_wall,dis)

        ! an_correction
        call correction_weight_function                  &
        &    (i,drr,rr0,nn,idim,x,grid_num,wall_n,grid_p, &
        &     dis,r2,nn_wall,g_wall,wall_anc,wall_nor)

        ! calculation_r**1
        call calc_ww(6,rr,r3,ww_num,ww_wall,dis)

        ! an_correction
        call correction_weight_function                  &
        &    (i,drr,rr0,nn,idim,x,grid_num,wall_n,grid_p, &
        &     dis,r3,nn_wall,g_wall,wall_anc,wall_nor)

        ! potential_function
        s_p = s_p-1.0/3.0*(r3*(dis**3)                                &
        &                 +r2*(dis**2)*(-1.5*ker_r(1)-1.5*dis)        &
        &                 +r1*(dis**1)*(3.0*dis*ker_r(1))             &
        &                 +rn*(-1.5*dis*(ker_r(1)**2)+0.5*(ker_r(1)**3)))

        if(s_p .ne. 0.0)then
          ! potential_coefficient
          l_c = l_p/s_p
          ! calculation_surface_tension_term(wall_part)
          do m=0,idim-1
            if(swi_solid(i) .eq. 0) then
              surf(i,m) = surf(i,m)-l_c*0.5*(1.0+cos(theta_sf))     &
              &           /mass*(r2*(dis**2)-(dis+ker_r(1))*r1*dis  &
              &                  +dis*ker_r(1)*rn)*drr(m)/rr
            end if
            fi(m) = fi(m) - (r2*(dis**2)-(dis+ker_r(1))*r1*dis  &
            &                  +dis*ker_r(1)*rn)*drr(m)/rr
          end do
        end if
      end if

      if(s_p .ne. 0.0)then
      ! surface_tension_calculation_for_particle_part ****************
        do l=1,neigh(i,0)
          j = neigh(i,l)
          if((itypep(j) .eq. type(2)) .or. &
          &  (itypep(j) .eq. type(4)))then
            ! interparticle_distance
            rr2 = 0.0
            do m=0,idim-1
              drr(m) = x(j,m)-x(i,m)
              rr2    = rr2+drr(m)**2
            end do
            rr = sqrt(rr2)

  !          if(itypep(j) .eq. type(4))then
  !            do m=0,idim-1
  !              surf(i,m) = surf(i,m)+l_c/mass*(rr-dis)*(rr-ker_r(1))*drr(m)/rr
  !            end do
  !          else if(itypep(j) .eq. type(2))then
  !            do m=0,idim-1
  !              surf(i,m) = surf(i,m)+l_c*0.5*(1.0+cos(theta_sf)) &
  !              &           /mass*(rr-dis)*(rr-ker_r(1))*drr(m)/rr
  !            end do
  !         end if

           if((itypep(j) .eq. type(2)) .or. &
          &  (itypep(j) .eq. type(4)))then
              do m=0,idim-1
                if(swi_solid(i) .eq. 0) then
                  surf(i,m) = surf(i,m)+l_c/mass*(rr-dis)*(rr-ker_r(1))*drr(m)/rr
                end if
                fi(m) = fi(m) + (rr-dis)*(rr-ker_r(1))*drr(m)/rr
              end do
           end if
          end if
        end do
      end if
      ff = 0.0
      do m = 0,idim-1
       ff = ff + fi(m)**2.0
      end do
      ff = sqrt(ff)
      if(ff/fflat .ge. 0.6) swi_hc(i) = 1 !0.32
    end if

  return
end subroutine surface_tension_pot2

!*********************************************************************
!*** Calculation F_i_flat for surface tension Sugii model          ***
!*********************************************************************
subroutine calculation_fi_flat(idim,dis,fflat)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in)  :: idim
  real   , intent(in)  :: dis
  real   , intent(out) :: fflat

  !subroutine_variable ***********************************************
  integer :: i,j,k,l,m
  real    :: rr2,rr
  real    :: fflat2
  real    :: dfflat
  integer :: ker_num
  real    :: r_kernel
  integer, parameter :: i_max = 11
  integer, parameter :: k_max = 11
  integer, parameter :: j_max = 4
  real   , parameter :: zero = 1.0e-10

  ! allocatable_variable *********************************************
  real, allocatable :: x_main(:,:)
  real, allocatable :: x_wall(:,:)
  real, allocatable :: drr(:)
  real, allocatable :: dfflatt(:)

  ! initial_setting **************************************************
  ker_num = 4
  fflat = 0.0
  r_kernel = 3.1 * dis

  ! allocation *******************************************************
  allocate(x_main(0:idim-1,1:ker_num))
  allocate(x_wall(0:idim-1,1:i_max*j_max*k_max))
  allocate(drr(0:idim-1))
  allocate(dfflatt(0:idim-1))

  ! main_positioning *************************************************
  if(idim .eq. 2)then
    do i=1,ker_num
      x_main(0,i) = 0.0
      x_main(1,i) = dis*0.5*(-1.0)
    end do
  else if(idim .eq. 3)then
    do i=1,ker_num
      x_main(0,i) = 0.0
      x_main(1,i) = dis*0.5*(-1.0)
      x_main(2,i) = 0.0
    end do
  end if

  ! wall_positioning *************************************************
  l = 0
  if(idim .eq. 3)then
    do i=1,i_max
      do j=1,j_max
        do k=1,k_max
          l = l+1
          x_wall(0,l) = dis*real(i-int((i_max+1)/2))
          x_wall(1,l) = dis*(real(j)-0.5)*(-1.0)
          x_wall(2,l) = dis*real(k-int((k_max+1)/2))
        end do
      end do
    end do
  else if(idim .eq. 2)then
    do i=1,i_max
      do j=1,j_max
        l = l+1
        x_wall(0,l) = dis*real(i-int((i_max+1)/2))
        x_wall(1,l) = dis*(real(j)-0.5)*(-1.0)
      end do
    end do
  end if

  dfflatt=0.0
  i = 1
  do j = 1,l
    ! interparticle_distance

    rr2 = 0.0
    do m=0,idim-1
     drr(m) = x_wall(m,j)-x_main(m,i)
     rr2   = rr2+drr(m)**2
    end do
    rr = sqrt(rr2)

    if((rr .le. r_kernel) .and. (rr .ge. zero))then
      fflat2=0.0
      ! potential_function(Kondo_model)
      do m=0,idim-1
!        dfflat = (rr-r_kernel)*(rr-dis)*drr(m)/rr
        dfflatt(m) =dfflatt(m) + (rr-r_kernel)*(rr-dis)*drr(m)/rr
!        fflat2 = fflat2 + (dfflat**2)
      end do
!      fflat2= sqrt(fflat2)
!      fflat = fflat+fflat2
    end if
  end do

  fflat2 = (dfflatt(0)**2)+(dfflatt(1)**2)
  fflat = sqrt(fflat2)
!  write(*,*) fflat

  ! array_deallocation ***********************************************
  deallocate(x_main)
  deallocate(x_wall)
  deallocate(drr)
  deallocate(dfflatt)

  return
  
end subroutine calculation_fi_flat

!*********************************************************************
!*** file_name_maker                                               ***
!*********************************************************************
subroutine file_name(dir,dat,datfile,fcount)
  implicit none
  ! mainroutine_variable *********************************************
  character(len = 10), intent(in)    :: dir
  character(len = 10), intent(in)    :: dat
  integer            , intent(in)    :: fcount
  character(len = 40), intent(inout) :: datfile

  ! subroutine_variable **********************************************
  integer            :: i
  character(len = 5) :: count

  ! fcount(int_to_char)***********************************************
  write(count,'(i5)')fcount

  ! replace space to zero ********************************************
  do i=1,len_trim(count)
    if(count(i:i) .eq. ' ')count(i:i) = '0'
  end do

  ! file_name ********************************************************
  write(datfile,*)'./',trim(adjustl(dir)),'/',trim(adjustl(dat)),  &
  &               count,'.vtk'

  ! erase space ******************************************************
  datfile = trim(adjustl(datfile))

  return
end subroutine file_name

!*********************************************************************
!*** nan_checker                                                   ***
!*********************************************************************
integer function nan(a)
  implicit none
  ! variable_definition **********************************************
  real, intent(in) :: a

  ! nan_check ********************************************************
  if(a*0.0 /= 0.0)then
    ! a_is_nan
    nan = 1
  else
    ! a_is_not_nan
    nan = 0
  end if

  return
end function nan

!*********************************************************************
!*** particle_collision                                            ***
!*********************************************************************
subroutine collision                         &
&          (i,nn,idim,ifluid,neighbor,nump,    &
&           itypep,type,cal_sw,              &
&           den,neigh,x,v,dist_min,coll_rat, &
&           dt,phase,dis,wall_num,x_temp,v_temp)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in)    :: i,nn,idim,ifluid,neighbor
  real   , intent(in)    :: dis
  integer, intent(in)    :: nump
  integer, intent(in)    :: itypep(0:nn-1)
  integer, intent(in)    :: type(0:ifluid)
  real   , intent(in)    :: den(0:ifluid)
  integer, intent(in)    :: neigh(0:nn-1,0:neighbor)
  real   , intent(in)    :: dt
  integer, intent(in)    :: cal_sw(0:nn-1,0:2)
  integer, intent(in)    :: phase
  real   , intent(in)    :: dist_min
  real   , intent(in)    :: coll_rat
  integer, intent(in)    :: wall_num
  real   , intent(inout) :: x(0:nn-1,0:idim-1)
  real   , intent(inout) :: v(0:nn-1,0:idim-1)
  real   , intent(in) :: x_temp(0:nn-1,0:idim-1)
  real   , intent(in) :: v_temp(0:nn-1,0:idim-1)

  ! subroutine_variable **********************************************
  integer :: j,l,m
  real    :: m1,m2,mm
  real    :: rr,rr2
  real    :: vabs,vrat
  real    :: col_vel(0:idim-1,1:2)
  real    :: distance(0:idim-1)
  real    :: vel_g(0:idim-1)
  real    :: vel_r(0:idim-1)
  real    :: vel_m(0:idim-1)

  ! particle_collision ***********************************************
    do l=1,neigh(i,0)
      j = neigh(i,l)
!      if(j .gt. i)then
        if(cal_sw(i,0) .eq. 1)then
          if(cal_sw(j,1) .eq. 1)then
            ! particle_type_check

            do m=0,idim-1
              if(cal_sw(i,1) .eq. 1)then
                col_vel(m,1) = v_temp(i,m)
              else
                col_vel(m,1) = 0.0
              end if

              if(cal_sw(j,1) .eq. 1)then
                col_vel(m,2) = v_temp(j,m)
              else
                col_vel(m,2) = 0.0
              end if
            end do

            rr2 = 0.0
            do m=0,idim-1
              distance(m) = x_temp(j,m)-x_temp(i,m)
              rr2   = rr2+distance(m)**2
            end do

            rr = sqrt(rr2)

            ! mass_setting
            if(rr .lt. (dis*dist_min))then

              m1 = den(itypep(i))
              m2 = den(itypep(i))

              if(itypep(i) .eq. type(2))m1 = den(phase)
              if(itypep(j) .eq. type(2))m2 = den(phase)

              ! total_mass
              mm = m1+m2

              ! collision_velocity
              do m=0,idim-1
                vel_g(m) = (m1*col_vel(m,1)+m2*col_vel(m,2))/mm
                vel_r(m) = m1*(col_vel(m,1)-vel_g(m))
              end do

              ! collision_coefficient
              vabs = 0.0
              do m=0,idim-1
                vabs = vabs+vel_r(m)*distance(m)/rr
              end do

              if(vabs .ge. 0.0)then
                vrat = 1.0+coll_rat
                do m=0,idim-1
                  vel_m(m) = vrat*vabs*distance(m)/rr
                end do

                if(cal_sw(i,0) .eq. 1)then
                  do m=0,idim-1
                    v(i,m) = v(i,m)-vel_m(m)/m1
                    x(i,m) = x(i,m)-dt*vel_m(m)/m1
                  end do
                end if

!                if(cal_sw(j,0) .eq. 1)then
!                  do m=0,idim-1
!                    v(j,m) = v(j,m)+vel_m(m)/m2
!                    x(j,m) = x(j,m)+dt*vel_m(m)/m2
!                   end do
!                end if
              end if
            end if
          end if
        end if
!      end if
    end do

  return
end subroutine collision

!*********************************************************************
!*** solid_surface_search                                          ***
!*********************************************************************
subroutine solid_surface                      &
&          (nn,idim,ifluid,neighbor,nump,dis, &
&           wall_num,type,itypep,x,neigh,solid_sur)
  implicit none
  ! mainroutine_variable *********************************************
  integer,intent(in) :: nn,idim,ifluid,neighbor
  integer,intent(in) :: nump
  integer,intent(in) :: wall_num
  integer,intent(in) :: type(0:ifluid)
  integer,intent(in) :: itypep(0:nn-1)
  real   ,intent(in) :: x(0:nn-1,0:idim-1)
  integer,intent(in) :: neigh(0:nn-1,0:neighbor)
  real   ,intent(in) :: dis
  integer,intent(inout) :: solid_sur(0:nn-1)

  ! subroutine_variable **********************************************
  integer :: i,k,j,m
  real    :: rr
  real    :: r_kernel
  real    :: dif_all
  real    :: sur_beta
  integer, parameter :: r_ker = 3.1
  real   , parameter :: sur_2d = 6.0
  real   , parameter :: sur_3d = 6.0
  real   , parameter :: zero = 1.0e-10
  real    :: dif(0:idim-1)

  ! array_initializing ***********************************************
  solid_sur(0:nn-1) = 0

  ! initial_setting **************************************************
  if(idim .eq. 2)then
    sur_beta = sur_2d
  else if(idim .eq. 3)then
    sur_beta = sur_3d
  else
    write(*,*)'error : subroutine solid_surface'
    stop
  end if

  ! calculation_solid_phase ******************************************
  do i=wall_num,nump-1
    if(itypep(i) .eq. type(3))then
      dif(0:idim-1) = 0.0

      do j=wall_num,nump-1
        if(itypep(i) .eq. itypep(j))then
          ! calculation_interparticle_distance
          rr = 0.0
          do k=0,idim-1
            rr = rr+(x(j,k)-x(i,k))**2
          end do

          rr = sqrt(rr)

          if(rr .le. r_ker*dis)then
            do k=0,idim-1
              dif(k) = dif(k)+(x(j,k)-x(i,k))/max(zero,rr)
            end do
          end if
        end if
      end do

      ! calculation_particle_position_unbalance
      dif_all = 0.0
      do k=0,idim-1
        dif_all = dif_all+dif(k)**2
      end do

      dif_all = sqrt(dif_all)

      if(dif_all .ge. sur_beta)then
        solid_sur(i) = 1
      end if
    end if
  end do

  return
end subroutine solid_surface

!*********************************************************************
!*** calculation_solid_phase                                       ***
!*********************************************************************
subroutine solid_phase                        &
&          (nn,idim,ifluid,neighbor,nump,dis, &
&           wall_num,type,itypep,rho,nu,x,v,  &
&           neigh,ker_r,accel,solid_sur,n0)
  implicit none
  ! mainroutine_variable *********************************************
  integer,intent(in) :: nn,idim,ifluid,neighbor
  integer,intent(in) :: nump
  integer,intent(in) :: wall_num
  integer,intent(in) :: type(0:ifluid)
  integer,intent(in) :: itypep(0:nn-1)
  real   ,intent(in) :: rho(0:nn-1)
  real   ,intent(in) :: nu(0:nn-1)
  real   ,intent(in) :: x(0:nn-1,0:idim-1)
  real   ,intent(in) :: v(0:nn-1,0:idim-1)
  integer,intent(in) :: neigh(0:nn-1,0:neighbor)
  real   ,intent(in) :: dis
  real   ,intent(in) :: ker_r(0:3)
  real   , intent(in) :: n0(0:nn-1,0:3)
  integer,intent(in) :: solid_sur(0:nn-1)
  real   , intent(inout) :: accel(0:nn-1,0:idim-1)

  ! subroutine_variable **********************************************
  integer :: i,k,j,m
  real    :: rn_p
  real    :: vel_p
  real    :: vel_f
  real    :: nu_f
  real    :: c_d
  real    :: rho_f
  real    :: rho_s
  real    :: rr_min
  real    :: vel2
  integer :: nearest_p
  real    :: rr2,rr,rrr
  real    :: ww,ww_sum
  real    :: r_kernel
  real   , parameter :: zero = 1.0e-10
  real    :: sol_vel(0:idim-1)
  real    :: drr(0:idim-1)

  ! calculation_solid_phase ******************************************
  do i=wall_num,nump-1
!    if(itypep(i) .eq. type(3))then
    if(solid_sur(i) .eq. 1)then
      ! calculation_inter_particle_distance
      rr_min  = 1.0e10
      ww_sum  = 0.0
      do m=0,idim-1
        sol_vel(m) = 0.0
      end do

      do j=1,neigh(i,0)
        k = neigh(i,j)
        if((itypep(k) .eq. type(5)).or.  &
        &  (itypep(k) .eq. type(4)))then
          rr2 = 0.0
          do m=0,idim-1
            drr(m) = x(k,m)-x(i,m)
            rr2   = rr2+drr(m)**2
          end do
          rr = sqrt(rr2)

          if(rr .lt. rr_min)then
            nearest_p = k
            rr_min = rr
          end if

          ! calculation_weight_function
          r_kernel = ker_r(0)
          ww = r_kernel/max(zero,rr)+rr/max(zero,r_kernel)-2.0

          if(rr .gt. r_kernel)ww = 0.0
          ww_sum = ww_sum+ww

          do m=0,idim-1
            sol_vel(m) = sol_vel(m)+ww*v(k,m)/n0(i,0)
          end do
        end if
      end do

      ! calculation_solid_average_velocity
      if(ww_sum .ne. 0.0)then
!        do m=0,idim-1
!          sol_vel(m) = sol_vel(m)/ww_sum
!        end do

        nu_f    = nu(nearest_p)
        rho_f   = rho(nearest_p)
        rho_s   = rho(i)

        vel2 = 0.0
        do m=0,idim-1
          vel2 = vel2+(sol_vel(m)-v(i,m))**2
        end do
        vel_f = sqrt(vel2)

        rn_p  = dis*vel_f/nu_f

        ! calculation_cd
        if(rn_p .le. 1.0e-10)then
          c_d   = 0.0
        else if(rn_p .lt. 1000.0)then
          c_d   = 24.0/rn_p*(1.0+0.15*(rn_p**0.687))
        else
          c_d   = 0.44
        end if

        ! calculation_external_force_term
        do m=0,idim-1
          accel(i,m) = accel(i,m)                   &
          &           +0.75*c_d*rho_f/rho_s/dis     &
          &            *(sol_vel(m)-v(i,m))*abs(sol_vel(m)-v(i,m))
        end do
      end if
    end if
  end do

  return
end subroutine solid_phase

!*********************************************************************
!*** solid_acceleration_correction                                 ***
!*********************************************************************
subroutine solid_correction             &
&          (nn,idim,ifluid,type,itypep, &
&           acc_sub,wall_num,nump,dt)
 implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: nn,idim,ifluid
  integer, intent(in) :: type(0:ifluid)
  integer, intent(in) :: itypep(0:nn-1)
  integer, intent(in) :: wall_num
  integer, intent(in) :: nump
  real   , intent(in) :: dt
  real   , intent(inout) :: acc_sub(0:idim-1,0:nn-1)

  ! subroutine_variable **********************************************
  integer              :: i,m
  integer              :: sol_count
  real                 :: accel_sum(0:idim-1)

  ! solid_acceleration_correction ************************************
  ! initial_setting
  sol_count           = 0
  do m=0,idim-1
    accel_sum(m) = 0.0
  end do

  ! sol_count, accel_sum calclation
  do i=wall_num,nump-1
    if(itypep(i) .eq. type(3))then
      sol_count = sol_count+1
      do m=0,idim-1
        accel_sum(m) = accel_sum(m)+acc_sub(m,i)
      end do
    end if
  end do

  ! acceleration_correction
  if(sol_count .ne. 0)then
    do i=wall_num,nump-1
     if(itypep(i) .eq. type(3))then
       do m=0,idim-1
         acc_sub(m,i) = accel_sum(m)/real(sol_count)
       end do
     end if
    end do
  end if

  return
end subroutine solid_correction

!*********************************************************************
!*** validation                                                    ***
!*********************************************************************
subroutine validation                    &
&          (nn,idim,ifluid,neighbor,     &
&           nump,x,v,p,time_sim,ite,     &
&           surf,visc,accel,type,itypep, &
&           p_out)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: nn,idim,ifluid,neighbor
  real   , intent(in) :: time_sim
  integer, intent(in) :: nump
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  real   , intent(in) :: v(0:nn-1,0:idim-1)
  real   , intent(in) :: p(0:nn-1)
  integer, intent(in) :: ite
  real   , intent(in) :: surf(0:nn-1,0:idim-1)
  real   , intent(in) :: visc(0:nn-1,0:idim-1)
  real   , intent(in) :: accel(0:nn-1,0:idim-1)
  integer, intent(in) :: type(0:ifluid)
  integer, intent(in) :: itypep(0:nn-1)
  real   , intent(in) :: p_out(0:nn-1)

  ! subroutine_variable **********************************************
  integer :: i
  integer :: nump_out

  ! initialize *******************************************************
  nump_out = 0

  ! output_file ******************************************************
  open(23,file='validation.dat')
  write(23,*)time_sim
  write(23,*)nump
  do i=0,nump-1
    if(itypep(i) .ge. type(4))then
      write(23,'(f15.5)')p_out(i)
      nump_out = nump_out+1
    end if
  end do

  write(23,*)nump_out
  close(23)

  ! validation_output ************************************************
  write(900,*)'stop in validation(mod_mps)'
  write(900,*)'output : validation file'

  ! calculation_stop *************************************************
  stop

  return
end subroutine validation

!*********************************************************************
!*** correction_term                                               ***
!*********************************************************************
subroutine correction                         &
&          (i,nn,idim,ifluid,neighbor,wall_num, &
&           nump,type,v,x,dt,cal_sw,press,itypep)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i,nn,idim,ifluid,neighbor
  integer, intent(in) :: wall_num
  integer, intent(in) :: nump
  integer, intent(in) :: type(0:ifluid)
  real   , intent(in) :: dt
  real   , intent(inout) :: press(0:nn-1,0:idim-1)
  integer, intent(in) :: cal_sw(0:nn-1,0:2)
  integer, intent(in) :: itypep(0:nn-1)
  real   , intent(inout) :: v(0:nn-1,0:idim-1)
  real   , intent(inout) :: x(0:nn-1,0:idim-1)

  ! subroutine_variable **********************************************
  integer :: l
  real    :: dv_cor

  ! correction_for_solid_phase ***************************************
  ! 中でiが回ってるから結果がいかれるので注意
!  call solid_correction             &
!  &    (nn,idim,ifluid,type,itypep, &
!  &     press,wall_num,nump,dt)

!  ! correction_term **************************************************
!  do i=wall_num,nump-1
!    if(itypep(i) .eq. type(3))then
!      do l=0,idim-1
!        dv_cor = press(i,l)
!        v(i,l) = v(i,l)+dv_cor*dt
!        x(i,l) = x(i,l)+dv_cor*(dt**2)
!      end do
!    else if(cal_sw(i,1) .eq. 1)then
!      do l=0,idim-1
!        dv_cor = press(i,l)
!        v(i,l) = v(i,l)+dv_cor*dt
!        x(i,l) = x(i,l)+dv_cor*(dt**2)
!      end do
!    end if
!  end do

  ! correction_term **************************************************

  if(cal_sw(i,1) .eq. 1)then
    do l=0,idim-1
      dv_cor = press(i,l)
      v(i,l) = v(i,l)+dv_cor*dt
      x(i,l) = x(i,l)+dv_cor*(dt**2)
    end do
  end if


  return
end subroutine correction

!*********************************************************************
!*** prediction term                                               ***
!*********************************************************************
subroutine prediction                         &
&          (i,nn,idim,ifluid,neighbor,wall_num, &
&           itypep,nump,type,v,x,dt,visc,     &
&           accel,cal_sw,surf,time_sim)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i,nn,idim,ifluid,neighbor
  integer, intent(in) :: wall_num
  integer, intent(in) :: itypep(0:nn-1)
  integer, intent(in) :: nump
  integer, intent(in) :: type(0:ifluid)
  real   , intent(in) :: dt
  real   , intent(in) :: visc(0:nn-1,0:idim-1)
  real   , intent(inout) :: accel(0:nn-1,0:idim-1)
  real   , intent(in) :: surf(0:nn-1,0:idim-1)
  real   , intent(in) :: time_sim
  integer, intent(in) :: cal_sw(0:nn-1,0:2)
  real   , intent(inout) :: v(0:nn-1,0:idim-1)
  real   , intent(inout) :: x(0:nn-1,0:idim-1)

  ! subroutine_variable **********************************************
  integer :: l
  real    :: dv_pre

  ! correction_for_solid_phase ***************************************
  ! 中でiが回ってるから結果がいかれるので注意
!  call solid_correction             &
!  &    (nn,idim,ifluid,type,itypep, &
!  &     accel,wall_num,nump,dt)

  ! prediction_term **************************************************
    if(itypep(i) .eq. type(3))then
      do l=0,idim-1
        dv_pre = accel(i,l)
        v(i,l) = v(i,l)+dv_pre*dt
        x(i,l) = x(i,l)+v(i,l)*dt
      end do
    else if(cal_sw(i,1) .eq. 1)then
      do l=0,idim-1
        dv_pre = visc(i,l)+accel(i,l)+surf(i,l)
        v(i,l) = v(i,l)+dv_pre*dt
        x(i,l) = x(i,l)+v(i,l)*dt
      end do
    else if(itypep(i) .eq. type(2))then
      do l=0,idim-1
        x(i,l) = x(i,l)+v(i,l)*dt
      end do
    end if


  return
end subroutine prediction

!*********************************************************************
!*** calculation_viscous_term                                      ***
!*********************************************************************
subroutine viscous_term                       &
&          (i,nn,idim,ifluid,neighbor,wall_num, &
&           nump,visc,itypep,type,cal_sw,     &
&           neigh,x_temp,v_temp,nu,ker_r,rlambda,n0,    &
&           grid_num,ww_num,grid,ww_wall,     &
&           wall_n,grid_p,grid_dist,dis,      &
&           near_wall,nn_wall,g_wall,wall_anc,&
&           wall_nor,solid_sur)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i,nn,idim,ifluid,neighbor
  integer, intent(in) :: wall_num
  integer, intent(in) :: nump
  real   , intent(in) :: nu(0:nn-1)
  integer, intent(in) :: itypep(0:nn-1)
  integer, intent(in) :: type(0:ifluid)
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  real   , intent(in) :: x_temp(0:nn-1,0:idim-1)
  real   , intent(in) :: v_temp(0:nn-1,0:idim-1)
  real   , intent(in) :: rlambda
  real   , intent(in) :: n0(0:nn-1,0:3)
  integer, intent(in) :: cal_sw(0:nn-1,0:2)
  real   , intent(in) :: ker_r(0:3)
  real   , intent(in) :: dis
  integer, intent(in) :: grid_num
  integer, intent(in) :: ww_num
  real   , intent(in) :: grid(1:grid_num,0:idim-1,0:1)
  real   , intent(in) :: ww_wall(0:7,0:ww_num)
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: grid_dist(1:grid_num,0:4*(idim-1)-1)
  integer, intent(in) :: near_wall(1:grid_num,0:4*(idim-1)-1)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  integer, intent(in) :: nn_wall
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  integer,intent(in) :: solid_sur(0:nn-1)
  real   , intent(inout) :: visc(0:nn-1,0:idim-1)

  ! subroutine_variable **********************************************
  integer :: j,l,m,k
  real    :: rr,rr2,rr0
  real    :: ww
  real    :: dvis0
  real    :: r_kernel
  integer :: l_ref
  real    :: ww_ref
  real    :: dvis(0:idim-1)
  real    :: drr0(0:idim-1)
  real    :: drr(0:idim-1)
  real    :: x0(0:idim-1)
  real, parameter   :: zero = 1.0e-10

  ! initital_setting *************************************************
  r_kernel = ker_r(1)

  ! array_initialize *************************************************
  do m=0,idim-1
    visc(i,m) = 0.0
  end do

  ! calculation_viscous_term *****************************************
!  do i=wall_num,nump-1
    if((cal_sw(i,1) .eq. 1) .and. &
    &  (itypep(i) .ne. type(3)))then
      dvis0 = 0.0
      do l=1,neigh(i,0)
        j = neigh(i,l)
        if((itypep(j) .eq. type(2)).or.  &
        &  (cal_sw(j,1) .eq. 1))then
          ! interparticle_distance
          rr2 = 0.0
          do m=0,idim-1
            drr(m) = x_temp(j,m)-x_temp(i,m)
            rr2   = rr2+drr(m)**2
          end do

          rr = sqrt(rr2)

          ! weight_function
          ww = r_kernel/max(zero,rr)+rr/max(zero,r_kernel)-2.0
          if(rr .gt. r_kernel)ww = 0.0

          ! viscous_term **********************************************
          if(itypep(j) .eq. type(3))then
            dvis0 = nu(i)*2.0*real(idim)*ww/rlambda/n0(i,1)
!            if(solid_sur(i) .eq. 1)dvis0 = nu(i)*2.0*real(idim)*ww/rlambda/n0(i,1)
!            if(solid_sur(i) .ne. 0)dvis0 = 0.0
          else
            dvis0 = nu(j)*2.0*real(idim)*ww/rlambda/n0(i,1)
          end if

          ! viscous_term_component ************************************
          do m=0,idim-1
            dvis(m)   = dvis0*(v_temp(j,m)-v_temp(i,m))
            visc(i,m) = visc(i,m)+dvis(m)
          end do
        end if
      end do

      ! wall_influence ************************************************************
      if(g_wall(grid_p(i),0) .ne. 0)then
        rr = r_kernel

        do k=1,g_wall(grid_p(i),0)
          l = g_wall(grid_p(i),k)
          do m=0,idim-1
            x0(m) = x_temp(i,m)
          end do

          call calc_dist                     &
          &    (x0,l,drr0,rr0,               &
          &     idim,wall_n,wall_anc,wall_nor)

          if(rr0 .lt. rr)then
            rr = rr0
            do m=0,idim-1
              drr(m) = drr0(m)
            end do
          end if
        end do

        ! calculation_particle_number_density
        call calc_ww(3,rr,ww_ref,ww_num,ww_wall,dis)

        ! an_correction
        call correction_weight_function_for_viscosity    &
        &    (i,drr,rr,nn,idim,x_temp,grid_num,wall_n,grid_p, &
        &     dis,ww_ref,nn_wall,g_wall,wall_anc,wall_nor)

        ! viscous_term
        dvis0 = nu(i)*2.0*real(idim)*ww_ref/rlambda/n0(i,1)

        ! viscous_term_component
        do m=0,idim-1
          dvis(m)   = dvis0*(0.0-v_temp(i,m))
          visc(i,m) = visc(i,m)+dvis(m)
        end do
      end if
    end if
!  end do

  return
end subroutine viscous_term

!*********************************************************************
!*** cal_standard_particle_number_density                          ***
!*********************************************************************
subroutine cal_par_num               &
&          (i,ineigh,ker_r,          &
&           nn,idim,ifluid,neighbor, &
&           an,neigh,x,itypep,grid,  &
&           grid_num,ww_wall,ww_num, &
&           grid_p,wall_n,wall_anc)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: ineigh
  integer, intent(in) :: i
  integer, intent(in) :: nn,idim,ifluid,neighbor
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  integer, intent(in) :: itypep(0:nn-1)
  real   , intent(in) :: ker_r(0:3)
  integer, intent(in) :: grid_num
  integer, intent(in) :: ww_num
  real   , intent(in) :: grid(1:grid_num,0:idim-1,0:1)
  real   , intent(in) :: ww_wall(0:7,0:ww_num)
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(inout) :: an(0:nn-1)

  ! subroutine_variable **********************************************
  integer :: l,m,j,k
  real    :: rr2,rr,ww,ww2
  real    :: r_kernel
  real    :: ww_w
  real    :: dis_x,dis_y,dis_z
  real    :: theta
  integer :: l_ref
  real   , parameter :: pi   = acos(-1.0)
  real   , parameter :: zero = 1.0e-10
  real    :: drr(0:idim-1)

  ! initial_setting **************************************************
  an(i)    = 0.0
  r_kernel = ker_r(ineigh)

  ! calculation_standard_particle_number_density *********************
  do l=1,neigh(i,0)
    j=neigh(i,l)
    k=itypep(j)
    ! interparticle_distance
    rr2 = 0.0
    do m=0,idim-1
      drr(m) = x(j,m)-x(i,m)
      rr2 = rr2+drr(m)**2
    end do

    rr = sqrt(rr2)

    ! weight_function
    if(rr .le. r_kernel)then
      ww = 1.0
    else if(rr .gt. r_kernel)then
      ww = 0.0
    end if

    ! particle_number_density
    an(i)=an(i)+ww
  end do

  ! count_up_own_particle ******************************************
  an(i) = an(i)+1.0

  return
end subroutine cal_par_num

!*********************************************************************
!*** count_particle_number_density                                 ***
!*********************************************************************
subroutine cal_n_div                        &
&           (i,ineigh,ker_r,cal_sw,         &
&            nn,idim,ifluid,neighbor,phase, &
&            an,neigh,x,itypep,type,ww_num, &
&            grid_num,grid,ww_wall,wall_n,  &
&            grid_p,grid_dist,wall_anc,dis, &
&            near_wall,nn_wall,g_wall,wall_nor)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: ineigh
  integer, intent(in) :: i
  integer, intent(in) :: nn,idim,ifluid,neighbor
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  integer, intent(in) :: itypep(0:nn-1)
  real   , intent(in) :: ker_r(0:3)
  integer, intent(in) :: cal_sw(0:nn-1,0:2)
  integer, intent(in) :: type(0:ifluid)
  integer, intent(in) :: phase
  real   , intent(in) :: dis
  integer, intent(in) :: grid_num
  integer, intent(in) :: ww_num
  real   , intent(in) :: grid(1:grid_num,0:idim-1,0:1)
  real   , intent(in) :: ww_wall(0:7,0:ww_num)
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: grid_dist(1:grid_num,0:4*(idim-1)-1)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  integer, intent(in) :: nn_wall
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  integer, intent(in) :: near_wall(1:grid_num,0:4*(idim-1)-1)
  real   , intent(inout) :: an(0:nn-1)

  ! subroutine_variable **********************************************
  integer :: l,m,j,k
  real    :: rr2,rr,rr0,ww
  real    :: r_kernel
  integer :: number
  real    :: theta
  real    :: ww_w
  integer :: l_ref
  real   , parameter :: pi   = acos(-1.0)
  real   , parameter :: zero = 1.0e-10
  real    :: drr(0:idim-1)
  real    :: drr0(0:idim-1)
  real    :: x0(0:idim-1)

  ! initial_setting ***************************************************
  an(i) = 0.0
  r_kernel = ker_r(ineigh)

  ! count_particle_number_density *************************************
  do l=1,neigh(i,0)
    j = neigh(i,l)
    k = itypep(j)
    if((itypep(j) .le. phase).or. &
    &  (itypep(j) .eq. type(2)))then
    ! interparticle_distance
      rr2 = 0.0
      do m=0,idim-1
        drr(m) = x(j,m)-x(i,m)
        rr2 = rr2+drr(m)**2
      end do

      rr = sqrt(rr2)

    ! weight_function
      if(ineigh .eq. 2)then
        ww = r_kernel/max(zero,rr)-rr/max(r_kernel,zero)
      else
        ww = r_kernel/max(zero,rr)+rr/max(zero,r_kernel)-2.0
      end if

      if(rr .gt. r_kernel)ww = 0.0

    ! particle_number_density
      an(i) = an(i)+ww
    end if
  end do

  ! particle_number_density(wall_part) **********************************
  ! search_near_wall
  if(g_wall(grid_p(i),0) .ne. 0)then
    if(ineigh .eq. 2)then
      number = 2
    else
      number = 3
    end if

    rr = r_kernel

    do k=1,g_wall(grid_p(i),0)
      l = g_wall(grid_p(i),k)
      do m=0,idim-1
       x0(m) = x(i,m)
      end do

      call calc_dist                     &
      &    (x0,l,drr0,rr0,               &
      &     idim,wall_n,wall_anc,wall_nor)

      if(rr0 .lt. rr)then
        rr = rr0
        do m=0,idim-1
          drr(m) = drr0(m)
        end do
      end if
    end do

    ! calculation_particle_number_density
    call calc_ww(number,rr,ww_w,ww_num,ww_wall,dis)

    ! an_correction
    call correction_weight_function                  &
    &    (i,drr,rr,nn,idim,x,grid_num,wall_n,grid_p, &
    &     dis,ww_w,nn_wall,g_wall,wall_anc,wall_nor)

    ! wall_particle_number_density
    an(i) = an(i)+ww_w
  end if

  return
end subroutine cal_n_div

!*********************************************************************
!*** count_particle_number_density(without_wall)                   ***
!*********************************************************************
subroutine cal_n_wall                     &
&           (i,ineigh,ker_r,              &
&            nn,idim,ifluid,neighbor,     &
&            grid_num,ww_num,grid,        &
&            ww_wall,wall_n,grid_p,       &
&            an,neigh,x,itypep,wall_anc)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: ineigh
  integer, intent(in) :: i
  integer, intent(in) :: nn,idim,ifluid,neighbor
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  integer, intent(in) :: itypep(0:nn-1)
  real   , intent(in) :: ker_r(0:3)
  integer, intent(in) :: grid_num
  integer, intent(in) :: ww_num
  real   , intent(in) :: grid(1:grid_num,0:idim-1,0:1)
  real   , intent(in) :: ww_wall(0:7,0:ww_num)
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(inout) :: an(0:nn-1)

  ! subroutine_variable **********************************************
  integer :: l,m,j,k
  real    :: rr2,rr,ww
  real    :: ww_w
  real    :: r_kernel
  integer :: number
  real    :: theta
  real   , parameter :: pi   = acos(-1.0)
  real   , parameter :: zero = 1.0e-10
  real    :: drr(0:idim-1)

  ! initial_setting ***************************************************
  an(i)=0.0
  r_kernel = ker_r(ineigh)

  ! particle_number_density *******************************************
  do l=1,neigh(i,0)
    j=neigh(i,l)
    k=itypep(j)
    ! interparticle_distance
    rr2 = 0.0
    do m=0,idim-1
      drr(m) = x(j,m)-x(i,m)
      rr2 = rr2+drr(m)**2
    end do

    rr = sqrt(rr2)

    ! waight_function
    if(ineigh .eq. 2)then
      ww = r_kernel/max(zero,rr)-rr/max(r_kernel,zero)
    else
      ww = r_kernel/max(zero,rr)+rr/max(zero,r_kernel)-2.0
    end if

    if(rr .gt. r_kernel)ww = 0.0

    ! particle_number_density
    an(i) = an(i)+ww
  end do

  return
end subroutine cal_n_wall

!*********************************************************************
!*** count_particle_number_density                                 ***
!*********************************************************************
subroutine cal_n                      &
&           (i,ineigh,ker_r,          &
&            nn,idim,ifluid,neighbor, &
&            grid_num,ww_num,         &
&            ww_wall,wall_n,grid_p,   &
&            an,neigh,x,itypep,       &
&            wall_anc,dis,            &
&            nn_wall,g_wall,wall_nor)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: ineigh
  integer, intent(in) :: i
  integer, intent(in) :: nn,idim,ifluid,neighbor
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  real   , intent(in) :: dis
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  integer, intent(in) :: itypep(0:nn-1)
  real   , intent(in) :: ker_r(0:3)
  integer, intent(in) :: grid_num
  integer, intent(in) :: ww_num
  real   , intent(in) :: ww_wall(0:7,0:ww_num)
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  integer, intent(in) :: nn_wall
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  real   , intent(inout) :: an(0:nn-1)

  ! subroutine_variable ***********************************************
  integer :: l,m,j,k
  real    :: rr2,rr,rr0
  real    :: ww,ww_w
  real    :: r_kernel
  integer :: number
  real   , parameter :: pi   = acos(-1.0)
  real   , parameter :: zero = 1.0e-10
  real    :: drr(0:idim-1)
  real    :: drr0(0:idim-1)
  real    :: x0(0:idim-1)

  ! initializing *******************************************************
  an(i) = 0.0
  r_kernel = ker_r(ineigh)

  ! count_particle_number_density **************************************
  do l=1,neigh(i,0)
    j = neigh(i,l)
    k = itypep(j)
    ! interparticle_distance
    rr2 = 0.0
    do m=0,idim-1
      drr(m) = x(j,m)-x(i,m)
      rr2 = rr2+drr(m)**2
    end do

    rr = sqrt(rr2)

    ! weight_function
    if(ineigh .eq. 2)then
      ww = r_kernel/max(zero,rr)-rr/max(r_kernel,zero)
    else
      ww = r_kernel/max(zero,rr)+rr/max(zero,r_kernel)-2.0
    end if

    if(rr .gt. r_kernel)ww = 0.0

    ! particle_number_density
    an(i) = an(i)+ww
  end do

  ! particle_number_density(wall_part) **********************************
  ! search_near_wall
  if(g_wall(grid_p(i),0) .ne. 0)then
    if(ineigh .eq. 2)then
      number = 2
    else
      number = 3
    end if

    rr = r_kernel

    do k=1,g_wall(grid_p(i),0)
      l = g_wall(grid_p(i),k)
      do m=0,idim-1
        x0(m) = x(i,m)
      end do

      call calc_dist                     &
      &    (x0,l,drr0,rr0,               &
      &     idim,wall_n,wall_anc,wall_nor)

      if(rr0 .lt. rr)then
        rr = rr0
        do m=0,idim-1
          drr(m) = drr0(m)
        end do
      end if
    end do

    ! calculation_particle_number_density
    call calc_ww(number,rr,ww_w,ww_num,ww_wall,dis)

    ! an_correction
    call correction_weight_function                  &
    &    (i,drr,rr,nn,idim,x,grid_num,wall_n,grid_p, &
    &     dis,ww_w,nn_wall,g_wall,wall_anc,wall_nor)

    ! wall_particle_number_density
    an(i) = an(i)+ww_w
  end if

  return
end subroutine cal_n

subroutine cal_n_2                    &
&           (i,ineigh,ker_r,          &
&            nn,idim,ifluid,neighbor, &
&            grid_num,ww_num,         &
&            ww_wall,wall_n,grid_p,   &
&            an,neigh,x,itypep,       &
&            wall_anc,dis,            &
&            nn_wall,g_wall,wall_nor, &
&            swi_hc,swi_solid)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: ineigh
  integer, intent(in) :: i
  integer, intent(in) :: nn,idim,ifluid,neighbor
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  real   , intent(in) :: dis
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  integer, intent(in) :: itypep(0:nn-1)
  real   , intent(in) :: ker_r(0:3)
  integer, intent(in) :: grid_num
  integer, intent(in) :: ww_num
  real   , intent(in) :: ww_wall(0:7,0:ww_num)
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  integer, intent(in) :: nn_wall
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  real   , intent(inout) :: an(0:nn-1)
  integer, intent(inout) :: swi_hc(0:nn-1)
  integer, intent(inout) :: swi_solid(0:nn-1)

  ! subroutine_variable ***********************************************
  integer :: l,m,j,k
  real    :: rr2,rr,rr0
  real    :: ww,ww_w
  real    :: r_kernel
  integer :: number
  real   , parameter :: pi   = acos(-1.0)
  real   , parameter :: zero = 1.0e-10
  real    :: drr(0:idim-1)
  real    :: drr0(0:idim-1)
  real    :: x0(0:idim-1)
  integer :: num_2,num_4

  ! initializing *******************************************************
  an(i) = 0.0
  r_kernel = ker_r(ineigh)

  num_2 = 0
  num_4 = 0

  ! count_particle_number_density **************************************
  do l=1,neigh(i,0)
    j = neigh(i,l)
    k = itypep(j)
    ! interparticle_distance
    rr2 = 0.0
    do m=0,idim-1
      drr(m) = x(j,m)-x(i,m)
      rr2 = rr2+drr(m)**2
    end do

    rr = sqrt(rr2)

    ! weight_function
    if(ineigh .eq. 2)then
      ww = r_kernel/max(zero,rr)-rr/max(r_kernel,zero)
    else
      ww = r_kernel/max(zero,rr)+rr/max(zero,r_kernel)-2.0
    end if

    if(rr .gt. r_kernel)ww = 0.0

    ! particle_number_density
    an(i) = an(i)+ww

    select case(itypep(j))
     case(2)
      num_2 = num_2 + 1
     case(4)
      num_4 = num_4 + 1
    end select

  end do

  if(itypep(i) .eq. 2 .and. num_2 .gt. num_4 .and. num_2 .gt. 5) then
   swi_solid(i) = 1
  else
   swi_solid(i) = 0
  end if

  ! particle_number_density(wall_part) **********************************
  ! search_near_wall
  if(g_wall(grid_p(i),0) .ne. 0)then
    if(ineigh .eq. 2)then
      number = 2
    else
      number = 3
    end if

    rr = r_kernel

    do k=1,g_wall(grid_p(i),0)
      l = g_wall(grid_p(i),k)
      do m=0,idim-1
        x0(m) = x(i,m)
      end do

      call calc_dist                     &
      &    (x0,l,drr0,rr0,               &
      &     idim,wall_n,wall_anc,wall_nor)

      if(rr0 .lt. rr)then
        rr = rr0
        do m=0,idim-1
          drr(m) = drr0(m)
        end do
        if(itypep(i) .eq. 2) swi_solid(i) = 1
      end if
    end do

    ! calculation_particle_number_density
    call calc_ww(number,rr,ww_w,ww_num,ww_wall,dis)

    ! an_correction
    call correction_weight_function                  &
    &    (i,drr,rr,nn,idim,x,grid_num,wall_n,grid_p, &
    &     dis,ww_w,nn_wall,g_wall,wall_anc,wall_nor)

    ! wall_particle_number_density
    an(i) = an(i)+ww_w
  end if

  return
end subroutine cal_n_2

!*********************************************************************
!*** cal_lambda                                                    ***
!*********************************************************************
subroutine cal_lambda                      &
&           (i,ineigh,ker_r,rlambda_ref,   &
&            nn,idim,ifluid,neighbor,      &
&            grid_num,ww_num,grid,ww_wall, &
&            wall_n,grid_p,neigh,x,n0,wall_anc)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: ineigh
  integer, intent(in) :: i
  integer, intent(in) :: nn,idim,ifluid,neighbor
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  real   , intent(in) :: ker_r(0:3)
  integer, intent(in) :: grid_num
  integer, intent(in) :: ww_num
  real   , intent(in) :: grid(1:grid_num,0:idim-1,0:1)
  real   , intent(in) :: ww_wall(0:7,0:ww_num)
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: n0(0:nn-1,0:3)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(inout) :: rlambda_ref

  ! subroutine_variable ***********************************************
  integer :: l,m,j,k
  real    :: rr2,rr,ww
  real    :: r_kernel
  real    :: ww_w
  integer :: number
  real    :: theta
  real   , parameter :: pi = acos(-1.0)
  real    :: drr(0:idim-1)
  integer :: wall_sw(0:1)
  real   , parameter :: zero = 1.0e-10

  ! initial_setting ****************************************************
  rlambda_ref = 0.0
  r_kernel    = ker_r(ineigh)

  ! calculation_lambda *************************************************
  do l=1,neigh(i,0)
    j=neigh(i,l)
    ! interparticle_distance
    rr2 = 0.0
    do m=0,idim-1
      drr(m) = x(j,m)-x(i,m)
      rr2 = rr2+drr(m)**2
    end do

    rr = sqrt(rr2)

    ! weight_function
    if(ineigh .eq. 2)then
      ww = r_kernel/max(zero,rr)-rr/max(zero,r_kernel)
    else
      ww = r_kernel/max(zero,rr)+rr/max(zero,r_kernel)-2.0
    end if

    if(rr .gt. r_kernel)ww = 0.0

    ! lambda
    rlambda_ref =rlambda_ref+rr2*ww/n0(i,ineigh)
  end do

  return
end subroutine cal_lambda

subroutine cal_pressure                            &
&           (i,nn,idim,ifluid,neighbor,nump,          &
&            wall_num,cal_sw,itypep,type,rho,p,     &
&            neigh,x,den,ker_r,n0,press,c_vel,phase,&
&            grid_num,ww_num,grid,ww_wall,wall_n,   &
&            grid_p,grid_dist,near_wall,wall_anc,   &
&            dis,nn_wall,g_wall,wall_nor)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i,nn,idim,ifluid,neighbor
  integer, intent(in) :: wall_num
  integer, intent(in) :: nump
  integer, intent(in) :: cal_sw(0:nn-1,0:2)
  integer, intent(in) :: itypep(0:nn-1)
  integer, intent(in) :: type(0:ifluid)
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  real   , intent(in) :: den(0:ifluid)
  real   , intent(in) :: ker_r(0:3)
  real   , intent(in) :: n0(0:nn-1,0:3)
  real   , intent(in) :: rho(0:nn-1)
  real   , intent(in) :: c_vel
  integer, intent(in) :: phase
  integer, intent(in) :: grid_num
  integer, intent(in) :: ww_num
  real   , intent(in) :: grid(1:grid_num,0:idim-1,0:1)
  real   , intent(in) :: ww_wall(0:7,0:ww_num)
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: grid_dist(1:grid_num,0:4*(idim-1)-1)
  integer, intent(in) :: near_wall(1:grid_num,0:4*(idim-1)-1)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: dis
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  integer, intent(in) :: nn_wall
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  real   , intent(inout) :: p(0:nn-1)
  real   , intent(inout) :: press(0:nn-1,0:idim-1)

  ! subroutine_variable **********************************************
  real   , parameter :: zero = 1.0e-10

  ! array_initialize *************************************************
  ! pressure_calculation *********************************************
  if(cal_sw(i,0) .eq. 1)then
    p(i) = 0.0
    press(i,0:idim-1) = 0.0
    p(i) = max((rho(i)-den(phase)),0.0)
  end if

  return
end subroutine cal_pressure
!*********************************************************************
!*** calculation_pressure_term                                     ***
!*********************************************************************
subroutine pressure_term                            &
&           (i,nn,idim,ifluid,neighbor,nump,          &
&            wall_num,cal_sw,itypep,type,rho,p,     &
&            neigh,x,den,ker_r,n0,press,c_vel,phase,&
&            grid_num,ww_num,grid,ww_wall,wall_n,   &
&            grid_p,grid_dist,near_wall,wall_anc,   &
&            dis,nn_wall,g_wall,wall_nor,swi_solid)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i,nn,idim,ifluid,neighbor
  integer, intent(in) :: wall_num
  integer, intent(in) :: nump
  integer, intent(in) :: cal_sw(0:nn-1,0:2)
  integer, intent(in) :: itypep(0:nn-1)
  integer, intent(in) :: type(0:ifluid)
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  real   , intent(in) :: den(0:ifluid)
  real   , intent(in) :: ker_r(0:3)
  real   , intent(in) :: n0(0:nn-1,0:3)
  real   , intent(in) :: rho(0:nn-1)
  real   , intent(in) :: c_vel
  integer, intent(in) :: phase
  integer, intent(in) :: grid_num
  integer, intent(in) :: ww_num
  real   , intent(in) :: grid(1:grid_num,0:idim-1,0:1)
  real   , intent(in) :: ww_wall(0:7,0:ww_num)
  integer, intent(in) :: wall_n
  integer, intent(in) :: grid_p(0:nn-1)
  real   , intent(in) :: grid_dist(1:grid_num,0:4*(idim-1)-1)
  integer, intent(in) :: near_wall(1:grid_num,0:4*(idim-1)-1)
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: dis
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  integer, intent(in) :: nn_wall
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  real   , intent(inout) :: p(0:nn-1)
  real   , intent(inout) :: press(0:nn-1,0:idim-1)
  integer, intent(in)    :: swi_solid(0:nn-1)

  ! subroutine_variable **********************************************
  integer :: j,l,k,m
  real    :: ddv0
  integer :: l_ref
  real    :: ww_ref
  real    :: rr,rr2,rr0
  real    :: r2pg
  real    :: ww
  real    :: ww_w
  real    :: r_kernel
  real    :: p_wall
  real   , parameter :: zero = 1.0e-10
  real    :: drr0(0:idim-1)
  real    :: drr(0:idim-1)
  real    :: ddv(0:idim-1)
  real    :: x0(0:idim-1)

  r_kernel = ker_r(2)

  ! pressure_term_calculation ****************************************
  if((cal_sw(i,0) .eq. 1) .and. (itypep(i) .ne. type(3)))then
    do l=1,neigh(i,0)
      j = neigh(i,l)
      if((itypep(j) .ge. type(2)) .or. &
      &  (itypep(j) .ne. type(3)) .or. &
      &  (itypep(j) .le. type(5)))then
        ! interparticle_distance
        rr2 = 0.0
        do m=0,idim-1
          drr(m) = x(i,m)-x(j,m)
          rr2   = rr2+drr(m)**2
        end do
        rr = sqrt(rr2)

        ! variable_setting
        r2pg     = max(zero,rho(i)/(c_vel**2))

        ! weight_function
        ww = r_kernel/max(rr,zero)-rr/max(r_kernel,zero)
        if(rr .gt. r_kernel)ww = 0.0

        ! pressure_term
        if(itypep(j) .eq. type(2) .and. swi_solid(j) .eq. 1)then
          ddv0 = real(idim)*(p(i))/max(zero,rr)*ww/r2pg/n0(i,2)
!        else if(itypep(j) .lt. phase)then
!!          ddv0 = real(idim)*(p(i))/max(zero,rr)*ww/r2pg/n0(i,2)
!          ddv0 = real(idim)*(2.0*p(i))/max(zero,rr)*ww/r2pg/n0(i,2)
        else if(itypep(j) .eq. phase .or. (itypep(j) .eq. type(2) .and. swi_solid(j) .eq. 0))then
          ddv0 = real(idim)*(p(i)+p(j))/max(zero,rr)*ww/r2pg/n0(i,2)
        else if(itypep(j) .gt. phase)then
!          ddv0 = real(idim)*(p(j))/max(zero,rr)*ww/r2pg/n0(i,2)
          ddv0 = real(idim)*(2.0*p(j))/max(zero,rr)*ww/r2pg/n0(i,2)
        else
          write(900,*)'error in pressure_term'
          stop
        end if

        do m=0,idim-1
          ddv(m) = ddv0*drr(m)/max(zero,rr)
        end do

        ! pressure_term_component
        do m=0,idim-1
          press(i,m) = press(i,m)+ddv(m)
        end do
      end if
    end do

    ! wall_influence ***********************************************
    if(g_wall(grid_p(i),0) .ne. 0)then
      rr = r_kernel

      do k=1,g_wall(grid_p(i),0)
        l = g_wall(grid_p(i),k)
        do m=0,idim-1
          x0(m) = x(i,m)
        end do

        call calc_dist                     &
        &    (x0,l,drr0,rr0,               &
        &     idim,wall_n,wall_anc,wall_nor)

        if(rr0 .lt. rr)then
          rr = rr0
          do m=0,idim-1
            drr(m) = drr0(m)
          end do
        end if
      end do

     ! data_copy
     do m=0,idim-1
       drr(m) = drr(m)/max(zero,rr)
       drr0(m) = drr(m)
     end do

      ! calculation_particle_number_density
      call calc_ww(7,rr,ww_ref,ww_num,ww_wall,dis)

      ! an_correction
      call correction_weight_function_for_pressure     &
      &    (i,drr,rr,nn,idim,x,grid_num,wall_n,grid_p, &
      &     dis,ww_ref,nn_wall,g_wall,wall_anc,wall_nor)

      ! variable_setting
      r2pg = max(zero,rho(i)/(c_vel**2))
      ddv0 = real(idim)*(p(i))/dis*ww_ref/r2pg/n0(i,2)

      do m=0,idim-1
        ddv(m) = ddv0*drr(m)
      end do

      ! pressure_term_component
      do m=0,idim-1
        press(i,m) = press(i,m)+ddv(m)
      end do
    end if
  end if

  return
end subroutine pressure_term

!*********************************************************************
!*** calculation_external_force_term                               ***
!*********************************************************************
subroutine external_force_term       &
&          (i,nn,idim,ifluid, &
&           accel,grav,nump,x,  &
&           num_surf,tau,ICMsurf, &
&           dis,den,phase, &
&           swi_hc)

  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i,nn,idim,ifluid
  integer, intent(in) :: nump
  real   , intent(in) :: grav,dis
  real   , intent(inout) :: accel(0:nn-1,0:idim-1)
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  integer, intent(in) :: num_surf
  real   , intent(in) :: tau(0:num_surf-1,0:1)
  real   , intent(in) :: ICMsurf(0:num_surf-1,0:1)
  real   , intent(in) :: den(0:ifluid)
  integer, intent(in) :: phase
  integer, intent(in) :: swi_hc(0:nn-1)

  ! subroutine_variable **********************************************
  integer :: m,i1,imin
  real    :: dista,min
  real,dimension(0:1) :: tau1

  ! array_initialize *************************************************
  do m=0,idim-1
   accel(i,m) = 0.0
  end do

  ! acceleration_calculation *****************************************
  if(swi_hc(i) .eq. 1) then

   min = 1.0e10
   do i1 = 0,num_surf-1 !粒子から最も近い格子点におけるtauを使う
    dista = (x(i,0)-ICMsurf(i1,0))**2+(x(i,1)-ICMsurf(i1,1))**2
    if(dista .lt. min) then
     min = dista
     imin = i1
    end if
   end do
   tau1(0) = tau(imin,0)/den(phase)/dis
   tau1(1) = tau(imin,1)/den(phase)/dis

   accel(i,0) = accel(i,0)+tau1(0)
   accel(i,1) = accel(i,1)+tau1(1)

  else
   accel(i,1) = accel(i,1)-grav
  end if

  return
end subroutine external_force_term

!*********************************************************************
!*** calculation_rho_&_nu                                          ***
!*********************************************************************
subroutine cal_rhonu                               &
&          (i,nn,idim,ifluid,neighbor,phase,num_ref, &
&           den,vis,rho,nu,nump,n0,an,wall_num)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i,nn,idim,ifluid,neighbor
  integer, intent(in) :: wall_num
  integer, intent(in) :: nump
  real   , intent(in) :: den(0:ifluid)
  real   , intent(in) :: vis(0:ifluid)
  real   , intent(in) :: n0(0:nn-1,0:3)
  real   , intent(in) :: an(0:nn-1)
  integer, intent(in) :: phase
  real   , intent(in) :: num_ref
  real   , intent(inout) :: rho(0:nn-1)
  real   , intent(inout) :: nu(0:nn-1)

  ! subroutine_variable **********************************************
  real, parameter   :: zero = 1.0e-10

  ! calculation_rho_&_nu *********************************************
  nu(i)  = vis(phase)
  rho(i) = max(an(i)/n0(i,1),1.0/num_ref)*den(phase)

  return
end subroutine cal_rhonu

!*********************************************************************
!*** set_rho_&_nu                                                  ***
!*********************************************************************
subroutine set_rhonu                    &
&           (i,nn,idim,ifluid,neighbor,   &
&            itypep,den,vis,rho,nu,nump)
  implicit none
  ! manroutine_variable **********************************************
  integer, intent(in) :: i,nn,idim,ifluid,neighbor
  integer, intent(in) :: nump
  integer, intent(in) :: itypep(0:nn-1)
  real   , intent(in) :: den(0:ifluid)
  real   , intent(in) :: vis(0:ifluid)
  real   , intent(inout) :: rho(0:nn-1)
  real   , intent(inout) :: nu(0:nn-1)

  ! subroutine_variable **********************************************

  ! set_rho_&_nu *****************************************************

    nu(i)  = vis(itypep(i))
    rho(i) = den(itypep(i))


  return
end subroutine set_rhonu

!*********************************************************************
!*** calculation_dt                                                ***
!*********************************************************************
subroutine cal_dt                          &
&          (nn,idim,ifluid,neighbor,       &
&           neigh,nump,itypep,type,cfl,x,  &
&           dt_fix,dt,time_sim,c_vel,dis,v,vv_max)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in)    :: nn,idim,ifluid,neighbor
  integer, intent(in)    :: neigh(0:nn-1,0:neighbor)
  integer, intent(in)    :: nump
  integer, intent(in)    :: itypep(0:nn-1)
  integer, intent(in)    :: type(0:ifluid)
  real   , intent(in)    :: cfl
  real   , intent(in)    :: dt_fix
  real   , intent(in)    :: x(0:nn-1,0:idim-1)
  real   , intent(in)    :: c_vel
  real   , intent(in)    :: dis
  real   , intent(in)    :: v(0:nn-1,0:idim-1)
  real   , intent(inout) :: dt
!  real   , intent(inout) :: time_sim
  double precision, intent(inout) :: time_sim
  real   , intent(inout)    :: vv_max

  ! subroutine_variable **********************************************
  integer :: i,j,l,m
  integer :: dt_sw
  real    :: vv_sca
  real    :: rr,rr2
  real    :: rr_min
  real    :: drr(0:idim-1)
  double precision :: time_sim_before

  ! switch_setting ***************************************************
  dt_sw = 1

  ! calculation_dt ***************************************************
  if(dt_sw .eq. 0)then
    dt = dt_fix
  else if(dt_sw .eq. 1)then
    rr_min = 1.0e10
    vv_max = 1.0e-10

    do i=0,nump-1
      if(itypep(i) .ge. type(3))then
        do l=1,neigh(i,0)
          j = neigh(i,l)
          if(itypep(j) .ge. type(2))then
            rr2 = 0.0
            do m=0,idim-1
             drr(m) = x(j,m)-x(i,m)
             rr2   = rr2+drr(m)**2
            end do
            rr = sqrt(rr2)
            rr_min = min(rr_min,rr)
          end if
        end do

         vv_sca = 0.0
         do m=0,idim-1
           vv_sca = vv_sca+v(i,m)**2
         end do

         vv_sca = sqrt(vv_sca)
         vv_max = max(vv_max,vv_sca)

      end if
    end do

    ! calculation_dt
    dt = rr_min/max(vv_max/cfl,c_vel)
  else
    dt = dis/c_vel
  end if

!  if(dt .gt. dt_fix)then
  if(dt .gt. 1.0)then
    dt = dt_fix
  end if

  ! simulation_time
  time_sim_before = time_sim
  time_sim = time_sim+dble(dt)

  if((time_sim - time_sim_before) .lt. 1.0d-15)then
    write(*,*) 'Error in subroutine cal_dt'
    write(*,*) dt,time_sim,time_sim_before
    stop
  end if

  ! delete **********************************************************
  if(dt .eq. 0.0)then
    write(*,*)'in subroutine cal_dt'
    write(*,*)dt,vv_max,cfl
    stop
  end if
  ! delete **********************************************************

  return
end subroutine cal_dt

!!*********************************************************************
!!*** paraview_output(binary)                                       ***
!!*********************************************************************
!subroutine output_para_kmgt           &
!&          (nn,idim,ifluid,fcount,t,  &
!&           x,v,itypep,type,an,p,surf)
!  implicit none
!  ! mainroutine_variable *********************************************
!  integer   , intent(in) :: nn,idim,ifluid
!  integer   , intent(in) :: fcount
!  real      , intent(in) :: x(0:nn-1,0:idim-1)
!  integer   , intent(in) :: itypep(0:nn-1)
!  integer   , intent(in) :: type(0:ifluid)
!  real      , intent(in) :: an(0:nn-1)
!  real      , intent(in) :: p(0:nn-1)
!  real      , intent(in) :: surf(0:nn-1,0:idim-1)
!  real      , intent(in) :: v(0:nn-1,0:idim-1)
!  real      , intent(in) :: t(0:nn-1)
!
!  ! subroutine_variable **********************************************
!  integer :: i,j
!  type(vtk_file_handle):: fd
!  real    , allocatable :: x_sub(:)
!  real    , allocatable :: y_sub(:)
!  real    , allocatable :: z_sub(:)
!  real    , allocatable :: u_sub(:)
!  real    , allocatable :: v_sub(:)
!  real    , allocatable :: w_sub(:)
!  real    , allocatable :: an_sub(:)
!  real    , allocatable :: p_sub(:)
!  real    , allocatable :: surx_sub(:)
!  real    , allocatable :: sury_sub(:)
!  real    , allocatable :: surz_sub(:)
!  integer :: output_type
!  integer :: out_count
!
!  ! initial_setting **************************************************
!  output_type = type(3)
!
!  ! array_allocation *************************************************
!  ! particles_count
!  out_count = 0
!
!   do i=0,nn-1
!     if(itypep(i) .ge. output_type)then
!       out_count = out_count+1
!     end if
!   end do
!
!   ! array_allocation
!   allocate(x_sub(0:out_count-1))
!   allocate(y_sub(0:out_count-1))
!   allocate(z_sub(0:out_count-1))
!   allocate(u_sub(0:out_count-1))
!   allocate(v_sub(0:out_count-1))
!   allocate(w_sub(0:out_count-1))
!   allocate(p_sub(0:out_count-1))
!   allocate(an_sub(0:out_count-1))
!   allocate(surx_sub(0:out_count-1))
!   allocate(sury_sub(0:out_count-1))
!   allocate(surz_sub(0:out_count-1))
!
!  ! output_vtk *******************************************************
!  ! data_copy
!  out_count = 0
!
!  do i=0,nn-1
!    if(itypep(i) .ge. output_type)then
!      x_sub(out_count) = x(i,0)
!      y_sub(out_count) = x(i,1)
!
!      if(idim .eq. 3)then
!        z_sub(out_count) = x(i,2)
!      else if(idim .eq. 2)then
!        z_sub(out_count) = 0.0
!      end if
!
!      u_sub(out_count) = v(i,0)
!      v_sub(out_count) = v(i,1)
!
!      if(idim .eq. 3)then
!        w_sub(out_count) = v(i,2)
!      else if(idim .eq. 2)then
!        w_sub(out_count) = 0.0
!      end if
!
!      an_sub(out_count) = an(i)
!      p_sub(out_count)  = p(i)
!
!      surx_sub(out_count) = surf(i,0)
!      sury_sub(out_count) = surf(i,1)
!
!      if(idim .eq. 3)then
!        surz_sub(out_count) = surf(i,2)
!      else if(idim .eq. 2)then
!        surz_sub(out_count) = 0.0
!      end if
!
!      out_count = out_count+1
!    end if
!  end do
!
!  ! output_vtk *******************************************************
!  call vtk_rgcy_ini( fd             = fd,             &
!  &                  output_format  = 'binary',  &
!  &                  output_endian  = 'big_endian',  &
!  &                  root_directory = './result/',  &
!  &                  prefix         = 'mps', &
!  &                  file_extension = '.vtk', &
!  &                  time_step      = fcount,         &
!  &                  hedder         = 'output',       &
!  &                  mesh_topology  = 'polydata' )
!  call vtk_rgcy_geo( fd    = fd,       &
!  &                  x     = x_sub(:), &
!  &                  y     = y_sub(:), &
!  &                  z     = z_sub(:))
!!  call vtk_rgcy_var( fd    = fd,                       &
!!  &                  name  = "particle_velocity[m/s]", &
!!  &                  vx    = u_sub(:),                 &
!!  &                  vy    = v_sub(:),                 &
!!  &                  vz    = w_sub(:) )
!
!  call vtk_rgcy_end( fd = fd )
!
!  ! array_deallocation ***********************************************
!  deallocate(x_sub)
!  deallocate(y_sub)
!  deallocate(z_sub)
!  deallocate(u_sub)
!  deallocate(v_sub)
!  deallocate(w_sub)
!  deallocate(p_sub)
!  deallocate(an_sub)
!  deallocate(surx_sub)
!  deallocate(sury_sub)
!  deallocate(surz_sub)
!
!  return
!end subroutine output_para_kmgt

subroutine output_para_bin            &
&          (nn,idim,ifluid,fcount,t,  &
&           x,v,itypep,type,an,p,surf,cal_type0,swi_koeki,ganma, &
&           swi_solid,swi_hc)
  implicit none
  ! mainroutine_variable *********************************************
  integer   , intent(in) :: nn,idim,ifluid
  integer   , intent(in) :: fcount
  real      , intent(in) :: x(0:nn-1,0:idim-1)
  integer   , intent(in) :: itypep(0:nn-1)
  integer   , intent(in) :: type(0:ifluid)
  real      , intent(in) :: an(0:nn-1)
  real      , intent(in) :: p(0:nn-1)
  real      , intent(in) :: surf(0:nn-1,0:idim-1)
  real      , intent(in) :: v(0:nn-1,0:idim-1)
  real      , intent(in) :: t(0:nn-1)
  integer, intent(in)    :: cal_type0(0:nn-1)
  integer, intent(in)    :: swi_koeki(0:nn-1)
  real, intent(in)    :: ganma(0:nn-1)
  integer, intent(in)    :: swi_solid(0:nn-1)
  integer, intent(in)    :: swi_hc(0:nn-1)

  ! subroutine_variable **********************************************
  integer             :: i,j,n
  integer             :: out_count
  integer             :: output_type
  real, allocatable   :: out_real_1(:,:)
  real, allocatable   :: out_real_2(:,:)
  real, allocatable   :: out_real_3(:,:)
  real, allocatable   :: out_real_4(:,:)
  real, allocatable   :: out_real_5(:,:)
  real, allocatable   :: out_real_6(:,:)
  real, allocatable   :: out_real_7(:,:)
  real, allocatable   :: out_real_8(:,:)
  real, allocatable   :: out_real_9(:,:)
  real, allocatable   :: out_real_10(:,:)
  real, allocatable   :: out_real_11(:,:)

  character(len =  40) :: datfile
  character(len = 200) :: number_com
  character(len =  10), parameter :: dir = 'result/MPS'
  character(len =  10), parameter :: dat = 'mps'
  character(len =   1), parameter :: newline = char(10)

  ! array_allocate ***************************************************
  allocate(out_real_1(0:2,0:nn-1))
  allocate(out_real_2(0:0,0:nn-1))
  allocate(out_real_3(0:0,0:nn-1))
  allocate(out_real_4(0:2,0:nn-1))
  allocate(out_real_5(0:2,0:nn-1))
  allocate(out_real_6(0:0,0:nn-1))
  allocate(out_real_7(0:0,0:nn-1))
  allocate(out_real_8(0:0,0:nn-1))
  allocate(out_real_9(0:0,0:nn-1))
  allocate(out_real_10(0:0,0:nn-1))
  allocate(out_real_11(0:0,0:nn-1))

  ! array_allocation *************************************************
  out_real_1(0:2,0:nn-1) = 0.0
  out_real_2(0:0,0:nn-1) = 0.0
  out_real_3(0:0,0:nn-1) = 0.0
  out_real_4(0:2,0:nn-1) = 0.0
  out_real_5(0:2,0:nn-1) = 0.0
  out_real_6(0:0,0:nn-1) = 0.0
  out_real_7(0:0,0:nn-1) = 0.0
  out_real_8(0:0,0:nn-1) = 0.0
  out_real_9(0:0,0:nn-1) = 0.0
  out_real_10(0:0,0:nn-1) = 0.0
  out_real_11(0:0,0:nn-1) = 0.0

  ! input_variable ***************************************************
  out_count   = 0
  output_type = type(2)

  ! output_particle_number_count *************************************
  do n=0,nn-1
    if(itypep(n) .ge. output_type)then
      ! position
      do i=0,2
        if(i <= idim-1)then
          out_real_1(i,out_count) = x(n,i)
        else
          out_real_1(i,out_count) = 0.0
        end if
      end do

      ! velocity
      do i=0,2
        if(i <= idim-1)then
          out_real_5(i,out_count) = v(n,i)
        else
          out_real_5(i,out_count) = 0.0
        end if
      end do

      ! particle_type
      out_real_6(0,out_count) = real(itypep(n))

      !particle_type_transient_mode（凝固粒子判定）
      if      (real(itypep(n)) .eq. type(2) )  then
        out_real_8(0,out_count) = 2.0
      else if ((real(itypep(n)) .eq. type(4) ) .and. ( swi_koeki(n) .eq. 0 ))then
        out_real_8(0,out_count) = 3.0
      else if ((real(itypep(n)) .eq. type(4) ) .and. ( swi_koeki(n) .eq. 1 ))then
        out_real_8(0,out_count) = 4.0
      end if

      !liquid_fraction
      out_real_9(0,out_count) = ganma(n)

      ! pressure
      out_real_2(0,out_count) = p(n)

      ! particle_numer_density
      out_real_3(0,out_count) = an(n)

      ! surface_tension
      do i=0,2
        if(i <= idim-1)then
          out_real_4(i,out_count) = surf(n,i)
        else
          out_real_4(i,out_count) = 0.0
        end if
      end do

      ! temperature
      out_real_7(0,out_count) = t(n)

      ! switch_solidification
      out_real_10(0,out_count) = swi_solid(n)

      ! switch_hc
      out_real_11(0,out_count) = swi_hc(n)

      ! count_up
      out_count = out_count+1
    end if

    if((itypep(n) .eq. type(0)) .and. &
    &  (cal_type0(n) .ge. 1))then
      ! position
      do i=0,2
        if(i <= idim-1)then
          out_real_1(i,out_count) = x(n,i)
        else
          out_real_1(i,out_count) = 0.0
        end if
      end do

      ! velocity
      do i=0,2
        if(i <= idim-1)then
          out_real_5(i,out_count) = v(n,i)
        else
          out_real_5(i,out_count) = 0.0
        end if
      end do

      ! particle_type
      out_real_6(0,out_count) = real(itypep(n))

      ! pressure
      out_real_2(0,out_count) = p(n)

      ! particle_numer_density
      out_real_3(0,out_count) = an(n)

      ! surface_tension
      do i=0,2
        if(i <= idim-1)then
          out_real_4(i,out_count) = surf(n,i)
        else
          out_real_4(i,out_count) = 0.0
        end if
      end do

      ! temperature
      out_real_7(0,out_count) = t(n)

      ! count_up
      out_count = out_count+1
    end if
  end do

  ! file_name_setting ************************************************
  call file_name(dir,dat,datfile,fcount)

  ! particle_number_check ********************************************
  if(out_count > 100000000)then
    write(900,'(a)')"Error in subroutine output_para"
    write(900,'(a)')"out_count is too large"
    stop
  end if

  ! outputfile *******************************************************
!BINARY
  open(unit       = 12,             &
  &    file       = datfile,        &
  &    form       = 'unformatted',  &
! &    access     = 'sequential',   &
  &    access     = 'stream',       &
  &    convert    = 'big_endian',   &
  &    status     = 'replace',      &
  &    action     = 'write')
!ASCII
!  open(unit       = 12,             &
!  &    file       = datfile,        &
!  &    form       = 'formatted',    &
!  &    status     = 'replace')

  write(number_com,*)out_count

!BINARY
  write(12)'# vtk DataFile Version 3.0'//newline
  write(12)'vtk output'//newline
  write(12)'BINARY'//newline
  write(12)'DATASET POLYDATA'//newline
  write(12)'POINTS '//trim(adjustl(number_com))//' float'//newline
!ASCII
!  write(12,'(a)')'# vtk DataFile Version 3.0'//newline
!  write(12,'(a)')'checkhedder'//newline
!  write(12,'(a)')'ASCII'//newline
!  write(12,'(a)')'DATASET POLYDATA'//newline
!  write(12,'(a)')'POINTS'
!  write(12,*)trim(adjustl(number_com))
!  write(12,'(a)')'double'

!BINARY
  ! particle_position
  write(12)(out_real_1(0,j),out_real_1(1,j),out_real_1(2,j),j=0,out_count-1)
  write(12)newline
  write(12)'POINT_DATA '//trim(adjustl(number_com))//newline

  ! particle_type
  write(12)'SCALARS particle_type float 1'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_real_6(0,j),j=0,out_count-1)
  write(12)newline

  !particle_type_transient_mode（凝固粒子判定）
  write(12)'SCALARS particle_type_transient_mode float 1'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_real_8(0,j),j=0,out_count-1)
  write(12)newline

  !liquid_fraction
  write(12)'SCALARS liquid_fraction float 1'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_real_9(0,j),j=0,out_count-1)
  write(12)newline

  ! velocity
  write(12)'SCALARS velocity float 3'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_real_5(0,j),out_real_5(1,j),out_real_5(2,j),j=0,out_count-1)
  write(12)newline

  ! pressure
  write(12)'SCALARS pressure float 1'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_real_2(0,j),j=0,out_count-1)
  write(12)newline

  ! particle_number_density
  write(12)'SCALARS num_density float 1'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_real_3(0,j),j=0,out_count-1)
  write(12)newline

  ! surface_tension
  write(12)'SCALARS surface_tension float 3'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_real_4(0,j),out_real_4(1,j),out_real_4(2,j),j=0,out_count-1)
  write(12)newline

  ! temperature
  write(12)'SCALARS temperature float 1'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_real_7(0,j),j=0,out_count-1)
  write(12)newline

  ! switch_solidification
  write(12)'SCALARS swi_solid float 1'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_real_10(0,j),j=0,out_count-1)
  write(12)newline

  ! switch_hc
  write(12)'SCALARS swi_hc float 1'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_real_11(0,j),j=0,out_count-1)
  write(12)newline

!ASCII
!  ! particle_position
!  do j=0,out_count-1
!    write(12,*)out_real_1(0,j),out_real_1(1,j),out_real_1(2,j)
!  end do
!  write(12,*)newline
!  write(12,*)'POINT_DATA '//trim(adjustl(number_com))//newline
!
!  ! particle_type
!  write(12,'(a)')'SCALARS particle_type float 1'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  do j=0,out_count-1
!  write(12,*)out_real_6(0,j)
!  end do
!  write(12,*)newline
!
!  !particle_type_transient_mode（凝固粒子判定）
!  write(12,'(a)')'SCALARS particle_type_transient_mode float 1'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  do j=0,out_count-1
!  write(12,*)out_real_8(0,j)
!  end do
!  write(12,*)newline
!
!  !liquid_fraction
!  write(12,'(a)')'SCALARS liquid_fraction float 1'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  do j=0,out_count-1
!  write(12,*)out_real_9(0,j)
!  end do
!  write(12,*)newline
!
!  ! velocity
!  write(12,'(a)')'SCALARS velocity float 3'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  do j=0,out_count-1
!  write(12,*)out_real_5(0,j),out_real_5(1,j),out_real_5(2,j)
!  end do
!  write(12,*)newline
!
!  ! pressure
!  write(12,'(a)')'SCALARS pressure float 1'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  do j=0,out_count-1
!  write(12,*)out_real_2(0,j)
!  end do
!  write(12,*)newline
!
!  ! particle_number_density
!  write(12,'(a)')'SCALARS num_density float 1'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  do j=0,out_count-1
!  write(12,*)out_real_3(0,j)
!  end do
!  write(12,*)newline
!
!  ! surface_tension
!  write(12,'(a)')'SCALARS surface_tension float 3'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  do j=0,out_count-1
!  write(12,*)out_real_4(0,j),out_real_4(1,j),out_real_4(2,j)
!  end do
!  write(12,*)newline
!
!  ! temperature
!  write(12,'(a)')'SCALARS temperature float 1'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  do j=0,out_count-1
!  write(12,*)out_real_7(0,j)
!  end do
!  write(12,*)newline

  close(12)

  ! array_deallocation **********************************************
  deallocate(out_real_1)
  deallocate(out_real_2)
  deallocate(out_real_3)
  deallocate(out_real_4)
  deallocate(out_real_5)
  deallocate(out_real_6)
  deallocate(out_real_7)
  deallocate(out_real_8)
  deallocate(out_real_9)
  deallocate(out_real_10)
  deallocate(out_real_11)

  return
end subroutine output_para_bin

!*********************************************************************
!*** PARAVIEW_output_for_terms(binary)                             ***
!*********************************************************************
subroutine output_term                   &
&          (nn,idim,ifluid,fcount,type,  &
&           x,itypep,surf,visc,accel,press)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: nn,idim,ifluid
  integer, intent(in) :: fcount
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  integer, intent(in) :: itypep(0:nn-1)
  integer, intent(in) :: type(0:ifluid)
  real   , intent(in) :: surf(0:nn-1,0:idim-1)
  real   , intent(in) :: visc(0:nn-1,0:idim-1)
  real   , intent(in) :: accel(0:nn-1,0:idim-1)
  real   , intent(in) :: press(0:nn-1,0:idim-1)

  ! subroutine_variable **********************************************
  integer           :: i,j,n
  integer           :: out_count
  integer           :: output_type
  real, allocatable :: out_posi(:,:)
  real, allocatable :: out_visc(:,:)
  real, allocatable :: out_surf(:,:)
  real, allocatable :: out_press(:,:)
  real, allocatable :: out_accel(:,:)
  real, allocatable :: out_type(:,:)
  character(len =  40) :: datfile
  character(len = 200) :: number_com
  character(len =  10), parameter :: dir = 'result'
  character(len =  10), parameter :: dat = 'term'
  character(len =   1), parameter :: newline = char(10)

  ! array_allocate ***************************************************
  allocate(out_posi(0:2,0:nn-1))
  allocate(out_visc(0:2,0:nn-1))
  allocate(out_surf(0:2,0:nn-1))
  allocate(out_press(0:2,0:nn-1))
  allocate(out_accel(0:2,0:nn-1))
  allocate(out_type(0:0,0:nn-1))

  ! array_allocation *************************************************
  out_posi(0:2,0:nn-1) = 0.0
  out_visc(0:2,0:nn-1) = 0.0
  out_surf(0:2,0:nn-1) = 0.0
  out_press(0:2,0:nn-1) = 0.0
  out_accel(0:2,0:nn-1) = 0.0
  out_type(0:0,0:nn-1)  = 0.0

  ! input_variable ***************************************************
  out_count   = 0
  output_type = type(2)

  ! output_particle_number_count *************************************
  do n=0,nn-1
    if(itypep(n) .ge. output_type)then
      ! position
      do i=0,2
        if(i <= idim-1)then
          out_posi(i,out_count) = x(n,i)
        else
          out_posi(i,out_count) = 0.0
        end if
      end do

      ! particle_type
      out_type(0,out_count) = real(itypep(n))

      ! pressure
      do i=0,2
        if(i <= idim-1)then
          out_press(i,out_count) = press(n,i)
        else
          out_press(i,out_count) = 0.0
        end if
      end do

      ! surface_tension
      do i=0,2
        if(i <= idim-1)then
          out_surf(i,out_count) = surf(n,i)
        else
          out_surf(i,out_count) = 0.0
        end if
      end do

      ! external_force
      do i=0,2
        if(i <= idim-1)then
          out_accel(i,out_count) = accel(n,i)
        else
          out_accel(i,out_count) = 0.0
        end if
      end do

      ! count_up
      out_count = out_count+1
    end if
  end do

  ! file_name_setting ************************************************
  call file_name(dir,dat,datfile,fcount)

  ! outputfile *******************************************************
!BINARY
  open(unit       = 12,             &
  &    file       = datfile,        &
  &    form       = 'unformatted',  &
! &    access     = 'sequential',   &
  &    access     = 'stream',       &
  &    convert    = 'big_endian',   &
  &    status     = 'replace',      &
  &    action     = 'write')
!ASCII
!  open(unit       = 12,             &
!  &    file       = datfile,        &
!  &    form       = 'formatted',    &
!  &    status     = 'replace')

  write(number_com,*)out_count

!BINARY
  write(12)'# vtk DataFile Version 3.0'//newline
  write(12)'vtk output'//newline
  write(12)'BINARY'//newline
  write(12)'DATASET POLYDATA'//newline
  write(12)'POINTS '//trim(adjustl(number_com))//' float'//newline

  ! particle_position
  write(12)(out_posi(0,j),out_posi(1,j),out_posi(2,j),j=0,out_count-1)
  write(12)newline
  write(12)'POINT_DATA '//trim(adjustl(number_com))//newline

  ! particle_type
  write(12)'SCALARS particle_type float 1'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_type(0,j),j=0,out_count-1)
  write(12)newline

   ! viscous_term
  write(12)'SCALARS viscous_term float 3'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_visc(0,j),out_visc(1,j),out_visc(2,j),j=0,out_count-1)
  write(12)newline

 ! viscous_term_magnitude
  write(12)'SCALARS viscous_term_mag float 1'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(sqrt(out_visc(0,j)**2+out_visc(1,j)**2+out_visc(2,j)),j=0,out_count-1)
  write(12)newline

  ! pressure_term
  write(12)'SCALARS pressure_term float 3'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_press(0,j),out_press(1,j),out_press(2,j),j=0,out_count-1)
  write(12)newline

  ! pressure_term_magnitude
  write(12)'SCALARS pressure_term_term_mag float 1'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(sqrt(out_press(0,j)**2+out_press(1,j)**2+out_press(2,j)**2),j=0,out_count-1)
  write(12)newline

  ! external_force_term
  write(12)'SCALARS external_force_term float 3'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_accel(0,j),out_accel(1,j),out_accel(2,j),j=0,out_count-1)
  write(12)newline

  ! external_force_term_magnitude
  write(12)'SCALARS external_force_term_mag float 1'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(sqrt(out_accel(0,j)**2+out_accel(1,j)**2+out_accel(2,j)**2),j=0,out_count-1)
  write(12)newline

  ! surface_tension_term
  write(12)'SCALARS surface_tension_term float 3'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(out_surf(0,j),out_surf(1,j),out_surf(2,j),j=0,out_count-1)
  write(12)newline

  ! surface_tension_term_magnitude
  write(12)'SCALARS surface_tension_term_mag float 1'//newline
  write(12)'LOOKUP_TABLE default'//newline
  write(12)(sqrt(out_surf(0,j)**2+out_surf(1,j)**2+out_surf(2,j)**2),j=0,out_count-1)
  write(12)newline

!ASCII
!  write(12,'(a)')'# vtk DataFile Version 3.0'//newline
!  write(12,'(a)')'vtk output'//newline
!  write(12,'(a)')'BINARY'//newline
!  write(12,'(a)')'DATASET POLYDATA'//newline
!  write(12,*)'POINTS '//trim(adjustl(number_com))//' float'//newline
!
!  ! particle_position
!  write(12,*)(out_posi(0,j),out_posi(1,j),out_posi(2,j),j=0,out_count-1)
!  write(12,*)newline
!  write(12,*)'POINT_DATA '//trim(adjustl(number_com))//newline
!
!  ! particle_type
!  write(12,'(a)')'SCALARS particle_type float 1'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  write(12,*)(out_type(0,j),j=0,out_count-1)
!  write(12,*)newline
!
!   ! viscous_term
!  write(12,'(a)')'SCALARS viscous_term float 3'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  write(12,*)(out_visc(0,j),out_visc(1,j),out_visc(2,j),j=0,out_count-1)
!  write(12,*)newline
!
! ! viscous_term_magnitude
!  write(12,'(a)')'SCALARS viscous_term_mag float 1'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  write(12,*)(sqrt(out_visc(0,j)**2+out_visc(1,j)**2+out_visc(2,j)),j=0,out_count-1)
!  write(12,*)newline
!
!  ! pressure_term
!  write(12,'(a)')'SCALARS pressure_term float 3'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  write(12,*)(out_press(0,j),out_press(1,j),out_press(2,j),j=0,out_count-1)
!  write(12,*)newline
!
!  ! pressure_term_magnitude
!  write(12,'(a)')'SCALARS pressure_term_term_mag float 1'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  write(12,*)(sqrt(out_press(0,j)**2+out_press(1,j)**2+out_press(2,j)**2),j=0,out_count-1)
!  write(12,*)newline
!
!  ! external_force_term
!  write(12,'(a)')'SCALARS external_force_term float 3'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  write(12,*)(out_accel(0,j),out_accel(1,j),out_accel(2,j),j=0,out_count-1)
!  write(12,*)newline
!
!  ! external_force_term_magnitude
!  write(12,'(a)')'SCALARS external_force_term_mag float 1'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  write(12,*)(sqrt(out_accel(0,j)**2+out_accel(1,j)**2+out_accel(2,j)**2),j=0,out_count-1)
!  write(12,*)newline
!
!  ! surface_tension_term
!  write(12,'(a)')'SCALARS surface_tension_term float 3'//newline
!  write(12,,'(a)')'LOOKUP_TABLE default'//newline
!  write(12,*)(out_surf(0,j),out_surf(1,j),out_surf(2,j),j=0,out_count-1)
!  write(12,*)newline
!
!  ! surface_tension_term_magnitude
!  write(12,'(a)')'SCALARS surface_tension_term_mag float 1'//newline
!  write(12,'(a)')'LOOKUP_TABLE default'//newline
!  write(12,*)(sqrt(out_surf(0,j)**2+out_surf(1,j)**2+out_surf(2,j)**2),j=0,out_count-1)
!  write(12,*)newline

  close(12)

  ! array_deallocation **********************************************
  deallocate(out_posi)
  deallocate(out_visc)
  deallocate(out_surf)
  deallocate(out_press)
  deallocate(out_accel)
  deallocate(out_type)

  return
end subroutine output_term

!*********************************************************************
!*** output_paraview                                               ***
!*********************************************************************
subroutine output_para                &
&          (nn,idim,ifluid,fcount,t,  &
&           x,v,itypep,type,an,p,surf)
  implicit none
  ! mainroutine_variable *********************************************
  integer   , intent(in) :: nn,idim,ifluid
  integer   , intent(in) :: fcount
  real      , intent(in) :: x(0:nn-1,0:idim-1)
  integer   , intent(in) :: itypep(0:nn-1)
  integer   , intent(in) :: type(0:ifluid)
  real      , intent(in) :: an(0:nn-1)
  real      , intent(in) :: p(0:nn-1)
  real      , intent(in) :: surf(0:nn-1,0:idim-1)
  real      , intent(in) :: v(0:nn-1,0:idim-1)
  real      , intent(in) :: t(0:nn-1)

  ! subroutine_variable **********************************************
  integer             :: i,n
  integer             :: out_count
  integer             :: output_type
  real, allocatable   :: out_real_1(:,:)
  real, allocatable   :: out_real_2(:,:)
  real, allocatable   :: out_real_3(:,:)
  real, allocatable   :: out_real_4(:,:)
  real, allocatable   :: out_real_5(:,:)
  real, allocatable   :: out_real_6(:,:)
  real, allocatable   :: out_real_7(:,:)
  character(len = 40) :: datfile
  character(len = 10), parameter :: dir = 'data'
  character(len = 10), parameter :: dat = 'mps'

  ! array_allocation *************************************************
  allocate(out_real_1(0:2,0:nn-1))
  allocate(out_real_2(0:0,0:nn-1))
  allocate(out_real_3(0:0,0:nn-1))
  allocate(out_real_4(0:2,0:nn-1))
  allocate(out_real_5(0:2,0:nn-1))
  allocate(out_real_6(0:0,0:nn-1))
  allocate(out_real_7(0:0,0:nn-1))

  ! array_initialize *************************************************
  out_real_1(0:2,0:nn-1) = 0.0
  out_real_2(0:0,0:nn-1) = 0.0
  out_real_3(0:0,0:nn-1) = 0.0
  out_real_4(0:2,0:nn-1) = 0.0
  out_real_5(0:2,0:nn-1) = 0.0
  out_real_6(0:0,0:nn-1) = 0.0
  out_real_7(0:0,0:nn-1) = 0.0

  ! initial_setting **************************************************
  out_count   = 0
  output_type = type(1)

  ! data_copy ********************************************************
  do n=0,nn-1
    if(itypep(n) >= output_type)then
      ! position
      do i=0,2
        if(i <= idim-1)then
          out_real_1(i,out_count) = x(n,i)
        else
          out_real_1(i,out_count) = 0.0
        end if
      end do

      ! velocity
      do i=0,2
        if(i <= idim-1)then
          out_real_5(i,out_count) = v(n,i)
        else
          out_real_5(i,out_count) = 0.0
        end if
      end do

      ! particle_type
      out_real_6(0,out_count) = real(itypep(n))

      ! pressure
      out_real_2(0,out_count) = p(n)

      ! particle_number_density
      out_real_3(0,out_count) = an(n)

      ! surface_tension
      do i=0,2
        if(i <= idim-1)then
          out_real_4(i,out_count) = surf(n,i)
        else
          out_real_4(i,out_count) = 0.0
        end if
      end do

      ! temperature
      out_real_7(0,out_count) = t(n)

      ! count_up
      out_count = out_count+1
    end if
  end do

  ! file_name ********************************************************
  call file_name(dir,dat,datfile,fcount)

  ! particle_number_check ********************************************
  if(out_count > 100000000)then
    write(900,'(a)')"Error in subroutine output_para"
    write(900,'(a)')"out_count is too large"
    stop
  end if

  ! output_file ******************************************************
  open(12,file = datfile)
  write(12,'(a)')'# vtk DataFile Version 3.0'
  write(12,'(a)')'vtk output'
  write(12,'(a)')'ASCII'
  write(12,'(a)')'DATASET POLYDATA'
  write(12,'(a,i8,a)')'POINTS',out_count,' float'

  do n=0,out_count-1
    write(12,*)out_real_1(0:2,n)
  end do

  write(12,'(a,i9)')'POINT_DATA',out_count
  write(12,'(a   )')'FIELD attributes 6'

  write(12,'(a,i9,a)')'particle_type  1',out_count,'float'
  do n=0,out_count-1
    write(12,*)out_real_6(0,n)
  end do

  write(12,'(a,i9,a)')'velocity  3',out_count,'float'
  do n=0,out_count-1
    write(12,*)out_real_5(0:2,n)
  end do

  write(12,'(a,i9,a)')'pressure  1',out_count,'float'
  do  n=0,out_count-1
    write(12,*)out_real_2(0,n)
  end do

  write(12,'(a,i9,a)')'particle_number_density  1',out_count,'float'
  do  n=0,out_count-1
    write(12,*)out_real_3(0,n)
  end do

  write(12,'(a,i9,a)')'surface_tension  3',out_count,'float'
  do  n=0,out_count-1
    write(12,*)out_real_4(0:2,n)
  end do

  write(12,'(a,i9,a)')'temperature  1',out_count,'float'
  do  n=0,out_count-1
    write(12,*)out_real_7(0,n)
  end do

  close(12)

  ! array_deallocation ***********************************************
  deallocate(out_real_1)
  deallocate(out_real_2)
  deallocate(out_real_3)
  deallocate(out_real_4)

  return
end subroutine output_para

!*********************************************************************
!*** output_an_paraview                                            ***
!*********************************************************************
subroutine output_an         &
&          (nn,idim,ifluid,  &
&           x,itypep,type,an)
  implicit none
  ! mainroutine_variable *********************************************
  integer   , intent(in) :: nn,idim,ifluid
  real      , intent(in) :: x(0:nn-1,0:idim-1)
  integer   , intent(in) :: itypep(0:nn-1)
  integer   , intent(in) :: type(0:ifluid)
  real      , intent(in) :: an(0:nn-1)

  ! subroutine_variable **********************************************
  integer             :: i,n
  integer             :: out_count
  integer             :: output_type
  real, allocatable   :: out_real_1(:,:)
  real, allocatable   :: out_real_2(:,:)
  character(len = 40) :: datfile

  ! array_allocation *************************************************
  allocate(out_real_1(0:2,1:nn-1))
  allocate(out_real_2(0:0,1:nn-1))

  ! array_initialize *************************************************
  out_real_1(0:2,1:nn-1) = 0.0
  out_real_2(0:0,1:nn-1) = 0.0

  ! variable_initialize **********************************************
  out_count   = 0
  output_type = type(2)

  ! data_copy ********************************************************
  do n=0,nn-1
    if(itypep(n) >= output_type)then
      out_count = out_count+1

      ! particle_poisition
      do i=0,2
        if(i <= idim-1)then
          out_real_1(i,out_count) = x(n,i)
        else
          out_real_1(i,out_count) = 0.0
        end if
      end do

      ! particle_number_density
      out_real_2(0,out_count) = an(n)

    end if
  end do

  ! output_file_name *************************************************
  write(datfile,'(a)')'./result/output_an.vtk'

  ! output_vtk_file **************************************************
  open(12,file = datfile)
  write(12,'(a)')'# vtk DataFile Version 3.0'
  write(12,'(a)')'vtk output'
  write(12,'(a)')'ASCII'
  write(12,'(a)')'DATASET POLYDATA'
  write(12,'(a,i8,a)')'POINTS',out_count,' float'

  do n=1,out_count
    write(12,*)out_real_1(0:2,n)
  end do

  write(12,'(a,i9)')'POINT_DATA',out_count
  write(12,'(a   )')'FIELD attributes 1'

  write(12,'(a,i9,a)')'pressure  1',out_count,'float'
  do n=1,out_count
    write(12,'(f10.3)')out_real_2(0,n)
  end do

  close(12)

  ! array_deallocation ***********************************************
  deallocate(out_real_1)
  deallocate(out_real_2)

  ! calculation_stop *************************************************
  write(*,*)'stop at output_an'
  stop

  return
end subroutine output_an

!*********************************************************************
!*** set_neighboring_particle                                      ***
!*********************************************************************
subroutine set_neigh                   &
&          (i,                         &
&           nn,idim,ifluid,neighbor,x, &
&           type,neigh,itypep,         &
&           ker_r,nn_grid,             &
&           grid_num,grid_p,g_neigh)
  implicit none
  ! mainroutine_variable ********************************************
  integer, intent(in)    :: i
  integer, intent(in)    :: nn,idim,ifluid,neighbor
  integer, intent(in)    :: type(0:ifluid)
  integer, intent(in)    :: itypep(0:nn-1)
  real   , intent(in)    :: ker_r(0:3)
  real   , intent(in)    :: x(0:nn-1,0:idim-1)
  integer, intent(in)    :: grid_num
  integer, intent(in)    :: grid_p(0:nn-1)
  integer, intent(in)    :: nn_grid
  integer, intent(in)    :: g_neigh(1:grid_num,0:nn_grid)
  integer, intent(inout) :: neigh(0:nn-1,0:neighbor)

  ! subroutine_variable *********************************************
  integer :: j,k,l,m
  real    :: rr,rr2
  real    :: drr(0:idim-1)

  ! array_initialize ************************************************
  neigh(i,0)= 0

  ! neighboring_particle_count **************************************
  ! particle_number_copy
  k = grid_p(i)

  if(itypep(i) .ge. type(2))then
    ! particle_count
    do l=1,g_neigh(k,0)
      j = g_neigh(k,l)
      ! particle_excepting
      if(i .ne. j)then
        if((itypep(j) .ne. type(0)))then
          ! interparticle_distance
          rr2 = 0.0
          do m=0,idim-1
            drr(m) = x(j,m)-x(i,m)
            rr2   = rr2+drr(m)**2
          end do

          rr = sqrt(rr2)

          ! neighboring_particle_count and numbering
          if(rr .le. ker_r(2))then
            neigh(i,0)= neigh(i,0)+1
            neigh(i,neigh(i,0)) = j
          end if
        end if
      end if
    end do
  end if

  ! neighboring_particle_number_check ******************************
  if(neigh(i,0) .ge. neighbor)then
    write(900,*) ' neighboring particle number is too much'
    stop
  end if

  return
end subroutine set_neigh


!*********************************************************************
!*** array_allocation                                              ***
!*********************************************************************
subroutine allocation                          &
&          (idim,nn,neighbor,grid_num,nn_grid, &
&           visc,accel,press,cal_sw,           &
&           p_out,grid_p,out,neigh,solid_sur,  &
&           g_neigh,rho,nu,sol_vel,ker_r)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in)               :: idim,nn,neighbor
  integer, intent(in)               :: grid_num,nn_grid
  real   , intent(out), allocatable :: visc(:,:)
  real   , intent(out), allocatable :: accel(:,:)
  real   , intent(out), allocatable :: press(:,:)
  integer, intent(out), allocatable :: cal_sw(:,:)
  real   , intent(out), allocatable :: p_out(:)
  integer, intent(out), allocatable :: grid_p(:)
  real   , intent(out), allocatable :: out(:,:)
  integer, intent(out), allocatable :: neigh(:,:)
  integer, intent(out), allocatable :: g_neigh(:,:)
  real   , intent(out), allocatable :: rho(:)
  real   , intent(out), allocatable :: nu(:)
  real   , intent(out), allocatable :: sol_vel(:)
  real   , intent(out), allocatable :: ker_r(:)
  integer, intent(out), allocatable :: solid_sur(:)

  ! array_allocation *************************************************
  allocate(visc(0:nn-1,0:idim-1))
  allocate(accel(0:nn-1,0:idim-1))
  allocate(press(0:nn-1,0:idim-1))
  allocate(cal_sw(0:nn-1,0:2))
  allocate(p_out(0:nn-1))
  allocate(grid_p(0:nn-1))
  allocate(out(1:16,0:nn-1))
  allocate(neigh(0:nn-1,0:neighbor))
  allocate(g_neigh(1:grid_num,0:nn_grid))
  allocate(rho(0:nn-1))
  allocate(nu(0:nn-1))
  allocate(sol_vel(0:idim-1))
  allocate(ker_r(0:3))
  allocate(solid_sur(0:nn-1))

  ! array_initialize *************************************************
  visc(0:nn-1,0:idim-1)             = 0.0
  accel(0:nn-1,0:idim-1)            = 0.0
  press(0:nn-1,0:idim-1)            = 0.0
  p_out(0:nn-1)                     = 0.0
  out(1:16,0:nn-1)                  = 0.0
  rho(0:nn-1)                       = 0.0
  nu(0:nn-1)                        = 0.0
  cal_sw(0:nn-1,0:2)                = 0
  grid_p(0:nn-1)                    = 1!0
  neigh(0:nn-1,0:neighbor)          = 0
  g_neigh(1:grid_num,0:nn_grid)     = 0
  solid_sur(0:nn-1)                 = 0

  return
end subroutine allocation

!*********************************************************************
!*** calc_ww                                                       ***
!*********************************************************************
subroutine calc_ww(number,dist,ww,ww_num,ww_wall,dis)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: number
  integer, intent(in) :: ww_num
  real   , intent(in) :: ww_wall(0:7,0:ww_num)
  real   , intent(in) :: dist
  real   , intent(in) :: dis
  real   , intent(inout) :: ww

  ! subroutine_variable  *********************************************
  integer :: n
  real    :: x1,x2
  real    :: y1,y2
  real    :: dist_c

  ! correction_distance **********************************************
  dist_c = dist

  ! calc_ww **********************************************************
  if(dist_c .le. ww_wall(0,0))then
    ww = ww_wall(number,0)
    return
  end if

  if(ww_wall(0,ww_num) .lt. dist_c)then
    ww = 0.0
    return
  end if

  do n=1,ww_num
    if(ww_wall(0,n-1) .lt. dist_c)then
      if(dist_c .le. ww_wall(0,n))then
        y1 = ww_wall(number,n-1)
        y2 = ww_wall(number,n)
        x1 = ww_wall(0,n-1)
        x2 = ww_wall(0,n)

        ww = y1+(dist_c-x1)/(x2-x1)*(y2-y1)
        return
      end if
    end if
  end do

  return
end subroutine calc_ww

!*********************************************************************
!*** input_initial_condition                                       ***
!*********************************************************************
subroutine input_file                  &
&    (nn,idim,ifluid,neighbor,type,    &
&     den,itypep,x,v,an,grav,vis,      &
&     time_max,ite_max,dt_fix,cfl,     &
&     dtout,dis,ker_c,            &
&     delta,sigma,n0,         &
&     ite,time_sim,nump,rlambda,       &
&     grid_n,grid_num,grid,near_wall,  &
&     non_cal,dis_rat,c_vel,           &
&     f_type,num_ref,p,           &
&     nn_grid,output_limit,t,          &
&     dist_min,coll_rat,surf,          &
&     wall_n,grid_dist,wall_anc,       &
&     wall_nor)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(out)              :: neighbor
  integer, intent(out)              :: nn_grid
  integer, intent(out)              :: output_limit
  integer, intent(out)              :: idim
  integer, intent(out)              :: nn,ifluid
  integer, intent(out)              :: grid_num
  real   , intent(out)              :: grav
  real   , intent(out)              :: time_max
  integer, intent(out)              :: ite_max
  real   , intent(out)              :: dt_fix,cfl
  real   , intent(out)              :: dtout
  real   , intent(out)              :: dis
  real   , intent(out)              :: delta,sigma
  integer, intent(out)              :: ite
  real   , intent(out)              :: time_sim
  integer, intent(out)              :: nump
  real   , intent(out)              :: dis_rat
  real   , intent(out)              :: c_vel
  real   , intent(out)              :: num_ref
  real   , intent(out)              :: rlambda
  real   , intent(out)              :: dist_min
  real   , intent(out)              :: coll_rat
  integer, intent(out)              :: wall_n
  real   , intent(out), allocatable :: x(:,:)
  real   , intent(out), allocatable :: v(:,:)
  real   , intent(out), allocatable :: an(:)
  real   , intent(out), allocatable :: n0(:,:)
  integer, intent(out), allocatable :: itypep(:)
  real   , intent(out), allocatable :: den(:)
  real   , intent(out), allocatable :: vis(:)
  integer, intent(out), allocatable :: type(:)
  integer, intent(out), allocatable :: grid_n(:)
  real   , intent(out), allocatable :: non_cal(:)
  integer, intent(out), allocatable :: f_type(:)
  real   , intent(out), allocatable :: ker_c(:)
  real   , intent(out), allocatable :: grid(:,:,:)
  real   , intent(out), allocatable :: t(:)
  real   , intent(out), allocatable :: surf(:,:)
  real   , intent(out), allocatable :: p(:)
  real   , intent(out), allocatable :: wall_anc(:,:,:)
  real   , intent(out), allocatable :: grid_dist(:,:)
  integer, intent(out), allocatable :: near_wall(:,:)
  real   , intent(out), allocatable :: wall_nor(:,:)

  ! subroutine_varible ***********************************************
  integer              :: i,j,n
  real                 :: dummy_real
  integer, allocatable :: check_pos(:)
  character(len = 40)  :: dummy

  ! input_dirctry_setting ********************************************
  open(11,file='./data/MPSini/initial_setting.dat')
  open(12,file='./data/MPSini/mps000.vtk')

  ! input_initial_condition ******************************************
  read(11,*)idim, nn, nump, ifluid
  read(11,*)neighbor, nn_grid, output_limit
  read(12,*)dummy
  read(12,*)dummy
  read(12,*)dummy
  read(12,*)dummy
  read(12,*)dummy,nump,dummy

  ! array_allocation
  allocate(x(0:nn-1,0:idim-1))
  allocate(v(0:nn-1,0:idim-1))
  allocate(an(0:nn-1))
  allocate(n0(0:nn-1,0:3))
  allocate(itypep(0:nn-1))
  allocate(den(0:ifluid))
  allocate(vis(0:ifluid))
  allocate(type(0:ifluid))
  allocate(ker_c(0:3))
  allocate(grid_n(0:idim-1))
  allocate(non_cal(0:idim-1))
  allocate(f_type(1:3))
  allocate(t(0:nn-1))
  allocate(surf(0:nn-1,0:idim-1))
  allocate(p(0:nn-1))

  ! array_initialize
  x(0:nn-1,0:idim-1)    = 0.0
  v(0:nn-1,0:idim-1)    = 0.0
  an(0:nn-1)            = 0.0
  n0(0:nn-1,0:3)        = 0.0
  den(0:ifluid)         = 0.0
  vis(0:ifluid)         = 0.0
  ker_c(0:3)            = 0.0
  non_cal(0:idim-1)     = 0.0
  t(0:nn-1)             = 0.0
  p(0:nn-1)             = 0.0
  surf(0:nn-1,0:idim-1) = 0.0
  itypep(0:nn-1)        = 0
  type(0:ifluid)        = 0
  f_type(1:3)           = 0
  grid_n(0:idim-1)      = 0

  do n=0,nump-1
    read(12,*)x(n,0:idim-1)
  end do

  read(12,*)dummy
  read(12,*)dummy
  read(12,*)dummy

  do n=0,nump-1
    read(12,*)dummy_real
    itypep(n) = int(dummy_real)
  end do

  read(12,*)dummy

  do n=0,nump-1
    read(12,*)v(n,0:idim-1)
  end do

  read(12,*)dummy

  do n=0,nump-1
    read(12,*)p(n)
  end do

  read(12,*)dummy

  do n=0,nump-1
    read(12,*)p(n)
  end do

  read(12,*)dummy

  do n=0,nump-1
    read(12,*)surf(n,0:idim-1)
  end do

  read(12,*)dummy

  do n=0,nump-1
    read(12,*)t(n)
  end do

  do i=0,ifluid
    read(11,*)type(i), den(i), vis(i)
  end do

  read(11,*)f_type(1), f_type(2), f_type(3)

  read(11,*)ite_max, time_max
  read(11,*)ite,     time_sim
  read(11,*)cfl,     dt_fix
  read(11,*)dtout
  read(11,*)ker_c(0), ker_c(1), ker_c(2), ker_c(3)
  read(11,*)n0(0,0),  n0(0,1),  n0(0,2),  n0(0,3)
  read(11,*)num_ref
  read(11,*)rlambda
  read(11,*)dis_rat
  read(11,*)c_vel
  read(11,*)dist_min,coll_rat

  do i=0,idim-1
    read(11,*)non_cal(i)
  end do

  do i=0,nn-1
    do j=0,3
      n0(i,j) = n0(0,j)
    end do
  end do

  read(11,*)grav
  read(11,*)delta, sigma
  read(11,*)dis

  do i=0,idim-1
    read(11,*)grid_n(i)
  end do

  read(11,*)wall_n
  read(11,*)grid_num

  allocate(grid(1:grid_num,0:idim-1,0:1))
  allocate(grid_dist(1:grid_num,0:4*(idim-1)-1))
  allocate(near_wall(1:grid_num,0:4*(idim-1)-1))
  allocate(wall_anc(0:idim-1,0:wall_n-1,0:idim-1))
  allocate(wall_nor(0:idim-1,0:wall_n-1))

  read(11,*)grid(1:grid_num,0:idim-1,0:1)
  read(11,*)wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  read(11,*)wall_nor(0:idim-1,0:wall_n-1)

  do i=1,grid_num
    do j=0,4*(idim-1)-1
      read(11,*)grid_dist(i,j),near_wall(i,j)
    end do
  end do

  close(11)
  close(12)

   allocate(check_pos(0:idim-1))
  do i=0,nump-1
    if(itypep(i) .ge. type(1))then
      do j=0,nump-1
        if(i .ne. j)then
          check_pos(0:idim-1) = 0

          do n=0,idim-1
            if(x(i,n) .eq. x(j,n))then
              check_pos(n) = 1
            end if
          end do

          do n=0,idim-1
            if(minval(check_pos(0:idim-1)) .ne. 0)then
              write(900,*)'alart in input_file(mod_mps)'
              write(900,*)'input data have same position'
            end if
          end do
        end if
      end do
    end if
  end do

  ! array_deallocation ***********************************************
  deallocate(check_pos)

  return
end subroutine input_file

!*********************************************************************
!*** renewal_initial_file                                          ***
!*********************************************************************
subroutine output_file                   &
&          (idim,nn,nump,ite,time_sim,   &
&           type,ifluid,itypep,          &
&           den,vis,grav,delta,sigma,    &
&           time_max,ite_max,            &
&           dt,dtout,cfl,       &
&           dis,ker_c,n0,           &
&           x,non_cal,num_ref,rlambda,   &
&           grid,grid_n,grid_num,        &
&           dis_rat,c_vel,          &
&           f_type,output_limit,         &
&           nn_grid,neighbor,t,wall_nor, &
&           dist_min,coll_rat,near_wall, &
&           wall_n,grid_dist,wall_anc)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: nn
  integer, intent(in) :: idim
  integer, intent(in) :: nump
  integer, intent(in) :: ifluid
  integer, intent(in) :: ite
  integer, intent(in) :: type(0:ifluid)
  real   , intent(in) :: time_sim
  real   , intent(in) :: den(0:ifluid)
  real   , intent(in) :: vis(0:ifluid)
  real   , intent(in) :: grav
  real   , intent(in) :: time_max
  real   , intent(in) :: dt
  real   , intent(in) :: cfl
  integer, intent(in) :: ite_max
  real   , intent(in) :: dtout
  real   , intent(in) :: dis
  real   , intent(in) :: delta
  real   , intent(in) :: sigma
  integer, intent(in) :: itypep(0:nn-1)
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  integer, intent(in) :: grid_num
  integer, intent(in) :: grid_n(0:idim-1)
  real   , intent(in) :: grid(1:grid_num,0:idim-1,0:1)
  real   , intent(in) :: non_cal(0:idim-1)
  real   , intent(in) :: dis_rat
  real   , intent(in) :: c_vel
  integer, intent(in) :: f_type(1:3)
  real   , intent(in) :: ker_c(0:3)
  integer, intent(in) :: nn_grid
  integer, intent(in) :: neighbor
  integer, intent(in) :: output_limit
  real   , intent(in) :: dist_min
  real   , intent(in) :: coll_rat
  real   , intent(in) :: t(0:nn-1)
  real   , intent(in) :: n0(0:nn-1,0:3)
  real   , intent(in) :: num_ref
  real   , intent(in) :: rlambda
  integer, intent(in) :: wall_n
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: grid_dist(1:grid_num,0:4*(idim-1)-1)
  integer, intent(in) :: near_wall(1:grid_num,0:4*(idim-1)-1)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)

  ! subroutine_variable **********************************************
  integer             :: n,i,j

  ! output_file_name *************************************************
  open(11,file='./data/MPS/initial_setting.dat')

  ! renewal_initial_file *********************************************
  write(11,*)idim, nn, nump, ifluid
  write(11,*)neighbor, nn_grid, output_limit

  do i=0,ifluid
    write(11,*)type(i), den(i), vis(i)
  end do

  write(11,*)f_type(1), f_type(2), f_type(3)

  write(11,*)ite_max, time_max
  write(11,*)ite,     time_sim
  write(11,*)cfl,     dt
  write(11,*)dtout
  write(11,*)ker_c(0), ker_c(1), ker_c(2), ker_c(3)
  write(11,*)n0(0,0),  n0(0,1),  n0(0,2),  n0(0,3)
  write(11,*)num_ref
  write(11,*)rlambda
  write(11,*)dis_rat
  write(11,*)c_vel
  write(11,*)dist_min,coll_rat

  do n=0,idim-1
    write(11,*)non_cal(n)
  end do

  write(11,*)grav
  write(11,*)delta, sigma
  write(11,*)dis

  do n=0,idim-1
    write(11,*)grid_n(n)
  end do

  write(11,*)wall_n
  write(11,*)grid_num
  write(11,*)grid(1:grid_num,0:idim-1,0:1)
  write(11,*)wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  write(11,*)wall_nor(0:idim-1,0:wall_n-1)

  do i=1,grid_num
    do j=0,4*(idim-1)-1
      write(11,*)grid_dist(i,j),near_wall(i,j)
    end do
  end do

  close(11)

  return
end subroutine output_file

!*********************************************************************
!*** calculation_temperature                                       ***
!*********************************************************************
subroutine cal_temperature               &
&          (nn,idim,ifluid,neighbor,     &
&           nump,itypep,type,cal_sw,den, &
&           neigh,x,ker_r,rlambda,n0,t,dt)
  implicit none
  ! mainroutine_variable *********************************************
  integer ,intent(in) :: nn,idim,ifluid,neighbor
  integer ,intent(in) :: nump
  integer ,intent(in) :: itypep(0:nn-1)
  integer ,intent(in) :: type(0:ifluid)
  integer ,intent(in) :: neigh(0:nn-1,0:neighbor)
  real    ,intent(in) :: x(0:nn-1,0:idim-1)
  real    ,intent(in) :: rlambda
  real    ,intent(in) :: n0(0:nn-1,0:3)
  integer ,intent(in) :: cal_sw(0:nn-1,0:2)
  real    ,intent(in) :: ker_r(0:3)
  real   , intent(inout) :: t(0:nn-1)
  real   , intent(in) :: den(0:ifluid)
  real   , intent(in) :: dt

  ! subroutine_variable **********************************************
  integer :: i,j,l,m
  real    :: rr,rr2
  real    :: ww
  real    :: dter
  real    :: r_kernel
  real    :: drr(0:idim-1)
  real, allocatable :: t_temp(:)
  real, parameter   :: zero = 1.0e-10
  real, parameter   :: c_v  = 4.2*1.0e3
  real, parameter   :: lam  = 0.552

  ! array_allocation *************************************************
  allocate(t_temp(0:nump-1))

  ! initial_setting **************************************************
  r_kernel = ker_r(1)
  t_temp(0:nump-1) = t(0:nump-1)

  ! calculation_temperature ******************************************
  do i=0,nump-1
    if((itypep(i) .ge. type(4)).and. &
    &  (itypep(i) .le. type(4)))then
      dter = 0.0
      do l=1,neigh(i,0)
        j = neigh(i,l)
        if((itypep(j) .ge. type(2)).and.  &
        &  (itypep(j) .le. type(4)))then
          ! interparticle_distance
          rr2 = 0.0
          do m=0,idim-1
            drr(m) = x(j,m)-x(i,m)
            rr2   = rr2+drr(m)**2
          end do

          rr = sqrt(rr2)

          ! weight_function
          ww = r_kernel/max(zero,rr)+rr/max(zero,r_kernel)-2.0
          if(rr .gt. r_kernel)ww = 0.0

          ! calculation_temperature
          dter  = 2.0*real(idim)*ww/rlambda/n0(i,1)
          dter  = lam/(den(itypep(i))*c_v)*dter*(t_temp(j)-t_temp(i))
          t(i) = t(i)+dter*dt
        end if
      end do
    end if
  end do

  ! array_deallocation ***********************************************
  deallocate(t_temp)

  return
end subroutine cal_temperature

!*********************************************************************
!*** calculation_temperature(without_wall)                         ***
!*********************************************************************
subroutine cal_temperature_wall_out      &
&          (nn,idim,ifluid,neighbor,     &
&           nump,itypep,type,cal_sw,den, &
&           neigh,x,ker_r,rlambda,n0,t,dt)
  implicit none
  ! maintourine_variable *********************************************
  integer ,intent(in) :: nn,idim,ifluid,neighbor
  integer ,intent(in) :: nump
  integer ,intent(in) :: itypep(0:nn-1)
  integer ,intent(in) :: type(0:ifluid)
  integer ,intent(in) :: neigh(0:nn-1,0:neighbor)
  real    ,intent(in) :: x(0:nn-1,0:idim-1)
  real    ,intent(in) :: rlambda
  real    ,intent(in) :: n0(0:nn-1,0:3)
  integer ,intent(in) :: cal_sw(0:nn-1,0:2)
  real    ,intent(in) :: ker_r(0:3)
  real   , intent(inout) :: t(0:nn-1)
  real   , intent(in) :: den(0:ifluid)
  real   , intent(in) :: dt

  ! subroutine_variable **********************************************
  integer :: i,j,l,m
  real    :: rr,rr2
  real    :: ww
  real    :: dter
  real    :: r_kernel
  real    :: drr(0:idim-1)
  real, allocatable :: t_temp(:)
  real, allocatable :: lam(:)
  real, allocatable :: c_v(:)
  real, parameter   :: zero = 1.0e-10

  ! array_allocation *************************************************
  allocate(t_temp(0:nump-1))
  allocate(lam(0:ifluid))
  allocate(c_v(0:ifluid))

  ! initial_setting **************************************************
  r_kernel = ker_r(1)
  t_temp(0:nump-1) = t(0:nump-1)

  ! parameter_setting ************************************************
  c_v(type(4)) = 4217.0
  c_v(type(5)) = 1005.0

  lam(type(4)) = 0.569
  lam(type(5)) = 0.0241

  ! calculation_temperature ******************************************
  do i=0,nump-1
    if((itypep(i) .ge. type(4)).and. &
    &  (itypep(i) .le. type(5)))then
      dter = 0.0
      do l=1,neigh(i,0)
        j = neigh(i,l)
        if((itypep(j) .ge. type(3)).and.  &
        &  (itypep(j) .le. type(5)))then
          ! interparticle_distance
          rr2 = 0.0
          do m=0,idim-1
            drr(m) = x(j,m)-x(i,m)
            rr2   = rr2+drr(m)**2
          end do

          rr = sqrt(rr2)

          ! weight_function
          ww = r_kernel/max(zero,rr)+rr/max(zero,r_kernel)-2.0
          if(rr .gt. r_kernel)ww = 0.0

          ! calculation_temperature
          dter  = 2.0*real(idim)*ww/rlambda/n0(i,1)
          dter  = lam(itypep(i))/c_v(itypep(i))*dter*(t_temp(j)-t_temp(i))
          t(i) = t(i)+dter*dt
        end if
      end do
    end if
  end do

  ! array_deallocation ***********************************************
  deallocate(t_temp)

  return
end subroutine cal_temperature_wall_out

!*********************************************************************
!*** calculation_temperature(with_icing)                           ***
!*********************************************************************
subroutine cal_temperature_icing         &
&          (nn,idim,ifluid,neighbor,     &
&           nump,itypep,type,cal_sw,den, &
&           neigh,x,ker_r,rlambda,n0,t)
  implicit none
  ! mainroutine_variable *********************************************
  integer ,intent(in) :: nn,idim,ifluid,neighbor
  integer ,intent(in) :: nump
  integer ,intent(in) :: itypep(0:nn-1)
  integer ,intent(in) :: type(0:ifluid)
  integer ,intent(in) :: neigh(0:nn-1,0:neighbor)
  real    ,intent(in) :: x(0:nn-1,0:idim-1)
  real    ,intent(in) :: rlambda
  real    ,intent(in) :: n0(0:nn-1,0:3)
  integer ,intent(in) :: cal_sw(0:nn-1,0:2)
  real    ,intent(in) :: ker_r(0:3)
  real   , intent(inout) :: t(0:nn-1)
  real   , intent(in) :: den(0:ifluid)

  ! subroutine_variable **********************************************
  integer :: i,j,l,m
  real    :: rr,rr2
  real    :: ww
  real    :: dter
  real    :: r_kernel
  real    :: drr(0:idim-1)
  real, allocatable :: t_temp(:)
  real, parameter   :: zero = 1.0e-10
  real, parameter   :: c_v  = 4.2*1.0e3
  real, parameter   :: lam  = 0.552

  ! array_allocation *************************************************
  allocate(t_temp(0:nump-1))

  ! initial_setting **************************************************
  r_kernel = ker_r(1)
  t_temp(0:nump-1) = t(0:nump-1)

  ! calculation_temperature ******************************************
  do i=0,nump-1
    if((itypep(i) .ge. type(3)).and. &
    &  (itypep(i) .le. type(4)))then
      dter = 0.0
      do l=1,neigh(i,0)
        j = neigh(i,l)
        if((itypep(j) .ge. type(2)).and.  &
        &  (itypep(j) .le. type(4)))then
          ! interparticle_distance
          rr2 = 0.0
          do m=0,idim-1
            drr(m) = x(j,m)-x(i,m)
            rr2   = rr2+drr(m)**2
          end do

          rr = sqrt(rr2)

          ! weight_function
          ww = r_kernel/max(zero,rr)+rr/max(zero,r_kernel)-2.0
          if(rr .gt. r_kernel)ww = 0.0

          ! calculation_temperature
          dter  = 2.0*real(idim)*ww/rlambda/n0(i,1)
          dter  = lam/(den(itypep(i))*c_v)*dter*(t_temp(j)-t_temp(i))
          t(i) = t(i)+dter
        end if
      end do

      if(itypep(i) .eq. type(3))t(i) = max(0.0,t(i))
    end if
  end do

  ! array_deallocation ***********************************************
  deallocate(t_temp)

  return
end subroutine cal_temperature_icing

!*********************************************************************
!*** calculation_temperature(with_supercooling)                    ***
!*********************************************************************
subroutine cal_temperature_supercooling  &
&          (nn,idim,ifluid,neighbor,     &
&           nump,itypep,type,cal_sw,den, &
&           neigh,x,ker_r,rlambda,n0,t)
  implicit none
  ! mainroutine_variable *********************************************
  integer ,intent(in) :: nn,idim,ifluid,neighbor
  integer ,intent(in) :: nump
  integer ,intent(in) :: itypep(0:nn-1)
  integer ,intent(in) :: type(0:ifluid)
  integer ,intent(in) :: neigh(0:nn-1,0:neighbor)
  real    ,intent(in) :: x(0:nn-1,0:idim-1)
  real    ,intent(in) :: rlambda
  real    ,intent(in) :: n0(0:nn-1,0:3)
  integer ,intent(in) :: cal_sw(0:nn-1,0:2)
  real    ,intent(in) :: ker_r(0:3)
  real   , intent(inout) :: t(0:nn-1)
  real   , intent(in) :: den(0:ifluid)

  ! subroutine_variable **********************************************
  integer :: i,j,l,m
  real    :: rr,rr2
  real    :: ww
  real    :: dter
  real    :: r_kernel
  real    :: drr(0:idim-1)
  real, allocatable :: t_temp(:)
  real, parameter   :: zero = 1.0e-10
  real, parameter   :: c_v  = 4.2*1.0e3
  real, parameter   :: lam  = 0.552

  ! array_allocation *************************************************
  allocate(t_temp(0:nump-1))

  ! initial_setting **************************************************
  r_kernel = ker_r(1)
  t_temp(0:nump-1) = t(0:nump-1)

  ! calculation_temperature ******************************************
  do i=0,nump-1
    if((itypep(i) .eq. type(4)))then
      dter = 0.0
      do l=1,neigh(i,0)
        j = neigh(i,l)
        if((itypep(j) .ge. type(3)).and.  &
        &  (itypep(j) .le. type(4)))then
          ! interparticle_distance
          rr2 = 0.0
          do m=0,idim-1
            drr(m) = x(j,m)-x(i,m)
            rr2   = rr2+drr(m)**2
          end do

          rr = sqrt(rr2)

          ! weight_function
          ww = r_kernel/max(zero,rr)+rr/max(zero,r_kernel)-2.0
          if(rr .gt. r_kernel)ww = 0.0

          ! calculation_temperature
          dter  = 2.0*real(idim)*ww/rlambda/n0(i,1)
          dter  = lam/(den(itypep(i))*c_v)*dter*(t_temp(j)-t_temp(i))
          t(i) = t(i)+dter
        end if
      end do
    end if

    if(itypep(i) .eq. type(3))t(i) = max(0.0,t(i))
  end do

  ! array_deallocation ***********************************************
  deallocate(t_temp)

  return
end subroutine cal_temperature_supercooling
end module mod_mps
