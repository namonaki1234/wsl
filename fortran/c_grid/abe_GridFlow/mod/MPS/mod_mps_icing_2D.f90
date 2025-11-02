module mod_icing_dim
  use mod_mps_dim
  use mod_mps
  implicit none
  ! subroutine_definition ********************************************

  ! function_definition ***********************************************

contains

! ********************************************************************
! ***********************  二次元対応可  *****************************
! ********************************************************************

!*********************************************************************
!*** count_icing particle_number_density                           ***
!*********************************************************************
subroutine cal_n_ice                      &
&           (i,ineigh,ker_r,an,nump,        &
&            nn,idim,ifluid,neighbor, &
&            type,grid_num,ww_num,         &
&            ww_wall,wall_n,grid_p,   &
&            neigh,x,itypep,       &
&            wall_anc,dis,            &
&            nn_wall,g_wall,wall_nor)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: ineigh
  integer, intent(in) :: i
  integer, intent(in) :: nn,idim,ifluid,neighbor
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  real   , intent(in) :: dis
  integer, intent(in)    :: nump
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  integer, intent(in) :: itypep(0:nn-1)
  integer, intent(in) :: type(0:ifluid)
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

  ! とりあえずカーネル半径を以前のに合わせて手打ち
  r_kernel = 2.1 * dis
  an(i) = 0.0

  ! count_particle_number_density **************************************
  do l=1,neigh(i,0)
    j = neigh(i,l)
    k = itypep(j)
    if(itypep(j) .eq. type(2))then
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
end subroutine cal_n_ice

!*********************************************************************
!*** icing                           ***
!*********************************************************************
subroutine icing                      &
&           (i,nn,idim,ifluid,nump,x,v, &
&            itypep,type,dis,ice_an,icing_col, &
&            ice_dist_min,ice_an_min,changep)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i,nn,idim,ifluid
  real   , intent(in) :: dis
  integer, intent(in)    :: nump
  real   , intent(inout) :: x(0:nn-1,0:idim-1)
  real   , intent(inout) :: v(0:nn-1,0:idim-1)
  integer, intent(in)    :: type(0:ifluid)
  integer, intent(inout) :: itypep(0:nn-1)
  real   , intent(in) :: ice_an(0:nn-1)
  integer, intent(in) :: icing_col(0:nn-1)
  real   , intent(in) :: ice_dist_min,ice_an_min
  integer, intent(inout) :: changep

  ! subroutine_variable ***********************************************
  integer :: l,m,j,k
  real   , parameter :: pi   = acos(-1.0)
  real   , parameter :: zero = 1.0e-10
  ! array_allocation ***************************************************

  ! initializing *******************************************************
  ! count_particle_number_density **************************************
  if((icing_col(i) .eq. 1)       .and. &
  &  (ice_an(i)    .ge. ice_an_min))then
    v(i,0) = 0.0
    v(i,1) = 0.0
    itypep(i) = type(2)
    changep = 1
  end if

  ! array_deallocation ***********************************************

  return
end subroutine icing

!*********************************************************************
!*** calculation_viscous_term (書き直した)                         ***
!*********************************************************************
subroutine temperature                       &
&          (i,nn,idim,ifluid,neighbor,wall_num, &
&           nump,visc,itypep,type,cal_sw,     &
&           neigh,x,v,nu,ker_r,rlambda,n0,    &
&           grid_num,ww_num,grid,ww_wall,     &
&           wall_n,grid_p,grid_dist,dis,      &
&           near_wall,nn_wall,g_wall,wall_anc,&
&           wall_nor,solid_sur,den,t,dt,      &
&           t_wall,t_temp,ganma,ganma_temp,swi_koeki,Qhc, &
&           swi_solid,flag_temp,wall,wall_dist)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: i,nn,idim,ifluid,neighbor
  integer, intent(in) :: wall_num
  integer, intent(in) :: nump
  real   , intent(in) :: nu(0:nn-1)
  integer, intent(inout) :: itypep(0:nn-1)
  integer, intent(in) :: type(0:ifluid)
  integer, intent(in) :: neigh(0:nn-1,0:neighbor)
  real   , intent(in)  :: x(0:nn-1,0:idim-1)
  real   , intent(out) :: v(0:nn-1,0:idim-1)
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
  real   , intent(in) :: den(0:ifluid)
  real   , intent(inout) :: visc(0:nn-1,0:idim-1)
  real   , intent(inout) :: t(0:nn-1)
  real   , intent(in) :: dt
  real   , intent(in) :: t_wall
  real   , intent(in) :: t_temp(0:nn-1)
  real   , intent(inout) :: ganma(0:nn-1)
  real   , intent(in) :: ganma_temp(0:nn-1)
  integer, intent(inout) :: swi_koeki(0:nn-1)
  real   , intent(in) :: Qhc(0:nn-1)
  integer, intent(in) :: swi_solid(0:nn-1)
  integer, intent(in) :: flag_temp(0:nn-1)
  real   , intent(in) :: wall(0:130,0:1)
  real   , intent(inout) :: wall_dist(0:nn-1)

  ! subroutine_variable **********************************************
  integer :: j,l,m,k,i1
  real    :: rr,rr2,rr0
  real    :: ww
  real    :: dvis0
  real    :: r_kernel
  integer :: l_ref
  real    :: ww_ref
!  real    :: lam_ij !lam_ijとden_ijを使った式は間違い
  real    :: delta_h
  real    :: h_ls
  real    :: t_melt
!  real    :: den_ij
  real    :: lam_trans,c_p_trans,den_trans
  real    :: dvis(0:idim-1)
  real    :: drr0(0:idim-1)
  real    :: drr(0:idim-1)
  real    :: x0(0:idim-1)
  real, parameter   :: zero = 1.0e-10
!  real, parameter   :: c_p  = 4.2*1.0e3
!  real, parameter   :: lam  = 0.552
  real    :: dter
  real    :: c_p(0:ifluid)
  real    :: lam(0:ifluid)
  real, allocatable :: t_s(:)
  real, allocatable :: t_s_wall(:)
  real              :: den_wall
  real              :: c_p_wall
  real              :: lam_wall

  ! array_allocation *************************************************
  allocate(t_s(0:nn-1))
  allocate(t_s_wall(0:nn-1))

  ! lam : 熱伝導率 [w/(m*k)]
  ! c_p : 定圧比熱 [j/(kg*k)]
    ! 液体
  c_p(type(4)) = 4.2*1.0e3
  lam(type(4)) = 0.552
  ! 気体
  c_p(type(3)) = 519.0
  lam(type(3)) = 17.0
  ! 壁面
  c_p(type(2)) = 2.100*1.0e3  !     !
  lam(type(2)) = 2.2  !14.9      !62.2
  !プレート
den_wall = 2790.0
  c_p_wall = 840.0!477.0
  lam_wall = 170.0 !14.9
  !den_wall = 2329.0!7900.0
  !c_p_wall = 740.0!519.0!477.0
  !lam_wall =156.0!17.0 !14.9
  ! 潜熱
  h_ls =3.344*1.0e5  !Sn(凝固潜熱)!58500.0

  ! 融点温度
  t_melt = 273.15

  ! initital_setting *************************************************
  r_kernel = ker_r(1)

  ! array_initialize *************************************************
  delta_h = 0.0

  ! calculation_heat_transfer*****************************************
  if(Qhc(i) .ne. 0.0) then
   t(i) = t(i) + Qhc(i)/(den(itypep(i))*c_p(itypep(i))) * dt
   delta_h = Qhc(i)/den(itypep(i)) * dt
   if((itypep(i) .eq. type(4)).and.(swi_koeki(i) .eq. 0)) then
    t(i) = t_melt
    ganma(i) = ganma(i)+(delta_h/h_ls)
    if (ganma(i) .le. 0.0)then
     itypep(i) = type(2)
     if(swi_solid(i) .eq. 1) then
      do m = 0,idim-1
       v(i,m) = 0.0
      end do
     end if
    end if
   end if
  end if

  if(t(i) .ne. t(i)) then
   write(*,*) 'Error : temperature is diverged'
   stop
  end if

  ! calculation_viscous_term *****************************************
  do l=1,neigh(i,0)
    j = neigh(i,l)
    if(((itypep(j) .ge. type(2)).and.  &
    &  (itypep(j) .le. type(4)))) then

      ! interparticle_distance

      dter = 0.0

      rr2 = 0.0
      do m=0,idim-1
        drr(m) = x(j,m)-x(i,m)
        rr2   = rr2+drr(m)**2
      end do

      rr = sqrt(rr2)

      ! weight_function
      ww = r_kernel/max(zero,rr)+rr/max(zero,r_kernel)-2.0
      if(rr .gt. r_kernel)ww = 0.0

     !jの粒子が内部粒子の場合はノイマン境界条件を適用
     if(flag_temp(j) .eq. 0) then
      if(wall_dist(i) .eq. -1.0) then
       call wall_dist_t(x(i,0),x(i,1),wall,wall_dist(i))
      end if
      if(wall_dist(j) .eq. -1.0) then
       call wall_dist_t(x(j,0),x(j,1),wall,wall_dist(j))
      end if
      !write(*,*) t(j),wall_dist(j)/wall_dist(i),(t(i) - t_wall) / wall_dist(i) * wall_dist(j) + t_wall
      if(wall_dist(j)/wall_dist(i) .ge. 10.0) then
       write(*,*) 'error: wall distance',wall_dist(j)/wall_dist(i)
       write(*,*) wall_dist(j),wall_dist(i)
       write(*,*) x(j,0),x(j,1),x(i,0),x(i,1)
       stop
      end if
      t(j) = min((t(i) - t_wall) / wall_dist(i) * wall_dist(j) + t_wall,t_melt)
      if(itypep(i) .eq. type(4)) wall_dist(i) = -1.0
     end if

      !iが固液でない かつ jも固液でない
      if (( (itypep(i) .ne. type(4)) .or.                                &
       &   ((itypep(i) .eq. type(4)) .and. (swi_koeki(i) .ne. 0))) .and. &
       &  ( (itypep(j) .ne. type(4)) .or.                                &
       &   ((itypep(j) .eq. type(4)) .and. (swi_koeki(j) .ne. 0))) ) then

        if (itypep(i) .eq. itypep(j)) then !i,j 同一物性の場合の温度変化
          dter    = 2.0*real(idim)*ww/rlambda/n0(i,1)
          dter    = lam(itypep(i))/(den(itypep(i))*c_p(itypep(i)))*dter*(t_temp(j)-t_temp(i))  !&
!                & + Qhc(i)/(den(itypep(i))*c_p(itypep(i)))
          t(i)    = t(i)+dter*dt

          delta_h = 2.0*real(idim)*ww/rlambda/n0(i,1)
          delta_h = delta_h*lam(itypep(i))/den(itypep(i))*(t_temp(j)-t_temp(i))   !&
!                & + Qhc(i)/den(itypep(i))
          delta_h = delta_h*dt
        else  !i,j 異なる物性の場合の温度変化
          !接触点温度
          t_s(j) =  (    (sqrt(den(itypep(i))*c_p(itypep(i))*lam(itypep(i)))*t_temp(i))  &
                   &  +  (sqrt(den(itypep(j))*c_p(itypep(j))*lam(itypep(j)))*t_temp(j)) )&
                   & / ( (sqrt(den(itypep(i))*c_p(itypep(i))*lam(itypep(i))))            &
                   &  +  (sqrt(den(itypep(j))*c_p(itypep(j))*lam(itypep(j))))           )

          dter     = 2.0*real(idim)*ww/rlambda/n0(i,1)
          dter     = lam(itypep(i))/(den(itypep(i))*c_p(itypep(i)))*dter*0.5*(t_s(j)-t_temp(i))  !&
!                & + Qhc(i)/(den(itypep(i))*c_p(itypep(i)))
          t(i)     = t(i)+dter*dt

          delta_h  = 2.0*real(idim)*ww/rlambda/n0(i,1)
          delta_h  = delta_h*lam(itypep(i))/den(itypep(i))*0.5*(t_s(j)-t_temp(i))  !&
!                & + Qhc(i)/den(itypep(i))
          delta_h  = delta_h*dt

        end if

      !iが固液である かつ jは固液でない
      else if ( ((itypep(i) .eq. type(4)) .and. (swi_koeki(i) .eq. 0)) .and. &
       &        ((itypep(j) .ne. type(4)) .or.                               &
       &        ((itypep(j) .eq. type(4)) .and. (swi_koeki(j) .ne. 0))) ) then

        !体積分率から物性を決定
        lam_trans = lam(2)*lam(4)/(ganma_temp(i)*lam(4)+(1.0-ganma_temp(i))*lam(2))
        den_trans = ganma_temp(i)*den(4)+(1.0-ganma_temp(i))*den(2)
        c_p_trans = ganma_temp(i)*c_p(4)+(1.0-ganma_temp(i))*c_p(2)

        !接触点温度
        t_s(j) =  (    (sqrt(den_trans     *c_p_trans     *lam_trans     )*t_temp(i))  &
                 &  +  (sqrt(den(itypep(j))*c_p(itypep(j))*lam(itypep(j)))*t_temp(j)) )&
                 & / ( (sqrt(den_trans     *c_p_trans     *lam_trans)     )            &
                 &  +  (sqrt(den(itypep(j))*c_p(itypep(j))*lam(itypep(j))))           )

        dter    = 2.0*real(idim)*ww/rlambda/n0(i,1)
        dter    = lam_trans/(den_trans*c_p_trans)*dter*0.5*(t_s(j)-t_temp(i))  !&
!              & + Qhc(i)/(den_trans*c_p_trans)
        t(i)    = t(i)+dter*dt

        delta_h = 2.0*real(idim)*ww/rlambda/n0(i,1)
        delta_h = delta_h*lam_trans/den_trans*0.5*(t_s(j)-t_temp(i))  !&
!              & + Qhc(i)/den_trans
        delta_h = delta_h*dt

      !iは固液でない かつ jは固液である
      else if (( (itypep(i) .ne. type(4)) .or.                              &
       &        ((itypep(i) .eq. type(4)) .and. (swi_koeki(i) .ne. 0))) .and.  &
       &       ( (itypep(j) .eq. type(4)) .and. (swi_koeki(j) .eq. 0)) ) then

        !体積分率から物性を決定
        lam_trans = lam(2)*lam(4)/(ganma_temp(j)*lam(4)+(1.0-ganma_temp(j))*lam(2))
        den_trans = ganma_temp(j)*den(4)+(1.0-ganma_temp(j))*den(2)
        c_p_trans = ganma_temp(j)*c_p(4)+(1.0-ganma_temp(j))*c_p(2)

        !接触点温度
        t_s(j) =  (    (sqrt(den(itypep(i))*c_p(itypep(i))*lam(itypep(i)))*t_temp(i))  &
                 &  +  (sqrt(den_trans     *c_p_trans     *lam_trans     )*t_temp(j)) )&
                 & / ( (sqrt(den(itypep(i))*c_p(itypep(i))*lam(itypep(i))))            &
                 &  +  (sqrt(den_trans     *c_p_trans     *lam_trans     ))           )


        dter    = 2.0*real(idim)*ww/rlambda/n0(i,1)
        dter    = lam(itypep(i))/(den(itypep(i))*c_p(itypep(i)))*dter*0.5*(t_s(j)-t_temp(i))  !&
!              & + Qhc(i)/(den(itypep(i))*c_p(itypep(i)))
        t(i)    = t(i)+dter*dt

        delta_h = 2.0*real(idim)*ww/rlambda/n0(i,1)
        delta_h = delta_h*lam(itypep(i))/den(itypep(i))*0.5*(t_s(j)-t_temp(i))  !&
!              & + Qhc(i)/den(itypep(i))
        delta_h = delta_h*dt

      !iもjも固液だと温度が融点で等しく熱の移動はない
      end if


      if((itypep(i) .eq. type(4)).and.(swi_koeki(i) .eq. 0)) then

        t(i) = t_melt

        ganma(i) = ganma(i)+(delta_h/h_ls)

        if (ganma(i) .le. 0.0)then
          itypep(i) = type(2)
          if(swi_solid(i) .eq. 1) then
            do m = 0,idim-1
             v(i,m) = 0.0
            end do
          end if
        end if
      end if
    end if
  end do

  ! wall_influence ************************************************************
  if(g_wall(grid_p(i),0) .ne. 0)then
    rr = r_kernel

    dter = 0.0
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
    call calc_ww(3,rr,ww_ref,ww_num,ww_wall,dis)

    ! an_correction
    call correction_weight_function_for_viscosity    &
    &    (i,drr,rr,nn,idim,x,grid_num,wall_n,grid_p, &
    &     dis,ww_ref,nn_wall,g_wall,wall_anc,wall_nor)

    !iが固液でない場合
    if (( itypep(i) .ne. type(4)) .or.  &
     &  ((itypep(i) .eq. type(4)) .and. (swi_koeki(i) .ne. 0)) ) then

      !接触点温度
      t_s_wall(i) =  (    (sqrt(den(itypep(i))*c_p(itypep(i))*lam(itypep(i)))*t_temp(i))  &
                    &  +  (sqrt(den_wall      *c_p_wall      *lam_wall      )*t_wall   ) )&
                    & / ( (sqrt(den(itypep(i))*c_p(itypep(i))*lam(itypep(i))))            &
                    &  +  (sqrt(den_wall      *c_p_wall      *lam_wall      ))           )

    !iが固液の場合
    else

      !体積分率から物性を決定
      lam_trans = lam(2)*lam(4)/(ganma_temp(i)*lam(4)+(1.0-ganma_temp(i))*lam(2))
      den_trans = ganma_temp(i)*den(4)+(1.0-ganma_temp(i))*den(2)
      c_p_trans = ganma_temp(i)*c_p(4)+(1.0-ganma_temp(i))*c_p(2)

      !接触点温度
      t_s_wall(i) =   ( (sqrt(den_trans     *c_p_trans     *lam_trans )*t_temp(i)) &
                    &  +(sqrt(den_wall      *c_p_wall      *lam_wall  )*t_wall   ))&
                    & / (sqrt(den_trans    *c_p_trans     *lam_trans)           &
                    &  + sqrt(den_wall     *c_p_wall      *lam_wall)             )
    end if

    dter   = 2.0*real(idim)*ww_ref/rlambda/n0(i,1)
    dter   = lam(itypep(i))/(den(itypep(i))*c_p(itypep(i)))*dter*0.5*(t_s_wall(i)-t_temp(i))  !&
!                & + Qhc(i)/(den(itypep(i))*c_p(itypep(i)))
    t(i)   = t(i)+dter*dt

    delta_h  = 2.0*real(idim)*ww_ref/rlambda/n0(i,1)
    delta_h  = delta_h*lam(itypep(i))/den(itypep(i))*0.5*(t_s_wall(i)-t_temp(i))  !&
!                & + Qhc(i)/den(itypep(i))
    delta_h  = delta_h*dt


    if((itypep(i) .eq. type(4)).and.(swi_koeki(i) .eq. 0)) then

      t(i) = t_melt

      ganma(i) = ganma(i)+(delta_h/h_ls)

      if (ganma(i) .le. 0.0)then
        itypep(i) = type(2)
        if(swi_solid(i) .eq. 1) then
          do m = 0,idim-1
            v(i,m) = 0.0
          end do
        end if
      end if
    end if
  end if

  ! array_deallocation **********************************************
  deallocate(t_s)
  deallocate(t_s_wall)
  return
end subroutine temperature

subroutine wall_dist_t(x,y,wall,wall_dist)
 implicit none
 ! mainroutine_variable
 real,intent(in) :: x,y
 real,intent(in) :: wall(0:130,0:1)
 real,intent(inout) :: wall_dist
 ! subroutine_variable
 integer	:: i
 real		:: dist
 real		:: dist_min
 integer	:: min_i
 integer	:: neigh_i

 dist_min = 1.0e10
 do i = 0,130
  dist = sqrt((x - wall(i,0))**2.0 + (y - wall(i,1))**2.0)
  if(dist .lt. dist_min) then
   dist_min = dist
   min_i = i
  end if
 end do

 dist_min = 1.0e10
 do i = 0,130
  if(i .eq. min_i) cycle
  dist = sqrt((x - wall(i,0))**2.0 + (y - wall(i,1))**2.0)
  if(dist .lt. dist_min) then
   dist_min = dist
   neigh_i = i
  end if
 end do

 wall_dist = abs((wall(min_i,0) - x) * (wall(neigh_i,1) - y) - (wall(neigh_i,0) - x) * (wall(min_i,1) - y)) &
 &           / sqrt((wall(min_i,0) - wall(neigh_i,0))**2.0 + (wall(min_i,1) - wall(neigh_i,1))**2.0)
!write(*,*) 'wall_dist:',wall_dist

end subroutine wall_dist_t


! ********************************************************************
! **********************  二次元対応不可  ****************************
! ********************************************************************

subroutine circle_positioning         &
&          (cent0,cent1,radius,v0,v1,       &
&           inter_l,inter_t,p_type,nump,    &
&           nn,idim,ifluid,x,v,type,itypep, &
&           dis,grid,grid_num,wall_num,     &
&           nn_grid,g_neigh,grid_n,         &
&           f_type,phase_min,phase_max,inptnum,&
&           incount,droplet_num,cal_type0)
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
  integer, intent(inout) :: inptnum
  integer, intent(in)    :: incount
  integer, intent(inout) :: droplet_num(0:nn-1)
  integer, intent(inout) :: cal_type0(0:nn-1)

  ! subroutine_variable **********************************************
  integer              :: i,j,n,k,l,g
  real   , parameter   :: pi = acos(0.0)*2.0
  real                 :: rr,rr_min
  integer              :: num(0:1)
  real                 :: x_pos(0:1)
  real :: theta
  integer :: switch

  inptnum = 0

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

      switch = 0

      do g=1,grid_num
        if(x_pos(0) .ge. grid(g,0,0))then
          if(x_pos(0) .le. grid(g,0,1))then
            if(x_pos(1) .ge. grid(g,1,0))then
              if(x_pos(1) .le. grid(g,1,1))then
                switch = 1
                exit
              end if
            end if
          end if
        end if
      end do

      ! particle distance calculation ********************************
      rr_min = 1.0e10

      if(switch .eq. 1)then
        do l=1,g_neigh(g,0)
          n      = g_neigh(g,l)
          if(itypep(n) .ne. type(0))then
            rr     = dist_2d(x(n,0),x(n,1),x_pos(0),x_pos(1))
            rr_min = min(rr,rr_min)
          end if
          if((itypep(n) .eq. type(0)) .and. &
          &  (cal_type0(n) .ge. 1))then
            rr     = dist_2d(x(n,0),x(n,1),x_pos(0),x_pos(1))
            rr_min = min(rr,rr_min)
          end if
        end do

        ! particle_positioning *****************************************
        if(rr_min > dis)then
          do n=wall_num,nn-1
            if((itypep(n) .eq. type(0)) .and. &
            &  (cal_type0(n) .eq. 0))then
              do k=0,idim-1
                x(n,k) = x_pos(k)
              end do

              v(n,0) = v0
              v(n,1) = v1
              itypep(n) = p_type
              nump = max(nump,n)+1
              inptnum = inptnum + 1
              droplet_num(n) = incount
              cal_type0(n)   = 1

              if(itypep(n) .eq. type(4))f_type(1) = type(4)
              if(itypep(n) .eq. type(5))f_type(2) = type(5)
              if(itypep(n) .eq. type(3))f_type(3) = type(3)

!              phase_min = min(phase_min,p_type)
!              phase_max = max(phase_max,p_type)

               ! g_neigh_setting ****************************************
              call g_neigh_setting                        &
              &    (n,g,nn_grid,idim,grid_num,grid_n,g_neigh)
              exit
            end if
          end do
        end if
      end if
    end do
  end do

  return
end subroutine circle_positioning

subroutine circle_positioning_txt         &
&          (cent0,cent1,radius,v0,v1,       &
&           inter_l,inter_t,p_type,nump,    &
&           nn,idim,ifluid,x,v,type,itypep, &
&           dis,grid,grid_num,wall_num,     &
&           nn_grid,g_neigh,grid_n,         &
&           f_type,phase_min,phase_max,inptnum, &
&           incount,droplet_num,cal_type0,swi_solid)
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
  integer, intent(inout) :: inptnum
  integer, intent(in)    :: incount
  integer, intent(inout) :: droplet_num(0:nn-1)
  integer, intent(inout) :: cal_type0(0:nn-1)
  integer, intent(inout) :: swi_solid(0:nn-1)

  ! subroutine_variable **********************************************
  integer              :: i,j,n,k,l,g
  real   , parameter   :: pi = acos(0.0)*2.0
  real                 :: rr,rr_min
  integer              :: num(0:1)
  real                 :: x_pos(0:1)
  real :: theta
  integer :: switch
  real   , allocatable :: x_temp(:,:)
  integer              :: ptnum_txt
  real                 :: dis_txt
  real                 :: mvd_txt

  inptnum = 0

  ! search_point_number **********************************************
  open(11,file='./data/MPS/circle_position.dat')
    READ(11,*)ptnum_txt
    READ(11,*)dis_txt
    READ(11,*)mvd_txt


    allocate(x_temp(0:ptnum_txt-1,0:1))

    do i = 0,ptnum_txt-1
      do j = 0,idim-1
        read(11,*)x_temp(i,j)
      end do
    end do
  close(11)

  if(dis.ne.dis_txt)then
    write(*,*)"ERROR : different dis in txt_positioning"
    stop
  end if
  if(radius.ne.mvd_txt)then
    write(*,*)"ERROR : different MVD in txt_positioning"
    stop
  end if

  ! particle_positioning *********************************************
  do i = 0,ptnum_txt-1

    switch = 0

    x_pos(0) = cent0 + x_temp(i,0)
    x_pos(1) = cent1 + x_temp(i,1)

    do g=1,grid_num
      if(x_pos(0) .ge. grid(g,0,0))then
        if(x_pos(0) .le. grid(g,0,1))then
          if(x_pos(1) .ge. grid(g,1,0))then
            if(x_pos(1) .le. grid(g,1,1))then
              switch = 1
              exit
            end if
          end if
        end if
      end if
    end do

    ! particle distance calculation ********************************
    rr_min = 1.0e10

    if(switch .eq. 1)then
      do l=1,g_neigh(g,0)
        n      = g_neigh(g,l)
        if(itypep(n) .ne. type(0))then
          rr     = dist_2d(x(n,0),x(n,1),x_pos(0),x_pos(1))
          rr_min = min(rr,rr_min)
        end if
        if((itypep(n) .eq. type(0)) .and. &
        &  (cal_type0(n) .ge. 1))then
            rr     = dist_2d(x(n,0),x(n,1),x_pos(0),x_pos(1))
            rr_min = min(rr,rr_min)
          end if
      end do

      ! particle_positioning *****************************************
  !    if(rr_min > dis)then
        do n=wall_num,nn-1
          if((itypep(n) .eq. type(0)) .and. &
          &  (cal_type0(n) .eq. 0))then
            do k=0,idim-1
              x(n,k) = x_pos(k)
            end do
            v(n,0) = v0
            v(n,1) = v1
            itypep(n) = p_type
            swi_solid(n) = -1
            nump = max(nump,n)+1
            inptnum = inptnum + 1
            droplet_num(n) = incount
            cal_type0(n)   = 1

            if(itypep(n) .eq. type(4))f_type(1) = type(4)
            if(itypep(n) .eq. type(5))f_type(2) = type(5)
            if(itypep(n) .eq. type(3))f_type(3) = type(3)

!            phase_min = min(phase_min,p_type)
!            phase_max = max(phase_max,p_type)

             ! g_neigh_setting ****************************************
            call g_neigh_setting                        &
            &    (n,g,nn_grid,idim,grid_num,grid_n,g_neigh)
            exit
          end if
        end do
  !    end if
    end if

  end do

  ! array_allocation *************************************************
  deallocate(x_temp)

  return
end subroutine circle_positioning_txt

!**********************************************************************
!*** 計算をリスタートするための諸条件の読み込み                     ***
!**********************************************************************

subroutine read_restart_file         &
&           (nump,droplet_nn,time_sim_dble,    &
&           nn,idim,ifluid,x,v,type,itypep, &
&           dis,grid,grid_num,wall_num,     &
&           nn_grid,g_neigh,grid_n,         &
&           f_type,phase_min,phase_max,inptnum, &
&           incount,droplet_num,cal_type0,  &
&           t,gamma,swi_koeki,time_sim,nextin, &
&           fcount,count_remesh, &
&           droplet_swi,droplet_x,droplet_y,droplet_z, &
&           droplet_u,droplet_v,droplet_w, &
&           flag_temp)

  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: nn,idim,ifluid
  integer, intent(in) :: droplet_nn
!  integer, intent(in) :: p_type
  real   , intent(in) :: dis
!  real   , intent(in) :: cent0
!  real   , intent(in) :: cent1
!  real   , intent(in) :: v0
!  real   , intent(in) :: v1
!  real   , intent(in) :: radius
!  real   , intent(in) :: inter_l
!  real   , intent(in) :: inter_t
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
  integer, intent(inout) :: inptnum
  integer, intent(inout) :: incount
  integer, intent(inout) :: droplet_num(0:nn-1)
  integer, intent(inout) :: cal_type0(0:nn-1)
  real   , intent(inout) :: t(0:nn-1)
  real   , intent(inout) :: gamma(0:nn-1)
  integer, intent(inout) :: swi_koeki(0:nn-1)
  real   , intent(inout) :: time_sim
  real   , intent(inout) :: nextin
  integer, intent(inout) :: fcount
  integer, intent(inout) :: count_remesh
  integer, intent(inout) :: droplet_swi(0:droplet_nn)
  real   , intent(inout) :: droplet_x(0:droplet_nn)
  real   , intent(inout) :: droplet_y(0:droplet_nn)
  real   , intent(inout) :: droplet_z(0:droplet_nn)
  real   , intent(inout) :: droplet_u(0:droplet_nn)
  real   , intent(inout) :: droplet_v(0:droplet_nn)
  real   , intent(inout) :: droplet_w(0:droplet_nn)
  double precision, intent(inout) :: time_sim_dble
  integer,intent(inout)	:: flag_temp(0:nn-1)

  ! subroutine_variable **********************************************
  integer              :: i,j,n,k,l,g
  real   , parameter   :: pi = acos(0.0)*2.0
  real                 :: rr,rr_min
  integer, allocatable :: num(:)
!  real   , allocatable :: x_pos(:)
  real   , allocatable :: drr(:)
  real :: theta
  integer :: switch
  real   , allocatable :: x_temp(:,:)
  INTEGER             :: ptnum_txt
  real                 :: dis_txt
  real                 :: mvd_txt
  real                 :: time_sim_txt
  real                 :: nextin_txt
  integer              :: incount_txt
  integer              :: fcount_txt
  integer              :: count_remesh_txt
  real   , allocatable :: v_temp(:,:)
  real   , allocatable :: t_temp(:)
  integer, allocatable :: type_temp(:)
  real   , allocatable :: gamma_temp(:)
  integer, allocatable :: swi_koeki_temp(:)
  
  
  inptnum = 0
  
  ! array_allocation *************************************************
  allocate(num(0:1))
!  allocate(x_pos(0:1))
  
  ! search_point_number **********************************************
  open(11,file='./result/ICM/restart_file.dat')
    READ(11,*)ptnum_txt                 !計算粒子数
    READ(11,*)dis_txt                   !計算粒子径
    READ(11,*)mvd_txt                   !MVD
    read(11,*)time_sim_txt              !シミュレーション時間
    read(11,*)nextin_txt                !次の液滴の投入時刻
    read(11,*)incount_txt               !液滴投入数
    read(11,*)fcount_txt                !出力ファイル数
    read(11,*)count_remesh_txt          !リメッシュした回数
    read(11,*)time_sim_dble             !シミュレーション時間（倍精度）
    
    
    
!    allocate(x_temp(0:1,0:ptnum_txt-1))
!    allocate(v_temp(0:1,0:ptnum_txt-1))
!    allocate(t_temp(0:ptnum_txt-1))
!    allocate(type_temp(0:ptnum_txt-1))
!    allocate(gamma_temp(0:ptnum_txt-1))
!    allocate(swi_koeki_temp(0:ptnum_txt-1))

  if(dis.ne.dis_txt)then
    write(*,*)"ERROR : different dis in txt_positioning"
    stop
  end if
!  if(radius.ne.mvd_txt)then
!    write(*,*)"ERROR : different MVD in txt_positioning"
!    stop
!  end if

  time_sim = time_sim_txt
  nextin = nextin_txt
  incount = incount_txt
  fcount = fcount_txt
  count_remesh = count_remesh_txt
  nump = ptnum_txt
    
    do i = 0,ptnum_txt-1
      do j = 0,idim-1
        read(11,*)x(i,j)
      end do
    end do
    do i = 0,ptnum_txt-1
      do j = 0,idim-1
        read(11,*)v(i,j)
      end do
    end do
    do i = 0,ptnum_txt-1
      read(11,*)t(i)
    end do
    do i = 0,ptnum_txt-1
      read(11,*)itypep(i)
    end do
    do i = 0,ptnum_txt-1
      read(11,*) flag_temp(i)
    end do
    do i = 0,ptnum_txt-1
      read(11,*)gamma(i)
    end do
    do i = 0,ptnum_txt-1
      read(11,*)swi_koeki(i)
    end do
    do i = 0,ptnum_txt-1
      read(11,*)cal_type0(i)
    end do
    do i = 0,ptnum_txt-1
      read(11,*)droplet_num(i)
    end do

    do i = 0,incount
      read(11,*)droplet_swi(i)
    end do
    do i = 0,incount
      read(11,*)droplet_x(i)
      read(11,*)droplet_y(i)
      read(11,*)droplet_z(i)
      read(11,*)droplet_u(i)
      read(11,*)droplet_v(i)
      read(11,*)droplet_w(i)
    end do
  close(11)
  
  do i=0,nump-1
    if(itypep(n) .eq. type(4))f_type(1) = type(4)
    if(itypep(n) .eq. type(5))f_type(2) = type(5)
    if(itypep(n) .eq. type(3))f_type(3) = type(3)
  end do
  
  
  ! particle_positioning *********************************************
!  do i = 0,ptnum_txt-1
    
!    switch = 1
    
  
!    do g=1,grid_num
!      if(x_temp(0) .ge. grid(g,0,0))then
!        if(x_temp(0) .le. grid(g,0,1))then
!          if(x_temp(1) .ge. grid(g,1,0))then
!            if(x_temp(1) .le. grid(g,1,1))then
!              switch = 1
!              exit
!            end if
!          end if
!        end if
!      end if
!    end do
!  
!    ! particle distance calculation ********************************
!    rr_min = 1.0e10
!    
!    if(switch .eq. 1)then
!      do l=1,g_neigh(g,0)
!        n      = g_neigh(g,l)
!        if(itypep(n) .ne. type(0))then
!          rr     = dist_2d(x(0,n),x(1,n),x_temp(0),x_temp(1))
!          rr_min = min(rr,rr_min)
!        end if
!        if((itypep(n) .eq. type(0)) .and. &
!        &  (cal_type0(n) .ge. 1))then
!            rr     = dist_2d(x(0,n),x(1,n),x_temp(0),x_temp(1))
!            rr_min = min(rr,rr_min)
!          end if
!      end do
    
      ! particle_positioning *****************************************
  !    if(rr_min > dis)then
!        do n=0,nump-1
!          if((itypep(n) .eq. type(0)) .and. &
!          &  (cal_type0(n) .eq. 0))then
!            do k=0,idim-1
!              x(k,n) = x_temp(k,i)
!              v(k,n) = v_temp(k,i)
!            end do
            
!            t(n) = t_temp(i)
!            itypep(n) = type_temp(i)
!            gamma(n) = gamma_temp(i)
!            swi_koeki = swi_koeki_temp(i)
!            nump = max(nump,n)+1
!            cal_type0(n)   = 2
!            droplet_num(n) = incount 
!            inptnum = inptnum + 1           
!            cal_type0(n)   = 1

!            v(0,n) = v0
!            v(1,n) = v1
!            itypep(n) = p_type
!

!            
!            if(itypep(n) .eq. type(4))f_type(1) = type(4)
!            if(itypep(n) .eq. type(5))f_type(2) = type(5)
!            if(itypep(n) .eq. type(3))f_type(3) = type(3)
            
!            phase_min = min(phase_min,p_type)
!            phase_max = max(phase_max,p_type)
   
             ! g_neigh_setting ****************************************
!            call g_neigh_setting                        &
!            &    (n,g,nn_grid,idim,grid_num,grid_n,g_neigh)
!            exit
!          end if
!        end do
  !    end if
  !  end if
      
!  end do

  ! array_allocation *************************************************
  deallocate(num)
!  deallocate(drr)
!  deallocate(x_temp)
!  deallocate(v_temp)
!  deallocate(t_temp)
!  deallocate(type_temp)
!  deallocate(gamma_temp)
!  deallocate(swi_koeki_temp)

  return
end subroutine read_restart_file

!*********************************************************************
!*** write_restart_file                                            ***
!*********************************************************************
subroutine write_restart_file(   &
    &    nn,idim,dis,mvd,nump,wall_num,droplet_nn,time_sim_dble, &
    &    time_sim,nextin,incount,fcount,count_remesh,           &
    &    x,v,t,itypep,gamma,swi_koeki,cal_type0,droplet_num,    &
    &    droplet_swi,droplet_x,droplet_y,droplet_z,             &
    &    droplet_u,droplet_v,droplet_w)
  
  integer, intent(in) :: nn,idim
  integer, intent(in) :: droplet_nn
  real   , intent(in) :: dis
  real   , intent(in) :: mvd
  integer, intent(in) :: nump
  integer, intent(in) :: wall_num
  real   , intent(in) :: time_sim
  real   , intent(in) :: nextin
  integer, intent(in) :: incount
  integer, intent(in) :: fcount
  integer, intent(in) :: count_remesh
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  real   , intent(in) :: v(0:nn-1,0:idim-1)
  real   , intent(in) :: t(0:nn-1)
  integer, intent(in) :: itypep(0:nn-1)
  real   , intent(in) :: gamma(0:nn-1)
  integer, intent(in) :: swi_koeki(0:nn-1)
  integer, intent(in) :: cal_type0(0:nn-1)
  integer, intent(in) :: droplet_num(0:nn-1)
  integer, intent(in) :: droplet_swi(0:droplet_nn)
  real   , intent(in) :: droplet_x(0:droplet_nn)
  real   , intent(in) :: droplet_y(0:droplet_nn)
  real   , intent(in) :: droplet_z(0:droplet_nn)
  real   , intent(in) :: droplet_u(0:droplet_nn)
  real   , intent(in) :: droplet_v(0:droplet_nn)
  real   , intent(in) :: droplet_w(0:droplet_nn)
  double precision, intent(in) :: time_sim_dble
  
  integer  :: i,j,k
  
  
  open(11,file='./result/MPS/restart_file.dat')
    write(11,*)nump-wall_num    !-temp_nump
    write(11,*)dis
    write(11,*)mvd
    write(11,*)time_sim
    write(11,*)nextin
    write(11,*)incount
    write(11,*)fcount
    write(11,*)count_remesh
    write(11,*)time_sim_dble

    do i = wall_num, nump-1
      do j = 0,idim-1
        write(11,*)x(i,j)
      end do
    end do
    do i = wall_num, nump-1
      do j = 0,idim-1
        write(11,*)v(i,j)
      end do 
    end do
    do i = wall_num,nump-1
      write(11,*)t(i)
    end do
    do i = wall_num,nump-1
      write(11,*)itypep(i)
    end do
    do i = wall_num,nump-1
      write(11,*)gamma(i)
    end do
    do i = wall_num,nump-1
      write(11,*)swi_koeki(i)
    end do
    do i = wall_num,nump-1
      write(11,*)cal_type0(i)
    end do
    do i = wall_num,nump-1
      write(11,*)droplet_num(i)
    end do
    do i = 0,incount
      write(11,*)droplet_swi(i)
    end do
    do i = 0,incount
      write(11,*)droplet_x(i)
      write(11,*)droplet_y(i)
      write(11,*)droplet_z(i)
      write(11,*)droplet_u(i)
      write(11,*)droplet_v(i)
      write(11,*)droplet_w(i)
    end do
  close(11)
  
  return
end subroutine write_restart_file

!*********************************************************************
!*** write_sphere_position                                         ***
!*********************************************************************
subroutine write_circle_position(   &
    &    nn,idim,dis,mvd,x,nump,temp_nump, &
    &    randomx,randomy,randomz           &
    &    )

  integer, intent(in) :: nn,idim
  real   , intent(in) :: dis
  real   , intent(in) :: mvd
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  integer, intent(in) :: nump
  integer, intent(in) :: temp_nump
  real   , intent(in) :: RANDOMX
  real   , intent(in) :: RANDOMY
  real   , intent(in) :: RANDOMZ

  integer  :: i,j,k
  real     :: temp_x(0:idim-1,0:nn-1)

  DO i = temp_nump, nump-1
    temp_x(0,i) = x(i,0) - RANDOMX
    temp_x(1,i) = x(i,1) - RANDOMY
!    temp_x(2,i) = x(i,2) - RANDOMZ
  end do


  open(11,file='./data/MPS/circle_position.dat')
    write(11,*)nump-temp_nump
    write(11,*)dis
    write(11,*)mvd

    DO i = temp_nump, nump-1
      DO j = 0,idim-1
        write(11,*)temp_x(j,i)
      END DO
    END DO

  close(11)

  return
end subroutine write_circle_position

!*********************************************************************
!*** box_positioning                                            ***
!*********************************************************************
subroutine box_positioning         &
&          (hor0,hor1,ver0,ver1,v0,v1, &
&           interval,dis_min, p_type,nump,    &
&           nn,idim,ifluid,x,v,type,itypep, &
&           dis,grid,grid_num,wall_num,     &
&           nn_grid,g_neigh,grid_n,         &
&           f_type,phase_min,phase_max)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: nn,idim,ifluid
  integer, intent(in) :: p_type
  real   , intent(in) :: dis
  real   , intent(in) :: hor0
  real   , intent(in) :: hor1
  real   , intent(in) :: ver0
  real   , intent(in) :: ver1
  real   , intent(in) :: v0
  real   , intent(in) :: v1
  real   , intent(in) :: interval
  real   , intent(in) :: dis_min
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
  integer         :: num(0:1)
  integer         :: pos_sw(0:1)
  real            :: x_pos(0:1)
  real            :: rr,rr_min
  real            :: x_min(0:1)
  real            :: x_max(0:1)

  ! calculation_start ************************************************
  ! data_copy ********************************************************
  x_min(0) = min(hor0,hor1)
  x_max(0) = max(hor0,hor1)
  x_min(1) = min(ver0,ver1)
  x_max(1) = max(ver0,ver1)

  ! reference_position ***********************************************
  ! x_direction
  if(x_min(0) == hor0)then
    pos_sw(0) = 0
  else
    pos_sw(0) = 1
  end if

  ! y_direction
  if(x_min(1) == ver0)then
    pos_sw(1) = 0
  else
    pos_sw(1) = 1
  end if
  if(interval <= 0.0)then
    write(900,'(a)')"error in subroutine box_positioning_2d"
    write(900,'(a)')"interval is negative or zero"
    stop
  end if

  ! number_of_searching_point ****************************************
  do k=0,1
    num(k) = int(abs(x_max(k)-x_min(k))/interval)
  end do

  do i=num(0)-1,0,-1
    do j=0,num(1)-1
      ! serch position
      ! searching_position
      if(pos_sw(0) == 0)then
        x_pos(0) = x_min(0)+real(i)*interval
      else
        x_pos(0) = x_max(0)-real(i)*interval
      end if

      if(pos_sw(1) == 0)then
        x_pos(1) = x_min(1)+real(j)*interval
      else
        x_pos(1) = x_max(1)-real(j)*interval
      end if

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
      rr_min = 1.0e10

      do l=1,g_neigh(g,0)
        n      = g_neigh(g,l)
        if(itypep(n) .ne. type(0))then
          rr     = dist_2d(x(n,0),x(n,1),x_pos(0),x_pos(1))
          rr_min = min(rr,rr_min)
        end if
      end do

      if(rr_min > dis_min*dis)then
        do n=wall_num,nn-1
          if(itypep(n) .eq. type(0))then
            do k=0,idim-1
              x(n,k) = x_pos(k)
            end do

            v(n,0) = v0
            v(n,1) = v1
            itypep(n) = p_type
            nump = max(nump,n)+1

            if(itypep(n) .eq. type(4))f_type(1) = type(4)
            if(itypep(n) .eq. type(5))f_type(2) = type(5)
            if(itypep(n) .eq. type(3))f_type(3) = type(3)

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
end subroutine box_positioning

!*********************************************************************
!*** collision                                            ***
!*********************************************************************
subroutine icing_judge_collision                         &
&          (i,nn,idim,ifluid,neighbor,nump,    &
&           itypep,type,cal_sw,              &
&           den,neigh,x,v,dist_min,coll_rat, &
&           dt,phase,dis,wall_num,wall_n,wall_anc, &
&           grid_p,nn_wall,grid_num,g_wall,wall_nor, &
&           ice_dist_min,icing_col,col_wall,j_temp,ker_r)
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
  integer, intent(in) :: wall_n
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  integer, intent(in) :: grid_p(0:nn-1)
  integer, intent(in) :: nn_wall
  integer, intent(in) :: grid_num
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  real   , intent(in)    :: ice_dist_min
  integer, intent(inout) :: icing_col(0:nn-1)
  real   , intent(inout) :: col_wall
  real   , intent(inout) :: j_temp
  real   , intent(in) :: ker_r(0:3)



  ! subroutine_variable **********************************************
  integer :: j,l,k,m
  real    :: m1,m2,mm
  real    :: rr,rr2
  real    :: vabs,vrat
  real    :: col_vel(0:idim-1,1:2)
  real    :: distance(0:idim-1)
  real    :: vel_g(0:idim-1)
  real    :: vel_r(0:idim-1)
  real    :: vel_m(0:idim-1)

  real    :: distance_wall
  real    :: vel_a
  real    :: vel_n
  real    :: cos_t
  real    :: theta
  real    :: n_a,h_a
  real    :: x1,x2,x3
  real    :: y1,y2,y3
  real    :: z1,z2,z3
  real    :: nd
  real    :: v_ref0,v_ref1
  real    :: nx(0:idim-1)
  real    :: hx(0:idim-1)
  real    :: v_sub(0:idim-1)
  real, parameter   :: zero = 1.0e-10
  real, parameter   :: pi   = acos(-1.0)
  real, parameter   :: col_rat = 0.0
  real    :: r_kernel

  r_kernel = maxval(ker_r)


  ! particle_collision ***********************************************
!  do l=1,neigh(i,0)
!    j = neigh(i,l)

!    rr2 = 0.0
!    if((itypep(j) .eq. type(2)))then
!      do m=0,idim-1
!        distance(m) = x(j,m)-x(i,m)
!        rr2   = rr2+distance(m)**2
!      end do
!      rr = sqrt(rr2)
!      if(rr .le. r_kernel)then
!        icing_col(i) = 1
!        j_temp = j
!      end if
!    end if
!  end do

  ! particle-wall_collision *******************************************
  if(icing_col(i).ne.1)then
    do k=1,g_wall(grid_p(i),0)
      j = g_wall(grid_p(i),k)
      ! wall_normal_vector
      x1 = wall_anc(0,j,0)
      x2 = wall_anc(0,j,1)
!      x3 = wall_anc(0,j,2)

      y1 = wall_anc(1,j,0)
      y2 = wall_anc(1,j,1)
!      y3 = wall_anc(1,j,2)

!      z1 = wall_anc(2,j,0)
!      z2 = wall_anc(2,j,1)
!      z3 = wall_anc(2,j,2)

      nx(0) = -wall_nor(0,j)
      nx(1) = -wall_nor(1,j)
!      nx(2) = -wall_nor(2,j)
!      nd    = -(nx(0)*x1+nx(1)*y1+nx(2)*z1)
      ND    = -(NX(0)*X1+NX(1)*Y1)

!      n_a = 1.0/max(zero,sqrt(nx(0)**2+nx(1)**2+nx(2)**2))
      N_A = 1.0/MAX(ZERO,SQRT(NX(0)**2+NX(1)**2))

      do l=0,idim-1
        nx(l) = nx(l)*n_a
      end do

      ! particle-wall_distance
!      distance_wall = abs(nx(0)*x(i,0)+nx(1)*x(i,1)+nx(2)*x(i,2)+nd) &
!      &               /sqrt(nx(0)**2+nx(1)**2+nx(2)**2)

       DISTANCE_WALL = ABS(NX(0)*x(i,0)+NX(1)*x(i,1)+ND) &
        &          /SQRT(NX(0)**2+NX(1)**2)


      ! collision_judgment
      if(distance_wall .le. dis*ice_dist_min)then
        icing_col(i) = 1
        col_wall = 1
      end if
    end do
  end if


  return
end subroutine icing_judge_collision

!*********************************************************************
!*** dom_search                                    ***
!*********************************************************************
subroutine dom_search                  &
&          (nn,nn_grid,idim,ifluid,nump,  &
&           grid_num,grid_n,grid,         &
&           g_neigh,x,itypep,type,grid_p, &
&           xmax,xmin,ymax,ymin,zmax,zmin)
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
  real, intent(inout) :: xmax
  real, intent(inout) :: xmin
  real, intent(inout) :: ymax
  real, intent(inout) :: ymin
  real, intent(inout) :: zmax
  real, intent(inout) :: zmin

  ! subroutine_variable **********************************************
  integer :: i,j,k
  integer :: switch

  ! grid_serch *******************************************************
  xmax = 0.0
  xmin = 0.0
  ymax = 0.0
  ymin = 0.0
  zmax = 0.0
  zmin = 0.0


  do j=1,grid_num
    if(xmin .ge. grid(j,0,0))then
      xmin = grid(j,0,0)
    end if
    if(xmax .le. grid(j,0,1))then
      xmax = grid(j,0,1)
    end if
    if(ymin .ge. grid(j,1,0))then
      ymin = grid(j,1,0)
    end if
    if(ymax .le. grid(j,1,1))then
      ymax = grid(j,1,1)
    end if
!    if(zmin .ge. grid(j,2,0))then
!      zmin = grid(j,2,0)
!    end if
!    if(zmax .le. grid(j,2,1))then
!      zmax = grid(j,2,1)
!    end if

  end do

  write(*,*)"xmin=",xmin,"xmax=",xmax
  write(*,*)"ymin=",ymin,"ymax=",ymax
!  write(*,*)"zmin=",zmin,"zmax=",zmax


  return
end subroutine dom_search

! ********************************************************************
! *****************  カップリング用 サブルーチン  ********************
! ********************************************************************


subroutine switch_grid                       &
&          (i,nn,idim,ifluid,neighbor,nump,    &
&           itypep,type,cal_sw,              &
&           den,neigh,x,v,dist_min,coll_rat, &
&           dt,phase,dis,wall_num,wall_n,wall_anc, &
&           grid_p,nn_wall,grid_num,g_wall,wall_nor, &
&           ker_r,mps_swi,droplet_num,change_incount)
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
  real   , intent(in) :: x(0:nn-1,0:idim-1)
  real   , intent(in) :: v(0:nn-1,0:idim-1)
  integer, intent(in) :: wall_n
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  integer, intent(in) :: grid_p(0:nn-1)
  integer, intent(in) :: nn_wall
  integer, intent(in) :: grid_num
  integer, intent(in) :: g_wall(1:grid_num,0:nn_wall)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)
  real   , intent(in) :: ker_r(0:3)
  integer, intent(in) :: droplet_num(0:nn-1)
  integer, intent(inout) :: mps_swi
  integer, intent(inout) :: change_incount

  ! subroutine_variable **********************************************
  integer :: j,l,k,m
  real    :: m1,m2,mm
  real    :: rr,rr2
  real    :: vabs,vrat
  real    :: col_vel(0:idim-1,1:2)
  real    :: distance(0:idim-1)
  real    :: vel_g(0:idim-1)
  real    :: vel_r(0:idim-1)
  real    :: vel_m(0:idim-1)

  real    :: distance_wall
  real    :: vel_a
  real    :: vel_n
  real    :: cos_t
  real    :: theta
  real    :: n_a,h_a
  real    :: x1,x2,x3
  real    :: y1,y2,y3
  real    :: z1,z2,z3
  real    :: nd
  real    :: v_ref0,v_ref1
  real    :: nx(0:idim-1)
  real    :: hx(0:idim-1)
  real    :: v_sub(0:idim-1)
  real, parameter   :: zero = 1.0e-10
  real, parameter   :: pi   = acos(-1.0)
  real    :: r_kernel

  r_kernel = maxval(ker_r)

  ! particle_collision ***********************************************
  do l=1,neigh(i,0)
    j = neigh(i,l)

    rr2 = 0.0
    if((itypep(j) .eq. type(2)) .or. &
    &  (itypep(j) .eq. type(4)))then
      do m=0,idim-1
        distance(m) = x(j,m)-x(i,m)
        rr2   = rr2+distance(m)**2
      end do
      rr = sqrt(rr2)
      if(rr .le. r_kernel)then
        mps_swi = 1
        change_incount = droplet_num(i)
      end if
    end if
  end do

  ! particle-wall_collision *******************************************
  do k=1,g_wall(grid_p(i),0)
    j = g_wall(grid_p(i),k)
    ! wall_normal_vector
    x1 = wall_anc(0,j,0)
    x2 = wall_anc(0,j,1)
!    x3 = wall_anc(0,j,2)

    y1 = wall_anc(1,j,0)
    y2 = wall_anc(1,j,1)
!    y3 = wall_anc(1,j,2)

!    z1 = wall_anc(2,j,0)
!    z2 = wall_anc(2,j,1)
!    z3 = wall_anc(2,j,2)

    nx(0) = -wall_nor(0,j)
    nx(1) = -wall_nor(1,j)
!    nx(2) = -wall_nor(2,j)
!    nd    = -(nx(0)*x1+nx(1)*y1+nx(2)*z1)
    ND    = -(NX(0)*X1+NX(1)*Y1)

!    n_a = 1.0/max(zero,sqrt(nx(0)**2+nx(1)**2+nx(2)**2))
    N_A = 1.0/MAX(ZERO,SQRT(NX(0)**2+NX(1)**2))

    do l=0,idim-1
      nx(l) = nx(l)*n_a
    end do

    ! particle-wall_distance
!    distance_wall = abs(nx(0)*x(i,0)+nx(1)*x(i,1)+nx(2)*x(i,2)+nd) &
!    &               /sqrt(nx(0)**2+nx(1)**2+nx(2)**2)

     DISTANCE_WALL = ABS(NX(0)*x(i,0)+NX(1)*x(i,1)+ND) &
      &          /SQRT(NX(0)**2+NX(1)**2)


    ! collision_judgment
    if(distance_wall .le. r_kernel)then
      mps_swi = 1
      change_incount = droplet_num(i)

    end if
  end do

  return
end subroutine switch_grid


subroutine grid_division_type0                  &
&          (nn,nn_grid,idim,ifluid,nump,  &
&           grid_num,grid_n,grid,         &
&           g_neigh,x,itypep,type,grid_p,cal_type0)
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
  integer, intent(in)    :: cal_type0(0:nn-1)
  ! subroutine_variable **********************************************
  integer :: ip,jg
  integer :: around

  ! この変数達をallocatableにするとバグる、、、なんで？？
  integer :: jj(0:nn-1)
  integer :: switch2(0:nn-1)
  integer :: shift(9)
  logical :: mask(0:nump-1)

  ! array_initialize *************************************************

  shift(1) = -grid_n(0)-1
  shift(2) = -grid_n(0)
  shift(3) = -grid_n(0)+1
  shift(4) =           -1
  shift(5) =            0
  shift(6) =           +1
  shift(7) = +grid_n(0)-1
  shift(8) = +grid_n(0)
  shift(9) = +grid_n(0)+1

  switch2(:) = 1

  do jg=1,grid_num
    g_neigh(jg,0) = 0
  end do

  do ip=0,nump-1
    mask(ip) = .false.
  end do
  do ip=0,nump-1
    if(((itypep(ip) .eq. type(0)) .and. &
    &  (cal_type0(ip) .eq. 2))    .or.  &
    &  (itypep(ip) .ge. type(2)))then
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
end subroutine grid_division_type0
! ***********************************************************************
! ****************  これは他次元も対応可能  *****************************
! ***********************************************************************

subroutine set_neigh_type0                   &
&          (i,                         &
&           nn,idim,ifluid,neighbor,x, &
&           type,neigh,itypep,         &
&           ker_r,nn_grid,             &
&           grid_num,grid_p,g_neigh,cal_type0)
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
  integer, intent(in)    :: cal_type0(0:nn-1)

  ! subroutine_variable *********************************************
  integer :: j,k,l,m
  real    :: rr,rr2
  real    :: drr(0:idim-1)

  ! array_initialize ************************************************
  neigh(i,0)= 0

  ! neighboring_particle_count **************************************
  ! particle_number_copy
  k = grid_p(i)

  if((itypep(i) .eq. type(0)) .and. &
  &  (cal_type0(i) .eq. 2))then
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
end subroutine set_neigh_type0
end module mod_icing_dim

