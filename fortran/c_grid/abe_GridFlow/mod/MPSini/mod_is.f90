!*********************************************************************
!****** MODULE_FOR_INITIALFILE_MAKER                               ***
!****** VER. 2014.08.13                                            ***
!*********************************************************************
MODULE MOD_IS
  USE MOD_IS_DIM
  IMPLICIT NONE
  ! SUBROUTINE_DEFINITION ********************************************
  PUBLIC :: ALLOCATION
  PUBLIC :: FLAG_SETTING
  PUBLIC :: OUTPUT_PARA
  PUBLIC :: OUTPUT_PARA_GRID
  PUBLIC :: OUTPUT_PARA_WALL
  PUBLIC :: OUTPUT_INITIAL_FILE
CONTAINS

!*********************************************************************
!*** ARRAY_ALLOCATION                                              ***
!*********************************************************************
SUBROUTINE ALLOCATION                   &
&          (NN,IDIM,NUMP,IFLUID,        &
&           TYPE,DEN,VIS,F_TYPE,ITYPEP, &
&           DOM_SIZE,GRID_N,NON_CAL,    &
&           KER_C,X,OUT_TYPE,T,WALL_N)
  IMPLICIT NONE
  ! MAINROUTINE_VARIABLE *********************************************
  INTEGER, INTENT(IN)               :: NN,IDIM,NUMP
  INTEGER, INTENT(IN)               :: IFLUID
  INTEGER, INTENT(IN)               :: WALL_N
  REAL   , INTENT(OUT), ALLOCATABLE :: DOM_SIZE(:,:)
  INTEGER, INTENT(OUT), ALLOCATABLE :: TYPE(:)
  REAL   , INTENT(OUT), ALLOCATABLE :: DEN(:)
  REAL   , INTENT(OUT), ALLOCATABLE :: VIS(:)
  INTEGER, INTENT(OUT), ALLOCATABLE :: F_TYPE(:)
  INTEGER, INTENT(OUT), ALLOCATABLE :: GRID_N(:)
  INTEGER, INTENT(OUT), ALLOCATABLE :: ITYPEP(:)
  REAL   , INTENT(OUT), ALLOCATABLE :: NON_CAL(:)
  REAL   , INTENT(OUT), ALLOCATABLE :: KER_C(:)
  REAL   , INTENT(OUT), ALLOCATABLE :: X(:,:)
  INTEGER, INTENT(OUT), ALLOCATABLE :: OUT_TYPE(:)
  REAL   , INTENT(OUT), ALLOCATABLE :: T(:)

  ! CALCULATION_START ************************************************
  WRITE(*,'(A)')"********** SUBROUTINE_ALLOCATION START"

  ! ARRAY_ALLOCATION *************************************************
  ALLOCATE(DOM_SIZE(0:IDIM-1,0:1))
  ALLOCATE(TYPE(0:IFLUID))
  ALLOCATE(DEN(0:IFLUID))
  ALLOCATE(VIS(0:IFLUID))
  ALLOCATE(F_TYPE(1:3))
  ALLOCATE(GRID_N(0:IDIM-1))
  ALLOCATE(ITYPEP(1:NN))
  ALLOCATE(NON_CAL(0:IDIM-1))
  ALLOCATE(KER_C(0:3))
  ALLOCATE(X(0:IDIM-1,1:NN))
  ALLOCATE(OUT_TYPE(1:NN))
  ALLOCATE(T(1:NN))

  ! ARRAY_INITIALIZE *************************************************
  DOM_SIZE(0:IDIM-1,0:1)            = 0.0
  DEN(0:IFLUID)                     = 0.0
  VIS(0:IFLUID)                     = 0.0
  NON_CAL(0:IDIM-1)                 = 0.0
  KER_C(0:3)                        = 0.0
  X(0:IDIM-1,1:NN)                  = 0.0
  T(1:NN)                           = 0.0
  OUT_TYPE(1:NN)                    = 0
  TYPE(0:IFLUID)                    = 0
  F_TYPE(1:3)                       = 0
  GRID_N(0:IDIM-1)                  = 0
  ITYPEP(1:NN)                      = 0

  ! CALCULATION_END **************************************************
  WRITE(*,'(A)')"********** SUBROUTINE_ALLOCATION END"
  WRITE(*,'(A)')""

  RETURN
END SUBROUTINE ALLOCATION

!*********************************************************************
!*** PARTICLE_FLAG_SETTING                                         ***
!*********************************************************************
SUBROUTINE FLAG_SETTING                       &
&          (NN,NUMP,IFLUID,ITYPEP,TYPE,F_TYPE)
  IMPLICIT NONE
  ! MAINROUTINE_VARIABLE *********************************************
  INTEGER, INTENT(IN)    :: NN,NUMP
  INTEGER, INTENT(IN)    :: IFLUID
  INTEGER, INTENT(IN)    :: TYPE(0:IFLUID)
  INTEGER, INTENT(INOUT) :: ITYPEP(1:NN)
  INTEGER, INTENT(INOUT) :: F_TYPE(1:3)

  ! SUBROUTINE_VARIABLE **********************************************
  INTEGER :: I

  ! CALCULATION_START ************************************************
  WRITE(*,'(A)')"********** SUBROUTINE_FLAG_SETTING START"

  ! FLAG_SETTING *****************************************************
  DO I=1,NUMP
    IF(ITYPEP(I) .EQ. TYPE(4))F_TYPE(1) = TYPE(4)
    IF(ITYPEP(I) .EQ. TYPE(5))F_TYPE(2) = TYPE(5)
    IF(ITYPEP(I) .EQ. TYPE(3))F_TYPE(3) = TYPE(3)
  END DO

  ! CALCULATION_END **************************************************
  WRITE(*,'(A)')"********** SUBROUTINE_FLAG_SETTING END"
  WRITE(*,'(A)')""

  RETURN
END SUBROUTINE FLAG_SETTING

!*********************************************************************
!*** output_paraview_file                                          ***
!*********************************************************************
subroutine output_para              &
&          (nn,idim,ifluid,nump,x,t,out_type,itypep,type)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: nn,idim,nump
  integer, intent(in) :: ifluid
  integer, intent(in) :: type(0:ifluid)
  real   , intent(in) :: x(0:idim-1,1:nn)
  integer, intent(in) :: out_type(1:nn)
  integer, intent(in) :: itypep(1:nn)
  real   , intent(in) :: t(1:nn)

  ! subroutine_variable **********************************************
  integer             :: i,n
  integer             :: out_count  = 0
  integer             :: real_num
  real, allocatable   :: out_pos(:,:)
  real, allocatable   :: out_real(:,:)
  character(len = 40) :: datfile

  ! calculation_start ************************************************
  write(*,'(a)')"********** subroutine_output_para start"

  ! data_number ******************************************************
  real_num = 2

  ! array_allocation *************************************************
  allocate(out_pos(0:2,1:nn))
  allocate(out_real(0:real_num-1,1:nn))

  ! array_initialize *************************************************
  out_pos(0:2,1:nn)           = 0.0
  out_real(0:real_num-1,1:nn) = 0.0

  ! data_copy ********************************************************
  do n=1,nump
    if(out_type(n) == 0)then
      if(itypep(n) .ge. type(3))then
        out_count = out_count+1
        ! position
        do i=0,2
          if(i <= idim-1)then
            out_pos(i,out_count) = x(i,n)
          else
            out_pos(i,out_count) = 0.0
          end if
        end do

        ! particle_type
        out_real(0,out_count) = real(itypep(n))

        ! temperature
        out_real(1,out_count) = t(n)
      end if
    end if
  end do

  ! file_name ********************************************************
  write(datfile,'(a)')'./data/MPSini/mps000.vtk'

  ! data_number_check ************************************************
  if(out_count > 100000000)then
    write(900,'(a)')"Error in subroutine_output_para(mod_is)"
    write(900,'(a)')"out_count : over 100000000"
    stop
  end if

  ! output_file ******************************************************
  open(11,file = datfile)
  write(11,'(a)')'# vtk DataFile Version 3.0'
  write(11,'(a)')'vtk output'
  write(11,'(a)')'ASCII'
  write(11,'(a)')'DATASET POLYDATA'
  write(11,'(a,i8,a)')'POINTS',out_count,' float'

  do n=1,out_count
    write(11,*)out_pos(0:2,n)
  end do

  write(11,'(a,i9)')'POINT_DATA',out_count
  write(11,'(a   )')'FIELD attributes 6'

  write(11,'(a,i9,a)')'particle_type  1',out_count,'float'
  do n=1,out_count
    write(11,*)out_real(0,n)
  end do

  write(11,'(a,i9,a)')'velocity  3',out_count,'float'
  do n=1,out_count
    write(11,*)0.0,0.0,0.0
  end do

  write(11,'(a,i9,a)')'pressure  1',out_count,'float'
  do n=1,out_count
    write(11,*)0.0
  end do

  write(11,'(a,i9,a)')'particle_number_density  1',out_count,'float'
  do n=1,out_count
    write(11,*)0.0
  end do

  write(11,'(a,i9,a)')'surface_tension  3',out_count,'float'
  do n=1,out_count
    write(11,*)0.0,0.0,0.0
  end do

  write(11,'(a,i9,a)')'temperature  1',out_count,'float'
  do n=1,out_count
    write(11,*)out_real(1,n)
  end do

  close(11)

  ! array_deallocation ***********************************************
  deallocate(out_pos)
  deallocate(out_real)

  ! calculation_end **************************************************
  write(*,'(a)')"********** subroutine_output_para end"
  write(*,'(a)')""

  return
end subroutine output_para

!*********************************************************************
!*** output_grid_paraview_file                                     ***
!*********************************************************************
subroutine output_para_grid                    &
&          (idim,grid_num,wall_n,grid,grid_dist,near_wall)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: idim
  integer, intent(in) :: grid_num
  integer, intent(in) :: wall_n
  real   , intent(in) :: grid(1:grid_num,0:idim-1,0:1)
  real   , intent(in) :: grid_dist(1:grid_num,0:4*(idim-1)-1)
  integer, intent(in) :: near_wall(1:grid_num,0:4*(idim-1)-1)

  ! subroutine_variable **********************************************
  integer             :: n,i,j,k
  integer             :: dnum
  integer              :: gnum
  real ,parameter     :: mar = 0.0
  real, allocatable   :: grid_sub(:,:,:)
  character(len = 40) :: datfile
  character(len = 200) :: d_num

  ! calculation_start ************************************************
  write(*,'(a)')"********** subroutine_output_para_grid start"


  ! array_allocation & data_copy *************************************
  allocate(grid_sub(1:grid_num,0:2,0:1))

  gnum = 0

  do i=1,grid_num
  	if(minval(near_wall(i,0:4*(idim-1)-1)) .ge. 0)then
    do j=0,2
      do k=0,1
        if(j .le. idim-1)then
          grid_sub(gnum,j,k) = grid(i,j,k)
        else
          grid_sub(gnum,j,k) = 0.0
        end if
      end do
    end do
    
    gnum = gnum+1
    
    end if
  end do
  ! data_number ******************************************************
  dnum = (gnum-1)*4*(idim-1)

  ! file_name ********************************************************
  write(datfile,'(a)')'./data/MPSini/grid.vtk'

  ! output_grid_paraview_file ****************************************
  open(11,file = datfile)
  write(11,'(a)')'# vtk DataFile Version 3.0'
  write(11,'(a)')'vtk output'
  write(11,'(a)')'ASCII'
  write(11,'(a)')'DATASET POLYDATA'
  write(11,'(a,i8,a)')'POINTS',dnum,' float'

  do n=1,gnum-1
    if(idim .eq. 2)then
      write(11,*)grid_sub(n,0,0),grid_sub(n,1,0),'0.0'
      write(11,*)grid_sub(n,0,1),grid_sub(n,1,0),'0.0'
      write(11,*)grid_sub(n,0,1),grid_sub(n,1,1),'0.0'
      write(11,*)grid_sub(n,0,0),grid_sub(n,1,1),'0.0'
    else
      write(11,*)grid_sub(n,0,0),grid_sub(n,1,0),grid_sub(n,2,0)
      write(11,*)grid_sub(n,0,1),grid_sub(n,1,0),grid_sub(n,2,0)
      write(11,*)grid_sub(n,0,0),grid_sub(n,1,1),grid_sub(n,2,0)
      write(11,*)grid_sub(n,0,1),grid_sub(n,1,1),grid_sub(n,2,0)
      write(11,*)grid_sub(n,0,0),grid_sub(n,1,0),grid_sub(n,2,1)
      write(11,*)grid_sub(n,0,1),grid_sub(n,1,0),grid_sub(n,2,1)
      write(11,*)grid_sub(n,0,0),grid_sub(n,1,1),grid_sub(n,2,1)
      write(11,*)grid_sub(n,0,1),grid_sub(n,1,1),grid_sub(n,2,1)
    end if
  end do

!  write(11,'(a,i9)')'POINT_DATA',dnum
!
!  write(11,'(a)')'FIELD attributes 2'
!
!  write(11,'(a,i9,a)')'wall_dist 1',dnum,'float'
!  do n=1,grid_num
!    do j=0,4*(idim-1)-1
!      write(11,*)grid_dist(n,j)
!    end do
!  end do
!
!  write(11,'(a,i9,a)')'nearest_wall 1',dnum,'float'
!  do n=1,grid_num
!    do j=0,4*(idim-1)-1
!      write(11,*)real(near_wall(n,j))
!    end do
!  end do

  close(11)

  ! deallocationg ****************************************************
  deallocate(grid_sub)

  ! calculation_end **************************************************
  write(*,'(a)')"********** subroutine_output_para_grid end"
  write(*,'(a)')""

  return
end subroutine output_para_grid

!*********************************************************************
!*** output_wall_paraview_file                                     ***
!*********************************************************************
subroutine output_para_wall              &
&          (idim,grid_num,wall_n,wall_anc)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in) :: idim
  integer, intent(in) :: grid_num
  integer, intent(in) :: wall_n
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)

  ! subroutine_variable **********************************************
  integer             :: n,i,j,k
  integer             :: dnum
  real, allocatable   :: grid_sub(:,:,:)
  real, allocatable   :: grid_vec(:,:,:)
  character(len = 40) :: datfile

  ! calculation_start ************************************************
  write(*,'(a)')"********** subroutine_output_para_wall start"

  ! data_number ******************************************************
  dnum = wall_n*(idim)

  ! array_alocation & data_copy **************************************
  allocate(grid_sub(0:2,0:wall_n-1,0:idim-1))
  allocate(grid_vec(0:2,0:wall_n-1,0:idim-1))

  do i=0,2
    do j=0,wall_n-1
      do k=0,idim-1
        if(i .le. idim-1)then
          grid_sub(i,j,k) = wall_anc(i,j,k)
        else
          grid_sub(i,j,k) = 0.0
        end if

        if(i .le. idim-1)then
          if(k .ne. idim-1)then
            grid_vec(i,j,k) = wall_anc(i,j,k+1)-wall_anc(i,j,k)
          else
            grid_vec(i,j,k) = wall_anc(i,j,0)-wall_anc(i,j,k)
          end if
        else
          grid_vec(i,j,k) = 0.0
        end if
      end do
    end do
  end do

  ! file_name ********************************************************
  write(datfile,'(a)')'./data/MPSini/wall.vtk'

  ! output_paraview_file *********************************************
  open(11,file = datfile)
  write(11,'(a)')'# vtk DataFile Version 3.0'
  write(11,'(a)')'vtk output'
  write(11,'(a)')'ASCII'
  write(11,'(a)')'DATASET POLYDATA'
  write(11,'(a,i8,a)')'POINTS',dnum,' float'

  do n=0,wall_n-1
    do j=0,idim-1
      write(11,*)grid_sub(0:2,n,j)
    end do
  end do

  write(11,'(a,i9)')'POINT_DATA',dnum

  write(11,'(a)')'FIELD attributes 1'

  write(11,'(a,i9,a)')'wall_dist 3',dnum,'float'
  do n=0,wall_n-1
    do j=0,idim-1
      write(11,*)grid_vec(0:2,n,j)
    end do
  end do
  close(11)

  ! array_deallocationg **********************************************
  deallocate(grid_sub)
  deallocate(grid_vec)

  ! calculation_end **************************************************
  write(*,'(a)')"********** subroutine_output_para_wall end"
  write(*,'(a)')""

  return
end subroutine output_para_wall

!*********************************************************************
!*** output_initial_setting_file                                   ***
!*********************************************************************
subroutine output_initial_file         &
&          (idim,nn,nump,ite,time_sim, &
&           type,ifluid,itypep,        &
&           den,vis,grav,delta,sigma,  &
&           time_max,ite_max,          &
&           dt,dtout,cfl,near_wall,    &
&           dis,ker_c,                 &
&           x,dom_size,non_cal,        &
&           grid,grid_n,grid_num,      &
&           dis_rat,c_vel,wall_nor,    &
&           f_type,output_limit,       &
&           nn_grid,neighbor,t,        &
&           dist_min,coll_rat,         &
&           wall_n,wall_anc,grid_dist)
  implicit none
  ! mainroutine_variable ********************************************
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
  integer, intent(in) :: itypep(1:nn)
  real   , intent(in) :: x(0:idim-1,1:nn)
  integer, intent(in) :: grid_num
  integer, intent(in) :: grid_n(0:idim-1)
  real   , intent(in) :: grid(1:grid_num,0:idim-1,0:1)
  real   , intent(in) :: non_cal(0:idim-1)
  real   , intent(in) :: dis_rat
  real   , intent(in) :: c_vel
  integer, intent(in) :: f_type(1:3)
  real   , intent(in) :: ker_c(0:3)
  real   , intent(in) :: dom_size(0:idim-1,0:1)
  integer, intent(in) :: nn_grid
  integer, intent(in) :: neighbor
  integer, intent(in) :: output_limit
  real   , intent(in) :: dist_min
  real   , intent(in) :: coll_rat
  real   , intent(in) :: t(1:nn)
  integer, intent(in) :: wall_n
  real   , intent(in) :: wall_anc(0:idim-1,0:wall_n-1,0:idim-1)
  real   , intent(in) :: grid_dist(1:grid_num,0:4*(idim-1)-1)
  integer, intent(in) :: near_wall(1:grid_num,0:4*(idim-1)-1)
  real   , intent(in) :: wall_nor(0:idim-1,0:wall_n-1)

  ! subroutine_variable *********************************************
  integer             :: n,i,j

  ! calculationstart ************************************************
  write(*,'(a)')"********** subroutine_output_initial_file start"

  ! file_name_*******************************************************
  open(11,file='./data/MPSini/initial_setting.dat')

  ! output_file *****************************************************
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
  write(11,*)0.0,      0.0,      0.0,      0.0
  write(11,*)0.0
  write(11,*)0.0
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

  ! calculation_end *************************************************
  write(*,'(a)')"********** subroutine_output_initial_file end"
  write(*,'(a)')""

  return
end subroutine output_initial_file

end module mod_is
