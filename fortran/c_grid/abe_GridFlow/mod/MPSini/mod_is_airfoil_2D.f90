!*********************************************************************
!****** module_for_initialfile_maker                               ***
!****** ver. 2014.08.13                                            ***
!*********************************************************************
module mod_is_airfoil_dim
  use mod_is_dim
  implicit none
  ! subroutine_definition ********************************************
contains

subroutine change_airfoil_param &
&          (data_n, xdata, ydata, &
&           wscale, xpos, ypos,  &
&           cent,idim, aoa)
  implicit none
  ! mainroutine_variable *********************************************
  integer, intent(in)  :: data_n
  real, intent(inout)  :: xdata(1:data_n)
  real, intent(inout)  :: ydata(1:data_n)
  real, intent(in)     :: xpos
  real, intent(in)     :: ypos
  real, intent(in)     :: wscale
  integer, intent(in)  :: idim
  real, intent(in)     :: aoa
  real, intent(inout)  :: cent(0:idim-1)
  
  ! subroutine_variable **********************************************
  integer         :: i,j,k,n,nn
  real, parameter :: pi = acos(0.0)*2.0
  real            :: tri(0:data_n*4,0:2,0:2)
  real            :: temp_a, temp_b, temp_c
  real            :: a(0:data_n*4)
  real            :: b(0:data_n*4)
  real            :: c(0:data_n*4)
  real            :: sca
  real            :: sx, sy, sz
  real            :: tx, ty, tz
  real            :: c_vector(0:2)
  real            :: inner
  real            :: base_xdata(1:data_n)
  real            :: base_ydata(1:data_n)
  real            :: temp_cent(0:idim-1)
  

  ! calculation_start ************************************************
  
  open(11, file = './data/MPSini/NACA0012.dat')
    do n = 1, data_n
      read(11,*) base_xdata(n), base_ydata(n)
    end do
  close(11)
  
  temp_cent(0) = 0.5*(maxval(base_xdata(1:data_n))+minval(base_xdata(1:data_n)))
  temp_cent(1) = 0.5*(maxval(base_ydata(1:data_n))+minval(base_ydata(1:data_n)))
  
  ! 翼の大きさ変更
  do n = 1, data_n
    base_xdata(n) = base_xdata(n) * wscale
    base_ydata(n) = base_ydata(n) * wscale
  end do
  temp_cent(0) = temp_cent(0) *  wscale
  temp_cent(1) = temp_cent(1) *  wscale
  
  
  ! 迎角の変更
  do n = 2, 66
    xdata(n) = sqrt(base_xdata(n)**2 + base_ydata(n)**2) * &
    &          cos(atan(base_ydata(n)/base_xdata(n)) - aoa/180*pi)
    ydata(n) = sqrt(base_xdata(n)**2 + base_ydata(n)**2) * &
    &          sin(atan(base_ydata(n)/base_xdata(n)) - aoa/180*pi)
  end do
  do n = 68, data_n
    xdata(n) = sqrt(base_xdata(n)**2 + base_ydata(n)**2) * &
    &          cos(atan(base_ydata(n)/base_xdata(n)) - aoa/180*pi)
    ydata(n) = sqrt(base_xdata(n)**2 + base_ydata(n)**2) * &
    &          sin(atan(base_ydata(n)/base_xdata(n)) - aoa/180*pi)
  end do
  cent(0) = sqrt(temp_cent(0)**2 + temp_cent(1)**2) * &
  &         cos(atan(temp_cent(1)/temp_cent(0)) - aoa/180*pi)
  cent(1) = sqrt(temp_cent(0)**2 + temp_cent(1)**2) * &
  &         sin(atan(temp_cent(1)/temp_cent(0)) - aoa/180*pi)
  
  ! 翼の位置の変更
  do n = 1, data_n
    xdata(n) = xdata(n) + xpos
    ydata(n) = ydata(n) + ypos
  end do
  cent(0) = cent(0) + xpos
  cent(1) = cent(1) + ypos
  
  
  return
end subroutine change_airfoil_param



subroutine make_airfoil_stl &
&          (data_n, xdata, ydata, &
&           zmin, zmax, &
&           cent, idim)
  implicit none
  ! mainroutine_variable *********************************************
  integer  :: data_n
  real     :: xdata(1:data_n)
  real     :: ydata(1:data_n)
  real     :: zmin
  real     :: zmax
  integer, intent(in)    :: idim
  real            :: cent(0:idim-1)
  
  ! subroutine_variable **********************************************
  integer         :: i,j,k,n,nn
  real, parameter :: pi = acos(0.0)*2.0
  real            :: temp_cent(0:idim-1)
  real            :: tri(0:data_n*4,0:2,0:2)
  real            :: temp_a, temp_b, temp_c
  real            :: a(0:data_n*4)
  real            :: b(0:data_n*4)
  real            :: c(0:data_n*4)
  real            :: sca
  real            :: sx, sy, sz
  real            :: tx, ty, tz
  real            :: c_vector(0:2)
  real            :: inner
  

  ! calculation_start ************************************************
  
  nn = 0
  
  ! 超ゴミコード注意
  ! 
  ! 手前の面の生成
  
  do n = 1, 65
    tri(nn, 0, 0) = cent(0)
    tri(nn, 0, 1) = cent(1)
    tri(nn, 0, 2) = zmin
    
    tri(nn, 1, 0) = xdata(n)
    tri(nn, 1, 1) = ydata(n)
    tri(nn, 1, 2) = zmin
    
    tri(nn, 2, 0) = xdata(n+1)
    tri(nn, 2, 1) = ydata(n+1)
    tri(nn, 2, 2) = zmin
    
    nn = nn + 1
  end do
  do n = 67, 131
    tri(nn, 0, 0) = cent(0)
    tri(nn, 0, 1) = cent(1)
    tri(nn, 0, 2) = zmin
    
    tri(nn, 1, 0) = xdata(n)
    tri(nn, 1, 1) = ydata(n)
    tri(nn, 1, 2) = zmin
    
    tri(nn, 2, 0) = xdata(n+1)
    tri(nn, 2, 1) = ydata(n+1)
    tri(nn, 2, 2) = zmin
    
    nn = nn + 1
  end do
  
!  ! 奥の面の生成
!  
!  do n = 1, 65
!    tri(nn, 0, 0) = cent(0)
!    tri(nn, 0, 1) = cent(1)
!    tri(nn, 0, 2) = zmax
!    
!    tri(nn, 1, 0) = xdata(n)
!    tri(nn, 1, 1) = ydata(n)
!    tri(nn, 1, 2) = zmax
!    
!    tri(nn, 2, 0) = xdata(n+1)
!    tri(nn, 2, 1) = ydata(n+1)
!    tri(nn, 2, 2) = zmax
!    
!    nn = nn + 1
!  end do
!  do n = 67, 131
!    tri(nn, 0, 0) = cent(0)
!    tri(nn, 0, 1) = cent(1)
!    tri(nn, 0, 2) = zmax
!    
!    tri(nn, 1, 0) = xdata(n)
!    tri(nn, 1, 1) = ydata(n)
!    tri(nn, 1, 2) = zmax
!    
!    tri(nn, 2, 0) = xdata(n+1)
!    tri(nn, 2, 1) = ydata(n+1)
!    tri(nn, 2, 2) = zmax
!    
!    nn = nn + 1
!  end do
!   
!  ! 上面の生成
!  do n = 1, 65
!    tri(nn, 0, 0) = xdata(n)
!    tri(nn, 0, 1) = ydata(n)
!    tri(nn, 0, 2) = zmax
!    
!    tri(nn, 1, 0) = xdata(n)
!    tri(nn, 1, 1) = ydata(n)
!    tri(nn, 1, 2) = zmin
!    
!    tri(nn, 2, 0) = xdata(n+1)
!    tri(nn, 2, 1) = ydata(n+1)
!    tri(nn, 2, 2) = zmin
!    
!    nn = nn + 1
!  end do
!  do n = 1, 65
!    tri(nn, 0, 0) = xdata(n)
!    tri(nn, 0, 1) = ydata(n)
!    tri(nn, 0, 2) = zmax
!    
!    tri(nn, 1, 0) = xdata(n+1)
!    tri(nn, 1, 1) = ydata(n+1)
!    tri(nn, 1, 2) = zmax
!    
!    tri(nn, 2, 0) = xdata(n+1)
!    tri(nn, 2, 1) = ydata(n+1)
!    tri(nn, 2, 2) = zmin
!    
!    nn = nn + 1
!  end do
!  
!  ! 下面の生成
!  do n = 67, 131
!    tri(nn, 0, 0) = xdata(n)
!    tri(nn, 0, 1) = ydata(n)
!    tri(nn, 0, 2) = zmax
!    
!    tri(nn, 1, 0) = xdata(n)
!    tri(nn, 1, 1) = ydata(n)
!    tri(nn, 1, 2) = zmin
!    
!    tri(nn, 2, 0) = xdata(n+1)
!    tri(nn, 2, 1) = ydata(n+1)
!    tri(nn, 2, 2) = zmin
!    
!    nn = nn + 1
!  end do
!  do n = 67, 131
!    tri(nn, 0, 0) = xdata(n)
!    tri(nn, 0, 1) = ydata(n)
!    tri(nn, 0, 2) = zmax
!    
!    tri(nn, 1, 0) = xdata(n+1)
!    tri(nn, 1, 1) = ydata(n+1)
!    tri(nn, 1, 2) = zmax
!    
!    tri(nn, 2, 0) = xdata(n+1)
!    tri(nn, 2, 1) = ydata(n+1)
!    tri(nn, 2, 2) = zmin
!    
!    nn = nn + 1
!  end do
  
  ! ケツが割れているとこの補間
  tri(nn, 0, 0) = cent(0)
  tri(nn, 0, 1) = cent(1)
  tri(nn, 0, 2) = zmin
  
  tri(nn, 1, 0) = xdata(66)
  tri(nn, 1, 1) = ydata(66)
  tri(nn, 1, 2) = zmin
  
  tri(nn, 2, 0) = xdata(132)
  tri(nn, 2, 1) = ydata(132)
  tri(nn, 2, 2) = zmin
    
  nn = nn + 1
  
!  tri(nn, 0, 0) = cent(0)
!  tri(nn, 0, 1) = cent(1)
!  tri(nn, 0, 2) = zmax
!  
!  tri(nn, 1, 0) = xdata(66)
!  tri(nn, 1, 1) = ydata(66)
!  tri(nn, 1, 2) = zmax
!  
!  tri(nn, 2, 0) = xdata(132)
!  tri(nn, 2, 1) = ydata(132)
!  tri(nn, 2, 2) = zmax
!    
!  nn = nn + 1
!  
!  tri(nn, 0, 0) = xdata(132)
!  tri(nn, 0, 1) = ydata(132)
!  tri(nn, 0, 2) = zmin
!  
!  tri(nn, 1, 0) = xdata(66)
!  tri(nn, 1, 1) = ydata(66)
!  tri(nn, 1, 2) = zmin
!  
!  tri(nn, 2, 0) = xdata(66)
!  tri(nn, 2, 1) = ydata(66)
!  tri(nn, 2, 2) = zmax
!  
!  nn = nn + 1
!  
!  tri(nn, 0, 0) = xdata(132)
!  tri(nn, 0, 1) = ydata(132)
!  tri(nn, 0, 2) = zmax
!  
!  tri(nn, 1, 0) = xdata(66)
!  tri(nn, 1, 1) = ydata(66)
!  tri(nn, 1, 2) = zmax
!  
!  tri(nn, 2, 0) = xdata(132)
!  tri(nn, 2, 1) = ydata(132)
!  tri(nn, 2, 2) = zmin
!  
!  nn = nn + 1
  
  do n = 0, nn-1
    sx = tri(n,1,0) - tri(n,0,0)
    sy = tri(n,1,1) - tri(n,0,1)
    sz = tri(n,1,2) - tri(n,0,2)
    
    tx = tri(n,2,0) - tri(n,0,0)
    ty = tri(n,2,1) - tri(n,0,1)
    tz = tri(n,2,2) - tri(n,0,2)
    
    temp_a = sy*tz - sz*ty
    temp_b = sz*tx - sx*tz
    temp_c = sx*ty - sy*tx
    
    sca = sqrt(temp_a**2 + temp_b**2 + temp_c**2)
    
    a(n) = temp_a / sca
    b(n) = temp_b / sca
    c(n) = temp_c / sca
    
    c_vector(0) = cent(0) - tri(n,0,0)
    c_vector(1) = cent(1) - tri(n,0,1)
    c_vector(2) = 0.5*(zmax - zmin) - tri(n,0,2)
    
    inner = a(n)*c_vector(0) + b(n)*c_vector(1) + c(n)*c_vector(2)
    
    if(inner .lt. 0.0)then
      a(n) = - a(n)
      b(n) = - b(n)
      c(n) = - c(n)
    end if
    
  end do
  
  
  open(100,file ='./data/MPSini/NACA0012.stl')
  ! header
  write(100,'(a)')'solid ascii'

  ! data
  do i=0,nn-1
    write(100,'(a,2e13.5,2e13.5,2e13.5)')'  facet normal',a(i),b(i),c(i)
!    write(100,'(a,2e13.5,2e13.5,2e13.5)')'  facet normal',0.0,0.0,1.0
    write(100,'(a)')'    outer loop'
    do  j = 0, 2
        write(100,'(a,2e13.5,2e13.5,2e13.5)')'      vertex',tri(i,j,0),tri(i,j,1),tri(i,j,2)
    end do
    write(100,'(a)')'    endloop'
    write(100,'(a)')'  endfacet'
  end do

  ! footer
  write(100,'(a)')'endsolid'

  close(100)
  
  
  open(11,file = './data/MPSini/normal_vector.vtk')
  write(11,'(a)')'# vtk datafile version 3.0'
  write(11,'(a)')'vtk output'
  write(11,'(a)')'ascii'
  write(11,'(a)')'dataset polydata'
  write(11,'(a,i8,a)')'points',nn*3,' float'

  do n=0,nn-1
    write(11,*)tri(n, 0, 0),tri(n, 0, 1),tri(n, 0, 2)
    write(11,*)tri(n, 1, 0),tri(n, 1, 1),tri(n, 1, 2)
    write(11,*)tri(n, 2, 0),tri(n, 2, 1),tri(n, 2, 2)
  end do

  write(11,'(a,i9)')'point_data',nn*3
  write(11,'(a   )')'field attributes 1'

  write(11,'(a,i9,a)')'normal_vector  3',nn*3,'float'
  do n=0,nn-1
    write(11,*)a(n),b(n),c(n)
    write(11,*)a(n),b(n),c(n)
    write(11,*)a(n),b(n),c(n)
  end do

  close(11)
  
  
  
  return
end subroutine make_airfoil_stl

subroutine make_airfoil_wall &
&          (data_n, xdata, ydata, aoa,idim, &
&           wall_n,wall_anc,wall_out,wall_nor,dis)
  implicit none
  ! mainroutine_variable *********************************************
  integer  :: data_n
  real     :: xdata(1:data_n)
  real     :: ydata(1:data_n)
  integer, intent(in)    :: idim
  integer, intent(in) :: wall_n
  real   , intent(in) :: dis
  real   , intent(inout) :: wall_out(0:idim-1,0:wall_n-1)
  real   , intent(inout) :: wall_nor(0:idim-1,0:wall_n-1)
  real   , intent(inout) :: wall_anc(0:idim-1,0:wall_n-1,0:1)
  real     :: aoa
  
  
  ! subroutine_variable **********************************************
  integer         :: i,j,k,n,nn
  real, parameter :: pi = acos(0.0)*2.0
  real            :: cent(0:idim-1)
  real            :: temp_cent(0:idim-1)
  real            :: tri(0:data_n*4,0:2,0:2)
  real            :: temp_a, temp_b, temp_c
  real            :: a(0:data_n*4)
  real            :: b(0:data_n*4)
  real            :: c(0:data_n*4)
  real            :: sca
  real            :: sx, sy, sz
  real            :: tx, ty, tz
  real            :: c_vector(0:2)
  real            :: inner
  

  ! calculation_start ************************************************
  
  ! 壁面境界条件の生成
  ! 上面
  do n = 1, 65
    call wall_setting        &
    &    (n-1, xdata(n), ydata(n), xdata(n+1), ydata(n+1), &
    &     idim,wall_n,wall_anc,wall_out,wall_nor,dis)
  end do
  ! 下面
  do n = 67, 131
    call wall_setting        &
    &    (n-2, xdata(n+1), ydata(n+1), xdata(n), ydata(n), &
    &     idim,wall_n,wall_anc,wall_out,wall_nor,dis)
  end do
  ! ケツ
  call wall_setting        &
  &    (130, xdata(66), ydata(66), xdata(132), ydata(132), &
  &     idim,wall_n,wall_anc,wall_out,wall_nor,dis)
  
  
  return
end subroutine make_airfoil_wall

end module mod_is_airfoil_dim
