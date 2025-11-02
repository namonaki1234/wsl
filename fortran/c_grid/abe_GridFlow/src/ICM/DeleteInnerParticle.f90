!粒子外挿コード
!"ifort DeleteInnerParticle.f90 -heap-arrays -qopenmp -O3"でコンパイル
program DeleteInnerParticle
implicit none
!variable
integer,parameter	:: nn = 50000000
integer,parameter	:: droplet_nn = 100000
real,allocatable	:: x(:,:)
real,allocatable	:: v(:,:)
real,allocatable	:: t(:)
integer,allocatable	:: itypep(:)
real,allocatable	:: gamma(:)
integer,allocatable	:: swi_koeki(:)
integer,allocatable	:: cal_type0(:)
integer,allocatable	:: droplet_num(:)
real,allocatable	:: droplet(:,:)

real,allocatable	:: x_type2(:,:)
real,allocatable	:: t_type2(:)
integer			:: nump_type4
integer			:: nump_type0

real,allocatable	:: dom(:,:)
real,allocatable	:: copy(:,:)
real,allocatable	:: x_over(:,:)

real,allocatable	:: x_set(:,:)

real,allocatable	:: an_set(:)
real,allocatable	:: wall(:,:)
real,allocatable	:: dist_w_set(:)
integer,allocatable	:: boundary_set(:)

real,parameter		:: pi = acos(-1.0)
real,parameter		:: zero = 1.0e-8
integer,parameter	:: idim = 2
real			:: chord
real			:: AOA
real			:: dis
real			:: mvd
real			:: dtout

real			:: dom_size
real			:: dom_xmin,dom_xmax,dom_ymin,dom_ymax
real			:: dom_j
integer			:: nump
integer			:: numw
integer			:: nump_preIP	!前計算で削除した内部粒子の総数
integer			:: nump_pre
integer			:: nump_dom	!計算領域の粒子数
integer			:: nump_over	!mps粒子と重なる配列粒子
integer			:: nump_set
integer			:: nump_restart
integer			:: nump_restart_set
integer			:: numdx,numdy
real			:: rr,rr_jdg
real			:: rr_min
real			:: an_jdg
real			:: inter_rat

character		:: char*200
real			:: dummy
integer			:: i,j,k,n
integer			:: ios
integer,allocatable	:: flag_over(:)
integer			:: flag_dom
integer,allocatable	:: flag_restart_set(:)
integer,allocatable	:: flag_restart(:)
integer,allocatable	:: flag_temp_set(:)
integer,allocatable	:: flag_temp(:)

integer			:: swi_restart

real			:: time_sim
real			:: nextin
integer			:: incount
integer			:: fcount
integer			:: count_remesh
double precision	:: time_sim_dble

integer			:: swi_dom
integer			:: swi_over

integer			:: step

!openmp
!integer		:: parallel_threads = 4
integer		:: start_time_c,end_time_c,CountPerSec,CountMax
real		:: process_time
!call omp_set_num_threads(parallel_threads)

!Calculation_Time
call system_clock(start_time_c,CountPerSec,CountMax)

!array_allocation
allocate(x(0:nn-1,0:2))
allocate(v(0:nn-1,0:2))
allocate(t(0:nn-1))
allocate(itypep(0:nn-1))
allocate(gamma(0:nn-1))
allocate(swi_koeki(0:nn-1))
allocate(cal_type0(0:nn-1))
allocate(droplet_num(0:nn-1))
allocate(droplet(0:droplet_nn-1,0:6))
x = 0
v = 0
t = 273.15
itypep = -1
gamma = -1.0
swi_koeki = -1
cal_type0 = -1
droplet_num = -1
droplet = -1.0

allocate(x_type2(0:nn-1,0:2))
allocate(t_type2(0:nn-1))

allocate(copy(0:nn-1,0:2))
allocate(dom(0:nn-1,0:2))
allocate(x_over(0:nn-1,0:2))

!リスタートファイル用スイッチ
swi_restart = 0 !0：リスタートファイルから読み込み（最初）
                !1：リスタートファイルから読み込み
swi_dom = 0	!0：計算領域を設定
		!1：テキストファイルから読み込み
swi_over = 0	!0：MPS粒子と重なる粒子を計算
		!1：テキストファイルから読み込み

open(1,file = './data/STEP.txt',status = 'old')
 read(1,*) step
close(1)
if(step .eq. 1) then
 swi_restart = 0
else
 swi_restart = 1
end if

if(swi_restart .ne. 0 .and. swi_restart .ne. 1) then
 write(*,*) 'error : swi_restart'
 stop
end if

!setting_parameter
an_jdg = 9.0	!粒子数密度の閾値

chord = 0.53 !0.533 !0.267
AOA = 4.0
dis = 100.0e-6 !2.0e-5 !1.35e-5
mvd = dis * 10.0
dtout = 0.01
inter_rat = 0.9

dom_xmin = -0.05 * chord !-0.015
dom_xmax = 0.15 * chord !0.04
dom_ymin = -0.075 * chord !-0.02
dom_ymax = 0.075 * chord !0.02
dom_j = abs(dom_xmin) * 2.0 / chord + 1.0

!calculation
write(*,*) '**************************************'
write(*,*) 'Computation Start'
write(*,*) '**************************************'

!input_mps_result
select case(swi_restart)
 case(0)
  nump_preIP = 0
 case(1)
  open(1,file = './result/ICM/preInnerParticle.dat',iostat = ios,status = 'old')
   i = 0
   do
    read(1,*,iostat = ios) x_type2(i,0),x_type2(i,1),x_type2(i,2)
    if(ios .ne. 0) exit
    i = i + 1
   end do
   nump_preIP = i
  close(1)
end select
call Input_restartFile(nn,droplet_nn,idim,x,v,t,itypep,gamma, &
     &                 swi_koeki,cal_type0,droplet_num,droplet, &
     &                 nump_preIP,nump_pre,time_sim,nextin,incount,fcount,count_remesh,time_sim_dble)


!remove_not_type2
nump_type4 = 0
nump_type0 = 0
j = nump_preIP
do i=0,nump_pre-1
 if(itypep(i) .eq. 2) then
  x_type2(j,0) = x(i,0)
  x_type2(j,1) = x(i,1)
  t_type2(j) = t(i)
  j = j + 1
 else if(itypep(i) .eq. 4) then
  nump_type4 = nump_type4 + 1
 else if(itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1) then
  nump_type0 = nump_type0 + 1
 end if
end do
nump = j

!array_copy
do i = 0,nump-1
 copy(i,0) = x_type2(i,0)
 copy(i,1) = x_type2(i,1)
 copy(i,2) = x_type2(i,2)
end do
deallocate(x_type2)
allocate(x_type2(0:nump-1,0:2))
do i = 0,nump-1
 x_type2(i,0) = copy(i,0)
 x_type2(i,1) = copy(i,1)
 x_type2(i,2) = copy(i,2)
end do

write(*,*) nump_preIP,nump_pre,nump_type4,nump_type0,nump

write(*,*) '--------------------------------------'
write(*,*) 'Input info'
write(*,*) '  number of type2 : ',nump - nump_preIP
write(*,*) '  number of type4 : ',nump_type4
write(*,*) '  number of type0 : ',nump_type0
write(*,*) '--------------------------------------'

!input_wall.vtk
open(1,file = './data/MPSini/wall.vtk',status = 'old')
 do i=0,3
  read(1,'(A)') char
 end do
 read(1,'(A,I3,A)') char(1:11),numw
 allocate(wall(0:numw-1,0:2))
 do i=0,numw-1
  read(1,*) wall(i,0),wall(i,1),wall(i,2)
 end do
close(1)

!setting_area
write(*,*) '--------------------------------------'
write(*,*) 'start setting computation domain'
select case(swi_dom)
 case(0)
  numdx = int((dom_xmax - dom_xmin) / (dis * inter_rat))
  numdy = int((dom_ymax - dom_ymin) / (dis * inter_rat))
  k = 0
  do i = 0,numdx-1
   do j = 0,numdy-1
    dom(k,0) = real(i) * dis * inter_rat + dom_xmin
    dom(k,1) = real(j) * dis * inter_rat + dom_ymin
    call judge_dom(dom(k,0),dom(k,1),numw,wall,dom_j,dis,chord,flag_dom)
    if(flag_dom .eq. 1) then
     k = k + 1
    else
     cycle
    end if
   end do
  end do
  nump_dom = k
  !output_data
  open(1,file = './data/ICM/domain.dat',status = 'replace')
   write(1,*) nump_dom
   do i = 0,nump_dom-1
    write(1,*) dom(i,0),dom(i,1)
   end do
  close(1)
  open(1,file = './data/ICM/domain.vtk',status = 'replace')
   write(1,'(A)') '# vtk DataFile Version 3.0'
   write(1,*) 'vtk output'
   write(1,*) 'ASCII'
   write(1,*) 'DATASET POLYDATA'
   write(1,*) 'POINTS ',nump_dom,' float'
   do i=0,nump_dom-1
    write(1,*) dom(i,0),dom(i,1),0.0
   end do
   write(1,*) 'POINT_DATA ',nump_dom
   write(1,*) 'SCALARS particle_type float 1'
   write(1,*) 'LOOKUP_TABLE default'
   do i=0,nump_dom-1
    write(1,*) 1
   end do
  close(1)
 case(1)
  open(1,file = './data/ICM/domain.dat',status = 'old')
   read(1,*) nump_dom
   do i = 0,nump_dom-1
    read(1,*) dom(i,0),dom(i,1)
   end do
  close(1)
end select
write(*,*) 'complete setting computation domain'

write(*,*) '--------------------------------------'
write(*,*) 'start setting overlap particles'
allocate(flag_over(0:nump_dom-1))
flag_over = 1

n = 0
if(swi_over .eq. 0) then
 k = 0
 !$omp parallel do default(shared) private(i)
 do i = 0,nump_dom-1
  k = k + 1
  if(mod(k,int(real(nump_dom-1)/100.0)) .eq. 0.0) write(*,*) real(k)/real(nump_dom)*100.0, '%'
  do j = 0,nump-1
   rr = sqrt((dom(i,0) - x_type2(j,0))**2.0 + (dom(i,1) - x_type2(j,1))**2.0)
   if(rr .le. dis*0.75) then
    flag_over(i) = 0
    exit
   end if
  end do
 end do
 !$omp end parallel do
 open(1,file = './data/ICM/flag_over.dat',status = 'replace')
  do i = 0,nump_dom-1
   write(1,*) flag_over(i)
  end do
 close(1)
end if
open(1,file = './data/ICM/flag_over.dat',status = 'old')
 do i = 0,nump_dom-1
  read(1,*) flag_over(i)
  !setting_overlap_particles
  if(flag_over(i) .eq. 0) then
   x_over(n,0) = dom(i,0)
   x_over(n,1) = dom(i,1)
   n = n + 1
  end if
 end do
close(1)
nump_over = n
write(*,*) ' over particle :',nump_over

open(1,file = './data/ICM/domain.vtk',status = 'replace')
 write(1,'(A)') '# vtk DataFile Version 3.0'
 write(1,*) 'vtk output'
 write(1,*) 'ASCII'
 write(1,*) 'DATASET POLYDATA'
 write(1,*) 'POINTS ',nump_dom,' float'
 do i=0,nump_dom-1
  write(1,*) dom(i,0),dom(i,1),0.0
 end do
 write(1,*) 'POINT_DATA ',nump_dom
 write(1,*) 'SCALARS particle_type float 1'
 write(1,*) 'LOOKUP_TABLE default'
 do i=0,nump_dom-1
  write(1,*) flag_over(i)
 end do
close(1)

write(*,*) 'complete setting overlap particles'

copy = 0.0
n = 0
do i = 0,nump_dom-1
 if(flag_over(i) .eq. 1) then
  copy(n,0) = dom(i,0)
  copy(n,1) = dom(i,1)
  n = n + 1
 end if
end do
nump_dom = n !nump_dom更新
dom = 0.0
do i = 0,nump_dom-1
 dom(i,0) = copy(i,0)
 dom(i,1) = copy(i,1)
end do

!calculate_icing_cell
write(*,*) '--------------------------------------'
write(*,*) 'start comute icing cell'
call read(nump,AOA,dis,x_type2)!,itypep)
write(*,*) 'complete compute icing cell'

!calculation_number_density_for_setting_particle
write(*,*) '--------------------------------------'
write(*,*) 'start compute setting particle'
nump_set = nump_over
allocate(x_set(0:nump_set-1,0:2))
allocate(an_set(0:nump_set-1))
allocate(dist_w_set(0:nump_set-1))
an_set = 0.0
j = 0
do i=0,nump_over-1
 x_set(j,0) = x_over(i,0)
 x_set(j,1) = x_over(i,1)
 j = j + 1
end do
call wall_dist(nump_set,numw,x_set,wall,dist_w_set)
call cal_an(nump_set,x_set,an_set,dist_w_set,dis,inter_rat, &
     &      dom_xmin,dom_xmax,dom_ymin,dom_ymax)

!output_setting_prticle_data
open(1,file = './result/ICM/ptcle_set.vtk',status = 'replace')
 write(1,'(A)') '# vtk DataFile Version 3.0'
 write(1,*) 'vtk output'
 write(1,*) 'ASCII'
 write(1,*) 'DATASET POLYDATA'
 write(1,*) 'POINTS ',nump_set,' float'
 do i=0,nump_set-1
  write(1,*) x_set(i,0),x_set(i,1),0.0
 end do
 write(1,*) 'POINT_DATA ',nump_set
 write(1,*) 'SCALARS an float 1'
 write(1,*) 'LOOKUP_TABLE default'
 do i=0,nump_set-1
  write(1,*) an_set(i)
 end do
close(1)
write(*,*) 'complete compute setting particle'

!judge_boundary
write(*,*) '--------------------------------------'
write(*,*) 'start judge boundary'
allocate(boundary_set(0:nump_set-1))
boundary_set = 0
do i = 0,nump_set-1
 if((an_set(i) .le. an_jdg) .and. (dist_w_set(i) .gt. dis * 0.5)) then
  boundary_set(i) = 1
 end if
end do
write(*,*) 'complete judge boundary'

!making_restart_vtk
allocate(flag_restart_set(0:nump_set-1))
allocate(flag_temp_set(0:nump_set-1))
flag_temp_set = 0
flag_restart_set = 0
do i=0,nump_set-1
 if(boundary_set(i) .eq. 0) cycle
 flag_restart_set(i) = 1
 flag_temp_set(i) = 1
 do j=0,nump_set-1
  if(j .eq. i) cycle
  if(boundary_set(j) .eq. 1) cycle
  rr = sqrt((x_set(j,0) - x_set(i,0))**2.0 + (x_set(j,1) - x_set(i,1))**2.0)
  if(rr .lt. 4.1 * dis) then
   flag_restart_set(j) = 1
  end if
  if(rr .lt. 1.1 * dis) then
   flag_temp_set(j) = 1
  end if
 end do
end do
nump_restart_set = 0
do i = 0,nump_set-1
 if(flag_restart_set(i) .eq. 1) then
  nump_restart_set = nump_restart_set + 1
 end if
end do

allocate(flag_restart(0:nump-1))
allocate(flag_temp(0:nump-1))
flag_restart = 0
flag_temp = 0
!$omp parallel do default(shared) private(i)
do i = 0,nump-1
 do j = 0,nump_set-1
  if(flag_restart_set(j) .eq. 0) cycle
  rr = sqrt((x_set(j,0) - x_type2(i,0))**2.0 + (x_set(j,1) - x_type2(i,1))**2.0)
  if(rr .lt. dis * 0.75) then
   flag_restart(i) = 1
   if(flag_temp_set(j) .eq. 1) then
    flag_temp(i) = 1
   end if
  end if
 end do
end do
!$omp end parallel do
nump_restart = 0
do i = 0,nump-1
 if(flag_restart(i) .eq. 1) nump_restart = nump_restart + 1
end do

open(1,file = './result/ICM/restart_particle.vtk',status = 'replace')
 write(1,'(A)') '# vtk DataFile Version 3.0'
 write(1,*) 'vtk output'
 write(1,*) 'ASCII'
 write(1,*) 'DATASET POLYDATA'
 write(1,*) 'POINTS ',nump_restart+nump_type4+nump_type0,' float'
 do i=0,nump-1
  if(flag_restart(i) .eq. 0) cycle
  write(1,*) x_type2(i,0),x_type2(i,1),x_type2(i,2)
 end do
 do i=0,nump_pre-1
  if(itypep(i) .eq. 4 .or. (itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1)) then
   write(1,*) x(i,0),x(i,1),x(i,2)
  end if
 end do
 write(1,*) 'POINT_DATA ',nump_restart+nump_type4+nump_type0
 write(1,*) 'SCALARS temperature float 1'
 write(1,*) 'LOOKUP_TABLE default'
 do i=0,nump-1
  if(flag_restart(i) .eq. 0) cycle
  write(1,*) t_type2(i)
 end do
 do i=0,nump_pre-1
  if(itypep(i) .eq. 4 .or. (itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1)) then
   write(1,*) t(i)
  end if
 end do
 write(1,*) 'SCALARS flag_temperature float 1'
 write(1,*) 'LOOKUP_TABLE default'
 do i=0,nump-1
  if(flag_restart(i) .eq. 0) cycle
  write(1,*) flag_temp(i)
 end do
 do i=0,nump_pre-1
  if(itypep(i) .eq. 4 .or. (itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1)) then
   write(1,*) 1
  end if
 end do
 write(1,*) 'SCALARS itypep float 1'
 write(1,*) 'LOOKUP_TABLE default'
 do i=0,nump-1
  if(flag_restart(i) .eq. 0) cycle
  write(1,*) 2
 end do
 do i=0,nump_pre-1
  if(itypep(i) .eq. 4 .or. (itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1)) then
   write(1,*) itypep(i)
  end if
 end do
close(1)

!output_restart_dat(for_coupling_code)
open(1,file = './result/ICM/restart_file.dat',status = 'replace')
 write(1,*) nump_restart + nump_type4 + nump_type0
 write(1,*) dis
 write(1,*) mvd
 write(1,*) time_sim
 write(1,*) nextin
 write(1,*) incount
 write(1,*) fcount
 write(1,*) count_remesh
 write(1,*) time_sim_dble
 !coordinate
 k = 0
 do i=0,nump-1
  if(flag_restart(i) .eq. 0) cycle
  do j = 0,idim-1
   write(1,*) x_type2(i,j)
  end do
  k = k + 1
 end do
 if(k .ne. nump_restart) then
  write(*,*) 'error:coord,type2'
  stop
 end if
 k = 0
 do i=0,nump_pre-1
  if(itypep(i) .eq. 4 .or. (itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1)) then
   do j = 0,idim-1
    write(1,*) x(i,j)
   end do
   k = k + 1
  end if
 end do
 if(k .ne. nump_type4+nump_type0) then
  write(*,*) 'error:coord,type4and0'
  stop
 end if
 !velocity
 do i=0,nump-1
  if(flag_restart(i) .eq. 0) cycle
  do j = 0,idim-1
   write(1,*) 0.0
  end do
 end do
 do i=0,nump_pre-1
  if(itypep(i) .eq. 4 .or. (itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1)) then
   do j = 0,idim-1
    write(1,*) v(i,j)
   end do
  end if
 end do
 !temperature
 do i=0,nump-1
  if(flag_restart(i) .eq. 0) cycle
  write(1,*) t_type2(i)
 end do
 do i=0,nump_pre-1
  if(itypep(i) .eq. 4 .or. (itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1)) then
   write(1,*) t(i)
  end if
 end do
 !itypep
 do i=0,nump-1
  if(flag_restart(i) .eq. 0) cycle
  write(1,*) 2
 end do
 do i=0,nump_pre-1
  if(itypep(i) .eq. 4 .or. (itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1)) then
   write(1,*) itypep(i)
  end if
 end do
 !flag_temp
 do i=0,nump-1
  if(flag_restart(i) .eq. 0) cycle
  write(1,*) flag_temp(i)
 end do
 do i=0,nump_pre-1
  if(itypep(i) .eq. 4 .or. (itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1)) then
   write(1,*) 1
  end if
 end do
 !gamma
 do i=0,nump-1
  if(flag_restart(i) .eq. 0) cycle
  write(1,*) 1.0
 end do
 do i=0,nump_pre-1
  if(itypep(i) .eq. 4 .or. (itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1)) then
   write(1,*) gamma(i)
  end if
 end do
 !swi_koeki
 do i=0,nump-1
  if(flag_restart(i) .eq. 0) cycle
  write(1,*) 0
 end do
 do i=0,nump_pre-1
  if(itypep(i) .eq. 4 .or. (itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1)) then
   write(1,*) swi_koeki(i)
  end if
 end do
 !swi_cal_type0
 do i=0,nump-1
  if(flag_restart(i) .eq. 0) cycle
  write(1,*) 0
 end do
 do i=0,nump_pre-1
  if(itypep(i) .eq. 4 .or. (itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1)) then
   write(1,*) cal_type0(i)
  end if
 end do
 !droplet_num
 do i=0,nump-1
  if(flag_restart(i) .eq. 0) cycle
  write(1,*) 0
 end do
 do i=0,nump_pre-1
  if(itypep(i) .eq. 4 .or. (itypep(i) .eq. 0 .and. cal_type0(i) .eq. 1)) then
   write(1,*) droplet_num(i)
  end if
 end do
 !droplet_swi
 do i=0,incount
  write(1,*) droplet(i,0)
 end do
 !droplet
 do i=0,incount
  write(1,*) droplet(i,1)
  write(1,*) droplet(i,2)
  write(1,*) droplet(i,3)
  write(1,*) droplet(i,4)
  write(1,*) droplet(i,5)
  write(1,*) droplet(i,6)
 end do
close(1)
write(*,*) 'complete output file'

open(1,file = './result/ICM/preInnerParticle.dat',status = 'replace')
 do i=0,nump-1
  if(flag_restart(i) .eq. 1) cycle
  write(1,*) x_type2(i,0),x_type2(i,1),x_type2(i,2)
 end do
close(1)

call system_clock(end_time_c)
process_time = real((end_time_c - start_time_c) / CountPerSec)

write(*,*) '**************************************'
write(*,*) 'Computation End'
write(*,*) 'TIME : ',int(process_time / 60.0),'min',mod(process_time,60.0),'sec'
write(*,*) '**************************************'

contains
!**************************************************************************
subroutine Input_restartFile( &
& nn,droplet_nn,idim,x,v,t,itypep,gamma, &
& swi_koeki,cal_type0,droplet_num,droplet, &
& nump_pre,nump,time_sim,nextin,incount,fcount,count_remesh,time_sim_dble)
 !mainroutine_variable
 integer,intent(in)	:: nn
 integer,intent(in)	:: droplet_nn
 integer,intent(in)	:: idim
 real,intent(inout)	:: x(0:nn-1,0:2)
 real,intent(inout)	:: v(0:nn-1,0:2)
 real,intent(inout)	:: t(0:nn-1)
 integer,intent(inout)	:: itypep(0:nn-1)
 real,intent(inout)	:: gamma(0:nn-1)
 integer,intent(inout)	:: swi_koeki(0:nn-1)
 integer,intent(inout)	:: cal_type0(0:nn-1)
 integer,intent(inout)	:: droplet_num(0:nn-1)
 real,intent(inout)	:: droplet(0:droplet_nn-1,0:6)
 integer,intent(in)	:: nump_pre
 integer,intent(inout)	:: nump
 real,intent(out)	:: time_sim
 real,intent(out)	:: nextin
 integer,intent(out)	:: incount
 integer,intent(out)	:: fcount
 integer,intent(out)	:: count_remesh
 double precision,intent(out)	:: time_sim_dble

 !subroutine_variable
 integer	:: i,j
 real		:: dis
 real		:: mvd
 real	:: dummy_real
 integer		:: dummy_int

 open(1,file = './result/MPS/restart_file.dat',status = 'old')
  read(1,*) nump
  read(1,*) dis
  read(1,*) mvd
  read(1,*) time_sim
  read(1,*) nextin
  read(1,*) incount
  read(1,*) fcount
  read(1,*) count_remesh
  read(1,*) time_sim_dble

  nump = nump_pre + nump
  !coordinate
  do i = nump_pre, nump-1
   do j = 0,idim-1
    read(1,*) x(i,j)
   end do
  end do
  !velocity
  do i = nump_pre, nump-1
   do j = 0,idim-1
    read(1,*) v(i,j)
   end do
  end do
  !temperature
  do i = nump_pre,nump-1
   read(1,*) t(i)
  end do
  !itypep
  do i = nump_pre,nump-1
   read(1,*) itypep(i)
  end do
  !gamma
  do i = nump_pre,nump-1
   read(1,*) gamma(i)
  end do
  !swi_koeki
  do i = nump_pre,nump-1
   read(1,*) swi_koeki(i)
  end do
  !cal_type0
  do i = nump_pre,nump-1
   read(1,*) cal_type0(i)
  end do
  !sroplet_num
  do i = nump_pre,nump-1
   read(1,*) droplet_num(i)
  end do
  !droplet_swi
  do i = 0,incount
   read(1,*) droplet(i,0)
  end do
  !droplet_coord_and_velocity
  do i = 0,incount
   read(1,*) droplet(i,1)
   read(1,*) droplet(i,2)
   read(1,*) droplet(i,3)
   read(1,*) droplet(i,4)
   read(1,*) droplet(i,5)
   read(1,*) droplet(i,6)
  end do
 close(1)
end subroutine Input_restartFile
!**************************************************************************
subroutine Median(nn,nump_pre,nump,neigh,neigh_med)
 implicit none
 !mainroutine_variable
 integer,intent(in)	:: nn
 integer,intent(in)	:: nump_pre
 integer,intent(in)	:: nump
 integer,intent(in)	:: neigh(0:nn-1)
 real,intent(inout)	:: neigh_med
 !subroutine_variable
 real,allocatable	:: array(:)
 real			:: a
 integer		:: i,j

 !array_allocation
 allocate(array(nump_pre:nump-1))

 do i=nump_pre,nump-1
  array(i) = neigh(i)
 end do

 do i=nump_pre,nump-1
  do j=nump_pre,nump-1
   if(array(i) .lt. array(j)) then
    a = array(i)
    array(i) = array(j)
    array(j) = a
   end if
  end do
 end do

 if(mod((nump-nump_pre),2) .eq. 0) then
  neigh_med = (array((nump - nump_pre) / 2 - 1 + nump_pre) + array((nump - nump_pre) / 2 + nump_pre)) / 2.0
 else
  neigh_med = array(int((nump - nump_pre) / 2) + nump_pre)
 end if
end subroutine Median
!**************************************************************************
subroutine judge_dom(x,y,numw,wall,dom_j,dis,chord,flag_dom)
 implicit none
 !mainroutine_variable
 real,intent(in)	:: x,y
 integer,intent(in)	:: numw
 real,intent(in)	:: wall(0:numw-1,0:2)
 real,intent(in)	:: dom_j
 real,intent(in)	:: dis
 real,intent(in)	:: chord
 integer,intent(inout)	:: flag_dom

 !subroutine_variable
 integer	:: i
 integer,dimension(0:1)	:: flag
 real		:: y_wall
 real,dimension(0:numw-1,0:2)	:: outside
 real		:: y_out

 flag = 0

 do i = 0,numw-1
  outside(i,0) = dom_j * (wall(i,0) - chord * 0.5) + chord * 0.5
  outside(i,1) = dom_j * wall(i,1) * 2.0
 end do

 !wall
 if(x .ge. 0.0) then
  do i = 0,numw-2 !iとi+1の区間なのでnum-2まで
   if(i .le. 129) then
    if(mod(i,2) .eq. 1) cycle
    if(x .ge. wall(i,0) .and. x .le. wall(i+1,0)) then
     y_wall = (wall(i+1,1) - wall(i,1)) / (wall(i+1,0) - wall(i,0)) &
            &  * (x - wall(i,0)) + wall(i,1)
     if(y .gt. y_wall + dis * 0.5) then
      flag(0) = 1
      exit
     end if
    end if
   else
    if(mod(i,2) .eq. 0) cycle
    if(x .ge. wall(i,0) .and. x .le. wall(i-1,0)) then
     y_wall = (wall(i-1,1) - wall(i,1)) / (wall(i-1,0) - wall(i,0)) &
            & * (x - wall(i,0)) + wall(i,1)
     if(y .lt. y_wall - dis * 0.5) then
      flag(0) = 1
      exit
     end if
    end if
   end if
  end do
 else
  flag(0) = 1
 end if

 !outside
 do i = 0,numw-2 !iとi+1の区間なのでnum-2まで
  if(i .le. 129) then
   if(mod(i,2) .eq. 1) cycle
   if(x .ge. outside(i,0) .and. x .le. outside(i+1,0)) then
    y_out = (outside(i+1,1) - outside(i,1)) / (outside(i+1,0) - outside(i,0)) &
          & * (x - outside(i,0)) + outside(i,1)
    if(y .lt. y_out .and. y .ge. 0.0) then
     flag(1) = 1
     exit
    end if
   end if
  else
   if(mod(i,2) .eq. 0) cycle
   if(x .ge. outside(i,0) .and. x .le. outside(i-1,0)) then
    y_out = (outside(i-1,1) - outside(i,1)) / (outside(i-1,0) - outside(i,0)) &
          & * (x - outside(i,0)) + outside(i,1)
    if(y .gt. y_out .and. y .le. 0.0) then
     flag(1) = 1
     exit
    end if
   end if
  end if
 end do

 if(flag(0) .eq. 1 .and. flag(1) .eq. 1) then
  flag_dom = 1
 else
  flag_dom = 0
 end if

end subroutine judge_dom
!**************************************************************************
subroutine wall_dist(nump,numw,x,wall,dist_w)
 !mainroutine_variable
 integer,intent(in)	:: nump,numw
 real,intent(in)	:: x(0:nump-1,0:2)
 real,intent(in)	:: wall(0:numw-1,0:2)
 real,intent(inout)	:: dist_w(0:nump-1)
 !subroutine_variable
 integer			:: i,j,k
 real,dimension(0:numw-1)	:: dist
 integer,dimension(0:1)		:: min_wall
 real				:: a,b,c
 real				:: sum
 real				:: S

 do i=0,nump-1
  !find_nearest_wall_point
  a = 10.0
  do j=0,numw-1
   dist(j) = sqrt((x(i,0) - wall(j,0))**2.0 + (x(i,1) - wall(j,1))**2.0)
   if(dist(j) .lt. a) then
    a = dist(j)
    min_wall(0) = j
   end if
  end do
  !find_second_nearest_wall_point
  b = 10.0
  do j = 0,numw-1
   if(wall(j,0) .eq. wall(min_wall(0),0)  &
   & .and. wall(j,1) .eq. wall(min_wall(0),1)) cycle
   if(dist(j) .lt. b) then
    b = dist(j)
    min_wall(1) = j
   end if
  end do
  c = sqrt((wall(min_wall(0),0) - wall(min_wall(1),0))**2.0  &
      &  + (wall(min_wall(0),1) - wall(min_wall(1),1))**2.0)
  sum = (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c)
  if(sum .le. 0.0) then
   S = 0.0
   dist_w(i) = 0.0
  else
   S = sqrt(sum) / 4.0
   dist_w(i) = 2.0 * S / c
  end if
 end do
end subroutine wall_dist
!**************************************************************************
subroutine cal_an(nump,x,an,dist_w,dis,inter_rat,dom_xmin,dom_xmax,dom_ymin,dom_ymax)
 !mainroutine_valiable
 integer,intent(in)	:: nump
 real,intent(in)	:: x(0:nump-1,0:2)
 real,intent(inout)	:: an(0:nump-1)
 real,intent(in)	:: dis
 real,intent(in)	:: dist_w(0:nump-1)
 real,intent(in)	:: inter_rat
 real,intent(in)	:: dom_xmin,dom_xmax,dom_ymin,dom_ymax
 !subroutine_variable
 integer		:: i,j,k,l
 real			:: dist
 real			:: re
 real,dimension(0:50000-1)	:: an_wall
 integer,allocatable	:: backet(:,:,:)
 integer		:: x_num,y_num
 integer,dimension(0:nump-1,0:1)	:: backet_num
 integer		:: ii,jj
 real,dimension(0:1)	:: x_neigh

 re = 3.1 * dis

 x_num = int((dom_xmax - dom_xmin) / re)
 y_num = int((dom_ymax - dom_ymin) / re)
 allocate(backet(0:x_num-1,0:y_num-1,-1:100))
 !第三成分(-1:バケット内粒子数,i:i番目の粒子番号)
 backet = -1
 do i = 0,nump-1
  ii = int((x(i,0) - dom_xmin) / re) - 1
  jj = int((x(i,1) - dom_ymin) / re) - 1
  backet_num(i,0) = ii
  backet_num(i,1) = jj

  if(ii .gt. x_num) then
   write(*,*) 'error: particle is not existed in backet',ii,x_num,x(i,0)
   write(*,*) i,x(i,0),x(i,1)
   write(*,*) 3.1 * dis
   write(*,*) dom_xmax - dom_xmin
   stop
  end if

  if(backet(ii,jj,-1) .eq. -1) then
   backet(ii,jj,-1) = 1
   backet(ii,jj,backet(ii,jj,-1)-1) = i
  else
   backet(ii,jj,-1) = backet(ii,jj,-1) + 1
   backet(ii,jj,backet(ii,jj,-1)-1) = i
  end if
 end do

 !set_wall_number_density
 an_wall = 0.0
 do l = 0,31
  do i = -3,3
   do j = -3,0
    dist = sqrt((real(i) * dis * inter_rat)**2.0 &
           & + ((real(l) * 0.1 - real(j) * inter_rat + 0.5) * dis)**2.0)
    an_wall(l) = an_wall(l) + w(dist,re)
   end do
  end do
 end do

 !calculate_number_density
 k = 0
 !$omp parallel do default(shared) private(i)
 do i=0,nump-1
  an(i) = an_wall(nint(dist_w(i) / dis * 10.0)) * 3.0 !壁の重みを大きくするため3倍
  k = k + 1
!  if(mod(k,int(real(nump)*0.01)) .eq. 0.0) write(*,*) k,real(k)/real(nump)*100.0
  do ii = -1,1
   do jj = -1,1
    if(backet(backet_num(i,0)+ii,backet_num(i,1)+jj,-1) .eq. -1) cycle
    do j = 0,backet(backet_num(i,0)+ii,backet_num(i,1)+jj,-1)-1
     if(backet(backet_num(i,0)+ii,backet_num(i,1)+jj,j) .eq. i) cycle
     x_neigh(0) = x(backet(backet_num(i,0)+ii,backet_num(i,1)+jj,j),0)
     x_neigh(1) = x(backet(backet_num(i,0)+ii,backet_num(i,1)+jj,j),1)
     dist = sqrt((x(i,0) - x_neigh(0))**2.0 + (x(i,1) - x_neigh(1))**2.0)
     an(i) = an(i) + w(dist,re)
    end do
   end do
  end do
 end do
 !$omp end parallel do

 deallocate(backet)

end subroutine cal_an

real function w(dist,R)
 real dist,R
 if(dist .ge. R) then
  w = 0.0
 else
  w = R / dist + dist / R - 2.0
 end if
end function
!**************************************************************************
subroutine read(nn,AOA,dis,xp)!,itypep)
 !mainroutine_variable
 integer,intent(in)	:: nn !nump
 real,intent(in)	:: AOA
 real,intent(in)	:: dis
 real,intent(in)	:: xp(0:nn-1,0:2) !x
! integer,intent(in)	:: itypep(0:nn-1)
 !subroutine_variable
 integer	:: i,j,k
 integer	:: n
 integer	:: ip,jp
 integer	:: is,ie,js,je,ks,ke
 real,allocatable	:: xi(:,:)
 real,allocatable	:: v(:,:)
 real,allocatable	:: x(:,:,:)
 real,allocatable	:: y(:,:,:)
 real,allocatable	:: z(:,:,:)
 real,allocatable	:: x2(:,:,:)
 real,allocatable	:: y2(:,:,:)
 real,allocatable	:: z2(:,:,:)
 real,allocatable	:: xsub(:,:,:)
 real,allocatable	:: ysub(:,:,:)
 real,allocatable	:: zsub(:,:,:)
 real,allocatable	:: icing_col(:,:,:)
 real,allocatable	:: cell_data(:,:)
 real,allocatable	:: subcell_data(:,:)
 real,allocatable	:: cell_area(:,:)
 real			:: theta
 real			:: pi = 3.141592
 real,allocatable	:: f(:,:,:)
 real,allocatable	:: fsub(:,:,:)
 real			:: IC_rat
 integer		:: count

 integer		:: is1,ie1

 IC_rat = 0.7

 theta = AOA / 180.0 * PI

 allocate(xi(0:nn-1,0:2))
 allocate(v(0:nn-1,0:2))

 do n=0,nn-1
  xi(n,0) = xp(n,0)
  xi(n,1) = xp(n,1)
 end do

 do n=0,nn-1!翼座標
  v(n,0) = xi(n,0)
  v(n,1) = xi(n,1)
  xi(n,0) = v(n,0)*cos(theta) - v(n,1)*sin(theta)
  xi(n,1) = v(n,0)*sin(theta) + v(n,1)*cos(theta)
 end do

 !-----------------MAIN GRID-----------------------------
 is=0
 ie=260
 js=0
 je=70
 ks=0
 ke=3

 allocate(f(is:ie,js:je,ks:ke))
 allocate(x(is:ie,js:je,ks:ke))
 allocate(y(is:ie,js:je,ks:ke))
 allocate(z(is:ie,js:je,ks:ke))
 allocate(x2(is:ie,js:je,ks:ke))
 allocate(y2(is:ie,js:je,ks:ke))
 allocate(z2(is:ie,js:je,ks:ke))
 allocate(icing_col(is:ie,js:je,ks:ke))
 allocate(cell_data(is:ie,js:je))
 allocate(cell_area(is:ie,js:je))

 open(1, FILE='./result/GridFlow/grid/clean/Main_ViewGrid.bin', form = 'unformatted', status = 'old')
  read(1) f
  read(1) x
  read(1) y
  read(1) z
 close(1)

 cell_data = 0.0

 k=1

 do n=0,nn-1
!  if(itypep(n) .eq. 2)then
   ip = -1
   jp = -1
   j = 0
   do i = is, (ie/2)+1
    if(v(n,1) .lt. 0.0)then
     if(xi(n,1) .ge. y(i,j,k)+((y(i,j+1,k)-y(i,j,k))* &
     & (xi(n,0)-x(i,j,k))/(x(i,j+1,k)-x(i,j,k))))then
      if(xi(n,1) .lt. y(i+1,j,k)+((y(i+1,j+1,k)-y(i+1,j,k))* &
      & (xi(n,0)-x(i+1,j,k))/(x(i+1,j+1,k)-x(i+1,j,k))))then
       ip = i
       exit
      end if
     end if
    end if
   end do
   do i = (ie/2)-1, ie-1
    if(v(n,1) .ge. 0.0)then
     if(xi(n,1) .ge. y(i,j,k)+((y(i,j+1,k)-y(i,j,k))* &
     & (xi(n,0)-x(i,j,k))/(x(i,j+1,k)-x(i,j,k))))then
      if(xi(n,1) .lt. y(i+1,j,k)+((y(i+1,j+1,k)-y(i+1,j,k))* &
      & (xi(n,0)-x(i+1,j,k))/(x(i+1,j+1,k)-x(i+1,j,k))))then
       ip = i
       exit
      end if
     end if
    end if
   end do

   if(ip .ne. -1)then
    do j = js,je-1
     if(ip .le. (ie/2)+1)then
      if(xi(n,1) .le. y(ip,j,k)+((y(ip+1,j,k)-y(ip,j,k))* &
      & (xi(n,0)-x(ip,j,k))/(x(ip+1,j,k)-x(ip,j,k))))then
       if(xi(n,1) .gt. y(ip,j+1,k)+((y(ip+1,j+1,k)-y(ip,j+1,k))* &
       & (xi(n,0)-x(ip,j+1,k))/(x(ip+1,j+1,k)-x(ip,j+1,k))))then
        jp = j
        exit
       end if
      end if
     end if
     if(ip .ge. (ie/2)-1)then
      if(xi(n,1) .ge. y(ip,j,k)+((y(ip+1,j,k)-y(ip,j,k))* &
      & (xi(n,0)-x(ip,j,k))/(x(ip+1,j,k)-x(ip,j,k))))then
       if(xi(n,1) .lt. y(ip,j+1,k)+((y(ip+1,j+1,k)-y(ip,j+1,k))* &
       & (xi(n,0)-x(ip,j+1,k))/(x(ip+1,j+1,k)-x(ip,j+1,k))))then
        jp = j
        exit
       end if
      end if
     end if
    end do
   end if
   if((ip .ne. -1) .and. (jp .ne. -1)) cell_data(ip,jp) = cell_data(ip,jp) + 1.0
!  end if
 end do

 k=1
 cell_area = 0.0
 do j = js,je-1
  do i = is,ie-1
   cell_area(i,j) = 0.5*abs(((x(i+1,j,k)-x(i,j,k))*(y(i,j+1,k)-y(i,j,k))) &
   &                -((x(i,j+1,k)-x(i,j,k))*(y(i+1,j,k)-y(i,j,k))))
   cell_area(i,j) = cell_area(i,j) &
   &               +0.5*abs(((x(i+1,j,k)-x(i+1,j+1,k))*(y(i,j+1,k)-y(i+1,j+1,k))) &
   &               -((x(i,j+1,k)-x(i+1,j+1,k))*(y(i+1,j,k)-y(i+1,j+1,k))))
   cell_data(i,j) = cell_data(i,j) * pi*(dis**2)/4.0
  end do
 end do

 count = 0
 do k = ks,ke
  do j = js,je-1
   do i = is,ie-1
    if(cell_data(i,j) .ge. cell_area(i,j) * IC_rat )then
     icing_col(i,j,k) = 1.0
     if(k .eq. 1) count = count + 1
    end if
   end do
  end do 
 end do
 write(*,*) count,'Particles are Judged as Icing Cell at Main Grid'

 open(1,file='./data/ICM/icing.txt',status = 'replace')
  do k = ks, ke
   do j = js, je
    do i = is, ie
     write(1,*)icing_col(i,j,k) 
    end do
   end do
  end do
 close(1)
	
 open(1,file='./data/ICM/grid.txt',status = 'replace')
  do k=ks,ke
   do j=js,je
    do i=is,ie
     write(1,*)x(i,j,k), y(i,j,k), z(i,j,k)
    end do
   end do
  end do
 close(1)

 open(1,file='./data/ICM/icepos.txt',status = 'replace')
  do k = ks, ke
   do j = js, je
    do i = is, ie
     if(icing_col(i,j,k) .eq. 1.0)then
      write(1,*)x(i,j,k), y(i,j,k), z(i,j,k)
     end if
    end do
   end do
  end do
 close(1)

 call output_icingdata(1,is,ie,js,je,ks,ke,icing_col)
 deallocate(icing_col)

 !-----------------SUB GRID-----------------------------
 is=0
 ie=300
 js=0
 je=50
 ks=0
 ke=3

 is1 = is + 90
 ie1 = ie - 90

 allocate(fsub(is:ie,js:je,ks:ke))
 allocate(xsub(is:ie,js:je,ks:ke))
 allocate(ysub(is:ie,js:je,ks:ke))
 allocate(zsub(is:ie,js:je,ks:ke))
 allocate(subcell_data(is:ie,js:je))
 allocate(icing_col(is:ie,js:je,ks:ke))

 icing_col= 0.0   !icing_colのリセット

! open(30, FILE='./grid_data/Sub_ND_Grid.bin', form = 'unformatted', status = 'old')
!  read(30) (((xsub(i,j,k), i = is, ie), j = js, je), k = ks, ke)
!  read(30) (((ysub(i,j,k), i = is, ie), j = js, je), k = ks, ke)
!  read(30) (((zsub(i,j,k), i = is, ie), j = js, je), k = ks, ke)
! close(30)

 open(1, FILE='./result/GridFlow/grid/clean/Sub_ViewGrid.bin', form = 'unformatted', status = 'old')
  read(1) fsub
  read(1) xsub
  read(1) ysub
  read(1) zsub
 close(1)

 subcell_data = 0.0

 k=1
 do n=0,nn-1
!  if(itypep(n) .eq. 2)then
   ip = -1
   jp = -1
   j = 0
   do i = is1,(ie/2)+1 !is, (ie/2)+1
    if(v(n,1) .le. 0.0)then
     if(xi(n,1) .ge. ysub(i,j,k)+((ysub(i,j+1,k)-ysub(i,j,k))* &
     & (xi(n,0)-xsub(i,j,k))/(xsub(i,j+1,k)-xsub(i,j,k))))then
      if(xi(n,1) .lt. ysub(i+1,j,k)+((ysub(i+1,j+1,k)-ysub(i+1,j,k))* &
      & (xi(n,0)-xsub(i+1,j,k))/(xsub(i+1,j+1,k)-xsub(i+1,j,k))))then
       ip = i
       exit
      end if
     end if
    end if
   end do
   do i = (ie/2)-1, ie1-1 !(ie/2)-1, ie-1
    if(v(n,1) .gt. 0.0)then
     if(xi(n,1) .ge. ysub(i,j,k)+((ysub(i,j+1,k)-ysub(i,j,k))* &
     & (xi(n,0)-xsub(i,j,k))/(xsub(i,j+1,k)-xsub(i,j,k))))then
      if(xi(n,1) .lt. ysub(i+1,j,k)+((ysub(i+1,j+1,k)-ysub(i+1,j,k))* &
      & (xi(n,0)-xsub(i+1,j,k))/(xsub(i+1,j+1,k)-xsub(i+1,j,k))))then
       ip = i
       exit
      end if
     end if
    end if
   end do
   if(ip .ne. -1)then
    do j = js,je-1
     if(ip .lt. (ie/2))then
      if(xi(n,1) .le. ysub(ip,j,k)+((ysub(ip+1,j,k)-ysub(ip,j,k))* &
      & (xi(n,0)-xsub(ip,j,k))/(xsub(ip+1,j,k)-xsub(ip,j,k))))then
       if(xi(n,1) .gt. ysub(ip,j+1,k)+((ysub(ip+1,j+1,k)-ysub(ip,j+1,k))* &
       & (xi(n,0)-xsub(ip,j+1,k))/(xsub(ip+1,j+1,k)-xsub(ip,j+1,k))))then
        jp = j
        exit
       end if
      end if
     end if
     if(ip .ge. (ie/2))then
      if(xi(n,1) .ge. ysub(ip,j,k)+((ysub(ip+1,j,k)-ysub(ip,j,k))* &
      & (xi(n,0)-xsub(ip,j,k))/(xsub(ip+1,j,k)-xsub(ip,j,k))))then
       if(xi(n,1) .lt. ysub(ip,j+1,k)+((ysub(ip+1,j+1,k)-ysub(ip,j+1,k))* &
       & (xi(n,0)-xsub(ip,j+1,k))/(xsub(ip+1,j+1,k)-xsub(ip,j+1,k))))then
        jp = j
        exit
       end if
      end if
     end if
    end do
   end if
   if((ip .ne. -1) .and. (jp .ne. -1)) cell_data(ip,jp) = cell_data(ip,jp) + 1.0
!  end if
 end do

  k=1
  cell_area = 0.0
  do j = js,je-1
    do i = is,ie-1
      cell_area(i,j) = 0.5*abs(((xsub(i+1,j,k)-xsub(i,j,k))*(ysub(i,j+1,k)-ysub(i,j,k))) &
      &                -((xsub(i,j+1,k)-xsub(i,j,k))*(ysub(i+1,j,k)-ysub(i,j,k))))
      cell_area(i,j) = cell_area(i,j) &
      &               +0.5*abs(((xsub(i+1,j,k)-xsub(i+1,j+1,k))*(ysub(i,j+1,k)-ysub(i+1,j+1,k))) &
      &               -((xsub(i,j+1,k)-xsub(i+1,j+1,k))*(ysub(i+1,j,k)-ysub(i+1,j+1,k))))
      cell_data(i,j) = cell_data(i,j) * pi*(dis**2)/4.0
    end do
  end do

 count = 0
 do k = ks,ke
  do j = js,je-1
   do i = is,ie-1
    if(cell_data(i,j) .ge. cell_area(i,j) * IC_rat )then
     icing_col(i,j,k) = 1.0
     if(k .eq. 1) count = count + 1
    end if
   end do
  end do 
 end do
 write(*,*) count,'Particles are Judged as Icing Cell at Sub Grid'

 open(1,file='./data/ICM/icingsub.txt',status = 'replace')
  do k=ks,ke
   do j=js,je
    do i=is,ie
     write(1,*)icing_col(i,j,k)
    end do
   end do
  end do
 close(1)

 open(1,file='./data/ICM/gridsub.txt',status = 'replace')
  do k=ks,ke
   do j=js,je
    do i=is,ie
     write(1,*)xsub(i,j,k), ysub(i,j,k), zsub(i,j,k)
    end do
   end do
  end do
 close(1)

 call output_icingdata(2,is,ie,js,je,ks,ke,icing_col)

end subroutine read
!**************************************************************************
subroutine output_icingdata(mm,is,ie,js,je,ks,ke,icing_col)
 implicit none
 !mainroutine_variale
 integer,intent(in)	:: mm !1:Main,2:Sub
 integer,intent(in)	:: is,ie,js,je,ks,ke
 real,intent(in)	:: icing_col(is:ie,js:je,ks:ke)

 !subroutine_variable
 integer		:: i1, i2
 integer		:: i_back
 character(5)		:: head
 character(40)		:: fn1,fn2
 integer,allocatable	:: ic(:, :)	!着氷セル
 integer,allocatable	:: ji(:)		!着氷高さ
 integer,allocatable	:: flag(:, :)	!着氷判定

 i_back = 50

 !file_name
 select case(mm)
  case (1)
   head = 'Main_'
  case (2)
   head = 'Sub_'
 end select
 fn1 = './result/ICM/'//trim(head)//'IcingIndex.dat'
 fn2 = './result/ICM/'//trim(head)//'IcingCell.dat'

 !array_allocation
 allocate(ic(is:ie,js:je))
 allocate(ji(is:ie))
 allocate(flag(is:ie,js:je))

 !calculate_j-index
 do i = is,ie
  do j = js + 1,je
   if(int(icing_col(i,j,0)) == 0) then
    ji(i) = j-1
    exit
   end if
  end do
 end do

 !remove_back_of_airfoil
 if(mm == 1) then
  i1 = is + i_back
  i2 = ie - i_back
 else
  i1 = is + i_back
  i2 = ie - i_back
 end if

 !calculate_icing-cell
  ! 場合分け --------------------------------------------------------------------------------------------
  ! ic =  0; 流体セル
  ! ic =  1; 着氷セル (完全に氷の中)
  ! ic =  2; i-1が流体セル，i+1が流体セルであるj方向に最大の氷層セル
  ! ic =  3; i-1が氷層セル，i+1が氷層セルであるj方向に最大の氷層セル
  ! ic =  4; i-1が流体セル，i+1が氷層セルであるj方向に最大の氷層セル
  ! ic =  5; i-1が氷層セル，i+1が流体セルであるj方向に最大の氷層セル
  ! ic =  6; i-1が氷層セル，i+1が氷層セルであるj方向に最小の氷層セル
  ! ic =  7; i-1が着氷セル，i+1が氷層セルであるj方向に最小の氷層セル	おそらく稲川さんのコードは7と8が逆
  !									最小セルは流れ場計算には使わないので問題なしなはず
  ! ic =  8; i-1が氷層セル，i+1が着氷セルであるj方向に最小の氷層セル
  ! ic =  9; i-1が流体セル，i+1が流体セルである氷層セル
  ! ic = 10; i-1が流体セル，i+1が着氷セルである氷層セル
  ! ic = 11; i-1が着氷セル，i+1が流体セルである氷層セル
  ! ic = 12; 壁面セル (着氷なし)
 ! ic = 1
 ic(:,:) = icing_col(:,:,0) !着氷セルは最大高さだけでなく全て1で初期設定で良いと思われ
 flag(:,:) = ic(:,:)
  ! j 方向最大 ------------------------------------------------------------------------------------------
  do i = i1 + 1,i2 - 1
   if(ji(i) == js .and. ic(i,ji(i)) == 0) cycle
   ! ic =  2
   if( flag(i-1,ji(i)) == 0 .and. flag(i+1,ji(i)) == 0) then
     ic(i,ji(i)) = 2
   ! ic =  3
   else if(flag(i-1,ji(i)) == 1 .and. flag(i+1,ji(i)) == 1) then
    ic(i,ji(i)) = 3
   ! ic =  4
   else if(flag(i-1,ji(i)) == 0 .and. flag(i+1,ji(i)) == 1) then
    ic(i,ji(i)) = 4
   ! ic =  5
   else
    ic(i,ji(i)) = 5
   endif
  enddo
  ! j 方向最小 ------------------------------------------------------------------------------------------
  do i = is + 1,ie - 1
   if(ji(i) == js) cycle
   if(ji(i) >= min(ji(i-1), ji(i+1)) &
   &   .and. (ji(i) /= ji(i-1) .or. ji(i) /= ji(i+1))) then
    ! ic = 6
    if(ji(i) > max(ji(i-1), ji(i+1))) then
     ic(i, min(ji(i-1), ji(i+1))) = 6
    ! ic = 7
    else if(ji(i-1) < ji(i)) then
     ic(i, ji(i-1)) = 7
    ! ic = 8
    else if(ji(i+1) < ji(i)) then
     ic(i, ji(i+1)) = 8
    end if
   endif
  end do
  ! j 方向その他 ----------------------------------------------------------------------------------------
  do i = is + 1, ie - 1
   if(ji(i) == js) cycle
   do j = ji(i) - 1, js, -1
    if(6 <= ic(i,j) .and. ic(i,j) <= 8) exit
    ! ic = 9
    if(ic(i-1,j) == 0 .and. ic(i+1,j) == 0) then
     ic(i,j) =  9
    ! ic = 10
    else if(ic(i-1,j) == 0) then
     ic(i,j) = 10
    ! ic = 11
    else if(ic(i+1,j) == 0) then
     ic(i,j) = 11
    end if
   end do
  end do
  ! 壁面セル --------------------------------------------------------------------------------------------
  do i = i1, i2
   if( ic(i,js) == 0 ) then
    ic(i,js) = 12
   endif
  enddo
 
 !output_Icingfile
 open(1,file = fn1,status = 'replace')
  do i = is,ie
   write(1,'(I2)') ji(i)
  end do
 close(1)
 
 open(1,file = fn2,status = 'replace')
  do j = js,je
   do i = is,ie
    write(1,'(I2)') ic(i,j)
   end do
  end do
 close(1)

 deallocate(ic)
 deallocate(ji)
 deallocate(flag)
 
end subroutine output_icingdata

end program DeleteInnerParticle








