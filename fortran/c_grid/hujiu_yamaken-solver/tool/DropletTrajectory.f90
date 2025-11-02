!液滴軌道fldファイル結合プログラム
!"Drop~~~~~_DropTrajectory.fld"と同じファイル内で実行
program DropletTrajectory
implicit none
!variable
integer,parameter	:: nn = 10000000
integer			:: n
integer			:: fnum
character(len=200)	:: filename
character(len=5)	:: fnum_char
real,allocatable	:: x_droplet(:,:)
integer,allocatable	:: n_droplet(:)
integer			:: count
integer			:: nn_droplet
character(len=10)	:: dummy
integer			:: interval

interval = 1 !何ファイルおきに読むか(1に近づけるとファイルサイズが大きくなるので注意)

!array_allocation
allocate(x_droplet(0:nn-1,0:1))
allocate(n_droplet(0:nn-1))

write(*,*) '***** Computation Strat *****'

!input_data
count = 0
fnum = 1
do
 write(fnum_char,'(I5.5)') fnum
 filename='Drop'//fnum_char//'_DropTrajectory.fld'
 if(access(trim(filename),'r') .gt. 0) exit
 write(*,*) filename
 open(1,file=trim(filename),status='old')
  !hedder
  do n = 0,1
   read(1,*) dummy
  end do
  read(1,*) dummy,dummy,nn_droplet
  do n = 0,10
   read(1,*) dummy
  end do
  !data
  do n = 0,nn_droplet-1
   read(1,*) x_droplet(count,0),x_droplet(count,1)
   n_droplet(count) = fnum
   count = count + 1
   if(count .ge. nn) then
    write(*,*) 'nn should be larger number'
    stop
   end if
  end do
 close(1)
 fnum = fnum + interval
end do

!output_data
filename='DropletTrajectory.fld'
open(1,file=trim(filename),status='replace')
 !hedder
 write(1,'(a)') '# AVS field file'
 write(1,'(a)') 'ndim = 1'
 write(1,'(a, i6)') 'dim1 = ', count
 write(1,'(a)') 'nspace = 3'
 write(1,'(a)') 'veclen = 1'
 write(1,'(a)') 'data = float'
 write(1,'(a)') 'field = irregular'
 write(1,'(a)') 'label = n'
 !data
 write(1,'(3a)') 'variable 1 file= ', trim(filename), &
 &               ' filetype = ascii skip = 13 offset = 2 stride = 3 close = 1'
 write(1,'(3a)') 'coord    1 file= ', trim(filename), &
 &               ' filetype = ascii skip = 13 offset = 0 stride = 3 close = 1'
 write(1,'(3a)') 'coord    2 file= ', trim(filename), &
 &               ' filetype = ascii skip = 13 offset = 1 stride = 3 close = 1'
 write(1,'(a)') 'eot'
 write(1,'(a)') ''
 do n = 0,count-1
  write(1,*) x_droplet(n,0),x_droplet(n,1),n_droplet(n)
 end do
close(1)

filename='DropletTrajectory.vtk'
open(1,file=trim(filename),status='replace')
 write(1,'(A)') '# vtk DataFile Version 3.0'
 write(1,*) 'vtk output'
 write(1,*) 'ASCII'
 write(1,*) 'DATASET POLYDATA'
 write(1,*) 'POINTS ',count,' float'
 do n = 0,count-1
  write(1,*) x_droplet(n,0),x_droplet(n,1),0.0
 end do
 write(1,*) 'POINT_DATA ', count
 write(1,*) 'SCALARS n float 1'
 write(1,*) 'LOOKUP_TABLE default'
 do n = 0,count-1
  write(1,*) n_droplet(n)
 end do
close(1)

write(*,*) '***** Computation End *****'

end program DropletTrajectory