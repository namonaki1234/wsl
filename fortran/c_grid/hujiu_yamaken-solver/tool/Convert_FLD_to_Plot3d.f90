!*****************************************
!*****   Convert_FLD_to_Plot3d.f90   *****
!*****************************************
!******2019.11.14 written by TomoyaYAGI***
!MicroAVS用のFLD形式ファイルからplot3d形式のファイルに変換するプログラム
!流れ場の可視化に対応

program grid

!変数宣言
implicit none
integer i,j,k,m,n,l					!変数
integer nb,imax,jmax,kmax				!マルチブロック数,各座標の最大数
integer dimensions					!fld形式のおける次元数
integer variables					!fld形式のおける変数の数
integer,allocatable::ni(:),nj(:),nk(:),nl(:)		!各座標の最大値+ファンクション数
real,allocatable::x(:,:,:,:),y(:,:,:,:),z(:,:,:,:)	!xyz座標
real,allocatable::flag(:,:,:,:)				!fld形式におけるブロック番号データ
real,allocatable::fun(:,:,:,:,:)			!物理量データ（ファンクションファイル用）
character(len=50) fn,ex_bin,ex_g,ex_f,fn1,fn2,fn3	!読み込むファイル名
character(len=10) q1,q2					!対話用

!対話===============================================================================================
print*,'*****************************************'
print*,'*****   Convert_FLD_to_Plot3d.f90   *****'
print*,'*****************************************'
print*,'---------------------------------------------------------'
print*,'Conversion for grid or flowfield?'
print*,'grid      ==> input "g"'
print*,'flowfield ==> input "f"'
read*,q1
if (q1 == 'g' .or. q1 == 'f') then
else
 stop
end if
print*,'---------------------------------------------------------'
print*,'Input the number of dimensions'
read*,dimensions
print*,'Input the number of grid points in i direction'
read*,imax
print*,'Input the number of grid points in j direction'
read*,jmax
if (dimensions == 3) then
 print*,'Input the number of grid points in k direction'
 read*,kmax
end if
if (q1 == 'f') then
 print*,'Input the number of variables (u,v,rho...etc )'
 read*,variables
end if
print*,'---------------------------------------------------------'
print*,'Input the file name before conversion'
print*,'ex) Main_ViewFlow.bin ==> input "Main_ViewFlow"'
read*,fn
ex_bin = '.bin'
ex_g   = '.g'
ex_f   = '.f'
fn1    = trim(fn)//trim(ex_bin)
fn2    = trim(fn)//trim(ex_g  )
fn3    = trim(fn)//trim(ex_f  )

!格子データの場合===================================================================================
if (q1 == 'g') then
 !ファイル入力
 if (dimensions == 2) then
 nb   = 1		!ブロック数は1
 kmax = 1		!k方向の格子点数は1
 allocate(ni(nb),nj(nb),nk(nb))
 ni(1) = imax
 nj(1) = jmax
 nk(1) = kmax
 allocate(flag(imax,jmax,kmax,nb))
 allocate(x(imax,jmax,kmax,nb),y(imax,jmax,kmax,nb),z(imax,jmax,kmax,nb))
 
 open (unit=7, form='unformatted', file=fn1)
 read(7) flag		!ブロック数読み込み（1ブロックを想定しているため、ここの値は読み込んでも使わない）
 read(7) x		!x座標読み込み
 read(7) y		!y座標読み込み
 close(7)
 
 end if
 
 if (dimensions == 3) then
 nb   = 1		!ブロック数は1
 allocate(ni(nb),nj(nb),nk(nb))
 ni(1) = imax
 nj(1) = jmax
 nk(1) = kmax
 allocate(flag(imax,jmax,kmax,nb))
 allocate(x(imax,jmax,kmax,nb),y(imax,jmax,kmax,nb),z(imax,jmax,kmax,nb))
 
 open (unit=7, form='unformatted', file=fn1)
 read(7) flag		!ブロック数読み込み（1ブロックを想定しているため、ここの値は読み込んでも使わない）
 read(7) x		!x座標読み込み
 read(7) y		!y座標読み込み
 read(7) z		!z座標読み込み
 close(7)
 
 end if
 
 !ファイル出力
 open (unit=7, form='unformatted', file=fn2)
 write(7) nb
 write(7) (ni(m),nj(m),nk(m),m=1,nb)
 do m=1, nb
  write(7) &
   ((( x(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)), &
   ((( y(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)), &
   ((( z(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m))
 end do
 close(8)
 print*,'Output "',trim(fn2),'" successfully!!'
end if

!流れ場データの場合=================================================================================
if (q1 == 'f') then
 !ファイル入力
 if (dimensions == 2) then
 nb   = 1		!ブロック数は1
 kmax = 1		!k方向の格子点数は1
 allocate(ni(nb),nj(nb),nk(nb),nl(nb))
 ni(1) = imax
 nj(1) = jmax
 nk(1) = kmax
 nl(1) = variables+1
 allocate(x(imax,jmax,kmax,nb),y(imax,jmax,kmax,nb),z(imax,jmax,kmax,nb))
 allocate(fun(imax,jmax,kmax,variables+1,nb))
 
 open (unit=7, form='unformatted', file=fn1)
 do l=1, variables
  read(7) fun(:,:,:,l,:)	!物理量データの読み込み
 end do
 read(7) x			!x座標読み込み
 read(7) y			!y座標読み込み
 close(7)
 
 fun(:,:,:,variables+1,:) = sqrt(fun(:,:,:,1,:)**2+fun(:,:,:,2,:)**2)	!速度場計算
 
 end if
 
 if (dimensions == 3) then
 nb   = 1		!ブロック数は1
 allocate(ni(nb),nj(nb),nk(nb),nl(nb))
 ni(1) = imax
 nj(1) = jmax
 nk(1) = kmax
 nl(1) = variables+1
 allocate(x(imax,jmax,kmax,nb),y(imax,jmax,kmax,nb),z(imax,jmax,kmax,nb))
 allocate(fun(imax,jmax,kmax,variables+1,nb))
 
 open (unit=7, form='unformatted', file=fn1)
 do l=1, variables
  read(7) fun(:,:,:,l,:)	!物理量データの読み込み
 end do
 read(7) x			!x座標読み込み
 read(7) y			!y座標読み込み
 read(7) z			!z座標読み込み
 close(7)
 
 fun(:,:,:,variables+1,:) = sqrt(fun(:,:,:,1,:)**2+fun(:,:,:,2,:)**2+fun(:,:,:,3,:)**2)	!速度場計算
 
 end if
 
 !ファイル出力
 open (unit=7, form='unformatted', file=fn2)
 write(7) nb
 write(7) (ni(m),nj(m),nk(m),m=1,nb)
 do m=1, nb
  write(7) &
   ((( x(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)), &
   ((( y(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)), &
   ((( z(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m))
 end do
 close(7)
 print*,'Output "',trim(fn2),'" successfully!!'
 
 open (unit=8, form='unformatted', file=fn3)
 write(8) nb
 write(8) (ni(m), nj(m), nk(m), nl(m), m=1,nb )
 do m=1, nb
 write(8) &
    (((( fun(i,j,k,l,m),  i=1,ni(m)), j=1,nj(m)), k=1,nk(m)), l=1,nl(m))
 end do
 close(8)
 print*,'Output "',trim(fn3),'" successfully!!'
end if

end program grid