!*******************************************************************************************************
!*******************************************************************************************************
!******** 着氷後の格子再構成プログラム								********
!******** (NACA翼，三次元，C-type，重合格子法)							********
!********					      2013.02.07  PROGRAMMED BY RYOSUKE HAYASHI ********
!******** ※ ダブルホーンに対応									********
!********					      2013.05.05     UPDATED BY RYOSUKE HAYASHI ********
!********					      2013.07.22     UPDATED BY RYOSUKE HAYASHI ********
!*******************************************************************************************************
!*******************************************************************************************************
program Remesh_NACA
 ! モジュール宣言 **************************************************************************************
 use Package_NACA
 use Package_FileIO
 use Package_Equation
 use Package_Grid
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! ファイル名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: ViewGrdFile *  8 = 'ViewGrid'
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, parameter :: mRef = 2
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: kRef
 ! 低レイノルズ数型格子の場合,格子を切ることが困難であると予想されるため,
 ! スイッチでステップを進めるか決める(0:進めない,1:進める)
 ! 0で切れるか試して,上手くいけば1で本番
 integer :: swi_IceStep = 1
 ! 処理開始 ********************************************************************************************
 ! 検証実験ケース選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(a)') "<< Exp. Case Selection >>"
 call SelectExpCase
 ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') '<< Initial Setting >>'
 call InitialSetting
 ! スムージング ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Smoothing >>"
 call SmoothIce
 ! 着氷後の格子 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Icing Grid Geometry >>"
 call IcingGrid
 ! 着氷後の流れ場初期値 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') '<< Initial Flow Condition >>'
 call InitialFlow
 ! 内部手続き ******************************************************************************************
 stop
contains
!*******************************************************************************************************
!******** 検証実験ケース選択 									********
!*******************************************************************************************************
subroutine SelectExpCase
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character :: fname * 20
 ! 処理開始 ********************************************************************************************
 ! 計算条件ファイル入力
 call Input_CalSetting( trim(ND_CalSetFile) // strtxt )
 Span = Span * lRef
 ! ディレクトリ設定
 if( IceStep == 0 ) then
   GrdInDir     = bckdir // 'grid//clean//'
   FlwCalInDir  = bckdir // 'flow//cal//clean//'
   FlwIniDir    = bckdir // 'flow//initial//clean//'
  else
   GrdInDir     = bckdir // 'grid//icing//'
   FlwCalInDir  = bckdir // 'flow//cal//icing//'
   FlwIniDir    = bckdir // 'flow//initial//icing//'
 endif
 IceCalInDir  = bckdir // 'icing//cal//'
 GrdOutDir    = bckdir // 'grid//icing//'
 FlwCalOutDir = bckdir // 'flow//initial//icing//'
 write(*, '(a)') '+++ Icing Step +++'
 write(*, '(a,i2)') '* Ice step      = ', IceStep
 write(*, '(a,i2)') '* Ice step max. = ', IceStepMax
 write(*, '(/,a)') '+++ Exp. Condition +++'
 write(*, '(a,e16.8e3)') '* Ts    = ', TsExp * aRef**2
 write(*, '(a,e16.8e3)') '* Ps    = ', PsExp * (rhoRef * aRef**2)
 write(*, '(a,e16.8e3)') '* V     = ', VelExp * aRef
 write(*, '(a,e16.8e3)') '* LWC   = ', LWC * RhoRef
 write(*, '(a,e16.8e3)') '* MVD   = ', MVD * LRef
 write(*, '(a,e16.8e3)') '* Rho   = ', Rhod * RhoRef
 write(*, '(a,e16.8e3)') '* Chord = ', Chord * LRef
 write(*, '(a,e16.8e3)') '* AOA   = ', AOA * 180.0 / pi
 ! 処理終了 ********************************************************************************************
 return
end subroutine SelectExpCase
!*******************************************************************************************************
!********* 初期設定										********
!*******************************************************************************************************
subroutine InitialSetting
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: m, i, k
integer :: j
 integer   :: j0
 character :: fname * 20
 ! 処理開始 ********************************************************************************************
 ! ブロック名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Set_BlockName
 ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( Flw(ms:me), Ice(ms:me) )
 ! 格子解像度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_Resolution1D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(IceRslFile), strtxt, &
  &      Ice(m)%is, Ice(m)%ie )
  call Input_Resolution3D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(RslFile), strtxt, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke )
! enddo
 ! C 型格子分割点 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! m = 1
 m = me
 call Input_CtypeGridPoint( &
 &      trim(GrdInDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
 &      Flw(m)%i1, Flw(m)%i2, Flw(m)%i3 )
 Ice(m)%i1 = Flw(m)%i1
 Ice(m)%i2 = Flw(m)%i2
 Ice(m)%i3 = Flw(m)%i3
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  allocate( Flw(m)%x  ( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%y  ( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%z  ( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%f  ( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%xix( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%xiy( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%xiz( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%etx( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%ety( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%etz( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%zex( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%zey( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%zez( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%jac( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%qh ( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le) )
  allocate( Ice(m)%x  ( Ice(m)%is: Ice(m)%ie ), &
  &         Ice(m)%y  ( Ice(m)%is: Ice(m)%ie ), &
  &         Ice(m)%z  ( Ice(m)%is: Ice(m)%ie ), &
  &         Ice(m)%f  ( Ice(m)%is: Ice(m)%ie ), &
  &         Ice(m)%Bi ( Ice(m)%is: Ice(m)%ie ), &
  &         Ice(m)%dBi( Ice(m)%is: Ice(m)%ie ), &
  &         Ice(m)%Ti ( Ice(m)%is: Ice(m)%ie ) )
! enddo
 ! 格子座標 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_Grid3D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(GrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
! enddo
 ! 氷層厚さ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
!  if(m == 1) cycle
  call Input_IceThickTem2D( &
  &      trim(IceCalInDir) // trim(BlkName(m)) // trim(IceThickTemFile), strdat, &
  &      Ice(m)%is, Ice(m)%ie, &
  &      Ice(m)%f, Ice(m)%Bi, Ice(m)%dBi, Ice(m)%Ti )
  call Input_IceBladeSurface2D( &
  &      trim(IceCalInDir) // trim(BlkName(m)) // trim(IceBladeFile), strdat, &
  &      Ice(m)%is, Ice(m)%ie, &
  &      Ice(m)%x, Ice(m)%y )
! enddo
 ! 流れ場初期解 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_Flux3D( &
  &      trim(FlwCalInDir) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      ls, le, Flw(m)%qh )
! enddo
 ! 対象ブロック ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 m    = mRef
 kRef = Flw(m)%ks + int( 0.5 * (Flw(m)%ke - Flw(m)%ks) )
 ! 着氷限界位置探索 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do i = Ice(m)%i1, Ice(m)%i3 !Flw(m)%is, Flw(m)%ie
  if(Ice(m)%f(i) /= 0) then
    Flw(m)%i1 = i
    exit
  endif
 enddo
 do i = Ice(m)%i3, Ice(m)%i1, -1 !Flw(m)%ie, Flw(m)%is, -1
  if(Ice(m)%f(i) /= 0) then
    Flw(m)%i2 = i
    exit
  endif
 enddo
 if( Flw(m)%i1 < Flw(m)%is + 5 .or. Flw(m)%i2 > Flw(m)%ie - 5 ) then
   write(*, '(a)') '!!!!! Error : Icing limit point !!!!!'
   write(*, '(2i4)') Flw(m)%i1, Flw(m)%i2
   stop
 endif
  ! 処理終了 ********************************************************************************************
 return
end subroutine InitialSetting
!*******************************************************************************************************
!******** 氷層厚さのスムージング 								********
!*******************************************************************************************************
subroutine SmoothIce
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, parameter :: nSmooth = 5 !20 !5
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: GrdCleanDir * 20 = '../grid/clean/'
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: Cx(:,:,:), Cy(:,:,:), Cz(:,:,:)
 real   , pointer :: Bi0(:)
 real   , pointer :: dBi0(:)
 integer :: i, k, m, n
 integer :: j0, iHorn
 integer :: iHorn0, iHorn1, iHorn2
 real    :: ax, ay, az, bx, by, bz, nx, ny, nz, na
 logical :: DoubleHorn
 ! 処理開始 ********************************************************************************************
 m = mRef; k = kRef; DoubleHorn = .false.
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( Cx( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
 &         Cy( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
 &         Cz( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ) )
 allocate( Bi0(Ice(m)%is: Ice(m)%ie), dBi0(Ice(m)%is: Ice(m)%ie) )
 ! 着氷前翼座標 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Input_Grid3D( &
 &      trim(GrdCleanDir) // trim(BlkName(m)) // trim(GrdFile), strbin, &
 &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
 &      Cx, Cy, Cz )
 ! 1ステップ前の氷層高さ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Bi0(:) = Ice(m)%Bi(:) - Ice(m)%dBi(:)
 ! ホーン位置の探索 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(a,i3,a,i4)') '* Ice limit index (i) = ', Flw(m)%i1, ',', Flw(m)%i2
 do i = Flw(m)%i1, Flw(m)%i2
  if( Ice(m)%dBi(i) >= maxval(Ice(m)%dBi(:)) ) then
    iHorn1 = i
    exit
  endif
 enddo
 do i = iHorn1 - 1, Flw(m)%i1, - 1
  if( Ice(m)%dBi(i) - Ice(m)%dBi(i-1) < 0.0 ) then
    iHorn0 = i
    DoubleHorn = .true.
    iHorn2 = iHorn1
    exit
  endif
 enddo
 if(DoubleHorn) then
   do i = iHorn0, Flw(m)%i1, - 1
    if( Ice(m)%dBi(i) - Ice(m)%dBi(i+1) < 0.0 ) exit
   enddo
   iHorn1 = i
 endif
 if(DoubleHorn) then
   write(*, '(a,i3)') '* Horn-1 index (i) = ', iHorn1
   write(*, '(a,i3)') '* Horn-2 index (i) = ', iHorn2
  else
   write(*, '(a,i3)') '* Horn-1 index (i) = ', iHorn1
 endif
 ! スムージング ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do n = 1, nSmooth
  dBi0(:) = Ice(m)%dBi(:)
  do i = Flw(m)%i1, Flw(m)%i2
   Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
  enddo

!  ! Single Horn ----------------------------------------------------------------------------------------
!  if(.not. DoubleHorn) then
!    ! 終端
!    if(Flw(m)%i2 - 1 <= iHorn1) then
!      do i = iHorn1 - 1, Flw(m)%i1, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    ! 始端
!     else if(iHorn1 <= Flw(m)%i1 + 1) then
!      do i = iHorn1 + 1, Flw(m)%i2
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    ! その他
!     else
!      do i = iHorn1 + 1, Flw(m)%i2
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn1 - 1, Flw(m)%i1, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    endif
!    write(*, '(i1,6(e16.8e3))') n, Ice(m)%dBi(iHorn1-2:iHorn1+2), sum(Ice(m)%dBi(:))
!  ! Double Horn ----------------------------------------------------------------------------------------
!   else
!    ! 終端かつ始端
!    if(iHorn1 <= Flw(m)%i1 + 1 .and. Flw(m)%i2 - 1 <= iHorn2) then
!      do i = iHorn1 + 1, iHorn0
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn2 - 1, iHorn0, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    ! 始端
!     else if(iHorn1 <= Flw(m)%i1 + 1) then
!      do i = iHorn1 + 1, iHorn0
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn2 + 1, Flw(m)%i2
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn2 - 1, iHorn0, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    ! 終端
!     else if(Flw(m)%i2 - 1 <= iHorn2) then
!      do i = iHorn2 - 1, iHorn0, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn1 + 1, iHorn0
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn1 - 1, Flw(m)%i1, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    ! その他
!     else
!      do i = iHorn1 + 1, iHorn0
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn1 - 1, Flw(m)%i1, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn2 + 1, Flw(m)%i2
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn2 - 1, iHorn0, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    endif
!    write(*, '(i1,6(e16.8e3))') n, Ice(m)%dBi(iHorn1-2:iHorn1+2), sum(Ice(m)%dBi(:))
!    write(*, '(i1,6(e16.8e3))') n, Ice(m)%dBi(iHorn2-2:iHorn2+2), sum(Ice(m)%dBi(:))
!  endif
 enddo
 ! スムージング後の翼座標 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 j0 = Flw(m)%js
 do i = Flw(m)%i1, Flw(m)%i2
  ! 単位法線ベクトル
  if(i == Flw(m)%is) then
    ax = 1.0 * ( - Cx(i,j0,k) + Cx(i+1,j0,k) )
    ay = 1.0 * ( - Cy(i,j0,k) + Cy(i+1,j0,k) )
    az = 1.0 * ( - Cz(i,j0,k) + Cz(i+1,j0,k) )
   else if(i == Flw(m)%ie) then
    ax = 1.0 * ( - Cx(i-1,j0,k) + Cx(i,j0,k) )
    ay = 1.0 * ( - Cy(i-1,j0,k) + Cy(i,j0,k) )
    az = 1.0 * ( - Cz(i-1,j0,k) + Cz(i,j0,k) )
   else
    ax = 0.5 * ( - Cx(i-1,j0,k) + Cx(i+1,j0,k) )
    ay = 0.5 * ( - Cy(i-1,j0,k) + Cy(i+1,j0,k) )
    az = 0.5 * ( - Cz(i-1,j0,k) + Cz(i+1,j0,k) )
  endif
  bx = 0.5 * ( + Cx(i,j0,k-1) - Cx(i,j0,k+1) )
  by = 0.5 * ( + Cy(i,j0,k-1) - Cy(i,j0,k+1) )
  bz = 0.5 * ( + Cz(i,j0,k-1) - Cz(i,j0,k+1) )
  nx = ay * bz - az * by
  ny = az * bx - ax * bz
  nz = ax * by - ay * bx
  na = sqrt(nx**2 + ny**2 + nz**2)
  nx = +1.0 * nx / na
  ny = +1.0 * ny / na
  nz = +1.0 * nz / na
  ! 翼法線方向に変化
!  Ice(m)%dBi(i) = Ice(m)%Bi(i) - Bi0(i)
  Ice(m)%x(i)   = Flw(m)%x(i,j0,k) + nx * Ice(m)%dBi(i)
  Ice(m)%y(i)   = Flw(m)%y(i,j0,k) + ny * Ice(m)%dBi(i)
!  Ice(m)%x(i)   = Cx(i,j0,k) + nx * Ice(m)%Bi(i)
!  Ice(m)%y(i)   = Cy(i,j0,k) + ny * Ice(m)%Bi(i)
 enddo
 ! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 if(swi_IceStep .eq. 1) then
! do m = ms, me
 m = me
!  if(m == 1) cycle
  call Output_IceThickTem3D( &
  &      trim(IceCalInDir) // trim(BlkName(m)) // trim(IceThickTemFile), strdat, &
  &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, &
  &      Ice(m)%f, Ice(m)%Bi, Ice(m)%dBi, Ice(m)%Ti )
! enddo
 end if
 deallocate(Cx, Cy, Cz, Bi0)
 ! 処理終了 ********************************************************************************************
 return
end subroutine SmoothIce
!*******************************************************************************************************
!********* 着氷後の格子再構成									********
!*******************************************************************************************************
subroutine IcingGrid
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: m, i
 integer   :: jp
 character :: fname1 * 20, fname2 * 20
 character :: fname * 20
 ! 処理開始 ********************************************************************************************
 ! 対象ブロック ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 m  = mRef
 jp = Flw(m)%je
 ! 格子生成 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call HtypeGridIceLE( &
! &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, Flw(m)%i1, Flw(m)%i2, &
! &      Ice(m)%x, Ice(m)%y, Flw(m)%x(:, jp, kRef), Flw(m)%y(:, jp, kRef), &
! &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
 call CtypeGridBlade( &
 &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
 &      Ice(m)%i1, Ice(m)%i2, Ice(m)%i3, &
 &      Ice(m)%x, Ice(m)%y, Flw(m)%x(:, jp, kRef), Flw(m)%y(:, jp, kRef), &
 &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
 ! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 if(swi_IceStep .eq. 1) then
 write(fname, '(a, 2(i2.2,a))') 'IceStep', IceStep, 'of', IceStepMax, '_'
! do m = ms, me
 m = me
!  select case(m)
!   case(1)
!    call Output_CtypeGridPoint( &
!    &      trim(GrdOutDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
!    &      Flw(m)%i1, Flw(m)%i2, Flw(m)%i3 )
!   case(2)
!    call Output_IceLimitPoint( &
!    &      trim(GrdOutDir) // trim(BlkName(m)) // trim(IceLimitPointFile) // strtxt, &
!    &      Flw(m)%i1, Flw(m)%i2 )
!  end select
  call Output_CtypeGridPoint( &
  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
  &      Ice(m)%i1, Ice(m)%i2, Ice(m)%i3 )
  call Output_IceLimitPoint( &
  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(IceLimitPointFile) // strtxt, &
  &      Flw(m)%i1, Flw(m)%i2 )

  call Output_Resolution1D( &
  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(IceRslFile), strtxt, &
  &      Ice(m)%is, Ice(m)%ie )
  call Output_Resolution3D( &
  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(RslFile), strtxt, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke )
  call Output_Grid3D( &
  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(GrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
  call Output_Grid3D( &
  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(fname) // trim(GrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
! enddo
 end if
 ! 可視化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  Flw(m)%f(:,:,:) = m
  call MakeMAVSFile3D( &
  &      trim(GrdOutDir), trim(BlkName(m)) // trim(ViewGrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%f , Flw(m)%x , Flw(m)%y , Flw(m)%z )
  call OutputPara_bin( &
   &     trim(GrdOutDir), trim(BlkName(m)) // trim(ViewGrdFile), &
   &     Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
   &     Flw(m)%x, Flw(m)%y, Flw(m)%z )
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine IcingGrid
!*******************************************************************************************************
!********* 翼前縁付近の着氷計算用の格子 (H-type)						********
!*******************************************************************************************************
subroutine HtypeGridIceLE( &
&            is, ie, js, je, ks, ke, i1, i2, xi, yi, xe, ye, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, js, je, ks, ke
 integer, intent(in)  :: i1, i2						! 着氷限界位置
 real   , intent(in)  :: xi(is:ie), yi(is:ie)				! 内部境界
 real   , intent(in)  :: xe(is:ie), ye(is:ie)				! 外部境界
 real   , intent(out) :: x(is:ie, js:je, ks:ke), &
 &                       y(is:ie, js:je, ks:ke), &
 &                       z(is:ie, js:je, ks:ke)
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: rs1 = 1.0 * 1.0e-3				! 翼表面境界格子幅
 real   , parameter :: re1 = 1.0 * 1.0e-3				! 翼遠方境界格子幅
 real   , parameter :: MGN = 1.5e-0					! 楕円-双曲の重みの許容範囲
 real   , parameter :: Rsd = 1.0e-5					! 楕円-双曲の収束判定値
 real   , parameter :: rs2 = 5.0 * 1.0e-1				! 翼表面境界格子幅
 real   , parameter :: re2 = 7.0 * 1.0e-2				! 翼遠方境界格子幅
 real   , parameter :: tb1 = 5.0 * 1.0e-2				! 内部境界直交性のパラメータ
 real   , parameter :: tb2 = 2.0 * 1.0e-1				! 外部境界直交性のパラメータ
 ! 処理開始 ********************************************************************************************
 ! 初期境界値 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call InitialBoundary( &
&       is, ie, js, je, xi, yi, xe, ye, x(:, :, ks), y(:, :, ks) )
 ! 側部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call HtypeSideBoundary( &
 &      is, ie, js, je, rs2, re2, x(:, :, ks), y(:, :, ks) )
 ! Transfinite 補間 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call TransfiniteInterpolation( &
 &      is, ie, js, je, x(:, :, ks), y(:, :, ks) )
 ! 楕円-双曲型偏微分法 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call HtypeGenerationEHPDE( &
 &      is, i1, i2, ie, js, je, rs1, re1, MGN, Rsd, x(:, :, ks), y(:, :, ks) )
! ! 二境界法 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call GenerationTwoBoundary( &
! &      is, ie, js, je, rs2, re2, tb1, tb2, x(:, :, ks), y(:, :, ks) )
 ! 三次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call ThreeDimensionalized( &
 &      is, ie, js, je, ks, ke, span, x, y, z )
 ! 処理終了 ********************************************************************************************
 return
end subroutine HtypeGridIceLE

!*******************************************************************************************************
!********* 翼周りの格子 (C-type)								********
!*******************************************************************************************************
subroutine CtypeGridBlade( &
&            is, ie, js, je, ks, ke, i1, i2, i3, xi, yi, xe, ye, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, js, je, ks, ke
 integer, intent(in)  :: i1, i2, i3
 real   , intent(in)  :: xi(i1:i3), yi(i1:i3)				! 内部境界
 real   , intent(in)  :: xe(is:ie), ye(is:ie)				! 外部境界
 real   , intent(inout) :: x(is:ie, js:je, ks:ke), &
 &                         y(is:ie, js:je, ks:ke), &
 &                         z(is:ie, js:je, ks:ke)
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: rs5 = 7.0 * 1.0e-2 !3.28649582292613 * 1.0e-3	! 内部境界格子幅
 real   , parameter :: re5 = 1.0 * 1.0e+2				! 外部境界格子幅
 real   , parameter :: tb1 = 1.0e-2 !2.6579 * 1.0e-1			! 内部境界直交性のパラメータ
 real   , parameter :: tb2 = 7.0 * 1.0e0				! 外部境界直交性のパラメータ
 ! 処理開始 ********************************************************************************************
 ! 初期境界値 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call InitialBoundaryCtype( &
 &      is, ie, i1, i3, js, je, xi, yi, xe, ye, x(:, :, ks), y(:, :, ks) )
 ! 二境界法 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call GenerationTwoBoundary( &
 &      is, ie, js, je, rs5, re5, tb1, tb2, x(:, :, ks), y(:, :, ks) )
 ! 三次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call ThreeDimensionalized( &
 &      is, ie, js, je, ks, ke, span, x, y, z )
 ! 処理終了 ********************************************************************************************
 return
end subroutine CtypeGridBlade

!*******************************************************************************************************
!******** 初期境界値										********
!*******************************************************************************************************
subroutine InitialBoundary( &
&            is, ie, js, je, xi, yi, xe, ye, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, js, je
 real   , intent(in)  :: xi(is:ie), yi(is:ie)				! 内部境界
 real   , intent(in)  :: xe(is:ie), ye(is:ie)				! 外部境界
 real   , intent(out) :: x(is:ie, js:je), y(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 do i = is, ie
  x(i,js) = xi(i)
  y(i,js) = yi(i)
  x(i,je) = xe(i)
  y(i,je) = ye(i)
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine InitialBoundary
!*******************************************************************************************************
!******** 初期境界値(C型格子用)									********
!*******************************************************************************************************
subroutine InitialBoundaryCtype( &
&            is, ie, i1, i3, js, je, xi, yi, xe, ye, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, js, je
 integer, intent(in)  :: i1, i3
 real   , intent(in)  :: xi(i1:i3), yi(i1:i3)				! 内部境界
 real   , intent(in)  :: xe(is:ie), ye(is:ie)				! 外部境界
 real   , intent(out) :: x(is:ie, js:je), y(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 do i = i1, i3
  x(i,js) = xi(i)
  y(i,js) = yi(i)
  x(i,je) = xe(i)
  y(i,je) = ye(i)
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine InitialBoundaryCtype

!*******************************************************************************************************
!******** 側部境界 (H-type)									********
!*******************************************************************************************************
subroutine HtypeSideBoundary( &
&            is, ie, js, je, rs, re, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je
 real   , intent(in)    :: rs, re
 real   , intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: t(:)
 integer :: j
 real    :: rr
 ! 処理開始 ********************************************************************************************
 allocate( t(js:je) )
 call GeometricInterpolation( rs / real(je-js+1), je-js+1, t, rr )
 do j = js + 1, je - 1
!  x(is,j) = x(is,js)
!  x(ie,j) = x(ie,js)
  x(is,j) = ( x(is,je) - x(is,js) ) * t(j) + x(is,js)
  x(ie,j) = ( x(ie,je) - x(ie,js) ) * t(j) + x(ie,js)
  y(is,j) = ( y(is,je) - y(is,js) ) * t(j) + y(is,js)
  y(ie,j) = ( y(ie,je) - y(ie,js) ) * t(j) + y(ie,js)
 enddo
 deallocate(t)
 ! 処理開始 ********************************************************************************************
 return
end subroutine HtypeSideBoundary
!*******************************************************************************************************
!******** Transfinite 補間									********
!*******************************************************************************************************
subroutine TransfiniteInterpolation( &
&            is, ie, js, je, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je
 real   , intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: alp1(:), alp2(:), bet1(:), bet2(:)
 integer :: i, j
 ! 処理開始 ********************************************************************************************
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( alp1(is:ie), alp2(is:ie), bet1(js:je), bet2(js:je) ) 
 ! Transfinite 補間 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 媒介変数
 do i = is, ie
  alp1(i) = real(ie-i) / real(ie)
  alp2(i) = 1.0 - alp1(i)
 enddo
 do j = js, je
  bet1(j) = real(je-j) / real(je)
  bet2(j) = 1.0 - bet1(j)
 enddo
 ! 補間
 call Transfinite2D(is, ie, js, je, alp1, alp2, bet1, bet2, x)
 call Transfinite2D(is, ie, js, je, alp1, alp2, bet1, bet2, y)
 ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 deallocate( alp1, alp2, bet1, bet2 )
 ! 処理終了 ********************************************************************************************
 return
end subroutine TransfiniteInterpolation
!*******************************************************************************************************
!******** 楕円-双曲型偏微分方程式法に基づく格子生成 (H-type)					********
!*******************************************************************************************************
subroutine HtypeGenerationEHPDE( &
&            is, i1, i2, ie, js, je, rs, re, margin, resi, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, i1, i2, ie, js, je
 real   , intent(in)    :: rs, re
 real   , intent(in)    :: margin					! 楕円-双曲型重み関数許容誤差
 real   , intent(in)    :: resi						! 計算ループ収束判定値
 real   , intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: omg  = 1.0					! 緩和係数
 integer, parameter :: nmax = 100000					! 最大計算回数
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: Cs(:, :)
 real   , pointer :: cv1m(:, :), cv1p(:, :), cv2m(:, :), cv2p(:, :)
 real   , pointer :: dx1m(:, :), dx1p(:, :), dy1m(:, :), dy1p(:, :)
 real   , pointer :: dx2m(:, :), dx2p(:, :), dy2m(:, :), dy2p(:, :)
 integer :: i, j, n
 real    :: dmax
 real    :: a1, a2, a3, b1, b2, b3, c1, c2, c3, hjs1, hjs2, hje1, hje2
 ! 処理開始 ********************************************************************************************
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( cs  (is:ie, js:je), &
 &         cv1m(is:ie, js:je), cv1p(is:ie, js:je), cv2m(is:ie, js:je), cv2p(is:ie, js:je), &
 &         dx1m(is:ie, js:je), dx1p(is:ie, js:je), dy1m(is:ie, js:je), dy1p(is:ie, js:je), &
 &         dx2m(is:ie, js:je), dx2p(is:ie, js:je), dy2m(is:ie, js:je), dy2p(is:ie, js:je) )
 ! パラメータ設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 双曲型重み関数 --------------------------------------------------------------------------------------
 do j = js+1, je-1
 do i = is+1, ie-1
  ! i - 方向
  if(i < i1) then
    cv1m(i,j) = 0.0
   else if(i <= i2) then
    cv1m(i,j) = real(j) / real(je)**2 * real(i2-i) / real(i2-i1)
   else
    cv1m(i,j) = 0.0
  endif	
  ! i + 方向
  if(i2 < i) then
    cv1p(i,j) = 0.0
   else if(i1 <= i) then
    cv1p(i,j) = real(j) / real(je)**2 * real(i-i1) / real(i2-i1)
   else
    cv1p(i,j) = 0.0
  endif
  ! j - 方向
  cv2m(i,j) = (real(je-j) / real(je))**2
  ! j + 方向
  cv2p(i,j) = 0.0
 enddo
 enddo
 ! 楕円型重み関数 --------------------------------------------------------------------------------------
 do j = js+1, je-1
 do i = is+1, ie-1
  ! 内部境界法線ベクトル
  a1 =  0.5 * ( - x(i-1,js) + x(i+1,js) )
  a2 =  0.5 * ( - y(i-1,js) + y(i+1,js) )
  a3 =  0.0
  b1 =  0.0
  b2 =  0.0
  b3 = -1.0
  c1 = a2 * b3 - a3 * b2
  c2 = a3 * b1 - a1 * b3
  hjs1 = c1 / sqrt(c1**2 + c2**2)
  hjs2 = c2 / sqrt(c1**2 + c2**2)
  ! 外部境界法線ベクトル
  a1 =  0.5 * ( - x(i-1,je) + x(i+1,je) )
  a2 =  0.5 * ( - y(i-1,je) + y(i+1,je) )
  a3 =  0.0
  b1 =  0.0
  b2 =  0.0
  b3 = -1.0
  c1 = a2 * b3 - a3 * b2
  c2 = a3 * b1 - a1 * b3
  hje1 = c1 / sqrt(c1**2 + c2**2)
  hje2 = c2 / sqrt(c1**2 + c2**2)
  ! i - 方向
  if(i < i1) then
    dx1m(i,j) = x(i,je) - x(i-1,je)
    dy1m(i,j) = y(i,je) - y(i-1,je)
   else if(i <= i2) then
    dx1m(i,j) = ( -x(i-1,js) + x(i,js)  ) *   real(je-j) / real(je) &
    &         + ( -x(i-1,je) + x(i,je)  ) * ( real(j-js) / real(je) )**2
    dy1m(i,j) = ( -y(i-1,js) + y(i,js)  ) *   real(je-j) / real(je) &
    &         + ( -y(i-1,je) + y(i,je)  ) * ( real(j-js) / real(je) )**2
   else
    dx1m(i,j) = 0.0
    dy1m(i,j) = 0.0
  endif
  ! i + 方向
  if(i2 < i) then
    dx1p(i,j) = x(i+1,je) - x(i,je)
    dy1p(i,j) = y(i+1,je) - y(i,je)
   else if(i1 <= i) then
    dx1p(i,j) = ( -x(i,js) + x(i+1,js) ) *   real(je-j) / real(je) &
    &         + ( -x(i,je) + x(i+1,je) ) * ( real(j-js) / real(je) )**2
    dy1p(i,j) = ( -y(i,js) + y(i+1,js) ) *   real(je-j) / real(je) &
    &         + ( -y(i,je) + y(i+1,je) ) * ( real(j-js) / real(je) )**2
   else
    dx1p(i,j) = 0.0
    dy1p(i,j) = 0.0
  endif
  ! j - 方向
  dx2m(i,j) = hjs1 * rs
  dy2m(i,j) = hjs2 * rs
  ! j + 方向
  dx2p(i,j) = hje1 * re
  dy2p(i,j) = hje2 * re
 enddo
 enddo
 ! 楕円-双曲型重み関数 ---------------------------------------------------------------------------------
 do j = js+1, je-1
 do i = is+1, ie-1
  Cs(i,j) = 1.0 - max( cv1m(i,j), cv1p(i,j), cv2m(i,j), cv2p(i,j) ) + margin
 enddo
 enddo
 ! 格子生成計算ループ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do n = 1, nmax
  ! ソルバー部 -----------------------------------------------------------------------------------------
  call GridEHPDE2D( &
  &      omg, is, ie, js, je, &
  &      cs, cv1m, cv1p, cv2m, cv2p, dx1m, dy1m, dx1p, dy1p, dx2m, dy2m, dx2p, dy2p, &
  &      x, y, dmax )
  if( mod(n, 1000) == 0.0 ) write(*, "(a,i5,e16.8e3)") "* Elliptic-Hyperbolic calculation...", n, dmax
  if( dmax < resi ) exit
  ! 境界条件 -------------------------------------------------------------------------------------------
   do j = js + 1, je - 1
    a1 = x(is+1,j) - x(is+1,j-1)
    a2 = y(is+1,j) - y(is+1,j-1)
    x(is,j) = x(is,j-1) + a1
    y(is,j) = y(is,j-1) + a2
    a1 = x(ie-1,j) - x(ie-1,j-1)
    a2 = y(ie-1,j) - y(ie-1,j-1)
    x(ie,j) = x(ie,j-1) + a1
    y(ie,j) = y(ie,j-1) + a2
   enddo
 enddo
 if( dmax > resi ) then
  write(*, '(a)') "!!!!! Elliptic-Hyperbolic calculation error !!!!!"
  stop
 endif
 ! 処理終了 ********************************************************************************************
 return
end subroutine HtypeGenerationEHPDE
!*******************************************************************************************************
!******** 二境界法に基づく格子生成								********
!*******************************************************************************************************
subroutine GenerationTwoBoundary( &
&            is, ie, js, je, rs, re, t1, t2, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je
 real   , intent(in)    :: rs, re, t1, t2
 real   , intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: etabar(:)
 real    :: rr
 ! 処理開始 ********************************************************************************************
 ! メモリ確保
 allocate( etabar(js:je) )
 ! 媒介変数
 call GeometricInterpolationRemesh( rs / real(je-js+1), je-js+1, etabar, rr )
 ! 二境界法
 call TwoBoundaryMethod2D( &
 &      is, ie, js, je, t1, t2, etabar, x, y )
 ! メモリ解放
 deallocate(etabar)
 ! 処理終了 ********************************************************************************************
 return
end subroutine GenerationTwoBoundary
!***********************************************************************
!**** 等比級数による一次元補間関数                                  ****
!**** 初項a、公比rの等比数列の和(等比級数)により、                  ****
!**** 区間0 <= x <= 1にn個の点を配置する                            ****
!***********************************************************************
SUBROUTINE GeometricInterpolationRemesh(a, n, x, r)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: a
  INTEGER, INTENT(IN)  :: n
  REAL,    INTENT(OUT) :: x(n)
  REAL,    INTENT(OUT) :: r
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i
  LOGICAL :: Check
  ! 処理開始 ***********************************************************
  ! 公比を決定 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  r = FUNCR(n - 1, a, 1.0)
  !r is calculated with Mathematica FindRoot
  ! 等比級数による補間関数を計算 +++++++++++++++++++++++++++++++++++++++
  x(1) = 0.0
  x(n) = 1.0
  Check = .FALSE.
  IF(r .NE. 1.0) THEN
    ! 等比級数を計算 ---------------------------------------------------
    DO i = 2, n - 1
      x(i) = a * (1.0 - r**(i - 1)) / (1.0 - r)
    ENDDO
    ! 単調性を検査 -----------------------------------------------------
    DO i = 2, n - 1
      IF(x(i-1) .GE. x(i) .OR. x(i) .GE. x(i+1)) THEN
        WRITE(*, '(A)') 'GeometricInterpolation -> No Monotone Error'
        Check = .TRUE.
        EXIT
      ENDIF
    ENDDO
  ELSE
    Check = .TRUE.
  ENDIF
  ! 公比が1の場合か単調性が保証されていない場合は等間隔の結果を返す ++++
  IF(Check) THEN
    DO i = 2, n - 1
      x(i) = REAL(i - 1) / REAL(n - 1)
    ENDDO
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
! 定義終了 *************************************************************
END SUBROUTINE GeometricInterpolationRemesh
!*******************************************************************************************************
!******** 三次元化 (スパン方向変化なし)								********
!*******************************************************************************************************
subroutine ThreeDimensionalized( &
 &      is, ie, js, je, ks, ke, dom, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je, ks, ke
 real   , intent(in)    :: dom
 real   , intent(inout) :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke)
 real   , intent(out)   :: z(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 do k = ks, ke
 do j = js, je
 do i = is, ie
  x(i,j,k) = x(i,j,ks)
  y(i,j,k) = y(i,j,ks)
  z(i,j,k) = -0.5 * dom + dom * real(k) / real(ke)
 enddo
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine ThreeDimensionalized
!*******************************************************************************************************
!******** 着氷後の流れ場初期値									********
!*******************************************************************************************************
subroutine InitialFlow
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m
 integer :: i,j,k,l
 real,allocatable :: jac2(:,:,:)
 ! 処理開始 ********************************************************************************************
 ! 無次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call NondimensionalizedCoord3D( &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
! enddo
 ! メトリックス ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_Metrics3D( &
  &      trim(FlwIniDir) // trim(BlkName(m)) // trim(ND_MetFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%jac, &
  &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
  &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
  &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez )

 allocate(jac2(Flw(m)%is:Flw(m)%ie,Flw(m)%js:Flw(m)%je,Flw(m)%ks:Flw(m)%ke))
 do k = Flw(m)%ks,Flw(m)%ke
 do j = Flw(m)%js,Flw(m)%je
 do i = Flw(m)%is,Flw(m)%ie
  jac2(i,j,k) = Flw(m)%jac(i,j,k)
 end do
 end do
 end do

 m = me
  call Metrics3D( &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z, &
  &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
  &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
  &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
  &      Flw(m)%jac )

 do k = Flw(m)%ks,Flw(m)%ke
 do j = Flw(m)%js,Flw(m)%je
 do i = Flw(m)%is,Flw(m)%ie
 do l = ls,le
  Flw(m)%qh(i,j,k,l) = Flw(m)%qh(i,j,k,l) * jac2(i,j,k) / Flw(m)%jac(i,j,k)
 end do
 end do
 end do
 end do

 deallocate(jac2)

! enddo
 ! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 if(swi_IceStep .eq. 1) then
 ! 格子座標 --------------------------------------------------------------------------------------------
! do m = ms, me
 m = me
  call Output_Grid3D( &
  &      trim(FlwCalOutDir) // trim(BlkName(m)) // trim(ND_GrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
! enddo
 ! メトリックス ----------------------------------------------------------------------------------------
! do m = ms, me
 m = me
  call Output_Metrics3D( &
  &      trim(FlwCalOutDir) // trim(BlkName(m)) // trim(ND_MetFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%jac, &
  &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
  &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
  &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez )
! enddo
 ! 流束関数 --------------------------------------------------------------------------------------------
! do m = ms, me
 m = me
  call Output_Flux3D( &
  &      trim(FlwCalOutDir) // trim(BlkName(m)) // trim(ND_IniFlxFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      ls, le, Flw(m)%qh )
! enddo
 end if
 ! 計算条件ファイル出力 --------------------------------------------------------------------------------
 nCount = 0; nDrop = 0
 IceStep = IceStep + 1
 Span = Span / lRef
 if(swi_IceStep .eq. 1) then
  call Output_CalSetting( trim(ND_CalSetFile) // strtxt )
 end if
 ! 処理終了 ********************************************************************************************
 return
end subroutine InitialFlow
! 定義終了 *********************************************************************************************
end program Remesh_NACA
!*******************************************************************************************************
!******** MicroAVSファイル作成（三次元）  		                                        ********
!*******************************************************************************************************
subroutine MakeMAVSFile3D( &
&            strdir, strname, ext, is, ie, js, je, ks, ke, f, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strdir*(*), strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 integer  , intent(in)  :: f(is:ie, js:je, ks:ke)
 real     , intent(in)  :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke), z(is:ie, js:je, ks:ke)
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: strbin * 4 = '.bin'
 character, parameter :: strfld * 4 = '.fld'
 integer  , parameter :: ndim   = 3
 integer  , parameter :: nspace = 3
 integer  , parameter :: veclen = 1
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: nrkind, nskip
 real      :: r
 integer   :: i, j, k, n
 ! 処理開始 ********************************************************************************************
 ! ヘッダファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 open(1, file = trim(strdir) // trim(strname) // trim(strfld), form = 'formatted')
  write(1, '(a)')     '# AVS field file'
  write(1, '(a,i1)') 'ndim   = ', ndim
  write(1, '(a,i4)') 'dim1   = ', ie - is + 1
  write(1, '(a,i4)') 'dim2   = ', je - js + 1
  write(1, '(a,i4)') 'dim3   = ', ke - ks + 1
  write(1, '(a)')    'label  = flag'
  write(1, '(a,i1)') 'nspace = ', nspace
  write(1, '(a,i2)') 'veclen = ', veclen
  write(1, '(a)')    'data   = float'
  write(1, '(a)')    'field  = irregular'
  select case(ext)
 ! Binary ----------------------------------------------------------------------------------------------
   case(strbin)
    nrkind = kind(r)
    nskip  = nrkind * (ie - is + 1) * (je - js + 1) * (ke - ks + 1) + 8
    do n = 1, veclen
     write(1, '((a,i2), (x,2a), (x,a), (x,a,i11), 2(x,a,i1))') &
     & 'variable ', n, &
     & 'file = ', trim(strname) // ext, &
     & 'filetype = binary', &
     & 'skip = ', 4 + nskip * (n - 1), &
     & 'stride = ', 1, &
     & 'close = ', 1
    enddo
    do n = 1, nspace
     write(1, '((a,i1), (x,2a), (x,a), (x,a,i11), 2(x,a,i1))') &
     & 'coord ', n, &
     & 'file = ', trim(strname) // ext, &
     & 'filetype = binary', &
     & 'skip = ', 4 + nskip * (veclen + n - 1), &
     & 'stride = ', 1, &
     & 'close = ', 1
    enddo
 ! Aschii ----------------------------------------------------------------------------------------------
   case default
    do n = 1, veclen
     write(1, '((a,i2), (1x,2a), (1x,a), (1x,a,i1), 2(1x,a,i2), (1x,a,i1))') &
     & 'variable ', n, &
     & 'file = ', trim(strname) // ext, &
     & 'filetype = ascii', &
     & 'skip = ', 0, &
     & 'offset = ', n - 1, &
     & 'stride = ', veclen + nspace, &
     & 'close = ', 1
    enddo
    do n = 1, nspace
     write(1, '((a,i1), (1x,2a), (1x,a), (1x,a,i1), 2(1x,a,i2), (1x,a,i1))') &
     & 'coord ', n, &
     & 'file = ', trim(strname) // ext, &
     & 'filetype = ascii', &
     & 'skip = ', 0, &
     & 'offset = ', veclen + n - 1, &
     & 'stride = ', veclen + nspace, &
     & 'close = ', 1
    enddo
  end select
 close(1)
 ! データファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 select case(ext)
 ! Binary ----------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strdir) // trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) real(f); write(1) x; write(1) y; write(1) z
    close(1)
   close(1)
 ! Aschii ----------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strdir) // trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     write(1, '(e16.8e3, 3(x,e16.8e3))') real(f(i,j,k)), x(i,j,k), y(i,j,k), z(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine MakeMAVSFile3D

!*******************************************************************************************************
!******** vtkファイル出力サブルーチン 								********
!*******************************************************************************************************
subroutine OutputPara_bin( &
&      strdir, strname, is, ie, js, je, ks, ke, &
&      x, y, z )
 implicit none
 !mainroutine_variable
 character, intent(in)  :: strdir*(*), strname*(*)
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 real     , intent(in)  :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke), z(is:ie, js:je, ks:ke)
 !subroutine_variable
 character(len = 4), parameter		:: strvtk = '.vtk'
 character(len = 1), parameter		:: newline = char(10)
 character(len = 200)	:: strnum
 integer		:: i,j,k,n
 integer		:: ni,nj,nk
 integer		:: npoint

 npoint = (ie - is) * (je - js) * (ke - ks)

 open(unit       = 1,             &
 &    file       = trim(strdir) // trim(strname) // trim(strvtk),        &
 &    form       = 'unformatted',  &
 &    access     = 'sequential',   &
 &    convert    = 'big_endian',   &
 &    recordtype = 'stream',       &
 &    action     = 'write')
  write(1) '# vtk DataFile Version 3.0'//newline
  write(1) 'vtk output'//newline
  write(1) 'BINARY'//newline
  write(1) 'DATASET UNSTRUCTURED_GRID'//newline
  write(strnum,*) npoint * 8
  write(1) 'POINTS'//trim(strnum)//' float'//newline
  do n = 0,npoint-1
   ni = mod(mod(n, (ie - is) * (je - js)), (ie - is))
   nj = int(mod(n, (ie - is) * (je - js)) / (ie - is))
   nk = int(n / ((ie - is) * (je - js)))
   write(1) x(ni,nj,nk)		,y(ni,nj,nk)		,z(ni,nj,nk)
   write(1) x(ni+1,nj,nk)	,y(ni+1,nj,nk)		,z(ni+1,nj,nk)
   write(1) x(ni+1,nj+1,nk)	,y(ni+1,nj+1,nk)	,z(ni+1,nj+1,nk)
   write(1) x(ni,nj+1,nk)	,y(ni,nj+1,nk)		,z(ni,nj+1,nk)
   write(1) x(ni,nj,nk+1)	,y(ni,nj,nk+1)		,z(ni,nj,nk+1)
   write(1) x(ni+1,nj,nk+1)	,y(ni+1,nj,nk+1)	,z(ni+1,nj,nk+1)
   write(1) x(ni+1,nj+1,nk+1)	,y(ni+1,nj+1,nk+1)	,z(ni+1,nj+1,nk+1)
   write(1) x(ni,nj+1,nk+1)	,y(ni,nj+1,nk+1)	,z(ni,nj+1,nk+1)
  end do
  write(1) newline
  write(strnum,*) npoint, npoint * 9
  write(1) 'CELLS'//trim(strnum)//newline
  do n = 0,npoint-1
   write(1) 8, n * 8 + 0, n * 8 + 1, n * 8 + 2, n * 8 + 3, n * 8 + 4, n * 8 + 5, n * 8 + 6, n * 8 + 7
  end do
  write(1) newline
  write(strnum,*) npoint
  write(1) 'CELL_TYPES'//trim(strnum)//newline
  do n = 0, npoint-1
   write(1) 12
  end do
  write(1) newline
  write(strnum,*) npoint * 8
  write(1) 'POINT_DATA'//trim(strnum)//newline
  write(1) 'SCALARS n float'//newline
  write(1) 'LOOKUP_TABLE default'//newline
  do n = 0, npoint-1
   write(1) real(n * 8 + 0)
   write(1) real(n * 8 + 1)
   write(1) real(n * 8 + 2)
   write(1) real(n * 8 + 3)
   write(1) real(n * 8 + 4)
   write(1) real(n * 8 + 5)
   write(1) real(n * 8 + 6)
   write(1) real(n * 8 + 7)
  end do
 close(1)

end subroutine OutputPara_bin
