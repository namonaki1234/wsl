!*******************************************************************************************************
!*******************************************************************************************************
!******** 重合格子法の補間係数係数プログラム							********
!******** (NACA翼，三次元)									********
!********					      2013.02.02  PROGRAMMED BY RYOSUKE HAYASHI ********
!********					      2013.10.23    MODIFIED BY RYOSUKE HAYASHI	********
!*******************************************************************************************************
!*******************************************************************************************************
program OversetInterpolation_NACA
 ! モジュール宣言 **************************************************************************************
 use Package_NACA
 use Package_FileIO
 use Package_Grid
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! ファイル名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: OverCheckFile * 13 = 'OversetCheck'
 ! 処理開始 ********************************************************************************************
 ! 検証実験ケース選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(a)') "<< Exp. Case Selection >>"
 call SelectExpCase
 ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') '<< Initial Setting >>'
 call InitialSetting
 ! 探索範囲限定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Searching Area Limitation >>"
 call LimitSearchingArea
 ! 補間係数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! サブからメインの補間 --------------------------------------------------------------------------------
 write(*, '(/,a)') "<< Sub to Main Computation >>"
 call CalInterpolationSGtoMG
 ! メインからサブの補間 --------------------------------------------------------------------------------
 write(*, '(/,a)') "<< Main to Sub Computation >>"
 call CalInterpolationMGtoSG
 ! 着氷後の計算除去点探索 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Computational Removal Point >>"
 call SearchRemovalPoint
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
 ! ディレクトリ設定
 if( IceStep == 0 ) then
   GrdInDir = bckdir // 'grid//clean//'
   OSGDir   = bckdir // 'overset//clean//'
  else
   GrdInDir = bckdir // 'grid//icing//'
   OSGDir   = bckdir // 'overset//icing//'
 endif
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
 integer :: m
 real    :: zmax
 ! 処理開始 ********************************************************************************************
 ! ブロック名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Set_BlockName
 ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( Grd(ms:me) )
 ! メモリ確保及び格子ファイル ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do m = ms, me
  ! 格子解像度
  call Input_Resolution3D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(RslFile), strtxt, &
  &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke )
  ! メモリ確保
  allocate( Grd(m)%x(Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke), &
  &         Grd(m)%y(Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke), &
  &         Grd(m)%z(Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke) )
  ! 格子座標
  call Input_Grid3D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(GrdFile), strbin, &
  &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
  &      Grd(m)%x, Grd(m)%y, Grd(m)%z )
  ! 格子分割点
  call Input_CtypeGridPoint( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
  &      Grd(m)%i1, Grd(m)%i2, Grd(m)%i3 )
 enddo
 ! 重合格子構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( OSG(ms:me) )
 ! 処理終了 ********************************************************************************************
 return
end subroutine InitialSetting
!*******************************************************************************************************
!******** 探索範囲限定										********
!*******************************************************************************************************
subroutine LimitSearchingArea
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k, m, n
 ! 処理開始 ********************************************************************************************
 ! サブからメインの補間 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 m = 1; n = 2
 ! i 方向 ----------------------------------------------------------------------------------------------
 do i = Grd(m)%is, Grd(m)%ie
  if( minval( Grd(m)%x(i,:,:) ) < maxval( Grd(n)%x(:,:,:) )  ) exit
 enddo
 OSG(m)%is = i - 1
 do i = Grd(m)%ie, Grd(m)%is, - 1
  if( minval( Grd(m)%x(i,:,:) ) < maxval( Grd(n)%x(:,:,:) )  ) exit
 enddo
 OSG(m)%ie = i + 1
 ! j 方向 ----------------------------------------------------------------------------------------------
 OSG(m)%js = Grd(m)%js
 do j = Grd(m)%js, Grd(m)%je
  if( maxval( Grd(m)%y(:,j,:) ) > maxval( Grd(n)%y(:,:,:) )  ) exit
 enddo
 OSG(m)%je = j + 1
 ! k 方向 ----------------------------------------------------------------------------------------------
 OSG(m)%ks = Grd(m)%ks + 1
 OSG(m)%ke = Grd(m)%ke - 1
 ! メインからサブの補間 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 m = 2
 OSG(m)%is = Grd(m)%is    ; OSG(m)%ie = Grd(m)%ie
 OSG(m)%js = Grd(m)%js    ; OSG(m)%je = Grd(m)%je
 OSG(m)%ks = Grd(m)%ks + 1; OSG(m)%ke = Grd(m)%ke - 1
 ! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do m = ms, me
  call Output_Resolution3D( &
  &      trim(OSGDir) // trim(BlkName(m)) // trim(OverAreaFile), strtxt, &
  &      OSG(m)%is, OSG(m)%ie, OSG(m)%js, OSG(m)%je, OSG(m)%ks, OSG(m)%ke )
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine LimitSearchingArea
!*******************************************************************************************************
!******** 補間係数計算 (サブからメイン)								********
!*******************************************************************************************************
subroutine CalInterpolationSGtoMG
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, pointer :: flag(:, :, :)
 integer :: i, j, k, m, n
 integer :: i0, i1, i2, i3, i4, i5, i6, i7, i8, &
 &          j0, j1, j2, j3, j4, j5, j6, j7, j8, &
 &          k0, k1, k2, k3, k4, k5, k6, k7, k8
 integer :: ip, jp, kp
 integer :: isp, iep, jsp, jep, ksp, kep
 real    :: alp, bet, gam
 logical :: fSearch
 ! 処理開始 ********************************************************************************************
 m = 1; n = 2
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( OSG(m)%ip   ( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%jp   ( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%kp   ( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%fOver( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term1( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term2( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term3( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term4( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term5( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term6( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term7( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term8( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ) )
 OSG(m)%ip(:,:,:) = 0; OSG(m)%jp(:,:,:) = 0; OSG(m)%kp(:,:,:) = 0; OSG(m)%fOver(:,:,:) = 0
 OSG(m)%term1(:,:,:) = 0.0; OSG(m)%term2(:,:,:) = 0.0
 OSG(m)%term3(:,:,:) = 0.0; OSG(m)%term4(:,:,:) = 0.0
 OSG(m)%term5(:,:,:) = 0.0; OSG(m)%term6(:,:,:) = 0.0
 OSG(m)%term7(:,:,:) = 0.0; OSG(m)%term7(:,:,:) = 0.0
 ! 補間成功判定のフラグ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( flag( Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke ) )
 flag(:,:,:) = 0
 ! 補間係数許容誤差  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 OSG(m)%MGN = 0.0
 ! 補間係数計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !$omp parallel do default(shared) &
 !$                private(i, j, k, i0, j0, k0, &
 !$                        i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4, &
 !$                        i5, j5, k5, i6, j6, k6, i7, j7, k7, i8, j8, k8, &
 !$                        fSearch, alp, bet, gam)
 do k0 = OSG(m)%ks, OSG(m)%ke
 write(*, "(a,i1,a)") "* SG to MG : k = ", k0, " section calculation..."
 do j0 = OSG(m)%js, OSG(m)%je
 do i0 = OSG(m)%is, OSG(m)%ie
  ! 探索範囲を前の補間点の近くに限定 -------------------------------------------------------------------
!  if( fSearch ) then
!    isp = max(ip - SR, Grd(m)%is); iep = min(ip + SR, Grd(m)%ie - 1)
!    jsp = max(jp - SR, Grd(m)%js); jep = min(jp + SR, Grd(m)%je - 1)
!    ksp = max(kp - SR, Grd(m)%ks); kep = min(kp + SR, Grd(m)%ke - 1)
!   else
    isp = Grd(n)%is; iep = Grd(n)%ie - 1
    jsp = Grd(n)%js; jep = Grd(n)%je - 1
    ksp = Grd(n)%ks; kep = Grd(n)%ke - 1
!  endif
  ! 探索開始 -------------------------------------------------------------------------------------------
  do k = ksp, kep
  do j = jsp, jep
  do i = isp, iep
   ! メインから補間された点の除去
!   if( i == Grd(n)%is .or. i == Grd(n)%ie - 1 .or. &
!   &                       j == Grd(n)%je - 1 .or. &
!   &   k == Grd(n)%ks .or. k == Grd(n)%ke - 1 ) cycle
   if( i == Grd(n)%is .or. i == Grd(n)%ie - 1 .or. &
   &                       j == Grd(n)%je - 1 ) cycle
   ! ワイド・サーチ
   i1 = i    ; j1 = j    ; k1 = k    ; i2 = i + 1; j2 = j    ; k2 = k
   i3 = i    ; j3 = j + 1; k3 = k    ; i4 = i + 1; j4 = j + 1; k4 = k
   i5 = i    ; j5 = j    ; k5 = k + 1; i6 = i + 1; j6 = j    ; k6 = k + 1
   i7 = i    ; j7 = j + 1; k7 = k + 1; i8 = i + 1; j8 = j + 1; k8 = k + 1
   call WideSearch3D8Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &      Grd(n)%x(i5,j5,k5), Grd(n)%y(i5,j5,k5), Grd(n)%z(i5,j5,k5), &
   &      Grd(n)%x(i6,j6,k6), Grd(n)%y(i6,j6,k6), Grd(n)%z(i6,j6,k6), &
   &      Grd(n)%x(i7,j7,k7), Grd(n)%y(i7,j7,k7), Grd(n)%z(i7,j7,k7), &
   &      Grd(n)%x(i8,j8,k8), Grd(n)%y(i8,j8,k8), Grd(n)%z(i8,j8,k8), &
   &      OSG(m)%MGN, fSearch )
   if(.not. fSearch) cycle
   ! 三重線形補間
   alp = 0.5; bet = 0.5; gam = 0.5
   call Interpolation3D8Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &      Grd(n)%x(i5,j5,k5), Grd(n)%y(i5,j5,k5), Grd(n)%z(i5,j5,k5), &
   &      Grd(n)%x(i6,j6,k6), Grd(n)%y(i6,j6,k6), Grd(n)%z(i6,j6,k6), &
   &      Grd(n)%x(i7,j7,k7), Grd(n)%y(i7,j7,k7), Grd(n)%z(i7,j7,k7), &
   &      Grd(n)%x(i8,j8,k8), Grd(n)%y(i8,j8,k8), Grd(n)%z(i8,j8,k8), &
   &  	  OSG(m)%MGN, alp, bet, gam, fSearch )
   if(fSearch) then
     OSG(m)%fOver(i0,j0,k0) = 1
     exit
   endif
  enddo; if(fSearch) exit
  enddo; if(fSearch) exit
  enddo
  ! 補間係数 -------------------------------------------------------------------------------------------
  if(fSearch) then
    OSG(m)%ip(i0,j0,k0) = i; OSG(m)%jp(i0,j0,k0) = j; OSG(m)%kp(i0,j0,k0) = k
    OSG(m)%term1(i0,j0,k0) = (1.0 - alp) * (1.0 - bet) * (1.0 - gam)
    OSG(m)%term2(i0,j0,k0) =        alp  * (1.0 - bet) * (1.0 - gam)
    OSG(m)%term3(i0,j0,k0) = (1.0 - alp) *        bet  * (1.0 - gam)
    OSG(m)%term4(i0,j0,k0) =        alp  *        bet  * (1.0 - gam)
    OSG(m)%term5(i0,j0,k0) = (1.0 - alp) * (1.0 - bet) *        gam
    OSG(m)%term6(i0,j0,k0) =        alp  * (1.0 - bet) *        gam
    OSG(m)%term7(i0,j0,k0) = (1.0 - alp) *        bet  *        gam
    OSG(m)%term8(i0,j0,k0) =        alp  *        bet  *        gam
    flag(i0,j0,k0) = 1
!    ip = i; jp = j; kp = k
    cycle
  endif
  ! 探索開始 -------------------------------------------------------------------------------------------
  do k = ksp, kep
  do j = jsp, jep
  do i = isp, iep
   ! メインから補間された点の除去
!   if( i == Grd(n)%is .or. i == Grd(n)%ie - 1 .or. &
!   &                       j == Grd(n)%je - 1 .or. &
!   &   k == Grd(n)%ks .or. k == Grd(n)%ke - 1 ) cycle
   if( i == Grd(n)%is .or. i == Grd(n)%ie - 1 .or. &
   &                       j == Grd(n)%je - 1 ) cycle
   ! ワイド・サーチ
   i1 = i    ; j1 = j    ; k1 = k    ; i2 = i + 1; j2 = j    ; k2 = k
   i3 = i    ; j3 = j + 1; k3 = k    ; i4 = i + 1; j4 = j + 1; k4 = k
   i5 = i    ; j5 = j    ; k5 = k + 1; i6 = i + 1; j6 = j    ; k6 = k + 1
   i7 = i    ; j7 = j + 1; k7 = k + 1; i8 = i + 1; j8 = j + 1; k8 = k + 1
   call WideSearch3D8Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &      Grd(n)%x(i5,j5,k5), Grd(n)%y(i5,j5,k5), Grd(n)%z(i5,j5,k5), &
   &      Grd(n)%x(i6,j6,k6), Grd(n)%y(i6,j6,k6), Grd(n)%z(i6,j6,k6), &
   &      Grd(n)%x(i7,j7,k7), Grd(n)%y(i7,j7,k7), Grd(n)%z(i7,j7,k7), &
   &      Grd(n)%x(i8,j8,k8), Grd(n)%y(i8,j8,k8), Grd(n)%z(i8,j8,k8), &
   &      OSG(m)%MGN, fSearch )
   if(.not. fSearch) cycle
   ! 三次元線形補間 (Pattern-1)
   i1 = i    ; j1 = j    ; k1 = k    ; i2 = i    ; j2 = j    ; k2 = k + 1
   i3 = i + 1; j3 = j    ; k3 = k + 1; i4 = i    ; j4 = j + 1; k4 = k + 1
   alp = 0.5; bet = 0.5; gam = 0.5
   call Interpolation3D4Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &	  OSG(m)%MGN, alp, bet, gam, fSearch )
   if(fSearch) then
     OSG(m)%fOver(i0,j0,k0) = 2
     exit
   endif
   ! 三次元線形補間 (Pattern-2)
   i1 = i    ; j1 = j    ; k1 = k    ; i2 = i    ; j2 = j + 1; k2 = k    
   i3 = i + 1; j3 = j + 1; k3 = k    ; i4 = i    ; j4 = j + 1; k4 = k + 1
   alp = 0.5; bet = 0.5; gam = 0.5
   call Interpolation3D4Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &	  OSG(m)%MGN, alp, bet, gam, fSearch )
   if(fSearch) then
     OSG(m)%fOver(i0,j0,k0) = 3
     exit
   endif
   ! 三次元線形補間 (Pattern-3)
   i1 = i    ; j1 = j    ; k1 = k    ; i2 = i + 1; j2 = j + 1; k2 = k    
   i3 = i + 1; j3 = j    ; k3 = k + 1; i4 = i    ; j4 = j + 1; k4 = k + 1
   alp = 0.5; bet = 0.5; gam = 0.5
   call Interpolation3D4Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &	  OSG(m)%MGN, alp, bet, gam, fSearch )
   if(fSearch) then
     OSG(m)%fOver(i0,j0,k0) = 4
     exit
   endif
   ! 三次元線形補間 (Pattern-4)
   i1 = i    ; j1 = j    ; k1 = k    ; i2 = i + 1; j2 = j    ; k2 = k    
   i3 = i + 1; j3 = j + 1; k3 = k    ; i4 = i + 1; j4 = j    ; k4 = k + 1
   alp = 0.5; bet = 0.5; gam = 0.5
   call Interpolation3D4Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &	  OSG(m)%MGN, alp, bet, gam, fSearch )
   if(fSearch) then
     OSG(m)%fOver(i0,j0,k0) = 5
     exit
   endif
   ! 三次元線形補間 (Pattern-5)
   i1 = i + 1; j1 = j + 1; k1 = k    ; i2 = i + 1; j2 = j    ; k2 = k + 1
   i3 = i    ; j3 = j + 1; k3 = k + 1; i4 = i + 1; j4 = j + 1; k4 = k + 1
   alp = 0.5; bet = 0.5; gam = 0.5
   call Interpolation3D4Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &	  OSG(m)%MGN, alp, bet, gam, fSearch )
   if(fSearch) then
     OSG(m)%fOver(i0,j0,k0) = 6
     exit
   endif
  enddo; if(fSearch) exit
  enddo; if(fSearch) exit
  enddo
  ! 補間係数決定 --------------------------------------------------------------------------------------
  if(fSearch) then
    OSG(m)%ip(i0,j0,k0) = i; OSG(m)%jp(i0,j0,k0) = j; OSG(m)%kp(i0,j0,k0) = k
    OSG(m)%term1(i0,j0,k0) = 1.0
    OSG(m)%term2(i0,j0,k0) = alp
    OSG(m)%term3(i0,j0,k0) = bet
    OSG(m)%term4(i0,j0,k0) = gam
    flag(i0,j0,k0) = 1
!    ip = i; jp = j; kp = k
  endif
 enddo
 enddo
 enddo
 !$omp end parallel do
 ! ファイル出力 ---------------------------------------------------------------------------------------
 ! 補間係数
 call Output_OversetCoe3D( &
 &      trim(OSGDir) // trim(BlkName(m)) // trim(OverCoeFile), strbin, &
 &      OSG(m)%is, OSG(m)%ie, OSG(m)%js, OSG(m)%je, OSG(m)%ks, OSG(m)%ke, &
 &      OSG(m)%fover, OSG(m)%ip, OSG(m)%jp, OSG(m)%kp, &
 &      OSG(m)%term1, OSG(m)%term2, OSG(m)%term3, OSG(m)%term4, &
 &      OSG(m)%term5, OSG(m)%term6, OSG(m)%term7, OSG(m)%term8 )
 ! 補間チェック
 call MakeMAVSFile3D( &
 &      trim(OSGDir), trim(BlkName(m)) // trim(OverCheckFile), strbin, &
 &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
 &      flag , Grd(m)%x , Grd(m)%y , Grd(m)%z )
 deallocate(flag)
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalInterpolationSGtoMG
!*******************************************************************************************************
!******** 補間係数計算 (メインからサブ)								********
!*******************************************************************************************************
subroutine CalInterpolationMGtoSG
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, pointer :: flag(:, :, :)
 integer :: i, j, k, m, n
 integer :: i0, i1, i2, i3, i4, i5, i6, i7, i8, &
 &          j0, j1, j2, j3, j4, j5, j6, j7, j8, &
 &          k0, k1, k2, k3, k4, k5, k6, k7, k8
 integer :: ip, jp, kp
 integer :: isp, iep, jsp, jep, ksp, kep
 real    :: alp, bet, gam
 logical :: fSearch
 ! 処理開始 ********************************************************************************************
 m = 2; n = 1
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( OSG(m)%ip   ( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%jp   ( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%kp   ( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%fOver( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term1( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term2( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term3( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term4( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term5( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term6( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term7( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
 &         OSG(m)%term8( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ) )
 OSG(m)%ip(:,:,:) = 0; OSG(m)%jp(:,:,:) = 0; OSG(m)%kp(:,:,:) = 0; OSG(m)%fOver(:,:,:) = 0
 OSG(m)%term1(:,:,:) = 0.0; OSG(m)%term2(:,:,:) = 0.0
 OSG(m)%term3(:,:,:) = 0.0; OSG(m)%term4(:,:,:) = 0.0
 OSG(m)%term5(:,:,:) = 0.0; OSG(m)%term6(:,:,:) = 0.0
 OSG(m)%term7(:,:,:) = 0.0; OSG(m)%term7(:,:,:) = 0.0
 ! 補間のフラグ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( flag( Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke ) )
 flag(:,:,:) = 0
 ! 補間係数許容誤差 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 OSG(m)%MGN = 0.0
 ! 補間係数計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !$omp parallel do default(shared) &
 !$                private(i, j, k, i0, j0, k0, &
 !$                        i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4, &
 !$                        i5, j5, k5, i6, j6, k6, i7, j7, k7, i8, j8, k8, &
 !$                        fSearch, alp, bet, gam)
 do k0 = OSG(m)%ks, OSG(m)%ke
 write(*, "(a,i1,a)") "* MG to SG : k = ", k0, " section calculation..."
 do j0 = OSG(m)%js, OSG(m)%je
 do i0 = OSG(m)%is, OSG(m)%ie
  ! 補間点される点を外部領域に限定 ---------------------------------------------------------------------
!  if( (i0 /= Grd(m)%is .and. i0 /= Grd(m)%ie) .and. &
!      (                      j0 /= Grd(m)%je) .and. &
!      (k0 /= Grd(m)%ks .and. k0 /= Grd(m)%ke) ) cycle
  if( (i0 /= Grd(m)%is .and. i0 /= Grd(m)%ie) .and. &
      (                      j0 /= Grd(m)%je) ) cycle
  ! 探索範囲を前の補間点の近くに限定 -------------------------------------------------------------------
!  if( fSearch ) then
!    isp = max(ip - SR, Grd(m)%is); iep = min(ip + SR, Grd(m)%ie - 1)
!    jsp = max(jp - SR, Grd(m)%js); jep = min(jp + SR, Grd(m)%je - 1)
!    ksp = max(kp - SR, Grd(m)%ks); kep = min(kp + SR, Grd(m)%ke - 1)
!   else
    isp = Grd(n)%is; iep = Grd(n)%ie - 1
    jsp = Grd(n)%js; jep = Grd(n)%je - 1
    ksp = Grd(n)%ks; kep = Grd(n)%ke - 1
!  endif
  ! 探索開始 ------------------------------------------------------------------------------------------
  do k = ksp, kep
  do j = jsp, jep
  do i = isp, iep
   ! ワイド・サーチ
   i1 = i    ; j1 = j    ; k1 = k    ; i2 = i + 1; j2 = j    ; k2 = k
   i3 = i    ; j3 = j + 1; k3 = k    ; i4 = i + 1; j4 = j + 1; k4 = k
   i5 = i    ; j5 = j    ; k5 = k + 1; i6 = i + 1; j6 = j    ; k6 = k + 1
   i7 = i    ; j7 = j + 1; k7 = k + 1; i8 = i + 1; j8 = j + 1; k8 = k + 1
   call WideSearch3D8Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &      Grd(n)%x(i5,j5,k5), Grd(n)%y(i5,j5,k5), Grd(n)%z(i5,j5,k5), &
   &      Grd(n)%x(i6,j6,k6), Grd(n)%y(i6,j6,k6), Grd(n)%z(i6,j6,k6), &
   &      Grd(n)%x(i7,j7,k7), Grd(n)%y(i7,j7,k7), Grd(n)%z(i7,j7,k7), &
   &      Grd(n)%x(i8,j8,k8), Grd(n)%y(i8,j8,k8), Grd(n)%z(i8,j8,k8), &
   &      OSG(m)%MGN, fSearch )
   if(.not. fSearch) cycle
   ! 三重線形補間
   alp = 0.5; bet = 0.5; gam = 0.5
   call Interpolation3D8Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &      Grd(n)%x(i5,j5,k5), Grd(n)%y(i5,j5,k5), Grd(n)%z(i5,j5,k5), &
   &      Grd(n)%x(i6,j6,k6), Grd(n)%y(i6,j6,k6), Grd(n)%z(i6,j6,k6), &
   &      Grd(n)%x(i7,j7,k7), Grd(n)%y(i7,j7,k7), Grd(n)%z(i7,j7,k7), &
   &      Grd(n)%x(i8,j8,k8), Grd(n)%y(i8,j8,k8), Grd(n)%z(i8,j8,k8), &
   &  	  OSG(m)%MGN, alp, bet, gam, fSearch )
   if(fSearch) then
     OSG(m)%fOver(i0,j0,k0) = 1
     exit
   endif
  enddo; if(fSearch) exit
  enddo; if(fSearch) exit
  enddo
  ! 補間係数 -------------------------------------------------------------------------------------------
  if(fSearch) then
    OSG(m)%ip(i0,j0,k0) = i; OSG(m)%jp(i0,j0,k0) = j; OSG(m)%kp(i0,j0,k0) = k
    OSG(m)%term1(i0,j0,k0) = (1.0 - alp) * (1.0 - bet) * (1.0 - gam)
    OSG(m)%term2(i0,j0,k0) =        alp  * (1.0 - bet) * (1.0 - gam)
    OSG(m)%term3(i0,j0,k0) = (1.0 - alp) *        bet  * (1.0 - gam)
    OSG(m)%term4(i0,j0,k0) =        alp  *        bet  * (1.0 - gam)
    OSG(m)%term5(i0,j0,k0) = (1.0 - alp) * (1.0 - bet) *        gam
    OSG(m)%term6(i0,j0,k0) =        alp  * (1.0 - bet) *        gam
    OSG(m)%term7(i0,j0,k0) = (1.0 - alp) *        bet  *        gam
    OSG(m)%term8(i0,j0,k0) =        alp  *        bet  *        gam
    flag(i0,j0,k0) = 1
!    ip = i; jp = j; kp = k
    cycle
  endif
  ! 探索開始 ------------------------------------------------------------------------------------------
  do k = ksp, kep
  do j = jsp, jep
  do i = isp, iep
   ! ワイド・サーチ
   i1 = i    ; j1 = j    ; k1 = k    ; i2 = i + 1; j2 = j    ; k2 = k
   i3 = i    ; j3 = j + 1; k3 = k    ; i4 = i + 1; j4 = j + 1; k4 = k
   i5 = i    ; j5 = j    ; k5 = k + 1; i6 = i + 1; j6 = j    ; k6 = k + 1
   i7 = i    ; j7 = j + 1; k7 = k + 1; i8 = i + 1; j8 = j + 1; k8 = k + 1
   call WideSearch3D8Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &      Grd(n)%x(i5,j5,k5), Grd(n)%y(i5,j5,k5), Grd(n)%z(i5,j5,k5), &
   &      Grd(n)%x(i6,j6,k6), Grd(n)%y(i6,j6,k6), Grd(n)%z(i6,j6,k6), &
   &      Grd(n)%x(i7,j7,k7), Grd(n)%y(i7,j7,k7), Grd(n)%z(i7,j7,k7), &
   &      Grd(n)%x(i8,j8,k8), Grd(n)%y(i8,j8,k8), Grd(n)%z(i8,j8,k8), &
   &      OSG(m)%MGN, fSearch )
   if(.not. fSearch) cycle
   ! 三次元線形補間 (Pattern-1)
   i1 = i    ; j1 = j    ; k1 = k    ; i2 = i    ; j2 = j    ; k2 = k + 1
   i3 = i + 1; j3 = j    ; k3 = k + 1; i4 = i    ; j4 = j + 1; k4 = k + 1
   alp = 0.5; bet = 0.5; gam = 0.5
   call Interpolation3D4Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &	  OSG(m)%MGN, alp, bet, gam, fSearch )
   if(fSearch) then
     OSG(m)%fOver(i0,j0,k0) = 2
     exit
   endif
   ! 三次元線形補間 (Pattern-2)
   i1 = i    ; j1 = j    ; k1 = k    ; i2 = i    ; j2 = j + 1; k2 = k    
   i3 = i + 1; j3 = j + 1; k3 = k    ; i4 = i    ; j4 = j + 1; k4 = k + 1
   alp = 0.5; bet = 0.5; gam = 0.5
   call Interpolation3D4Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &	  OSG(m)%MGN, alp, bet, gam, fSearch )
   if(fSearch) then
     OSG(m)%fOver(i0,j0,k0) = 3
     exit
   endif
   ! 三次元線形補間 (Pattern-3)
   i1 = i    ; j1 = j    ; k1 = k    ; i2 = i + 1; j2 = j + 1; k2 = k    
   i3 = i + 1; j3 = j    ; k3 = k + 1; i4 = i    ; j4 = j + 1; k4 = k + 1
   alp = 0.5; bet = 0.5; gam = 0.5
   call Interpolation3D4Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &	  OSG(m)%MGN, alp, bet, gam, fSearch )
   if(fSearch) then
     OSG(m)%fOver(i0,j0,k0) = 4
     exit
   endif
   ! 三次元線形補間 (Pattern-4)
   i1 = i    ; j1 = j    ; k1 = k    ; i2 = i + 1; j2 = j    ; k2 = k    
   i3 = i + 1; j3 = j + 1; k3 = k    ; i4 = i + 1; j4 = j    ; k4 = k + 1
   alp = 0.5; bet = 0.5; gam = 0.5
   call Interpolation3D4Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &	  OSG(m)%MGN, alp, bet, gam, fSearch )
   if(fSearch) then
     OSG(m)%fOver(i0,j0,k0) = 5
     exit
   endif
   ! 三次元線形補間 (Pattern-5)
   i1 = i + 1; j1 = j + 1; k1 = k    ; i2 = i + 1; j2 = j    ; k2 = k + 1
   i3 = i    ; j3 = j + 1; k3 = k + 1; i4 = i + 1; j4 = j + 1; k4 = k + 1
   alp = 0.5; bet = 0.5; gam = 0.5
   call Interpolation3D4Point( &
   &      Grd(m)%x(i0,j0,k0), Grd(m)%y(i0,j0,k0), Grd(m)%z(i0,j0,k0), &
   &      Grd(n)%x(i1,j1,k1), Grd(n)%y(i1,j1,k1), Grd(n)%z(i1,j1,k1), &
   &      Grd(n)%x(i2,j2,k2), Grd(n)%y(i2,j2,k2), Grd(n)%z(i2,j2,k2), &
   &      Grd(n)%x(i3,j3,k3), Grd(n)%y(i3,j3,k3), Grd(n)%z(i3,j3,k3), &
   &      Grd(n)%x(i4,j4,k4), Grd(n)%y(i4,j4,k4), Grd(n)%z(i4,j4,k4), &
   &	  OSG(m)%MGN, alp, bet, gam, fSearch )
   if(fSearch) then
     OSG(m)%fOver(i0,j0,k0) = 6
     exit
   endif
  enddo; if(fSearch) exit
  enddo; if(fSearch) exit
  enddo
  ! 補間係数決定 --------------------------------------------------------------------------------------
  if(fSearch) then
    OSG(m)%ip(i0,j0,k0) = i; OSG(m)%jp(i0,j0,k0) = j; OSG(m)%kp(i0,j0,k0) = k
    OSG(m)%term1(i0,j0,k0) = 1.0
    OSG(m)%term2(i0,j0,k0) = alp
    OSG(m)%term3(i0,j0,k0) = bet
    OSG(m)%term4(i0,j0,k0) = gam
    flag(i0,j0,k0) = 1
!    ip = i; jp = j; kp = k
    else
     write(*, '(a)') '!!!! Error : Interpolation !!!!!'
     write(*, '(a,3(i3,a))') '@ (i,j,k) = (', i0, ',', j0, ',', k0, ')'
  endif
 enddo
 enddo
 enddo
 !$omp end parallel do
 ! ファイル出力 ---------------------------------------------------------------------------------------
 ! 補間係数
 call Output_OversetCoe3D( &
 &      trim(OSGDir) // trim(BlkName(m)) // trim(OverCoeFile), strbin, &
 &      OSG(m)%is, OSG(m)%ie, OSG(m)%js, OSG(m)%je, OSG(m)%ks, OSG(m)%ke, &
 &      OSG(m)%fover, OSG(m)%ip, OSG(m)%jp, OSG(m)%kp, &
 &      OSG(m)%term1, OSG(m)%term2, OSG(m)%term3, OSG(m)%term4, &
 &      OSG(m)%term5, OSG(m)%term6, OSG(m)%term7, OSG(m)%term8 )
 ! 補間チェック
 call MakeMAVSFile3D( &
 &      trim(OSGDir), trim(BlkName(m)) // trim(OverCheckFile), strbin, &
 &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
 &      flag , Grd(m)%x , Grd(m)%y , Grd(m)%z )
 deallocate(flag)
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalInterpolationMGtoSG
!*******************************************************************************************************
!******** 着氷後の計算除去点探索								********
!*******************************************************************************************************
subroutine SearchRemovalPoint
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, m
 integer :: kp
 ! 処理開始 ********************************************************************************************
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 m = 1
 allocate( IceIn( Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke ) )
 IceIn(:, :, :) = 0
 kp = nint( 0.5 * Grd(m)%ke )
 ! 氷の中の点探索 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 if( IceStep == 0 ) return
 do j = OSG(m)%js     , OSG(m)%je - 10
 do i = OSG(m)%is + 10, OSG(m)%ie - 10
  if( OSG(m)%fOver(i,j,kp) == 0 ) then
    IceIn(i,j,:) = 1
  endif
 enddo
 enddo
 ! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Output_ArrayInt3D( &
 &      trim(OSGDir) // trim(BlkName(m)) // trim(IceInPointFile), strbin, &
 &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
 &      IceIn )
 ! 処理終了 ********************************************************************************************
 return
end subroutine SearchRemovalPoint
! 定義終了 *********************************************************************************************
end program OversetInterpolation_NACA
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
