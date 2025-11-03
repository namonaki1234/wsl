!*******************************************************************************************************
!*******************************************************************************************************
!******** 液滴軌道計算プログラム								********
!******** (NACA翼，三次元圧縮性乱流場，重合格子法，液滴径分布なし)				********
!********					      2013.06.18  PROGRAMMED BY RYOSUKE HAYASHI ********
!********					      2013.07.18     UPDATED BY RYOSUKE HAYASHI ********
!********					      2014.04.15     UPDATED BY RYOSUKE HAYASHI ********
!******** 液滴変形モデル			      2014.12.08     UPDATED BY MIKI	SHIMURA	********
!********					      2015.04.29     UPDATED BY MIKI	SHIMURA ********
!********					      2015.05.04     UPDATED BY MIKI	SHIMURA	********
!******** 二次液滴の衝突判定入り		      2015.07.27     UPDATED BY MIKI    SHIMURA ********
!*******************************************************************************************************
!*******************************************************************************************************
program DropletTrajectory_NACA
 ! モジュール宣言 **************************************************************************************
 use Package_NACA
 use Package_FileIO
 use Package_Droplet
 use mt95
 use Package_Grid
 use Package_Flow
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 共有定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, parameter :: mini     =   me !1					! 液滴初期ブロック
 integer, parameter :: CalNum   =   3					! 1 : 衝突　2 : 可視化　3 : 両方
 integer, parameter :: nDrpView = 500					! 可視化液滴数
 real   , parameter :: SplLim   = 40.0e-6				! スプラッシュ判定
 ! 共有変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 logical, pointer :: fOver(:, :, :)					! 重合格子部のフラグ
 integer :: DrpNum							! 液滴ファイル番号
 real    :: xmin, ymin, ymax, zmin, zmax				! 液滴投入範囲
 real    :: uini, vini, wini						! 液滴初期速度
 real    :: xp, yp, zp                              			! 液滴位置 (物理空間)
 real    :: up, vp, wp                         				! 液滴速度 (物理空間)
 real    :: xip, etp, zep                           			! 液滴位置 (計算空間)
 real    :: uxp, vep, wzp                           			! 液滴速度 (計算空間)
 real    :: ypini							! 液滴初期位置 (物理空間)
 real    :: dp								! 液滴直径
 real    :: fx, fy, fz                              			! 液滴に働く力
 real    :: Ohw, Rew
 real    :: Cd								! 抗力係数
 real    :: AsRatio							! 液滴のアスペクト比
 real    :: ImpMass							! 液滴の衝突質量
 real    :: Time0							! 時間 (時間ごとの様子を見たい時)
 integer :: Nlost							! 液滴ロスト数
 logical :: fSearch							! 液滴探索成功フラグ
 logical :: fSplash							! 液滴スプラッシュフラグ
 logical :: fBounce							! 液滴バウンドフラグ
 logical :: fImpi1							! 1回目の液滴衝突フラグ
 logical :: fImpi2							! 2回目の液滴衝突フラグ
 logical :: fBreakup							! 液滴分裂フラグ
 ! 処理開始 ********************************************************************************************
 ! 検証実験ケース選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(a)'  ) "<< Exp. Case Selection >>"
 call SelectExpCase
 ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Initial Condition Set >>"
 call InitialSetting
 ! 液滴軌道計算初期条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Inlet Droplet Set >>"
 call InitialDropletCondition
 ! 衝突面積計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Impingemenr Area Set >>"
 call CalImpingementArea
 ! 液滴軌道計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Droplet Trajectory calculation Computation >>"
 call CalDropletTrajectory
 ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Memory Deallocation >>"
 call Deallocating
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
   GrdInDir    = bckdir // 'grid/clean/'
   OSGDir      = bckdir // 'overset/clean/'
   FlwIniDir   = bckdir // 'flow/initial/clean/'
   FlwCalInDir = bckdir // 'flow/cal/clean/'
   DrpImpDir   = bckdir // 'droplet/impingement/clean/'
   DrpTraDir   = bckdir // 'droplet/trajectory/clean/'
  else
   GrdInDir    = bckdir // 'grid/icing/'
   OSGDir      = bckdir // 'overset/icing/'
   FlwIniDir   = bckdir // 'flow/initial//icing/'
   FlwCalInDir = bckdir // 'flow/cal/icing/'
   DrpImpDir   = bckdir // 'droplet/impingement/icing/'
   DrpTraDir   = bckdir // 'droplet/trajectory/icing/'
 endif
 write(*, '(a,i2)') '* Ice step      = ', IceStep
 write(*, '(a,i2)') '* Ice step max. = ', IceStepMax
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
!******** 初期設定										********
!*******************************************************************************************************
subroutine InitialSetting
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, j, k, m
 integer   :: jp
 character :: fname * 30
 ! 処理開始 ********************************************************************************************
 ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! allocate( Flw(ms:me), OSG(ms:me), Drp(1:nDrpFile), Ice(ms:me) )
 allocate( Flw(ms:me), Drp(1:nDrpFile), Ice(ms:me) )
 ! ブロック名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Set_BlockName
 ! 格子解像度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_Resolution3D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(RslFile), strtxt, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke )
! enddo
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  allocate( Flw(m)%x   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%y   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%z   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%rho (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%u   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%v   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%w   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%p   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%t   (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%mu  (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%kin (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%eps (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%xix (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%xiy (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%xiz (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%etx (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%ety (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%etz (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%zex (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%zey (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%zez (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%jac (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%qh  (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le) )
! enddo
 ! 格子座標 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_Grid3D( &
  &      trim(FlwIniDir) // trim(BlkName(m)) // trim(ND_GrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
! enddo
 ! C 型格子分割点 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! m = 1
 m = me
 call Input_CtypeGridPoint( &
 &      trim(GrdInDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
 &      Flw(m)%i1, Flw(m)%i2, Flw(m)%i3 )
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
! enddo
 ! 流束関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_Flux3D( &
  &      trim(FlwCalInDir) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      ls, le, Flw(m)%qh )
!  call Input_Flux3D( &
!  &      trim(FlwCalInDir) // trim(fname) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
!  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
!  &      ls, le, Flw(m)%qh )
! enddo
 ! 物理量 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call SetPhysics3DKEM( &
  &      Rg, gamma, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
  &      Flw(m)%qh, Flw(m)%jac, &
  &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps )
  call ViscosityCoefficient3D( &
  &      muSth, TsSth, s1, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%t, &
  &      Flw(m)%mu )
  Flw(m)%w(:,:,:) = 0.0
! enddo
! ! 重合格子補間係数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
!  ! 補間探索範囲 ---------------------------------------------------------------------------------------
!  call Input_Resolution3D( &
!  &      trim(OSGDir) // trim(BlkName(m)) // trim(OverAreaFile), strtxt, &
!  &      OSG(m)%is, OSG(m)%ie, OSG(m)%js, OSG(m)%je, OSG(m)%ks, OSG(m)%ke )
!  ! メモリ確保 -----------------------------------------------------------------------------------------
!  allocate( OSG(m)%ip   (OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
!  &         OSG(m)%jp   (OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
!  &         OSG(m)%kp   (OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
!  &         OSG(m)%fOver(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
!  &         OSG(m)%term1(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
!  &         OSG(m)%term2(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
!  &         OSG(m)%term3(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
!  &         OSG(m)%term4(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
!  &         OSG(m)%term5(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
!  &         OSG(m)%term6(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
!  &         OSG(m)%term7(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke), &
!  &         OSG(m)%term8(OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke) )
!  ! 補間係数 -------------------------------------------------------------------------------------------
!  call Input_OversetCoe3D( &
!  &      trim(OSGDir) // trim(BlkName(m)) // trim(OverCoeFile), strbin, &
!  &      OSG(m)%is, OSG(m)%ie, OSG(m)%js, OSG(m)%je, OSG(m)%ks, OSG(m)%ke, &
!  &      OSG(m)%fOver, OSG(m)%ip, OSG(m)%jp, OSG(m)%kp, &
!  &      OSG(m)%term1, OSG(m)%term2, OSG(m)%term3, OSG(m)%term4, &
!  &      OSG(m)%term5, OSG(m)%term6, OSG(m)%term7, OSG(m)%term8 )
! enddo
! ! 重合格子部のフラグ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! m = 1
! allocate( fOver(Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke) )
! fOver(:,:,:) = .false.
! do k = OSG(m)%ks, OSG(m)%ke
! do j = OSG(m)%js, OSG(m)%je
! do i = OSG(m)%is, OSG(m)%ie
!  if( OSG(m)%fOver(i,j,k) == 0 ) cycle
!  fOver(i,j,k) = .true.
! enddo
! enddo
! enddo
! ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
!  deallocate( Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps, Flw(m)%jac, Flw(m)%qh )
!  deallocate( OSG(m)%fOver, OSG(m)%ip   , OSG(m)%jp   , OSG(m)%kp   , &
!  &           OSG(m)%term1, OSG(m)%term2, OSG(m)%term3, OSG(m)%term4, &
!  &           OSG(m)%term5, OSG(m)%term6, OSG(m)%term7, OSG(m)%term8 )
! enddo
 ! 液滴軌道 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do m = 1, nDrpFile
  allocate( Drp(m)%m(1:nCalMax), &
  &         Drp(m)%x(1:nCalMax), Drp(m)%y(1:nCalMax), Drp(m)%z(1:nCalMax), &
  &         Drp(m)%u(1:nCalMax), Drp(m)%v(1:nCalMax), Drp(m)%w(1:nCalMax), &
  &         Drp(m)%mr(1:nCalMax), Drp(m)%Cd(1:nCalMax), Drp(m)%AsRatio(1:nCalMax)  )
  Drp(m)%m(:) = 1
  Drp(m)%u(:) = 0.0; Drp(m)%v(:) = 0.0; Drp(m)%w(:) = 0.0
  Drp(m)%mr(:) = 1.0; Drp(m)%Cd(:) = 0.0; Drp(m)%AsRatio(:) = 0.0
 enddo
 ! 着氷計算領域 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_Resolution1D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(IceRslFile), strtxt, &
  &      Ice(m)%is, Ice(m)%ie )
! enddo
 ! 衝突セル ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  ! メモリ確保
  allocate( Flw(m)%nimp(Ice(m)%is: Ice(m)%ie, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%uimp(Ice(m)%is: Ice(m)%ie, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%vimp(Ice(m)%is: Ice(m)%ie, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%wimp(Ice(m)%is: Ice(m)%ie, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%simp(Ice(m)%is: Ice(m)%ie, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%CE  (Ice(m)%is: Ice(m)%ie, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%ImpMass(Ice(m)%is: Ice(m)%ie, Flw(m)%ks: Flw(m)%ke) )
  Flw(m)%nimp(:,:) = 0
  Flw(m)%uimp(:,:) = 0.0; Flw(m)%vimp(:,:) = 0.0; Flw(m)%wimp(:,:) = 0.0; Flw(m)%simp(:,:) = 0
  Flw(m)%ImpMass(:,:) =0.0
! enddo
 ! 液滴計算途中解 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 if( nDrop > nDrpOutputCount .and. CalNum /= 2 ) then
   write(fname, '(a,i8.8,a)') 'Drop', nDrop, '_'
!   do m = ms, me
 m = me
    call Input_Impingement3D( &
    &      trim(DrpImpDir) // trim(BlkName(m)) // trim(ND_DrpImpiDataFile), strbin, &
    &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, &
    &      Flw(m)%nimp, Flw(m)%uimp, Flw(m)%vimp, Flw(m)%wimp )
!    call Input_Impingement3D( &
!    &      trim(DrpImpDir) // trim(fname) // trim(BlkName(m)) // trim(ND_DrpImpiDataFile), strbin, &
!    &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, &
!    &      Flw(m)%nimp, Flw(m)%uimp, Flw(m)%vimp, Flw(m)%wimp )
    nStart = nDrop
!   enddo
  else
   nStart = 1
 endif
 ! 処理終了 ********************************************************************************************
 return
end subroutine InitialSetting
!*******************************************************************************************************
!******** 液滴軌道計算の初期条件								********
!*******************************************************************************************************
subroutine InitialDropletCondition
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: rxMin = -7.0
 real   , parameter :: ryMin = -0.9!-0.2 !-0.9 !4deg用　!-0.1 0deg用
 real   , parameter :: ryMax = +0.2!+0.2 !-0.3 !4deg用　!+0.1 0deg用
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m
 integer :: i
 integer :: ip, jp, kp
 integer :: imin, imax
 real    :: dis
 ! 処理開始 ********************************************************************************************
 ! 液滴投入ブロック ------------------------------------------------------------------------------------
 m = mini
 ip = Flw(m)%i2; jp = Flw(m)%js; kp = int(0.5 * Flw(m)%ke)
 ! 液滴投入範囲 ----------------------------------------------------------------------------------------
 ! x 方向
 xmin = Flw(m)%x(ip,jp,kp) + rxMin * chord
 ! y 方向
 ymin = minval( Flw(m)%y(Flw(m)%i1:Flw(m)%i3,jp,kp) ) + ryMin * chord
 ymax = maxval( Flw(m)%y(Flw(m)%i1:Flw(m)%i3,jp,kp) ) + ryMax * chord
 write(*,'(a,e16.8e3)')'* ymin = ',ymin
 write(*,'(a,e16.8e3)')'* ymax = ',ymax
 ! z 方向
! zmin = Flw(m)%z(ip, jp, Flw(m)%ks + 1)
! zmax = Flw(m)%z(ip, jp, Flw(m)%ke - 1)
!Mainグリッドはz方向に4分割,Subグリッドは5分割なので,Mainに合わせる
 zmin = Flw(m)%z(ip, jp, Flw(m)%ks) + (Flw(m)%z(ip, jp, Flw(m)%ke) - Flw(m)%z(ip, jp, Flw(m)%ks)) / 3.0
 zmax = Flw(m)%z(ip, jp, Flw(m)%ke) - (Flw(m)%z(ip, jp, Flw(m)%ke) - Flw(m)%z(ip, jp, Flw(m)%ks)) / 3.0
 ! 液滴投入速度 ----------------------------------------------------------------------------------------
 uini = VelExp * cos(AOA)
 vini = VelExp * sin(AOA)
 wini = 0.0
 DrpInVel = sqrt(uini**2 + vini**2 + wini**2)
 ! 液滴投入面積 ----------------------------------------------------------------------------------------
 DrpInArea = (ymax - ymin) * (zmax - zmin)
 ! ファイル出力 ----------------------------------------------------------------------------------------
 if(nDrop == 0) then
   call Output_DropInCondition( &
   &      trim(DrpImpDir) // trim(DrpInConditionFile), strtxt, &
   &      DrpInArea * lRef**2,  DrpInVel * aRef )
 endif
 ! 処理終了 ********************************************************************************************
 return
end subroutine InitialDropletCondition
!*******************************************************************************************************
!******** 液滴軌道計算										********
!*******************************************************************************************************
subroutine CalDropletTrajectory
 ! 変数宣言 ********************************************************************************************
 implicit none
 character :: fname * 50
 ! 処理開始 ********************************************************************************************
 ! ログファイル設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 open(9, file = trim(DrpImpDir) // 'Log.txt')
 write(9, '(a)')         '+++ Numerical Condition +++'
 write(9, '(a,i2)')      '* Exp. case        = ', ExpCaseNum
 write(9, '(a,e10.4e1)') '* Cn               = ', Cnd
 write(9, '(a,i9)')      '* Droplet Count    = ', nDrpIn
 write(9, '(a,i9)')      '* Max. Calculation = ', nCalMax
 write(*, '(a)')          '+++ Numerical Condition +++'
 write(*, '(a,i2)')       '* Exp. case        = ', ExpCaseNum
 write(*, '(a,e10.4e1)')  '* Cn               = ', Cnd
 write(*, '(a,i9)')       '* Droplet Count    = ', nDrpIn
 write(*, '(a,i9)')       '* Max. Calculation = ', nCalMax
 ! 初期値 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 DrpNum = 0
 Nlost  = 0
 if(CalNum == 2) then
   nDrpIn       = nDrpView
   nDrpFile     = nDrpView
   nOutputCount = nDrpView
 endif
 call genrand_init( nStart )
 ! 液滴投入ループ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do nDrop = nStart, nDrpIn
  ! 液滴ファイル番号
  DrpNum = mod(DrpNum, nDrpFile) + 1
  ! 液滴軌道計算
  call DropletTimeStep
  ! 計算途中解出力
  if( mod(nDrop, nDrpOutputCount) == 0.0 ) then
    select case(CalNum)
     case(1)
      call CalCollectionEfficiency
     case(2)
      call OutputViewFile
     case(3)
      call CalCollectionEfficiency
      call OutputViewFile
    end select
  endif
!!!  if(fSplash .or. fBounce) then
!!!    call ViewSLDPhenomena( DrpNum )
!!!  endif
  ! 計算ログ出力
  if( mod(nDrop, nOutputLog) == 0.0 .or. &
  &   nDrop == 1 .or. CalNum == 2 ) then
!  open(15,file='alp1.txt')
!  write(15,'(e16.8e3,2(5X,e16.8e3))') Ohw, Rew
    write(9, '(/, a)')    '+++ Computation Progress +++'
    write(9, '(a, i11)')  '* Droplet Number      = ', nDrop
!    write(9, '(a, i11)')  '* Impingement on Main = ', nint( sum(Flw(1)%nimp(:,:)) )
    write(9, '(a, i11)')  '* Impingement on Sub1 = ', nint( sum(Flw(2)%nimp(:,Flw(2)%ks+1)) )
    write(9, '(a, i11)')  '* Impingement on Sub2 = ', nint( sum(Flw(2)%nimp(:,Flw(2)%ks+2)) )
    write(9, '(a, i11)')  '* Impingement on Sub3 = ', nint( sum(Flw(2)%nimp(:,Flw(2)%ke-1)) )
    write(9, '(a, i11)')  '* Missing Number      = ', Nlost
    write(9, '(a, i11)')  '* Splash Number       = ', nDrpSpl
    write(9, '(a, i11)')  '* Bounce Number       = ', nDrpBou
    write(9, '(a, i11)')  '* 2 impinge Number    = ', nDrp2Imp
    write(9, '(a, i11)')  '* Breakup Number      = ', nDrpBrk
    write(*, '(/, a)')    '+++ Computation Progress +++'
    write(*, '(a, i11)')  '* Droplet Number      = ', nDrop
!    write(*, '(a, i11)')  '* Impingement on Main = ', nint( sum(Flw(1)%nimp(:,:)) )
    write(*, '(a, i11)')  '* Impingement on Sub1 = ', nint( sum(Flw(2)%nimp(:,Flw(2)%ks+1)) )
    write(*, '(a, i11)')  '* Impingement on Sub2 = ', nint( sum(Flw(2)%nimp(:,Flw(2)%ks+2)) )
    write(*, '(a, i11)')  '* Impingement on Sub3 = ', nint( sum(Flw(2)%nimp(:,Flw(2)%ke-1)) )
    write(*, '(a, i11)')  '* Missing Number      = ', Nlost
    write(*, '(a, i11)')  '* Splash Number       = ', nDrpSpl
    write(*, '(a, i11)')  '* Bounce Number       = ', nDrpBou
    write(*, '(a, i11)')  '* 2 impinge Number    = ', nDrp2Imp
    write(*, '(a, i11)')  '* Breakup Number      = ', nDrpBrk
  endif
 enddo
 call OutputViewFile
 close(9)
! close(15)
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalDropletTrajectory
!*******************************************************************************************************
!******** 液滴軌道の時間進行									********
!*******************************************************************************************************
subroutine DropletTimeStep
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, parameter :: nRetry = 100					! 液滴探索最大再試行回数
 integer, parameter :: mRetry = 10					! 液滴補間最大再試行回数
 integer, parameter :: SR     = 2					! 液滴探索範囲
 real   , parameter :: MGN    = 1.0e-3					! 補間経緯数許容誤差
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: n, m
 integer :: ip, jp, kp, mp
 integer :: ip0, jp0, kp0
 real(8) :: rr
 real    :: xp0, yp0, zp0, up0, vp0, wp0
 real    :: dt
 real    :: mr
 logical :: fMove, fImpi, fExit
 character :: fname * 20
 ! 処理開始 ********************************************************************************************
 ! 液滴初期位置及び速度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 mp   = mini
 xp   = xmin
 up   = uini; vp = vini; wp = wini
 dp   = MVD
 mr   = 1.0
 fSplash = .false.
 fBounce = .false.
 fImpi1 = .false.
 fImpi2 = .false.
 fBreakup = .false.
 do n = 1, nRetry
  select case(CalNum)
   case(2)
    ! y 方向
    yp = ymin + (ymax - ymin) * real(nDrop - 1) / real(nDrpView - 1)
    ypini = yp
    ! z 方向
    zp = 0.5 *  (zmax + zmin)
   case default
    ! y 方向 (乱数)
    call genrand_real1(rr)
    yp = ymin + (ymax - ymin) * real(rr)
    ypini = yp
    ! z 方向 (乱数)
    call genrand_real1(rr)
    zp = zmin + (zmax - zmin) * real(rr)
  end select
  ! 初期位置探索
  do m = 1, mRetry
   call SearchDroplet( &
   &      mp, &
   &      Flw(mp)%is, Flw(mp)%ie - 1, &
   &      Flw(mp)%js, Flw(mp)%je - 1, &
   &      Flw(mp)%ks, Flw(mp)%ke - 1, &
   &      MGN * real(m - 1), ip, jp, kp )
   if(fSearch) exit
  enddo
  if(fSearch) exit
 enddo
 if(.not. fSearch) then
   write(9, '(a)') '!!!!! Error : Droplet initial position !!!!!'
   write(9, '(a, i11)') '* Droplet number   = ', nDrop
   write(9, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
   write(*, '(a)') '!!!!! Error : Droplet initial position !!!!!'
   write(*, '(a, i11)') '* Droplet number   = ', nDrop
   write(*, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
   return
 endif
 ! 液滴のアスペクト比の出力ファイル ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(fname, '(a, 2(i2.2,a))') 'IceStep', IceStep, 'of', IceStepMax, '_'
 open(12, file = trim(DrpImpDir) // trim(fname) // 'AspectRatio.txt')
 Time0 = 0.0
 ! 液滴の時間進行 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do nCount = 1, nCalMax
  ! 時間刻み -------------------------------------------------------------------------------------------
  call dt3D( uxp, vep, wzp, Cnd, dt )
  if(CalNum == 2) then
   write(12,'(e16.8e3,3x,e16.8e3)')Time0,AsRatio
  Time0 = Time0 + dt
  endif
  ! 液滴位置探索 ---------------------------------------------------------------------------------------
  ip0 = ip; jp0 = jp; kp0 = kp
  do n = 1, nRetry
   xp0 = xp; yp0 = yp; zp0 = zp
   up0 = up; vp0 = vp; wp0 = wp
   ! 時間進行
   call TimeEuler3D( &
   &      dt, fx, fy, fz, xp0, yp0, zp0, up0, vp0, wp0, &
   &      xp, yp, zp, up, vp, wp )
   ! 液滴探索
   do m = 1, mRetry
    call SearchDroplet( &
    &      mp, &
    &      max(Flw(mp)%is, ip0 - SR), min(Flw(mp)%ie - 1, ip0 + SR), &
    &      max(Flw(mp)%js, jp0 - SR), min(Flw(mp)%je - 1, jp0 + SR), &
    &      max(Flw(mp)%ks, kp0 - SR), min(Flw(mp)%ke - 1, kp0 + SR), &
    &      MGN * real(m - 1), ip, jp, kp )
    if(fSearch) exit
   enddo
   if(fSearch) then
     exit
    else
     call BoundaryBladeSurface( ip0, jp0, kp0, Flw(mp)%js, mp, fImpi, up, vp, wp, mr )
     if(fImpi) return
   endif
  enddo
  if(.not. fSearch) then
    write(9, '(a)') '!!!!! Error : Droplet search !!!!!'
    write(9, '(a, i6)')  '* Iteration number = ', nCount
    write(9, '(a, i11)') '* Droplet number   = ', nDrop
    write(9, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip0, ',', jp0, ',', kp0, ')'
    write(*, '(a)') '!!!!! Error : Droplet search !!!!!'
    write(*, '(a, i6)')  '* Iteration number = ', nCount
    write(*, '(a, i11)') '* Droplet number   = ', nDrop
    write(*, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip0, ',', jp0, ',', kp0, ')'
    Nlost = Nlost + 1
    return								! ロスト．次の液滴へ
  endif
  ! 境界条件 -------------------------------------------------------------------------------------------
  fMove = .false.; fImpi = .false.; fExit = .false.
  select case(mp)
   ! Main Grid
   case(1)
!    ! 流出境界
!    if(.not. fMove) call BoundaryOutlet( ip, jp, mp, fExit )
!    if(fExit) return
!    ! Sub Grid 境界
!    if(.not. fMove) call BoundaryMGtoSG( ip, jp, kp, mp, fMove )
!    ! C 型格子ブランチ・カット
!    if(.not. fMove) call BoundaryCtypeBranchCut( ip, jp, mp, fMove )
!    ! 翼表面
!    if(.not. fMove) call BoundaryBladeSurface( ip, jp, kp, Flw(mp)%js, mp, fImpi, up, vp, wp, mr )
!    if(fImpi) return
!    ! 周期境界
!    if(.not. fMove) call BoundaryPeriodic( ip, jp, kp, mp, fMove )
    write(*,*) 'error:boundary condition in droplet trajectory'
    stop
   ! Sub Grid
   case(2)
!    ! Main Grid 境界
!    if(.not. fMove) call BoundarySGtoMG( ip, jp, kp, mp, fMove )
    ! 流出境界
    if(.not. fMove) call BoundaryOutlet( ip, jp, mp, fExit )
    if(fExit) return
    ! C 型格子ブランチ・カット
    if(.not. fMove) call BoundaryCtypeBranchCut( ip, jp, mp, fMove )
    ! 翼表面
    if(.not. fMove) call BoundaryBladeSurface( ip, jp, kp, Flw(mp)%js, mp, fImpi, up, vp, wp, mr )
    if(fImpi) return
    ! 周期境界
    if(.not. fMove) call BoundaryPeriodic( ip, jp, kp, mp, fMove )
   case default
    write(*, '(a)') '!!!!! Error : Block index !!!!!'
    stop
  end select
  ! 格子移動後の液滴位置探索 ---------------------------------------------------------------------------
  if(fMove) then
    do n = 1, nRetry
     ! 液滴位置探索
     do m = 1, mRetry
      call SearchDroplet( &
      &      mp, &
      &      Flw(mp)%is, Flw(mp)%ie - 1, &
      &      Flw(mp)%js, Flw(mp)%je - 1, &
      &      Flw(mp)%ks, Flw(mp)%ke - 1, &
      &      MGN * real(m - 1), ip, jp, kp )
      if(fSearch) exit
     enddo
     if(fSearch) exit
     ! 時間進行 (探索失敗時)
     xp0 = xp; yp0 = yp; zp0 = zp
     up0 = up; vp0 = vp; wp0 = wp
     call TimeEuler3D( &
     &      0.1 * dt, fx, fy, fz, xp0, yp0, zp0, up0, vp0, wp0, &
     &      xp, yp, zp, up, vp, wp )
     ! 格子を突き抜けた場合の処理
     call BoundaryBladeSurface( ip, jp, kp, Flw(mp)%js, mp, fImpi, up, vp, wp, mr )
     if(fImpi) return
    enddo
    if(.not. fSearch) then
      write(9, '(a)') '!!!!! Error : Droplet search !!!!!'
      write(9, '(a, i6)')  '* Iteration number = ', nCount
      write(9, '(a, i11)') '* Droplet number   = ', nDrop
      write(9, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip0, ',', jp0, ',', kp0, ')'
      write(*, '(a)') '!!!!! Error : Droplet search !!!!!'
      write(*, '(a, i6)')  '* Iteration number = ', nCount
      write(*, '(a, i11)') '* Droplet number   = ', nDrop
      write(*, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip0, ',', jp0, ',', kp0, ')'
      Nlost = Nlost + 1
      return								! ロスト．次の液滴へ
    endif
  endif
  ! 液滴軌道データログ ---------------------------------------------------------------------------------
  Drp(DrpNum)%step       = nCount
  Drp(DrpNum)%m (nCount) = mp
  Drp(DrpNum)%x (nCount) = xp
  Drp(DrpNum)%y (nCount) = yp
  Drp(DrpNum)%z (nCount) = zp
  Drp(DrpNum)%u (nCount) = up
  Drp(DrpNum)%v (nCount) = vp
  Drp(DrpNum)%w (nCount) = wp
  Drp(DrpNum)%mr(nCount) = mr
  Drp(DrpNum)%Cd(nCount) = Cd
  Drp(DrpNum)%AsRatio(nCount) = AsRatio
 enddo
 close(12)
 ! 処理終了 ********************************************************************************************
 return
end subroutine DropletTimeStep
!*******************************************************************************************************
!******** 液滴探索										********
!*******************************************************************************************************
subroutine SearchDroplet( &
&            mp, isp, iep, jsp, jep, ksp, kep, MGN, ip, jp, kp )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: mp						! 補間点探索ブロック
 integer, intent(in)  :: isp, iep, jsp, jep, ksp, kep			! 補間点探索範囲
 real   , intent(in)  :: MGN						! 補間係数許容誤差
 integer, intent(out) :: ip, jp, kp					! 補間点格子番号
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i1, i2, i3, i4, i5, i6, i7, i8, &
 &          j1, j2, j3, j4, j5, j6, j7, j8, &
 &          k1, k2, k3, k4, k5, k6, k7, k8				! 補間点
 real    :: alp, bet, gam						! 補間係数
 integer :: fInter							! 補間パターン
 real    :: term1, term2, term3, term4, term5, term6, term7, term8	! 補間係数
 real    :: rhof, uf, vf, wf, muf                                       ! 周囲流体の物理量
 real    :: fdx, fdy, fdz						! 液滴に働く抗力
 real    :: fgx, fgy, fgz						! 液滴に働く重力
 real    :: gx, gy, gz							! 重力係数
 real    :: xixp, xiyp, xizp, etxp, etyp, etzp, zexp, zeyp, zezp
 real    :: dp1, dp2
 ! 処理開始 ********************************************************************************************
 ! 探索開始 (物理空間) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 fInter = 0
 do kp = ksp, kep
 do jp = jsp, jep
 do ip = isp, iep
  ! ワイド・サーチ -------------------------------------------------------------------------------------
  i1 = ip    ; j1 = jp    ; k1 = kp    ; i2 = ip + 1; j2 = jp    ; k2 = kp
  i3 = ip    ; j3 = jp + 1; k3 = kp    ; i4 = ip + 1; j4 = jp + 1; k4 = kp
  i5 = ip    ; j5 = jp    ; k5 = kp + 1; i6 = ip + 1; j6 = jp    ; k6 = kp + 1
  i7 = ip    ; j7 = jp + 1; k7 = kp + 1; i8 = ip + 1; j8 = jp + 1; k8 = kp + 1
  call WideSearch3D8Point( &
  &      xp, yp, zp, &
  &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
  &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
  &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
  &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
  &      Flw(mp)%x(i5,j5,k5), Flw(mp)%y(i5,j5,k5), Flw(mp)%z(i5,j5,k5), &
  &      Flw(mp)%x(i6,j6,k6), Flw(mp)%y(i6,j6,k6), Flw(mp)%z(i6,j6,k6), &
  &      Flw(mp)%x(i7,j7,k7), Flw(mp)%y(i7,j7,k7), Flw(mp)%z(i7,j7,k7), &
  &      Flw(mp)%x(i8,j8,k8), Flw(mp)%y(i8,j8,k8), Flw(mp)%z(i8,j8,k8), &
  &      MGN, fSearch )
  if(.not. fSearch) then
   cycle
  endif
  ! 三重線形補間 ---------------------------------------------------------------------------------------
  alp = 0.5; bet = 0.5; gam = 0.5
  call Interpolation3D8Point( &
  &      xp, yp, zp, &
  &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
  &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
  &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
  &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
  &      Flw(mp)%x(i5,j5,k5), Flw(mp)%y(i5,j5,k5), Flw(mp)%z(i5,j5,k5), &
  &      Flw(mp)%x(i6,j6,k6), Flw(mp)%y(i6,j6,k6), Flw(mp)%z(i6,j6,k6), &
  &      Flw(mp)%x(i7,j7,k7), Flw(mp)%y(i7,j7,k7), Flw(mp)%z(i7,j7,k7), &
  &      Flw(mp)%x(i8,j8,k8), Flw(mp)%y(i8,j8,k8), Flw(mp)%z(i8,j8,k8), &
  &  	 MGN, alp, bet, gam, fSearch )
  if(fSearch) then
    fInter = 1
    exit
  endif
  ! 三次元線形補間 -------------------------------------------------------------------------------------
  ! Pattern-1
  i1 = ip    ; j1 = jp    ; k1 = kp    ; i2 = ip    ; j2 = jp    ; k2 = kp + 1
  i3 = ip + 1; j3 = jp    ; k3 = kp + 1; i4 = ip    ; j4 = jp + 1; k4 = kp + 1
  alp = 0.5; bet = 0.5; gam = 0.5
  call Interpolation3D4Point( &
  &      xp, yp, zp, &
  &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
  &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
  &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
  &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
  &	 MGN, alp, bet, gam, fSearch )
  if(fSearch) then
    fInter = 2
    exit
  endif
  ! Pattern-2
  i1 = ip    ; j1 = jp    ; k1 = kp    ; i2 = ip    ; j2 = jp + 1; k2 = kp
  i3 = ip + 1; j3 = jp + 1; k3 = kp    ; i4 = ip    ; j4 = jp + 1; k4 = kp + 1
  alp = 0.5; bet = 0.5; gam = 0.5
  call Interpolation3D4Point( &
  &      xp, yp, zp, &
  &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
  &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
  &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
  &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
  &	 MGN, alp, bet, gam, fSearch )
  if(fSearch) then
    fInter = 3
    exit
  endif
  ! Pattern-3
  i1 = ip    ; j1 = jp    ; k1 = kp    ; i2 = ip + 1; j2 = jp + 1; k2 = kp
  i3 = ip + 1; j3 = jp    ; k3 = kp + 1; i4 = ip    ; j4 = jp + 1; k4 = kp + 1
  alp = 0.5; bet = 0.5; gam = 0.5
  call Interpolation3D4Point( &
  &      xp, yp, zp, &
  &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
  &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
  &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
  &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
  &	 MGN, alp, bet, gam, fSearch )
  if(fSearch) then
    fInter = 4
    exit
  endif
  ! Pattern-4
  i1 = ip    ; j1 = jp    ; k1 = kp    ; i2 = ip + 1; j2 = jp    ; k2 = kp
  i3 = ip + 1; j3 = jp + 1; k3 = kp    ; i4 = ip + 1; j4 = jp    ; k4 = kp + 1
  alp = 0.5; bet = 0.5; gam = 0.5
  call Interpolation3D4Point( &
  &      xp, yp, zp, &
  &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
  &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
  &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
  &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
  &	 MGN, alp, bet, gam, fSearch )
  if(fSearch) then
    fInter = 5
    exit
  endif
  ! Pattern-5
  i1 = ip + 1; j1 = jp + 1; k1 = kp    ; i2 = ip + 1; j2 = jp    ; k2 = kp + 1
  i3 = ip    ; j3 = jp + 1; k3 = kp + 1; i4 = ip + 1; j4 = jp + 1; k4 = kp + 1
  alp = 0.5; bet = 0.5; gam = 0.5
  call Interpolation3D4Point( &
  &      xp, yp, zp, &
  &      Flw(mp)%x(i1,j1,k1), Flw(mp)%y(i1,j1,k1), Flw(mp)%z(i1,j1,k1), &
  &      Flw(mp)%x(i2,j2,k2), Flw(mp)%y(i2,j2,k2), Flw(mp)%z(i2,j2,k2), &
  &      Flw(mp)%x(i3,j3,k3), Flw(mp)%y(i3,j3,k3), Flw(mp)%z(i3,j3,k3), &
  &      Flw(mp)%x(i4,j4,k4), Flw(mp)%y(i4,j4,k4), Flw(mp)%z(i4,j4,k4), &
  &	 MGN, alp, bet, gam, fSearch )
  if(fSearch) then
    fInter = 6
    exit
  endif
 enddo
 if(fInter /= 0) exit
 enddo
 if(fInter /= 0) exit
 enddo
 ! 補間係数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 select case(fInter)
  case(0)
   return
  case(1)
   term1 = (1.0 - alp) * (1.0 - bet) * (1.0 - gam)
   term2 =        alp  * (1.0 - bet) * (1.0 - gam)
   term3 = (1.0 - alp) *        bet  * (1.0 - gam)
   term4 =        alp  *        bet  * (1.0 - gam)
   term5 = (1.0 - alp) * (1.0 - bet) *        gam
   term6 =        alp  * (1.0 - bet) *        gam
   term7 = (1.0 - alp) *        bet  *        gam
   term8 =        alp  *        bet  *        gam
  case(2:6)
   term1 = 1.0
   term2 = alp
   term3 = bet
   term4 = gam
  case default
   write(*, '(a)') '!!!!! Error : fInter number !!!!!'
   stop
 end select
 ! 周囲流体の情報を補間 (計算空間) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 select case(fInter)
  case(1)
   xixp = term1 * Flw(mp)%xix(i1,j1,k1) + term2 * Flw(mp)%xix(i2,j2,k2) &
   &    + term3 * Flw(mp)%xix(i3,j3,k3) + term4 * Flw(mp)%xix(i4,j4,k4) &
   &    + term5 * Flw(mp)%xix(i5,j5,k5) + term6 * Flw(mp)%xix(i6,j6,k6) &
   &    + term7 * Flw(mp)%xix(i7,j7,k7) + term8 * Flw(mp)%xix(i8,j8,k8)
   xiyp = term1 * Flw(mp)%xiy(i1,j1,k1) + term2 * Flw(mp)%xiy(i2,j2,k2) &
   &    + term3 * Flw(mp)%xiy(i3,j3,k3) + term4 * Flw(mp)%xiy(i4,j4,k4) &
   &    + term5 * Flw(mp)%xiy(i5,j5,k5) + term6 * Flw(mp)%xiy(i6,j6,k6) &
   &    + term7 * Flw(mp)%xiy(i7,j7,k7) + term8 * Flw(mp)%xiy(i8,j8,k8)
   xizp = term1 * Flw(mp)%xiz(i1,j1,k1) + term2 * Flw(mp)%xiz(i2,j2,k2) &
   &    + term3 * Flw(mp)%xiz(i3,j3,k3) + term4 * Flw(mp)%xiz(i4,j4,k4) &
   &    + term5 * Flw(mp)%xiz(i5,j5,k5) + term6 * Flw(mp)%xiz(i6,j6,k6) &
   &    + term7 * Flw(mp)%xiz(i7,j7,k7) + term8 * Flw(mp)%xiz(i8,j8,k8)
   etxp = term1 * Flw(mp)%etx(i1,j1,k1) + term2 * Flw(mp)%etx(i2,j2,k2) &
   &    + term3 * Flw(mp)%etx(i3,j3,k3) + term4 * Flw(mp)%etx(i4,j4,k4) &
   &    + term5 * Flw(mp)%etx(i5,j5,k5) + term6 * Flw(mp)%etx(i6,j6,k6) &
   &    + term7 * Flw(mp)%etx(i7,j7,k7) + term8 * Flw(mp)%etx(i8,j8,k8)
   etyp = term1 * Flw(mp)%ety(i1,j1,k1) + term2 * Flw(mp)%ety(i2,j2,k2) &
   &    + term3 * Flw(mp)%ety(i3,j3,k3) + term4 * Flw(mp)%ety(i4,j4,k4) &
   &    + term5 * Flw(mp)%ety(i5,j5,k5) + term6 * Flw(mp)%ety(i6,j6,k6) &
   &    + term7 * Flw(mp)%ety(i7,j7,k7) + term8 * Flw(mp)%ety(i8,j8,k8)
   etzp = term1 * Flw(mp)%etz(i1,j1,k1) + term2 * Flw(mp)%etz(i2,j2,k2) &
   &    + term3 * Flw(mp)%etz(i3,j3,k3) + term4 * Flw(mp)%etz(i4,j4,k4) &
   &    + term5 * Flw(mp)%etz(i5,j5,k5) + term6 * Flw(mp)%etz(i6,j6,k6) &
   &    + term7 * Flw(mp)%etz(i7,j7,k7) + term8 * Flw(mp)%etz(i8,j8,k8)
   zexp = term1 * Flw(mp)%zex(i1,j1,k1) + term2 * Flw(mp)%zex(i2,j2,k2) &
   &    + term3 * Flw(mp)%zex(i3,j3,k3) + term4 * Flw(mp)%zex(i4,j4,k4) &
   &    + term5 * Flw(mp)%zex(i5,j5,k5) + term6 * Flw(mp)%zex(i6,j6,k6) &
   &    + term7 * Flw(mp)%zex(i7,j7,k7) + term8 * Flw(mp)%zex(i8,j8,k8)
   zeyp = term1 * Flw(mp)%zey(i1,j1,k1) + term2 * Flw(mp)%zey(i2,j2,k2) &
   &    + term3 * Flw(mp)%zey(i3,j3,k3) + term4 * Flw(mp)%zey(i4,j4,k4) &
   &    + term5 * Flw(mp)%zey(i5,j5,k5) + term6 * Flw(mp)%zey(i6,j6,k6) &
   &    + term7 * Flw(mp)%zey(i7,j7,k7) + term8 * Flw(mp)%zey(i8,j8,k8)
   zezp = term1 * Flw(mp)%zez(i1,j1,k1) + term2 * Flw(mp)%zez(i2,j2,k2) &
   &    + term3 * Flw(mp)%zez(i3,j3,k3) + term4 * Flw(mp)%zez(i4,j4,k4) &
   &    + term5 * Flw(mp)%zez(i5,j5,k5) + term6 * Flw(mp)%zez(i6,j6,k6) &
   &    + term7 * Flw(mp)%zez(i7,j7,k7) + term8 * Flw(mp)%zez(i8,j8,k8)
   rhof = term1 * Flw(mp)%rho(i1,j1,k1) + term2 * Flw(mp)%rho(i2,j2,k2) &
   &    + term3 * Flw(mp)%rho(i3,j3,k3) + term4 * Flw(mp)%rho(i4,j4,k4) &
   &    + term5 * Flw(mp)%rho(i5,j5,k5) + term6 * Flw(mp)%rho(i6,j6,k6) &
   &    + term7 * Flw(mp)%rho(i7,j7,k7) + term8 * Flw(mp)%rho(i8,j8,k8)
   uf   = term1 * Flw(mp)%u  (i1,j1,k1) + term2 * Flw(mp)%u  (i2,j2,k2) &
   &    + term3 * Flw(mp)%u  (i3,j3,k3) + term4 * Flw(mp)%u  (i4,j4,k4) &
   &    + term5 * Flw(mp)%u  (i5,j5,k5) + term6 * Flw(mp)%u  (i6,j6,k6) &
   &    + term7 * Flw(mp)%u  (i7,j7,k7) + term8 * Flw(mp)%u  (i8,j8,k8)
   vf   = term1 * Flw(mp)%v  (i1,j1,k1) + term2 * Flw(mp)%v  (i2,j2,k2) &
   &    + term3 * Flw(mp)%v  (i3,j3,k3) + term4 * Flw(mp)%v  (i4,j4,k4) &
   &    + term5 * Flw(mp)%v  (i5,j5,k5) + term6 * Flw(mp)%v  (i6,j6,k6) &
   &    + term7 * Flw(mp)%v  (i7,j7,k7) + term8 * Flw(mp)%v  (i8,j8,k8)
   wf   = term1 * Flw(mp)%w  (i1,j1,k1) + term2 * Flw(mp)%w  (i2,j2,k2) &
   &    + term3 * Flw(mp)%w  (i3,j3,k3) + term4 * Flw(mp)%w  (i4,j4,k4) &
   &    + term5 * Flw(mp)%w  (i5,j5,k5) + term6 * Flw(mp)%w  (i6,j6,k6) &
   &    + term7 * Flw(mp)%w  (i7,j7,k7) + term8 * Flw(mp)%w  (i8,j8,k8)
   muf  = term1 * Flw(mp)%mu (i1,j1,k1) + term2 * Flw(mp)%mu (i2,j2,k2) &
   &    + term3 * Flw(mp)%mu (i3,j3,k3) + term4 * Flw(mp)%mu (i4,j4,k4) &
   &    + term5 * Flw(mp)%mu (i5,j5,k5) + term6 * Flw(mp)%mu (i6,j6,k6) &
   &    + term7 * Flw(mp)%mu (i7,j7,k7) + term8 * Flw(mp)%mu (i8,j8,k8)
  case(2:6)
   xixp = term1 * Flw(mp)%xix(i1,j1,k1) + term2 * (Flw(mp)%xix(i2,j2,k2) - Flw(mp)%xix(i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%xix(i3,j3,k3) - Flw(mp)%xix(i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%xix(i4,j4,k4) - Flw(mp)%xix(i1,j1,k1))
   xiyp = term1 * Flw(mp)%xiy(i1,j1,k1) + term2 * (Flw(mp)%xiy(i2,j2,k2) - Flw(mp)%xiy(i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%xiy(i3,j3,k3) - Flw(mp)%xiy(i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%xiy(i4,j4,k4) - Flw(mp)%xiy(i1,j1,k1))
   xizp = term1 * Flw(mp)%xiz(i1,j1,k1) + term2 * (Flw(mp)%xiz(i2,j2,k2) - Flw(mp)%xiz(i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%xiz(i3,j3,k3) - Flw(mp)%xiz(i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%xiz(i4,j4,k4) - Flw(mp)%xiz(i1,j1,k1))
   etxp = term1 * Flw(mp)%etx(i1,j1,k1) + term2 * (Flw(mp)%etx(i2,j2,k2) - Flw(mp)%etx(i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%etx(i3,j3,k3) - Flw(mp)%etx(i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%etx(i4,j4,k4) - Flw(mp)%etx(i1,j1,k1))
   etyp = term1 * Flw(mp)%ety(i1,j1,k1) + term2 * (Flw(mp)%ety(i2,j2,k2) - Flw(mp)%ety(i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%ety(i3,j3,k3) - Flw(mp)%ety(i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%ety(i4,j4,k4) - Flw(mp)%ety(i1,j1,k1))
   etzp = term1 * Flw(mp)%etz(i1,j1,k1) + term2 * (Flw(mp)%etz(i2,j2,k2) - Flw(mp)%etz(i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%etz(i3,j3,k3) - Flw(mp)%etz(i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%etz(i4,j4,k4) - Flw(mp)%etz(i1,j1,k1))
   zexp = term1 * Flw(mp)%zex(i1,j1,k1) + term2 * (Flw(mp)%zex(i2,j2,k2) - Flw(mp)%zex(i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%zex(i3,j3,k3) - Flw(mp)%zex(i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%zex(i4,j4,k4) - Flw(mp)%zex(i1,j1,k1))
   zeyp = term1 * Flw(mp)%zey(i1,j1,k1) + term2 * (Flw(mp)%zey(i2,j2,k2) - Flw(mp)%zey(i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%zey(i3,j3,k3) - Flw(mp)%zey(i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%zey(i4,j4,k4) - Flw(mp)%zey(i1,j1,k1))
   zezp = term1 * Flw(mp)%zez(i1,j1,k1) + term2 * (Flw(mp)%zez(i2,j2,k2) - Flw(mp)%zez(i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%zez(i3,j3,k3) - Flw(mp)%zez(i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%zez(i4,j4,k4) - Flw(mp)%zez(i1,j1,k1))
   rhof = term1 * Flw(mp)%rho(i1,j1,k1) + term2 * (Flw(mp)%rho(i2,j2,k2) - Flw(mp)%rho(i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%rho(i3,j3,k3) - Flw(mp)%rho(i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%rho(i4,j4,k4) - Flw(mp)%rho(i1,j1,k1))
   uf   = term1 * Flw(mp)%u  (i1,j1,k1) + term2 * (Flw(mp)%u  (i2,j2,k2) - Flw(mp)%u  (i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%u  (i3,j3,k3) - Flw(mp)%u  (i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%u  (i4,j4,k4) - Flw(mp)%u  (i1,j1,k1))
   vf   = term1 * Flw(mp)%v  (i1,j1,k1) + term2 * (Flw(mp)%v  (i2,j2,k2) - Flw(mp)%v  (i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%v  (i3,j3,k3) - Flw(mp)%v  (i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%v  (i4,j4,k4) - Flw(mp)%v  (i1,j1,k1))
   wf   = term1 * Flw(mp)%w  (i1,j1,k1) + term2 * (Flw(mp)%w  (i2,j2,k2) - Flw(mp)%w  (i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%w  (i3,j3,k3) - Flw(mp)%w  (i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%w  (i4,j4,k4) - Flw(mp)%w  (i1,j1,k1))
   muf  = term1 * Flw(mp)%mu (i1,j1,k1) + term2 * (Flw(mp)%mu (i2,j2,k2) - Flw(mp)%mu (i1,j1,k1)) &
   &                                    + term3 * (Flw(mp)%mu (i3,j3,k3) - Flw(mp)%mu (i1,j1,k1)) &
   &                                    + term4 * (Flw(mp)%mu (i4,j4,k4) - Flw(mp)%mu (i1,j1,k1))
  case default
   write(*, '(a)') '!!!!! Error : fInter number !!!!!'
   stop
 end select
 ! 液滴速度 (計算空間) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Velocity3D( &
 &      xixp, xiyp, xizp, etxp, etyp, etzp, zexp, zeyp, zezp, up, vp, wp, &
 &      uxp, vep, wzp )
! ! 液滴のアスペクト比 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call CalAspectRatio( &
! &       dp, rhod, up, vp, wp, rhof, uf, vf, wf, mud, sigd, muf, &
! &       AsRatio)
 ! 液滴の分裂 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 dp1 = dp
 call DropletBreakup( &
  &      dp1, Rhod, up, vp, wp, Rhof, uf, vf, wf, mud, sigd, fBreakup, nDrpBrk, &
  &      dp2 )
 dp = dp2
 ! 液滴に働く抗力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 select case(DragNum)
 case(1)
  call Drag3D( &
  &      dp, Rhod, up, vp, wp, Rhof, uf, vf, wf, muf, &
  &      fdx, fdy, fdz, Cd )
 case(2)
  call DragGennaro( &
  &      dp, Rhod, up, vp, wp, Rhof, uf, vf, wf, muf, &
  &      fdx, fdy, fdz, Cd )
 case(3)
  call WiegandDrag( &
  &      dp, Rhod, up, vp, wp, Rhof, uf, vf, wf, muf, sigd, &
  &      fdx, fdy, fdz, Cd )
 case(4)
  call WebernumberDrag( &
  &      dp, Rhod, up, vp, wp, Rhof, uf, vf, wf, muf, sigd, &
  &      fdx, fdy, fdz, Cd )
end select
 ! 液滴に働く重力と浮力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 gx = +0.0
 gy = -0.0
 gz = 0.0
 gx = gx / (aRef**2 / lRef)
 gy = gy / (aRef**2 / lRef)
 gz = gz / (aRef**2 / lRef)
 call GravityBuoyancy3D( &
 &      Rhod, Rhof, gx, gy, gz, &
 &      fgx, fgy, fgz )
 ! 合力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 fx = fdx + fgx; fy = fdy + fgy; fz = fdz + fgz
 ! 処理終了 ********************************************************************************************
 return
end subroutine SearchDroplet
!*******************************************************************************************************
!******** Main Grid から Sub Grid への移動							********
!*******************************************************************************************************
subroutine BoundaryMGtoSG( ip, jp, kp, mp, fMove )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: ip, jp, kp
 integer, intent(inout) :: mp
 logical, intent(out)   :: fMove
 ! 処理開始 ********************************************************************************************
 ! 例外処理 --------------------------------------------------------------------------------------------
 if( ip == Flw(mp)%ie .or. jp == Flw(mp)%je .or. kp == Flw(mp)%ke ) then
  fMove = .false.
  return
 endif
 ! ブロック移動判定 ------------------------------------------------------------------------------------
 if( fOver(ip  , jp  , kp) .and. fOver(ip+1, jp  , kp) .and. &
 &   fOver(ip  , jp+1, kp) .and. fOver(ip+1, jp+1, kp) ) then
   fMove = .true.
   mp    = 2
  else
   fMove = .false.
   mp    = 1
 endif
 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryMGtoSG
!*******************************************************************************************************
!******** Sub Grid から Main Grid への移動							********
!*******************************************************************************************************
subroutine BoundarySGtoMG( ip, jp, kp, mp, fMove )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: ip, jp, kp
 integer, intent(inout) :: mp
 logical, intent(out)   :: fMove
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 logical :: isOut, ieOut, jeOut
 ! 処理開始 ********************************************************************************************
 isOut = ip <  Flw(mp)%is + 1 .and. uxp < 0.0
 ieOut = ip >= Flw(mp)%ie - 1 .and. uxp > 0.0
 jeOut = jp >= Flw(mp)%je - 1 .and. vep > 0.0
 if( isOut .or. ieOut .or. jeOut ) then
   fMove = .true.
   mp    = 1
  else
   fMove = .false.
   mp    = 2
 endif
 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundarySGtoMG
!*******************************************************************************************************
!******** 翼表面の衝突判定									********
!*******************************************************************************************************
subroutine BoundaryBladeSurface( ip, jp, kp, j0, mp, fImpi, up, vp, wp, mr )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: ip, jp, kp
 integer, intent(in)    :: j0
 integer, intent(in)    :: mp
 logical, intent(out)   :: fImpi
 real   , intent(inout) :: up, vp, wp, mr
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, parameter :: jMin   = 20
 integer, parameter :: nRetry = 100
 real   , parameter :: MGN    = 0.01
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real    :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
 real    :: alp, bet, gam
 integer :: n
 integer :: i1, i2, i3, i4, j1, j2, j3, j4, k1, k2, k3, k4
 logical :: fDivi
 real    :: wgh1, wgh2, wgh3, wgh4, wgh5, wgh6, wgh7, wgh8
 real    :: r1, r2, r3, r4
 real    :: u1, v1, w1, d1, m1, n1, u2, v2, w2, d2, m2, n2
 real    :: SLDLim
 integer,parameter	:: swi_splash_bound = 0 !0:全量付着,1:スプラッシュ・バウンドモデル(wright)
 ! 処理開始 ********************************************************************************************
 ! 例外処理 --------------------------------------------------------------------------------------------
 if( ip == Flw(mp)%ie .or. jp == Flw(mp)%je .or. kp == Flw(mp)%ke ) then
   fImpi = .false.
   return
 endif
 if( (mp == 1) .and. (ip < Flw(mp)%i1 .or. ip >= Flw(mp)%i3) ) then
   fImpi = .false.
   return
 endif
 if( jp > jMin ) then
   fImpi = .false.
   return
 endif
 ! 衝突判定 --------------------------------------------------------------------------------------------
 ! 周囲点のインデックス
 i1 = ip; k1 = kp; i2 = ip; k2 = kp + 1; i3 = ip + 1; k3 = kp
 ! 周囲 3点の座標
 x1 = Flw(mp)%x(i1,j0,k1); y1 = Flw(mp)%y(i1,j0,k1); z1 = Flw(mp)%z(i1,j0,k1)
 x2 = Flw(mp)%x(i2,j0,k2); y2 = Flw(mp)%y(i2,j0,k2); z2 = Flw(mp)%z(i2,j0,k2)
 x3 = Flw(mp)%x(i3,j0,k3); y3 = Flw(mp)%y(i3,j0,k3); z3 = Flw(mp)%z(i3,j0,k3)
 ! 衝突前の物理量
 u1 = up; v1 = vp; w1 = wp; d1 = dp; m1 = mr
 n1 = nDrpSpl + nDrpBou
 SLDLim = SplLim / lRef
 ! 衝突判定
 select case(swi_splash_bound)
  case(0)
   call ImpingementJudge3D( &
   &      xp, yp, zp, x1, y1, z1, x2, y2, z2, x3, y3, z3, MVD, &
   &      fImpi, Rhod, ImpMass )
   u2 = 0.0; v2 = 0.0; w2 = 0.0
   d2 = 0.0; m2 = 0.0; m2 = 0.0
  case(1)
   call SplashLEWICE3D( &
   &      xp, yp, zp, u1, v1, w1, d1, m1, x1, y1, z1, x2, y2, z2, x3, y3, z3, &
   &      LWC, Rhod, Sigd, mud, SLDLim, &
   &      fImpi, u2, v2, w2, d2, m2, nDrpSpl, nDrpBou, nDrp2Imp, ImpMass, Ohw, Rew, fSplash, fBounce, fImpi2 )
  ! call SplashSamenfink( &
  ! &      xp, yp, zp, u1, v1, w1, d1, m1, x1, y1, z1, x2, y2, z2, x3, y3, z3, &
  ! &      LWC, Rhod, Sigd, mud, SLDLim, Chord, &
  ! &      fImpi, u2, v2, w2, d2, m2, nDrpSpl, nDrpBou )
   ! 衝突後の物理量
   up = u2; vp = v2; wp = w2; dp = d2; mr = m2
   n2 = nDrpSpl + nDrpBou
   ! スプラッシュ・バウンドのログ
   if(n1 < n2 ) then
     if(d1 > d2 )then
  !     write(*, '(a)') '----- Splash -----'
  !     write(*, '(a, i11)') '* Droplet number   = ', nDrop
  !     write(*, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
  !     write(9, '(a)') '----- Splash -----'
  !     write(9, '(a, i11)') '* Droplet number   = ', nDrop
  !     write(9, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
  !     fSplash = .true.
      else
  !     write(*, '(a)') '----- Bounce -----'
  !     write(*, '(a, i11)') '* Droplet number   = ', nDrop
  !     write(*, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
  !     write(9, '(a)') '----- Bounce -----'
  !     write(9, '(a, i11)') '* Droplet number   = ', nDrop
  !     write(9, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
  !     fBounce = .true.
     endif
    else
     if(fImpi) then
  !     write(*, '(a)') '----- Impingement -----'
  !     write(*, '(a, i11)') '* Droplet number   = ', nDrop
  !     write(*, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
  !     write(9, '(a)') '----- Impingement -----'
  !     write(9, '(a, i11)') '* Droplet number   = ', nDrop
  !     write(9, '(a, 4(i4, a))') '* Cell index (', mp, ',', ip, ',', jp, ',', kp, ')'
     endif
   endif
 end select
 if(.not. fImpi) return
 ! 衝突量の分配 ----------------------------------------------------------------------------------------
 ! 周囲 4点と液滴の距離
 i1 = ip; k1 = kp; i2 = ip; k2 = kp + 1; i3 = ip + 1; k3 = kp; i4 = ip + 1; k4 = kp + 1
 x1 = Flw(mp)%x(i1,j0,k1); y1 = Flw(mp)%y(i1,j0,k1); z1 = Flw(mp)%z(i1,j0,k1)
 x2 = Flw(mp)%x(i2,j0,k2); y2 = Flw(mp)%y(i2,j0,k2); z2 = Flw(mp)%z(i2,j0,k2)
 x3 = Flw(mp)%x(i3,j0,k3); y3 = Flw(mp)%y(i3,j0,k3); z3 = Flw(mp)%z(i3,j0,k3)
 x4 = Flw(mp)%x(i4,j0,k4); y4 = Flw(mp)%y(i4,j0,k4); z4 = Flw(mp)%z(i4,j0,k4)
 r1 = sqrt( (x1 - xp)**2 + (y1 - yp)**2 + (z1 - zp)**2 )
 r2 = sqrt( (x2 - xp)**2 + (y2 - yp)**2 + (z2 - zp)**2 )
 r3 = sqrt( (x3 - xp)**2 + (y3 - yp)**2 + (z3 - zp)**2 )
 r4 = sqrt( (x4 - xp)**2 + (y4 - yp)**2 + (z4 - zp)**2 )
 ! 重み関数
 if( r1 == 0.0 ) then
   wgh1 = 1.0; wgh2 = 0.0; wgh3 = 0.0; wgh4 = 0.0
  else if( r2 == 0.0 ) then 
   wgh1 = 0.0; wgh2 = 1.0; wgh3 = 0.0; wgh4 = 0.0
  else if( r3 == 0.0 ) then 
   wgh1 = 0.0; wgh2 = 0.0; wgh3 = 1.0; wgh4 = 0.0
  else if( r4 == 0.0 ) then 
   wgh1 = 0.0; wgh2 = 0.0; wgh3 = 0.0; wgh4 = 1.0
  else
   alp  = (r2 * r3 * r4 + r1 * r3 * r4 + r1 * r2 * r4 + r1 * r2 * r3) / (r1 * r2 * r3 * r4)
   wgh1 = 1.0 / (r1 * alp)
   wgh2 = 1.0 / (r2 * alp)
   wgh3 = 1.0 / (r3 * alp)
   wgh4 = 1.0 / (r4 * alp)
 endif
 ! 衝突個数
 Flw(mp)%nimp(i1,k1) = Flw(mp)%nimp(i1,k1) + (m1 - m2) * wgh1
 Flw(mp)%uimp(i1,k1) = Flw(mp)%uimp(i1,k1) + u1        * wgh1
 Flw(mp)%vimp(i1,k1) = Flw(mp)%vimp(i1,k1) + v1        * wgh1
 Flw(mp)%wimp(i1,k1) = Flw(mp)%wimp(i1,k1) + w1        * wgh1
 Flw(mp)%nimp(i2,k2) = Flw(mp)%nimp(i2,k2) + (m1 - m2) * wgh2
 Flw(mp)%uimp(i2,k2) = Flw(mp)%uimp(i2,k2) + u1        * wgh2
 Flw(mp)%vimp(i2,k2) = Flw(mp)%vimp(i2,k2) + v1        * wgh2
 Flw(mp)%wimp(i2,k2) = Flw(mp)%wimp(i2,k2) + w1        * wgh2
 Flw(mp)%nimp(i3,k3) = Flw(mp)%nimp(i3,k3) + (m1 - m2) * wgh3
 Flw(mp)%uimp(i3,k3) = Flw(mp)%uimp(i3,k3) + u1        * wgh3
 Flw(mp)%vimp(i3,k3) = Flw(mp)%vimp(i3,k3) + v1        * wgh3
 Flw(mp)%wimp(i3,k3) = Flw(mp)%wimp(i3,k3) + w1        * wgh3
 Flw(mp)%nimp(i4,k4) = Flw(mp)%nimp(i4,k4) + (m1 - m2) * wgh4
 Flw(mp)%uimp(i4,k4) = Flw(mp)%uimp(i4,k4) + u1        * wgh4
 Flw(mp)%vimp(i4,k4) = Flw(mp)%vimp(i4,k4) + v1        * wgh4
 Flw(mp)%wimp(i4,k4) = Flw(mp)%wimp(i4,k4) + w1        * wgh4
 ! 衝突量
 Flw(mp)%ImpMass(i1,k1) = Flw(mp)%ImpMass(i1,k1) + ImpMass * wgh1
 Flw(mp)%ImpMass(i2,k2) = Flw(mp)%ImpMass(i2,k2) + ImpMass * wgh2
 Flw(mp)%ImpMass(i3,k3) = Flw(mp)%ImpMass(i3,k3) + ImpMass * wgh3
 Flw(mp)%ImpMass(i4,k4) = Flw(mp)%ImpMass(i4,k4) + ImpMass * wgh4
 ! スプラッシュの場合は計算続行
 if(n1 < n2 ) fImpi = .false.
 ! 2回めの衝突の場合 ***********************************************************************************
 if(fImpi2)then
  write(*,'(a)')               '***** droplet impinges twice *****'
  write(*,'(a, i11)')          '* Droplet Number =', nDrop
  write(*,'(a, e16.8e3)')      '* Initial yp     =', ypini
  write(*,'(a, 2(i4, a))')     '* ( ip, jp )     = ( ', ip, ', ', jp, ' )'
  write(*,'(a, 2(e16.8e3, a))')'* ( xp, yp )     = ( ', xp, ', ', yp, ' )'
 endif
 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryBladeSurface
!*******************************************************************************************************
!******** C 型格子ブランチ・カット								********
!*******************************************************************************************************
subroutine BoundaryCtypeBranchCut( ip, jp, mp, fMove )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: ip, jp
 integer, intent(in)  :: mp
 logical, intent(out) :: fMove
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 logical :: jsOut
 ! 処理開始 ********************************************************************************************
 ! 例外処理 --------------------------------------------------------------------------------------------
 if( ip >= Flw(mp)%i1 .and. ip < Flw(mp)%i3 ) return
 ! ブロック移動判定 ------------------------------------------------------------------------------------
 jsOut = jp < Flw(mp)%js + 1 .and. vep < 0.0
 if(jsOut) then
   xp    = xp + up / abs(vep) * 2.0
   yp    = yp + vp / abs(vep) * 2.0
   zp    = zp + wp / abs(vep) * 2.0
   fMove = .true.
  else
   fMove = .false.
 endif
 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryCtypeBranchCut
!*******************************************************************************************************
!******** 周期境界										********
!*******************************************************************************************************
subroutine BoundaryPeriodic( ip, jp, kp, mp, fMove )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: ip, jp, kp
 integer, intent(in)  :: mp
 logical, intent(out) :: fMove
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real    :: dz
 logical :: ksOut, keOut
 ! 処理開始 ********************************************************************************************
 ksOut = kp <  Flw(mp)%ks + 1 .and. wzp < 0.0
 keOut = kp >= Flw(mp)%ke - 1 .and. wzp > 0.0
 dz    = maxval( Flw(mp)%z(ip,jp,:Flw(mp)%ke-1) )  - minval( Flw(mp)%z(ip,jp,:) )
 if(ksOut) then
   zp = zp + dz
   fMove = .true.
   return
 endif
 if(keOut) then
   zp = zp - dz
   fMove = .true.
   return
 endif
 fMove = .false.
 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryPeriodic
!*******************************************************************************************************
!******** 流出境界										********
!*******************************************************************************************************
subroutine BoundaryOutlet( ip, jp, mp, fExit )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: ip, jp
 integer, intent(in)  :: mp
 logical, intent(out) :: fExit
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: j0
 integer :: ibs, ibe
 logical :: TEOut, isOut, ieOut, jeOut
 ! 処理開始 ********************************************************************************************
 ! ξ方向 ----------------------------------------------------------------------------------------------
! j0 = Flw(1)%js
 j0 = Flw(2)%js
 ! T.E. 通過
! ibs = Flw(1)%i1 + 1; ibe = Flw(1)%i3 - 1
! TEOut = xp > maxval( Flw(1)%x(ibs:ibe, j0, :) )
 ibs = Flw(2)%i1 + 1; ibe = Flw(2)%i3 - 1
 TEOut = xp > maxval( Flw(2)%x(ibs:ibe, j0, :) )
 ! 計算領域
 isOut = ip <  Flw(mp)%is + 1 .and. uxp < 0.0
 ieOut = ip >= Flw(mp)%ie - 1 .and. uxp > 0.0
 ! 流出判定
 if( TEOut .or. isOut .or. ieOut ) then
   fExit = .true.
   return
  else
   fExit = .false.
 endif
 ! η方向 ----------------------------------------------------------------------------------------------
 jeOut = jp >= Flw(mp)%je - 1 .and. vep > 0.0
 ! 流出判定
 if( jeOut ) then
   fExit = .true.
   return
  else
   fExit = .false.
 endif
 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryOutlet
!*******************************************************************************************************
!******** 衝突面積計算										********
!*******************************************************************************************************
subroutine CalImpingementArea
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m
 integer :: j0
 ! 処理開始 ********************************************************************************************
! do m = ms, me
 m = me
  j0 = Flw(m)%js
  ! 衝突面積
  call ImpingementArea3D( &
  &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x(Ice(m)%is:Ice(m)%ie, j0, :), &
  &      Flw(m)%y(Ice(m)%is:Ice(m)%ie, j0, :), &
  &      Flw(m)%z(Ice(m)%is:Ice(m)%ie, j0, :), &
  &      Flw(m)%simp )
  ! ファイル出力
  if( nDrop == 0 ) then
    call Output_ArrayReal2D( &
    &      trim(DrpImpDir) // trim(BlkName(m)) // trim(DrpImpiAreaFile), strbin, &
    &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, Flw(m)%simp * lRef**2 )
  endif
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalImpingementArea
!*******************************************************************************************************
!******** 収集効率計算										********
!*******************************************************************************************************
subroutine CalCollectionEfficiency
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: ui(:, :), vi(:, :), wi(:, :)
 integer   :: m
 integer   :: j0
 character :: fname * 20
 ! 処理開始 ********************************************************************************************
 write(fname, '(a, 2(i2.2,a))') 'IceStep', IceStep, 'of', IceStepMax, '_'
 ! 計算ログ出力 ----------------------------------------------------------------------------------------
 call Output_CalSetting( trim(ND_CalSetFile) // strtxt )
 ! 衝突速度及び収集効率 --------------------------------------------------------------------------------
! do m = ms, me
 m = me
  ! 初期設定
  allocate( ui(Ice(m)%is:Ice(m)%ie, Flw(m)%ks:Flw(m)%ke), &
  &         vi(Ice(m)%is:Ice(m)%ie, Flw(m)%ks:Flw(m)%ke), &
  &         wi(Ice(m)%is:Ice(m)%ie, Flw(m)%ks:Flw(m)%ke) )
  ui(:, :) = Flw(m)%uimp(:, :)
  vi(:, :) = Flw(m)%vimp(:, :)
  wi(:, :) = Flw(m)%wimp(:, :)
  ! 衝突速度
  call ImpingementVelocity3D( &
  &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, Flw(m)%nimp, &
  &      ui, vi, wi )
  ! 収集効率
  call CollectionEfficiency3D( &
  &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, &
  &      nDrop, MVD, LWC, Rhod, DrpInArea, DrpInVel, &
  &      Flw(m)%nimp, Flw(m)%simp, Flw(m)%ImpMass, Flw(m)%CE )
  ! ファイル出力
  call Output_Impingement3D( &
  &      trim(DrpImpDir) // trim(BlkName(m)) // trim(ND_DrpImpiDataFile), strbin, &
  &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%nimp, Flw(m)%uimp, Flw(m)%vimp, Flw(m)%wimp )
  call Output_CollectionEfficiency3D( &
  &      trim(DrpimpDir) // trim(BlkName(m)) // trim(CollectionFile), strbin, &
  &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%CE * aRef / lRef**3, ui * aRef, vi * aRef, wi * aRef )
  call Output_Impingement3D( &
  &      trim(DrpImpDir) // trim(fname) // trim(BlkName(m)) // trim(ND_DrpImpiDataFile), strbin, &
  &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%nimp, Flw(m)%uimp, Flw(m)%vimp, Flw(m)%wimp )
  call Output_CollectionEfficiency3D( &
  &      trim(DrpimpDir) // trim(fname) // trim(BlkName(m)) // trim(CollectionFile), strbin, &
  &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%CE * aRef / lRef**3, ui * aRef, vi * aRef, wi * aRef )
  ! メモリ解法
  deallocate( ui, vi, wi )
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalCollectionEfficiency
!*******************************************************************************************************
!******** 液滴軌道可視化ファイル出力								********
!*******************************************************************************************************
subroutine OutputViewFile
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: n
 character :: fname * 20
 ! 処理開始 ********************************************************************************************
! do n = 1, nDrpFile
!  write(fname, '(a,i5.5,a)') 'Drop', n, '_'
!  call MakeMAVSFileDroplet3D( &
!  &      trim(DrpTraDir), trim(fname) // trim(DrpTrajectoryFile), strfld, &
!  &      lRef, aRef, Drp(n)%step, &
!  &      Drp(n)%x(1:Drp(n)%step), Drp(n)%y(1:Drp(n)%step), Drp(n)%z(1:Drp(n)%step), &
!  &      Drp(n)%u(1:Drp(n)%step), Drp(n)%v(1:Drp(n)%step), Drp(n)%w(1:Drp(n)%step) )
! enddo
 do n = 1, nDrpFile
  write(fname, '(a,i5.5,a)') 'Drop', n, '_'
   call MakeMAVSFileDropletWe3D( &
  &      trim(DrpTradir), trim(fname) // trim(DrpTrajectoryFile), strfld, &
  &      Drp(n)%step, lRef, Drp(n)%x(1:Drp(n)%step), Drp(n)%y(1:Drp(n)%step), Drp(n)%z(1:Drp(n)%step), &
  &      Drp(n)%cd(1:Drp(n)%step) )
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine OutputViewFile
!*******************************************************************************************************
!******** SLD現象の可視化ファイル								********
!*******************************************************************************************************
subroutine ViewSLDPhenomena(n)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer  , intent(in) :: n
 character :: fname * 20
 ! 処理開始 ********************************************************************************************
 if(fSplash .and. fBounce) then
!   write(fname, '(a,i5.5,a)') 'SLD', nDrpSpl + nDrpBou, '_'
!   call MakeMAVSFileSLD3D( &
!   &      trim(DrpTraDir), trim(fname) // trim(DrpTrajectoryFile), strfld, &
!   &      lRef, aRef, Drp(n)%step, &
!   &      Drp(n)%x(1:Drp(n)%step), Drp(n)%y(1:Drp(n)%step), Drp(n)%z(1:Drp(n)%step), &
!   &      Drp(n)%mr(1:Drp(n)%step) )
  else if(fSplash) then
   write(fname, '(a,i5.5,a)') 'Splash', nDrpSpl, '_'
   call MakeMAVSFileSLD3D( &
   &      trim(DrpTraDir), trim(fname) // trim(DrpTrajectoryFile), strfld, &
   &      lRef, aRef, Drp(n)%step, &
   &      Drp(n)%x(1:Drp(n)%step), Drp(n)%y(1:Drp(n)%step), Drp(n)%z(1:Drp(n)%step), &
   &      Drp(n)%mr(1:Drp(n)%step) )
  else
!   write(fname, '(a,i5.5,a)') 'Bounce', nDrpBou, '_'
!   call MakeMAVSFileDroplet3D( &
!   &      trim(DrpTraDir), trim(fname) // trim(DrpTrajectoryFile), strfld, &
!   &      lRef, aRef, Drp(n)%step, &
!   &      Drp(n)%x(1:Drp(n)%step), Drp(n)%y(1:Drp(n)%step), Drp(n)%z(1:Drp(n)%step), &
!   &      Drp(n)%u(1:Drp(n)%step), Drp(n)%v(1:Drp(n)%step), Drp(n)%w(1:Drp(n)%step) )
 endif
 ! 処理終了 ********************************************************************************************
 return
end subroutine ViewSLDPhenomena
!*******************************************************************************************************
!******** メモリ解放										********
!*******************************************************************************************************
subroutine Deallocating
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m
 ! 処理開始 ********************************************************************************************
! do m = ms, me
 m = me
  deallocate( Flw(m)%x, Flw(m)%y, Flw(m)%z, &
  &           Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%mu, &
  &           Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
  &           Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
  &           Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
  &           Flw(m)%nimp, Flw(m)%uimp, Flw(m)%vimp, Flw(m)%wimp, Flw(m)%CE )
  deallocate( Drp(m)%m, Drp(m)%x, Drp(m)%y, Drp(m)%z, Drp(m)%u, Drp(m)%v, Drp(m)%w, &
  &           Drp(m)%mr, Drp(m)%Cd, Drp(m)%AsRatio )
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine Deallocating
! 定義終了 *********************************************************************************************
end program DropletTrajectory_NACA
