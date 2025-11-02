!*******************************************************************************************************
!*******************************************************************************************************
!******** 計算結果検証プログラム								********
!******** (NACA翼，三次元圧縮性乱流場，重合格子法，高レイノルズ数型 k-eモデル)		  	********
!********					      2013.01.27  PROGRAMMED BY RYOSUKE HAYASHI ********
!*******************************************************************************************************
!*******************************************************************************************************
program Validation_NACA
 ! モジュール宣言 **************************************************************************************
 use Package_NACA
 use Package_FileIO
 use Package_Flow
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! ファイル名 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: YpFile         *  2 = 'Yp'
 character, parameter :: TsFile         *  2 = 'Ts'
 character, parameter :: T0File         *  2 = 'T0'
 character, parameter :: CEFile         * 20 = 'CollectionEfficiency'
 character, parameter :: IceThickFile   * 12 = 'IceThickness'
 character, parameter :: CleanBladeFile * 10 = 'CleanBlade'
 character, parameter :: IcingBladeFile * 10 = 'IcingBlade'
 ! 共有変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 type Type_Validation
  real   , pointer :: DisLE(:)						! L.E.からの距離
  real   , pointer :: yp(:)						! y+
  real   , pointer :: Ts(:)						! 静温
  real   , pointer :: CE(:)						! 収集効率
  real   , pointer :: Bi(:)						! 氷層厚さ
  real   , pointer :: Ti(:)						! 表面温度
  real   , pointer :: Bw(:)						! 水膜厚さ
  real   , pointer :: Tw(:)						! 水膜温度
  real   , pointer :: xi(:), yi(:)					! 着氷翼座標 (氷層)
  integer :: is, ie
  integer :: ks, ke
  integer :: i0, j0, k0							! 前縁の格子番号
 end type Type_Validation
 type(Type_Validation), pointer :: Vld(:)
 character :: IceStepName * 20
 ! 処理開始 ********************************************************************************************
 ! 検証実験ケース選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call SelectExpCase
 write(*, '(a)') "+++++ Select Exp. case complete. +++++"
 ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call InitialSetting
 write(*, '(a)') "+++++ Initial setting complete. +++++"
 ! 翼面に沿った前縁からの距離 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call DistanceLE
 write(*, '(a)') "+++++ Distance from L.E. complete. +++++"
 ! y+ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call CalYp
 write(*, '(a)') "+++++ y+ complete. +++++"
 ! 静温分布及び淀み点温度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call CalStaticTemperature
 write(*, '(a)') "+++++ Static temperature complete. +++++"
 ! 収集効率 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call CalCollectionEfficiency
 write(*, '(a)') "+++++ Collection efficiency complete. +++++"
 ! 氷層厚さ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call CalIceThickness
 write(*, '(a)') "+++++ Icing surface complete. +++++"
 ! 着氷翼の可視化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call ViewIcingBlade
 write(*, '(a)') "+++++ Visualized icing blade complete. +++++"
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
 ! ディレクトリ設定
 if(IceStep == 0 .or. IceStep == 1) then
   FlwCalInDir = bckdir // 'flow//cal//clean//'
   DrpImpDir   = bckdir // 'droplet//impingement//clean//'
  else
   FlwCalInDir = bckdir // 'flow//cal//icing//'
   DrpImpDir   = bckdir // 'droplet//impingement//icing//'
 endif
 GrdInDir    = bckdir // 'grid//icing//'
 OSGDir      = bckdir // 'overset//icing//'
 FlwIniDir   = bckdir // 'flow//initial//icing//'
! FlwCalInDir = bckdir // 'flow//cal//icing//'
! DrpImpDir   = bckdir // 'droplet//impingement//icing//'
 IceCalInDir = bckdir // 'icing//cal//'
 VldDir      = bckdir // 'validation//'
 write(*, '(a)') '<< Icing Step >>'
 write(*, '(a,i2)') '* Ice step      = ', IceStep
 write(*, '(a,i2)') '* Ice step max. = ', IceStepMax
 write(*, '(/,a)') '<< Exp. Condition >>'
 write(*, '(a,e16.8e3)') '* Ts    = ', TsExp * aRef**2
 write(*, '(a,e16.8e3)') '* Ps    = ', PsExp * (rhoRef * aRef**2)
 write(*, '(a,e16.8e3)') '* V     = ', VelExp * aRef
 write(*, '(a,e16.8e3)') '* LWC   = ', LWC * RhoRef
 write(*, '(a,e16.8e3)') '* MVD   = ', MVD * LRef
 write(*, '(a,e16.8e3)') '* Rho   = ', Rhod * RhoRef
 write(*, '(a,e16.8e3)') '* Chord = ', Chord * LRef
 write(*, '(a,e16.8e3)') '* AOA   = ', AOA * 180.0 / pi
 write(IceStepName, '(a, 2(i2.2,a))') 'IceStep', IceStep, 'of', IceStepMax, '_'
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
 integer   :: m
 character :: fname * 30
 ! 処理開始 ********************************************************************************************
 ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( Flw(ms:me), Vld(ms:me) )
 ! ブロック名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Set_BlockName
 ! 格子解像度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
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
  &         Flw(m)%mut (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
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
! enddo
 ! 物理量に変換 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call SetPhysics3DKEM( &
  &      Rg, gamma, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
  &      Flw(m)%qh, Flw(m)%jac, &
  &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps )
! enddo
 ! 着氷領域 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_Resolution1D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(IceRslFile), strtxt, &
  &      Vld(m)%is, Vld(m)%ie )
  Vld(m)%ks = Flw(m)%ks; Vld(m)%ke = Flw(m)%ke
! enddo
 ! 前縁の格子番号 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  Vld(m)%i0 = 0.5 * (Vld(m)%is + Vld(m)%ie)
  Vld(m)%j0 = Flw(m)%js; Vld(m)%k0 = nint(0.5 * Flw(m)%ke)
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine InitialSetting
!*******************************************************************************************************
!******** 翼面に沿った前縁からの距離 								********
!*******************************************************************************************************
subroutine DistanceLE
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, m
 real    :: disMax
 ! 処理開始 ********************************************************************************************
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  allocate( Vld(m)%DisLE( Vld(m)%is: Vld(m)%ie ) )
! enddo
 ! 前縁からの距離 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  ! S.S.
  Vld(m)%DisLE(Vld(m)%i0) = 0.0
  do i = Vld(m)%i0 + 1, Vld(m)%ie
   Vld(m)%DisLE(i) = Vld(m)%DisLE(i-1) &
   &               + sqrt( ( Flw(m)%x(i, Vld(m)%j0, Vld(m)%k0) - Flw(m)%x(i-1, Vld(m)%j0, Vld(m)%k0) )**2 &
   &               +       ( Flw(m)%y(i, Vld(m)%j0, Vld(m)%k0) - Flw(m)%y(i-1, Vld(m)%j0, Vld(m)%k0) )**2 )
  enddo
  ! P.S.
  do i = Vld(m)%i0 - 1, Vld(m)%is, -1
   Vld(m)%DisLE(i) = Vld(m)%DisLE(i+1) &
   &               + sqrt( ( Flw(m)%x(i, Vld(m)%j0, Vld(m)%k0) - Flw(m)%x(i+1, Vld(m)%j0, Vld(m)%k0) )**2 &
   &               +       ( Flw(m)%y(i, Vld(m)%j0, Vld(m)%k0) - Flw(m)%y(i+1, Vld(m)%j0, Vld(m)%k0) )**2 )
  enddo
! enddo
 ! 無次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 disMax = max( maxval(Vld(1)%DisLE(:)), maxval(Vld(2)%DisLE(:)) )
! do m = ms, me
 m = me
  do i = Vld(m)%is, Vld(m)%ie
   Vld(m)%DisLE(i) = Vld(m)%DisLE(i) / disMax
  enddo
  do i = Vld(m)%is, Vld(m)%i0
   Vld(m)%DisLE(i) = Vld(m)%DisLE(i) * (- 1.0)
  enddo
! enddo
 ! 処理開始 ********************************************************************************************
 return
end subroutine DistanceLE
!*******************************************************************************************************
!******** y+											********
!*******************************************************************************************************
subroutine CalYp
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: kappa = 0.40, B = 5.5
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, m
 integer :: j0, j1, k0
 real    :: dx, dy, dz, disp, mup, nup, utau
 ! 処理開始 ********************************************************************************************
! do m = ms, me
 m = me
  ! メモリ確保 -----------------------------------------------------------------------------------------
  allocate( Vld(m)%yp( Vld(m)%is:Vld(m)%ie ) )
  ! y+ -------------------------------------------------------------------------------------------------
  j0 = Flw(m)%js; j1 = Flw(m)%je
  k0 = int( 0.5 * Flw(m)%ke )
  do i = Vld(m)%is, Vld(m)%ie
   dx   = Flw(m)%x(i,j1,k0) - Flw(m)%x(i,j0,k0)
   dy   = Flw(m)%y(i,j1,k0) - Flw(m)%y(i,j0,k0)
   dz   = Flw(m)%z(i,j1,k0) - Flw(m)%z(i,j0,k0)
   disp = sqrt(dx**2 + dy**2 + dz**2)
   mup  = muSth * (Flw(m)%t(i,j1,k0) / TsSth)**1.5 * (TsSth + s1) / (Flw(m)%t(i,j1,k0) + s1)
   nup  = mup / Flw(m)%rho(i,j1,k0)
   utau = cmu / (kappa * disp) * Flw(m)%kin(i,j1,k0)**2 / Flw(m)%eps(i,j1,k0)
   Vld(m)%yp(i) = disp * utau / nup * ReRef
  enddo
  ! ファイル出力 ---------------------------------------------------------------------------------------
  open(1, file = trim(VldDir) // trim(BlkName(m)) // trim(YpFile) // strtxt)
   do i = Vld(m)%is, Vld(m)%ie
    write(1, '(e16.8e3, x, e16.8e3)') Vld(m)%DisLE(i), Vld(m)%yp(i)
   enddo
  close(1)
  ! メモリ解放 -----------------------------------------------------------------------------------------
  deallocate(Vld(m)%yp)
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalYp
!*******************************************************************************************************
!******** 静温分布 										********
!*******************************************************************************************************
subroutine CalStaticTemperature
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, j, m
 real      :: disMax
 real      :: t0, mach, c
 character :: fname * 20
 ! 処理開始 ********************************************************************************************
 ! Ngraph ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  ! メモリ確保 -----------------------------------------------------------------------------------------
  allocate( Vld(m)%Ts( Vld(m)%is: Vld(m)%ie ) )
  ! 有次元化 -------------------------------------------------------------------------------------------
  do i = Vld(m)%is, Vld(m)%ie
   Vld(m)%Ts(i) = Flw(m)%t(i, Vld(m)%j0, Vld(m)%k0) * aRef**2 - IceTs
  enddo
  ! ファイル出力 ---------------------------------------------------------------------------------------
  open(1, file = trim(VldDir) // trim(BlkName(m)) // trim(TsFile) // strtxt)
   do i = Vld(m)%is, Vld(m)%ie
    write(1, '(e16.8e3, 2x, e16.8e3)') Vld(m)%DisLE(i), Vld(m)%Ts(i)
   enddo
  close(1)
  ! メモリ解放 -----------------------------------------------------------------------------------------
! enddo
 ! 淀み点温度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! マッハ数より導出 ------------------------------------------------------------------------------------
 c    = sqrt( gamma * Rg * TsExp )
 mach = VelExp / c
 t0   = TsExp * (1.0 + 0.5 * (gamma - 1.0) * mach**2) - IceTs
 ! ファイル出力 ---------------------------------------------------------------------------------------
 open(1, file = trim(VldDir) // trim(T0File) // strtxt)
  write(1, *) t0, maxval(Vld(2)%Ts(:))
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalStaticTemperature
!*******************************************************************************************************
!******** 収集効率分布										*******
!*******************************************************************************************************
subroutine CalCollectionEfficiency
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real     , pointer :: CE(:, :), ui(:, :), vi(:, :), wi(:, :)
 integer   :: m, i
 character :: fname * 30
 ! 処理開始 ********************************************************************************************
 ! Ngraph ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
!  if(m == 1) cycle
  ! メモリ確保 -----------------------------------------------------------------------------------------
  allocate( Vld(m)%CE( Vld(m)%is: Vld(m)%ie ) )
  allocate( CE( Vld(m)%is: Vld(m)%ie, Vld(m)%ks: Vld(m)%ke ), &
  &         ui( Vld(m)%is: Vld(m)%ie, Vld(m)%ks: Vld(m)%ke ), &
  &         vi( Vld(m)%is: Vld(m)%ie, Vld(m)%ks: Vld(m)%ke ), &
  &         wi( Vld(m)%is: Vld(m)%ie, Vld(m)%ks: Vld(m)%ke ) )
  ! ファイル入力 ---------------------------------------------------------------------------------------
  call Input_CollectionEfficiency3D( &
  &      trim(DrpImpDir) // trim(BlkName(m)) // trim(CollectionFile), strbin, &
  &      Vld(m)%is, Vld(m)%ie, Flw(m)%ks, Flw(m)%ke, &
  &      CE, ui, vi, wi  )
!  call Input_DropCalLog( &
!  &      trim(DrpImpDir) // trim(DrpCalLogFile) // strtxt )
!  write(fname, '(a,i8.8,a)') 'Drop', nDrop, '_'
!  call Input_CollectionEfficiency3D( &
!  &      trim(DrpImpDir) // trim(fname) // trim(BlkName(m)) // trim(CollectionFile), strbin, &
!  &      Vld(m)%is, Vld(m)%ie, Flw(m)%ks, Flw(m)%ke, &
!  &      CE, ui, vi, wi  )
  ! 配列整理 -------------------------------------------------------------------------------------------
  do i = Vld(m)%is, Vld(m)%ie
   Vld(m)%CE(i) = CE(i, Vld(m)%k0)
  enddo
  ! ファイル出力 ---------------------------------------------------------------------------------------
  open(1, file = trim(VldDir) // trim(CEFile) // strtxt)
  open(2, file = trim(VldDir) // trim(IceStepName) // trim(CEFile) // strtxt)
   do i = Vld(m)%is, Vld(m)%ie
    write(1, '(e16.8e3, 2x, e16.8e3)') Vld(m)%DisLE(i), Vld(m)%CE(i)
    write(2, '(e16.8e3, 2x, e16.8e3)') Vld(m)%DisLE(i), Vld(m)%CE(i)
   enddo
  close(1); close(2)
  ! メモリ解放 -----------------------------------------------------------------------------------------
  deallocate(CE, ui, vi, wi)
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalCollectionEfficiency
!*******************************************************************************************************
!******** 氷層厚さ 										********
!*******************************************************************************************************
subroutine CalIceThickness
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, pointer :: f(:)
 real   , pointer :: dBi(:)
 integer   :: m, i
 real      :: aa
 ! 処理開始 ********************************************************************************************
 ! Ngraph ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
!  if(m == 1) cycle
  ! メモリ確保 -----------------------------------------------------------------------------------------
  allocate( Vld(m)%Bi( Vld(m)%is: Vld(m)%ie ), Vld(m)%Ti( Vld(m)%is: Vld(m)%ie ) )
  allocate( f(Vld(m)%is: Vld(m)%ie), dBi(Vld(m)%is: Vld(m)%ie) )
  ! ファイル入力 ---------------------------------------------------------------------------------------
  call Input_IceThickTem2D( &
  &      trim(IceCalInDir) // trim(BlkName(m)) // trim(IceThickTemFile), strdat, &
  &      Vld(m)%is, Vld(m)%ie, &
  &      f, Vld(m)%Bi, dBi, Vld(m)%Ti )
  ! ファイル出力 ---------------------------------------------------------------------------------------
  open(1, file = trim(VldDir) // trim(IceThickFile) // strtxt)
  open(2, file = trim(VldDir) // trim(IceStepName) // trim(IceThickFile) // strtxt)
   do i = Vld(m)%is, Vld(m)%ie
    write(1, '(e16.8e3, x, e16.8e3)') Vld(m)%DisLE(i), Vld(m)%Bi(i) * 1.0e3
    write(2, '(e16.8e3, x, e16.8e3)') Vld(m)%DisLE(i), Vld(m)%Bi(i) * 1.0e3
   enddo
  close(1); close(2)
  ! メモリ解放 -----------------------------------------------------------------------------------------
  deallocate(f, dBi)
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalIceThickness
!*******************************************************************************************************
!******** 着氷翼の可視化 									********
!*******************************************************************************************************
subroutine ViewIcingBlade
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: xi(:, :), yi(:, :), zi(:, :)
 integer :: i, m
 ! 処理開始 ********************************************************************************************
 ! Ngraph ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
!  if(m == 1) cycle
  ! メモリ確保 -----------------------------------------------------------------------------------------
  allocate( Vld(m)%xi( Vld(m)%is: Vld(m)%ie ), Vld(m)%yi( Vld(m)%is: Vld(m)%ie ) )
  allocate( xi( Vld(m)%is: Vld(m)%ie, Vld(m)%ks: Vld(m)%ke ), &
  &         yi( Vld(m)%is: Vld(m)%ie, Vld(m)%ks: Vld(m)%ke ), &
  &         zi( Vld(m)%is: Vld(m)%ie, Vld(m)%ks: Vld(m)%ke ) )
  ! 配列整理 -------------------------------------------------------------------------------------------
  do i = Vld(m)%is, Vld(m)%ie
    Vld(m)%xi(i) = Flw(m)%x(i,Vld(m)%j0,Vld(m)%k0)
    Vld(m)%yi(i) = Flw(m)%y(i,Vld(m)%j0,Vld(m)%k0)
  enddo
  ! ファイル出力 ---------------------------------------------------------------------------------------
  open(1, file = trim(VldDir) // trim(IcingBladeFile) // strtxt)
  open(2, file = trim(VldDir) // trim(IceStepName) // trim(IcingBladeFile) // strtxt)
   do i = Vld(m)%is, Vld(m)%ie
    write(1, '(e16.8e3, x, e16.8e3)') Vld(m)%xi(i) / Chord, Vld(m)%yi(i) / Chord
    write(2, '(e16.8e3, x, e16.8e3)') Vld(m)%xi(i) / Chord, Vld(m)%yi(i) / Chord
   enddo
  close(1); close(2)
  ! メモリ解放 -----------------------------------------------------------------------------------------
  deallocate( xi, yi, zi )
  deallocate(Vld(m)%xi, Vld(m)%yi)
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine ViewIcingBlade
! 定義終了 *********************************************************************************************
end program Validation_NACA
