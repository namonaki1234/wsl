!*******************************************************************************************************
!*******************************************************************************************************
!******** 流れ場初期解プログラム								********
!******** (NACA翼，三次元圧縮性乱流場)								********
!********					      2013.02.02  PROGRAMMED BY RYOSUKE HAYASHI ********
!********					      2014.04.15     UPDATED BY RYOSUKE HAYASHI ********
!*******************************************************************************************************
!*******************************************************************************************************
program InitialFlow_NACA
 ! モジュール宣言 **************************************************************************************
 use Package_NACA
 use Package_FileIO
 use Package_Grid
 use Package_Flow
 use Package_ViewFlow
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! ファイル名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: ViewFlowFile  *  8 = 'ViewFlow'
 ! 共有定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, parameter :: nSmooth = 100					! スムージング回数
 ! 処理開始 ********************************************************************************************
 ! 検証実験ケース選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call SelectExpCase
 write(*, '(a)') "+++++ Select Exp. case complete. +++++"
 ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call InitialSetting
 write(*, '(a)') "+++++ Initial setting complete. +++++"
 ! メトリックス ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call CalMetrics
 write(*, '(a)') "+++++ Metrics clculation complete. +++++"
 ! 初期条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call InitialCondition
 write(*, '(a)') "+++++ Initial condition complete. +++++"
 ! 境界条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call BoundaryCondition
 write(*, '(a)') "+++++ Boundary condition complete. +++++"
! ! 重合格子の補間 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call InterpolationOversetGrid
 write(*, '(a)') "+++++ Interpolation overset condition complete. +++++"
 ! スムージング ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Smoothing
 write(*, '(a)') "+++++ Smoothing complete. +++++"
 ! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call OutputFile
 write(*, '(a)') "+++++ Output file complete. +++++"
 ! 無次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call CalNondimensionalized
 write(*, '(a)') "+++++ Nondimensionalized complete. +++++"
 ! 無次元化ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call OutputNDFile
 write(*, '(a)') "+++++ Output ND file complete. +++++"
 ! 可視化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call VisualizedFlow
 write(*, '(a)') "+++++ Visualized flow complete. +++++"
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
 ! ディレクトリ設定
 GrdInDir    = bckdir // 'grid/clean/'
 OSGDir      = bckdir // 'overset/clean/'
 FlwIniDir   = bckdir // 'flow/initial//clean/'
 ! 計算条件ファイル入力
 call Input_CalSetting( trim(ND_CalSetFile) // strtxt )
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
 integer :: m
 ! 処理開始 ********************************************************************************************
 ! ブロック名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Set_BlockName
 ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! allocate( Flw(ms:me), OSG(ms:me) )
 allocate( Flw(ms:me) )
 ! メモリ確保及び格子ファイル ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  ! 格子解像度 -----------------------------------------------------------------------------------------
  call Input_Resolution3D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(RslFile), strtxt, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke )
  ! メモリ確保 -----------------------------------------------------------------------------------------
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
  ! 格子座標 -------------------------------------------------------------------------------------------
  call Input_Grid3D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(GrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
  ! C 型格子分割点 -------------------------------------------------------------------------------------
  call Input_CtypeGridPoint( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
  &      Flw(m)%i1, Flw(m)%i2, Flw(m)%i3 )
! enddo
! ! 重合格子補間係数ファイル ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
!  ! 補間探索範囲 ---------------------------------------------------------------------------------------
!  call Input_Resolution3D( &
!  &      trim(OSGDir) // trim(BlkName(m)) // trim(OverAreaFile), strtxt, &
!  &      OSG(m)%is, OSG(m)%ie, OSG(m)%js, OSG(m)%je, OSG(m)%ks, OSG(m)%ke )
!  ! メモリ確保 -----------------------------------------------------------------------------------------
!  allocate( OSG(m)%ip   ( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
!  &         OSG(m)%jp   ( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
!  &         OSG(m)%kp   ( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
!  &         OSG(m)%fOver( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
!  &         OSG(m)%term1( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
!  &         OSG(m)%term2( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
!  &         OSG(m)%term3( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
!  &         OSG(m)%term4( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
!  &         OSG(m)%term5( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
!  &         OSG(m)%term6( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
!  &         OSG(m)%term7( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ), &
!  &         OSG(m)%term8( OSG(m)%is: OSG(m)%ie, OSG(m)%js: OSG(m)%je, OSG(m)%ks: OSG(m)%ke ) )
!  ! 補間係数 -------------------------------------------------------------------------------------------
!  call Input_OversetCoe3D( &
!  &      trim(OSGDir) // trim(BlkName(m)) // trim(OverCoeFile), strbin, &
!  &      OSG(m)%is, OSG(m)%ie, OSG(m)%js, OSG(m)%je, OSG(m)%ks, OSG(m)%ke, &
!  &      OSG(m)%fOver, OSG(m)%ip, OSG(m)%jp, OSG(m)%kp, &
!  &      OSG(m)%term1, OSG(m)%term2, OSG(m)%term3, OSG(m)%term4, &
!  &      OSG(m)%term5, OSG(m)%term6, OSG(m)%term7, OSG(m)%term8 )
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine InitialSetting
!*******************************************************************************************************
!******** メトリックス計算									********
!*******************************************************************************************************
subroutine CalMetrics
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m
 ! 処理開始 ********************************************************************************************
! do m = ms, me
 m = me
  call Metrics3D( &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z, &
  &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
  &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
  &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
  &      Flw(m)%jac )
! enddo
 ! 処理終了 ********************************************************************************************
return
end subroutine CalMetrics
!*******************************************************************************************************
!******** 初期条件設定 										********
!*******************************************************************************************************
subroutine InitialCondition
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k, m
 integer :: j0
 real    :: RhoIn, Mach
 real    :: BL
 ! 処理開始 ********************************************************************************************
 ! 流入条件及び流出条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Mach  = VelExp / sqrt(gamma * Rg * TsExp)
 VelIn = VelExp
 RhoIn = PsExp / (Rg * TsExp)
 TtIn  = TsExp * (1.0 + 0.5 * (gamma - 1.0) * Mach**2)
 PtIn  = PsExp * (1.0 + 0.5 * (gamma - 1.0) * Mach**2)**(gamma / (gamma - 1.0))
 TsOut = TsExp
 PsOut = PsExp
 ! 初期条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  ! 静温，静圧
  do k = Flw(m)%ks, Flw(m)%ke
  do j = Flw(m)%js, Flw(m)%je
  do i = Flw(m)%is, Flw(m)%ie
   Flw(m)%t  (i,j,k) = TsExp
   Flw(m)%p  (i,j,k) = PsExp
  enddo
  enddo
  enddo
  ! 粘性係数
  call ViscosityCoefficient3D( &
  &      muSth, TsSth, s1, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%t, Flw(m)%mu )
  ! 密度，速度，乱流量
  do k = Flw(m)%ks, Flw(m)%ke
  do j = Flw(m)%js, Flw(m)%je
  do i = Flw(m)%is, Flw(m)%ie
   Flw(m)%rho(i,j,k) = PsExp / (Rg * TsExp)
   Flw(m)%u  (i,j,k) = VelExp * cos(AOA)
   Flw(m)%v  (i,j,k) = VelExp * sin(AOA)
   Flw(m)%w  (i,j,k) = 0.0
   Flw(m)%kin(i,j,k) = 1.5 * ( 0.01 * sqrt( Flw(m)%u(i,j,k)**2 &
   &                                      + Flw(m)%v(i,j,k)**2 &
   &                                      + Flw(m)%w(i,j,k)**2) )**2
   Flw(m)%eps(i,j,k) = Flw(m)%kin(i,j,k)**2 * Flw(m)%rho(i,j,k) / (Flw(m)%mu(i,j,k) * Ret)
   Flw(m)%mut(i,j,k) = Cmu * Flw(m)%rho(i,j,k) * Flw(m)%kin(i,j,k)**2 / Flw(m)%eps(i,j,k)
  enddo
  enddo
  enddo
! enddo
 ! 1/7 乗則 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 境界層厚さ ------------------------------------------------------------------------------------------
 BL = 0.01 * chord
 ! 速度分布 --------------------------------------------------------------------------------------------
! do m = ms, me
 m = me
  j0 = Flw(m)%js
  do k = Flw(m)%ks, Flw(m)%ke
  do j = Flw(m)%js, Flw(m)%je
  do i = Flw(m)%is, Flw(m)%ie
    if( Flw(m)%i1 <= i .and. i <= Flw(m)%i3 ) then
      Flw(m)%u(i,j,k) = Flw(m)%u(i,j,k) &
      &               * min( 1.0, ( sqrt( (Flw(m)%x(i,j,k) - Flw(m)%x(i,j0,k))**2 &
      &                                 + (Flw(m)%y(i,j,k) - Flw(m)%y(i,j0,k))**2 &
      &                                 + (Flw(m)%z(i,j,k) - Flw(m)%z(i,j0,k))**2 ) &
      &                           / BL )**(1.0 / 7.0) )
    else
      Flw(m)%u(i,j,k) = Flw(m)%u(i,j,k) &
      &               * min( 1.0, ( sqrt( (Flw(m)%x(i,j,k) - Flw(m)%x(Flw(m)%i1,j0,k))**2 &
      &                                 + (Flw(m)%y(i,j,k) - Flw(m)%y(Flw(m)%i1,j0,k))**2 &
      &                                 + (Flw(m)%z(i,j,k) - Flw(m)%z(Flw(m)%i1,j0,k))**2 ) &
      &                           / BL )**(1.0 / 7.0) )
    endif
  enddo
  enddo
  enddo
! enddo
 ! 計算量に変換 ----------------------------------------------------------------------------------------
! do m = ms, me
 m = me
  call SetFlux3DKEM( &
  &      gamma, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
  &      Flw(m)%qh, Flw(m)%jac, &
  &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%kin, Flw(m)%eps )
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine InitialCondition
!*******************************************************************************************************
!******** 境界条件 										********
!*******************************************************************************************************
subroutine BoundaryCondition
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k, l, m
 ! 処理開始 ********************************************************************************************
 ! Main Grid, Sub Grid +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms,me
 m = me
 ! 流入境界 --------------------------------------------------------------------------------------------
 call BoundaryInlet( &
 &      m, Flw(m)%is + 1, Flw(m)%ie - 1, Flw(m)%ks + 1, Flw(m)%ke - 1, &
 &      Flw(m)%je - 0, Flw(m)%je - 1 )
 ! C 型格子ブランチ・カット ----------------------------------------------------------------------------
 call BoundaryCtypeBranch( &
 &      m, Flw(m)%is + 1, Flw(m)%i1 - 1, Flw(m)%ks + 1, Flw(m)%ke - 1, &
 &      Flw(m)%js + 0, Flw(m)%js + 1, &
 &      Flw(m)%ie )
 call BoundaryCtypeBranch( &
 &      m, Flw(m)%i3 + 1, Flw(m)%ie - 1, Flw(m)%ks + 1, Flw(m)%ke - 1, &
 &      Flw(m)%js + 0, Flw(m)%js + 1, &
 &      Flw(m)%ie )
 ! 翼壁面 ----------------------------------------------------------------------------------------------
 call BoundaryBladeSurface( &
 &      m, Flw(m)%i1 + 0, Flw(m)%i3 - 0, Flw(m)%ks + 1, Flw(m)%ke - 1, &
 &      Flw(m)%js + 0, Flw(m)%js + 1, Flw(m)%js + 2, TurbNum )
 ! 流出境界 --------------------------------------------------------------------------------------------
 call BoundaryOutlet( &
 &      m, Flw(m)%js + 0, Flw(m)%je - 0, Flw(m)%ks + 1, Flw(m)%ke - 1, &
 &      Flw(m)%is + 0, Flw(m)%is + 1 )
 call BoundaryOutlet( &
 &      m, Flw(m)%js + 0, Flw(m)%je - 0, Flw(m)%ks + 1, Flw(m)%ke - 1, &
 &      Flw(m)%ie - 0, Flw(m)%ie - 1 )
 ! 周期境界 --------------------------------------------------------------------------------------------
 call BoundaryPeriodic( &
 &      m, Flw(m)%is + 0, Flw(m)%ie - 0, Flw(m)%js + 0, Flw(m)%je - 0, &
 &      Flw(m)%ks + 0, Flw(m)%ke - 1 )
 call BoundaryPeriodic( &
 &      m, Flw(m)%is + 0, Flw(m)%ie - 0, Flw(m)%js + 0, Flw(m)%je - 0, &
 &      Flw(m)%ke - 0, Flw(m)%ks + 1 )
! end do
 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryCondition
!*******************************************************************************************************
!******** 流入境界										********
!******** ⅰ．密度外挿，その他固定　			ⅱ．角度・全温・全圧固定，マッハ数外挿	********
!******** ⅲ．角度・体積流量・全温固定，密度外挿  						********
!*******************************************************************************************************
subroutine BoundaryInlet( &
&            m, is, ie, ks, ke, j0, j1 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 integer, intent(in) :: is, ie, ks, ke
 integer, intent(in) :: j0, j1
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k, l
 real    :: u, v, w, p, t
 real    :: Tt_Ts, Ts, Ps, mac, vel, rho, kin, eps
 ! 処理開始 ********************************************************************************************
 select case(BCNum)
  case(1)
   do k = ks, ke
   do i = is, ie
    ! 密度外挿
    rho = Flw(m)%qh(i,j1,k,1) * Flw(m)%jac(i,j1,k)
    ! その他固定
    vel = VelIn
    u   = vel * cos(AOA)
    v   = vel * sin(AOA)
    w   = 0.0
    t   = TsExp
    kin = 0.5 * 3.0 * ( 0.01 * vel )**2
    eps = rho * kin**2 / (Flw(m)%mu(i,j0,k) * Ret)
    ! 流束関数
    Flw(m)%qh(i,j0,k,1) = rho / Flw(m)%jac(i,j0,k)
    Flw(m)%qh(i,j0,k,2) = u * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,3) = v * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,4) = w * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,5) = ( Rg * t / (gamma - 1.0) + 0.5 * vel**2 + kin ) * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,6) =  kin * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,7) =  eps * Flw(m)%qh(i,j0,k,1)
   enddo
   enddo
  case(2)
   do k = ks, ke
   do i = is, ie
    ! 外挿する速度からマッハ数導出
    u = Flw(m)%qh(i,j1,k,2) / Flw(m)%qh(i,j1,k,1)
    v = Flw(m)%qh(i,j1,k,3) / Flw(m)%qh(i,j1,k,1)
    w = Flw(m)%qh(i,j1,k,4) / Flw(m)%qh(i,j1,k,1)
    p = ( Flw(m)%qh(i,j1,k,5) - Flw(m)%qh(i,j1,k,6) &
    &     - 0.5 * ( Flw(m)%qh(i,j1,k,2)**2 &
    &             + Flw(m)%qh(i,j1,k,3)**2 &
    &             + Flw(m)%qh(i,j1,k,4)**2 ) / Flw(m)%qh(i,j1,k,1) &
    &   ) * (gamma - 1.0) * Flw(m)%jac(i,j1,k)
    t = p / (Flw(m)%qh(i,j1,k,1) * Flw(m)%jac(i,j1,k) * Rg)
    mac = sqrt( u**2 + v**2 + w**2 ) / sqrt( gamma * Rg * t )
    ! 外挿するマッハ数と固定する全温・全圧から静温・静圧導出
    Tt_Ts = 1.0 + 0.5 * (gamma - 1.0) * mac**2
    Ts = TtIn / Tt_Ts
    Ps = PtIn / Tt_Ts**( gamma / (gamma - 1.0) )
    ! その他導出
    rho = Ps / ( Rg * Ts )
    vel = mac * sqrt(gamma * Rg * Ts)
    u   = vel * cos(AOA)
    v   = vel * sin(AOA)
    w   = 0.0
    kin = 0.5 * 3.0 * ( 0.01 * vel )**2
    eps = rho * kin**2 / (Flw(m)%mu(i,j0,k) * Ret)
    ! 流束関数
    Flw(m)%qh(i,j0,k,1) = rho / Flw(m)%jac(i,j0,k)
    Flw(m)%qh(i,j0,k,2) = u   * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,3) = v   * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,4) = w   * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,5) = ( Rg * Ts / (gamma - 1.0) + 0.5 * vel**2 + kin ) * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,6) =  kin * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,7) =  eps * Flw(m)%qh(i,j0,k,1)
   enddo
   enddo
  case(3)
   do k = ks, ke
   do i = is, ie
    ! 密度外挿
    rho = Flw(m)%qh(i,j1,k,1) * Flw(m)%jac(i,j1,k)
    ! 体積流量固定
    vel = VelIn
    ! 全温固定
    Ts  = TtIn - (gamma - 1.0) / (2.0 * gamma * Rg) * vel**2
    Ps  = rho * Rg * Ts
    ! 速度から乱流量導出
    kin = 0.5 * 3.0 * ( 0.01 * vel )**2
    eps = rho * kin**2 / (Flw(m)%mu(i,j0,k) * Ret)
    ! 流束関数
    Flw(m)%qh(i,j0,k,1) = rho / Flw(m)%jac(i,j0,k)
    Flw(m)%qh(i,j0,k,2) = vel * cos(AOA) * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,3) = vel * sin(AOA) * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,4) = 0.0
    Flw(m)%qh(i,j0,k,5) = ( Rg * Ts / (gamma - 1.0) + 0.5 * vel**2 + kin ) * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,6) =  kin * Flw(m)%qh(i,j0,k,1)
    Flw(m)%qh(i,j0,k,7) =  eps * Flw(m)%qh(i,j0,k,1)
   enddo
   enddo
  case default
   write(*, '(a)') '!!!! Error : Inlet boundary condition pattern !!!!'
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryInlet
!*******************************************************************************************************
!******** 流出境界 (静圧固定，その他外挿) 							********
!*******************************************************************************************************
subroutine BoundaryOutlet( &
&            m, js, je, ks, ke, i0, i1 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 integer, intent(in) :: js, je, ks, ke
 integer, intent(in) :: i0, i1
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: j, k, l
 real    :: Ts, Ps
 ! 処理開始 ********************************************************************************************
 do k = ks, ke
 do j = js, je
  Ps = ( Flw(m)%qh(i1,j,k,5) - Flw(m)%qh(i1,j,k,6) &
  &      - 0.5 * ( Flw(m)%qh(i1,j,k,2)**2 &
  &              + Flw(m)%qh(i1,j,k,3)**2 &
  &              + Flw(m)%qh(i1,j,k,4)**2 ) / Flw(m)%qh(i1,j,k,1) &
  &    ) * (gamma - 1.0) * Flw(m)%jac(i1,j,k)
  Ts = Ps / (Flw(m)%qh(i1,j,k,1) * Flw(m)%jac(i1,j,k) * Rg)
  do l = ls + 1, le
   Flw(m)%qh(i0,j,k,l) = Flw(m)%qh(i1,j,k,l) / Flw(m)%qh(i1,j,k,1)
  enddo
  Flw(m)%qh(i0,j,k,1) = PsOut / ( Rg * Ts * Flw(m)%jac(i0,j,k) )
  do l = ls + 1, le
    Flw(m)%qh(i0,j,k,l) = Flw(m)%qh(i0,j,k,l) * Flw(m)%qh(i0,j,k,1)
  enddo
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryOutlet
!*******************************************************************************************************
!******** C 型格子ブランチ・カット境界 (平均値外挿)                           			********
!*******************************************************************************************************
subroutine BoundaryCtypeBranch( &
&            m, is, ie, ks, ke, j0, j1, imax )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 integer, intent(in) :: is, ie, ks, ke
 integer, intent(in) :: j0, j1
 integer, intent(in) :: imax
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k, l
 real    :: qAve
 ! 処理開始 ********************************************************************************************
 do l = ls, le
 do k = ks, ke
 do i = is, ie
  qAve = 0.5 * ( Flw(m)%qh(i     ,j1,k,l) * Flw(m)%jac(i     ,j1,k) &
  &            + Flw(m)%qh(imax-i,j1,k,l) * Flw(m)%jac(imax-i,j1,k) )
  Flw(m)%qh(i,j0,k,l) = qAve / Flw(m)%jac(i,j0,k)
 enddo
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryCtypeBranch
!*******************************************************************************************************
!******** 翼壁面境界										********
!*******************************************************************************************************
subroutine BoundaryBladeSurface( &
&            m, is, ie, ks, ke, j0, j1, j2, TurbNum )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 integer, intent(in) :: is, ie, ks, ke
 integer, intent(in) :: j0, j1, j2
 integer, intent(in) :: TurbNum
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 ! 処理開始 ********************************************************************************************
 do k = ks, ke
 do i = is, ie
  if(TurbNum .eq. 4) then
    call WallNoSlipLRe(m, i, j0, k, i, j1, k)
  else
    if(fSlip) then
      if( Flw(m)%qh(i,j1,k,1) <= 0.0 .or. Flw(m)%qh(i,j2,k,1) <= 0.0 ) cycle
      call WallSlipWF(m, i, j0, k, i, j1, k, i, j2, k)
    else
      if( Flw(m)%qh(i,j1,k,1) <= 0.0 ) cycle
      call WallNoSlipWF(m, i, j0, k, i, j1, k)
    endif
  endif
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryBladeSurface
!*******************************************************************************************************
!******** 壁境界 (滑りなし, 断熱, 壁関数)                                			********
!*******************************************************************************************************
subroutine WallNoslipWF( &
&            m, i0, j0, k0, i1, j1, k1 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 integer, intent(in) :: i0, j0, k0
 integer, intent(in) :: i1, j1, k1
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real :: dx, dy, dz
 real :: uc, vc, wc
 real :: yp, up, nup, kp, epsp, utau, dudy1, dudy2
 ! 処理開始 ********************************************************************************************
 ! 壁関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁近傍の分布
 dx = Flw(m)%x(i1,j1,k1) - Flw(m)%x(i0,j0,k0)
 dy = Flw(m)%y(i1,j1,k1) - Flw(m)%y(i0,j0,k0)
 dz = Flw(m)%z(i1,j1,k1) - Flw(m)%z(i0,j0,k0)
 uc = Flw(m)%qh(i1,j1,k1,2) / Flw(m)%qh(i1,j1,k1,1)
 vc = Flw(m)%qh(i1,j1,k1,3) / Flw(m)%qh(i1,j1,k1,1)
 wc = Flw(m)%qh(i1,j1,k1,4) / Flw(m)%qh(i1,j1,k1,1)
 ! 壁関数
 yp   = sqrt(dx**2 + dy**2 + dz**2)
 up   = sqrt(uc**2 + vc**2 + wc**2)
 up   = max(zero, up)
 nup  = Flw(m)%mu(i1,j1,k1) / (Flw(m)%qh(i1,j1,k1,1) * Flw(m)%jac(i1,j1,k1))
 kp   = Flw(m)%kin(i1,j1,k1)
 epsp = Flw(m)%eps(i1,j1,k1)
! call WallFunctionKEM1S( &
! &      yp, up, nup, kp, epsp, utau, dudy1, dudy2, kp, epsp )
 call WallFunctionKEM2S( &
 &      yp, up, nup, utau, dudy1, dudy2, kp, epsp )
 ! 速度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁面
 uc = 0.0
 vc = 0.0
 wc = 0.0
 ! 流束関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁面から一点目
 Flw(m)%qh(i1,j1,k1,5) = Flw(m)%qh(i1,j1,k1,5) - Flw(m)%qh(i1,j1,k1,6)
 Flw(m)%qh(i1,j1,k1,6) = Flw(m)%qh(i1,j1,k1,1) * kp
 Flw(m)%qh(i1,j1,k1,7) = Flw(m)%qh(i1,j1,k1,1) * epsp
 Flw(m)%qh(i1,j1,k1,5) = Flw(m)%qh(i1,j1,k1,5) + Flw(m)%qh(i1,j1,k1,6)
 ! 壁面
 Flw(m)%qh(i0,j0,k0,1) = Flw(m)%jac(i1,j1,k1) / Flw(m)%jac(i0,j0,k0) * Flw(m)%qh(i1,j1,k1,1)
 Flw(m)%qh(i0,j0,k0,2) = Flw(m)%qh(i0,j0,k0,1) * uc
 Flw(m)%qh(i0,j0,k0,3) = Flw(m)%qh(i0,j0,k0,1) * vc
 Flw(m)%qh(i0,j0,k0,4) = Flw(m)%qh(i0,j0,k0,1) * wc
 Flw(m)%qh(i0,j0,k0,5) = Flw(m)%jac(i1,j1,k1) / Flw(m)%jac(i0,j0,k0) &
 &                     * ( Flw(m)%qh(i1,j1,k1,5) &
 &                       - 0.5 * ( Flw(m)%qh(i1,j1,k1,2)**2 &
 &                               + Flw(m)%qh(i1,j1,k1,3)**2 &
 &                               + Flw(m)%qh(i1,j1,k1,4)**2 ) / Flw(m)%qh(i1,j1,k1,1) ) &
 &                     + 0.5 * ( Flw(m)%qh(i0,j0,k0,2)**2 &
 &                             + Flw(m)%qh(i0,j0,k0,3)**2 &
 &                             + Flw(m)%qh(i0,j0,k0,4)**2 ) / Flw(m)%qh(i0,j0,k0,1)
 Flw(m)%qh(i0,j0,k0,6) = Flw(m)%jac(i1,j1,k1) / Flw(m)%jac(i0,j0,k0) * Flw(m)%qh(i1,j1,k1,6)
 Flw(m)%qh(i0,j0,k0,7) = Flw(m)%jac(i1,j1,k1) / Flw(m)%jac(i0,j0,k0) * Flw(m)%qh(i1,j1,k1,7)
 ! 処理終了 ********************************************************************************************
 return
end subroutine WallNoSlipWF
!*******************************************************************************************************
!******** 壁境界 (滑りあり, 断熱, 壁関数)                                			********
!*******************************************************************************************************
subroutine WallSlipWF( &
&            m, i0, j0, k0, i1, j1, k1, i2, j2, k2 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 integer, intent(in) :: i0, j0, k0
 integer, intent(in) :: i1, j1, k1
 integer, intent(in) :: i2, j2, k2
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real :: dx2, dy2, dz2
 real :: uc2, vc2, wc2
 real :: uc1, vc1, wc1
 real :: uc0, vc0, wc0
 real :: yp, up, nup, kp, epsp, utau, dudy1, dudy2
 ! 処理開始 ********************************************************************************************
 ! 壁関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁近傍の分布
 dx2 = Flw(m)%x(i2,j2,k2) - Flw(m)%x(i0,j0,k0)
 dy2 = Flw(m)%y(i2,j2,k2) - Flw(m)%y(i0,j0,k0)
 dz2 = Flw(m)%z(i2,j2,k2) - Flw(m)%z(i0,j0,k0)
 uc2 = Flw(m)%qh(i2,j2,k2,2) / Flw(m)%qh(i2,j2,k2,1)
 vc2 = Flw(m)%qh(i2,j2,k2,3) / Flw(m)%qh(i2,j2,k2,1)
 wc2 = Flw(m)%qh(i2,j2,k2,4) / Flw(m)%qh(i2,j2,k2,1)
 ! 壁関数
 yp   = sqrt(dx2**2 + dy2**2 + dz2**2)
 up   = sqrt(uc2**2 + vc2**2 + wc2**2)
 up   = max(zero, up)
 nup  = Flw(m)%mu(i2,j2,k2) / (Flw(m)%qh(i2,j2,k2,1) * Flw(m)%jac(i2,j2,k2))
 kp   = Flw(m)%kin(i2,j2,k2)
 epsp = Flw(m)%eps(i2,j2,k2)
! call WallFunctionKEM1S( &
! &      yp, up, nup, kp, epsp, utau, dudy1, dudy2, kp, epsp )
 call WallFunctionKEM2S( &
 &      yp, up, nup, utau, dudy1, dudy2, kp, epsp )
 ! 速度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁面から一点目
 uc1 = Flw(m)%qh(i1,j1,k1,2) / Flw(m)%qh(i1,j1,k1,1)
 vc1 = Flw(m)%qh(i1,j1,k1,3) / Flw(m)%qh(i1,j1,k1,1)
 wc1 = Flw(m)%qh(i1,j1,k1,4) / Flw(m)%qh(i1,j1,k1,1)
 yp  = sqrt( (Flw(m)%x(i2,j2,k2) - Flw(m)%x(i1,j1,k1))**2 &
 &         + (Flw(m)%y(i2,j2,k2) - Flw(m)%y(i1,j1,k1))**2 &
 &         + (Flw(m)%z(i2,j2,k2) - Flw(m)%z(i1,j1,k1))**2 )
 up  = max(zero, up - dudy1 * yp) / max(zero, sqrt(uc1**2 + vc1**2 + wc1**2))
 uc1 = uc1 * up
 vc1 = vc1 * up
 wc1 = wc1 * up
 ! 壁面
 uc0 = 0.0
 vc0 = 0.0
 wc0 = 0.0
 ! 流束関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁面から二点目
 Flw(m)%qh(i2,j2,k2,5) = Flw(m)%qh(i2,j2,k2,5) - Flw(m)%qh(i2,j2,k2,6)
 Flw(m)%qh(i2,j2,k2,6) = Flw(m)%qh(i2,j2,k2,1) * kp
 Flw(m)%qh(i2,j2,k2,7) = Flw(m)%qh(i2,j2,k2,1) * epsp
 Flw(m)%qh(i2,j2,k2,5) = Flw(m)%qh(i2,j2,k2,5) + Flw(m)%qh(i2,j2,k2,6)
 ! 壁面から一点目
 Flw(m)%qh(i1,j1,k1,5) = Flw(m)%qh(i1,j1,k1,5) &
 &                    - 0.5 * ( Flw(m)%qh(i1,j1,k1,2)**2 &
 &                            + Flw(m)%qh(i1,j1,k1,3)**2 &
 &                            + Flw(m)%qh(i1,j1,k1,4)**2 ) / Flw(m)%qh(i1,j1,k1,1) &
 &                    - Flw(m)%qh(i1,j1,k1,6)
 Flw(m)%qh(i1,j1,k1,2) = Flw(m)%qh(i1,j1,k1,1) * uc1
 Flw(m)%qh(i1,j1,k1,3) = Flw(m)%qh(i1,j1,k1,1) * vc1
 Flw(m)%qh(i1,j1,k1,4) = Flw(m)%qh(i1,j1,k1,1) * wc1
 Flw(m)%qh(i1,j1,k1,6) = Flw(m)%qh(i1,j1,k1,1) * kp
 Flw(m)%qh(i1,j1,k1,7) = Flw(m)%qh(i1,j1,k1,1) * epsp
 Flw(m)%qh(i1,j1,k1,5) = Flw(m)%qh(i1,j1,k1,5) &
 &                    + 0.5 * ( Flw(m)%qh(i1,j1,k1,2)**2 &
 &                            + Flw(m)%qh(i1,j1,k1,3)**2 &
 &                            + Flw(m)%qh(i1,j1,k1,4)**2 ) / Flw(m)%qh(i1,j1,k1,1) &
 &                    + Flw(m)%qh(i1,j1,k1,6)
 ! 壁面
 Flw(m)%qh(i0,j0,k0,1) = Flw(m)%jac(i1,j1,k1) / Flw(m)%jac(i0,j0,k0) * Flw(m)%qh(i1,j1,k1,1)
 Flw(m)%qh(i0,j0,k0,2) = Flw(m)%qh(i0,j0,k0,1) * uc0
 Flw(m)%qh(i0,j0,k0,3) = Flw(m)%qh(i0,j0,k0,1) * vc0
 Flw(m)%qh(i0,j0,k0,4) = Flw(m)%qh(i0,j0,k0,1) * wc0
 Flw(m)%qh(i0,j0,k0,5) = Flw(m)%jac(i1,j1,k1) / Flw(m)%jac(i0,j0,k0) * ( &
 &                      Flw(m)%qh(i1,j1,k1,5) &
 &                    - 0.5 * ( Flw(m)%qh(i1,j1,k1,2)**2 &
 &                            + Flw(m)%qh(i1,j1,k1,3)**2 &
 &                            + Flw(m)%qh(i1,j1,k1,4)**2 ) / Flw(m)%qh(i1,j1,k1,1) &
 &                    ) &
 &                    + 0.5 * ( Flw(m)%qh(i0,j0,k0,2)**2 &
 &                            + Flw(m)%qh(i0,j0,k0,3)**2 &
 &                            + Flw(m)%qh(i0,j0,k0,4)**2 ) / Flw(m)%qh(i0,j0,k0,1)
 Flw(m)%qh(i0,j0,k0,6) = Flw(m)%jac(i1,j1,k1) / Flw(m)%jac(i0,j0,k0) * Flw(m)%qh(i1,j1,k1,6)
 Flw(m)%qh(i0,j0,k0,7) = Flw(m)%jac(i1,j1,k1) / Flw(m)%jac(i0,j0,k0) * Flw(m)%qh(i1,j1,k1,7)
 ! 処理終了 ********************************************************************************************
 return
end subroutine WallslipWF
!*******************************************************************************************************
!******** 壁境界 (滑りなし, 断熱, 低レイノルズ数型AKNモデル)                  			********
!*******************************************************************************************************
subroutine WallNoslipLRe( &
&            m, i0, j0, k0, i1, j1, k1 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 integer, intent(in) :: i0, j0, k0
 integer, intent(in) :: i1, j1, k1
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real :: dx, dy, dz, yp
 double precision :: epsp
 real :: uc, vc, wc
 ! 処理開始 ********************************************************************************************
 ! 速度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁面 ------------------------------------------------------------------------------------------------
 uc = 0.0
 vc = 0.0
 wc = 0.0
 ! 流束関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁面から一点目 --------------------------------------------------------------------------------------
 dx = Flw(m)%x(i1,j1,k1) - Flw(m)%x(i0,j0,k0)
 dy = Flw(m)%y(i1,j1,k1) - Flw(m)%y(i0,j0,k0)
 dz = Flw(m)%z(i1,j1,k1) - Flw(m)%z(i0,j0,k0)
 yp   = sqrt(dx**2 + dy**2 + dz**2)
 epsp = 2.0d0 * dble(Flw(m)%mu(i1,j1,k1) * Flw(m)%kin(i1,j1,k1)) / &
 &      (dble(Flw(m)%qh(i1,j1,k1,1) * Flw(m)%jac(i1,j1,k1)) * dble(yp)**2.0)
 if(epsp .ne. epsp) then
  write(*,*) 'eps_p is diverged'
  stop
 end if
 ! epsp
 Flw(m)%qh(i1,j1,k1,7) = real(dble(Flw(m)%qh(i1,j1,k1,1)) * epsp)
 ! 壁面 ------------------------------------------------------------------------------------------------
 ! 全て外挿
 Flw(m)%qh(i0,j0,k0,1) = Flw(m)%qh(i1,j1,k1,1) * Flw(m)%jac(i1,j1,k1) / Flw(m)%jac(i0,j0,k0)
 Flw(m)%qh(i0,j0,k0,2) = uc * Flw(m)%qh(i0,j0,k0,1)
 Flw(m)%qh(i0,j0,k0,3) = vc * Flw(m)%qh(i0,j0,k0,1)
 Flw(m)%qh(i0,j0,k0,4) = wc * Flw(m)%qh(i0,j0,k0,1)
 Flw(m)%qh(i0,j0,k0,5) = Flw(m)%jac(i1,j1,k1) / Flw(m)%jac(i0,j0,k0) &
 &                     * ( Flw(m)%qh(i1,j1,k1,5) &
 &                       - 0.5 * ( Flw(m)%qh(i1,j1,k1,2)**2 &
 &                               + Flw(m)%qh(i1,j1,k1,3)**2 &
 &                               + Flw(m)%qh(i1,j1,k1,4)**2 ) / Flw(m)%qh(i1,j1,k1,1) ) &
 &                     + 0.5 * ( Flw(m)%qh(i0,j0,k0,2)**2 &
 &                             + Flw(m)%qh(i0,j0,k0,3)**2 &
 &                             + Flw(m)%qh(i0,j0,k0,4)**2 ) / Flw(m)%qh(i0,j0,k0,1)
 Flw(m)%qh(i0,j0,k0,6) = 0.0
 Flw(m)%qh(i0,j0,k0,7) = Flw(m)%qh(i1,j1,k1,7) * Flw(m)%jac(i1,j1,k1) / Flw(m)%jac(i0,j0,k0)
 ! 処理終了 ********************************************************************************************
 return
end subroutine WallNoSlipLRe
!*******************************************************************************************************
!******** 周期境界                 				                                 *******
!*******************************************************************************************************
subroutine BoundaryPeriodic( &
&            m, is, ie, js, je, k0, k1 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 integer, intent(in) :: is, ie, js, je
 integer, intent(in) :: k0, k1
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, l
 ! 処理開始 ********************************************************************************************
 do l = ls, le
 do j = js, je
 do i = is, ie
  Flw(m)%qh(i,j,k0,l) = Flw(m)%qh(i,j,k1,l) * Flw(m)%jac(i,j,k1) / Flw(m)%jac(i,j,k0)
 enddo
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryPeriodic
!*******************************************************************************************************
!******** 重合格子の補間									********
!*******************************************************************************************************
subroutine InterpolationOversetGrid
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k, l, m, n
 integer :: i0, j0, k0, &
 &          i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4, &
 &          i5, j5, k5, i6, j6, k6, i7, j7, k7, i8, j8, k8
 ! 処理開始 ********************************************************************************************
 do m = ms, me
  select case(m)
   case(1); n = 2
   case(2); n = 1
  end select
  !$omp parallel do default(shared) &
  !$                private(i, j, k, l, i0, j0, k0, &
  !$                        i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4, &
  !$                        i5, j5, k5, i6, j6, k6, i7, j7, k7, i8, j8, k8)
  do k0 = OSG(m)%ks, OSG(m)%ke
  do j0 = OSG(m)%js, OSG(m)%je
  do i0 = OSG(m)%is, OSG(m)%ie
   i = OSG(m)%ip(i0,j0,k0); j = OSG(m)%jp(i0,j0,k0); k = OSG(m)%kp(i0,j0,k0)
   select case( OSG(m)%fOver(i0,j0,k0) )
    ! 補間しない点 -------------------------------------------------------------------------------------
    case(0)
     cycle
    ! 三重線形補間点 -----------------------------------------------------------------------------------
    case(1)
     i1 = i    ; j1 = j    ; k1 = k    ; i2 = i + 1; j2 = j    ; k2 = k    
     i3 = i    ; j3 = j + 1; k3 = k    ; i4 = i + 1; j4 = j + 1; k4 = k    
     i5 = i    ; j5 = j    ; k5 = k + 1; i6 = i + 1; j6 = j    ; k6 = k + 1
     i7 = i    ; j7 = j + 1; k7 = k + 1; i8 = i + 1; j8 = j + 1; k8 = k + 1
     do l = ls, le
      Flw(m)%qh(i0,j0,k0,l) = ( OSG(m)%term1(i0,j0,k0) * Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) &
      &                       + OSG(m)%term2(i0,j0,k0) * Flw(n)%qh(i2,j2,k2,l) * Flw(n)%jac(i2,j2,k2) &
      &                       + OSG(m)%term3(i0,j0,k0) * Flw(n)%qh(i3,j3,k3,l) * Flw(n)%jac(i3,j3,k3) &
      &                       + OSG(m)%term4(i0,j0,k0) * Flw(n)%qh(i4,j4,k4,l) * Flw(n)%jac(i4,j4,k4) &
      &                       + OSG(m)%term5(i0,j0,k0) * Flw(n)%qh(i5,j5,k5,l) * Flw(n)%jac(i5,j5,k5) &
      &                       + OSG(m)%term6(i0,j0,k0) * Flw(n)%qh(i6,j6,k6,l) * Flw(n)%jac(i6,j6,k6) &
      &                       + OSG(m)%term7(i0,j0,k0) * Flw(n)%qh(i7,j7,k7,l) * Flw(n)%jac(i7,j7,k7) &
      &                       + OSG(m)%term8(i0,j0,k0) * Flw(n)%qh(i8,j8,k8,l) * Flw(n)%jac(i8,j8,k8) &
      &                         ) / Flw(m)%jac(i0,j0,k0)
     enddo
    ! 三次元線形補間点 ---------------------------------------------------------------------------------
    case(2)
     i1 = i    ; j1 = j    ; k1 = k    ; i2 = i    ; j2 = j    ; k2 = k + 1
     i3 = i + 1; j3 = j    ; k3 = k + 1; i4 = i    ; j4 = j + 1; k4 = k + 1
     do l = ls, le
      Flw(m)%qh(i0,j0,k0,l) = ( OSG(m)%term1(i0,j0,k0) * Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) &
      &                       + OSG(m)%term2(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i2,j2,k2,l) * Flw(n)%jac(i2,j2,k2) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                       + OSG(m)%term3(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i3,j3,k3,l) * Flw(n)%jac(i3,j3,k3) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                       + OSG(m)%term4(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i4,j4,k4,l) * Flw(n)%jac(i4,j4,k4) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                         ) / Flw(m)%jac(i0,j0,k0)
     enddo
    case(3)
     i1 = i    ; j1 = j    ; k1 = k    ; i2 = i    ; j2 = j + 1; k2 = k    
     i3 = i + 1; j3 = j + 1; k3 = k    ; i4 = i    ; j4 = j + 1; k4 = k + 1
     do l = ls, le
      Flw(m)%qh(i0,j0,k0,l) = ( OSG(m)%term1(i0,j0,k0) * Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) &
      &                       + OSG(m)%term2(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i2,j2,k2,l) * Flw(n)%jac(i2,j2,k2) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                       + OSG(m)%term3(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i3,j3,k3,l) * Flw(n)%jac(i3,j3,k3) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                       + OSG(m)%term4(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i4,j4,k4,l) * Flw(n)%jac(i4,j4,k4) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                         ) / Flw(m)%jac(i0,j0,k0)
     enddo
    case(4)
     i1 = i    ; j1 = j    ; k1 = k    ; i2 = i + 1; j2 = j + 1; k2 = k    
     i3 = i + 1; j3 = j    ; k3 = k + 1; i4 = i    ; j4 = j + 1; k4 = k + 1
     do l = ls, le
      Flw(m)%qh(i0,j0,k0,l) = ( OSG(m)%term1(i0,j0,k0) * Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) &
      &                       + OSG(m)%term2(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i2,j2,k2,l) * Flw(n)%jac(i2,j2,k2) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                       + OSG(m)%term3(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i3,j3,k3,l) * Flw(n)%jac(i3,j3,k3) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                       + OSG(m)%term4(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i4,j4,k4,l) * Flw(n)%jac(i4,j4,k4) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                         ) / Flw(m)%jac(i0,j0,k0)
     enddo
    case(5)
     i1 = i    ; j1 = j    ; k1 = k    ; i2 = i + 1; j2 = j    ; k2 = k    
     i3 = i + 1; j3 = j + 1; k3 = k    ; i4 = i + 1; j4 = j    ; k4 = k + 1
     do l = ls, le
      Flw(m)%qh(i0,j0,k0,l) = ( OSG(m)%term1(i0,j0,k0) * Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) &
      &                       + OSG(m)%term2(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i2,j2,k2,l) * Flw(n)%jac(i2,j2,k2) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                       + OSG(m)%term3(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i3,j3,k3,l) * Flw(n)%jac(i3,j3,k3) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                       + OSG(m)%term4(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i4,j4,k4,l) * Flw(n)%jac(i4,j4,k4) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                         ) / Flw(m)%jac(i0,j0,k0)
     enddo
    case(6)
     i1 = i + 1; j1 = j + 1; k1 = k    ; i2 = i + 1; j2 = j    ; k2 = k + 1
     i3 = i    ; j3 = j + 1; k3 = k + 1; i4 = i + 1; j4 = j + 1; k4 = k + 1
     do l = ls, le
      Flw(m)%qh(i0,j0,k0,l) = ( OSG(m)%term1(i0,j0,k0) * Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) &
      &                       + OSG(m)%term2(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i2,j2,k2,l) * Flw(n)%jac(i2,j2,k2) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                       + OSG(m)%term3(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i3,j3,k3,l) * Flw(n)%jac(i3,j3,k3) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                       + OSG(m)%term4(i0,j0,k0) &
      &                         * ( Flw(n)%qh(i4,j4,k4,l) * Flw(n)%jac(i4,j4,k4) &
      &                           - Flw(n)%qh(i1,j1,k1,l) * Flw(n)%jac(i1,j1,k1) ) &
      &                         ) / Flw(m)%jac(i0,j0,k0)
     enddo
   end select
  enddo
  enddo
  enddo
  !$omp end parallel do
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine InterpolationOversetGrid
!*******************************************************************************************************
!******** スムージング										********
!*******************************************************************************************************
subroutine Smoothing
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m, n
 ! 処理開始 ********************************************************************************************
 do n  = 1, nSmooth
!  do m = ms, me
  m = me
   call SmoothingFlux3D( &
   &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
   &      Flw(m)%jac, Flw(m)%qh )
!  enddo
!  call InterpolationOversetGrid
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine Smoothing
!*******************************************************************************************************
!******** ファイル出力										********
!*******************************************************************************************************
subroutine OutputFile
 ! 変数定義 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m
 ! 処理開始 ********************************************************************************************
 ! 計算条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Output_CalSetting( trim(CalSetFile) // strtxt )
 ! メトリックス ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Output_Metrics3D( &
  &      trim(FlwIniDir) // trim(BlkName(m)) // trim(MetFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%jac, &
  &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
  &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
  &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez )
! enddo
 ! 流束関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Output_Flux3D( &
  &      trim(FlwIniDir) // trim(BlkName(m)) // trim(IniFlxFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      ls, le, Flw(m)%qh )
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine OutputFile
!*******************************************************************************************************
!******** 無次元化										********
!*******************************************************************************************************
subroutine CalNondimensionalized
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: j, k, m
 ! 処理開始 ********************************************************************************************
 ! 無次元化参照値 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! m = 1
 m = me
 aRef   = sqrt( gamma * Rg * maxval(Flw(m)%t(:,:,:)) )
 LRef   = maxval( sqrt( Flw(m)%x(:,:,:)**2 + Flw(m)%y(:,:,:)**2  + Flw(m)%z(:,:,:)**2 ) )
 RhoRef = maxval(Flw(m)%rho(:,:,:))
 ! 無次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  ! 座標系 ---------------------------------------------------------------------------------------------
  call NondimensionalizedCoord3D( &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
  ! 物理量 ---------------------------------------------------------------------------------------------
  call SetPhysics3DKEM( &
  &      Rg, gamma, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
  &      Flw(m)%qh, Flw(m)%jac, &
  &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps )
  call ViscosityCoefficient3D( &
  &      muSth, TsSth, s1, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%t, Flw(m)%mu )
  call EddyViscosityCoefficient3D( &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%rho, Flw(m)%kin, Flw(m)%eps, cmu, Flw(m)%mut )
  call NondimensionalizedPhysics3D( &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%mu, &
  &      Flw(m)%kin, Flw(m)%eps, Flw(m)%mut )
  ! メトリックス ---------------------------------------------------------------------------------------
  call Metrics3D( &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z, &
  &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
  &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
  &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
  &      Flw(m)%jac )
  ! 流束関数 -------------------------------------------------------------------------------------------
  call SetFlux3DKEM( &
  &      gamma, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
  &      Flw(m)%qh, Flw(m)%jac, &
  &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%kin, Flw(m)%eps )
! enddo
 ! サザーランドの式の定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 s1    = s1 / aRef**2
 muSth = muSth / (RhoRef * aRef * lRef)
 TsSth = TsSth / aRef**2
 ! 流入条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 VelIn  = VelIn / aRef
 PtIn   = PtIn / (rhoRef * aRef**2)
 TtIn   = TtIn / aRef**2
 PsOut  = PsOut / (rhoRef * aRef**2)
 TsOut  = TsOut / aRef**2
 LWC    = LWC  / RhoRef
 MVD    = MVD  / LRef
 RhoD   = RhoD / RhoRef
 SigD   = SigD / (rhoRef * aRef**2 * lRef)
 muD    = muD / (RhoRef * aRef * lRef)
 Span   = Span / LRef
 Chord  = Chord / LRef
 VelExp = VelExp / aRef
 PsExp  = PsExp / (rhoRef * aRef**2)
 TsExp  = TsExp / aRef**2
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalNondimensionalized
!*******************************************************************************************************
!******** 無次元化ファイル出力									********
!*******************************************************************************************************
subroutine OutputNDFile
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m
 ! 処理開始 ********************************************************************************************
 ! 計算条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Output_CalSetting( trim(ND_CalSetFile) // strtxt )
 ! 格子座標 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Output_Grid3D( &
  &      trim(FlwIniDir) // trim(BlkName(m)) // trim(ND_GrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
! enddo
 ! メトリックス ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Output_Metrics3D( &
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
  call Output_Flux3D( &
  &      trim(FlwIniDir) // trim(BlkName(m)) // trim(ND_IniFlxFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      ls, le, Flw(m)%qh )
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine OutputNDFile
!*******************************************************************************************************
!******** 可視化										********
!*******************************************************************************************************
subroutine VisualizedFlow
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: mach(:, :, :), Pt(:, :, :), Tt(:, :, :)
 integer :: m
 ! 処理開始 ********************************************************************************************
! do m = ms, me
 m = me
  ! 物理量に変換
  call SetPhysics3DKEM( &
  &      Rg, gamma, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
  &      Flw(m)%qh, Flw(m)%jac, &
  &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps )
  ! 粘性係数
  call ViscosityCoefficient3D( &
  &      muSth, TsSth, s1, Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%t, Flw(m)%mu )
  ! 渦粘性係数
  call EddyViscosityCoefficient3D( &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%rho, Flw(m)%kin, Flw(m)%eps, cmu, Flw(m)%mut )
  ! マッハ数，全温，全圧
  allocate( mach(Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Pt  (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Tt  (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke) )
  call CalMachTtPt3D( &
  &      Rg, gamma, Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, &
  &      mach, Pt, Tt )
  ! MicroAVSファイル
  call MakeMAVSFile3D( &
  &      trim(FlwIniDir), trim(BlkName(m)) // trim(ViewFlowFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      rhoRef, aRef, lRef, &
  &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%mu, &
  &      Flw(m)%kin, Flw(m)%eps, Flw(m)%mut, mach, Pt, Tt, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
  deallocate( mach, Pt, Tt )
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine VisualizedFlow
! 定義終了 *********************************************************************************************
end program InitialFlow_NACA
