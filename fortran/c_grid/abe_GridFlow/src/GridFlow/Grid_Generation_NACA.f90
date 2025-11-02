!*******************************************************************************************************
!*******************************************************************************************************
!******** 格子生成プログラム									********
!******** (ACA0012翼，三次元，C-type，重合格子法)						********
!******** 					　　　2014.04.10  PROGRAMMED BY RYOSUKE HAYASHI ********
!********					      2014.04.15     UPDATED BY RYOSUKE HAYASHI ********
!******** 可視化用のDATデータ出力		      2016.12.15     UPDATED BY SHO     URANAI  ********
!*******************************************************************************************************
!*******************************************************************************************************
program GridGeneration_NACA
 ! モジュール宣言 **************************************************************************************
 use Package_NACA
 use Package_FileIO
 use Package_Equation
 use Package_Grid
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! ファイル名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: ViewGrdFile * 8 = 'ViewGrid'
 ! 共有定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, parameter :: NACA1  =  0					! NACA 翼の反り比
 integer, parameter :: NACA2  =  0					! NACA 翼の最大反り位置
 integer, parameter :: NACA34 =  12					! NACA 翼の最大翼厚
 real   , parameter :: a0     =  0.2969					! 翼厚分布式の定数
 real   , parameter :: a1     = -0.1260					! 翼厚分布式の定数
 real   , parameter :: a2     = -0.3516					! 翼厚分布式の定数
 real   , parameter :: a3     =  0.2843					! 翼厚分布式の定数
 real   , parameter :: a4     = -0.1015					! 翼厚分布式の定数
 ! 共有変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   :: thick, ytmax, ycmax, xc
 ! 処理開始 ********************************************************************************************
 ! 検証実験ケース選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(a)') "<< Computational Condition >>"
 call SelectExpCase
 ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(a)') '<< Initial Setting >>'
 call InitialSetting
 ! Main Grid +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call MainGrid
 write(*, '(a)') '+++++ Main grid complete. +++++'
 ! Sub Grid +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call SubGrid
 write(*, '(a)') '+++++ Sub grid complete. +++++'
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
 ! 実験条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 実験ケース入力
 write(*, '(a)') '* Select Exp. Case...'
 read(*, *) ExpCaseNum
 write(fname, '(i2.2)') ExpCaseNum
 ! 実験条件入力
 call Input_ExpCondition( &
 &      './data/GridFlow/' // trim(ExpConditionFile) // trim(fname) // strtxt )
 ! 着氷モデル選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(a)') '* Select Icing Model...'
 read(*, *) ThermoNum
 write(fname, '(i2.2)') ThermoNum  !これって意味あるの?
 ! ディレクトリ設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(fname, '(a)') 'clean/'
 GrdOutDir = resultdir // 'grid/' // trim(fname)
 ! 計算条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 計算条件設定
 nCount = 0; nDrop = 0; IceStep = 0; IceStepMax = 4
 IceTime = 0.0; dti = 1.0e-3
 Cn = 1.0; nRunge = 4; nCalMax = 20000; nOutputLog = 1000; nOutputCount = 1000
 nTVD = 2; eTVD = 0.1; LmtPro = 0.1; LmtAve = 100.0
 TurbNum = 2; fSteady = .true.; fTime = .false.
 BCNum = 2; fSlip = .false.
 RoughNum = 1; RunbackNum = 1; DropShedNum = 0; IceShedNum = 0; DragNum = 1
 Cnd = 0.03
 nDrpIn = 1000000; nDrpCalMax = 100000; nDrpFile = 100; nDrpOutputLog = 1000; nDrpOutputCount = 50000
 Rg = 287.1; gamma = 1.4; Pr = 0.72; Prt = 0.90; Ret = 500.0; Cmu = 0.09
 muSth = 1.82e-5; TsSth = 273.15; s1 = 117.0
 ! 計算条件ファイル出力
 call Output_CalSetting( './data/GridFlow/'//trim(ND_CalSetFile) // strtxt )
 ! 処理終了 ********************************************************************************************
 return
end subroutine SelectExpCase
!*******************************************************************************************************
!********* 初期設定										********
!*******************************************************************************************************
subroutine InitialSetting
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 処理開始 ********************************************************************************************
 ! ブロック名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Set_BlockName
 ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( Grd(ms:me) )
 ! 翼形状 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ycmax = chord * real(NACA1)  * 1.0e-2
 xc    = chord * real(NACA2)  * 1.0e-1
 thick = chord * real(NACA34) * 1.0e-2
 ytmax = 0.5 * thick
 ! 処理終了 ********************************************************************************************
 return
end subroutine InitialSetting
!*******************************************************************************************************
!********* Main Grid 生成									********
!*******************************************************************************************************
subroutine MainGrid
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m
 real    :: dom1, dom2, dom3
 ! 処理開始 ********************************************************************************************
 ! ブロック番号 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 m = 1
 ! 解像度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Grd(m)%is =  0
 Grd(m)%i1 = 50 + Grd(m)%is !50+ !kouen
 Grd(m)%i2 = 80 + Grd(m)%i1 !60+ !C
 Grd(m)%i3 = 80 + Grd(m)%i2 !60+ !C
 Grd(m)%ie = 50 + Grd(m)%i3 !50+ !kouen
 Grd(m)%js =  0
 Grd(m)%je = 70
 Grd(m)%ks =  0
 Grd(m)%ke =  3
 ! 計算領域 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 dom1 = 15.0 * chord	!10*						! 後縁から流出境界まで
 dom2 = 15.0 * chord	!10*						! 後縁から外部境界まで
 dom3 = 15.0 * chord	!10*						! 最大翼厚位置から流入境界まで
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( Grd(m)%x( Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke ), &
 &         Grd(m)%y( Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke ), &
 &         Grd(m)%z( Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke ), &
 &         Grd(m)%f( Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke ) )
 ! 格子生成 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(a)') '+++++ C-typr grid generation start. +++++'
 call CtypeGridBlade( &
 &      m, Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
 &      Grd(m)%i1, Grd(m)%i2, Grd(m)%i3, dom1, dom2, dom3, &
 &      Grd(m)%x , Grd(m)%y , Grd(m)%z )
 ! 翼のフラグ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call ViewBlade( &
 &      m, Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
 &      Grd(m)%i1, Grd(m)%i3, Grd(m)%js, Grd(m)%f )
 ! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Output_Resolution3D( &
 &      trim(GrdOutDir) // trim(BlkName(m)) // trim(RslFile), strtxt, &
 &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke )
 call Output_Resolution1D( &
 &      trim(GrdOutDir) // trim(BlkName(m)) // trim(IceRslFile), strtxt, &
 &      Grd(m)%i1, Grd(m)%i3 )
 call Output_CtypeGridPoint( &
 &      trim(GrdOutDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
 &      Grd(m)%i1, Grd(m)%i2, Grd(m)%i3 )
 call Output_Grid3D( &
 &      trim(GrdOutDir) // trim(BlkName(m)) // trim(GrdFile), strbin, &
 &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
 &      Grd(m)%x, Grd(m)%y, Grd(m)%z )
! call Output_Grid3D( &
! &      trim(GrdOutDir) // trim(BlkName(m)) // trim(GrdFile), strdat, &
! &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
! &      Grd(m)%x, Grd(m)%y, Grd(m)%z )
 ! 可視化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call MakeMAVSFile3D( &
 &      trim(GrdOutDir), trim(BlkName(m)) // trim(ViewGrdFile), strbin, &
 &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
 &      Grd(m)%f , Grd(m)%x , Grd(m)%y , Grd(m)%z )
 call OutputPara_bin( &
    &      trim(GrdOutDir), trim(BlkName(m)) // trim(ViewGrdFile), &
    &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
    &      Grd(m)%x, Grd(m)%y, Grd(m)%z )

 ! 処理終了 ********************************************************************************************
 return
end subroutine MainGrid
!*******************************************************************************************************
!********* Sub Grid 生成									********
!*******************************************************************************************************
subroutine SubGrid
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m
 real    :: dom1, dom2
 real    :: dom3, dom4
 ! 処理開始 ********************************************************************************************
 ! ブロック番号 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 m = me
 ! 解像度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 select case(swi_subgrid)
  case(1) ! H_type
   Grd(m)%is =   0
   Grd(m)%ie = 500 !250!200
   Grd(m)%js =   0
   Grd(m)%je =  50
   Grd(m)%ks =   0
   Grd(m)%ke =   3
  case(2) ! C_type
   Grd(m)%is =  0
   Grd(m)%i1 = 50 + Grd(m)%is !50+ !kouen
   Grd(m)%i2 = 100 + Grd(m)%i1 !60+ !C
   Grd(m)%i3 = 100 + Grd(m)%i2 !60+ !C
   Grd(m)%ie = 50 + Grd(m)%i3 !50+ !kouen
   Grd(m)%js =  0
   Grd(m)%je = 50
   Grd(m)%ks =  0
   Grd(m)%ke =  3
  case(3) ! O_type
   Grd(m)%is =   0
   Grd(m)%i1 = 150 + Grd(m)%is !110
   Grd(m)%i2 = 100 + Grd(m)%i1 !50
   Grd(m)%i3 = 100 + Grd(m)%i2 !50
   Grd(m)%ie = 150 + Grd(m)%i3 !110	! 300
   Grd(m)%js =   0
   Grd(m)%j1 =  30 + Grd(m)%js !30
   Grd(m)%je =  10 + Grd(m)%j1 !100	! 30
   Grd(m)%ks =  0
   Grd(m)%ke =  3
 end select

 ! 計算領域 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 select case(swi_subgrid)
  case(1) ! H_type
   dom1 = 0.8 * chord			! 前縁から着氷限界位置付近まで
   dom2 = 0.3 * chord			! 最大氷層厚さ以上
  case(2) ! C_type
   dom1 = 1.3 * chord	!10*		! 後縁から流出境界まで
   dom2 = 0.3 * chord	!10*		! 後縁から外部境界まで
   dom3 = 0.3 * chord	!10*		! 最大翼厚位置から流入境界まで
  case(3) ! O_type
   dom1 = 0.3 * chord 			! 翼法線方向計算領域
   dom2 = 0.4 * chord 			! 前縁から着氷するとこまで (i)内
   dom3 = 1.0 * chord 			! 前縁から着氷するとこまで (i)外
   dom4 = 0.1 * chord 			! 前縁から着氷するとこまで (j)外
 end select

 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( Grd(m)%x( Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke ), &
 &         Grd(m)%y( Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke ), &
 &         Grd(m)%z( Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke ), &
 &         Grd(m)%f( Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je, Grd(m)%ks: Grd(m)%ke ) )
 ! 格子生成 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 select case(swi_subgrid)
  case(1) ! H_type
   call HtypeGridIceLE( &
   &      m, Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, dom1, dom2, &
   &      Grd(m)%x , Grd(m)%y , Grd(m)%z )
  case(2) ! C_type
   call CtypeGridBlade_sub( &
   &      m, Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
   &      Grd(m)%i1, Grd(m)%i2, Grd(m)%i3, dom1, dom2, dom3, &
   &      Grd(m)%x , Grd(m)%y , Grd(m)%z )
  case(3) ! O_type
   call OtypeStatorVaneGrid( &
   &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
   &      Grd(m)%i1, Grd(m)%i2, Grd(m)%i3, Grd(m)%j1, dom1, dom2, dom3, dom4, &
   &      Grd(m)%x , Grd(m)%y , Grd(m)%z )
 end select

 ! 翼のフラグ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call ViewBlade( &
 &      m, Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
 &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%f )
 ! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Output_Resolution3D( &
 &      trim(GrdOutDir) // trim(BlkName(m)) // trim(RslFile), strtxt, &
 &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke )
 call Output_Resolution1D( &
 &      trim(GrdOutDir) // trim(BlkName(m)) // trim(IceRslFile), strtxt, &
 &      Grd(m)%is, Grd(m)%ie )
 call Output_Grid3D( &
 &      trim(GrdOutDir) // trim(BlkName(m)) // trim(GrdFile), strbin, &
 &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, Grd(m)%x, Grd(m)%y, Grd(m)%z )
 call Output_Grid3D( & !全体用
 &      trim(GrdOutDir) // trim(BlkName(m)) // trim(GrdFile), strdat, &
 &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, Grd(m)%x, Grd(m)%y, Grd(m)%z )
 call Output_Grid3D( & !翼形状のみ
 &      trim(GrdOutDir) // trim(BlkName(m)) // trim(GrdFile)//'_clean', strdat, &
 &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%js, Grd(m)%ks, Grd(m)%ks, Grd(m)%x, Grd(m)%y, Grd(m)%z )

 if(swi_subgrid .eq. 2) then ! C_type
  call Output_CtypeGridPoint( &
  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
  &      Grd(m)%i1, Grd(m)%i2, Grd(m)%i3 )
 end if

 ! 可視化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call MakeMAVSFile3D( &
 &      trim(GrdOutDir), trim(BlkName(m)) // trim(ViewGrdFile), strbin, &
 &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
 &      Grd(m)%f , Grd(m)%x , Grd(m)%y , Grd(m)%z )
 call OutputPara_bin( &
    &      trim(GrdOutDir), trim(BlkName(m)) // trim(ViewGrdFile), &
    &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
    &      Grd(m)%x, Grd(m)%y, Grd(m)%z )
 ! 処理終了 ********************************************************************************************
 return
end subroutine SubGrid
!*******************************************************************************************************
!********* 翼周りの格子 (C-type)								********
!*******************************************************************************************************
subroutine CtypeGridBlade( &
&            m, is, ie, js, je, ks, ke, i1, i2, i3, dom1, dom2, dom3, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: m
 integer, intent(in)  :: is, ie, js, je, ks, ke
 integer, intent(in)  :: i1, i2, i3
 real   , intent(in)  :: dom1, dom2, dom3
 real   , intent(out) :: x(is:ie, js:je, ks:ke), &
 &                       y(is:ie, js:je, ks:ke), &
 &                       z(is:ie, js:je, ks:ke)
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: rs1 = 0.0 * 1.0e-1				! 領域ⅰ起点格子幅
 real   , parameter :: re1 = 6.0 * 1.0e-2				! 領域ⅰ終点格子幅
 real   , parameter :: rs2 = 7.0 * 1.0e-1				! 領域ⅱ起点格子幅
 real   , parameter :: re2 = 1.0 * 1.0e-1                               ! 領域ⅱ終点格子幅
 real   , parameter :: rs3 = 1.0 * 1.0e-1				! 領域ⅲ起点格子幅
 real   , parameter :: re3 = 7.0 * 1.0e-1                               ! 領域ⅲ終点格子幅
 real   , parameter :: rs4 = 6.0 * 1.0e-2				! 領域ⅳ起点格子幅
 real   , parameter :: re4 = 0.0 * 1.0e0                                ! 領域ⅳ終点格子幅
 real   , parameter :: rs5 = 3.0 * 1.0e-1				! 内部境界格子幅 !j-第一格子点
 real   , parameter :: re5 = 0.0 * 1.0e0				! 外部境界格子幅
! real   , parameter :: rs5 = 1.0 * 1.0e-3				! 内部境界格子幅
! real   , parameter :: re5 = 0.0 * 1.0e-2				! 外部境界格子幅
 real   , parameter :: tb1 = 1.0 * 3.0e-1				! 内部境界直交性のパラメータ
 real   , parameter :: tb2 = 1.0 * 1.0e+1				! 外部境界直交性のパラメータ
 real   , parameter :: mgn = 3.0					! 楕円-双曲型重み関数の許容誤差
 real   , parameter :: rsd = 1.0e-5					! 楕円-双曲型最大計算回数
 ! 処理開始 ********************************************************************************************
 ! 内部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call CtypeInternalBoundary( &
 &      is, i1, i2, i3, ie, dom1, rs1, re1, rs2, re2, rs3, re3, rs4 ,re4, x(:, js, ks), y(:, js, ks) )
 ! 外部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call CtypeExternalBoundary( &
 &      is, i1, i2, i3, ie, js, je, dom2, dom3, x(:, :, ks), y(:, :, ks) )
! ! 側部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call CtypeSideBoundary( &
! &      is, ie, js, je, x(:, :, ks), y(:, :, ks) )
! ! Transfinite 補間 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call TransfiniteInterpolation( &
! &      is, ie, js, je, x(:, :, ks), y(:, :, ks) )
! ! 楕円-双曲型偏微分法 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call CtypeGenerationEHPDE( &
! &      is, i1, i2, i3, ie, js, je, rs5, re5, mgn, rsd, x(:, :, ks), y(:, :, ks) )
 ! 二境界法 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call GenerationTwoBoundary( &
 &      m, is, ie, i1, i3, js, je, rs5, re5, tb1, tb2, x(:, :, ks), y(:, :, ks) )
 ! 三次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call ThreeDimensionalized( &
 &      is, ie, js, je, ks, ke, span, x, y, z )
 ! 処理終了 ********************************************************************************************
 return
end subroutine CtypeGridBlade

subroutine CtypeGridBlade_sub( &
&            m, is, ie, js, je, ks, ke, i1, i2, i3, dom1, dom2, dom3, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: m
 integer, intent(in)  :: is, ie, js, je, ks, ke
 integer, intent(in)  :: i1, i2, i3
 real   , intent(in)  :: dom1, dom2, dom3
 real   , intent(out) :: x(is:ie, js:je, ks:ke), &
 &                       y(is:ie, js:je, ks:ke), &
 &                       z(is:ie, js:je, ks:ke)
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: rs1 = 0.0 * 1.0e-1				! 領域ⅰ起点格子幅
 real   , parameter :: re1 = 2.0e-1 !1.0e-3 !6.0 * 1.0e-2				! 領域ⅰ終点格子幅
 real   , parameter :: rs2 = 7.0 * 1.0e-1				! 領域ⅱ起点格子幅
 real   , parameter :: re2 = 1.0 * 1.0e-1                               ! 領域ⅱ終点格子幅
 real   , parameter :: rs3 = 1.0 * 1.0e-1				! 領域ⅲ起点格子幅
 real   , parameter :: re3 = 7.0 * 1.0e-1                               ! 領域ⅲ終点格子幅
 real   , parameter :: rs4 = 2.0e-1 !1.0e-3 !6.0 * 1.0e-2				! 領域ⅳ起点格子幅
 real   , parameter :: re4 = 0.0 * 1.0e0                                ! 領域ⅳ終点格子幅
 real   , parameter :: rs5 = 3.0e-1 !6.0e-1 !5.0e-2 !1.0e-0 !3.0 * 1.0e-1			! 内部境界格子幅 !j-第一格子点
 real   , parameter :: re5 = 0.0 * 1.0e0				! 外部境界格子幅
! real   , parameter :: rs5 = 1.0 * 1.0e-3				! 内部境界格子幅
! real   , parameter :: re5 = 0.0 * 1.0e-2				! 外部境界格子幅
 real   , parameter :: tb1 = 5.0e-2 !1.0e-1 !3.0e-2!1.0 * 3.0e-1		! 内部境界直交性のパラメータ
 real   , parameter :: tb2 = 1.0e-2 !2.0e-1 !1.5e-1!1.0 * 1.0e+1				! 外部境界直交性のパラメータ
 real   , parameter :: mgn = 3.0					! 楕円-双曲型重み関数の許容誤差
 real   , parameter :: rsd = 1.0e-5					! 楕円-双曲型最大計算回数
 ! 処理開始 ********************************************************************************************
 ! 内部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call CtypeInternalBoundary( &
 &      is, i1, i2, i3, ie, dom1, rs1, re1, rs2, re2, rs3, re3, rs4 ,re4, x(:, js, ks), y(:, js, ks) )
 ! 外部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call CtypeExternalBoundary( &
 &      is, i1, i2, i3, ie, js, je, dom2, dom3, x(:, :, ks), y(:, :, ks) )
! ! 側部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call CtypeSideBoundary( &
! &      is, ie, js, je, x(:, :, ks), y(:, :, ks) )
! ! Transfinite 補間 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call TransfiniteInterpolation( &
! &      is, ie, js, je, x(:, :, ks), y(:, :, ks) )
! ! 楕円-双曲型偏微分法 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call CtypeGenerationEHPDE( &
! &      is, i1, i2, i3, ie, js, je, rs5, re5, mgn, rsd, x(:, :, ks), y(:, :, ks) )
 ! 二境界法 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call GenerationTwoBoundary( &
 &      m, is, ie, i1, i3, js, je, rs5, re5, tb1, tb2, x(:, :, ks), y(:, :, ks) )
 ! スムージング
 call Smoothing( &
 &      is, ie, js, je, x(:, :, ks), y(:, :, ks) )
 ! 三次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call ThreeDimensionalized( &
 &      is, ie, js, je, ks, ke, span, x, y, z )
 ! 処理終了 ********************************************************************************************
 return
end subroutine CtypeGridBlade_sub

!*******************************************************************************************************
!******** 静翼周り生成										********
!*******************************************************************************************************
SUBROUTINE OtypeStatorVaneGrid( &
&            is, ie, js, je, ks, ke, i1, i2, i3, j1, dom1, dom2, dom3,dom4, x, y, z )
! 変数型宣言 ********************************************************************************************
 implicit none
! 変数宣言 **********************************************************************************************
 integer, intent(in)  :: is, ie, js, je, ks, ke
 integer, intent(in)  :: i1, i2, i3
 integer, intent(in)  :: j1
 real   , intent(in)  :: dom1, dom2, dom3, dom4
 real   , intent(out) :: x(is:ie, js:je, ks:ke), &
 &                       y(is:ie, js:je, ks:ke), &
 &                       z(is:ie, js:je, ks:ke)
  ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   real   , parameter	:: rs1  = 1.2e-1, re1 = 2.6e-1			! 翼周方向S.S.格子幅（内側境界）
!    real   , parameter	:: rs2  = 5.5e-1, re2 = 5.5e-1
!    real   , parameter	:: rs3  = 2.6e-1, re3 = 1.2e-1			! 翼周方向P.S.格子幅（内側境界）
!    real   , parameter	:: rs4  = 4.2e-1, re4 = 4.8e-1			! 翼周方向S.S.格子幅（外側境界）
!    real   , parameter	:: rs5  = 4.8e-1, re5 = 4.2e-1			! 翼周方向P.S.格子幅（外側境界）
!    real   , parameter	:: rs6  = 9.0e-1, re6 = 0.0e+0			! 翼法線方向格子幅（側部）
!    real   , parameter	:: rs7  = 5.0e-2, re7 = 7.0e-2
 real   , parameter	:: rs1  = 5.0e-1, re1 = 4.0e-1			! 翼周方向S.S.格子幅（内側境界）
 real   , parameter	:: rs2  = 3.0e-1, re2 = 3.0e-1
 real   , parameter	:: rs3  = 4.0e-1, re3 = 5.0e-1			! 翼周方向P.S.格子幅（内側境界）
 real   , parameter	:: rs4  = 3.0e-0, re4 = 1.0e-0			! 翼周方向S.S.格子幅（外側境界）
 real   , parameter	:: rs5  = 1.0e-0, re5 = 3.0e-0			! 翼周方向P.S.格子幅（外側境界）
 real   , parameter	:: rs6  = 4.0e-2, re6 = 0.0e-0			! 翼法線方向格子幅（側部）
 real   , parameter	:: rs7  = 3.0e-2, re7 = 1.0e-1			! 翼法線方向格子幅（楕円双曲）
 real   , parameter	:: MGN  = 5.0 !1.0e-3				! 楕円-双曲型偏微分法許容誤差
 real   , parameter	:: Res  = 1.0e-5				! 楕円-双曲型偏微分法許容誤差

! real   , parameter	:: rs1  = 1.2e-1, re1 = 2.6e-1			! 翼周方向S.S.格子幅（内側境界）
! real   , parameter	:: rs2  = 5.5e-1, re2 = 5.5e-1
! real   , parameter	:: rs3  = 2.6e-1, re3 = 1.2e-1			! 翼周方向P.S.格子幅（内側境界）
! real   , parameter	:: rs4  = 4.2e-1, re4 = 4.8e-1			! 翼周方向S.S.格子幅（外側境界）
! real   , parameter	:: rs5  = 4.8e-1, re5 = 4.2e-1			! 翼周方向P.S.格子幅（外側境界）
! real   , parameter	:: rs6  = 9.0e-1, re6 = 0.0e+0			! 翼法線方向格子幅（側部）
! real   , parameter	:: rs7  = 5.0e-2, re7 = 7.0e-2			! 翼法線方向格子幅（楕円双曲）
! real   , parameter	:: MGN  = 2.0					! 楕円-双曲型偏微分法許容誤差
! real   , parameter	:: Res  = 1.0e-5				! 楕円-双曲型偏微分法許容誤差

  ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer :: jp


integer :: i,j

! 処理開始 ********************************************************************************************
!! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  allocate( Grd(m)%x( Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je ), &
!  &         Grd(m)%y( Grd(m)%is: Grd(m)%ie, Grd(m)%js: Grd(m)%je )  )
! 格子解像度入力 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! iは周方向の格子点数，jは翼と垂直方向の格子点数
!  Grd(m)%is = SVOis
!  Grd(m)%i1 = SVOi1 + Grd(m)%is
!  Grd(m)%i2 = SVOi2 + Grd(m)%i1
!  Grd(m)%i3 = SVOi3 + Grd(m)%i2
!  Grd(m)%ie = SVOie + Grd(m)%i3	! 300
!  Grd(m)%js = SVOjs
!  Grd(m)%j1 = SVOj1 + Grd(m)%js
!  Grd(m)%je = SVOje + Grd(m)%j1	! 30
!! スプライン補間 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  write(*, '(a)') '* Spline interpolation start...'
!  jp = js
!  call SplineInterpolation( &
!  &      nSV, SVX, SVY, &
!  &      is, ie, &
!  &      x(:,jp), y(:,jp) )
! 内部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  write(*, '(a)') '* Inner boundary start...'
  call OtypeInnerBoundary( &
  &      is, ie, i1, i2, i3, &
  &      dom2, rs1, re1, rs2, re2, rs3, re3, &
  &      x(:,js,ks), y(:,js,ks) )
! 外部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  write(*, '(a)') '* Outer boundary start...'
  call OtypeOuterBoundary( &
  &      is, ie, js, je, i1, i2, i3, &
  &      dom1, dom3, rs4, re4, rs5, re5, &
  &      x, y )
! 側部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  write(*, '(a)') '* Side boundary start...'
  call OtypeSideBoundary( &
  &      is, ie, js, je, j1, &
  &      dom4, rs6, re6, &
  &      x, y )
! Transfinite 補間 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  write(*, '(a)') '* Transfinite interpolation start...'
  call TransfiniteInterpolation( &
  &      is, ie, js, je, &
  &      x, y )
! 楕円-双曲型偏微分法 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  write(*, '(a)') '* EHPDE grid generation start...'
  call OtypeEHPDEGeneration( &
  &      is, ie, js, je, rs7, re7, MGN, Res, &
  &      x, y )
 ! 三次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call ThreeDimensionalized( &
 &      is, ie, js, je, ks, ke, span, x, y, z )

!! 回転座標 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  write(*, '(a)') '* Rotation matrix start...'
!  call RotationMatrix( &
!  &      AOA, Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, &
!  &      Grd(m)%x, Grd(m)%y )
!! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  write(*, '(a)') '* Check grid start...'
!  call Output_Resolution2D( &
!  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(RslFile), strtxt, &
!  &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je )
! 処理終了 ********************************************************************************************
  return
END SUBROUTINE OtypeStatorVaneGrid


!*******************************************************************************************************
!******** 内部境界 (C-type)									********
!*******************************************************************************************************
subroutine CtypeInternalBoundary( &
&            is, i1, i2, i3, ie, dom, rs1, re1, rs2, re2, rs3, re3, rs4 ,re4, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, i1, i2, i3, ie
 real   , intent(in)  :: dom
 real   , intent(in)  :: rs1, re1, rs2, re2, rs3, re3, rs4 ,re4
 real   , intent(out) :: x(is:ie), y(is:ie)
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, parameter :: nmax = 10000
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: xu(:), xl(:), yu(:), yl(:)
 real   , pointer :: tu(:), tl(:), lu(:), ll(:)
 real   , pointer :: x2(:), y2(:)
 real   , pointer :: tb(:)
 real   , pointer :: sy(:)
 integer :: n, i
 integer :: n1, n2, n3, n4, n5
 integer :: si
 real    :: xx, yy, yt, yc, theta
 real    :: rr
 ! 処理開始 ********************************************************************************************
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( xu(0:nmax), yu(0:nmax), xl(0:nmax), yl(0:nmax) )
 allocate( tu(i2:i3), tl(i1:i2), lu(0:nmax), ll(0:nmax) )
 allocate( x2(i1:i2), y2(i1:i2) )
 allocate( tb(is:i1) )
 ! 翼周り ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 翼形状 ----------------------------------------------------------------------------------------------!check
 do n = 1, nmax - 1
  ! コード方向位置
  xx = chord * real(n) / real(nmax)
  ! 翼厚分布
  yt = 5.0 * thick * ( a0 * (xx / chord)**0.5 + a1 * (xx / chord)**1 + a2 * (xx / chord)**2 &
  &                  + a3 * (xx / chord)**3   + a4 * (xx / chord)**4 )
  ! 翼反り分布
  if(xc == 0.0) then
    yc = 0.0
   else
    if(xx < xc) then
      yc = ycmax * chord**2 / xc**2 * (2.0 * xc * xx / chord**2 - xx**2 / chord**2)
     else
      yc = ycmax * chord**2 / (chord - xc)**2 &
      &  * (1.0 - 2.0 * xc / chord + 2.0 * xc * xx / chord**2 - xx**2 / chord**2)
    endif
  endif
  ! 翼座標
  theta = atan(yc / xx)
  xu(n) = xx - yt * sin(theta)
  yu(n) = yc + yt * cos(theta)
  xl(n) = xx + yt * sin(theta)
  yl(n) = yc - yt * cos(theta)
 enddo
 xu(0) = 0.0; xu(nmax) = chord
 yu(0) = 0.0; yu(nmax) = 0.0
 xl(0) = 0.0; xl(nmax) = chord
 yl(0) = 0.0; yl(nmax) = 0.0
 ! 格子幅調整 ------------------------------------------------------------------------------------------
 ! 媒介変数
 call VinokurInterpolation( rs3 / real(i3-i2+1), re3 / real(i3-i2+1), i3-i2+1, tu )
 call VinokurInterpolation( re2 / real(i2-i1+1), rs2 / real(i2-i1+1), i2-i1+1, tl )
 ! 前縁からの距離
 lu(0) = 0.0; ll(0) = 0.0
 do n = 1, nmax
  lu(n) = sqrt( (xu(n) - xu(n-1))**2 + (yu(n) - yu(n-1))**2 ) + lu(n-1)
  ll(n) = sqrt( (xl(n) - xl(n-1))**2 + (yl(n) - yl(n-1))**2 ) + ll(n-1)
 enddo
 ! 無次元化
 do n = 0, nmax
  lu(n) = lu(n) / lu(nmax)
  ll(n) = ll(n) / ll(nmax)
 enddo
 ! ラグランジュ補間で探索 (上面)
 do i = i2 + 1, i3 - 1
  do n = 0, nmax - 1
   if( lu(n) <= tu(i) .and. tu(i) < lu(n+1) ) then
     select case(n)
      case default
       n1 = n-2; n2 = n-1; n3 = n; n4 = n+1; n5 = n+2
       call LagrangeInterpolation( 5, lu(n1:n5), xu(n1:n5), tu(i), x(i) )
       call LagrangeInterpolation( 5, lu(n1:n5), yu(n1:n5), tu(i), y(i) )
      case(0)
       n1 = n; n2 = n+1; n3 = n+2
       call LagrangeInterpolation( 3, lu(n1:n3), xu(n1:n3), tu(i), x(i) )
       call LagrangeInterpolation( 3, lu(n1:n3), yu(n1:n3), tu(i), y(i) )
      case(1)
       n1 = n-1; n2 = n; n3 = n+1; n4 = n+2
       call LagrangeInterpolation( 4, lu(n1:n4), xu(n1:n4), tu(i), x(i) )
       call LagrangeInterpolation( 4, lu(n1:n4), yu(n1:n4), tu(i), y(i) )
      case(nmax-1)
       n1 = n-2; n2 = n-1; n3 = n; n4 = n+1
       call LagrangeInterpolation( 4, lu(n1:n4), xu(n1:n4), tu(i), x(i) )
       call LagrangeInterpolation( 4, lu(n1:n4), yu(n1:n4), tu(i), y(i) )
     end select
     exit
   endif
  enddo
 enddo
 x(i2) = xu(0); x(i3) = xu(nmax)
 y(i2) = yu(0); y(i3) = yu(nmax)
 ! ラグランジュ補間で探索 (下面)
 do i = i1 + 1, i2 - 1
  do n = 0, nmax - 1
   if( ll(n) <= tl(i) .and. tl(i) < ll(n+1) ) then
     select case(n)
      case default
       n1 = n-2; n2 = n-1; n3 = n; n4 = n+1; n5 = n+2
       call LagrangeInterpolation( 5, ll(n1:n5), xl(n1:n5), tl(i), x(i) )
       call LagrangeInterpolation( 5, ll(n1:n5), yl(n1:n5), tl(i), y(i) )
      case(0)
       n1 = n; n2 = n+1; n3 = n+2
       call LagrangeInterpolation( 3, ll(n1:n3), xl(n1:n3), tl(i), x(i) )
       call LagrangeInterpolation( 3, ll(n1:n3), yl(n1:n3), tl(i), y(i) )
      case(1)
       n1 = n-1; n2 = n; n3 = n+1; n4 = n+2
       call LagrangeInterpolation( 4, ll(n1:n4), xl(n1:n4), tl(i), x(i) )
       call LagrangeInterpolation( 4, ll(n1:n4), yl(n1:n4), tl(i), y(i) )
      case(nmax-1)
       n1 = n-2; n2 = n-1; n3 = n; n4 = n+1
       call LagrangeInterpolation( 4, ll(n1:n4), xl(n1:n4), tl(i), x(i) )
       call LagrangeInterpolation( 4, ll(n1:n4), yl(n1:n4), tl(i), y(i) )
     end select
     exit
   endif
  enddo
 enddo
 x(i1) = xl(0); x(i2) = xl(nmax)
 y(i1) = yl(0); y(i2) = yl(nmax)
 ! 格子番号整理
 do i = i1, i2
  x2(i) = x(i1+i2-i)
  y2(i) = y(i1+i2-i)
 enddo
 do i = i1, i2
  x(i) = x2(i)
  y(i) = y2(i)
 enddo
 ! 後縁をスムージング ----------------------------------------------------------------------------------
 allocate( sy(is:ie) )
 si = 3
 do n = 1, 5
  do i = i1 + 1, i1 + si
   sy(i) = ( y(i-1) + y(i) + y(i+1) ) / 3.0
  enddo
  do i = i1 + 1, i1 + si
   y(i) = sy(i)
  enddo
  do i = i3 - si, i3 - 1
   sy(i) = ( y(i-1) + y(i) + y(i+1) ) / 3.0
  enddo
  do i = i3 - si, i3 - 1
   y(i) = sy(i)
  enddo
 enddo
 deallocate(sy)
 ! C 型格子ブランチ・カット ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call GeometricInterpolationInv( re1 / real(i1-is+1), i1-is+1, tb, rr )
 do i = is, i1 - 1
  x(i) = dom - (dom - x(i1)) * tb(i)
  y(i) = y(i1)
 enddo
 do i = i3, ie
  x(i) = x(ie-i)
  y(i) = y(ie-i)
 enddo
 ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 deallocate(xu, yu, tu, lu, xl, yl, tl, ll, x2, y2, tb)
 ! 処理終了 ********************************************************************************************
 return
end subroutine CtypeInternalBoundary
!*******************************************************************************************************
!******** 外部境界 (C-type)									********
!*******************************************************************************************************
subroutine CtypeExternalBoundary( &
&            is, i1, i2, i3, ie, js, je, dom1, dom2, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, i1, i2, i3, ie, js, je
 real   , intent(in)    :: dom1, dom2
 real   , intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 integer :: ic
 real    :: a1, a2, a3, b1, b2, b3, c1, c2, c3
 ! 処理開始 ********************************************************************************************
 ! C 型の湾曲部 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do i = i1 + 1, i3 - 1
  ! 翼周り方向ベクトル
  a1 = 0.5 * (- x(i-1,js) + x(i+1,js))
  a2 = 0.5 * (- y(i-1,js) + y(i+1,js))
  a3 = 0.0
  ! スパン方向ベクトル
  b1 =  0.0
  b2 =  0.0
  b3 = -1.0
  ! 法線方向ベクトル
  c1 = a2 * b3 - a3 * b2
  c2 = a3 * b1 - a1 * b3
  c3 = a1 * b2 - a2 * b1
  ! 翼周り部外部領域
  x(i,je) = c1 / sqrt(c1**2 + c2**2 + c3**2) * dom1 + x(i,js)
  y(i,je) = c2 / sqrt(c1**2 + c2**2 + c3**2) * dom1 + y(i,js)
 enddo
 ! C 型の直線部 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! x 方向 ----------------------------------------------------------------------------------------------
 ! 下面
 do i = i1 + 1, i2
  if( x(i,je) < x(i1,js) ) exit
 enddo
 ic = i - 1
 do i = is, ic
  x(i,je) = x(ic,je) + (x(is,js) - x(ic,je)) * real(ic-i) / real(ic)
 enddo
 ! 上面
 do i = i2, i3 - 1
  if( x(i,je) > x(i1,js) ) exit
 enddo
 ic = i
 do i = ic, ie
  x(i,je) = x(ic,je) + (x(is,js) - x(ic,je)) * real(i-ic) / real(ie-ic)
 enddo
 ! y 方向 ----------------------------------------------------------------------------------------------
 ! 下面
 do i = i1 + 1, i2
  if( y(i,je) == minval(y(i1+1:i2,je)) ) exit
 enddo
 ic = i
 do i = is, ic - 1
  y(i,je) = minval(y(ic:i2,je))
 enddo
 ! 上面
 do i = i2, i3 - 1
  if( y(i,je) == maxval(y(i2:i3-1,je)) ) exit
 enddo
 ic = i
 do i = ic + 1, ie
  y(i,je) = maxval(y(i2:ic,je))
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine CtypeExternalBoundary
!*******************************************************************************************************
!******** 側部境界 (C-type)									********
!*******************************************************************************************************
subroutine CtypeSideBoundary( &
&            is, ie, js, je, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je
 real   , intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: j
 ! 処理開始 ********************************************************************************************
 do j = js + 1, je - 1
  x(is,j) = x(is,js)
  x(ie,j) = x(ie,js)
  y(is,j) = ( y(is,je) - y(is,js) ) * real(j) / real(je) + y(is,js)
  y(ie,j) = ( y(ie,je) - y(ie,js) ) * real(j) / real(je) + y(ie,js)
 enddo
 ! 処理開始 ********************************************************************************************
 return
end subroutine CtypeSideBoundary
!*******************************************************************************************************
!******** 楕円-双曲型偏微分方程式法に基づく格子生成						********
!*******************************************************************************************************
subroutine CtypeGenerationEHPDE( &
&            is, i1, i2, i3, ie, js, je, rs, re, margin, resi, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, i1, i2, i3, ie, js, je
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
    cv1m(i,j) = 1.0
   else if(i <= i3) then
    cv1m(i,j) = real(j) / real(je) * real(i3-i) / real(i3-i1)
   else
    cv1m(i,j) = 0.0
  endif	
  ! i + 方向
  if(i3 < i) then
    cv1p(i,j) = 1.0
   else if(i1 <= i) then
    cv1p(i,j) = real(j) / real(je) * real(i-i1) / real(i3-i1)
   else
    cv1p(i,j) = 0.0
  endif
  ! j - 方向
  cv2m(i,j) = (real(je-j) / real(je))**2
  ! j + 方向
  if(i2 < i) then
    cv2p(i,j) = (real(ie-i) / real(ie))**2
   else if(i == i2) then
    cv2p(i,j) = 0.0
   else
    cv2p(i,j) = (real(i) / real(ie))**2
  endif
!  cv2p(i,j) = (real(j-js) / real(je))**1 * 0.3
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
   else if(i <= i3) then
    dx1m(i,j) = ( -x(i-1,js) + x(i,js)  ) *   real(je-j) / real(je) &
    &         + ( -x(i-1,je) + x(i,je)  ) * ( real(j-js) / real(je) )**2
    dy1m(i,j) = ( -y(i-1,js) + y(i,js)  ) *   real(je-j) / real(je) &
    &         + ( -y(i-1,je) + y(i,je)  ) * ( real(j-js) / real(je) )**2
   else
    dx1m(i,j) = 0.0
    dy1m(i,j) = 0.0
  endif
  ! i + 方向
  if(i3 < i) then
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
  do j = js+1, je-1
   y(is,j) = y(is+1,j)
   y(ie,j) = y(ie-1,j)
  enddo
 enddo
 if( dmax > resi ) then
  write(*, '(a)') "!!!!! Elliptic-Hyperbolic calculation error !!!!!"
  stop
 endif
 ! 処理終了 ********************************************************************************************
 return
end subroutine CtypeGenerationEHPDE
!*******************************************************************************************************
!********* 翼前縁付近の着氷計算用の格子 (H-type)						********
!*******************************************************************************************************
subroutine HtypeGridIceLE( &
&            m, is, ie, js, je, ks, ke, dom1, dom2, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer , intent(in)  :: m
 integer, intent(in)  :: is, ie, js, je, ks, ke
 real   , intent(in)  :: dom1, dom2
 real   , intent(out) :: x(is:ie, js:je, ks:ke), &
 &                       y(is:ie, js:je, ks:ke), &
 &                       z(is:ie, js:je, ks:ke)
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: rs1 = 1.0e-1 !7.0 * 1.0e-2				! 前縁付近格子幅
 real   , parameter :: re1 = 8.0e-1 !2.0 * 1.0e-0				! 着氷限界位置付近格子幅
 real   , parameter :: rs2 = 2.0e0 !5.0 * 1.0e-1				! 翼表面境界格子幅
 real   , parameter :: re2 = 1.0e0 !7.0 * 1.0e-2				! 翼遠方境界格子幅
 real   , parameter :: tb1 = 1.0e-2 !1.0 * 1.0e-1				! 内部境界直交性のパラメータ
 real   , parameter :: tb2 = 1.0e-1 !2.0 * 1.0e-1				! 外部境界直交性のパラメータ
 ! 処理開始 ********************************************************************************************
 ! 内部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call HtypeInternalBoundary( &
 &      is, ie, dom1, rs1, re1, x(:, js, ks), y(:, js, ks) )
 ! 外部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call HtypeExternalBoundary( &
 &      is, ie, js, je, dom2, x(:, :, ks), y(:, :, ks) )
! ! 側部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call HtypeSideBoundary( &
! &      is, ie, js, je, rs2, re2, x(:, :, ks), y(:, :, ks) )
! ! Transfinite 補間 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call TransfiniteInterpolation( &
! &      is, ie, js, je, x(:, :, ks), y(:, :, ks) )
 ! 二境界法 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call GenerationTwoBoundary( &
 &      m, is, ie, -1, -1, js, je, rs2, re2, tb1, tb2, x(:, :, ks), y(:, :, ks) )
 ! 三次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call ThreeDimensionalized( &
 &      is, ie, js, je, ks, ke, span, x, y, z )
 ! 処理終了 ********************************************************************************************
 return
end subroutine HtypeGridIceLE
!*******************************************************************************************************
!******** 内部境界 (H-type)									********
!*******************************************************************************************************
subroutine HtypeInternalBoundary( &
&            is, ie, dom, rs, re, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: dom
 real   , intent(in)  :: rs, re
 real   , intent(out) :: x(is:ie), y(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: xu(:), xl(:), yu(:), yl(:)
 real   , pointer :: tu(:), tl(:), lu(:), ll(:)
 real   , pointer :: x0(:), y0(:)
 real   , pointer :: sy(:)
 integer :: n, i
 integer :: n1, n2, n3, n4, n5
 integer :: si
 integer :: i1
 integer :: nmax
 real    :: xx, yy, yt, yc, theta
 real    :: rr
 ! 処理開始 ********************************************************************************************
 ! 前縁格子番号 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 i1 = int(0.5 * ie)
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 nmax = 10000
 allocate( xu(0:nmax), yu(0:nmax), xl(0:nmax), yl(0:nmax) )
 allocate( tu(i1:ie), tl(is:i1), lu(0:nmax), ll(0:nmax) )
 allocate( x0(is:i1), y0(is:i1) )
 ! 翼周り ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 翼形状 ----------------------------------------------------------------------------------------------
 do n = 1, nmax - 1
  ! コード方向位置
  xx = chord * real(n) / real(nmax)
  ! 翼厚分布
  yt = 5.0 * thick * ( a0 * (xx / chord)**0.5 + a1 * (xx / chord)**1 + a2 * (xx / chord)**2 &
  &                  + a3 * (xx / chord)**3   + a4 * (xx / chord)**4 )
  ! 翼反り分布
  if(xc == 0.0) then
    yc = 0.0
   else
    if(xx < xc) then
      yc = ycmax * chord**2 / xc**2 * (2.0 * xc * xx / chord**2 - xx**2 / chord**2)
     else
      yc = ycmax * chord**2 / (chord - xc)**2 &
      &  * (1.0 - 2.0 * xc / chord + 2.0 * xc * xx / chord**2 - xx**2 / chord**2)
    endif
  endif
  ! 翼座標
  theta = atan(yc / xx)
  xu(n) = xx - yt * sin(theta)
  yu(n) = yc + yt * cos(theta)
  xl(n) = xx + yt * sin(theta)
  yl(n) = yc - yt * cos(theta)
 enddo
 xu(0) = 0.0; xu(nmax) = chord
 yu(0) = 0.0; yu(nmax) = 0.0
 xl(0) = 0.0; xl(nmax) = chord
 yl(0) = 0.0; yl(nmax) = 0.0
 ! 計算領域 --------------------------------------------------------------------------------------------
 do n = 0, nmax
  if( xu(n) - xu(0) > dom ) exit
 enddo
 nmax = n
 ! 格子幅調整 ------------------------------------------------------------------------------------------
 ! 媒介変数
 call VinokurInterpolation( rs / real(ie-i1+1), re / real(ie-i1+1), ie-i1+1, tu )
 call VinokurInterpolation( rs / real(i1-is+1), re / real(i1-is+1), i1-is+1, tl )
 ! 前縁からの距離
 lu(0) = 0.0; ll(0) = 0.0
 do n = 1, nmax
  lu(n) = sqrt( (xu(n) - xu(n-1))**2 + (yu(n) - yu(n-1))**2 ) + lu(n-1)
  ll(n) = sqrt( (xl(n) - xl(n-1))**2 + (yl(n) - yl(n-1))**2 ) + ll(n-1)
 enddo
 ! 無次元化
 do n = 0, nmax
  lu(n) = lu(n) / lu(nmax)
  ll(n) = ll(n) / ll(nmax)
 enddo
 ! ラグランジュ補間で探索 (上面)
 do i = i1 + 1, ie - 1
  do n = 0, nmax - 1
   if( lu(n) <= tu(i) .and. tu(i) < lu(n+1) ) then
     if(n <= 0) then
       n1 = n; n2 = n+1; n3 = n+2
       call LagrangeInterpolation( 3, lu(n1:n3), xu(n1:n3), tu(i), x(i) )
       call LagrangeInterpolation( 3, lu(n1:n3), yu(n1:n3), tu(i), y(i) )
      else if(n <= 1) then
       n1 = n-1; n2 = n; n3 = n+1; n4 = n+2
       call LagrangeInterpolation( 4, lu(n1:n4), xu(n1:n4), tu(i), x(i) )
       call LagrangeInterpolation( 4, lu(n1:n4), yu(n1:n4), tu(i), y(i) )
      else if(n < nmax-1) then
       n1 = n-2; n2 = n-1; n3 = n; n4 = n+1; n5 = n+2
       call LagrangeInterpolation( 5, lu(n1:n5), xu(n1:n5), tu(i), x(i) )
       call LagrangeInterpolation( 5, lu(n1:n5), yu(n1:n5), tu(i), y(i) )
      else
       n1 = n-2; n2 = n-1; n3 = n; n4 = n+1
       call LagrangeInterpolation( 4, lu(n1:n4), xu(n1:n4), tu(i), x(i) )
       call LagrangeInterpolation( 4, lu(n1:n4), yu(n1:n4), tu(i), y(i) )
    endif
    exit
   endif
  enddo
 enddo
 x(i1) = xu(0); x(ie) = xu(nmax)
 y(i1) = yu(0); y(ie) = yu(nmax)
 ! ラグランジュ補間で探索 (下面)
 do i = is + 1, i1 - 1
  do n = 0, nmax - 1
   if( ll(n) <= tl(i) .and. tl(i) < ll(n+1) ) then
     if(n <= 0) then
       n1 = n; n2 = n+1; n3 = n+2
       call LagrangeInterpolation( 3, ll(n1:n3), xl(n1:n3), tl(i), x(i) )
       call LagrangeInterpolation( 3, ll(n1:n3), yl(n1:n3), tl(i), y(i) )
      else if(n <= 1) then
       n1 = n-1; n2 = n; n3 = n+1; n4 = n+2
       call LagrangeInterpolation( 4, ll(n1:n4), xl(n1:n4), tl(i), x(i) )
       call LagrangeInterpolation( 4, ll(n1:n4), yl(n1:n4), tl(i), y(i) )
      else if(n < nmax-1) then
       n1 = n-2; n2 = n-1; n3 = n; n4 = n+1; n5 = n+2
       call LagrangeInterpolation( 5, ll(n1:n5), xl(n1:n5), tl(i), x(i) )
       call LagrangeInterpolation( 5, ll(n1:n5), yl(n1:n5), tl(i), y(i) )
      else
       n1 = n-2; n2 = n-1; n3 = n; n4 = n+1
       call LagrangeInterpolation( 4, ll(n1:n4), xl(n1:n4), tl(i), x(i) )
       call LagrangeInterpolation( 4, ll(n1:n4), yl(n1:n4), tl(i), y(i) )
    endif
    exit
   endif
  enddo
 enddo
 x(is) = xu(0); x(i1) = xl(nmax)
 y(is) = yu(0); y(i1) = yl(nmax)
 ! 格子番号整理
 do i = is, i1
  x0(i) = x(is+i1-i)
  y0(i) = y(is+i1-i)
 enddo
 do i = is, i1
  x(i) = x0(i)
  y(i) = y0(i)
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine HtypeInternalBoundary
!*******************************************************************************************************
!******** 外部境界 (H-type)									********
!*******************************************************************************************************
subroutine HtypeExternalBoundary( &
&            is, ie, js, je, dom, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je
 real   , intent(in)    :: dom
 real   , intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 integer :: ic
 real    :: a1, a2, a3, b1, b2, b3, c1, c2, c3
 ! 処理開始 ********************************************************************************************
 do i = is, ie
  ! 翼周り方向ベクトル
  if(i == is) then
    a1 = - x(i  ,js) + x(i+1,js)
    a2 = - y(i  ,js) + y(i+1,js)
   else if(i == ie) then
    a1 = - x(i-1,js) + x(i  ,js)
    a2 = - y(i-1,js) + y(i  ,js)
   else
    a1 = 0.5 * (- x(i-1,js) + x(i+1,js))
    a2 = 0.5 * (- y(i-1,js) + y(i+1,js))
  endif
  a3 = 0.0
  ! スパン方向ベクトル
  b1 =  0.0
  b2 =  0.0
  b3 = -1.0
  ! 法線方向ベクトル
  c1 = a2 * b3 - a3 * b2
  c2 = a3 * b1 - a1 * b3
  c3 = a1 * b2 - a2 * b1
  ! 翼周り部外部領域
  x(i,je) = c1 / sqrt(c1**2 + c2**2 + c3**2) * dom + x(i,js)
  y(i,je) = c2 / sqrt(c1**2 + c2**2 + c3**2) * dom + y(i,js)
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine HtypeExternalBoundary
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
!******** O 型格子内部境界									********
!*******************************************************************************************************
subroutine OtypeInnerBoundary( &
&            is, ie, i1, i2, i3, dom2, rs1, re1, rs2, re2, rs3, re3, x, y )
! 変数型宣言 ********************************************************************************************
  implicit none
! 変数宣言 ***********************************************************************************************
  ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, intent(in)  :: is, ie							! 格子数
    integer, intent(in)  :: i1, i2, i3						! 前縁付近の格子番号
    real   , intent(in)  :: dom2							! 氷着くとこ
    real   , intent(in)  :: rs1, re1, rs2, re2, rs3, re3				! 格子幅
    real   , intent(inout) :: x(is:ie), y(is:ie)					! 座標
  ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, parameter :: nmax = 10000
  ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real   , pointer :: xs(:), ys(:), xp(:), yp(:), ts0(:), ts(:), tp0(:), tp(:), ls(:), lp(:)
    integer :: i
    integer :: n, n1, n2, n3, n4, n5
    integer :: smax, ip
    real    :: tt, rr
    real    :: xx, yt, yc, theta
! 処理開始 ********************************************************************************************

 allocate( xs(0:nmax), ys(0:nmax), xp(0:nmax), yp(0:nmax) )
 ! 翼周り ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 翼形状 ----------------------------------------------------------------------------------------------
 do n = 1, nmax - 1
  ! コード方向位置
  xx = chord * real(n) / real(nmax)
  ! 翼厚分布
  yt = 5.0 * thick * ( a0 * (xx / chord)**0.5 + a1 * (xx / chord)**1 + a2 * (xx / chord)**2 &
  &                  + a3 * (xx / chord)**3   + a4 * (xx / chord)**4 )
  ! 翼反り分布
  if(xc == 0.0) then
    yc = 0.0
   else
    if(xx < xc) then
      yc = ycmax * chord**2 / xc**2 * (2.0 * xc * xx / chord**2 - xx**2 / chord**2)
     else
      yc = ycmax * chord**2 / (chord - xc)**2 &
      &  * (1.0 - 2.0 * xc / chord + 2.0 * xc * xx / chord**2 - xx**2 / chord**2)
    endif
  endif
  ! 翼座標
  theta = atan(yc / xx)
  xp(n) = xx - yt * sin(theta)
  yp(n) = yc + yt * cos(theta)
  xs(nmax-n) = xx + yt * sin(theta)
  ys(nmax-n) = yc - yt * cos(theta)
 enddo
 xp(0) = 0.0; xp(nmax) = chord
 yp(0) = 0.0; yp(nmax) = 0.0
 xs(nmax) = 0.0; xs(0) = chord
 ys(nmax) = 0.0; ys(0) = 0.0
!! 格子点めっちゃ増やす --------------------------------------------------------------------------------
!  ! 構造体メモリ確保
!    allocate( xs(0:nmax), ys(0:nmax), xp(0:nmax), yp(0:nmax) )
!    DO n = 0, nmax
!      tt = real(i2 - is - 4 + 2) * real(n) / real(nmax)
!      ! S.S.
!      call BSplineCurve( &
!      &      i2-is, x(is:i2), y(is:i2), 4, tt, xs(n), ys(n) )
!      tt = real(ie - i2 - 4 + 2) * real(n) / real(nmax)
!      ! P.S.
!      call BSplineCurve( &
!      &      ie-i2, x(i2:ie), y(i2:ie), 4, tt, xp(n), yp(n) )
!    END DO
! 媒介変数 --------------------------------------------------------------------------------------------
  ! 構造体メモリ確保
    allocate( ts(is:i2), tp(i2:ie) )
  ! S.S.
    allocate( ts0(is:i1) )
    call VinokurInterpolation( rs1 / real(i1-is+1), re1 / real(i1-is+1), i1-is+1, ts0 )
    DO i = is, i1
      ts(i) = (1.0 - dom2) * ts0(i)
    END DO
    deallocate( ts0 )
    allocate( ts0(i1:i2) )
    call GeometricInterpolationInv(re2 / real(i2-i1+1), i2-i1+1, ts0, rr) 
    DO i = i1 + 1, i2
      ts(i) = ts(i1) + dom2 * ts0(i)
    END DO
    deallocate( ts0 )
  ! P.S.
    allocate( tp0(i2:i3) )
    call GeometricInterpolation(rs2 / real(i3-i2+1), i3-i2+1, tp0, rr) 
    DO i = i2, i3
      tp(i) = dom2 * tp0(i)
    END DO
    deallocate( tp0 )
    allocate( tp0(i3:ie) )
    call VinokurInterpolation( rs3 / real(ie-i3+1), re3 / real(ie-i3+1), ie-i3+1, tp0 )
    DO i = i3, ie
      tp(i) = dom2 + (1.0 - dom2) * tp0(i)
    END DO
! 起点からの距離 -------------------------------------------------------------------------------------
  allocate( ls(0:nmax), lp(0:nmax) )
  ls(0) = 0.0; lp(0) = 0.0
  DO n = 1, nmax
   ! S.S.
   ls(n) = sqrt( (xs(n) - xs(n-1))**2 + (ys(n) - ys(n-1))**2 ) + ls(n-1)
   ! P.S.
   lp(n) = sqrt( (xp(n) - xp(n-1))**2 + (yp(n) - yp(n-1))**2 ) + lp(n-1)
  END DO
  ! 無次元化
  DO n = 0, nmax
   ls(n) = ls(n) / ls(nmax)
   lp(n) = lp(n) / lp(nmax)
  END DO
! ラグランジュ補間 -----------------------------------------------------------------------------------
  ! S.S.
  DO i = is + 1, i2 - 1
    DO n = 0, nmax - 1
    IF( ls(n) <= ts(i) .and. ts(i) < ls(n+1) ) THEN
    select case(n)
    case default
     n1 = n-2; n2 = n-1; n3 = n; n4 = n+1; n5 = n+2
     call LagrangeInterpolation( 5, ls(n1:n5), xs(n1:n5), ts(i), x(i) )
     call LagrangeInterpolation( 5, ls(n1:n5), ys(n1:n5), ts(i), y(i) )
    case(0)
     n1 = n; n2 = n+1; n3 = n+2
     call LagrangeInterpolation( 3, ls(n1:n3), xs(n1:n3), ts(i), x(i) )
     call LagrangeInterpolation( 3, ls(n1:n3), ys(n1:n3), ts(i), y(i) )
    case(1)
     n1 = n-1; n2 = n; n3 = n+1; n4 = n+2
     call LagrangeInterpolation( 4, ls(n1:n4), xs(n1:n4), ts(i), x(i) )
     call LagrangeInterpolation( 4, ls(n1:n4), ys(n1:n4), ts(i), y(i) )
    case(nmax-1)
     n1 = n-2; n2 = n-1; n3 = n; n4 = n+1
     call LagrangeInterpolation( 4, ls(n1:n4), xs(n1:n4), ts(i), x(i) )
     call LagrangeInterpolation( 4, ls(n1:n4), ys(n1:n4), ts(i), y(i) )
    end select
    exit
    END IF
    END DO
  END DO
  x(is) = xs(0); x(i2) = xs(nmax)
  y(is) = ys(0); y(i2) = ys(nmax)
  ! P.S.
  DO i = i2 + 1, ie - 1
    DO n = 0, nmax - 1
    IF( lp(n) <= tp(i) .and. tp(i) < lp(n+1) ) THEN
      select case(n)
       case default
        n1 = n-2; n2 = n-1; n3 = n; n4 = n+1; n5 = n+2
        call LagrangeInterpolation( 5, lp(n1:n5), xp(n1:n5), tp(i), x(i) )
        call LagrangeInterpolation( 5, lp(n1:n5), yp(n1:n5), tp(i), y(i) )
       case(0)
        n1 = n; n2 = n+1; n3 = n+2
        call LagrangeInterpolation( 3, lp(n1:n3), xp(n1:n3), tp(i), x(i) )
        call LagrangeInterpolation( 3, lp(n1:n3), yp(n1:n3), tp(i), y(i) )
       case(1)
        n1 = n-1; n2 = n; n3 = n+1; n4 = n+2
        call LagrangeInterpolation( 4, lp(n1:n4), xp(n1:n4), tp(i), x(i) )
        call LagrangeInterpolation( 4, lp(n1:n4), yp(n1:n4), tp(i), y(i) )
       case(nmax-1)
        n1 = n-2; n2 = n-1; n3 = n; n4 = n+1
        call LagrangeInterpolation( 4, lp(n1:n4), xp(n1:n4), tp(i), x(i) )
        call LagrangeInterpolation( 4, lp(n1:n4), yp(n1:n4), tp(i), y(i) )
      end select
      exit
    END IF
    END DO
  END DO
  x(i2) = xp(0); x(ie) = xp(nmax)
  y(i2) = yp(0); y(ie) = yp(nmax)
! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  deallocate( xs, ys, xp, yp, ts, tp, ls, lp )
! 端点のスムージング ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  smax = 50; ip = 4
  allocate( xs(is:ie), ys(is:ie) )
  DO n = 1, smax
  xs(:) = x(:); ys(:) = y(:)
    DO i = is, is + ip
    IF(i == is) THEN
     x(i) = (xs(ie-1) + xs(i) + xs(i+1)) / 3.0
     y(i) = (ys(ie-1) + ys(i) + ys(i+1)) / 3.0
    ELSE
     x(i) = (xs(i-1) + xs(i) + xs(i+1)) / 3.0
     y(i) = (ys(i-1) + ys(i) + ys(i+1)) / 3.0
    END IF
    END DO
    DO i = ie - ip, ie
    IF(i == ie) THEN
      x(i) = (xs(i-1) + xs(i) + xs(is+1)) / 3.0
      y(i) = (ys(i-1) + ys(i) + ys(is+1)) / 3.0
    ELSE
      x(i) = (xs(i-1) + xs(i) + xs(i+1)) / 3.0
      y(i) = (ys(i-1) + ys(i) + ys(i+1)) / 3.0
    END IF
    END DO
  END DO
! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  deallocate(xs, ys)
!!!!!!!
!allocate(xs(is:ie), ys(is:ie))
!do i = is, ie
!xs(i) = sqrt((x(i) - x(i2))**2 + (y(i) - y(i2))**2)
!enddo
!tt = 0.0
!do i = is, i2 - 1
!tt = tt + xs(i) - xs(i+1)
!enddo
!do i = is, ie
!if(i < i2) then
!  ys(i) = -xs(i) / tt
! else
!  ys(i) = xs(i) / tt
!endif
!enddo
!do i = is, ie
!write(*,*) ys(i)
!enddo
!stop
!!!!!!!
! 処理終了 ********************************************************************************************
  return
END SUBROUTINE OtypeInnerBoundary

!*******************************************************************************************************
!******** O 型格子外部境界									********
!*******************************************************************************************************
SUBROUTINE OtypeOuterBoundary( &
&            is, ie, js, je, i1, i2, i3, dom1, dom2, rs1, re1, rs2, re2, x, y )
! 変数型宣言 ********************************************************************************************
  implicit none
! 変数宣言 ***********************************************************************************************
  ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, intent(in)		:: is, ie, js, je				! 格子数
    integer, intent(in)		:: i1, i2, i3					! 前縁格子番号
    real   , intent(in)		:: dom1						! 計算領域
    real   , intent(in)		:: dom2						! 氷着くとこ
    real   , intent(in)		:: rs1, re1, rs2, re2				! 格子幅
    real   , intent(inout)	:: x(is:ie, js:je), y(is:ie, js:je)		! 格子座標
  ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, parameter		:: nmax = 10000
  ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real   , pointer		:: xs(:), ys(:), xp(:), yp(:), ts(:), tp0(:), tp(:), ls(:), lp(:)
    real   , pointer		:: h1(:), h2(:), h3(:), h4(:)
    integer			:: i
    integer			:: n, n1, n2, n3, n4, n5
    real			:: tt
    real			:: a1, a2, c1, c2
    integer			:: smax, ip
! 処理開始 ********************************************************************************************
! 外部境界座標 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  allocate( h1(is:ie), h2(is:ie), h3(is:ie), h4(is:ie) )
  DO i = is, ie - 1
    ! 翼周り方向ベクトル
      IF(i == is) THEN
        a1 = 0.5 * (- x(is,js) + x(is+1,js) - x(ie-1,js) + x(ie,js))
        a2 = 0.5 * (- y(is,js) + y(is+1,js) - y(ie-1,js) + y(ie,js))
      ELSE
        a1 = 0.5 * (- x(i-1,js) + x(i+1,js))
        a2 = 0.5 * (- y(i-1,js) + y(i+1,js))
      END IF
    ! 法線方向ベクトル
      c1 = -a2
      c2 =  a1
      h1(i) = c1 / sqrt(c1**2 + c2**2)
      h2(i) = c2 / sqrt(c1**2 + c2**2)
  END DO
  h1(ie) = h1(is)
  h2(ie) = h2(is)
! スムージング +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO n = 1, 1
  h3(:) = h1(:); h4(:) = h2(:)
    DO i = is, ie
      IF(i == is) THEN
        h1(i) = (h3(ie-1) + h3(i) + h3(i+1)) / 3.0
        h2(i) = (h4(ie-1) + h4(i) + h4(i+1)) / 3.0
      ELSE IF(i == ie) THEN
        h1(i) = (h3(i-1) + h3(i) + h3(is+1)) / 3.0
        h2(i) = (h4(i-1) + h4(i) + h4(is+1)) / 3.0
      ELSE
        h1(i) = (h3(i-1) + h3(i) + h3(i+1)) / 3.0
        h2(i) = (h4(i-1) + h4(i) + h4(i+1)) / 3.0
      END IF
    END DO
  END DO
! 座標設定 ---------------------------------------------------------------------------------------------
  DO i = is, ie
    x(i,je) = h1(i) * dom1 + x(i,js)
    y(i,je) = h2(i) * dom1 + y(i,js)

!write(*,*) x(i,je),y(i,je)

  END DO
! 格子幅 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! 格子点めっちゃ増やす --------------------------------------------------------------------------------
    allocate( xs(0:nmax), ys(0:nmax), xp(0:nmax), yp(0:nmax) )
    DO n = 0, nmax
    tt = real(i2 - is - 4 + 2) * real(n) / real(nmax)
    ! S.S.
    call BSplineCurve( &
    &      i2-is, x(is:i2,je), y(is:i2,je), 4, tt, xs(n), ys(n) )
    tt = real(ie - i2 - 4 + 2) * real(n) / real(nmax)
    ! P.S.
    call BSplineCurve( &
    &      ie-i2, x(i2:ie,je), y(i2:ie,je), 4, tt, xp(n), yp(n) )
    END DO
  ! 媒介変数 -------------------------------------------------------------------------------------------
    allocate( ts(is:i2), tp0(i2:ie), tp(i2:ie) )
    ! S.S.
    IF(rs1 /= 0.0 .and. re1 /= 0.0) THEN
      call VinokurInterpolation( rs1 / real(i1-is+1), re1 / real(i1-is+1), i1-is+1, ts )
      DO i = is, i1
      ts(i) = (1.0 - dom2) * ts(i)
      END DO
      DO i = i1 + 1, i2
      ts(i) = ts(i1) + dom2 * real(i - i1) / real(i2 - i1)
      END DO
    ELSE
      DO i = is, i2
      ts(i) = real(i-is) / real(i2-is)
      END DO
    END IF
    ! P.S.
    IF(rs2 /= 0.0 .and. re2 /= 0.0) THEN
      call VinokurInterpolation( rs2 / real(ie-i3+1), re2 / real(ie-i3+1), ie-i3+1, tp0 )
      DO i = i3, ie
        tp(i) = dom2 + (1.0 - dom2) * tp0(i-(i3-i2))
      END DO
      DO i = i2, i3 - 1
        tp(i) = dom2 * real(i - i2) / real(i3 - i2)
      END DO
    ELSE
      DO i = i2, ie
        tp(i) = real(i-i2) / real(ie-i2)
      END DO
    END IF
  ! 起点からの距離 -------------------------------------------------------------------------------------
    allocate( ls(0:nmax), lp(0:nmax) )
    ls(0) = 0.0; lp(0) = 0.0
    DO n = 1, nmax
      ! S.S.
      ls(n) = sqrt( (xs(n) - xs(n-1))**2 + (ys(n) - ys(n-1))**2 ) + ls(n-1)
      ! P.S.
      lp(n) = sqrt( (xp(n) - xp(n-1))**2 + (yp(n) - yp(n-1))**2 ) + lp(n-1)
    END DO
  ! 無次元化
    DO n = 0, nmax
      ls(n) = ls(n) / ls(nmax)
      lp(n) = lp(n) / lp(nmax)
    END DO
  ! ラグランジュ補間 -----------------------------------------------------------------------------------
    ! S.S.
      DO i = is + 1, i2 - 1
        DO n = 0, nmax - 1
        IF( ls(n) <= ts(i) .and. ts(i) < ls(n+1) ) THEN
          select case(n)
           case default
            n1 = n-2; n2 = n-1; n3 = n; n4 = n+1; n5 = n+2
            call LagrangeInterpolation( 5, ls(n1:n5), xs(n1:n5), ts(i), x(i,je) )
            call LagrangeInterpolation( 5, ls(n1:n5), ys(n1:n5), ts(i), y(i,je) )
           case(0)
            n1 = n; n2 = n+1; n3 = n+2
            call LagrangeInterpolation( 3, ls(n1:n3), xs(n1:n3), ts(i), x(i,je) )
            call LagrangeInterpolation( 3, ls(n1:n3), ys(n1:n3), ts(i), y(i,je) )
           case(1)
            n1 = n-1; n2 = n; n3 = n+1; n4 = n+2
            call LagrangeInterpolation( 4, ls(n1:n4), xs(n1:n4), ts(i), x(i,je) )
            call LagrangeInterpolation( 4, ls(n1:n4), ys(n1:n4), ts(i), y(i,je) )
           case(nmax-1)
            n1 = n-2; n2 = n-1; n3 = n; n4 = n+1
            call LagrangeInterpolation( 4, ls(n1:n4), xs(n1:n4), ts(i), x(i,je) )
            call LagrangeInterpolation( 4, ls(n1:n4), ys(n1:n4), ts(i), y(i,je) )
          END select
          exit
        END IF
        END DO
      END DO
      x(is,je) = xs(0); x(i2,je) = xs(nmax)
      y(is,je) = ys(0); y(i2,je) = ys(nmax)
    ! P.S.
      DO i = i2 + 1, ie - 1
        DO n = 0, nmax - 1
        IF( lp(n) <= tp(i) .and. tp(i) < lp(n+1) ) then
        select case(n)
           case default
            n1 = n-2; n2 = n-1; n3 = n; n4 = n+1; n5 = n+2
            call LagrangeInterpolation( 5, lp(n1:n5), xp(n1:n5), tp(i), x(i,je) )
            call LagrangeInterpolation( 5, lp(n1:n5), yp(n1:n5), tp(i), y(i,je) )
           case(0)
            n1 = n; n2 = n+1; n3 = n+2
            call LagrangeInterpolation( 3, lp(n1:n3), xp(n1:n3), tp(i), x(i,je) )
            call LagrangeInterpolation( 3, lp(n1:n3), yp(n1:n3), tp(i), y(i,je) )
           case(1)
            n1 = n-1; n2 = n; n3 = n+1; n4 = n+2
            call LagrangeInterpolation( 4, lp(n1:n4), xp(n1:n4), tp(i), x(i,je) )
            call LagrangeInterpolation( 4, lp(n1:n4), yp(n1:n4), tp(i), y(i,je) )
           case(nmax-1)
            n1 = n-2; n2 = n-1; n3 = n; n4 = n+1
            call LagrangeInterpolation( 4, lp(n1:n4), xp(n1:n4), tp(i), x(i,je) )
            call LagrangeInterpolation( 4, lp(n1:n4), yp(n1:n4), tp(i), y(i,je) )
        END SELECT
        EXIT
        END IF
        END DO
      END DO
      x(i2,je) = xp(0); x(ie,je) = xp(nmax)
      y(i2,je) = yp(0); y(ie,je) = yp(nmax)
      deallocate( xs, ys, xp, yp, ts, tp0, tp, ls, lp )
! 端点のスムージング ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  smax = 10; ip = 30
  ! 構造体メモリ確保
    allocate( xs(is:ie), ys(is:ie) )
  DO n = 1, smax
  xs(:) = x(:,je); ys(:) = y(:,je)
    DO i = i2 - ip, i2 + ip
      x(i,je) = (xs(i-1) + xs(i) + xs(i+1)) / 3.0
      y(i,je) = (ys(i-1) + ys(i) + ys(i+1)) / 3.0
    END DO
  END DO
  ! 構造体メモリ解法
    deallocate( xs, ys )
  smax = 5; ip = 20
  ! 構造体メモリ確保
    allocate( xs(is:ie), ys(is:ie) )
  DO n = 1, smax
  xs(:) = x(:,je); ys(:) = y(:,je)
    DO i = i2 - ip, i2 + ip
      x(i,je) = (xs(i-1) + xs(i) + xs(i+1)) / 3.0
      y(i,je) = (ys(i-1) + ys(i) + ys(i+1)) / 3.0
    END DO
    DO i = is, is + ip
      IF(i == is) THEN
         x(i,je) = (xs(ie-1) + xs(i) + xs(i+1)) / 3.0
         y(i,je) = (ys(ie-1) + ys(i) + ys(i+1)) / 3.0
      ELSE
         x(i,je) = (xs(i-1) + xs(i) + xs(i+1)) / 3.0
         y(i,je) = (ys(i-1) + ys(i) + ys(i+1)) / 3.0
      END IF
    END DO
    DO i = ie - ip, ie
      IF(i == ie) THEN
        x(i,je) = (xs(i-1) + xs(i) + xs(is+1)) / 3.0
        y(i,je) = (ys(i-1) + ys(i) + ys(is+1)) / 3.0
      ELSE
        x(i,je) = (xs(i-1) + xs(i) + xs(i+1)) / 3.0
        y(i,je) = (ys(i-1) + ys(i) + ys(i+1)) / 3.0
      END IF
    END DO
  END DO
  ! メモリ解放
    deallocate(xs, ys)
! 処理終了 ********************************************************************************************
  return
END SUBROUTINE OtypeOuterBoundary

!*******************************************************************************************************
!******** O 型格子側部境界									********
!*******************************************************************************************************
SUBROUTINE OtypeSideBoundary( &
&            is, ie, js, je, j1, dom3, rs, re, x, y )
! 変数型宣言 ********************************************************************************************
  implicit none
! 変数宣言 ***********************************************************************************************
  ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, intent(in)		:: is, ie, js, je				! 格子数
    integer, intent(in)		:: j1						! 氷着くとこまで
    real   , intent(in)		:: dom3						! 氷着くとこ
    real   , intent(in)		:: rs, re					! 格子幅
    real   , intent(inout)	:: x(is:ie, js:je), y(is:ie, js:je)		! 格子座標
  ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real   , pointer :: t0(:), t(:)
    integer :: j
    real    :: rr
! 処理開始 ********************************************************************************************
  allocate( t0(j1:je), t(js:je) )
  IF(rs /= 0.0) THEN
  call GeometricInterpolation( rs / real(je-j1+1), je-j1+1, t0, rr )
    DO j = j1, je
      t(j) = dom3 + (1.0 - dom3) * t0(j)
    END DO
    DO j = js, j1 - 1
      t(j) = dom3 * real(j) / real(j1)
    END DO
  ELSE
    DO j = js, je
      t(j) = real(j) / real(je)
    END DO
  END IF
  DO j = js, je
    x(is,j) = ( x(is,je) - x(is,js) ) * t(j) + x(is,js)
    x(ie,j) = ( x(ie,je) - x(ie,js) ) * t(j) + x(ie,js)
    y(is,j) = ( y(is,je) - y(is,js) ) * t(j) + y(is,js)
    y(ie,j) = ( y(ie,je) - y(ie,js) ) * t(j) + y(ie,js)
  END DO
! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 deallocate(t0, t)
! 処理開始 ********************************************************************************************
  return
END SUBROUTINE OtypeSideBoundary


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
!******** O 型格子生成 (楕円-双曲型偏微分法)							********
!*******************************************************************************************************
SUBROUTINE OtypeEHPDEGeneration( &
&            is, ie, js, je, rs, re, margin, resi, x, y )
! 変数型宣言 ********************************************************************************************
  implicit none
! 変数宣言 **********************************************************************************************
  ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, intent(in)		:: is, ie, js, je				! 格子数
    real   , intent(in)		:: rs, re					! 境界格子幅
    real   , intent(in)		:: margin					! 重み関数許容誤差
    real   , intent(in)		:: resi						! 収束判定
    real   , intent(inout)	:: x(is:ie, js:je), y(is:ie, js:je)		! 格子座標
  ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real   , parameter		:: omg  = 1.0					! 緩和係数
    integer, parameter		:: nmax = 100000				! 最大計算回数
  ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    real   , pointer		:: Cs(:, :)
    real   , pointer		:: cv1m(:, :), cv1p(:, :), cv2m(:, :), cv2p(:, :)
    real   , pointer		:: dx1m(:, :), dx1p(:, :), dy1m(:, :), dy1p(:, :)
    real   , pointer		:: dx2m(:, :), dx2p(:, :), dy2m(:, :), dy2p(:, :)
    real   , pointer		:: x2(:, :), y2(:, :)
    integer			:: i, j, n, im
    real			:: dmax
    real			:: a1, a2, c1, c2, hjs1, hjs2, hje1, hje2
! 処理開始 ********************************************************************************************
  ! 格子生成計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! 二点重ね -------------------------------------------------------------------------------------------
 im = ie + 1
   allocate( x2(is:im, js:je), y2(is:im, js:je) )
   x2(is:im-1,:) = x(:,:); x2(im,:) = x(is+1,:)
   y2(is:im-1,:) = y(:,:); y2(im,:) = y(is+1,:)
 ! 双曲型重み関数 -------------------------------------------------------------------------------------
 allocate( cs  (is:im, js:je), &
 &         cv1m(is:im, js:je), cv1p(is:im, js:je), cv2m(is:im, js:je), cv2p(is:im, js:je), &
 &         dx1m(is:im, js:je), dx1p(is:im, js:je), dy1m(is:im, js:je), dy1p(is:im, js:je), &
 &         dx2m(is:im, js:je), dx2p(is:im, js:je), dy2m(is:im, js:je), dy2p(is:im, js:je) )
 DO j = js + 1, je - 1
 DO i = is + 1, im - 1
  ! i- 方向
!  cv1m(i,j) = real(j) / real(je)**2 * real(im-i) / real(im-is)
  cv1m(i,j) = 0.0
  ! i+ 方向
!  cv1p(i,j) = real(j) / real(je)**2 * real(i-is) / real(im-is)
  cv1p(i,j) = 0.0
  ! j- 方向
  cv2m(i,j) = (real(je-j) / real(je))**2
  ! j+ 方向
  cv2p(i,j) = (real(j-js) / real(je))**2
!  cv2p(i,j) = 0.0
  END DO
  END DO
  ! 楕円型重み関数 -------------------------------------------------------------------------------------
  DO j = js + 1, je - 1
  DO i = is + 1, im - 1
  ! i- 方向
!  dx1m(i,j) = ( - x2(i-1,js) + x2(i,js)  ) *   real(je-j) / real(je) &
!  &         + ( - x2(i-1,je) + x2(i,je)  ) * ( real(j-js) / real(je) )**2
!  dy1m(i,j) = ( - y2(i-1,js) + y2(i,js)  ) *   real(je-j) / real(je) &
!  &         + ( - y2(i-1,je) + y2(i,je)  ) * ( real(j-js) / real(je) )**2
  dx1m(i,j) = 0.0
  dy1m(i,j) = 0.0
  ! i+ 方向
!  dx1p(i,j) = ( - x2(i,js) + x2(i+1,js) ) *   real(je-j) / real(je) &
!  &         + ( - x2(i,je) + x2(i+1,je) ) * ( real(j-js) / real(je) )**2
!  dy1p(i,j) = ( - y2(i,js) + y2(i+1,js) ) *   real(je-j) / real(je) &
!  &         + ( - y2(i,je) + y2(i+1,je) ) * ( real(j-js) / real(je) )**2
  dx1p(i,j) = 0.0
  dy1p(i,j) = 0.0
  ! j- 方向
  a1 =  0.5 * ( - x2(i-1,js) + x2(i+1,js) )
  a2 =  0.5 * ( - y2(i-1,js) + y2(i+1,js) )
  c1 = -a2
  c2 =  a1
  hjs1 = c1 / sqrt(c1**2 + c2**2)
  hjs2 = c2 / sqrt(c1**2 + c2**2)
  dx2m(i,j) = hjs1 * rs / real(je - js)
  dy2m(i,j) = hjs2 * rs / real(je - js)
  ! j+ 方向
  a1 =  0.5 * ( - x2(i-1,je) + x2(i+1,je) )
  a2 =  0.5 * ( - y2(i-1,je) + y2(i+1,je) )
  c1 = -a2
  c2 =  a1
  hje1 = c1 / sqrt(c1**2 + c2**2)
  hje2 = c2 / sqrt(c1**2 + c2**2)
  dx2p(i,j) = hje1 * re / real(je - js)
  dy2p(i,j) = hje2 * re / real(je - js)
  END DO
  END DO
! 重み関数
  DO j = js + 1, je - 1
  DO i = is + 1, im - 1
    Cs(i,j) = 1.0 - max( cv1m(i,j), cv1p(i,j), cv2m(i,j), cv2p(i,j) ) + margin
  END DO
  END DO
! 格子生成計算 ---------------------------------------------------------------------------------------
  DO n = 1, nmax
    ! ソルバー部
    call GridEHPDE2D( &
    &      omg, is, im, js, je, &
    &      cs, cv1m, cv1p, cv2m, cv2p, dx1m, dy1m, dx1p, dy1p, dx2m, dy2m, dx2p, dy2p, &
    &      x2, y2, dmax )
!    call GridEPDE2D( &
!    &      omg, is, im, js, je, &
!    &      x2, y2, dmax )
    IF( mod(n, 1000) == 0.0 ) write(*, "(a,i5,e16.8e3)") "/ Elliptic-Hyperbolic calculation...", n, dmax
    IF( dmax < resi ) THEN
    write(*, "(a,i5,e16.8e3)") "/ Elliptic-Hyperbolic calculation...", n, dmax
    EXIT
    END IF
    ! 境界条件
    DO j = js + 1, je - 1
     x2(is,j) = x2(im-1,j)
     y2(is,j) = y2(im-1,j)
     x2(im,j) = x2(is+1,j)
     y2(im,j) = y2(is+1,j)
    END DO
  END DO
  IF( dmax > resi ) THEN
  write(*, '(a)') "!!!!! Elliptic-Hyperbolic calculation error !!!!!"
  stop
  END IF
  x(:,:) = x2(is:im-1,:)
  y(:,:) = y2(is:im-1,:)
  deallocate(cs, cv1m, cv1p, cv2m, cv2p, dx1m, dx1p, dy1m, dy1p, dx2m, dx2p, dy2m, dy2p, x2, y2)
! 処理終了 ********************************************************************************************
  return
END SUBROUTINE OtypeEHPDEGeneration

!*******************************************************************************************************
!******** 二境界法に基づく格子生成								********
!*******************************************************************************************************
subroutine GenerationTwoBoundary( &
&            m, is, ie, i1, i3, js, je, rs, re, t1, t2, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: m
 integer, intent(in)    :: is, ie, js, je
 integer, intent(in)    :: i1, i3
 real   , intent(in)    :: rs, re, t1, t2
 real   , intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: etabar(:)
 real    :: rr
 ! 処理開始 ********************************************************************************************
 ! メモリ確保
 allocate( etabar(js:je) )
 ! 媒介変数
 call GeometricInterpolation( rs / real(je-js+1), je-js+1, etabar, rr )
 ! 二境界法
 call TwoBoundaryMethod2D( &
 &      m, is, ie, i1, i3, js, je, t1, t2, etabar, x, y, swi_subgrid )
 ! メモリ解放
 deallocate(etabar)
 ! 処理終了 ********************************************************************************************
 return
end subroutine GenerationTwoBoundary
!*******************************************************************************************************
!******** 三次元化 (スパン方向変化なし)								********
!*******************************************************************************************************
subroutine ThreeDimensionalized( &
&            is, ie, js, je, ks, ke, dom, x, y, z )
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
! 定義終了 *********************************************************************************************
end program GridGeneration_NACA
!*******************************************************************************************************
!******** 翼のフラグ										********
!*******************************************************************************************************
subroutine ViewBlade( &
&            m, is, ie, js, je, ks, ke, ibs, ibe, j0, f )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: m
 integer, intent(in)  :: is, ie, js, je, ks, ke
 integer, intent(in)  :: ibs, ibe, j0
 integer, intent(out) :: f(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 do k = ks, ke
 do j = js, je
 do i = is, ie
  if( ibs <= i .and. i <= ibe .and. j == j0 ) then
    f(i,j,k) = 1
   else
    f(i,j,k) = 0
!    f(i,j,k) = m
  endif
 enddo
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
end subroutine ViewBlade
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

subroutine Smoothing(is, ie, js, je, x, y )
!mainroutine variable
 integer,intent(in) :: is
 integer,intent(in) :: ie
 integer,intent(in) :: js
 integer,intent(in) :: je
 real,intent(inout) :: x(is:ie,js:je)
 real,intent(inout) :: y(is:ie,js:je)

!subroutine variable
 integer :: i,j
 integer :: n
 integer,parameter :: nmax = 2

 do n = 0,nmax-1 !add
  do i = is+1,is+80 !ie-1
   do j = js+1,je-1
    x(i,j) = 0.25 * (x(i-1,j) + x(i+1,j) + x(i,j-1) + x(i,j+1))
    y(i,j) = 0.25 * (y(i-1,j) + y(i+1,j) + y(i,j-1) + y(i,j+1))
   end do
  end do
  do i = ie-80,ie-1
   do j = js+1,je-1
    x(i,j) = 0.25 * (x(i-1,j) + x(i+1,j) + x(i,j-1) + x(i,j+1))
    y(i,j) = 0.25 * (y(i-1,j) + y(i+1,j) + y(i,j-1) + y(i,j+1))
   end do
  end do
 end do

end subroutine Smoothing
