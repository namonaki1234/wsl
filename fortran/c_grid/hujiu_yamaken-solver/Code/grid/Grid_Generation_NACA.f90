!*******************************************************************************************************
!*******************************************************************************************************
!******** 格子生成プログラム                                                                        ********
!******** (NACA0012翼，三次元，C-type，重合格子法)                                                ********
!********                                         　　　2014.04.10  PROGRAMMED BY RYOSUKE HAYASHI ********
!********                                              2014.04.15     UPDATED BY RYOSUKE HAYASHI ********
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
   character, parameter :: ViewGrdFile*8 = 'ViewGrid'
   ! 共有定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer, parameter :: NACA1 = 0                                        ! NACA 翼の反り比
   integer, parameter :: NACA2 = 0                                        ! NACA 翼の最大反り位置
   integer, parameter :: NACA34 = 12                                        ! NACA 翼の最大翼厚
   real, parameter :: a0 = 0.2969                                        ! 翼厚分布式の定数
   real, parameter :: a1 = -0.1260                                        ! 翼厚分布式の定数
   real, parameter :: a2 = -0.3516                                        ! 翼厚分布式の定数
   real, parameter :: a3 = 0.2843                                        ! 翼厚分布式の定数
   real, parameter :: a4 = -0.1015                                        ! 翼厚分布式の定数
   ! 共有変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   real   :: thick, ytmax, ycmax, xc
   ! 処理開始 ********************************************************************************************
   ! 検証実験ケース選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   write (*, '(a)') "<< Computational Condition >>"
   call SelectExpCase
   ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   write (*, '(a)') '<< Initial Setting >>'
   call InitialSetting
! ! Main Grid +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call MainGrid
! write(*, '(a)') '+++++ Main grid complete. +++++'
   ! Sub Grid +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   call SubGrid
   write (*, '(a)') '+++++ Sub grid complete. +++++'
   ! 内部手続き ******************************************************************************************
   stop
contains
!*******************************************************************************************************
!******** 検証実験ケース選択                                                                         ********
!*******************************************************************************************************
   subroutine SelectExpCase
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      character :: fname*20
      ! 処理開始 ********************************************************************************************
      ! 実験条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 実験ケース入力
      write (*, '(a)') '* Select Exp. Case...'
      read (*, *) ExpCaseNum
      write (fname, '(i2.2)') ExpCaseNum
      ! 実験条件入力
      call Input_ExpCondition( &
      &      trim(ExpConditionFile)//trim(fname)//strtxt)
      ! 着氷モデル選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write (*, '(a)') '* Select Icing Model...'
      read (*, *) ThermoNum
      write (fname, '(i2.2)') ThermoNum
      ! ディレクトリ設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write (fname, '(a)') 'clean/'
      GrdOutDir = bckdir//'grid/'//trim(fname)
      ! 計算条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 計算条件設定
      nCount = 0; nDrop = 0; IceStep = 0; IceStepMax = 5
      IceTime = 0.0; dti = 1.0e-3
      Cn = 1.0; nRunge = 4; nCalMax = 20000; nOutputLog = 1000; nOutputCount = 50000
      nTVD = 2; eTVD = 0.1; LmtPro = 0.1; LmtAve = 100.0
      TurbNum = 4; fSteady = .true.; fTime = .false.
      BCNum = 2; fSlip = .false.
      RoughNum = 1; RunbackNum = 1; DropShedNum = 0; IceShedNum = 0; DragNum = 1
      Cnd = 0.03
      nDrpIn = 1000000; nDrpCalMax = 100000; nDrpFile = 100; nDrpOutputLog = 100000; nDrpOutputCount = 50000
      Rg = 287.1; gamma = 1.4; Pr = 0.72; Prt = 0.90; Ret = 500.0; Cmu = 0.09
      muSth = 1.82e-5; TsSth = 273.15; s1 = 117.0
      ! 計算条件ファイル出力
      call Output_CalSetting(trim(ND_CalSetFile)//strtxt)
      ! 処理終了 ********************************************************************************************
      return
   end subroutine SelectExpCase
!*******************************************************************************************************
!********* 初期設定                                                                                ********
!*******************************************************************************************************
   subroutine InitialSetting
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 処理開始 ********************************************************************************************
      ! ブロック名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call Set_BlockName
      ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      allocate (Grd(ms:me))
      ! 翼形状 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ycmax = chord*real(NACA1)*1.0e-2
      xc = chord*real(NACA2)*1.0e-1
      thick = chord*real(NACA34)*1.0e-2
      ytmax = 0.5*thick
      ! 処理終了 ********************************************************************************************
      return
   end subroutine InitialSetting
!*******************************************************************************************************
!********* Main Grid 生成                                                                        ********
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
      Grd(m)%is = 0
      Grd(m)%i1 = 192 + Grd(m)%is
      Grd(m)%i2 = 256 + Grd(m)%i1
      Grd(m)%i3 = 256 + Grd(m)%i2
      Grd(m)%ie = 192 + Grd(m)%i3
      Grd(m)%js = 0
      Grd(m)%je = 256
      Grd(m)%ks = 0
      Grd(m)%ke = 3
      ! 計算領域 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      dom1 = 10.0*chord                                                        ! 後縁から流出境界まで
      dom2 = 10.0*chord                                                        ! 後縁から外部境界まで
      dom3 = 10.0*chord                                                        ! 最大翼厚位置から流入境界まで
      ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      allocate (Grd(m)%x(Grd(m)%is:Grd(m)%ie, Grd(m)%js:Grd(m)%je, Grd(m)%ks:Grd(m)%ke), &
      &         Grd(m)%y(Grd(m)%is:Grd(m)%ie, Grd(m)%js:Grd(m)%je, Grd(m)%ks:Grd(m)%ke), &
      &         Grd(m)%z(Grd(m)%is:Grd(m)%ie, Grd(m)%js:Grd(m)%je, Grd(m)%ks:Grd(m)%ke), &
      &         Grd(m)%f(Grd(m)%is:Grd(m)%ie, Grd(m)%js:Grd(m)%je, Grd(m)%ks:Grd(m)%ke))
      ! 格子生成 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write (*, '(a)') '+++++ C-typr grid generation start. +++++'
      call CtypeGridBlade( &
      &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
      &      Grd(m)%i1, Grd(m)%i2, Grd(m)%i3, dom1, dom2, dom3, &
      &      Grd(m)%x, Grd(m)%y, Grd(m)%z)
      ! 翼のフラグ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call ViewBlade( &
      &      m, Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
      &      Grd(m)%i1, Grd(m)%i3, Grd(m)%js, Grd(m)%f)
      ! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call Output_Resolution3D( &
      &      trim(GrdOutDir)//trim(BlkName(m))//trim(RslFile), strtxt, &
      &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke)
      call Output_Resolution1D( &
      &      trim(GrdOutDir)//trim(BlkName(m))//trim(IceRslFile), strtxt, &
      &      Grd(m)%i1, Grd(m)%i3)
      call Output_CtypeGridPoint( &
      &      trim(GrdOutDir)//trim(BlkName(m))//trim(CtypePointFile)//strtxt, &
      &      Grd(m)%i1, Grd(m)%i2, Grd(m)%i3)
      call Output_Grid3D( &
      &      trim(GrdOutDir)//trim(BlkName(m))//trim(GrdFile), strbin, &
      &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
      &      Grd(m)%x, Grd(m)%y, Grd(m)%z)
      ! 可視化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call MakeMAVSFile3D( &
      &      trim(GrdOutDir), trim(BlkName(m))//trim(ViewGrdFile), strbin, &
      &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
      &      Grd(m)%f, Grd(m)%x, Grd(m)%y, Grd(m)%z)
      call OutputPara_bin( &
       &      trim(GrdOutDir), trim(BlkName(m))//trim(ViewGrdFile), &
       &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
       &      Grd(m)%x, Grd(m)%y, Grd(m)%z)
      ! 処理終了 ********************************************************************************************
      return
   end subroutine MainGrid
!*******************************************************************************************************
!********* Sub Grid 生成                                                                        ********
!*******************************************************************************************************
   subroutine SubGrid
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer :: m
      real    :: dom1, dom2, dom3
      ! 処理開始 ********************************************************************************************
      ! ブロック番号 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      m = me
      ! 解像度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Grd(m)%is = 0
      Grd(m)%i1 = 192 + Grd(m)%is
      Grd(m)%i2 = 256 + Grd(m)%i1
      Grd(m)%i3 = 256 + Grd(m)%i2
      Grd(m)%ie = 192 + Grd(m)%i3
      Grd(m)%js = 0
      Grd(m)%je = 256
      Grd(m)%ks = 0
      Grd(m)%ke = 4
      ! 計算領域 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      dom1 = 10.0*chord                                                        ! 後縁から流出境界まで
      dom2 = 10.0*chord                                                        ! 後縁から外部境界まで
      dom3 = 10.0*chord                                                        ! 最大翼厚位置から流入境界まで
      ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      allocate (Grd(m)%x(Grd(m)%is:Grd(m)%ie, Grd(m)%js:Grd(m)%je, Grd(m)%ks:Grd(m)%ke), &
      &         Grd(m)%y(Grd(m)%is:Grd(m)%ie, Grd(m)%js:Grd(m)%je, Grd(m)%ks:Grd(m)%ke), &
      &         Grd(m)%z(Grd(m)%is:Grd(m)%ie, Grd(m)%js:Grd(m)%je, Grd(m)%ks:Grd(m)%ke), &
      &         Grd(m)%f(Grd(m)%is:Grd(m)%ie, Grd(m)%js:Grd(m)%je, Grd(m)%ks:Grd(m)%ke))
      ! 格子生成 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write (*, '(a)') '+++++ C-typr grid generation start. +++++'
      call CtypeGridBlade( &
      &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
      &      Grd(m)%i1, Grd(m)%i2, Grd(m)%i3, dom1, dom2, dom3, &
      &      Grd(m)%x, Grd(m)%y, Grd(m)%z)
      ! 翼のフラグ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call ViewBlade( &
      &      m, Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
      &      Grd(m)%i1, Grd(m)%i3, Grd(m)%js, Grd(m)%f)
      ! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call Output_Resolution3D( &
      &      trim(GrdOutDir)//trim(BlkName(m))//trim(RslFile), strtxt, &
      &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke)
      call Output_Resolution1D( &
      &      trim(GrdOutDir)//trim(BlkName(m))//trim(IceRslFile), strtxt, &
      &      Grd(m)%i1, Grd(m)%i3)
      call Output_CtypeGridPoint( &
      &      trim(GrdOutDir)//trim(BlkName(m))//trim(CtypePointFile)//strtxt, &
      &      Grd(m)%i1, Grd(m)%i2, Grd(m)%i3)
      call Output_Grid3D( &
      &      trim(GrdOutDir)//trim(BlkName(m))//trim(GrdFile), strbin, &
      &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
      &      Grd(m)%x, Grd(m)%y, Grd(m)%z)
      ! 可視化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call MakeMAVSFile3D( &
      &      trim(GrdOutDir), trim(BlkName(m))//trim(ViewGrdFile), strbin, &
      &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
      &      Grd(m)%f, Grd(m)%x, Grd(m)%y, Grd(m)%z)
      call OutputPara_bin( &
       &      trim(GrdOutDir), trim(BlkName(m))//trim(ViewGrdFile), &
       &      Grd(m)%is, Grd(m)%ie, Grd(m)%js, Grd(m)%je, Grd(m)%ks, Grd(m)%ke, &
       &      Grd(m)%x, Grd(m)%y, Grd(m)%z)
      ! 処理終了 ********************************************************************************************
      return
   end subroutine SubGrid
!*******************************************************************************************************
!********* 翼周りの格子 (C-type)                                                                ********
!*******************************************************************************************************
   subroutine CtypeGridBlade( &
   &            is, ie, js, je, ks, ke, i1, i2, i3, dom1, dom2, dom3, x, y, z)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in)  :: is, ie, js, je, ks, ke
      integer, intent(in)  :: i1, i2, i3
      real, intent(in)  :: dom1, dom2, dom3
      real, intent(out) :: x(is:ie, js:je, ks:ke), &
      &                       y(is:ie, js:je, ks:ke), &
      &                       z(is:ie, js:je, ks:ke)
      ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, parameter :: rs1 = 1.0*1.0e-1                                ! 領域ⅰ起点格子幅
      real, parameter :: re1 = 4.79512133220952*1.0e-1                        ! 領域ⅰ終点格子幅
      real, parameter :: rs2 = 6.0*1.0e-0                                ! 領域ⅱ起点格子幅
      real, parameter :: re2 = 5.0*1.0e-2                               ! 領域ⅱ終点格子幅
      real, parameter :: rs3 = 5.0*1.0e-2                                ! 領域ⅲ起点格子幅
      real, parameter :: re3 = 6.0*1.0e-0                               ! 領域ⅲ終点格子幅
      real, parameter :: rs4 = 4.79512133220952*1.0e-1                        ! 領域ⅳ起点格子幅
      real, parameter :: re4 = 1.0*1.0e-1                               ! 領域ⅳ終点格子幅
      real, parameter :: rs5 = 3.28649582292613*1.0e-3                        ! 内部境界格子幅
      real, parameter :: re5 = 1.0*1.0e+2                                ! 外部境界格子幅
      real, parameter :: tb1 = 2.6579*1.0e-1                                ! 内部境界直交性のパラメータ
      real, parameter :: tb2 = 7.0*1.0e0                                ! 外部境界直交性のパラメータ
      real, parameter :: mgn = 3.0                                        ! 楕円-双曲型重み関数の許容誤差
      real, parameter :: rsd = 1.0e-5                                        ! 楕円-双曲型最大計算回数
      ! 処理開始 ********************************************************************************************
      ! 内部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call CtypeInternalBoundary( &
      &      is, i1, i2, i3, ie, dom1, rs1, re1, rs2, re2, rs3, re3, rs4, re4, x(:, js, ks), y(:, js, ks))
      ! 外部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call CtypeExternalBoundary( &
      &      is, i1, i2, i3, ie, js, je, dom2, dom3, x(:, :, ks), y(:, :, ks))
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
      &      is, ie, js, je, rs5, re5, tb1, tb2, x(:, :, ks), y(:, :, ks))
      ! 三次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call ThreeDimensionalized( &
      &      is, ie, js, je, ks, ke, span, x, y, z)
      ! 処理終了 ********************************************************************************************
      return
   end subroutine CtypeGridBlade
!*******************************************************************************************************
!******** 内部境界 (C-type)                                                                        ********
!*******************************************************************************************************
   subroutine CtypeInternalBoundary( &
   &            is, i1, i2, i3, ie, dom, rs1, re1, rs2, re2, rs3, re3, rs4, re4, x, y)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in)  :: is, i1, i2, i3, ie
      real, intent(in)  :: dom
      real, intent(in)  :: rs1, re1, rs2, re2, rs3, re3, rs4, re4
      real, intent(out) :: x(is:ie), y(is:ie)
      ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, parameter :: nmax = 10000
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, pointer :: xu(:), xl(:), yu(:), yl(:)
      real, pointer :: tu(:), tl(:), lu(:), ll(:)
      real, pointer :: x2(:), y2(:)
      real, pointer :: tb(:)
      real, pointer :: sy(:)
      integer :: n, i
      integer :: n1, n2, n3, n4, n5
      integer :: si
      real    :: xx, yy, yt, yc, theta
      real    :: rr
      ! 処理開始 ********************************************************************************************
      ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      allocate (xu(0:nmax), yu(0:nmax), xl(0:nmax), yl(0:nmax))
      allocate (tu(i2:i3), tl(i1:i2), lu(0:nmax), ll(0:nmax))
      allocate (x2(i1:i2), y2(i1:i2))
      allocate (tb(is:i1))
      ! 翼周り ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 翼形状 ----------------------------------------------------------------------------------------------
      do n = 1, nmax - 1
         ! コード方向位置
         xx = chord*real(n)/real(nmax)
         ! 翼厚分布
         yt = 5.0*thick*(a0*(xx/chord)**0.5 + a1*(xx/chord)**1 + a2*(xx/chord)**2 &
         &                  + a3*(xx/chord)**3 + a4*(xx/chord)**4)
         ! 翼反り分布
         if (xc == 0.0) then
            yc = 0.0
         else
            if (xx < xc) then
               yc = ycmax*chord**2/xc**2*(2.0*xc*xx/chord**2 - xx**2/chord**2)
            else
               yc = ycmax*chord**2/(chord - xc)**2 &
                   &  *(1.0 - 2.0*xc/chord + 2.0*xc*xx/chord**2 - xx**2/chord**2)
            end if
         end if
         ! 翼座標
         theta = atan(yc/xx)
         xu(n) = xx - yt*sin(theta)
         yu(n) = yc + yt*cos(theta)
         xl(n) = xx + yt*sin(theta)
         yl(n) = yc - yt*cos(theta)
      end do
      xu(0) = 0.0; xu(nmax) = chord
      yu(0) = 0.0; yu(nmax) = 0.0
      xl(0) = 0.0; xl(nmax) = chord
      yl(0) = 0.0; yl(nmax) = 0.0
      ! 格子幅調整 ------------------------------------------------------------------------------------------
      ! 媒介変数
      call VinokurInterpolation(rs3/real(i3 - i2 + 1), re3/real(i3 - i2 + 1), i3 - i2 + 1, tu)
      call VinokurInterpolation(re2/real(i2 - i1 + 1), rs2/real(i2 - i1 + 1), i2 - i1 + 1, tl)
      ! 前縁からの距離
      lu(0) = 0.0; ll(0) = 0.0
      do n = 1, nmax
         lu(n) = sqrt((xu(n) - xu(n - 1))**2 + (yu(n) - yu(n - 1))**2) + lu(n - 1)
         ll(n) = sqrt((xl(n) - xl(n - 1))**2 + (yl(n) - yl(n - 1))**2) + ll(n - 1)
      end do
      ! 無次元化
      do n = 0, nmax
         lu(n) = lu(n)/lu(nmax)
         ll(n) = ll(n)/ll(nmax)
      end do
      ! ラグランジュ補間で探索 (上面)
      do i = i2 + 1, i3 - 1
         do n = 0, nmax - 1
            if (lu(n) <= tu(i) .and. tu(i) < lu(n + 1)) then
               select case (n)
               case default
                  n1 = n - 2; n2 = n - 1; n3 = n; n4 = n + 1; n5 = n + 2
                  call LagrangeInterpolation(5, lu(n1:n5), xu(n1:n5), tu(i), x(i))
                  call LagrangeInterpolation(5, lu(n1:n5), yu(n1:n5), tu(i), y(i))
               case (0)
                  n1 = n; n2 = n + 1; n3 = n + 2
                  call LagrangeInterpolation(3, lu(n1:n3), xu(n1:n3), tu(i), x(i))
                  call LagrangeInterpolation(3, lu(n1:n3), yu(n1:n3), tu(i), y(i))
               case (1)
                  n1 = n - 1; n2 = n; n3 = n + 1; n4 = n + 2
                  call LagrangeInterpolation(4, lu(n1:n4), xu(n1:n4), tu(i), x(i))
                  call LagrangeInterpolation(4, lu(n1:n4), yu(n1:n4), tu(i), y(i))
               case (nmax - 1)
                  n1 = n - 2; n2 = n - 1; n3 = n; n4 = n + 1
                  call LagrangeInterpolation(4, lu(n1:n4), xu(n1:n4), tu(i), x(i))
                  call LagrangeInterpolation(4, lu(n1:n4), yu(n1:n4), tu(i), y(i))
               end select
               exit
            end if
         end do
      end do
      x(i2) = xu(0); x(i3) = xu(nmax)
      y(i2) = yu(0); y(i3) = yu(nmax)
      ! ラグランジュ補間で探索 (下面)
      do i = i1 + 1, i2 - 1
         do n = 0, nmax - 1
            if (ll(n) <= tl(i) .and. tl(i) < ll(n + 1)) then
               select case (n)
               case default
                  n1 = n - 2; n2 = n - 1; n3 = n; n4 = n + 1; n5 = n + 2
                  call LagrangeInterpolation(5, ll(n1:n5), xl(n1:n5), tl(i), x(i))
                  call LagrangeInterpolation(5, ll(n1:n5), yl(n1:n5), tl(i), y(i))
               case (0)
                  n1 = n; n2 = n + 1; n3 = n + 2
                  call LagrangeInterpolation(3, ll(n1:n3), xl(n1:n3), tl(i), x(i))
                  call LagrangeInterpolation(3, ll(n1:n3), yl(n1:n3), tl(i), y(i))
               case (1)
                  n1 = n - 1; n2 = n; n3 = n + 1; n4 = n + 2
                  call LagrangeInterpolation(4, ll(n1:n4), xl(n1:n4), tl(i), x(i))
                  call LagrangeInterpolation(4, ll(n1:n4), yl(n1:n4), tl(i), y(i))
               case (nmax - 1)
                  n1 = n - 2; n2 = n - 1; n3 = n; n4 = n + 1
                  call LagrangeInterpolation(4, ll(n1:n4), xl(n1:n4), tl(i), x(i))
                  call LagrangeInterpolation(4, ll(n1:n4), yl(n1:n4), tl(i), y(i))
               end select
               exit
            end if
         end do
      end do
      x(i1) = xl(0); x(i2) = xl(nmax)
      y(i1) = yl(0); y(i2) = yl(nmax)
      ! 格子番号整理
      do i = i1, i2
         x2(i) = x(i1 + i2 - i)
         y2(i) = y(i1 + i2 - i)
      end do
      do i = i1, i2
         x(i) = x2(i)
         y(i) = y2(i)
      end do
      ! 後縁をスムージング ----------------------------------------------------------------------------------
      allocate (sy(is:ie))
      si = 3
      do n = 1, 5
         do i = i1 + 1, i1 + si
            sy(i) = (y(i - 1) + y(i) + y(i + 1))/3.0
         end do
         do i = i1 + 1, i1 + si
            y(i) = sy(i)
         end do
         do i = i3 - si, i3 - 1
            sy(i) = (y(i - 1) + y(i) + y(i + 1))/3.0
         end do
         do i = i3 - si, i3 - 1
            y(i) = sy(i)
         end do
      end do
      deallocate (sy)
      ! C 型格子ブランチ・カット ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call GeometricInterpolationInv(re1/real(i1 - is + 1), i1 - is + 1, tb, rr)
      do i = is, i1 - 1
         x(i) = dom - (dom - x(i1))*tb(i)
         y(i) = y(i1)
      end do
      do i = i3, ie
         x(i) = x(ie - i)
         y(i) = y(ie - i)
      end do
      ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      deallocate (xu, yu, tu, lu, xl, yl, tl, ll, x2, y2, tb)
      ! 処理終了 ********************************************************************************************
      return
   end subroutine CtypeInternalBoundary
!*******************************************************************************************************
!******** 外部境界 (C-type)                                                                        ********
!*******************************************************************************************************
   subroutine CtypeExternalBoundary( &
   &            is, i1, i2, i3, ie, js, je, dom1, dom2, x, y)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in)    :: is, i1, i2, i3, ie, js, je
      real, intent(in)    :: dom1, dom2
      real, intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer :: i
      integer :: ic
      real    :: a1, a2, a3, b1, b2, b3, c1, c2, c3
      ! 処理開始 ********************************************************************************************
      ! C 型の湾曲部 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do i = i1 + 1, i3 - 1
         ! 翼周り方向ベクトル
         a1 = 0.5*(-x(i - 1, js) + x(i + 1, js))
         a2 = 0.5*(-y(i - 1, js) + y(i + 1, js))
         a3 = 0.0
         ! スパン方向ベクトル
         b1 = 0.0
         b2 = 0.0
         b3 = -1.0
         ! 法線方向ベクトル
         c1 = a2*b3 - a3*b2
         c2 = a3*b1 - a1*b3
         c3 = a1*b2 - a2*b1
         ! 翼周り部外部領域
         x(i, je) = c1/sqrt(c1**2 + c2**2 + c3**2)*dom1 + x(i, js)
         y(i, je) = c2/sqrt(c1**2 + c2**2 + c3**2)*dom1 + y(i, js)
      end do
      ! C 型の直線部 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! x 方向 ----------------------------------------------------------------------------------------------
      ! 下面
      do i = i1 + 1, i2
         if (x(i, je) < x(i1, js)) exit
      end do
      ic = i - 1
      do i = is, ic
         x(i, je) = x(ic, je) + (x(is, js) - x(ic, je))*real(ic - i)/real(ic)
      end do
      ! 上面
      do i = i2, i3 - 1
         if (x(i, je) > x(i1, js)) exit
      end do
      ic = i
      do i = ic, ie
         x(i, je) = x(ic, je) + (x(is, js) - x(ic, je))*real(i - ic)/real(ie - ic)
      end do
      ! y 方向 ----------------------------------------------------------------------------------------------
      ! 下面
      do i = i1 + 1, i2
         if (y(i, je) == minval(y(i1 + 1:i2, je))) exit
      end do
      ic = i
      do i = is, ic - 1
         y(i, je) = minval(y(ic:i2, je))
      end do
      ! 上面
      do i = i2, i3 - 1
         if (y(i, je) == maxval(y(i2:i3 - 1, je))) exit
      end do
      ic = i
      do i = ic + 1, ie
         y(i, je) = maxval(y(i2:ic, je))
      end do
      ! 処理終了 ********************************************************************************************
      return
   end subroutine CtypeExternalBoundary
!*******************************************************************************************************
!******** 側部境界 (C-type)                                                                        ********
!*******************************************************************************************************
   subroutine CtypeSideBoundary( &
   &            is, ie, js, je, x, y)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in)    :: is, ie, js, je
      real, intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer :: j
      ! 処理開始 ********************************************************************************************
      do j = js + 1, je - 1
         x(is, j) = x(is, js)
         x(ie, j) = x(ie, js)
         y(is, j) = (y(is, je) - y(is, js))*real(j)/real(je) + y(is, js)
         y(ie, j) = (y(ie, je) - y(ie, js))*real(j)/real(je) + y(ie, js)
      end do
      ! 処理開始 ********************************************************************************************
      return
   end subroutine CtypeSideBoundary
!*******************************************************************************************************
!******** 楕円-双曲型偏微分方程式法に基づく格子生成                                                ********
!*******************************************************************************************************
   subroutine CtypeGenerationEHPDE( &
   &            is, i1, i2, i3, ie, js, je, rs, re, margin, resi, x, y)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in)    :: is, i1, i2, i3, ie, js, je
      real, intent(in)    :: rs, re
      real, intent(in)    :: margin                                        ! 楕円-双曲型重み関数許容誤差
      real, intent(in)    :: resi                                                ! 計算ループ収束判定値
      real, intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
      ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, parameter :: omg = 1.0                                        ! 緩和係数
      integer, parameter :: nmax = 100000                                        ! 最大計算回数
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, pointer :: Cs(:, :)
      real, pointer :: cv1m(:, :), cv1p(:, :), cv2m(:, :), cv2p(:, :)
      real, pointer :: dx1m(:, :), dx1p(:, :), dy1m(:, :), dy1p(:, :)
      real, pointer :: dx2m(:, :), dx2p(:, :), dy2m(:, :), dy2p(:, :)
      integer :: i, j, n
      real    :: dmax
      real    :: a1, a2, a3, b1, b2, b3, c1, c2, c3, hjs1, hjs2, hje1, hje2
      ! 処理開始 ********************************************************************************************
      ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      allocate (cs(is:ie, js:je), &
      &         cv1m(is:ie, js:je), cv1p(is:ie, js:je), cv2m(is:ie, js:je), cv2p(is:ie, js:je), &
      &         dx1m(is:ie, js:je), dx1p(is:ie, js:je), dy1m(is:ie, js:je), dy1p(is:ie, js:je), &
      &         dx2m(is:ie, js:je), dx2p(is:ie, js:je), dy2m(is:ie, js:je), dy2p(is:ie, js:je))
      ! パラメータ設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 双曲型重み関数 --------------------------------------------------------------------------------------
      do j = js + 1, je - 1
      do i = is + 1, ie - 1
         ! i - 方向
         if (i < i1) then
            cv1m(i, j) = 1.0
         else if (i <= i3) then
            cv1m(i, j) = real(j)/real(je)*real(i3 - i)/real(i3 - i1)
         else
            cv1m(i, j) = 0.0
         end if
         ! i + 方向
         if (i3 < i) then
            cv1p(i, j) = 1.0
         else if (i1 <= i) then
            cv1p(i, j) = real(j)/real(je)*real(i - i1)/real(i3 - i1)
         else
            cv1p(i, j) = 0.0
         end if
         ! j - 方向
         cv2m(i, j) = (real(je - j)/real(je))**2
         ! j + 方向
         if (i2 < i) then
            cv2p(i, j) = (real(ie - i)/real(ie))**2
         else if (i == i2) then
            cv2p(i, j) = 0.0
         else
            cv2p(i, j) = (real(i)/real(ie))**2
         end if
!  cv2p(i,j) = (real(j-js) / real(je))**1 * 0.3
      end do
      end do
      ! 楕円型重み関数 --------------------------------------------------------------------------------------
      do j = js + 1, je - 1
      do i = is + 1, ie - 1
         ! 内部境界法線ベクトル
         a1 = 0.5*(-x(i - 1, js) + x(i + 1, js))
         a2 = 0.5*(-y(i - 1, js) + y(i + 1, js))
         a3 = 0.0
         b1 = 0.0
         b2 = 0.0
         b3 = -1.0
         c1 = a2*b3 - a3*b2
         c2 = a3*b1 - a1*b3
         hjs1 = c1/sqrt(c1**2 + c2**2)
         hjs2 = c2/sqrt(c1**2 + c2**2)
         ! 外部境界法線ベクトル
         a1 = 0.5*(-x(i - 1, je) + x(i + 1, je))
         a2 = 0.5*(-y(i - 1, je) + y(i + 1, je))
         a3 = 0.0
         b1 = 0.0
         b2 = 0.0
         b3 = -1.0
         c1 = a2*b3 - a3*b2
         c2 = a3*b1 - a1*b3
         hje1 = c1/sqrt(c1**2 + c2**2)
         hje2 = c2/sqrt(c1**2 + c2**2)
         ! i - 方向
         if (i < i1) then
            dx1m(i, j) = x(i, je) - x(i - 1, je)
            dy1m(i, j) = y(i, je) - y(i - 1, je)
         else if (i <= i3) then
            dx1m(i, j) = (-x(i - 1, js) + x(i, js))*real(je - j)/real(je) &
            &         + (-x(i - 1, je) + x(i, je))*(real(j - js)/real(je))**2
            dy1m(i, j) = (-y(i - 1, js) + y(i, js))*real(je - j)/real(je) &
            &         + (-y(i - 1, je) + y(i, je))*(real(j - js)/real(je))**2
         else
            dx1m(i, j) = 0.0
            dy1m(i, j) = 0.0
         end if
         ! i + 方向
         if (i3 < i) then
            dx1p(i, j) = x(i + 1, je) - x(i, je)
            dy1p(i, j) = y(i + 1, je) - y(i, je)
         else if (i1 <= i) then
            dx1p(i, j) = (-x(i, js) + x(i + 1, js))*real(je - j)/real(je) &
            &         + (-x(i, je) + x(i + 1, je))*(real(j - js)/real(je))**2
            dy1p(i, j) = (-y(i, js) + y(i + 1, js))*real(je - j)/real(je) &
            &         + (-y(i, je) + y(i + 1, je))*(real(j - js)/real(je))**2
         else
            dx1p(i, j) = 0.0
            dy1p(i, j) = 0.0
         end if
         ! j - 方向
         dx2m(i, j) = hjs1*rs
         dy2m(i, j) = hjs2*rs
         ! j + 方向
         dx2p(i, j) = hje1*re
         dy2p(i, j) = hje2*re
      end do
      end do
      ! 楕円-双曲型重み関数 ---------------------------------------------------------------------------------
      do j = js + 1, je - 1
      do i = is + 1, ie - 1
         Cs(i, j) = 1.0 - max(cv1m(i, j), cv1p(i, j), cv2m(i, j), cv2p(i, j)) + margin
      end do
      end do
      ! 格子生成計算ループ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do n = 1, nmax
         ! ソルバー部 -----------------------------------------------------------------------------------------
         call GridEHPDE2D( &
         &      omg, is, ie, js, je, &
         &      cs, cv1m, cv1p, cv2m, cv2p, dx1m, dy1m, dx1p, dy1p, dx2m, dy2m, dx2p, dy2p, &
         &      x, y, dmax)
         if (mod(n, 1000) == 0.0) write (*, "(a,i5,e16.8e3)") "* Elliptic-Hyperbolic calculation...", n, dmax
         if (dmax < resi) exit
         ! 境界条件 -------------------------------------------------------------------------------------------
         do j = js + 1, je - 1
            y(is, j) = y(is + 1, j)
            y(ie, j) = y(ie - 1, j)
         end do
      end do
      if (dmax > resi) then
         write (*, '(a)') "!!!!! Elliptic-Hyperbolic calculation error !!!!!"
         stop
      end if
      ! 処理終了 ********************************************************************************************
      return
   end subroutine CtypeGenerationEHPDE
!*******************************************************************************************************
!********* 翼前縁付近の着氷計算用の格子 (H-type)                                                ********
!*******************************************************************************************************
   subroutine HtypeGridIceLE( &
   &            is, ie, js, je, ks, ke, dom1, dom2, x, y, z)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in)  :: is, ie, js, je, ks, ke
      real, intent(in)  :: dom1, dom2
      real, intent(out) :: x(is:ie, js:je, ks:ke), &
      &                       y(is:ie, js:je, ks:ke), &
      &                       z(is:ie, js:je, ks:ke)
      ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, parameter :: rs1 = 7.0*1.0e-2                                ! 前縁付近格子幅
      real, parameter :: re1 = 2.0*1.0e-0                                ! 着氷限界位置付近格子幅
      real, parameter :: rs2 = 5.0*1.0e-1                                ! 翼表面境界格子幅
      real, parameter :: re2 = 7.0*1.0e-2                                ! 翼遠方境界格子幅
      real, parameter :: tb1 = 1.0*1.0e-1                                ! 内部境界直交性のパラメータ
      real, parameter :: tb2 = 2.0*1.0e-1                                ! 外部境界直交性のパラメータ
      ! 処理開始 ********************************************************************************************
      ! 内部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call HtypeInternalBoundary( &
      &      is, ie, dom1, rs1, re1, x(:, js, ks), y(:, js, ks))
      ! 外部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call HtypeExternalBoundary( &
      &      is, ie, js, je, dom2, x(:, :, ks), y(:, :, ks))
! ! 側部境界 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call HtypeSideBoundary( &
! &      is, ie, js, je, rs2, re2, x(:, :, ks), y(:, :, ks) )
! ! Transfinite 補間 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call TransfiniteInterpolation( &
! &      is, ie, js, je, x(:, :, ks), y(:, :, ks) )
      ! 二境界法 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call GenerationTwoBoundary( &
      &      is, ie, js, je, rs2, re2, tb1, tb2, x(:, :, ks), y(:, :, ks))
      ! 三次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call ThreeDimensionalized( &
      &      is, ie, js, je, ks, ke, span, x, y, z)
      ! 処理終了 ********************************************************************************************
      return
   end subroutine HtypeGridIceLE
!*******************************************************************************************************
!******** 内部境界 (H-type)                                                                        ********
!*******************************************************************************************************
   subroutine HtypeInternalBoundary( &
   &            is, ie, dom, rs, re, x, y)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in)  :: is, ie
      real, intent(in)  :: dom
      real, intent(in)  :: rs, re
      real, intent(out) :: x(is:ie), y(is:ie)
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, pointer :: xu(:), xl(:), yu(:), yl(:)
      real, pointer :: tu(:), tl(:), lu(:), ll(:)
      real, pointer :: x0(:), y0(:)
      real, pointer :: sy(:)
      integer :: n, i
      integer :: n1, n2, n3, n4, n5
      integer :: si
      integer :: i1
      integer :: nmax
      real    :: xx, yy, yt, yc, theta
      real    :: rr
      ! 処理開始 ********************************************************************************************
      ! 前縁格子番号 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      i1 = int(0.5*ie)
      ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      nmax = 10000
      allocate (xu(0:nmax), yu(0:nmax), xl(0:nmax), yl(0:nmax))
      allocate (tu(i1:ie), tl(is:i1), lu(0:nmax), ll(0:nmax))
      allocate (x0(is:i1), y0(is:i1))
      ! 翼周り ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 翼形状 ----------------------------------------------------------------------------------------------
      do n = 1, nmax - 1
         ! コード方向位置
         xx = chord*real(n)/real(nmax)
         ! 翼厚分布
         yt = 5.0*thick*(a0*(xx/chord)**0.5 + a1*(xx/chord)**1 + a2*(xx/chord)**2 &
         &                  + a3*(xx/chord)**3 + a4*(xx/chord)**4)
         ! 翼反り分布
         if (xc == 0.0) then
            yc = 0.0
         else
            if (xx < xc) then
               yc = ycmax*chord**2/xc**2*(2.0*xc*xx/chord**2 - xx**2/chord**2)
            else
               yc = ycmax*chord**2/(chord - xc)**2 &
                   &  *(1.0 - 2.0*xc/chord + 2.0*xc*xx/chord**2 - xx**2/chord**2)
            end if
         end if
         ! 翼座標
         theta = atan(yc/xx)
         xu(n) = xx - yt*sin(theta)
         yu(n) = yc + yt*cos(theta)
         xl(n) = xx + yt*sin(theta)
         yl(n) = yc - yt*cos(theta)
      end do
      xu(0) = 0.0; xu(nmax) = chord
      yu(0) = 0.0; yu(nmax) = 0.0
      xl(0) = 0.0; xl(nmax) = chord
      yl(0) = 0.0; yl(nmax) = 0.0
      ! 計算領域 --------------------------------------------------------------------------------------------
      do n = 0, nmax
         if (xu(n) - xu(0) > dom) exit
      end do
      nmax = n
      ! 格子幅調整 ------------------------------------------------------------------------------------------
      ! 媒介変数
      call VinokurInterpolation(rs/real(ie - i1 + 1), re/real(ie - i1 + 1), ie - i1 + 1, tu)
      call VinokurInterpolation(rs/real(i1 - is + 1), re/real(i1 - is + 1), i1 - is + 1, tl)
      ! 前縁からの距離
      lu(0) = 0.0; ll(0) = 0.0
      do n = 1, nmax
         lu(n) = sqrt((xu(n) - xu(n - 1))**2 + (yu(n) - yu(n - 1))**2) + lu(n - 1)
         ll(n) = sqrt((xl(n) - xl(n - 1))**2 + (yl(n) - yl(n - 1))**2) + ll(n - 1)
      end do
      ! 無次元化
      do n = 0, nmax
         lu(n) = lu(n)/lu(nmax)
         ll(n) = ll(n)/ll(nmax)
      end do
      ! ラグランジュ補間で探索 (上面)
      do i = i1 + 1, ie - 1
         do n = 0, nmax - 1
            if (lu(n) <= tu(i) .and. tu(i) < lu(n + 1)) then
               if (n <= 0) then
                  n1 = n; n2 = n + 1; n3 = n + 2
                  call LagrangeInterpolation(3, lu(n1:n3), xu(n1:n3), tu(i), x(i))
                  call LagrangeInterpolation(3, lu(n1:n3), yu(n1:n3), tu(i), y(i))
               else if (n <= 1) then
                  n1 = n - 1; n2 = n; n3 = n + 1; n4 = n + 2
                  call LagrangeInterpolation(4, lu(n1:n4), xu(n1:n4), tu(i), x(i))
                  call LagrangeInterpolation(4, lu(n1:n4), yu(n1:n4), tu(i), y(i))
               else if (n < nmax - 1) then
                  n1 = n - 2; n2 = n - 1; n3 = n; n4 = n + 1; n5 = n + 2
                  call LagrangeInterpolation(5, lu(n1:n5), xu(n1:n5), tu(i), x(i))
                  call LagrangeInterpolation(5, lu(n1:n5), yu(n1:n5), tu(i), y(i))
               else
                  n1 = n - 2; n2 = n - 1; n3 = n; n4 = n + 1
                  call LagrangeInterpolation(4, lu(n1:n4), xu(n1:n4), tu(i), x(i))
                  call LagrangeInterpolation(4, lu(n1:n4), yu(n1:n4), tu(i), y(i))
               end if
               exit
            end if
         end do
      end do
      x(i1) = xu(0); x(ie) = xu(nmax)
      y(i1) = yu(0); y(ie) = yu(nmax)
      ! ラグランジュ補間で探索 (下面)
      do i = is + 1, i1 - 1
         do n = 0, nmax - 1
            if (ll(n) <= tl(i) .and. tl(i) < ll(n + 1)) then
               if (n <= 0) then
                  n1 = n; n2 = n + 1; n3 = n + 2
                  call LagrangeInterpolation(3, ll(n1:n3), xl(n1:n3), tl(i), x(i))
                  call LagrangeInterpolation(3, ll(n1:n3), yl(n1:n3), tl(i), y(i))
               else if (n <= 1) then
                  n1 = n - 1; n2 = n; n3 = n + 1; n4 = n + 2
                  call LagrangeInterpolation(4, ll(n1:n4), xl(n1:n4), tl(i), x(i))
                  call LagrangeInterpolation(4, ll(n1:n4), yl(n1:n4), tl(i), y(i))
               else if (n < nmax - 1) then
                  n1 = n - 2; n2 = n - 1; n3 = n; n4 = n + 1; n5 = n + 2
                  call LagrangeInterpolation(5, ll(n1:n5), xl(n1:n5), tl(i), x(i))
                  call LagrangeInterpolation(5, ll(n1:n5), yl(n1:n5), tl(i), y(i))
               else
                  n1 = n - 2; n2 = n - 1; n3 = n; n4 = n + 1
                  call LagrangeInterpolation(4, ll(n1:n4), xl(n1:n4), tl(i), x(i))
                  call LagrangeInterpolation(4, ll(n1:n4), yl(n1:n4), tl(i), y(i))
               end if
               exit
            end if
         end do
      end do
      x(is) = xu(0); x(i1) = xl(nmax)
      y(is) = yu(0); y(i1) = yl(nmax)
      ! 格子番号整理
      do i = is, i1
         x0(i) = x(is + i1 - i)
         y0(i) = y(is + i1 - i)
      end do
      do i = is, i1
         x(i) = x0(i)
         y(i) = y0(i)
      end do
      ! 処理終了 ********************************************************************************************
      return
   end subroutine HtypeInternalBoundary
!*******************************************************************************************************
!******** 外部境界 (H-type)                                                                        ********
!*******************************************************************************************************
   subroutine HtypeExternalBoundary( &
   &            is, ie, js, je, dom, x, y)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in)    :: is, ie, js, je
      real, intent(in)    :: dom
      real, intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer :: i
      integer :: ic
      real    :: a1, a2, a3, b1, b2, b3, c1, c2, c3
      ! 処理開始 ********************************************************************************************
      do i = is, ie
         ! 翼周り方向ベクトル
         if (i == is) then
            a1 = -x(i, js) + x(i + 1, js)
            a2 = -y(i, js) + y(i + 1, js)
         else if (i == ie) then
            a1 = -x(i - 1, js) + x(i, js)
            a2 = -y(i - 1, js) + y(i, js)
         else
            a1 = 0.5*(-x(i - 1, js) + x(i + 1, js))
            a2 = 0.5*(-y(i - 1, js) + y(i + 1, js))
         end if
         a3 = 0.0
         ! スパン方向ベクトル
         b1 = 0.0
         b2 = 0.0
         b3 = -1.0
         ! 法線方向ベクトル
         c1 = a2*b3 - a3*b2
         c2 = a3*b1 - a1*b3
         c3 = a1*b2 - a2*b1
         ! 翼周り部外部領域
         x(i, je) = c1/sqrt(c1**2 + c2**2 + c3**2)*dom + x(i, js)
         y(i, je) = c2/sqrt(c1**2 + c2**2 + c3**2)*dom + y(i, js)
      end do
      ! 処理終了 ********************************************************************************************
      return
   end subroutine HtypeExternalBoundary
!*******************************************************************************************************
!******** 側部境界 (H-type)                                                                        ********
!*******************************************************************************************************
   subroutine HtypeSideBoundary( &
   &            is, ie, js, je, rs, re, x, y)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in)    :: is, ie, js, je
      real, intent(in)    :: rs, re
      real, intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, pointer :: t(:)
      integer :: j
      real    :: rr
      ! 処理開始 ********************************************************************************************
      allocate (t(js:je))
      call GeometricInterpolation(rs/real(je - js + 1), je - js + 1, t, rr)
      do j = js + 1, je - 1
!  x(is,j) = x(is,js)
!  x(ie,j) = x(ie,js)
         x(is, j) = (x(is, je) - x(is, js))*t(j) + x(is, js)
         x(ie, j) = (x(ie, je) - x(ie, js))*t(j) + x(ie, js)
         y(is, j) = (y(is, je) - y(is, js))*t(j) + y(is, js)
         y(ie, j) = (y(ie, je) - y(ie, js))*t(j) + y(ie, js)
      end do
      deallocate (t)
      ! 処理開始 ********************************************************************************************
      return
   end subroutine HtypeSideBoundary
!*******************************************************************************************************
!******** Transfinite 補間                                                                        ********
!*******************************************************************************************************
   subroutine TransfiniteInterpolation( &
   &            is, ie, js, je, x, y)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in)    :: is, ie, js, je
      real, intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, pointer :: alp1(:), alp2(:), bet1(:), bet2(:)
      integer :: i, j
      ! 処理開始 ********************************************************************************************
      ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      allocate (alp1(is:ie), alp2(is:ie), bet1(js:je), bet2(js:je))
      ! Transfinite 補間 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 媒介変数
      do i = is, ie
         alp1(i) = real(ie - i)/real(ie)
         alp2(i) = 1.0 - alp1(i)
      end do
      do j = js, je
         bet1(j) = real(je - j)/real(je)
         bet2(j) = 1.0 - bet1(j)
      end do
      ! 補間
      call Transfinite2D(is, ie, js, je, alp1, alp2, bet1, bet2, x)
      call Transfinite2D(is, ie, js, je, alp1, alp2, bet1, bet2, y)
      ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      deallocate (alp1, alp2, bet1, bet2)
      ! 処理終了 ********************************************************************************************
      return
   end subroutine TransfiniteInterpolation
!*******************************************************************************************************
!******** 二境界法に基づく格子生成                                                                ********
!*******************************************************************************************************
   subroutine GenerationTwoBoundary( &
   &            is, ie, js, je, rs, re, t1, t2, x, y)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in)    :: is, ie, js, je
      real, intent(in)    :: rs, re, t1, t2
      real, intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, pointer :: etabar(:)
      real    :: rr
      ! 処理開始 ********************************************************************************************
      ! メモリ確保
      allocate (etabar(js:je))
      ! 媒介変数
      call GeometricInterpolation(rs/real(je - js + 1), je - js + 1, etabar, rr)
      ! 二境界法
      call TwoBoundaryMethod2D( &
      &      is, ie, js, je, t1, t2, etabar, x, y)
      ! メモリ解放
      deallocate (etabar)
      ! 処理終了 ********************************************************************************************
      return
   end subroutine GenerationTwoBoundary
!*******************************************************************************************************
!******** 三次元化 (スパン方向変化なし)                                                                ********
!*******************************************************************************************************
   subroutine ThreeDimensionalized( &
   &            is, ie, js, je, ks, ke, dom, x, y, z)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in)    :: is, ie, js, je, ks, ke
      real, intent(in)    :: dom
      real, intent(inout) :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke)
      real, intent(out)   :: z(is:ie, js:je, ks:ke)
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer :: i, j, k
      ! 処理開始 ********************************************************************************************
      do k = ks, ke
      do j = js, je
      do i = is, ie
         x(i, j, k) = x(i, j, ks)
         y(i, j, k) = y(i, j, ks)
         z(i, j, k) = -0.5*dom + dom*real(k)/real(ke)
      end do
      end do
      end do
      ! 処理終了 ********************************************************************************************
      return
   end subroutine ThreeDimensionalized
!*******************************************************************************************************
!******** 翼のフラグ                                                                                ********
!*******************************************************************************************************
   subroutine ViewBlade( &
   &            m, is, ie, js, je, ks, ke, ibs, ibe, j0, f)
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
         if (ibs <= i .and. i <= ibe .and. j == j0) then
            f(i, j, k) = 1
         else
            f(i, j, k) = 0
!    f(i,j,k) = m
         end if
      end do
      end do
      end do
      ! 処理終了 ********************************************************************************************
   end subroutine ViewBlade
!*******************************************************************************************************
!******** MicroAVSファイル作成（三次元）                                                          ********
!*******************************************************************************************************
   subroutine MakeMAVSFile3D( &
   &            strdir, strname, ext, is, ie, js, je, ks, ke, f, x, y, z)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      character, intent(in)  :: strdir*(*), strname*(*), ext*4
      integer, intent(in)  :: is, ie, js, je, ks, ke
      integer, intent(in)  :: f(is:ie, js:je, ks:ke)
      real, intent(in)  :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke), z(is:ie, js:je, ks:ke)
      ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      character, parameter :: strbin*4 = '.bin'
      character, parameter :: strfld*4 = '.fld'
      integer, parameter :: ndim = 3
      integer, parameter :: nspace = 3
      integer, parameter :: veclen = 1
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer   :: nrkind, nskip
      real      :: r
      integer   :: i, j, k, n
      ! 処理開始 ********************************************************************************************
      ! ヘッダファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      open (1, file=trim(strdir)//trim(strname)//trim(strfld), form='formatted')
      write (1, '(a)') '# AVS field file'
      write (1, '(a,i1)') 'ndim   = ', ndim
      write (1, '(a,i4)') 'dim1   = ', ie - is + 1
      write (1, '(a,i4)') 'dim2   = ', je - js + 1
      write (1, '(a,i4)') 'dim3   = ', ke - ks + 1
      write (1, '(a)') 'label  = flag'
      write (1, '(a,i1)') 'nspace = ', nspace
      write (1, '(a,i2)') 'veclen = ', veclen
      write (1, '(a)') 'data   = float'
      write (1, '(a)') 'field  = irregular'
      select case (ext)
         ! Binary ----------------------------------------------------------------------------------------------
      case (strbin)
         nrkind = kind(r)
         nskip = nrkind*(ie - is + 1)*(je - js + 1)*(ke - ks + 1) + 8
         do n = 1, veclen
            write (1, '((a,i2), (x,2a), (x,a), (x,a,i11), 2(x,a,i1))') &
            & 'variable ', n, &
            & 'file = ', trim(strname)//ext, &
            & 'filetype = binary', &
            & 'skip = ', 4 + nskip*(n - 1), &
            & 'stride = ', 1, &
            & 'close = ', 1
         end do
         do n = 1, nspace
            write (1, '((a,i1), (x,2a), (x,a), (x,a,i11), 2(x,a,i1))') &
            & 'coord ', n, &
            & 'file = ', trim(strname)//ext, &
            & 'filetype = binary', &
            & 'skip = ', 4 + nskip*(veclen + n - 1), &
            & 'stride = ', 1, &
            & 'close = ', 1
         end do
         ! Aschii ----------------------------------------------------------------------------------------------
      case default
         do n = 1, veclen
            write (1, '((a,i2), (1x,2a), (1x,a), (1x,a,i1), 2(1x,a,i2), (1x,a,i1))') &
            & 'variable ', n, &
            & 'file = ', trim(strname)//ext, &
            & 'filetype = ascii', &
            & 'skip = ', 0, &
            & 'offset = ', n - 1, &
            & 'stride = ', veclen + nspace, &
            & 'close = ', 1
         end do
         do n = 1, nspace
            write (1, '((a,i1), (1x,2a), (1x,a), (1x,a,i1), 2(1x,a,i2), (1x,a,i1))') &
            & 'coord ', n, &
            & 'file = ', trim(strname)//ext, &
            & 'filetype = ascii', &
            & 'skip = ', 0, &
            & 'offset = ', veclen + n - 1, &
            & 'stride = ', veclen + nspace, &
            & 'close = ', 1
         end do
      end select
      close (1)
      ! データファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      select case (ext)
         ! Binary ----------------------------------------------------------------------------------------------
      case (strbin)
         open (1, file=trim(strdir)//trim(strname)//ext, form='unformatted', status='replace')
         write (1) real(f); write (1) x; write (1) y; write (1) z
         close (1)
         close (1)
         ! Aschii ----------------------------------------------------------------------------------------------
      case default
         open (1, file=trim(strdir)//trim(strname)//ext, form='formatted', status='replace')
         do k = ks, ke
         do j = js, je
         do i = is, ie
            write (1, '(e16.8e3, 3(x,e16.8e3))') real(f(i, j, k)), x(i, j, k), y(i, j, k), z(i, j, k)
         end do
         end do
         end do
         close (1)
      end select
      ! 処理終了 ********************************************************************************************
      return
   end subroutine MakeMAVSFile3D
!*******************************************************************************************************
!******** vtkファイル出力サブルーチン                                                                 ********
!*******************************************************************************************************
   subroutine OutputPara_bin( &
   &      strdir, strname, is, ie, js, je, ks, ke, &
   &      x, y, z)
      implicit none
      !mainroutine_variable
      character, intent(in)  :: strdir*(*), strname*(*)
      integer, intent(in)  :: is, ie, js, je, ks, ke
      real, intent(in)  :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke), z(is:ie, js:je, ks:ke)
      !subroutine_variable
      character(len=4), parameter                :: strvtk = '.vtk'
      character(len=1), parameter                :: newline = char(10)
      character(len=200)        :: strnum
      integer                :: i, j, k, n
      integer                :: ni, nj, nk
      integer                :: npoint

      npoint = (ie - is)*(je - js)*(ke - ks)
      open (unit=1, &
            file=trim(strdir)//trim(strname)//trim(strvtk), &
            form='unformatted', &
            access='stream', &          ! ← ここを stream に
            convert='big_endian', &      ! GNU拡張（OK）。気になるなら -fconvert で全体指定でも可
            action='write')
      write (1) '# vtk DataFile Version 3.0'//newline
      write (1) 'vtk output'//newline
      write (1) 'BINARY'//newline
      write (1) 'DATASET UNSTRUCTURED_GRID'//newline
      write (strnum, *) npoint*8
      write (1) 'POINTS'//trim(strnum)//' float'//newline
      do n = 0, npoint - 1
         ni = mod(mod(n, (ie - is)*(je - js)), (ie - is))
         nj = int(mod(n, (ie - is)*(je - js))/(ie - is))
         nk = int(n/((ie - is)*(je - js)))
         write (1) x(ni, nj, nk), y(ni, nj, nk), z(ni, nj, nk)
         write (1) x(ni + 1, nj, nk), y(ni + 1, nj, nk), z(ni + 1, nj, nk)
         write (1) x(ni + 1, nj + 1, nk), y(ni + 1, nj + 1, nk), z(ni + 1, nj + 1, nk)
         write (1) x(ni, nj + 1, nk), y(ni, nj + 1, nk), z(ni, nj + 1, nk)
         write (1) x(ni, nj, nk + 1), y(ni, nj, nk + 1), z(ni, nj, nk + 1)
         write (1) x(ni + 1, nj, nk + 1), y(ni + 1, nj, nk + 1), z(ni + 1, nj, nk + 1)
         write (1) x(ni + 1, nj + 1, nk + 1), y(ni + 1, nj + 1, nk + 1), z(ni + 1, nj + 1, nk + 1)
         write (1) x(ni, nj + 1, nk + 1), y(ni, nj + 1, nk + 1), z(ni, nj + 1, nk + 1)
      end do
      write (1) newline
      write (strnum, *) npoint, npoint*9
      write (1) 'CELLS'//trim(strnum)//newline
      do n = 0, npoint - 1
         write (1) 8, n*8 + 0, n*8 + 1, n*8 + 2, n*8 + 3, n*8 + 4, n*8 + 5, n*8 + 6, n*8 + 7
      end do
      write (1) newline
      write (strnum, *) npoint
      write (1) 'CELL_TYPES'//trim(strnum)//newline
      do n = 0, npoint - 1
         write (1) 12
      end do
      write (1) newline
      write (strnum, *) npoint*8
      write (1) 'POINT_DATA'//trim(strnum)//newline
      write (1) 'SCALARS n float'//newline
      write (1) 'LOOKUP_TABLE default'//newline
      do n = 0, npoint - 1
         write (1) real(n*8 + 0)
         write (1) real(n*8 + 1)
         write (1) real(n*8 + 2)
         write (1) real(n*8 + 3)
         write (1) real(n*8 + 4)
         write (1) real(n*8 + 5)
         write (1) real(n*8 + 6)
         write (1) real(n*8 + 7)
      end do
      close (1)

   end subroutine OutputPara_bin
! 定義終了 *********************************************************************************************
end program GridGeneration_NACA
