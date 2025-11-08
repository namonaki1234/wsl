!*******************************************************************************************************
!*******************************************************************************************************
!******** NACA翼用の着氷モデル検証計算共通モジュール						********
!********					      2014.04.15  PROGRAMMED BY RYOSUKE HAYASHI ********
!*******************************************************************************************************
!*******************************************************************************************************
module Package_NACA
 ! 変数宣言 ********************************************************************************************
 implicit none
 public
 ! 構造体宣言 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 格子 ------------------------------------------------------------------------------------------------
 type Type_Grid
  integer, pointer :: f(:, :, :)
  real   , pointer :: x(:, :, :), y(:, :, :), z(:, :, :)
  integer :: is, ie, js, je, ks, ke
  integer :: i1, i2, i3, i4, i5, &
  &          j1, j2, j3, j4, j5, &
  &          k1, k2, k3, k4, k5
  real    :: xDom, yDom, zDom
 end type Type_Grid
 public :: Type_Grid
 type(Type_Grid), pointer :: Grd(:)
 ! 重合格子法補間係数 ----------------------------------------------------------------------------------
 type Type_Overset
  integer, pointer :: ip(:, :, :), jp(:, :, :), kp(:, :, :)		! 補間点
  integer, pointer :: fOver(:, :, :)					! 補間点のフラグ
  real   , pointer :: term1(:, :, :), term2(:, :, :), &
  &		      term3(:, :, :), term4(:, :, :), &
  &                   term5(:, :, :), term6(:, :, :), &
  &		      term7(:, :, :), term8(:, :, :)			! 補間係数
  integer :: is, ie, js, je, ks, ke					! 補間範囲
  integer :: m1								! 被補間ブロック
  integer :: m2								! 加補間ブロック
  real    :: MGN							! 補間係数許容誤差
 end type Type_Overset
 type(Type_Overset), pointer :: OSG(:)
 ! 流れ場 ----------------------------------------------------------------------------------------------
 type Type_Flow
  real   , pointer :: x(:, :, :), y(:, :, :), z(:, :, :)		! 格子座標
  integer, pointer :: f(:, :, :)
  real   , pointer :: xix(:, :, :), xiy(:, :, :), xiz(:, :, :), &
  &                   etx(:, :, :), ety(:, :, :), etz(:, :, :), &
  &                   zex(:, :, :), zey(:, :, :), zez(:, :, :), &
  &	              jac(:, :, :)					! メトリックス
  real   , pointer :: rho(:, :, :), u(:, :, :), v(:, :, :), w(:, :, :), &
  &                   p(:, :, :), t(:, :, :), &
  &                   kin(:, :, :), eps(:, :, :), &
  &                   mu(:, :, :), mut(:, :, :)				! 物理量
  real   , pointer :: qh(:, :, :, :), dqh(:, :, :, :), &
  &                   dqc(:, :, :, :), dqd(:, :, :, :), &
  &                   dqp(:, :, :, :), dqr(:, :, :, :), &
  &                   qh0(:, :, :, :), dqh0(:, :, :, :)			! 流束関数
  real   , pointer :: dt(:, :, :)					! 時間刻み
  real   , pointer :: Res(:)						! 残差
  real   , pointer :: nimp(:, :)					! 衝突個数
  real   , pointer :: uimp(:, :), vimp(:, :), wimp(:, :)		! 衝突速度の和
  real   , pointer :: simp(:, :) 					! 衝突面積
  real   , pointer :: CE(:, :)						! 収集効率
  real   , pointer :: ImpMass(:, :)					! 衝突質量
  integer, pointer :: fRough(:, :)					! 壁面粗さのフラグ
  real   , pointer :: RH(:, :)						! 局所粗さ
  integer :: is, ie, js, je, ks, ke					! 格子数
  integer :: i1, i2, i3, i4, i5, &
  &          j1, j2, j3, j4, j5, &
  &          k1, k2, k3, k4, k5						! 境界格子番号
 end type Type_Flow
 type(Type_Flow), pointer :: Flw(:)
 ! 液滴軌道 --------------------------------------------------------------------------------------------
 type Type_Drop
  integer, pointer :: m(:)						! 液滴軌道の格子番号
  real   , pointer :: x(:), y(:), z(:)					! 液滴軌道の位置
  real   , pointer :: u(:), v(:), w(:)					! 液滴の速度
  real   , pointer :: mr(:)						! 液滴質量比
  real   , pointer :: cd(:)						! 抗力係数
  real   , pointer :: Re(:)						! レイノルズ数
  real   , pointer :: AsRatio(:)					! 液滴のアスペクト比
  integer  :: step							! 液滴軌道の計算ステップ
 end type Type_Drop
 type(Type_Drop), pointer :: Drp(:)
 ! 着氷 ------------------------------------------------------------------------------------------------
 type Type_Ice
  real   , pointer :: x(:), y(:), z(:)					! 格子座標
  integer, pointer :: f(:)
  real   , pointer :: Uimp(:), Vimp(:)					! 液滴衝突速度
  real   , pointer :: CE(:)						! 液滴収集効率
  real   , pointer :: SA(:)						! 表面セル面積
  real   , pointer :: Yp(:)						! 表面から一点目の距離
  real   , pointer :: Up(:), Vp(:), Wp(:)				! 表面から一点目の速度成分
  real   , pointer :: Velp(:)						! 表面から一点目の速度
  real   , pointer :: nup(:)						! 表面から一点目の動粘性係数
  real   , pointer :: Ptp(:)						! 表面から一点目の全圧
  real   , pointer :: Rho0(:)						! 翼表面密度
  real   , pointer :: mu0(:)						! 翼表面粘性係数
  real   , pointer :: Ts0(:)						! 翼表面温度
  real   , pointer :: Ta(:)						! 周囲流体の温度
  real   , pointer :: Pa(:)						! 周囲流体の圧力
  real   , pointer :: Ue(:)						! 境界層外端の速度
  integer, pointer :: fRough(:)						! 壁面粗さのフラグ
  integer, pointer :: nPhase(:)						! 表面の状態
  real   , pointer :: Mes(:)						! 蒸発・昇華の質量
  real   , pointer :: Mim(:)						! 液滴衝突の質量
  real   , pointer :: Mrin(:)						! Runback-in の質量
  real   , pointer :: Mrout(:)						! Runback-out の質量
  real   , pointer :: Mst(:)						! 検査体積内に残る質量
  real   , pointer :: Mac(:)						! 堆積する氷の質量
  real   , pointer :: Bi(:)						! 氷層厚さ (絶対量)
  real   , pointer :: Bi0(:)						! 氷層厚さ (絶対量)
  real   , pointer :: dBi(:)						! 氷層厚さ (変化量)
  real   , pointer :: Ti(:)						! 氷層温度
  real   , pointer :: xi(:), yi(:)					! 着氷翼座標
  integer :: is, ie
  integer :: i1, i2, i3
 end type Type_Ice
 type(Type_Ice), pointer :: Ice(:)
 ! 共有定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: pi    = 3.1415926535897932			! 円周率
 real   , parameter :: zero  = 1.0e-20					! ゼロ割防止
 real   , parameter :: InvSix = 1.0 / 6.0
 real   , parameter :: IceTs = 273.15					! 氷結温度 [K]
 real   , parameter :: knot  = 0.514444444				! 1 [knot] = 0.514444444 [m/s]
 integer, parameter :: ls = 1
 integer, parameter :: le = 7						! 流れ場配列数
 integer, parameter :: ms = 1						! Main Grid 構造体番号
 integer, parameter :: me = 2						! Sub Grid 構造体番号
 ! 共有変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, pointer :: IceIn(:, :, :)					! 氷の中の点
 integer :: BCNum, TurbNum						! 流れ場モデルの選択
 integer :: ThermoNum, RoughNum, RunbackNum, DropShedNum, IceShedNum	! 着氷モデル選択
 integer :: DragNum
 real    :: span, chord, AOA						! スパン長，コード長，迎角
 real    :: VelExp, TsExp, PsExp					! 実験条件
 real    :: VelIn, PtIn, TtIn, PsOut, TsOut				! 流入出境界固定値
 real    :: TsSth, muSth						! サザーランドの式の定数
 real    :: ReRef, aRef, LRef, RhoRef					! 無次元化参照値
 integer :: nRunge							! ルンゲ・クッタ法の階層
 integer :: nTVD							! TVDスキームの精度
 real    :: Cn, Cnd, Cmu, eTVD						! 計算条件
 real    :: Rg, gamma, Pr, Prt, Ret, s1					! 物性値
 real    :: dtMin, dti							! 時間刻み
 real    :: LmtPro, LmtAve						! k-eリミッター定数
 integer :: nStart, nCalMax, nDrpCalMax, nCount				! 計算回数
 integer :: nOutputLog, nDrpOutputLog, nOutputCount, nDrpOutputCount	! ファイル出力周期
 real    :: time							! 経過時間
 integer :: nDrop							! 液滴数のカウンタ
 integer :: nDrpIn							! 投入液滴数
 integer :: nDrpFile							! 液滴軌道ファイル出力数
 real    :: LWC								! 水分含有率
 real    :: MVD								! 平均液滴径
 real    :: RhoD							! 液滴密度
 real    :: SigD							! 液滴表面張力
 real    :: muD								! 液滴粘性係数
 real    :: DrpInArea							! 液滴投入面積
 real    :: DrpInVel							! 液滴投入速度
 real    :: IceTime, IceTimeMax						! 着氷時間
 integer :: IceStep, IceStepMax						! 計算ステップ
 integer :: ExpCaseNum							! 実験ケース
 integer :: nDrpSpl, nDrpBou						! スプラッシュ・バウンドの液滴数
 integer :: nDrp2Imp							! 2回衝突した液滴数
 integer :: nDrpBrk							! 分裂した液滴数
 logical :: fTime, fSteady, fSlip
 ! ディレクトリ名 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character :: GrdInDir * 50, GrdOutDir * 50, OSGDir * 50, &
 &            FlwIniDir * 50, FlwCalInDir * 50, FlwCalOutDir * 50, FlwViewDir * 50, IceViewDir * 50, &
 &            DrpTraDir * 50, DrpImpDir * 50, IceCalInDir * 50, IceCalOutDir * 50, VldDir * 50
 ! ファイル名 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: ExpCaseFile          * 19 = 'Input_ExpCaseNumber'
 character, parameter :: ExpConditionFile     * 19 = 'Input_ExpCondition'
 character, parameter :: CalSetFile           * 10 = 'CalSetting'
 character, parameter :: ND_CalSetFile        * 13 = 'ND_CalSetting'
 character, parameter :: RslFile              * 10 = 'Resolution'
 character, parameter :: IceRslFile           * 13 = 'IceResolution'
 character, parameter :: CtypePointFile	      * 14 = 'CtypeGridPoint'
 character, parameter :: GrdFile              *  4 = 'Grid'
 character, parameter :: OverAreaFile         * 11 = 'OversetArea'
 character, parameter :: OverCoeFile          * 18 = 'OversetCoefficient'
 character, parameter :: MetFile              *  7 = 'Metrics'
 character, parameter :: IniFlxFile           * 11 = 'InitialFlux'
 character, parameter :: FlxFile              *  4 = 'Flux'
 character, parameter :: PhyFile              *  7 = 'Physics'
 character, parameter :: ND_RefFile           *  6 = 'ND_Ref'
 character, parameter :: ND_SutherlandFile    * 13 = 'ND_Sutherland'
 character, parameter :: ND_InflowFile        * 18 = 'ND_InflowCondition'
 character, parameter :: ND_GrdFile           *  7 = 'ND_Grid'
 character, parameter :: ND_MetFile           * 10 = 'ND_Metrics'
 character, parameter :: ND_IniFlxFile        * 14 = 'ND_InitialFlux'
 character, parameter :: ND_FlxFile           *  7 = 'ND_Flux'
 character, parameter :: ND_DrpImpiDataFile   * 18 = 'ND_ImpingementData'
 character, parameter :: ND_CollectionFile    * 20 = 'ND_CollectionEfficiency'
 character, parameter :: FlwCalCountFile      * 12 = 'FlowCalCount'
 character, parameter :: FlwCalLogFile        *  7 = 'FlowLog'
 character, parameter :: DrpCalLogFile        *  7 = 'DropLog'
 character, parameter :: DrpTrajectoryFile    * 14 = 'DropTrajectory'
 character, parameter :: DrpImpiAreaFile      * 15 = 'ImpingementArea'
 character, parameter :: DrpImpiDataFile      * 15 = 'ImpingementData'
 character, parameter :: DrpInConditionFile   * 15 = 'DropInCondition'
 character, parameter :: CollectionFile       * 20 = 'CollectionEfficiency'
 character, parameter :: IceCalLogFile        *  6 = 'IceLog'
 character, parameter :: IceThickTemFile      * 11 = 'IceThickTem'
 character, parameter :: IceMassFlxFile       * 11 = 'IceMassFlux'
 character, parameter :: IceBladeFile         *  8 = 'IceBlade'
 character, parameter :: IceStepFile          * 17 = 'Input_IceTimeStep'
 character, parameter :: WaterBladeFile       * 10 = 'WaterBlade'
 character, parameter :: MessingerFile        * 13 = 'MessingerPara'
 character, parameter :: ExMessingerFile      * 15 = 'ExMessingerPara'
 character, parameter :: IceLimitPointFile    * 13 = 'IceLimitPoint'
 character, parameter :: IceInPointFile       * 10 = 'IceInPoint'
 character, parameter :: RoughFlagFile        * 13 = 'RoughnessFlag'
 character :: BlkName(ms:me) * 5
! 内部手続き開始 ***************************************************************************************
contains
!*******************************************************************************************************
!******** C 型格子の分割点入力									********
!*******************************************************************************************************
subroutine Input_CtypeGridPoint(strname, i1, i2, i3)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*)
 integer  , intent(out) :: i1, i2, i3
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname), form = 'formatted', status = 'old')
  read(1, '(i4,2(x,i4))') i1, i2, i3
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_CtypeGridPoint
!*******************************************************************************************************
!******** C 型格子の分割点出力									********
!*******************************************************************************************************
subroutine Output_CtypeGridPoint(strname, i1, i2, i3)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*)
 integer  , intent(in) :: i1, i2, i3
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname), form = 'formatted', status = 'replace')
  write(1, '(i4,2(x,i4))') i1, i2, i3
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_CtypeGridPoint
!*******************************************************************************************************
!******** 着氷限界割点入力									********
!*******************************************************************************************************
subroutine Input_IceLimitPoint(strname, i1, i2)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*)
 integer  , intent(out) :: i1, i2
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname), form = 'formatted', status = 'old')
  read(1, *) i1
  read(1, *) i2
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_IceLimitPoint
!*******************************************************************************************************
!******** 着氷限界点出力									********
!*******************************************************************************************************
subroutine Output_IceLimitPoint(strname, i1, i2)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*)
 integer  , intent(in) :: i1, i2
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname), form = 'formatted', status = 'replace')
  write(1, *) i1
  write(1, *) i2
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_IceLimitPoint
!*******************************************************************************************************
!******** 計算条件設定ファイルの入力								********
!*******************************************************************************************************
subroutine Input_CalSetting(strname)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*)
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname), form = 'formatted', status = 'old')
  read(1, '(i2)') ExpCaseNum
  read(1, '(i2,x,i2,2(x,i10))') IceStep, IceStepMax, nCount, nDrop
  read(1, '(e16.8e3,x,e16.8e3)') IceTime, dti
  read(1, '(e16.8e3,2(x,e16.8e3))') Span, Chord, AOA
  read(1, '(e16.8e3,2(x,e16.8e3))') VelExp, PsExp, TsExp
  read(1, '(e16.8e3,5(x,e16.8e3))') LWC, MVD, Rhod, Sigd, mud, IceTimeMax
  read(1, '(e16.8e3,4(x,i10))') Cn, nRunge, nCalMax, nOutputLog, nOutputCount
  read(1, '(i3,3(x,e16.8e3))') nTVD, eTVD, LmtPro, LmtAve
  read(1, '(i3,2(x,l1))') TurbNum, fSteady, fTime
  read(1, '(i3,x,l1)') BCNum, fSlip
  read(1, '(i3,5(x,i2))') ThermoNum, RoughNum, RunbackNum, DropShedNum, IceShedNum, DragNum
  read(1, '(e16.8e3,5(x,e16.8e3))') Rg, gamma, Pr, Prt, Ret, Cmu
  read(1, '(e16.8e3,2(x,e16.8e3))') muSth, TsSth, s1
  read(1, '(e16.8e3,2(x,e16.8e3))') aRef, LRef, RhoRef
  read(1, '(e16.8e3,5(x,i10))') Cnd, nDrpIn, nDrpCalMax, nDrpFile, nDrpOutputLog, nDrpOutputCount
  read(1, '(4i10)') nDrpSpl, nDrpBou, nDrp2Imp, nDrpBrk
  read(1, '(e16.8e3,4(x,e16.8e3))') VelIn, TtIn, PtIn, TsOut, PsOut
 close(1)
 ! 処理終了 ********************************************************************************************
end subroutine Input_CalSetting
!*******************************************************************************************************
!******** 計算条件設定ファイルの出力								********
!*******************************************************************************************************
subroutine Output_CalSetting(strname)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*)
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname), form = 'formatted', status = 'replace')
  write(1, '(i2)') ExpCaseNum
  write(1, '(i2,x,i2,2(x,i10))') IceStep, IceStepMax, nCount, nDrop
  write(1, '(e16.8e3,x,e16.8e3)') IceTime, dti
  write(1, '(e16.8e3,2(x,e16.8e3))') Span, Chord, AOA
  write(1, '(e16.8e3,2(x,e16.8e3))') VelExp, PsExp, TsExp
  write(1, '(e16.8e3,5(x,e16.8e3))') LWC, MVD, Rhod, Sigd, mud, IceTimeMax
  write(1, '(e16.8e3,4(x,i10))') Cn, nRunge, nCalMax, nOutputLog, nOutputCount
  write(1, '(i3,3(x,e16.8e3))') nTVD, eTVD, LmtPro, LmtAve
  write(1, '(i3,2(x,l1))') TurbNum, fSteady, fTime
  write(1, '(i3,x,l1)') BCNum, fSlip
  write(1, '(i3,5(x,i2))') ThermoNum, RoughNum, RunbackNum, DropShedNum, IceShedNum, DragNum
  write(1, '(e16.8e3,5(x,e16.8e3))') Rg, gamma, Pr, Prt, Ret, Cmu
  write(1, '(e16.8e3,2(x,e16.8e3))') muSth, TsSth, s1
  write(1, '(e16.8e3,2(x,e16.8e3))') aRef, LRef, RhoRef
  write(1, '(e16.8e3,5(x,i10))') Cnd, nDrpIn, nDrpCalMax, nDrpFile, nDrpOutputLog, nDrpOutputCount
  write(1, '(4i10)') nDrpSpl, nDrpBou, nDrp2Imp, nDrpBrk
  write(1, '(e16.8e3,4(x,e16.8e3))') VelIn, TtIn, PtIn, TsOut, PsOut
 close(1)
 ! 処理終了 ********************************************************************************************
end subroutine Output_CalSetting
!*******************************************************************************************************
!******** 実験条件入力										********
!*******************************************************************************************************
subroutine Input_ExpCondition(strname)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*)
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname), form = 'formatted', status = 'old')
  read(1, *) Chord							! [m]
  read(1, *) TsExp							! [℃]
  read(1, *) VelExp							! [m/s]
  read(1, *) PsExp							! [kPa]
  read(1, *) MVD							! [μm]
  read(1, *) LWC							! [g/m^3]
  read(1, *) IceTimeMax							! [min.]
  read(1, *) Span							! [m]
  read(1, *) AOA							! [deg.]
  read(1, *) RhoD							! [g/cm^3]
  read(1, *) SigD							! [N/m]
  read(1, *) muD							! [N・s/m^2]
 close(1)
 TsExp      = TsExp + IceTs
 PsExp      = PsExp * 1.0e+3
 MVD        = MVD * 1.0e-6
 LWC        = LWC * 1.0e-3
 IceTimeMax = IceTimeMax * 60.0
 AOA        = AOA * pi / 180.0
 RhoD       = RhoD * 1.0e+3
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_ExpCondition
!*******************************************************************************************************
!******** ブロック名設定 									********
!*******************************************************************************************************
subroutine Set_BlockName
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m
 ! 処理開始 ********************************************************************************************
 BlkName(1) = 'Main_'
 BlkName(2) = 'Sub_'
 do m = ms, me
  BlkName(m) = trim(BlkName(m))
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine Set_BlockName
!*******************************************************************************************************
!******** 座標系の無次元化                                              			********
!*******************************************************************************************************
subroutine NondimensionalizedCoord3D( &
&            is, ie, js, je, ks, ke, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je, ks, ke
 real   , intent(inout) :: x(is:ie, js:je, ks:ke)
 real   , intent(inout) :: y(is:ie, js:je, ks:ke)
 real   , intent(inout) :: z(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 !$omp parallel do default(shared) private(i, j, k)
 do k = ks, ke
 do j = js, je
 do i = is, ie
  x(i,j,k) = x(i,j,k) / LRef
  y(i,j,k) = y(i,j,k) / LRef
  z(i,j,k) = z(i,j,k) / LRef
 enddo
 enddo
 enddo
 !$omp end parallel do
 ! 処理終了 ********************************************************************************************
 return
end subroutine NondimensionalizedCoord3D
!*******************************************************************************************************
!******** 物理量の無次元化                                              			********
!*******************************************************************************************************
subroutine NondimensionalizedPhysics3D( &
&            is, ie, js, je, ks, ke, rho, u, v, w, p, t, mu, kin, eps, mut )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je, ks, ke
 real   , intent(inout) :: rho(is:ie, js:je, ks:ke)
 real   , intent(inout) :: u(is:ie, js:je, ks:ke)
 real   , intent(inout) :: v(is:ie, js:je, ks:ke)
 real   , intent(inout) :: w(is:ie, js:je, ks:ke)
 real   , intent(inout) :: p(is:ie, js:je, ks:ke)
 real   , intent(inout) :: t(is:ie, js:je, ks:ke)
 real   , intent(inout) :: kin(is:ie, js:je, ks:ke)
 real   , intent(inout) :: eps(is:ie, js:je, ks:ke)
 real   , intent(inout) :: mu(is:ie, js:je, ks:ke)
 real   , intent(inout) :: mut(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 !$omp parallel do default(shared) private(i, j, k)
 do k = ks, ke
 do j = js, je
 do i = is, ie
  rho(i,j,k) = rho(i,j,k) / RhoRef
  u(i,j,k)   = u(i,j,k)   / aRef
  v(i,j,k)   = v(i,j,k)   / aRef
  w(i,j,k)   = w(i,j,k)   / aRef
  p(i,j,k)   = p(i,j,k)   / (RhoRef * aRef**2)
  t(i,j,k)   = t(i,j,k)   / aRef**2
  mu(i,j,k)  = mu(i,j,k)  / (RhoRef * aRef * LRef)
  kin(i,j,k) = kin(i,j,k) / aRef**2
  eps(i,j,k) = eps(i,j,k) / (aRef**3 / LRef)
  mut(i,j,k) = mut(i,j,k) / (RhoRef * aRef * LRef)
 enddo
 enddo
 enddo
 !$omp end parallel do
 ! 処理終了 ********************************************************************************************
 return
end subroutine NondimensionalizedPhysics3D
!*******************************************************************************************************
!******** 座標系の有次元化                                              			********
!*******************************************************************************************************
subroutine DimensionalizedCoord3D( &
&            is, ie, js, je, ks, ke, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je, ks, ke
 real   , intent(inout) :: x(is:ie, js:je, ks:ke)
 real   , intent(inout) :: y(is:ie, js:je, ks:ke)
 real   , intent(inout) :: z(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 !$omp parallel do default(shared) private(i, j, k)
 do k = ks, ke
 do j = js, je
 do i = is, ie
  x(i,j,k) = x(i,j,k) * LRef
  y(i,j,k) = y(i,j,k) * LRef
  z(i,j,k) = z(i,j,k) * LRef
 enddo
 enddo
 enddo
 !$omp end parallel do
 ! 処理終了 ********************************************************************************************
 return
end subroutine DimensionalizedCoord3D
!*******************************************************************************************************
!******** 物理量の有次元化                                              			********
!*******************************************************************************************************
subroutine DimensionalizedPhysics3D( &
&            is, ie, js, je, ks, ke, rho, u, v, w, p, t, mu, kin, eps, mut )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je, ks, ke
 real   , intent(inout) :: rho(is:ie, js:je, ks:ke)
 real   , intent(inout) :: u(is:ie, js:je, ks:ke)
 real   , intent(inout) :: v(is:ie, js:je, ks:ke)
 real   , intent(inout) :: w(is:ie, js:je, ks:ke)
 real   , intent(inout) :: p(is:ie, js:je, ks:ke)
 real   , intent(inout) :: t(is:ie, js:je, ks:ke)
 real   , intent(inout) :: kin(is:ie, js:je, ks:ke)
 real   , intent(inout) :: eps(is:ie, js:je, ks:ke)
 real   , intent(inout) :: mu(is:ie, js:je, ks:ke)
 real   , intent(inout) :: mut(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 !$omp parallel do default(shared) private(i, j, k)
 do k = ks, ke
 do j = js, je
 do i = is, ie
  rho(i,j,k) = rho(i,j,k) * RhoRef
  u(i,j,k)   = u(i,j,k)   * aRef
  v(i,j,k)   = v(i,j,k)   * aRef
  w(i,j,k)   = w(i,j,k)   * aRef
  p(i,j,k)   = p(i,j,k)   * (RhoRef * aRef**2)
  t(i,j,k)   = t(i,j,k)   * aRef**2
  mu(i,j,k)  = mu(i,j,k)  * (RhoRef * aRef * LRef)
  kin(i,j,k) = kin(i,j,k) * aRef**2
  eps(i,j,k) = eps(i,j,k) * (aRef**3 / LRef)
  mut(i,j,k) = mut(i,j,k) * (RhoRef * aRef * LRef)
 enddo
 enddo
 enddo
 !$omp end parallel do
 ! 処理終了 ********************************************************************************************
 return
end subroutine DimensionalizedPhysics3D
!*******************************************************************************************************
!******** 数値流束のスムージング 								********
!*******************************************************************************************************
subroutine SmoothingFlux3D( &
 &           is, ie, js, je, ks, ke, ls, le, jac, qh )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je, ks, ke, ls, le
 real   , intent(in)    :: jac(is:ie, js:je, ks:ke)
 real   , intent(inout) :: qh (is:ie, js:je, ks:ke, ls:le)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: q0(:, :, :, :)
 integer :: i, j, k, l
 ! 処理開始 ********************************************************************************************
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( q0(is:ie, js:je, ks:ke, ls:le) )
 ! 元の物理量を保存 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do k = ks, ke
 do j = js, je
 do i = is, ie
  do l = ls + 1, le
   q0(i,j,k,l) = qh(i,j,k,l) / qh(i,j,k,1)
  enddo
  q0(i,j,k,1) = qh(i,j,k,1) * jac(i,j,k)
 enddo
 enddo
 enddo
 ! 周囲点の平均値を外挿 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do k = ks + 1, ke - 1
 do j = js + 1, je - 1
 do i = is + 1, ie - 1
  qh(i,j,k,1) = ( q0(i  ,j  ,k-1,1) + q0(i  ,j-1,k  ,1) + q0(i-1,j  ,k  ,1) &
  &             + q0(i+1,j  ,k  ,1) + q0(i  ,j+1,k  ,1) + q0(i  ,j  ,k+1,1) &
  &             ) * InvSix / jac(i,j,k)
  do l = ls + 1, le
   qh(i,j,k,l) = ( q0(i  ,j  ,k-1,l) + q0(i  ,j-1,k  ,l) + q0(i-1,j  ,k  ,l) &
   &             + q0(i+1,j  ,k  ,l) + q0(i  ,j+1,k  ,l) + q0(i  ,j  ,k+1,l) &
   &             ) * InvSix * qh(i,j,k,1)
  enddo
 enddo
 enddo
 enddo
 ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 deallocate(q0)
 ! 処理終了 ********************************************************************************************
 return
end subroutine SmoothingFlux3D
! 定義終了 *********************************************************************************************
end module Package_NACA
