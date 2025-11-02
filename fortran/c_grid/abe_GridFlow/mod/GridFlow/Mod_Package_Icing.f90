!*******************************************************************************************************
!*******************************************************************************************************
!******** パッケージ型モジュール								********
!******** 着氷モデル										********
!******** 			                      2013.07.08  PROGRAMMED BY RYOSUKE HAYASHI ********
!******** EMM 改良                                    2015.06.11     UPDATED BY MIKI    SHIMURA ********
!******** Ice Density                                 2015.07.17     UPDATED BY MIKI    SHIMURA ********
!******** エネルギー流束の修正                        2016.12.09     UPDATED BY SHO     URANAI  ********
!******** 霧氷・雨氷の厚さと温度の修正                2016.12.10     UPDATED BY SHO     URANAI  ********
!******** 氷の成長 修正			              2016.12.21     UPDATED BY SHO     URANAI  ********
!******** 加熱も考慮したモデル		              2018.09.26     UPDATED BY SHO     URANAI  ********
!*******************************************************************************************************
!*******************************************************************************************************
module Package_Icing
 ! 変数宣言 ********************************************************************************************
 implicit none
 private
 ! サブルーチン宣言 ************************************************************************************
 ! 共有サブルーチン ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 public :: Smooth2P1D, Smooth3P1D
 public :: Smooth4P2D, Smooth5P2D
 public :: WallFrictionVelocity2D, WallFrictionVelocity3D
 public :: HeatTransferCoefficient2D, HeatTransferCoefficient3D
 public :: RoughnessShinBond2D, RoughnessShinBond3D
 public :: RoughnessCIRAMIL2D, RoughnessCIRAMIL3D
 public :: RoughnessWright2D, RoughnessWright3D
 public :: RoughnessLEWICE2D, RoughnessLEWICE3D
 public :: IceDensity2D, IceDensity3D
 public :: SurfaceTemperature2D, SurfaceTemperature3D
 PUBLIC :: RenewSurfaceTemperature2D
 public :: DropImpingementOrg2D, DropImpingementOrg3D
 public :: DropImpingementExt2D, DropImpingementExt3D
 public :: EvapolationOrg2D, EvapolationOrg3D
 public :: EvapolationExt2D, EvapolationExt3D
 public :: HeatBalanceOrg2D, HeatBalanceOrg3D
 public :: HeatBalanceExt2D, HeatBalanceExt3D
 public :: HeatBalanceExt2DverUR
 public :: FreezingFractionOrg2D, FreezingFractionOrg3D
 public :: FreezingFractionExt2D, FreezingFractionExt3D
 public :: RunbackOutOrg2D, RunbackOutOrg3D
 public :: RunbackOutExt2D, RunbackOutExt3D
 public :: WaterSheddingOrg2D, WaterSheddingOrg3D
 public :: WaterSheddingExt2D, WaterSheddingExt3D
 public :: RunbackInOrg2D, RunbackInOrg3D
 public :: RunbackInExt2D, RunbackInExt3D
 public :: IceThicknessOrg2D, IceThicknessOrg3D
 public :: IceThicknessExt2D, IceThicknessExt3D
 public :: IceThicknessExt2DverUR
 public :: PhaseChangeExt2D, PhaseChangeExt3D
 public :: PhaseChangeExt2DverUR
 public :: IceTimeGrowth2D, IceTimeGrowth3D
 public :: SaveWaterTemperature2D
 public :: SaveIceThickness2D, SaveIceThickness3D
 public :: IceSheddingTotal
 public :: IceSheddingCell
 ! 共有定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: zero          = 1.0e-15
 real   , parameter :: pi            = 3.1415926535897932
 real   , parameter :: OneThird      = 0.33333333
 real   , parameter :: gamma         = 1.4				! 気体の比熱比
 real   , parameter :: Rg            = 287.1				! 気体のガス定数
 real   , parameter :: Cpa           = 1.006 * 1.0e3			! 空気の比熱 [J/(kg*K)]
 real   , parameter :: Cpi           = 2.050 * 1.0e3			! 氷の比熱   [J/(kg*K)]
 real   , parameter :: Cpw           = 4.218 * 1.0e3			! 水の比熱   [J/(kg*K)]
 real   , parameter :: e0            = 27.03				! 飽和蒸気圧
 real   , parameter :: ka            = 24.2e-3				! 空気の熱伝導率 [W/m*K]
 real   , parameter :: ki            = 2.18				! 氷の熱伝導率 [W/m*K]
 real   , parameter :: kw            = 0.571				! 水の熱伝導率 [W/m*K]
 real   , parameter :: LF            = 3.344 * 1.0e5			! 凝固の潜熱 [J/kg]
 real   , parameter :: LE            = 2.50  * 1.0e6			! 蒸発の潜熱 [J/kg]
 real   , parameter :: LS            = 2.8344 * 1.0e6			! 昇華の潜熱 [J/kg]
 real   , parameter :: Pr            = 0.72				! 層流プランルトル数
 real   , parameter :: Prt           = 0.9				! 乱流プランルトル数
! real   , parameter :: rr            = Pr**0.5				! 回復係数（層流）
 real   , parameter :: rr            = Pr**(1.0 / 3.0)			! 回復係数（乱流）
 real   , parameter :: eps           = 0.5				! 氷の放射率（0.5-0.8）
 real   , parameter :: muw           = 1.795 * 1.0e-3			! 水の粘性係数 [Pa*s]
 real   , parameter :: Rhor          = 880.0				! 霧氷の密度 [kg/m^3]
 real   , parameter :: Rhog          = 917.0				! 雨氷の密度 [kg/m^3]
 real   , parameter :: Rhow          = 999.0				! 水の密度 [kg/m^3]
 real   , parameter :: sigr          = 5.6704 * 1.0e-8			! ステファン・ボルツマン定数
 real   , parameter :: sigw          = 0.072				! 水の表面張力 [N/m]
 real   , parameter :: Tf            = 273.15				! 水の氷結温度 [K]
 real   , parameter :: AdhAve_AlRa06 = 286.94912 * 1.0e3		! アルミニウムの付着限界値 [Pa]
 real   , parameter :: AdhMax_AlRa06 = 399.51268 * 1.0e3		! アルミニウムの付着限界値 [Pa]
 real   , parameter :: AdhMin_AlRa06 = 155.94802 * 1.0e3		! アルミニウムの付着限界値 [Pa]
 real   , parameter :: AdhAve_TiRa06 = 274.67861 * 1.0e3		! チタンの付着限界値 [Pa]
 real   , parameter :: AdhMax_TiRa06 = 569.76729 * 1.0e3		! チタンの付着限界値 [Pa]
 real   , parameter :: AdhMin_TiRa06 = 197.44101 * 1.0e3		! チタンの付着限界値 [Pa]
 real   , parameter :: Thick_Wing    = 3.0e-6				! 翼材(Ti)の厚さ[m]
 real   , parameter :: k_wing        = 17.0				! 翼材(Ti)の熱伝導率[W/(m・K)]
 real   , parameter :: c_wing        = 519.0				! 翼材(Ti)の比熱[J/(kg・K)]
 real   , parameter :: rho_wing      = 4510.0				! 翼材(Ti)の密度[kg/(m^3)]
! real   , parameter :: Thick_Wing    = 2.0e-3				! 翼材(Al)の厚さ[m]
! real   , parameter :: k_wing        = 256.05				! 翼材(Al)の熱伝導率[W/(m・K)]
! real   , parameter :: c_wing        = 913.0				! 翼材(Al)の比熱[J/(kg・K)]
! real   , parameter :: rho_wing      = 2700.0				! 翼材(Al)の密度[kg/(m^3)]
 ! 内部手続き ******************************************************************************************
contains
!*******************************************************************************************************
!******** スムージング (一次元，対象点を含まない周囲2点の平均を外挿)				********
!*******************************************************************************************************
subroutine Smooth2P1D( &
&             nSmooth, is, ie, f )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: nSmooth
 integer, intent(in)    :: is, ie
 real   , intent(inout) :: f(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: f0(:)
 integer :: i, n
 ! 処理開始 ********************************************************************************************
 if(nSmooth < 1) return
 allocate( f0(is:ie) )
 do n = 1, nSmooth
  f0(:) = f(:)
  do i = is + 1, ie - 1
   f(i) = ( f0(i-1) + f0(i+1) ) * 0.5
  enddo
 enddo
 deallocate(f0)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Smooth2P1D
!*******************************************************************************************************
!******** スムージング (一次元，対象点を含む周囲3点の平均を外挿)				********
!*******************************************************************************************************
subroutine Smooth3P1D( &
&             nSmooth, is, ie, f )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: nSmooth
 integer, intent(in)    :: is, ie
 real   , intent(inout) :: f(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: f0(:)
 integer :: i, n
 ! 処理開始 ********************************************************************************************
 if(nSmooth < 1) return
 allocate( f0(is:ie) )
 do n = 1, nSmooth
  f0(:) = f(:)
  do i = is + 1, ie - 1
   f(i) = ( f0(i-1) + f0(i) + f0(i+1) ) * OneThird
  enddo
 enddo
 deallocate(f0)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Smooth3P1D
!*******************************************************************************************************
!******** スムージング (二次元，対象点を含まない周囲4点の平均を外挿)				********
!*******************************************************************************************************
subroutine Smooth4P2D( &
&             nSmooth, is, ie, ks, ke, f )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: nSmooth
 integer, intent(in)    :: is, ie, ks, ke
 real   , intent(inout) :: f(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: f0(:, :)
 integer :: i, k, n
 ! 処理開始 ********************************************************************************************
 if(nSmooth < 1) return
 allocate( f0(is:ie, ks:ke) )
 do n = 1, nSmooth
  f0(:,:) = f(:,:)
  do k = ks + 1, ke - 1
  do i = is + 1, ie - 1
   f(i,k) = ( f0(i,k-1) + f0(i-1,k) + f0(i+1,k) + f0(i,k+1) ) * 0.25
  enddo
  enddo
 enddo
 deallocate(f0)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Smooth4P2D
!*******************************************************************************************************
!******** スムージング (二次元，対象点を含む周囲5点の平均を外挿)				********
!*******************************************************************************************************
subroutine Smooth5P2D( &
&             nSmooth, is, ie, ks, ke, f )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: nSmooth
 integer, intent(in)    :: is, ie, ks, ke
 real   , intent(inout) :: f(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: f0(:, :)
 integer :: i, k, n
 ! 処理開始 ********************************************************************************************
 if(nSmooth < 1) return
 allocate( f0(is:ie, ks:ke) )
 do n = 1, nSmooth
  f0(:,:) = f(:,:)
  do k = ks + 1, ke - 1
  do i = is + 1, ie - 1
   f(i,k) = ( f0(i,k-1) + f0(i-1,k) + f0(i,k) + f0(i+1,k) + f0(i,k+1) ) * 0.2
  enddo
  enddo
 enddo
 deallocate(f0)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Smooth5P2D
!*******************************************************************************************************
!******** 壁面摩擦速度 (二次元)									********
!*******************************************************************************************************
subroutine WallFrictionVelocity2D( &
 &           is, ie, fRough, yp, up, nup, RH, utau )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 integer, intent(in)  :: fRough(is:ie)					! 粗さのフラグ
 real   , intent(in)  :: yp(is:ie)      				! 壁からの距離   (yp > 0)
 real   , intent(in)  :: up(is:ie)   					! 壁水平方向速度 (up > 0)
 real   , intent(in)  :: nup(is:ie)   					! 動粘度         (nup > 0)
 real   , intent(in)  :: RH(is:ie)					! 粗さ高さ       (RH  > 0)
 real   , intent(out) :: utau(is:ie)  					! 摩擦速度       (utau > 0)
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: kappa  = 0.40, b = 5.5				! 対数則のモデル定数
 integer, parameter :: nmax   = 100					! 繰返し計算数
 real   , parameter :: ResMin = 1.0e-8					! 収束判定値
 real   , parameter :: ksp1 = 4.0, ksp2 = 15.0, ksp3 = 55.0
 real   , parameter :: At = 9.5, Ar = 8.5
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 integer :: n
 real    :: f, futau, dutau
 real    :: ksp
 real    :: lnksp1, lnksp2, lnksp3
 real    :: alpha, beta, gamma
 real    :: db, dbutau
 ! 処理開始 ********************************************************************************************
 utau(:) = 0.0
 do i = is + 1, ie - 1
  select case( fRough(i) )
   ! 滑面 ----------------------------------------------------------------------------------------------
   case(0)
    utau(i) = sqrt( nup(i) * up(i) / yp(i))
    do n = 1, nmax
     f     = up(i) / utau(i) - log(yp(i) * utau(i) / nup(i)) / kappa - b
     futau = - (up(i) / utau(i) + 1 / kappa) / utau(i)
     if(abs(futau) .le. zero) futau = sign(zero, futau)
     dutau = - f / futau
     utau(i)  = max(zero, utau(i) + dutau)
     if( abs(dutau) .lt. utau(i) * ResMin ) exit
    enddo
    if(yp(i) * utau(i) / nup(i) .le. 11.635) then
      utau(i) = sqrt(nup(i) * up(i) / yp(i))
    endif
   ! 粗面 ----------------------------------------------------------------------------------------------
   case(1)
    lnksp1 = log(ksp1)
    lnksp2 = log(ksp2)
    lnksp3 = log(ksp3)
    beta   = 0.5 &
    &      * ( (b - Ar + lnksp3 / kappa) * (lnksp2**2 - lnksp1**2) &
    &        - (b - At + lnksp2 / kappa) * (lnksp3**2 - lnksp1**2) ) &
    &      / ( (b - At + lnksp2 / kappa) * (lnksp1 - lnksp3) &
    &        - (b - Ar + lnksp3 / kappa) * (lnksp1 - lnksp2) )
    alpha  = (b - At + lnksp2 / kappa) &
    &      / ((lnksp2 - beta)**2 - (lnksp1 - beta)**2)
    gamma  = - alpha * (lnksp1 - beta)**2
    utau(i) = sqrt(nup(i) * up(i) / yp(i))
    do n = 1, nmax
     ksp  = RH(i) * utau(i) / nup(i)
     if(ksp .lt. ksp1) then
       db     = 0.0
       dbutau = 0.0
      else if(ksp .gt. ksp3) then
       db     = b - Ar + log(ksp) / kappa
       dbutau = 1.0 / (kappa * ksp)
      else
       db     = alpha * (log(ksp) - beta)**2 + gamma
       dbutau = 2.0 * alpha * (log(ksp) - beta) / ksp
     endif
     f     = up(i) / utau(i) - log(yp(i) * utau(i) / nup(i)) / kappa - b + db
     futau = - (up(i) / utau(i) + 1 / kappa) / utau(i) + dbutau
     if(abs(futau) .le. zero) futau = sign(zero, futau)
     dutau = - f / futau
     utau(i)  = max(zero, utau(i) + dutau)
     if(abs(dutau) .lt. utau(i) * ResMin) exit
    enddo
    if(yp(i) * utau(i) / nup(i) .le. 11.635) then
      utau(i) = sqrt(nup(i) * up(i) / yp(i))
    endif
   ! 例外処理 ------------------------------------------------------------------------------------------
   case default; write(*, '(a)') '!!!!! Error : Roughness Number !!!!!'
    stop
  end select
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine WallFrictionVelocity2D
!*******************************************************************************************************
!******** 壁面摩擦速度 (三次元)									********
!*******************************************************************************************************
subroutine WallFrictionVelocity3D( &
 &           is, ie, ks, ke, fRough, yp, up, nup, RH, utau )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 integer, intent(in)  :: fRough(is:ie, ks:ke)				! 粗さのフラグ
 real   , intent(in)  :: yp(is:ie, ks:ke)      				! 壁からの距離   (yp > 0)
 real   , intent(in)  :: up(is:ie, ks:ke)   				! 壁水平方向速度 (up > 0)
 real   , intent(in)  :: nup(is:ie, ks:ke)   				! 動粘度         (nup > 0)
 real   , intent(in)  :: RH(is:ie, ks:ke)				! 粗さ高さ       (RH  > 0)
 real   , intent(out) :: utau(is:ie, ks:ke)  				! 摩擦速度       (utau > 0)
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: kappa  = 0.40, b = 5.5				! 対数則のモデル定数
 integer, parameter :: nmax   = 100					! 繰返し計算数
 real   , parameter :: ResMin = 1.0e-8					! 収束判定値
 real   , parameter :: ksp1 = 4.0, ksp2 = 15.0, ksp3 = 55.0
 real   , parameter :: At = 9.5, Ar = 8.5
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 integer :: n
 real    :: f, futau, dutau
 real    :: ksp
 real    :: lnksp1, lnksp2, lnksp3
 real    :: alpha, beta, gamma
 real    :: db, dbutau
 ! 処理開始 ********************************************************************************************
 utau(:, :) = 0.0
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  select case( fRough(i,k) )
   ! 滑面 ----------------------------------------------------------------------------------------------
   case(0)
    utau(i,k) = sqrt( nup(i,k) * up(i,k) / yp(i,k))
    do n = 1, nmax
     f     = up(i,k) / utau(i,k) - log(yp(i,k) * utau(i,k) / nup(i,k)) / kappa - b
     futau = - (up(i,k) / utau(i,k) + 1 / kappa) / utau(i,k)
     if(abs(futau) .le. zero) futau = sign(zero, futau)
     dutau = - f / futau
     utau(i,k)  = max(zero, utau(i,k) + dutau)
     if( abs(dutau) .lt. utau(i,k) * ResMin ) exit
    enddo
    if(yp(i,k) * utau(i,k) / nup(i,k) .le. 11.635) then
      utau(i,k) = sqrt(nup(i,k) * up(i,k) / yp(i,k))
    endif
   ! 粗面 ----------------------------------------------------------------------------------------------
   case(1)
    lnksp1 = log(ksp1)
    lnksp2 = log(ksp2)
    lnksp3 = log(ksp3)
    beta   = 0.5 &
    &      * ( (b - Ar + lnksp3 / kappa) * (lnksp2**2 - lnksp1**2) &
    &        - (b - At + lnksp2 / kappa) * (lnksp3**2 - lnksp1**2) ) &
    &      / ( (b - At + lnksp2 / kappa) * (lnksp1 - lnksp3) &
    &        - (b - Ar + lnksp3 / kappa) * (lnksp1 - lnksp2) )
    alpha  = (b - At + lnksp2 / kappa) &
    &      / ((lnksp2 - beta)**2 - (lnksp1 - beta)**2)
    gamma  = - alpha * (lnksp1 - beta)**2
    utau(i,k) = sqrt(nup(i,k) * up(i,k) / yp(i,k))
    do n = 1, nmax
     ksp  = RH(i,k) * utau(i,k) / nup(i,k)
     if(ksp .lt. ksp1) then
       db     = 0.0
       dbutau = 0.0
      else if(ksp .gt. ksp3) then
       db     = b - Ar + log(ksp) / kappa
       dbutau = 1.0 / (kappa * ksp)
      else
       db     = alpha * (log(ksp) - beta)**2 + gamma
       dbutau = 2.0 * alpha * (log(ksp) - beta) / ksp
     endif
     f     = up(i,k) / utau(i,k) - log(yp(i,k) * utau(i,k) / nup(i,k)) / kappa - b + db
     futau = - (up(i,k) / utau(i,k) + 1 / kappa) / utau(i,k) + dbutau
     if(abs(futau) .le. zero) futau = sign(zero, futau)
     dutau = - f / futau
     utau(i,k)  = max(zero, utau(i,k) + dutau)
     if(abs(dutau) .lt. utau(i,k) * ResMin) exit
    enddo
    if(yp(i,k) * utau(i,k) / nup(i,k) .le. 11.635) then
      utau(i,k) = sqrt(nup(i,k) * up(i,k) / yp(i,k))
    endif
   ! 例外処理 ------------------------------------------------------------------------------------------
   case default; write(*, '(a)') '!!!!! Error : Roughness Number !!!!!'
    stop
  end select
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine WallFrictionVelocity3D
!*******************************************************************************************************
!******** 熱伝達率 (二次元)									********
!*******************************************************************************************************
subroutine HeatTransferCoefficient2D( &
&            is, ie, istg, x, y, Ue, mu, rho, utau, RH, hc )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, istg
 real   , intent(in)  :: x(is:ie), y(is:ie)					! 
 real   , intent(in)  :: Ue(is:ie)					! 境界層外端の速度
 real   , intent(in)  :: mu(is:ie)					! 粘性係数
 real   , intent(in)  :: rho(is:ie)					! 密度
 real   , intent(in)  :: utau(is:ie)					! 壁面摩擦速度
 real   , intent(in)  :: RH(is:ie)					! 表面粗さ
 real   , intent(out) :: hc(is:ie)					! 熱伝達率 [W/(m^2*K)]
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, iht1, iht2
 real    :: Rek								! 粗さレイノルズ数
 real    :: Stk								! 粗さスタントン数
 real    :: St								! スタントン数
 real    :: cf								! 壁面摩擦係数
 real    :: Vds, AA(2,2), AB(2), r, xa, yb
 ! 処理開始 ********************************************************************************************
 hc(:) = 0.0
 iht1 = int(( is + ie ) / 2.0) - 2; iht2 = int(( is + ie ) / 2.0) + 2
 iht1 = min( iht2, istg-2 )       ; iht2 = max( iht1, istg+2 )
 AA(1,1) = x(istg) - x(istg+1); AA(1,2) = y(istg) - y(istg+1)
 AA(2,1) = x(istg) - x(istg-1); AA(2,2) = y(istg) - y(istg-1)
 AB(1) = 0.5 * ( x(istg)**2 - x(istg+1)**2 + y(istg)**2 - y(istg+1)**2 )
 AB(2) = 0.5 * ( x(istg)**2 - x(istg-1)**2 + y(istg)**2 - y(istg-1)**2 )
 xa = ( AA(2,2) * AB(1) - AA(1,2) * AB(2) ) / ( AA(1,1) * AA(2,2) - AA(1,2) * AA(2,1) )
 yb = ( AA(1,1) * AB(2) - AA(2,1) * AB(1) ) / ( AA(1,1) * AA(2,2) - AA(1,2) * AA(2,1) )
 r = sqrt( (x(istg) - xa)**2 + (y(istg) - yb)**2 )
 do i = is + 1, ie - 1
    Rek     = rho(i) * Ue(i) * RH(i) / mu(i)
!  if( Rek < 600.0 .and. i /= istg )then
!    ! 層流
!    iht1 = min( i, istg )
!    iht2 = max( i, istg )
!    Vds  = sum( Ue(iht1:iht2)**(1.87) * ds(iht1:iht2) )
!    hc(i) = 0.296 * ka * Ue(i)**(1.435) &
!    &       / sqrt( mu(i) / rho(i) * Vds )				! [J/(m^2*K*s)]
!   else
    ! 乱流
    Stk   = 1.92 * (Rek)**(-0.45) * Pr**(-0.8)
    cf    = 2.0 * (utau(i) / Ue(i))**2
    St    = 0.5 * cf / ( Prt + sqrt(0.5 * cf ) / Stk )
    hc(i) = St * rho(i) * Ue(i) * Cpa					! [J/(m^2*K*s)]
!  end if
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine HeatTransferCoefficient2D
!*******************************************************************************************************
!******** 熱伝達率 (三次元)									********
!*******************************************************************************************************
subroutine HeatTransferCoefficient3D( &
&            is, ie, ks, ke, Ue, mu, rho, utau, RH, hc )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: Ue(is:ie, ks:ke)				! 境界層外端の速度
 real   , intent(in)  :: mu(is:ie, ks:ke)
 real   , intent(in)  :: rho(is:ie, ks:ke)
 real   , intent(in)  :: utau(is:ie, ks:ke)				! 壁面摩擦速度
 real   , intent(in)  :: RH(is:ie, ks:ke)				! 表面粗さ
 real   , intent(out) :: hc(is:ie, ks:ke)				! 熱伝達率 [W/(m^2*K)]
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 real    :: Rek								! 粗さレイノルズ数
 real    :: Stk								! 粗さスタントン数
 real    :: St								! スタントン数
 real    :: cf								! 壁面摩擦係数
 ! 処理開始 ********************************************************************************************
 hc(:, :) = 0.0
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  Rek     = rho(i,k) * utau(i,k) * RH(i,k) / mu(i,k)
  Stk     = 1.156 * (Rek)**(-0.2)
  cf      = 2.0 * (utau(i,k) / Ue(i,k))**2
  St      = 0.5 * cf / ( Prt + sqrt(0.5 * cf / Stk) )
  hc(i,k) = St * rho(i,k) * Ue(i,k) * Cpa				! [J/(m^2*K*s)]
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine HeatTransferCoefficient3D
!*******************************************************************************************************
!******** 壁面粗さモデル (二次元，ShinとBondの経験式, 1992)					********
!*******************************************************************************************************
subroutine RoughnessShinBond2D( &
&            is, ie, Ts, Chord, Vfree, MVD, LWC, RH )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: Ts(is:ie)					! 表面温度
 real   , intent(in)  :: Chord, Vfree, MVD, LWC
 real   , intent(out) :: RH(is:ie)					! 粗さ高さ
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: k0 = 1.177 * 1.0e-3
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 RH(:) = 0.0
 do i = is + 1, ie - 1
  if( MVD <= 20.0e-6 )then
    RH(i) = 0.6839 * (0.571 + 0.246 * LWC + 1.257 * LWC**2) * (0.047 * Ts(i) - 11.27) &
    &       * ( 0.4286 + 0.0044139 * Vfree * 2.237 ) * k0 * Chord * 1.0
   else
    RH(i) = 0.6839 * (0.571 + 0.246 * LWC + 1.257 * LWC**2) * (0.047 * Ts(i) - 11.27) &
    &       * ( 0.4286 + 0.0044139 * Vfree * 2.237 ) * k0 * Chord * ( 1.667 - 0.00333 * MVD )
  end if
  if( RH(i) <= 0.0 ) RH(i) = k0
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine RoughnessShinBond2D
!*******************************************************************************************************
!******** 壁面粗さモデル (三次元，ShinとBondの経験式, 1992)					********
!*******************************************************************************************************
subroutine RoughnessShinBond3D( &
&            is, ie, ks, ke, Ts, LWC, RH )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: Ts(is:ie, ks:ke)				! 表面温度
 real   , intent(in)  :: LWC
 real   , intent(out) :: RH(is:ie, ks:ke)				! 粗さ高さ
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: k0 = 0.628e-3					! 物体の粗さ
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 ! 処理開始 ********************************************************************************************
 RH(:, :) = 0.0
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  RH(i,k) = k0 * 0.6839 * (0.047 * Ts(i,k) - 11.27) * (0.571 + 0.246 * LWC + 1.257 * LWC**2)
  if( RH(i,k) <= 0.0 ) RH(i,k) = k0
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine RoughnessShinBond3D
!*******************************************************************************************************
!******** 壁面粗さモデル (二次元，CIRAMIL code)							********
!*******************************************************************************************************
subroutine RoughnessCIRAMIL2D( &
&            is, ie, Ts, RH )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: Ts(is:ie)					! 表面温度
 real   , intent(out) :: RH(is:ie)					! 粗さ高さ
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 RH(:) = 0.0
 do i = is + 1, ie - 1
  if( Tf > Ts(i) ) then
    RH(i) = 2.82 * (Tf - Ts(i))**(-0.6) * 1.0e-3
   else
    RH(i) = 0.0
  endif
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine RoughnessCIRAMIL2D
!*******************************************************************************************************
!******** 壁面粗さモデル (三次元，CIRAMIL code)							********
!*******************************************************************************************************
subroutine RoughnessCIRAMIL3D( &
&            is, ie, ks, ke, Ts, RH )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: Ts(is:ie, ks:ke)				! 表面温度
 real   , intent(out) :: RH(is:ie, ks:ke)				! 粗さ高さ
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 ! 処理開始 ********************************************************************************************
 RH(:, :) = 0.0
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  if( Tf > Ts(i,k) ) then
    RH(i,k) = 2.82 * (Tf - Ts(i,k))**(-0.6) * 1.0e-3
   else
    RH(i,k) = 0.0
  endif
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine RoughnessCIRAMIL3D
!*******************************************************************************************************
!******** 壁面粗さモデル (二次元，Wrightらの式, 1997)						********
!*******************************************************************************************************
subroutine RoughnessWright2D( &
&            is, ie, x, y, SA, Up, Vp, Utau, Rho0, F, RH )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: x(is:ie), y(is:ie)				! 表面座標
 real   , intent(in)  :: SA(is:ie)					! 表面セル面積
 real   , intent(in)  :: Up(is:ie), Vp(is:ie)				! 壁面から一点目の速度
 real   , intent(in)  :: Utau(is:ie)					! 壁面摩擦速度
 real   , intent(in)  :: Rho0(is:ie)					! 翼表面密度
 real   , intent(in)  :: F						! 液滴の湿り率 (= LWC？)
 real   , intent(out) :: RH(is:ie)					! 粗さ高さ
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 real    :: ax, ay, az, bx, by, bz, nx, ny, na
 real    :: Uw, Vw, Velw
 real    :: SFx, SFy
 real    :: tau
 ! 処理開始 ********************************************************************************************
 RH(:) = 0.0
 do i = is + 1, ie - 1
  ! 壁面せん断力
  ax =  0.5 * ( - x(i-1) + x(i+1) )
  ay =  0.5 * ( - y(i-1) + y(i+1) )
  az =  0.0
  bx =  0.0
  by =  0.0
  bz = -1.0
  nx = ay * bz - az * by
  ny = az * bx - ax * bz
  na = sqrt(nx**2 + ny**2)
  nx = +1.0 * nx / na
  ny = +1.0 * ny / na
  Uw = Up(i) - (Up(i) * nx + Vp(i) * ny) * nx
  Vw = Vp(i) - (Up(i) * nx + Vp(i) * ny) * ny
  Velw = sqrt(Uw**2 + Vw**2)
  if(Velw /= 0) then
    SFx = Utau(i)**2 * Rho0(i) * SA(i) * Uw / Velw
    SFy = Utau(i)**2 * Rho0(i) * SA(i) * Vw / Velw
   else
    SFx = 0.0
    SFy = 0.0
  endif
  tau = sqrt(SFx**2 + SFy**2)
  ! 粗さ高さ
  if( tau > 0.0 ) then
    RH(i) = ( (4.0 * sigw * muw**2) / (Rhow * F * tau) )**(1.0/3.0) * 1.0e-3
   else
    RH(i) = 0.0
  endif
 enddo
 RH(is) = RH(is+1)
 RH(ie) = RH(ie-1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine RoughnessWright2D
!*******************************************************************************************************
!******** 壁面粗さモデル (三次元，Wrightらの式, 1997)						********
!*******************************************************************************************************
subroutine RoughnessWright3D( &
&            is, ie, ks, ke, x, y, z, SA, Up, Vp, Wp, Utau, Rho0, F, RH )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: x(is:ie, ks:ke), &
 &                       y(is:ie, ks:ke), &
 &                       z(is:ie, ks:ke)				! 表面座標
 real   , intent(in)  :: SA(is:ie, ks:ke)				! 表面セル面積
 real   , intent(in)  :: Up(is:ie, ks:ke), &
 &                       Vp(is:ie, ks:ke), &
 &                       Wp(is:ie, ks:ke)				! 壁面から一点目の速度
 real   , intent(in)  :: Utau(is:ie, ks:ke)				! 壁面摩擦速度
 real   , intent(in)  :: Rho0(is:ie, ks:ke)				! 翼表面密度
 real   , intent(in)  :: F						! 液滴の湿り率 (= LWC？)
 real   , intent(out) :: RH(is:ie, ks:ke)				! 粗さ高さ
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 real    :: ax, ay, az, bx, by, bz, nx, ny, nz, na
 real    :: Uw, Vw, Ww, Velw
 real    :: SFx, SFy, SFz
 real    :: tau
 ! 処理開始 ********************************************************************************************
 RH(:, :) = 0.0
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  ! 壁面せん断力
  ax = 0.5 * ( - x(i-1,k  ) + x(i+1,k  ) )
  ay = 0.5 * ( - y(i-1,k  ) + y(i+1,k  ) )
  az = 0.5 * ( - z(i-1,k  ) + z(i+1,k  ) )
  bx = 0.5 * ( + x(i  ,k-1) - x(i  ,k+1) )
  by = 0.5 * ( + y(i  ,k-1) - y(i  ,k+1) )
  bz = 0.5 * ( + z(i  ,k-1) - z(i  ,k+1) )
  nx = ay * bz - az * by
  ny = az * bx - ax * bz
  nz = ax * by - ay * bx
  na = sqrt(nx**2 + ny**2 + nz**2)
  nx = +1.0 * nx / na
  ny = +1.0 * ny / na
  nz = +1.0 * nz / na
  Uw = Up(i,k) - (Up(i,k) * nx + Vp(i,k) * ny + Wp(i,k) * nz) * nx
  Vw = Vp(i,k) - (Up(i,k) * nx + Vp(i,k) * ny + Wp(i,k) * nz) * ny
  Ww = Wp(i,k) - (Up(i,k) * nx + Vp(i,k) * ny + Wp(i,k) * nz) * nz
  Velw = sqrt(Uw**2 + Vw**2 + Ww**2)
  if(Velw /= 0) then
    SFx = Utau(i,k)**2 * Rho0(i,k) * SA(i,k) * Uw / Velw
    SFy = Utau(i,k)**2 * Rho0(i,k) * SA(i,k) * Vw / Velw
    SFz = Utau(i,k)**2 * Rho0(i,k) * SA(i,k) * Ww / Velw
   else
    SFx = 0.0
    SFy = 0.0
    SFz = 0.0
  endif
  tau = sqrt(SFx**2 + SFy**2 + SFz**2)
  ! 粗さ高さ
  if( tau > 0.0 ) then
    RH(i,k) = ( (4.0 * sigw * muw**2) / (Rhow * F * tau) )**(1.0/3.0) * 1.0e-3
   else
    RH(i,k) = 0.0
  endif
 enddo
 enddo
 do k = ks + 1, ke - 1
  RH(is,k) = RH(is+1,k)
  RH(ie,k) = RH(ie-1,k)
 enddo
 do i = is + 0, ie - 0
  RH(i,ks) = RH(i,ks+1)
  RH(i,ke) = RH(i,ke-1)
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine RoughnessWright3D
!*******************************************************************************************************
!******** 壁面粗さモデル (二次元，LEWICE 3.0, 1997)						********
!*******************************************************************************************************
subroutine RoughnessLEWICE2D( &
&            is, ie, FF, RH )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: FF(is:ie)					! 氷結率
 real   , intent(out) :: RH(is:ie)					! 粗さ高さ
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 do i = is, ie
!  RH(i) = 0.5 * sqrt( 0.15 + 0.3 / max(FF(i), zero) )
  RH(i) = 0.5 * sqrt( 0.15 + 0.3 / max(FF(i), zero) )
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine RoughnessLEWICE2D
!*******************************************************************************************************
!******** 壁面粗さモデル (三次元，LEWICE 3.0, 1997)						********
!*******************************************************************************************************
subroutine RoughnessLEWICE3D( &
&            is, ie, ks, ke, FF, RH )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: FF(is:ie, ks:ke)				! 氷結率
 real   , intent(out) :: RH(is:ie, ks:ke)				! 粗さ高さ
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 ! 処理開始 ********************************************************************************************
 do k = ks, ke
 do i = is, ie
  RH(i,k) = 0.5 * sqrt( 0.15 + 0.3 / max(FF(i,k), zero) )
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine RoughnessLEWICE3D
!*******************************************************************************************************
!******** 氷の密度 (二次元)									********
!******** Macklin 1962										********
!*******************************************************************************************************
subroutine IceDensity2D( &
&            is, ie, MVD, Uim, Vim, Ts, Rhoi )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: MVD
 real   , intent(in)  :: Uim, Vim
 real   , intent(in)  :: Ts
 real   , intent(out) :: Rhoi(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 real    :: R
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
  R = 0.5 * MVD * sqrt(Uim**2 + Vim**2) / Ts
  if(R <= 10.0)then
   Rhoi(i) = 110.0 * R**0.76
  else if(R > 60.0)then
   Rhoi(i) = 917.0
  else
   Rhoi(i) = R / (R + 5.61) * 10**3
  endif
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine IceDensity2D
!!*******************************************************************************************************
!!******** 氷の密度 (二次元)									********
!!******** !!! 密度モデル考える !!!								********
!!*******************************************************************************************************
!subroutine IceDensity2D( &
!&            is, ie, Tfree, Rhoi )
! ! 変数宣言 ********************************************************************************************
! implicit none
! ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! integer, intent(in)  :: is, ie
! real   , intent(in)  :: Tfree
! real   , intent(out) :: Rhoi(is:ie)
! ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! integer :: i
! ! 処理開始 ********************************************************************************************
! do i = is + 1, ie - 1
!  if( Tfree < Tf - 15.0 ) then
!    Rhoi(i) = Rhor
!   else
!    Rhoi(i) = Rhog
!  endif
! enddo
! ! 処理終了 ********************************************************************************************
! return
!end subroutine IceDensity2D
!*******************************************************************************************************
!******** 氷の密度 (三次元)									********
!******** !!! 密度モデル考える !!!								********
!*******************************************************************************************************
subroutine IceDensity3D( &
&            is, ie, ks, ke, Tfree, Rhoi )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: Tfree
 real   , intent(out) :: Rhoi(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  if( Tfree < Tf - 15.0 ) then
    Rhoi(i,k) = Rhor
   else
    Rhoi(i,k) = Rhog
  endif
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine IceDensity3D
!*******************************************************************************************************
!******** 表面温度 (二次元)									********
!*******************************************************************************************************
subroutine SurfaceTemperature2D( &
&            is, ie, Ufree, Tfree, Ue, Ta, Ts )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: Ufree, Tfree
 real   , intent(in)  :: Ue(is:ie)
 real   , intent(out) :: Ta(is:ie)					! 流体温度 (境界層外端)
 real   , intent(out) :: Ts(is:ie)					! 表面温度
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 real    :: Mach
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
  Ta(i) = Tfree
  Mach  = Ufree / sqrt(gamma * Rg * Ta(i))
  Ts(i) = Ta(i) + 0.5 * (Ufree**2 - Ue(i)**2) / Cpa &
  &             * (1.0 + 0.2 * rr * Mach**2) / (1.0 + 0.2 * Mach**2)
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine SurfaceTemperature2D
!*******************************************************************************************************
!******** 表面温度 (三次元)									********
!*******************************************************************************************************
subroutine SurfaceTemperature3D( &
&            is, ie, ks, ke, Ufree, Tfree, Ue, Ta, Ts )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: Ufree, Tfree
 real   , intent(in)  :: Ue(is:ie, ks:ke)
 real   , intent(out) :: Ta(is:ie, ks:ke)				! 流体温度 (境界層外端)
 real   , intent(out) :: Ts(is:ie, ks:ke)				! 表面温度
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 real    :: Mach
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  Ta(i,k) = Tfree
  Mach    = Ufree / sqrt(gamma * Rg * Ta(i,k))
  Ts(i,k) = Ta(i,k) + 0.5 * (Ufree**2 - Ue(i,k)**2) / Cpa &
  &                 * (1.0 + 0.2 * rr * Mach**2) / (1.0 + 0.2 * Mach**2)
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine SurfaceTemperature3D
!*******************************************************************************************************
!******** 表面温度の更新 (二次元)									********
!*******************************************************************************************************
subroutine RenewSurfaceTemperature2D( &
&            is, ie, Ts, Bi, Ti, Twall )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie
 real   , intent(in)    :: Bi(is:ie)					! 氷層厚さ
 real   , intent(in)    :: Ti(is:ie)					! 氷層温度
 real   , intent(inout) :: Ts(is:ie)					! 表面温度
 real   , intent(in)    :: Twall(is:ie)					! 翼面温度
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 ! 着氷無し⇒翼面温度(氷層が既に存在する場合は更新しない); 氷成長あり⇒氷層温度
 do i = is, ie
  if ( Bi(i) > 0.0 ) then
    Ts(i) = Ti(i)
   else
    Ts(i) = Twall(i)
  end if
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine RenewSurfaceTemperature2D
!*******************************************************************************************************
!******** 衝突液滴の質量・エンタルピ (二次元, Original Messinger)				********
!*******************************************************************************************************
subroutine DropImpingementOrg2D( &
&            is, ie, MVD, Rhod, Tfree, SA, CE, Uim, Vim, Mim, Him )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: MVD, Rhod, Tfree
 real   , intent(in)  :: SA(is:ie)
 real   , intent(in)  :: CE(is:ie)
 real   , intent(in)  :: Uim(is:ie), Vim(is:ie)
 real   , intent(out) :: Mim(is:ie), Him(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
!  Mim(i) = Rhod * CE(i) * 0.5 * pi * (0.5 * MVD)**2 * SA(i)
  Mim(i) = Rhod * CE(i) * 4.0 / 3.0 * pi * (0.5 * MVD)**3 * SA(i)
  Him(i) = Cpw * (Tfree - Tf) + 0.5 * (Uim(i)**2 + Vim(i)**2)
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine DropImpingementOrg2D
!*******************************************************************************************************
!******** 衝突液滴の質量・エンタルピ (三次元, Original Messinger)				********
!*******************************************************************************************************
subroutine DropImpingementOrg3D( &
&            is, ie, ks, ke, MVD, Rhod, Tfree, SA, CE, Uim, Vim, Wim, Mim, Him )
 ! 変数宣言 ********************************************************************************************
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: MVD, Rhod, Tfree
 real   , intent(in)  :: SA(is:ie, ks:ke)
 real   , intent(in)  :: CE(is:ie, ks:ke)
 real   , intent(in)  :: Uim(is:ie, ks:ke), Vim(is:ie, ks:ke), Wim(is:ie, ks:ke)
 real   , intent(out) :: Mim(is:ie, ks:ke), Him(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  Mim(i,k) = Rhod * CE(i,k) * 4.0 / 3.0 * pi * (0.5 * MVD)**3 * SA(i,k)
  Him(i,k) = Cpw * (Tfree - Tf) + 0.5 * (Uim(i,k)**2 + Vim(i,k)**2 + Wim(i,k)**2)
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine DropImpingementOrg3D
!*******************************************************************************************************
!******** 衝突液滴の質量流束 (二次元, Extended Messinger)					********
!*******************************************************************************************************
subroutine DropImpingementExt2D( &
&            is, ie, MVD, Rhod, CE, Mim )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: MVD, Rhod
 real   , intent(in)  :: CE(is:ie)
 real   , intent(out) :: Mim(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
!  Mim(i) = Rhod * CE(i) * 0.5 * pi * (0.5 * MVD)**2
  Mim(i) = Rhod * CE(i) * 4.0 / 3.0 * pi * (0.5 * MVD)**3
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine DropImpingementExt2D
!*******************************************************************************************************
!******** 衝突液滴の質量流束 (三次元, Extended Messinger)					********
!*******************************************************************************************************
subroutine DropImpingementExt3D( &
&            is, ie, ks, ke, MVD, Rhod, CE, Mim )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: MVD, Rhod
 real   , intent(in)  :: CE(is:ie, ks:ke)
 real   , intent(out) :: Mim(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  Mim(i,k) = Rhod * CE(i,k) * 4.0 / 3.0 * pi * (0.5 * MVD)**3
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine DropImpingementExt3D
!*******************************************************************************************************
!******** 蒸発・昇華の質量・エンタルピ (二次元, Original Messinger)				********
!*******************************************************************************************************
subroutine EvapolationOrg2D( &
&            is, ie, Tfree, Ts, Pt, hc, SA, dBi, Mes, Hes )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: Tfree
 real   , intent(in)  :: Ts(is:ie), Pt(is:ie)
 real   , intent(in)  :: hc(is:ie), SA(is:ie)
 real   , intent(in)  :: dBi(is:ie)
 real   , intent(out) :: Mes(is:ie), Hes(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 real    :: Tv, Tsur,  Pvsur, Pvinf
 real    :: kai, Qes
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
  Tsur = Ts(i)
  kai  = 0.622 * hc(i) * LS / ( Cpa * Pt(i) * (1.0 / Pr)**(2.0 / 3.0) )
  Qes  = kai * e0 * (Tsur - Tfree)
  Mes(i) = Qes / LS * SA(i)
  if( dBi(i) > 0.0 ) then
    Hes(i) = LS + Cpi * (Tsur - Tf) - LF
   else
    Hes(i) = LE + Cpw * (Tsur - Tf)
  endif
  if( Mes(i) < 0.0 ) then
    Mes(i) = 0.0
    Hes(i) = 0.0
  endif
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine EvapolationOrg2D
!*******************************************************************************************************
!******** 蒸発・昇華の質量・エンタルピ (三次元, Original Messinger)				********
!*******************************************************************************************************
subroutine EvapolationOrg3D( &
&            is, ie, ks, ke, Tfree, Ts, Pt, hc, SA, dBi, Mes, Hes )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: Tfree
 real   , intent(in)  :: Ts(is:ie, ks:ke), Pt(is:ie, ks:ke)
 real   , intent(in)  :: hc(is:ie, ks:ke), SA(is:ie, ks:ke)
 real   , intent(in)  :: dBi(is:ie, ks:ke)
 real   , intent(out) :: Mes(is:ie, ks:ke), Hes(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 real    :: Tv, Tsur,  Pvsur, Pvinf
 real    :: kai, Qes
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  Tsur = Ts(i,k)
  kai  = 0.622 * hc(i,k) * LS / ( Cpa * Pt(i,k) * (1.0 / Pr)**(2.0 / 3.0) )
  Qes  = kai * e0 * (Tsur - Tfree)
  Mes(i,k) = Qes / LS * SA(i,k)
  if( dBi(i,k) > 0.0 ) then
    Hes(i,k) = LS + Cpi * (Tsur - Tf) - LF
   else
    Hes(i,k) = LE + Cpw * (Tsur - Tf)
  endif
  if( Mes(i,k) < 0.0 ) then
    Mes(i,k) = 0.0
    Hes(i,k) = 0.0
  endif
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine EvapolationOrg3D
!*******************************************************************************************************
!******** 蒸発・昇華の質量流束 (二次元, Extended Messinger)					********
!*******************************************************************************************************
subroutine EvapolationExt2D( &
&            is, ie, Tfree, Pfree, Ta, Ts, Ti, Tw, Pt, hc, nPhase, Mes )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: Tfree, Pfree
 real   , intent(in)  :: Ts(is:ie), Ta(is:ie), Ti(is:ie), Tw(is:ie), Pt(is:ie)
 real   , intent(in)  :: hc(is:ie)
 integer, intent(in)  :: nPhase(is:ie)
 real   , intent(out) :: Mes(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 real    :: Tsur, Tbsur, Tbinf						! 表面温度，Tbar
 real    :: Pvsur, Pvinf						! 表面と周囲の蒸気圧
 real    :: kai, Qes
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
  select case(nPhase(i))
 ! 着氷なし --------------------------------------------------------------------------------------------
   case(0)
    Tsur = 0.0; Mes(i) = 0.0
 ! 霧氷 ------------------------------------------------------------------------------------------------
 ! 昇華
   case(1)
    Tsur = Ti(i)
 ! 雨氷 or 防氷 ----------------------------------------------------------------------------------------
 ! 蒸発
   case(2,3)
    Tsur = Tw(i)
   case default; write(*, '(a)') '!!!!! Error : Ice phase !!!!!'
  end select
  ! 表面温度，周囲温度から蒸気圧を計算
  Tbsur = 72.0 + 1.8 * (Tsur - 273.15)
  Pvsur = 3386.0 * (0.0039 + 6.8096 * 10**(-6.0) * Tbsur**2 + 3.5579 * 10**(-7.0) * Tbsur**3)
  Tbinf = 72.0 + 1.8 * (Tfree - 273.15)
  Pvinf = 3386.0 * (0.0039 + 6.8096 * 10**(-6.0) * Tbinf**2 + 3.5579 * 10**(-7.0) * Tbinf**3)
  if(Pvsur > Pvinf) then
   Mes(i) = 0.7 / Cpa * hc(i) * (Pvsur - Pvinf) / Pfree
  else
   Mes(i) = 0.0
  endif
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine EvapolationExt2D
!*******************************************************************************************************
!******** 蒸発・昇華の質量流束 (三次元, Extended Messinger)					********
!*******************************************************************************************************
subroutine EvapolationExt3D( &
&            is, ie, ks, ke, Ta, Ti, Tw, Pt, hc, nPhase, Mes )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: Ta(is:ie, ks:ke), Ti(is:ie, ks:ke), Tw(is:ie, ks:ke), Pt(is:ie, ks:ke)
 real   , intent(in)  :: hc(is:ie, ks:ke)
 integer, intent(in)  :: nPhase(is:ie, ks:ke)
 real   , intent(out) :: Mes(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 real    :: Tsur, kai, Qes
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  select case(nPhase(i,k))
   case(0,1)
    Tsur = Ti(i,k)
    kai  = 0.622 * hc(i,k) * LS / ( Cpa * Pt(i,k) * (1.0 / Pr)**(2.0 / 3.0) )
    Qes  = kai * e0 * (Tsur - Ta(i,k))
    Mes(i,k) = Qes / LS
    if( Mes(i,k) < 0.0 ) Mes(i,k) = 0.0
   case(2)
    Tsur = Tw(i,k)
    kai  = 0.622 * hc(i,k) * LE / ( Cpa * Pt(i,k) * (1.0 / Pr)**(2.0 / 3.0) )
    Qes  = kai * e0 * (Tsur - Ta(i,k))
    Mes(i,k) = Qes / LE
    if( Mes(i,k) < 0.0 ) Mes(i,k) = 0.0
   case default; write(*, '(a)') '!!!!! Error : Ice phase !!!!!'
  end select
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine EvapolationExt3D
!*******************************************************************************************************
!******** 熱の収支 (二次元，Original Messinger, 摩擦熱・対流熱，他は無視)			********
!*******************************************************************************************************
subroutine HeatBalanceOrg2D( &
&            is, ie, Tfree, Ts, hc, SA, utau, Q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: Tfree
 real   , intent(in)  :: Ts(is:ie)
 real   , intent(in)  :: hc(is:ie), SA(is:ie), utau(is:ie)
 real   , intent(out) :: Q(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 real    :: Qf, Qc
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
  ! 摩擦熱
  Qf = hc(i) * SA(i) * rr * 0.5 * Utau(i)**2 / Cpa
  ! 対流熱
  Qc = hc(i) * SA(i) * (Ts(i) - Tfree)
  ! 考慮する熱の項
  Q(i) = Qf - Qc
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine HeatBalanceOrg2D
!*******************************************************************************************************
!******** 熱の収支 (三次元，Original Messinger, 摩擦熱・対流熱，他は無視)			********
!*******************************************************************************************************
subroutine HeatBalanceOrg3D( &
&            is, ie, ks, ke, Tfree, Ts, hc, SA, utau, Q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: Tfree
 real   , intent(in)  :: Ts(is:ie, ks:ke)
 real   , intent(in)  :: hc(is:ie, ks:ke), SA(is:ie, ks:ke), utau(is:ie, ks:ke)
 real   , intent(out) :: Q(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 real    :: Qf, Qc
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  ! 摩擦熱
  Qf = hc(i,k) * SA(i,k) * rr * 0.5 * Utau(i,k)**2 / Cpa
  ! 対流熱
  Qc = hc(i,k) * SA(i,k) * (Ts(i,k) - Tfree)
  ! 考慮する熱の項
  Q(i,k) = Qf - Qc
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine HeatBalanceOrg3D
!*******************************************************************************************************
!******** 熱の収支 (二次元，Extended Messinger)							********
!*******************************************************************************************************
subroutine HeatBalanceExt2D( &
&            is, ie, Tfree, Ta, Tw, Tw0, Pt, Ufree, hc, dBi, Mim, Mrin, Mst, Mes, nPhase, Q0, Q1 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: Tfree, Ufree
 real   , intent(in)  :: Ta(is:ie), Tw(is:ie), Tw0(is:ie), Pt(is:ie)
 real   , intent(in)  :: hc(is:ie), dBi(is:ie)
 real   , intent(in)  :: Mim(is:ie), Mrin(is:ie), Mst(is:ie), Mes(is:ie)
 integer, intent(in)  :: nPhase(is:ie)
 real   , intent(inout) :: Q0(is:ie), Q1(is:ie) !Q0:熱量IN,Q1:熱量OUT
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 real    :: kai, T_rin
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
!write(*,*)'* Mim(', i, ') = ', Mim(i)
!write(*,*)'* Mrin(', i, ') = ', Mrin(i)
!write(*,*)'* Mst(', i, ') = ', Mst(i)
!write(*,*)'* Mes(', i, ') = ', Mes(i)
  select case(nPhase(i))
 ! 着氷なし or 霧氷 ------------------------------------------------------------------------------------
 ! 昇華
   case(0,1)
!    Q0(i) = Q0(i) + rhor * LF * dBi(i)
!    Q1(i) = Q1(i)
    kai  = 0.622 * hc(i) * LS / ( Cpa * Pt(i) * (1.0 / Pr)**(2.0 / 3.0) )
    Q0(i) = rhor * LF * dBi(i) &
    &     + 0.5 * Mim(i) * Ufree**2 &
    &     + 0.5 * rr * hc(i) * Ufree**2 / Cpa &
    &     + Mim(i) * Cpw * Ta(i) &
    &     + hc(i) * Ta(i) & 
    &     + 4.0 * eps * sigr * Ta(i)**4 &
    &     + kai * e0 * Ta(i) & 
    &     + ( Mrin(i) + Mst(i) ) * Cpw * Tf
    Q1(i) = Mim(i) * Cpw &
    &     + hc(i) &
    &     + 4.0 * eps * sigr * Tfree**3 &
    &     + kai * e0 &
    &     + ( Mrin(i) + Mst(i) ) * Cpw
 ! 雨氷 ------------------------------------------------------------------------------------------------
 ! 蒸発
   case(2)
    Q0(i) = Q0(i)
    Q1(i) = Q1(i)
!    kai  = 0.622 * hc(i) * LE / ( Cpa * Pt(i) * (1.0 / Pr)**(2.0 / 3.0) )
!    Q0(i) = 0.5 * Mim(i) * Ufree**2 &
!    &     + 0.5 * rr * hc(i) * Ufree**2 / Cpa &
!    &     + Mim(i) * Cpw * Ta(i) &
!    &     + hc(i) * Ta(i) & 
!    &     + 4.0 * eps * sigr * Ta(i)**4 &
!    &     + kai * e0 * Ta(i) & 
!    &     + ( Mrin(i) + Mst(i) ) * Cpw * Tf
!    Q1(i) = Mim(i) * Cpw &
!    &     + hc(i) &
!    &     + 4.0 * eps * sigr * Tfree**3 &
!    &     + kai * e0 &
!    &     + ( Mrin(i) + Mst(i) ) * Cpw
   case default; write(*, '(a)') '!!!!! Error : Ice phase !!!!!'
  end select
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine HeatBalanceExt2D
!*******************************************************************************************************
!******** 熱の収支 (二次元，Extended Messinger, 翼面加熱を考慮)					********
!*******************************************************************************************************
subroutine HeatBalanceExt2DverUR( &
&          is, ie, nPhase, Tfree, Ufree, Pt, Tc, Tw, hc, Mim, Mrin, Mst, Mes, Qin, Qout )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, nPhase(is:ie)
 real   , intent(in)    :: Tfree, Ufree
 real   , intent(in)    :: Pt(is:ie), Tc(is:ie), Tw(is:ie), hc(is:ie)
 real   , intent(in)    :: Mim(is:ie), Mrin(is:ie), Mst(is:ie), Mes(is:ie)
 real   , intent(inout) :: Qin(is:ie), Qout(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 real    :: kai, T_rin
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
  select case(nPhase(i))
   ! 着氷なし ------------------------------------------------------------------------------------------
   case(0)
   Qin(i) = 0.0; Qout(i) = 0.0
   ! 霧氷 ----------------------------------------------------------------------------------------------
   case(1)
    kai     = 0.622 * hc(i) * LS / ( Cpa * Pt(i) * (1.0 / Pr)**(2.0 / 3.0) )
    Qin(i)  = Qin(i) + ( Mrin(i) + Mim(i) + Mst(i) ) * LF ! Q_latentを追加
    Qout(i) =  hc(i) * ( Tc(i) - Tfree )            &
    &        + kai * e0 * ( Tc(i) - Tfree )         & ! 昇華熱で再計算
    &        + Mim(i) * Cpw * ( Tc(i) - Tfree )    &
    &        + eps * sigr * ( Tc(i)**4 - Tfree**4 )
   ! 雨氷  ---------------------------------------------------------------------------------------------
   case(2)
    kai     = 0.622 * hc(i) * LE / ( Cpa * Pt(i) * (1.0 / Pr)**(2.0 / 3.0) )
    Qin(i) =  0.5 * Mim(i) * Ufree**2 &
    &       + 0.5 * rr * hc(i) * Ufree**2 / Cpa &
    &       + ( Mrin(i) + Mim(i) + Mst(i) ) * ( Tc(i) - Tw(i) )
    Qout(i) =  hc(i) * ( Tc(i) - Tfree )            &
    &        + kai * e0 * ( Tc(i) - Tfree )         &
    &        + eps * sigr * ( Tc(i)**4 - Tfree**4 )
   ! 防氷 ----------------------------------------------------------------------------------------------
   case(3)
    ! 値の更新はしない
    Qin(i)  = Qin(i)
    Qout(i) = Qout(i)
   case default; write(*, '(a)') '!!!!! Error : Ice phase !!!!!'
  end select
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine HeatBalanceExt2DverUR
!*******************************************************************************************************
!******** 熱の収支 (三次元，Extended Messinger)							********
!*******************************************************************************************************
subroutine HeatBalanceExt3D( &
&            is, ie, ks, ke, Tfree, Ta, Pt, utau, Uim, Vim, Wim, &
&            hc, Mim, Mrin, Mst, nPhase, Q0, Q1 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: Tfree
 real   , intent(in)  :: Ta(is:ie, ks:ke), Pt(is:ie, ks:ke)
 real   , intent(in)  :: utau(is:ie, ks:ke)
 real   , intent(in)  :: Uim(is:ie, ks:ke), Vim(is:ie, ks:ke), Wim(is:ie, ks:ke)
 real   , intent(in)  :: hc(is:ie, ks:ke)
 real   , intent(in)  :: Mim(is:ie, ks:ke), Mrin(is:ie, ks:ke), Mst(is:ie, ks:ke)
 integer, intent(in)  :: nPhase(is:ie, ks:ke)
 real   , intent(out) :: Q0(is:ie, ks:ke), Q1(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 real    :: kai
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  select case(nPhase(i,k))
 ! 着氷なし or 霧氷 ------------------------------------------------------------------------------------
 ! 昇華
   case(0,1)
    kai  = 0.622 * hc(i,k) * LS / ( Cpa * Pt(i,k) * (1.0 / Pr)**(2.0 / 3.0) )
 ! 雨氷 ------------------------------------------------------------------------------------------------
 ! 蒸発
   case(2)
    kai  = 0.622 * hc(i,k) * LE / ( Cpa * Pt(i,k) * (1.0 / Pr)**(2.0 / 3.0) )
   case default; write(*, '(a)') '!!!!! Error : Ice phase !!!!!'
  end select
  Q0(i,k) = 0.5 * Mim(i,k) * (Uim(i,k)**2 + Vim(i,k)**2 + Wim(i,k)**2) &
  &       + 0.5 * rr * hc(i,k) * Utau(i,k)**2 / Cpa &
  &       + Mim(i,k) * Cpw * Ta(i,k) + hc(i,k) * Ta(i,k) + 4.0 * eps * sigr * Ta(i,k)**4 &
  &       + kai * e0 * Ta(i,k) + Mrin(i,k) * Cpw * Tf
!  &       + kai * e0 * Tfree + Mrin(i,k) * Cpw * Tf			! 結構変わる
!  &       + kai * e0 * Tfree + (Mrin(i,k) + Mst(i,k)) * Cpw * Tf
!  Q1(i,k) = Mim(i,k) * Cpw + hc(i,k) + 4.0 * eps * sigr * Ta(i,k)**3 &	! 若干変わる
  Q1(i,k) = Mim(i,k) * Cpw + hc(i,k) + 4.0 * eps * sigr * Tfree**3 &
  &       + kai * e0 + Mrin(i,k) * Cpw
!  &       + kai * e0 + (Mrin(i,k) + Mst(i,k)) * Cpw
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine HeatBalanceExt3D
!*******************************************************************************************************
!******** 氷結率 (二次元，Original Messinger)							********
!*******************************************************************************************************
subroutine FreezingFractionOrg2D( &
&            is, ie, Mim, Mes, Mrin, Mst, Him, Hes, Hrin, Hrout, Hac, Q, FF )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: Mim(is:ie), Mes(is:ie), Mrin(is:ie), Mst(is:ie)
 real   , intent(in)  :: Him(is:ie), Hes(is:ie), Hrin(is:ie), Hrout(is:ie), Hac(is:ie)
 real   , intent(in)  :: Q(is:ie)
 real   , intent(out) :: FF(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: ff0(:, :)
 integer :: i, l
 ! 処理開始 ********************************************************************************************
 ! 氷結率 (0 ; 氷結なし = 全てRunback，1 ; 全て氷結 = Runbackなし) +++++++++++++++++++++++++++++++++++++
 do i = is + 1, ie - 1
!  if( (Mim(i) + Mrin(i)) * (Hac(i) - Hrout(i)) /= 0.0 ) then
!    ff(i) = ( Mim(i) * Him(i) + Mrin(i) * Hrin(i) + Q(i) &
!    &         - Mes(i) * Hes(i) - (Mim(i) + Mrin(i) - Mes(i)) * Hrout(i)  ) &
!    &       / ( (Mim(i) + Mrin(i)) * (Hac(i) - Hrout(i)) )
  if( (Mim(i) + Mrin(i) + Mst(i)) * (Hac(i) - Hrout(i)) /= 0.0 ) then
    ff(i) = ( Mim(i) * Him(i) + (Mrin(i) + Mst(i)) * Hrin(i) + Q(i) &
    &         - Mes(i) * Hes(i) - (Mim(i) + Mrin(i) + Mst(i) - Mes(i)) * Hrout(i)  ) &
    &       / ( (Mim(i) + Mrin(i) + Mst(i)) * (Hac(i) - Hrout(i)) )
   else
    ff(i) = 1.0
  endif
 enddo
! ! スムージング ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! allocate( ff0(is:ie) )
! do l = 1, 5
!  ff0(:, :) = ff(:, :)
!  do i = is + 1, ie - 1
!   ff(i) = (ff0(i-1) + 2.0 * ff0(i) + ff0(i+1)) * 0.25
!  enddo
! enddo
! deallocate(ff0)
 ! リミッタ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do i = is + 1, ie - 1
  if( ff(i) < 0.0 ) then
    ff(i) = 0.0
   else if( ff(i) <= 1.0 ) then
    cycle
   else 
    ff(i) = 1.0
  endif
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine FreezingFractionOrg2D
!*******************************************************************************************************
!******** 氷結率 (三次元，Original Messinger)							********
!*******************************************************************************************************
subroutine FreezingFractionOrg3D( &
&            is, ie, ks, ke, Mim, Mes, Mrin, Mst, Him, Hes, Hrin, Hrout, Hac, Q, FF )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: Mim(is:ie, ks:ke), Mes(is:ie, ks:ke), Mrin(is:ie, ks:ke), Mst(is:ie, ks:ke)
 real   , intent(in)  :: Him(is:ie, ks:ke), Hes(is:ie, ks:ke), &
 &                       Hrin(is:ie, ks:ke), Hrout(is:ie, ks:ke), Hac(is:ie, ks:ke)
 real   , intent(in)  :: Q(is:ie, ks:ke)
 real   , intent(out) :: FF(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: ff0(:, :)
 integer :: i, k, l
 ! 処理開始 ********************************************************************************************
 ! 氷結率 (0 : 氷結なし = 全てRunback，1 : 全て氷結 = Runbackなし) +++++++++++++++++++++++++++++++++++++
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
!  if( (Mim(i,k) + Mrin(i,k)) * (Hac(i,k) - Hrout(i,k)) /= 0.0 ) then
!    ff(i,k) = ( Mim(i,k) * Him(i,k) + Mrin(i,k) * Hrin(i,k) + Q(i,k) &
!    &         - Mes(i,k) * Hes(i,k) - (Mim(i,k) + Mrin(i,k) - Mes(i,k)) * Hrout(i,k)  ) &
!    &       / ( (Mim(i,k) + Mrin(i,k)) * (Hac(i,k) - Hrout(i,k)) )
  if( (Mim(i,k) + Mrin(i,k) + Mst(i,k)) * (Hac(i,k) - Hrout(i,k)) /= 0.0 ) then
    ff(i,k) = ( Mim(i,k) * Him(i,k) + (Mrin(i,k) + Mst(i,k)) * Hrin(i,k) + Q(i,k) &
    &         - Mes(i,k) * Hes(i,k) - (Mim(i,k) + Mrin(i,k) + Mst(i,k) - Mes(i,k)) * Hrout(i,k)  ) &
    &       / ( (Mim(i,k) + Mrin(i,k) + Mst(i,k)) * (Hac(i,k) - Hrout(i,k)) )
   else
    ff(i,k) = 1.0
  endif
 enddo
 enddo
! ! スムージング ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! allocate( ff0(is:ie, ks:ke) )
! do l = 1, 5
!  ff0(:, :) = ff(:, :)
!  do k = ks + 1, ke - 1
!  do i = is + 1, ie - 1
!   ff(i,k) = (ff0(i-1,k) + 2.0 * ff0(i,k) + ff0(i+1,k)) * 0.25
!  enddo
!  enddo
! enddo
! deallocate(ff0)
 ! リミッタ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  if( ff(i,k) < 0.0 ) then
    ff(i,k) = 0.0
   else if( ff(i,k) <= 1.0 ) then
    cycle
   else 
    ff(i,k) = 1.0
  endif
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine FreezingFractionOrg3D
!*******************************************************************************************************
!******** 氷結率 (二次元，Extended Messinger)							********
!*******************************************************************************************************
subroutine FreezingFractionExt2D( &
&            is, ie, time, Mim, Mes, Mrin, Mst, nPhase, Bi, Bw, Bg, dBi, FF )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: time
 real   , intent(in)  :: Mim(is:ie), Mes(is:ie), Mrin(is:ie), Mst(is:ie)
 integer, intent(in)  :: nPhase(is:ie)
 real   , intent(in)  :: Bi(is:ie), Bw(is:ie), Bg(is:ie)
 real   , intent(in)  :: dBi(is:ie)
 real   , intent(out) :: FF(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 real    :: kai
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
  ! 氷結率 (0 = 氷結なし，1 = 全て氷結) ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  select case(nPhase(i))
   ! 着氷なし ------------------------------------------------------------------------------------------
   case(0)
    FF(i) = 0.0
   ! 霧氷 ----------------------------------------------------------------------------------------------
   case(1)
    if( ((Mim(i) + Mrin(i) + Mst(i)) * time) > 0.0 ) then
!      FF(i) = Rhor * Bi(i) / ( (Mim(i) + Mrin(i) + Mst(i)) * time )
      FF(i) = Rhor * dBi(i) / ( Mim(i) + Mrin(i) + Mst(i) )
     else
      FF(i) = 0.0
    endif
   ! 雨氷 ----------------------------------------------------------------------------------------------
   case(2)
    if( ((Mim(i) + Mrin(i) + Mst(i)) * time) > 0.0 ) then
!      FF(i) = (rhor * Bg(i) + rhog * (Bi(i) - Bg(i))) / ((Mim(i) + Mrin(i) + Mst(i)) * time)
      FF(i) = rhog * dBi(i) / ( Mim(i) + Mrin(i) + Mst(i) )
     else
      FF(i) = 0.0
    endif
   ! 防氷 ----------------------------------------------------------------------------------------------
   case(3)
    FF(i) = 0.0
   ! 例外処理 ------------------------------------------------------------------------------------------
   case default
    write(*, '(a)') '!!!!! Error : Ice phase !!!!!'
  end select
  ! 範囲外 (0～1)の処理 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( FF(i) < 0.0 ) then
    FF(i) = 0.0
   else if( FF(i) >= 1.0 ) then
    FF(i) = 1.0
  endif
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine FreezingFractionExt2D
!*******************************************************************************************************
!******** 氷結率 (三次元，Extended Messinger)							********
!*******************************************************************************************************
subroutine FreezingFractionExt3D( &
&            is, ie, ks, ke, time, Mim, Mes, Mrin, Mst, nPhase, Bi, Bw, FF )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: time
 real   , intent(in)  :: Mim(is:ie, ks:ke), Mes(is:ie, ks:ke), Mrin(is:ie, ks:ke), Mst(is:ie, ks:ke)
 integer, intent(in)  :: nPhase(is:ie, ks:ke)
 real   , intent(in)  :: Bi(is:ie, ks:ke), Bw(is:ie, ks:ke)
 real   , intent(out) :: FF(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 real    :: kai
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  ! 氷結率 (0 = 氷結なし，1 = 全て氷結) ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  select case(nPhase(i,k))
   ! 着氷なし ------------------------------------------------------------------------------------------
   case(0)
    FF(i,k) = 1.0
   ! 霧氷 ----------------------------------------------------------------------------------------------
   case(1)
    if( (Mim(i,k) + Mrin(i,k) + Mst(i,k)) * time > 0.0 ) then
      FF(i,k) = Rhor * Bi(i,k) / ( (Mim(i,k) + Mrin(i,k) + Mst(i,k)) * time )
!K    if( (Mim(i,k) + Mrin(i,k) + Mst(i,k)) > 0.0 ) then
!K      FF(i,k) = Rhor * dBi(i,k) / (Mim(i,k) + Mrin(i,k) + Mst(i,k))
     else
      FF(i,k) = 1.0
    endif
   ! 雨氷 ----------------------------------------------------------------------------------------------
   case(2)
    if( (Mim(i,k) + Mrin(i,k) + Mst(i,k)) * time > 0.0 ) then
      FF(i,k) = (Rhog * Bi(i,k) + Rhow * Bw(i,k)) / ( (Mim(i,k) + Mrin(i,k) + Mst(i,k)) * time )
!      FF(i,k) = (Rhog * Bi(i,k) + Rhow * Bw(i,k) * FF(I,K)) / ( (Mim(i,k) + Mrin(i,k) + Mst(i,k)) * time )
!      FF(i,k) = Rhog * Bi(i,k) / ( (Mim(i,k) + Mrin(i,k) + Mst(i,k)) * time )
!    if( Mim(i,k) + Mrin(i,k) + Mst(i,k) > 0.0 ) then
!      FF(i,k) = Rhog * dBi(i,k) / ( Mim(i,k) + Mrin(i,k) + Mst(i,k) )
!    if( (Mim(i,k) + Mrin(i,k) + Mst(i,k)) * time > 0.0 ) then
!      FF(i,k) = (Rhor * Bg(i,k) + Rhog * (Bi(i,k) - Bg(i,k))) &
!      &       / ( (Mim(i,k) + Mrin(i,k) + Mst(i,k)) * time )
!K    if( (Mim(i,k) + Mrin(i,k) + Mst(i,k)) > 0.0 ) then
!K      FF(i,k) = Rhog * dBi(i,k) / (Mim(i,k) + Mrin(i,k) + Mst(i,k))
     else
      FF(i,k) = 1.0
    endif
   ! 例外処理 ------------------------------------------------------------------------------------------
   case default
    write(*, '(a)') '!!!!! Error : Ice phase !!!!!'
  end select
  ! 範囲外 (0～1)の処理 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( FF(i,k) < 0.0 ) then
    FF(i,k) = 0.0
   else if( FF(i,k) <= 1.0 ) then
    cycle
   else 
    FF(i,k) = 1.0
  endif
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine FreezingFractionExt3D
!*******************************************************************************************************
!******** Runback-out の質量・エンタルピ (二次元，Original Messinger)				********
!*******************************************************************************************************
subroutine RunbackOutOrg2D( &
&            is, ie, RunBackNum, dt, Ts, SA, utau, RH, Mim, Mes, Mrin, FF, Mrout, Hrout, Mst )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie
 integer, intent(in)    :: RunbackNum
 real   , intent(in)    :: dt
 real   , intent(in)    :: Ts(is:ie)
 real   , intent(in)    :: SA(is:ie), utau(is:ie), RH(is:ie)
 real   , intent(in)    :: Mim(is:ie), Mes(is:ie), Mrin(is:ie)
 real   , intent(in)    :: FF(is:ie)
 real   , intent(out)   :: Mrout(is:ie), Hrout(is:ie)
 real   , intent(inout) ::  Mst(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: Mout0(:), Bw(:)
 integer :: i
 real    :: MFR
 real    :: Aout, Ast
 ! 処理開始 ********************************************************************************************
 allocate( Bw(is:ie) )
 Mrout(:) = 0.0; Bw(:) = 0.0
 do i = is + 1, ie - 1
  ! Runback-out 質量流束 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Mrout(i) = (1.0 - ff(i)) * (Mim(i) + Mrin(i) + Mst(i)) - Mes(i)
  if( Mrout(i) <= 0.0 ) then
    Mrout(i) = 0.0
    Mst(i)   = 0.0
   else
    select case( RunbackNum )
    ! 全てランバック -----------------------------------------------------------------------------------
     case(0)
      Mst(i) = 0.0
     case(1)
      ! 水膜の厚さ
      Bw(i) = ( Mrout(i) + Mst(i) ) / Rhow / SA(i) * dt
      ! 面積充填率 (氷表面を半円と仮定)
      MFR  = 0.5 * 0.5 * pi
      Ast  = (1.0 - MFR) * 2.0 * RH(i)**2 / dt
      Aout = Mrout(i) / Rhow * Bw(i) / SA(i)
      if( Bw(i) <= 0.0 ) then
        Mst  (i) = Mrout(i)
        Mrout(i) = 0.0
       else
        if( Ast < Aout ) then
          Mst  (i) = Ast * Rhow / Bw(i) * SA(i)
          Mrout(i) = Mrout(i) - Mst(i)
         else
          Mst  (i) = Mrout(i)
          Mrout(i) = 0.0
        endif
       endif
     case default; write(*, '(a)') '!!!!! Error : Runback model number !!!!!'
    end select
  endif
  if( Mrout(i) > 0.0 ) then
    Hrout(i) = Cpw * (Ts(i) - Tf) + 0.5 * Utau(i)**2
   else
    Hrout(i) = 0.0
  endif
 enddo
 deallocate(Bw)
 ! 処理終了 ********************************************************************************************
 return
end subroutine RunbackOutOrg2D
!*******************************************************************************************************
!******** Runback-out の質量・エンタルピ (三次元，Original Messinger)				********
!*******************************************************************************************************
subroutine RunbackOutOrg3D( &
&            is, ie, ks, ke, RunbackNum, dt, Ts, SA, utau, RH, Mim, Mes, Mrin, FF, Mrout, Hrout, Mst )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, ks, ke
 integer, intent(in)    :: RunbackNum
 real   , intent(in)    :: dt
 real   , intent(in)    :: Ts(is:ie, ks:ke)
 real   , intent(in)    :: SA(is:ie, ks:ke), utau(is:ie, ks:ke), RH(is:ie, ks:ke)
 real   , intent(in)    :: Mim(is:ie, ks:ke), Mes(is:ie, ks:ke), Mrin(is:ie, ks:ke)
 real   , intent(in)    :: FF(is:ie, ks:ke)
 real   , intent(out)   :: Mrout(is:ie, ks:ke), Hrout(is:ie, ks:ke)
 real   , intent(inout) :: Mst(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: Mout0(:, :), Bw(:, :)
 integer :: i, k
 real    :: MFR
 real    :: Aout, Ast
 ! 処理開始 ********************************************************************************************
 allocate( Bw(is:ie, ks:ke) )
 Mrout(:,:) = 0.0; Bw(:,:) = 0.0
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  ! Runback-out 質量流束 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Mrout(i,k) = (1.0 - ff(i,k)) * (Mim(i,k) + Mrin(i,k) + Mst(i,k)) - Mes(i,k)
  if( Mrout(i,k) <= 0.0 ) then
    Mrout(i,k) = 0.0
    Mst(i,k)   = 0.0
   else
    select case( RunbackNum )
    ! 全てランバック -----------------------------------------------------------------------------------
     case(0)
      Mst(i,k) = 0.0
     case(1)
      ! 水膜の厚さ
      Bw(i,k) = ( Mrout(i,k) + Mst(i,k) ) / Rhow / SA(i,k)
      ! 面積充填率 (氷表面を半球と仮定)
      MFR  = 0.25 * 0.5 * 4.0 / 3.0 * pi
      Ast  = (1.0 - MFR) * 4.0 * RH(i,k)**3
      Aout = Mrout(i,k) / Rhow * Bw(i,k) / SA(i,k)
      if( Bw(i,k) <= 0.0 ) then
        Mst  (i,k) = Mrout(i,k)
        Mrout(i,k) = 0.0
       else
        if( Ast < Aout ) then
          Mst  (i,k) = Ast * Rhow / Bw(i,k) * SA(i,k)
          Mrout(i,k) = Mrout(i,k) - Mst(i,k)
         else
          Mst  (i,k) = Mrout(i,k)
          Mrout(i,k) = 0.0
        endif
       endif
     case default; write(*, '(a)') '!!!!! Error : Runback model number !!!!!'
    end select
  endif
  if( Mrout(i,k) > 0.0 ) then
    Hrout(i,k) = Cpw * (Ts(i,k) - Tf) + 0.5 * Utau(i,k)**2
   else
    Hrout(i,k) = 0.0
  endif
 enddo
 enddo
 deallocate(Bw)
 ! 処理終了 ********************************************************************************************
 return
end subroutine RunbackOutOrg3D
!*******************************************************************************************************
!******** Runback-out の質量流束 (二次元，Extended Messinger)					********
!*******************************************************************************************************
subroutine RunbackOutExt2D( &
&            is, ie, RunbackNum, dt, time, RH, Mim, Mes, Mrin, Bw, FF, nPhase, Mrout, Mst )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie
 integer, intent(in)    :: RunbackNum
 real   , intent(in)    :: dt, time
 real   , intent(in)    :: RH(is:ie)
 real   , intent(in)    :: Mim(is:ie), Mes(is:ie), Mrin(is:ie)
 real   , intent(in)    :: Bw(is:ie)
 real   , intent(in)    :: FF(is:ie)
 integer, intent(in)    :: nPhase(is:ie)
 real   , intent(out)   :: Mrout(is:ie)
 real   , intent(inout) :: Mst(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: Mout0(:)
 integer :: i
 real    :: MFR
 real    :: Aout, Ast
 real    :: Vwt, Vst
 real   , parameter :: A1 = 0.524
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
 ! Runback-out 質量流束 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 select case(nPhase(i))
  ! 着氷なし and 霧氷 (ランバックなし) -----------------------------------------------------------------
  case(0,1)
!  case(0)
   Mrout(i) = 0.0
   Mst(i)   = 0.0
  ! 雨氷 (ランバック発生) ------------------------------------------------------------------------------
  case(2)
!  case(1,2)
   Mrout(i) = ( 1.0 - FF(i) ) * ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) )
!   Mrout(i) = (1.0 - FF(i)) * (Mim(i) + Mrin(i) + Mst(i)) - Mes(i)
    if( Mrout(i) < 0.0 ) then
      Mrout(i) = 0.0
      Mst(i)   = 0.0
     else
     select case( RunbackNum )
      ! 全てランバック ---------------------------------------------------------------------------------
      case(1)
       Mst(i) = 0.0
       ! 粗さにひっかかって残る質量流束 ----------------------------------------------------------------
       case(0)
        ! 面積充填率 (氷表面を半球と仮定)
        MFR  = 0.5 * 0.5 * pi
        Ast  = (1.0 - MFR) * 2.0 * RH(i)**2
        Aout = Mrout(i) / Rhow * Bw(i)
!        Aout = Mrout(i) / Rhow * Bw(i) * (1.0 - FF(I)) * TIME
        if( Ast < Aout ) then
          Mst  (i) = Ast * Rhow / Bw(i)
!          Mst  (i) = Ast * Rhow / Bw(i) / (1.0 - FF(I)) / TIME
          Mrout(i) = Mrout(i) - Mst(i)
         else
          Mst  (i) = Mrout(i)
          Mrout(i) = 0.0
        endif
!        ! 磯部さんバージョン --------------------------------------------------------------------------
!        Vwt = Mrout(i) / rhow
!        Vst = (1.0 - A1) * RH(i)
!        if( Vst >= Vwt ) then
!         Mst(i) = Vwt * rhow
!         Mrout(i) = 0.0
!        else
!         Mst(i) = Vst * rhow
!         Mrout(i) = ( Vwt - Vst ) * rhow
!        endif
       case default; write(*, '(a)') '!!!!! Error : Runback model number !!!!!'
     end select
    endif
  ! 防氷 (ランバック発生) ------------------------------------------------------------------------------
  case(3)
   Mrout(i) = ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) )
 end select
 enddo
 ! リミッタ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! allocate( Mout0(is:ie) )
! Mout0(:,:) = Mrout(:,:)
! do i = is + 2, ie - 2
!  if( Mout0(i-2) == 0.0 .and. Mout0(i-1) == 0.0 .and. &
!  &   Mout0(i+1) == 0.0 .and. Mout0(i+2) == 0.0 ) then
!    Mst  (i) = Mst  (i) + Mrout(i)
!    Mrout(i) = 0.0
!  endif
! enddo
! deallocate( Mout0 )
 ! 処理終了 ********************************************************************************************
 return
end subroutine RunbackOutExt2D
!*******************************************************************************************************
!******** Runback-out の質量流束 (三次元，Extended Messinger)					********
!*******************************************************************************************************
subroutine RunbackOutExt3D( &
&            is, ie, ks, ke, RunbackNum, dt, time, RH, Mim, Mes, Mrin, Bw, FF, nPhase, Mrout, Mst )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, ks, ke
 integer, intent(in)    :: RunbackNum
 real   , intent(in)    :: dt, time
 real   , intent(in)    :: RH(is:ie, ks:ke)
 real   , intent(in)    :: Mim(is:ie, ks:ke), Mes(is:ie, ks:ke), Mrin(is:ie, ks:ke)
 real   , intent(in)    :: Bw(is:ie, ks:ke)
 real   , intent(in)    :: FF(is:ie, ks:ke)
 integer, intent(in)    :: nPhase(is:ie, ks:ke)
 real   , intent(out)   :: Mrout(is:ie, ks:ke)
 real   , intent(inout) :: Mst(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: Mout0(:, :)
 integer :: i, k
 real    :: MFR
 real    :: Aout, Ast
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  ! Runback-out 質量流束 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( nPhase(i,k) == 2 ) then
    Mrout(i,k) = (1.0 - FF(i,k)) * (Mim(i,k) + Mrin(i,k) + Mst(i,k)) - Mes(i,k)
    if( Mrout(i,k) < 0.0 ) then
      Mrout(i,k) = 0.0
      Mst(i,k)   = 0.0
     else
      select case( RunbackNum )
       ! 全てランバック -------------------------------------------------------------------------------
       case(0)
        Mst(i,k) = 0.0
       ! 粗さにひっかかって残る質量流束 ---------------------------------------------------------------
       case(1)
        ! 面積充填率 (氷表面を半球と仮定)
        MFR  = 0.25 * 0.5 * 4.0 / 3.0 * pi
        Ast  = (1.0 - MFR) * 4.0 * RH(i,k)**3
        Aout = Mrout(i,k) / Rhow * Bw(i,k)
!        Aout = Mrout(i,k) / Rhow * Bw(i,k) * (1.0 - FF(I,K)) * TIME
        if( Ast < Aout ) then
          Mst  (i,k) = Ast * Rhow / Bw(i,k)
!          Mst  (i,k) = Ast * Rhow / Bw(i,k) / (1.0 - FF(I,K)) / TIME
          Mrout(i,k) = Mrout(i,k) - Mst(i,k)
         else
          Mst  (i,k) = Mrout(i,k)
          Mrout(i,k) = 0.0
        endif
       case default; write(*, '(a)') '!!!!! Error : Runback model number !!!!!'
      end select
     endif
    else
     Mrout(i,k) = 0.0
     Mst  (i,k) = 0.0
   endif
 enddo
 enddo
 ! リミッタ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! allocate( Mout0(is:ie, ks:ke) )
! Mout0(:,:) = Mrout(:,:)
! do k = ks + 1, ke - 1
! do i = is + 1, ie - 1
!  if( Mout0(i-1,k  ) == 0.0 .and. Mout0(i-1,k  ) == 0.0 .and. &
!  &   Mout0(i  ,k-1) == 0.0 .and. Mout0(i  ,k-1) == 0.0 .and. &
!  &   Mout0(i+1,k  ) == 0.0 .and. Mout0(i+1,k  ) == 0.0 .and. &
!  &   Mout0(i  ,k+1) == 0.0 .and. Mout0(i  ,k+1) == 0.0 ) then
!!    Mst  (i,k) = Mst  (i,k) + Mrout(i,k)
!    Mrout(i,k) = 0.0
!  endif
! enddo
! enddo
! deallocate( Mout0 )
 ! 処理終了 ********************************************************************************************
 return
end subroutine RunbackOutExt3D

!*******************************************************************************************************
!******** 水膜逸脱 (二次元，Original Messinger, せん断力，遠心力のみ考慮)			********
!*******************************************************************************************************
subroutine WaterSheddingOrg2D( &
&            is, ie, dt, GravityX, GravityY, x, y, Rhof, up, vp, SA, utau, Mrout, Hrout )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie
 real   , intent(in)    :: dt, GravityX, GravityY
 real   , intent(in)    :: x(is:ie), y(is:ie)
 real   , intent(in)    :: Rhof(is:ie)
 real   , intent(in)    :: up(is:ie), vp(is:ie)
 real   , intent(in)    :: SA(is:ie), utau(is:ie)
 real   , intent(inout) :: Mrout(is:ie), Hrout(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 real    :: ax, ay, az, bx, by, bz, nx, ny, na
 real    :: Fx, Fy, FN
 real    :: SFx, SFy, GFx, GFy
 real    :: Velw, Uw, Vw
 real    :: Vd, Rd, ST
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
  if( Mrout(i) <= 0.0 ) cycle
  ! 壁面せん断力 ---------------------------------------------------------------------------------------
  nx = up(i)
  ny = vp(i)
  na = sqrt(nx**2 + ny**2)
  if(na /= 0) then
    nx = nx / na
    ny = ny / na
    SFx = Utau(i)**2 * Rhof(i) * SA(i) * nx
    SFy = Utau(i)**2 * Rhof(i) * SA(i) * ny
   else
    SFx = 0.0
    SFy = 0.0
  endif
  ! 重力 -----------------------------------------------------------------------------------------------
  GFx = Mrout(i) * dt * GravityX
  GFy = Mrout(i) * dt * GravityY
  ! 壁面法線方向に働く力 -------------------------------------------------------------------------------
  ax =  0.5 * ( - x(i-1) + x(i+1) )
  ay =  0.5 * ( - y(i-1) + y(i+1) )
  az =  0.0
  bx =  0.0
  by =  0.0
  bz = -1.0
  nx = ay * bz - az * by
  ny = az * bx - ax * bz
  na = sqrt(nx**2 + ny**2)
  nx = +1.0 * nx / na
  ny = +1.0 * ny / na
  Fx = SFx + GFx; Fy = SFy + GFy
  FN = Fx * nx + Fy * ny
  ! 表面張力 -------------------------------------------------------------------------------------------
  ! 液滴体積
  Vd = Mrout(i) * dt / Rhow
  ! 液滴半径
  Rd = (1.5 / pi * Vd)**(1.0 / 3.0)
  ! 表面張力
  ST = sigw * 2.0 * pi * Rd
  ! 逸脱判定 -------------------------------------------------------------------------------------------
  if( FN > ST ) then
    Mrout(i) = 0.0
    Hrout(i) = 0.0
  endif
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine WaterSheddingOrg2D
!*******************************************************************************************************
!******** 水膜逸脱 (三次元，Original Messinger, せん断力，遠心力のみ考慮)			********
!*******************************************************************************************************
subroutine WaterSheddingOrg3D( &
&            is, ie, ks, ke, dt, GravityX, GravityY, GravityZ, OmegaX, OmegaY, OmegaZ, &
&            x, y, z, Rhof, up, vp, wp, SA, utau, Mrout, Hrout )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, ks, ke
 real   , intent(in)    :: dt, GravityX, GravityY, GravityZ, OmegaX, OmegaY, OmegaZ
 real   , intent(in)    :: x(is:ie, ks:ke), y(is:ie, ks:ke), z(is:ie, ks:ke)
 real   , intent(in)    :: Rhof(is:ie, ks:ke)
 real   , intent(in)    :: up(is:ie, ks:ke), vp(is:ie, ks:ke), wp(is:ie, ks:ke)
 real   , intent(in)    :: SA(is:ie, ks:ke), utau(is:ie, ks:ke)
 real   , intent(inout) :: Mrout(is:ie, ks:ke), Hrout(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 real    :: ax, ay, az, bx, by, bz, nx, ny, nz, na
 real    :: Fx, Fy, Fz, FN
 real    :: SFx, SFy, SFz, GFx, GFy, GFz, CFx, CFy, CFz
 real    :: Velw, Uw, Vw, Ww
 real    :: Vd, Rd, ST
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  if( Mrout(i,k) <= 0.0 ) cycle
  ! 壁面せん断力 ---------------------------------------------------------------------------------------
  nx = up(i,k)
  ny = vp(i,k)
  nz = wp(i,k)
  na = sqrt(nx**2 + ny**2 + nz**2)
  if(na /= 0) then
    nx = nx / na
    ny = ny / na
    nz = nz / na
    SFx = Utau(i,k)**2 * Rhof(i,k) * SA(i,k) * nx
    SFy = Utau(i,k)**2 * Rhof(i,k) * SA(i,k) * ny
    SFz = Utau(i,k)**2 * Rhof(i,k) * SA(i,k) * nz
   else
    SFx = 0.0
    SFy = 0.0
    SFz = 0.0
  endif
  ! 遠心力 ---------------------------------------------------------------------------------------------
  CFx = Omegay * (Omegax * y(i,k) - Omegay * x(i,k)) - Omegaz * (Omegaz * x(i,k) - Omegax * z(i,k))
  CFy = Omegaz * (Omegay * z(i,k) - Omegaz * y(i,k)) - Omegax * (Omegax * y(i,k) - Omegay * x(i,k))
  CFz = Omegax * (Omegaz * x(i,k) - Omegax * z(i,k)) - Omegay * (Omegay * z(i,k) - Omegaz * y(i,k))
  CFx = -1.0 * CFx * Mrout(i,k) * dt
  CFy = -1.0 * CFy * Mrout(i,k) * dt
  CFz = -1.0 * CFz * Mrout(i,k) * dt
  ! 重力 -----------------------------------------------------------------------------------------------
  GFx = Mrout(i,k) * dt * GravityX
  GFy = Mrout(i,k) * dt * GravityY
  GFz = Mrout(i,k) * dt * GravityZ
  ! 壁面法線方向に働く力 -------------------------------------------------------------------------------
  ax = 0.5 * ( - x(i-1,k  ) + x(i+1,k  ) )
  ay = 0.5 * ( - y(i-1,k  ) + y(i+1,k  ) )
  az = 0.5 * ( - z(i-1,k  ) + z(i+1,k  ) )
  bx = 0.5 * ( + x(i  ,k-1) - x(i  ,k+1) )
  by = 0.5 * ( + y(i  ,k-1) - y(i  ,k+1) )
  bz = 0.5 * ( + z(i  ,k-1) - z(i  ,k+1) )
  nx = ay * bz - az * by
  ny = az * bx - ax * bz
  nz = ax * by - ay * bx
  na = sqrt(nx**2 + ny**2 + nz**2)
  nx = +1.0 * nx / na
  ny = +1.0 * ny / na
  nz = +1.0 * nz / na
  Fx = SFx + CFx + GFx; Fy = SFy + CFy + GFy; Fz = SFz + CFz + GFz
  FN = Fx * nx + Fy * ny + Fz * nz
  ! 表面張力 -------------------------------------------------------------------------------------------
  ! 液滴体積
  Vd = Mrout(i,k) * dt / Rhow
  ! 液滴半径
  Rd = (1.5 / pi * Vd)**(1.0 / 3.0)
  ! 表面張力
  ST = sigw * 2.0 * pi * Rd
  ! 逸脱判定 -------------------------------------------------------------------------------------------
  if( FN > ST ) then
    Mrout(i,k) = 0.0
    Hrout(i,k) = 0.0
  endif
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine WaterSheddingOrg3D
!*******************************************************************************************************
!******** 水膜逸脱 (二次元，Extended Messinger, せん断力，遠心力のみ考慮)			********
!*******************************************************************************************************
subroutine WaterSheddingExt2D( &
&            is, ie, dt, GravityX, GravityY, x, y, Rhof, up, vp, SA, utau, Mrout )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie
 real   , intent(in)    :: dt, GravityX, GravityY
 real   , intent(in)    :: x(is:ie), y(is:ie)
 real   , intent(in)    :: Rhof(is:ie)
 real   , intent(in)    :: up(is:ie), vp(is:ie)
 real   , intent(in)    :: SA(is:ie), utau(is:ie)
 real   , intent(inout) :: Mrout(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 real    :: ax, ay, az, bx, by, bz, nx, ny, na
 real    :: Fx, Fy, FN
 real    :: SFx, SFy, GFx, GFy
 real    :: Velw, Uw, Vw
 real    :: Vd, Rd, ST
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
  if( Mrout(i) == 0.0 ) cycle
  ! 壁面せん断力 ---------------------------------------------------------------------------------------
  nx = up(i)
  ny = vp(i)
  na = sqrt(nx**2 + ny**2)
  if(na /= 0) then
    nx = nx / na
    ny = ny / na
    SFx = Utau(i)**2 * Rhof(i) * SA(i) * nx
    SFy = Utau(i)**2 * Rhof(i) * SA(i) * ny
   else
    SFx = 0.0
    SFy = 0.0
  endif
  ! 重力 -----------------------------------------------------------------------------------------------
  GFx = Mrout(i) * SA(i) * dt * GravityX
  GFy = Mrout(i) * SA(i) * dt * GravityY
  ! 壁面法線方向に働く力 -------------------------------------------------------------------------------
  ax =  0.5 * ( - x(i-1) + x(i+1) )
  ay =  0.5 * ( - y(i-1) + y(i+1) )
  az =  0.0
  bx =  0.0
  by =  0.0
  bz = -1.0
  nx = ay * bz - az * by
  ny = az * bx - ax * bz
  na = sqrt(nx**2 + ny**2)
  nx = +1.0 * nx / na
  ny = +1.0 * ny / na
  Fx = SFx + GFx; Fy = SFy + GFy
  FN = Fx * nx + Fy * ny
  ! 表面張力 -------------------------------------------------------------------------------------------
  ! 液滴体積
  Vd = Mrout(i) * SA(i) * dt / Rhow
  ! 液滴半径
  Rd = (1.5 / pi * Vd)**(1.0 / 3.0)
  ! 表面張力
  ST = sigw * 2.0 * pi * Rd
  ! 逸脱判定 -------------------------------------------------------------------------------------------
  if( FN > ST ) Mrout(i) = 0.0
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine WaterSheddingExt2D
!*******************************************************************************************************
!******** 水膜逸脱 (三次元，Extended Messinger, せん断力，遠心力のみ考慮)			********
!*******************************************************************************************************
subroutine WaterSheddingExt3D( &
&            is, ie, ks, ke, dt, GravityX, GravityY, GravityZ, OmegaX, OmegaY, OmegaZ, &
&            x, y, z, Rhof, up, vp, wp, SA, utau, Mrout )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, ks, ke
 real   , intent(in)    :: dt, GravityX, GravityY, GravityZ, OmegaX, OmegaY, OmegaZ
 real   , intent(in)    :: x(is:ie, ks:ke), y(is:ie, ks:ke), z(is:ie, ks:ke)
 real   , intent(in)    :: Rhof(is:ie, ks:ke)
 real   , intent(in)    :: up(is:ie, ks:ke), vp(is:ie, ks:ke), wp(is:ie, ks:ke)
 real   , intent(in)    :: SA(is:ie, ks:ke), utau(is:ie, ks:ke)
 real   , intent(inout) :: Mrout(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 real    :: ax, ay, az, bx, by, bz, nx, ny, nz, na
 real    :: Fx, Fy, Fz, FN
 real    :: SFx, SFy, SFz, CFx, CFy, CFz, GFx, GFy, GFz
 real    :: Velw, Uw, Vw, Ww
 real    :: Vd, Rd, ST
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  if( Mrout(i,k) == 0.0 ) cycle
  ! 壁面せん断力 ---------------------------------------------------------------------------------------
  nx = up(i,k)
  ny = vp(i,k)
  nz = wp(i,k)
  na = sqrt(nx**2 + ny**2 + nz**2)
  if(na /= 0) then
    nx = nx / na
    ny = ny / na
    nz = nz / na
    SFx = Utau(i,k)**2 * Rhof(i,k) * SA(i,k) * nx
    SFy = Utau(i,k)**2 * Rhof(i,k) * SA(i,k) * ny
    SFz = Utau(i,k)**2 * Rhof(i,k) * SA(i,k) * nz
   else
    SFx = 0.0
    SFy = 0.0
    SFz = 0.0
  endif
  ! 遠心力 ---------------------------------------------------------------------------------------------
  CFx = Omegay * (Omegax * y(i,k) - Omegay * x(i,k)) - Omegaz * (Omegaz * x(i,k) - Omegax * z(i,k))
  CFy = Omegaz * (Omegay * z(i,k) - Omegaz * y(i,k)) - Omegax * (Omegax * y(i,k) - Omegay * x(i,k))
  CFz = Omegax * (Omegaz * x(i,k) - Omegax * z(i,k)) - Omegay * (Omegay * z(i,k) - Omegaz * y(i,k))
  CFx = -1.0 * CFx * Mrout(i,k) * SA(i,k) * dt
  CFy = -1.0 * CFy * Mrout(i,k) * SA(i,k) * dt
  CFz = -1.0 * CFz * Mrout(i,k) * SA(i,k) * dt
  ! 重力 -----------------------------------------------------------------------------------------------
  GFx = Mrout(i,k) * SA(i,k) * dt * GravityX
  GFy = Mrout(i,k) * SA(i,k) * dt * GravityY
  GFz = Mrout(i,k) * SA(i,k) * dt * GravityZ
  ! 壁面法線方向に働く力 -------------------------------------------------------------------------------
  ax = 0.5 * ( - x(i-1,k  ) + x(i+1,k  ) )
  ay = 0.5 * ( - y(i-1,k  ) + y(i+1,k  ) )
  az = 0.5 * ( - z(i-1,k  ) + z(i+1,k  ) )
  bx = 0.5 * ( + x(i  ,k-1) - x(i  ,k+1) )
  by = 0.5 * ( + y(i  ,k-1) - y(i  ,k+1) )
  bz = 0.5 * ( + z(i  ,k-1) - z(i  ,k+1) )
  nx = ay * bz - az * by
  ny = az * bx - ax * bz
  nz = ax * by - ay * bx
  na = sqrt(nx**2 + ny**2 + nz**2)
  nx = +1.0 * nx / na
  ny = +1.0 * ny / na
  nz = +1.0 * nz / na
  Fx = SFx + CFx + GFx; Fy = SFy + CFy + GFy; Fz = SFz + CFz + GFz
  FN = Fx * nx + Fy * ny + Fz * nz
  ! 表面張力 -------------------------------------------------------------------------------------------
  ! 液滴体積
  Vd = Mrout(i,k) * SA(i,k) * dt / Rhow
  ! 液滴半径
  Rd = (1.5 / pi * Vd)**(1.0 / 3.0)
  ! 表面張力
  ST = sigw * 2.0 * pi * Rd
  ! 逸脱判定 -------------------------------------------------------------------------------------------
  if( FN > ST ) Mrout(i,k) = 0.0
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine WaterSheddingExt3D
!*******************************************************************************************************
!******** Runback-in の質量・エンタルピ (二次元，Original Messinger)				********
!*******************************************************************************************************
subroutine RunbackInOrg2D( &
&            is, ie, dt, GravityX, GravityY, x, y, Rhof, up, vp, SA, utau, Mrout, Hrout, Mrin, Hrin )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: dt, GravityX, GravityY
 real   , intent(in)  :: x(is:ie), y(is:ie)
 real   , intent(in)  :: Rhof(is:ie)
 real   , intent(in)  :: up(is:ie), vp(is:ie)
 real   , intent(in)  :: SA(is:ie), utau(is:ie)
 real   , intent(in)  :: Mrout(is:ie), Hrout(is:ie)
 real   , intent(out) :: Mrin(is:ie), Hrin(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, pointer :: ip(:)
 real   , pointer :: px(:), py(:)
 real   , pointer :: pro(:)
 integer :: i, m
 integer :: mp
 real    :: nx, ny, na
 real    :: Fx, Fy, FA
 real    :: SFx, SFy, GFx, GFy
 real    :: Velw, Uw, Vw
 real    :: pa, ProMax
 ! 処理開始 ********************************************************************************************
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( ip(1:2), px(1:2), py(1:2), pro(1:2) )
 ! 初期化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Mrin(:) = 0.0; Hrin(:) = 0.0
 ! 周囲 2点の内最も力が作用する 1点のみに Runback-in すると仮定 ++++++++++++++++++++++++++++++++++++++++
 do i = is + 1, ie - 1
  if( Mrout(i) <= 0.0 ) cycle
  ! 壁面せん断力 ---------------------------------------------------------------------------------------
  nx = up(i)
  ny = vp(i)
  na = sqrt(nx**2 + ny**2)
  if(na /= 0) then
    nx = nx / na
    ny = ny / na
    SFx = Utau(i)**2 * Rhof(i) * SA(i) * nx
    SFy = Utau(i)**2 * Rhof(i) * SA(i) * ny
   else
    SFx = 0.0
    SFy = 0.0
  endif
  ! 重力 -----------------------------------------------------------------------------------------------
  GFx = Mrout(i) * dt * GravityX
  GFy = Mrout(i) * dt * GravityY
  ! 壁面法線方向に働く力 -------------------------------------------------------------------------------
  Fx = SFx + GFx; Fy = SFy + GFy
  ! 周囲 2点の単位位置ベクトル -------------------------------------------------------------------------
  ip(1) = i-1; ip(2) = i+1
  do m = 1, 2
   px(m) = x(ip(m)) - x(i)
   py(m) = y(ip(m)) - y(i)
   pa    = sqrt( px(m)**2 + py(m)**2 )
   px(m) = px(m) / pa
   py(m) = py(m) / pa
  enddo
  ! 液滴に働く力のベクトルと位置ベクトルの内積 ---------------------------------------------------------
  Fa = sqrt(Fx**2 + Fy**2)
  if( Fa /= 0.0 ) then
    Fx = Fx / Fa
    Fy = Fy / Fa
   else
    cycle
  endif
  do m = 1, 2
   pro(m) = Fx * px(m) + Fy * py(m)
  enddo
  ! 液滴に働く力が傾いている方向にRunback (ベクトルの成す角が小さい方向) -------------------------------
  mp = maxval( maxloc(pro(:)) )
  Mrin(ip(mp)) = Mrin(ip(mp)) + Mrout(i)
  Hrin(ip(mp)) = Hrin(ip(mp)) + Hrout(i)
 enddo
 ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 deallocate( ip, px, py, pro)
 ! 処理終了 ********************************************************************************************
 return
end subroutine RunbackInOrg2D
!*******************************************************************************************************
!******** Runback-in の質量・エンタルピ (三次元，Original Messinger)				********
!*******************************************************************************************************
subroutine RunbackInOrg3D( &
&            is, ie, ks, ke, dt, GravityX, GravityY, GravityZ, OmegaX, OmegaY, OmegaZ, &
&            x, y, z, Rhof, up, vp, wp, SA, utau, Mrout, Hrout, Mrin, Hrin )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: dt, GravityX, GravityY, GravityZ, OmegaX, OmegaY, OmegaZ
 real   , intent(in)  :: x(is:ie, ks:ke), y(is:ie, ks:ke), z(is:ie, ks:ke)
 real   , intent(in)  :: Rhof(is:ie, ks:ke)
 real   , intent(in)  :: up(is:ie, ks:ke), vp(is:ie, ks:ke), wp(is:ie, ks:ke)
 real   , intent(in)  :: SA(is:ie, ks:ke), utau(is:ie, ks:ke)
 real   , intent(in)  :: Mrout(is:ie, ks:ke), Hrout(is:ie, ks:ke)
 real   , intent(out) :: Mrin(is:ie, ks:ke), Hrin(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, pointer :: ip(:), kp(:), mp(:)
 real   , pointer :: px(:), py(:), pz(:)
 real   , pointer :: pro(:), wr(:)
 integer :: i, k, m
 real    :: ax, ay, az, bx, by, bz, nx, ny, nz, na
 real    :: Fx, Fy, Fz, FA
 real    :: SFx, SFy, SFz, CFx, CFy, CFz, GFx, GFy, GFz
 real    :: Velw, Uw, Vw, Ww
 real    :: pa, ProMax
 ! 処理開始 ********************************************************************************************
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( ip(1:8), kp(1:8), px(1:8), py(1:8), pz(1:8), pro(1:8), mp(1:3), wr(1:3) )
 ! 初期化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Mrin(:, :) = 0.0; Hrin(:, :) = 0.0
 ! 周囲 8点の内で力が作用する 3点のみに Runback-in すると仮定 ++++++++++++++++++++++++++++++++++++++++++
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  if( Mrout(i,k) <= 0.0 ) cycle
  ! 壁面せん断力 ---------------------------------------------------------------------------------------
  nx = up(i,k)
  ny = vp(i,k)
  nz = wp(i,k)
  na = sqrt(nx**2 + ny**2 + nz**2)
  if(na /= 0) then
    nx = nx / na
    ny = ny / na
    nz = nz / na
    SFx = Utau(i,k)**2 * Rhof(i,k) * SA(i,k) * nx
    SFy = Utau(i,k)**2 * Rhof(i,k) * SA(i,k) * ny
    SFz = Utau(i,k)**2 * Rhof(i,k) * SA(i,k) * nz
   else
    SFx = 0.0
    SFy = 0.0
    SFz = 0.0
  endif
  ! 遠心力 ---------------------------------------------------------------------------------------------
  CFx = Omegay * (Omegax * y(i,k) - Omegay * x(i,k)) - Omegaz * (Omegaz * x(i,k) - Omegax * z(i,k))
  CFy = Omegaz * (Omegay * z(i,k) - Omegaz * y(i,k)) - Omegax * (Omegax * y(i,k) - Omegay * x(i,k))
  CFz = Omegax * (Omegaz * x(i,k) - Omegax * z(i,k)) - Omegay * (Omegay * z(i,k) - Omegaz * y(i,k))
  CFx = -1.0 * CFx * Mrout(i,k) * dt
  CFy = -1.0 * CFy * Mrout(i,k) * dt
  CFz = -1.0 * CFz * Mrout(i,k) * dt
  ! 重力 -----------------------------------------------------------------------------------------------
  GFx = Mrout(i,k) * dt * GravityX
  GFy = Mrout(i,k) * dt * GravityY
  GFz = Mrout(i,k) * dt * GravityZ
  ! 壁面法線方向に働く力 -------------------------------------------------------------------------------
  Fx = SFx + CFx + GFx; Fy = SFy + CFy + GFy; Fz = SFz + CFz + GFz
  ! 周囲 8点の単位位置ベクトル -------------------------------------------------------------------------
  ip(1) = i+1; kp(1) = k+1; ip(2) = i  ; kp(2) = k+1; ip(3) = i-1; kp(3) = k+1
  ip(4) = i+1; kp(4) = k  ;                           ip(5) = i-1; kp(5) = k
  ip(6) = i+1; kp(6) = k-1; ip(7) = i  ; kp(7) = k-1; ip(8) = i-1; kp(8) = k-1
  do m = 1, 8
   px(m) = x(ip(m), kp(m)) - x(i,k)
   py(m) = y(ip(m), kp(m)) - y(i,k)
   pz(m) = z(ip(m), kp(m)) - z(i,k)
   pa    = sqrt( px(m)**2 + py(m)**2 + pz(m)**2 )
   px(m) = px(m) / pa
   py(m) = py(m) / pa
   pz(m) = pz(m) / pa
  enddo
  ! 液滴に働く力のベクトルと位置ベクトルの内積 ---------------------------------------------------------
  Fa = sqrt(Fx**2 + Fy**2 + Fz**2)
  if( Fa /= 0.0 ) then
    Fx = Fx / Fa
    Fy = Fy / Fa
    Fz = Fz / Fa
   else
    cycle
  endif
  do m = 1, 8
   pro(m) = Fx * px(m) + Fy * py(m) + Fz * pz(m)
  enddo
  ! 液滴に働く力が傾いている方向にRunback (ベクトルの成す角が小さい方向) -------------------------------
  mp(1) = maxval( maxloc(pro(:)) )
  mp(2) = maxval( maxloc(pro(:), mask = pro(:) < pro(mp(1))) )
  mp(3) = maxval( maxloc(pro(:), mask = pro(:) < pro(mp(2))) )
  if( pro(mp(1)) >= 0.0 .and. pro(mp(2)) >= 0.0 .and. pro(mp(3)) >= 0.0 ) then
    do m = 1, 3
     wr(m) = pro(mp(m)) / ( pro(mp(1)) + pro(mp(2)) + pro(mp(3)) )
     Mrin(ip(mp(m)), kp(mp(m))) = Mrin(ip(mp(m)), kp(mp(m))) + Mrout(i,k) * wr(m)
     Hrin(ip(mp(m)), kp(mp(m))) = Hrin(ip(mp(m)), kp(mp(m))) + Hrout(i,k) * wr(m)
    enddo
   else if( pro(mp(1)) >= 0.0 .and. pro(mp(2)) >= 0.0 ) then
    do m = 1, 2
     wr(m) = pro(mp(m)) / ( pro(mp(1)) + pro(mp(2)) )
     Mrin(ip(mp(m)), kp(mp(m))) = Mrin(ip(mp(m)), kp(mp(m))) + Mrout(i,k) * wr(m)
     Hrin(ip(mp(m)), kp(mp(m))) = Hrin(ip(mp(m)), kp(mp(m))) + Hrout(i,k) * wr(m)
    enddo
   else if( pro(mp(1)) >= 0.0 ) then
    m = 1
    wr(m) = pro(mp(m)) / pro(mp(1))
    Mrin(ip(mp(m)), kp(mp(m))) = Mrin(ip(mp(m)), kp(mp(m))) + Mrout(i,k) * wr(m)
    Hrin(ip(mp(m)), kp(mp(m))) = Hrin(ip(mp(m)), kp(mp(m))) + Hrout(i,k) * wr(m)
   else
    cycle
  endif
 enddo
 enddo
 ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 deallocate( ip, kp, px, py, pz, pro, mp, wr )
 ! 処理終了 ********************************************************************************************
 return
end subroutine RunbackInOrg3D
!*******************************************************************************************************
!******** Runback-in の質量流束 (二次元，Extended Messinger)					********
!*******************************************************************************************************
subroutine RunbackInExt2D( &
&            is, ie, dt, GravityX, GravityY, x, y, Rhof, up, vp, SA, utau, Mrout, Mrin )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: dt, GravityX, GravityY
 real   , intent(in)  :: x(is:ie), y(is:ie)
 real   , intent(in)  :: Rhof(is:ie)
 real   , intent(in)  :: up(is:ie), vp(is:ie)
 real   , intent(in)  :: SA(is:ie), utau(is:ie)
 real   , intent(in)  :: Mrout(is:ie)
 real   , intent(out) :: Mrin(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, pointer :: ip(:), kp(:)
 real   , pointer :: px(:), py(:)
 real   , pointer :: pro(:)
 integer :: i, m
 integer :: mp
 real    :: nx, ny, na
 real    :: Fx, Fy, Fz, FA
 real    :: SFx, SFy, SFz, CFx, CFy, CFz, GFx, GFy, GFz
 real    :: Velw, Uw, Vw, Ww
 real    :: pa, ProMax
 ! 処理開始 ********************************************************************************************
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( ip(1:2), px(1:2), py(1:2), pro(1:2) )
 ! 初期化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Mrin(:) = 0.0
 ! 周囲 2点の内最も力が作用する 1点のみに Runback-in すると仮定 ++++++++++++++++++++++++++++++++++++++++
 do i = is + 1, ie - 1
  if( Mrout(i) == 0.0 ) cycle
  ! 壁面せん断力 ---------------------------------------------------------------------------------------
  nx = up(i)
  ny = vp(i)
  na = sqrt(nx**2 + ny**2)
  if(na /= 0) then
    nx = nx / na
    ny = ny / na
    SFx = Utau(i)**2 * Rhof(i) * SA(i) * nx
    SFy = Utau(i)**2 * Rhof(i) * SA(i) * ny
   else
    SFx = 0.0
    SFy = 0.0
  endif
  ! 重力 -----------------------------------------------------------------------------------------------
  GFx = Mrout(i) * SA(i) * dt * GravityX
  GFy = Mrout(i) * SA(i) * dt * GravityY
  ! 壁面法線方向に働く力 -------------------------------------------------------------------------------
  Fx = SFx + GFx; Fy = SFy + GFy
  ! 周囲 2点の単位位置ベクトル -------------------------------------------------------------------------
  ip(1) = i-1; ip(2) = i+1
  do m = 1, 2
   px(m) = x(ip(m)) - x(i)
   py(m) = y(ip(m)) - y(i)
   pa    = sqrt( px(m)**2 + py(m)**2 )
   if( pa /= 0.0 ) then
     px(m) = px(m) / pa
     py(m) = py(m) / pa
    else
     px(m) = 0.0
     py(m) = 0.0
   endif
  enddo
  ! 液滴に働く力のベクトルと位置ベクトルの内積 ---------------------------------------------------------
  Fa = sqrt(Fx**2 + Fy**2)
  if( Fa /= 0.0 ) then
    Fx = Fx / Fa
    Fy = Fy / Fa
   else
    cycle
  endif
  do m = 1, 2
   pro(m) = Fx * px(m) + Fy * py(m)
  enddo
  ! 液滴に働く力が傾いている方向にRunback (ベクトルの成す角が小さい方向) -------------------------------
  mp = maxval( maxloc(pro(:)) )
  Mrin(ip(mp)) = Mrin(ip(mp)) + Mrout(i)
 enddo
 ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 deallocate( ip, px, py, pro )
 ! 処理終了 ********************************************************************************************
 return
end subroutine RunbackInExt2D
!!*******************************************************************************************************
!!******** Runback-in の質量流束 (二次元，Extended Messinger)					********
!!*******************************************************************************************************
!subroutine RunbackInExt2D( &
!&            is, ie, dt, GravityX, GravityY, x, y, Rhof, up, vp, SA, utau, Mrout, Mrin )
! ! 変数宣言 ********************************************************************************************
! implicit none
! ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! integer, intent(in)  :: is, ie
! real   , intent(in)  :: dt, GravityX, GravityY
! real   , intent(in)  :: x(is:ie), y(is:ie)
! real   , intent(in)  :: Rhof(is:ie)
! real   , intent(in)  :: up(is:ie), vp(is:ie)
! real   , intent(in)  :: SA(is:ie), utau(is:ie)
! real   , intent(in)  :: Mrout(is:ie)
! real   , intent(out) :: Mrin(is:ie)
! ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! integer, pointer :: ip(:), kp(:)
! real   , pointer :: px(:), py(:)
! real   , pointer :: pro(:)
! integer :: i, m
! integer :: mp
! real    :: ax, ay, az
! real    :: bx, by, bz
! real    :: nx, ny, na
! real    :: Fx, Fy, Fz, FA
! real    :: SFx, SFy, SFz, CFx, CFy, CFz, GFx, GFy, GFz
! real    :: Velw, Uw, Vw, Ww
! real    :: pa, ProMax
! ! 処理開始 ********************************************************************************************
! ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! allocate( ip(1:2), px(1:2), py(1:2), pro(1:2) )
! ! 初期化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Mrin(:) = 0.0
! ! 周囲 2点の内最も力が作用する 1点のみに Runback-in すると仮定 ++++++++++++++++++++++++++++++++++++++++
! do i = is + 1, ie - 1
!  if( Mrout(i) == 0.0 ) cycle
!  ! 壁面法線ベクトル -----------------------------------------------------------------------------------
!  ax =  0.5 * ( - x(i-1) + x(i+1) )
!  ay =  0.5 * ( - y(i-1) + y(i+1) )
!  az =  0.0
!  bx =  0.0
!  by =  0.0
!  bz = -1.0
!  nx = ay * bz - az * by
!  ny = az * bx - ax * bz
!  na = sqrt(nx**2 + ny**2)
!  nx = +1.0 * nx / na
!  ny = +1.0 * ny / na
!  ! 壁面せん断力 ---------------------------------------------------------------------------------------
!  nx = up(i)
!  ny = vp(i)
!  na = sqrt(nx**2 + ny**2)
!  if(na /= 0) then
!    nx = nx / na
!    ny = ny / na
!    SFx = Utau(i)**2 * Rhof(i) * SA(i) * nx
!    SFy = Utau(i)**2 * Rhof(i) * SA(i) * ny
!   else
!    SFx = 0.0
!    SFy = 0.0
!  endif
!  ! 重力 -----------------------------------------------------------------------------------------------
!  GFx = Mrout(i) * SA(i) * dt * GravityX
!  GFy = Mrout(i) * SA(i) * dt * GravityY
!  ! 壁面法線方向に働く力 -------------------------------------------------------------------------------
!  Fx = SFx + GFx; Fy = SFy + GFy
!  ! 周囲 2点の単位位置ベクトル -------------------------------------------------------------------------
!  ip(1) = i-1; ip(2) = i+1
!  do m = 1, 2
!   px(m) = x(ip(m)) - x(i)
!   py(m) = y(ip(m)) - y(i)
!   pa    = sqrt( px(m)**2 + py(m)**2 )
!   if( pa /= 0.0 ) then
!     px(m) = px(m) / pa
!     py(m) = py(m) / pa
!    else
!     px(m) = 0.0
!     py(m) = 0.0
!   endif
!  enddo
!  ! 液滴に働く力のベクトルと位置ベクトルの内積 ---------------------------------------------------------
!  Fa = sqrt(Fx**2 + Fy**2)
!  if( Fa /= 0.0 ) then
!    Fx = Fx / Fa
!    Fy = Fy / Fa
!   else
!    cycle
!  endif
!  do m = 1, 2
!   pro(m) = Fx * px(m) + Fy * py(m)
!  enddo
!  ! 液滴に働く力が傾いている方向にRunback (ベクトルの成す角が小さい方向) -------------------------------
!  mp = maxval( maxloc(pro(:)) )
!  Mrin(ip(mp)) = Mrin(ip(mp)) + Mrout(i)
! enddo
! ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! deallocate( ip, px, py, pro )
! ! 処理終了 ********************************************************************************************
! return
!end subroutine RunbackInExt2D
!*******************************************************************************************************
!******** Runback-in の質量流束 (三次元，Extended Messinger)					********
!*******************************************************************************************************
subroutine RunbackInExt3D( &
&            is, ie, ks, ke, dt, GravityX, GravityY, GravityZ, OmegaX, OmegaY, OmegaZ, &
&            x, y, z, Rhof, up, vp, wp, SA, utau, Mrout, Mrin )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: dt, GravityX, GravityY, GravityZ, OmegaX, OmegaY, OmegaZ
 real   , intent(in)  :: x(is:ie, ks:ke), y(is:ie, ks:ke), z(is:ie, ks:ke)
 real   , intent(in)  :: Rhof(is:ie, ks:ke)
 real   , intent(in)  :: up(is:ie, ks:ke), vp(is:ie, ks:ke), wp(is:ie, ks:ke)
 real   , intent(in)  :: SA(is:ie, ks:ke), utau(is:ie, ks:ke)
 real   , intent(in)  :: Mrout(is:ie, ks:ke)
 real   , intent(out) :: Mrin(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, pointer :: ip(:), kp(:), mp(:)
 real   , pointer :: px(:), py(:), pz(:)
 real   , pointer :: pro(:), wr(:)
 integer :: i, k, m
 real    :: ax, ay, az, bx, by, bz, nx, ny, nz, na
 real    :: Fx, Fy, Fz, FA
 real    :: SFx, SFy, SFz, CFx, CFy, CFz, GFx, GFy, GFz
 real    :: Velw, Uw, Vw, Ww
 real    :: pa, ProMax
 ! 処理開始 ********************************************************************************************
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( ip(1:8), kp(1:8), px(1:8), py(1:8), pz(1:8), pro(1:8), mp(1:3), wr(1:3) )
 ! 初期化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Mrin(:, :) = 0.0
 ! 周囲 8点の内で力が作用する 3点のみに Runback-in すると仮定 ++++++++++++++++++++++++++++++++++++++++++
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  if( Mrout(i,k) == 0.0 ) cycle
  ! 壁面せん断力 ---------------------------------------------------------------------------------------
  nx = up(i,k)
  ny = vp(i,k)
  nz = wp(i,k)
  na = sqrt(nx**2 + ny**2 + nz**2)
  if(na /= 0) then
    nx = nx / na
    ny = ny / na
    nz = nz / na
    SFx = Utau(i,k)**2 * Rhof(i,k) * SA(i,k) * nx
    SFy = Utau(i,k)**2 * Rhof(i,k) * SA(i,k) * ny
    SFz = Utau(i,k)**2 * Rhof(i,k) * SA(i,k) * nz
   else
    SFx = 0.0
    SFy = 0.0
    SFz = 0.0
  endif
  ! 遠心力 ---------------------------------------------------------------------------------------------
  CFx = Omegay * (Omegax * y(i,k) - Omegay * x(i,k)) - Omegaz * (Omegaz * x(i,k) - Omegax * z(i,k))
  CFy = Omegaz * (Omegay * z(i,k) - Omegaz * y(i,k)) - Omegax * (Omegax * y(i,k) - Omegay * x(i,k))
  CFz = Omegax * (Omegaz * x(i,k) - Omegax * z(i,k)) - Omegay * (Omegay * z(i,k) - Omegaz * y(i,k))
  CFx = -1.0 * CFx * Mrout(i,k) * SA(i,k) * dt
  CFy = -1.0 * CFy * Mrout(i,k) * SA(i,k) * dt
  CFz = -1.0 * CFz * Mrout(i,k) * SA(i,k) * dt
  ! 重力 -----------------------------------------------------------------------------------------------
  GFx = Mrout(i,k) * SA(i,k) * dt * GravityX
  GFy = Mrout(i,k) * SA(i,k) * dt * GravityY
  GFz = Mrout(i,k) * SA(i,k) * dt * GravityZ
  ! 壁面法線方向に働く力 -------------------------------------------------------------------------------
  Fx = SFx + CFx + GFx; Fy = SFy + CFy + GFy; Fz = SFz + CFz + GFz
  ! 周囲 8点の単位位置ベクトル -------------------------------------------------------------------------
  ip(1) = i+1; kp(1) = k+1; ip(2) = i  ; kp(2) = k+1; ip(3) = i-1; kp(3) = k+1
  ip(4) = i+1; kp(4) = k  ;                           ip(5) = i-1; kp(5) = k
  ip(6) = i+1; kp(6) = k-1; ip(7) = i  ; kp(7) = k-1; ip(8) = i-1; kp(8) = k-1
  do m = 1, 8
   px(m) = x(ip(m), kp(m)) - x(i,k)
   py(m) = y(ip(m), kp(m)) - y(i,k)
   pz(m) = z(ip(m), kp(m)) - z(i,k)
   pa    = sqrt( px(m)**2 + py(m)**2 + pz(m)**2 )
   if( pa /= 0.0 ) then
     px(m) = px(m) / pa
     py(m) = py(m) / pa
     pz(m) = pz(m) / pa
    else
     px(m) = 0.0
     py(m) = 0.0
     pz(m) = 0.0
   endif
  enddo
  ! 液滴に働く力のベクトルと位置ベクトルの内積 ---------------------------------------------------------
  Fa = sqrt(Fx**2 + Fy**2 + Fz**2)
  if( Fa /= 0.0 ) then
    Fx = Fx / Fa
    Fy = Fy / Fa
    Fz = Fz / Fa
   else
    cycle
  endif
  do m = 1, 8
   pro(m) = Fx * px(m) + Fy * py(m) + Fz * pz(m)
  enddo
  ! 液滴に働く力が傾いている方向にRunback (ベクトルの成す角が小さい方向) -------------------------------
  mp(1) = maxval( maxloc(pro(:)) )
  mp(2) = maxval( maxloc(pro(:), mask = pro(:) < pro(mp(1))) )
  mp(3) = maxval( maxloc(pro(:), mask = pro(:) < pro(mp(2))) )
  if( pro(mp(1)) >= 0.0 .and. pro(mp(2)) >= 0.0 .and. pro(mp(3)) >= 0.0 ) then
    do m = 1, 3
     wr(m) = pro(mp(m)) / ( pro(mp(1)) + pro(mp(2)) + pro(mp(3)) )
     Mrin(ip(mp(m)), kp(mp(m))) = Mrin(ip(mp(m)), kp(mp(m))) + Mrout(i,k) * wr(m)
    enddo
   else if( pro(mp(1)) >= 0.0 .and. pro(mp(2)) >= 0.0 ) then
    do m = 1, 2
     wr(m) = pro(mp(m)) / ( pro(mp(1)) + pro(mp(2)) )
     Mrin(ip(mp(m)), kp(mp(m))) = Mrin(ip(mp(m)), kp(mp(m))) + Mrout(i,k) * wr(m)
    enddo
   else if( pro(mp(1)) >= 0.0 ) then
    m = 1
    wr(m) = pro(mp(m)) / pro(mp(1))
    Mrin(ip(mp(m)), kp(mp(m))) = Mrin(ip(mp(m)), kp(mp(m))) + Mrout(i,k) * wr(m)
   else
    cycle
  endif
 enddo
 enddo
 ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 deallocate( ip, kp, px, py, pz, pro, mp, wr )
 ! 処理終了 ********************************************************************************************
 return
end subroutine RunbackInExt3D
!*******************************************************************************************************
!******** 氷層厚さ (二次元，Original Messinger)							********
!*******************************************************************************************************
subroutine IceThicknessOrg2D( &
&            is, ie, Ts, SA, Rhoi, FF, Mim, Mrin, Mst, Mac, Hac, dBi )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: Ts(is:ie)
 real   , intent(in)  :: SA(is:ie)
 real   , intent(in)  :: Rhoi(is:ie)
 real   , intent(in)  :: FF(is:ie)
 real   , intent(in)  :: Mim(is:ie), Mrin(is:ie), Mst(is:ie)
 real   , intent(out) :: Mac(is:ie), Hac(is:ie), dBi(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
  Mac(i) = ff(i) * (Mim(i) + Mrin(i) + Mst(i))
  if( Mac(i) > 0.0 ) then
    Hac(i) = Cpi * (Ts(i) - Tf)
   else if( Mac(i) /= 0.0 ) then
    write(*, '(a)') '!!!!! Error : Mac is negative !!!!!'
    write(*, '(i4,5(x,e16.8e3))') i, ff(i), Mim(i), Mrin(i), Mst(i), Mac(i)
    stop
  endif
  dBi(i) = Mac(i) / Rhoi(i) / SA(i)
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine IceThicknessOrg2D
!*******************************************************************************************************
!******** 氷層厚さ (三次元，Original Messinger)							********
!*******************************************************************************************************
subroutine IceThicknessOrg3D( &
&            is, ie, ks, ke, Ts, SA, Rhoi, FF, Mim, Mrin, Mst, Mac, Hac, dBi )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: Ts(is:ie, ks:ke)
 real   , intent(in)  :: SA(is:ie, ks:ke)
 real   , intent(in)  :: Rhoi(is:ie, ks:ke)
 real   , intent(in)  :: FF(is:ie, ks:ke)
 real   , intent(in)  :: Mim(is:ie, ks:ke), Mrin(is:ie, ks:ke), Mst(is:ie, ks:ke)
 real   , intent(out) :: Mac(is:ie, ks:ke), Hac(is:ie, ks:ke), dBi(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  Mac(i,k) = ff(i,k) * (Mim(i,k) + Mrin(i,k) + Mst(i,k))
  if( Mac(i,k) > 0.0 ) then
    Hac(i,k) = Cpi * (Ts(i,k) - Tf)
   else if( Mac(i,k) /= 0.0 ) then
    write(*, '(a)') '!!!!! Error : mac is negative !!!!!'
    write(*, '(2i4,5(x,e16.8e3))') i, k, ff(i,k), Mim(i,k), Mrin(i,k), Mst(i,k), Mac(i,k)
    stop
  endif
  dBi(i,k) = Mac(i,k) / Rhoi(i,k) / SA(i,k)
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine IceThicknessOrg3D
!*******************************************************************************************************
!******** 氷層厚さ (二次元，Extended Messinger)							********
!*******************************************************************************************************
subroutine IceThicknessExt2D( &
&            is, ie, Ts, Mim, Mes, Mrin, Mst, Q0, Q1, Bi, nPhase, time, dt, Bg, Tg, dBi, Bw, Ti, Tw )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie
 real   , intent(in)    :: Ts(is:ie)
 real   , intent(in)    :: Mim(is:ie), Mes(is:ie), Mrin(is:ie), Mst(is:ie)
 real   , intent(in)    :: Q1(is:ie)
 real   , intent(in)    :: Bi(is:ie)
 real   , intent(in)    :: time, dt
 real   , intent(in)    :: Bg(is:ie), Tg(is:ie)
 real   , intent(out)   :: dBi(is:ie), Bw(is:ie)
 real   , intent(out)   :: Ti(is:ie), Tw(is:ie)
 real   , intent(inout) :: Q0(is:ie)
 integer, intent(inout) :: nPhase(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 do i = is + 1, ie - 1
  select case(nPhase(i))
   ! 着氷なし ------------------------------------------------------------------------------------------
   case(0)
    dBi(i) = 0.0
   ! 霧氷 ----------------------------------------------------------------------------------------------
   case(1)
    ! 氷層厚さ
    dBi(i) = ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) ) / Rhor
    if(dBi(i) < 0.0) dBi(i) = 0.0
    ! 氷層温度
!    Q0(i) = Rhor * LF * dBi(i) + Q0(i)
!    Ti(i) = (Q0(i) - Q1(i) * Ts(i)) / (ki + Q1(i) * Bi(i)) * Bi(i) + Ts(i) ! #ver_Hayashi san
    Ti(i) = (Q0(i) - Q1(i)) / ki  * (Bi(i) + dBi(i) * dt) + Ts(i)  !#ver_U
    if(Ti(i) > Tf) dBi(i) = 0.0 !#ver_Uranai 今までの→if(Ti(i) > Ts(i)) Ti(i) = Ts(i)
    ! 水膜厚さ
    if(dBi(i) <= 0.0)THEN
      Bw(i) = ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) ) / Rhow !; nPhase(i) = 2 !#ver_Uranai
     ELSE
      Bw(i) = 0.0
    END IF
   ! 雨氷 ----------------------------------------------------------------------------------------------
   case(2)
!#これまでのCode
!    ! 水膜厚さ
!!    Bw(i) = (Mim(i) + Mrin(i) + Mst(i) - Mes(i)) / rhow * (time - tg(i)) &
!!    &     - rhog / rhow * (Bi(i) - Bg(i))					! 式(32)
!    Bw(i) = ( Mrin(i) + Mst(i) - Mes(i) ) / Rhow * dt				! 林さん
!    if(Bw(i) < 0.0) Bw(i) = 0.0
!    ! 水膜温度
!    Tw(i) = (Q0(i) - Q1(i) * Tf) / (kw + Q1(i) * Bw(i)) * Bw(i) + Tf  !#ver_Hayashi san
!    ! if(Tw(i) > Tf) Tw(i) = Tf 元々は入っていた(2016/12/09)
!    ! 氷層厚さ
!    if( Bi(i) > 0.0 ) then
!!#ver_Hayashi san
!!     dBi(i) = (ki * (Tf - Ts(i)) / Bi(i) - kw * (Q0(i) - Q1(i) * Tf) / (kw + Q1(i) * Bw(i))) / (rhog * LF)
!     else
!      dBi(i) = 0.0
!    endif
!    if(dBi(i) < 0.0) dBi(i) = 0.0
!#Uranai version
    ! 氷層厚さ
    IF( Bi(i) <= 0.0 )THEN
      dBi(i) = ( Q1(i) - Q0(i) ) / (rhog * LF)
     ELSE
      dBi(i) = ( ki * (Tf - Ts(i)) / Bi(i) + (Q1(i) - Q0(i)) ) / (rhog * LF)
    END IF

    if(dBi(i) < 0.0) dBi(i) = 0.0
    ! 水膜厚さ
    Bw(i) = (Mim(i) + Mrin(i) + Mst(i) - Mes(i)) / rhow * (dt - tg(i)) &
    &     - rhog / rhow * ( dBi(i) * dt - Bg(i))		        	! 式(32)
    if(Bw(i) < 0.0) Bw(i) = 0.0
    ! 水膜温度
    if(Bi(i) <= 0.0 .and. dBi(i) <= 0.0) THEN
      Tw(i) = (Q0(i) - Q1(i)) / kw * Bw(i) + Ts(i)
     ELSE
      Tw(i) = (Q0(i) - Q1(i)) / kw * Bw(i) + Tf
    END IF

   ! 例外処理 ------------------------------------------------------------------------------------------
   case default
    write(*, '(a)') '!!!!! Error : Ice phase !!!!!'
  end select
!  if( dBi(i) <= 0.0 .and. Bw(i) <= 0.0 ) nPhase(i) = 0
 enddo
 ! 処理開始 ********************************************************************************************
 return
end subroutine IceThicknessExt2D
!*******************************************************************************************************
!******** 氷層厚さ (二次元，Extended Messinger, 翼面加熱を考慮)					********
!*******************************************************************************************************
subroutine IceThicknessExt2DverUR( &
&            is, ie, Ts, Tfree, Mim, Mes, Mrin, Mst, Mshed, &
&            Q0, Q1, Qcond, Bi, nPhase, dt, dtRK, Bg, Tg, dBi, Bw, Ti, Tw, Tw0, hf, Tb, Twall )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie
 real   , intent(in)    :: Ts(is:ie), hf(is:ie)
 real   , intent(in)    :: Mim(is:ie), Mrin(is:ie), Mst(is:ie)
 real   , intent(in)    :: Tfree
 real   , intent(in)    :: Q0(is:ie), Q1(is:ie), Bi(is:ie)
 real   , intent(in)    :: dt, dtRK
 real   , intent(inout) :: Mes(is:ie), dBi(is:ie), Bg(is:ie), Tg(is:ie), Tb(is:ie), Twall(is:ie)
 real   , intent(out)   :: Bw(is:ie), Mshed(is:ie)
 real   , intent(out)   :: Ti(is:ie)
 real   , intent(out)   :: Qcond(is:ie)
 real   , intent(inout) :: Tw(is:ie), Tw0(is:ie)	! Q0=Q_in, Q1=Q_out
 integer, intent(inout) :: nPhase(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 Real    :: hg, R_cond, ICMNum, Bio, Qs, Ql, Qd, Tws, IceThick, BiRK, dBi_Max
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
  Mshed(i) = 0.0; Qcond(i) = 0.0
  IF( ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) ) <= 0.0 )THEN
    nPhase(i) = 0; Mes(i) = Mim(i) + Mrin(i) + Mst(i)
  END IF
  select case(nPhase(i))
   ! 着氷なし(No Ice) ----------------------------------------------------------------------------------
   case(0)
    dBi(i) = 0.0; Bg(i) = 0.0; Tg(i) = 0.0; Bw(i) = 0.0; Tw(i) = 0.0
    ! 氷層温度(既に氷層が形成されている場合は表面温度を代入)
    if( Bi(i) > 0.0 )then
      Ti(i) = Ts(i)
     else
      Ti(i) = 0.0
    end if
   ! 霧氷(Rime Ice) ------------------------------------------------------------------------------------
   case(1)
    ! 氷層厚さ
    dBi(i) = ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) ) / Rhor
    BiRK = Bi(i) + dBi(i) * dtRK
    ! 氷層温度
    Ti(i) = Tb(i) ! + ( Q0(i) - Q1(i) ) / ki * BiRK
    ! 水膜厚さ
    Bw(i) = 0.0
    ! 水膜温度
    Tw(i) = 0.0
    ! 翼面への伝熱量(冷却を正)
    Qcond(i) = Q1(i) - Q0(i)
   ! 雨氷(Glaze Ice) -----------------------------------------------------------------------------------
   case(2)
    BiRK = Bi(i) + dBi(i) * dtRK
    ! 氷層厚さ
    IF( BiRK <= 0.0 )THEN
      dBi(i) = ( Q1(i) - Q0(i) ) / (Rhog * LF)
     ELSE
      dBi(i) = ( ki * ( Tf - Tb(i) ) / BiRK + ( Q1(i) - Q0(i) ) ) / ( Rhog * LF )
    END IF
    ! 氷層温度
    Ti(i) = Tb(i) !0.5 * ( Tf + Tb(i) )
    dBi_Max = ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) ) / Rhog
    if(dBi(i) < 0.0)then
      dBi(i) = 0.0; Ti(i) = 0.0
     else if( dBi(i) > dBi_Max )then
      dBi(i) = ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) ) / Rhog
    end if
    ! 水膜厚さ
    if( dBi(i) > 0.0 )then
      Bw(i) = ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) ) / Rhow * dt &
      &       - Rhog / Rhow * dBi(i) * dt ! 式(32)
     else
      Bw(i) = ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) ) / Rhow * dtRK
    end if
    ! 水膜温度
    if( dBi(i) <= 0.0 )THEN
      Tw(i) = (Q0(i) - Q1(i)) / kw * Bw(i) + Tb(i)
     ELSE
      Tw(i) = (Q0(i) - Q1(i)) / kw * Bw(i) + Tf
    END IF
    ! 翼面への伝熱量(冷却を正) 
    if( dBi(i) > 0.0 )then
      Qcond(i) = Q1(i) - Q0(i) - Rhog * LF * dBi(i)
     else
      Qcond(i) = Q1(i) - Q0(i)
    end if
   ! 防氷 ----------------------------------------------------------------------------------------------
   case(3)
    dBi(i) = 0.0; Bg(i) = 0.0; Tg(i) = 0.0
    ! 氷層温度(既に氷層が形成されている場合は表面温度を代入)
    if( Bi(i) > 0.0 )then
      Ti(i) = Ts(i)
     else
      Ti(i) = 0.0
    end if
    Bw(i) = ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) ) / Rhow * dt
    Qcond(i) = Q1(i) - Q0(i)
    if( Tw(i) > Tb(i) .and. Qcond(i) > 0.0 )then
      Qcond(i) = ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) ) * Cpw * ( Tb(i) - Tw(i) )
    end if
    Tw(i) = Tb(i)
   ! 例外処理 ------------------------------------------------------------------------------------------
   case default
    write(*, '(a)') '!!!!! Error : Ice phase !!!!!'
  end select
 enddo
 ! 処理開始 ********************************************************************************************
 return
end subroutine IceThicknessExt2DverUR
!*******************************************************************************************************
!******** 氷層厚さ (三次元，Extended Messinger)							********
!*******************************************************************************************************
subroutine IceThicknessExt3D( &
&            is, ie, ks, ke, Ts, Mim, Mes, Mrin, Mst, Q0, Q1, Bi, nPhase, time, dt, dBi, Bw, Ti, Tw )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, ks, ke
 real   , intent(in)    :: Ts(is:ie, ks:ke)
 real   , intent(in)    :: Mim(is:ie, ks:ke), Mes(is:ie, ks:ke), Mrin(is:ie, ks:ke), Mst(is:ie, ks:ke)
 real   , intent(in)    :: Q1(is:ie, ks:ke)
 real   , intent(in)    :: Bi(is:ie, ks:ke)
 real   , intent(in)    :: time, dt
 real   , intent(out)   :: dBi(is:ie, ks:ke), Bw(is:ie, ks:ke)
 real   , intent(out)   :: Ti(is:ie, ks:ke), Tw(is:ie, ks:ke)
 real   , intent(inout) :: Q0(is:ie, ks:ke)
 integer, intent(inout) :: nPhase(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  select case(nPhase(i,k))
   ! 着氷なし ------------------------------------------------------------------------------------------
   case(0)
    dBi(i,k) = 0.0
   ! 霧氷 ----------------------------------------------------------------------------------------------
   case(1)
    ! 氷層厚さ
    dBi(i,k) = ( Mim(i,k) + Mrin(i,k) + Mst(i,k) - Mes(i,k) ) / Rhor
    if(dBi(i,k) < 0.0) dBi(i,k) = 0.0
    ! 氷層温度
    Q0(i,k) = Rhor * LF * dBi(i,k) + Q0(i,k)
    Ti(i,k) = (Q0(i,k) - Q1(i,k) * Ts(i,k)) / (ki + Q1(i,k) * Bi(i,k)) * Bi(i,k) + Ts(i,k)
    if(Ti(i,k) > Ts(i,k)) Ti(i,k) = Ts(i,k)
!    Ti(i,k) = Ts(i,k)
    ! 水膜厚さ
    Bw(i,k) = 0.0
   ! 雨氷 ----------------------------------------------------------------------------------------------
   case(2)
    ! 水膜厚さ
!    Bw(i,k) = ( Mim(i,k) + Mrin(i,k) + Mst(i,k) - Mes(i,k) ) / Rhow *  (time - tg(i,k)) &
!    &       - Rhog / Rhow * (Bi(i,k) - Bg(i,k))
!    Bw(i,k) = ( Mrin(i,k) + Mst(i,k) ) / Rhow * dt
    Bw(i,k) = ( Mrin(i,k) + Mst(i,k) - Mes(i,k) ) / Rhow * dt
    if(Bw(i,k) < 0.0) Bw(i,k) = 0.0
    ! 水膜温度
    Tw(i,k) = (Q0(i,k) - Q1(i,k) * Tf) / (kw + Q1(i,k) * Bw(i,k)) * Bw(i,k) + Tf
    if(Tw(i,k) > Tf) Tw(i,k) = Tf
!    Tw(i,k) = Tf
    ! 氷層厚さ
    if( Bi(i,k) > 0.0 ) then
      dBi(i,k) = ( ki * (Tf - Ts(i,k)) / Bi(i,k) &
!      &          - kw * (Q0(i,k) - Q1(i,k) * Tf) / (kw + Q1(i,k) * Bw(i,k) * FF(I,K)) &
      &          - kw * (Q0(i,k) - Q1(i,k) * Tf) / (kw + Q1(i,k) * Bw(i,k)) &
      &            ) / (Rhog * LF)
     else
      dBi(i,k) = 0.0
    endif
    if(dBi(i,k) < 0.0) dBi(i,k) = 0.0
   ! 例外処理 ------------------------------------------------------------------------------------------
   case default
    write(*, '(a)') '!!!!! Error : Ice phase !!!!!'
  end select
  if( dBi(i,k) <= 0.0 .and. Bw(i,k) <= 0.0 ) nPhase(i,k) = 0
 enddo
 enddo
 ! 処理開始 ********************************************************************************************
 return
end subroutine IceThicknessExt3D
!*******************************************************************************************************
!******** 霧氷 / 雨氷の相変化 (二次元，Extended Messinger)					********
!*******************************************************************************************************
subroutine PhaseChangeExt2D( &
&  is, ie, hc, Tfree, Ts, Ta, Tw, Tw0, Pt, Ufree, &
&  Mim, Mrin, Mst, Mes, time, dt, Bi, Q0, Q1, Bg, tg, nPhase ,i_Pt)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, i_Pt
 real   , intent(in)  :: hc(is:ie)
 real   , intent(in)  :: Tfree, Ufree
 real   , intent(in)  :: Ts(is:ie), Ta(is:ie), Pt(is:ie)
 real   , intent(in)  :: Mim(is:ie), Mrin(is:ie), Mst(is:ie), Mes(is:ie)
 real   , intent(in)  :: time, dt
 real   , intent(inout)  :: Bi(is:ie), Tw(is:ie), Tw0(is:ie)
 real   , intent(out) :: Q0(is:ie), Q1(is:ie)
 real   , intent(out) :: Bg(is:ie)
 real   , intent(out) :: tg(is:ie)
 integer, intent(out) :: nPhase(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 real    :: kai
 real    :: term1, term2, term3
 real    :: term11, term22
 real    :: B_water, C1, C2, Tm !水膜の厚さ 計算用の係数C1,C2 翼面と外気温の平均温度
 real    :: Tws(is:ie)          !水膜の表面温度
 real    :: T_rin               !ランバックしてくる液滴の温度
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
   kai  = 0.622 * hc(i) * LE / ( Cpa * Pt(i) * (1.0 / Pr)**(2.0 / 3.0) )
   Q0(i) = 0.5 * Mim(i) * Ufree**2 &
   &     + 0.5 * rr * hc(i) * Ufree**2 / Cpa &
   &     + Mim(i) * Cpw * Ta(i) &
   &     + hc(i) * Ta(i) & 
   &     + 4.0 * eps * sigr * Ta(i)**4 &
   &     + kai * e0 * Ta(i) & 
   &     + ( Mrin(i) + Mst(i) ) * Cpw * Tf
   Q1(i) = Mim(i) * Cpw &
   &     + hc(i) &
   &     + 4.0 * eps * sigr * Tfree**3 &
   &     + kai * e0 &
   &     + ( Mrin(i) + Mst(i) ) * Cpw
  ! 着氷判定 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( (Mim(i) + Mrin(i) + Mst(i)) > 0.0 ) then
      ! 雨氷が現れる氷層厚さ及ぶ時間 -------------------------------------------------------------------
      term1 = ki * ( Tf - Ts(i) )
      term2 = Rhog * LF * ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) ) / Rhor
!      term2 =  LF * ( Mim(i) + Mrin(i) + Mst(i) - Mes(i) )
      term3 = Q0(i) - Q1(i) * Tf  !#ver_Hayashi san
!      term3 = Q0(i) - Q1(i)  !#ver_U
      Bg(i) = term1 / (term2 + term3)
      Tg(i) = rhor * Bg(i) / (Mim(i) + Mrin(i) + Mst(i) - Mes(i))
!! Kamagata san Version ---------------------------------------------------------------------------------
!   if(Bg(i) >= 0.0)then
!     if(term1 < 0.0)then
!      nPhase(i) = 0	! 着氷なし
!     else if(time >= Tg(i) .and. Bi(i) >= Bg(i))then
!      nPhase(i) = 2	! 雨氷
!     else
!      nPhase(i) = 1	! 霧氷
!     endif
!   else if(term1 < 0.0)then
!     nPhase(i) = 0	! 着氷なし
!   else if((term2 + term3) < 0.0)then
!     nPhase(i) = 1	! 霧氷
!   endif
!! -------------------------------------------------------------------------------- Kamagata san Version
! Hayashi san Version ----------------------------------------------------------------------------------
      if( (term2 + term3) > 0.0 ) then
        Bg(i) = term1 / (term2 + term3)
       else
        nPhase(i) = 1
        cycle
      endif
      if( Bg(i) > 0.0 ) then
        term11 = Rhor * Bg(i)
        term22 = Mim(i) + Mrin(i) + Mst(i) - Mes(i)
        if( term22 >  0.0 ) then
          tg(i) = term11 / term22
         else
          nPhase(i) = 1
          cycle
        endif
        else
        nPhase(i) = 1
        cycle
      endif
      ! 霧氷 / 雨氷の判定 ------------------------------------------------------------------------------
      if(time <= tg(i)) then
        nPhase(i) = 1
       else
        nPhase(i) = 2
      endif
   else
    ! 着氷なし
    nPhase(i) = 0 ; Tws(i) = 0.0
  endif
 enddo

 !求めた水膜温度を代入
 DO i = is + 1, ie - 1
   Tw(i) = Tws(i)
 END DO
 ! 処理終了 ********************************************************************************************
 return
end subroutine PhaseChangeExt2D
!*******************************************************************************************************
!******** 霧氷 / 雨氷の相変化 (二次元，Extended Messinger, 翼面加熱を考慮)			********
!*******************************************************************************************************
subroutine PhaseChangeExt2DverUR( &
&  is, ie, x, y, hc, Tfree, Ts, Ta, Tw, Tw0, Pt, Ufree, Pfree, &
&  Mim, Mrin, Mst, Mes, time, dt, Bi, Qin, Qout, Bg, tg, nPhase ,i_Pt, Tc)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, i_Pt
 real   , intent(in)    :: x(is:ie), y(is:ie), hc(is:ie)
 real   , intent(in)    :: Tfree, Ufree, Pfree
 real   , intent(in)    :: Ts(is:ie), Ta(is:ie), Pt(is:ie)
 real   , intent(in)    :: Mim(is:ie), Mrin(is:ie), Mst(is:ie)
 real   , intent(in)    :: time, dt
 real   , intent(inout) :: Bi(is:ie), Tw(is:ie), Tw0(is:ie), Mes(is:ie)
 real   , intent(out)   :: Qin(is:ie), Qout(is:ie)				! Q0=Q_in, Q1=Q_out
 real   , intent(out)   :: Bg(is:ie), Tc(is:ie)
 real   , intent(out)   :: tg(is:ie)
 integer, intent(out)   :: nPhase(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 real    :: kai, R_cond, ds, Q_latent, a1, a2, Qd, Bice, Tice, mass, Bimax, dtheta1, dtheta2
 real    :: term1, term2, term3
 real    :: term11, term22
 real    :: Twref(is:ie)							! 水膜の代表温度
 real    :: B_water, T_rin							! Runbackした水膜の温度
 real    :: Tsur, Tbsur, Tbinf							! 表面温度，Tbar
 real    :: Pvsur, Pvinf							! 表面と周囲の蒸気圧
 real    :: ICMNum						! 集中熱容量モデルの判定値(Bio数の応用)
 ! 処理開始 ********************************************************************************************
 Qin(i) = 0.0; Qout(i) = 0.0; Bg(i) = 0.0; tg(i) = 0.0
 !$omp parallel do
 do i = is + 1, ie - 1
  ! 着氷判定 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( ( Mim(i) + Mrin(i) + Mst(i) ) > 0.0 ) then !流入する質量流束があるかどうか
    ! RunbackInしてくる水膜の温度
    IF( i < i_Pt )THEN
      T_rin = Tw0(i+1)
     ELSE
      T_rin = Tw0(i-1)
    END IF
    Twref(i) = ( Mrin(i) * T_rin + Mim(i) * Tfree + Mst(i) * Tf ) / ( Mrin(i) + Mim(i) + Mst(i) )
    Tc(i) = ( SQRT( Rhow * Cpw * kw ) * Twref(i) + SQRT( rho_wing * c_wing * k_wing ) * Ts(i) ) &
    &       / ( SQRT( Rhow * Cpw * kw ) + SQRT( rho_wing * c_wing * k_wing ) )
    if( Tc(i) > Tf )then
      nPhase(i) = 3 ! Anti-Icing Mode
      Qin(i) =  0.5 * Mim(i) * Ufree**2 &
      &       + 0.5 * rr * hc(i) * Ufree**2 / Cpa
      Qout(i) =  hc(i) * ( Tc(i) - Tfree ) &
      &        + kai * e0 * ( Tc(i) - Tfree ) &
      &        + ( Mrin(i) + Mim(i) + Mst(i)   ) * Cpw * ( Tc(i) - Twref(i) )
     else
      ! 霧氷 / 雨氷の判定 ------------------------------------------------------------------------------
      kai  = 0.622 * hc(i) * LS / ( Cpa * Pt(i) * (1.0 / Pr)**(2.0 / 3.0) )
      Qin(i) =  0.5 * Mim(i) * Ufree**2 &
      &       + 0.5 * rr * hc(i) * Ufree**2 / Cpa
      Qout(i) =  hc(i) * ( Tc(i) - Tfree ) &
      &        + kai * e0 * ( Tc(i) - Tfree ) &
      &        + eps * sigr * ( Tc(i)**4 - Tfree**4 ) &
      &        + ( Mrin(i) + Mim(i) + Mst(i) ) * Cpw * ( Tc(i) - Twref(i) )
      term1 = ki * ( Tf - Tc(i) )
      term2 = ( Mrin(i) + Mim(i) + Mst(i) ) * LF
      term3 = Qin(i) - Qout(i)
      Bg(i) = term1 / ( term2 + term3 )
      tg(i) = Rhor * Bg(i) / ( Mrin(i) + Mim(i) + Mst(i) )
      Bimax = ( Mrin(i) + Mim(i) + Mst(i) ) * dt / Rhor
      if( Qout(i) - Qin(i) > term2 )then
        if( Bg(i) >= Bimax .and. tg(i) >= dt )then
          nPhase(i) = 1 ! Dry Mode
         else
          nPhase(i) = 2 ! Wet Mode
        end if
       else
        nPhase(i) = 2 ! Wet Mode
      end if
    end if
   else
    ! 着氷なし
    nPhase(i) = 0; Twref(i) = 0.0; Tc(i) = 0.0
  endif
 enddo
 !$omp end parallel do
 !求めた水膜温度を代入
 do i = is + 1, ie - 1
   Tw(i) = Twref(i)
 end do
 ! 処理終了 ********************************************************************************************
 return
end subroutine PhaseChangeExt2DverUR
!*******************************************************************************************************
!******** 霧氷 / 雨氷の相変化 (三次元，Extended Messinger)					********
!*******************************************************************************************************
subroutine PhaseChangeExt3D( &
&            is, ie, ks, ke, Ts, Mim, Mes, Mrin, Mst, Q0, Q1, time, dt, Bg, tg, nPhase )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: Ts(is:ie, ks:ke)
 real   , intent(in)  :: Mim(is:ie, ks:ke), Mes(is:ie, ks:ke), Mrin(is:ie, ks:ke), Mst(is:ie, ks:ke)
 real   , intent(in)  :: Q0(is:ie, ks:ke), Q1(is:ie, ks:ke)
 real   , intent(in)  :: time, dt
 real   , intent(out) :: Bg(is:ie, ks:ke)
 real   , intent(out) :: tg(is:ie, ks:ke)
 integer, intent(out) :: nPhase(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, pointer :: nPK(:)
 integer :: i, k
 real    :: term1, term2, term3
 ! 処理開始 ********************************************************************************************
 ! 着氷判定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  if( Mim(i,k) + Mrin(i,k) + Mst(i,k) > 0.0 ) then
      ! 雨氷が現れる氷層厚さ及ぶ時間 -------------------------------------------------------------------
      term1 = ki * ( Tf - Ts(i,k) )
      term2 = Rhog * LF * ( Mim(i,k) + Mrin(i,k) + Mst(i,k) - Mes(i,k) ) / Rhor
      term3 = Q0(i,k) - Q1(i,k) * Tf
      if( term2 + term3 /= 0.0 ) then
        Bg(i,k) = term1 / (term2 + term3)
       else
        nPhase(i,k) = 1
        cycle
      endif
      if( Bg(i,k) > 0.0 ) then
        term1 = Rhor * Bg(i,k)
        term2 = Mim(i,k) + Mrin(i,k) + Mst(i,k) - Mes(i,k)
        if( term2 >  0.0 ) then
          tg(i,k) = term1 / term2
         else
          nPhase(i,k) = 1
          cycle
        endif
       else
        nPhase(i,k) = 1
        cycle
      endif
      ! 霧氷 / 雨氷の判定 ------------------------------------------------------------------------------
      if( time <= tg(i,k) ) then
        nPhase(i,k) = 1
       else
        nPhase(i,k) = 2
      endif
   else
    ! 着氷なし
    nPhase(i,k) = 0
  endif
 enddo
 enddo
 ! リミッタ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( nPK(ks:ke) )
 do i = is + 1, ie - 1
  nPK(:) = nPhase(i,:)
  do k = ks + 1, ke - 1
   if(nPK(k-1) == 1 .and. nPK(k) == 2 .and. nPK(k+1) == 1) then
     nPhase(i,k) = 1
   endif
  enddo
 enddo
 deallocate(nPK)
 ! 処理終了 ********************************************************************************************
 return
end subroutine PhaseChangeExt3D
!*******************************************************************************************************
!******** 氷の時間成長 (二次元)									********
!*******************************************************************************************************
subroutine IceTimeGrowth2D( &
&            is, ie, dt, alp, dBi, Bi0, Bi )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: dt
 real   , intent(in)  :: alp
 real   , intent(in)  :: dBi(is:ie)
 real   , intent(in)  :: Bi0(is:ie)
 real   , intent(out) :: Bi(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 real    :: dBiAve
 ! 処理開始 ********************************************************************************************
 do i = is + 1, ie - 1
  if(alp == 0.0) then
    Bi(i) = Bi0(i) + dt * dBi(i)
   else
    dBiAve  = ( dBi(i-1) + 2.0 * dBi(i) + dBi(i+1) ) / 4.0
    Bi(i) = Bi0(i) + dt * dBi(i)
  endif
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine IceTimeGrowth2D
!*******************************************************************************************************
!******** 氷の時間成長 (三次元)									********
!*******************************************************************************************************
subroutine IceTimeGrowth3D( &
&            is, ie, ks, ke, dt, alp, dBi, Bi0, Bi )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: dt
 real   , intent(in)  :: alp
 real   , intent(in)  :: dBi(is:ie, ks:ke)
 real   , intent(in)  :: Bi0(is:ie, ks:ke)
 real   , intent(out) :: Bi(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: B0(:, :)
 integer :: i, k
 real    :: dBiAve
 ! 処理開始 ********************************************************************************************
 do k = ks + 1, ke - 1
 do i = is + 1, ie - 1
  if(alp == 0.0) then
    Bi(i,k) = Bi0(i,k) + dt * dBi(i,k)
   else
    dBiAve  = 0.25 * (dBi(i-1,k) + dBi(i+1,k) + dBi(i,k-1) + dBi(i,k+1)) 
!    dBiAve  = 0.5 * (dBi(i,k-1) + dBi(i,k+1))
!    Bi(i,k) = Bi0(i,k) + dt * min(dBi(i,k), dBiAve * alp)
!    Bi(i,k) = Bi0(i,k) + dt * max( dBiAve / alp, dBi(i,k) )
    Bi(i,k) = Bi0(i,k) + dt * max( dBiAve / alp, (min(dBi(i,k), dBiAve * alp)) )
  endif
 enddo
 enddo
! do k = ks + 2, ke - 2
! do i = is + 2, ie - 2
!  if(alp == 0.0) then
!    Bi(i,k) = Bi0(i,k) + dt * dBi(i,k)
!   else
!    dBiAve  = 0.125 * ( dBi(i-1,k) + dBi(i+1,k) + dBi(i,k-1) + dBi(i,k+1) &
!    &                 + dBi(i-2,k) + dBi(i+2,k) + dBi(i,k-2) + dBi(i,k+2) ) 
!!    dBiAve  = 0.5 * (dBi(i,k-1) + dBi(i,k+1))
!!    Bi(i,k) = Bi0(i,k) + dt * min(dBi(i,k), dBiAve * alp)
!!    Bi(i,k) = Bi0(i,k) + dt * max( dBiAve / alp, dBi(i,k) )
!    Bi(i,k) = Bi0(i,k) + dt * max( dBiAve / alp, (min(dBi(i,k), dBiAve * alp)) )
!  endif
! enddo
! enddo
!!! do k = ks + 1, ke - 1
!!! do i = is + 1, ie - 1
!!!  if(alp == 0.0) then
!!!    Bi(i,k) = Bi0(i,k) + dt * dBi(i,k)
!!!   else
!!!    dBiAve  = 0.125 * ( dBi(i-1,k  ) + dBi(i+1,k  ) + dBi(i,k+1) + dBi(i-1,k+1) &
!!!    &                 + dBi(i-1,k-1) + dBi(i+1,k-1) + dBi(i,k-1) + dBi(i+1,k+1) ) 
!!!!    dBiAve  = 0.5 * (dBi(i,k-1) + dBi(i,k+1))
!!!!    Bi(i,k) = Bi0(i,k) + dt * min(dBi(i,k), dBiAve * alp)
!!!!    Bi(i,k) = Bi0(i,k) + dt * max( dBiAve / alp, dBi(i,k) )
!!!    Bi(i,k) = Bi0(i,k) + dt * max( dBiAve / alp, (min(dBi(i,k), dBiAve * alp)) )
!!!  endif
!!! enddo
!!! enddo
! if(alp /= 0.0) then
!   allocate( B0(is:ie, ks:ke) )
!   B0(:, :) = Bi(:, :)
!   do k = ks + 2, ke - 2
!   do i = is + 2, ie - 2
!    dBiAve  = 0.25 * (B0(i,k-1) + B0(i,k+1) + B0(i,k-2) + B0(i,k+2))
!    Bi(i,k) = max( dBiAve / alp, (min(B0(i,k), dBiAve * alp)) )
!   enddo
!   enddo
!   deallocate(B0)
! endif
 ! 処理終了 ********************************************************************************************
 return
end subroutine IceTimeGrowth3D
!*******************************************************************************************************
!******** 過去の水膜温度保存 (二次元)								********
!*******************************************************************************************************
subroutine SaveWaterTemperature2D( &
&            is, ie, Tw, Tw0 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: Tw(is:ie)
 real   , intent(out) :: Tw0(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 do i = is, ie
  Tw0(i) = Tw(i)
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine SaveWaterTemperature2D
!*******************************************************************************************************
!******** 過去の氷層厚さ保存 (二次元)								********
!*******************************************************************************************************
subroutine SaveIceThickness2D( &
&            is, ie, Bi, Bi0 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie
 real   , intent(in)  :: Bi(is:ie)
 real   , intent(out) :: Bi0(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 do i = is, ie
  Bi0(i) = Bi(i)
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine SaveIceThickness2D
!*******************************************************************************************************
!******** 過去の氷層厚さ保存 (三次元)								********
!*******************************************************************************************************
subroutine SaveIceThickness3D( &
&            is, ie, ks, ke, Bi, Bi0 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: Bi(is:ie, ks:ke)
 real   , intent(out) :: Bi0(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 ! 処理開始 ********************************************************************************************
 do k = ks, ke
 do i = is, ie
  Bi0(i,k) = Bi(i,k)
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine SaveIceThickness3D
!*******************************************************************************************************
!******** アイス・シェディング (全堆積量による判定)						********
!*******************************************************************************************************
subroutine IceSheddingTotal( &
&            is, ie, ks, ke, OmegaX, OmegaY, OmegaZ, x, y, z, up, vp, wp, Utau, Rhop, SA, Bi,  &
&            CFa, SFa, AFa, fShed )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, ks, ke
 real   , intent(in)  :: OmegaX, OmegaY, OmegaZ
 real   , intent(in)  :: x(is:ie, ks:ke), y(is:ie, ks:ke), z(is:ie, ks:ke)
 real   , intent(in)  :: up(is:ie, ks:ke), vp(is:ie, ks:ke), wp(is:ie, ks:ke)
 real   , intent(in)  :: utau(is:ie, ks:ke), Rhop(is:ie, ks:ke)
 real   , intent(in)  :: SA(is:ie, ks:ke), Bi(is:ie, ks:ke)
 real   , intent(out) :: CFa, SFa, AFa
 integer, intent(out) :: fShed
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 real    :: Rhoi, Mi
 real    :: nx, ny, nz, na
 real    :: CFx, CFy, CFz
 real    :: SFx, SFy, SFz
 ! 処理開始 ********************************************************************************************
 ! 氷に作用する力
 CFa = 0.0; SFa = 0.0; AFa = 0.0
 Rhoi = 0.5 * (Rhog + Rhor)
 do k = ks, ke
 do i = is, ie
  if(Bi(i,k) > 0.0) then
    ! 質量
    Mi = Rhoi * SA(i,k) * Bi(i,k)
    ! 遠心力
    CFx = - ( OmegaY * (OmegaX * y(i,k) - OmegaY * x(i,k)) &
    &       - OmegaZ * (OmegaZ * x(i,k) - OmegaX * z(i,k)) )
    CFy = - ( OmegaZ * (OmegaY * z(i,k) - OmegaZ * y(i,k)) &
    &       - OmegaX * (OmegaX * y(i,k) - OmegaY * x(i,k)) )
    CFz = - ( OmegaX * (OmegaZ * x(i,k) - OmegaX * z(i,k)) &
    &       - OmegaY * (OmegaY * z(i,k) - OmegaZ * y(i,k)) )
    CFa = CFa + Mi * sqrt(CFx**2 + CFy**2 + CFz**2)
    ! せん断応力
    nx = up(i,k)
    ny = vp(i,k)
    nz = wp(i,k)
    na = sqrt(nx**2 + ny**2 + nz**2)
    if(na /= 0) then
      nx = nx / na
      ny = ny / na
      nz = nz / na
      SFx = Utau(i,k)**2 * Rhop(i,k) * SA(i,k) * nx
      SFy = Utau(i,k)**2 * Rhop(i,k) * SA(i,k) * ny
      SFz = Utau(i,k)**2 * Rhop(i,k) * SA(i,k) * nz
     else
      SFx = 0.0
      SFy = 0.0
      SFz = 0.0
    endif
    SFa = SFa + sqrt(SFx**2 + SFy**2 + SFz**2)
    ! 付着力
    AFa = AFa + AdhMin_AlRa06 * SA(i,k)
!    AFa = AFa + AdhMax_AlRa06 * SA(i,k)
!    AFa = AFa + AdhAve_AlRa06 * SA(i,k)
   else
    cycle
  endif
 enddo
 enddo
 ! シェディング判定
 if( sqrt(CFa**2 + SFa**2) > AFa ) then
   fShed = 2
  else
   fShed = 0
 endif
 ! 処理終了 ********************************************************************************************
 return
end subroutine IceSheddingTotal
!*******************************************************************************************************
!******** アイス・シェディング (各セルの堆積量による判定)					********
!*******************************************************************************************************
subroutine IceSheddingCell( &
&            is, ie, ks, ke, OmegaX, OmegaY, OmegaZ, x, y, z, up, vp, wp, Utau, Rhop, SA, Bi, &
&            CFa, SFa, AFa, fShed )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, ks, ke
 real   , intent(in)    :: OmegaX, OmegaY, OmegaZ
 real   , intent(in)    :: x (is:ie, ks:ke), y (is:ie, ks:ke), z (is:ie, ks:ke)
 real   , intent(in)    :: up(is:ie, ks:ke), vp(is:ie, ks:ke), wp(is:ie, ks:ke)
 real   , intent(in)    :: utau(is:ie, ks:ke), Rhop(is:ie, ks:ke)
 real   , intent(in)    :: SA(is:ie, ks:ke)
 real   , intent(out)   :: CFa(is:ie, ks:ke), SFa(is:ie, ks:ke), AFa(is:ie, ks:ke)
 integer, intent(inout) :: fShed(is:ie, ks:ke)
 real   , intent(inout) :: Bi(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, pointer :: fShed0(:, :)
 real   , pointer :: Bi0(:, :)
 integer :: i, k
 real    :: Rhoi, Mi
 real    :: nx, ny, nz, na
 real    :: CFx, CFy, CFz
 real    :: SFx, SFy, SFz
 ! 処理開始 ********************************************************************************************
 ! 氷に作用する力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( fShed0(is:ie, ks:ke) )
 CFa(:,:) = 0.0; SFa(:,:) = 0.0; AFa(:,:) = 0.0
 Rhoi = 0.5 * (Rhog + Rhor)
 fShed0(:, :) =  fShed(:, :)
 do k = ks, ke
 do i = is, ie
  if(Bi(i,k) > 0.0) then
!  if(Bi(i,k) > 1.0e-3) then
    ! 質量
    Mi = Rhoi * SA(i,k) * Bi(i,k)
    ! 遠心力
    CFx = - ( OmegaY * (OmegaX * y(i,k) - OmegaY * x(i,k)) &
    &       - OmegaZ * (OmegaZ * x(i,k) - OmegaX * z(i,k)) )
    CFy = - ( OmegaZ * (OmegaY * z(i,k) - OmegaZ * y(i,k)) &
    &       - OmegaX * (OmegaX * y(i,k) - OmegaY * x(i,k)) )
    CFz = - ( OmegaX * (OmegaZ * x(i,k) - OmegaX * z(i,k)) &
    &       - OmegaY * (OmegaY * z(i,k) - OmegaZ * y(i,k)) )
    CFa(i,k) = Mi * sqrt(CFx**2 + CFy**2 + CFz**2)
    ! せん断応力
    nx = up(i,k)
    ny = vp(i,k)
    nz = wp(i,k)
    na = sqrt(nx**2 + ny**2 + nz**2)
    if(na /= 0) then
      nx = nx / na
      ny = ny / na
      nz = nz / na
      SFx = Utau(i,k)**2 * Rhop(i,k) * SA(i,k) * nx
      SFy = Utau(i,k)**2 * Rhop(i,k) * SA(i,k) * ny
      SFz = Utau(i,k)**2 * Rhop(i,k) * SA(i,k) * nz
     else
      SFx = 0.0
      SFy = 0.0
      SFz = 0.0
    endif
    SFa(i,k) = sqrt(SFx**2 + SFy**2 + SFz**2)
    ! 付着力
    AFa(i,k) = AdhMin_AlRa06 * SA(i,k)
!    AFa(i,k) = AFa(i,k) + AdhMax_AlRa06 * SA(i,k)
!    AFa(i,k) = AFa(i,k) + AdhAve_AlRa06 * SA(i,k)
    ! シェディング判定
    if( sqrt(CFa(i,k)**2 + SFa(i,k)**2) > AFa(i,k) ) then
      fShed(i,k) = fShed(i,k) + 1
      Bi(i,k)    = 0.0
    endif
  endif
 enddo
 enddo
 ! 周りのセルの処理 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! allocate( Bi0(is:ie, ks:ke) )
! Bi0(:, :) = Bi(:, :)
 do k = ks + 2, ke - 2
 do i = is + 2, ie - 2
  if(fShed(i,k) > fShed0(i,k)) then
    Bi(i-1,k  ) = 0.0; Bi(i+1,k  ) = 0.0
    Bi(i  ,k-1) = 0.0; Bi(i  ,k+1) = 0.0
    Bi(i-1,k-1) = 0.25 * ( Bi(i-2,k-1) + Bi(i-0,k-1) + Bi(i-1,k-2) + Bi(i-1,k-0) )
    Bi(i-1,k+1) = 0.25 * ( Bi(i-2,k+1) + Bi(i-0,k+1) + Bi(i-1,k+2) + Bi(i-1,k+0) )
    Bi(i+1,k-1) = 0.25 * ( Bi(i+2,k-1) + Bi(i+0,k-1) + Bi(i+1,k-2) + Bi(i+1,k-0) )
    Bi(i+1,k+1) = 0.25 * ( Bi(i+2,k+1) + Bi(i+0,k+1) + Bi(i+1,k+2) + Bi(i+1,k+0) )
  endif
 enddo
 enddo
! deallocate(Bi0)
 deallocate(fShed0)
 ! 処理終了 ********************************************************************************************
 return
end subroutine IceSheddingCell
! 定義終了 *********************************************************************************************
end module Package_Icing
