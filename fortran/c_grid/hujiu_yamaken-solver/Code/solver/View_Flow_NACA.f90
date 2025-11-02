!*******************************************************************************************************
!*******************************************************************************************************
!******** 流れ場可視化プログラム								********
!******** (NACA翼，三次元圧縮性乱流場，重合格子法，高レイノルズ数型 k-eモデル)		  	********
!********					      2013.01.27  PROGRAMMED BY RYOSUKE HAYASHI ********
!********					      2013.07.18     UPDATED BY RYOSUKE HAYASHI ********
!*******************************************************************************************************
!*******************************************************************************************************
program ViewFlow_NACA
 ! モジュール宣言 **************************************************************************************
 use Package_NACA
 use Package_FileIO
 use Package_Flow
 use Package_ViewFlow
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! ファイル名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: ViewFlowFile * 10 = 'ViewFlow'
 ! 共有定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, parameter :: nViewType = 1					! 0 : 途中解, 1 : 定常解
 ! 処理開始 ********************************************************************************************
 ! 検証実験ケース選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call SelectExpCase
 write(*, '(a)') "+++++ Select Exp. case complete. +++++"
 ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call InitialSetting
 write(*, '(a)') "+++++ Initial setting complete. +++++"
 ! 可視化計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call VisualizedFlow
 write(*, '(a)') "+++++ Flow visualized complete. +++++"
 ! 空力計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call airodynamic_performance
 write(*, '(a)') "+++++ Airodynamic performance complete. +++++"
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
   GrdInDir    = bckdir // 'grid//clean//'
   OSGDir      = bckdir // 'overset//clean//'
   FlwIniDir   = bckdir // 'flow//initial//clean//'
   FlwCalInDir = bckdir // 'flow//cal//clean//'
   FlwViewDir  = bckdir // 'flow//view//clean//'
  else
   GrdInDir    = bckdir // 'grid//icing//'
   OSGDir      = bckdir // 'overset//icing//'
   FlwIniDir   = bckdir // 'flow//initial//icing//'
   FlwCalInDir = bckdir // 'flow//cal//icing//'
   FlwViewDir  = bckdir // 'flow//view//icing//'
 endif
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
 integer   :: kp
 character :: fname * 30
 ! 処理開始 ********************************************************************************************
 ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( Flw(ms:me) )
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
 select case(nViewType)
  case(0)
!   do m = ms, me
 m = me
    write(fname, '(a,i6.6,a)') 'Count', nCount, '_'
    call Input_Flux3D( &
    &      trim(FlwCalInDir) // trim(fname) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      ls, le, Flw(m)%qh )
!   enddo
  case(1)
!   do m = ms, me
 m = me
    call Input_Flux3D( &
    &      trim(FlwCalInDir) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      ls, le, Flw(m)%qh )
!   enddo
 end select
 ! C 型格子分割点 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_CtypeGridPoint( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
  &      Flw(m)%i1, Flw(m)%i2, Flw(m)%i3 )
 ! 計算除去点 (氷の中の点) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! m = 1
 m = me
 allocate( IceIn( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ) )
 if( IceStep == 0 ) then
   IceIn(:, :, :) = 0
  else
   call Input_ArrayInt3D( &
   &      trim(OSGDir) // trim(BlkName(m)) // trim(IceInPointFile), strbin, &
   &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
   &      IceIn )
 endif
 kp = nint( 0.5 * Flw(m)%ke )
 do k = Flw(m)%ks, Flw(m)%ke
 do j = Flw(m)%js, Flw(m)%je
 do i = Flw(m)%is, Flw(m)%ie
  if(  k == kp ) cycle
  IceIn(i,j,k) = IceIn(i,j,kp)
 enddo
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine InitialSetting
!*******************************************************************************************************
!******** 可視化計算 										********
!*******************************************************************************************************
subroutine VisualizedFlow
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real  , pointer :: mach(:, :, :), Tt(:, :, :), Pt(:, :, :)
 character :: fname * 20
 integer   :: m
 ! 処理開始 ********************************************************************************************
! do m = ms, me
 m = me
  ! メモリ確保 -----------------------------------------------------------------------------------------
  allocate( mach( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Tt  ( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Pt  ( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ) )
  ! 物理量に変換 ---------------------------------------------------------------------------------------
  call SetPhysics3DKEM( &
  &      Rg, gamma, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
  &      Flw(m)%qh, Flw(m)%jac, &
  &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps )
  ! 粘性係数及び渦粘性係数 -----------------------------------------------------------------------------
  call ViscosityCoefficient3D( &
  &      muSth, TsSth, s1, Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%t, Flw(m)%mu )
  call EddyViscosityCoefficient3D( &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%rho, Flw(m)%kin, Flw(m)%eps, cmu, Flw(m)%mut )
  ! マッハ数，全温，全圧 -------------------------------------------------------------------------------
  call CalMachTtPt3D( &
  &      Rg, gamma, Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, &
  &      mach, Pt, Tt )
  ! 氷の中の点の処理 -----------------------------------------------------------------------------------
  if(m == 1) then
    call ViewIceIn3D( &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      IceIn, Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%mu, &
    &      Flw(m)%kin, Flw(m)%eps, Flw(m)%mut, mach, Pt, Tt )
  endif
  ! ファイル出力 ---------------------------------------------------------------------------------------
  select case(nViewType)
   case(0)
    write(fname, '(a,i6.6,a)') 'Count', nCount, '_'
    call MakeMAVSFile3D( &
    &      trim(FlwViewDir), trim(fname) // trim(BlkName(m)) // trim(ViewFlowFile), strbin, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      rhoRef, aRef, lRef, &
    &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%mu, &
    &      Flw(m)%kin, Flw(m)%eps, Flw(m)%mut, mach, Pt, Tt, &
    &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
   case(1)
    call MakeMAVSFile3D( &
    &      trim(FlwViewDir), trim(BlkName(m)) // trim(ViewFlowFile), strbin, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      rhoRef, aRef, lRef, &
    &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%mu, &
    &      Flw(m)%kin, Flw(m)%eps, Flw(m)%mut, mach, Pt, Tt, &
    &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
  end select
  ! メモリ解放 -----------------------------------------------------------------------------------------
  deallocate( mach, Tt, Pt )
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine VisualizedFlow
!*******************************************************************************************************
!******** 空力性能計算 										********
!*******************************************************************************************************
subroutine airodynamic_performance
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer	:: m
 integer	:: i,j,k
 integer	:: istart,iend
 integer	:: j0,k0
 real		:: nup
 real		:: yp
 real		:: up
 real		:: x1,x2,y1,y2
 real		:: Dp,Df,Lp,Lf
 real		:: dA
 real		:: dx,dy
 real		:: theta
 real		:: press
 real		:: RhoIn,MuIn
 real,allocatable	:: Cp(:)
 real,allocatable	:: utau(:),tau(:)
 real,dimension(0:1)	:: sgn_p,sgn_tau
 real		:: swi_tau
 real		:: Drag,Lift
 real		:: Cd,Cl
 real		:: Re
 real,allocatable	:: yplus(:)

 ! 処理開始 ********************************************************************************************
 m = me

 open(123,file = 'PressureDistribution.dat',status = 'replace')
 open(124,file = './CdCl.dat',status = 'replace')
 open(125,file = './Grid_NACA0012.csv',status = 'replace')
 open(126,file = './CdCl_sgn.dat',status = 'replace')
 open(127,file = './grid.txt',status = 'replace')
 open(128,file = './yplus.csv',status = 'replace')

 istart = Flw(m)%i1
 iend = Flw(m)%i3-1
 j0 = Flw(m)%js
 k0 = Flw(m)%ks

 !array_allocation
 allocate(utau(istart:iend))
 allocate(tau(istart:iend))
 allocate(yplus(istart:iend))
 allocate(Cp(istart:iend))

 !流入条件
 RhoIn = PsExp / (Rg * TsExp)
 MuIn = muSth * (TsExp / TsSth)**1.5 * (TsSth + s1) / (TsExp + s1)

 !摩擦速度
 do i = istart,iend
  nup = Flw(m)%mu(i,j0+1,k) / Flw(m)%rho(i,j0+1,k)
  up = sqrt((Flw(m)%u(i,j0+1,k0) - Flw(m)%u(i,j0,k0))**2.0 &
    &     + (Flw(m)%v(i,j0+1,k0) - Flw(m)%v(i,j0,k0))**2.0)
  yp = sqrt((Flw(m)%x(i,j0+1,k0) - Flw(m)%x(i,j0,k0))**2.0 &
    &     + (Flw(m)%y(i,j0+1,k0) - Flw(m)%y(i,j0,k0))**2.0)
  utau(i) = sqrt(nup * up / yp)
!write(*,*) i,nup
 end do

 !圧力係数
 do i = istart,iend
  Cp(i) = (Flw(m)%p(i,j0,k0) - PsExp) / (0.5 * RhoIn * VelExp**2.0)
  write(123,*) Flw(m)%x(i,j0,k0)/chord,Cp(i)
 end do

 !抗力,揚力の計算
 !初期値設定
 Dp = 0.0
 Df = 0.0
 Lp = 0.0
 Lf = 0.0

 do i = istart,iend
  x1 = (Flw(m)%x(i,j0,k0) - 0.5 * chord) * cos(AOA) &
    & + Flw(m)%y(i,j0,k0) * sin(AOA)
  y1 = -(Flw(m)%x(i,j0,k0) - 0.5 * chord) * sin(AOA) &
    & + Flw(m)%y(i,j0,k0) * cos(AOA)
  x2 = (Flw(m)%x(i+1,j0,k0) - 0.5 * chord) * cos(AOA) &
    & + Flw(m)%y(i+1,j0,k0) * sin(AOA)
  y2 = -(Flw(m)%x(i+1,j0,k0) - 0.5 * chord) * sin(AOA) &
    & + Flw(m)%y(i+1,j0,k0) * cos(AOA)

  write(125,*) Flw(m)%x(i,j0,k0), Flw(m)%y(i,j0,k0), x1, y1
  write(127,*) x1,y1

  yplus(i) = yp * utau(i) / Flw(m)%mu(i,j0,k0) * Flw(m)%rho(i,j0,k0)
!write(*,*) i,yp,utau(i),Flw(m)%mu(i,j0,k0),Flw(m)%rho(i,j0,k0)
  write(128,*) i, yplus(i)

  dx = x2 - x1
  dy = y2 - y1
  dx = dx * lRef
  dy = dy * lRef

  if(dx .ge. 0.0 .and. dy .ge. 0) then
   sgn_p(0) = 1.0
   sgn_p(1) = -1.0
   sgn_tau(0) = 1.0
   sgn_tau(1) = 1.0
  else if(dx .ge. 0.0 .and. dy .le. 0.0) then
   sgn_p(0) = -1.0
   sgn_p(1) = -1.0
   sgn_tau(0) = 1.0
   sgn_tau(1) = -1.0
  else if(dx .le. 0.0 .and. dy .ge. 0.0) then
   sgn_p(0) = 1.0
   sgn_p(1) = 1.0
   sgn_tau(0) = 1.0
   sgn_tau(1) = -1.0
  else if(dx .le. 0.0 .and. dy .le. 0.0) then
   sgn_p(0) = -1.0
   sgn_p(1) = 1.0
   sgn_tau(0) = 1.0
   sgn_tau(1) = 1.0
  end if

  swi_tau = Flw(m)%u(i,j0+1,k0) / abs(Flw(m)%u(i,j0+1,k0))

  dA = sqrt(dx**2.0 + dy**2.0)
  theta = atan(abs(dx) / abs(dy))

  write(126,*) sgn_p(0),swi_tau * sgn_tau(0),sgn_p(1),swi_tau * sgn_tau(1),&
             & sgn_p(0) * cos(theta),sgn_p(1) * sin(theta)

  press = 0.5 * (Flw(m)%P(i,j0,k0) + Flw(m)%P(i+1,j0,k0)) * rhoRef * aRef**2
!write(124,*) i,press

  tau(i) = ((0.5 * (utau(i) + utau(i+1))) * aRef)**2 * &
          & (0.5 * RhoRef * (Flw(m)%rho(i,j0,k0) + Flw(m)%rho(i+1,j0,k0)))
!write(*,*) i,utau(i),aRef,RhoRef,tau(i)
! write(124,*) i,tau(i)

  Dp = Dp + sgn_p(0) * press * cos(theta) * dA
  Df = Df + swi_tau * sgn_tau(0) * tau(i) * sin(theta) * dA
!write(*,*) i,Df,swi_tau,sgn_tau(0),tau(i),sin(theta),dA
!write(*,*) swi_tau,sgn_tau(0),tau(i),sin(theta),dA
!write(*,*) Df

  Lp = Lp + sgn_p(1) * press * sin(theta) * dA
  Lf = Lf + swi_tau * sgn_tau(1) * tau(i) * cos(theta) * dA

 end do

 Drag = Dp + Df
 Lift = Lp + Lf

 Cd = Drag / (0.5 * (RhoIn * RhoRef) * (VelExp * aRef)**2.0 * (chord * lRef))
 Cl = Lift / (0.5 * (RhoIn * RhoRef) * (VelExp * aRef)**2.0 * (chord * lRef))

 Re = (RhoIn * RhoRef) * (VelExp * aRef) * (chord * lRef) / (muIn * rhoRef * aRef * lRef)

!output_data
!write(124,*) rhoin*RhoRef,velexp*aref,chord*lref
  write(124,*) '============================================'
  write(124,*) 'Initial value'
  write(124,*) 'chord length [m]    : ', chord * lRef,chord,lRef
  write(124,*) 'AOA [deg]           : ', AOA / pi * 180.0
  write(124,*) 'velocity [m/s]      : ', VelExp * aRef
  write(124,*) 'density [g/m^3]     : ', RhoIn * RhoRef
  write(124,*) 'viscosity [Pa s]    : ', muIn * rhoRef * aRef * lRef
  write(124,*) 'Reynolds number [-] : ', Re
  write(124,*) '============================================'
  write(124,*) 'Calculation Result'
  write(124,*) '--------------------------------------------'
  write(124,*) 'Drag Coefficient'
  write(124,*) ' Drag by Pressure : ',Dp
  write(124,*) ' Drag by Friction : ',Df
  write(124,*) ' Total Drag : ',Drag
  write(124,*) ' Cd = ',Cd
  write(124,*) '--------------------------------------------'
  write(124,*) 'Lift Coefficient'
  write(124,*) ' Lift by Pressure : ',Lp
  write(124,*) ' Lift by Friction : ',Lf
  write(124,*) ' Total Lift : ',Lift
  write(124,*) ' Cl = ',Cl
  write(124,*) '--------------------------------------------'

 close(123)
 close(124)
 close(125)
 close(126)
 close(127)
 close(128)

end subroutine airodynamic_performance
! 定義終了 *********************************************************************************************
end program ViewFlow_NACA
