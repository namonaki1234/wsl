!*******************************************************************************************************
!*******************************************************************************************************
!******** 流れ場可視化プログラム								********
!******** (NACA翼，三次元圧縮性乱流場，重合格子法，高レイノルズ数型 k-eモデル)		  	********
!********					      2013.01.27  PROGRAMMED BY RYOSUKE HAYASHI ********
!********					      2013.07.18     UPDATED BY RYOSUKE HAYASHI ********
!******** 抗力・揚力の計算			      2017.10.12     UPDATED BY SHO     URANAI  ********
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
 character, parameter :: ViewIcingFile * 9 = 'ViewIcing'
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
 ! メモリ解法 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Deallocating
 write(*, '(a)') "+++++ Calculation complete. +++++"
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
 call Input_CalSetting( './data/GridFlow/' // trim(ND_CalSetFile) // strtxt )
 ! ディレクトリ設定
 if( IceStep == 0 ) then
   write( fname, '(a)') 'clean//'
  else
   write( fname, '(a)') 'icing//'
   IceCalInDir  = resultdir // 'icing//cal//'
 end if
 GrdInDir    = resultdir // 'grid//' // Trim(fname)
 OSGDir      = resultdir // 'overset//' // Trim(fname)
 FlwIniDir   = resultdir // 'flow//initial//' // Trim(fname)
 FlwCalInDir = resultdir // 'flow//cal//' // Trim(fname)
 FlwViewDir  = resultdir // 'flow//view//' // Trim(fname)
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
 character :: head * 5
 character :: fn1 * 50
 character :: fn2 * 50
 integer :: step
 ! 処理開始 ********************************************************************************************
 ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( Flw(ms:me) )
 ! ブロック名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Set_BlockName
 ! 格子解像度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do m = ms, me
  call Input_Resolution3D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(RslFile), strtxt, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke )
 enddo
 ! C 型格子分割点(抗力・揚力の計算のために) ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do m = ms,me
  if(m .eq. me .and. swi_subgrid .ne. 2) cycle
  call Input_CtypeGridPoint( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
  &      Flw(m)%i1, Flw(m)%i2, Flw(m)%i3 )
 end do
 ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do m = ms, me
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
  &         Flw(m)%qh  (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le), &
  &         Flw(m)%ic  (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je), &
  &         Flw(m)%ji  (Flw(m)%is: Flw(m)%ie) )
 enddo
 ! 格子座標 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do m = ms, me
  call Input_Grid3D( &
  &      trim(FlwIniDir) // trim(BlkName(m)) // trim(ND_GrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
 enddo
 ! メトリックス ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do m = ms, me
  call Input_Metrics3D( &
  &      trim(FlwIniDir) // trim(BlkName(m)) // trim(ND_MetFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%jac, &
  &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
  &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
  &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez )
 enddo
 ! 流束関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 select case(nViewType)
  case(0)
   call Input_FlowCalLog( &
   &      trim(FlwCalInDir) // trim(FlwCalLogFile) // strtxt )
   do m = ms, me
    write(fname, '(a,i6.6,a)') 'Count', nCount, '_'
    call Input_Flux3D( &
    &      trim(FlwCalInDir) // trim(fname) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      ls, le, Flw(m)%qh )
   enddo
  case(1)
   do m = ms, me
    call Input_Flux3D( &
    &      trim(FlwCalInDir) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      ls, le, Flw(m)%qh )
   enddo
 end select
 ! 着氷セルのフラグ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do m = ms, me
  if(IceStep == 0) then
!    Flw(m)%ji(:,:) = Flw(m)%js
!    Flw(m)%ic(:,:,:) = 0; Flw(m)%ic(Flw(m)%i1:Flw(m)%i3,Flw(m)%js,:) = 12
!    if(m == 2) then
!      Flw(m)%ic(:,Flw(m)%js,:) = 12
!    endif
   else
!    call Input_ArrayInt3D( &
!    &      trim(IceCalInDir) // trim(BlkName(m)) // trim(IcingCellFile), strdat, &
!    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
!    &      Flw(m)%ic )
!    call Input_ArrayInt2D( &
!    &      trim(IceCalInDir) // trim(BlkName(m)) // trim(IcingIndexFile), strdat, &
!    &      Flw(m)%is, Flw(m)%ie, Flw(m)%ks, Flw(m)%ke, &
!    &      Flw(m)%ji )
  endif
 enddo

 if(swi_ICM .eq. 1) then
  do m = ms,me
   select case (m)
    case (1)
     head = 'Main_'
    case (2)
     head = 'Sub_'
   end select

   open(1,file = './data/STEP.txt',status = 'old')
    read(1,*) step
   close(1)
   if(step .eq. 1) then
    fn1 = './data/GridFlow/IcingCellData/'//trim(head)//'IcingCell'
    fn2 = './data/GridFlow/IcingCellData/'//trim(head)//'IcingIndex'
   else
    fn1 = './result/ICM/'//trim(head)//'IcingCell'
    fn2 = './result/ICM/'//trim(head)//'IcingIndex'
   end if

   call Input_ArrayInt2D( &
   &      trim(fn1), strdat, &
   &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, &
   &      Flw(m)%ic )
   call Input_ArrayInt1D( &
   &      trim(fn2), strdat, &
   &      Flw(m)%is, Flw(m)%ie, &
   &      Flw(m)%ji )
  end do
 end if

 ! 計算除去点 (氷の中の点) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 m = 1
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
 integer   :: i, j, k, kp
 ! 処理開始 ********************************************************************************************
 do m = ms, me
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

   
!    call CdCl

  endif
  if(m == 2) then
    allocate( IceIn( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ) )
    if( IceStep == 0 ) then
      IceIn(:, :, :) = 0
     else
      IceIn(:, :, :) = 0
      do k = Flw(m)%ks, Flw(m)%ke
      do j = Flw(m)%js, Flw(m)%je
      do i = Flw(m)%is, Flw(m)%ie
!      if( Flw(m)%ic(i,j,k) == 1 ) then
!         IceIn(i,j,k) = 1
!        else
!         IceIn(i,j,k) = 0
!       endif
      enddo
      enddo
      enddo
    endif
    kp = nint( 0.5 * Flw(m)%ke )
    do k = Flw(m)%ks, Flw(m)%ke
    do j = Flw(m)%js, Flw(m)%je
    do i = Flw(m)%is, Flw(m)%ie
     if( k == kp ) cycle
     IceIn(i,j,k) = IceIn(i,j,kp)
    enddo
    enddo
    enddo
    call ViewIceIn3D( &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      IceIn, Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%mu, &
    &      Flw(m)%kin, Flw(m)%eps, Flw(m)%mut, mach, Pt, Tt )
   !抗力の計算
    Call forcek
    ! Cp分布の出力
    Call Output_Cp( &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      Flw(m)%x, Flw(m)%y, Flw(m)%u, Flw(m)%p )

  endif

  if(swi_ICM .eq. 0) then
   if(m .eq. ms .and. swi_subgrid .eq. 1) then
    call CdCl(m)
   else if(m .eq. me .and. swi_subgrid .ne. 1) then
    call CdCl(m)
   end if
  else
   if(m .eq. ms .and. swi_subgrid .eq. 1) then
    call CdCl_ICM(m)
   else if(m .eq. me .and. swi_subgrid .ne. 1) then
    call CdCl_ICM(m)
   end if
  end if

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
    ! MicroAVSファイル
    write(fname, '(a)') trim(ViewIcingFile)
    call ViewIceAccretion( trim(fname) )
  end select

 ! Paraviewファイルの出力
 call OutputPara_bin( &
    &      trim(FlwViewDir), trim(BlkName(m)) // trim(ViewFlowFile), &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      rhoRef, aRef, lRef, &
    &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%mu, &
    &      Flw(m)%kin, Flw(m)%eps, Flw(m)%mut, mach, Pt, Tt, &
    &      Flw(m)%x, Flw(m)%y, Flw(m)%z )

  ! メモリ解放 -----------------------------------------------------------------------------------------
  deallocate( mach, Tt, Pt )
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine VisualizedFlow
!*******************************************************************************************************
!******** 抗力・揚力の計算 									********
!*******************************************************************************************************
subroutine forcek
 !変数宣言
 implicit none
 !局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer             :: i1, i2, is, ie, js, je, ks, ke
 integer	     :: mf, i_f, j_f, k_f
 Integer             :: mf1, mf2, osg_is, osg_ie, osg_js, osg_je, osg_ks, osg_ke, AreaNum, i_fs, i_fe
 Integer             :: i_f1, i_f2
 Real,    parameter  :: B_f = 5.5, kappa_f = 0.40, cmu_f = 0.09
 Real,    pointer    :: utau_f(:)
 real :: ip, pavg, CL, cd1, cd2, Rho_air, alpha, area, ip1, ip2, Cd
 real :: F_drag, F_lift, F_DragViscous, delta_x, delta_y, norm_s, norm_x, cos_alpha, sin_alpha
 Real :: dx_j, dy_j, uc_j, vc_j, yp_j, up_j, nup_j, Tau_wall, Force_drag, Force_lift, u_j
 Real :: dudy1_f, dudy2_f !!今回の計算では意味をもたない変数
 character :: fname * 20, HeaterArea * 6
! !配列定義
! !! メモリ確保
! Allocate( utau_f(i1:i2) )
!処理開始
 k_f = 2
 Write(HeaterArea, '(f0.2)') HeatRange
 Open(2, file = trim(FlwViewDir)//'Cd_'//Trim(HeaterArea)//'%.dat', &
 &       form = 'formatted', status = 'replace')
 Open(3, file = trim(FlwViewDir)//'CL_'//Trim(HeaterArea)//'%.dat', &
 &       form = 'formatted', status = 'replace')

 ! 補間位置の入力 --------------------------------------------------------------------------------------
 mf = ms
 call Input_Resolution3D( &
 &      trim(OSGDir) // trim(BlkName(mf)) // trim(OverAreaFile), strtxt, &
 &      osg_is, osg_ie, osg_ks, osg_ke, osg_ks, osg_ke )

! 抗力係数の計算 ---------------------------------------------------------------------------------------
 Do j_f = Flw(mf)%js, Flw(mf)%je - 1

 ! 下面側の算出位置(x座標)の探索
 Do i_f = Flw(mf)%is,  Flw(mf)%i2
  delta_x = ABS( Flw(mf)%x(i_f,j_f,k_f) - Flw(mf)%x(Flw(mf)%i1,js,k_f) )
  ip2 = i_f
  If( delta_x < chord ) EXIT
 End Do
 ! 上面側の算出位置(x座標)の探索
 Do i_f = Flw(mf)%ie, Flw(mf)%i2, -1
  delta_x = ABS( Flw(mf)%x(i_f,j_f,k_f) - Flw(mf)%x(Flw(mf)%i3,js,k_f) )
  ip1 = i_f
  If( delta_x < chord ) EXIT
 End Do

 !上面側
 delta_y = ABS( Flw(mf)%y(ip1,j_f+1,k_f) - Flw(mf)%y(ip1,j_f,k_f) )
 uc_j = Flw(mf)%qh(ip1,j_f,k_f,2) / Flw(mf)%qh(ip1,j_f,k_f,1)
 up_j = Flw(mf)%qh(ip1+1,j_f,k_f,2) / Flw(mf)%qh(ip1+1,j_f,k_f,1)
 u_j = ( up_j - uc_j) / (Flw(mf)%x(ip1+1,j_f,k_f) - Flw(mf)%x(ip1,j_f,k_f) ) &
 &      * ( 2.0 * chord - Flw(mf)%x(ip1,j_f,k_f) ) + uc_j
 If( u_j >= VelExp ) u_j = VelExp
 cd1 = u_j / VelExp * ( 1.0 - u_j / VelExp ) * delta_y
  !下面側
  delta_y = ABS( Flw(mf)%y(ip2,j_f+1,k_f) - Flw(mf)%y(ip2,j_f,k_f) )
  uc_j = Flw(mf)%qh(ip2,j_f,k_f,2)   / Flw(mf)%qh(ip2,j_f,k_f,1)
  up_j = Flw(mf)%qh(ip2-1,j_f,k_f,2) / Flw(mf)%qh(ip2-1,j_f,k_f,1)
  u_j = ( up_j - uc_j) / ( Flw(mf)%x(ip2-1,j_f,k_f) - Flw(mf)%x(ip2,j_f,k_f) ) &
  &      * ( 2.0 * chord - Flw(mf)%x(ip2,j_f,k_f) ) + uc_j
  If( u_j >= VelExp ) u_j = VelExp
  cd2 = u_j / VelExp * ( 1.0 - u_j / VelExp ) * delta_y
  CD = CD + ( cd1 + cd2 )
 End Do

 CD = 2.0 / chord * CD

 F_drag = 0.0; F_lift = 0.0
 Rho_air = PsExp / ( Rg * TsExp )

 Do AreaNum = 1, 5
  !!! 計算範囲の設定 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !!! 1: Mainから補間領域まで(下面側), 2: 補間領域の境目(下面側), 3: Sub領域 +
  !!! 4: 補間領域の境目(上面側), 5: 補間領域からMainまで(上面側)             +
  !!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Select Case(AreaNum)
   Case(1)
    mf1 = ms; mf2 = ms; i_fs = Flw(mf1)%i1; i_fe = osg_is     ; j_f = Flw(mf1)%js
   Case(2)
    mf1 = ms; mf2 = me; i_fs = Flw(mf2)%is; i_fe = Flw(mf2)%is; j_f = Flw(mf2)%js
    i_f1 = osg_is; i_f2 = Flw(mf2)%is
   Case(3)
    mf1 = me; mf2 = me; i_fs = Flw(mf1)%is; i_fe = Flw(mf1)%ie; j_f = Flw(mf1)%js
   Case(4)
    mf1 = me; mf2 = ms; i_fs = Flw(mf1)%ie; i_fe = Flw(mf1)%ie; j_f = Flw(mf1)%js
    i_f1 = Flw(mf1)%ie; i_f2 = osg_ie
   Case(5)
    mf1 = ms; mf2 = ms; i_fs = osg_ie     ; i_fe = Flw(mf1)%i3; j_f = Flw(mf1)%js
   Case Default
    Write(*, '(a)') '****    Error AreaNum    ****'
  End Select

  Do i_f = i_fs, i_fe - 1
   Select Case(AreaNum)
    Case(1,3,5)
     delta_x = Flw(mf1)%x(i_f+1,j_f,k_f) - Flw(mf1)%x(i_f,j_f,k_f)
     delta_y = Flw(mf1)%y(i_f+1,j_f,k_f) - Flw(mf1)%y(i_f,j_f,k_f)
     norm_s = SQRT( delta_x**2 + delta_y**2 )
     norm_x = ABS( delta_x )
     cos_alpha = ( delta_x**2 )/ ( norm_s * norm_x )
     sin_alpha = SQRT( 1.0 - ( cos_alpha**2 ) )
     pavg = 0.5 * ( Flw(mf1)%p(i_f,j_f,k_f) + Flw(mf1)%p(i_f+1,j_f,k_f) )
    Case(2,4)
     delta_x = Flw(mf2)%x(i_f2,j_f,k_f) - Flw(mf1)%x(i_f1,j_f,k_f)
     delta_y = Flw(mf2)%y(i_f2,j_f,k_f) - Flw(mf1)%y(i_f1,j_f,k_f)
     norm_s = SQRT( delta_x**2 + delta_y**2 )
     norm_x = ABS( delta_x )
     cos_alpha = ( delta_x**2 )/ ( norm_s * norm_x )
     sin_alpha = SQRT( 1.0 - ( cos_alpha**2 ) )
     pavg = 0.5 * ( Flw(mf1)%p(i_f1,j_f,k_f) + Flw(mf2)%p(i_f2,j_f,k_f) )
    Case Default
     Write(*, '(a)') '****    Error AreaNum    ****'
   End Select

   F_drag = F_drag + pavg * sin_alpha * sign( 1.0, delta_y ) * norm_s
   F_lift = F_lift - pavg * cos_alpha * sign( 1.0, delta_x ) * norm_s

  End Do
 End Do

  Force_drag = F_lift * sin(AOA) + F_drag * cos(AOA)
  Force_lift = F_lift * cos(AOA) - F_drag * sin(AOA)
  CL = Force_lift / ( 0.5 * Rho_air * (VelExp**2) * chord )
  Cd = Force_drag / ( 0.5 * Rho_air * (VelExp**2) * chord )

! ip = i2 - 1
! Do i_f = i1, ip
!  j_f = j_Ref(i_f,k_f)
!  !++  圧力抗力  ++
!  delta_x = x(i_f+1,j_f,k_f) - x(i_f,j_f,k_f)
!  delta_y = y(i_f+1,j_f,k_f) - y(i_f,j_f,k_f)
!  norm_s = SQRT( delta_x**2 + delta_y**2 )
!  norm_x = ABS( delta_x )
!  cos_alpha = ( delta_x**2 )/ ( norm_s * norm_x )
!  sin_alpha = SQRT( 1.0 - ( cos_alpha**2 ) )
!  alpha = acos( cos_alpha ) * 180.0 / pi
!  pavg = 0.5 * ( p(i_f,j_f,k_f) + p(i_f + 1,j_f,k_f) )
!  area = area + norm_s
!  !++  摩擦抗力  ++
!  !!++  壁近傍の分布  ++
!  dx_j = x(i_f,j_f+1,k_f) - x(i_f,j_f,k_f)
!  dy_j = y(i_f,j_f+1,k_f) - y(i_f,j_f,k_f)
!  uc_j = Flw(mf)%qh(i_f,j_f+1,k_f,2) / Flw(mf)%qh(i_f,j_f+1,k_f,1)
!  vc_j = Flw(mf)%qh(i_f,j_f+1,k_f,3) / Flw(mf)%qh(i_f,j_f+1,k_f,1)
!  !!++  壁関数  ++
!  yp_j   = sqrt(dx_j**2 + dy_j**2)
!  up_j   = sqrt(uc_j**2 + vc_j**2)
!  nup_j  = Flw(mf)%mu(i_f,j_f+1,k_f) / (Flw(mf)%qh(i_f,j_f+1,k_f,1) * Flw(mf)%jac(i_f,j_f+1,k_f))
!  !!++  摩擦速度の計算  ++
!  Call CalUtauS( &
!  &            kappa_f, B_f, yp_j, up_j, nup_j, &
!  &            utau_f(i_f), dudy1_f, dudy2_f &
!  &          )
!  !!++  壁面せん断応力の計算  ++
!  Tau_wall = (Flw(mf)%qh(i_f,j_f+1,k_f,1) * Flw(mf)%jac(i_f,j_f+1,k_f)) * utau_f(i_f)**2
!
!  F_DragViscous = F_DragViscous + Tau_wall * cos_alpha * sign( 1.0, delta_y ) * norm_s
!  F_drag = F_drag + pavg * sin_alpha * sign( 1.0, delta_y ) * norm_s
!  F_lift = F_lift - pavg * cos_alpha * sign( 1.0, delta_x ) * norm_s
!  Write(*, '(i4,f6.2,e11.3)') j_f, alpha, pavg * sin_alpha * sign( 1.0, delta_y ) * norm_s
! End Do
!
! Force_drag = ( F_drag * cos(AOA) + F_lift * sin(AOA) ) / area
! Force_lift = ( F_drag * -1.0 * sin(AOA) + F_lift * cos(AOA) ) / area
!
! cdp1= Force_drag / ( 0.5 * Rho_air * (VelExp**2) )
! cdp2= ( Force_drag + F_DragViscous * cos(AOA) / area ) / ( 0.5 * Rho_air * VelExp**2 )
! clp = Force_lift / ( 0.5 * Rho_air * (VelExp**2) )

 Write(2, '(e16.8e3,x,e16.8e3)') HeatRange, CD
 Write(3, '(e16.8e3,x,e16.8e3)') HeatRange, CL

 Close(2); Close(3)

! !! メモリ解放
! deallocate( utau_f )
end subroutine forcek
!*******************************************************************************************************
!******** Cp分布の出力(速度分布も)								********
!*******************************************************************************************************
subroutine Output_Cp( is, ie, js, je, ks, ke, x, y, u, p)
 ! 変数宣言 ********************************************************************************************
 implicit none
 !局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: is, ie, js, je, ks, ke
 integer	     :: i_f, j_f, k_f, ip
 real,    intent(in) :: x(is:ie,js:je,ks:ke), y(is:ie,js:je,ks:ke)
 real,    intent(in) :: u(is:ie,js:je,ks:ke), p(is:ie,js:je,ks:ke)
 real :: Rho_air, cp, fp, x_min, xp
 character :: HeaterArea * 6
 ! 処理開始 ********************************************************************************************
 Write(HeaterArea, '(f0.2)') HeatRange
 Open(1, file = trim(FlwViewDir)//'Cp_'//Trim(HeaterArea)//'%.dat', &
 &       form = 'formatted', status = 'replace')
 Open(2, file = trim(FlwViewDir)//'U(0.8c)_'//Trim(HeaterArea)//'%.dat', &
 &       form = 'formatted', status = 'replace')
 Open(3, file = trim(FlwViewDir)//'U(0.5c)_'//Trim(HeaterArea)//'%.dat', &
 &       form = 'formatted', status = 'replace')
 Open(4, file = trim(FlwViewDir)//'U(0.1c)_'//Trim(HeaterArea)//'%.dat', &
 &       form = 'formatted', status = 'replace')
 Open(7, file = trim(FlwViewDir)//'U(0.05c)_'//Trim(HeaterArea)//'%.dat', &
 &       form = 'formatted', status = 'replace')
! 圧力係数の計算 ---------------------------------------------------------------------------------------
 j_f = js; k_f = 2
 Rho_air = PsExp / ( Rg * TsExp )
 x_min = minval( x(:,js,k_f) )
 Do i_f = is, ie
  cp = ( p(i_f,j_f,k_f) - PsExp ) / ( 0.5 * Rho_air * VelExp**2 )
  xp = ( x(i_f,j_f,k_f) - x_min ) / ( Chord - x_min )
  Write(1, '(e16.8e3,x,e16.8e3)') xp, cp
 End Do
! 速度分布 ---------------------------------------------------------------------------------------------
 Do j_f = je, js, -1 !! 上面側
  Write(2, '(e16.8e3,x,e16.8e3)') u(ie,j_f,k_f) * aRef, ( y(ie,j_f,k_f) - y(ie,js,k_f) ) * LRef
 End Do

 ! 0.5chordの位置i(x座標)の探索
 fp = 0.5 * chord
 Do i_f = ie, is, -1
  ip = i_f
  If( x(i_f,js,k_f) < fp ) EXIT
 End Do
 Do j_f = je, js, -1 !! 上面側
  Write(3, '(e16.8e3,x,e16.8e3)') u(ip,j_f,k_f) * aRef, ( y(ip,j_f,k_f) - y(ip,js,k_f) ) * LRef
 End Do

 ! 0.1chordの位置i(x座標)の探索
 fp = 0.1 * chord
 Do i_f = ie, is, -1
  ip = i_f
  If( x(i_f,js,k_f) < fp ) EXIT
 End Do
 Do j_f = je, js, -1 !! 上面側
  Write(4, '(e16.8e3,x,e16.8e3)') u(ip,j_f,k_f) * aRef, ( y(ip,j_f,k_f) - y(ip,js,k_f) ) * LRef
 End Do

 ! 0.05chordの位置i(x座標)の探索
 fp = 0.05 * chord
 Do i_f = ie, is, -1
  ip = i_f
  If( x(i_f,js,k_f) < fp ) EXIT
 End Do
 Do j_f = je, js, -1 !! 上面側
  Write(7, '(e16.8e3,x,e16.8e3)') u(ip,j_f,k_f) * aRef, ( y(ip,j_f,k_f) - y(ip,js,k_f) ) * LRef
 End Do

 Close(1); Close(2); Close(3); Close(4); Close(7)
end subroutine Output_Cp
!*******************************************************************************************************
!******** 着氷パラメータの可視化								********
!*******************************************************************************************************
subroutine ViewIceAccretion( fname )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent (in) :: fname*(*)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: x (:, :, :), y (:, :, :), z (:, :, :)
 real   , pointer :: x0(:, :, :), y0(:, :, :), z0(:, :, :)
 real   , pointer :: Ts(:, :, :), CE(:, :, :), Bi(:, :, :), f(:, :, :)
 integer :: m, k_VIA, kRef
 Character :: FlwCleanDir * 40
 ! 処理開始 ********************************************************************************************
 m = 2
 allocate( x (Flw(m)%is: Flw(m)%ie, 0: 1, Flw(m)%ks: Flw(m)%ke), &
 &         y (Flw(m)%is: Flw(m)%ie, 0: 1, Flw(m)%ks: Flw(m)%ke), &
 &         z (Flw(m)%is: Flw(m)%ie, 0: 1, Flw(m)%ks: Flw(m)%ke), &
 &         Ts(Flw(m)%is: Flw(m)%ie, 0: 1, Flw(m)%ks: Flw(m)%ke), &
 &         CE(Flw(m)%is: Flw(m)%ie, 0: 1, Flw(m)%ks: Flw(m)%ke), &
 &         Bi(Flw(m)%is: Flw(m)%ie, 0: 1, Flw(m)%ks: Flw(m)%ke), &
 &         f (Flw(m)%is: Flw(m)%ie, 0: 1, Flw(m)%ks: Flw(m)%ke) )
 allocate( x0 (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
 &         y0 (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
 &         z0 (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke) )
 ! (着氷前の座標)
 m = 2
 FlwCleanDir  = resultdir // 'flow//initial//clean//'
  call Input_Grid3D( &
  &      trim(FlwCleanDir) // trim(BlkName(m)) // trim(ND_GrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      x0, y0, z0 )
 kRef = Flw(m)%ks + int( 0.5 * (Flw(m)%ke - Flw(m)%ks) )
 Do k_VIA = Flw(m)%ks, Flw(m)%ke 
  x (:,0,k_VIA) = x0(:,0,k_VIA) * lRef; x (:,1,k_VIA) = Flw(m)%x(:,0,k_VIA) * lRef
  y (:,0,k_VIA) = y0(:,0,k_VIA) * lRef; y (:,1,k_VIA) = Flw(m)%y(:,0,k_VIA) * lRef
  z (:,0,k_VIA) = z0(:,0,k_VIA) * lRef; z (:,1,k_VIA) = Flw(m)%z(:,0,k_VIA) * lRef
  Ts(:,0,k_VIA) = 0.0 ; CE(:,0,k_VIA) = 0.0 ; Bi(:,0,k_VIA) = 0.0
  Ts(:,1,k_VIA) = 1.0 ; CE(:,1,k_VIA) = 1.0 ; Bi(:,1,k_VIA) = 1.0
  f (:,0,k_VIA) = 0.0
  f (:,1,k_VIA) = 1.0
 End Do
 call MakeMAVSFile3D4Para( &
 &      trim(FlwViewDir), trim(BlkName(m)) // trim(fname), strbin, &
 &      Flw(m)%is, Flw(m)%ie, 0, 1, Flw(m)%ks, Flw(m)%ke, &
 &      Ts, CE, Bi, f, x, y, z )
 deallocate( Ts, CE, Bi, f, x, y, z, x0, y0, z0 )
 ! 処理終了 ********************************************************************************************
 return
end subroutine ViewIceAccretion
!*******************************************************************************************************
!******** MicroAVSファイル作成（三次元）  		                                        ********
!*******************************************************************************************************
subroutine MakeMAVSFile3D4Para( &
&            strdir, strname, ext, is, ie, js, je, ks, ke, &
&            Ts, CE, Bi, f, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strdir*(*), strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 real     , intent(in)  :: Ts(is:ie, js:je, ks:ke), CE(is:ie, js:je, ks:ke), &
 &                         Bi(is:ie, js:je, ks:ke), f(is:ie, js:je, ks:ke)
 real     , intent(in)  :: x (is:ie, js:je, ks:ke), y (is:ie, js:je, ks:ke), z (is:ie, js:je, ks:ke)
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: strbin * 4 = '.bin'
 character, parameter :: strfld * 4 = '.fld'
 integer  , parameter :: ndim   = 3
 integer  , parameter :: nspace = 3
 integer  , parameter :: veclen = 4
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: nrkind, nskip
 real      :: r
 integer   :: i, j, k, n
 ! 処理開始 ********************************************************************************************
 ! ヘッダファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 open(1, file = trim(strdir) // trim(strname) // trim(strfld), form = 'formatted')
  write(1, '(a)')    '# AVS field file'
  write(1, '(a,i1)') 'ndim   = ', ndim
  write(1, '(a,i4)') 'dim1   = ', ie - is + 1
  write(1, '(a,i4)') 'dim2   = ', je - js + 1
  write(1, '(a,i4)') 'dim3   = ', ke - ks + 1
  write(1, '(a)')    'label  = Ts CE Bi flag'
  write(1, '(a,i1)') 'nspace = ', nspace
  write(1, '(a,i2)') 'veclen = ', veclen
  write(1, '(a)')    'data   = float'
  write(1, '(a)')    'field  = irregular'
  select case(ext)
   ! Binary --------------------------------------------------------------------------------------------
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
   ! Aschii --------------------------------------------------------------------------------------------
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
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strdir) // trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) Ts; write(1) CE; write(1) Bi; write(1) f
    write(1) x ; write(1) y ; write(1) z
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strdir) // trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     write(1, '(e16.8e3,6(x,e16.8e3))') Ts(i,j,k), CE(i,j,k), Bi(i,j,k), f(i,j,k), &
     &                                  x(i,j,k), y(i,j,k), z(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine MakeMAVSFile3D4Para
!*******************************************************************************************************
!******** メモリ解放 										********
!*******************************************************************************************************
subroutine Deallocating
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m
 ! 処理開始 ********************************************************************************************
 do m = ms, me
  deallocate( Flw(m)%x  , Flw(m)%y  , Flw(m)%z  , &
  &           Flw(m)%rho, Flw(m)%u  , Flw(m)%v  , Flw(m)%w  , &
  &           Flw(m)%p  , Flw(m)%t  , Flw(m)%mu , &
  &           Flw(m)%kin, Flw(m)%eps, Flw(m)%mut, &
  &           Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
  &           Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
  &           Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
  &           Flw(m)%jac, Flw(m)%qh )
 enddo
 ! 処理終了 ********************************************************************************************
return
end subroutine Deallocating

subroutine CdCl(m)
implicit none
integer,intent(in) :: m
integer	:: i,j,k
integer :: istart,iend
integer	:: n
integer	:: j_cal
real	:: x1,x2,y1,y2
real	:: Dp,Df,Lp,Lf
real	:: A,dA
real	:: dx,dy
real	:: theta
real	:: press
real	:: RhoIn,MuIn
real,allocatable	:: utau(:),tau(:)
real,dimension(0:1)	:: sgn_p,sgn_tau
real	:: swi_tau
real	:: Drag,Lift
real	:: Cd,Cl
real	:: Re

real,allocatable	:: yplus(:)

open(124,file = './data/GridFlow/' // './CdCl_ICM.dat',status = 'replace')
open(125,file = './data/GridFlow/' // './Grid_NACA0012.csv',status = 'replace')
open(126,file = './data/GridFlow/' // './CdCl_sgn.dat',status = 'replace')
open(127,file = './data/GridFlow/' // './grid.txt',status = 'replace')
open(128,file = './data/GridFlow/' // './yplus.csv',status = 'replace')

!input_utau_data
allocate(utau(Flw(m)%is:Flw(m)%ie))
allocate(tau(Flw(m)%is:Flw(m)%ie))

allocate(yplus(Flw(m)%is:Flw(m)%ie))

if(m .eq. ms) then
 open(123, file = './data/GridFlow/' // 'Main_utau.txt', status = 'old') !check
else if(m .eq. me) then
 open(123, file = './data/GridFlow/' // 'utau.txt', status = 'old')
end if
 do i = Flw(m)%is,Flw(m)%ie
  read(123,*) utau(i)
 end do
close(123)

!calculate_inlet_value
RhoIn = PsExp / (Rg * TsExp)
MuIn = muSth * (TsExp / TsSth)**1.5 * (TsSth + s1) / (TsExp + s1)

!setting_initial_variable
Dp = 0.0
Df = 0.0
Lp = 0.0
Lf = 0.0
A = 0.0

if(m .eq. me .or. swi_subgrid .eq. 2) then
 istart = Flw(m)%i1
 iend = Flw(m)%i3-1
else
 istart = Flw(m)%is
 iend = Flw(m)%ie-1
end if

!calculation_drag_and_lift
do i = istart,iend
 select case(swi_ICM)
  case(0)
   j_cal = Flw(m)%js
  case(1)
  if(Flw(m)%ic(i,Flw(m)%js) .eq. 12) then
   j_cal = Flw(m)%js
  else
   j_cal = Flw(m)%ji(i) + 1
  end if
 end select

 x1 = (Flw(m)%x(i,j_cal,Flw(m)%ks) - 0.5 * chord) * cos(AOA) &
   & + Flw(m)%y(i,j_cal,Flw(m)%ks) * sin(AOA)
 y1 = -(Flw(m)%x(i,j_cal,Flw(m)%ks) - 0.5 * chord) * sin(AOA) &
   & + Flw(m)%y(i,j_cal,Flw(m)%ks) * cos(AOA)
 x2 = (Flw(m)%x(i+1,j_cal,Flw(m)%ks) - 0.5 * chord) * cos(AOA) &
   & + Flw(m)%y(i+1,j_cal,Flw(m)%ks) * sin(AOA)
 y2 = -(Flw(m)%x(i+1,j_cal,Flw(m)%ks) - 0.5 * chord) * sin(AOA) &
   & + Flw(m)%y(i+1,j_cal,Flw(m)%ks) * cos(AOA)
 
 write(125,*) Flw(m)%x(i,j_cal,Flw(m)%ks), Flw(m)%y(i,j_cal,Flw(m)%ks), x1, y1
 write(127,*) x1,y1

 dx = x2 - x1 !Flw(m)%x(i+1,j_cal,Flw(m)%ks) - Flw(m)%x(i,j_cal,Flw(m)%ks)
 dy = y2 - y1 !Flw(m)%y(i+1,j_cal,Flw(m)%ks) - Flw(m)%y(i,j_cal,Flw(m)%ks)

 yplus(i) = sqrt((Flw(m)%x(i,j_cal+1,Flw(m)%ks) - Flw(m)%x(i,j_cal,Flw(m)%ks))**2.0 &
          &    + (Flw(m)%y(i,j_cal+1,Flw(m)%ks) - Flw(m)%y(i,j_cal,Flw(m)%ks))**2.0) &
          &    * utau(i) / Flw(m)%mu(i,j_cal,Flw(m)%ks) * Flw(m)%rho(i,j_cal,Flw(m)%ks)
 write(128,*) i, yplus(i)

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

! write(*,*) Flw(m)%u(i,Flw(m)%js+1,Flw(m)%ks)

 n = 0
 do
  if(Flw(m)%u(i,j_cal+n,Flw(m)%ks) .ne. 0) then
   swi_tau = Flw(m)%u(i,j_cal+n,Flw(m)%ks) / abs(Flw(m)%u(i,j_cal+n,Flw(m)%ks))
   exit
  else
   n = n + 1
!   if(n .ge. Flw(m)%je) then
!    write(*,*) i,Flw(m)%u(i,j_cal+n-1,Flw(m)%ks),Flw(m)%u(i,j_cal+n-2,Flw(m)%ks)
!   end if
  end if
 end do

! if(swi_tau .ne. swi_tau) write(*,*) i,j_cal+n, Flw(m)%u(i,j_cal+n,Flw(m)%ks)


!write(*,*) '---------'
!write(*,*) 'sgnp0',sgn_p(0)
!write(*,*) 'sgnp1',sgn_p(1)
!write(*,*) 'sgntau0',sgn_tau(0)
!write(*,*) 'sgntau1',sgn_tau(1)
!write(*,*) 'swi_tau',swi_tau

! dA = sqrt((Flw(m)%x(i+1,j_cal,Flw(m)%ks) - Flw(m)%x(i,j_cal,Flw(m)%ks))**2.0 &
!       & + (Flw(m)%y(i+1,j_cal,Flw(m)%ks) - Flw(m)%y(i,j_cal,Flw(m)%ks))**2.0)

 dA = sqrt(dx**2.0 + dy**2.0)

 A = A + dA
!write(124,*) A

 theta = atan(abs(dx) / abs(dy)) !0.5 * pi - atan(abs(dy) / abs(dx))
!write(*,*) theta * 180.0 / pi
! write(*,*) dy,dx,atan(abs(dy) / abs(dx))

 write(126,*) sgn_p(0),swi_tau * sgn_tau(0),sgn_p(1),swi_tau * sgn_tau(1),&
            & sgn_p(0) * cos(theta),sgn_p(1) * sin(theta)

 press = 0.5 * (Flw(m)%P(i,j_cal,Flw(m)%ks) + Flw(m)%P(i+1,j_cal,Flw(M)%ks)) * rhoRef * aRef**2
!write(124,*) i,press

! write(*,*) i,i-i1,sgn_p(0),press
 tau(i) = ((0.5 * (utau(i) + utau(i+1))) * aRef)**2 * &
         & (0.5 * RhoRef * (Flw(m)%rho(i,j_cal,Flw(m)%ks) + Flw(m)%rho(i+1,j_cal,Flw(m)%ks)))
! write(124,*) i,tau(i)

 Dp = Dp + sgn_p(0) * press * cos(theta) * dA
 Df = Df + swi_tau * sgn_tau(0) * tau(i) * sin(theta) * dA
!write(*,*) swi_tau,sgn_tau(0),tau(i),sin(theta),dA
!write(*,*) Df

 Lp = Lp + sgn_p(1) * press * sin(theta) * dA
 Lf = Lf + swi_tau * sgn_tau(1) * tau(i) * cos(theta) * dA

end do

!do i = 0,(i2-i1)/2
! write(*,*) P(i1+i,j_cal,ks)-P(i2-i,j_cal,ks)
!end do

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
!write(124,*) muIn,A !1.7608199E-05
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

close(124)
close(125)
close(126)
close(127)
close(128)

end subroutine CdCl


subroutine CdCl_ICM(m)
 implicit none
 !main_routine_variable
 integer,intent(in)	:: m
 !subroutine_variable
 integer		:: nn
 integer		:: i,j,n
 integer		:: ks
 real,allocatable	:: flag(:,:)
 integer		:: start_i,end_i
 integer		:: temp_i,temp_j
 integer,allocatable	:: surf_ij(:,:)
 integer,allocatable	:: normal(:)
 integer		:: num_surf
 integer		:: pi,pj
 real			:: x1,x2,y1,y2
 real			:: Dp,Df,Lp,Lf
 real			:: A,dA
 real			:: dx,dy
 real			:: theta
 real			:: press
 real			:: RhoIn,MuIn
 integer		:: tau_i,tau_j
 real,allocatable	:: utau(:),tau(:)
 real,allocatable	:: yplus(:)
 real,dimension(0:1)	:: sgn_p,sgn_tau
 real			:: swi_tau
 real			:: Drag,Lift
 real			:: Cd,Cl
 real			:: Re

 real			:: yp
 real,parameter		:: cmu = 0.09
 real,parameter		:: kapp = 0.42

 !setting_initial_variable
 Dp = 0.0
 Df = 0.0
 Lp = 0.0
 Lf = 0.0
 A = 0.0

 !calculate_inlet_value
 RhoIn = PsExp / (Rg * TsExp)
 MuIn = muSth * (TsExp / TsSth)**1.5 * (TsSth + s1) / (TsExp + s1)

 nn = (Flw(m)%ie - Flw(m)%is) * 2
 allocate(flag(Flw(m)%is:Flw(m)%ie,Flw(m)%js:Flw(m)%ie))
 allocate(surf_ij(0:nn-1,0:1))
 allocate(normal(0:nn-1))

 if(swi_subgrid .eq. 2) then
  start_i = Flw(m)%i1
  end_i = Flw(m)%i3
 else
  start_i = Flw(m)%is
  end_i = Flw(m)%ie
 end if

 flag = 0
 do i = start_i,end_i-1
  do j = Flw(m)%js,Flw(m)%ji(i)
   flag(i,j) = 1
   flag(i+1,j) = 1
   if(Flw(m)%ic(i,Flw(m)%js) .eq. 12) exit
   flag(i,j+1) = 1
   flag(i+1,j+1) = 1
  end do
 end do

 n = 0
 surf_ij(n,0) = start_i
 surf_ij(n,1) = Flw(m)%js
 temp_i = surf_ij(n,0)
 temp_j = surf_ij(n,1)
 normal(n) = 0
 n = n + 1

 do
  if(flag(temp_i,temp_j+1) .eq. 1) then
   surf_ij(n,0) = temp_i
   surf_ij(n,1) = temp_j + 1
   normal(n) = -1
  else if(flag(temp_i+1,temp_j) .eq. 1) then
   surf_ij(n,0) = temp_i + 1
   surf_ij(n,1) = temp_j
   normal(n) = 0
  else if(flag(temp_i,temp_j-1) .eq. 1) then
   surf_ij(n,0) = temp_i
   surf_ij(n,1) = temp_j - 1
   normal(n) = 1
  else
   write(*,*) 'error: fail to track surface'
   stop
  end if
!write(*,*) n,surf_ij(n,0),surf_ij(n,1)
  flag(temp_i,temp_j) = -1
  temp_i = surf_ij(n,0)
  temp_j = surf_ij(n,1)
  n = n + 1
  if(surf_ij(n-1,0) .eq. end_i) exit
 end do

 num_surf = n

 !output_vtk_file
 ks = Flw(m)%ks
 open(1,file = './data/GridFlow/' // 'ICM_Surface.vtk',status = 'replace')
  write(1,'(A)') '# vtk DataFile Version 3.0'
  write(1,*) 'vtk output'
  write(1,*) 'ASCII'
  write(1,*) 'DATASET POLYDATA'
  write(1,*) 'POINTS ',num_surf,' float'
  do n = 0,num_surf-1
   write(1,*) Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks)*lRef,Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks)*lRef,0.0
  end do
  write(1,*) 'POINT_DATA ', num_surf
  write(1,*) 'SCALARS i float 1'
  write(1,*) 'LOOKUP_TABLE default'
  do n = 0,num_surf-1
   write(1,*) n
  end do
 close(1)

! !input_utau_file
! if(m .eq. ms) then
!  open(1, file = 'Main_utau.txt', status = 'old') !check
! else if(m .eq. me) then
!  open(1, file = 'utau.txt', status = 'old')
! end if
!  do i = Flw(m)%is,Flw(m)%ie
!   read(1,*) utau(i)
!  end do
! close(1)

 !calculate_utau
 allocate(utau(0:num_surf-1))
 allocate(tau(0:num_surf-1))
 allocate(yplus(0:num_surf-1))
 do n = 0,num_surf-1
  select case(normal(n))
   case(0)
    pi = surf_ij(n,0)
    pj = surf_ij(n,1) + 1
   case(-1)
    pi = surf_ij(n,0) - 1
    pj = surf_ij(n,1)
   case(1)
    pi = surf_ij(n,0) + 1
    pj = surf_ij(n,1)
  end select
  yp = sqrt((Flw(m)%x(pi,pj,ks) - Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks))**2.0 + &
       &    (Flw(m)%y(pi,pj,ks) - Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks))**2.0)
  utau(n) = cmu * Flw(m)%kin(pi,pj,ks)**2.0 / (kapp * yp * Flw(m)%eps(pi,pj,ks))
  if(utau(n) .ne. utau(n)) utau(n) = 0.0
  yplus(n) = yp * utau(n) / Flw(m)%mu(surf_ij(n,0),surf_ij(n,1),ks) * Flw(m)%rho(surf_ij(n,0),surf_ij(n,1),ks)
 end do


 !calculation_drag_and_lift
 do n = 0,num_surf-1-1
  x1 = (Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks) - 0.5 * chord) * cos(AOA) &
    & + Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks) * sin(AOA)
  y1 = -(Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks) - 0.5 * chord) * sin(AOA) &
    & + Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks) * cos(AOA)
  x2 = (Flw(m)%x(surf_ij(n+1,0),surf_ij(n+1,1),ks) - 0.5 * chord) * cos(AOA) &
    & + Flw(m)%y(surf_ij(n+1,0),surf_ij(n+1,1),ks) * sin(AOA)
  y2 = -(Flw(m)%x(surf_ij(n+1,0),surf_ij(n+1,1),ks) - 0.5 * chord) * sin(AOA) &
    & + Flw(m)%y(surf_ij(n+1,0),surf_ij(n+1,1),ks) * cos(AOA)

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

  if(normal(n) .eq. -1) then
   tau_i = surf_ij(n,0) - 1
   tau_j = surf_ij(n,1)
  else if(normal(n) .eq. 1) then
   tau_i = surf_ij(n,0) + 1
   tau_j = surf_ij(n,1)
  else if(normal(n) .eq. 0) then
   tau_i = surf_ij(n,0)
   tau_j = surf_ij(n,1) + 1
  end if

  if(Flw(m)%u(tau_i,tau_j,ks) .ne. 0.0) then
   swi_tau = Flw(m)%u(tau_i,tau_j,ks) / abs(Flw(m)%u(tau_i,tau_j,ks))
  else
   swi_tau = 0.0
  end if

  dA = sqrt(dx**2.0 + dy**2.0)

  A = A + dA
  theta = atan(abs(dx) / abs(dy))
  press = 0.5 * (Flw(m)%P(surf_ij(n,0),surf_ij(n,1),ks) + Flw(m)%P(surf_ij(n+1,0),surf_ij(n+1,1),ks)) * rhoRef * aRef**2


  tau(n) = ((0.5 * (utau(n) + utau(n+1))) * aRef)**2 * &
          & (0.5 * RhoRef * (Flw(m)%rho(surf_ij(n,0),surf_ij(n,1),ks) + Flw(m)%rho(surf_ij(n+1,0),surf_ij(n+1,1),ks)))

  Dp = Dp + sgn_p(0) * press * cos(theta) * dA
  Df = Df + swi_tau * sgn_tau(0) * tau(n) * sin(theta) * dA
  Lp = Lp + sgn_p(1) * press * sin(theta) * dA
  Lf = Lf + swi_tau * sgn_tau(1) * tau(n) * cos(theta) * dA

 end do

 Drag = Dp + Df
 Lift = Lp + Lf

 Cd = Drag / (0.5 * (RhoIn * RhoRef) * (VelExp * aRef)**2.0 * (chord * lRef))
 Cl = Lift / (0.5 * (RhoIn * RhoRef) * (VelExp * aRef)**2.0 * (chord * lRef))

 Re = (RhoIn * RhoRef) * (VelExp * aRef) * (chord * lRef) / (muIn * rhoRef * aRef * lRef)

 open(1,file = './data/GridFlow/' // './CdCl_ICM.dat',status = 'replace')
  write(1,*) '============================================'
  write(1,*) 'Initial value'
  write(1,*) 'chord length [m]    : ', chord * lRef,chord,lRef
  write(1,*) 'AOA [deg]           : ', AOA / pi * 180.0
  write(1,*) 'velocity [m/s]      : ', VelExp * aRef
  write(1,*) 'density [g/m^3]     : ', RhoIn * RhoRef
  write(1,*) 'viscosity [Pa s]    : ', muIn * rhoRef * aRef * lRef
  write(1,*) 'Reynolds number [-] : ', Re
 !write(124,*) muIn,A !1.7608199E-05
  write(1,*) '============================================'
  write(1,*) 'Calculation Result'
  write(1,*) '--------------------------------------------'
  write(1,*) 'Drag Coefficient'
  write(1,*) ' Drag by Pressure : ',Dp
  write(1,*) ' Drag by Friction : ',Df
  write(1,*) ' Total Drag : ',Drag
  write(1,*) ' Cd = ',Cd
  write(1,*) '--------------------------------------------'
  write(1,*) 'Lift Coefficient'
  write(1,*) ' Lift by Pressure : ',Lp
  write(1,*) ' Lift by Friction : ',Lf
  write(1,*) ' Total Lift : ',Lift
  write(1,*) ' Cl = ',Cl
  write(1,*) '--------------------------------------------'
 close(1)

end subroutine CdCl_ICM


subroutine OutputPara_bin( &
&      strdir, strname, is, ie, js, je, ks, ke, &
&      rhoRef, aRef, lRef, rho, u, v, w, Ps, Ts, mu, &
&      kin, eps, mut, mac, Pt, Tt, &
&      x, y, z )
 implicit none
 !mainroutine_variable
 character, intent(in)  :: strdir*(*), strname*(*)
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 real     , intent(in)  :: rhoRef, aRef, lRef
 real     , intent(in)  :: rho(is:ie, js:je, ks:ke), &
 &                         u  (is:ie, js:je, ks:ke), v (is:ie, js:je, ks:ke), w (is:ie, js:je, ks:ke), &
 &                         Ps (is:ie, js:je, ks:ke), Ts(is:ie ,js:je, ks:ke), mu(is:ie ,js:je, ks:ke), &
 &                         kin(is:ie, js:je, ks:ke), eps(is:ie, js:je, ks:ke), &
 &                         mut(is:ie, js:je, ks:ke), &
 &                         mac(is:ie, js:je, ks:ke), Pt(is:ie, js:je, ks:ke), Tt(is:ie, js:je, ks:ke)
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
   write(1) x(ni,nj,nk)		* lRef,	y(ni,nj,nk)		* lRef,	z(ni,nj,nk)		* lRef
   write(1) x(ni+1,nj,nk)	* lRef,	y(ni+1,nj,nk)		* lRef,	z(ni+1,nj,nk)		* lRef
   write(1) x(ni+1,nj+1,nk)	* lRef,	y(ni+1,nj+1,nk)		* lRef,	z(ni+1,nj+1,nk)		* lRef
   write(1) x(ni,nj+1,nk)	* lRef,	y(ni,nj+1,nk)		* lRef,	z(ni,nj+1,nk)		* lRef
   write(1) x(ni,nj,nk+1)	* lRef,	y(ni,nj,nk+1)		* lRef,	z(ni,nj,nk+1)		* lRef
   write(1) x(ni+1,nj,nk+1)	* lRef,	y(ni+1,nj,nk+1)		* lRef,	z(ni+1,nj,nk+1)		* lRef
   write(1) x(ni+1,nj+1,nk+1)	* lRef,	y(ni+1,nj+1,nk+1)	* lRef,	z(ni+1,nj+1,nk+1)	* lRef
   write(1) x(ni,nj+1,nk+1)	* lRef,	y(ni,nj+1,nk+1)		* lRef,	z(ni,nj+1,nk+1)		* lRef
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
  write(1) 'VECTORS velocity float'//newline
  do n = 0, npoint-1
   ni = mod(mod(n, (ie - is) * (je - js)), (ie - is))
   nj = int(mod(n, (ie - is) * (je - js)) / (ie - is))
   nk = int(n / ((ie - is) * (je - js)))
   write(1) u(ni,nj,nk)		* aRef,	v(ni,nj,nk)		* aRef,	w(ni,nj,nk)		* aRef
   write(1) u(ni+1,nj,nk)	* aRef,	v(ni+1,nj,nk)		* aRef,	w(ni+1,nj,nk)		* aRef
   write(1) u(ni+1,nj+1,nk)	* aRef,	v(ni+1,nj+1,nk)		* aRef,	w(ni+1,nj+1,nk)		* aRef
   write(1) u(ni,nj+1,nk)	* aRef,	v(ni,nj+1,nk)		* aRef,	w(ni,nj+1,nk)		* aRef
   write(1) u(ni,nj,nk+1)	* aRef,	v(ni,nj,nk+1)		* aRef,	w(ni,nj,nk+1)		* aRef
   write(1) u(ni+1,nj,nk+1)	* aRef,	v(ni+1,nj,nk+1)		* aRef,	w(ni+1,nj,nk+1)		* aRef
   write(1) u(ni+1,nj+1,nk+1)	* aRef,	v(ni+1,nj+1,nk+1)	* aRef,	w(ni+1,nj+1,nk+1)	* aRef
   write(1) u(ni,nj+1,nk+1)	* aRef,	v(ni,nj+1,nk+1)		* aRef,	w(ni,nj+1,nk+1)		* aRef
  end do
  write(1) newline
  write(1) 'SCALARS rho float'//newline
  write(1) 'LOOKUP_TABLE default'//newline
  do n = 0, npoint-1
   ni = mod(mod(n, (ie - is) * (je - js)), (ie - is))
   nj = int(mod(n, (ie - is) * (je - js)) / (ie - is))
   nk = int(n / ((ie - is) * (je - js)))
   write(1) rho(ni,nj,nk)	* rhoRef
   write(1) rho(ni+1,nj,nk)	* rhoRef
   write(1) rho(ni+1,nj+1,nk)	* rhoRef
   write(1) rho(ni,nj+1,nk)	* rhoRef
   write(1) rho(ni,nj,nk+1)	* rhoRef
   write(1) rho(ni+1,nj,nk+1)	* rhoRef
   write(1) rho(ni+1,nj+1,nk+1)	* rhoRef
   write(1) rho(ni,nj+1,nk+1)	* rhoRef
  end do
  write(1) newline
  write(1) 'SCALARS Ps float'//newline
  write(1) 'LOOKUP_TABLE default'//newline
  do n = 0, npoint-1
   ni = mod(mod(n, (ie - is) * (je - js)), (ie - is))
   nj = int(mod(n, (ie - is) * (je - js)) / (ie - is))
   nk = int(n / ((ie - is) * (je - js)))
   write(1) Ps(ni,nj,nk)	* rhoRef * aRef**2.0
   write(1) Ps(ni+1,nj,nk)	* rhoRef * aRef**2.0
   write(1) Ps(ni+1,nj+1,nk)	* rhoRef * aRef**2.0
   write(1) Ps(ni,nj+1,nk)	* rhoRef * aRef**2.0
   write(1) Ps(ni,nj,nk+1)	* rhoRef * aRef**2.0
   write(1) Ps(ni+1,nj,nk+1)	* rhoRef * aRef**2.0
   write(1) Ps(ni+1,nj+1,nk+1)	* rhoRef * aRef**2.0
   write(1) Ps(ni,nj+1,nk+1)	* rhoRef * aRef**2.0
  end do
  write(1) newline
  write(1) 'SCALARS Ts float'//newline
  write(1) 'LOOKUP_TABLE default'//newline
  do n = 0, npoint-1
   ni = mod(mod(n, (ie - is) * (je - js)), (ie - is))
   nj = int(mod(n, (ie - is) * (je - js)) / (ie - is))
   nk = int(n / ((ie - is) * (je - js)))
   write(1) Ts(ni,nj,nk)	* aRef**2.0
   write(1) Ts(ni+1,nj,nk)	* aRef**2.0
   write(1) Ts(ni+1,nj+1,nk)	* aRef**2.0
   write(1) Ts(ni,nj+1,nk) 	* aRef**2.0
   write(1) Ts(ni,nj,nk+1)	* aRef**2.0
   write(1) Ts(ni+1,nj,nk+1)	* aRef**2.0
   write(1) Ts(ni+1,nj+1,nk+1)	* aRef**2.0
   write(1) Ts(ni,nj+1,nk+1)	* aRef**2.0
  end do
  write(1) newline
  write(1) 'SCALARS mu float'//newline
  write(1) 'LOOKUP_TABLE default'//newline
  do n = 0, npoint-1
   ni = mod(mod(n, (ie - is) * (je - js)), (ie - is))
   nj = int(mod(n, (ie - is) * (je - js)) / (ie - is))
   nk = int(n / ((ie - is) * (je - js)))
   write(1) mu(ni,nj,nk)	* rhoRef * aRef * lRef
   write(1) mu(ni+1,nj,nk)	* rhoRef * aRef * lRef
   write(1) mu(ni+1,nj+1,nk)	* rhoRef * aRef * lRef
   write(1) mu(ni,nj+1,nk) 	* rhoRef * aRef * lRef
   write(1) mu(ni,nj,nk+1)	* rhoRef * aRef * lRef
   write(1) mu(ni+1,nj,nk+1)	* rhoRef * aRef * lRef
   write(1) mu(ni+1,nj+1,nk+1)	* rhoRef * aRef * lRef
   write(1) mu(ni,nj+1,nk+1)	* rhoRef * aRef * lRef
  end do
  write(1) newline
  write(1) 'SCALARS k float'//newline
  write(1) 'LOOKUP_TABLE default'//newline
  do n = 0, npoint-1
   ni = mod(mod(n, (ie - is) * (je - js)), (ie - is))
   nj = int(mod(n, (ie - is) * (je - js)) / (ie - is))
   nk = int(n / ((ie - is) * (je - js)))
   write(1) kin(ni,nj,nk)	* aRef**2.0
   write(1) kin(ni+1,nj,nk)	* aRef**2.0
   write(1) kin(ni+1,nj+1,nk)	* aRef**2.0
   write(1) kin(ni,nj+1,nk) 	* aRef**2.0
   write(1) kin(ni,nj,nk+1)	* aRef**2.0
   write(1) kin(ni+1,nj,nk+1)	* aRef**2.0
   write(1) kin(ni+1,nj+1,nk+1)	* aRef**2.0
   write(1) kin(ni,nj+1,nk+1)	* aRef**2.0
  end do
  write(1) newline
  write(1) 'SCALARS eps float'//newline
  write(1) 'LOOKUP_TABLE default'//newline
  do n = 0, npoint-1
   ni = mod(mod(n, (ie - is) * (je - js)), (ie - is))
   nj = int(mod(n, (ie - is) * (je - js)) / (ie - is))
   nk = int(n / ((ie - is) * (je - js)))
   write(1) eps(ni,nj,nk)	* aRef**3.0 / lRef
   write(1) eps(ni+1,nj,nk)	* aRef**3.0 / lRef
   write(1) eps(ni+1,nj+1,nk)	* aRef**3.0 / lRef
   write(1) eps(ni,nj+1,nk) 	* aRef**3.0 / lRef
   write(1) eps(ni,nj,nk+1)	* aRef**3.0 / lRef
   write(1) eps(ni+1,nj,nk+1)	* aRef**3.0 / lRef
   write(1) eps(ni+1,nj+1,nk+1)	* aRef**3.0 / lRef
   write(1) eps(ni,nj+1,nk+1)	* aRef**3.0 / lRef
  end do
  write(1) newline
  write(1) 'SCALARS mut float'//newline
  write(1) 'LOOKUP_TABLE default'//newline
  do n = 0, npoint-1
   ni = mod(mod(n, (ie - is) * (je - js)), (ie - is))
   nj = int(mod(n, (ie - is) * (je - js)) / (ie - is))
   nk = int(n / ((ie - is) * (je - js)))
   write(1) mut(ni,nj,nk)	* rhoRef * aRef * lRef
   write(1) mut(ni+1,nj,nk)	* rhoRef * aRef * lRef
   write(1) mut(ni+1,nj+1,nk)	* rhoRef * aRef * lRef
   write(1) mut(ni,nj+1,nk) 	* rhoRef * aRef * lRef
   write(1) mut(ni,nj,nk+1)	* rhoRef * aRef * lRef
   write(1) mut(ni+1,nj,nk+1)	* rhoRef * aRef * lRef
   write(1) mut(ni+1,nj+1,nk+1)	* rhoRef * aRef * lRef
   write(1) mut(ni,nj+1,nk+1)	* rhoRef * aRef * lRef
  end do
  write(1)newline
  write(1) 'SCALARS mach float'//newline
  write(1) 'LOOKUP_TABLE default'//newline
  do n = 0, npoint-1
   ni = mod(mod(n, (ie - is) * (je - js)), (ie - is))
   nj = int(mod(n, (ie - is) * (je - js)) / (ie - is))
   nk = int(n / ((ie - is) * (je - js)))
   write(1) mac(ni,nj,nk)
   write(1) mac(ni+1,nj,nk)
   write(1) mac(ni+1,nj+1,nk)
   write(1) mac(ni,nj+1,nk)
   write(1) mac(ni,nj,nk+1)
   write(1) mac(ni+1,nj,nk+1)
   write(1) mac(ni+1,nj+1,nk+1)
   write(1) mac(ni,nj+1,nk+1)
  end do
  write(1) newline
  write(1) 'SCALARS Pt float'//newline
  write(1) 'LOOKUP_TABLE default'//newline
  do n = 0, npoint-1
   ni = mod(mod(n, (ie - is) * (je - js)), (ie - is))
   nj = int(mod(n, (ie - is) * (je - js)) / (ie - is))
   nk = int(n / ((ie - is) * (je - js)))
   write(1) Pt(ni,nj,nk)	* rhoRef * aRef**2.0
   write(1) Pt(ni+1,nj,nk)	* rhoRef * aRef**2.0
   write(1) Pt(ni+1,nj+1,nk)	* rhoRef * aRef**2.0
   write(1) Pt(ni,nj+1,nk)	* rhoRef * aRef**2.0
   write(1) Pt(ni,nj,nk+1)	* rhoRef * aRef**2.0
   write(1) Pt(ni+1,nj,nk+1)	* rhoRef * aRef**2.0
   write(1) Pt(ni+1,nj+1,nk+1)	* rhoRef * aRef**2.0
   write(1) Pt(ni,nj+1,nk+1)	* rhoRef * aRef**2.0
  end do
  write(1) 'SCALARS Tt float'//newline
  write(1) 'LOOKUP_TABLE default'//newline
  do n = 0, npoint-1
   ni = mod(mod(n, (ie - is) * (je - js)), (ie - is))
   nj = int(mod(n, (ie - is) * (je - js)) / (ie - is))
   nk = int(n / ((ie - is) * (je - js)))
   write(1) Tt(ni,nj,nk)	* aRef**2.0
   write(1) Tt(ni+1,nj,nk)	* aRef**2.0
   write(1) Tt(ni+1,nj+1,nk)	* aRef**2.0
   write(1) Tt(ni,nj+1,nk) 	* aRef**2.0
   write(1) Tt(ni,nj,nk+1)	* aRef**2.0
   write(1) Tt(ni+1,nj,nk+1)	* aRef**2.0
   write(1) Tt(ni+1,nj+1,nk+1)	* aRef**2.0
   write(1) Tt(ni,nj+1,nk+1)	* aRef**2.0
  end do
 close(1)

end subroutine OutputPara_bin


! 定義終了 *********************************************************************************************
end program ViewFlow_NACA
