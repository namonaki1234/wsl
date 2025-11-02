!*******************************************************************************************************
!*******************************************************************************************************
!******** 流れ場計算プログラム									********
!******** (NACA翼，三次元圧縮性乱流場，重合格子法，高レイノルズ数型 k-eモデル)		  	********
!********					      2012.02.02  PROGRAMMED BY RYOSUKE HAYASHI ********
!********					      2013.07.18     UPDATED BY RYOSUKE HAYASHI ********
!********					      2014.04.15     UPDATED BY RYOSUKE HAYASHI ********
!******** 翼表面に温度分布あり			      2016.10.26     UPDATED BY SHO     URANAI  ********
!******** L2残差を適用				      2018.03.06     UPDATED BY SHO     URANAI  ********
!******** 翼内部の熱伝導計算あり		      2018.10.16     UPDATED BY SHO     URANAI  ********
!*******************************************************************************************************
!*******************************************************************************************************
program FlowField_NACA
 ! モジュール宣言 **************************************************************************************
 use Package_NACA
 use Package_FileIO
 use Package_Flow
 use Package_Icing
 ! 変数宣言 ********************************************************************************************
 implicit none
 real, pointer :: q0_main(:,:,:,:), q0_sub(:,:,:,:)
 ! ファイル名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: ResidualFile * 8 = 'Residual'
 ! 処理開始 ********************************************************************************************
 write(*,*) '==============================================='
 select case(swi_ICM)
 case(0)
  write(*,*) 'Seitch ICM : OFF'
 case(1)
  write(*,*) 'Switch ICM : ON'
 end select
 write(*,*) '==============================================='
 ! 検証実験ケース選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(a)') "<< Exp. Case Selecation >>"
 call SelectExpCase
 ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Initial Setting >>"
 call InitialSetting
 ! 流れ場計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Flow Field Computation >>"
 call CalFlowField
 ! メモリ解法 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Computation Finallized >>"
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
   GrdInDir     = resultdir // 'grid//' // Trim(fname)
   OSGDir       = resultdir // 'overset//' // Trim(fname)
   FlwIniDir    = resultdir // 'flow//initial//' // Trim(fname)
   FlwCalInDir  = resultdir // 'flow//cal//' // Trim(fname)
   FlwCalOutDir = resultdir // 'flow//cal//' // Trim(fname)
 write(*, '(a)') '+++ Icing Step +++'
 write(*, '(a,i2)') '* Ice step      = ', IceStep
 write(*, '(a,i2)') '* Ice step max. = ', IceStepMax
 write(*, '(/,a)') '+++ Exp. Condition +++'
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
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, parameter :: nSmooth = 100
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: m, n, ih
 character :: fname * 30
 character :: head * 5
 character :: fn1 * 50
 character :: fn2 * 50
 integer :: step
 ! 処理開始 ********************************************************************************************
 ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( Flw(ms:me), OSG(ms:me) )
 ! ブロック名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Set_BlockName
 ! 格子解像度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do m = ms, me
  call Input_Resolution3D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(RslFile), strtxt, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke )
 enddo
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
  &         Flw(m)%qh0 (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le), &
  &         Flw(m)%dqh (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le), &
  &         Flw(m)%dqh0(Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le), &
  &         Flw(m)%dqc (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le), &
  &         Flw(m)%dqd (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le), &
  &         Flw(m)%dqp (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le), &
  &         Flw(m)%dqr (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le), &
  &         Flw(m)%dt  (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%HeatFlag (Flw(m)%is: Flw(m)%ie), &
  &         Flw(m)%Res (ls: le), &
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
 nStart = nCount; time = 0.0
 if(nCount == 0) then
   do m = ms, me
    call Input_Flux3D( &
    &      trim(FlwIniDir) // trim(BlkName(m)) // trim(ND_IniFlxFile), strbin, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      ls, le, Flw(m)%qh )
   enddo
  else
   do m = ms, me
    call Input_Flux3D( &
    &      trim(FlwCalInDir) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      ls, le, Flw(m)%qh )
   enddo
 endif
 ! 物理量 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do m = ms, me
  call SetPhysics3DKEM( &
  &      Rg, gamma, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
  &      Flw(m)%qh, Flw(m)%jac, &
  &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps )
 enddo
 ! C 型格子分割点 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do m = ms,me
  if(m .eq. me .and. swi_subgrid .ne. 2) cycle
  call Input_CtypeGridPoint( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
  &      Flw(m)%i1, Flw(m)%i2, Flw(m)%i3 )
 end do

 ! 重合格子補間係数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do m = ms, me
  ! 補間探索範囲 ---------------------------------------------------------------------------------------
  call Input_Resolution3D( &
  &      trim(OSGDir) // trim(BlkName(m)) // trim(OverAreaFile), strtxt, &
  &      OSG(m)%is, OSG(m)%ie, OSG(m)%js, OSG(m)%je, OSG(m)%ks, OSG(m)%ke )
  ! メモリ確保 -----------------------------------------------------------------------------------------
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
  ! 補間係数 -------------------------------------------------------------------------------------------
  call Input_OversetCoe3D( &
  &      trim(OSGDir) // trim(BlkName(m)) // trim(OverCoeFile), strbin, &
  &      OSG(m)%is, OSG(m)%ie, OSG(m)%js, OSG(m)%je, OSG(m)%ks, OSG(m)%ke, &
  &      OSG(m)%fOver, OSG(m)%ip, OSG(m)%jp, OSG(m)%kp, &
  &      OSG(m)%term1, OSG(m)%term2, OSG(m)%term3, OSG(m)%term4, &
  &      OSG(m)%term5, OSG(m)%term6, OSG(m)%term7, OSG(m)%term8 )
 enddo
 ! 加熱位置 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 FlwIniClDir = resultdir // 'flow/initial//clean/'
 DO m=ms,me
  OPEN(1, FILE = trim(FlwIniClDir) // trim(BlkName(m)) // 'HeatFlag.dat', &
   &      FORM = 'formatted', STATUS = 'old')
   DO ih = Flw(m)%is, Flw(m)%ie
    READ(1,'(i4)') Flw(m)%HeatFlag(ih)
   END DO
  CLOSE(1)
 END DO
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
 ! スムージング ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 if(IceStep /= 0) then
  do n  = 1, nSmooth
   do m = ms, me
    call SmoothingFlux3D( &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
    &      Flw(m)%jac, Flw(m)%qh )
   enddo
   call InterpolationOversetGrid
  enddo
 endif
 ! 粗さのフラグ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 m = 1
 allocate( Flw(m)%fRough(Flw(m)%is: Flw(m)%ie, Flw(m)%ks: Flw(m)%ke) )
 Flw(m)%fRough(:, :) = 0
 m = 2
 allocate( Flw(m)%fRough(Flw(m)%is: Flw(m)%ie, Flw(m)%ks: Flw(m)%ke), &
 &         Flw(m)%RH    (Flw(m)%is: Flw(m)%ie, Flw(m)%ks: Flw(m)%ke) )
 if(IceStep > 0) then
   call Input_ArrayInt2D( &
   &      trim(IceCalInDir) // trim(BlkName(m)) // trim(RoughFlagFile), strdat, &
   &      Flw(m)%is, Flw(m)%ie, Flw(m)%ks, Flw(m)%ke, &
   &      Flw(m)%fRough )
  else
   Flw(m)%fRough(:, :) = 0
 endif
 ! 着氷セルのフラグ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
 ! 処理終了 ********************************************************************************************
 return
end subroutine InitialSetting
!*******************************************************************************************************
!******** 流れ場計算 										********
!*******************************************************************************************************
subroutine CalFlowField
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所定数 ********************************************************************************************
 real    , parameter :: dtMax = 1.0e0
 real    , parameter :: dtMin = 1.0e-20
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m, nn, l, j, i
 integer :: nCompCount
 real, pointer :: Order(:,:)
 real    :: dt, dtTmp
 real    :: RK(1:nRunge)
 logical :: fOutputCount, fOutputLog, CalEnd
 integer   :: fnum = 0
 character :: num*4
 ! 処理開始 ********************************************************************************************
 Allocate( Order(ms:me,ls:le), &
 &         q0_main(Flw(ms)%is: Flw(ms)%ie, Flw(ms)%js: Flw(ms)%je, Flw(ms)%ks: Flw(ms)%ke, ls: le), &
 &         q0_sub (Flw(me)%is: Flw(me)%ie, Flw(me)%js: Flw(me)%je, Flw(me)%ks: Flw(me)%ke, ls: le) )
 ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Order(:,:) = 0.0; q0_main(:,:,:,:) = Flw(ms)%qh(:,:,:,:); q0_sub(:,:,:,:) = Flw(me)%qh(:,:,:,:)
 ! 収束判定ファイル出力 --------------------------------------------------------------------------------
 open(20, file = trim(FlwCalOutDir) // trim(BlkName(1)) // 'L2Norm' //strtxt)
 open(21, file = trim(FlwCalOutDir) // trim(BlkName(2)) // 'L2Norm' //strtxt)
 open(22, file = trim(FlwCalOutDir) // trim(BlkName(1)) // 'Res_WallTemperature' //strtxt)
 open(23, file = trim(FlwCalOutDir) // trim(BlkName(2)) // 'Res_WallTemperature' //strtxt)
 write(*, '(a)') '+++ Computational loop start +++'
 ! 流れ場計算ループ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 if(IceStep == 0) then
   nCompCount = nCalMax
  else
   nCompCount = nCalMax * 2.0
 endif
 nStart = 1
 do nCount = nStart, nCompCount
  ! 時間刻み幅 -----------------------------------------------------------------------------------------
  do m = ms, me
   ! 局所時間刻み
   call CalLocalDt3D( &
   &      Cn, Rg, gamma, &
   &      Flw(m)%is , Flw(m)%ie , Flw(m)%js , Flw(m)%je , Flw(m)%ks , Flw(m)%ke , &
   &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
   &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
   &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
   &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%t, Flw(m)%mu, &
   &      Flw(m)%dt )
  enddo
  ! 非定常計算 -----------------------------------------------------------------------------------------
  if(.not. fSteady) then !非定常の場合、メインとサブ両方での最小のdtを出し、全てのdtに代入してる
    dt = dtMax
    Call Calculation_HeatConductionTimeStep(dt)
    do m = ms, me
     dt = min( minval( Flw(m)%dt(:,:,:) ), dt )
    enddo
    do m = ms, me
     Flw(m)%dt(:,:,:) = dt
    enddo
  endif
  ! 過去の流束を保存 -----------------------------------------------------------------------------------
  do m = ms, me !ルンゲクッタ法の最初のQの値のこと
   call SaveFlux3DKEM( &
   &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
   &      Flw(m)%qh, Flw(m)%dqh, Flw(m)%qh0, Flw(m)%dqh0 )
  enddo
  ! 内部反復 -------------------------------------------------------------------------------------------
  do nn = 1, nRunge
   do m = ms, me
    ! 粘性係数
    call ViscosityCoefficient3D( &
    &      muSth, TsSth, s1, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
    &      Flw(m)%t, &
    &      Flw(m)%mu )
    ! 対流項
    call Convection3D( &
    &      nTVD, eTVD, Rg, gamma, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
    &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
    &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
    &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
    &      Flw(m)%jac, Flw(m)%qh, &
    &      Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, &
    &      Flw(m)%dqc )
    ! 拡散項 & 生産項
    select case(TurbNum)
     ! 層流
     case(0)
      call Viscosity3D( &
      &      Rg, gamma, Pr, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
      &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
      &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
      &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
      &      Flw(m)%jac, Flw(m)%qh, &
      &      Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%t, Flw(m)%mu, &
      &      Flw(m)%dqd &
      &          )
      Flw(m)%dqp(:,:,:,:) = 0.0
     ! Launder-Spalding
     case(1)
      call Turbulence3DEvmStd( &
      &	     LmtPro, Rg, gamma, Pr, Prt, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
      &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
      &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
      &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
      &      Flw(m)%jac, Flw(m)%qh, &
      &      Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%t, Flw(m)%kin, Flw(m)%eps, Flw(m)%mu, &
      &      Flw(m)%mut, Flw(m)%dqd, Flw(m)%dqp )
     ! Kato-Launder
     case(2)
      call Turbulence3DEvmStdKL( &
      &	     LmtPro, Rg, gamma, Pr, Prt, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
      &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
      &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
      &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
      &      Flw(m)%jac, Flw(m)%qh, &
      &      Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%t, Flw(m)%kin, Flw(m)%eps, Flw(m)%mu, &
      &      Flw(m)%mut, Flw(m)%dqd, Flw(m)%dqp )
     ! Launder-Sharma
     case(3)
      call Turbulence3DEvmLS( &
      &	     LmtPro, Rg, gamma, Pr, Prt, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
      &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
      &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
      &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
      &      Flw(m)%jac, &
      &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%t, Flw(m)%kin, Flw(m)%eps, Flw(m)%mu, &
      &      Flw(m)%mut, Flw(m)%dqd, Flw(m)%dqp )
     case(4)
      call Turbulence3DEvmANKlRe( &
      &	     LmtPro, Rg, gamma, Pr, Prt, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%i1, Flw(m)%i3, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
      &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
      &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
      &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
      &      Flw(m)%jac, Flw(m)%qh, &
      &      Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%t, Flw(m)%kin, Flw(m)%eps, Flw(m)%mu, &
      &      Flw(m)%mut, Flw(m)%dqd, Flw(m)%dqp, &
      &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
     
     case default; write(*, '(a)') '!!!!! Error : TurbNum !!!!!' 
    end select
    ! 外力項
      Flw(m)%dqr(:,:,:,:) = 0.0
    ! 各項の和
    call SumDQH3D( &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
    &      Flw(m)%dqc, Flw(m)%dqd, Flw(m)%dqp, Flw(m)%dqr, &
    &      Flw(m)%dqh )
    ! 未来の流束
    if(fTime) then
      call RungeKutta3D( &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
      &      nn, nRunge, Flw(m)%dt, &
      &      Flw(m)%dqh, Flw(m)%qh0, Flw(m)%qh )
     else
      call LUADI3D( &
      &      0.5, rg, gamma, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
      &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
      &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
      &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
      &      Flw(m)%jac, &
      &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%t, Flw(m)%mu + Flw(m)%mut, &
      &      Flw(m)%dt, &
      &      Flw(m)%dqh0, Flw(m)%dqh, Flw(m)%qh0, Flw(m)%qh )
    endif
    ! 乱流量のリミッター
    if(LmtAve /= 0.0) then
      call Limiter3DKEM( &
      &      LmtAve, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
      &      Flw(m)%jac, &
      &      Flw(m)%qh )
    endif
   enddo
   ! 境界条件
   call BoundaryCondition
   if(swi_ICM .eq. 1) then
    do m = ms,me
     call Icepoint(m)
    end do
   end if
   ! 物理量 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do m = ms, me
    call SetPhysics3DKEM( &
    &      Rg, gamma, &
    &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
    &      Flw(m)%qh, Flw(m)%jac, &
    &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps )
   enddo
   ! 重合格子の補間
   call InterpolationOversetGrid
  enddo
  ! 翼面における熱移動の計算 ----------------------------------------------------------------------------
!  if( IceStep == 0 ) then
!    !! +++ Main Grid +++
!    Call Calculation_HeatConduction( ms, Flw(ms)%i1, Flw(ms)%i3, dt, nCount )
!    !! +++ Sub Grid +++
!    Call Calculation_HeatConduction( me, Flw(me)%is, Flw(me)%ie, dt, nCount )
!    ! 重合格子の補間
!    call InterpolationOversetGrid
!  end if
  ! 経過時間 -------------------------------------------------------------------------------------------
  if(.not. fSteady) then
    time = time + dt
  endif
  ! ファイル出力 ---------------------------------------------------------------------------------------
  fOutputCount = mod(nCount, nOutputCount) .eq. 0.0 .or. nCount == nStart .or. nCount == nCompCount
  fOutputLog   = mod(nCount, nOutputLog)   .eq. 0.0 .or. nCount == nStart .or. nCount == nCompCount
  if( fOutputLog .or. fOutputCount ) then
    if( nCount == nStart ) then
      write(*, '(a)'        ) '+++ Numerical Condition +++'
      write(*, '(a,e10.4e1)') '* Cn                  = ', Cn
      write(*, '(a,i2)'     ) '* TVD Order           = ', nTVD
      write(*, '(a,e10.4e1)') '* TVD Entropy         = ', eTVD
      write(*, '(a,e10.4e1)') '* Limitter k-e Ave    = ', LmtAve
      write(*, '(a,e10.4e1)') '* Limitter k-e Pro    = ', LmtPro
      write(*, '(a,i7)'     ) '* Computational Count = ', nCompCount
    endif
    if(fOutputCount)then
      call OutputFileCount
      call CalResidual( nCount, CalEnd, Order )
    endif
    write(*, '(a)'     ) '+++ Calculation progress +++'
    write(*, '(a,i6,a)') '* Calculation Count = ', nCount
    if(fSteady) then
      write(*, '(a,e16.8e3)') '* Main Min. dt = ', minval(Flw(1)%dt(:,:,:))
      write(*, '(a,e16.8e3)') '* Main Max. dt = ', maxval(Flw(1)%dt(:,:,:))
      write(*, '(a,e16.8e3)') '* Sub  Min. dt = ', minval(Flw(2)%dt(:,:,:))
      write(*, '(a,e16.8e3)') '* Sub  Max. dt = ', maxval(Flw(2)%dt(:,:,:))
     else
      write(*, '(a,e16.8e3)') '* Calculation time  = ', time
      write(*, '(a,e16.8e3)') '* Calculation dt    = ', dt
    endif

   fnum = fnum + 1
   write(num,'(I4.4)') fnum
   do m = ms, me
    write(*,*) 'Output File : ', trim(BlkName(m)) // 'Flow' // trim(num)
    call OutputPara_bin( &
       &      resultdir // 'flow//view//clean//overtime//', trim(BlkName(m)) // 'Flow' // trim(num), &
       &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
       &      rhoRef, aRef, lRef, &
       &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%mu, &
       &      Flw(m)%kin, Flw(m)%eps, Flw(m)%mut, &
       &      Flw(m)%x, Flw(m)%y, Flw(m)%z )

    if(swi_ICM .eq.0) then
     if(m .eq. me) then
      call calyplus(m)
     end if
     if(m .eq. ms .and. swi_subgrid .eq. 1) then
      call CdCl(m)
     else if(m .eq. me .and. swi_subgrid .ne. 1) then
      call CdCl(m)
     end if
    else
     if(m .eq. ms .and. swi_subgrid .eq. 1) then
      call ICM_data(m)
     else if(m .eq. me .and. swi_subgrid .ne. 1) then
      call ICM_data(m)
     end if
    end if

   end do

  endif
!  if( CalEnd ) then
!    write(*, '(a)')              '<< Flow Calculation Converged >>'
!    write(*, '(a,i8,a)')         '* Calculation Count = ', nCount
!    if(fSteady) then
!      write(*, '(a,e16.8e3)')      '* Main Min. dt = ', minval(Flw(1)%dt(:,:,:))
!      write(*, '(a,e16.8e3)')      '* Main Max. dt = ', maxval(Flw(1)%dt(:,:,:))
!      write(*, '(a,e16.8e3)')      '* Sub  Min. dt = ', minval(Flw(2)%dt(:,:,:))
!      write(*, '(a,e16.8e3)')      '* Sub  Max. dt = ', maxval(Flw(2)%dt(:,:,:))
!     else
!      write(*, '(a,e16.8e3)')      '* Calculation time  = ', time
!      write(*, '(a,e16.8e3)')      '* Calculation dt    = ', dt
!    endif
!    EXIT !
!  end if
 enddo
! if( .not. CalEnd ) write(*, '(a)') '!!! Not Convergence !!!'
 close(20); close(21)
 deallocate( Order, q0_main, q0_sub )
 ! 処理終了 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 return
end subroutine CalFlowField
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

open(700,file='./data/GridFlow/utau.txt',status='replace')
open(701,file = './data/GridFlow/Main_utau.txt',status = 'replace')


 ! Main Grid +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 m = 1
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
 select case(swi_ICM)
  case(0)
   call BoundaryBladeSurface( &
   &      m, Flw(m)%i1 + 0, Flw(m)%i3 - 0, Flw(m)%ks + 1, Flw(m)%ks + 1, &
   &      Flw(m)%js + 0, Flw(m)%js + 1, Flw(m)%js + 2 )
  case(1)
   call BoundaryIcingBladeSurface(m)
 end select
 ! 流出境界 --------------------------------------------------------------------------------------------
 call BoundaryOutlet( &
 &      m, Flw(m)%js + 0, Flw(m)%je - 0, Flw(m)%ks + 1, Flw(m)%ke - 1, &
 &      Flw(m)%is + 0, Flw(m)%is + 1 )
 call BoundaryOutlet( &
 &      m, Flw(m)%js + 0, Flw(m)%je - 0, Flw(m)%ks + 1, Flw(m)%ke - 1, &
 &      Flw(m)%ie - 0, Flw(m)%ie - 1 )
 ! 氷の中の点 ------------------------------------------------------------------------------------------
 call BoundaryIceIn( &
 &      1, 2, Flw(1)%i1, Flw(1)%i3, Flw(1)%js, Flw(1)%je, Flw(1)%ks + 1, Flw(1)%ke - 1, &
 &      Flw(2)%i1, Flw(2)%i2, Flw(2)%js, int(0.5 * Flw(2)%ke) )
 ! 周期境界 --------------------------------------------------------------------------------------------
 call BoundaryPeriodic( &
 &      m, Flw(m)%is + 0, Flw(m)%ie - 0, Flw(m)%js + 0, Flw(m)%je - 0, &
 &      Flw(m)%ks + 0, Flw(m)%ks + 1 )
 call BoundaryPeriodic( &
 &      m, Flw(m)%is + 0, Flw(m)%ie - 0, Flw(m)%js + 0, Flw(m)%je - 0, &
 &      Flw(m)%ke - 1, Flw(m)%ks + 1 )
 call BoundaryPeriodic( &
 &      m, Flw(m)%is + 0, Flw(m)%ie - 0, Flw(m)%js + 0, Flw(m)%je - 0, &
 &      Flw(m)%ke - 0, Flw(m)%ks + 1 )
 ! Sub Grid ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 m = 2
 if(swi_subgrid .eq. 2) then
  ! C 型格子ブランチ・カット ----------------------------------------------------------------------------
  call BoundaryCtypeBranch( &
  &      m, Flw(m)%is + 1, Flw(m)%i1 - 1, Flw(m)%ks + 1, Flw(m)%ke - 1, &
  &      Flw(m)%js + 0, Flw(m)%js + 1, &
  &      Flw(m)%ie )
  call BoundaryCtypeBranch( &
  &      m, Flw(m)%i3 + 1, Flw(m)%ie - 1, Flw(m)%ks + 1, Flw(m)%ke - 1, &
  &      Flw(m)%js + 0, Flw(m)%js + 1, &
  &      Flw(m)%ie )
 end if
 ! 翼壁面 ----------------------------------------------------------------------------------------------
 select case( RoughNum )
  case(1)
   call RoughnessShinBond3D( &
   &      Flw(m)%is, Flw(m)%ie, Flw(m)%ks, Flw(m)%ke, Flw(m)%t * aRef**2, LWC * RhoRef, Flw(m)%RH )
  case default; write(*, '(a)') '!!!!! Error : Rougness model number !!!!!'
 end select
 Flw(m)%RH(:, :) = Flw(m)%RH(:, :) / lRef
 select case (swi_ICM)
  case(0)
   if(swi_subgrid .eq. 2) then
    if(TurbNum .eq. 4) then
     call BoundaryBladeSurface_lRe( &
     &      m, Flw(m)%i1 + 0, Flw(m)%i3 - 0, Flw(m)%ks + 1, Flw(m)%ks + 1, &
     &      Flw(m)%js + 0, Flw(m)%js + 1, Flw(m)%js + 2 )
    else
     call BoundaryBladeSurface( &
     &      m, Flw(m)%i1 + 0, Flw(m)%i3 - 0, Flw(m)%ks + 1, Flw(m)%ks + 1, &
     &      Flw(m)%js + 0, Flw(m)%js + 1, Flw(m)%js + 2 )
    end if
   else
    if(TurbNum .eq. 4) then
     call BoundaryBladeSurface_lRe( &
     &      m, Flw(m)%is + 0, Flw(m)%ie - 0, Flw(m)%ks + 1, Flw(m)%ks + 1, &
     &      Flw(m)%js + 0, Flw(m)%js + 1, Flw(m)%js + 2 )
    else
     call BoundaryBladeSurface( &
     &      m, Flw(m)%is + 0, Flw(m)%ie - 0, Flw(m)%ks + 1, Flw(m)%ks + 1, &
     &      Flw(m)%js + 0, Flw(m)%js + 1, Flw(m)%js + 2 )
    end if
   end if
  case(1)
   if(TurbNum .eq. 4) then
    call BoundaryIcingBladeSurface_lRe(m)
   else
    call BoundaryIcingBladeSurface(m)
   end if
 end select
 ! 周期境界 --------------------------------------------------------------------------------------------
 call BoundaryPeriodic( &
 &      m, Flw(m)%is + 0, Flw(m)%ie - 0, Flw(m)%js + 0, Flw(m)%je - 0, &
 &      Flw(m)%ks + 0, Flw(m)%ks + 1 )
 call BoundaryPeriodic( &
 &      m, Flw(m)%is + 0, Flw(m)%ie - 0, Flw(m)%js + 0, Flw(m)%je - 0, &
 &      Flw(m)%ke - 1, Flw(m)%ks + 1 )
 call BoundaryPeriodic( &
 &      m, Flw(m)%is + 0, Flw(m)%ie - 0, Flw(m)%js + 0, Flw(m)%je - 0, &
 &      Flw(m)%ke - 0, Flw(m)%ks + 1 )

close(700)
close(701)


 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryCondition
!*******************************************************************************************************
!******** 流入境界										********
!******** @．密度外挿，その他固定　			A．角度・全温・全圧固定，マッハ数外挿	********
!******** B．角度・体積流量・全温固定，密度外挿  						********
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
    t   = TsExp / aRef**2
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
!******** 流出境界										********
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
&            m, is, ie, ks, ke, j0, j1, j2 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 integer, intent(in) :: is, ie, ks, ke
 integer, intent(in) :: j0, j1, j2
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k, heat
 real :: utau
 ! 処理開始 ********************************************************************************************
 if(m .eq. 1) then
  do i = Flw(m)%is,Flw(m)%i1-1
   write(701,*) 0.0
  end do
 else if(m .eq. 2 .and. swi_subgrid .eq. 2) then
  do i = Flw(m)%is,Flw(m)%i1-1
   write(700,*) 0.0
  end do
 end if

 do k = ks, ke
 do i = is, ie
  if(fSlip) then
    if( Flw(m)%qh(i,j1,k,1) <= 0.0 .or. Flw(m)%qh(i,j2,k,1) <= 0.0 ) cycle
    call WallSlipWF(m, i, j0, k, i, j1, k, i, j2, k)
   else
    if( Flw(m)%qh(i,j1,k,1) <= 0.0 ) cycle
    call WallNoslipWF(m, i, j0, k, i, j1, k ,utau)
  endif

  if( m .eq. 1 ) then
   !write(700,*)'k     i     utau '
   write(701,*) utau
  else if( m .eq. 2 ) then
   write(700,*) utau
  end if

 enddo
 enddo

 if(m .eq. 1) then
  do i = Flw(m)%i3+1,Flw(m)%ie
   write(701,*) 0.0
  end do
 else if(m .eq. 2 .and. swi_subgrid .eq. 2) then
  do i = Flw(m)%i3+1,Flw(m)%ie
   write(700,*) 0.0
  end do
 end if

 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryBladeSurface

subroutine BoundaryBladeSurface_lRe( &
&            m, is, ie, ks, ke, j0, j1, j2 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 integer, intent(in) :: is, ie, ks, ke
 integer, intent(in) :: j0, j1, j2
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k, heat
 real :: utau
 ! 処理開始 ********************************************************************************************
 if(m .eq. 1) then
  do i = Flw(m)%is,Flw(m)%i1-1
   write(701,*) 0.0
  end do
 else if(m .eq. 2 .and. swi_subgrid .eq. 2) then
  do i = Flw(m)%is,Flw(m)%i1-1
   write(700,*) 0.0
  end do
 end if

 do k = ks, ke
 do i = is, ie
  if(fSlip) then
    if( Flw(m)%qh(i,j1,k,1) <= 0.0 .or. Flw(m)%qh(i,j2,k,1) <= 0.0 ) cycle
    call WallSlipWF(m, i, j0, k, i, j1, k, i, j2, k)
   else
    if( Flw(m)%qh(i,j1,k,1) <= 0.0 ) cycle
!    call WallNoslipWF(m, i, j0, k, i, j1, k ,utau)
    call WallNoSlipLRe(m, i, j0, k, i, j1, k ,utau)
  endif

  if( m .eq. 1 ) then
   !write(700,*)'k     i     utau '
   write(701,*) utau
  else if( m .eq. 2 ) then
   write(700,*) utau
  end if

 enddo
 enddo

 if(m .eq. 1) then
  do i = Flw(m)%i3+1,Flw(m)%ie
   write(701,*) 0.0
  end do
 else if(m .eq. 2 .and. swi_subgrid .eq. 2) then
  do i = Flw(m)%i3+1,Flw(m)%ie
   write(700,*) 0.0
  end do
 end if

 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryBladeSurface_lRe

subroutine BoundaryIcingBladeSurface(m)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: is,ie
 real   , pointer :: qq(:, :)
 integer :: i, j, n, nmax
 integer :: ii, jj, l, k
 real :: ut
 real,dimension(Flw(2)%is:Flw(2)%ie) :: utau
 real,dimension(Flw(1)%is:Flw(1)%ie) :: utau_main
 ! 処理開始 ********************************************************************************************
 ! ic =  0; 流体セル
 ! ic =  1; 着氷セル (完全に氷の中)
 ! ic =  2; i-1が流体セル，i+1が流体セルであるj方向に最大の氷層セル
 ! ic =  3; i-1が氷層セル，i+1が氷層セルであるj方向に最大の氷層セル
 ! ic =  4; i-1が流体セル，i+1が氷層セルであるj方向に最大の氷層セル
 ! ic =  5; i-1が氷層セル，i+1が流体セルであるj方向に最大の氷層セル
 ! ic =  6; i-1が氷層セル，i+1が氷層セルであるj方向に最小の氷層セル
 ! ic =  7; i-1が氷層セル，i+1が着氷セルであるj方向に最小の氷層セル
 ! ic =  8; i-1が着氷セル，i+1が氷層セルであるj方向に最小の氷層セル
 ! ic =  9; i-1が流体セル，i+1が流体セルである氷層セル
 ! ic = 10; i-1が流体セル，i+1が着氷セルである氷層セル
 ! ic = 11; i-1が着氷セル，i+1が流体セルである氷層セル
 ! ic = 12; 壁面セル (着氷なし)

 utau(:) = 0.0
 utau_main(:) = 0.0

! select case(m)
! case(ms)
!  is = Flw(m)%i1
!  ie = Flw(m)%i3
! case(me)
!  is = Flw(m)%is + 1
!  ie = Flw(m)%ie - 1
! end select

 if(m .eq. ms .or. swi_subgrid .eq. 2) then
  is = Flw(m)%i1
  ie = Flw(m)%i3
 else
  is = Flw(m)%is + 1
  ie = Flw(m)%ie - 1
 end if


 j = Flw(m)%js
 do k = Flw(m)%ks,Flw(m)%ke
  do i = is,ie

   if(Flw(m)%qh(i,j,k,1) .ne. Flw(m)%qh(i,j,k,1)) then
    write(*,*) 'rho = NaN @BC',i,j,k,nCount
    stop
   end if

   if(Flw(m)%ic(i-1,j) .eq. 12 .or. Flw(m)%ic(i,j) .eq. 12) then
    call WallNoSlipWF(m, i, j, k, i, j+1, k, ut)
    ! utau出力
    if( m .eq. 2 .and. k .eq. Flw(m)%ks ) utau(i) = ut
    if( m .eq. 1 .and. k .eq. Flw(m)%ks ) utau_main(i) = ut
   end if
  end do
 end do

 do k = Flw(m)%ks,Flw(m)%ke
  do i = is,ie
   do j = Flw(m)%js+1,amax0(Flw(m)%ji(i),Flw(m)%ji(i-1),Flw(m)%ji(i+1))+1
    if(Flw(m)%ic(i-1,j-1) .eq. 3 .or. Flw(m)%ic(i,j-1) .eq. 3) then
     call WallNoSlipWF(m, i, j, k, i, j+1, k, ut)
     ! utau出力
     if( m .eq. 2 .and. k .eq. Flw(m)%ks ) utau(i) = ut
     if( m .eq. 1 .and. k .eq. Flw(m)%ks ) utau_main(i) = ut
    else if(Flw(m)%ic(i,j-1) .eq. 2 .or. Flw(m)%ic(i,j-1) .eq. 4) then
     nmax = 2
     allocate( qq(1:nmax, ls:le) )
     n = 1
     call WallNoSlipWF(m, i, j, k, i-1, j, k, ut)
     do l = ls, le
      qq(n,l) = Flw(m)%qh(i,j+1,k,l) * Flw(m)%jac(i,j+1,k)
     end do
     n = 2
     call WallNoSlipWF(m, i, j, k, i, j+1, k, ut)
     do l = ls, le
      qq(n,l) = Flw(m)%qh(i,j+1,k,l) * Flw(m)%jac(i,j+1,k)
     end do
     ! utau出力
     if( m .eq. 2 .and. k .eq. Flw(m)%ks ) utau(i) = ut
     if( m .eq. 1 .and. k .eq. Flw(m)%ks ) utau_main(i) = ut
     ! 平均値外挿
     do l = ls, le
      Flw(m)%qh(i,j+1,k,l) = sum(qq(:,l)) * 0.5 / Flw(m)%jac(i,j+1,k)
     end do
     deallocate(qq)
    else if(Flw(m)%ic(i-1,j-1) .eq. 2 .or. Flw(m)%ic(i-1,j-1) .eq. 5) then
     nmax = 2
     allocate( qq(1:nmax, ls:le) )
     n = 1
     call WallNoSlipWF(m, i, j, k, i+1, j, k, ut)
     do l = ls, le
      qq(n,l) = Flw(m)%qh(i,j+1,k,l) * Flw(m)%jac(i,j+1,k)
     end do
     n = 2
     call WallNoSlipWF(m, i, j, k, i, j+1, k, ut)
     do l = ls, le
      qq(n,l) = Flw(m)%qh(i,j+1,k,l) * Flw(m)%jac(i,j+1,k)
     end do
     ! utau出力
     if( m .eq. 2 .and. k .eq. Flw(m)%ks ) utau(i) = ut
     if( m .eq. 1 .and. k .eq. Flw(m)%ks ) utau_main(i) = ut
     ! 平均値外挿
     do l = ls, le
      Flw(m)%qh(i,j+1,k,l) = sum(qq(:,l)) * 0.5 / Flw(m)%jac(i,j+1,k)
     end do
     deallocate(qq)
    else if(Flw(m)%ic(i,j-1) .eq. 9 .or. Flw(m)%ic(i,j-1) .eq. 10 .or. &
         &  Flw(m)%ic(i,j) .eq. 10) then
     call WallNoSlipWF(m, i, j, k, i-1, j, k, ut)
    else if(Flw(m)%ic(i-1,j-1) .eq. 9 .or. Flw(m)%ic(i-1,j-1) .eq. 11 .or. &
         &  Flw(m)%ic(i-1,j) .eq. 11) then
     call WallNoSlipWF(m, i, j, k, i+1, j, k, ut)
    end if
   end do
  end do
 end do

 if(m .eq. 2) then
  write(700,*) utau(Flw(m)%is+1)
  do i = Flw(m)%is+1,Flw(m)%ie-1
   write(700,*) utau(i)
  end do
  write(700,*) utau(Flw(m)%ie-1)
 end if

 if(m .eq. 1) then
  write(701,*) utau_main(Flw(m)%is+1)
  do i = Flw(m)%is+1,Flw(m)%ie-1
   write(701,*) utau_main(i)
  end do
  write(701,*) utau_main(Flw(m)%ie-1)
 end if

 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryIcingBladeSurface

subroutine BoundaryIcingBladeSurface_lRe(m)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: is,ie
 real   , pointer :: qq(:, :)
 integer :: i, j, n, nmax
 integer :: ii, jj, l, k
 real :: ut
 real,dimension(Flw(2)%is:Flw(2)%ie) :: utau
 real,dimension(Flw(1)%is:Flw(1)%ie) :: utau_main
 ! 処理開始 ********************************************************************************************
 ! ic =  0; 流体セル
 ! ic =  1; 着氷セル (完全に氷の中)
 ! ic =  2; i-1が流体セル，i+1が流体セルであるj方向に最大の氷層セル
 ! ic =  3; i-1が氷層セル，i+1が氷層セルであるj方向に最大の氷層セル
 ! ic =  4; i-1が流体セル，i+1が氷層セルであるj方向に最大の氷層セル
 ! ic =  5; i-1が氷層セル，i+1が流体セルであるj方向に最大の氷層セル
 ! ic =  6; i-1が氷層セル，i+1が氷層セルであるj方向に最小の氷層セル
 ! ic =  7; i-1が氷層セル，i+1が着氷セルであるj方向に最小の氷層セル
 ! ic =  8; i-1が着氷セル，i+1が氷層セルであるj方向に最小の氷層セル
 ! ic =  9; i-1が流体セル，i+1が流体セルである氷層セル
 ! ic = 10; i-1が流体セル，i+1が着氷セルである氷層セル
 ! ic = 11; i-1が着氷セル，i+1が流体セルである氷層セル
 ! ic = 12; 壁面セル (着氷なし)

 utau(:) = 0.0
 utau_main(:) = 0.0

 select case(m)
 case(ms)
  is = Flw(m)%i1
  ie = Flw(m)%i3
 case(me)
  is = Flw(m)%is + 1
  ie = Flw(m)%ie - 1
 end select

 j = Flw(m)%js
 do k = Flw(m)%ks,Flw(m)%ke
  do i = is,ie

   if(Flw(m)%qh(i,j,k,1) .ne. Flw(m)%qh(i,j,k,1)) then
    write(*,*) 'rho = NaN @BC',i,j,k,nCount
    stop
   end if

   if(Flw(m)%ic(i-1,j) .eq. 12 .or. Flw(m)%ic(i,j) .eq. 12) then
    call WallNoSlipLRe(m, i, j, k, i, j+1, k, ut)
    ! utau出力
    if( m .eq. 2 .and. k .eq. Flw(m)%ks ) utau(i) = ut
    if( m .eq. 1 .and. k .eq. Flw(m)%ks ) utau_main(i) = ut
   end if
  end do
 end do

 do k = Flw(m)%ks,Flw(m)%ke
  do i = is,ie
   do j = Flw(m)%js+1,amax0(Flw(m)%ji(i),Flw(m)%ji(i-1),Flw(m)%ji(i+1))+1
    if(Flw(m)%ic(i-1,j-1) .eq. 3 .or. Flw(m)%ic(i,j-1) .eq. 3) then
     call WallNoSlipLRe(m, i, j, k, i, j+1, k, ut)
     ! utau出力
     if( m .eq. 2 .and. k .eq. Flw(m)%ks ) utau(i) = ut
     if( m .eq. 1 .and. k .eq. Flw(m)%ks ) utau_main(i) = ut
    else if(Flw(m)%ic(i,j-1) .eq. 2 .or. Flw(m)%ic(i,j-1) .eq. 4) then
     nmax = 2
     allocate( qq(1:nmax, ls:le) )
     n = 1
     call WallNoSlipLRe(m, i, j, k, i-1, j, k, ut)
     do l = ls, le
      qq(n,l) = Flw(m)%qh(i,j+1,k,l) * Flw(m)%jac(i,j+1,k)
     end do
     n = 2
     call WallNoSlipLRe(m, i, j, k, i, j+1, k, ut)
     do l = ls, le
      qq(n,l) = Flw(m)%qh(i,j+1,k,l) * Flw(m)%jac(i,j+1,k)
     end do
     ! utau出力
     if( m .eq. 2 .and. k .eq. Flw(m)%ks ) utau(i) = ut
     if( m .eq. 1 .and. k .eq. Flw(m)%ks ) utau_main(i) = ut
     ! 平均値外挿
     do l = ls, le
      Flw(m)%qh(i,j+1,k,l) = sum(qq(:,l)) * 0.5 / Flw(m)%jac(i,j+1,k)
     end do
     deallocate(qq)
    else if(Flw(m)%ic(i-1,j-1) .eq. 2 .or. Flw(m)%ic(i-1,j-1) .eq. 5) then
     nmax = 2
     allocate( qq(1:nmax, ls:le) )
     n = 1
     call WallNoSlipLRe(m, i, j, k, i+1, j, k, ut)
     do l = ls, le
      qq(n,l) = Flw(m)%qh(i,j+1,k,l) * Flw(m)%jac(i,j+1,k)
     end do
     n = 2
     call WallNoSlipLRe(m, i, j, k, i, j+1, k, ut)
     do l = ls, le
      qq(n,l) = Flw(m)%qh(i,j+1,k,l) * Flw(m)%jac(i,j+1,k)
     end do
     ! utau出力
     if( m .eq. 2 .and. k .eq. Flw(m)%ks ) utau(i) = ut
     if( m .eq. 1 .and. k .eq. Flw(m)%ks ) utau_main(i) = ut
     ! 平均値外挿
     do l = ls, le
      Flw(m)%qh(i,j+1,k,l) = sum(qq(:,l)) * 0.5 / Flw(m)%jac(i,j+1,k)
     end do
     deallocate(qq)
    else if(Flw(m)%ic(i,j-1) .eq. 9 .or. Flw(m)%ic(i,j-1) .eq. 10 .or. &
         &  Flw(m)%ic(i,j) .eq. 10) then
     call WallNoSlipLRe(m, i, j, k, i-1, j, k, ut)
    else if(Flw(m)%ic(i-1,j-1) .eq. 9 .or. Flw(m)%ic(i-1,j-1) .eq. 11 .or. &
         &  Flw(m)%ic(i-1,j) .eq. 11) then
     call WallNoSlipLRe(m, i, j, k, i+1, j, k, ut)
    end if
   end do
  end do
 end do

 if(m .eq. 2) then
  write(700,*) utau(Flw(m)%is+1)
  do i = Flw(m)%is+1,Flw(m)%ie-1
   write(700,*) utau(i)
  end do
  write(700,*) utau(Flw(m)%ie-1)
 end if

 if(m .eq. 1) then
  write(701,*) utau_main(Flw(m)%is+1)
  do i = Flw(m)%is+1,Flw(m)%ie-1
   write(701,*) utau_main(i)
  end do
  write(701,*) utau_main(Flw(m)%ie-1)
 end if

 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryIcingBladeSurface_lRe

!*******************************************************************************************************
!******** 壁境界 (滑りなし, 温度固定, 壁関数)							********
!*******************************************************************************************************
subroutine WallNoslipWF( &
&            m, i0, j0, k0, i1, j1, k1,utau)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 integer, intent(in) :: i0, j0, k0
 integer, intent(in) :: i1, j1, k1
 real,intent(inout)  :: utau
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real :: dx, dy, dz
 real :: uc, vc, wc
 real :: yp, up, nup, kp, epsp, dudy1, dudy2, pp, Tp
 ! 処理開始 ********************************************************************************************
 ! 壁関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁近傍の分布 ----------------------------------------------------------------------------------------
 dx = Flw(m)%x(i1,j1,k1) - Flw(m)%x(i0,j0,k0)
 dy = Flw(m)%y(i1,j1,k1) - Flw(m)%y(i0,j0,k0)
 dz = Flw(m)%z(i1,j1,k1) - Flw(m)%z(i0,j0,k0)
 uc = Flw(m)%qh(i1,j1,k1,2) / Flw(m)%qh(i1,j1,k1,1)
 vc = Flw(m)%qh(i1,j1,k1,3) / Flw(m)%qh(i1,j1,k1,1)
 wc = Flw(m)%qh(i1,j1,k1,4) / Flw(m)%qh(i1,j1,k1,1)
 ! 壁関数 ----------------------------------------------------------------------------------------------
 yp   = sqrt(dx**2 + dy**2 + dz**2)
 up   = sqrt(uc**2 + vc**2 + wc**2)
 up   = max(zero, up)
 nup  = Flw(m)%mu(i1,j1,k1) / (Flw(m)%qh(i1,j1,k1,1) * Flw(m)%jac(i1,j1,k1))
 kp   = Flw(m)%kin(i1,j1,k1)
 epsp = Flw(m)%eps(i1,j1,k1)
 nup  = max(zero, nup)
 kp   = max(zero, kp)
 epsp = max(zero, epsp)
 select case( Flw(m)%fRough(i0,k0) )
  case(0)
   call WallFunctionKEM2S( &
   &      yp, up, nup, utau, dudy1, dudy2, kp, epsp )
  case(1)
   call WallFunctionKEM2R( &
   &      Flw(m)%RH(i0,k0), yp, up, nup, utau, dudy1, dudy2, kp, epsp )
  case default
   write(*, '(a)') '!!!!! Error : fRough number !!!!!'
 end select
 ! 速度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁面 ------------------------------------------------------------------------------------------------
 uc = 0.0
 vc = 0.0
 wc = 0.0
 ! 流束関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁面から一点目 --------------------------------------------------------------------------------------
 ! kp, epsp
 Flw(m)%qh(i1,j1,k1,5) = Flw(m)%qh(i1,j1,k1,5) - Flw(m)%qh(i1,j1,k1,6)
 Flw(m)%qh(i1,j1,k1,6) = Flw(m)%qh(i1,j1,k1,1) * kp
 Flw(m)%qh(i1,j1,k1,7) = Flw(m)%qh(i1,j1,k1,1) * epsp
 Flw(m)%qh(i1,j1,k1,5) = Flw(m)%qh(i1,j1,k1,5) + Flw(m)%qh(i1,j1,k1,6)
 ! 壁面 ------------------------------------------------------------------------------------------------
 ! 温度固定,圧力・k・εは外挿
 pp = ( Flw(m)%qh(i1,j1,k1,5) - Flw(m)%qh(i1,j1,k1,6) &
 &       - 0.5 * ( Flw(m)%qh(i1,j1,k1,2)**2 &
 &               + Flw(m)%qh(i1,j1,k1,3)**2 &
 &               + Flw(m)%qh(i1,j1,k1,4)**2 &
 &                ) / Flw(m)%qh(i1,j1,k1,1) ) * (gamma - 1.0) * Flw(m)%jac(i1,j1,k1)
 Tp = ( Flw(m)%qh0(i0,j0,k0,5) - Flw(m)%qh0(i0,j0,k0,6) &
 &       - 0.5 * ( Flw(m)%qh0(i0,j0,k0,2)**2 &
 &               + Flw(m)%qh0(i0,j0,k0,3)**2 &
 &               + Flw(m)%qh0(i0,j0,k0,4)**2 &
 &                ) / Flw(m)%qh0(i0,j0,k0,1) ) * (gamma - 1.0) / ( Flw(m)%qh0(i0,j0,k0,1) * Rg )
 Flw(m)%qh(i0,j0,k0,1) = pp / ( Rg * Tp ) / Flw(m)%jac(i0,j0,k0)
 Flw(m)%qh(i0,j0,k0,2) = Flw(m)%qh(i0,j0,k0,1) * uc
 Flw(m)%qh(i0,j0,k0,3) = Flw(m)%qh(i0,j0,k0,1) * vc
 Flw(m)%qh(i0,j0,k0,4) = Flw(m)%qh(i0,j0,k0,1) * wc
 Flw(m)%qh(i0,j0,k0,6) = Flw(m)%qh(i1,j1,k1,6) / Flw(m)%qh(i1,j1,k1,1) * Flw(m)%qh(i0,j0,k0,1)
 Flw(m)%qh(i0,j0,k0,7) = Flw(m)%qh(i1,j1,k1,7) / Flw(m)%qh(i1,j1,k1,1) * Flw(m)%qh(i0,j0,k0,1)
 Flw(m)%qh(i0,j0,k0,5) = Rg * Tp / ( gamma - 1.0 ) * Flw(m)%qh(i0,j0,k0,1) &
 &                      + 0.5 * ( Flw(m)%qh(i0,j0,k0,2)**2 &
 &                              + Flw(m)%qh(i0,j0,k0,3)**2 &
 &                              + Flw(m)%qh(i0,j0,k0,4)**2 ) / Flw(m)%qh(i0,j0,k0,1) &
 &                      + Flw(m)%qh(i0,j0,k0,6)
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
 real :: yp, up, nup, kp, epsp, utau, dudy1, dudy2, pp, Tp
 ! 処理開始 ********************************************************************************************
 ! 壁関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁近傍の分布 ----------------------------------------------------------------------------------------
 dx2 = Flw(m)%x(i2,j2,k2) - Flw(m)%x(i0,j0,k0)
 dy2 = Flw(m)%y(i2,j2,k2) - Flw(m)%y(i0,j0,k0)
 dz2 = Flw(m)%z(i2,j2,k2) - Flw(m)%z(i0,j0,k0)
 uc2 = Flw(m)%qh(i2,j2,k2,2) / Flw(m)%qh(i2,j2,k2,1)
 vc2 = Flw(m)%qh(i2,j2,k2,3) / Flw(m)%qh(i2,j2,k2,1)
 wc2 = Flw(m)%qh(i2,j2,k2,4) / Flw(m)%qh(i2,j2,k2,1)
 ! 壁関数 ----------------------------------------------------------------------------------------------
 yp   = sqrt(dx2**2 + dy2**2 + dz2**2)
 up   = sqrt(uc2**2 + vc2**2 + wc2**2)
 nup  = Flw(m)%mu(i2,j2,k2) / (Flw(m)%qh(i2,j2,k2,1) * Flw(m)%jac(i2,j2,k2))
 kp   = Flw(m)%kin(i2,j2,k2)
 epsp = Flw(m)%eps(i2,j2,k2)
 up   = max(zero, up)
 nup  = max(zero, nup)
 kp   = max(zero, kp)
 epsp = max(zero, epsp)
 select case( Flw(m)%fRough(i0,k0) )
  case(0)
   call WallFunctionKEM2S( &
   &      yp, up, nup, utau, dudy1, dudy2, kp, epsp )
  case(1)
   call WallFunctionKEM2R( &
   &      Flw(m)%RH(i0,k0), yp, up, nup, utau, dudy1, dudy2, kp, epsp )
  case default
   write(*, '(a)') '!!!!! Error : fRough number !!!!!'
 end select
 ! 速度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁面から一点目 --------------------------------------------------------------------------------------
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
 ! 壁面 ------------------------------------------------------------------------------------------------
 uc0 = 0.0
 vc0 = 0.0
 wc0 = 0.0
 ! 流束関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁面から二点目 --------------------------------------------------------------------------------------
 ! kp, epsp 外挿
 Flw(m)%qh(i2,j2,k2,5) = Flw(m)%qh(i2,j2,k2,5) - Flw(m)%qh(i2,j2,k2,6)
 Flw(m)%qh(i2,j2,k2,6) = Flw(m)%qh(i2,j2,k2,1) * kp
 Flw(m)%qh(i2,j2,k2,7) = Flw(m)%qh(i2,j2,k2,1) * epsp
 Flw(m)%qh(i2,j2,k2,5) = Flw(m)%qh(i2,j2,k2,5) + Flw(m)%qh(i2,j2,k2,6)
 ! 壁面から一点目 --------------------------------------------------------------------------------------
 ! up, vp, wp, kp, epsp 外挿
 Flw(m)%qh(i1,j1,k1,5) = Flw(m)%qh(i1,j1,k1,5) &
 &                    - 0.5 * ( Flw(m)%qh(i1,j1,k1,2)**2 &
 &                            + Flw(m)%qh(i1,j1,k1,3)**2 &
 &                            + Flw(m)%qh(i1,j1,k1,4)**2 ) / Flw(m)%qh(i1,j1,k1,1) &
 &                    - Flw(m)%qh(i1,j1,k1,6)
 Flw(m)%qh(i1,j1,k1,2) = uc1  * Flw(m)%qh(i1,j1,k1,1)
 Flw(m)%qh(i1,j1,k1,3) = vc1  * Flw(m)%qh(i1,j1,k1,1)
 Flw(m)%qh(i1,j1,k1,4) = wc1  * Flw(m)%qh(i1,j1,k1,1)
 Flw(m)%qh(i1,j1,k1,6) = kp   * Flw(m)%qh(i1,j1,k1,1)
 Flw(m)%qh(i1,j1,k1,7) = epsp * Flw(m)%qh(i1,j1,k1,1)
 Flw(m)%qh(i1,j1,k1,5) = Flw(m)%qh(i1,j1,k1,5) &
 &                    + 0.5 * ( Flw(m)%qh(i1,j1,k1,2)**2 &
 &                            + Flw(m)%qh(i1,j1,k1,3)**2 &
 &                            + Flw(m)%qh(i1,j1,k1,4)**2 ) / Flw(m)%qh(i1,j1,k1,1) &
 &                    + Flw(m)%qh(i1,j1,k1,6)
 ! 壁面 ------------------------------------------------------------------------------------------------
 ! 全て外挿
 ! 温度固定,その他の物理量は外挿
 pp = ( Flw(m)%qh(i1,j1,k1,5) - Flw(m)%qh(i1,j1,k1,6) &
 &       - 0.5 * ( Flw(m)%qh(i1,j1,k1,2)**2 &
 &               + Flw(m)%qh(i1,j1,k1,3)**2 &
 &               + Flw(m)%qh(i1,j1,k1,4)**2 &
 &                ) / Flw(m)%qh(i1,j1,k1,1) ) * (gamma - 1.0) * Flw(m)%jac(i1,j1,k1)
 Tp = ( Flw(m)%qh0(i0,j0,k0,5) - Flw(m)%qh0(i0,j0,k0,6) &
 &       - 0.5 * ( Flw(m)%qh0(i0,j0,k0,2)**2 &
 &               + Flw(m)%qh0(i0,j0,k0,3)**2 &
 &               + Flw(m)%qh0(i0,j0,k0,4)**2 &
 &                ) / Flw(m)%qh0(i0,j0,k0,1) ) * (gamma - 1.0) / ( Flw(m)%qh0(i0,j0,k0,1) * Rg )
 Flw(m)%qh(i0,j0,k0,1) = pp / ( Rg * Tp ) / Flw(m)%jac(i0,j0,k0)
 Flw(m)%qh(i0,j0,k0,2) = Flw(m)%qh(i0,j0,k0,1) * uc0
 Flw(m)%qh(i0,j0,k0,3) = Flw(m)%qh(i0,j0,k0,1) * vc0
 Flw(m)%qh(i0,j0,k0,4) = Flw(m)%qh(i0,j0,k0,1) * wc0
 Flw(m)%qh(i0,j0,k0,6) = Flw(m)%qh(i1,j1,k1,6) / Flw(m)%qh(i1,j1,k1,1) * Flw(m)%qh(i0,j0,k0,1)
 Flw(m)%qh(i0,j0,k0,7) = Flw(m)%qh(i1,j1,k1,7) / Flw(m)%qh(i1,j1,k1,1) * Flw(m)%qh(i0,j0,k0,1)
 Flw(m)%qh(i0,j0,k0,5) = Rg * Tp / ( gamma - 1.0 ) * Flw(m)%qh(i0,j0,k0,1) &
 &                      + 0.5 * ( Flw(m)%qh(i0,j0,k0,2)**2 &
 &                              + Flw(m)%qh(i0,j0,k0,3)**2 &
 &                              + Flw(m)%qh(i0,j0,k0,4)**2 ) / Flw(m)%qh(i0,j0,k0,1) &
 &                      + Flw(m)%qh(i0,j0,k0,6)
 ! 処理終了 ********************************************************************************************
 return
end subroutine WallslipWF

!*******************************************************************************************************
!******** 壁境界 (滑りなし, 温度固定, 低レイノルズ数型モデル用)					********
!*******************************************************************************************************
subroutine WallNoslipLRe( &
&            m, i0, j0, k0, i1, j1, k1,utau)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m
 integer, intent(in) :: i0, j0, k0
 integer, intent(in) :: i1, j1, k1
 real,intent(inout)  :: utau
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real :: dx, dy, dz
 real :: uc, vc, wc
 real :: yp, pp, Tp
 double precision :: epsp
 real :: cross, uwall, uwx, uwy, uwz
 ! 処理開始 ********************************************************************************************
 ! 壁近傍の分布 ----------------------------------------------------------------------------------------
 dx = Flw(m)%x(i1,j1,k1) - Flw(m)%x(i0,j0,k0)
 dy = Flw(m)%y(i1,j1,k1) - Flw(m)%y(i0,j0,k0)
 dz = Flw(m)%z(i1,j1,k1) - Flw(m)%z(i0,j0,k0)
 yp   = sqrt(dx**2 + dy**2 + dz**2)
 ! 摩擦速度 --------------------------------------------------------------------------------------------
 cross = Flw(m)%u(i1,j1,k1) * dx + Flw(m)%v(i1,j1,k1) * dy + Flw(m)%w(i1,j1,k1) * dz
 uwx = Flw(m)%u(i1,j1,k1) - cross * dx / yp**2.0
 uwy = Flw(m)%v(i1,j1,k1) - cross * dy / yp**2.0
 uwz = Flw(m)%w(i1,j1,k1) - cross * dz / yp**2.0
 uwall = sqrt(uwx**2.0 + uwy**2.0 + uwz**2.0)
 utau = sqrt((Flw(m)%mu(i1,j1,k1) * uwall / yp) / (Flw(m)%qh(i1,j1,k1,1) * Flw(m)%jac(i1,j1,k1)))
 ! 速度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! 壁面 ------------------------------------------------------------------------------------------------
 uc = 0.0
 vc = 0.0
 wc = 0.0
 ! 流束関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! kp, epsp
 epsp = 2.0 * dble(Flw(m)%mu(i1,j1,k1) * Flw(m)%kin(i1,j1,k1)) &
      & / (dble(Flw(m)%qh(i1,j1,k1,1) * Flw(m)%jac(i1,j1,k1)) * dble(yp)**2.0)
 if(epsp .ne. epsp) then
  write(*,*) 'eps_p is diverged'
  write(*,*) Flw(m)%mu(i1,j1,k1),Flw(m)%qh(i1,j1,k1,6)*Flw(m)%jac(i1,j1,k1), &
  &          Flw(m)%qh(i1,j1,k1,1)*Flw(m)%jac(i1,j1,k1),yp
  stop
 end if
 Flw(m)%qh(i1,j1,k1,7) = Flw(m)%qh(i1,j1,k1,1) * real(epsp) / Flw(m)%jac(i1,j1,k1)
 ! 壁面 ------------------------------------------------------------------------------------------------
 ! 温度固定,圧力・k・εは外挿
 pp = ( Flw(m)%qh(i1,j1,k1,5) - Flw(m)%qh(i1,j1,k1,6) &
 &       - 0.5 * ( Flw(m)%qh(i1,j1,k1,2)**2 &
 &               + Flw(m)%qh(i1,j1,k1,3)**2 &
 &               + Flw(m)%qh(i1,j1,k1,4)**2 &
 &                ) / Flw(m)%qh(i1,j1,k1,1) ) * (gamma - 1.0) * Flw(m)%jac(i1,j1,k1)
 Tp = ( Flw(m)%qh0(i0,j0,k0,5) - Flw(m)%qh0(i0,j0,k0,6) &
 &       - 0.5 * ( Flw(m)%qh0(i0,j0,k0,2)**2 &
 &               + Flw(m)%qh0(i0,j0,k0,3)**2 &
 &               + Flw(m)%qh0(i0,j0,k0,4)**2 &
 &                ) / Flw(m)%qh0(i0,j0,k0,1) ) * (gamma - 1.0) / ( Flw(m)%qh0(i0,j0,k0,1) * Rg )
 Flw(m)%qh(i0,j0,k0,1) = pp / ( Rg * Tp ) / Flw(m)%jac(i0,j0,k0)
 Flw(m)%qh(i0,j0,k0,2) = Flw(m)%qh(i0,j0,k0,1) * uc
 Flw(m)%qh(i0,j0,k0,3) = Flw(m)%qh(i0,j0,k0,1) * vc
 Flw(m)%qh(i0,j0,k0,4) = Flw(m)%qh(i0,j0,k0,1) * wc
 Flw(m)%qh(i0,j0,k0,6) = 0.0
 Flw(m)%qh(i0,j0,k0,7) = Flw(m)%qh(i1,j1,k1,7) / Flw(m)%qh(i1,j1,k1,1) * Flw(m)%qh(i0,j0,k0,1)
 Flw(m)%qh(i0,j0,k0,5) = Rg * Tp / ( gamma - 1.0 ) * Flw(m)%qh(i0,j0,k0,1) &
 &                      + 0.5 * ( Flw(m)%qh(i0,j0,k0,2)**2 &
 &                              + Flw(m)%qh(i0,j0,k0,3)**2 &
 &                              + Flw(m)%qh(i0,j0,k0,4)**2 ) / Flw(m)%qh(i0,j0,k0,1) &
 &                      + Flw(m)%qh(i0,j0,k0,6)
 ! 処理終了 ********************************************************************************************
 return
end subroutine WallNoSlipLRe

!*******************************************************************************************************
!******** 氷の中の点の処理 (氷層の平均値を外挿)		 	                                *******
!*******************************************************************************************************
subroutine BoundaryIceIn( &
&            m1, m2, is, ie, js, je, ks, ke, i1, i2, j0, k0 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: m1, m2
 integer, intent(in) :: is, ie, js, je, ks, ke
 integer, intent(in) :: i1, i2, j0, k0
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: qAve(:)
 integer :: i, j, k, l
 ! 処理開始 ********************************************************************************************
 ! メモリ確保 ------------------------------------------------------------------------------------------
 allocate( qAve(ls:le) )
 ! 氷層の平均値 ----------------------------------------------------------------------------------------
 do l = ls, le
  qAve(l) = sum( Flw(m2)%qh(i1:i2,j0,k0,l) * Flw(m2)%jac(i1:i2,j0,k0) ) &
  &       / real(i2 - i1 + 1)
 enddo
 ! 平均値を外挿 ----------------------------------------------------------------------------------------
 do k = ks, ke
 do j = js, je
 do i = is, ie
  if( IceIn(i,j,k) == 0 ) cycle
  do l = ls, le
    Flw(m1)%qh(i,j,k,l) = qAve(l) / Flw(m1)%jac(i,j,k)
  enddo
 enddo
 enddo
 enddo
 ! メモリ解放 ------------------------------------------------------------------------------------------
 deallocate(qAve)
 ! 処理終了 ********************************************************************************************
 return
end subroutine BoundaryIceIn
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
!******** 残差計算 										********
!*******************************************************************************************************
subroutine CalResidual( CalculateNumber, Convergence, Order0 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Integer, pointer       :: OrderNum(:)
 real   , pointer       :: res0(:,:,:,:), L2_Norm(:)
 integer                :: l, m, ir, jr, kr, NodeNumber
 Integer, Intent(in)    :: CalculateNumber				!計算回数
 Real, Intent(inout)    :: Order0(ms:me,ls:le)				!各成分のL2残差
 Real                   :: Res_Local, eps, Order_L2Norm, EpsValue
 Logical, Intent(out)   :: Convergence 
 ! 処理開始 ********************************************************************************************
 Convergence = .true.
 eps = 2.0
 do m = ms, me
 ! メモリ確保 -----------------------------------------------------------------------------------------
  allocate( res0(Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le), &
  &         L2_Norm(ls: le), OrderNum(ls: le) )
 ! 残差計算 -------------------------------------------------------------------------------------------
  Select Case (m) !realの有効桁数は約7桁(10進法)
  Case(1) !Main
  call CalResidual3D( &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &     ls, le, Flw(m)%dt, q0_main, Flw(m)%qh, &
  &     res0 )
  Case(2) !Sub
  call CalResidual3D( &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &     ls, le, Flw(m)%dt, q0_sub, Flw(m)%qh, &
  &     res0 )
  Case Default
   Write(*,'(a)') '//  ERROR Grid Number //'
  End Select
 ! L2ノルムの算出 -----------------------------------------------------------------------------------------
  L2_Norm(:) = 0.0
  Do l = ls, le
   NodeNumber = 0
   Do kr = Flw(m)%ks + 1, Flw(m)%ke - 1
    Do jr = Flw(m)%js + 1, Flw(m)%je - 1
     Do ir = Flw(m)%is + 1, Flw(m)%ie - 1
       Res_Local = res0(ir,jr,kr,l) * Flw(m)%qh(ir,jr,kr,1)**2 * Flw(m)%dt(ir,jr,kr)
       L2_Norm(l)  = L2_Norm(l)  + Res_Local * Flw(m)%jac(ir,jr,kr)**2
       NodeNumber = NodeNumber + 1
     End Do
    End Do
   End Do
   L2_Norm(l)  = SQRT( L2_Norm(l) / Real(NodeNumber) )
   If( L2_Norm(l) == 0.0 ) L2_Norm(l) = 1.0e-20 !桁落ちした場合
   EpsValue = MAX( L2_Norm(l), Order0(m,l) ) / MIN( L2_Norm(l), Order0(m,l) )
   ! 流れ場の収束判定
   if ( Convergence .and. ( EpsValue < eps ) ) then
     Convergence = .true.
    else
     Convergence = .false.
   end if
   Order0(m,l) = L2_Norm(l)
  End Do
  ! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Select Case (m) !realの有効桁数は約7桁(10進法)
  Case(1) !Main
   write(20, '(i8, 7e16.7e3)') nCount, L2_Norm(ls:le)
  Case(2) !Sub
   write(21, '(i8, 7e16.7e3)') nCount, L2_Norm(ls:le)
  Case Default
   Write(*,'(a)') '//  ERROR Grid Number //'
  End Select
  ! メモリ解放 -----------------------------------------------------------------------------------------
  deallocate(res0, L2_Norm, OrderNum)
 enddo
 q0_main(:,:,:,:) = Flw(ms)%qh(:,:,:,:); q0_sub(:,:,:,:) = Flw(me)%qh(:,:,:,:)
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalResidual
!*******************************************************************************************************
!******** 熱伝導計算の時間刻み									********
!*******************************************************************************************************
Subroutine Calculation_HeatConductionTimeStep(dt_Conduction)
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer       :: i_h, j_h, k_h, m_h
 real, pointer :: d_l(:)						! 要素の中心点間の長さ
 Real, Intent(inout) :: dt_Conduction					!熱伝導計算の最小時間刻み
 Real          :: d_lmin
 Real          :: k_wing, c_wing, rho_wing				! 翼の物性値(熱伝導率,比熱,密度)
 ! 処理開始 ********************************************************************************************
 m_h = me  !! SubGrid
 j_h = Flw(m_h)%js
 !!++++++++++++++++++++++++++++++++++++++++++
 k_h = Flw(m_h)%ks + 1 !!!k方向も後で拡張する
 !!++++++++++++++++++++++++++++++++++++++++++
 ! 翼の物性値(アルミ)
 k_wing = 256.05; c_wing = 913.0; rho_wing = 2700.0 !![W/(m・k)] [J/(kg・K)] [kg/m^3]!!
 ! メモリ確保 ------------------------------------------------------------------------------------------
 allocate( d_l( Flw(m_h)%is: Flw(m_h)%ie ) )
 Do i_h = Flw(m_h)%is + 1, Flw(m_h)%ie
  d_l(i_h)    = SQRT( ( Flw(m_h)%x(i_h ,j_h,k_h) - Flw(m_h)%x(i_h - 1,j_h,k_h) )**2 &
  &                 + ( Flw(m_h)%y(i_h ,j_h,k_h) - Flw(m_h)%y(i_h - 1,j_h,k_h) )**2 ) * LRef
 End Do
 !! 要素点間の最小幅
 d_lmin = minval( d_l(Flw(m_h)%is + 1 : Flw(m_h)%ie) )
 !! TimeStepを算出
 dt_Conduction = 0.3 * c_wing * rho_wing / k_wing * d_lmin**2
 !! TimeStepの無次元化
 dt_Conduction = dt_Conduction / ( LRef / aRef )
 ! メモリ解放 ------------------------------------------------------------------------------------------
  deallocate( d_l )
End Subroutine Calculation_HeatConductionTimeStep
!*******************************************************************************************************
!******** 翼表面の熱伝導計算(有次元)								********
!*******************************************************************************************************
Subroutine Calculation_HeatConduction( mNum, ihs, ihe, d_t, nnNum )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in) :: nnNum, mNum, ihs, ihe
 real   , intent(in) :: d_t					! 時間刻み幅(非定常計算用)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer       :: i_h, j_h, k_h, m_h, ipt
 integer       :: HeatType						! 加熱方式
 real, pointer :: Rh_h(:), yp_h(:), Velp_h(:), Tph(:), Theta(:), d_th(:), q_input(:)
 real, pointer :: T0_wall(:), T_wall(:)				! 温度(時間発展前・後)
 real, pointer :: Utau_h(:)					! 摩擦速度
 real, pointer :: hc_h(:)					! 熱伝達係数
 real, pointer :: d_s(:)					! 要素の長さ
 real, pointer :: d_l(:)					! 要素の中心点間の長さ
 real, pointer :: LT(:), UT1(:), UT2(:), Th(:), bh(:)		! LU分解用の配列
 Real             :: A0, A1, A2, AA1, AA2, Pressurep		! LU分解用の変数,第1格子点上の圧力
 Real             :: uch, vch, wch, uph, vph, wph		! 境界上と第1格子点上の速度
 Real             :: ElectricPower, S_HA, InputPower, P_Max
! Real             :: k_wing, c_wing, rho_wing			! 翼の物性値(熱伝導率,比熱,密度)
! Real             :: thick					! 翼の厚さ
 Logical          :: fOutputLog_T				! ファイル出力の判定
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real, parameter :: dn = 3.0					! 拡散数
 real, parameter :: k_air    = 24.2e-3				! 空気の熱伝導率 [W/m*K]
 real, parameter :: Thick    = 3.0e-6				! 翼材(Ti)の厚さ[m]
 real, parameter :: k_wing        = 17.0			! 翼材(Ti)の熱伝導率[W/(m・K)]
 real, parameter :: c_wing        = 519.0			! 翼材(Ti)の比熱[J/(kg・K)]
 real, parameter :: rho_wing      = 4510.0			! 翼材(Ti)の密度[kg/(m^3)]
 ! 処理開始 ********************************************************************************************
 ! 加熱方式の選択( 1⇒温度一定; 2⇒等熱流束 )
 HeatType = 1
 m_h = mNum
 j_h = Flw(m_h)%js
 !!++++++++++++++++++++++++++++++++++++++++++
 k_h = Flw(m_h)%ks + 1 !!!k方向も後で拡張する
 !!++++++++++++++++++++++++++++++++++++++++++
! ! 翼の物性値(アルミ)
! k_wing = 256.05; c_wing = 913.0; rho_wing = 2700.0 !![W/(m・k)] [J/(kg・K)] [kg/m^3]!!
! ! 翼の厚さ(2mmとする)
! thick = 0.002
 ! メモリ確保 ------------------------------------------------------------------------------------------
 allocate( T_wall    ( ihs: ihe ), T0_wall   ( ihs: ihe ), Rh_h      ( ihs: ihe ), &
 &         Utau_h    ( ihs: ihe ), yp_h      ( ihs: ihe ), Velp_h    ( ihs: ihe ), &
 &         hc_h      ( ihs: ihe ), d_s       ( ihs: ihe ), d_l       ( ihs: ihe ), &
 &         bh        ( ihs: ihe ), LT        ( ihs: ihe ), UT1       ( ihs: ihe ), &
 &         d_th      ( ihs: ihe ), UT2       ( ihs: ihe ), Th        ( ihs: ihe ), &
 &         Theta     ( ihs: ihe ), q_input   ( ihs: ihe ), Tph       ( ihs: ihe ) )
 ! 温度の初期値
 Do i_h = ihs, ihe
  T0_wall(i_h) = ( Flw(m_h)%qh0(i_h,j_h,k_h,5) - Flw(m_h)%qh0(i_h,j_h,k_h,6) &
  &             - 0.5 * (  Flw(m_h)%qh0(i_h,j_h,k_h,2)**2 &
  &                      + Flw(m_h)%qh0(i_h,j_h,k_h,3)**2 &
  &                      + Flw(m_h)%qh0(i_h,j_h,k_h,4)**2 ) / Flw(m_h)%qh0(i_h,j_h,k_h,1) &
  &                      ) * (gamma - 1.0)  / ( Flw(m_h)%qh0(i_h,j_h,k_h,1)  * Rg )
  T0_wall(i_h) = T0_wall(i_h) * aRef**2
 End Do
 !境界からの第一格子点の物理量,要素間距離,
 Do i_h = ihs + 1, ihe - 1
  yp_h(i_h) = sqrt(  ( Flw(m_h)%x(i_h,j_h + 1,k_h) - Flw(m_h)%x(i_h,j_h,k_h) )**2 &
  &                + ( Flw(m_h)%y(i_h,j_h + 1,k_h) - Flw(m_h)%y(i_h,j_h,k_h) )**2 &
  &                + ( Flw(m_h)%z(i_h,j_h + 1,k_h) - Flw(m_h)%z(i_h,j_h,k_h) )**2 ) * LRef
  uph =  Flw(m_h)%qh(i_h,j_h+1,k_h,2) / Flw(m_h)%qh(i_h,j_h+1,k_h,1) &
  &    - Flw(m_h)%qh(i_h,j_h  ,k_h,2) / Flw(m_h)%qh(i_h,j_h  ,k_h,1)
  vph =  Flw(m_h)%qh(i_h,j_h+1,k_h,3) / Flw(m_h)%qh(i_h,j_h+1,k_h,1) &
  &    - Flw(m_h)%qh(i_h,j_h  ,k_h,3) / Flw(m_h)%qh(i_h,j_h  ,k_h,1)
  wph =  Flw(m_h)%qh(i_h,j_h+1,k_h,4) / Flw(m_h)%qh(i_h,j_h+1,k_h,1) &
  &    - Flw(m_h)%qh(i_h,j_h  ,k_h,4) / Flw(m_h)%qh(i_h,j_h  ,k_h,1)
  Velp_h(i_h) = SQRT( uph**2 + vph**2 + wph**2 ) * aRef
  Tph(i_h) = ( Flw(m_h)%qh0(i_h,j_h+1,k_h,5) - Flw(m_h)%qh0(i_h,j_h+1,k_h,6) &
  &             - 0.5 * (  Flw(m_h)%qh0(i_h,j_h+1,k_h,2)**2 &
  &                      + Flw(m_h)%qh0(i_h,j_h+1,k_h,3)**2 &
  &                      + Flw(m_h)%qh0(i_h,j_h+1,k_h,4)**2 ) / Flw(m_h)%qh0(i_h,j_h+1,k_h,1) &
  &                      ) * (gamma - 1.0)  / ( Flw(m_h)%qh0(i_h,j_h+1,k_h,1)  * Rg )
  d_s(i_h)    = SQRT( ( Flw(m_h)%x(i_h + 1 ,j_h,k_h) - Flw(m_h)%x(i_h - 1,j_h,k_h) )**2 &
  &                 + ( Flw(m_h)%y(i_h + 1 ,j_h,k_h) - Flw(m_h)%y(i_h - 1,j_h,k_h) )**2 ) * 0.5 * LRef
  d_l(i_h)    = SQRT( ( Flw(m_h)%x(i_h ,j_h,k_h) - Flw(m_h)%x(i_h - 1,j_h,k_h) )**2 &
  &                 + ( Flw(m_h)%y(i_h ,j_h,k_h) - Flw(m_h)%y(i_h - 1,j_h,k_h) )**2 ) * LRef
 End Do
 d_l(ihs) &
 & = SQRT( ( Flw(m_h)%x(ihs,j_h,k_h) - Flw(m_h)%x(ihs + 1,j_h,k_h) )**2 &
 &       + ( Flw(m_h)%y(ihs,j_h,k_h) - Flw(m_h)%y(ihs + 1,j_h,k_h) )**2 ) * LRef
 d_l(ihe) &
 & = SQRT( ( Flw(m_h)%x(ihe,j_h,k_h) - Flw(m_h)%x(ihe - 1,j_h,k_h) )**2 &
 &       + ( Flw(m_h)%y(ihe,j_h,k_h) - Flw(m_h)%y(ihe - 1,j_h,k_h) )**2 ) * LRef
 d_s(ihs) &
 & = SQRT( ( Flw(m_h)%x(ihs,j_h,k_h) - Flw(m_h)%x(ihs + 1,j_h,k_h) )**2 &
 &       + ( Flw(m_h)%y(ihs,j_h,k_h) - Flw(m_h)%y(ihs + 1,j_h,k_h) )**2 ) * 0.5 * LRef
 d_s(ihe) &
 & = SQRT( ( Flw(m_h)%x(ihe,j_h,k_h) - Flw(m_h)%x(ihe - 1,j_h,k_h) )**2 &
 &       + ( Flw(m_h)%y(ihe,j_h,k_h) - Flw(m_h)%y(ihe - 1,j_h,k_h) )**2 ) * 0.5 * LRef
 ! 加熱面積の算出
 S_HA = 0.0
 do i_h =  ihs+1, ihe
  if( Flw(m_h)%HeatFlag(i_h) .eq. 1 )S_HA = S_HA + d_l(i_h)
 end do
 !時間ステップの設定
 if(.not. fSteady) then
   d_th(:) = d_t * ( LRef / aRef ) !時間ステップの有次元化
  else
   d_th(:) = dn * c_wing * rho_wing / k_wing * d_l(:)**2
 end if
 ! 表面粗さ --------------------------------------------------------------------------------------------
 select case( RoughNum )
  case(1)
   call RoughnessShinBond2D( &
   &      ihs, ihe, T_wall, Chord * LRef, VelExp * aRef, MVD * LRef, LWC * RhoRef, Rh_h )
  case(2)
   call RoughnessCIRAMIL2D( &
   &      ihs, ihe, T_wall, Rh_h )
  case(3)
   write(*, '(a)') '!!!!! Error : Rougness model number !!!!!'
   Stop  !** 必要になったら修正する('17/11/02時点では必要無し)
!   call WallFrictionVelocity2D( &
!   &      Flw(m)%is, Flw(m)%ie, T_wall, LWC, &
!   &      Rh_h )
!   &      Ice(m)%is, Ice(m)%ie, &
!   &      Ice(m)%fRough, Ice(m)%yp, Ice(m)%Velp, Ice(m)%nup, RH, &
!   &      utau )
!   call RoughnessWright2D( &
!   &      Ice(m)%is, Ice(m)%ie, Ice(m)%x, Ice(m)%y, &
!   &      Ice(m)%SA, Ice(m)%Up, Ice(m)%Vp, Utau, Ice(m)%Rho0, LWC, &
!   &      RH )
  case default; write(*, '(a)') '!!!!! Error : Rougness model number !!!!!'
 end select
 ! 摩擦速度 --------------------------------------------------------------------------------------------
 call WallFrictionVelocity2D( &
 &      ihs, ihe, Flw(m_h)%fRough, yp_h, Velp_h, &
 &      Flw(m_h)%mu(:,j_h + 1,Flw(m_h)%ks) * (RhoRef * aRef * LRef), Rh_h, utau_h )
 ! 熱伝達率 --------------------------------------------------------------------------------------------
 !よどみ点の位置の座標ipt
 P_Max = maxval( Flw(m_h)%p(:,j_h,k_h) )
 do i_h = ihs, ihe
  if( Flw(m_h)%p(i_h,j_h,k_h) == P_Max )then
   ipt = i_h
  end if
 end do
 call HeatTransferCoefficient2D( &
 &      ihs, ihe, ipt, Flw(m_h)%x(:,j_h,k_h) * LRef, Flw(m_h)%y(:,j_h,k_h) * LRef, &
 &      ( Flw(m_h)%qh0(:,j_h + 30,k_h,2) / Flw(m_h)%qh0(:,j_h + 30,k_h,1) ) * aRef, &
 &      ( Flw(m_h)%mu(:,j_h,k_h) + Flw(m_h)%mut(:,j_h,k_h) )* (RhoRef * aRef * LRef), &
 &      Flw(m_h)%qh0(:,j_h,k_h,1) * Flw(m_h)%jac(:,j_h,k_h) * RhoRef, utau_h, Rh_h, &
 &      hc_h )
 ! T(n+1)の係数行列AのLU分解&既知の値による列ベクトルの計算
 UT1(ihs) = 1.0; UT2(ihs) = -1.0; bh(ihs) = 0.0
 Do i_h = ihs + 1, ihe - 1
  select case(HeatType)
   case(1)
    q_input(i_h) = 0.0
   case(2)
    InputPower = 13.0 ! 出力電力[W]
    if( Flw(m_h)%HeatFlag(i_h) .eq. 1 )then
      q_input(i_h) = InputPower / ( S_HA * Span * LRef )! [W/m^2]
     else
      q_input(i_h) = 0.0 ! [W]
    end if
   case default
    write(*,'(a)') "!!! Error : HeatType Number"
    stop
  end select
  A0 = 0.5 * k_wing / ( d_s(i_h) * d_l(i_h) ) / ( rho_wing * c_wing )
  A1 = 0.5 * ( k_wing / d_s(i_h) * ( 1.0 / d_l(i_h) + 1.0 / d_l(i_h + 1) ) &
  &               + hc_h(i_h) / thick ) / ( rho_wing * c_wing )
  A2 = 0.5 * k_wing / ( d_s(i_h) * d_l(i_h + 1) ) / ( rho_wing * c_wing )
  bh(i_h) = A0 * T0_wall(i_h -1) &
  &         + ( 1.0 / d_th(i_h) - A1 ) * T0_wall(i_h) &
  &         + A2 * T0_wall(i_h + 1) &
  &         + ( q_input(i_h) + hc_h(i_h) * TsExp * aRef**2 ) / ( rho_wing * c_wing * thick )
  LT (i_h) = -1.0 * A0 / UT1(i_h - 1)
  UT1(i_h) = ( 1.0 / d_th(i_h) + A1 ) - LT (i_h) * UT2(i_h - 1)
  UT2(i_h) = -1.0 * A2
 End Do
 LT (ihe) = -1.0 / UT1(ihe - 1); UT1(ihe) = 1.0 - LT(ihe) * UT2(ihe - 1); bh(ihe) = 0.0
 ! 前進代入
 Th(ihs) = bh(ihs)
 Do i_h = ihs + 1, ihe
  Th(i_h) = bh(i_h) - LT(i_h) * Th(i_h - 1)
 End do
 ! 後退代入
 Theta(ihe) = Th(ihe) / UT1(ihe)
 Do i_h = ihe - 1, ihs, -1
  Theta(i_h) = ( Th(i_h) - UT2(i_h) * Theta(i_h + 1) ) / UT1(i_h)
 End do
 ! 値の更新(非加熱部のみ)
 Do i_h = ihs, ihe
 Select Case( Flw(m_h)%HeatFlag(i_h) )
  Case(0)
   T_wall(i_h) = Theta(i_h)
  Case(1)
   select case( HeatType )
    case(1)
     T_wall(i_h) = T0_wall(i_h)
    case(2)
     T_wall(i_h) = Theta(i_h)
    case default
     write(*,'(a)') "!!! Error : HeatType Number"
     stop
   end select
  Case Default
   Write(*,'(a, i4)') "!!! Error : HeatFlag Number", i_h
   Stop
  End Select
 End Do
 !! 計算後の翼面温度の出力
 Open(1,File = trim(FlwCalOutDir) // trim(BlkName(m_h)) // 'Temperature_wing.txt',&
 &      Form = 'formatted', Status = 'replace')
 Do i_h = ihs, ihe
  Write(1,'(e16.8e2, x, e16.8e2)') Flw(m_h)%x(i_h,j_h,k_h) / chord, T_wall(i_h)
 End Do
 Close(1)
 ! 無次元化
 T_wall(:) = T_wall(:) / aRef**2
 !流れ場の境界条件の更新
 Do k_h = Flw(m_h)%ks, Flw(m_h)%ke
  Do i_h = ihs, ihe
   ! 圧力は0次外挿
   Pressurep = ( Flw(m_h)%qh(i_h,j_h+1,k_h,5) - Flw(m_h)%qh(i_h,j_h+1,k_h,6) &
   &             - 0.5 * (  Flw(m_h)%qh(i_h,j_h+1,k_h,2)**2 &
   &                      + Flw(m_h)%qh(i_h,j_h+1,k_h,3)**2 &
   &                      + Flw(m_h)%qh(i_h,j_h+1,k_h,4)**2 ) / Flw(m_h)%qh(i_h,j_h+1,k_h,1) &
   &            ) * (gamma - 1.0) * Flw(m_h)%jac(i_h,j_h+1,k_h)
   uch = Flw(m_h)%qh(i_h,j_h,k_h,2) / Flw(m_h)%qh(i_h,j_h,k_h,1)
   vch = Flw(m_h)%qh(i_h,j_h,k_h,3) / Flw(m_h)%qh(i_h,j_h,k_h,1)
   wch = Flw(m_h)%qh(i_h,j_h,k_h,4) / Flw(m_h)%qh(i_h,j_h,k_h,1)
   ! 保存ベクトルの更新
   Flw(m_h)%qh(i_h,j_h,k_h,1) = Pressurep &
   &                            / ( Rg * T_wall(i_h) * Flw(m_h)%jac(i_h,j_h,k_h) )
   Flw(m_h)%qh(i_h,j_h,k_h,2) = uch * Flw(m_h)%qh(i_h,j_h,k_h,1)
   Flw(m_h)%qh(i_h,j_h,k_h,3) = vch * Flw(m_h)%qh(i_h,j_h,k_h,1)
   Flw(m_h)%qh(i_h,j_h,k_h,4) = wch * Flw(m_h)%qh(i_h,j_h,k_h,1)
   Flw(m_h)%qh(i_h,j_h,k_h,6) =  Flw(m_h)%qh(i_h,j_h+1,k_h,6) / Flw(m_h)%qh(i_h,j_h+1,k_h,1) &
   &                           * Flw(m_h)%qh(i_h,j_h,k_h,1)
   Flw(m_h)%qh(i_h,j_h,k_h,7) =  Flw(m_h)%qh(i_h,j_h+1,k_h,7) / Flw(m_h)%qh(i_h,j_h+1,k_h,1) &
   &                           * Flw(m_h)%qh(i_h,j_h,k_h,1)
   Flw(m_h)%qh(i_h,j_h,k_h,5) =  Rg * T_wall(i_h) / ( gamma - 1.0 ) * Flw(m_h)%qh(i_h,j_h,k_h,1) &
   &                           + 0.5 * ( Flw(m_h)%qh(i_h,j_h,k_h,2)**2 &
   &                                   + Flw(m_h)%qh(i_h,j_h,k_h,3)**2 &
   &                                   + Flw(m_h)%qh(i_h,j_h,k_h,4)**2 ) / Flw(m_h)%qh(i_h,j_h,k_h,1) &
   &                           + Flw(m_h)%qh(i_h,j_h,k_h,6)
  End Do
 End Do
 ! ファイル出力
 fOutputLog_T = mod(nnNum, nOutputLog) .eq. 0.0 .or. nnNum == nstart
 If(fOutputLog_T)Then
  select case(mNum)
   case(ms)
   write(22, '(i7, e12.4e2)') &
   &     nnNum, SUM( ABS( T_wall(:) * aRef**2 - T0_wall(:) ) )/ Real(ihe - 1)
   case(me)
   write(23, '(i7, e12.4e2)') &
   &     nnNum, SUM( ABS( T_wall(:) * aRef**2 - T0_wall(:) ) )/ Real(ihe - 1)
  end select
 End If
 ! メモリ解放
 deallocate( T_wall, Rh_h, Utau_h, yp_h, Velp_h, hc_h, d_s, d_l, Theta, bh, LT, UT1, UT2, d_th, Th, &
 &           T0_wall, q_input, Tph )
 ! 処理終了 ********************************************************************************************
 return
End Subroutine Calculation_HeatConduction
!*******************************************************************************************************
!******** 流れ場ファイル出力 (計算回数毎) 							********
!*******************************************************************************************************
subroutine OutputFileCount
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character :: fname * 30
 integer   :: m, i
 ! 処理開始 ********************************************************************************************
 ! 計算条件 --------------------------------------------------------------------------------------------
 call OUtput_CalSetting( './data/GridFlow/' // trim(ND_CalSetFile) // strtxt )
 ! 流束関数 --------------------------------------------------------------------------------------------
 write(fname, '(a, 2(i2.2,a))') 'IceStep', IceStep, 'of', IceStepMax, '_'
 do m = ms, me
  call Output_Flux3D( &
  &      trim(FlwCalOutDir) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      ls, le, Flw(m)%qh )
  call Output_Flux3D( &
  &      trim(FlwCalOutDir) // trim(fname) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      ls, le, Flw(m)%qh )
  ! 翼面温度の出力
  Open(11,File = trim(FlwCalOutDir) // trim(BlkName(m)) // 'Temperature_wing.txt',&
  &      Form = 'formatted', Status = 'replace')
  Do i = Flw(m)%is, Flw(m)%ie
   Write(11,'(e16.8e2, x, e16.8e2)') Flw(m)%x(i,Flw(m)%js,Flw(m)%ks+1) * LRef, &
   &                                 Flw(m)%t(i,Flw(m)%js,Flw(m)%ks+1) * aRef**2
  End Do
  Close(11)
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine OutputFileCount
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
  &           Flw(m)%jac, Flw(m)%qh , Flw(m)%qh0, &
  &           Flw(m)%dqc, Flw(m)%dqd, Flw(m)%dqp, Flw(m)%dqr, &
  &           Flw(m)%dqh, Flw(m)%dqh0, &
  &           Flw(m)%dt , Flw(m)%Res, Flw(m)%HeatFlag )
  deallocate( OSG(m)%ip, OSG(m)%jp, OSG(m)%kp, &
  &           OSG(m)%term1, OSG(m)%term2, OSG(m)%term3, OSG(m)%term4, &
  &           OSG(m)%term5, OSG(m)%term6, OSG(m)%term7, OSG(m)%term8, &
  &           OSG(m)%fOver )
 enddo
 ! 処理終了 ********************************************************************************************
return
end subroutine Deallocating

subroutine IcePoint(m)
 !mainroutine
 integer	:: m
 !subroutine_variable
 integer	:: i,j,k,l,ii,jj
 real		:: uc,vc,wc

 do i = Flw(m)%is,Flw(m)%ie
  do j = Flw(m)%js,Flw(m)%je
   do k = Flw(m)%ks,Flw(m)%ke
    if(Flw(m)%qh(i,j,k,1) .ne. Flw(m)%qh(i,j,k,1)) then
     write(*,*) 'rho = NaN @IcePoint',nCount
     stop
    else if(Flw(m)%qh(i,j,k,1) .eq. 0.0) then
     write(*,*) 'rho = 0',nCount
     stop
    end if
   end do
  end do
 end do

 do i = Flw(m)%is,Flw(m)%ie
  do j = Flw(m)%js,Flw(m)%ji(i) !Flw(m)%ji(i),Flw(m)%js,-1
   do k = Flw(m)%ks,Flw(m)%ke
    if((Flw(m)%ic(i,j) .ne. 12) .and. (Flw(m)%ic(i,j) .ne. 0)) then
     do ii = 0, 1
      do jj = 0, 1
       Flw(m)%u(i+ii,j+jj,k) = 0.0
       Flw(m)%v(i+ii,j+jj,k) = 0.0
       do l = 2,4
        Flw(m)%qh(i+ii,j+jj,k,l) = 0.0
       end do
       if(j+jj .ge. Flw(m)%ji(i)+1) cycle
    !   do l = ls,le
    !    Flw(m)%qh(i+ii,j+jj,k,l) = Flw(m)%qh(i+ii,Flw(m)%ji(i)+1,k,l)
    !   end do
      end do
     end do
    end if
   end do
  end do
 end do

end subroutine IcePoint

subroutine OutputPara_bin( &
&      strdir, strname, is, ie, js, je, ks, ke, &
&      rhoRef, aRef, lRef, rho, u, v, w, Ps, Ts, mu, &
&      kin, eps, mut, &
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
 &                         mut(is:ie, js:je, ks:ke)
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
 close(1)

end subroutine OutputPara_bin

subroutine calyplus(m)
 implicit none
 !mainroutine_variable
 integer,intent(in)	:: m
 !subroutine_variable
 integer		:: i,j,k
 integer		:: istart,iend
 real,allocatable	:: utau(:)
 real,allocatable	:: yplus(:)
 integer		:: j_cal
 real			:: ypmax,ypmin

 !array_allocation
 allocate(utau(Flw(m)%is:Flw(m)%ie))

 if(swi_subgrid .eq. 2) then
  allocate(yplus(Flw(m)%i1:Flw(m)%i3))
 else
  allocate(yplus(Flw(m)%is:Flw(m)%ie))
 end if 

 open(1,file = './data/GridFlow/' // 'utau.txt',status = 'old')
  do i = Flw(m)%is,Flw(m)%ie
   read(1,*) utau(i)
  end do
 close(1)

 if(swi_subgrid .eq. 2) then
  istart = Flw(m)%i1
  iend = Flw(m)%i3
 else
  istart = Flw(m)%is
  iend = Flw(m)%ie
 end if 

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
  yplus(i) = sqrt((Flw(m)%x(i,j_cal+1,Flw(m)%ks) - Flw(m)%x(i,j_cal,Flw(m)%ks))**2.0 &
           &    + (Flw(m)%y(i,j_cal+1,Flw(m)%ks) - Flw(m)%y(i,j_cal,Flw(m)%ks))**2.0) &
           &    * utau(i) / Flw(m)%mu(i,j_cal,Flw(m)%ks) * Flw(m)%rho(i,j_cal,Flw(m)%ks)
 end do

 ypmax = maxval(yplus)
 ypmin = minval(yplus)

 write(*,*) 'Range of yplus:',ypmin,'to',ypmax

end subroutine calyplus

subroutine CdCl(m)
 implicit none
 integer,intent(in)	:: m
 integer	:: i,j,k
 integer	:: istart,iend
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
 
 
 !input_utau_data
 allocate(utau(Flw(m)%is:Flw(m)%ie))
 allocate(tau(Flw(m)%is:Flw(m)%ie))
 
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
 
 !calculation_drag_and_lift
 if(swi_subgrid .eq. 2) then
  istart = Flw(m)%i1
  iend = Flw(m)%i3-1
 else
  istart = Flw(m)%is
  iend = Flw(m)%ie-1
 end if 

 do i = istart,iend
! do i = Flw(m)%i1,Flw(m)%i3-1 !check
! do i = Flw(m)%is,Flw(m)%ie-1
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
   end if
  end do
 
  dA = sqrt(dx**2.0 + dy**2.0)
 
  A = A + dA
 
  theta = atan(abs(dx) / abs(dy)) !0.5 * pi - atan(abs(dy) / abs(dx))
 
  press = 0.5 * (Flw(m)%P(i,j_cal,Flw(m)%ks) + Flw(m)%P(i+1,j_cal,Flw(M)%ks)) * rhoRef * aRef**2

!write(*,*) press
 
  tau(i) = ((0.5 * (utau(i) + utau(i+1))) * aRef)**2 * &
          & (0.5 * RhoRef * (Flw(m)%rho(i,j_cal,Flw(m)%ks) + Flw(m)%rho(i+1,j_cal,Flw(m)%ks)))
 
  Dp = Dp + sgn_p(0) * press * cos(theta) * dA
  Df = Df + swi_tau * sgn_tau(0) * tau(i) * sin(theta) * dA
 
  Lp = Lp + sgn_p(1) * press * sin(theta) * dA
  Lf = Lf + swi_tau * sgn_tau(1) * tau(i) * cos(theta) * dA
 
 end do
 
 Drag = Dp + Df
 Lift = Lp + Lf
 
 Cd = Drag / (0.5 * (RhoIn * RhoRef) * (VelExp * aRef)**2.0 * (chord * lRef))
 Cl = Lift / (0.5 * (RhoIn * RhoRef) * (VelExp * aRef)**2.0 * (chord * lRef))
 
 Re = (RhoIn * RhoRef) * (VelExp * aRef) * (chord * lRef) / (muIn * rhoRef * aRef * lRef)
 
 !output_data
 write(*,*) '--------------------------------------------'
 write(*,*) 'Reynolds number [-] : ', Re
 write(*,*) '--------------------------------------------'
 write(*,*) 'Drag Coefficient'
 write(*,*) ' Drag by Pressure : ',Dp
 write(*,*) ' Drag by Friction : ',Df
 write(*,*) ' Total Drag : ',Drag
 write(*,*) ' Cd = ',Cd
 write(*,*) '--------------------------------------------'
 write(*,*) 'Lift Coefficient'
 write(*,*) ' Lift by Pressure : ',Lp
 write(*,*) ' Lift by Friction : ',Lf
 write(*,*) ' Total Lift : ',Lift
 write(*,*) ' Cl = ',Cl
 write(*,*) '--------------------------------------------'
 
end subroutine CdCl

subroutine ICM_data(m)
 implicit none
 !main_routine_variable
 integer,intent(in)	:: m
 !subroutine_variable
 integer		:: nn
 integer		:: i,j,n,l
 integer		:: ks
 real,allocatable	:: flag(:,:)
 integer		:: start_i,end_i
 integer		:: temp_i,temp_j
 integer,allocatable	:: surf_ij(:,:)
 integer,allocatable	:: normal(:)
 integer		:: num_surf
 integer,allocatable	:: n_i(:)
 integer,allocatable	:: n_j(:)
 real			:: x1,x2,y1,y2
 real			:: Dp,Df,Lp,Lf
 real			:: A,dA
 real			:: dx,dy
 real			:: theta
 real			:: press
 real			:: RhoIn,MuIn
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

 integer	:: output_num
 real,dimension(0:1000)	:: pre_cd,pre_cl

 real,allocatable	:: RH(:)
 real,allocatable	:: hc(:)
 real,allocatable	:: Nuc(:)
 real,allocatable	:: Frc(:)
 real			:: k0
 real			:: Rek
 real			:: Nmd
 real			:: Pr
 real			:: Nu
 real,parameter		:: k_air = 24.2e-3 !空気の熱伝導率 [W/(m K)]
 real,parameter		:: Cpa = 1.006 * 1.0e3

 real			:: dist
 real			:: dist_min
 integer		:: s0
 real,allocatable	:: s(:)
 real,allocatable	:: t_vec(:,:)
 real,allocatable	:: n_vec(:,:)
 real,allocatable	:: output_xy(:,:)
 real,allocatable	:: output_tau(:,:)
 real,allocatable	:: output_t(:)
 real			:: aa,bb

 integer,parameter	:: smoothing_num = 10
 real,allocatable	:: temp_hc(:)
 real,allocatable	:: temp_t(:)
 real,allocatable	:: temp_utau(:)

 !calculate_inlet_value
 RhoIn = PsExp / (Rg * TsExp)
 MuIn = muSth * (TsExp / TsSth)**1.5 * (TsSth + s1) / (TsExp + s1)

 nn = (Flw(m)%ie - Flw(m)%is) * 10
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

 allocate(n_i(0:num_surf-1))
 allocate(n_j(0:num_surf-1))
 do n = 0,num_surf-1
  select case(normal(n))
   case(0)
    n_i(n) = surf_ij(n,0)
    n_j(n) = surf_ij(n,1) + 1
   case(-1)
    n_i(n) = surf_ij(n,0) - 1
    n_j(n) = surf_ij(n,1)
   case(1)
    n_i(n) = surf_ij(n,0) + 1
    n_j(n) = surf_ij(n,1)
  end select
 end do

 !output_vtk_file
 ks = Flw(m)%ks+1
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

 !CdCl********************************************************************
 !setting_initial_variable
 Dp = 0.0
 Df = 0.0
 Lp = 0.0
 Lf = 0.0
 A = 0.0
 !calculate_utau
 allocate(utau(0:num_surf-1))
 allocate(tau(0:num_surf-1))
 allocate(yplus(0:num_surf-1))
 do n = 0,num_surf-1
  yp = sqrt((Flw(m)%x(n_i(n),n_j(n),ks) - Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks))**2.0 + &
       &    (Flw(m)%y(n_i(n),n_j(n),ks) - Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks))**2.0)
  utau(n) = cmu * Flw(m)%kin(n_i(n),n_j(n),ks)**2.0 / (kapp * yp * Flw(m)%eps(n_i(n),n_j(n),ks))
  if(utau(n) .ne. utau(n)) utau(n) = 0.0
  yplus(n) = yp * utau(n) / Flw(m)%mu(surf_ij(n,0),surf_ij(n,1),ks) * Flw(m)%rho(surf_ij(n,0),surf_ij(n,1),ks)
 end do

 !smoothing
 allocate(temp_utau(0:num_surf-1))
 do i = 0,smoothing_num-1
  do n = 0,num_surf-1
   temp_utau(n) = utau(n)
  end do
  do n = 0+1,num_surf-1-1
   utau(n) = (temp_utau(n-1) + temp_utau(n) + temp_utau(n+1)) / 3.0
  end do
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

  if(Flw(m)%u(n_i(n),n_j(n),ks) .ne. 0.0) then
   swi_tau = Flw(m)%u(n_i(n),n_j(n),ks) / abs(Flw(m)%u(n_i(n),n_j(n),ks))
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

 write(*,*) '--------------------------------------------'
 write(*,*) 'Drag Coefficient'
 write(*,*) ' Drag by Pressure : ',Dp
 write(*,*) ' Drag by Friction : ',Df
 write(*,*) ' Total Drag : ',Drag
 write(*,*) ' Cd = ',Cd
 write(*,*) '--------------------------------------------'
 write(*,*) 'Lift Coefficient'
 write(*,*) ' Lift by Pressure : ',Lp
 write(*,*) ' Lift by Friction : ',Lf
 write(*,*) ' Total Lift : ',Lift
 write(*,*) ' Cl = ',Cl
 write(*,*) '--------------------------------------------'

 if(ncount .eq. 1) then
  output_num = 1
  open(1,file = './data/GridFlow/' // 'CdCl_log.csv',status = 'replace')
   write(1,*) output_num
   write(1,*) Cd,Cl
  close(1)
 else
  open(1,file = './data/GridFlow/' // 'CdCl_log.csv',status = 'old')
   read(1,*) output_num
   do n = 0,output_num-1
    read(1,*) pre_cd(n),pre_cl(n)
   end do
  close(1)
  open(1,file = './data/GridFlow/' // 'CdCl_log.csv',status = 'replace')
   write(1,*) output_num + 1
   do n = 0,output_num-1
    write(1,*) pre_cd(n),pre_cl(n)
   end do
   write(1,*) Cd,Cl
  close(1)
 end if

 !hc********************************************************************
 allocate(RH(0:num_surf-1))
 allocate(hc(0:num_surf-1))
 allocate(Nuc(0:num_surf-1))
 allocate(Frc(0:num_surf-1))
 allocate(s(0:num_surf-1))
 allocate(t_vec(0:num_surf-1,0:1))
 allocate(n_vec(0:num_surf-1,0:1))

 !ShinBond !修正2023.12
 k0 = 0.001177!0.628e-3
 do n = 0,num_surf-1
!  RH(n) = 0.6839 * (0.571 + 0.246 * (LWC * rhoRef * 1000.0) &
!  &       + 1.257 * (LWC * rhoRef * 1000.0)**2) * (0.047 * Flw(m)%t(surf_ij(n,0),surf_ij(n,1),ks) * aRef**2.0 - 11.27) &
!  &       * (0.4286 + 0.0044139 * VelExp * aRef * 2.237) * k0 * Chord * lRef * (1.667 - 0.00333 * MVD * lRef)
  RH(n) = 0.6839 * (0.5714 + 0.2457 * (LWC * rhoRef * 1000.0) &
  &       + 1.257 * (LWC * rhoRef * 1000.0)**2.0) * (0.047 * Flw(m)%t(surf_ij(n,0),surf_ij(n,1),ks) * aRef**2.0 - 11.27) &
  &       * k0 * Chord * lRef
  if(MVD*lRef .gt.  20.0) RH(n) = RH(n) * (1.667 - 0.00333 * MVD * lRef)
 end do
 !Fallast
 do n = 0,num_surf-1
  Rek = Flw(m)%rho(surf_ij(n,0),surf_ij(n,1),ks) * rhoRef * utau(n) * aRef * RH(n) &
  &     / (Flw(m)%mu(surf_ij(n,0),surf_ij(n,1),ks) * rhoRef * aRef * lRef)
  Nmd = 7.655 * 10.0**(-5.0) * TsExp * aRef**2.0 + 0.024139
  Pr = Flw(m)%mu(surf_ij(n,0),surf_ij(n,1),ks) * rhoRef * aRef * lRef * Cpa / Nmd
  Nu = 2.0 + 0.66 / (1.0 + (0.84 * Pr**(1.0 / 6.0))**3.0)**(1.0 / 3.0) * (Rek * Pr)**1.7 &
  &    / (1.0 + (Rek * Pr)**1.2)
  hc(n) = Nu * Nmd / (MVD * lRef) ! [J/(m^2*K*s)]
 end do

 !smoothing
 allocate(temp_hc(0:num_surf-1))
 do i = 0,smoothing_num-1
  do n = 0,num_surf-1
   temp_hc(n) = hc(n)
  end do
  do n = 0+1,num_surf-1-1
   hc(n) = (temp_hc(n-1) + temp_hc(n) + temp_hc(n+1)) / 3.0
  end do
 end do

 do n = 0,num_surf-1
  Nuc(n) = hc(n) * (chord * lRef) / k_air
  Frc(n) = Nuc(n) / sqrt(Re)
 end do

 !output********************************************************************
 dist_min = 1.0e10
 do n = 0,num_surf-1
  dist = Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks)
  if(Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks) .ge. 0.0 .and. dist .lt. dist_min) then
   dist_min = dist
   s0 = n
  end if
 end do
 s = 0.0
 do n = s0+1,num_surf-1
  s(n) = s(n-1) + sqrt((Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks) - Flw(m)%x(surf_ij(n-1,0),surf_ij(n-1,1),ks))**2.0 + &
  &                    (Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks) - Flw(m)%y(surf_ij(n-1,0),surf_ij(n-1,1),ks))**2.0)
 end do
 do n = s0-1,0,-1
  s(n) = s(n+1) + sqrt((Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks) - Flw(m)%x(surf_ij(n+1,0),surf_ij(n+1,1),ks))**2.0 + &
  &                    (Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks) - Flw(m)%y(surf_ij(n+1,0),surf_ij(n+1,1),ks))**2.0)
 end do

 open(1,file = './data/GridFlow/' // 'FrosslingNumber.txt',status = 'replace')
  do n = 0,num_surf-1
   write(1,*) s(n) / chord,Frc(n)
  end do
 close(1)

 do n = 0,num_surf-1
  if(n .eq. 0) then
   t_vec(n,0) = (Flw(m)%x(surf_ij(n+1,0),surf_ij(n+1,1),ks) - Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks)) &
   &            / sqrt((Flw(m)%x(surf_ij(n+1,0),surf_ij(n+1,1),ks) - Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks))**2.0 &
   &                 + (Flw(m)%y(surf_ij(n+1,0),surf_ij(n+1,1),ks) - Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks))**2.0)
   t_vec(n,1) = (Flw(m)%y(surf_ij(n+1,0),surf_ij(n+1,1),ks) - Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks)) &
   &            / sqrt((Flw(m)%x(surf_ij(n+1,0),surf_ij(n+1,1),ks) - Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks))**2.0 &
   &                 + (Flw(m)%y(surf_ij(n+1,0),surf_ij(n+1,1),ks) - Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks))**2.0)
  else if(n .eq. num_surf-1) then
   t_vec(n,0) = (Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks) - Flw(m)%x(surf_ij(n-1,0),surf_ij(n-1,1),ks)) &
   &            / sqrt((Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks) - Flw(m)%x(surf_ij(n-1,0),surf_ij(n-1,1),ks))**2.0 &
   &                 + (Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks) - Flw(m)%y(surf_ij(n-1,0),surf_ij(n-1,1),ks))**2.0)
   t_vec(n,1) = (Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks) - Flw(m)%y(surf_ij(n-1,0),surf_ij(n-1,1),ks)) &
   &            / sqrt((Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks) - Flw(m)%x(surf_ij(n-1,0),surf_ij(n-1,1),ks))**2.0 &
   &                 + (Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks) - Flw(m)%y(surf_ij(n-1,0),surf_ij(n-1,1),ks))**2.0)
  else
   t_vec(n,0) = (Flw(m)%x(surf_ij(n+1,0),surf_ij(n+1,1),ks) - Flw(m)%x(surf_ij(n-1,0),surf_ij(n-1,1),ks)) &
   &            / sqrt((Flw(m)%x(surf_ij(n+1,0),surf_ij(n+1,1),ks) - Flw(m)%x(surf_ij(n-1,0),surf_ij(n-1,1),ks))**2.0 &
   &                 + (Flw(m)%y(surf_ij(n+1,0),surf_ij(n+1,1),ks) - Flw(m)%y(surf_ij(n-1,0),surf_ij(n-1,1),ks))**2.0)
   t_vec(n,1) = (Flw(m)%y(surf_ij(n+1,0),surf_ij(n+1,1),ks) - Flw(m)%y(surf_ij(n-1,0),surf_ij(n-1,1),ks)) &
   &            / sqrt((Flw(m)%x(surf_ij(n+1,0),surf_ij(n+1,1),ks) - Flw(m)%x(surf_ij(n-1,0),surf_ij(n-1,1),ks))**2.0 &
   &                 + (Flw(m)%y(surf_ij(n+1,0),surf_ij(n+1,1),ks) - Flw(m)%y(surf_ij(n-1,0),surf_ij(n-1,1),ks))**2.0)
  end if
  l = 0
  do
   if(Flw(m)%u(surf_ij(n,0),surf_ij(n,1)+l,ks) .ne. 0.0) then
    t_vec(n,0) = sign(t_vec(n,0),Flw(m)%u(surf_ij(n,0),surf_ij(n,1)+l,ks))
    t_vec(n,1) = sign(t_vec(n,1),Flw(m)%v(surf_ij(n,0),surf_ij(n,1)+l,ks))
    exit
   else
    l = l + 1
   end if
  end do
!  if(sqrt(Flw(m)%u(n_i(n),n_j(n),ks)**2.0 + Flw(m)%v(n_i(n),n_j(n),ks)**2.0) .ne. 0.0) then
!   t_vec(n,0) = Flw(m)%u(n_i(n),n_j(n),ks) / sqrt(Flw(m)%u(n_i(n),n_j(n),ks)**2.0 + Flw(m)%v(n_i(n),n_j(n),ks)**2.0)
!   t_vec(n,1) = Flw(m)%v(n_i(n),n_j(n),ks) / sqrt(Flw(m)%u(n_i(n),n_j(n),ks)**2.0 + Flw(m)%v(n_i(n),n_j(n),ks)**2.0)
!  else
!   t_vec(n,0) = 0.0
!   t_vec(n,1) = 0.0
!  end if
  n_vec(n,0) = (Flw(m)%x(n_i(n),n_j(n),ks) - Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks)) &
  &            / sqrt((Flw(m)%x(n_i(n),n_j(n),ks) - Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks))**2.0 &
  &                 + (Flw(m)%y(n_i(n),n_j(n),ks) - Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks))**2.0)
  n_vec(n,1) = (Flw(m)%y(n_i(n),n_j(n),ks) - Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks)) &
  &            / sqrt((Flw(m)%x(n_i(n),n_j(n),ks) - Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks))**2.0 &
  &                 + (Flw(m)%y(n_i(n),n_j(n),ks) - Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks))**2.0)
 end do

 allocate(output_xy(0:num_surf-1,0:1))
 allocate(output_tau(0:num_surf-1,0:1))
 allocate(output_t(0:num_surf-1))
 do n = 0,num_surf-1
  aa = Flw(m)%x(surf_ij(n,0),surf_ij(n,1),ks)*lRef
  bb = Flw(m)%y(surf_ij(n,0),surf_ij(n,1),ks)*lRef
  output_xy(n,0) = aa * cos(AOA) + bb * sin(AOA)
  output_xy(n,1) = -aa * sin(AOA) + bb * cos(AOA)
  aa = (utau(n)*aRef)**2.0*Flw(m)%rho(surf_ij(n,0),surf_ij(n,1),ks)*rhoRef*t_vec(n,0)
  bb = (utau(n)*aRef)**2.0*Flw(m)%rho(surf_ij(n,0),surf_ij(n,1),ks)*rhoRef*t_vec(n,1)
  output_tau(n,0) = aa * cos(AOA) + bb * sin(AOA)
  output_tau(n,1) = -aa * sin(AOA) + bb * cos(AOA)
 end do

 !smoothing
 allocate(temp_t(0:num_surf-1))
 do i = 0,smoothing_num-1
  do n = 0,num_surf-1
   temp_t(n) = Flw(m)%t(n_i(n),n_j(n),ks)*aRef**2.0
   output_t(n) = temp_t(n)
  end do
  do n = 0+1,num_surf-1-1
   output_t(n) = (temp_t(n-1) + temp_t(n) + temp_t(n+1)) / 3.0
  end do
 end do

 open(1,file = './data/GridFlow/' // 'GridFlowData.dat',status = 'replace')
  write(1,*) num_surf
  do n = 0,num_surf-1
   write(1,'(f10.6,f10.6,f10.6,f10.6,f10.4,f10.4)') &
   &          output_xy(n,0),output_xy(n,1),output_tau(n,0),output_tau(n,1), &
   &          hc(n),output_t(n)
  end do
 close(1)

 open(1,file = './data/GridFlow/' // 'GridFlowData.vtk',status = 'replace')
  write(1,'(A)') '# vtk DataFile Version 3.0'
  write(1,*) 'vtk output'
  write(1,*) 'ASCII'
  write(1,*) 'DATASET POLYDATA'
  write(1,*) 'POINTS ',num_surf,' float'
  do n = 0,num_surf-1
   write(1,*) output_xy(n,0),output_xy(n,1),0.0
  end do
  write(1,*) 'POINT_DATA ', num_surf
  write(1,*) 'VECTORS tau float'
  do n = 0,num_surf-1
   write(1,*) output_tau(n,0),output_tau(n,1),0.0
  end do
  write(1,*) 'VECTORS hc float'
  do n = 0,num_surf-1
   write(1,*) hc(n)*n_vec(n,0),hc(n)*n_vec(n,1),0.0
  end do
  write(1,*) 'VECTORS surface_temperature float'
  do n = 0,num_surf-1
   write(1,*) output_t(n)*n_vec(n,0),output_t(n)*n_vec(n,1),0.0
  end do
 close(1)

 deallocate(flag)
 deallocate(surf_ij)
 deallocate(normal)
 deallocate(utau)
 deallocate(tau)
 deallocate(yplus)
 deallocate(RH)
 deallocate(hc)
 deallocate(Nuc)
 deallocate(Frc)
 deallocate(s)
 deallocate(t_vec)
 deallocate(n_vec)
 deallocate(temp_hc)
 deallocate(temp_t)
 deallocate(output_xy)
 deallocate(output_tau)
 deallocate(output_t)

end subroutine ICM_data

! 定義終了 *********************************************************************************************
end program FlowField_NACA
