!*******************************************************************************************************
!*******************************************************************************************************
!******** 流れ場計算プログラム                                                                        ********
!******** (NACA翼，三次元圧縮性乱流場，重合格子法，高レイノルズ数型 k-eモデル)                          ********
!********                                              2012.02.02  PROGRAMMED BY RYOSUKE HAYASHI ********
!********                                              2013.07.18     UPDATED BY RYOSUKE HAYASHI ********
!********                                              2014.04.15     UPDATED BY RYOSUKE HAYASHI ********
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
   ! ファイル名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   character, parameter :: ResidualFile*8 = 'Residual'
   ! 処理開始 ********************************************************************************************
   ! 検証実験ケース選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   write (*, '(a)') "<< Exp. Case Selecation >>"
   call SelectExpCase
   ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   write (*, '(/,a)') "<< Initial Setting >>"
   call InitialSetting
   ! 流れ場計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   write (*, '(/,a)') "<< Flow Field Computation >>"
   call CalFlowField
   ! メモリ解法 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   write (*, '(/,a)') "<< Computation Finallized >>"
   call Deallocating
   ! 内部手続き ******************************************************************************************
   stop
contains
!*******************************************************************************************************
!******** 検証実験ケース選択                                                                         ********
!*******************************************************************************************************
   subroutine SelectExpCase
      ! 変数宣言 ********************************************************************************************
      implicit none
      character :: fname*20
      ! 処理開始 ********************************************************************************************
      ! 計算条件ファイル入力
      call Input_CalSetting(trim(ND_CalSetFile)//strtxt)
      ! ディレクトリ設定
      if (IceStep == 0) then
         GrdInDir = bckdir//'grid//clean//'
         OSGDir = bckdir//'overset//clean//'
         FlwIniDir = bckdir//'flow//initial//clean//'
         FlwCalInDir = bckdir//'flow//cal//clean//'
         FlwCalOutDir = bckdir//'flow//cal//clean//'
      else
         GrdInDir = bckdir//'grid//icing//'
         OSGDir = bckdir//'overset//icing//'
         FlwIniDir = bckdir//'flow//initial//icing//'
         FlwCalInDir = bckdir//'flow//cal//icing//'
         FlwCalOutDir = bckdir//'flow//cal//icing//'
         IceCalInDir = bckdir//'icing//cal//'
      end if
      write (*, '(a)') '+++ Icing Step +++'
      write (*, '(a,i2)') '* Ice step      = ', IceStep
      write (*, '(a,i2)') '* Ice step max. = ', IceStepMax
      write (*, '(/,a)') '+++ Exp. Condition +++'
      write (*, '(a,e16.8e3)') '* Ts    = ', TsExp*aRef**2
      write (*, '(a,e16.8e3)') '* Ps    = ', PsExp*(rhoRef*aRef**2)
      write (*, '(a,e16.8e3)') '* V     = ', VelExp*aRef
      write (*, '(a,e16.8e3)') '* LWC   = ', LWC*RhoRef
      write (*, '(a,e16.8e3)') '* MVD   = ', MVD*LRef
      write (*, '(a,e16.8e3)') '* Rho   = ', Rhod*RhoRef
      write (*, '(a,e16.8e3)') '* Chord = ', Chord*LRef
      write (*, '(a,e16.8e3)') '* AOA   = ', AOA*180.0/pi
      ! 処理終了 ********************************************************************************************
      return
   end subroutine SelectExpCase
!*******************************************************************************************************
!******** 初期設定                                                                                ********
!*******************************************************************************************************
   subroutine InitialSetting
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, parameter :: nSmooth = 100
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer   :: m, n
      integer        :: j
      character :: fname*30
      ! 処理開始 ********************************************************************************************
      ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! allocate( Flw(ms:me), OSG(ms:me) )
      allocate (Flw(ms:me))
      ! ブロック名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call Set_BlockName
      ! 格子解像度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
      m = me
      call Input_Resolution3D( &
      &      trim(GrdInDir)//trim(BlkName(m))//trim(RslFile), strtxt, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke)
! enddo
      ! メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
      m = me
      allocate (Flw(m)%x(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%y(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%z(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%rho(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%u(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%v(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%w(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%p(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%t(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%mu(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%kin(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%eps(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%mut(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%xix(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%xiy(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%xiz(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%etx(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%ety(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%etz(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%zex(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%zey(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%zez(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%jac(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%qh(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke, ls:le), &
      &         Flw(m)%qh0(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke, ls:le), &
      &         Flw(m)%dqh(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke, ls:le), &
      &         Flw(m)%dqh0(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke, ls:le), &
      &         Flw(m)%dqc(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke, ls:le), &
      &         Flw(m)%dqd(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke, ls:le), &
      &         Flw(m)%dqp(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke, ls:le), &
      &         Flw(m)%dqr(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke, ls:le), &
      &         Flw(m)%dt(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%Res(ls:le))
! enddo
      ! 格子座標 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
      m = me
      call Input_Grid3D( &
      &      trim(FlwIniDir)//trim(BlkName(m))//trim(ND_GrdFile), strbin, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
      &      Flw(m)%x, Flw(m)%y, Flw(m)%z)
! enddo
      ! メトリックス ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
      m = me
      call Input_Metrics3D( &
      &      trim(FlwIniDir)//trim(BlkName(m))//trim(ND_MetFile), strbin, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
      &      Flw(m)%jac, &
      &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
      &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
      &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez)
! enddo
      ! 流束関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      nStart = nCount; time = 0.0
      if (nCount == 0) then
!   do m = ms, me
         m = me
         call Input_Flux3D( &
         &      trim(FlwIniDir)//trim(BlkName(m))//trim(ND_IniFlxFile), strbin, &
         &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
         &      ls, le, Flw(m)%qh)
!   enddo
      else
!   do m = ms, me
         m = me
         call Input_Flux3D( &
         &      trim(FlwCalInDir)//trim(BlkName(m))//trim(ND_FlxFile), strbin, &
         &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
         &      ls, le, Flw(m)%qh)
!   enddo
      end if
      ! 物理量 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
      m = me
      call SetPhysics3DKEM( &
      &      Rg, gamma, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
      &      Flw(m)%qh, Flw(m)%jac, &
      &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps)
! enddo
      ! C 型格子分割点 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
      m = me
      call Input_CtypeGridPoint( &
      &      trim(GrdInDir)//trim(BlkName(m))//trim(CtypePointFile)//strtxt, &
      &      Flw(m)%i1, Flw(m)%i2, Flw(m)%i3)
! end do
      ! 重合格子補間係数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
      ! 計算除去点 (氷の中の点) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! m = 1
! allocate( IceIn( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ) )
! if( IceStep == 0 ) then
!   IceIn(:, :, :) = 0
!  else
!   call Input_ArrayInt3D( &
!   &      trim(OSGDir) // trim(BlkName(m)) // trim(IceInPointFile), strbin, &
!   &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
!   &      IceIn )
!   IceIn(:, :, :) = 0
! endif
      ! スムージング ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (IceStep /= 0) then
         do n = 1, nSmooth
!   do m = ms, me
            m = me
            call SmoothingFlux3D( &
            &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
            &      Flw(m)%jac, Flw(m)%qh)
!   enddo
!   call InterpolationOversetGrid
         end do
      end if
      ! 粗さのフラグ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! m = 1
      m = me
      allocate (Flw(m)%fRough(Flw(m)%is:Flw(m)%ie, Flw(m)%ks:Flw(m)%ke))
      Flw(m)%fRough(:, :) = 0
      m = 2
      allocate (Flw(m)%fRough(Flw(m)%is:Flw(m)%ie, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%RH(Flw(m)%is:Flw(m)%ie, Flw(m)%ks:Flw(m)%ke))
      if (IceStep > 0) then
         call Input_ArrayInt2D( &
         &      trim(IceCalInDir)//trim(BlkName(m))//trim(RoughFlagFile), strdat, &
         &      Flw(m)%is, Flw(m)%ie, Flw(m)%ks, Flw(m)%ke, &
         &      Flw(m)%fRough)
      else
         Flw(m)%fRough(:, :) = 0
      end if
      ! 処理終了 ********************************************************************************************
      return
   end subroutine InitialSetting
!*******************************************************************************************************
!******** 流れ場計算                                                                                 ********
!*******************************************************************************************************
   subroutine CalFlowField
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 局所定数 ********************************************************************************************
      real, parameter :: dtMax = 1.0e0
      real, parameter :: dtMin = 1.0e-20
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer :: m, nn
      integer        :: i
      integer :: nCompCount
      real    :: dt, dtTmp
      real    :: RK(1:nRunge)
      logical :: fOutputCount, fOutputLog
      integer   :: fnum = 0
      character :: num*4
      character :: fstep*2
      character :: fdir*100
      character :: fn*100
      ! 処理開始 ********************************************************************************************
      write (*, '(a)') '+++ Computational loop start +++'
      ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 収束判定ファイル出力 --------------------------------------------------------------------------------
      open (20, file=trim(FlwCalOutDir)//trim(BlkName(1))//trim(ResidualFile)//strtxt)
      open (21, file=trim(FlwCalOutDir)//trim(BlkName(2))//trim(ResidualFile)//strtxt)
      ! 流れ場計算ループ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (IceStep == 0) then
         nCompCount = nCalMax !* 2.0
      else
         nCompCount = nCalMax
      end if
      do nCount = nStart, nCompCount
         ! 時間刻み幅 -----------------------------------------------------------------------------------------
!  do m = ms, me
         m = me
         ! 局所時間刻み
         call CalLocalDt3D( &
         &      cn, Rg, gamma, &
         &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
         &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
         &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
         &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
         &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%t, Flw(m)%mu, &
         &      Flw(m)%dt)
!  enddo
         ! 非定常計算 -----------------------------------------------------------------------------------------
         if (.not. fSteady) then
            dt = dtMax
!    do m = ms, me
            m = me
            dt = min(minval(Flw(m)%dt(:, :, :)), dt)
!    enddo
!    do m = ms, me
            Flw(m)%dt(:, :, :) = dt
!    enddo
         end if
         ! 過去の流束を保存 -----------------------------------------------------------------------------------
!  do m = ms, me
         m = me
         call SaveFlux3DKEM( &
         &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
         &      Flw(m)%qh, Flw(m)%dqh, Flw(m)%qh0, Flw(m)%dqh0)
!  enddo
         ! 内部反復 -------------------------------------------------------------------------------------------
         do nn = 1, nRunge
!   do m = ms, me
            m = me
            ! 粘性係数
            call ViscosityCoefficient3D( &
            &      muSth, TsSth, s1, &
            &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
            &      Flw(m)%t, &
            &      Flw(m)%mu)
            ! 対流項
            call Convection3D( &
            &      nTVD, eTVD, Rg, gamma, &
            &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
            &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
            &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
            &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
            &      Flw(m)%jac, Flw(m)%qh, &
            &      Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, &
            &      Flw(m)%dqc)
            ! 拡散項 & 生産項
            select case (TurbNum)
               ! 層流
            case (0)
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
               Flw(m)%dqp(:, :, :, :) = 0.0
               ! Launder-Spalding
            case (1)
               call Turbulence3DEvmStd( &
               &             LmtPro, Rg, gamma, Pr, Prt, &
               &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
               &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
               &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
               &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
               &      Flw(m)%jac, Flw(m)%qh, &
               &      Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%t, Flw(m)%kin, Flw(m)%eps, Flw(m)%mu, &
               &      Flw(m)%mut, Flw(m)%dqd, Flw(m)%dqp)
               ! Kato-Launder
            case (2)
               call Turbulence3DEvmStdKL( &
               &             LmtPro, Rg, gamma, Pr, Prt, &
               &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
               &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
               &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
               &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
               &      Flw(m)%jac, Flw(m)%qh, &
               &      Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%t, Flw(m)%kin, Flw(m)%eps, Flw(m)%mu, &
               &      Flw(m)%mut, Flw(m)%dqd, Flw(m)%dqp)
               ! Launder-Sharma
            case (3)
               call Turbulence3DEvmLS( &
               &             LmtPro, Rg, gamma, Pr, Prt, &
               &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
               &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
               &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
               &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
               &      Flw(m)%jac, &
               &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%t, Flw(m)%kin, Flw(m)%eps, Flw(m)%mu, &
               &      Flw(m)%mut, Flw(m)%dqd, Flw(m)%dqp)
               ! AKN k-e
            case (4)
               call Turbulence3DEvmAKN( &
               &             LmtPro, Rg, gamma, Pr, Prt, &
               &      Flw(m)%is, Flw(m)%ie, Flw(m)%i1, Flw(m)%i3, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
               &      Flw(m)%x, Flw(m)%y, Flw(m)%z, &
               &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
               &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
               &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
               &      Flw(m)%jac, &
               &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%t, Flw(m)%kin, Flw(m)%eps, Flw(m)%mu, &
               &      Flw(m)%mut, Flw(m)%dqd, Flw(m)%dqp)
            case default; write (*, '(a)') '!!!!! Error : TurbNum !!!!!'
            end select
            ! 外力項
            Flw(m)%dqr(:, :, :, :) = 0.0
            ! 各項の和
            call SumDQH3D( &
            &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
            &      Flw(m)%dqc, Flw(m)%dqd, Flw(m)%dqp, Flw(m)%dqr, &
            &      Flw(m)%dqh)
            ! 未来の流束
            if (fTime) then
               call RungeKutta3D( &
               &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
               &      nn, nRunge, Flw(m)%dt, &
               &      Flw(m)%dqh, Flw(m)%qh0, Flw(m)%qh)
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
               &      Flw(m)%dqh0, Flw(m)%dqh, Flw(m)%qh0, Flw(m)%qh)
            end if
            ! 乱流量のリミッター
            if (LmtAve /= 0.0) then
               call Limiter3DKEM( &
               &      LmtAve, &
               &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
               &      Flw(m)%jac, &
               &      Flw(m)%qh)
            end if
!   enddo
            ! 境界条件
            call BoundaryCondition
!   ! 重合格子の補間
!   call InterpolationOversetGrid
            ! 物理量
!   do m = ms, me
            call SetPhysics3DKEM( &
            &      Rg, gamma, &
            &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
            &      Flw(m)%qh, Flw(m)%jac, &
            &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps)
!   enddo
         end do
         ! 経過時間 -------------------------------------------------------------------------------------------
         if (.not. fSteady) then
            time = time + dt
         end if
         ! ファイル出力 ---------------------------------------------------------------------------------------
         fOutputCount = mod(nCount, nOutputCount) .eq. 0.0 .or. nCount == nCompCount
         fOutputLog = mod(nCount, nOutputLog) .eq. 0.0 .or. nCount == nCompCount
         if (fOutputLog .or. fOutputCount) then
            if (nCount == nStart) then
               write (*, '(a)') '+++ Numerical Condition +++'
               write (*, '(a,e10.4e1)') '* Cn                  = ', Cn
               write (*, '(a,i2)') '* TVD Order           = ', nTVD
               write (*, '(a,e10.4e1)') '* TVD Entropy         = ', eTVD
               write (*, '(a,e10.4e1)') '* Limitter k-e Ave    = ', LmtAve
               write (*, '(a,e10.4e1)') '* Limitter k-e Pro    = ', LmtPro
               write (*, '(a,i7)') '* Computational Count = ', nCompCount
            end if
            if (fOutputLog) then
               call CalResidual
            end if
            if (fOutputCount) call OutputFileCount
            write (*, '(a)') '+++ Calculation progress +++'
            write (*, '(a,i6,a)') '* Calculation Count = ', nCount
            if (fSteady) then
!      write(*, '(a,e16.8e3)') '* Main Min. dt = ', minval(Flw(1)%dt(:,:,:))
!      write(*, '(a,e16.8e3)') '* Main Max. dt = ', maxval(Flw(1)%dt(:,:,:))
               write (*, '(a,e16.8e3)') '* Sub  Min. dt = ', minval(Flw(2)%dt(:, :, :))
               write (*, '(a,e16.8e3)') '* Sub  Max. dt = ', maxval(Flw(2)%dt(:, :, :))
            else
               write (*, '(a,e16.8e3)') '* Calculation time  = ', time
               write (*, '(a,e16.8e3)') '* Calculation dt    = ', dt
            end if

            fnum = fnum + 1
            write (num, '(I4.4)') fnum
            write (fstep, '(I2.2)') IceStep
!   do m = ms, me
            m = me
            write (*, *) 'Output File : ', trim(BlkName(m))//'Flow'//trim(num)
            if (IceStep .eq. 0) then
               fdir = bckdir//'flow//view//clean//overtime//'
            else
               fdir = bckdir//'flow//view//icing//overtime//'
            end if

!    !ParaViewファイル出力
!    fn = trim(BlkName(m)) // 'Flow' // trim(fstep) // '_' // trim(num)
!    call OutputPara_bin( &
!       &      trim(fdir), trim(fn), &
!       &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
!       &      rhoRef, aRef, lRef, &
!       &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%mu, &
!       &      Flw(m)%kin, Flw(m)%eps, Flw(m)%mut, &
!       &      Flw(m)%x, Flw(m)%y, Flw(m)%z )

            if (m .eq. me) then
               call calyplus(m)
            end if
!   end do

         end if
      end do
      close (20); close (21); close (22)
      ! 処理終了 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
   end subroutine CalFlowField
!*******************************************************************************************************
!******** 境界条件                                                                                 ********
!*******************************************************************************************************
   subroutine BoundaryCondition
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer :: i, j, k, l, m
      ! 処理開始 ********************************************************************************************
      ! Main Grid, Sub Grid++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms,me
      m = me
      ! 流入境界 --------------------------------------------------------------------------------------------
      call BoundaryInlet( &
      &      m, Flw(m)%is + 1, Flw(m)%ie - 1, Flw(m)%ks + 1, Flw(m)%ke - 1, &
      &      Flw(m)%je - 0, Flw(m)%je - 1)
      ! C 型格子ブランチ・カット ----------------------------------------------------------------------------
      call BoundaryCtypeBranch( &
      &      m, Flw(m)%is + 1, Flw(m)%i1 - 1, Flw(m)%ks + 1, Flw(m)%ke - 1, &
      &      Flw(m)%js + 0, Flw(m)%js + 1, &
      &      Flw(m)%ie)
! call BoundaryCtypeBranch( &
! &      m, Flw(m)%i3 + 1, Flw(m)%ie - 1, Flw(m)%ks + 1, Flw(m)%ke - 1, &
! &      Flw(m)%js + 0, Flw(m)%js + 1, &
! &      Flw(m)%ie )
      ! 翼壁面 ----------------------------------------------------------------------------------------------
      call BoundaryBladeSurface( &
      &      m, Flw(m)%i1 + 0, Flw(m)%i3 - 0, Flw(m)%ks + 1, Flw(m)%ke - 1, &
      &      Flw(m)%js + 0, Flw(m)%js + 1, Flw(m)%js + 2, TurbNum)
      ! 流出境界 --------------------------------------------------------------------------------------------
      call BoundaryOutlet( &
      &      m, Flw(m)%js + 0, Flw(m)%je - 0, Flw(m)%ks + 1, Flw(m)%ke - 1, &
      &      Flw(m)%is + 0, Flw(m)%is + 1)
      call BoundaryOutlet( &
      &      m, Flw(m)%js + 0, Flw(m)%je - 0, Flw(m)%ks + 1, Flw(m)%ke - 1, &
      &      Flw(m)%ie - 0, Flw(m)%ie - 1)
      ! 氷の中の点 ------------------------------------------------------------------------------------------
! call BoundaryIceIn( &
! &      1, 2, Flw(1)%i1, Flw(1)%i3, Flw(1)%js, Flw(1)%je, Flw(1)%ks + 1, Flw(1)%ke - 1, &
! &      Flw(2)%i1, Flw(2)%i2, Flw(2)%js, int(0.5 * Flw(2)%ke) )
      ! 周期境界 --------------------------------------------------------------------------------------------
      call BoundaryPeriodic( &
      &      m, Flw(m)%is + 0, Flw(m)%ie - 0, Flw(m)%js + 0, Flw(m)%je - 0, &
      &      Flw(m)%ks + 0, Flw(m)%ke - 1)
      call BoundaryPeriodic( &
      &      m, Flw(m)%is + 0, Flw(m)%ie - 0, Flw(m)%js + 0, Flw(m)%je - 0, &
      &      Flw(m)%ke - 0, Flw(m)%ks + 1)
! end do
      ! 処理終了 ********************************************************************************************
      return
   end subroutine BoundaryCondition
!*******************************************************************************************************
!******** 流入境界                                                                                ********
!******** ⅰ．密度外挿，その他固定　                        ⅱ．角度・全温・全圧固定，マッハ数外挿        ********
!******** ⅲ．角度・体積流量・全温固定，密度外挿                                                  ********
!*******************************************************************************************************
   subroutine BoundaryInlet( &
   &            m, is, ie, ks, ke, j0, j1)
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
      select case (BCNum)
      case (1)
         do k = ks, ke
         do i = is, ie
            ! 密度外挿
            rho = Flw(m)%qh(i, j1, k, 1)*Flw(m)%jac(i, j1, k)
            ! その他固定
            vel = VelIn
            u = vel*cos(AOA)
            v = vel*sin(AOA)
            w = 0.0
            t = TsExp/aRef**2
            kin = 0.5*3.0*(0.01*vel)**2
            eps = rho*kin**2/(Flw(m)%mu(i, j0, k)*Ret)
            ! 流束関数
            Flw(m)%qh(i, j0, k, 1) = rho/Flw(m)%jac(i, j0, k)
            Flw(m)%qh(i, j0, k, 2) = u*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 3) = v*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 4) = w*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 5) = (Rg*t/(gamma - 1.0) + 0.5*vel**2 + kin)*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 6) = kin*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 7) = eps*Flw(m)%qh(i, j0, k, 1)
         end do
         end do
      case (2)
         do k = ks, ke
         do i = is, ie
            ! 外挿する速度からマッハ数導出
            u = Flw(m)%qh(i, j1, k, 2)/Flw(m)%qh(i, j1, k, 1)
            v = Flw(m)%qh(i, j1, k, 3)/Flw(m)%qh(i, j1, k, 1)
            w = Flw(m)%qh(i, j1, k, 4)/Flw(m)%qh(i, j1, k, 1)
            p = (Flw(m)%qh(i, j1, k, 5) - Flw(m)%qh(i, j1, k, 6) &
            &     - 0.5*(Flw(m)%qh(i, j1, k, 2)**2 &
            &             + Flw(m)%qh(i, j1, k, 3)**2 &
            &             + Flw(m)%qh(i, j1, k, 4)**2)/Flw(m)%qh(i, j1, k, 1) &
            &   )*(gamma - 1.0)*Flw(m)%jac(i, j1, k)
            t = p/(Flw(m)%qh(i, j1, k, 1)*Flw(m)%jac(i, j1, k)*Rg)
            mac = sqrt(u**2 + v**2 + w**2)/sqrt(gamma*Rg*t)
            ! 外挿するマッハ数と固定する全温・全圧から静温・静圧導出
            Tt_Ts = 1.0 + 0.5*(gamma - 1.0)*mac**2
            Ts = TtIn/Tt_Ts
            Ps = PtIn/Tt_Ts**(gamma/(gamma - 1.0))
            ! その他導出
            rho = Ps/(Rg*Ts)
            vel = mac*sqrt(gamma*Rg*Ts)
            u = vel*cos(AOA)
            v = vel*sin(AOA)
            w = 0.0
            kin = 0.5*3.0*(0.01*vel)**2
            eps = rho*kin**2/(Flw(m)%mu(i, j0, k)*Ret)
            ! 流束関数
            Flw(m)%qh(i, j0, k, 1) = rho/Flw(m)%jac(i, j0, k)
            Flw(m)%qh(i, j0, k, 2) = u*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 3) = v*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 4) = w*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 5) = (Rg*Ts/(gamma - 1.0) + 0.5*vel**2 + kin)*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 6) = kin*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 7) = eps*Flw(m)%qh(i, j0, k, 1)
         end do
         end do
      case (3)
         do k = ks, ke
         do i = is, ie
            ! 密度外挿
            rho = Flw(m)%qh(i, j1, k, 1)*Flw(m)%jac(i, j1, k)
            ! 体積流量固定
            vel = VelIn
            ! 全温固定
            Ts = TtIn - (gamma - 1.0)/(2.0*gamma*Rg)*vel**2
            Ps = rho*Rg*Ts
            ! 速度から乱流量導出
            kin = 0.5*3.0*(0.01*vel)**2
            eps = rho*kin**2/(Flw(m)%mu(i, j0, k)*Ret)
            ! 流束関数
            Flw(m)%qh(i, j0, k, 1) = rho/Flw(m)%jac(i, j0, k)
            Flw(m)%qh(i, j0, k, 2) = vel*cos(AOA)*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 3) = vel*sin(AOA)*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 4) = 0.0
            Flw(m)%qh(i, j0, k, 5) = (Rg*Ts/(gamma - 1.0) + 0.5*vel**2 + kin)*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 6) = kin*Flw(m)%qh(i, j0, k, 1)
            Flw(m)%qh(i, j0, k, 7) = eps*Flw(m)%qh(i, j0, k, 1)
         end do
         end do
      case default
         write (*, '(a)') '!!!! Error : Inlet boundary condition pattern !!!!'
      end select
      ! 処理終了 ********************************************************************************************
      return
   end subroutine BoundaryInlet
!*******************************************************************************************************
!******** 流出境界                                                                                ********
!*******************************************************************************************************
   subroutine BoundaryOutlet( &
   &            m, js, je, ks, ke, i0, i1)
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
         Ps = (Flw(m)%qh(i1, j, k, 5) - Flw(m)%qh(i1, j, k, 6) &
         &      - 0.5*(Flw(m)%qh(i1, j, k, 2)**2 &
         &              + Flw(m)%qh(i1, j, k, 3)**2 &
         &              + Flw(m)%qh(i1, j, k, 4)**2)/Flw(m)%qh(i1, j, k, 1) &
         &    )*(gamma - 1.0)*Flw(m)%jac(i1, j, k)
         Ts = Ps/(Flw(m)%qh(i1, j, k, 1)*Flw(m)%jac(i1, j, k)*Rg)
         do l = ls + 1, le
            Flw(m)%qh(i0, j, k, l) = Flw(m)%qh(i1, j, k, l)/Flw(m)%qh(i1, j, k, 1)
         end do
         Flw(m)%qh(i0, j, k, 1) = PsOut/(Rg*Ts*Flw(m)%jac(i0, j, k))
         do l = ls + 1, le
            Flw(m)%qh(i0, j, k, l) = Flw(m)%qh(i0, j, k, l)*Flw(m)%qh(i0, j, k, 1)
         end do
      end do
      end do
      ! 処理終了 ********************************************************************************************
      return
   end subroutine BoundaryOutlet
!*******************************************************************************************************
!******** C 型格子ブランチ・カット境界 (平均値外挿)                                                   ********
!*******************************************************************************************************
   subroutine BoundaryCtypeBranch( &
   &            m, is, ie, ks, ke, j0, j1, imax)
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
         qAve = 0.5*(Flw(m)%qh(i, j1, k, l)*Flw(m)%jac(i, j1, k) &
         &            + Flw(m)%qh(imax - i, j1, k, l)*Flw(m)%jac(imax - i, j1, k))
         Flw(m)%qh(i, j0, k, l) = qAve/Flw(m)%jac(i, j0, k)
         Flw(m)%qh(imax - i, j0, k, l) = qAve/Flw(m)%jac(imax - i, j0, k)
      end do
      end do
      end do
      ! 処理終了 ********************************************************************************************
      return
   end subroutine BoundaryCtypeBranch
!*******************************************************************************************************
!******** 翼壁面境界                                                                                ********
!*******************************************************************************************************
   subroutine BoundaryBladeSurface( &
   &            m, is, ie, ks, ke, j0, j1, j2, TurbNum)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in) :: m
      integer, intent(in) :: is, ie, ks, ke
      integer, intent(in) :: j0, j1, j2
      integer, intent(in) :: TurbNum !低Re数型乱流モデル用
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer :: i, k
      ! 処理開始 ********************************************************************************************
      do k = ks, ke
      do i = is, ie
         if (TurbNum .eq. 4) then
            call WallNoSlipLRe(m, i, j0, k, i, j1, k)
         else
            if (fSlip) then
               if (Flw(m)%qh(i, j1, k, 1) <= 0.0 .or. Flw(m)%qh(i, j2, k, 1) <= 0.0) cycle
               call WallSlipWF(m, i, j0, k, i, j1, k, i, j2, k)
            else
               if (Flw(m)%qh(i, j1, k, 1) <= 0.0) cycle
               call WallNoSlipWF(m, i, j0, k, i, j1, k)
            end if
         end if
      end do
      end do
      ! 処理終了 ********************************************************************************************
      return
   end subroutine BoundaryBladeSurface
!*******************************************************************************************************
!******** 壁境界 (滑りなし, 断熱, 壁関数)                                                        ********
!*******************************************************************************************************
   subroutine WallNoslipWF( &
   &            m, i0, j0, k0, i1, j1, k1)
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
      ! 壁近傍の分布 ----------------------------------------------------------------------------------------
      dx = Flw(m)%x(i1, j1, k1) - Flw(m)%x(i0, j0, k0)
      dy = Flw(m)%y(i1, j1, k1) - Flw(m)%y(i0, j0, k0)
      dz = Flw(m)%z(i1, j1, k1) - Flw(m)%z(i0, j0, k0)
      uc = Flw(m)%qh(i1, j1, k1, 2)/Flw(m)%qh(i1, j1, k1, 1)
      vc = Flw(m)%qh(i1, j1, k1, 3)/Flw(m)%qh(i1, j1, k1, 1)
      wc = Flw(m)%qh(i1, j1, k1, 4)/Flw(m)%qh(i1, j1, k1, 1)
      ! 壁関数 ----------------------------------------------------------------------------------------------
      yp = sqrt(dx**2 + dy**2 + dz**2)
      up = sqrt(uc**2 + vc**2 + wc**2)
      up = max(zero, up)
      nup = Flw(m)%mu(i1, j1, k1)/(Flw(m)%qh(i1, j1, k1, 1)*Flw(m)%jac(i1, j1, k1))
      kp = Flw(m)%kin(i1, j1, k1)
      epsp = Flw(m)%eps(i1, j1, k1)
      nup = max(zero, nup)
      kp = max(zero, kp)
      epsp = max(zero, epsp)
      select case (Flw(m)%fRough(i0, k0))
      case (0)
         call WallFunctionKEM2S( &
         &      yp, up, nup, utau, dudy1, dudy2, kp, epsp)
      case (1)
         call WallFunctionKEM2R( &
         &      Flw(m)%RH(i0, k0), yp, up, nup, utau, dudy1, dudy2, kp, epsp)
      case default
         write (*, '(a)') '!!!!! Error : fRough number !!!!!'
      end select
      ! 速度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 壁面 ------------------------------------------------------------------------------------------------
      uc = 0.0
      vc = 0.0
      wc = 0.0
      ! 流束関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 壁面から一点目 --------------------------------------------------------------------------------------
      ! kp, epsp
      Flw(m)%qh(i1, j1, k1, 5) = Flw(m)%qh(i1, j1, k1, 5) - Flw(m)%qh(i1, j1, k1, 6)
      Flw(m)%qh(i1, j1, k1, 6) = Flw(m)%qh(i1, j1, k1, 1)*kp
      Flw(m)%qh(i1, j1, k1, 7) = Flw(m)%qh(i1, j1, k1, 1)*epsp
      Flw(m)%qh(i1, j1, k1, 5) = Flw(m)%qh(i1, j1, k1, 5) + Flw(m)%qh(i1, j1, k1, 6)
      ! 壁面 ------------------------------------------------------------------------------------------------
      ! 全て外挿
      Flw(m)%qh(i0, j0, k0, 1) = Flw(m)%qh(i1, j1, k1, 1)*Flw(m)%jac(i1, j1, k1)/Flw(m)%jac(i0, j0, k0)
      Flw(m)%qh(i0, j0, k0, 2) = uc*Flw(m)%qh(i0, j0, k0, 1)
      Flw(m)%qh(i0, j0, k0, 3) = vc*Flw(m)%qh(i0, j0, k0, 1)
      Flw(m)%qh(i0, j0, k0, 4) = wc*Flw(m)%qh(i0, j0, k0, 1)
      Flw(m)%qh(i0, j0, k0, 5) = Flw(m)%jac(i1, j1, k1)/Flw(m)%jac(i0, j0, k0) &
      &                     *(Flw(m)%qh(i1, j1, k1, 5) &
      &                       - 0.5*(Flw(m)%qh(i1, j1, k1, 2)**2 &
      &                               + Flw(m)%qh(i1, j1, k1, 3)**2 &
      &                               + Flw(m)%qh(i1, j1, k1, 4)**2)/Flw(m)%qh(i1, j1, k1, 1)) &
      &                     + 0.5*(Flw(m)%qh(i0, j0, k0, 2)**2 &
      &                             + Flw(m)%qh(i0, j0, k0, 3)**2 &
      &                             + Flw(m)%qh(i0, j0, k0, 4)**2)/Flw(m)%qh(i0, j0, k0, 1)
      Flw(m)%qh(i0, j0, k0, 6) = Flw(m)%qh(i1, j1, k1, 6)*Flw(m)%jac(i1, j1, k1)/Flw(m)%jac(i0, j0, k0)
      Flw(m)%qh(i0, j0, k0, 7) = Flw(m)%qh(i1, j1, k1, 7)*Flw(m)%jac(i1, j1, k1)/Flw(m)%jac(i0, j0, k0)
      ! 処理終了 ********************************************************************************************
      return
   end subroutine WallNoSlipWF
!*******************************************************************************************************
!******** 壁境界 (滑りあり, 断熱, 壁関数)                                                        ********
!*******************************************************************************************************
   subroutine WallSlipWF( &
   &            m, i0, j0, k0, i1, j1, k1, i2, j2, k2)
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
      ! 壁近傍の分布 ----------------------------------------------------------------------------------------
      dx2 = Flw(m)%x(i2, j2, k2) - Flw(m)%x(i0, j0, k0)
      dy2 = Flw(m)%y(i2, j2, k2) - Flw(m)%y(i0, j0, k0)
      dz2 = Flw(m)%z(i2, j2, k2) - Flw(m)%z(i0, j0, k0)
      uc2 = Flw(m)%qh(i2, j2, k2, 2)/Flw(m)%qh(i2, j2, k2, 1)
      vc2 = Flw(m)%qh(i2, j2, k2, 3)/Flw(m)%qh(i2, j2, k2, 1)
      wc2 = Flw(m)%qh(i2, j2, k2, 4)/Flw(m)%qh(i2, j2, k2, 1)
      ! 壁関数 ----------------------------------------------------------------------------------------------
      yp = sqrt(dx2**2 + dy2**2 + dz2**2)
      up = sqrt(uc2**2 + vc2**2 + wc2**2)
      nup = Flw(m)%mu(i2, j2, k2)/(Flw(m)%qh(i2, j2, k2, 1)*Flw(m)%jac(i2, j2, k2))
      kp = Flw(m)%kin(i2, j2, k2)
      epsp = Flw(m)%eps(i2, j2, k2)
      up = max(zero, up)
      nup = max(zero, nup)
      kp = max(zero, kp)
      epsp = max(zero, epsp)
      select case (Flw(m)%fRough(i0, k0))
      case (0)
         call WallFunctionKEM2S( &
         &      yp, up, nup, utau, dudy1, dudy2, kp, epsp)
      case (1)
         call WallFunctionKEM2R( &
         &      Flw(m)%RH(i0, k0), yp, up, nup, utau, dudy1, dudy2, kp, epsp)
      case default
         write (*, '(a)') '!!!!! Error : fRough number !!!!!'
      end select
      ! 速度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 壁面から一点目 --------------------------------------------------------------------------------------
      uc1 = Flw(m)%qh(i1, j1, k1, 2)/Flw(m)%qh(i1, j1, k1, 1)
      vc1 = Flw(m)%qh(i1, j1, k1, 3)/Flw(m)%qh(i1, j1, k1, 1)
      wc1 = Flw(m)%qh(i1, j1, k1, 4)/Flw(m)%qh(i1, j1, k1, 1)
      yp = sqrt((Flw(m)%x(i2, j2, k2) - Flw(m)%x(i1, j1, k1))**2 &
      &         + (Flw(m)%y(i2, j2, k2) - Flw(m)%y(i1, j1, k1))**2 &
      &         + (Flw(m)%z(i2, j2, k2) - Flw(m)%z(i1, j1, k1))**2)
      up = max(zero, up - dudy1*yp)/max(zero, sqrt(uc1**2 + vc1**2 + wc1**2))
      uc1 = uc1*up
      vc1 = vc1*up
      wc1 = wc1*up
      ! 壁面 ------------------------------------------------------------------------------------------------
      uc0 = 0.0
      vc0 = 0.0
      wc0 = 0.0
      ! 流束関数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 壁面から二点目 --------------------------------------------------------------------------------------
      ! kp, epsp 外挿
      Flw(m)%qh(i2, j2, k2, 5) = Flw(m)%qh(i2, j2, k2, 5) - Flw(m)%qh(i2, j2, k2, 6)
      Flw(m)%qh(i2, j2, k2, 6) = Flw(m)%qh(i2, j2, k2, 1)*kp
      Flw(m)%qh(i2, j2, k2, 7) = Flw(m)%qh(i2, j2, k2, 1)*epsp
      Flw(m)%qh(i2, j2, k2, 5) = Flw(m)%qh(i2, j2, k2, 5) + Flw(m)%qh(i2, j2, k2, 6)
      ! 壁面から一点目 --------------------------------------------------------------------------------------
      ! up, vp, wp, kp, epsp 外挿
      Flw(m)%qh(i1, j1, k1, 5) = Flw(m)%qh(i1, j1, k1, 5) &
      &                    - 0.5*(Flw(m)%qh(i1, j1, k1, 2)**2 &
      &                            + Flw(m)%qh(i1, j1, k1, 3)**2 &
      &                            + Flw(m)%qh(i1, j1, k1, 4)**2)/Flw(m)%qh(i1, j1, k1, 1) &
      &                    - Flw(m)%qh(i1, j1, k1, 6)
      Flw(m)%qh(i1, j1, k1, 2) = uc1*Flw(m)%qh(i1, j1, k1, 1)
      Flw(m)%qh(i1, j1, k1, 3) = vc1*Flw(m)%qh(i1, j1, k1, 1)
      Flw(m)%qh(i1, j1, k1, 4) = wc1*Flw(m)%qh(i1, j1, k1, 1)
      Flw(m)%qh(i1, j1, k1, 6) = kp*Flw(m)%qh(i1, j1, k1, 1)
      Flw(m)%qh(i1, j1, k1, 7) = epsp*Flw(m)%qh(i1, j1, k1, 1)
      Flw(m)%qh(i1, j1, k1, 5) = Flw(m)%qh(i1, j1, k1, 5) &
      &                    + 0.5*(Flw(m)%qh(i1, j1, k1, 2)**2 &
      &                            + Flw(m)%qh(i1, j1, k1, 3)**2 &
      &                            + Flw(m)%qh(i1, j1, k1, 4)**2)/Flw(m)%qh(i1, j1, k1, 1) &
      &                    + Flw(m)%qh(i1, j1, k1, 6)
      ! 壁面 ------------------------------------------------------------------------------------------------
      ! 全て外挿
      Flw(m)%qh(i0, j0, k0, 1) = Flw(m)%qh(i1, j1, k1, 1)*Flw(m)%jac(i1, j1, k1)/Flw(m)%jac(i0, j0, k0)
      Flw(m)%qh(i0, j0, k0, 2) = uc0*Flw(m)%qh(i0, j0, k0, 1)
      Flw(m)%qh(i0, j0, k0, 3) = vc0*Flw(m)%qh(i0, j0, k0, 1)
      Flw(m)%qh(i0, j0, k0, 4) = wc0*Flw(m)%qh(i0, j0, k0, 1)
      Flw(m)%qh(i0, j0, k0, 5) = Flw(m)%jac(i1, j1, k1)/Flw(m)%jac(i0, j0, k0) &
      &                     *(Flw(m)%qh(i1, j1, k1, 5) &
      &                       - 0.5*(Flw(m)%qh(i1, j1, k1, 2)**2 &
      &                               + Flw(m)%qh(i1, j1, k1, 3)**2 &
      &                               + Flw(m)%qh(i1, j1, k1, 4)**2)/Flw(m)%qh(i1, j1, k1, 1)) &
      &                     + 0.5*(Flw(m)%qh(i0, j0, k0, 2)**2 &
      &                             + Flw(m)%qh(i0, j0, k0, 3)**2 &
      &                             + Flw(m)%qh(i0, j0, k0, 4)**2)/Flw(m)%qh(i0, j0, k0, 1)
      Flw(m)%qh(i0, j0, k0, 6) = Flw(m)%qh(i1, j1, k1, 6)*Flw(m)%jac(i1, j1, k1)/Flw(m)%jac(i0, j0, k0)
      Flw(m)%qh(i0, j0, k0, 7) = Flw(m)%qh(i1, j1, k1, 7)*Flw(m)%jac(i1, j1, k1)/Flw(m)%jac(i0, j0, k0)
      ! 処理終了 ********************************************************************************************
      return
   end subroutine WallslipWF
!*******************************************************************************************************
!******** 壁境界 (滑りなし, 断熱, 低レイノルズ数型AKNモデル)                                          ********
!*******************************************************************************************************
   subroutine WallNoslipLRe( &
   &            m, i0, j0, k0, i1, j1, k1)
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
      dx = Flw(m)%x(i1, j1, k1) - Flw(m)%x(i0, j0, k0)
      dy = Flw(m)%y(i1, j1, k1) - Flw(m)%y(i0, j0, k0)
      dz = Flw(m)%z(i1, j1, k1) - Flw(m)%z(i0, j0, k0)
      yp = sqrt(dx**2 + dy**2 + dz**2)
      epsp = 2.0d0*dble(Flw(m)%mu(i1, j1, k1)*Flw(m)%kin(i1, j1, k1))/ &
      &      ((dble(Flw(m)%qh(i1, j1, k1, 1))*dble(Flw(m)%jac(i1, j1, k1)))*dble(yp)**2.0)
!write(*,*) i1,j1,k1,epsp,Flw(m)%kin(i1,j1,k1),Flw(m)%qh(i1,j1,k1,1),Flw(m)%jac(i1,j1,k1)
      if (epsp .ne. epsp) then
         write (*, *) 'eps_p is diverged'
         write (*, *) i1, j1, Flw(m)%mu(i1, j1, k1), Flw(m)%qh(i1, j1, k1, 6), Flw(m)%jac(i1, j1, k1), &
         &          Flw(m)%qh(i1, j1, k1, 1), Flw(m)%jac(i1, j1, k1), yp
         stop
      end if
      ! epsp
      Flw(m)%qh(i1, j1, k1, 7) = real(dble(Flw(m)%qh(i1, j1, k1, 1))*epsp)
      ! 壁面 ------------------------------------------------------------------------------------------------
      ! 全て外挿
      Flw(m)%qh(i0, j0, k0, 1) = Flw(m)%qh(i1, j1, k1, 1)*Flw(m)%jac(i1, j1, k1)/Flw(m)%jac(i0, j0, k0)
      Flw(m)%qh(i0, j0, k0, 2) = uc*Flw(m)%qh(i0, j0, k0, 1)
      Flw(m)%qh(i0, j0, k0, 3) = vc*Flw(m)%qh(i0, j0, k0, 1)
      Flw(m)%qh(i0, j0, k0, 4) = wc*Flw(m)%qh(i0, j0, k0, 1)
      Flw(m)%qh(i0, j0, k0, 5) = Flw(m)%jac(i1, j1, k1)/Flw(m)%jac(i0, j0, k0) &
      &                     *(Flw(m)%qh(i1, j1, k1, 5) &
      &                       - 0.5*(Flw(m)%qh(i1, j1, k1, 2)**2 &
      &                               + Flw(m)%qh(i1, j1, k1, 3)**2 &
      &                               + Flw(m)%qh(i1, j1, k1, 4)**2)/Flw(m)%qh(i1, j1, k1, 1)) &
      &                     + 0.5*(Flw(m)%qh(i0, j0, k0, 2)**2 &
      &                             + Flw(m)%qh(i0, j0, k0, 3)**2 &
      &                             + Flw(m)%qh(i0, j0, k0, 4)**2)/Flw(m)%qh(i0, j0, k0, 1)
      Flw(m)%qh(i0, j0, k0, 6) = 0.0
      Flw(m)%qh(i0, j0, k0, 7) = Flw(m)%qh(i1, j1, k1, 7)*Flw(m)%jac(i1, j1, k1)/Flw(m)%jac(i0, j0, k0)
      ! 処理終了 ********************************************************************************************
      return
   end subroutine WallNoSlipLRe
!*******************************************************************************************************
!******** 氷の中の点の処理 (氷層の平均値を外挿)                                                         *******
!*******************************************************************************************************
   subroutine BoundaryIceIn( &
   &            m1, m2, is, ie, js, je, ks, ke, i1, i2, j0, k0)
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, intent(in) :: m1, m2
      integer, intent(in) :: is, ie, js, je, ks, ke
      integer, intent(in) :: i1, i2, j0, k0
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, pointer :: qAve(:)
      integer :: i, j, k, l
      ! 処理開始 ********************************************************************************************
      ! メモリ確保 ------------------------------------------------------------------------------------------
      allocate (qAve(ls:le))
      ! 氷層の平均値 ----------------------------------------------------------------------------------------
      do l = ls, le
         qAve(l) = sum(Flw(m2)%qh(i1:i2, j0, k0, l)*Flw(m2)%jac(i1:i2, j0, k0)) &
         &       /real(i2 - i1 + 1)
      end do
      ! 平均値を外挿 ----------------------------------------------------------------------------------------
      do k = ks, ke
      do j = js, je
      do i = is, ie
         if (IceIn(i, j, k) == 0) cycle
         do l = ls, le
            Flw(m1)%qh(i, j, k, l) = qAve(l)/Flw(m1)%jac(i, j, k)
         end do
      end do
      end do
      end do
      ! メモリ解放 ------------------------------------------------------------------------------------------
      deallocate (qAve)
      ! 処理終了 ********************************************************************************************
      return
   end subroutine BoundaryIceIn
!*******************************************************************************************************
!******** 周期境界                                                                                  *******
!*******************************************************************************************************
   subroutine BoundaryPeriodic( &
   &            m, is, ie, js, je, k0, k1)
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
         Flw(m)%qh(i, j, k0, l) = Flw(m)%qh(i, j, k1, l)*Flw(m)%jac(i, j, k1)/Flw(m)%jac(i, j, k0)
      end do
      end do
      end do
      ! 処理終了 ********************************************************************************************
      return
   end subroutine BoundaryPeriodic
!*******************************************************************************************************
!******** 重合格子の補間                                                                        ********
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
         select case (m)
         case (1); n = 2
         case (2); n = 1
         end select
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k, l, i0, j0, k0, &
!$OMP&  i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4, &
!$OMP&  i5, j5, k5, i6, j6, k6, i7, j7, k7, i8, j8, k8)
         do k0 = OSG(m)%ks, OSG(m)%ke
         do j0 = OSG(m)%js, OSG(m)%je
         do i0 = OSG(m)%is, OSG(m)%ie
            i = OSG(m)%ip(i0, j0, k0); j = OSG(m)%jp(i0, j0, k0); k = OSG(m)%kp(i0, j0, k0)
            select case (OSG(m)%fOver(i0, j0, k0))
               ! 補間しない点 -------------------------------------------------------------------------------------
            case (0)
               cycle
               ! 三重線形補間点 -----------------------------------------------------------------------------------
            case (1)
               i1 = i; j1 = j; k1 = k; i2 = i + 1; j2 = j; k2 = k
               i3 = i; j3 = j + 1; k3 = k; i4 = i + 1; j4 = j + 1; k4 = k
               i5 = i; j5 = j; k5 = k + 1; i6 = i + 1; j6 = j; k6 = k + 1
               i7 = i; j7 = j + 1; k7 = k + 1; i8 = i + 1; j8 = j + 1; k8 = k + 1
               do l = ls, le
                  Flw(m)%qh(i0, j0, k0, l) = (OSG(m)%term1(i0, j0, k0)*Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1) &
                  &                       + OSG(m)%term2(i0, j0, k0)*Flw(n)%qh(i2, j2, k2, l)*Flw(n)%jac(i2, j2, k2) &
                  &                       + OSG(m)%term3(i0, j0, k0)*Flw(n)%qh(i3, j3, k3, l)*Flw(n)%jac(i3, j3, k3) &
                  &                       + OSG(m)%term4(i0, j0, k0)*Flw(n)%qh(i4, j4, k4, l)*Flw(n)%jac(i4, j4, k4) &
                  &                       + OSG(m)%term5(i0, j0, k0)*Flw(n)%qh(i5, j5, k5, l)*Flw(n)%jac(i5, j5, k5) &
                  &                       + OSG(m)%term6(i0, j0, k0)*Flw(n)%qh(i6, j6, k6, l)*Flw(n)%jac(i6, j6, k6) &
                  &                       + OSG(m)%term7(i0, j0, k0)*Flw(n)%qh(i7, j7, k7, l)*Flw(n)%jac(i7, j7, k7) &
                  &                       + OSG(m)%term8(i0, j0, k0)*Flw(n)%qh(i8, j8, k8, l)*Flw(n)%jac(i8, j8, k8) &
                  &                         )/Flw(m)%jac(i0, j0, k0)
               end do
               ! 三次元線形補間点 ---------------------------------------------------------------------------------
            case (2)
               i1 = i; j1 = j; k1 = k; i2 = i; j2 = j; k2 = k + 1
               i3 = i + 1; j3 = j; k3 = k + 1; i4 = i; j4 = j + 1; k4 = k + 1
               do l = ls, le
                  Flw(m)%qh(i0, j0, k0, l) = (OSG(m)%term1(i0, j0, k0)*Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1) &
                  &                       + OSG(m)%term2(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i2, j2, k2, l)*Flw(n)%jac(i2, j2, k2) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                       + OSG(m)%term3(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i3, j3, k3, l)*Flw(n)%jac(i3, j3, k3) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                       + OSG(m)%term4(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i4, j4, k4, l)*Flw(n)%jac(i4, j4, k4) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                         )/Flw(m)%jac(i0, j0, k0)
               end do
            case (3)
               i1 = i; j1 = j; k1 = k; i2 = i; j2 = j + 1; k2 = k
               i3 = i + 1; j3 = j + 1; k3 = k; i4 = i; j4 = j + 1; k4 = k + 1
               do l = ls, le
                  Flw(m)%qh(i0, j0, k0, l) = (OSG(m)%term1(i0, j0, k0)*Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1) &
                  &                       + OSG(m)%term2(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i2, j2, k2, l)*Flw(n)%jac(i2, j2, k2) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                       + OSG(m)%term3(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i3, j3, k3, l)*Flw(n)%jac(i3, j3, k3) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                       + OSG(m)%term4(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i4, j4, k4, l)*Flw(n)%jac(i4, j4, k4) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                         )/Flw(m)%jac(i0, j0, k0)
               end do
            case (4)
               i1 = i; j1 = j; k1 = k; i2 = i + 1; j2 = j + 1; k2 = k
               i3 = i + 1; j3 = j; k3 = k + 1; i4 = i; j4 = j + 1; k4 = k + 1
               do l = ls, le
                  Flw(m)%qh(i0, j0, k0, l) = (OSG(m)%term1(i0, j0, k0)*Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1) &
                  &                       + OSG(m)%term2(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i2, j2, k2, l)*Flw(n)%jac(i2, j2, k2) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                       + OSG(m)%term3(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i3, j3, k3, l)*Flw(n)%jac(i3, j3, k3) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                       + OSG(m)%term4(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i4, j4, k4, l)*Flw(n)%jac(i4, j4, k4) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                         )/Flw(m)%jac(i0, j0, k0)
               end do
            case (5)
               i1 = i; j1 = j; k1 = k; i2 = i + 1; j2 = j; k2 = k
               i3 = i + 1; j3 = j + 1; k3 = k; i4 = i + 1; j4 = j; k4 = k + 1
               do l = ls, le
                  Flw(m)%qh(i0, j0, k0, l) = (OSG(m)%term1(i0, j0, k0)*Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1) &
                  &                       + OSG(m)%term2(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i2, j2, k2, l)*Flw(n)%jac(i2, j2, k2) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                       + OSG(m)%term3(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i3, j3, k3, l)*Flw(n)%jac(i3, j3, k3) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                       + OSG(m)%term4(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i4, j4, k4, l)*Flw(n)%jac(i4, j4, k4) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                         )/Flw(m)%jac(i0, j0, k0)
               end do
            case (6)
               i1 = i + 1; j1 = j + 1; k1 = k; i2 = i + 1; j2 = j; k2 = k + 1
               i3 = i; j3 = j + 1; k3 = k + 1; i4 = i + 1; j4 = j + 1; k4 = k + 1
               do l = ls, le
                  Flw(m)%qh(i0, j0, k0, l) = (OSG(m)%term1(i0, j0, k0)*Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1) &
                  &                       + OSG(m)%term2(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i2, j2, k2, l)*Flw(n)%jac(i2, j2, k2) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                       + OSG(m)%term3(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i3, j3, k3, l)*Flw(n)%jac(i3, j3, k3) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                       + OSG(m)%term4(i0, j0, k0) &
                  &                         *(Flw(n)%qh(i4, j4, k4, l)*Flw(n)%jac(i4, j4, k4) &
                  &                           - Flw(n)%qh(i1, j1, k1, l)*Flw(n)%jac(i1, j1, k1)) &
                  &                         )/Flw(m)%jac(i0, j0, k0)
               end do
            end select
         end do
         end do
         end do
         !$OMP END PARALLEL DO
      end do
      ! 処理終了 ********************************************************************************************
      return
   end subroutine InterpolationOversetGrid
!*******************************************************************************************************
!******** 残差計算                                                                                 ********
!*******************************************************************************************************
   subroutine CalResidual
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, pointer :: res0(:, :, :, :)
      integer :: l, m
      ! 処理開始 ********************************************************************************************
! do m = ms, me
      m = me
      ! メモリ確保 -----------------------------------------------------------------------------------------
      allocate (res0(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke, ls:le))
      ! 残差計算 -------------------------------------------------------------------------------------------
      call CalResidual3D( &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
      &     ls, le, Flw(m)%dt, Flw(m)%qh0, Flw(m)%qh, &
      &     res0)
      ! 最大値算出 -----------------------------------------------------------------------------------------
      do l = ls, le
         Flw(m)%Res(l) = maxval(res0(:, :, :, l))
      end do
      ! メモリ解放 -----------------------------------------------------------------------------------------
      deallocate (res0)
! enddo
      ! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! write(20, '(i6, 7e16.8e3)') nCount, Flw(1)%Res(ls:le)
      write (21, '(i6, 7e16.8e3)') nCount, Flw(2)%Res(ls:le)
      ! 処理終了 ********************************************************************************************
      return
   end subroutine CalResidual
!*******************************************************************************************************
!******** 流れ場ファイル出力 (計算回数毎)                                                         ********
!*******************************************************************************************************
   subroutine OutputFileCount
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      character :: fname*30
      integer   :: m
      ! 処理開始 ********************************************************************************************
      ! 計算条件 --------------------------------------------------------------------------------------------
      call OUtput_CalSetting(trim(ND_CalSetFile)//strtxt)
      ! 流束関数 --------------------------------------------------------------------------------------------
      write (fname, '(a, 2(i2.2,a))') 'IceStep', IceStep, 'of', IceStepMax, '_'
! do m = ms, me
      m = me
      call Output_Flux3D( &
      &      trim(FlwCalOutDir)//trim(BlkName(m))//trim(ND_FlxFile), strbin, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
      &      ls, le, Flw(m)%qh)
      call Output_Flux3D( &
      &      trim(FlwCalOutDir)//trim(fname)//trim(BlkName(m))//trim(ND_FlxFile), strbin, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
      &      ls, le, Flw(m)%qh)
! enddo
      ! 処理終了 ********************************************************************************************
      return
   end subroutine OutputFileCount
!*******************************************************************************************************
!******** メモリ解放                                                                                 ********
!*******************************************************************************************************
   subroutine Deallocating
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer :: m
      ! 処理開始 ********************************************************************************************
! do m = ms, me
      m = me
      deallocate (Flw(m)%x, Flw(m)%y, Flw(m)%z, &
      &           Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, &
      &           Flw(m)%p, Flw(m)%t, Flw(m)%mu, &
      &           Flw(m)%kin, Flw(m)%eps, Flw(m)%mut, &
      &           Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
      &           Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
      &           Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
      &           Flw(m)%jac, Flw(m)%qh, Flw(m)%qh0, &
      &           Flw(m)%dqc, Flw(m)%dqd, Flw(m)%dqp, Flw(m)%dqr, &
      &           Flw(m)%dqh, Flw(m)%dqh0, &
      &           Flw(m)%dt, Flw(m)%Res)
!  deallocate( OSG(m)%ip, OSG(m)%jp, OSG(m)%kp, &
!  &           OSG(m)%term1, OSG(m)%term2, OSG(m)%term3, OSG(m)%term4, &
!  &           OSG(m)%term5, OSG(m)%term6, OSG(m)%term7, OSG(m)%term8, &
!  &           OSG(m)%fOver )
! enddo
      ! 処理終了 ********************************************************************************************
      return
   end subroutine Deallocating
!*******************************************************************************************************
!******** y+計算                                                                                 ********
!*******************************************************************************************************
   subroutine calyplus(m)
      implicit none
      !mainroutine_variable
      integer, intent(in)        :: m
      !subroutine_variable
      integer                :: i, j, k
      real, allocatable        :: utau(:)
      real, allocatable        :: yplus(:)
      real                        :: yp, up
      real                        :: ypmax, ypmin

      !array_allocation
      allocate (utau(Flw(m)%is:Flw(m)%ie))
      allocate (yplus(Flw(m)%i1:Flw(m)%i3))

      do i = Flw(m)%i1, Flw(m)%i3
         yp = sqrt((Flw(m)%x(i, Flw(m)%js + 1, Flw(m)%ks) - Flw(m)%x(i, Flw(m)%js, Flw(m)%ks))**2.0 + &
              &    (Flw(m)%y(i, Flw(m)%js + 1, Flw(m)%ks) - Flw(m)%y(i, Flw(m)%js, Flw(m)%ks))**2.0)
         up = sqrt(Flw(m)%u(i, Flw(m)%js + 1, Flw(m)%ks)**2.0 + Flw(m)%v(i, Flw(m)%js + 1, Flw(m)%ks)**2.0)
         yplus(i) = sqrt(Flw(m)%rho(i, Flw(m)%js, Flw(m)%ks)/ &
         &              (Flw(m)%mut(i, Flw(m)%js, Flw(m)%ks) + Flw(m)%mu(i, Flw(m)%js, Flw(m)%ks))* &
         &              yp*up)
      end do
      ypmax = maxval(yplus)
      ypmin = minval(yplus)

      write (*, *) 'Range of yplus:', ypmin, 'to', ypmax

      deallocate (utau)
      deallocate (yplus)
   end subroutine calyplus

!*******************************************************************************************************
!******** ParaViewファイル出力                                                                        ********
!*******************************************************************************************************
   subroutine OutputPara_bin( &
   &      strdir, strname, is, ie, js, je, ks, ke, &
   &      rhoRef, aRef, lRef, rho, u, v, w, Ps, Ts, mu, &
   &      kin, eps, mut, &
   &      x, y, z)
      implicit none
      !mainroutine_variable
      character, intent(in)  :: strdir*(*), strname*(*)
      integer, intent(in)  :: is, ie, js, je, ks, ke
      real, intent(in)  :: rhoRef, aRef, lRef
      real, intent(in)  :: rho(is:ie, js:je, ks:ke), &
      &                         u(is:ie, js:je, ks:ke), v(is:ie, js:je, ks:ke), w(is:ie, js:je, ks:ke), &
      &                         Ps(is:ie, js:je, ks:ke), Ts(is:ie, js:je, ks:ke), mu(is:ie, js:je, ks:ke), &
      &                         kin(is:ie, js:je, ks:ke), eps(is:ie, js:je, ks:ke), &
      &                         mut(is:ie, js:je, ks:ke)
      real, intent(in)  :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke), z(is:ie, js:je, ks:ke)
      !subroutine_variable
      character(len=4), parameter                :: strvtk = '.vtk'
      character(len=1), parameter                :: newline = char(10)
      character(len=200)        :: strnum
      integer                :: i, j, k, n
      integer                :: ni, nj, nk
      integer                :: npoint

      npoint = (ie - is)*(je - js)*(ke - ks)

open(unit       = 1, &
     file       = trim(strdir)//trim(strname)//trim(strvtk), &
     form       = 'unformatted', &
     access     = 'stream', &          ! ← ここを stream に
     convert    = 'big_endian', &      ! GNU拡張（OK）。気になるなら -fconvert で全体指定でも可
     action     = 'write')
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
         write (1) x(ni, nj, nk)*lRef, y(ni, nj, nk)*lRef, z(ni, nj, nk)*lRef
         write (1) x(ni + 1, nj, nk)*lRef, y(ni + 1, nj, nk)*lRef, z(ni + 1, nj, nk)*lRef
         write (1) x(ni + 1, nj + 1, nk)*lRef, y(ni + 1, nj + 1, nk)*lRef, z(ni + 1, nj + 1, nk)*lRef
         write (1) x(ni, nj + 1, nk)*lRef, y(ni, nj + 1, nk)*lRef, z(ni, nj + 1, nk)*lRef
         write (1) x(ni, nj, nk + 1)*lRef, y(ni, nj, nk + 1)*lRef, z(ni, nj, nk + 1)*lRef
         write (1) x(ni + 1, nj, nk + 1)*lRef, y(ni + 1, nj, nk + 1)*lRef, z(ni + 1, nj, nk + 1)*lRef
         write (1) x(ni + 1, nj + 1, nk + 1)*lRef, y(ni + 1, nj + 1, nk + 1)*lRef, z(ni + 1, nj + 1, nk + 1)*lRef
         write (1) x(ni, nj + 1, nk + 1)*lRef, y(ni, nj + 1, nk + 1)*lRef, z(ni, nj + 1, nk + 1)*lRef
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
      write (1) 'VECTORS velocity float'//newline
      do n = 0, npoint - 1
         ni = mod(mod(n, (ie - is)*(je - js)), (ie - is))
         nj = int(mod(n, (ie - is)*(je - js))/(ie - is))
         nk = int(n/((ie - is)*(je - js)))
         write (1) u(ni, nj, nk)*aRef, v(ni, nj, nk)*aRef, w(ni, nj, nk)*aRef
         write (1) u(ni + 1, nj, nk)*aRef, v(ni + 1, nj, nk)*aRef, w(ni + 1, nj, nk)*aRef
         write (1) u(ni + 1, nj + 1, nk)*aRef, v(ni + 1, nj + 1, nk)*aRef, w(ni + 1, nj + 1, nk)*aRef
         write (1) u(ni, nj + 1, nk)*aRef, v(ni, nj + 1, nk)*aRef, w(ni, nj + 1, nk)*aRef
         write (1) u(ni, nj, nk + 1)*aRef, v(ni, nj, nk + 1)*aRef, w(ni, nj, nk + 1)*aRef
         write (1) u(ni + 1, nj, nk + 1)*aRef, v(ni + 1, nj, nk + 1)*aRef, w(ni + 1, nj, nk + 1)*aRef
         write (1) u(ni + 1, nj + 1, nk + 1)*aRef, v(ni + 1, nj + 1, nk + 1)*aRef, w(ni + 1, nj + 1, nk + 1)*aRef
         write (1) u(ni, nj + 1, nk + 1)*aRef, v(ni, nj + 1, nk + 1)*aRef, w(ni, nj + 1, nk + 1)*aRef
      end do
      write (1) newline
      write (1) 'SCALARS rho float'//newline
      write (1) 'LOOKUP_TABLE default'//newline
      do n = 0, npoint - 1
         ni = mod(mod(n, (ie - is)*(je - js)), (ie - is))
         nj = int(mod(n, (ie - is)*(je - js))/(ie - is))
         nk = int(n/((ie - is)*(je - js)))
         write (1) rho(ni, nj, nk)*rhoRef
         write (1) rho(ni + 1, nj, nk)*rhoRef
         write (1) rho(ni + 1, nj + 1, nk)*rhoRef
         write (1) rho(ni, nj + 1, nk)*rhoRef
         write (1) rho(ni, nj, nk + 1)*rhoRef
         write (1) rho(ni + 1, nj, nk + 1)*rhoRef
         write (1) rho(ni + 1, nj + 1, nk + 1)*rhoRef
         write (1) rho(ni, nj + 1, nk + 1)*rhoRef
      end do
      write (1) newline
      write (1) 'SCALARS Ps float'//newline
      write (1) 'LOOKUP_TABLE default'//newline
      do n = 0, npoint - 1
         ni = mod(mod(n, (ie - is)*(je - js)), (ie - is))
         nj = int(mod(n, (ie - is)*(je - js))/(ie - is))
         nk = int(n/((ie - is)*(je - js)))
         write (1) Ps(ni, nj, nk)*rhoRef*aRef**2.0
         write (1) Ps(ni + 1, nj, nk)*rhoRef*aRef**2.0
         write (1) Ps(ni + 1, nj + 1, nk)*rhoRef*aRef**2.0
         write (1) Ps(ni, nj + 1, nk)*rhoRef*aRef**2.0
         write (1) Ps(ni, nj, nk + 1)*rhoRef*aRef**2.0
         write (1) Ps(ni + 1, nj, nk + 1)*rhoRef*aRef**2.0
         write (1) Ps(ni + 1, nj + 1, nk + 1)*rhoRef*aRef**2.0
         write (1) Ps(ni, nj + 1, nk + 1)*rhoRef*aRef**2.0
      end do
      write (1) newline
      write (1) 'SCALARS Ts float'//newline
      write (1) 'LOOKUP_TABLE default'//newline
      do n = 0, npoint - 1
         ni = mod(mod(n, (ie - is)*(je - js)), (ie - is))
         nj = int(mod(n, (ie - is)*(je - js))/(ie - is))
         nk = int(n/((ie - is)*(je - js)))
         write (1) Ts(ni, nj, nk)*aRef**2.0
         write (1) Ts(ni + 1, nj, nk)*aRef**2.0
         write (1) Ts(ni + 1, nj + 1, nk)*aRef**2.0
         write (1) Ts(ni, nj + 1, nk)*aRef**2.0
         write (1) Ts(ni, nj, nk + 1)*aRef**2.0
         write (1) Ts(ni + 1, nj, nk + 1)*aRef**2.0
         write (1) Ts(ni + 1, nj + 1, nk + 1)*aRef**2.0
         write (1) Ts(ni, nj + 1, nk + 1)*aRef**2.0
      end do
      write (1) newline
      write (1) 'SCALARS mu float'//newline
      write (1) 'LOOKUP_TABLE default'//newline
      do n = 0, npoint - 1
         ni = mod(mod(n, (ie - is)*(je - js)), (ie - is))
         nj = int(mod(n, (ie - is)*(je - js))/(ie - is))
         nk = int(n/((ie - is)*(je - js)))
         write (1) mu(ni, nj, nk)*rhoRef*aRef*lRef
         write (1) mu(ni + 1, nj, nk)*rhoRef*aRef*lRef
         write (1) mu(ni + 1, nj + 1, nk)*rhoRef*aRef*lRef
         write (1) mu(ni, nj + 1, nk)*rhoRef*aRef*lRef
         write (1) mu(ni, nj, nk + 1)*rhoRef*aRef*lRef
         write (1) mu(ni + 1, nj, nk + 1)*rhoRef*aRef*lRef
         write (1) mu(ni + 1, nj + 1, nk + 1)*rhoRef*aRef*lRef
         write (1) mu(ni, nj + 1, nk + 1)*rhoRef*aRef*lRef
      end do
      write (1) newline
      write (1) 'SCALARS k float'//newline
      write (1) 'LOOKUP_TABLE default'//newline
      do n = 0, npoint - 1
         ni = mod(mod(n, (ie - is)*(je - js)), (ie - is))
         nj = int(mod(n, (ie - is)*(je - js))/(ie - is))
         nk = int(n/((ie - is)*(je - js)))
         write (1) kin(ni, nj, nk)*aRef**2.0
         write (1) kin(ni + 1, nj, nk)*aRef**2.0
         write (1) kin(ni + 1, nj + 1, nk)*aRef**2.0
         write (1) kin(ni, nj + 1, nk)*aRef**2.0
         write (1) kin(ni, nj, nk + 1)*aRef**2.0
         write (1) kin(ni + 1, nj, nk + 1)*aRef**2.0
         write (1) kin(ni + 1, nj + 1, nk + 1)*aRef**2.0
         write (1) kin(ni, nj + 1, nk + 1)*aRef**2.0
      end do
      write (1) newline
      write (1) 'SCALARS eps float'//newline
      write (1) 'LOOKUP_TABLE default'//newline
      do n = 0, npoint - 1
         ni = mod(mod(n, (ie - is)*(je - js)), (ie - is))
         nj = int(mod(n, (ie - is)*(je - js))/(ie - is))
         nk = int(n/((ie - is)*(je - js)))
         write (1) eps(ni, nj, nk)*aRef**3.0/lRef
         write (1) eps(ni + 1, nj, nk)*aRef**3.0/lRef
         write (1) eps(ni + 1, nj + 1, nk)*aRef**3.0/lRef
         write (1) eps(ni, nj + 1, nk)*aRef**3.0/lRef
         write (1) eps(ni, nj, nk + 1)*aRef**3.0/lRef
         write (1) eps(ni + 1, nj, nk + 1)*aRef**3.0/lRef
         write (1) eps(ni + 1, nj + 1, nk + 1)*aRef**3.0/lRef
         write (1) eps(ni, nj + 1, nk + 1)*aRef**3.0/lRef
      end do
      write (1) newline
      write (1) 'SCALARS mut float'//newline
      write (1) 'LOOKUP_TABLE default'//newline
      do n = 0, npoint - 1
         ni = mod(mod(n, (ie - is)*(je - js)), (ie - is))
         nj = int(mod(n, (ie - is)*(je - js))/(ie - is))
         nk = int(n/((ie - is)*(je - js)))
         write (1) mut(ni, nj, nk)*rhoRef*aRef*lRef
         write (1) mut(ni + 1, nj, nk)*rhoRef*aRef*lRef
         write (1) mut(ni + 1, nj + 1, nk)*rhoRef*aRef*lRef
         write (1) mut(ni, nj + 1, nk)*rhoRef*aRef*lRef
         write (1) mut(ni, nj, nk + 1)*rhoRef*aRef*lRef
         write (1) mut(ni + 1, nj, nk + 1)*rhoRef*aRef*lRef
         write (1) mut(ni + 1, nj + 1, nk + 1)*rhoRef*aRef*lRef
         write (1) mut(ni, nj + 1, nk + 1)*rhoRef*aRef*lRef
      end do
      close (1)

   end subroutine OutputPara_bin
! 定義終了 *********************************************************************************************
end program FlowField_NACA
