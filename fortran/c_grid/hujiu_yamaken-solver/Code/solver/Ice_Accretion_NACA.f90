!*******************************************************************************************************
!*******************************************************************************************************
!******** 着氷計算プログラム                                                                        ********
!******** (NACA翼，三次元圧縮性乱流場，重合格子法)                                                ********
!******** 有次元                                                                                ********
!********                                              2013.07.06  PROGRAMMED BY RYOSUKE HAYASHI ********
!********                                              2013.07.18     UPDATED BY RYOSUKE HAYASHI ********
!********                                              2014.04.15     UPDATED BY RYOSUKE HAYASHI ********
!********                                              2015.06.11     UPDATED BY MIKI    SHIMURA ********
!*******************************************************************************************************
!*******************************************************************************************************
program IceAccretion_NACA
   ! モジュール宣言 **************************************************************************************
   use Package_NACA
   use Package_FileIO
   use Package_Flow
   use Package_Icing
   ! 変数宣言 ********************************************************************************************
   implicit none
   ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer, parameter :: mRef = 2
   ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer :: kRef
   ! 処理開始 ********************************************************************************************
   ! 検証実験ケース選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   write (*, '(a)') "<< Exp. Case Selection >>"
   call SelectExpCase
   ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   write (*, '(/,a)') "<< Initial Setting >>"
   call InitialSetting
   ! 着氷計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   write (*, '(/,a)') "<< Icing Computation >>"
   call CalIceAccretion
   ! 着氷翼座標 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   write (*, '(/,a)') "<< Icing Blade Geometry >>"
   call CalIcingBlade
   ! ファイル出力 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   write (*, '(/,a)') "<< File Output >>"
   call OutputFile
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
      ! 計算条件ファイル入力
      call Input_CalSetting(trim(ND_CalSetFile)//strtxt)
      ! 有次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      VelIn = VelIn*aRef
      PtIn = PtIn*(rhoRef*aRef**2)
      TtIn = TtIn*aRef**2
      PsOut = PsOut*(rhoRef*aRef**2)
      TsOut = TsOut*aRef**2
      LWC = LWC*RhoRef
      MVD = MVD*LRef
      RhoD = RhoD*RhoRef
      Span = Span*LRef
      Chord = Chord*LRef
      VelExp = VelExp*aRef
      PsExp = PsExp*(rhoRef*aRef**2)
      TsExp = TsExp*aRef**2
      ! ディレクトリ設定
      if (IceStep == 0) then
         GrdInDir = bckdir//'grid//clean//'
         OSGDir = bckdir//'overset//clean//'
         FlwIniDir = bckdir//'flow//initial//clean//'
         FlwCalInDir = bckdir//'flow//cal//clean//'
         DrpImpDir = bckdir//'droplet//impingement//clean//'
      else
         GrdInDir = bckdir//'grid//icing//'
         OSGDir = bckdir//'overset//icing//'
         FlwIniDir = bckdir//'flow//initial//icing//'
         FlwCalInDir = bckdir//'flow//cal//icing//'
         DrpImpDir = bckdir//'droplet//impingement//icing//'
         IceCalInDir = bckdir//'icing//cal//'
      end if
      IceCalOutDir = bckdir//'icing//cal//'
      write (*, '(a)') '+++ Icing Step +++'
      write (*, '(a,i2)') '* Ice step      = ', IceStep
      write (*, '(a,i2)') '* Ice step max. = ', IceStepMax
      write (*, '(/,a)') '+++ Exp. Condition +++'
      write (*, '(a,e16.8e3)') '* Ts    = ', TsExp
      write (*, '(a,e16.8e3)') '* Ps    = ', PsExp
      write (*, '(a,e16.8e3)') '* V     = ', VelExp
      write (*, '(a,e16.8e3)') '* LWC   = ', LWC
      write (*, '(a,e16.8e3)') '* MVD   = ', MVD
      write (*, '(a,e16.8e3)') '* Rho   = ', Rhod
      write (*, '(a,e16.8e3)') '* Chord = ', Chord
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
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, pointer :: CE2D(:)
      integer   :: i, k, m
      integer   :: j0, j1, jb
      character :: fname*30
      ! 処理開始 ********************************************************************************************
      ! 構造体メモリ確保 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      allocate (Flw(ms:me), Ice(ms:me))
      ! ブロック名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call Set_BlockName
      ! 格子解像度 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
      m = me
      call Input_Resolution3D( &
      &      trim(GrdInDir)//trim(BlkName(m))//trim(RslFile), strtxt, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke)
! enddo
      ! C 型格子分割点 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! m = 1
      m = me
      call Input_CtypeGridPoint( &
      &      trim(GrdInDir)//trim(BlkName(m))//trim(CtypePointFile)//strtxt, &
      &      Flw(m)%i1, Flw(m)%i2, Flw(m)%i3)
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
      &         Flw(m)%qh(Flw(m)%is:Flw(m)%ie, Flw(m)%js:Flw(m)%je, Flw(m)%ks:Flw(m)%ke, ls:le))
! enddo
      ! 格子座標 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
      m = me
      call Input_Grid3D( &
      &      trim(GrdInDir)//trim(BlkName(m))//trim(GrdFile), strbin, &
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
! do m = ms, me
      m = me
      call Input_Flux3D( &
      &      trim(FlwCalInDir)//trim(BlkName(m))//trim(ND_FlxFile), strbin, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
      &      ls, le, Flw(m)%qh)
! enddo
      ! 物理量 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
      m = me
      call SetPhysics3DKEM( &
      &      Rg, gamma, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
      &      Flw(m)%qh, Flw(m)%jac, &
      &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps)
      call ViscosityCoefficient3D( &
      &      muSth, TsSth, s1, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
      &      Flw(m)%t, &
      &      Flw(m)%mu)
      Flw(m)%w(:, :, :) = 0.0
! enddo
      ! 有次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
      m = me
      call DimensionalizedPhysics3D( &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
      &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%mu, &
      &      Flw(m)%kin, Flw(m)%eps, Flw(m)%mut)
! enddo
      ! 着氷計算領域 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
      m = me
      call Input_Resolution1D( &
      &      trim(GrdInDir)//trim(BlkName(m))//trim(IceRslFile), strtxt, &
      &      Ice(m)%is, Ice(m)%ie)
! enddo
      ! 計算対象断面設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      m = mRef; kRef = Flw(m)%ks + int(0.5*(Flw(m)%ke - Flw(m)%ks))
      ! 液滴投入条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call Input_DropInCondition( &
      &      trim(DrpImpDir)//trim(DrpInConditionFile), strtxt, &
      &      DrpInArea, DrpInVel)
      ! 液滴衝突データ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
      m = me
      ! メモリ確保
      allocate (Flw(m)%Nimp(Ice(m)%is:Ice(m)%ie, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%Uimp(Ice(m)%is:Ice(m)%ie, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%Vimp(Ice(m)%is:Ice(m)%ie, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%Wimp(Ice(m)%is:Ice(m)%ie, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%Simp(Ice(m)%is:Ice(m)%ie, Flw(m)%ks:Flw(m)%ke), &
      &         Flw(m)%CE(Ice(m)%is:Ice(m)%ie, Flw(m)%ks:Flw(m)%ke))
      ! 衝突面積
      call Input_ArrayReal2D( &
      &      trim(DrpImpDir)//trim(BlkName(m))//trim(DrpImpiAreaFile), strbin, &
      &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, Flw(m)%simp)
      ! 液滴収集効率及び衝突速度
      call Input_CollectionEfficiency3D( &
      &      trim(DrpImpDir)//trim(BlkName(m))//trim(CollectionFile), strbin, &
      &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, &
      &      Flw(m)%CE, Flw(m)%Uimp, Flw(m)%Vimp, Flw(m)%Wimp)
! enddo
      ! 着氷計算初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! メモリ確保
      m = mRef; k = kRef
      allocate (Ice(m)%x(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%y(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%z(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%f(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Rho0(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%mu0(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Ts0(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Yp(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Up(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Vp(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Wp(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Velp(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%nup(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Ptp(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%SA(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%fRough(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Uimp(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Vimp(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%CE(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Ue(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Ta(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%nPhase(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Mes(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Mim(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Mrin(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Mrout(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Mst(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Mac(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Bi(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Bi0(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%dBi(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%Ti(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%xi(Ice(m)%is:Ice(m)%ie), &
      &         Ice(m)%yi(Ice(m)%is:Ice(m)%ie))
      ! 物理量
      j0 = Flw(m)%js; j1 = Flw(m)%js + 1
! jb = nint(0.3 * Flw(m)%je)
! jb = 10
      jb = 100
      do i = Ice(m)%is, Ice(m)%ie
         Ice(m)%x(i) = Flw(m)%x(i, j0, k)
         Ice(m)%y(i) = Flw(m)%y(i, j0, k)
         Ice(m)%z(i) = Flw(m)%z(i, j0, k)
         Ice(m)%Rho0(i) = Flw(m)%rho(i, j0, k)
         Ice(m)%mu0(i) = Flw(m)%mu(i, j0, k)
         Ice(m)%Ts0(i) = Flw(m)%t(i, j0, k)
         Ice(m)%Yp(i) = sqrt((Flw(m)%x(i, j1, k) - Flw(m)%x(i, j0, k))**2 &
         &                    + (Flw(m)%y(i, j1, k) - Flw(m)%y(i, j0, k))**2 &
         &                    + (Flw(m)%z(i, j1, k) - Flw(m)%z(i, j0, k))**2)
         Ice(m)%Up(i) = Flw(m)%u(i, j1, k) - Flw(m)%u(i, j0, k)
         Ice(m)%Vp(i) = Flw(m)%v(i, j1, k) - Flw(m)%v(i, j0, k)
         Ice(m)%Wp(i) = Flw(m)%w(i, j1, k) - Flw(m)%w(i, j0, k)
         Ice(m)%Velp(i) = sqrt(Ice(m)%Up(i)**2 + Ice(m)%Vp(i)**2 + Ice(m)%Wp(i)**2)
         Ice(m)%Ptp(i) = Flw(m)%p(i, j1, k) + 0.5*Flw(m)%rho(i, j1, k)*Ice(m)%Velp(i)**2
         Ice(m)%Uimp(i) = Flw(m)%Uimp(i, k)
         Ice(m)%Vimp(i) = Flw(m)%Vimp(i, k)
         Ice(m)%SA(i) = Flw(m)%Simp(i, k)
         Ice(m)%CE(i) = Flw(m)%CE(i, k)
         Ice(m)%Ue(i) = maxval(sqrt(Flw(m)%u(i, :jb, k)**2 &
         &                            + Flw(m)%v(i, :jb, k)**2 &
         &                            + Flw(m)%w(i, :jb, k)**2))                !!! 境界層外端の速度に変更 !!!
         Ice(m)%Ta(i) = Flw(m)%t(i, jb, k)                                        !!! 周囲流体の温度に変更 !!!
         if (Flw(m)%rho(i, j1, k) > 0.0) then
            Ice(m)%nup(i) = Flw(m)%mu(i, j1, k)/Flw(m)%rho(i, j1, k)
         else
            Ice(m)%nup(i) = Flw(m)%mu(i, j1, k)
            write (*, '(a)') '!!!!! Error : density is negative !!!!!'
         end if
      end do
      ! 着氷計算途中解
      if (IceStep == 0) then
         Ice(m)%fRough(:) = 0
         Ice(m)%nPhase(:) = 0
         Ice(m)%Mes(:) = 0.0
         Ice(m)%Mim(:) = 0.0
         Ice(m)%Mrin(:) = 0.0
         Ice(m)%Mrout(:) = 0.0
         Ice(m)%Mst(:) = 0.0
         Ice(m)%Mac(:) = 0.0
         Ice(m)%Bi(:) = 0.0
         Ice(m)%dBi(:) = 0.0
      else
         ! 表面粗さのフラグ
         Ice(m)%fRough(:) = 1
         ! 氷層厚さ・氷層温度
         call Input_IceThickTem2D( &
         &      trim(IceCalInDir)//trim(BlkName(m))//trim(IceThickTemFile), strdat, &
         &      Ice(m)%is, Ice(m)%ie, &
         &      Ice(m)%f, Ice(m)%Bi, Ice(m)%dBi, Ice(m)%Ti)
         ! 熱力学モデルパラメータ
         select case (ThermoNum)
         case (1)
            call Input_OrgMessingerPara2D( &
            &      trim(IceCalInDir)//trim(BlkName(m))//trim(MessingerFile), strdat, &
            &      Ice(m)%is, Ice(m)%ie, &
            &      Ice(m)%Mac, Ice(m)%Mes, Ice(m)%Mim, Ice(m)%Mrin, Ice(m)%Mrout, Ice(m)%Mst)
         case (2)
            call Input_ExtMessingerPara2D( &
            &      trim(IceCalInDir)//trim(BlkName(m))//trim(ExMessingerFile), strdat, &
            &      Ice(m)%is, Ice(m)%ie, &
            &      Ice(m)%nPhase, &
            &      Ice(m)%Mes, Ice(m)%Mim, Ice(m)%Mrin, Ice(m)%Mrout, Ice(m)%Mst)
         end select
      end if
      ! メモリ解放 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
      m = me
      deallocate (Flw(m)%x, Flw(m)%y, Flw(m)%z, Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, &
      &           Flw(m)%p, Flw(m)%t, Flw(m)%mu, Flw(m)%kin, Flw(m)%eps, &
      &           Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
      &           Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
      &           Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
      &           Flw(m)%jac, Flw(m)%qh, &
      &           Flw(m)%Nimp, Flw(m)%Uimp, Flw(m)%Vimp, Flw(m)%Wimp, Flw(m)%Simp, Flw(m)%CE)
! enddo
      ! 処理終了 ********************************************************************************************
      return
   end subroutine InitialSetting
!*******************************************************************************************************
!******** 着氷計算                                                                                ********
!*******************************************************************************************************
   subroutine CalIceAccretion
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer, parameter :: nStart = 0
      integer, parameter :: nIceLog = 1000
      real, parameter :: LmtBi = 2.0
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real, pointer :: Utau(:)                                                ! 摩擦速度
      real, pointer :: RH(:)                                                ! 表面粗さ
      real, pointer :: hc(:)                                                ! 熱伝達率 [J/(s*m^2*K)]
      real, pointer :: ff(:)                                                ! 氷結率
      real, pointer :: Rhoi(:)                                                ! 氷の密度
      real, pointer :: Q(:)                                                ! 熱 [J]
      real, pointer :: Him(:), Hes(:), Hac(:), &
      &                   Hrin(:), Hrout(:)                                        ! 比エンタルピ [J/kg]
      real, pointer :: Bw(:), Tw(:)                                        ! 水膜厚さ及び温度
      real, pointer :: Bg(:), tg(:)                                        ! 雨氷が現れる氷層厚さ及び時間
      real, pointer :: Q0(:), Q1(:)                                        ! 熱流束
      real, pointer :: IceCc(:), IceSc(:), IceAc(:)                        ! 氷に作用する力
      real, pointer :: B0(:)
      integer   :: i, m, n, nn
      integer   :: nCountMax
      real      :: CalTime
      real      :: Ufree, Pfree, Tfree, Rfree                                ! 主流速度・圧力・温度
      real      :: GravityX, GravityY
      real      :: IceCt, IceSt, IceAt
      real      :: PhaseRK, timeRK, dtRK
      character :: fname*50
      ! 処理開始 ********************************************************************************************
      ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 計算回数 --------------------------------------------------------------------------------------------
      CalTime = IceTimeMax/real(IceStepMax)
      nCountMax = nint(CalTime/dti)
      ! 対象ブロック ----------------------------------------------------------------------------------------
      m = mRef
      ! メモリ確保 ------------------------------------------------------------------------------------------
      allocate (Utau(Ice(m)%is:Ice(m)%ie), &
      &         RH(Ice(m)%is:Ice(m)%ie), &
      &         hc(Ice(m)%is:Ice(m)%ie), &
      &         FF(Ice(m)%is:Ice(m)%ie), &
      &         B0(Ice(m)%is:Ice(m)%ie))
      if (IceStep == 0) FF(:) = 0.0
      ! 主流 ------------------------------------------------------------------------------------------------
! Ufree = VelExp; Pfree = PsExp; Tfree = TsExp
      Ufree = DrpInVel; Pfree = PsExp; Tfree = TsExp
      ! 重力 ------------------------------------------------------------------------------------------------
! GravityX = +9.806650 * sin(AOA)
! GravityY = -9.806650 * cos(AOA)
      GravityX = 0.0
      GravityY = 0.0
      ! 表面粗さ --------------------------------------------------------------------------------------------
      select case (RoughNum)
      case (1)
         call RoughnessShinBond2D( &
         &      Ice(m)%is, Ice(m)%ie, Ice(m)%Ts0, LWC, &
         &      RH)
      case (2)
         call RoughnessCIRAMIL2D( &
         &      Ice(m)%is, Ice(m)%ie, Ice(m)%Ts0, &
         &      RH)
      case (3)
         call WallFrictionVelocity2D( &
         &      Ice(m)%is, Ice(m)%ie, &
         &      Ice(m)%fRough, Ice(m)%yp, Ice(m)%Velp, Ice(m)%nup, RH, &
         &      utau)
         call RoughnessWright2D( &
         &      Ice(m)%is, Ice(m)%ie, Ice(m)%x, Ice(m)%y, &
         &      Ice(m)%SA, Ice(m)%Up, Ice(m)%Vp, Utau, Ice(m)%Rho0, LWC, &
         &      RH)
      case default; write (*, '(a)') '!!!!! Error : Rougness model number !!!!!'
      end select
      ! 摩擦速度 --------------------------------------------------------------------------------------------
      select case (TurbNum)
      case (4) !低Re数型モデル
         call WallFrictionVelocityLRe2D( &
         &      Ice(m)%is, Ice(m)%ie, &
         &      Ice(m)%yp, Ice(m)%Velp, Ice(m)%nup, utau)
      case default
         call WallFrictionVelocity2D( &
         &      Ice(m)%is, Ice(m)%ie, &
         &      Ice(m)%fRough, Ice(m)%yp, Ice(m)%Velp, Ice(m)%nup, RH, &
         &      utau)
      end select
      ! 熱伝達率 --------------------------------------------------------------------------------------------
      call HeatTransferCoefficient2D( &
      &      Ice(m)%is, Ice(m)%ie, &
      &      Ice(m)%Ue, Ice(m)%mu0, Ice(m)%rho0, utau, RH, &
      &      hc)
      ! 着氷計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      select case (ThermoNum)
         ! Original Messinger ---------------------------------------------------------------------------------
      case (1)
         write (*, '(/,a)') '+++ Original Messinger +++'
!   ! メモリ確保
!   allocate( Q    ( Ice(m)%is: Ice(m)%ie ), &
!   &         Him  ( Ice(m)%is: Ice(m)%ie ), &
!   &         Hes  ( Ice(m)%is: Ice(m)%ie ), &
!   &         Hac  ( Ice(m)%is: Ice(m)%ie ), &
!   &         Hrin ( Ice(m)%is: Ice(m)%ie ), &
!   &         Hrout( Ice(m)%is: Ice(m)%ie ), &
!   &         Rhoi ( Ice(m)%is: Ice(m)%ie ) )
!   ! 初期条件
!   Him(:) = 0.0; Hes(:) = 0.0; Hac(:) = 0.0; Hrin(:) = 0.0; Hrout(:) = 0.0
!   ! 氷の密度
!   call IceDensity2D( &
!   &      Ice(m)%is, Ice(m)%ie, Tfree, Rhoi )
!   ! 衝突液滴の質量・エンタルピ
!   call DropImpingementOrg2D( &
!   &      Ice(m)%is, Ice(m)%ie, MVD, Rhod, Tfree, &
!   &      Ice(m)%SA, Ice(m)%CE, Ice(m)%Uimp, Ice(m)%Vimp, &
!   &      Ice(m)%Mim, Him )
!   ! 昇華・蒸発の質量・エンタルピ
!   call EvapolationOrg2D( &
!   &      Ice(m)%is, Ice(m)%ie, &
!   &      Tfree, Ice(m)%Ts0, Ice(m)%Ptp, hc, Ice(m)%SA, Ice(m)%dBi, &
!   &      Ice(m)%Mes, Hes )
!   ! 熱の収支
!   call HeatBalanceOrg2D( &
!   &      Ice(m)%is, Ice(m)%ie, &
!   &      Tfree, Ice(m)%Ts0, hc, Ice(m)%SA, utau, &
!   &      Q )
         ! Extended Messinger ---------------------------------------------------------------------------------
      case (2)
         write (*, '(/,a)') '+++ Extended Messinger +++'
         ! メモリ確保
         allocate (Bw(Ice(m)%is:Ice(m)%ie), &
         &         Tw(Ice(m)%is:Ice(m)%ie), &
         &         Bg(Ice(m)%is:Ice(m)%ie), &
         &         tg(Ice(m)%is:Ice(m)%ie), &
         &         Q0(Ice(m)%is:Ice(m)%ie), &
         &         Q1(Ice(m)%is:Ice(m)%ie))
         ! 初期条件
         Ice(m)%Ti(:) = Ice(m)%Ts0(:); Tw(:) = 273.15
         Bw(:) = 0.0; Tw(:) = 0.0
         Bg(:) = 0.0; Tg(:) = 0.0
         Q0(:) = 0.0; Q1(:) = 0.0
!   ! 表面温度
!   call SurfaceTemperature2D( &
!   &      Ice(m)%is, Ice(m)%ie, Ufree, Tfree, Ice(m)%Ue, &
!   &      Ice(m)%Ta, Ice(m)%Ts0 )
         ! 衝突液滴の質量流束
         call DropImpingementExt2D( &
         &      Ice(m)%is, Ice(m)%ie, &
         &      MVD, Rhod, Ice(m)%CE, &
         &      Ice(m)%Mim)
      case default
         write (*, '(a)') '!!!!! Error : Thermodynamics model number !!!!!'
         stop
      end select
      ! 着氷時間進行 ========================================================================================
      write (*, '(a,e16.8e3)') '* Start icing time = ', IceTime
      B0(:) = Ice(m)%Bi(:)
      do n = nStart, nCountMax
         ! 1 ステップ前の氷層厚さ保存
         call SaveIceThickness2D( &
         &      Ice(m)%is, Ice(m)%ie, &
         &      Ice(m)%Bi, Ice(m)%Bi0)
         select case (ThermoNum)
            ! Original Messinger --------------------------------------------------------------------------------
         case (1)
!    ! 氷結率
!    call FreezingFractionOrg2D( &
!    &      Ice(m)%is, Ice(m)%ie, &
!    &      Ice(m)%Mim, Ice(m)%Mes, Ice(m)%Mrin, Ice(m)%Mst, Him, Hes, Hrin, Hrout, Hac, Q, &
!    &      FF )
!    ! Runback-out の質量・エンタルピ
!    call RunbackOutOrg2D( &
!    &      Ice(m)%is, Ice(m)%ie, RunbackNum, dti, &
!    &      Ice(m)%Ts0, Ice(m)%SA, utau, RH, Ice(m)%Mim, Ice(m)%Mes, Ice(m)%Mrin, FF, &
!    &      Ice(m)%Mrout, Hrout, Ice(m)%Mst )
!    ! 水膜の逸脱
!    if(DropShedNum == 1) then
!      call WaterSheddingOrg2D( &
!      &      Ice(m)%is, Ice(m)%ie, dti, GravityX, GravityY, &
!      &      Ice(m)%x, Ice(m)%y, Ice(m)%Rho0, Ice(m)%up, Ice(m)%vp, &
!      &      Ice(m)%SA, utau, &
!      &      Ice(m)%Mrout, Hrout )
!    endif
!    ! Runback-in の質量・エンタルピ
!    call RunbackInOrg2D( &
!    &      Ice(m)%is, Ice(m)%ie, dti, &
!    &      GravityX, GravityY, &
!    &      Ice(m)%x, Ice(m)%y, Ice(m)%Rho0, Ice(m)%up, Ice(m)%vp, &
!    &      Ice(m)%SA, utau, Ice(m)%Mrout, Hrout, &
!    &      Ice(m)%Mrin, Hrin )
!    ! 氷層厚さ
!    call IceThicknessOrg2D( &
!    &      Ice(m)%is, Ice(m)%ie, &
!    &      Ice(m)%Ts0, Ice(m)%SA, Rhoi, FF, Ice(m)%Mim, Ice(m)%Mrin, Ice(m)%Mst, &
!    &      Ice(m)%Mac, Hac, Ice(m)%dBi )
!    ! 氷成長の時間進行
!    call IceTimeGrowth2D( &
!    &      Ice(m)%is, Ice(m)%ie, dti, &
!    &      LmtBi, Ice(m)%dBi, Ice(m)%Bi0, &
!    &      Ice(m)%Bi )
            ! Extended Messinger ---------------------------------------------------------------------------------
         case (2)
            do nn = 1, nRunge
               PhaseRK = 1.0/real(nRunge + 1 - nn)
               dtRK = PhaseRK*dti
               timeRK = IceTime + dtRK
!     ! 熱の収支
!     call HeatBalanceExt2D( &
!     &      Ice(m)%is, Ice(m)%ie, Tfree, &
!     &      Ice(m)%Ta, Ice(m)%Ptp, Ufree, hc, Ice(m)%dBi, &
!     &      Ice(m)%Mim, Ice(m)%Mrin, Ice(m)%Mst, Ice(m)%Mes, Ice(m)%nPhase, &
!     &      Q0, Q1 )
               ! 霧氷 / 雨氷の相変化
               call PhaseChangeExt2D( &
               &      Ice(m)%is, Ice(m)%ie, &
               &      hc, Tfree, Ice(m)%Ts0, Ice(m)%Ta, Ice(m)%Ptp, Ufree, &
               &      Ice(m)%Mim, Ice(m)%Mrin, Ice(m)%Mst, Ice(m)%Mes, timeRK, dtRK, Ice(m)%Bi, &
               &      Q0, Q1, Bg, tg, Ice(m)%nPhase)
               ! 熱の収支
               call HeatBalanceExt2D( &
               &      Ice(m)%is, Ice(m)%ie, Tfree, &
               &      Ice(m)%Ta, Ice(m)%Ptp, Ufree, hc, Ice(m)%dBi, &
               &      Ice(m)%Mim, Ice(m)%Mrin, Ice(m)%Mst, Ice(m)%Mes, Ice(m)%nPhase, &
               &      Q0, Q1)
               ! 蒸発・昇華の質量流束
               call EvapolationExt2D( &
               &      Ice(m)%is, Ice(m)%ie, &
               &      Tfree, Pfree, Ice(m)%Ts0, Ice(m)%Ta, Ice(m)%Ti, Tw, Ice(m)%Ptp, hc, Ice(m)%nPhase, &
               &      Ice(m)%Mes)
!     call EvapolationExt2D( &
!     &      Ice(m)%is, Ice(m)%ie, &
!     &      Ice(m)%Ta, Ice(m)%Ti, Tw, Ice(m)%Ptp, hc, Ice(m)%nPhase, &
!     &      Ice(m)%Mes )
!     ! 氷結率
!     call FreezingFractionExt2D( &
!     &      Ice(m)%is, Ice(m)%ie, TimeRK, &
!     &      Ice(m)%Mim, Ice(m)%Mes, Ice(m)%Mrin, Ice(m)%Mst, Ice(m)%nPhase, Ice(m)%Bi, Bw, Bg, Ice(m)%dBi, &
!     &      FF )
               ! 氷層厚さ・氷層温度・水膜厚さ・水膜温度
               call IceThicknessExt2D( &
               &      Ice(m)%is, Ice(m)%ie, &
               &      Ice(m)%Ts0, Ice(m)%Mim, Ice(m)%Mes, Ice(m)%Mrin, Ice(m)%Mst, &
               &      Q0, Q1, Ice(m)%Bi, Ice(m)%nPhase, timeRK, dtRK, &
               &      Bg, Tg, &
               &      Ice(m)%dBi, Bw, Ice(m)%Ti, Tw)
               ! 氷結率
               call FreezingFractionExt2D( &
               &      Ice(m)%is, Ice(m)%ie, TimeRK, &
               &      Ice(m)%Mim, Ice(m)%Mes, Ice(m)%Mrin, Ice(m)%Mst, Ice(m)%nPhase, Ice(m)%Bi, Bw, Bg, Ice(m)%dBi, &
               &      FF)
               ! Runback-out の質量流束
               call RunbackOutExt2D( &
               &      Ice(m)%is, Ice(m)%ie, RunbackNum, dtRK, timeRK, &
               &      RH, Ice(m)%Mim, Ice(m)%Mes, Ice(m)%Mrin, Bw, FF, Ice(m)%nPhase, &
               &      Ice(m)%Mrout, Ice(m)%Mst)
               ! 水膜の逸脱
               if (DropShedNum == 1) then
                  call WaterSheddingExt2D( &
                  &      Ice(m)%is, Ice(m)%ie, dtRK, GravityX, GravityY, &
                  &      Ice(m)%x, Ice(m)%y, Ice(m)%Rho0, &
                  &      Ice(m)%up, Ice(m)%vp, Ice(m)%SA, utau, &
                  &      Ice(m)%Mrout)
               end if
               ! Runback-in の質量流束
               call RunbackInExt2D( &
               &      Ice(m)%is, Ice(m)%ie, dtRK, GravityX, GravityY, &
               &      Ice(m)%x, Ice(m)%y, Ice(m)%Rho0, &
               &      Ice(m)%up, Ice(m)%vp, Ice(m)%SA, utau, Ice(m)%Mrout, &
               &      Ice(m)%Mrin)
               ! 氷成長の時間進行
               call IceTimeGrowth2D( &
               &      Ice(m)%is, Ice(m)%ie, dtRK, &
               &      PhaseRK, Ice(m)%dBi, Ice(m)%Bi0, &
               &      Ice(m)%Bi)
            end do
         end select
         ! 境界条件 -------------------------------------------------------------------------------------------
         Ice(m)%Bi(Ice(m)%is) = 0.0
         Ice(m)%Bi(Ice(m)%ie) = 0.0
         ! 時間進行 -------------------------------------------------------------------------------------------
         IceTime = IceTime + dti
         ! ログ出力 -------------------------------------------------------------------------------------------
         if (mod(n, nIceLog) == 0 .or. n == nStart .or. n == nCountMax) then
            write (*, '(/,a)') '// Calculation progress //'
            write (*, '(2(a,i9))') '* Calculation count      = ', n, ' of ', nCountMax
            write (*, '(a,e16.8e3)') '* Accretion mass         = ', maxval(Ice(m)%Mac(:))
            write (*, '(a,e16.8e3)') '* Impingement mass       = ', maxval(Ice(m)%Mim(:))
            write (*, '(a,e16.8e3)') '* Runback mass           = ', maxval(Ice(m)%Mrin(:))
            write (*, '(a,e16.8e3)') '* Stay mass              = ', maxval(Ice(m)%Mst(:))
            write (*, '(a,e16.8e3)') '* Evaporation mass       = ', maxval(Ice(m)%Mes(:))
            write (*, '(a,e16.8e3)') '* Max. ice layer         = ', maxval(Ice(m)%Bi(:))
            write (*, '(a,e16.8e3)') '* Ice volume             = ', sum(Ice(m)%Bi(:)*Ice(m)%SA(:))
            write (*, '(a,e16.8e3)') '* Exposure time          = ', IceTime
         end if
      end do
      write (*, '(a,e16.8e3)') '* Finish icing time = ', IceTime
      Ice(m)%dBi(:) = Ice(m)%Bi(:) - B0(:)
      ! 着氷箇所のフラグ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 新しく着氷した点
      do i = Ice(m)%is, Ice(m)%ie
         if (Ice(m)%Bi(i) > 1.0e-5) then
            Ice(m)%f(i) = 1
            Ice(m)%fRough(i) = 1
         else
            Ice(m)%f(i) = 0
            Ice(m)%fRough(i) = 0
         end if
      end do
      ! 処理終了 ********************************************************************************************
      return
   end subroutine CalIceAccretion
!*******************************************************************************************************
!******** 着氷翼の座標                                                                                ********
!*******************************************************************************************************
   subroutine CalIcingBlade
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer :: m, i, k
      real    :: ax, ay, az, bx, by, bz, nx, ny, nz, na
      ! 処理開始 ********************************************************************************************
      m = mRef
      ! 着氷後の翼座標を算出 -------------------------------------------------------------------------------
      do k = Flw(m)%ks, Flw(m)%ke
      do i = Ice(m)%is, Ice(m)%ie
         if (i == Ice(m)%is) then
            ax = 1.0*(-Ice(m)%x(i) + Ice(m)%x(i + 1))
            ay = 1.0*(-Ice(m)%y(i) + Ice(m)%y(i + 1))
            az = 0.0
         else if (i == Ice(m)%ie) then
            ax = 1.0*(-Ice(m)%x(i - 1) + Ice(m)%x(i))
            ay = 1.0*(-Ice(m)%y(i - 1) + Ice(m)%y(i))
            az = 0.0
         else
            ax = 0.5*(-Ice(m)%x(i - 1) + Ice(m)%x(i + 1))
            ay = 0.5*(-Ice(m)%y(i - 1) + Ice(m)%y(i + 1))
            az = 0.0
         end if
         bx = 0.0
         by = 0.0
         bz = -1.0
         nx = ay*bz - az*by
         ny = az*bx - ax*bz
         nz = ax*by - ay*bx
         na = sqrt(nx**2 + ny**2 + nz**2)
         nx = +1.0*nx/na
         ny = +1.0*ny/na
         nz = +1.0*nz/na
         Ice(m)%xi(i) = Ice(m)%x(i) + nx*Ice(m)%dBi(i)
         Ice(m)%yi(i) = Ice(m)%y(i) + ny*Ice(m)%dBi(i)
      end do
      end do
      ! 処理終了 ********************************************************************************************
      return
   end subroutine CalIcingBlade
!*******************************************************************************************************
!******** 着氷計算結果ファイル出力                                                                ********
!*******************************************************************************************************
   subroutine OutputFile
      ! 変数宣言 ********************************************************************************************
      implicit none
      ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer   :: m, i, k
      character :: fname*20
      ! 処理開始 ********************************************************************************************
      ! 計算条件ファイル出力
      ! 無次元化 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      VelIn = VelIn/aRef
      PtIn = PtIn/(rhoRef*aRef**2)
      TtIn = TtIn/aRef**2
      PsOut = PsOut/(rhoRef*aRef**2)
      TsOut = TsOut/aRef**2
      LWC = LWC/RhoRef
      MVD = MVD/LRef
      RhoD = RhoD/RhoRef
      Span = Span/LRef
      Chord = Chord/LRef
      VelExp = VelExp/aRef
      PsExp = PsExp/(rhoRef*aRef**2)
      TsExp = TsExp/aRef**2
      call Output_CalSetting(trim(ND_CalSetFile)//strtxt)
      ! 着氷データ ------------------------------------------------------------------------------------------
      m = mRef
      write (fname, '(a, 2(i2.2,a))') 'IceStep', IceStep, 'of', IceStepMax, '_'
      ! 粗さのフラグ
      allocate (Flw(m)%fRough(Flw(m)%is:Flw(m)%ie, Flw(m)%ks:Flw(m)%ke))
      Flw(m)%fRough = 0
      do k = Flw(m)%ks, Flw(m)%ke
      do i = Ice(m)%is, Ice(m)%ie
         Flw(m)%fRough(i, k) = Ice(m)%fRough(i)
      end do
      end do
      call Output_ArrayInt2D( &
      &      trim(IceCalOutDir)//trim(BlkName(m))//trim(RoughFlagFile), strdat, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%ks, Flw(m)%ke, &
      &      Flw(m)%fRough)
      call Output_ArrayInt2D( &
      &      trim(IceCalOutDir)//trim(fname)//trim(BlkName(m))//trim(RoughFlagFile), strdat, &
      &      Flw(m)%is, Flw(m)%ie, Flw(m)%ks, Flw(m)%ke, &
      &      Flw(m)%fRough)
      ! 氷層厚さ・氷層温度
      call Output_IceThickTem2D( &
      &      trim(IceCalOutDir)//trim(BlkName(m))//trim(IceThickTemFile), strdat, &
      &      Ice(m)%is, Ice(m)%ie, &
      &      Ice(m)%f, Ice(m)%Bi, Ice(m)%dBi, Ice(m)%Ti)
      call Output_IceThickTem2D( &
      &      trim(IceCalOutDir)//trim(fname)//trim(BlkName(m))//trim(IceThickTemFile), strdat, &
      &      Ice(m)%is, Ice(m)%ie, &
      &      Ice(m)%f, Ice(m)%Bi, Ice(m)%dBi, Ice(m)%Ti)
      ! 熱力学モデルパラメータ
      select case (ThermoNum)
      case (1)
         call Output_OrgMessingerPara2D( &
         &      trim(IceCalOutDir)//trim(BlkName(m))//trim(MessingerFile), strdat, &
         &      Ice(m)%is, Ice(m)%ie, &
         &      Ice(m)%Mac, Ice(m)%Mes, Ice(m)%Mim, Ice(m)%Mrin, Ice(m)%Mrout, Ice(m)%Mst)
         call Output_OrgMessingerPara2D( &
         &      trim(IceCalOutDir)//trim(fname)//trim(BlkName(m))//trim(MessingerFile), strdat, &
         &      Ice(m)%is, Ice(m)%ie, &
         &      Ice(m)%Mac, Ice(m)%Mes, Ice(m)%Mim, Ice(m)%Mrin, Ice(m)%Mrout, Ice(m)%Mst)
      case (2)
         call Output_ExtMessingerPara2D( &
         &      trim(IceCalOutDir)//trim(BlkName(m))//trim(ExMessingerFile), strdat, &
         &      Ice(m)%is, Ice(m)%ie, &
         &      Ice(m)%nPhase, &
         &      Ice(m)%Mes, Ice(m)%Mim, Ice(m)%Mrin, Ice(m)%Mrout, Ice(m)%Mst)
         call Output_ExtMessingerPara2D( &
         &      trim(IceCalOutDir)//trim(fname)//trim(BlkName(m))//trim(ExMessingerFile), strdat, &
         &      Ice(m)%is, Ice(m)%ie, &
         &      Ice(m)%nPhase, &
         &      Ice(m)%Mes, Ice(m)%Mim, Ice(m)%Mrin, Ice(m)%Mrout, Ice(m)%Mst)
      end select
      ! 着氷翼座標
      call Output_IceBladeSurface2D( &
      &      trim(IceCalOutDir)//trim(BlkName(m))//trim(IceBladeFile), strdat, &
      &      Ice(m)%is, Ice(m)%ie, &
      &      Ice(m)%xi, Ice(m)%yi)
      call Output_IceBladeSurface2D( &
      &      trim(IceCalOutDir)//trim(fname)//trim(BlkName(m))//trim(IceBladeFile), strdat, &
      &      Ice(m)%is, Ice(m)%ie, &
      &      Ice(m)%xi, Ice(m)%yi)

      ! 処理終了 ********************************************************************************************
      return
   end subroutine OutputFile

end program IceAccretion_NACA
