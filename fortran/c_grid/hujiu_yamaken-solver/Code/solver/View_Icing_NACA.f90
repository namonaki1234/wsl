!*******************************************************************************************************
!*******************************************************************************************************
!******** 着氷現象可視化プログラム								********
!******** (NACA翼，三次元圧縮性乱流場，重合格子法，高レイノルズ数型 k-eモデル)		  	********
!********					      2012.10.01  PROGRAMMED BY RYOSUKE HAYASHI ********
!*******************************************************************************************************
!*******************************************************************************************************
program ViewIcing_NACA
 ! モジュール宣言 **************************************************************************************
 use Package_NACA
 use Package_FileIO
 use Package_Flow
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! ファイル名設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: ViewIcingFile * 10 = 'ViewIcing'
 ! 処理開始 ********************************************************************************************
 ! 検証実験ケース選択 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call SelectExpCase
 write(*, '(a)') "+++++ Select Exp. case complete. +++++"
 ! 初期設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call InitialSetting
 write(*, '(a)') "+++++ Initial setting complete. +++++"
 ! 可視化計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call VisualizedIcing
 write(*, '(a)') "+++++ Icing visualized complete. +++++"
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
 character :: fname1 * 20, fname2 * 20
 ! 処理開始 ********************************************************************************************
 ! 実験条件入力
 call Input_ExpCaseNumber( &
 &      trim(ExpCaseFile) // strtxt, ExpCaseNum )
 write(fname1, '(i2.2)') ExpCaseNum
 call Input_ExpCondition( &
 &      trim(ExpConditionFile) // trim(fname1) // strtxt )
 ! 着氷経過時間ステップ入力
 call Input_IceCalStep( &
 &      trim(IceStepFile), strtxt, &
 &      nIceStep, nIceStepMax )
 ! ディレクトリ設定
 write(fname1, '(a,i2.2,a,i2.2,a)') 'step', nIceStep - 1, 'of', nIceStepMax, '//'
 write(fname2, '(a,i2.2,a)') 'case', ExpCaseNum, '//'
 if( nIceStep == 1 ) then
   GrdInDir    = bckdir // 'grid//clean//'
   OSGDir      = bckdir // 'overset//clean//'
   FlwIniDir   = bckdir // 'flow//initial//clean//' // trim(fname2)
   FlwCalInDir = bckdir // 'flow//cal//clean//' // trim(fname2)
   DrpImpDir   = bckdir // 'droplet//impingement//clean//' // trim(fname2)
  else
   GrdInDir    = bckdir // 'grid//' // trim(fname1)  // trim(fname2)
   OSGDir      = bckdir // 'overset//' // trim(fname1)  // trim(fname2)
   FlwIniDir   = bckdir // 'flow//initial//' // trim(fname1)  // trim(fname2)
   FlwCalInDir = bckdir // 'flow//cal//' // trim(fname1)  // trim(fname2)
   DrpImpDir   = bckdir // 'droplet//impingement//' // trim(fname1)  // trim(fname2)
 endif
 write(fname1  , '(a,i2.2,a,i2.2,a)') 'step', nIceStep, 'of', nIceStepMax, '//'
 IceCalInDir = bckdir // 'icing//cal//' // trim(fname1) // trim(fname2)
 IceViewDir  = bckdir // 'icing//view//' // trim(fname1) // trim(fname2)
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
 integer   :: m
 integer   :: is, ie, ks, ke
 character :: fname * 30
 ! 処理開始 ********************************************************************************************
 ! 計算条件 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Input_FlowCalCondition( &
 &      trim(FlwCalConditionFile), strtxt, &
 &      Rg, gamma, Pr, Prt, Ret, nRunge, nTVD, eTVD, Cn, Cmu, nCalMax, nOutputLog, nOutputCount )
 ! サザーランドの式の定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Input_Sutherland( &
 &      trim(FlwIniDir) // trim(ND_SutherlandFile), strtxt, s1, TsSth, muSth )
 ! 無次元化参照値 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Input_NDRef( &
 &      trim(FlwIniDir) // trim(ND_RefFile), strtxt, aRef, LRef, RhoRef )
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
 ! C 型格子分割点 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! m = 1
 m = me
 call Input_CtypeGridPoint( &
 &      trim(GrdInDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
 &      Flw(m)%i1, Flw(m)%i2, Flw(m)%i3 )
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
  &         Flw(m)%qh  (Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le), &
  &         Flw(m)%nimp(Flw(m)%is: Flw(m)%ie, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%uimp(Flw(m)%is: Flw(m)%ie, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%vimp(Flw(m)%is: Flw(m)%ie, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%wimp(Flw(m)%is: Flw(m)%ie, Flw(m)%ks: Flw(m)%ke), &
  &         Flw(m)%CE  (Flw(m)%is: Flw(m)%ie, Flw(m)%ks: Flw(m)%ke) )
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
! do m = ms, me
 m = me
  call Input_Flux3D( &
  &      trim(FlwCalInDir) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      ls, le, Flw(m)%qh )
! enddo
 ! 収集効率 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  Flw(m)%CE(:,:) = 0; Flw(m)%uimp(:,:) = 0.0; Flw(m)%vimp(:,:) = 0.0; Flw(m)%wimp(:,:) = 0.0
  call Input_Resolution2D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(IceRslFile), strtxt, &
  &      is, ie, ks, ke )
  call Input_CollectionEfficiency3D( &
  &      trim(DrpImpDir) // trim(BlkName(m)) // trim(CollectionFile), strbin, &
  &      is, ie, ks, ke, &
  &      Flw(m)%CE(is:ie,:), Flw(m)%uimp(is:ie,:), Flw(m)%vimp(is:ie,:), Flw(m)%wimp(is:ie,:) )
! enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine InitialSetting
!*******************************************************************************************************
!******** 可視化計算 										********
!*******************************************************************************************************
subroutine VisualizedIcing
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real     , pointer :: x(:, :, :), y(:, :, :), z(:, :, :), &
 &                     Tsur(:, :, :), Vimp(:, :, :), CE(:, :, :)
 character :: fname * 20
 integer   :: i, j, k, m
 integer   :: j0
 ! 処理開始 ********************************************************************************************
! do m = ms, me
 m = me
  ! メモリ確保 -----------------------------------------------------------------------------------------
  j0 = Flw(m)%js
  allocate( x   (Flw(m)%is: Flw(m)%ie, 0:1, Flw(m)%ks: Flw(m)%ke), &
  &         y   (Flw(m)%is: Flw(m)%ie, 0:1, Flw(m)%ks: Flw(m)%ke), &
  &         z   (Flw(m)%is: Flw(m)%ie, 0:1, Flw(m)%ks: Flw(m)%ke), &
  &         Tsur(Flw(m)%is: Flw(m)%ie, 0:1, Flw(m)%ks: Flw(m)%ke), &
  &         Vimp(Flw(m)%is: Flw(m)%ie, 0:1, Flw(m)%ks: Flw(m)%ke), &
  &         CE  (Flw(m)%is: Flw(m)%ie, 0:1, Flw(m)%ks: Flw(m)%ke) )
  ! 物理量に変換 ---------------------------------------------------------------------------------------
  call SetPhysics3DKEM( &
  &      Rg, gamma, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, ls, le, &
  &      Flw(m)%qh, Flw(m)%jac, &
  &      Flw(m)%rho, Flw(m)%u, Flw(m)%v, Flw(m)%w, Flw(m)%p, Flw(m)%t, Flw(m)%kin, Flw(m)%eps )
  ! 翼座標 ---------------------------------------------------------------------------------------------
  do k = Flw(m)%ks, Flw(m)%ke
  do i = Flw(m)%is, Flw(m)%ie
   x(i,0,k) = Flw(m)%x(i,j0,k) * lRef
   y(i,0,k) = Flw(m)%y(i,j0,k) * lRef
   z(i,0,k) = Flw(m)%z(i,j0,k) * lRef
  enddo
  enddo
  ! 表面温度 -------------------------------------------------------------------------------------------
  do k = Flw(m)%ks, Flw(m)%ke
  do i = Flw(m)%is, Flw(m)%ie
   Tsur(i,0,k) = Flw(m)%t(i,j0,k) * aRef**2
  enddo
  enddo
  ! 液滴衝突速度 ----------------------------------------------------------------------------------------
  do k = Flw(m)%ks, Flw(m)%ke
  do i = Flw(m)%is, Flw(m)%ie
   Vimp(i,0,k) = sqrt( Flw(m)%uimp(i,k)**2 + Flw(m)%vimp(i,k)**2 + Flw(m)%wimp(i,k)**2 )
  enddo
  enddo
  ! 液滴収集効率 ---------------------------------------------------------------------------------------
  do k = Flw(m)%ks, Flw(m)%ke
  do i = Flw(m)%is, Flw(m)%ie
   CE(i,0,k) = Flw(m)%CE(i,k)
  enddo
  enddo
  ! 三次元化 -------------------------------------------------------------------------------------------
  do k = Flw(m)%ks, Flw(m)%ke
  do i = Flw(m)%is, Flw(m)%ie
   x   (i,1,k) = x   (i,0,k)
   y   (i,1,k) = y   (i,0,k)
   z   (i,1,k) = z   (i,0,k)
   Tsur(i,1,k) = Tsur(i,0,k)
   Vimp(i,1,k) = Vimp(i,0,k)
   CE  (i,1,k) = CE  (i,0,k)
  enddo
  enddo
  ! ファイル出力 ---------------------------------------------------------------------------------------
  call MakeMAVSFile3D( &
  &      trim(IceViewDir), trim(BlkName(m)) // trim(ViewIcingFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, 0, 1, Flw(m)%ks, Flw(m)%ke, &
  &      Tsur, CE, Vimp,  &
  &      x, y, z )
  deallocate( x, y, z, Tsur, Vimp, CE )
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine VisualizedIcing
!*******************************************************************************************************
!******** MicroAVSファイル作成（三次元）  		                                        ********
!*******************************************************************************************************
subroutine MakeMAVSFile3D( &
&            strdir, strname, ext, is, ie, js, je, ks, ke, &
&            Tsur, CE, Vimp, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strdir*(*), strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 real     , intent(in)  :: Tsur(is:ie, js:je, ks:ke), CE(is:ie, js:je, ks:ke), Vimp(is:ie, js:je, ks:ke)
 real     , intent(in)  :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke), z(is:ie, js:je, ks:ke)
 ! 局所定数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: strbin * 4 = '.bin'
 character, parameter :: strfld * 4 = '.fld'
 integer  , parameter :: ndim   =  3
 integer  , parameter :: nspace =  3
 integer  , parameter :: veclen =  3
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
  write(1, '(a)')    'label  = Tsur, CE, Vimp'
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
    write(1) Tsur; write(1) CE  ; write(1) Vimp
    write(1) x   ; write(1) y   ; write(1) z
   close(1)
 ! Aschii ----------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strdir) // trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     write(1, '(e16.5e3,15(x,e16.8e3))') Tsur(i,j,k), CE(i,j,k), Vimp(i,j,k), &
     &                                   x(i,j,k), y(i,j,k), z(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine MakeMAVSFile3D
! 定義終了 *********************************************************************************************
end program ViewIcing_NACA
