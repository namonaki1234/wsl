!*******************************************************************************************************
!*******************************************************************************************************
!******** パッケージ型モジュール								********
!******** ファイル入出力									********
!******** 			                      2013.01.27  PROGRAMMED BY RYOSUKE HAYASHI ********
!******** 			                      2013.06.17     UPDATED BY RYOSUKE HAYASHI ********
!*******************************************************************************************************
!*******************************************************************************************************
module Package_FileIO
 ! 変数宣言 ********************************************************************************************
 implicit none
 private
 ! 拡張子設定 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter, public :: bckdir * 3 = '../'
 character, parameter, public :: strdir * 2 = './'
 character, parameter, public :: strbin * 4 = '.bin'
 character, parameter, public :: strdat * 4 = '.dat'
 character, parameter, public :: strfld * 4 = '.fld'
 character, parameter, public :: strtxt * 4 = '.txt'
 ! 共有サブルーチン ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 public :: Input_Resolution1D, Output_Resolution1D
 public :: Input_Resolution2D, Output_Resolution2D, Input_Resolution3D, Output_Resolution3D
 public :: Input_Grid2D      , Output_Grid2D      , Input_Grid3D      , Output_Grid3D
 public :: Input_Metrics2D   , Output_Metrics2D   , Input_Metrics3D   , Output_Metrics3D
 public :: Input_Flux2D      , Output_Flux2D      , Input_Flux3D      , Output_Flux3D
 public :: Input_OversetCoe2D, Output_OversetCoe2D, Input_OversetCoe3D, Output_OversetCoe3D
 public :: Input_Physics2D   , Output_Physics2D   , Input_Physics3D   , Output_Physics3D
 public :: Input_ArrayReal1D , Output_ArrayReal1D , Input_ArrayReal2D , Output_ArrayReal2D, &
 &         Input_ArrayReal3D , Output_ArrayReal3D
 public :: Input_ArrayInt1D  , Output_ArrayInt1D  , Input_ArrayInt2D  , Output_ArrayInt2D, &
 &         Input_ArrayInt3D  , Output_ArrayInt3D
 public :: Input_GridPoint , Output_GridPoint
 public :: Input_DropInCondition, Output_DropInCondition
 public :: Input_Impingement3D  , Output_Impingement3D
 public :: Input_CollectionEfficiency3D, Output_CollectionEfficiency3D
 public :: Input_IceThickTem2D    , Output_IceThickTem2D
 public :: Input_IceThickTem3D    , Output_IceThickTem3D
 public :: Input_OrgMessingerPara2D, Output_OrgMessingerPara2D
 public :: Input_ExtMessingerPara2D, Output_ExtMessingerPara2D
 public :: Input_OrgMessingerPara3D, Output_OrgMessingerPara3D
 public :: Input_ExtMessingerPara3D, Output_ExtMessingerPara3D
 public :: Input_IceBladeSurface2D, Output_IceBladeSurface2D
 public :: Input_IceBladeSurface3D, Output_IceBladeSurface3D
 ! 内部手続き ******************************************************************************************
 contains
!*******************************************************************************************************
!******** 解像度入力 (一次元) 									********
!*******************************************************************************************************
subroutine Input_Resolution1D( &
&            strname, ext, is, ie )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(out) :: is, ie
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
  read(1, '(i4,x,i4)') is, ie
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_Resolution1D
!*******************************************************************************************************
!******** 解像度出力 (一次元) 									********
!*******************************************************************************************************
subroutine Output_Resolution1D( &
&            strname, ext, is, ie )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
  write(1, '(i4,x,i4)') is, ie
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_Resolution1D
!*******************************************************************************************************
!******** 解像度入力 (二次元) 									********
!*******************************************************************************************************
subroutine Input_Resolution2D( &
&            strname, ext, is, ie, js, je )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(out) :: is, ie, js, je
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
  read(1, '(i4,3(x,i4))') is, ie, js, je
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_Resolution2D
!*******************************************************************************************************
!******** 解像度出力 (二次元) 									********
!*******************************************************************************************************
subroutine Output_Resolution2D( &
&            strname, ext, is, ie, js, je )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
  write(1, '(i4,3(x,i4))') is, ie, js, je
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_Resolution2D
!*******************************************************************************************************
!******** 解像度入力 (三次元) 									********
!*******************************************************************************************************
subroutine Input_Resolution3D( &
&            strname, ext, is, ie, js, je, ks, ke )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(out) :: is, ie, js, je, ks, ke
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
  read(1, '(i4,5(x,i4))') is, ie, js, je, ks, ke
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_Resolution3D
!*******************************************************************************************************
!******** 解像度出力 (三次元) 									********
!*******************************************************************************************************
subroutine Output_Resolution3D( &
&            strname, ext, is, ie, js, je, ks, ke )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je, ks, ke
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
  write(1, '(i4,5(x,i4))') is, ie, js, je, ks, ke
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_Resolution3D
!*******************************************************************************************************
!******** 格子座標入力 (二次元) 								********
!*******************************************************************************************************
subroutine Input_Grid2D( &
&            strname, ext, is, ie, js, je, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je
 real     , intent(out) :: x(is:ie, js:je), y(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) ((x(i,j), i = is, ie), j = js, je)
    read(1) ((y(i,j), i = is, ie), j = js, je)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do j = js, je
    do i = is, ie
     read(1, '(e16.8e3,1(x,e16.8e3))') x(i,j), y(i,j)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_Grid2D
!*******************************************************************************************************
!******** 格子座標出力 (二次元) 								********
!*******************************************************************************************************
subroutine Output_Grid2D( &
&            strname, ext, is, ie, js, je, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je
 real     , intent(in) :: x(is:ie, js:je), y(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) ((x(i,j), i = is, ie), j = js, je)
    write(1) ((y(i,j), i = is, ie), j = js, je)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do j = js, je
    do i = is, ie
     write(1, '(e16.8e3,1(x,e16.8e3))') x(i,j), y(i,j)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_Grid2D
!*******************************************************************************************************
!******** 格子座標入力  (三次元)								********
!*******************************************************************************************************
subroutine Input_Grid3D( &
&            strname, ext, is, ie, js, je, ks, ke, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 real     , intent(out) :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke), z(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) (((x(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((y(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((z(i,j,k), i = is, ie), j = js, je), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     read(1, '(e16.8e3,2(x,e16.8e3))') x(i,j,k), y(i,j,k), z(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_Grid3D
!*******************************************************************************************************
!******** 格子座標出力 (三次元)									********
!*******************************************************************************************************
subroutine Output_Grid3D( &
&            strname, ext, is, ie, js, je, ks, ke, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je, ks, ke
 real     , intent(in) :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke), z(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) (((x(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((y(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((z(i,j,k), i = is, ie), j = js, je), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     write(1, '(e16.8e3,2(x,e16.8e3))') x(i,j,k), y(i,j,k), z(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_Grid3D
!*******************************************************************************************************
!******** メトリックス入力 (二次元) 								********
!*******************************************************************************************************
subroutine Input_Metrics2D( &
&            strname, ext, is, ie, js, je, jac, xix, xiy, etx, ety )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je
 real     , intent(out) :: jac(is:ie, js:je)
 real     , intent(out) :: xix(is:ie, js:je), xiy(is:ie, js:je), &
 &                         etx(is:ie, js:je), ety(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) ((jac(i,j), i = is, ie), j = js, je)
    read(1) ((xix(i,j), i = is, ie), j = js, je)
    read(1) ((xiy(i,j), i = is, ie), j = js, je)
    read(1) ((etx(i,j), i = is, ie), j = js, je)
    read(1) ((ety(i,j), i = is, ie), j = js, je)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do j = js, je
    do i = is, ie
     read(1, '(e16.8e3,4(x,e16.8e3))') jac(i,j), xix(i,j), xiy(i,j), etx(i,j), ety(i,j)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_Metrics2D
!*******************************************************************************************************
!******** メトリックス出力 (二次元) 								********
!*******************************************************************************************************
subroutine Output_Metrics2D( &
&            strname, ext, is, ie, js, je, jac, xix, xiy, etx, ety )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je
 real     , intent(in) :: jac(is:ie, js:je)
 real     , intent(in) :: xix(is:ie, js:je), xiy(is:ie, js:je), &
 &                        etx(is:ie, js:je), ety(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) ((jac(i,j), i = is, ie), j = js, je)
    write(1) ((xix(i,j), i = is, ie), j = js, je)
    write(1) ((xiy(i,j), i = is, ie), j = js, je)
    write(1) ((etx(i,j), i = is, ie), j = js, je)
    write(1) ((ety(i,j), i = is, ie), j = js, je)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do j = js, je
    do i = is, ie
     write(1, '(e16.8e3,4(x,e16.8e3))') jac(i,j), xix(i,j), xiy(i,j), etx(i,j), ety(i,j)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_Metrics2D
!*******************************************************************************************************
!******** メトリックス入力 (三次元) 								********
!*******************************************************************************************************
subroutine Input_Metrics3D( &
&            strname, ext, is, ie, js, je, ks, ke, jac, xix, xiy, xiz, etx, ety, etz, zex, zey, zez )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 real     , intent(out) :: jac(is:ie, js:je, ks:ke)
 real     , intent(out) :: xix(is:ie, js:je, ks:ke), &
 &                         xiy(is:ie, js:je, ks:ke), &
 &                         xiz(is:ie, js:je, ks:ke), &
 &                         etx(is:ie, js:je, ks:ke), &
 &                         ety(is:ie, js:je, ks:ke), &
 &                         etz(is:ie, js:je, ks:ke), &
 &                         zex(is:ie, js:je, ks:ke), &
 &                         zey(is:ie, js:je, ks:ke), &
 &                         zez(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) (((jac(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((xix(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((xiy(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((xiz(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((etx(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((ety(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((etz(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((zex(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((zey(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((zez(i,j,k), i = is, ie), j = js, je), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     read(1, '(e16.8e3,9(x,e16.8e3))') jac(i,j,k), &
     &                                 xix(i,j,k), xiy(i,j,k), xiz(i,j,k), &
     &                                 etx(i,j,k), ety(i,j,k), etz(i,j,k), &
     &                                 zex(i,j,k), zey(i,j,k), zez(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_Metrics3D
!*******************************************************************************************************
!******** メトリックス出力 (三次元) 								********
!*******************************************************************************************************
subroutine Output_Metrics3D( &
&            strname, ext, is, ie, js, je, ks, ke, jac, xix, xiy, xiz, etx, ety, etz, zex, zey, zez )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je, ks, ke
 real     , intent(in) :: jac(is:ie, js:je, ks:ke)
 real     , intent(in) :: xix(is:ie, js:je, ks:ke), &
 &                        xiy(is:ie, js:je, ks:ke), &
 &                        xiz(is:ie, js:je, ks:ke), &
 &                        etx(is:ie, js:je, ks:ke), &
 &                        ety(is:ie, js:je, ks:ke), &
 &                        etz(is:ie, js:je, ks:ke), &
 &                        zex(is:ie, js:je, ks:ke), &
 &                        zey(is:ie, js:je, ks:ke), &
 &                        zez(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) (((jac(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((xix(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((xiy(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((xiz(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((etx(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((ety(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((etz(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((zex(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((zey(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((zez(i,j,k), i = is, ie), j = js, je), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     write(1, '(e16.8e3,9(x,e16.8e3))') jac(i,j,k), &
     &                                  xix(i,j,k), xiy(i,j,k), xiz(i,j,k), &
     &                                  etx(i,j,k), ety(i,j,k), etz(i,j,k), &
     &                                  zex(i,j,k), zey(i,j,k), zez(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_Metrics3D
!*******************************************************************************************************
!******** 流束入力 (二次元) 									********
!*******************************************************************************************************
subroutine Input_Flux2D( &
&            strname, ext, is, ie, js, je, ls, le, qh )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je
 integer  , intent(in)  :: ls, le
 real     , intent(out) :: qh(is:ie, js:je, ls:le)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character :: strnum
 integer   :: i, j, l
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    do l = ls, le
     read(1) ((qh(i,j,l), i = is, ie), j = js, je)
    enddo
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   write(strnum, '(i1)') le-ls
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do j = js, je
    do i = is, ie
     read(1, '(e16.8e3, ' // strnum // '(x,e16.8e3))') qh(i,j,ls:le)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_Flux2D
!*******************************************************************************************************
!******** 流束出力 (二次元) 									********
!*******************************************************************************************************
subroutine Output_Flux2D( &
&            strname, ext, is, ie, js, je, ls, le, qh )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je
 integer  , intent(in) :: ls, le
 real     , intent(in) :: qh(is:ie, js:je, ls:le)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character :: strnum
 integer   :: i, j, l
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    do l = ls, le
     write(1) ((qh(i,j,l), i = is, ie), j = js, je)
    enddo
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   write(strnum, '(i1)') le-ls
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do j = js, je
    do i = is, ie
     write(1, '(e16.8e3, ' // strnum // '(x,e16.8e3))') qh(i,j,ls:le)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_Flux2D
!*******************************************************************************************************
!******** 流束入力 (三次元) 									********
!*******************************************************************************************************
subroutine Input_Flux3D( &
&            strname, ext, is, ie, js, je, ks, ke, ls, le, qh )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 integer  , intent(in)  :: ls, le
 real     , intent(out) :: qh(is:ie, js:je, ks:ke, ls:le)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character :: strnum
 integer   :: i, j, k, l
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    do l = ls, le
     read(1) (((qh(i,j,k,l), i = is, ie), j = js, je), k = ks, ke)
    enddo
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   write(strnum, '(i1)') le-ls
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     read(1, '(e16.8e3, ' // strnum // '(x,e16.8e3))') qh(i,j,k,ls:le)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_Flux3D
!*******************************************************************************************************
!******** 流束出力 (三次元) 									********
!*******************************************************************************************************
subroutine Output_Flux3D( &
&            strname, ext, is, ie, js, je, ks, ke, ls, le, qh )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je, ks, ke
 integer  , intent(in) :: ls, le
 real     , intent(in) :: qh(is:ie, js:je, ks:ke, ls:le)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character :: strnum
 integer   :: i, j, k, l
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    do l = ls, le
     write(1) (((qh(i,j,k,l), i = is, ie), j = js, je), k = ks, ke)
    enddo
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   write(strnum, '(i1)') le-ls
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     write(1, '(e16.8e3, ' // strnum // '(x,e16.8e3))') qh(i,j,k,ls:le)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_Flux3D
!*******************************************************************************************************
!******** 重合格子法補間係数入力 (二次元)							********
!*******************************************************************************************************
subroutine Input_OversetCoe2D( &
&            strname, ext, is, ie, js, je, fover, ip, jp, term1, term2, term3, term4 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je
 integer  , intent(out) :: fover(is:ie, js:je)
 integer  , intent(out) :: ip(is:ie, js:je), jp(is:ie, js:je)
 real     , intent(out) :: term1(is:ie, js:je), term2(is:ie, js:je), &
 &                         term3(is:ie, js:je), term4(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, j
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) ((fover(i,j), i = is, ie), j = js, je)
    read(1) ((   ip(i,j), i = is, ie), j = js, je)
    read(1) ((   jp(i,j), i = is, ie), j = js, je)
    read(1) ((term1(i,j), i = is, ie), j = js, je)
    read(1) ((term2(i,j), i = is, ie), j = js, je)
    read(1) ((term3(i,j), i = is, ie), j = js, je)
    read(1) ((term4(i,j), i = is, ie), j = js, je)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do j = js, je
    do i = is, ie
     read(1, '(i1,2(x,i4),4(x,e16.8e3))') fover(i,j), ip(i,j), jp(i,j), &
     &                                    term1(i,j), term2(i,j), term3(i,j), term4(i,j)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_OversetCoe2D
!*******************************************************************************************************
!******** 重合格子法補間係数出力 (二次元)							********
!*******************************************************************************************************
subroutine Output_OversetCoe2D( &
&            strname, ext, is, ie, js, je, fover, ip, jp, term1, term2, term3, term4 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je
 integer  , intent(in) :: fover(is:ie, js:je)
 integer  , intent(in) :: ip(is:ie, js:je), jp(is:ie, js:je)
 real     , intent(in) :: term1(is:ie, js:je), term2(is:ie, js:je), &
 &                        term3(is:ie, js:je), term4(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, j
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) ((fover(i,j), i = is, ie), j = js, je)
    write(1) ((   ip(i,j), i = is, ie), j = js, je)
    write(1) ((   jp(i,j), i = is, ie), j = js, je)
    write(1) ((term1(i,j), i = is, ie), j = js, je)
    write(1) ((term2(i,j), i = is, ie), j = js, je)
    write(1) ((term3(i,j), i = is, ie), j = js, je)
    write(1) ((term4(i,j), i = is, ie), j = js, je)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do j = js, je
    do i = is, ie
     write(1, '(i1,2(x,i4),4(x,e16.8e3))') fover(i,j), ip(i,j), jp(i,j), &
     &                                     term1(i,j), term2(i,j), term3(i,j), term4(i,j)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_OversetCoe2D
!*******************************************************************************************************
!******** 重合格子法補間係数入力 (三次元)							********
!*******************************************************************************************************
subroutine Input_OversetCoe3D( &
&            strname, ext, is, ie, js, je, ks, ke, fover, ip, jp, kp, &
&            term1, term2, term3, term4, term5, term6, term7, term8 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 integer  , intent(out) :: fover(is:ie, js:je, ks:ke)
 integer  , intent(out) :: ip(is:ie, js:je, ks:ke), jp(is:ie, js:je, ks:ke), kp(is:ie, js:je, ks:ke)
 real     , intent(out) :: term1(is:ie, js:je, ks:ke), term2(is:ie, js:je, ks:ke), &
 &                         term3(is:ie, js:je, ks:ke), term4(is:ie, js:je, ks:ke), &
 &                         term5(is:ie, js:je, ks:ke), term6(is:ie, js:je, ks:ke), &
 &                         term7(is:ie, js:je, ks:ke), term8(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, j, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) (((fover(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((   ip(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((   jp(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((   kp(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((term1(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((term2(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((term3(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((term4(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((term5(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((term6(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((term7(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((term8(i,j,k), i = is, ie), j = js, je), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     read(1, '(i1,3(x,i4),8(x,e16.8e3))') fover(i,j,k), ip(i,j,k), jp(i,j,k), kp(i,j,k), &
     &                                    term1(i,j,k), term2(i,j,k), term3(i,j,k), term4(i,j,k), &
     &                                    term5(i,j,k), term6(i,j,k), term7(i,j,k), term8(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_OversetCoe3D
!*******************************************************************************************************
!******** 重合格子法補間係数出力 (三次元)							********
!*******************************************************************************************************
subroutine Output_OversetCoe3D( &
&            strname, ext, is, ie, js, je, ks, ke, fover, ip, jp, kp, &
&            term1, term2, term3, term4, term5, term6, term7, term8 )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je, ks, ke
 integer  , intent(in) :: fover(is:ie, js:je, ks:ke)
 integer  , intent(in) :: ip(is:ie, js:je, ks:ke), jp(is:ie, js:je, ks:ke), kp(is:ie, js:je, ks:ke)
 real     , intent(in) :: term1(is:ie, js:je, ks:ke), term2(is:ie, js:je, ks:ke), &
 &                        term3(is:ie, js:je, ks:ke), term4(is:ie, js:je, ks:ke), &
 &                        term5(is:ie, js:je, ks:ke), term6(is:ie, js:je, ks:ke), &
 &                        term7(is:ie, js:je, ks:ke), term8(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, j, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) (((fover(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((   ip(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((   jp(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((   kp(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((term1(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((term2(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((term3(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((term4(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((term5(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((term6(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((term7(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((term8(i,j,k), i = is, ie), j = js, je), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     write(1, '(i1,3(x,i4),8(x,e16.8e3))') fover(i,j,k), ip(i,j,k), jp(i,j,k), kp(i,j,k), &
     &                                     term1(i,j,k), term2(i,j,k), term3(i,j,k), term4(i,j,k), &
     &                                     term5(i,j,k), term6(i,j,k), term7(i,j,k), term8(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_OversetCoe3D
!*******************************************************************************************************
!******** 物理量入力 (二次元)									********
!*******************************************************************************************************
subroutine Input_Physics2D( &
&            strname, ext, is, ie, js, je, rho, u, v, ps, ts, mu, kin, eps, mut )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je
 real     , intent(out) :: rho(is:ie, js:je), u(is:ie, js:je), v(is:ie, js:je), &
 &                         ps(is:ie, js:je), ts(is:ie, js:je), mu(is:ie, js:je), &
 &                         kin(is:ie, js:je), eps(is:ie, js:je), mut(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, j
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) ((rho(i,j), i = is, ie), j = js, je)
    read(1) ((  u(i,j), i = is, ie), j = js, je)
    read(1) ((  v(i,j), i = is, ie), j = js, je)
    read(1) (( ps(i,j), i = is, ie), j = js, je)
    read(1) (( ts(i,j), i = is, ie), j = js, je)
    read(1) (( mu(i,j), i = is, ie), j = js, je)
    read(1) ((kin(i,j), i = is, ie), j = js, je)
    read(1) ((eps(i,j), i = is, ie), j = js, je)
    read(1) ((mut(i,j), i = is, ie), j = js, je)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do j = js, je
    do i = is, ie
     read(1, '(e16.8e3,8(x,e16.8e3))') rho(i,j), u(i,j), v(i,j), ps(i,j), ts(i,j), mu(i,j), &
     &                                 kin(i,j), eps(i,j), mut(i,j)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_Physics2D
!*******************************************************************************************************
!******** 物理量出力 (二次元)									********
!*******************************************************************************************************
subroutine Output_Physics2D( &
&            strname, ext, is, ie, js, je, rho, u, v, ps, ts, mu, kin, eps, mut )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je
 real     , intent(in) :: rho(is:ie, js:je), u(is:ie, js:je), v(is:ie, js:je), &
 &                        ps(is:ie, js:je), ts(is:ie, js:je), mu(is:ie, js:je), &
 &                        kin(is:ie, js:je), eps(is:ie, js:je), mut(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, j
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) ((rho(i,j), i = is, ie), j = js, je)
    write(1) ((  u(i,j), i = is, ie), j = js, je)
    write(1) ((  v(i,j), i = is, ie), j = js, je)
    write(1) (( ps(i,j), i = is, ie), j = js, je)
    write(1) (( ts(i,j), i = is, ie), j = js, je)
    write(1) (( mu(i,j), i = is, ie), j = js, je)
    write(1) ((kin(i,j), i = is, ie), j = js, je)
    write(1) ((eps(i,j), i = is, ie), j = js, je)
    write(1) ((mut(i,j), i = is, ie), j = js, je)
    write(1) ((mut(i,j), i = is, ie), j = js, je)
    write(1) ((mut(i,j), i = is, ie), j = js, je)
    write(1) ((mut(i,j), i = is, ie), j = js, je)


   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do j = js, je
    do i = is, ie
     write(1, '(e16.8e3,8(x,e16.8e3))') rho(i,j), u(i,j), v(i,j), ps(i,j), ts(i,j), mu(i,j), &
     &                                  kin(i,j), eps(i,j), mut(i,j)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_Physics2D
!*******************************************************************************************************
!******** 物理量入力 (三次元)									********
!*******************************************************************************************************
subroutine Input_Physics3D( &
&            strname, ext, is, ie, js, je, ks, ke, rho, u, v, w, ps, ts, mu, kin, eps, mut )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 real     , intent(out) :: rho(is:ie, js:je, ks:ke), &
 &                         u(is:ie, js:je, ks:ke), v(is:ie, js:je, ks:ke), w(is:ie, js:je, ks:ke), &
 &                         ps(is:ie, js:je, ks:ke), ts(is:ie, js:je, ks:ke), mu(is:ie, js:je, ks:ke), &
 &                         kin(is:ie, js:je, ks:ke), eps(is:ie, js:je, ks:ke), mut(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, j, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) (((rho(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((  u(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((  v(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((  w(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) ((( ps(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) ((( ts(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) ((( mu(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((kin(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((eps(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    read(1) (((mut(i,j,k), i = is, ie), j = js, je), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     read(1, '(e16.8e3,9(x,e16.8e3))') rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), &
     &                                 ps(i,j,k), ts(i,j,k), mu(i,j,k), &
     &                                 kin(i,j,k), eps(i,j,k), mut(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_Physics3D
!*******************************************************************************************************
!******** 物理量出力 (三次元)									********
!*******************************************************************************************************
subroutine Output_Physics3D( &
&            strname, ext, is, ie, js, je, ks, ke, rho, u, v, w, ps, ts, mu, kin, eps, mut )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je, ks, ke
 real     , intent(in) :: rho(is:ie, js:je, ks:ke), &
 &                        u(is:ie, js:je, ks:ke), v(is:ie, js:je, ks:ke), w(is:ie, js:je, ks:ke), &
 &                        ps(is:ie, js:je, ks:ke), ts(is:ie, js:je, ks:ke), mu(is:ie, js:je, ks:ke), &
 &                        kin(is:ie, js:je, ks:ke), eps(is:ie, js:je, ks:ke), mut(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, j, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) (((rho(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((  u(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((  v(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((  w(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) ((( ps(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) ((( ts(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) ((( mu(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((kin(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((eps(i,j,k), i = is, ie), j = js, je), k = ks, ke)
    write(1) (((mut(i,j,k), i = is, ie), j = js, je), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     write(1, '(e16.8e3,9(x,e16.8e3))') rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), &
     &                                  ps(i,j,k), ts(i,j,k), mu(i,j,k), &
     &                                  kin(i,j,k), eps(i,j,k), mut(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_Physics3D
!*******************************************************************************************************
!******** 一次元配列ファイル入力（実数）                                                        ********
!*******************************************************************************************************
subroutine Input_ArrayReal1D( &
&            strname, ext, is, ie, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie
 real     , intent(out) :: q(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) (q(i), i = is, ie)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do i = is, ie
     read(1, '(e16.8e3)') q(i)
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_ArrayReal1D
!*******************************************************************************************************
!******** 一次元配列ファイル出力（実数）                                                        ********
!*******************************************************************************************************
subroutine Output_ArrayReal1D( &
&            strname, ext, is, ie, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie
 real     , intent(in) :: q(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) (q(i), i = is, ie)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do i = is, ie
     write(1, '(e16.8e3)') q(i)
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_ArrayReal1D
!*******************************************************************************************************
!******** 二次元配列ファイル入力（実数）                                                        ********
!*******************************************************************************************************
subroutine Input_ArrayReal2D( &
&            strname, ext, is, ie, js, je, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je
 real     , intent(out) :: q(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) ((q(i,j), i = is, ie), j = js, je)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do j = js, je
    do i = is, ie
     read(1, '(e16.8e3)') q(i,j)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_ArrayReal2D
!*******************************************************************************************************
!******** 二次元配列ファイル出力（実数）                                                        ********
!*******************************************************************************************************
subroutine Output_ArrayReal2D( &
&            strname, ext, is, ie, js, je, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je
 real     , intent(in) :: q(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) ((q(i,j), i = is, ie), j = js, je)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do j = js, je
    do i = is, ie
     write(1, '(e16.8e3)') q(i,j)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_ArrayReal2D
!*******************************************************************************************************
!******** 三次元配列ファイル入力（実数）                                                        ********
!*******************************************************************************************************
subroutine Input_ArrayReal3D( &
&            strname, ext, is, ie, js, je, ks, ke, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 real     , intent(out) :: q(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) (((q(i,j,k), i = is, ie), j = js, je), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     read(1, '(e16.8e3)') q(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_ArrayReal3D
!*******************************************************************************************************
!******** 三次元配列ファイル出力（実数）                                                        ********
!*******************************************************************************************************
subroutine Output_ArrayReal3D( &
&            strname, ext, is, ie, js, je, ks, ke, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je, ks, ke
 real     , intent(in) :: q(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) (((q(i,j,k), i = is, ie), j = js, je), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     write(1, '(e16.8e3)') q(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_ArrayReal3D
!*******************************************************************************************************
!******** 一次元配列ファイル入力（整数）                                                        ********
!*******************************************************************************************************
subroutine Input_ArrayInt1D( &
&            strname, ext, is, ie, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie
 integer  , intent(out) :: q(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) (q(i), i = is, ie)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do i = is, ie
     read(1, '(i16)') q(i)
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_ArrayInt1D
!*******************************************************************************************************
!******** 一次元配列ファイル出力（整数）                                                        ********
!*******************************************************************************************************
subroutine Output_ArrayInt1D( &
&            strname, ext, is, ie, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie
 integer  , intent(in) :: q(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) (q(i), i = is, ie)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do i = is, ie
     write(1, '(i16)') q(i)
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_ArrayInt1D
!*******************************************************************************************************
!******** 二次元配列ファイル入力（整数）                                                        ********
!*******************************************************************************************************
subroutine Input_ArrayInt2D( &
&            strname, ext, is, ie, js, je, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je
 integer  , intent(out) :: q(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) ((q(i,j), i = is, ie), j = js, je)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do j = js, je
    do i = is, ie
     read(1, '(i16)') q(i,j)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_ArrayInt2D
!*******************************************************************************************************
!******** 二次元配列ファイル出力（整数）                                                        ********
!*******************************************************************************************************
subroutine Output_ArrayInt2D( &
&            strname, ext, is, ie, js, je, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je
 integer  , intent(in) :: q(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) ((q(i,j), i = is, ie), j = js, je)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do j = js, je
    do i = is, ie
     write(1, '(i16)') q(i,j)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_ArrayInt2D
!*******************************************************************************************************
!******** 三次元配列ファイル入力（整数）                                                        ********
!*******************************************************************************************************
subroutine Input_ArrayInt3D( &
&            strname, ext, is, ie, js, je, ks, ke, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 integer  , intent(out) :: q(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) (((q(i,j,k), i = is, ie), j = js, je), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     read(1, '(i16)') q(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_ArrayInt3D
!*******************************************************************************************************
!******** 三次元配列ファイル出力（整数）                                                        ********
!*******************************************************************************************************
subroutine Output_ArrayInt3D( &
&            strname, ext, is, ie, js, je, ks, ke, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, js, je, ks, ke
 integer  , intent(in) :: q(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) (((q(i,j,k), i = is, ie), j = js, je), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     write(1, '(i16)') q(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_ArrayInt3D
!*******************************************************************************************************
!******** 計算格子分割点入力 									********
!*******************************************************************************************************
subroutine Input_GridPoint( &
&            strname, ext, ip )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(out) :: ip
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
  read(1, '(i4)') ip
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_GridPoint
!*******************************************************************************************************
!******** 計算格子分割点出力 									********
!*******************************************************************************************************
subroutine Output_GridPoint( &
&            strname, ext, ip )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: ip
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
  write(1, '(i4)') ip
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_GridPoint
!*******************************************************************************************************
!******** 液滴投入条件入力 									********
!*******************************************************************************************************
subroutine Input_DropInCondition( &
&            strname, ext, a, u )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 real     , intent(out) :: a, u
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
  read(1, '(e16.8e3)') a
  read(1, '(e16.8e3)') u
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_DropInCondition
!*******************************************************************************************************
!******** 液滴投入条件出力 									********
!*******************************************************************************************************
subroutine Output_DropInCondition( &
&            strname, ext, a, u )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 real     , intent(in) :: a, u
 ! 処理開始 ********************************************************************************************
 open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
  write(1, '(e16.8e3)') a
  write(1, '(e16.8e3)') u
 close(1)
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_DropInCondition
!*******************************************************************************************************
!******** 液滴衝突データ入力 (三次元)								********
!*******************************************************************************************************
subroutine Input_Impingement3D( &
&            strname, ext, is, ie, ks, ke, Nimp, Uimp, Vimp, Wimp )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, ks, ke
 real     , intent(out) :: Nimp(is:ie, ks:ke)				! 衝突個数
 real     , intent(out) :: Uimp(is:ie, ks:ke), &
 &                         Vimp(is:ie, ks:ke), &
 &                         Wimp(is:ie, ks:ke)				! 衝突速度の和
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) ((Nimp(i,k), i = is, ie), k = ks, ke)
    read(1) ((Uimp(i,k), i = is, ie), k = ks, ke)
    read(1) ((Vimp(i,k), i = is, ie), k = ks, ke)
    read(1) ((Wimp(i,k), i = is, ie), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do k = ks, ke
    do i = is, ie
     read(1, '(e16.8e3,3(x,e16.8e3))') Nimp(i,k), Uimp(i,k), Vimp(i,k), Wimp(i,k)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_Impingement3D
!*******************************************************************************************************
!******** 液滴衝突データ出力 (三次元)								********
!*******************************************************************************************************
subroutine Output_Impingement3D( &
&            strname, ext, is, ie, ks, ke, Nimp, Uimp, Vimp, Wimp )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, ks, ke
 real     , intent(in) :: Nimp(is:ie, ks:ke)				! 衝突個数
 real     , intent(in) :: Uimp(is:ie, ks:ke), &
 &                        Vimp(is:ie, ks:ke), &
 &                        Wimp(is:ie, ks:ke)				! 衝突速度の和
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) ((Nimp(i,k), i = is, ie), k = ks, ke)
    write(1) ((Uimp(i,k), i = is, ie), k = ks, ke)
    write(1) ((Vimp(i,k), i = is, ie), k = ks, ke)
    write(1) ((Wimp(i,k), i = is, ie), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do i = is, ie
     write(1, '(e16.8e3,3(x,e16.8e3))') Nimp(i,k), Uimp(i,k), Vimp(i,k), Wimp(i,k)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_Impingement3D
!*******************************************************************************************************
!******** 収集効率入力 (三次元)									********
!*******************************************************************************************************
subroutine Input_CollectionEfficiency3D( &
&            strname, ext, is, ie, ks, ke, CE, u, v, w )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, ks, ke
 real     , intent(out) :: CE(is:ie, ks:ke)
 real     , intent(out) :: u(is:ie, ks:ke), v(is:ie, ks:ke), w(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) ((CE(i,k), i = is, ie), k = ks, ke)
    read(1) ((u (i,k), i = is, ie), k = ks, ke)
    read(1) ((v (i,k), i = is, ie), k = ks, ke)
    read(1) ((w (i,k), i = is, ie), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do k = ks, ke
    do i = is, ie
     read(1, '(e16.8e3,3(x,e16.8e3))') CE(i,k), u(i,k), v(i,k), w(i,k)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_CollectionEfficiency3D
!*******************************************************************************************************
!******** 収集効率出力 (三次元)									********
!*******************************************************************************************************
subroutine Output_CollectionEfficiency3D( &
&            strname, ext, is, ie, ks, ke, CE, u, v, w )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, ks, ke
 real     , intent(in) :: CE(is:ie, ks:ke)
 real     , intent(in) :: u(is:ie, ks:ke), v(is:ie, ks:ke), w(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) ((CE(i,k), i = is, ie), k = ks, ke)
    write(1) ((u (i,k), i = is, ie), k = ks, ke)
    write(1) ((v (i,k), i = is, ie), k = ks, ke)
    write(1) ((w (i,k), i = is, ie), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do i = is, ie
     write(1, '(e16.8e3,3(x,e16.8e3))') CE(i,k), u(i,k), v(i,k), w(i,k)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_CollectionEfficiency3D
!*******************************************************************************************************
!******** 氷層厚さ・温度入力 (二次元)								********
!*******************************************************************************************************
subroutine Input_IceThickTem2D( &
&            strname, ext, is, ie, f, Bi, dBi, Ti )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie
 integer  , intent(out) :: f  (is:ie)					! 着氷箇所のフラグ
 real     , intent(out) :: Bi (is:ie)					! 氷層厚さ (絶対量)
 real     , intent(out) :: dBi(is:ie)					! 氷層厚さ (変化量)
 real     , intent(out) :: Ti (is:ie)					! 氷層温度
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) (f  (i), i = is, ie)
    read(1) (Bi (i), i = is, ie)
    read(1) (dBi(i), i = is, ie)
    read(1) (Ti (i), i = is, ie)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do i = is, ie
     read(1, '(i1,3(x,e16.8e3))') f(i), Bi(i), dBi(i), Ti(i)
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_IceThickTem2D
!*******************************************************************************************************
!******** 氷層厚さ・温度出力 (二次元)								********
!*******************************************************************************************************
subroutine Output_IceThickTem2D( &
&            strname, ext, is, ie, f, Bi, dBi, Ti )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie
 integer  , intent(in) :: f  (is:ie)					! 着氷箇所のフラグ
 real     , intent(in) :: Bi (is:ie)					! 氷層厚さ (絶対量)
 real     , intent(in) :: dBi(is:ie)					! 氷層厚さ (変化量)
 real     , intent(in) :: Ti (is:ie)					! 氷層温度
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) (f  (i), i = is, ie)
    write(1) (Bi (i), i = is, ie)
    write(1) (dBi(i), i = is, ie)
    write(1) (Ti (i), i = is, ie)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do i = is, ie
     write(1, '(i1,3(x,e16.8e3))') f(i), Bi(i), dBi(i), Ti(i)
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_IceThickTem2D
!*******************************************************************************************************
!******** 氷層厚さ・温度入力 (三次元)								********
!*******************************************************************************************************
subroutine Input_IceThickTem3D( &
&            strname, ext, is, ie, ks, ke, f, Bi, dBi, Ti )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, ks, ke
 integer  , intent(out) :: f  (is:ie, ks:ke)				! 着氷箇所のフラグ
 real     , intent(out) :: Bi (is:ie, ks:ke)				! 氷層厚さ (絶対量)
 real     , intent(out) :: dBi(is:ie, ks:ke)				! 氷層厚さ (変化量)
 real     , intent(out) :: Ti (is:ie, ks:ke)				! 氷層温度
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) ((f  (i,k), i = is, ie), k = ks, ke)
    read(1) ((Bi (i,k), i = is, ie), k = ks, ke)
    read(1) ((dBi(i,k), i = is, ie), k = ks, ke)
    read(1) ((Ti (i,k), i = is, ie), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do k = ks, ke
    do i = is, ie
     read(1, '(i1,3(x,e16.8e3))') f(i,k), Bi(i,k), dBi(i,k), Ti(i,k)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_IceThickTem3D
!*******************************************************************************************************
!******** 氷層厚さ・温度出力 (三次元)								********
!*******************************************************************************************************
subroutine Output_IceThickTem3D( &
&            strname, ext, is, ie, ks, ke, f, Bi, dBi, Ti )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, ks, ke
 integer  , intent(in) :: f (is:ie, ks:ke)				! 着氷箇所のフラグ
 real     , intent(in) :: Bi(is:ie, ks:ke)				! 氷層厚さ (絶対量)
 real     , intent(in) :: dBi(is:ie, ks:ke)				! 氷層厚さ (変化量)
 real     , intent(in) :: Ti(is:ie, ks:ke)				! 氷層温度
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) ((f (i,k), i = is, ie), k = ks, ke)
    write(1) ((Bi(i,k), i = is, ie), k = ks, ke)
    write(1) ((dBi(i,k), i = is, ie), k = ks, ke)
    write(1) ((Ti(i,k), i = is, ie), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do i = is, ie
     write(1, '(i1,3(x,e16.8e3))') f(i,k), Bi(i,k), dBi(i,k), Ti(i,k)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_IceThickTem3D
!*******************************************************************************************************
!******** Original Messinger モデルの変数入力 (二次元)						********
!*******************************************************************************************************
subroutine Input_OrgMessingerPara2D( &
&            strname, ext, is, ie, Mac, Mes, Mim, Min, Mout, Mre )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie
 real     , intent(out) :: Mac (is:ie)					! 堆積する氷の質量
 real     , intent(out) :: Mes (is:ie)					! 蒸発・昇華する質量
 real     , intent(out) :: Mim (is:ie)					! 衝突する質量
 real     , intent(out) :: Min (is:ie)					! Runback-in する質量
 real     , intent(out) :: Mout(is:ie)					! Runback-out する質量
 real     , intent(out) :: Mre (is:ie)					! 検査体積内に残る質量
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) (Mac   (i), i = is, ie)
    read(1) (Mes   (i), i = is, ie)
    read(1) (Mim   (i), i = is, ie)
    read(1) (Min   (i), i = is, ie)
    read(1) (Mout  (i), i = is, ie)
    read(1) (Mre   (i), i = is, ie)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do i = is, ie
     read(1, '(e16.8e3,5(x,e16.8e3))') Mac(i), Mes(i), Mim(i), Min(i), Mout(i), Mre(i)
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_OrgMessingerPara2D
!*******************************************************************************************************
!******** Original Messinger モデルの変数出力 (二次元)						********
!*******************************************************************************************************
subroutine Output_OrgMessingerPara2D( &
&            strname, ext, is, ie, Mac, Mes, Mim, Min, Mout, Mre )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie
 real     , intent(in) :: Mac (is:ie)					! 堆積する氷の質量
 real     , intent(in) :: Mes (is:ie)					! 蒸発・昇華する質量
 real     , intent(in) :: Mim (is:ie)					! 衝突する質量
 real     , intent(in) :: Min (is:ie)					! Runback-in する質量
 real     , intent(in) :: Mout(is:ie)					! Runback-out する質量
 real     , intent(in) :: Mre (is:ie)					! 検査体積内に残る質量
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) (Mac   (i), i = is, ie)
    write(1) (Mes   (i), i = is, ie)
    write(1) (Mim   (i), i = is, ie)
    write(1) (Min   (i), i = is, ie)
    write(1) (Mout  (i), i = is, ie)
    write(1) (Mre   (i), i = is, ie)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do i = is, ie
     write(1, '(e16.8e3,5(x,e16.8e3))') Mac(i), Mes(i), Mim(i), Min(i), Mout(i), Mre(i)
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_OrgMessingerPara2D
!*******************************************************************************************************
!******** Extended-Messinger モデルの変数入力 (二次元)						********
!*******************************************************************************************************
subroutine Input_ExtMessingerPara2D( &
&            strname, ext, is, ie, nPhase, Mes, Mim, Min, Mout, Mre )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie
 integer  , intent(out) :: nPhase(is:ie)				! 壁面の状態
 real     , intent(out) :: Mes   (is:ie)				! 蒸発・昇華する質量流束
 real     , intent(out) :: Mim   (is:ie)				! 衝突する質量流束
 real     , intent(out) :: Min   (is:ie)				! Runback-in する質量流束
 real     , intent(out) :: Mout  (is:ie)				! Runback-out する質量流束
 real     , intent(out) :: Mre   (is:ie)				! 検査体積内に残る質量流束
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) (nPhase(i), i = is, ie)
    read(1) (Mes   (i), i = is, ie)
    read(1) (Mim   (i), i = is, ie)
    read(1) (Min   (i), i = is, ie)
    read(1) (Mout  (i), i = is, ie)
    read(1) (Mre   (i), i = is, ie)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do i = is, ie
     read(1, '(i1,5(x,e16.8e3))') nPhase(i), Mes(i), Mim(i), Min(i), Mout(i), Mre(i)
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_ExtMessingerPara2D
!*******************************************************************************************************
!******** Extended-Messinger モデルの変数出力 (二次元)						********
!*******************************************************************************************************
subroutine Output_ExtMessingerPara2D( &
&            strname, ext, is, ie, nPhase, Mes, Mim, Min, Mout, Mre )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie
 integer  , intent(in) :: nPhase(is:ie)					! 壁面の状態
 real     , intent(in) :: Mes (is:ie)					! 蒸発・昇華する質量流束
 real     , intent(in) :: Mim (is:ie)					! 衝突する質量流束
 real     , intent(in) :: Min (is:ie)					! Runback-in する質量流束
 real     , intent(in) :: Mout(is:ie)					! Runback-out する質量流束
 real     , intent(in) :: Mre (is:ie)					! 検査体積内に残る質量流束
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) (nPhase(i), i = is, ie)
    write(1) (Mes   (i), i = is, ie)
    write(1) (Mim   (i), i = is, ie)
    write(1) (Min   (i), i = is, ie)
    write(1) (Mout  (i), i = is, ie)
    write(1) (Mre   (i), i = is, ie)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do i = is, ie
     write(1, '(i1,5(x,e16.8e3))') nPhase(i), Mes(i), Mim(i), Min(i), Mout(i), Mre(i)
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_ExtMessingerPara2D
!*******************************************************************************************************
!******** Original Messinger モデルの変数入力 (三次元)						********
!*******************************************************************************************************
subroutine Input_OrgMessingerPara3D( &
&            strname, ext, is, ie, ks, ke, Mac, Mes, Mim, Min, Mout, Mre )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, ks, ke
 real     , intent(out) :: Mac (is:ie, ks:ke)				! 堆積する氷の質量
 real     , intent(out) :: Mes (is:ie, ks:ke)				! 蒸発・昇華する質量
 real     , intent(out) :: Mim (is:ie, ks:ke)				! 衝突する質量
 real     , intent(out) :: Min (is:ie, ks:ke)				! Runback-in する質量
 real     , intent(out) :: Mout(is:ie, ks:ke)				! Runback-out する質量
 real     , intent(out) :: Mre (is:ie, ks:ke)				! 検査体積内に残る質量
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) ((Mac   (i,k), i = is, ie), k = ks, ke)
    read(1) ((Mes   (i,k), i = is, ie), k = ks, ke)
    read(1) ((Mim   (i,k), i = is, ie), k = ks, ke)
    read(1) ((Min   (i,k), i = is, ie), k = ks, ke)
    read(1) ((Mout  (i,k), i = is, ie), k = ks, ke)
    read(1) ((Mre   (i,k), i = is, ie), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do k = ks, ke
    do i = is, ie
     read(1, '(e16.8e3,5(x,e16.8e3))') Mac(i,k), Mes(i,k), Mim(i,k), Min(i,k), Mout(i,k), Mre(i,k)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_OrgMessingerPara3D
!*******************************************************************************************************
!******** Original Messinger モデルの変数出力 (三次元)						********
!*******************************************************************************************************
subroutine Output_OrgMessingerPara3D( &
&            strname, ext, is, ie, ks, ke, Mac, Mes, Mim, Min, Mout, Mre )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, ks, ke
 real     , intent(in) :: Mac (is:ie, ks:ke)				! 堆積する氷の質量
 real     , intent(in) :: Mes (is:ie, ks:ke)				! 蒸発・昇華する質量
 real     , intent(in) :: Mim (is:ie, ks:ke)				! 衝突する質量
 real     , intent(in) :: Min (is:ie, ks:ke)				! Runback-in する質量
 real     , intent(in) :: Mout(is:ie, ks:ke)				! Runback-out する質量
 real     , intent(in) :: Mre (is:ie, ks:ke)				! 検査体積内に残る質量
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) ((Mac   (i,k), i = is, ie), k = ks, ke)
    write(1) ((Mes   (i,k), i = is, ie), k = ks, ke)
    write(1) ((Mim   (i,k), i = is, ie), k = ks, ke)
    write(1) ((Min   (i,k), i = is, ie), k = ks, ke)
    write(1) ((Mout  (i,k), i = is, ie), k = ks, ke)
    write(1) ((Mre   (i,k), i = is, ie), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do i = is, ie
     write(1, '(e16.8e3,5(x,e16.8e3))') Mac(i,k), Mes(i,k), Mim(i,k), Min(i,k), Mout(i,k), Mre(i,k)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_OrgMessingerPara3D
!*******************************************************************************************************
!******** Extended-Messinger モデルの変数入力 (三次元)						********
!*******************************************************************************************************
subroutine Input_ExtMessingerPara3D( &
&            strname, ext, is, ie, ks, ke, nPhase, Mes, Mim, Min, Mout, Mre )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, ks, ke
 integer  , intent(out) :: nPhase(is:ie, ks:ke)				! 壁面の状態
 real     , intent(out) :: Mes   (is:ie, ks:ke)				! 蒸発・昇華する質量流束
 real     , intent(out) :: Mim   (is:ie, ks:ke)				! 衝突する質量流束
 real     , intent(out) :: Min   (is:ie, ks:ke)				! Runback-in する質量流束
 real     , intent(out) :: Mout  (is:ie, ks:ke)				! Runback-out する質量流束
 real     , intent(out) :: Mre   (is:ie, ks:ke)				! 検査体積内に残る質量流束
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) ((nPhase(i,k), i = is, ie), k = ks, ke)
    read(1) ((Mes   (i,k), i = is, ie), k = ks, ke)
    read(1) ((Mim   (i,k), i = is, ie), k = ks, ke)
    read(1) ((Min   (i,k), i = is, ie), k = ks, ke)
    read(1) ((Mout  (i,k), i = is, ie), k = ks, ke)
    read(1) ((Mre   (i,k), i = is, ie), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do k = ks, ke
    do i = is, ie
     read(1, '(i1,5(x,e16.8e3))') nPhase(i,k), Mes(i,k), Mim(i,k), Min(i,k), Mout(i,k), Mre(i,k)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_ExtMessingerPara3D
!*******************************************************************************************************
!******** Extended-Messinger モデルの変数出力 (三次元)						********
!*******************************************************************************************************
subroutine Output_ExtMessingerPara3D( &
&            strname, ext, is, ie, ks, ke, nPhase, Mes, Mim, Min, Mout, Mre )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, ks, ke
 integer  , intent(in) :: nPhase(is:ie, ks:ke)				! 壁面の状態
 real     , intent(in) :: Mes (is:ie, ks:ke)				! 蒸発・昇華する質量流束
 real     , intent(in) :: Mim (is:ie, ks:ke)				! 衝突する質量流束
 real     , intent(in) :: Min (is:ie, ks:ke)				! Runback-in する質量流束
 real     , intent(in) :: Mout(is:ie, ks:ke)				! Runback-out する質量流束
 real     , intent(in) :: Mre (is:ie, ks:ke)				! 検査体積内に残る質量流束
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: i, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) ((nPhase(i,k), i = is, ie), k = ks, ke)
    write(1) ((Mes   (i,k), i = is, ie), k = ks, ke)
    write(1) ((Mim   (i,k), i = is, ie), k = ks, ke)
    write(1) ((Min   (i,k), i = is, ie), k = ks, ke)
    write(1) ((Mout  (i,k), i = is, ie), k = ks, ke)
    write(1) ((Mre   (i,k), i = is, ie), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do i = is, ie
     write(1, '(i1,5(x,e16.8e3))') nPhase(i,k), Mes(i,k), Mim(i,k), Min(i,k), Mout(i,k), Mre(i,k)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_ExtMessingerPara3D
!*******************************************************************************************************
!******** 着氷翼表面座標入力  (二次元)								********
!*******************************************************************************************************
subroutine Input_IceBladeSurface2D( &
&            strname, ext, is, ie, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie
 real     , intent(out) :: x(is:ie), y(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) (x(i), i = is, ie)
    read(1) (y(i), i = is, ie)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do i = is, ie
     read(1, '(e16.8e3,x,e16.8e3)') x(i), y(i)
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_IceBladeSurface2D
!*******************************************************************************************************
!******** 着氷翼表面座標出力  (二次元)								********
!*******************************************************************************************************
subroutine Output_IceBladeSurface2D( &
&            strname, ext, is, ie, x, y )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie
 real     , intent(in) :: x(is:ie), y(is:ie)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) (x(i), i = is, ie)
    write(1) (y(i), i = is, ie)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do i = is, ie
     write(1, '(e16.8e3,x,e16.8e3)') x(i), y(i)
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_IceBladeSurface2D
!*******************************************************************************************************
!******** 着氷翼表面座標入力  (三次元)								********
!*******************************************************************************************************
subroutine Input_IceBladeSurface3D( &
&            strname, ext, is, ie, ks, ke, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strname*(*), ext*4
 integer  , intent(in)  :: is, ie, ks, ke
 real     , intent(out) :: x(is:ie, ks:ke), y(is:ie, ks:ke), z(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'old')
    read(1) ((x(i,k), i = is, ie), k = ks, ke)
    read(1) ((y(i,k), i = is, ie), k = ks, ke)
    read(1) ((z(i,k), i = is, ie), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'old')
    do k = ks, ke
    do i = is, ie
     read(1, '(e16.8e3,2(x,e16.8e3))') x(i,k), y(i,k), z(i,k)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Input_IceBladeSurface3D
!*******************************************************************************************************
!******** 着氷翼表面座標出力  (三次元)								********
!*******************************************************************************************************
subroutine Output_IceBladeSurface3D( &
&            strname, ext, is, ie, ks, ke, x, y, z )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in) :: strname*(*), ext*4
 integer  , intent(in) :: is, ie, ks, ke
 real     , intent(in) :: x(is:ie, ks:ke), y(is:ie, ks:ke), z(is:ie, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, k
 ! 処理開始 ********************************************************************************************
 select case(ext)
  ! Binary ---------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) ((x(i,k), i = is, ie), k = ks, ke)
    write(1) ((y(i,k), i = is, ie), k = ks, ke)
    write(1) ((z(i,k), i = is, ie), k = ks, ke)
   close(1)
  ! Aschii ---------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do i = is, ie
     write(1, '(e16.8e3,2(x,e16.8e3))') x(i,k), y(i,k), z(i,k)
    enddo
    enddo
   close(1)
 end select
 ! 処理終了 ********************************************************************************************
 return
end subroutine Output_IceBladeSurface3D
! 定義終了 *********************************************************************************************
end module Package_FileIO
