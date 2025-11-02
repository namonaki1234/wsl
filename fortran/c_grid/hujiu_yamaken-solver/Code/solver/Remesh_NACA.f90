!*******************************************************************************************************
!*******************************************************************************************************
!******** ���X��̊i�q�č\���v���O����								********
!******** (NACA���C�O�����CC-type�C�d���i�q�@)							********
!********					      2013.02.07  PROGRAMMED BY RYOSUKE HAYASHI ********
!******** �� �_�u���z�[���ɑΉ�									********
!********					      2013.05.05     UPDATED BY RYOSUKE HAYASHI ********
!********					      2013.07.22     UPDATED BY RYOSUKE HAYASHI ********
!*******************************************************************************************************
!*******************************************************************************************************
program Remesh_NACA
 ! ���W���[���錾 **************************************************************************************
 use Package_NACA
 use Package_FileIO
 use Package_Equation
 use Package_Grid
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �t�@�C�����ݒ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: ViewGrdFile *  8 = 'ViewGrid'
 ! �Ǐ��萔 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, parameter :: mRef = 2
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: kRef
 ! �჌�C�m���Y���^�i�q�̏ꍇ,�i�q��؂邱�Ƃ�����ł���Ɨ\�z����邽��,
 ! �X�C�b�`�ŃX�e�b�v��i�߂邩���߂�(0:�i�߂Ȃ�,1:�i�߂�)
 ! 0�Ő؂�邩������,��肭������1�Ŗ{��
 integer :: swi_IceStep = 1
 ! �����J�n ********************************************************************************************
 ! ���؎����P�[�X�I�� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(a)') "<< Exp. Case Selection >>"
 call SelectExpCase
 ! �����ݒ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') '<< Initial Setting >>'
 call InitialSetting
 ! �X���[�W���O ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Smoothing >>"
 call SmoothIce
 ! ���X��̊i�q ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') "<< Icing Grid Geometry >>"
 call IcingGrid
 ! ���X��̗���ꏉ���l ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(/,a)') '<< Initial Flow Condition >>'
 call InitialFlow
 ! �����葱�� ******************************************************************************************
 stop
contains
!*******************************************************************************************************
!******** ���؎����P�[�X�I�� 									********
!*******************************************************************************************************
subroutine SelectExpCase
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character :: fname * 20
 ! �����J�n ********************************************************************************************
 ! �v�Z�����t�@�C������
 call Input_CalSetting( trim(ND_CalSetFile) // strtxt )
 Span = Span * lRef
 ! �f�B���N�g���ݒ�
 if( IceStep == 0 ) then
   GrdInDir     = bckdir // 'grid//clean//'
   FlwCalInDir  = bckdir // 'flow//cal//clean//'
   FlwIniDir    = bckdir // 'flow//initial//clean//'
  else
   GrdInDir     = bckdir // 'grid//icing//'
   FlwCalInDir  = bckdir // 'flow//cal//icing//'
   FlwIniDir    = bckdir // 'flow//initial//icing//'
 endif
 IceCalInDir  = bckdir // 'icing//cal//'
 GrdOutDir    = bckdir // 'grid//icing//'
 FlwCalOutDir = bckdir // 'flow//initial//icing//'
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
 ! �����I�� ********************************************************************************************
 return
end subroutine SelectExpCase
!*******************************************************************************************************
!********* �����ݒ�										********
!*******************************************************************************************************
subroutine InitialSetting
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: m, i, k
integer :: j
 integer   :: j0
 character :: fname * 20
 ! �����J�n ********************************************************************************************
 ! �u���b�N���ݒ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Set_BlockName
 ! �\���̃������m�� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( Flw(ms:me), Ice(ms:me) )
 ! �i�q�𑜓x ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_Resolution1D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(IceRslFile), strtxt, &
  &      Ice(m)%is, Ice(m)%ie )
  call Input_Resolution3D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(RslFile), strtxt, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke )
! enddo
 ! C �^�i�q�����_ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! m = 1
 m = me
 call Input_CtypeGridPoint( &
 &      trim(GrdInDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
 &      Flw(m)%i1, Flw(m)%i2, Flw(m)%i3 )
 Ice(m)%i1 = Flw(m)%i1
 Ice(m)%i2 = Flw(m)%i2
 Ice(m)%i3 = Flw(m)%i3
 ! �������m�� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  allocate( Flw(m)%x  ( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%y  ( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%z  ( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%f  ( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%xix( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%xiy( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%xiz( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%etx( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%ety( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%etz( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%zex( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%zey( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%zez( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%jac( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
  &         Flw(m)%qh ( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke, ls: le) )
  allocate( Ice(m)%x  ( Ice(m)%is: Ice(m)%ie ), &
  &         Ice(m)%y  ( Ice(m)%is: Ice(m)%ie ), &
  &         Ice(m)%z  ( Ice(m)%is: Ice(m)%ie ), &
  &         Ice(m)%f  ( Ice(m)%is: Ice(m)%ie ), &
  &         Ice(m)%Bi ( Ice(m)%is: Ice(m)%ie ), &
  &         Ice(m)%dBi( Ice(m)%is: Ice(m)%ie ), &
  &         Ice(m)%Ti ( Ice(m)%is: Ice(m)%ie ) )
! enddo
 ! �i�q���W ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_Grid3D( &
  &      trim(GrdInDir) // trim(BlkName(m)) // trim(GrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
! enddo
 ! �X�w���� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
!  if(m == 1) cycle
  call Input_IceThickTem2D( &
  &      trim(IceCalInDir) // trim(BlkName(m)) // trim(IceThickTemFile), strdat, &
  &      Ice(m)%is, Ice(m)%ie, &
  &      Ice(m)%f, Ice(m)%Bi, Ice(m)%dBi, Ice(m)%Ti )
  call Input_IceBladeSurface2D( &
  &      trim(IceCalInDir) // trim(BlkName(m)) // trim(IceBladeFile), strdat, &
  &      Ice(m)%is, Ice(m)%ie, &
  &      Ice(m)%x, Ice(m)%y )
! enddo
 ! ����ꏉ���� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_Flux3D( &
  &      trim(FlwCalInDir) // trim(BlkName(m)) // trim(ND_FlxFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      ls, le, Flw(m)%qh )
! enddo
 ! �Ώۃu���b�N ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 m    = mRef
 kRef = Flw(m)%ks + int( 0.5 * (Flw(m)%ke - Flw(m)%ks) )
 ! ���X���E�ʒu�T�� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do i = Ice(m)%i1, Ice(m)%i3 !Flw(m)%is, Flw(m)%ie
  if(Ice(m)%f(i) /= 0) then
    Flw(m)%i1 = i
    exit
  endif
 enddo
 do i = Ice(m)%i3, Ice(m)%i1, -1 !Flw(m)%ie, Flw(m)%is, -1
  if(Ice(m)%f(i) /= 0) then
    Flw(m)%i2 = i
    exit
  endif
 enddo
 if( Flw(m)%i1 < Flw(m)%is + 5 .or. Flw(m)%i2 > Flw(m)%ie - 5 ) then
   write(*, '(a)') '!!!!! Error : Icing limit point !!!!!'
   write(*, '(2i4)') Flw(m)%i1, Flw(m)%i2
   stop
 endif
  ! �����I�� ********************************************************************************************
 return
end subroutine InitialSetting
!*******************************************************************************************************
!******** �X�w�����̃X���[�W���O 								********
!*******************************************************************************************************
subroutine SmoothIce
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �Ǐ��萔 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, parameter :: nSmooth = 5 !20 !5
 ! �Ǐ��萔 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: GrdCleanDir * 20 = '../grid/clean/'
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: Cx(:,:,:), Cy(:,:,:), Cz(:,:,:)
 real   , pointer :: Bi0(:)
 real   , pointer :: dBi0(:)
 integer :: i, k, m, n
 integer :: j0, iHorn
 integer :: iHorn0, iHorn1, iHorn2
 real    :: ax, ay, az, bx, by, bz, nx, ny, nz, na
 logical :: DoubleHorn
 ! �����J�n ********************************************************************************************
 m = mRef; k = kRef; DoubleHorn = .false.
 ! �������m�� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( Cx( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
 &         Cy( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ), &
 &         Cz( Flw(m)%is: Flw(m)%ie, Flw(m)%js: Flw(m)%je, Flw(m)%ks: Flw(m)%ke ) )
 allocate( Bi0(Ice(m)%is: Ice(m)%ie), dBi0(Ice(m)%is: Ice(m)%ie) )
 ! ���X�O�����W ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call Input_Grid3D( &
 &      trim(GrdCleanDir) // trim(BlkName(m)) // trim(GrdFile), strbin, &
 &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
 &      Cx, Cy, Cz )
 ! 1�X�e�b�v�O�̕X�w���� +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Bi0(:) = Ice(m)%Bi(:) - Ice(m)%dBi(:)
 ! �z�[���ʒu�̒T�� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 write(*, '(a,i3,a,i4)') '* Ice limit index (i) = ', Flw(m)%i1, ',', Flw(m)%i2
 do i = Flw(m)%i1, Flw(m)%i2
  if( Ice(m)%dBi(i) >= maxval(Ice(m)%dBi(:)) ) then
    iHorn1 = i
    exit
  endif
 enddo
 do i = iHorn1 - 1, Flw(m)%i1, - 1
  if( Ice(m)%dBi(i) - Ice(m)%dBi(i-1) < 0.0 ) then
    iHorn0 = i
    DoubleHorn = .true.
    iHorn2 = iHorn1
    exit
  endif
 enddo
 if(DoubleHorn) then
   do i = iHorn0, Flw(m)%i1, - 1
    if( Ice(m)%dBi(i) - Ice(m)%dBi(i+1) < 0.0 ) exit
   enddo
   iHorn1 = i
 endif
 if(DoubleHorn) then
   write(*, '(a,i3)') '* Horn-1 index (i) = ', iHorn1
   write(*, '(a,i3)') '* Horn-2 index (i) = ', iHorn2
  else
   write(*, '(a,i3)') '* Horn-1 index (i) = ', iHorn1
 endif
 ! �X���[�W���O ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do n = 1, nSmooth
  dBi0(:) = Ice(m)%dBi(:)
  do i = Flw(m)%i1, Flw(m)%i2
   Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
  enddo

!  ! Single Horn ----------------------------------------------------------------------------------------
!  if(.not. DoubleHorn) then
!    ! �I�[
!    if(Flw(m)%i2 - 1 <= iHorn1) then
!      do i = iHorn1 - 1, Flw(m)%i1, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    ! �n�[
!     else if(iHorn1 <= Flw(m)%i1 + 1) then
!      do i = iHorn1 + 1, Flw(m)%i2
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    ! ���̑�
!     else
!      do i = iHorn1 + 1, Flw(m)%i2
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn1 - 1, Flw(m)%i1, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    endif
!    write(*, '(i1,6(e16.8e3))') n, Ice(m)%dBi(iHorn1-2:iHorn1+2), sum(Ice(m)%dBi(:))
!  ! Double Horn ----------------------------------------------------------------------------------------
!   else
!    ! �I�[���n�[
!    if(iHorn1 <= Flw(m)%i1 + 1 .and. Flw(m)%i2 - 1 <= iHorn2) then
!      do i = iHorn1 + 1, iHorn0
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn2 - 1, iHorn0, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    ! �n�[
!     else if(iHorn1 <= Flw(m)%i1 + 1) then
!      do i = iHorn1 + 1, iHorn0
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn2 + 1, Flw(m)%i2
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn2 - 1, iHorn0, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    ! �I�[
!     else if(Flw(m)%i2 - 1 <= iHorn2) then
!      do i = iHorn2 - 1, iHorn0, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn1 + 1, iHorn0
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn1 - 1, Flw(m)%i1, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    ! ���̑�
!     else
!      do i = iHorn1 + 1, iHorn0
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn1 - 1, Flw(m)%i1, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn2 + 1, Flw(m)%i2
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!      do i = iHorn2 - 1, iHorn0, -1
!       Ice(m)%dBi(i) = ( dBi0(i-1) + dBi0(i) + dBi0(i+1) ) / 3.0
!      enddo
!    endif
!    write(*, '(i1,6(e16.8e3))') n, Ice(m)%dBi(iHorn1-2:iHorn1+2), sum(Ice(m)%dBi(:))
!    write(*, '(i1,6(e16.8e3))') n, Ice(m)%dBi(iHorn2-2:iHorn2+2), sum(Ice(m)%dBi(:))
!  endif
 enddo
 ! �X���[�W���O��̗����W ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 j0 = Flw(m)%js
 do i = Flw(m)%i1, Flw(m)%i2
  ! �P�ʖ@���x�N�g��
  if(i == Flw(m)%is) then
    ax = 1.0 * ( - Cx(i,j0,k) + Cx(i+1,j0,k) )
    ay = 1.0 * ( - Cy(i,j0,k) + Cy(i+1,j0,k) )
    az = 1.0 * ( - Cz(i,j0,k) + Cz(i+1,j0,k) )
   else if(i == Flw(m)%ie) then
    ax = 1.0 * ( - Cx(i-1,j0,k) + Cx(i,j0,k) )
    ay = 1.0 * ( - Cy(i-1,j0,k) + Cy(i,j0,k) )
    az = 1.0 * ( - Cz(i-1,j0,k) + Cz(i,j0,k) )
   else
    ax = 0.5 * ( - Cx(i-1,j0,k) + Cx(i+1,j0,k) )
    ay = 0.5 * ( - Cy(i-1,j0,k) + Cy(i+1,j0,k) )
    az = 0.5 * ( - Cz(i-1,j0,k) + Cz(i+1,j0,k) )
  endif
  bx = 0.5 * ( + Cx(i,j0,k-1) - Cx(i,j0,k+1) )
  by = 0.5 * ( + Cy(i,j0,k-1) - Cy(i,j0,k+1) )
  bz = 0.5 * ( + Cz(i,j0,k-1) - Cz(i,j0,k+1) )
  nx = ay * bz - az * by
  ny = az * bx - ax * bz
  nz = ax * by - ay * bx
  na = sqrt(nx**2 + ny**2 + nz**2)
  nx = +1.0 * nx / na
  ny = +1.0 * ny / na
  nz = +1.0 * nz / na
  ! ���@�������ɕω�
!  Ice(m)%dBi(i) = Ice(m)%Bi(i) - Bi0(i)
  Ice(m)%x(i)   = Flw(m)%x(i,j0,k) + nx * Ice(m)%dBi(i)
  Ice(m)%y(i)   = Flw(m)%y(i,j0,k) + ny * Ice(m)%dBi(i)
!  Ice(m)%x(i)   = Cx(i,j0,k) + nx * Ice(m)%Bi(i)
!  Ice(m)%y(i)   = Cy(i,j0,k) + ny * Ice(m)%Bi(i)
 enddo
 ! �t�@�C���o�� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 if(swi_IceStep .eq. 1) then
! do m = ms, me
 m = me
!  if(m == 1) cycle
  call Output_IceThickTem3D( &
  &      trim(IceCalInDir) // trim(BlkName(m)) // trim(IceThickTemFile), strdat, &
  &      Ice(m)%is, Ice(m)%ie, Flw(m)%ks, Flw(m)%ke, &
  &      Ice(m)%f, Ice(m)%Bi, Ice(m)%dBi, Ice(m)%Ti )
! enddo
 end if
 deallocate(Cx, Cy, Cz, Bi0)
 ! �����I�� ********************************************************************************************
 return
end subroutine SmoothIce
!*******************************************************************************************************
!********* ���X��̊i�q�č\��									********
!*******************************************************************************************************
subroutine IcingGrid
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: m, i
 integer   :: jp
 character :: fname1 * 20, fname2 * 20
 character :: fname * 20
 ! �����J�n ********************************************************************************************
 ! �Ώۃu���b�N ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 m  = mRef
 jp = Flw(m)%je
 ! �i�q���� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call HtypeGridIceLE( &
! &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, Flw(m)%i1, Flw(m)%i2, &
! &      Ice(m)%x, Ice(m)%y, Flw(m)%x(:, jp, kRef), Flw(m)%y(:, jp, kRef), &
! &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
 call CtypeGridBlade( &
 &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
 &      Ice(m)%i1, Ice(m)%i2, Ice(m)%i3, &
 &      Ice(m)%x, Ice(m)%y, Flw(m)%x(:, jp, kRef), Flw(m)%y(:, jp, kRef), &
 &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
 ! �t�@�C���o�� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 if(swi_IceStep .eq. 1) then
 write(fname, '(a, 2(i2.2,a))') 'IceStep', IceStep, 'of', IceStepMax, '_'
! do m = ms, me
 m = me
!  select case(m)
!   case(1)
!    call Output_CtypeGridPoint( &
!    &      trim(GrdOutDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
!    &      Flw(m)%i1, Flw(m)%i2, Flw(m)%i3 )
!   case(2)
!    call Output_IceLimitPoint( &
!    &      trim(GrdOutDir) // trim(BlkName(m)) // trim(IceLimitPointFile) // strtxt, &
!    &      Flw(m)%i1, Flw(m)%i2 )
!  end select
  call Output_CtypeGridPoint( &
  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(CtypePointFile) // strtxt, &
  &      Ice(m)%i1, Ice(m)%i2, Ice(m)%i3 )
  call Output_IceLimitPoint( &
  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(IceLimitPointFile) // strtxt, &
  &      Flw(m)%i1, Flw(m)%i2 )

  call Output_Resolution1D( &
  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(IceRslFile), strtxt, &
  &      Ice(m)%is, Ice(m)%ie )
  call Output_Resolution3D( &
  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(RslFile), strtxt, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke )
  call Output_Grid3D( &
  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(GrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
  call Output_Grid3D( &
  &      trim(GrdOutDir) // trim(BlkName(m)) // trim(fname) // trim(GrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
! enddo
 end if
 ! ���� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  Flw(m)%f(:,:,:) = m
  call MakeMAVSFile3D( &
  &      trim(GrdOutDir), trim(BlkName(m)) // trim(ViewGrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%f , Flw(m)%x , Flw(m)%y , Flw(m)%z )
  call OutputPara_bin( &
   &     trim(GrdOutDir), trim(BlkName(m)) // trim(ViewGrdFile), &
   &     Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
   &     Flw(m)%x, Flw(m)%y, Flw(m)%z )
! enddo
 ! �����I�� ********************************************************************************************
 return
end subroutine IcingGrid
!*******************************************************************************************************
!********* ���O���t�߂̒��X�v�Z�p�̊i�q (H-type)						********
!*******************************************************************************************************
subroutine HtypeGridIceLE( &
&            is, ie, js, je, ks, ke, i1, i2, xi, yi, xe, ye, x, y, z )
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �����ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, js, je, ks, ke
 integer, intent(in)  :: i1, i2						! ���X���E�ʒu
 real   , intent(in)  :: xi(is:ie), yi(is:ie)				! �������E
 real   , intent(in)  :: xe(is:ie), ye(is:ie)				! �O�����E
 real   , intent(out) :: x(is:ie, js:je, ks:ke), &
 &                       y(is:ie, js:je, ks:ke), &
 &                       z(is:ie, js:je, ks:ke)
 ! �Ǐ��萔 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: rs1 = 1.0 * 1.0e-3				! ���\�ʋ��E�i�q��
 real   , parameter :: re1 = 1.0 * 1.0e-3				! ���������E�i�q��
 real   , parameter :: MGN = 1.5e-0					! �ȉ~-�o�Ȃ̏d�݂̋��e�͈�
 real   , parameter :: Rsd = 1.0e-5					! �ȉ~-�o�Ȃ̎�������l
 real   , parameter :: rs2 = 5.0 * 1.0e-1				! ���\�ʋ��E�i�q��
 real   , parameter :: re2 = 7.0 * 1.0e-2				! ���������E�i�q��
 real   , parameter :: tb1 = 5.0 * 1.0e-2				! �������E���𐫂̃p�����[�^
 real   , parameter :: tb2 = 2.0 * 1.0e-1				! �O�����E���𐫂̃p�����[�^
 ! �����J�n ********************************************************************************************
 ! �������E�l ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call InitialBoundary( &
&       is, ie, js, je, xi, yi, xe, ye, x(:, :, ks), y(:, :, ks) )
 ! �������E ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call HtypeSideBoundary( &
 &      is, ie, js, je, rs2, re2, x(:, :, ks), y(:, :, ks) )
 ! Transfinite ��� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call TransfiniteInterpolation( &
 &      is, ie, js, je, x(:, :, ks), y(:, :, ks) )
 ! �ȉ~-�o�Ȍ^�Δ����@ +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call HtypeGenerationEHPDE( &
 &      is, i1, i2, ie, js, je, rs1, re1, MGN, Rsd, x(:, :, ks), y(:, :, ks) )
! ! �񋫊E�@ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! call GenerationTwoBoundary( &
! &      is, ie, js, je, rs2, re2, tb1, tb2, x(:, :, ks), y(:, :, ks) )
 ! �O������ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call ThreeDimensionalized( &
 &      is, ie, js, je, ks, ke, span, x, y, z )
 ! �����I�� ********************************************************************************************
 return
end subroutine HtypeGridIceLE

!*******************************************************************************************************
!********* ������̊i�q (C-type)								********
!*******************************************************************************************************
subroutine CtypeGridBlade( &
&            is, ie, js, je, ks, ke, i1, i2, i3, xi, yi, xe, ye, x, y, z )
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �����ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, js, je, ks, ke
 integer, intent(in)  :: i1, i2, i3
 real   , intent(in)  :: xi(i1:i3), yi(i1:i3)				! �������E
 real   , intent(in)  :: xe(is:ie), ye(is:ie)				! �O�����E
 real   , intent(inout) :: x(is:ie, js:je, ks:ke), &
 &                         y(is:ie, js:je, ks:ke), &
 &                         z(is:ie, js:je, ks:ke)
 ! �Ǐ��萔 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: rs5 = 7.0 * 1.0e-2 !3.28649582292613 * 1.0e-3	! �������E�i�q��
 real   , parameter :: re5 = 1.0 * 1.0e+2				! �O�����E�i�q��
 real   , parameter :: tb1 = 1.0e-2 !2.6579 * 1.0e-1			! �������E���𐫂̃p�����[�^
 real   , parameter :: tb2 = 7.0 * 1.0e0				! �O�����E���𐫂̃p�����[�^
 ! �����J�n ********************************************************************************************
 ! �������E�l ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call InitialBoundaryCtype( &
 &      is, ie, i1, i3, js, je, xi, yi, xe, ye, x(:, :, ks), y(:, :, ks) )
 ! �񋫊E�@ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call GenerationTwoBoundary( &
 &      is, ie, js, je, rs5, re5, tb1, tb2, x(:, :, ks), y(:, :, ks) )
 ! �O������ ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 call ThreeDimensionalized( &
 &      is, ie, js, je, ks, ke, span, x, y, z )
 ! �����I�� ********************************************************************************************
 return
end subroutine CtypeGridBlade

!*******************************************************************************************************
!******** �������E�l										********
!*******************************************************************************************************
subroutine InitialBoundary( &
&            is, ie, js, je, xi, yi, xe, ye, x, y )
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �����ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, js, je
 real   , intent(in)  :: xi(is:ie), yi(is:ie)				! �������E
 real   , intent(in)  :: xe(is:ie), ye(is:ie)				! �O�����E
 real   , intent(out) :: x(is:ie, js:je), y(is:ie, js:je)
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! �����J�n ********************************************************************************************
 do i = is, ie
  x(i,js) = xi(i)
  y(i,js) = yi(i)
  x(i,je) = xe(i)
  y(i,je) = ye(i)
 enddo
 ! �����I�� ********************************************************************************************
 return
end subroutine InitialBoundary
!*******************************************************************************************************
!******** �������E�l(C�^�i�q�p)									********
!*******************************************************************************************************
subroutine InitialBoundaryCtype( &
&            is, ie, i1, i3, js, je, xi, yi, xe, ye, x, y )
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �����ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, js, je
 integer, intent(in)  :: i1, i3
 real   , intent(in)  :: xi(i1:i3), yi(i1:i3)				! �������E
 real   , intent(in)  :: xe(is:ie), ye(is:ie)				! �O�����E
 real   , intent(out) :: x(is:ie, js:je), y(is:ie, js:je)
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i
 ! �����J�n ********************************************************************************************
 do i = i1, i3
  x(i,js) = xi(i)
  y(i,js) = yi(i)
  x(i,je) = xe(i)
  y(i,je) = ye(i)
 enddo
 ! �����I�� ********************************************************************************************
 return
end subroutine InitialBoundaryCtype

!*******************************************************************************************************
!******** �������E (H-type)									********
!*******************************************************************************************************
subroutine HtypeSideBoundary( &
&            is, ie, js, je, rs, re, x, y )
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �����ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je
 real   , intent(in)    :: rs, re
 real   , intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: t(:)
 integer :: j
 real    :: rr
 ! �����J�n ********************************************************************************************
 allocate( t(js:je) )
 call GeometricInterpolation( rs / real(je-js+1), je-js+1, t, rr )
 do j = js + 1, je - 1
!  x(is,j) = x(is,js)
!  x(ie,j) = x(ie,js)
  x(is,j) = ( x(is,je) - x(is,js) ) * t(j) + x(is,js)
  x(ie,j) = ( x(ie,je) - x(ie,js) ) * t(j) + x(ie,js)
  y(is,j) = ( y(is,je) - y(is,js) ) * t(j) + y(is,js)
  y(ie,j) = ( y(ie,je) - y(ie,js) ) * t(j) + y(ie,js)
 enddo
 deallocate(t)
 ! �����J�n ********************************************************************************************
 return
end subroutine HtypeSideBoundary
!*******************************************************************************************************
!******** Transfinite ���									********
!*******************************************************************************************************
subroutine TransfiniteInterpolation( &
&            is, ie, js, je, x, y )
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �����ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je
 real   , intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: alp1(:), alp2(:), bet1(:), bet2(:)
 integer :: i, j
 ! �����J�n ********************************************************************************************
 ! �������m�� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( alp1(is:ie), alp2(is:ie), bet1(js:je), bet2(js:je) ) 
 ! Transfinite ��� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! �}��ϐ�
 do i = is, ie
  alp1(i) = real(ie-i) / real(ie)
  alp2(i) = 1.0 - alp1(i)
 enddo
 do j = js, je
  bet1(j) = real(je-j) / real(je)
  bet2(j) = 1.0 - bet1(j)
 enddo
 ! ���
 call Transfinite2D(is, ie, js, je, alp1, alp2, bet1, bet2, x)
 call Transfinite2D(is, ie, js, je, alp1, alp2, bet1, bet2, y)
 ! ��������� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 deallocate( alp1, alp2, bet1, bet2 )
 ! �����I�� ********************************************************************************************
 return
end subroutine TransfiniteInterpolation
!*******************************************************************************************************
!******** �ȉ~-�o�Ȍ^�Δ����������@�Ɋ�Â��i�q���� (H-type)					********
!*******************************************************************************************************
subroutine HtypeGenerationEHPDE( &
&            is, i1, i2, ie, js, je, rs, re, margin, resi, x, y )
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �����ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, i1, i2, ie, js, je
 real   , intent(in)    :: rs, re
 real   , intent(in)    :: margin					! �ȉ~-�o�Ȍ^�d�݊֐����e�덷
 real   , intent(in)    :: resi						! �v�Z���[�v��������l
 real   , intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
 ! �Ǐ��萔 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , parameter :: omg  = 1.0					! �ɘa�W��
 integer, parameter :: nmax = 100000					! �ő�v�Z��
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: Cs(:, :)
 real   , pointer :: cv1m(:, :), cv1p(:, :), cv2m(:, :), cv2p(:, :)
 real   , pointer :: dx1m(:, :), dx1p(:, :), dy1m(:, :), dy1p(:, :)
 real   , pointer :: dx2m(:, :), dx2p(:, :), dy2m(:, :), dy2p(:, :)
 integer :: i, j, n
 real    :: dmax
 real    :: a1, a2, a3, b1, b2, b3, c1, c2, c3, hjs1, hjs2, hje1, hje2
 ! �����J�n ********************************************************************************************
 ! �������m�� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 allocate( cs  (is:ie, js:je), &
 &         cv1m(is:ie, js:je), cv1p(is:ie, js:je), cv2m(is:ie, js:je), cv2p(is:ie, js:je), &
 &         dx1m(is:ie, js:je), dx1p(is:ie, js:je), dy1m(is:ie, js:je), dy1p(is:ie, js:je), &
 &         dx2m(is:ie, js:je), dx2p(is:ie, js:je), dy2m(is:ie, js:je), dy2p(is:ie, js:je) )
 ! �p�����[�^�ݒ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! �o�Ȍ^�d�݊֐� --------------------------------------------------------------------------------------
 do j = js+1, je-1
 do i = is+1, ie-1
  ! i - ����
  if(i < i1) then
    cv1m(i,j) = 0.0
   else if(i <= i2) then
    cv1m(i,j) = real(j) / real(je)**2 * real(i2-i) / real(i2-i1)
   else
    cv1m(i,j) = 0.0
  endif	
  ! i + ����
  if(i2 < i) then
    cv1p(i,j) = 0.0
   else if(i1 <= i) then
    cv1p(i,j) = real(j) / real(je)**2 * real(i-i1) / real(i2-i1)
   else
    cv1p(i,j) = 0.0
  endif
  ! j - ����
  cv2m(i,j) = (real(je-j) / real(je))**2
  ! j + ����
  cv2p(i,j) = 0.0
 enddo
 enddo
 ! �ȉ~�^�d�݊֐� --------------------------------------------------------------------------------------
 do j = js+1, je-1
 do i = is+1, ie-1
  ! �������E�@���x�N�g��
  a1 =  0.5 * ( - x(i-1,js) + x(i+1,js) )
  a2 =  0.5 * ( - y(i-1,js) + y(i+1,js) )
  a3 =  0.0
  b1 =  0.0
  b2 =  0.0
  b3 = -1.0
  c1 = a2 * b3 - a3 * b2
  c2 = a3 * b1 - a1 * b3
  hjs1 = c1 / sqrt(c1**2 + c2**2)
  hjs2 = c2 / sqrt(c1**2 + c2**2)
  ! �O�����E�@���x�N�g��
  a1 =  0.5 * ( - x(i-1,je) + x(i+1,je) )
  a2 =  0.5 * ( - y(i-1,je) + y(i+1,je) )
  a3 =  0.0
  b1 =  0.0
  b2 =  0.0
  b3 = -1.0
  c1 = a2 * b3 - a3 * b2
  c2 = a3 * b1 - a1 * b3
  hje1 = c1 / sqrt(c1**2 + c2**2)
  hje2 = c2 / sqrt(c1**2 + c2**2)
  ! i - ����
  if(i < i1) then
    dx1m(i,j) = x(i,je) - x(i-1,je)
    dy1m(i,j) = y(i,je) - y(i-1,je)
   else if(i <= i2) then
    dx1m(i,j) = ( -x(i-1,js) + x(i,js)  ) *   real(je-j) / real(je) &
    &         + ( -x(i-1,je) + x(i,je)  ) * ( real(j-js) / real(je) )**2
    dy1m(i,j) = ( -y(i-1,js) + y(i,js)  ) *   real(je-j) / real(je) &
    &         + ( -y(i-1,je) + y(i,je)  ) * ( real(j-js) / real(je) )**2
   else
    dx1m(i,j) = 0.0
    dy1m(i,j) = 0.0
  endif
  ! i + ����
  if(i2 < i) then
    dx1p(i,j) = x(i+1,je) - x(i,je)
    dy1p(i,j) = y(i+1,je) - y(i,je)
   else if(i1 <= i) then
    dx1p(i,j) = ( -x(i,js) + x(i+1,js) ) *   real(je-j) / real(je) &
    &         + ( -x(i,je) + x(i+1,je) ) * ( real(j-js) / real(je) )**2
    dy1p(i,j) = ( -y(i,js) + y(i+1,js) ) *   real(je-j) / real(je) &
    &         + ( -y(i,je) + y(i+1,je) ) * ( real(j-js) / real(je) )**2
   else
    dx1p(i,j) = 0.0
    dy1p(i,j) = 0.0
  endif
  ! j - ����
  dx2m(i,j) = hjs1 * rs
  dy2m(i,j) = hjs2 * rs
  ! j + ����
  dx2p(i,j) = hje1 * re
  dy2p(i,j) = hje2 * re
 enddo
 enddo
 ! �ȉ~-�o�Ȍ^�d�݊֐� ---------------------------------------------------------------------------------
 do j = js+1, je-1
 do i = is+1, ie-1
  Cs(i,j) = 1.0 - max( cv1m(i,j), cv1p(i,j), cv2m(i,j), cv2p(i,j) ) + margin
 enddo
 enddo
 ! �i�q�����v�Z���[�v ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 do n = 1, nmax
  ! �\���o�[�� -----------------------------------------------------------------------------------------
  call GridEHPDE2D( &
  &      omg, is, ie, js, je, &
  &      cs, cv1m, cv1p, cv2m, cv2p, dx1m, dy1m, dx1p, dy1p, dx2m, dy2m, dx2p, dy2p, &
  &      x, y, dmax )
  if( mod(n, 1000) == 0.0 ) write(*, "(a,i5,e16.8e3)") "* Elliptic-Hyperbolic calculation...", n, dmax
  if( dmax < resi ) exit
  ! ���E���� -------------------------------------------------------------------------------------------
   do j = js + 1, je - 1
    a1 = x(is+1,j) - x(is+1,j-1)
    a2 = y(is+1,j) - y(is+1,j-1)
    x(is,j) = x(is,j-1) + a1
    y(is,j) = y(is,j-1) + a2
    a1 = x(ie-1,j) - x(ie-1,j-1)
    a2 = y(ie-1,j) - y(ie-1,j-1)
    x(ie,j) = x(ie,j-1) + a1
    y(ie,j) = y(ie,j-1) + a2
   enddo
 enddo
 if( dmax > resi ) then
  write(*, '(a)') "!!!!! Elliptic-Hyperbolic calculation error !!!!!"
  stop
 endif
 ! �����I�� ********************************************************************************************
 return
end subroutine HtypeGenerationEHPDE
!*******************************************************************************************************
!******** �񋫊E�@�Ɋ�Â��i�q����								********
!*******************************************************************************************************
subroutine GenerationTwoBoundary( &
&            is, ie, js, je, rs, re, t1, t2, x, y )
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �����ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je
 real   , intent(in)    :: rs, re, t1, t2
 real   , intent(inout) :: x(is:ie, js:je), y(is:ie, js:je)
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , pointer :: etabar(:)
 real    :: rr
 ! �����J�n ********************************************************************************************
 ! �������m��
 allocate( etabar(js:je) )
 ! �}��ϐ�
 call GeometricInterpolationRemesh( rs / real(je-js+1), je-js+1, etabar, rr )
 ! �񋫊E�@
 call TwoBoundaryMethod2D( &
 &      is, ie, js, je, t1, t2, etabar, x, y )
 ! ���������
 deallocate(etabar)
 ! �����I�� ********************************************************************************************
 return
end subroutine GenerationTwoBoundary
!***********************************************************************
!**** ���䋉���ɂ��ꎟ����Ԋ֐�                                  ****
!**** ����a�A����r�̓��䐔��̘a(���䋉��)�ɂ��A                  ****
!**** ���0 <= x <= 1��n�̓_��z�u����                            ****
!***********************************************************************
SUBROUTINE GeometricInterpolationRemesh(a, n, x, r)
  ! �ϐ��錾 ***********************************************************
  IMPLICIT NONE
  ! �����ϐ� +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: a
  INTEGER, INTENT(IN)  :: n
  REAL,    INTENT(OUT) :: x(n)
  REAL,    INTENT(OUT) :: r
  ! �Ǐ��ϐ� +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i
  LOGICAL :: Check
  ! �����J�n ***********************************************************
  ! ��������� +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  r = FUNCR(n - 1, a, 1.0)
  !r is calculated with Mathematica FindRoot
  ! ���䋉���ɂ���Ԋ֐����v�Z +++++++++++++++++++++++++++++++++++++++
  x(1) = 0.0
  x(n) = 1.0
  Check = .FALSE.
  IF(r .NE. 1.0) THEN
    ! ���䋉�����v�Z ---------------------------------------------------
    DO i = 2, n - 1
      x(i) = a * (1.0 - r**(i - 1)) / (1.0 - r)
    ENDDO
    ! �P���������� -----------------------------------------------------
    DO i = 2, n - 1
      IF(x(i-1) .GE. x(i) .OR. x(i) .GE. x(i+1)) THEN
        WRITE(*, '(A)') 'GeometricInterpolation -> No Monotone Error'
        Check = .TRUE.
        EXIT
      ENDIF
    ENDDO
  ELSE
    Check = .TRUE.
  ENDIF
  ! ���䂪1�̏ꍇ���P�������ۏ؂���Ă��Ȃ��ꍇ�͓��Ԋu�̌��ʂ�Ԃ� ++++
  IF(Check) THEN
    DO i = 2, n - 1
      x(i) = REAL(i - 1) / REAL(n - 1)
    ENDDO
  ENDIF
  ! �����I�� ***********************************************************
  RETURN
! ��`�I�� *************************************************************
END SUBROUTINE GeometricInterpolationRemesh
!*******************************************************************************************************
!******** �O������ (�X�p�������ω��Ȃ�)								********
!*******************************************************************************************************
subroutine ThreeDimensionalized( &
 &      is, ie, js, je, ks, ke, dom, x, y, z )
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �����ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)    :: is, ie, js, je, ks, ke
 real   , intent(in)    :: dom
 real   , intent(inout) :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke)
 real   , intent(out)   :: z(is:ie, js:je, ks:ke)
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! �����J�n ********************************************************************************************
 do k = ks, ke
 do j = js, je
 do i = is, ie
  x(i,j,k) = x(i,j,ks)
  y(i,j,k) = y(i,j,ks)
  z(i,j,k) = -0.5 * dom + dom * real(k) / real(ke)
 enddo
 enddo
 enddo
 ! �����I�� ********************************************************************************************
 return
end subroutine ThreeDimensionalized
!*******************************************************************************************************
!******** ���X��̗���ꏉ���l									********
!*******************************************************************************************************
subroutine InitialFlow
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: m
 integer :: i,j,k,l
 real,allocatable :: jac2(:,:,:)
 ! �����J�n ********************************************************************************************
 ! �������� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call NondimensionalizedCoord3D( &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
! enddo
 ! ���g���b�N�X ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! do m = ms, me
 m = me
  call Input_Metrics3D( &
  &      trim(FlwIniDir) // trim(BlkName(m)) // trim(ND_MetFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%jac, &
  &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
  &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
  &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez )

 allocate(jac2(Flw(m)%is:Flw(m)%ie,Flw(m)%js:Flw(m)%je,Flw(m)%ks:Flw(m)%ke))
 do k = Flw(m)%ks,Flw(m)%ke
 do j = Flw(m)%js,Flw(m)%je
 do i = Flw(m)%is,Flw(m)%ie
  jac2(i,j,k) = Flw(m)%jac(i,j,k)
 end do
 end do
 end do

 m = me
  call Metrics3D( &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z, &
  &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
  &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
  &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez, &
  &      Flw(m)%jac )

 do k = Flw(m)%ks,Flw(m)%ke
 do j = Flw(m)%js,Flw(m)%je
 do i = Flw(m)%is,Flw(m)%ie
 do l = ls,le
  Flw(m)%qh(i,j,k,l) = Flw(m)%qh(i,j,k,l) * jac2(i,j,k) / Flw(m)%jac(i,j,k)
 end do
 end do
 end do
 end do

 deallocate(jac2)

! enddo
 ! �t�@�C���o�� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 if(swi_IceStep .eq. 1) then
 ! �i�q���W --------------------------------------------------------------------------------------------
! do m = ms, me
 m = me
  call Output_Grid3D( &
  &      trim(FlwCalOutDir) // trim(BlkName(m)) // trim(ND_GrdFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%x, Flw(m)%y, Flw(m)%z )
! enddo
 ! ���g���b�N�X ----------------------------------------------------------------------------------------
! do m = ms, me
 m = me
  call Output_Metrics3D( &
  &      trim(FlwCalOutDir) // trim(BlkName(m)) // trim(ND_MetFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      Flw(m)%jac, &
  &      Flw(m)%xix, Flw(m)%xiy, Flw(m)%xiz, &
  &      Flw(m)%etx, Flw(m)%ety, Flw(m)%etz, &
  &      Flw(m)%zex, Flw(m)%zey, Flw(m)%zez )
! enddo
 ! �����֐� --------------------------------------------------------------------------------------------
! do m = ms, me
 m = me
  call Output_Flux3D( &
  &      trim(FlwCalOutDir) // trim(BlkName(m)) // trim(ND_IniFlxFile), strbin, &
  &      Flw(m)%is, Flw(m)%ie, Flw(m)%js, Flw(m)%je, Flw(m)%ks, Flw(m)%ke, &
  &      ls, le, Flw(m)%qh )
! enddo
 end if
 ! �v�Z�����t�@�C���o�� --------------------------------------------------------------------------------
 nCount = 0; nDrop = 0
 IceStep = IceStep + 1
 Span = Span / lRef
 if(swi_IceStep .eq. 1) then
  call Output_CalSetting( trim(ND_CalSetFile) // strtxt )
 end if
 ! �����I�� ********************************************************************************************
 return
end subroutine InitialFlow
! ��`�I�� *********************************************************************************************
end program Remesh_NACA
!*******************************************************************************************************
!******** MicroAVS�t�@�C���쐬�i�O�����j  		                                        ********
!*******************************************************************************************************
subroutine MakeMAVSFile3D( &
&            strdir, strname, ext, is, ie, js, je, ks, ke, f, x, y, z )
 ! �ϐ��錾 ********************************************************************************************
 implicit none
 ! �����ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, intent(in)  :: strdir*(*), strname*(*), ext*4
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 integer  , intent(in)  :: f(is:ie, js:je, ks:ke)
 real     , intent(in)  :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke), z(is:ie, js:je, ks:ke)
 ! �Ǐ��萔 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 character, parameter :: strbin * 4 = '.bin'
 character, parameter :: strfld * 4 = '.fld'
 integer  , parameter :: ndim   = 3
 integer  , parameter :: nspace = 3
 integer  , parameter :: veclen = 1
 ! �Ǐ��ϐ� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer   :: nrkind, nskip
 real      :: r
 integer   :: i, j, k, n
 ! �����J�n ********************************************************************************************
 ! �w�b�_�t�@�C���o�� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 open(1, file = trim(strdir) // trim(strname) // trim(strfld), form = 'formatted')
  write(1, '(a)')     '# AVS field file'
  write(1, '(a,i1)') 'ndim   = ', ndim
  write(1, '(a,i4)') 'dim1   = ', ie - is + 1
  write(1, '(a,i4)') 'dim2   = ', je - js + 1
  write(1, '(a,i4)') 'dim3   = ', ke - ks + 1
  write(1, '(a)')    'label  = flag'
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
 ! �f�[�^�t�@�C���o�� ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 select case(ext)
 ! Binary ----------------------------------------------------------------------------------------------
  case(strbin)
   open(1, file = trim(strdir) // trim(strname) // ext, form = 'unformatted', status = 'replace')
    write(1) real(f); write(1) x; write(1) y; write(1) z
    close(1)
   close(1)
 ! Aschii ----------------------------------------------------------------------------------------------
  case default
   open(1, file = trim(strdir) // trim(strname) // ext, form = 'formatted', status = 'replace')
    do k = ks, ke
    do j = js, je
    do i = is, ie
     write(1, '(e16.8e3, 3(x,e16.8e3))') real(f(i,j,k)), x(i,j,k), y(i,j,k), z(i,j,k)
    enddo
    enddo
    enddo
   close(1)
 end select
 ! �����I�� ********************************************************************************************
 return
end subroutine MakeMAVSFile3D

!*******************************************************************************************************
!******** vtk�t�@�C���o�̓T�u���[�`�� 								********
!*******************************************************************************************************
subroutine OutputPara_bin( &
&      strdir, strname, is, ie, js, je, ks, ke, &
&      x, y, z )
 implicit none
 !mainroutine_variable
 character, intent(in)  :: strdir*(*), strname*(*)
 integer  , intent(in)  :: is, ie, js, je, ks, ke
 real     , intent(in)  :: x(is:ie, js:je, ks:ke), y(is:ie, js:je, ks:ke), z(is:ie, js:je, ks:ke)
 !subroutine_variable
 character(len = 4), parameter		:: strvtk = '.vtk'
 character(len = 1), parameter		:: newline = char(10)
 character(len = 200)	:: strnum
 integer		:: i,j,k,n
 integer		:: ni,nj,nk
 integer		:: npoint

 npoint = (ie - is) * (je - js) * (ke - ks)

open(unit       = 1, &
     file       = trim(strdir)//trim(strname)//trim(strvtk), &
     form       = 'unformatted', &
     access     = 'stream', &          ! ← ここを stream に
     convert    = 'big_endian', &      ! GNU拡張（OK）。気になるなら -fconvert で全体指定でも可
     action     = 'write')

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
   write(1) x(ni,nj,nk)		,y(ni,nj,nk)		,z(ni,nj,nk)
   write(1) x(ni+1,nj,nk)	,y(ni+1,nj,nk)		,z(ni+1,nj,nk)
   write(1) x(ni+1,nj+1,nk)	,y(ni+1,nj+1,nk)	,z(ni+1,nj+1,nk)
   write(1) x(ni,nj+1,nk)	,y(ni,nj+1,nk)		,z(ni,nj+1,nk)
   write(1) x(ni,nj,nk+1)	,y(ni,nj,nk+1)		,z(ni,nj,nk+1)
   write(1) x(ni+1,nj,nk+1)	,y(ni+1,nj,nk+1)	,z(ni+1,nj,nk+1)
   write(1) x(ni+1,nj+1,nk+1)	,y(ni+1,nj+1,nk+1)	,z(ni+1,nj+1,nk+1)
   write(1) x(ni,nj+1,nk+1)	,y(ni,nj+1,nk+1)	,z(ni,nj+1,nk+1)
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
  write(1) 'SCALARS n float'//newline
  write(1) 'LOOKUP_TABLE default'//newline
  do n = 0, npoint-1
   write(1) real(n * 8 + 0)
   write(1) real(n * 8 + 1)
   write(1) real(n * 8 + 2)
   write(1) real(n * 8 + 3)
   write(1) real(n * 8 + 4)
   write(1) real(n * 8 + 5)
   write(1) real(n * 8 + 6)
   write(1) real(n * 8 + 7)
  end do
 close(1)

end subroutine OutputPara_bin
