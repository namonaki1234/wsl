!*******************************************************************************************************
!*******************************************************************************************************
!******** パッケージ型モジュール								********
!******** 計算結果検証用サブルーチン群								********
!********				              2012.05.24  PROGRAMMED BY RYOSUKE HAYASHI	********
!*******************************************************************************************************
!*******************************************************************************************************
module Package_Validation
 ! 変数宣言 ********************************************************************************************
 implicit none
 private
 ! サブルーチン宣言 ************************************************************************************
 public :: SetPhysics2DKEM, SetPhysics3DKEM
 public :: CalMachTtPt2D  , CalMachTtPt3D
 public :: CalVortex2D    , CalVortex3D
 public :: ViscosityCoefficient3D
 public :: EddyViscosityCoefficient3D
! 内部手続き *******************************************************************************************
contains
!*******************************************************************************************************
!******** 保存量流束から物理量を計算 (二次元)                                  			********
!*******************************************************************************************************
subroutine SetPhysics2DKEM( &
&            Rg, gamma, is, ie, js, je, ls, le, jac, qh, rho, u, v, p, t, kin, eps )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real,    intent(in)  :: Rg, gamma
 integer, intent(in)  :: is, ie, js, je
 integer, intent(in)  :: ls, le
 real,    intent(in)  :: jac(is:ie, js:je)
 real,    intent(in)  :: qh(is:ie, js:je, ls:le)
 real,    intent(out) :: rho(is:ie, js:je), u(is:ie, js:je), v(is:ie, js:je), &
 &                       p(is:ie, js:je), t(is:ie, js:je), kin(is:ie, js:je), eps(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j
 ! 処理開始 ********************************************************************************************
 ! 例外処理(計算対象が適合していない場合は終了) ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 if(ls .ne. 1 .or. le .ne. 6) then
   write(*,'(a)') 'error : ls or le irregal'
   stop
 endif
 ! 保存量流束から物理量を計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !$omp parallel do default(shared) private(i, j)
 do j = js, je
 do i = is, ie
  if(qh(i,j,1) .gt. 0.0) then
    rho(i,j) = jac(i,j) * qh(i,j,1)
    u(i,j)   = qh(i,j,2) / qh(i,j,1)
    v(i,j)   = qh(i,j,3) / qh(i,j,1)
    p(i,j)   = ( qh(i,j,4) - qh(i,j,5) - 0.5 * (qh(i,j,2)**2 + qh(i,j,3)**2) / qh(i,j,1) ) &
             * (gamma - 1.0) * jac(i,j)
    t(i,j)   = p(i,j) / (rho(i,j) * rg)
    kin(i,j) = qh(i,j,5) / qh(i,j,1)
    eps(i,j) = qh(i,j,6) / qh(i,j,1)
   else
    rho(i,j) = 0.0
    u(i,j)   = 0.0
    v(i,j)   = 0.0
    p(i,j)   = 0.0
    t(i,j)   = 0.0
    kin(i,j) = 0.0
    eps(i,j) = 0.0
  endif
 enddo
 enddo
 !$omp end parallel do
 ! 処理終了 ********************************************************************************************
 return
end subroutine SetPhysics2DKEM
!*******************************************************************************************************
!******** 保存量流束から物理量を計算 (三次元)                                  			********
!*******************************************************************************************************
subroutine SetPhysics3DKEM( &
&            Rg, gamma, is, ie, js, je, ks, ke, ls, le, jac, qh, rho, u, v, w, p, t, kin, eps )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real,    intent(in)  :: Rg, gamma
 integer, intent(in)  :: is, ie, js, je, ks, ke
 integer, intent(in)  :: ls, le
 real,    intent(in)  :: jac(is:ie, js:je, ks:ke)
 real,    intent(in)  :: qh(is:ie, js:je, ks:ke, ls:le)
 real,    intent(out) :: rho(is:ie, js:je, ks:ke), &
 &                       u(is:ie, js:je, ks:ke), v(is:ie, js:je, ks:ke), w(is:ie, js:je, ks:ke), &
 &                       p(is:ie, js:je, ks:ke), t(is:ie, js:je, ks:ke), &
 &                       kin(is:ie, js:je, ks:ke), eps(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 ! 例外処理 (計算対象が適合していない場合は終了) +++++++++++++++++++++++++++++++++++++++++++++++++++++++
 if(ls .ne. 1 .or. le .ne. 7) then
   write(*,'(a)') 'error : ls or le irregal'
   stop
 endif
 ! 保存量流束から物理量を計算 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !$omp parallel do default(shared) private(i, j, k)
 do k = ks, ke
 do j = js, je
 do i = is, ie
  if( qh(i,j,k,1) .gt. 0.0 ) then
    rho(i,j,k) = jac(i,j,k) * qh(i,j,k,1)
    u(i,j,k)   = qh(i,j,k,2) / qh(i,j,k,1)
    v(i,j,k)   = qh(i,j,k,3) / qh(i,j,k,1)
    w(i,j,k)   = qh(i,j,k,4) / qh(i,j,k,1)
    p(i,j,k)   = ( qh(i,j,k,5) - qh(i,j,k,6) &
    &            - 0.5 * ( qh(i,j,k,2)**2 + qh(i,j,k,3)**2 + qh(i,j,k,4)**2 ) / qh(i,j,k,1) ) &
    &          * (gamma - 1.0) * jac(i,j,k)
    t(i,j,k)   = p(i,j,k) / (rho(i,j,k) * rg)
    kin(i,j,k) = qh(i,j,k,6) / qh(i,j,k,1)
    eps(i,j,k) = qh(i,j,k,7) / qh(i,j,k,1)
   else
    rho(i,j,k) = 0.0
    u(i,j,k)   = 0.0
    v(i,j,k)   = 0.0
    w(i,j,k)   = 0.0
    p(i,j,k)   = 0.0
    t(i,j,k)   = 0.0
    kin(i,j,k) = 0.0
    eps(i,j,k) = 0.0
  endif
 enddo
 enddo
 enddo
 !$omp end parallel do
 ! 処理終了 ********************************************************************************************
 return
end subroutine SetPhysics3DKEM
!*******************************************************************************************************
!******** マッハ数，全温，全圧を計算 (二次元)                            			********
!*******************************************************************************************************
subroutine CalMachTtPt2D( &
&            Rg, gamma, is, ie, js, je, u, v, Ps, Ts, mac, Pt, Tt )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real,    intent(in)  :: rg, gamma
 integer, intent(in)  :: is, ie, js, je
 real,    intent(in)  :: u(is:ie, js:je), v(is:ie, js:je), Ps(is:ie, js:je), Ts(is:ie, js:je)
 real,    intent(out) :: mac(is:ie, js:je), Pt(is:ie, js:je), Tt(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j
 ! 処理開始 ********************************************************************************************
 do j = js, je
 do i = is, ie
  if(Ts(i,j) .gt. 0.0) then
    mac(i,j) = sqrt( (u(i,j)**2 + v(i,j)**2) / (gamma * Rg * Ts(i,j)) )
    Tt(i,j)   = Ts(i,j) &
    &         * (1.0 + (gamma - 1.0) * 0.5 * mac(i,j)**2)
    Pt(i,j)   = Ps(i,j) &
    &         * ( (1.0 + (gamma - 1.0) * 0.5 * mac(i,j)**2)**(gamma / (gamma - 1.0)) )
   else
    mac(i,j) = 0.0
    Pt(i,j)   = 0.0
    Tt(i,j)   = 0.0
  endif
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalMachTtPt2D
!*******************************************************************************************************
!******** マッハ数，全温，全圧 (三次元) 						        ********
!*******************************************************************************************************
subroutine CalMachTtPt3D( &
&            Rg, gamma,is, ie, js, je, ks, ke, u, v, w, Ps, Ts, mac, Pt, Tt )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real,    intent(in)  :: Rg, gamma
 integer, intent(in)  :: is, ie, js, je, ks, ke
 real,    intent(in)  :: u(is:ie, js:je, ks:ke), v(is:ie, js:je, ks:ke), w(is:ie, js:je, ks:ke), &
 &                       Ps(is:ie, js:je, ks:ke), Ts(is:ie, js:je, ks:ke)
 real,    intent(out) :: mac(is:ie, js:je, ks:ke)
 real,    intent(out) :: Pt(is:ie, js:je, ks:ke), &
 &                       Tt(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 !$omp parallel do default(shared) private(i, j, k)
 do k = ks, ke
 do j = js, je
 do i = is, ie
  if(Ts(i,j,k) > 0.0) then
    mac(i,j,k) = sqrt( (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2) &
    &                   / (gamma * rg * Ts(i,j,k)) )
    Tt(i,j,k) = Ts(i,j,k) &
    &         * ( (1.0 + (gamma - 1.0) * 0.5 * mac(i,j,k)**2) )
    Pt(i,j,k) = Ps(i,j,k) &
    &         * ( (1.0 + (gamma - 1.0) * 0.5 * mac(i,j,k)**2)**(gamma / (gamma - 1.0)) )
   else
    mac(i,j,k) = 0.0
    Pt(i,j,k)   = 0.0
    Tt(i,j,k)   = 0.0
  endif
 enddo
 enddo
 enddo
 !$omp end parallel do
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalMachTtPt3D
!*******************************************************************************************************
!******** 渦構造の可視化 (二次元)                                        			********
!******** 渦度ベクトル, エンストロフィ, ヘリシティ, 度勾配テンソルの第二不変量			********
!*******************************************************************************************************
subroutine CalVortex2D( &
&            OmgZ, is, ie, js, je, rho, u, v, xix, xiy, etx, ety, oz, phi, psi, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real,    intent(in)  :: OmgZ
 integer, intent(in)  :: is, ie, js, je
 real,    intent(in)  :: rho(is:ie, js:je)
 real,    intent(in)  :: u(is:ie, js:je), v(is:ie, js:je)
 real,    intent(in)  :: xix(is:ie, js:je), xiy(is:ie, js:je), &
 &                       etx(is:ie, js:je), ety(is:ie, js:je)
 real,    intent(out) :: oz(is:ie, js:je)
 real,    intent(out) :: phi(is:ie, js:je)
 real,    intent(out) :: psi(is:ie, js:je)
 real,    intent(out) :: q(is:ie, js:je)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j
 real    :: uxi, vxi, uet, vet
 real    :: ux, vx, uy, vy
 real    :: s11, s22, s33, s12, ss
 real    :: o12
 ! 処理開始 ********************************************************************************************
 do j = js, je
 do i = is, ie
  if(rho(i,j) .gt. 0.0) then
    ! 速度勾配
    if(i .eq. is) then
      if(rho(i+1,j) .gt. 0.0 .and. rho(i+2,j) .gt. 0.0) then
        uxi = -1.5 * u(i,j) + 2.0 * u(i+1,j) - 0.5 * u(i+2,j)
        vxi = -1.5 * v(i,j) + 2.0 * v(i+1,j) - 0.5 * v(i+2,j)
       else
        uxi = 0.0
        vxi = 0.0
      endif
     else if(i .eq. ie) then
      if(rho(i-2,j) .gt. 0.0 .and. rho(i-1,j) .gt. 0.0) then
        uxi = 0.5 * u(i-2,j) - 2.0 * u(i-1,j) + 1.5 * u(i,j)
        vxi = 0.5 * v(i-2,j) - 2.0 * v(i-1,j) + 1.5 * v(i,j)
       else
        uxi = 0.0
        vxi = 0.0
      endif
     else
      if(rho(i-1,j) .gt. 0.0 .and. rho(i+1,j) .gt. 0.0) then
        uxi = -0.5 * (u(i-1,j) - u(i+1,j))
        vxi = -0.5 * (v(i-1,j) - v(i+1,j))
       else
        uxi = 0.0
        vxi = 0.0
      endif
    endif
    if(j .eq. js)then
      if(rho(i,j+1) .gt. 0.0 .and. rho(i,j+2) .gt. 0.0) then
        uet =-1.5 * u(i,j) + 2.0 * u(i,j+1) - 0.5 * u(i,j+2)
        vet =-1.5 * v(i,j) + 2.0 * v(i,j+1) - 0.5 * v(i,j+2)
       else
        uet = 0.0
        vet = 0.0
      endif
     else if(j .eq. je)then
      if(rho(i,j-2) .gt. 0.0 .and. rho(i,j-1) .gt. 0.0) then
        uet = 0.5 * u(i,j-2) - 2.0 * u(i,j-1) + 1.5 * u(i,j)
        vet = 0.5 * v(i,j-2) - 2.0 * v(i,j-1) + 1.5 * v(i,j)
       else
        uet = 0.0
        vet = 0.0
      endif
     else
      if(rho(i,j-1) .gt. 0.0 .and. rho(i,j+1) .gt. 0.0) then
        uet =-0.5 * (u(i,j-1) - u(i,j+1))
        vet =-0.5 * (v(i,j-1) - v(i,j+1))
       else
        uet = 0.0
        vet = 0.0
      endif
    endif
    ux = xix(i,j) * uxi + etx(i,j) * uet
    uy = xiy(i,j) * uxi + ety(i,j) * uet
    vx = xix(i,j) * vxi + etx(i,j) * vet
    vy = xiy(i,j) * vxi + ety(i,j) * vet
    ! 歪み速度テンソル, 角速度テンソル
    ss  = (ux + vy) / 3.0
    s11 = ux - ss
    s22 = vy - ss
    s33 =    - ss
    s12 = 0.5 * (uy + vx)
    o12 = 0.5 * (uy - vx) - OmgZ
    ! 渦度ベクトル
    oz(i,j)  =-2.0 * o12
    ! エンストロフィ
    phi(i,j) = 0.5 * oz(i,j)**2
    ! ヘリシティ
    psi(i,j) = 0.0
    ! 速度勾配の第二不変量
    q(i,j)   = o12 * o12 - s12 * s12 &
    &        - 0.5 * (s11 * s11 + s22 * s22 + s33 * s33)
   else
    oz(i,j)  = 0.0
    phi(i,j) = 0.0
    psi(i,j) = 0.0
    q(i,j)   = 0.0
  endif
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalVortex2D
!*******************************************************************************************************
!******** 渦構造の可視化 (三次元)                                        			********
!******** 渦度ベクトル, エンストロフィ, ヘリシティ, 速度勾配テンソルの第二不変量                ********
!*******************************************************************************************************
subroutine CalVortex3D( &
&            OmgX, OmgY, OmgZ, is, ie, js, je, ks, ke, rho, u, v, w, &
&            xix, xiy, xiz, etx, ety, etz, zex, zey, zez, ox, oy, oz, phi, psi, q )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real   , intent(in)  :: OmgX, OmgY, OmgZ
 integer, intent(in)  :: is, ie, js, je, ks, ke
 real   , intent(in)  :: rho(is:ie, js:je, ks:ke), &
 &                       u(is:ie, js:je, ks:ke), v(is:ie, js:je, ks:ke), w(is:ie, js:je, ks:ke)
 real   , intent(in)  :: xix(is:ie, js:je, ks:ke), xiy(is:ie, js:je, ks:ke), xiz(is:ie, js:je, ks:ke), &
 &                       etx(is:ie, js:je, ks:ke), ety(is:ie, js:je, ks:ke), etz(is:ie, js:je, ks:ke), &
 &                       zex(is:ie, js:je, ks:ke), zey(is:ie, js:je, ks:ke), zez(is:ie, js:je, ks:ke)
 real   , intent(out) :: ox(is:ie, js:je, ks:ke), oy(is:ie, js:je, ks:ke), oz(is:ie, js:je, ks:ke)
 real   , intent(out) :: phi(is:ie, js:je, ks:ke)
 real   , intent(out) :: psi(is:ie, js:je, ks:ke)
 real   , intent(out) :: q(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 real    :: uxi, vxi, wxi, uet, vet, wet, uze, vze, wze
 real    :: ux, vx, wx, uy, vy, wy, uz, vz, wz
 real    :: s11, s22, s33, s12, s23, s31, ss
 real    :: o12, o23, o31
 real    :: psiabs
 ! 処理開始 ********************************************************************************************
 do k = ks, ke
 do j = js, je
 do i = is, ie
  if(rho(i,j,k) .gt. 0.0) then
    ! 速度勾配
    if(i .eq. is)then
      if(rho(i+1,j,k) .gt. 0.0 .and. rho(i+2,j,k) .gt. 0.0) then
        uxi = -1.5 * u(i,j,k) + 2.0 * u(i+1,j,k) - 0.5 * u(i+2,j,k)
        vxi = -1.5 * v(i,j,k) + 2.0 * v(i+1,j,k) - 0.5 * v(i+2,j,k)
        wxi = -1.5 * w(i,j,k) + 2.0 * w(i+1,j,k) - 0.5 * w(i+2,j,k)
       else
        uxi = 0.0
        vxi = 0.0
        wxi = 0.0
      endif
     else if(i .eq. ie)then
      if(rho(i-2,j,k) .gt. 0.0 .and. rho(i-1,j,k) .gt. 0.0) then
        uxi = 0.5 * u(i-2,j,k) - 2.0 * u(i-1,j,k) + 1.5 * u(i,j,k)
        vxi = 0.5 * v(i-2,j,k) - 2.0 * v(i-1,j,k) + 1.5 * v(i,j,k)
        wxi = 0.5 * w(i-2,j,k) - 2.0 * w(i-1,j,k) + 1.5 * w(i,j,k)
       else
        uxi = 0.0
        vxi = 0.0
        wxi = 0.0
      endif
     else
      if(rho(i-1,j,k) .gt. 0.0 .and. rho(i+1,j,k) .gt. 0.0) then
        uxi = -0.5 * (u(i-1,j,k) - u(i+1,j,k))
        vxi = -0.5 * (v(i-1,j,k) - v(i+1,j,k))
        wxi = -0.5 * (w(i-1,j,k) - w(i+1,j,k))
       else
        uxi = 0.0
        vxi = 0.0
        wxi = 0.0
      endif
    endif
    if(j .eq. js)then
      if(rho(i,j+1,k) .gt. 0.0 .and. rho(i,j+2,k) .gt. 0.0) then
        uet = -1.5 * u(i,j,k) + 2.0 * u(i,j+1,k) - 0.5 * u(i,j+2,k)
        vet = -1.5 * v(i,j,k) + 2.0 * v(i,j+1,k) - 0.5 * v(i,j+2,k)
        wet = -1.5 * w(i,j,k) + 2.0 * w(i,j+1,k) - 0.5 * w(i,j+2,k)
       else
        uet = 0.0
        vet = 0.0
        wet = 0.0
      endif
     else if(j .eq. je)then
      if(rho(i,j-2,k) .gt. 0.0 .and. rho(i,j-1,k) .gt. 0.0) then
        uet = 0.5 * u(i,j-2,k) - 2.0 * u(i,j-1,k) + 1.5 * u(i,j,k)
        vet = 0.5 * v(i,j-2,k) - 2.0 * v(i,j-1,k) + 1.5 * v(i,j,k)
        wet = 0.5 * w(i,j-2,k) - 2.0 * w(i,j-1,k) + 1.5 * w(i,j,k)
      else
        uet = 0.0
        vet = 0.0
        wet = 0.0
      endif
     else
      if(rho(i,j-1,k) .gt. 0.0 .and. rho(i,j+1,k) .gt. 0.0) then
        uet = -0.5 * (u(i,j-1,k) - u(i,j+1,k))
        vet = -0.5 * (v(i,j-1,k) - v(i,j+1,k))
        wet = -0.5 * (w(i,j-1,k) - w(i,j+1,k))
       else
        uet = 0.0
        vet = 0.0
        wet = 0.0
      endif
    endif
    if(k .eq. ks)then
      if(rho(i,j,k+1) .gt. 0.0 .and. rho(i,j,k+2) .gt. 0.0) then
        uze = -1.5 * u(i,j,k) + 2.0 * u(i,j,k+1) - 0.5 * u(i,j,k+2)
        vze = -1.5 * v(i,j,k) + 2.0 * v(i,j,k+1) - 0.5 * v(i,j,k+2)
        wze = -1.5 * w(i,j,k) + 2.0 * w(i,j,k+1) - 0.5 * w(i,j,k+2)
       else
        uze = 0.0
        vze = 0.0
        wze = 0.0
      endif
     else if(k .eq. ke)then
      if(rho(i,j,k-2) .gt. 0.0 .and. rho(i,j,k-1) .gt. 0.0) then
        uze = 0.5 * u(i,j,k-2) - 2.0 * u(i,j,k-1) + 1.5 * u(i,j,k)
        vze = 0.5 * v(i,j,k-2) - 2.0 * v(i,j,k-1) + 1.5 * v(i,j,k)
        wze = 0.5 * w(i,j,k-2) - 2.0 * w(i,j,k-1) + 1.5 * w(i,j,k)
       else
        uze = 0.0
        vze = 0.0
        wze = 0.0
      endif
     else
      if(rho(i,j,k-1) .gt. 0.0 .and. rho(i,j,k+1) .gt. 0.0) then
        uze = -0.5 * (u(i,j,k-1) - u(i,j,k+1))
        vze = -0.5 * (v(i,j,k-1) - v(i,j,k+1))
        wze = -0.5 * (w(i,j,k-1) - w(i,j,k+1))
       else
        uze = 0.0
        vze = 0.0
        wze = 0.0
      endif
    endif
    ux = xix(i,j,k) * uxi + etx(i,j,k) * uet + zex(i,j,k) * uze
    uy = xiy(i,j,k) * uxi + ety(i,j,k) * uet + zey(i,j,k) * uze
    uz = xiz(i,j,k) * uxi + etz(i,j,k) * uet + zez(i,j,k) * uze
    vx = xix(i,j,k) * vxi + etx(i,j,k) * vet + zex(i,j,k) * vze
    vy = xiy(i,j,k) * vxi + ety(i,j,k) * vet + zey(i,j,k) * vze
    vz = xiz(i,j,k) * vxi + etz(i,j,k) * vet + zez(i,j,k) * vze
    wx = xix(i,j,k) * wxi + etx(i,j,k) * wet + zex(i,j,k) * wze
    wy = xiy(i,j,k) * wxi + ety(i,j,k) * wet + zey(i,j,k) * wze
    wz = xiz(i,j,k) * wxi + etz(i,j,k) * wet + zez(i,j,k) * wze
    ! 歪み速度テンソル, 角速度テンソル
    ss  = (ux + vy + wz) / 3.0
    s11 = ux - ss
    s22 = vy - ss
    s33 = wz - ss
    s12 = 0.5 * (uy + vx)
    s23 = 0.5 * (vz + wy)
    s31 = 0.5 * (wx + uz)
    o12 = 0.5 * (uy - vx) - OmgZ
    o23 = 0.5 * (vz - wy) - OmgX
    o31 = 0.5 * (wx - uz) - OmgY
    ! 渦度ベクトル
    ox(i,j,k)  = -2.0 * o23
    oy(i,j,k)  = -2.0 * o31
    oz(i,j,k)  = -2.0 * o12
    ! エンストロフィ
    phi(i,j,k) = 0.5 * (ox(i,j,k)**2 + oy(i,j,k)**2 + oz(i,j,k)**2)
    ! ヘリシティ
    psi(i,j,k) = u(i,j,k) * ox(i,j,k) &
    &          + v(i,j,k) * oy(i,j,k) &
    &          + w(i,j,k) * oz(i,j,k)
    ! 無次元ヘリシティ
    psiabs = sqrt( ( u(i,j,k)**2 +  v(i,j,k)**2 +  w(i,j,k)**2) &
    &            * (ox(i,j,k)**2 + oy(i,j,k)**2 + oz(i,j,k)**2) )
    if(psiabs .gt. 0.0) psi(i,j,k) = psi(i,j,k) / psiabs
    ! 速度勾配の第二不変量
    q(i,j,k)   = o12 * o12 - s12 * s12 &
    &          + o23 * o23 - s23 * s23 &
    &          + o31 * o31 - s31 * s31 &
    &          - 0.5 * (s11 * s11 + s22 * s22 + s33 * s33)
   else
    ox(i,j,k)  = 0.0
    oy(i,j,k)  = 0.0
    oz(i,j,k)  = 0.0
    phi(i,j,k) = 0.0
    psi(i,j,k) = 0.0
    q(i,j,k)   = 0.0
  endif
 enddo
 enddo
 enddo
 ! 処理終了 ********************************************************************************************
 return
end subroutine CalVortex3D
!*******************************************************************************************************
!******** 粘性係数計算                                                  			********
!******** 計算対象 : 単相, 三次元, 圧縮性                               			********
!*******************************************************************************************************
subroutine ViscosityCoefficient3D( &
&            muIn, TsIn, s1, is, ie, js, je, ks, ke, t, mu )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 real,    intent(in)  :: muIn, TsIn, s1
 integer, intent(in)  :: is, ie, js, je, ks, ke
 real,    intent(in)  :: t(is:ie, js:je, ks:ke)
 real,    intent(out) :: mu(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 !$omp parallel do default(shared) private(i, j, k)
 do k = ks, ke
 do j = js, je
 do i = is, ie
  if(t(i,j,k) .gt. 0.0) then
    mu(i,j,k) = muIn * (t(i,j,k) / TsIn)**1.5 &
    &          * (TsIn + s1) / (t(i,j,k) + s1)
   else
!    mu(i,j,k) = 0.0
    mu(i,j,k) = muIn
  endif
 enddo
 enddo
 enddo
 !$omp end parallel do
 ! 処理終了 ********************************************************************************************
 return
end subroutine ViscosityCoefficient3D
!*******************************************************************************************************
!******** 渦粘性係数計算									********
!*******************************************************************************************************
subroutine EddyViscosityCoefficient3D( &
&            is, ie, js, je, ks, ke, rho, kin, eps, cmu, mut )
 ! 変数宣言 ********************************************************************************************
 implicit none
 ! 引数変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer, intent(in)  :: is, ie, js, je, ks, ke
 real   , intent(in)  :: rho(is:ie, js:je, ks:ke)
 real   , intent(in)  :: kin(is:ie, js:je, ks:ke)
 real   , intent(in)  :: eps(is:ie, js:je, ks:ke)
 real   , intent(in)  :: cmu
 real   , intent(out) :: mut(is:ie, js:je, ks:ke)
 ! 局所変数 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 integer :: i, j, k
 ! 処理開始 ********************************************************************************************
 !$omp parallel do default(shared) private(i, j, k)
 do k = ks, ke
 do j = js, je
 do i = is, ie
  if( rho(i,j,k) > 0.0 .and. eps(i,j,k) > 0.0 ) then
    mut(i,j,k) = cmu * kin(i,j,k)**2 / eps(i,j,k)
  endif
 enddo
 enddo
 enddo
 !$ end parallel do
 ! 処理終了 ********************************************************************************************
 return
end subroutine EddyViscosityCoefficient3D
! 定義終了 *********************************************************************************************
end module Package_Validation



