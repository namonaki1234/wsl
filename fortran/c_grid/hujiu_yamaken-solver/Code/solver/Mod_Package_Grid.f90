!***********************************************************************
!***********************************************************************
!**** パッケージ型モジュール                                        ****
!**** 格子生成用サブルーチン群                                      ****
!****                         2010.02.04 PROGRAMED BY SUZUKI MASAYA ****
!****                         2011.10.15 UPDATED BY RYOSUKE HAYASHI ****
!***********************************************************************
!***********************************************************************
MODULE Package_Grid
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  PRIVATE
  ! サブルーチン宣言 ***************************************************
  ! 共有サブルーチン (正也さん作) ++++++++++++++++++++++++++++++++++++++
  PUBLIC :: CheckGrid2D, CheckGrid3D
  PUBLIC :: FUNCR, GeometricInterpolation, GeometricInterpolationInv
  PUBLIC :: VinokurInterpolation
  PUBLIC :: GridEPDE1D,  GridEPDE2D,  GridEPDE3D
  PUBLIC :: GridEHPDE1D, GridEHPDE2D, GridEHPDE3D
  PUBLIC :: TwoBoundaryMethod2D
  PUBLIC :: Transfinite1D, Transfinite2D, Transfinite3D
  PUBLIC :: Metrics2D, Metrics3D
  PUBLIC :: WallDistance2D, WallDistance3D
  PUBLIC :: WideSearch2D4Point, WideSearch3D8Point
  PUBLIC :: Interpolation2D3Point, Interpolation3D4Point
  PUBLIC :: Interpolation2D4Point, Interpolation3D8Point
  ! 共有サブルーチン (自作) ++++++++++++++++++++++++++++++++++++++++++++
  PUBLIC :: BladeInOut
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 二次元格子間隔のチェック                                      ****
!***********************************************************************
SUBROUTINE CheckGrid2D( &
&            IS, IE, JS, JE, X, Y &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN) :: IS, IE, JS, JE
  REAL,    INTENT(IN) :: X(IS:IE, JS:JE), Y(IS:IE, JS:JE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: RATE_MAX, RATE_MIN
  REAL    :: DL1, DL2, RATE
  ! 処理開始 ***********************************************************
  ! 格子点間隔の比に問題がないかを判定 +++++++++++++++++++++++++++++++++
  RATE_MAX = 1.2
  RATE_MIN = RATE_MAX**(-1.0)
  OPEN(16,FILE='Check.txt')
  ! I方向
  DO J=JS,JE
    DO I=IS+1,IE-1
      DL1 = SQRT( (X(I  ,J) - X(I-1,J))**2 &
      &         + (Y(I  ,J) - Y(I-1,J))**2 )
      DL2 = SQRT( (X(I+1,J) - X(I  ,J))**2 &
      &         + (Y(I+1,J) - Y(I  ,J))**2 )
      RATE = DL2 / DL1
      IF((RATE .GT. RATE_MAX) .OR. (RATE .LT. RATE_MIN)) THEN
        WRITE(16,'(A,F12.8,A,I4,A,I4,A)') &
        & 'Bad Rate : ', RATE, ' Grid ', I, '-',J, ' I Direction'
      ENDIF
    ENDDO
  ENDDO
  ! J方向
  DO I=IS,IE
    DO J=JS+1,JE-1
      DL1 = SQRT( (X(I,J  ) - X(I,J-1))**2 &
      &         + (Y(I,J  ) - Y(I,J-1))**2 )
      DL2 = SQRT( (X(I,J+1) - X(I,J  ))**2 &
      &         + (Y(I,J+1) - Y(I,J  ))**2 )
      RATE = DL2 / DL1
      IF((RATE .GT. RATE_MAX ) .OR. (RATE .LT. RATE_MIN)) THEN
        WRITE(16,'(A,F12.8,A,I4,A,I4,A)') &
        & 'Bad Rate : ', RATE, ' Grid ', I, '-', J, ' J Direction'
      ENDIF
    ENDDO
  ENDDO
  CLOSE(16)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE CheckGrid2D
!***********************************************************************
!**** 三次元格子間隔のチェック                                      ****
!***********************************************************************
SUBROUTINE CheckGrid3D( &
&            IS, IE, JS, JE, KS, KE, X, Y, Z &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN) :: IS, IE, JS, JE, KS, KE
  REAL,    INTENT(IN) :: X(IS:IE, JS:JE, KS:KE), &
  &                      Y(IS:IE, JS:JE, KS:KE), &
  &                      Z(IS:IE, JS:JE, KS:KE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: RATE_MAX, RATE_MIN
  REAL    :: DL1, DL2, RATE
  ! 処理開始 ***********************************************************
  ! 格子点間隔の比に問題がないかを判定 +++++++++++++++++++++++++++++++++
  RATE_MAX = 1.2
  RATE_MIN = RATE_MAX**(-1.0)
  OPEN(16,FILE='Check.txt')
  ! I方向
  DO K=KS,KE
  DO J=JS,JE
    DO I=IS+1,IE-1
      DL1 = SQRT( ( X(I  ,J,K)-X(I-1,J,K) )**2 &
      &         + ( Y(I  ,J,K)-Y(I-1,J,K) )**2 &
      &         + ( Z(I  ,J,K)-Z(I-1,J,K) )**2 )
      DL2 = SQRT( ( X(I+1,J,K)-X(I  ,J,K) )**2 &
      &         + ( Y(I+1,J,K)-Y(I  ,J,K) )**2 &
      &         + ( Z(I+1,J,K)-Z(I  ,J,K) )**2 )
      RATE = DL2/DL1
      IF( ( RATE .GT. RATE_MAX ).OR.( RATE .LT. RATE_MIN ) ) THEN
        WRITE(16,'(A,F12.8,A,I4,A,I4,A,I4,A)') &
        & 'Bad Rate : ',RATE,' Grid ',I,'-',J,'-',K, &
        & ' I Direction'
      ENDIF
    ENDDO
  ENDDO
  ENDDO
  ! J方向
  DO K=KS,KE
  DO I=IS,IE
    DO J=JS+1,JE-1
      DL1 = SQRT( ( X(I,J  ,K)-X(I,J-1,K) )**2 &
      &         + ( Y(I,J  ,K)-Y(I,J-1,K) )**2 &
      &         + ( Z(I,J  ,K)-Z(I,J-1,K) )**2 )
      DL2 = SQRT( ( X(I,J+1,K)-X(I,J  ,K) )**2 &
      &         + ( Y(I,J+1,K)-Y(I,J  ,K) )**2 &
      &         + ( Z(I,J+1,K)-Z(I,J  ,K) )**2 )
      RATE = DL2/DL1
      IF( ( RATE .GT. RATE_MAX ).OR.( RATE .LT. RATE_MIN ) ) THEN
        WRITE(16,'(A,F12.8,A,I4,A,I4,A,I4,A)') &
        & 'Bad Rate : ',RATE,' Grid ',I,'-',J,'-',K, &
        & ' J Direction'
      ENDIF
    ENDDO
  ENDDO
  ENDDO
  ! K方向
  DO J=JS,JE
  DO I=IS,IE
    DO K=KS+1,KE-1
      DL1 = SQRT( ( X(I,J,K  )-X(I,J,K-1) )**2 &
      &         + ( Y(I,J,K  )-Y(I,J,K-1) )**2 &
      &         + ( Z(I,J,K  )-Z(I,J,K-1) )**2 )
      DL2 = SQRT( ( X(I,J,K+1)-X(I,J,K  ) )**2 &
      &         + ( Y(I,J,K+1)-Y(I,J,K  ) )**2 &
      &         + ( Z(I,J,K+1)-Z(I,J,K  ) )**2 )
      RATE = DL2/DL1
      IF( ( RATE .GT. RATE_MAX ).OR.( RATE .LT. RATE_MIN ) ) THEN
        WRITE(16,'(A,F12.8,A,I4,A,I4,A,I4,A)') &
        & 'Bad Rate : ',RATE,' Grid ',I,'-',J,'-',K, &
        & ' K Direction'
      ENDIF
    ENDDO
  ENDDO
  ENDDO
  CLOSE(16)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE CheckGrid3D
!***********************************************************************
!**** 最終項数n、初項a、等比級数Sから公比FUNCRを計算する            ****
!***********************************************************************
REAL FUNCTION FUNCR(n, a, S)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, PARAMETER :: MMAX = 1000 !100
  REAL(8), PARAMETER :: RMIN = 0.5D+0, RMAX = 2.0D+0
  REAL(8), PARAMETER :: ZERO = 1.0D-20, RESMIN = 1.0D-8
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN) :: n
  REAL,    INTENT(IN) :: a, S
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: m
  REAL(8) :: S1, Sn
  REAL(8) :: r, dr, f, fr
  ! 処理開始 ***********************************************************
  ! 初期値 -------------------------------------------------------------
  S1 = DBLE(a)
  Sn = DBLE(S)
  IF(S1 * DBLE(n) .LT. Sn) THEN
   r = RMAX
  ELSEIF(S1 * DBLE(n) .GT. Sn) THEN
   r = RMIN
  ELSE
   FUNCR = 1.0
   RETURN
  END IF
  ! Newton-Raphson Method ----------------------------------------------
  DO m = 1, MMAX
    f  = S1 * (1.0D+0 - r**DBLE(n)) / (1.0D+0 - r) - Sn
    fr = S1 * ( - DBLE(n) * r**DBLE(n - 1) * (1.0D+0 - r) &
    &         + 1.0D+0 - r**DBLE(n) &
    &    ) / (1.0D+0 - r)**2
    IF(ABS(fr) .LE. ZERO) fr = SIGN(ZERO, fr)
    dr =-f / fr
    r  = r + dr
    r  = MAX(RMIN, MIN(RMAX, r))
    IF(ABS(dr) .LE. ABS(r) * RESMIN .OR. r .EQ. 1.0D+0) EXIT
  ENDDO
  FUNCR = REAL(r)
  ! WRITE(*,*) 'Common ratio of geometric progression = ', FUNCR
  ! 処理終了 ***********************************************************
  RETURN
END FUNCTION FUNCR
!***********************************************************************
!**** 等比級数による一次元補間関数                                  ****
!**** 初項a、公比rの等比数列の和(等比級数)により、                  ****
!**** 区間0 <= x <= 1にn個の点を配置する                            ****
!***********************************************************************
SUBROUTINE GeometricInterpolation(a, n, x, r)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: a
  INTEGER, INTENT(IN)  :: n
  REAL,    INTENT(OUT) :: x(n)
  REAL,    INTENT(OUT) :: r
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i
  LOGICAL :: Check
  ! 処理開始 ***********************************************************
  ! 公比を決定 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  r = 1.030898977259116 !FUNCR(n - 1, a, 1.0)
  !r is calculated with Mathematica FindRoot
  ! 等比級数による補間関数を計算 +++++++++++++++++++++++++++++++++++++++
  x(1) = 0.0
  x(n) = 1.0
  Check = .FALSE.
  IF(r .NE. 1.0) THEN
    ! 等比級数を計算 ---------------------------------------------------
    DO i = 2, n - 1
      x(i) = a * (1.0 - r**(i - 1)) / (1.0 - r)
    ENDDO
    ! 単調性を検査 -----------------------------------------------------
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
  ! 公比が1の場合か単調性が保証されていない場合は等間隔の結果を返す ++++
  IF(Check) THEN
    DO i = 2, n - 1
      x(i) = REAL(i - 1) / REAL(n - 1)
    ENDDO
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
! 定義終了 *************************************************************
END SUBROUTINE GeometricInterpolation
!***********************************************************************
!**** 等比級数による一次元補間関数                                  ****
!**** 初項a、公比rの等比数列の和(等比級数)を逆順に並べて            ****
!**** 区間0 <= x <= 1にn個の点を配置する(終点の幅がaとなる)         ****
!***********************************************************************
SUBROUTINE GeometricInterpolationInv(a, n, x, r)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: a
  INTEGER, INTENT(IN)  :: n
  REAL,    INTENT(OUT) :: x(n)
  REAL,    INTENT(OUT) :: r
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i
  LOGICAL :: Check
  ! 処理開始 ***********************************************************
  ! 公比を決定 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  r = 1.00700309717312 !FUNCR(n - 1, a, 1.0)
  ! 等比級数による補間関数を計算 +++++++++++++++++++++++++++++++++++++++
  x(1) = 0.0
  x(n) = 1.0
  Check = .FALSE.
  IF(r .NE. 1.0) THEN
    ! 等比級数を計算 ---------------------------------------------------
    DO i = 2, n - 1
      x(i) = 1.0 - a * (1.0 - r**(n - i)) / (1.0 - r)
    ENDDO
    ! 単調性を検査 -----------------------------------------------------
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
  ! 公比が1の場合か単調性が保証されていない場合は等間隔の結果を返す ++++
  IF(Check) THEN
    DO i = 2, n - 1
      x(i) = REAL(i - 1) / REAL(n - 1)
    ENDDO
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
! 定義終了 *************************************************************
END SUBROUTINE GeometricInterpolationInv
!***********************************************************************
!**** Vinokurの一次元補間関数(1983, J. Comp. Phys., 50, 215-234)    ****
!****   区間0 <= x <= 1にN個の点を配置する                          ****
!****   Dx1 = x(2) - x(1), DxN = x(N) - x(N-1)                      ****
!***********************************************************************
SUBROUTINE VinokurInterpolation(Dx1, DxN, N, x)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, PARAMETER :: JMAX = 100
  REAL(8), PARAMETER :: ZERO = 1.0D-20, RESMIN = 1.0D-12
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: Dx1, DxN
  INTEGER, INTENT(IN)  :: N
  REAL,    INTENT(OUT) :: x(N)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i, j
  REAL(8) :: S0, S1, B, A, xi, DZ, U
  REAL(8) :: f, fdz, ddz
  LOGICAL :: Check
  ! 処理開始 ***********************************************************
  ! Vinokurの一次元補間関数を計算 ++++++++++++++++++++++++++++++++++++++
  S0 = 1.0D+0 / (DBLE(N - 1) * DBLE(Dx1))
  S1 = 1.0D+0 / (DBLE(N - 1) * DBLE(DxN))
  B  = SQRT(S0 * S1)
  A  = SQRT(S0 / S1)
  IF(B .LT. 1.0D+0) THEN
    DZ = B
    DO j = 1, JMAX
      f   = SIN(DZ) / DZ - B
      fdz = (DZ * COS(DZ) - SIN(DZ)) / DZ**2
      IF(ABS(fdz) .LE. ZERO) fdz = SIGN(ZERO, fdz)
      ddz = - f / fdz
      DZ  = DZ + ddZ
      IF(ABS(DZ) .LE. ZERO) DZ = SIGN(ZERO, DZ)
      IF(ABS(ddz) .LT. DZ * RESMIN) EXIT
    ENDDO
    DO i = 1, N
      xi = DBLE(i - 1) / DBLE(N - 1)
      x(i) = REAL( &
      &      TAN(xi * DZ) &
      &    / (A * SIN(DZ) + (1.0D+0 - A * COS(DZ)) * TAN(xi * DZ)) &
      &    )
    ENDDO
  ELSEIF(B .GT. 1.0D+0) THEN
    DZ = B
    DO j = 1, JMAX
      f   = SINH(DZ) / DZ - B
      fdz = (DZ * COSH(DZ) - SINH(DZ)) / DZ**2
      IF(ABS(fdz) .LE. ZERO) fdz = SIGN(ZERO, fdz)
      ddz = - f / fdz
      DZ  = DZ + ddZ
      IF(ABS(DZ) .LE. ZERO) DZ = SIGN(ZERO, DZ)
      IF(ABS(ddz) .LT. DZ * RESMIN) EXIT
    ENDDO
    DO i = 1, N
      xi = DBLE(i - 1) / DBLE(N - 1)
      x(i) = REAL( &
      &      TANH(xi * DZ) &
      &    / (A * SINH(DZ) + (1.0D+0 - A * COSH(DZ)) * TANH(xi * DZ)) &
      &    )
    ENDDO
  ELSE
    DO i = 1, N
      xi = DBLE(i - 1) / DBLE(N - 1)
      U  = xi * ( 1.0D+0 &
      &  + 2.0D+0 * (B - 1.0D+0) * (xi - 0.5D+0) * (1.0D+0 - xi) &
      &  )
      x(i) = REAL(U / (A + (1.0D+0 - A) * U))
    ENDDO
  ENDIF
  ! 始点と終点を補正 +++++++++++++++++++++++++++++++++++++++++++++++++++
  x(1) = 0.0
  x(n) = 1.0
  ! 単調性を検査 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Check = .FALSE.
  DO i = 2, n - 1
    IF(x(i-1) .GE. x(i) .OR. x(i) .GE. x(i+1)) THEN
      WRITE(*, '(A)') 'VinokurInterpolation -> No Monotone Error'
      Check = .TRUE.
      EXIT
    ENDIF
  ENDDO
  ! 公比が1の場合か単調性が保証されていない場合は等間隔の結果を返す ++++
  IF(Check) THEN
    DO i = 2, n - 1
      x(i) = REAL(i - 1) / REAL(n - 1)
    ENDDO
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
! 定義終了 *************************************************************
END SUBROUTINE VinokurInterpolation
!***********************************************************************
!**** 楕円型偏微分方程式(Elliptic Partial Differential Equation)    ****
!**** による格子生成法(一次元ラプラス方程式)                        ****
!***********************************************************************
SUBROUTINE GridEPDE1D( &
&            OMG, &
&            IS, IE, &
&            X, DMAX &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: OMG                             ! 緩和係数
  INTEGER, INTENT(IN)    :: IS, IE                          ! 格子数
  ! 格子点座標 ---------------------------------------------------------
  REAL,    INTENT(INOUT) :: X(IS: IE)                       ! x
  ! 残差 ---------------------------------------------------------------
  REAL,    INTENT(OUT)   :: DMAX                            ! 最大値
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER           :: I
  REAL, ALLOCATABLE :: DX(:)
  ! 処理開始 ***********************************************************
  ALLOCATE(DX(IS: IE)); DX(:) = 0.0
  DO I = IS + 1, IE - 1
    DX(I) = (X(I+1) + X(I-1)) / 2.0 - X(I)
  ENDDO
  DO I = IS + 1, IE - 1
    X(I) = X(I) + OMG * DX(I)
  ENDDO
  DMAX = MAXVAL(SQRT(DX(:)**2))
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE GridEPDE1D
!***********************************************************************
!**** 楕円型偏微分方程式(Elliptic Partial Differential Equation)    ****
!**** による格子生成法(二次元ラプラス方程式)                        ****
!***********************************************************************
SUBROUTINE GridEPDE2D( &
&            OMG, &
&            IS, IE, JS, JE, &
&            X, Y, DMAX &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: OMG                             ! 緩和係数
  INTEGER, INTENT(IN)    :: IS, IE, JS, JE                  ! 格子数
  ! 格子点座標 ---------------------------------------------------------
  REAL,    INTENT(INOUT) :: X(IS: IE, JS: JE), &            ! x
  &                         Y(IS: IE, JS: JE)               ! y
  ! 残差 ---------------------------------------------------------------
  REAL,    INTENT(OUT)   :: DMAX                            ! 最大値
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER           :: I, J
  REAL              :: XXI, XET, YXI, YET
  REAL              :: AJA
  REAL              :: AJA11, AJA12, AJA21, AJA22
  REAL              :: ALP11, ALP12, ALP21, ALP22
  REAL              :: DXC, DXM, DYC, DYM
  REAL, ALLOCATABLE :: DX(:, :), DY(:, :)
  ! 処理開始 ***********************************************************
  ALLOCATE(DX(IS: IE, JS: JE)); DX(:, :) = 0.0
  ALLOCATE(DY(IS: IE, JS: JE)); DY(:, :) = 0.0
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    ! x_xi (x = x, y; xi = xi, eta) ------------------------------------
    XXI =-0.5 * (X(I-1,J) - X(I+1,J))
    YXI =-0.5 * (Y(I-1,J) - Y(I+1,J))
    XET =-0.5 * (X(I,J-1) - X(I,J+1))
    YET =-0.5 * (Y(I,J-1) - Y(I,J+1))
    ! J (Jacobian) -----------------------------------------------------
    AJA = XXI * YET - XET * YXI
    ! J_ij (i = 1, 2; j = 1, 2) ----------------------------------------
    AJA11 = YET
    AJA12 =-YXI
    AJA21 =-XET
    AJA22 = XXI
    ! alpha_ij (i = 1, 2; j = 1, 2) ------------------------------------
    ALP11 = AJA11 * AJA11 + AJA21 * AJA21
    ALP12 = AJA11 * AJA12 + AJA21 * AJA22
    ALP21 = AJA12 * AJA11 + AJA22 * AJA21
    ALP22 = AJA12 * AJA12 + AJA22 * AJA22
    ! x, y -------------------------------------------------------------
    DXC = ALP11 * (X(I+1,J) + X(I-1,J)) &
    &   + ALP22 * (X(I,J+1) + X(I,J-1)) &
    &   + 0.5 * ALP12 * ( &
    &     X(I-1,J-1) - X(I+1,J-1) - X(I-1,J+1) + X(I+1,J+1) &
    &   )
    DXM = 2.0 * (ALP11 + ALP22)
    DYC = ALP11 * (Y(I+1,J) + Y(I-1,J)) &
    &   + ALP22 * (Y(I,J+1) + Y(I,J-1)) &
    &   + 0.5 * ALP12 * ( &
    &     Y(I-1,J-1) - Y(I+1,J-1) - Y(I-1,J+1) + Y(I+1,J+1) &
    &   )
    DYM = 2.0 * (ALP11 + ALP22)
    IF(DXM .NE. 0.0 .AND. DYM .NE. 0.0) THEN
      DX(I,J) = DXC / DXM - X(I,J)
      DY(I,J) = DYC / DYM - Y(I,J)
    ENDIF
  ENDDO
  ENDDO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    X(I,J) = X(I,J) + OMG * DX(I,J)
    Y(I,J) = Y(I,J) + OMG * DY(I,J)
  ENDDO
  ENDDO
  DMAX = MAXVAL(SQRT(DX(:, :)**2 + DY(:, :)**2))
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE GridEPDE2D
!***********************************************************************
!**** 楕円型偏微分方程式(Elliptic Partial Differential Equation)    ****
!**** による格子生成法(三次元ラプラス方程式)                        ****
!***********************************************************************
SUBROUTINE GridEPDE3D( &
&            OMG, &
&            IS, IE, JS, JE, KS, KE, &
&            X, Y, Z, DMAX &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: OMG                             ! 緩和係数
  INTEGER, INTENT(IN)    :: IS, IE, JS, JE, KS, KE          ! 格子数
  ! 格子点座標 ---------------------------------------------------------
  REAL,    INTENT(INOUT) :: X(IS: IE, JS: JE, KS: KE), &    ! x
  &                         Y(IS: IE, JS: JE, KS: KE), &    ! y
  &                         Z(IS: IE, JS: JE, KS: KE)       ! z
  ! 残差 ---------------------------------------------------------------
  REAL,    INTENT(OUT)   :: DMAX                            ! 最大値
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER           :: I, J, K
  REAL              :: XXI, XET, XZE, YXI, YET, YZE, ZXI, ZET, ZZE
  REAL              :: AJA
  REAL              :: AJA11, AJA12, AJA13, &
  &                    AJA21, AJA22, AJA23, &
  &                    AJA31, AJA32, AJA33
  REAL              :: ALP11, ALP12, ALP13, &
  &                    ALP21, ALP22, ALP23, &
  &                    ALP31, ALP32, ALP33
  REAL              :: DXC, DXM, DYC, DYM, DZC, DZM
  REAL, ALLOCATABLE :: DX(:, :, :), DY(:, :, :), DZ(:, :, :)
  ! 処理開始 ***********************************************************
  ALLOCATE(DX(IS: IE, JS: JE, KS: KE)); DX(:, :, :) = 0.0
  ALLOCATE(DY(IS: IE, JS: JE, KS: KE)); DY(:, :, :) = 0.0
  ALLOCATE(DZ(IS: IE, JS: JE, KS: KE)); DZ(:, :, :) = 0.0
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    ! x_xi (x = x, y, z; xi = xi, eta, zeta) ---------------------------
    XXI =-0.5 * (X(I-1,J,K) - X(I+1,J,K))
    YXI =-0.5 * (Y(I-1,J,K) - Y(I+1,J,K))
    ZXI =-0.5 * (Z(I-1,J,K) - Z(I+1,J,K))
    XET =-0.5 * (X(I,J-1,K) - X(I,J+1,K))
    YET =-0.5 * (Y(I,J-1,K) - Y(I,J+1,K))
    ZET =-0.5 * (Z(I,J-1,K) - Z(I,J+1,K))
    XZE =-0.5 * (X(I,J,K-1) - X(I,J,K+1))
    YZE =-0.5 * (Y(I,J,K-1) - Y(I,J,K+1))
    ZZE =-0.5 * (Z(I,J,K-1) - Z(I,J,K+1))
    ! J (Jacobian) -----------------------------------------------------
    AJA = XXI * (YET * ZZE - YZE * ZET) &
    &   + XET * (YZE * ZXI - YXI * ZZE) &
    &   + XZE * (YXI * ZET - YET * ZXI)
    ! J_ij (i = 1, 3; j = 1, 3) ----------------------------------------
    AJA11 = YET * ZZE - YZE * ZET
    AJA12 = YZE * ZXI - YXI * ZZE
    AJA13 = YXI * ZET - YET * ZXI
    AJA21 = XZE * ZET - XET * ZZE
    AJA22 = XXI * ZZE - XZE * ZXI
    AJA23 = XET * ZXI - XXI * ZET
    AJA31 = XET * YZE - XZE * YET
    AJA32 = XZE * YXI - XXI * YZE
    AJA33 = XXI * YET - XET * YXI
    ! alpha_ij (i = 1, 3; j = 1, 3) ------------------------------------
    ALP11 = AJA11 * AJA11 + AJA21 * AJA21 + AJA31 * AJA31
    ALP12 = AJA11 * AJA12 + AJA21 * AJA22 + AJA31 * AJA32
    ALP13 = AJA11 * AJA13 + AJA21 * AJA23 + AJA31 * AJA33
    ALP21 = AJA12 * AJA11 + AJA22 * AJA21 + AJA32 * AJA31
    ALP22 = AJA12 * AJA12 + AJA22 * AJA22 + AJA32 * AJA32
    ALP23 = AJA12 * AJA13 + AJA22 * AJA23 + AJA32 * AJA33
    ALP31 = AJA13 * AJA11 + AJA23 * AJA21 + AJA33 * AJA31
    ALP32 = AJA13 * AJA12 + AJA23 * AJA22 + AJA33 * AJA32
    ALP33 = AJA13 * AJA13 + AJA23 * AJA23 + AJA33 * AJA33
    ! x, y, z ----------------------------------------------------------
    DXC = ALP11 * (X(I+1,J,K) + X(I-1,J,K)) &
    &   + ALP22 * (X(I,J+1,K) + X(I,J-1,K)) &
    &   + ALP33 * (X(I,J,K+1) + X(I,J,K-1)) &
    &   + 0.5 * ALP12 * ( &
    &     X(I-1,J-1,K) - X(I+1,J-1,K) - X(I-1,J+1,K) + X(I+1,J+1,K) &
    &   ) &
    &   + 0.5 * ALP23 * ( &
    &     X(I,J-1,K-1) - X(I,J+1,K-1) - X(I,J-1,K+1) + X(I,J+1,K+1) &
    &   ) &
    &   + 0.5 * ALP31 * ( &
    &     X(I-1,J,K-1) - X(I+1,J,K-1) - X(I-1,J,K+1) + X(I+1,J,K+1) &
    &   )
    DXM = 2.0 * (ALP11 + ALP22 + ALP33)
    DYC = ALP11 * (Y(I+1,J,K) + Y(I-1,J,K)) &
    &   + ALP22 * (Y(I,J+1,K) + Y(I,J-1,K)) &
    &   + ALP33 * (Y(I,J,K+1) + Y(I,J,K-1)) &
    &   + 0.5 * ALP12 * ( &
    &     Y(I-1,J-1,K) - Y(I+1,J-1,K) - Y(I-1,J+1,K) + Y(I+1,J+1,K) &
    &   ) &
    &   + 0.5 * ALP23 * ( &
    &     Y(I,J-1,K-1) - Y(I,J+1,K-1) - Y(I,J-1,K+1) + Y(I,J+1,K+1) &
    &   ) &
    &   + 0.5 * ALP31 * ( &
    &     Y(I-1,J,K-1) - Y(I+1,J,K-1) - Y(I-1,J,K+1) + Y(I+1,J,K+1) &
    &   )
    DYM = 2.0 * (ALP11 + ALP22 + ALP33)
    DZC = ALP11 * (Z(I+1,J,K) + Z(I-1,J,K)) &
    &   + ALP22 * (Z(I,J+1,K) + Z(I,J-1,K)) &
    &   + ALP33 * (Z(I,J,K+1) + Z(I,J,K-1)) &
    &   + 0.5 * ALP12 * ( &
    &     Z(I-1,J-1,K) - Z(I+1,J-1,K) - Z(I-1,J+1,K) + Z(I+1,J+1,K) &
    &   ) &
    &   + 0.5 * ALP23 * ( &
    &     Z(I,J-1,K-1) - Z(I,J+1,K-1) - Z(I,J-1,K+1) + Z(I,J+1,K+1) &
    &   ) &
    &   + 0.5 * ALP31 * ( &
    &     Z(I-1,J,K-1) - Z(I+1,J,K-1) - Z(I-1,J,K+1) + Z(I+1,J,K+1) &
    &   )
    DZM = 2.0 * (ALP11 + ALP22 + ALP33)
    IF(DXM .NE. 0.0 .AND. DYM .NE. 0.0 .AND. DZM .NE. 0.0) THEN
      DX(I,J,K) = DXC / DXM - X(I,J,K)
      DY(I,J,K) = DYC / DYM - Y(I,J,K)
      DZ(I,J,K) = DZC / DZM - Z(I,J,K)
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    X(I,J,K) = X(I,J,K) + OMG * DX(I,J,K)
    Y(I,J,K) = Y(I,J,K) + OMG * DY(I,J,K)
    Z(I,J,K) = Z(I,J,K) + OMG * DZ(I,J,K)
  ENDDO
  ENDDO
  ENDDO
  DMAX = MAXVAL(SQRT(DX(:, :, :)**2 + DY(:, :, :)**2 + DZ(:, :, :)**2))
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE GridEPDE3D
!***********************************************************************
!**** 一次元楕円-双曲型格子生成法(Yamamoto Kazuomi, 1992)           ****
!**** (1D Elliptic-Hyperbolic Grid Generation Method)               ****
!***********************************************************************
SUBROUTINE GridEHPDE1D( &
&            OMG, &
&            IS, IE, &
&            Cs, &
&            Cv1m, Cv1p, &
&            DX1m, DX1p, &
&            X, DMAX &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: OMG                             ! 緩和係数
  INTEGER, INTENT(IN)    :: IS, IE                          ! 格子数
  ! 楕円型重み関数(0以上) ----------------------------------------------
  REAL,    INTENT(IN)    :: Cs(IS: IE)                      ! 1-Cvが基本
  ! 双曲型重み関数(0以上) ----------------------------------------------
  REAL,    INTENT(IN)    :: Cv1m(IS: IE), &                 ! I-方向
  &                         Cv1p(IS: IE)                    ! I+方向
  ! 境界の傾きと幅(partial x_i / partial xi_i) -------------------------
  REAL,    INTENT(IN)    :: DX1m(IS: IE), &                 ! I-方向(x)
  &                         DX1p(IS: IE)                    ! I+方向(x)
  ! 格子点座標 ---------------------------------------------------------
  REAL,    INTENT(INOUT) :: X(IS: IE)                       ! x
  ! 残差 ---------------------------------------------------------------
  REAL,    INTENT(OUT)   :: DMAX                            ! 最大値
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  LOGICAL           :: FLAG
  INTEGER           :: I
  REAL              :: DXC, DXM
  REAL, ALLOCATABLE :: DX(:)
  ! 処理開始 ***********************************************************
  ! 例外処理(重み関数が負の場合は終了) +++++++++++++++++++++++++++++++++
  FLAG = .FALSE.
  DO I = IS + 1, IE - 1
    FLAG = FLAG .OR. Cs(I) .LT. 0.0 .OR. &
    &      Cv1m(I) .LT. 0.0 .OR. Cv1p(I) .LT. 0.0
  ENDDO
  IF(FLAG) THEN
    WRITE(*, '(A)') 'Error : Weight function is negative...'
    WRITE(*, '(A)') 'GridEHPDE1D is abnormally Terminated !'
    RETURN
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(DX(IS: IE)); DX(:) = 0.0
  ! ソルバー部 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO I = IS + 1, IE - 1
    DXC = Cs(I) * (X(I+1) + X(I-1)) &
    &   + Cv1p(I) * (X(I+1) - DX1p(I)) &
    &   + Cv1m(I) * (X(I-1) + DX1m(I))
    DXM = 2.0 * Cs(I) + Cv1p(I) + Cv1m(I)
    IF(DXM .NE. 0.0) THEN
      DX(I) = DXC / DXM - X(I)
    ENDIF
  ENDDO
  DO I = IS + 1, IE - 1
    X(I) = X(I) + OMG * DX(I)
  ENDDO
  DMAX = MAXVAL(SQRT(DX(:)**2))
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE GridEHPDE1D
!***********************************************************************
!**** 二次元楕円-双曲型格子生成法(Yamamoto Kazuomi, 1992)           ****
!**** (2D Elliptic-Hyperbolic Grid Generation Method)               ****
!***********************************************************************
SUBROUTINE GridEHPDE2D( &
&            OMG, &
&            IS, IE, JS, JE, &
&            Cs, &
&            Cv1m, Cv1p, Cv2m, Cv2p, &
&            DX1m, DY1m, DX1p, DY1p, DX2m, DY2m, DX2p, DY2p, &
&            X, Y, DMAX &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: OMG                             ! 緩和係数
  INTEGER, INTENT(IN)    :: IS, IE, JS, JE                  ! 格子数
  ! 楕円型重み関数(0以上) ----------------------------------------------
  REAL,    INTENT(IN)    :: Cs(IS: IE, JS: JE)              ! 1-Cvが基本
  ! 双曲型重み関数(0以上) ----------------------------------------------
  REAL,    INTENT(IN)    :: Cv1m(IS: IE, JS: JE), &         ! I-方向
  &                         Cv1p(IS: IE, JS: JE), &         ! I+方向
  &                         Cv2m(IS: IE, JS: JE), &         ! J-方向
  &                         Cv2p(IS: IE, JS: JE)            ! J+方向
  ! 境界の傾きと幅(partial x_i / partial xi_i) -------------------------
  REAL,    INTENT(IN)    :: DX1m(IS: IE, JS: JE), &         ! I-方向(x)
  &                         DY1m(IS: IE, JS: JE), &         ! I-方向(y)
  &                         DX1p(IS: IE, JS: JE), &         ! I+方向(x)
  &                         DY1p(IS: IE, JS: JE), &         ! I+方向(y)
  &                         DX2m(IS: IE, JS: JE), &         ! J-方向(x)
  &                         DY2m(IS: IE, JS: JE), &         ! J-方向(y)
  &                         DX2p(IS: IE, JS: JE), &         ! J+方向(x)
  &                         DY2p(IS: IE, JS: JE)            ! J+方向(y)
  ! 格子点座標 ---------------------------------------------------------
  REAL,    INTENT(INOUT) :: X(IS: IE, JS: JE), &            ! x
  &                         Y(IS: IE, JS: JE)               ! y
  ! 残差 ---------------------------------------------------------------
  REAL,    INTENT(OUT)   :: DMAX                            ! 最大値
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  LOGICAL           :: FLAG
  INTEGER           :: I, J
  REAL              :: XXI, XET, YXI, YET
  REAL              :: AJA
  REAL              :: AJA11, AJA12, AJA21, AJA22
  REAL              :: ALP11, ALP12, ALP21, ALP22
  REAL              :: DXC, DXM, DYC, DYM
  REAL, ALLOCATABLE :: DX(:, :), DY(:, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(重み関数が負の場合は終了) +++++++++++++++++++++++++++++++++
  FLAG = .FALSE.
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    FLAG = FLAG .OR. Cs(I,J) .LT. 0.0 .OR. &
    &      Cv1m(I,J) .LT. 0.0 .OR. Cv1p(I,J) .LT. 0.0 .OR. &
    &      Cv2m(I,J) .LT. 0.0 .OR. Cv2p(I,J) .LT. 0.0
  ENDDO
  ENDDO
  IF(FLAG) THEN
    WRITE(*, '(A)') 'Error : Weight function is negative...'
    WRITE(*, '(A)') 'GridEHPDE2D is abnormally Terminated !'
    RETURN
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(DX(IS: IE, JS: JE)); DX(:, :) = 0.0
  ALLOCATE(DY(IS: IE, JS: JE)); DY(:, :) = 0.0
  ! ソルバー部 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    ! x_xi (x = x, y; xi = xi, eta) ------------------------------------
    XXI =-0.5 * (X(I-1,J) - X(I+1,J))
    YXI =-0.5 * (Y(I-1,J) - Y(I+1,J))
    XET =-0.5 * (X(I,J-1) - X(I,J+1))
    YET =-0.5 * (Y(I,J-1) - Y(I,J+1))
    ! J (Jacobian) -----------------------------------------------------
    AJA = XXI * YET - XET * YXI
    ! J_ij (i = 1, 2; j = 1, 2) ----------------------------------------
    AJA11 = YET
    AJA12 =-YXI
    AJA21 =-XET
    AJA22 = XXI
    ! alpha_ij (i = 1, 2; j = 1, 2) ------------------------------------
    ALP11 = AJA11 * AJA11 + AJA21 * AJA21
    ALP12 = AJA11 * AJA12 + AJA21 * AJA22
    ALP21 = AJA12 * AJA11 + AJA22 * AJA21
    ALP22 = AJA12 * AJA12 + AJA22 * AJA22
    ! x, y -------------------------------------------------------------
    DXC = Cs(I,J) * ( &
    &       ALP11 * (X(I+1,J) + X(I-1,J)) &
    &     + ALP22 * (X(I,J+1) + X(I,J-1)) &
    &     + 0.5 * ALP12 * ( &
    &       X(I-1,J-1) - X(I+1,J-1) - X(I-1,J+1) + X(I+1,J+1) &
    &     ) &
    &   ) &
    &   + ALP11 * ( Cv1p(I,J) * (X(I+1,J) - DX1p(I,J)) &
    &             + Cv1m(I,J) * (X(I-1,J) + DX1m(I,J)) ) &
    &   + ALP22 * ( Cv2p(I,J) * (X(I,J+1) - DX2p(I,J)) &
    &             + Cv2m(I,J) * (X(I,J-1) + DX2m(I,J)) )
    DXM = 2.0 * Cs(I,J) * (ALP11 + ALP22) &
    &   + ALP11 * (Cv1p(I,J) + Cv1m(I,J)) &
    &   + ALP22 * (Cv2p(I,J) + Cv2m(I,J))
    DYC = Cs(I,J) * ( &
    &       ALP11 * (Y(I+1,J) + Y(I-1,J)) &
    &     + ALP22 * (Y(I,J+1) + Y(I,J-1)) &
    &     + 0.5 * ALP12 * ( &
    &       Y(I-1,J-1) - Y(I+1,J-1) - Y(I-1,J+1) + Y(I+1,J+1) &
    &     ) &
    &   ) &
    &   + ALP11 * ( Cv1p(I,J) * (Y(I+1,J) - DY1p(I,J)) &
    &             + Cv1m(I,J) * (Y(I-1,J) + DY1m(I,J)) ) &
    &   + ALP22 * ( Cv2p(I,J) * (Y(I,J+1) - DY2p(I,J)) &
    &             + Cv2m(I,J) * (Y(I,J-1) + DY2m(I,J)) )
    DYM = 2.0 * Cs(I,J) * (ALP11 + ALP22) &
    &   + ALP11 * (Cv1p(I,J) + Cv1m(I,J)) &
    &   + ALP22 * (Cv2p(I,J) + Cv2m(I,J))
    IF(DXM .NE. 0.0 .AND. DYM .NE. 0.0) THEN
      DX(I,J) = DXC / DXM - X(I,J)
      DY(I,J) = DYC / DYM - Y(I,J)
    ENDIF
  ENDDO
  ENDDO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    X(I,J) = X(I,J) + OMG * DX(I,J)
    Y(I,J) = Y(I,J) + OMG * DY(I,J)
  ENDDO
  ENDDO
  DMAX = MAXVAL(SQRT(DX(:, :)**2 + DY(:, :)**2))
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE GridEHPDE2D
!***********************************************************************
!**** 三次元楕円-双曲型格子生成法(Yamamoto Kazuomi, 1992)           ****
!**** (3D Elliptic-Hyperbolic Grid Generation Method)               ****
!***********************************************************************
SUBROUTINE GridEHPDE3D( &
&            OMG, &
&            IS, IE, JS, JE, KS, KE, &
&            Cs, &
&            Cv1m, Cv1p, Cv2m, Cv2p, Cv3m, Cv3p, &
&            DX1m, DY1m, DZ1m, DX1p, DY1p, DZ1p, &
&            DX2m, DY2m, DZ2m, DX2p, DY2p, DZ2p, &
&            DX3m, DY3m, DZ3m, DX3p, DY3p, DZ3p, &
&            X, Y, Z, DMAX &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: OMG                             ! 緩和係数
  INTEGER, INTENT(IN)    :: IS, IE, JS, JE, KS, KE          ! 格子数
  ! 楕円型重み関数(0以上) ----------------------------------------------
  REAL,    INTENT(IN)    :: Cs(IS: IE, JS: JE, KS: KE)      ! 1-Cvが基本
  ! 双曲型重み関数(0以上) ----------------------------------------------
  REAL,    INTENT(IN)    :: Cv1m(IS: IE, JS: JE, KS: KE), & ! I-方向
  &                         Cv1p(IS: IE, JS: JE, KS: KE), & ! I+方向
  &                         Cv2m(IS: IE, JS: JE, KS: KE), & ! J-方向
  &                         Cv2p(IS: IE, JS: JE, KS: KE), & ! J+方向
  &                         Cv3m(IS: IE, JS: JE, KS: KE), & ! K-方向
  &                         Cv3p(IS: IE, JS: JE, KS: KE)    ! K+方向
  ! 境界の傾きと幅(partial x_i / partial xi_i) -------------------------
  REAL,    INTENT(IN)    :: DX1m(IS: IE, JS: JE, KS: KE), & ! I方向-(x)
  &                         DY1m(IS: IE, JS: JE, KS: KE), & ! I方向-(y)
  &                         DZ1m(IS: IE, JS: JE, KS: KE), & ! I方向-(z)
  &                         DX1p(IS: IE, JS: JE, KS: KE), & ! I方向+(x)
  &                         DY1p(IS: IE, JS: JE, KS: KE), & ! I方向+(y)
  &                         DZ1p(IS: IE, JS: JE, KS: KE), & ! I方向+(z)
  &                         DX2m(IS: IE, JS: JE, KS: KE), & ! J方向-(x)
  &                         DY2m(IS: IE, JS: JE, KS: KE), & ! J方向-(y)
  &                         DZ2m(IS: IE, JS: JE, KS: KE), & ! J方向-(z)
  &                         DX2p(IS: IE, JS: JE, KS: KE), & ! J方向+(x)
  &                         DY2p(IS: IE, JS: JE, KS: KE), & ! J方向+(y)
  &                         DZ2p(IS: IE, JS: JE, KS: KE), & ! J方向+(z)
  &                         DX3m(IS: IE, JS: JE, KS: KE), & ! K方向-(x)
  &                         DY3m(IS: IE, JS: JE, KS: KE), & ! K方向-(y)
  &                         DZ3m(IS: IE, JS: JE, KS: KE), & ! K方向-(z)
  &                         DX3p(IS: IE, JS: JE, KS: KE), & ! K方向+(x)
  &                         DY3p(IS: IE, JS: JE, KS: KE), & ! K方向+(y)
  &                         DZ3p(IS: IE, JS: JE, KS: KE)    ! K方向+(z)
  ! 格子点座標 ---------------------------------------------------------
  REAL,    INTENT(INOUT) :: X(IS: IE, JS: JE, KS: KE), &    ! x
  &                         Y(IS: IE, JS: JE, KS: KE), &    ! y
  &                         Z(IS: IE, JS: JE, KS: KE)       ! z
  ! 残差 ---------------------------------------------------------------
  REAL,    INTENT(OUT)   :: DMAX                            ! 最大値
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  LOGICAL           :: FLAG
  INTEGER           :: I, J, K
  REAL              :: XXI, XET, XZE, YXI, YET, YZE, ZXI, ZET, ZZE
  REAL              :: AJA
  REAL              :: AJA11, AJA12, AJA13, &
  &                    AJA21, AJA22, AJA23, &
  &                    AJA31, AJA32, AJA33
  REAL              :: ALP11, ALP12, ALP13, &
  &                    ALP21, ALP22, ALP23, &
  &                    ALP31, ALP32, ALP33
  REAL              :: DXC, DXM, DYC, DYM, DZC, DZM
  REAL, ALLOCATABLE :: DX(:, :, :), DY(:, :, :), DZ(:, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(重み関数が負の場合は終了) +++++++++++++++++++++++++++++++++
  FLAG = .FALSE.
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    FLAG = FLAG .OR. Cs(I,J,K) .LT. 0.0 .OR. &
    &      Cv1m(I,J,K) .LT. 0.0 .OR. Cv1p(I,J,K) .LT. 0.0 .OR. &
    &      Cv2m(I,J,K) .LT. 0.0 .OR. Cv2p(I,J,K) .LT. 0.0 .OR. &
    &      Cv3m(I,J,K) .LT. 0.0 .OR. Cv3p(I,J,K) .LT. 0.0
  ENDDO
  ENDDO
  ENDDO
  IF(FLAG) THEN
    WRITE(*, '(A)') 'Error : Weight function is negative...'
    WRITE(*, '(A)') 'GridEHPDE3D is abnormally Terminated !'
    RETURN
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(DX(IS: IE, JS: JE, KS: KE)); DX(:, :, :) = 0.0
  ALLOCATE(DY(IS: IE, JS: JE, KS: KE)); DY(:, :, :) = 0.0
  ALLOCATE(DZ(IS: IE, JS: JE, KS: KE)); DZ(:, :, :) = 0.0
  ! ソルバー部 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    ! x_xi (x = x, y, z; xi = xi, eta, zeta) ---------------------------
    XXI =-0.5 * (X(I-1,J,K) - X(I+1,J,K))
    YXI =-0.5 * (Y(I-1,J,K) - Y(I+1,J,K))
    ZXI =-0.5 * (Z(I-1,J,K) - Z(I+1,J,K))
    XET =-0.5 * (X(I,J-1,K) - X(I,J+1,K))
    YET =-0.5 * (Y(I,J-1,K) - Y(I,J+1,K))
    ZET =-0.5 * (Z(I,J-1,K) - Z(I,J+1,K))
    XZE =-0.5 * (X(I,J,K-1) - X(I,J,K+1))
    YZE =-0.5 * (Y(I,J,K-1) - Y(I,J,K+1))
    ZZE =-0.5 * (Z(I,J,K-1) - Z(I,J,K+1))
    ! J (Jacobian) -----------------------------------------------------
    AJA = XXI * (YET * ZZE - YZE * ZET) &
    &   + XET * (YZE * ZXI - YXI * ZZE) &
    &   + XZE * (YXI * ZET - YET * ZXI)
    ! J_ij (i = 1, 3; j = 1, 3) ----------------------------------------
    AJA11 = YET * ZZE - YZE * ZET
    AJA12 = YZE * ZXI - YXI * ZZE
    AJA13 = YXI * ZET - YET * ZXI
    AJA21 = XZE * ZET - XET * ZZE
    AJA22 = XXI * ZZE - XZE * ZXI
    AJA23 = XET * ZXI - XXI * ZET
    AJA31 = XET * YZE - XZE * YET
    AJA32 = XZE * YXI - XXI * YZE
    AJA33 = XXI * YET - XET * YXI
    ! alpha_ij (i = 1, 3; j = 1, 3) ------------------------------------
    ALP11 = AJA11 * AJA11 + AJA21 * AJA21 + AJA31 * AJA31
    ALP12 = AJA11 * AJA12 + AJA21 * AJA22 + AJA31 * AJA32
    ALP13 = AJA11 * AJA13 + AJA21 * AJA23 + AJA31 * AJA33
    ALP21 = AJA12 * AJA11 + AJA22 * AJA21 + AJA32 * AJA31
    ALP22 = AJA12 * AJA12 + AJA22 * AJA22 + AJA32 * AJA32
    ALP23 = AJA12 * AJA13 + AJA22 * AJA23 + AJA32 * AJA33
    ALP31 = AJA13 * AJA11 + AJA23 * AJA21 + AJA33 * AJA31
    ALP32 = AJA13 * AJA12 + AJA23 * AJA22 + AJA33 * AJA32
    ALP33 = AJA13 * AJA13 + AJA23 * AJA23 + AJA33 * AJA33
    ! x, y, z ----------------------------------------------------------
    DXC = Cs(I,J,K) * ( &
    &       ALP11 * (X(I+1,J,K) + X(I-1,J,K)) &
    &     + ALP22 * (X(I,J+1,K) + X(I,J-1,K)) &
    &     + ALP33 * (X(I,J,K+1) + X(I,J,K-1)) &
    &     + 0.5 * ALP12 * ( &
    &       X(I-1,J-1,K) - X(I+1,J-1,K) - X(I-1,J+1,K) + X(I+1,J+1,K) &
    &     ) &
    &     + 0.5 * ALP23 * ( &
    &       X(I,J-1,K-1) - X(I,J+1,K-1) - X(I,J-1,K+1) + X(I,J+1,K+1) &
    &     ) &
    &     + 0.5 * ALP31 * ( &
    &       X(I-1,J,K-1) - X(I+1,J,K-1) - X(I-1,J,K+1) + X(I+1,J,K+1) &
    &     ) &
    &   ) &
    &   + ALP11 * ( Cv1p(I,J,K) * (X(I+1,J,K) - DX1p(I,J,K)) &
    &             + Cv1m(I,J,K) * (X(I-1,J,K) + DX1m(I,J,K)) ) &
    &   + ALP22 * ( Cv2p(I,J,K) * (X(I,J+1,K) - DX2p(I,J,K)) &
    &             + Cv2m(I,J,K) * (X(I,J-1,K) + DX2m(I,J,K)) ) &
    &   + ALP33 * ( Cv3p(I,J,K) * (X(I,J,K+1) - DX3p(I,J,K)) &
    &             + Cv3m(I,J,K) * (X(I,J,K-1) + DX3m(I,J,K)) )
    DXM = 2.0 * Cs(I,J,K) * (ALP11 + ALP22 + ALP33) &
    &   + ALP11 * (Cv1p(I,J,K) + Cv1m(I,J,K)) &
    &   + ALP22 * (Cv2p(I,J,K) + Cv2m(I,J,K)) &
    &   + ALP33 * (Cv3p(I,J,K) + Cv3m(I,J,K))
    DYC = Cs(I,J,K) * ( &
    &       ALP11 * (Y(I+1,J,K) + Y(I-1,J,K)) &
    &     + ALP22 * (Y(I,J+1,K) + Y(I,J-1,K)) &
    &     + ALP33 * (Y(I,J,K+1) + Y(I,J,K-1)) &
    &     + 0.5 * ALP12 * ( &
    &       Y(I-1,J-1,K) - Y(I+1,J-1,K) - Y(I-1,J+1,K) + Y(I+1,J+1,K) &
    &     ) &
    &     + 0.5 * ALP23 * ( &
    &       Y(I,J-1,K-1) - Y(I,J+1,K-1) - Y(I,J-1,K+1) + Y(I,J+1,K+1) &
    &     ) &
    &     + 0.5 * ALP31 * ( &
    &       Y(I-1,J,K-1) - Y(I+1,J,K-1) - Y(I-1,J,K+1) + Y(I+1,J,K+1) &
    &     ) &
    &   ) &
    &   + ALP11 * ( Cv1p(I,J,K) * (Y(I+1,J,K) - DY1p(I,J,K)) &
    &             + Cv1m(I,J,K) * (Y(I-1,J,K) + DY1m(I,J,K)) ) &
    &   + ALP22 * ( Cv2p(I,J,K) * (Y(I,J+1,K) - DY2p(I,J,K)) &
    &             + Cv2m(I,J,K) * (Y(I,J-1,K) + DY2m(I,J,K)) ) &
    &   + ALP33 * ( Cv3p(I,J,K) * (Y(I,J,K+1) - DY3p(I,J,K)) &
    &             + Cv3m(I,J,K) * (Y(I,J,K-1) + DY3m(I,J,K)) )
    DYM = 2.0 * Cs(I,J,K) * (ALP11 + ALP22 + ALP33) &
    &   + ALP11 * (Cv1p(I,J,K) + Cv1m(I,J,K)) &
    &   + ALP22 * (Cv2p(I,J,K) + Cv2m(I,J,K)) &
    &   + ALP33 * (Cv3p(I,J,K) + Cv3m(I,J,K))
    DZC = Cs(I,J,K) * ( &
    &       ALP11 * (Z(I+1,J,K) + Z(I-1,J,K)) &
    &     + ALP22 * (Z(I,J+1,K) + Z(I,J-1,K)) &
    &     + ALP33 * (Z(I,J,K+1) + Z(I,J,K-1)) &
    &     + 0.5 * ALP12 * ( &
    &       Z(I-1,J-1,K) - Z(I+1,J-1,K) - Z(I-1,J+1,K) + Z(I+1,J+1,K) &
    &     ) &
    &     + 0.5 * ALP23 * ( &
    &       Z(I,J-1,K-1) - Z(I,J+1,K-1) - Z(I,J-1,K+1) + Z(I,J+1,K+1) &
    &     ) &
    &     + 0.5 * ALP31 * ( &
    &       Z(I-1,J,K-1) - Z(I+1,J,K-1) - Z(I-1,J,K+1) + Z(I+1,J,K+1) &
    &     ) &
    &   ) &
    &   + ALP11 * ( Cv1p(I,J,K) * (Z(I+1,J,K) - DZ1p(I,J,K)) &
    &             + Cv1m(I,J,K) * (Z(I-1,J,K) + DZ1m(I,J,K)) ) &
    &   + ALP22 * ( Cv2p(I,J,K) * (Z(I,J+1,K) - DZ2p(I,J,K)) &
    &             + Cv2m(I,J,K) * (Z(I,J-1,K) + DZ2m(I,J,K)) ) &
    &   + ALP33 * ( Cv3p(I,J,K) * (Z(I,J,K+1) - DZ3p(I,J,K)) &
    &             + Cv3m(I,J,K) * (Z(I,J,K-1) + DZ3m(I,J,K)) )
    DZM = 2.0 * Cs(I,J,K) * (ALP11 + ALP22 + ALP33) &
    &   + ALP11 * (Cv1p(I,J,K) + Cv1m(I,J,K)) &
    &   + ALP22 * (Cv2p(I,J,K) + Cv2m(I,J,K)) &
    &   + ALP33 * (Cv3p(I,J,K) + Cv3m(I,J,K))
    IF(DXM .NE. 0.0 .AND. DYM .NE. 0.0 .AND. DZM .NE. 0.0) THEN
      DX(I,J,K) = DXC / DXM - X(I,J,K)
      DY(I,J,K) = DYC / DYM - Y(I,J,K)
      DZ(I,J,K) = DZC / DZM - Z(I,J,K)
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    X(I,J,K) = X(I,J,K) + OMG * DX(I,J,K)
    Y(I,J,K) = Y(I,J,K) + OMG * DY(I,J,K)
    Z(I,J,K) = Z(I,J,K) + OMG * DZ(I,J,K)
  ENDDO
  ENDDO
  ENDDO
  DMAX = MAXVAL(SQRT(DX(:, :, :)**2 + DY(:, :, :)**2 + DZ(:, :, :)**2))
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE GridEHPDE3D
!***********************************************************************
!**** 二次元二境界法(エルミート三次補間)                            ****
!****   T1, T2 : 格子の直交性を規定する(値が大きいほど直交性が強い) ****
!****   ETABAR : 格子の集中を規定する関数 ETABAR = ETABAR(ETA)      ****
!****            D ETABAR / D ETA が小さいほど格子が集中する        ****
!****            ただし 0 <= ETA <= 1                               ****
!***********************************************************************
SUBROUTINE TwoBoundaryMethod2D( &
&            IS, IE, JS, JE, T1, T2, ETABAR, X, Y &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)    :: IS, IE, JS, JE
  REAL,    INTENT(IN)    :: T1,T2
  REAL,    INTENT(IN)    :: ETABAR(JS:JE)
  REAL,    INTENT(INOUT) :: X(IS:IE,JS:JE), Y(IS:IE,JS:JE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, PARAMETER :: DIMNUM = 2
  REAL :: S1(DIMNUM, IS:IE), S2(DIMNUM, IS:IE), &
  &       G1(DIMNUM, IS:IE), G2(DIMNUM, IS:IE), MU(1:4)
  REAL :: DXDXI1, DYDXI1, DXDXI2, DYDXI2
  INTEGER :: I,J
  ! 処理開始 ***********************************************************
  ! 二境界のベクトル S_K(XI,ZETA) K=1~DIMNUM を定義 ++++++++++++++++++++
  DO I = IS, IE
    S1(1, I) = X(I, JS)
    S1(2, I) = Y(I, JS)
    S2(1, I) = X(I, JE)
    S2(2, I) = Y(I, JE)
  ENDDO
  ! S_K の法線ベクトル G_K(S_K(XI,ZETA),XI) K=1~DIMNUM を定義 ++++++++++
  DO I = IS, IE
    IF( I .EQ. IS ) THEN
      DXDXI1 = S1(1, I+1) - S1(1, I)
      DYDXI1 = S1(2, I+1) - S1(2, I)
      DXDXI2 = S2(1, I+1) - S2(1, I)
      DYDXI2 = S2(2, I+1) - S2(2, I)
    ENDIF
    IF( I .EQ. IE ) THEN
      DXDXI1 = S1(1, I) - S1(1, I-1)
      DYDXI1 = S1(2, I) - S1(2, I-1)
      DXDXI2 = S2(1, I) - S2(1, I-1)
      DYDXI2 = S2(2, I) - S2(2, I-1)
    ENDIF
    IF( ( I .NE. IS ) .AND. ( I .NE. IE ) ) THEN
      DXDXI1 = 0.5 * ( S1(1, I+1) - S1(1, I-1) )
      DYDXI1 = 0.5 * ( S1(2, I+1) - S1(2, I-1) )
      DXDXI2 = 0.5 * ( S2(1, I+1) - S2(1, I-1) )
      DYDXI2 = 0.5 * ( S2(2, I+1) - S2(2, I-1) )
    ENDIF
    G1(1, I) = -DYDXI1 / SQRT( DXDXI1**2 + DYDXI1**2 )
    G1(2, I) =  DXDXI1 / SQRT( DXDXI1**2 + DYDXI1**2 )
    G2(1, I) = -DYDXI2 / SQRT( DXDXI2**2 + DYDXI2**2 )
    G2(2, I) =  DXDXI2 / SQRT( DXDXI2**2 + DYDXI2**2 )
  ENDDO
  ! 補間関数 MU を定義し、座標ベクトル F を計算 ++++++++++++++++++++++++
  DO I = IS, IE
  DO J = JS+1, JE-1
    ! 補間関数の定義 MU_K(ETABAR)  K=1~4 -------------------------------
    MU(1) = 2.0 * ETABAR(J)**3 - 3.0 * ETABAR(J)**2             + 1.0
    MU(2) =-2.0 * ETABAR(J)**3 + 3.0 * ETABAR(J)**2
    MU(3) =       ETABAR(J)**3 - 2.0 * ETABAR(J)**2 + ETABAR(J)
    MU(4) =       ETABAR(J)**3 -       ETABAR(J)**2
    ! 定義 F(XI,MU,ZETA) -----------------------------------------------
    X(I,J) = MU(1) * S1(1,I)           + MU(2) * S2(1,I) &
    &      + MU(3) * T1      * G1(1,I) + MU(4) * T2      * G2(1,I)
    Y(I,J) = MU(1) * S1(2,I)           + MU(2) * S2(2,I) &
    &      + MU(3) * T1      * G1(2,I) + MU(4) * T2      * G2(2,I)
  ENDDO
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
! 定義終了 *************************************************************
END SUBROUTINE TwoBoundaryMethod2D
!***********************************************************************
!**** Transfinite補間法(Eriksson, 1982)                             ****
!****   ・Fの入力(境界上の値は与えられていなくてはならない)         ****
!****   ・Fの出力(境界上以外の値が補間により計算される)             ****
!****   ・混合関数の制限(以下の関係を満たし、単調に変化すること)    ****
!****     1.0 = ALP1(IS), ALP1(IE) = 0.0                            ****
!****     0.0 = ALP2(IS), ALP2(IE) = 1.0                            ****
!***********************************************************************
SUBROUTINE Transfinite1D( &
&            IS, IE, ALP1, ALP2, F &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)    :: IS, IE
  REAL,    INTENT(IN)    :: ALP1(IS: IE), ALP2(IS: IE)
  REAL,    INTENT(INOUT) :: F(IS: IE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I
  ! 処理開始 ***********************************************************
  DO I = IS + 1, IE - 1
    F(I) = ALP1(I) * F(IS) + ALP2(I) * F(IE)
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
! 定義終了 *************************************************************
END SUBROUTINE Transfinite1D
!***********************************************************************
!**** Transfinite補間法(Eriksson, 1982)                             ****
!****   ・Fの入力(境界上の値は与えられていなくてはならない)         ****
!****   ・Fの出力(境界上以外の値が補間により計算される)             ****
!****   ・混合関数の制限(以下の関係を満たし、単調に変化すること)    ****
!****     1.0 = ALP1(IS), ALP1(IE) = 0.0                            ****
!****     0.0 = ALP2(IS), ALP2(IE) = 1.0                            ****
!****     1.0 = BET1(JS), BET1(JE) = 0.0                            ****
!****     0.0 = BET2(JS), BET2(JE) = 1.0                            ****
!***********************************************************************
SUBROUTINE Transfinite2D( &
&            IS, IE, JS, JE, ALP1, ALP2, BET1, BET2, F &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)    :: IS, IE, JS, JE
  REAL,    INTENT(IN)    :: ALP1(IS: IE), ALP2(IS: IE)
  REAL,    INTENT(IN)    :: BET1(JS: JE), BET2(JS: JE)
  REAL,    INTENT(INOUT) :: F(IS: IE, JS: JE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER           :: I, J
  REAL, ALLOCATABLE :: F1(:, :)
  ! 処理開始 ***********************************************************
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(F1(IS: IE, JS: JE))
  ! 一段目 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO J = JS, JE
  DO I = IS, IE
    F1(I,J) = ALP1(I) * F(IS,J) + ALP2(I) * F(IE,J)
  ENDDO
  ENDDO
  ! 二段目 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    F(I,J) = F1(I,J) + BET1(J) * (F(I,JS) - F1(I,JS)) &
    &                + BET2(J) * (F(I,JE) - F1(I,JE))
  ENDDO
  ENDDO
  ! メモリ解放 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DEALLOCATE(F1)
  ! 処理終了 ***********************************************************
  RETURN
! 定義終了 *************************************************************
END SUBROUTINE Transfinite2D
!***********************************************************************
!**** Transfinite補間法(Eriksson, 1982)                             ****
!****   ・Fの入力(境界上の値は与えられていなくてはならない)         ****
!****   ・Fの出力(境界上以外の値が補間により計算される)             ****
!****   ・混合関数の制限(以下の関係を満たし、単調に変化すること)    ****
!****     1.0 = ALP1(IS), ALP1(IE) = 0.0                            ****
!****     0.0 = ALP2(IS), ALP2(IE) = 1.0                            ****
!****     1.0 = BET1(JS), BET1(JE) = 0.0                            ****
!****     0.0 = BET2(JS), BET2(JE) = 1.0                            ****
!****     1.0 = GAM1(KS), GAM1(KE) = 0.0                            ****
!****     0.0 = GAM2(KS), GAM2(KE) = 1.0                            ****
!***********************************************************************
SUBROUTINE Transfinite3D( &
&            IS, IE, JS, JE, KS, KE, &
&            ALP1, ALP2, BET1, BET2, GAM1, GAM2, F &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)    :: IS, IE, JS, JE, KS, KE
  REAL,    INTENT(IN)    :: ALP1(IS: IE), ALP2(IS: IE)
  REAL,    INTENT(IN)    :: BET1(JS: JE), BET2(JS: JE)
  REAL,    INTENT(IN)    :: GAM1(KS: KE), GAM2(KS: KE)
  REAL,    INTENT(INOUT) :: F(IS: IE, JS: JE, KS: KE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER           :: I, J, K
  REAL, ALLOCATABLE :: F1(:, :, :), F2(:, :, :)
  ! 処理開始 ***********************************************************
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(F1(IS: IE, JS: JE, KS: KE), F2(IS: IE, JS: JE, KS: KE))
  ! 一段目 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
    F1(I,J,K) = ALP1(I) * F(IS,J,K) + ALP2(I) * F(IE,J,K)
  ENDDO
  ENDDO
  ENDDO
  ! 二段目 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
    F2(I,J,K) = F1(I,J,K) + BET1(J) * (F(I,JS,K) - F1(I,JS,K)) &
    &                     + BET2(J) * (F(I,JE,K) - F1(I,JE,K))
  ENDDO
  ENDDO
  ENDDO
  ! 三段目 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    F(I,J,K) = F2(I,J,K) + GAM1(K) * (F(I,J,KS) - F2(I,J,KS)) &
    &                    + GAM2(K) * (F(I,J,KE) - F2(I,J,KE))
  ENDDO
  ENDDO
  ENDDO
  ! メモリ解放 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DEALLOCATE(F1, F2)
  ! 処理終了 ***********************************************************
  RETURN
! 定義終了 *************************************************************
END SUBROUTINE Transfinite3D
!***********************************************************************
!**** 二次元一般座標変換メトリックス計算                            ****
!***********************************************************************
SUBROUTINE Metrics2D( &
&            IS, IE, JS, JE, &
&            X, Y, XIX, XIY, ETX, ETY, AJA &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  REAL,    INTENT(IN)  :: X(IS:IE, JS:JE), Y(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: XIX(IS:IE, JS:JE), XIY(IS:IE, JS:JE), &
  &                       ETX(IS:IE, JS:JE), ETY(IS:IE, JS:JE), &
  &                       AJA(IS:IE, JS:JE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I,J
  REAL    :: XXI(IS:IE,JS:JE), XET(IS:IE,JS:JE), &
  &          YXI(IS:IE,JS:JE), YET(IS:IE,JS:JE)
  REAL    :: AJAINV
  ! 処理開始 ***********************************************************
  ! X_XI, Y_XI の計算 --------------------------------------------------
  DO J=JS,JE
    I=IS
      XXI(I,J) = (-3.0 * X(I,J) + 4.0 * X(I+1,J) - X(I+2,J)) * 0.5
      YXI(I,J) = (-3.0 * Y(I,J) + 4.0 * Y(I+1,J) - Y(I+2,J)) * 0.5
    DO I=IS+1,IE-1
      XXI(I,J) = (-X(I-1,J) + X(I+1,J)) * 0.5
      YXI(I,J) = (-Y(I-1,J) + Y(I+1,J)) * 0.5
    ENDDO
    I=IE
      XXI(I,J) = (X(I-2,J) - 4.0 * X(I-1,J) + 3.0 * X(I,J)) * 0.5
      YXI(I,J) = (Y(I-2,J) - 4.0 * Y(I-1,J) + 3.0 * Y(I,J)) * 0.5
  ENDDO
  ! X_ET, Y_ET の計算 --------------------------------------------------
  DO I=IS,IE
    J=JS
      XET(I,J) = (-3.0 * X(I,J) + 4.0 * X(I,J+1) - X(I,J+2)) * 0.5
      YET(I,J) = (-3.0 * Y(I,J) + 4.0 * Y(I,J+1) - Y(I,J+2)) * 0.5
    DO J=JS+1,JE-1
      XET(I,J) = (-X(I,J-1) + X(I,J+1)) * 0.5
      YET(I,J) = (-Y(I,J-1) + Y(I,J+1)) * 0.5
    ENDDO
    J=JE
      XET(I,J) = (X(I,J-2) - 4.0 * X(I,J-1) + 3.0 * X(I,J)) * 0.5
      YET(I,J) = (Y(I,J-2) - 4.0 * Y(I,J-1) + 3.0 * Y(I,J)) * 0.5
  ENDDO
  ! Jacobian , metrics の計算 ------------------------------------------
  DO J=JS,JE
  DO I=IS,IE
    AJAINV = XXI(I,J) * YET(I,J) - XET(I,J) * YXI(I,J)
    IF(AJAINV .GT. 0.0) THEN
      AJA(I,J) = 1.0 / AJAINV
    ELSE
      AJA(I,J) = 0.0
    ENDIF
    XIX(I,J) = AJA(I,J) * YET(I,J)
    XIY(I,J) =-AJA(I,J) * XET(I,J)
    ETX(I,J) =-AJA(I,J) * YXI(I,J)
    ETY(I,J) = AJA(I,J) * XXI(I,J)
  ENDDO
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE Metrics2D
!***********************************************************************
!**** 三次元一般座標変換メトリックス計算                            ****
!***********************************************************************
SUBROUTINE Metrics3D( &
&            IS, IE, JS, JE, KS, KE, &
&            X, Y, Z, XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, AJA &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  REAL,    INTENT(IN)  :: X(IS:IE, JS:JE, KS:KE), &
  &                       Y(IS:IE, JS:JE, KS:KE), &
  &                       Z(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: XIX(IS:IE, JS:JE, KS:KE), &
  &                       XIY(IS:IE, JS:JE, KS:KE), &
  &                       XIZ(IS:IE, JS:JE, KS:KE), &
  &                       ETX(IS:IE, JS:JE, KS:KE), &
  &                       ETY(IS:IE, JS:JE, KS:KE), &
  &                       ETZ(IS:IE, JS:JE, KS:KE), &
  &                       ZEX(IS:IE, JS:JE, KS:KE), &
  &                       ZEY(IS:IE, JS:JE, KS:KE), &
  &                       ZEZ(IS:IE, JS:JE, KS:KE), &
  &                       AJA(IS:IE, JS:JE, KS:KE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I,J,K
  REAL    :: XXI(IS:IE,JS:JE,KS:KE), XET(IS:IE,JS:JE,KS:KE), &
  &          XZE(IS:IE,JS:JE,KS:KE), YXI(IS:IE,JS:JE,KS:KE), &
  &          YET(IS:IE,JS:JE,KS:KE), YZE(IS:IE,JS:JE,KS:KE), &
  &          ZXI(IS:IE,JS:JE,KS:KE), ZET(IS:IE,JS:JE,KS:KE), &
  &          ZZE(IS:IE,JS:JE,KS:KE)
  REAL    :: AJAINV
  ! 処理開始 ***********************************************************
  ! X_XI, Y_XI, Z_XI の計算 --------------------------------------------
  DO K=KS,KE
  DO J=JS,JE
    I=IS
      XXI(I,J,K) = ( -3.0*X(I,J,K)+4.0*X(I+1,J,K)-X(I+2,J,K) )*0.5
      YXI(I,J,K) = ( -3.0*Y(I,J,K)+4.0*Y(I+1,J,K)-Y(I+2,J,K) )*0.5
      ZXI(I,J,K) = ( -3.0*Z(I,J,K)+4.0*Z(I+1,J,K)-Z(I+2,J,K) )*0.5
    DO I=IS+1,IE-1
      XXI(I,J,K) = ( -X(I-1,J,K)+X(I+1,J,K) )*0.5
      YXI(I,J,K) = ( -Y(I-1,J,K)+Y(I+1,J,K) )*0.5
      ZXI(I,J,K) = ( -Z(I-1,J,K)+Z(I+1,J,K) )*0.5
    ENDDO
    I=IE
      XXI(I,J,K) = ( X(I-2,J,K)-4.0*X(I-1,J,K)+3.0*X(I,J,K) )*0.5
      YXI(I,J,K) = ( Y(I-2,J,K)-4.0*Y(I-1,J,K)+3.0*Y(I,J,K) )*0.5
      ZXI(I,J,K) = ( Z(I-2,J,K)-4.0*Z(I-1,J,K)+3.0*Z(I,J,K) )*0.5
  ENDDO
  ENDDO
  ! X_ET, Y_ET, Z_ET の計算 --------------------------------------------
  DO K=KS,KE
  DO I=IS,IE
    J=JS
      XET(I,J,K) = ( -3.0*X(I,J,K)+4.0*X(I,J+1,K)-X(I,J+2,K) )*0.5
      YET(I,J,K) = ( -3.0*Y(I,J,K)+4.0*Y(I,J+1,K)-Y(I,J+2,K) )*0.5
      ZET(I,J,K) = ( -3.0*Z(I,J,K)+4.0*Z(I,J+1,K)-Z(I,J+2,K) )*0.5
    DO J=JS+1,JE-1
      XET(I,J,K) = ( -X(I,J-1,K)+X(I,J+1,K) )*0.5
      YET(I,J,K) = ( -Y(I,J-1,K)+Y(I,J+1,K) )*0.5
      ZET(I,J,K) = ( -Z(I,J-1,K)+Z(I,J+1,K) )*0.5
    ENDDO
    J=JE
      XET(I,J,K) = ( X(I,J-2,K)-4.0*X(I,J-1,K)+3.0*X(I,J,K) )*0.5
      YET(I,J,K) = ( Y(I,J-2,K)-4.0*Y(I,J-1,K)+3.0*Y(I,J,K) )*0.5
      ZET(I,J,K) = ( Z(I,J-2,K)-4.0*Z(I,J-1,K)+3.0*Z(I,J,K) )*0.5
  ENDDO
  ENDDO
  ! X_ZE, Y_ZE, Z_ZE の計算 --------------------------------------------
  DO J=JS,JE
  DO I=IS,IE
    K=KS
      XZE(I,J,K) = ( -3.0*X(I,J,K)+4.0*X(I,J,K+1)-X(I,J,K+2) )*0.5
      YZE(I,J,K) = ( -3.0*Y(I,J,K)+4.0*Y(I,J,K+1)-Y(I,J,K+2) )*0.5
      ZZE(I,J,K) = ( -3.0*Z(I,J,K)+4.0*Z(I,J,K+1)-Z(I,J,K+2) )*0.5
    DO K=KS+1,KE-1
      XZE(I,J,K) = ( -X(I,J,K-1)+X(I,J,K+1) )*0.5
      YZE(I,J,K) = ( -Y(I,J,K-1)+Y(I,J,K+1) )*0.5
      ZZE(I,J,K) = ( -Z(I,J,K-1)+Z(I,J,K+1) )*0.5
    ENDDO
    K=KE
      XZE(I,J,K) = ( X(I,J,K-2)-4.0*X(I,J,K-1)+3.0*X(I,J,K) )*0.5
      YZE(I,J,K) = ( Y(I,J,K-2)-4.0*Y(I,J,K-1)+3.0*Y(I,J,K) )*0.5
      ZZE(I,J,K) = ( Z(I,J,K-2)-4.0*Z(I,J,K-1)+3.0*Z(I,J,K) )*0.5
  ENDDO
  ENDDO
  ! Jacobian , metrics の計算 ------------------------------------------
  DO K=KS,KE
  DO J=JS,JE
  DO I=IS,IE
    AJAINV = XXI(I,J,K) * ( YET(I,J,K) * ZZE(I,J,K) &
    &                     - YZE(I,J,K) * ZET(I,J,K) ) &
    &      + XET(I,J,K) * ( YZE(I,J,K) * ZXI(I,J,K) &
    &                     - YXI(I,J,K) * ZZE(I,J,K) ) &
    &      + XZE(I,J,K) * ( YXI(I,J,K) * ZET(I,J,K) &
    &                     - YET(I,J,K) * ZXI(I,J,K) )
    IF(AJAINV .GT. 0.0) THEN
      AJA(I,J,K) = 1.0 / AJAINV
    ELSE
      AJA(I,J,K) = 0.0
    ENDIF
    XIX(I,J,K) = AJA(I,J,K) * ( YET(I,J,K) * ZZE(I,J,K) &
    &                         - YZE(I,J,K) * ZET(I,J,K) )
    XIY(I,J,K) = AJA(I,J,K) * ( ZET(I,J,K) * XZE(I,J,K) &
    &                         - ZZE(I,J,K) * XET(I,J,K) )
    XIZ(I,J,K) = AJA(I,J,K) * ( XET(I,J,K) * YZE(I,J,K) &
    &                         - XZE(I,J,K) * YET(I,J,K) )
    ETX(I,J,K) = AJA(I,J,K) * ( YZE(I,J,K) * ZXI(I,J,K) &
    &                         - YXI(I,J,K) * ZZE(I,J,K) )
    ETY(I,J,K) = AJA(I,J,K) * ( ZZE(I,J,K) * XXI(I,J,K) &
    &                         - ZXI(I,J,K) * XZE(I,J,K) )
    ETZ(I,J,K) = AJA(I,J,K) * ( XZE(I,J,K) * YXI(I,J,K) &
    &                         - XXI(I,J,K) * YZE(I,J,K) )
    ZEX(I,J,K) = AJA(I,J,K) * ( YXI(I,J,K) * ZET(I,J,K) &
    &                         - YET(I,J,K) * ZXI(I,J,K) )
    ZEY(I,J,K) = AJA(I,J,K) * ( ZXI(I,J,K) * XET(I,J,K) &
    &                         - ZET(I,J,K) * XXI(I,J,K) )
    ZEZ(I,J,K) = AJA(I,J,K) * ( XXI(I,J,K) * YET(I,J,K) &
    &                         - XET(I,J,K) * YXI(I,J,K) )
  ENDDO
  ENDDO
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE Metrics3D
!***********************************************************************
!**** 壁からの距離を計算(二次元)                                    ****
!***********************************************************************
SUBROUTINE WallDistance2D( &
&            IS, IE, JS, JE, &
&            WallNum, WallX, WallY, X, Y, &
&            YWALL &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: WallNum
  REAL,    INTENT(IN)  :: WallX(WallNum), WallY(WallNum)
  REAL,    INTENT(IN)  :: X(IS:IE ,JS:JE), Y(IS:IE ,JS:JE)
  REAL,    INTENT(OUT) :: YWALL(IS:IE ,JS:JE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, N
  ! 処理開始 ***********************************************************
  DO J=JS,JE
  DO I=IS,IE
    N=1
      YWALL(I,J) = SQRT( (X(I,J) - WallX(N) )**2 &
      &                + (Y(I,J) - WallY(N) )**2 )
    DO N=1,WallNum
      YWALL(I,J) = MIN( &
      &            YWALL(I,J), &
      &            SQRT( (X(I,J) - WallX(N) )**2 &
      &                + (Y(I,J) - WallY(N) )**2 ) &
      &          )
    ENDDO
  ENDDO
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE WallDistance2D
!***********************************************************************
!**** 壁からの距離を計算(三次元)                                    ****
!***********************************************************************
SUBROUTINE WallDistance3D( &
&            IS, IE, JS, JE, KS, KE, &
&            WallNum, WallX, WallY, WallZ, X, Y, Z, &
&            YWALL &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: WallNum
  REAL,    INTENT(IN)  :: WallX(WallNum), WallY(WallNum), WallZ(WallNum)
  REAL,    INTENT(IN)  :: X(IS:IE ,JS:JE, KS:KE), &
  &                       Y(IS:IE ,JS:JE, KS:KE), &
  &                       Z(IS:IE ,JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: YWALL(IS:IE ,JS:JE, KS:KE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, N
  ! 処理開始 ***********************************************************
  DO K=KS,KE
  DO J=JS,JE
  DO I=IS,IE
    N=1
      YWALL(I,J,K) = SQRT( (X(I,J,K) - WallX(N) )**2 &
      &                  + (Y(I,J,K) - WallY(N) )**2 &
      &                  + (Z(I,J,K) - WallZ(N) )**2 )
    DO N=1,WallNum
      YWALL(I,J,K) = MIN( YWALL(I,J,K) , &
      &              SQRT( (X(I,J,K) - WallX(N) )**2 &
      &                  + (Y(I,J,K) - WallY(N) )**2 &
      &                  + (Z(I,J,K) - WallZ(N) )**2 ) &
      &            )
    ENDDO
  ENDDO
  ENDDO
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE WallDistance3D
!***********************************************************************
!**** ワイドサーチ(二次元, 四点)                                    ****
!***********************************************************************
SUBROUTINE WideSearch2D4Point( &
&            X0, Y0, X1, Y1, X2, Y2, X3, Y3, X4, Y4, &
&            MARGIN, &
&            FEXIST &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: X0, Y0
  REAL,    INTENT(IN)  :: X1, Y1
  REAL,    INTENT(IN)  :: X2, Y2
  REAL,    INTENT(IN)  :: X3, Y3
  REAL,    INTENT(IN)  :: X4, Y4
  REAL,    INTENT(IN)  :: MARGIN
  LOGICAL, INTENT(OUT) :: FEXIST
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL :: XMIN, XMAX, YMIN, YMAX
  REAL :: MARGINX, MARGINY
  ! 処理開始 ***********************************************************
  XMIN = MIN( X1, X2, X3, X4 )
  XMAX = MAX( X1, X2, X3, X4 )
  YMIN = MIN( Y1, Y2, Y3, Y4 )
  YMAX = MAX( Y1, Y2, Y3, Y4 )
  MARGINX = (XMAX - XMIN) * MARGIN
  MARGINY = (YMAX - YMIN) * MARGIN
  FEXIST = (XMIN - MARGINX .LE. X0) .AND. (X0 .LE. XMAX + MARGINX) &
  &  .AND. (YMIN - MARGINY .LE. Y0) .AND. (Y0 .LE. YMAX + MARGINY)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE WideSearch2D4Point
!***********************************************************************
!**** ワイドサーチ(三次元, 八点)                                    ****
!***********************************************************************
SUBROUTINE WideSearch3D8Point( &
&            X0, Y0, Z0, &
&            X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, &
&            X5, Y5, Z5, X6, Y6, Z6, X7, Y7, Z7, X8, Y8, Z8, &
&            MARGIN, &
&            FEXIST &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: X0, Y0, Z0
  REAL,    INTENT(IN)  :: X1, Y1, Z1
  REAL,    INTENT(IN)  :: X2, Y2, Z2
  REAL,    INTENT(IN)  :: X3, Y3, Z3
  REAL,    INTENT(IN)  :: X4, Y4, Z4
  REAL,    INTENT(IN)  :: X5, Y5, Z5
  REAL,    INTENT(IN)  :: X6, Y6, Z6
  REAL,    INTENT(IN)  :: X7, Y7, Z7
  REAL,    INTENT(IN)  :: X8, Y8, Z8
  REAL,    INTENT(IN)  :: MARGIN
  LOGICAL, INTENT(OUT) :: FEXIST
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
  REAL :: MARGINX, MARGINY, MARGINZ
  ! 処理開始 ***********************************************************
  XMIN = MIN( X1, X2, X3, X4, X5, X6, X7, X8 )
  XMAX = MAX( X1, X2, X3, X4, X5, X6, X7, X8 )
  YMIN = MIN( Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8 )
  YMAX = MAX( Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8 )
  ZMIN = MIN( Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8 )
  ZMAX = MAX( Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8 )
  MARGINX = (XMAX - XMIN) * MARGIN
  MARGINY = (YMAX - YMIN) * MARGIN
  MARGINZ = (ZMAX - ZMIN) * MARGIN
  FEXIST = (XMIN - MARGINX .LE. X0) .AND. (X0 .LE. XMAX + MARGINX) &
  &  .AND. (YMIN - MARGINY .LE. Y0) .AND. (Y0 .LE. YMAX + MARGINY) &
  &  .AND. (ZMIN - MARGINZ .LE. Z0) .AND. (Z0 .LE. ZMAX + MARGINZ)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE WideSearch3D8Point
!***********************************************************************
!**** 線形補間(二次元, 三点)                                        ****
!****   X0 = X1 + ALP (X2 - X1) + BET (X3 - X1)                     ****
!****   Y0 = Y1 + ALP (Y2 - Y1) + BET (Y3 - Y1)                     ****
!***********************************************************************
SUBROUTINE Interpolation2D3Point( &
&            X0, Y0, X1, Y1, X2, Y2, X3, Y3, &
&            MARGIN, &
&            ALP, BET, FCOMP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: X0, Y0
  REAL,    INTENT(IN)  :: X1, Y1
  REAL,    INTENT(IN)  :: X2, Y2
  REAL,    INTENT(IN)  :: X3, Y3
  REAL,    INTENT(IN)  :: MARGIN
  REAL,    INTENT(OUT) :: ALP, BET
  LOGICAL, INTENT(OUT) :: FCOMP
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL :: deti, detx, dety
  ! 処理開始 ***********************************************************
  deti = (X2 - X1) * (Y3 - Y1) - (Y2 - Y1) * (X3 - X1)
  detx = (X0 - X1) * (Y3 - Y1) - (Y0 - Y1) * (X3 - X1)
  dety = (X2 - X1) * (Y0 - Y1) - (Y2 - Y1) * (X0 - X1)
  IF(deti .NE. 0.0) THEN
    ALP = detx / deti
    BET = dety / deti
    FCOMP = 0.0 - MARGIN .LE. ALP .AND. ALP .LE. 1.0 + MARGIN .AND. &
    &       0.0 - MARGIN .LE. BET .AND. BET .LE. 1.0 + MARGIN .AND. &
    &       ALP + BET .LE. 1.0 + MARGIN
  ELSE
    ALP = 0.0
    BET = 0.0
    FCOMP = .FALSE.
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE Interpolation2D3Point
!***********************************************************************
!**** 線形補間(三次元, 四点)                                        ****
!****   X0 = X1 + ALP (X2 - X1) + BET (X3 - X1) + GAM (X4 - X1)     ****
!****   Y0 = Y1 + ALP (Y2 - Y1) + BET (Y3 - Y1) + GAM (Y4 - Y1)     ****
!****   Z0 = Z1 + ALP (Z2 - Z1) + BET (Z3 - Z1) + GAM (Z4 - Z1)     ****
!***********************************************************************
SUBROUTINE Interpolation3D4Point( &
&            X0, Y0, Z0, &
&            X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, &
&            MARGIN, &
&            ALP, BET, GAM, FCOMP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: X0, Y0, Z0
  REAL,    INTENT(IN)  :: X1, Y1, Z1
  REAL,    INTENT(IN)  :: X2, Y2, Z2
  REAL,    INTENT(IN)  :: X3, Y3, Z3
  REAL,    INTENT(IN)  :: X4, Y4, Z4
  REAL,    INTENT(IN)  :: MARGIN
  REAL,    INTENT(OUT) :: ALP, BET, GAM
  LOGICAL, INTENT(OUT) :: FCOMP
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL :: dx0, dx2, dx3, dx4
  REAL :: dy0, dy2, dy3, dy4
  REAL :: dz0, dz2, dz3, dz4
  REAL :: deti, detx, dety, detz
  ! 処理開始 ***********************************************************
  dx0 = X0 - X1; dx2 = X2 - X1; dx3 = X3 - X1; dx4 = X4 - X1
  dy0 = Y0 - Y1; dy2 = Y2 - Y1; dy3 = Y3 - Y1; dy4 = Y4 - Y1
  dz0 = Z0 - Z1; dz2 = Z2 - Z1; dz3 = Z3 - Z1; dz4 = Z4 - Z1
  deti = dx2 * dy3 * dz4 + dx3 * dy4 * dz2 + dx4 * dy2 * dz3 &
  &    - dz2 * dy3 * dx4 - dz3 * dy4 * dx2 - dz4 * dy2 * dx3
  detx = dx0 * dy3 * dz4 + dx3 * dy4 * dz0 + dx4 * dy0 * dz3 &
  &    - dz0 * dy3 * dx4 - dz3 * dy4 * dx0 - dz4 * dy0 * dx3
  dety = dx2 * dy0 * dz4 + dx0 * dy4 * dz2 + dx4 * dy2 * dz0 &
  &    - dz2 * dy0 * dx4 - dz0 * dy4 * dx2 - dz4 * dy2 * dx0
  detz = dx2 * dy3 * dz0 + dx3 * dy0 * dz2 + dx0 * dy2 * dz3 &
  &    - dz2 * dy3 * dx0 - dz3 * dy0 * dx2 - dz0 * dy2 * dx3
  IF(deti .NE. 0.0) THEN
    ALP = detx / deti
    BET = dety / deti
    GAM = detz / deti
    FCOMP = 0.0 - MARGIN .LE. ALP .AND. ALP .LE. 1.0 + MARGIN .AND. &
    &       0.0 - MARGIN .LE. BET .AND. BET .LE. 1.0 + MARGIN .AND. &
    &       0.0 - MARGIN .LE. GAM .AND. GAM .LE. 1.0 + MARGIN .AND. &
    &       ALP + BET + GAM .LE. 1.0 + MARGIN
  ELSE
    ALP = 0.0
    BET = 0.0
    GAM = 0.0
    FCOMP = .FALSE.
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE Interpolation3D4Point
!***********************************************************************
!**** 双一次補間(二次元, 四点)                                      ****
!****   X0 = (1 - ALP) (1 - BET) X1 + ALP (1 - BET) X2              ****
!****      + (1 - ALP) BET       X3 + ALP BET       X4              ****
!****   Y0 = (1 - ALP) (1 - BET) Y1 + ALP (1 - BET) Y2              ****
!****      + (1 - ALP) BET       Y3 + ALP BET       Y4              ****
!***********************************************************************
SUBROUTINE Interpolation2D4Point( &
&            X0, Y0, X1, Y1, X2, Y2, X3, Y3, X4, Y4, &
&            MARGIN, &
&            ALP0, BET0, FCOMP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, PARAMETER :: NMAX = 100
  REAL(8), PARAMETER :: DMIN = 1.0D-8
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: X0, Y0
  REAL,    INTENT(IN)    :: X1, Y1
  REAL,    INTENT(IN)    :: X2, Y2
  REAL,    INTENT(IN)    :: X3, Y3
  REAL,    INTENT(IN)    :: X4, Y4
  REAL,    INTENT(IN)    :: MARGIN
  REAL,    INTENT(INOUT) :: ALP0, BET0
  LOGICAL, INTENT(OUT)   :: FCOMP
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL(8) :: XX1, XX2, XX4, XX8, YY1, YY2, YY4, YY8
  REAL(8) :: ALP, BET
  INTEGER :: N
  REAL(8) :: FX, FY
  REAL(8) :: FXA, FXB, FYA, FYB
  REAL(8) :: DELTA, DA, DB, DABMIN, DABMAX, DMAX
  REAL(8) :: CRLX
  REAL(8) :: MARGIND
  ! 処理開始 ***********************************************************
  ! 初期値 -------------------------------------------------------------
  ALP   = DBLE(ALP0)
  BET   = DBLE(BET0)
  FCOMP = .FALSE.
  ! 補間関数用定数設定 -------------------------------------------------
  XX1 = DBLE(X2) - DBLE(X1)
  XX2 = DBLE(X3) - DBLE(X1)
  XX4 = DBLE(X1) - DBLE(X2) - DBLE(X3) + DBLE(X4)
  XX8 = DBLE(X1) - DBLE(X0)
  YY1 = DBLE(Y2) - DBLE(Y1)
  YY2 = DBLE(Y3) - DBLE(Y1)
  YY4 = DBLE(Y1) - DBLE(Y2) - DBLE(Y3) + DBLE(Y4)
  YY8 = DBLE(Y1) - DBLE(Y0)
  ! Newton-Raphson Method で、補間パラメータを計算 ---------------------
  CRLX = 1.0
  DO N = 1, NMAX
    FX = XX1*ALP + XX2*BET + XX4*ALP*BET + XX8
    FY = YY1*ALP + YY2*BET + YY4*ALP*BET + YY8
    FXA = XX1 + XX4*BET
    FXB = XX2 + XX4*ALP
    FYA = YY1 + YY4*BET
    FYB = YY2 + YY4*ALP
    DELTA = FXA*FYB - FXB*FYA
    IF(ABS(DELTA) .LE. 1.0D-20) DELTA = SIGN(1.0D-20, DELTA)
    DA =-( FX*FYB - FXB*FY ) / DELTA
    DB =-( FXA*FY - FX*FYA ) / DELTA
    ALP = ALP + DA * CRLX
    BET = BET + DB * CRLX
    DABMIN = MIN(ABS(DA), ABS(DB))
    DABMAX = MAX(ABS(DA), ABS(DB))
    SELECT CASE(N)
      CASE(1)
        DMAX = DABMIN
      CASE DEFAULT
        IF(DMAX .LT. DABMIN) THEN
          DMAX = DABMIN
          CRLX = CRLX * 0.5
        ENDIF
    END SELECT
    IF(DABMAX .LE. DMIN) EXIT
  ENDDO
  ! 正当性検査 ---------------------------------------------------------
  MARGIND = DBLE(MARGIN)
  IF( DABMAX .LE. DMIN .AND. &
  &   0.0D+0 - MARGIND .LE. ALP .AND. ALP .LT. 1.0D+0 + MARGIND .AND. &
  &   0.0D+0 - MARGIND .LE. BET .AND. BET .LT. 1.0D+0 + MARGIND ) THEN
    FCOMP = .TRUE.
  ENDIF
  ALP0 = REAL(ALP)
  BET0 = REAL(BET)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE Interpolation2D4Point
!***********************************************************************
!**** 三重線形補間(三次元, 八点)                                    ****
!****   X0 = (1 - ALP) (1 - BET) (1 - GAM) X1                       ****
!****      + ALP       (1 - BET) (1 - GAM) X2                       ****
!****      + (1 - ALP) BET       (1 - GAM) X3                       ****
!****      + ALP       BET       (1 - GAM) X4                       ****
!****      + (1 - ALP) (1 - BET) GAM       X5                       ****
!****      + ALP       (1 - BET) GAM       X6                       ****
!****      + (1 - ALP) BET       GAM       X7                       ****
!****      + ALP       BET       GAM       X8                       ****
!****   Y0 = (1 - ALP) (1 - BET) (1 - GAM) Y1                       ****
!****      + ALP       (1 - BET) (1 - GAM) Y2                       ****
!****      + (1 - ALP) BET       (1 - GAM) Y3                       ****
!****      + ALP       BET       (1 - GAM) Y4                       ****
!****      + (1 - ALP) (1 - BET) GAM       Y5                       ****
!****      + ALP       (1 - BET) GAM       Y6                       ****
!****      + (1 - ALP) BET       GAM       Y7                       ****
!****      + ALP       BET       GAM       Y8                       ****
!****   Z0 = (1 - ALP) (1 - BET) (1 - GAM) Z1                       ****
!****      + ALP       (1 - BET) (1 - GAM) Z2                       ****
!****      + (1 - ALP) BET       (1 - GAM) Z3                       ****
!****      + ALP       BET       (1 - GAM) Z4                       ****
!****      + (1 - ALP) (1 - BET) GAM       Z5                       ****
!****      + ALP       (1 - BET) GAM       Z6                       ****
!****      + (1 - ALP) BET       GAM       Z7                       ****
!****      + ALP       BET       GAM       Z8                       ****
!***********************************************************************
SUBROUTINE Interpolation3D8Point( &
&            X0, Y0, Z0, &
&            X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, &
&            X5, Y5, Z5, X6, Y6, Z6, X7, Y7, Z7, X8, Y8, Z8, &
&            MARGIN, &
&            ALP0, BET0, GAM0, FCOMP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, PARAMETER :: NMAX = 100
  REAL(8), PARAMETER :: DMIN = 1.0D-8
!  REAL(8), PARAMETER :: DMIN = 1.0D-10
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: X0, Y0, Z0
  REAL,    INTENT(IN)    :: X1, Y1, Z1
  REAL,    INTENT(IN)    :: X2, Y2, Z2
  REAL,    INTENT(IN)    :: X3, Y3, Z3
  REAL,    INTENT(IN)    :: X4, Y4, Z4
  REAL,    INTENT(IN)    :: X5, Y5, Z5
  REAL,    INTENT(IN)    :: X6, Y6, Z6
  REAL,    INTENT(IN)    :: X7, Y7, Z7
  REAL,    INTENT(IN)    :: X8, Y8, Z8
  REAL,    INTENT(IN)    :: MARGIN
  REAL,    INTENT(INOUT) :: ALP0, BET0, GAM0
  LOGICAL, INTENT(OUT)   :: FCOMP
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL(8) :: XX1, XX2, XX3, XX4, XX5, XX6, XX7, XX8, &
  &          YY1, YY2, YY3, YY4, YY5, YY6, YY7, YY8, &
  &          ZZ1, ZZ2, ZZ3, ZZ4, ZZ5, ZZ6, ZZ7, ZZ8
  REAL(8) :: ALP, BET, GAM
  INTEGER :: N
  REAL(8) :: FX, FY, FZ
  REAL(8) :: FXA, FXB, FXG, FYA, FYB, FYG, FZA, FZB, FZG
  REAL(8) :: DELTA, DA, DB, DG, DABGMIN, DABGMAX, DMAX
  REAL(8) :: CRLX
  REAL(8) :: MARGIND
  ! 処理開始 ***********************************************************
  ! 初期値 -------------------------------------------------------------
  ALP   = DBLE(ALP0)
  BET   = DBLE(BET0)
  GAM   = DBLE(GAM0)
  FCOMP = .FALSE.
  ! 補間関数用定数設定 -------------------------------------------------
  XX1 = DBLE(X2) - DBLE(X1)
  XX2 = DBLE(X3) - DBLE(X1)
  XX3 = DBLE(X5) - DBLE(X1)
  XX4 = DBLE(X1) - DBLE(X2) - DBLE(X3) + DBLE(X4)
  XX5 = DBLE(X1) - DBLE(X3) - DBLE(X5) + DBLE(X7)
  XX6 = DBLE(X1) - DBLE(X2) - DBLE(X5) + DBLE(X6)
  XX7 =-DBLE(X1) + DBLE(X2) + DBLE(X3) - DBLE(X4) &
  &   + DBLE(X5) - DBLE(X6) - DBLE(X7) + DBLE(X8)
  XX8 = DBLE(X1) - DBLE(X0)
  YY1 = DBLE(Y2) - DBLE(Y1)
  YY2 = DBLE(Y3) - DBLE(Y1)
  YY3 = DBLE(Y5) - DBLE(Y1)
  YY4 = DBLE(Y1) - DBLE(Y2) - DBLE(Y3) + DBLE(Y4)
  YY5 = DBLE(Y1) - DBLE(Y3) - DBLE(Y5) + DBLE(Y7)
  YY6 = DBLE(Y1) - DBLE(Y2) - DBLE(Y5) + DBLE(Y6)
  YY7 =-DBLE(Y1) + DBLE(Y2) + DBLE(Y3) - DBLE(Y4) &
  &   + DBLE(Y5) - DBLE(Y6) - DBLE(Y7) + DBLE(Y8)
  YY8 = DBLE(Y1) - DBLE(Y0)
  ZZ1 = DBLE(Z2) - DBLE(Z1)
  ZZ2 = DBLE(Z3) - DBLE(Z1)
  ZZ3 = DBLE(Z5) - DBLE(Z1)
  ZZ4 = DBLE(Z1) - DBLE(Z2) - DBLE(Z3) + DBLE(Z4)
  ZZ5 = DBLE(Z1) - DBLE(Z3) - DBLE(Z5) + DBLE(Z7)
  ZZ6 = DBLE(Z1) - DBLE(Z2) - DBLE(Z5) + DBLE(Z6)
  ZZ7 =-DBLE(Z1) + DBLE(Z2) + DBLE(Z3) - DBLE(Z4) &
  &   + DBLE(Z5) - DBLE(Z6) - DBLE(Z7) + DBLE(Z8)
  ZZ8 = DBLE(Z1) - DBLE(Z0)
  ! Newton-Raphson Method で、補間パラメータを計算 ---------------------
  CRLX = 1.0
  DO N = 1, NMAX
    FX = XX1*ALP + XX2*BET + XX3*GAM + XX4*ALP*BET &
    &  + XX5*BET*GAM + XX6*ALP*GAM + XX7*ALP*BET*GAM + XX8
    FY = YY1*ALP + YY2*BET + YY3*GAM + YY4*ALP*BET &
    &  + YY5*BET*GAM + YY6*ALP*GAM + YY7*ALP*BET*GAM + YY8
    FZ = ZZ1*ALP + ZZ2*BET + ZZ3*GAM + ZZ4*ALP*BET &
    &  + ZZ5*BET*GAM + ZZ6*ALP*GAM + ZZ7*ALP*BET*GAM + ZZ8
!    FXA = XX1 + XX4*BET + XX6*GAM + XX8*BET*GAM
!    FXB = XX2 + XX4*ALP + XX5*GAM + XX8*ALP*GAM
!    FXG = XX3 + XX5*BET + XX6*ALP + XX8*ALP*BET
!    FYA = YY1 + YY4*BET + YY6*GAM + YY8*BET*GAM
!    FYB = YY2 + YY4*ALP + YY5*GAM + YY8*ALP*GAM
!    FYG = YY3 + YY5*BET + YY6*ALP + YY8*ALP*BET
!    FZA = ZZ1 + ZZ4*BET + ZZ6*GAM + ZZ8*BET*GAM
!    FZB = ZZ2 + ZZ4*ALP + ZZ5*GAM + ZZ8*ALP*GAM
!    FZG = ZZ3 + ZZ5*BET + ZZ6*ALP + ZZ8*ALP*BET
    FXA = XX1 + XX4*BET + XX6*GAM + XX7*BET*GAM
    FXB = XX2 + XX4*ALP + XX5*GAM + XX7*ALP*GAM
    FXG = XX3 + XX5*BET + XX6*ALP + XX7*ALP*BET
    FYA = YY1 + YY4*BET + YY6*GAM + YY7*BET*GAM
    FYB = YY2 + YY4*ALP + YY5*GAM + YY7*ALP*GAM
    FYG = YY3 + YY5*BET + YY6*ALP + YY7*ALP*BET
    FZA = ZZ1 + ZZ4*BET + ZZ6*GAM + ZZ7*BET*GAM
    FZB = ZZ2 + ZZ4*ALP + ZZ5*GAM + ZZ7*ALP*GAM
    FZG = ZZ3 + ZZ5*BET + ZZ6*ALP + ZZ7*ALP*BET
    DELTA = FXA*(FYB*FZG-FYG*FZB) - FXB*(FYA*FZG-FYG*FZA) &
    &     + FXG*(FYA*FZB-FYB*FZA)
    IF(ABS(DELTA) .LE. 1.0D-20) DELTA = SIGN(1.0D-20, DELTA)
    IF(ABS(DELTA) .GE. 1.0D+5) EXIT
    DA =-( FX*(FYB*FZG-FYG*FZB) - FXB*(FY*FZG-FYG*FZ) &
    &    + FXG*(FY*FZB-FYB*FZ) ) / DELTA
    DB =-( FXA*(FY*FZG-FYG*FZ) - FX*(FYA*FZG-FYG*FZA) &
    &    + FXG*(FYA*FZ-FY*FZA) ) / DELTA
    DG =-( FXA*(FYB*FZ-FY*FZB) - FXB*(FYA*FZ-FY*FZA) &
    &    + FX*(FYA*FZB-FYB*FZA) ) / DELTA
    ALP = ALP + DA * CRLX
    BET = BET + DB * CRLX
    GAM = GAM + DG * CRLX
    DABGMIN = MIN(ABS(DA), ABS(DB), ABS(DG))
    DABGMAX = MAX(ABS(DA), ABS(DB), ABS(DG))
    IF( (ABS(DA) .LE. 1.0D-10) .AND. &
    &   (ABS(DB) .LE. 1.0D-10) .AND. &
    &   (ABS(DG) .LE. 1.0D-10) ) EXIT
    IF( (ABS(DA) .GE. 1.0D+5 ) .AND. &
    &   (ABS(DB) .GE. 1.0D+5 ) .AND. &
    &   (ABS(DG) .GE. 1.0D+5 ) ) EXIT
    SELECT CASE(N)
      CASE(1)
        DMAX = DABGMIN
      CASE DEFAULT
        IF(DMAX .LT. DABGMIN) THEN
          DMAX = DABGMIN
          CRLX = CRLX * 0.5
        ENDIF
    END SELECT
    IF(DABGMAX .LE. DMIN) EXIT
  ENDDO
  ! 正当性検査 ---------------------------------------------------------
  MARGIND = DBLE(MARGIN)
  IF( DABGMAX .LE. DMIN .AND. &
  &   0.0D+0 - MARGIND .LE. ALP .AND. ALP .LT. 1.0D+0 + MARGIND .AND. &
  &   0.0D+0 - MARGIND .LE. BET .AND. BET .LT. 1.0D+0 + MARGIND .AND. &
  &   0.0D+0 - MARGIND .LE. GAM .AND. GAM .LT. 1.0D+0 + MARGIND ) THEN
    FCOMP = .TRUE.
  ENDIF
  ALP0 = REAL(ALP)
  BET0 = REAL(BET)
  GAM0 = REAL(GAM)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE Interpolation3D8Point
!***********************************************************************
!**** 翼の中の点を探索                                              ****
!***********************************************************************
SUBROUTINE BladeInOut( &
&            IS1, IE1, JS1, JE1, KS1, KE1, X1, Y1, Z1,  &
&            IS2, IE2,           KS2, KE2, X2, Y2, Z2,  &
&            BIN )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! メイン・グリッド ---------------------------------------------------
  INTEGER, INTENT(IN)  :: IS1, IE1, JS1, JE1, KS1, KE1
  REAL   , INTENT(IN)  :: X1(IS1:IE1, JS1:JE1, KS1:KE1), &
&                         Y1(IS1:IE1, JS1:JE1, KS1:KE1), &
&                         Z1(IS1:IE1, JS1:JE1, KS1:KE1)
  ! サブ・グリッド -----------------------------------------------------
  INTEGER, INTENT(IN)  :: IS2, IE2, KS2, KE2
  REAL   , INTENT(IN)  :: X2(IS2:IE2, KS2:KE2), &
&                         Y2(IS2:IE2, KS2:KE2), &
&                         Z2(IS2:IE2, KS2:KE2)
  ! 翼の中か外の判定値 -------------------------------------------------
  INTEGER, INTENT(OUT) :: BIN(IS1:IE1, JS1:JE1, KS1:KE1)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, II, KK, IP, KP
  REAL    :: DL, DLMIN
  REAL    :: AX, AY, AZ, BX, BY, BZ, CX, CY, CZ, DX, DY, DZ
  REAL    :: HX1, HY1, HZ1, HX2, HY2, HZ2, PRO
  ! 処理開始 ***********************************************************
  DO K = KS1, KE1
  DO J = JS1, JE1
  DO I = IS1, IE1
  ! 翼面で最も近い点の算出 ---------------------------------------------
    DLMIN = 1.0E+10
    DO KK = KS2, KE2
    DO II = IS2, IE2
      DL = SQRT( (X2(II,KK) - X1(I,J,K))**2 &
      &        + (Y2(II,KK) - Y1(I,J,K))**2 &
      &        + (Z2(II,KK) - Z1(I,J,K))**2 )
      IF( DL .LT. DLMIN ) THEN
        DLMIN = DL
        IP = II
        KP = KK
      ENDIF
    END DO
    END DO
  ! 翼面の法線ベクトルを求め内積をとる ---------------------------------
  ! 翼表面からの単位法線ベクトル
    IF( KP .EQ. KS2 ) THEN
      AX = X2(IP,KP+1) - X2(IP,KP  )
      AY = Y2(IP,KP+1) - Y2(IP,KP  )
      AZ = Z2(IP,KP+1) - Z2(IP,KP  )
     ELSE IF( KP .EQ. KE2 ) THEN
      AX = X2(IP,KP  ) - X2(IP,KP-1)
      AY = Y2(IP,KP  ) - Y2(IP,KP-1)
      AZ = Z2(IP,KP  ) - Z2(IP,KP-1)
     ELSE
      AX = X2(IP,KP+1) - X2(IP,KP-1)
      AY = Y2(IP,KP+1) - Y2(IP,KP-1)
      AZ = Z2(IP,KP+1) - Z2(IP,KP-1)
    END IF
    IF( IP .EQ. IS2 )THEN
      BX = X2(IP+1,KP) - X2(IP  ,KP)
      BY = Y2(IP+1,KP) - Y2(IP  ,KP)
      BZ = Z2(IP+1,KP) - Z2(IP  ,KP)
     ELSE IF( IP .EQ. IE2 ) THEN
      BX = X2(IP  ,KP) - X2(IP-1,KP)
      BY = Y2(IP  ,KP) - Y2(IP-1,KP)
      BZ = Z2(IP  ,KP) - Z2(IP-1,KP)
     ELSE
      BX = X2(IP+1,KP) - X2(IP-1,KP)
      BY = Y2(IP+1,KP) - Y2(IP-1,KP)
      BZ = Z2(IP+1,KP) - Z2(IP-1,KP)
    END IF
    CX  = AY * BZ - AZ * BY
    CY  = AZ * BX - AX * BZ
    CZ  = AX * BY - AY * BX
    DL  = SQRT( CX**2 + CY**2 + CZ**2 )
    HX2 = CX / DL
    HY2 = CY / DL
    HZ2 = CZ / DL
  ! 翼近傍への単位ベクトル ---------------------------------------------
    DX  = X1(I,J,K) - X2(IP,KP)
    DY  = Y1(I,J,K) - Y2(IP,KP)
    DZ  = Z1(I,J,K) - Z2(IP,KP)
    DL  = SQRT( DX**2 + DY**2 + DZ**2 )
    HX1 = DX / DL
    HY1 = DY / DL
    HZ1 = DZ / DL
  ! 内積が正なら翼外:0 負なら翼内:1 ------------------------------------
    PRO = HX1 * HX2 + HY1 * HY2 + HZ1 * HZ2
    IF( PRO .GT. 0.0 )THEN
      BIN(I,J,K) = 0
     ELSE
      BIN(I,J,K) = 1
    END IF
  END DO
  END DO
  END DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE BladeInOut
! 定義終了 *************************************************************
END MODULE Package_Grid
