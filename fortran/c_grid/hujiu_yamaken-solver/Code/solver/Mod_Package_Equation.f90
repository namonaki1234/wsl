!***********************************************************************
!***********************************************************************
!**** パッケージ型モジュール                                        ****
!**** 方程式の解法及び近似用サブルーチン群                          ****
!****                         2008.12.01 PROGRAMED BY SUZUKI MASAYA ****
!****                         2011.10.15 UPDATED BY RYOSUKE HAYASHI ****
!***********************************************************************
!***********************************************************************
MODULE Package_Equation
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  PRIVATE
  ! サブルーチン宣言 ***************************************************
  ! 共有サブルーチン(正也さん作) +++++++++++++++++++++++++++++++++++++++
  PUBLIC :: Gauss
  PUBLIC :: NRM_SM
  PUBLIC :: LeastSquaresMethod
  PUBLIC :: LagrangeInterpolation
  PUBLIC :: DividedDifference, NewtonDividedDifference
  PUBLIC :: BSpline1, BSpline2, BSpline2dBC
  PUBLIC :: NormalizedBSpline, NormalizedBSplinedBC
  PUBLIC :: NSpline
  PUBLIC :: NaturalSpline
  PUBLIC :: SplineSecondInterpolationSet, SplineSecondInterpolationCalc
  PUBLIC :: SplineThirdInterpolationSet, SplineThirdInterpolationCalc
  PUBLIC :: BezierCurve, BSplineCurve
  PUBLIC :: BezierSurface, BSplineSurface
  ! 共有サブルーチン(自作) +++++++++++++++++++++++++++++++++++++++
  PUBLIC :: LinearInterpolation
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** Gauss の消去法                                                ****
!***********************************************************************
SUBROUTINE Gauss(DIMS, DIME, A, B, X)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: DIMS, DIME
  REAL,    INTENT(IN)  :: A(DIMS:DIME, DIMS:DIME), B(DIMS:DIME)
  REAL,    INTENT(OUT) :: X(DIMS:DIME)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL :: ML(DIMS:DIME, DIMS:DIME), MR(DIMS:DIME)
  INTEGER :: N, M, NP, MP, NN
  REAL :: P, SWAP
  REAL, PARAMETER :: ZeroLimit = 1.0E-10
  ! 処理開始 ***********************************************************
  ! 初期化 -------------------------------------------------------------
  DO M = DIMS, DIME
    DO N = DIMS, DIME
      ML(N, M) = A(N, M)
    ENDDO
  ENDDO
  DO N = DIMS, DIME
    MR(N) = B(N)
  ENDDO
  ! 前進消去 -----------------------------------------------------------
  DO N = DIMS, DIME
    ! ピボット選択
    P = ML(N, N)
    DO NP = N + 1, DIME
      IF( ABS(P) .LT. ABS(ML(NP, N)) ) THEN
        P = ML(NP, N)
        DO MP = DIMS, DIME
          SWAP       = ML(N,  MP)
          ML(N,  MP) = ML(NP, MP)
          ML(NP, MP) = SWAP
        ENDDO
          SWAP   = MR(N)
          MR(N ) = MR(NP)
          MR(NP) = SWAP
      ENDIF
    ENDDO
    IF(ABS(P) .LE. ZeroLimit) THEN
      WRITE(*, '(A)') 'Cant solve a simulaneous equation...'
      RETURN
    ENDIF
    ! ピボット行の標準化
    DO M = DIMS, DIME
      ML(N, M) = ML(N, M) / P
    ENDDO
    MR(N) = MR(N) / P
    ! ピボット行以下の掃き出し
    DO NN = N + 1, DIME
      IF(NN .EQ. N) CYCLE
      P = ML(NN, N)
      IF(ABS(P) .LE. ZeroLimit) CYCLE
      DO M = DIMS, DIME
        ML(NN, M) = ML(NN, M) - ML(N, M) * P
      ENDDO
      MR(NN) = MR(NN) - MR(N) * P
    ENDDO
  ENDDO
  ! 後進代入 -----------------------------------------------------------
  DO N = DIME, DIMS, -1
    DO M = N + 1, DIME
      MR(N) = MR(N) - ML(N, M) * MR(M)
    ENDDO
    MR(N) = MR(N) / ML(N, N)
    X(N)  = MR(N)
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE Gauss
!***********************************************************************
!**** Newton-Raphson Method (Secant Method)                         ****
!***********************************************************************
SUBROUTINE NRM_SM(Func, NMAX, EPS, x_init, x_solv)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    EXTERNAL    :: Func         ! 関数 Func(x) = 0
  INTEGER, INTENT(IN)  :: NMAX         ! 計算回数の最大値
  REAL,    INTENT(IN)  :: EPS          ! 収束判定値
  REAL,    INTENT(IN)  :: x_init       ! 初期値
  REAL,    INTENT(OUT) :: x_solv       ! Func(x) = 0 を満たす x の収束解
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: n
  REAL(8) :: x0, x1, x, dx, df
  ! 処理開始 ***********************************************************
  ! WRITE(*,'(A)') '------------ Newton-Raphson Method log ------------'
  x0 = DBLE(x_init)
  x  = DBLE(x_init) + 0.01D+0 * DBLE(x_init)
  DO n = 1, NMAX
    df = DBLE(Func(REAL(x)) - Func(REAL(x0)))
    IF(ABS(df) .LE. DBLE(EPS)) THEN
      EXIT
    ENDIF
    x1 = (x0 * DBLE(Func(REAL(x))) - x * DBLE(Func(REAL(x0)))) / df
    dx = ABS((x1 - x) / x1)
    x0 = x
    x  = x1
    ! WRITE(*,'(A, I4, 2X, A, 2X, A, E16.8E3, 2X, A, 2X, A, E16.8E3)') &
    ! &  'Iteration = ', n, '|', 'x = ', x, '|', 'dx = ', dx
    IF(dx .LE. DBLE(EPS)) THEN
      EXIT
    ENDIF
  ENDDO
  x_solv = REAL(x)
  ! WRITE(*,'(A)') '---------------------------------------------------'
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NRM_SM
!***********************************************************************
!**** 最小2乗法(Least Squares Method)による関数作成                 ****
!***********************************************************************
SUBROUTINE LeastSquaresMethod( &
&            FuncDIM, DATAS, DATAE, X, F, Weight, Coefficient &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: FuncDIM, DATAS, DATAE
  REAL,    INTENT(IN)  :: X(DATAS:DATAE), F(DATAS:DATAE), &
  &                       Weight(DATAS:DATAE)
  REAL,    INTENT(OUT) :: Coefficient(0:FuncDIM)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: L, N, M
  REAL :: MatrixLeft(0:FuncDIM, 0:FuncDIM), MatrixRight(0:FuncDIM)
  ! 処理開始 ***********************************************************
  ! 左辺、右辺配列初期化 +++++++++++++++++++++++++++++++++++++++++++++++
  DO M = 0, FuncDIM
    DO N = 0, FuncDIM
      MatrixLeft(N, M) = 0.0
    ENDDO
  ENDDO
  DO N = 0, FuncDIM
    MatrixRight(N) = 0.0
  ENDDO
  ! 左辺、右辺配列初期値 +++++++++++++++++++++++++++++++++++++++++++++++
  DO M = 0, FuncDIM
    DO N = 0, FuncDIM
      DO L = DATAS, DATAE
        MatrixLeft(N, M) = MatrixLeft(N, M) &
        &                + Weight(L) * PHI(X(L), N) * PHI(X(L), M)
      ENDDO
    ENDDO
  ENDDO
  DO N=0,FuncDIM
    DO L=DATAS,DATAE
      MatrixRight(N) = MatrixRight(N) &
      &              + Weight(L) * F(L) * PHI(X(L), N)
    ENDDO
  ENDDO
  ! ガウス法の消去法 +++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL Gauss(0, FuncDIM, MatrixLeft, MatrixRight, Coefficient)
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 関数定義                                                      ****
!***********************************************************************
FUNCTION PHI(X, J)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN) :: X
  INTEGER, INTENT(IN) :: J
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL PHI
  ! 処理開始 ***********************************************************
  IF( J .EQ. 0 ) THEN
    PHI = 1.0
  ELSE
    PHI = X**J
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
END FUNCTION PHI
! 定義終了 *************************************************************
END SUBROUTINE LeastSquaresMethod
!***********************************************************************
!**** Lagrangeの補間多項式                                          ****
!**** (節点数が多すぎるとルンゲの現象が生じるので注意)              ****
!***********************************************************************
SUBROUTINE LagrangeInterpolation( &
&            n, xn, yn, x, L &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: n             ! 節点数
  REAL,    INTENT(IN)  :: xn(n), yn(n)  ! 節点座標
  REAL,    INTENT(IN)  :: x             ! 補間点の x 座標
  REAL,    INTENT(OUT) :: L             ! 補間点の y 座標
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i, j
  REAL    :: P
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(n .LE. 1) THEN
    WRITE(*,'(A)') 'LagrangeInterpolation - Error : n <= 1'
    STOP
  ENDIF
  DO i = 1, n - 1
    IF(xn(i) .GE. xn(i + 1)) THEN
      WRITE(*,'(A)') 'LagrangeInterpolation - Error : xn(i) >= xn(i+1)'
      STOP
    ENDIF
  ENDDO
  ! 公式の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  L = 0.0
  DO i = 1, n
    P = 1.0
    DO j = 1, i - 1
      P = P * (x - xn(j)) / (xn(i) - xn(j))
    ENDDO
    DO j = i + 1, n
      P = P * (x - xn(j)) / (xn(i) - xn(j))
    ENDDO
    L = L + P * yn(i)
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE LagrangeInterpolation
!***********************************************************************
!**** 差分商                                                        ****
!***********************************************************************
REAL FUNCTION DividedDifference(n, xn, yn) RESULT(F)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: n             ! データ点数 (n - 1階差分商)
  REAL,    INTENT(IN)  :: xn(n), yn(n)  ! 関数 y = f(x)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(n .LE. 1) THEN
    WRITE(*,'(A)') 'DividedDifference - Error : n <= 1'
    STOP
  ENDIF
  DO i = 1, n - 1
    IF(xn(i) .GE. xn(i + 1)) THEN
      WRITE(*,'(A)') 'DividedDifference - Error : xn(i) >= xn(i+1)'
      STOP
    ENDIF
  ENDDO
  ! 公式の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  F = 0.0
  DO i = 1, n
    F = F + yn(i) / P(i, xn(i))
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 関数定義                                                      ****
!***********************************************************************
REAL FUNCTION P(j, t)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN) :: j
  REAL,    INTENT(IN) :: t
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i
  ! 処理開始 ***********************************************************
  P = 1.0
  DO i = 1, j - 1
    P = P * (t - xn(i))
  ENDDO
  DO i = j + 1, n
    P = P * (t - xn(i))
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END FUNCTION P
! 定義終了 *************************************************************
END FUNCTION DividedDifference
!***********************************************************************
!**** Newtonの差分商補間公式                                        ****
!***********************************************************************
REAL FUNCTION NewtonDividedDifference(n, xn, yn, x) RESULT(f)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: n             ! データ点数
  REAL,    INTENT(IN)  :: xn(n), yn(n)  ! 関数 y = f(x)
  REAL,    INTENT(IN)  :: x             ! 補間点の x
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: k, i
  REAL    :: Fk
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(n .LE. 1) THEN
    WRITE(*,'(A)') 'NewtonDividedDifference - Error : n <= 1'
    STOP
  ENDIF
  DO i = 1, n - 1
    IF(xn(i) .GE. xn(i + 1)) THEN
      WRITE(*,'(A)') 'NewtonDividedDifference - ' // &
      &              'Error : xn(i) >= xn(i+1)'
      STOP
    ENDIF
  ENDDO
  ! 公式の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  f = yn(1)
  DO k = 1, n - 1
    Fk = DividedDifference(k + 1, xn(1: k + 1), yn(1: k + 1))
    DO i = 1, k
      Fk = Fk * (x - xn(i))
    ENDDO
    f = f + Fk
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END FUNCTION NewtonDividedDifference
!***********************************************************************
!**** B-スプライン(Schoenbergの定義)                                ****
!***********************************************************************
REAL FUNCTION BSpline1(m, xm, x) RESULT(Mi)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: m             ! 位数 (m-1次)
  REAL,    INTENT(IN)  :: xm(m + 1)     ! 節点
  REAL,    INTENT(IN)  :: x             ! 補間点の x
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i
  REAL    :: ym(m + 1)
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(m .LE. 1) THEN
    WRITE(*,'(A)') 'BSpline1 - Error : m <= 1'
    STOP
  ENDIF
  DO i = 1, m
    IF(xm(i) .GE. xm(i + 1)) THEN
      WRITE(*,'(A)') 'BSpline1 - Error : xm(i) >= xm(i+1)'
      STOP
    ENDIF
  ENDDO
  ! 公式の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO i = 1, m + 1
    ym(i) = MM(x, xm(i))
  ENDDO
  Mi = DividedDifference(m + 1, xm, ym)
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 関数定義                                                      ****
!***********************************************************************
REAL FUNCTION MM(x, t)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN) :: x, t
  ! 処理開始 ***********************************************************
  IF(t - x .GT. 0.0) THEN
    MM = REAL(m) * (t - x)**(m - 1)
  ELSE
    MM = 0.0
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
END FUNCTION MM
! 定義終了 *************************************************************
END FUNCTION BSpline1
!***********************************************************************
!**** B-スプライン(de Boorの定義)                                   ****
!***********************************************************************
REAL FUNCTION BSpline2(k, xi_k, x) RESULT(M_ik)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: k             ! 位数 (k-1次)
  REAL,    INTENT(IN)  :: xi_k(k + 1)   ! 節点
  REAL,    INTENT(IN)  :: x             ! 補間点の x
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i
  REAL    :: M_k(k + 1)
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(k .LE. 1) THEN
    WRITE(*,'(A)') 'BSpline2 - Error : k <= 1'
    STOP
  ENDIF
  DO i = 1, k
    IF(xi_k(i) .GE. xi_k(i + 1)) THEN
      WRITE(*,'(A)') 'BSpline2 - Error : xi_k(i) >= xi_k(i+1)'
      STOP
    ENDIF
  ENDDO
  ! 公式の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO i = 1, k + 1
    M_k(i) = MM(x, xi_k(i))
  ENDDO
  M_ik = DividedDifference(k + 1, xi_k, M_k)
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 関数定義                                                      ****
!***********************************************************************
REAL FUNCTION MM(x, t)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN) :: x, t
  ! 処理開始 ***********************************************************
  IF(t - x .GT. 0.0) THEN
    MM = (t - x)**(k - 1)
  ELSE
    MM = 0.0
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
END FUNCTION MM
! 定義終了 *************************************************************
END FUNCTION BSpline2
!***********************************************************************
!**** Bスプライン(de Boorの定義, de Boor-Coxの算法)                 ****
!***********************************************************************
REAL FUNCTION BSpline2dBC(k, xi_k, x) RESULT(M_k)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: k             ! 位数 (k-1次)
  REAL,    INTENT(IN)  :: xi_k(k + 1)   ! 節点
  REAL,    INTENT(IN)  :: x             ! 補間点の x
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i, j
  REAL    :: M_ik(k, k), M_ik1, M_ik2
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(k .LE. 1) THEN
    WRITE(*,'(A)') 'BSpline2dBC - Error : k <= 1'
    STOP
  ENDIF
  DO i = 1, k
    IF(xi_k(i) .GT. xi_k(i + 1)) THEN
      WRITE(*,'(A)') 'BSpline2dBC - Error : xi_k(i) >= xi_k(i+1)'
      STOP
    ENDIF
  ENDDO
  ! 0次BスプラインM_i,1を計算 ++++++++++++++++++++++++++++++++++++++++++
  ! 多重節点の場合、計算は破綻しないが不正確な値しか得られない
  DO i = 1, k
    IF(xi_k(i) .LE. x .AND. x .LT. xi_k(i + 1)) THEN
      M_ik(i, 1) = 1.0 / (xi_k(i + 1) - xi_k(i))
    ELSE
      M_ik(i, 1) = 0.0
    ENDIF
  ENDDO
  ! k-1次Bスプラインを計算 +++++++++++++++++++++++++++++++++++++++++++++
  DO j = 2, k
    DO i = 1, k - j + 1
      IF(xi_k(i + j) - xi_k(i + 1) .NE. 0.0) THEN
        M_ik1 = (xi_k(i + j) - x) / (xi_k(i + j) - xi_k(i)) &
        &     * M_ik(i + 1, j - 1)
      ELSE
        M_ik1 = 0.0
      ENDIF
      IF(xi_k(i + j - 1) - xi_k(i) .NE. 0.0) THEN
        M_ik2 = (xi_k(i) - x) / (xi_k(i + j) - xi_k(i)) &
        &     * M_ik(i, j - 1)
      ELSE
        M_ik2 = 0.0
      ENDIF
      M_ik(i, j) = M_ik1 - M_ik2
    ENDDO
  ENDDO
  M_k = M_ik(1, k)
  ! 処理終了 ***********************************************************
  RETURN
END FUNCTION BSpline2dBC
!***********************************************************************
!**** 正規化されたB-スプライン                                      ****
!***********************************************************************
REAL FUNCTION NormalizedBSpline(k, xi_k, x) RESULT(N_ik)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: k             ! 位数 (k-1次)
  REAL,    INTENT(IN)  :: xi_k(k + 1)   ! 節点
  REAL,    INTENT(IN)  :: x             ! 補間点の x
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i
  REAL    :: M_ik, M_k(k + 1)
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(k .LE. 1) THEN
    WRITE(*,'(A)') 'NormalizedBSpline - Error : k <= 1'
    STOP
  ENDIF
  DO i = 1, k
    IF(xi_k(i) .GE. xi_k(i + 1)) THEN
      WRITE(*,'(A)') 'NormalizedBSpline - Error : xi_k(i) >= xi_k(i+1)'
      STOP
    ENDIF
  ENDDO
  ! Bスプライン(de Boorの定義)を計算 +++++++++++++++++++++++++++++++++++
  DO i = 1, k + 1
    M_k(i) = MM(x, xi_k(i))
  ENDDO
  M_ik = DividedDifference(k + 1, xi_k, M_k)
  ! 正規化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  N_ik = (xi_k(k + 1) - xi_k(1)) * M_ik
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 関数定義                                                      ****
!***********************************************************************
REAL FUNCTION MM(x, t)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN) :: x, t
  ! 処理開始 ***********************************************************
  IF(t - x .GT. 0.0) THEN
    MM = (t - x)**(k - 1)
  ELSE
    MM = 0.0
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
END FUNCTION MM
! 定義終了 *************************************************************
END FUNCTION NormalizedBSpline
!***********************************************************************
!**** 正規化されたBスプライン(de Boor-Coxの算法)                    ****
!***********************************************************************
REAL FUNCTION NormalizedBSplinedBC(k, xi_k, x) RESULT(N_k)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: k             ! 位数 (k-1次)
  REAL,    INTENT(IN)  :: xi_k(k + 1)   ! 節点
  REAL,    INTENT(IN)  :: x             ! 補間点の x
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i, j
  REAL    :: N_ik(k, k), N_ik1, N_ik2
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(k .LE. 1) THEN
    WRITE(*,'(A)') 'NormalizedBSplinedBC - Error : k <= 1'
    STOP
  ENDIF
  DO i = 1, k
    IF(xi_k(i) .GT. xi_k(i + 1)) THEN
      WRITE(*,'(A)') 'NormalizedBSplinedBC - ' // &
      &              'Error : xi_k(i) >= xi_k(i+1)'
      STOP
    ENDIF
  ENDDO
  ! 正規化された0次BスプラインN_i,1を計算 ++++++++++++++++++++++++++++++
  DO i = 1, k
    IF(xi_k(i) .LE. x .AND. x .LT. xi_k(i + 1)) THEN
      N_ik(i, 1) = 1.0
    ELSE
      N_ik(i, 1) = 0.0
    ENDIF
  ENDDO
  DO i = k, 1, - 1
    IF(ABS(xi_k(i) - xi_k(i + 1)) .LT. 1.0E-8 .AND. &
    &  ABS(xi_k(i) - x) .LT. 1.0E-8) THEN
      N_ik(i, 1) = 1.0
    ELSE
      EXIT
    ENDIF
  ENDDO
  ! 正規化されたk-1次Bスプラインを計算 +++++++++++++++++++++++++++++++++
  DO j = 2, k
    DO i = 1, k - j + 1
      IF(xi_k(i + j) - xi_k(i + 1) .NE. 0.0) THEN
        N_ik1 = (xi_k(i + j) - x) / (xi_k(i + j) - xi_k(i + 1)) &
        &     * N_ik(i + 1, j - 1)
      ELSE
        N_ik1 = N_ik(i + 1, j - 1)
      ENDIF
      IF(xi_k(i + j - 1) - xi_k(i) .NE. 0.0) THEN
        N_ik2 = (xi_k(i) - x) / (xi_k(i + j - 1) - xi_k(i)) &
        &     * N_ik(i, j - 1)
      ELSE
        N_ik2 = N_ik(i, j - 1)
      ENDIF
      N_ik(i, j) = N_ik1 - N_ik2
    ENDDO
  ENDDO
  N_k = N_ik(1, k)
  ! 処理終了 ***********************************************************
  RETURN
END FUNCTION NormalizedBSplinedBC
!***********************************************************************
!**** N-スプライン                                                  ****
!***********************************************************************
REAL FUNCTION NSpline(k, xk, x) RESULT(Ni)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: k             ! 位数 (2k-1次の自然スプライン)
  REAL,    INTENT(IN)  :: xk(k + 1)     ! 節点
  REAL,    INTENT(IN)  :: x             ! 補間点の x
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i
  REAL    :: yk(k + 1)
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(k .LE. 1) THEN
    WRITE(*,'(A)') 'NSpline - Error : k <= 1'
    STOP
  ENDIF
  DO i = 1, k
    IF(xk(i) .GE. xk(i + 1)) THEN
      WRITE(*,'(A)') 'NSpline - Error : xk(i) >= xk(i+1)'
      STOP
    ENDIF
  ENDDO
  ! 公式の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO i = 1, k + 1
    yk(i) = N(x, xk(i))
  ENDDO
  Ni = DividedDifference(k + 1, xk, yk)
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 関数定義                                                      ****
!***********************************************************************
REAL FUNCTION N(x, t)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN) :: x, t
  ! 処理開始 ***********************************************************
  IF(x - t .GT. 0.0) THEN
    N = (x - t)**(2 * k - 1)
  ELSE
    N = 0.0
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
END FUNCTION N
! 定義終了 *************************************************************
END FUNCTION NSpline
!***********************************************************************
!**** 自然スプライン                                                ****
!***********************************************************************
SUBROUTINE NaturalSpline(n, xn, yn, k, x, y)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: n             ! 節点数
  REAL,    INTENT(IN)  :: xn(n), yn(n)  ! 節点座標
  INTEGER, INTENT(IN)  :: k             ! 位数 (2k-1次)
  REAL,    INTENT(IN)  :: x             ! 補間値に対応する x の値
  REAL,    INTENT(OUT) :: y             ! 補間値 y = s(x)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i, j
  REAL    :: Ni1(n, n), Ni2(n, n), Nij(1: n - k, 1: n - k)
  REAL    :: delta(1: n - k)
  REAL    :: b(1: n - k)
  REAL    :: Nsp(1: n - k)
  REAL    :: P, Pk(1: k)
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(n .LE. 1) THEN
    WRITE(*,'(A)') 'NaturalSpline - Error : n <= 1'
    STOP
  ENDIF
  DO i = 1, n - 1
    IF(xn(i) .GE. xn(i + 1)) THEN
      WRITE(*,'(A)') 'NaturalSpline - Error : xn(i) >= xn(i+1)'
      STOP
    ENDIF
  ENDDO
  ! スプライン係数の計算 +++++++++++++++++++++++++++++++++++++++++++++++
  DO j = 1, n
  DO i = 1, n
    Ni1(i, j) = NN(xn(j), xn(i))
  ENDDO
  ENDDO
  ! DO j = 1, n
  !   WRITE(*,*) (Ni1(i, j), i = 1, n)
  ! ENDDO
  DO j = 1, n
  DO i = 1, n - k
    Ni2(i, j) = DividedDifference(k + 1, xn(i: i + k), Ni1(i: i + k, j))
  ENDDO
  ENDDO
  ! DO j = 1, n
  !   WRITE(*,*) (Ni2(i, j), i = 1, n - k)
  ! ENDDO
  DO j = 1, n - k
  DO i = 1, n - k
    Nij(i, j) = DividedDifference(k + 1, xn(j: j + k), Ni2(i, j: j + k))
  ENDDO
  ENDDO
  ! DO j = 1, n - k
  !   WRITE(*,*) (Nij(i, j), i = 1, n - k)
  ! ENDDO
  DO j = 1, n - k
    delta(j) = DividedDifference(k + 1, xn(j: j + k), yn(j: j + k))
  ENDDO
  ! DO j = 1, n - k
  !   WRITE(*,*) delta(j)
  ! ENDDO
  CALL Gauss(1, n - k, Nij, delta, b)
  ! DO i = 1, n - k
  !   WRITE(*,*) b(i)
  ! ENDDO
  ! N-スプラインの計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  DO i = 1, n - k
    Nsp(i) = NSpline(k, xn(i: i + k), x)
  ENDDO
  ! DO i = 1, n - k
  !   WRITE(*,*) Nsp(i)
  ! ENDDO
  ! k-1次多項式の計算 ++++++++++++++++++++++++++++++++++++++++++++++++++
  Pk(1) = yn(1)
  DO j = 2, k
    Pk(j) = DividedDifference(j, xn(1: j), yn(1: j))
    DO i = 1, n - k
      Pk(j) = Pk(j) - b(i) &
      &     * DividedDifference(k, xn(1: k), Ni2(i, 1: k))
    ENDDO
  ENDDO
  P = Pk(1)
  DO j = 2, k
    P = P + (x - xn(j - 1)) * Pk(j)
  ENDDO
  ! WRITE(*,*) P
  ! 自然スプラインの値を計算 +++++++++++++++++++++++++++++++++++++++++++
  y = P
  DO i = 1, n - k
    y = y + b(i) * NSp(i)
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 関数定義                                                      ****
!***********************************************************************
REAL FUNCTION NN(x, t)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN) :: x, t
  ! 処理開始 ***********************************************************
  IF(x - t .GT. 0.0) THEN
    NN = (x - t)**(2 * k - 1)
  ELSE
    NN = 0.0
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
END FUNCTION NN
! 定義終了 *************************************************************
END SUBROUTINE NaturalSpline
!***********************************************************************
!**** スプライン二次補間関数の導出                                  ****
!***********************************************************************
SUBROUTINE SplineSecondInterpolationSet( &
&            IS, IE, X, Y, SP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE
  REAL,    INTENT(IN)  :: X(IS: IE), Y(IS: IE)
  REAL,    INTENT(OUT) :: SP(1: 3, IS: IE - 1)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL    :: DX(IS: IE - 1)
  INTEGER :: I
  ! 処理開始 ***********************************************************
  ! Xの変化量を求める(DXにおけるIはXにおけるI+1/2) ---------------------
  DO I = IS, IE - 1
    DX(I) =-X(I) + X(I+1)
  ENDDO
  ! スプライン関数を求める ---------------------------------------------
  ! SP1
  DO I = IS, IE - 1
    SP(1, I) = Y(I)
  ENDDO
  ! SP2
  DO I = IS + 1, IE - 1
    SP(2, I) =-DX(I) / (DX(I-1) * (DX(I-1) + DX(I))) * Y(I-1) &
    &        - (DX(I-1) - DX(I)) / (DX(I-1) * DX(I)) * Y(I  ) &
    &        + DX(I-1) / (DX(I) * (DX(I-1) + DX(I))) * Y(I+1)
  ENDDO
  I = IS
    SP(2, I) = (Y(I+1) - Y(I)) / DX(I)
    ! 理論的には正しいが DX(I) = 0.5 のときゼロ割が発生する
    ! SP(2, I) = (2.0 * (Y(I+1) - Y(I)) - SP(2, I+1)) &
    ! &        / (2.0 * DX(I) - 1.0)
  ! SP3
  DO I = IS + 1, IE - 2
    SP(3, I) = (SP(2, I+1) - SP(2, I)) / (2.0 * DX(I))
  ENDDO
  I = IS
    SP(3, I) = 0.0
    ! 理論的には正しいが DX(I) = 0.5 のとき使えない
    ! SP(3, I) = (SP(2, I+1) - SP(2, I)) / (2.0 * DX(I))
  I = IE - 1
    SP(3, I) = (-Y(I) + Y(I+1)) / DX(I-1)**2 &
    &        - SP(2, I) / DX(I-1)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SplineSecondInterpolationSet
!***********************************************************************
!**** スプライン二次補間による補間値の計算                          ****
!***********************************************************************
SUBROUTINE SplineSecondInterpolationCalc( &
&            IS, IE, X, Y, SP, X0, Y0 &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE
  REAL,    INTENT(IN)  :: X(IS: IE), Y(IS: IE)
  REAL,    INTENT(IN)  :: SP(1: 3, IS: IE - 1)
  REAL,    INTENT(IN)  :: X0
  REAL,    INTENT(OUT) :: Y0
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I
  ! 処理開始 ***********************************************************
  Y0 = Y(IE)
  DO I = IS, IE - 1
    IF(X0 .GE. X(I) .AND. X0 .LT. X(I + 1)) THEN
      Y0 = SP(1, I) + (X0 - X(I)) * ( &
      &    SP(2, I) + (X0 - X(I)) * SP(3, I))
      EXIT
    ENDIF
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SplineSecondInterpolationCalc
!***********************************************************************
!**** スプライン三次補間関数の導出                                  ****
!***********************************************************************
SUBROUTINE SplineThirdInterpolationSet( &
&            IS, IE, X, Y, SP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! INTEGER, PARAMETER :: NMAX = 100
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE
  REAL,    INTENT(IN)  :: X(IS: IE), Y(IS: IE)
  REAL,    INTENT(OUT) :: SP(1: 4, IS: IE - 1)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL    :: DX(IS: IE), DY(IS: IE)
  REAL    :: A(IS: IE, IS: IE), B(IS: IE)
  INTEGER :: I, J
  ! INTEGER :: N
  ! 処理開始 ***********************************************************
  ! 連立方程式の各項を計算 ---------------------------------------------
  DO I = IS + 1, IE
    DX(I) = X(I) - X(I - 1)
  ENDDO
  I = IS
    DX(I) = DX(I + 1)
  DO I = IS, IE
    DO J = IS, IE
      A(J, I) = 0.0
    ENDDO
    B(I) = 0.0
  ENDDO
  DO I = IS + 1, IE - 1
    A(I, I - 1) = DX(I + 1)
    A(I, I    ) = 2.0 * (DX(I) + DX(I - 1))
    A(I, I + 1) = DX(I)
    B(I)        = 3.0 * ( &
    &               (Y(I) - Y(I - 1)) / DX(I) * DX(I + 1) &
    &             + (Y(I + 1) - Y(I)) / DX(I + 1) * DX(I) &
    &           )
  ENDDO
  I = IS
    A(I, I) = 1.0
    B(I)    = (Y(I + 1) - Y(I    )) / DX(I + 1)
  I = IE
    A(I, I) = 1.0
    B(I)    = (Y(I    ) - Y(I - 1)) / DX(I    )
  ! Gauss の消去法で連立方程式を解く -----------------------------------
  CALL Gauss(IS, IE, A, B, DY)
  ! 反復法で Y の勾配を再計算(Gauss の消去法の解を初期条件とする) ------
  ! DO N = 1, NMAX
  !   DO I = IS + 1, IE - 1
  !     DY(I) = ( &
  !     &         3.0 * ( &
  !     &           (Y(I  ) - Y(I-1)) * DX(I+1) / DX(I  ) &
  !     &         + (Y(I+1) - Y(I  )) * DX(I  ) / DX(I+1) &
  !     &       ) &
  !     &       - DX(I+1) * DY(I-1) &
  !     &       - DX(I  ) * DY(I+1) &
  !     &     ) * 0.5 / (DX(I) + DX(I-1))
  !   ENDDO
  !   DY(IS) = DY(IS + 1)
  !   DY(IE) = DY(IE - 1)
  ! ENDDO
  ! スプライン関数を求める ---------------------------------------------
  DO I = IS, IE - 1
    SP(1, I) = Y(I)
    SP(2, I) = DY(I)
    SP(4, I) = ( (DY(I + 1) + DY(I)) &
    &          - 2.0 * DBLE(Y(I + 1) - Y(I)) / DX(I + 1) &
    &        ) / DX(I + 1)**2
    SP(3, I) = ((Y(I + 1) - Y(I)) / DX(I + 1) - DY(I)) / DX(I + 1) &
    &          - SP(4, I) * DX(I + 1)
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SplineThirdInterpolationSet
!***********************************************************************
!**** スプライン三次補間による補間値の計算                          ****
!***********************************************************************
SUBROUTINE SplineThirdInterpolationCalc( &
&            IS, IE, X, Y, SP, X0, Y0 &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE
  REAL,    INTENT(IN)  :: X(IS: IE), Y(IS: IE)
  REAL,    INTENT(IN)  :: SP(1: 4, IS: IE - 1)
  REAL,    INTENT(IN)  :: X0
  REAL,    INTENT(OUT) :: Y0
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I
  ! 処理開始 ***********************************************************
  Y0 = Y(IE)
  DO I = IS, IE - 1
    IF(X0 .GE. X(I) .AND. X0 .LT. X(I + 1)) THEN
      Y0 = SP(1, I) + (X0 - X(I)) * ( &
      &    SP(2, I) + (X0 - X(I)) * ( &
      &    SP(3, I) + (X0 - X(I)) * SP(4, I)))
      EXIT
    ENDIF
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SplineThirdInterpolationCalc
!***********************************************************************
!**** Bezier曲線                                                    ****
!***********************************************************************
SUBROUTINE BezierCurve( &
&            n, Px, Py, t, x, Y &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: n                  ! 制御点数
  REAL,    INTENT(IN)  :: Px(0: n), Py(0: n) ! 制御点座標
  REAL,    INTENT(IN)  :: t                  ! 媒介変数 (0 <= t <= 1)
  REAL,    INTENT(OUT) :: x, y               ! Bezier曲線座標
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i, l
  INTEGER :: A, B, C
  REAL    :: J
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(n .LE. 0) THEN
    WRITE(*,'(A)') 'BezierCurve - Error : n <= 0'
    STOP
  ENDIF
  IF(t .LT. 0.0 .OR. 1.0 .LT. t) THEN
    WRITE(*,'(A)') 'BezierCurve - Error : t < 0 or 1 < t'
    STOP
  ENDIF
  ! Bezier 曲線の計算 ++++++++++++++++++++++++++++++++++++++++++++++++++
  x = 0.0
  y = 0.0
  DO i = 0, n
    A = 1
    DO l = 1, n
      A = A * l
    ENDDO
    B = 1
    DO l = 1, i
      B = B * l
    ENDDO
    C = 1
    DO l = 1, n - i
      C = C * l
    ENDDO
    J = REAL(A) / REAL(B * C) * t**i * (1.0 - t)**(n - i)
    x = x + Px(i) * J
    y = y + Py(i) * J
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE BezierCurve
!***********************************************************************
!**** B-スプライン曲線(正規化されたB-スプラインを基底として使用)    ****
!***********************************************************************
SUBROUTINE BSplineCurve(l, P_ix, P_iy, k, t, x, y)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: l                     ! 節点数
  REAL,    INTENT(IN)  :: P_ix(0: l), P_iy(0: l)! 節点座標
  INTEGER, INTENT(IN)  :: k                     ! 位数 (k-1次)
  REAL,    INTENT(IN)  :: t           ! 媒介変数 (0 <= t <= l - k + 2)
  REAL,    INTENT(OUT) :: x, y        ! B-スプライン曲線座標
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i
  REAL    :: xi_i(- k + 1: l + 1)
  REAL    :: N_ik
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(l .LE. 0) THEN
    WRITE(*,'(A)') 'BSplineCurve - Error : l <= 0'
    STOP
  ENDIF
  IF(k .GE. l + 2) THEN
    WRITE(*,'(A)') 'BSplineCurve - Error : k >= l + 2'
    STOP
  ENDIF
  IF(t .LT. 0.0 .OR. REAL(l - k + 2) .LT. t) THEN
    WRITE(*,'(A)') 'BSplineCurve - Error : t < 0 or l - k + 2 < t'
    STOP
  ENDIF
  ! 節点 xi_i の計算(内部節点+付加節点) ++++++++++++++++++++++++++++++++
  DO i = - k + 1, 0
    xi_i(i) = 0.0
  ENDDO
  DO i = 0, l - k + 2
    xi_i(i) = REAL(i)
  ENDDO
  DO i = l - k + 2, l + 1
    xi_i(i) = REAL(l - k + 2)
  ENDDO
  ! B-スプライン曲線の計算 +++++++++++++++++++++++++++++++++++++++++++++
  x = 0.0
  y = 0.0
  DO i = 0, l
    N_ik = NormalizedBSplinedBC(k, xi_i(i - k + 1: i + 1), t)
    x = x + P_ix(i) * N_ik
    y = y + P_iy(i) * N_ik
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE BSplineCurve
!***********************************************************************
!**** Bezier曲面                                                    ****
!***********************************************************************
SUBROUTINE BezierSurface( &
&            n, m, Px, Py, Pz, t, u, x, y, z &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: n, m               ! 制御点数
  REAL,    INTENT(IN)  :: Px(0: n, 0: m), &  ! 制御点座標 (x)
  &                       Py(0: n, 0: m), &  !            (y)
  &                       Pz(0: n, 0: m)     !            (z)
  REAL,    INTENT(IN)  :: t, &               ! 媒介変数 (0 <= t <= 1)
  &                       u                  !          (0 <= u <= 1)
  REAL,    INTENT(OUT) :: x, y, z            ! Bezier曲線座標 (x, y, z)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i, j, l
  INTEGER :: A, B, C
  REAL    :: Jn, Km
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(n .LE. 0) THEN
    WRITE(*,'(A)') 'BezierSurface - Error : n <= 0'
    STOP
  ENDIF
  IF(m .LE. 0) THEN
    WRITE(*,'(A)') 'BezierSurface - Error : m <= 0'
    STOP
  ENDIF
  IF(t .LT. 0.0 .OR. 1.0 .LT. t) THEN
    WRITE(*,'(A)') 'BezierSurface - Error : t < 0 or 1 < t'
    STOP
  ENDIF
  IF(u .LT. 0.0 .OR. 1.0 .LT. u) THEN
    WRITE(*,'(A)') 'BezierSurface - Error : u < 0 or 1 < u'
    STOP
  ENDIF
  ! Bezier 曲線の計算 ++++++++++++++++++++++++++++++++++++++++++++++++++
  x = 0.0
  y = 0.0
  z = 0.0
  DO j = 0, m
  DO i = 0, n
    ! J_n,i (t) --------------------------------------------------------
    A = 1
    DO l = 1, n
      A = A * l
    ENDDO
    B = 1
    DO l = 1, i
      B = B * l
    ENDDO
    C = 1
    DO l = 1, n - i
      C = C * l
    ENDDO
    Jn = REAL(A) / REAL(B * C) * t**i * (1.0 - t)**(n - i)
    ! K_m,j (u) --------------------------------------------------------
    A = 1
    DO l = 1, m
      A = A * l
    ENDDO
    B = 1
    DO l = 1, j
      B = B * l
    ENDDO
    C = 1
    DO l = 1, m - j
      C = C * l
    ENDDO
    Km = REAL(A) / REAL(B * C) * u**j * (1.0 - u)**(m - j)
    ! x, y, z ----------------------------------------------------------
    x = x + Px(i,j) * Jn * Km
    y = y + Py(i,j) * Jn * Km
    z = z + Pz(i,j) * Jn * Km
  ENDDO
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE BezierSurface
!***********************************************************************
!**** B-スプライン曲面(正規化されたB-スプラインを基底として使用)    ****
!***********************************************************************
SUBROUTINE BSplineSurface(n, m, Px, Py, Pz, k, l, t, u, x, y, z)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! 制御点 -------------------------------------------------------------
  INTEGER, INTENT(IN)  :: n, m              ! 制御点数
  REAL,    INTENT(IN)  :: Px(0: n, 0: m), & ! 座標 x
  &                       Py(0: n, 0: m), & ! 座標 y
  &                       Pz(0: n, 0: m)    ! 座標 z
  ! 位数 ---------------------------------------------------------------
  INTEGER, INTENT(IN)  :: k, l              ! k - 1次, l - 1次
  ! 媒介変数 -----------------------------------------------------------
  REAL,    INTENT(IN)  :: t, &              ! 0 <= t <= n - k + 2
  &                       u                 ! 0 <= u <= m - l + 2
  ! B-スプライン曲線座標 -----------------------------------------------
  REAL,    INTENT(OUT) :: x, y, z           ! x, y, z
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: i, j
  REAL    :: xi_i(- k + 1: n + 1), xi_j(- l + 1: m + 1)
  REAL    :: N_ik, N_jl
  ! 処理開始 ***********************************************************
  ! 例外処理 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF(n .LE. 0) THEN
    WRITE(*,'(A)') 'BSplineSurface - Error : n <= 0'
    STOP
  ENDIF
  IF(m .LE. 0) THEN
    WRITE(*,'(A)') 'BSplineSurface - Error : m <= 0'
    STOP
  ENDIF
  IF(k .GE. n + 2) THEN
    WRITE(*,'(A)') 'BSplineSurface - Error : k >= n + 2'
    STOP
  ENDIF
  IF(l .GE. m + 2) THEN
    WRITE(*,'(A)') 'BSplineSurface - Error : l >= m + 2'
    STOP
  ENDIF
  IF(t .LT. 0.0 .OR. REAL(n - k + 2) .LT. t) THEN
    WRITE(*,'(A)') 'BSplineSurface - Error : t < 0 or n - k + 2 < t'
    STOP
  ENDIF
  IF(u .LT. 0.0 .OR. REAL(m - l + 2) .LT. u) THEN
    WRITE(*,'(A)') 'BSplineSurface - Error : u < 0 or m - l + 2 < u'
    STOP
  ENDIF
  ! 節点 xi_i の計算(内部節点+付加節点) ++++++++++++++++++++++++++++++++
  DO i = - k + 1, 0
    xi_i(i) = 0.0
  ENDDO
  DO i = 0, n - k + 2
    xi_i(i) = REAL(i)
  ENDDO
  DO i = n - k + 2, n + 1
    xi_i(i) = REAL(n - k + 2)
  ENDDO
  ! 節点 xi_j の計算(内部節点+付加節点) ++++++++++++++++++++++++++++++++
  DO j = - l + 1, 0
    xi_j(j) = 0.0
  ENDDO
  DO j = 0, m - l + 2
    xi_j(j) = REAL(j)
  ENDDO
  DO j = m - l + 2, m + 1
    xi_j(j) = REAL(m - l + 2)
  ENDDO
  ! B-スプライン曲線の計算 +++++++++++++++++++++++++++++++++++++++++++++
  x = 0.0
  y = 0.0
  z = 0.0
  DO j = 0, m
  DO i = 0, n
    N_ik = NormalizedBSplinedBC(k, xi_i(i - k + 1: i + 1), t)
    N_jl = NormalizedBSplinedBC(l, xi_j(j - l + 1: j + 1), u)
    x = x + Px(i,j) * N_ik * N_jl
    y = y + Py(i,j) * N_ik * N_jl
    z = z + Pz(i,j) * N_ik * N_jl
  ENDDO
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE BSplineSurface
!***********************************************************************
!**** 線形補間                                                      ****
!***********************************************************************
SUBROUTINE LinearInterpolation(NS, NE, X, Y, XP, YP)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! 制御点 -------------------------------------------------------------
  INTEGER, INTENT(IN)  :: NS, NE
  REAL   , INTENT(IN)  :: X(NS:NE), Y(NS:NE)
  ! 補間点 -------------------------------------------------------------
  REAL   , INTENT(IN)  :: XP
  REAL   , INTENT(OUT) :: YP
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: N
  REAL    :: mm, nn
  ! 処理開始 ***********************************************************
  DO N = NS, NE
    IF( XP < X(N) ) EXIT
  ENDDO
  IF(N == NS) THEN
    YP = Y(NS)
   ELSE IF(N == NE+1) THEN
    YP = Y(NE)
   ELSE
    mm = XP - X(N-1)
    nn = X(N) - XP
    YP = (nn * Y(N-1) + mm * Y(N)) / (mm + nn)
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE LinearInterpolation
! 定義終了 *************************************************************
END MODULE Package_Equation
