!***********************************************************************
!***********************************************************************
!**** パッケージ型モジュール                                        ****
!**** 流れ場計算用サブルーチン群                                    ****
!****                         2011.09.05 PROGRAMED BY SUZUKI MASAYA ****
!***********************************************************************
!***********************************************************************
MODULE Package_Flow
  ! モジュール宣言 *****************************************************
  !$USE omp_lib
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  PRIVATE
  ! サブルーチン宣言 ***************************************************
  ! 共有サブルーチン +++++++++++++++++++++++++++++++++++++++++++++++++++
  PUBLIC :: SetFlux2D, SetFlux3D
  PUBLIC :: SetFlux2DKEM, SetFlux3DKEM
  PUBLIC :: SetFlux2DRSM, SetFlux3DRSM
  PUBLIC :: SetPhysics2D, SetPhysics3D
  PUBLIC :: SetPhysics2DKEM, SetPhysics3DKEM
  PUBLIC :: SetPhysics2DRSM, SetPhysics3DRSM
  PUBLIC :: ViscosityCoefficient2D, ViscosityCoefficient3D
  PUBLIC :: CalLocalDt2D, CalLocalDt3D
  PUBLIC :: CalDt2D, CalDt3D
  PUBLIC :: RungeKutta2D, RungeKutta3D
  PUBLIC :: LUADI2D, LUADI3D
  PUBLIC :: Convection2D, Convection3D
  PUBLIC :: RotationForce2D, RotationForce3D
  PUBLIC :: Viscosity2D, Viscosity3D
  PUBLIC :: Turbulence2DEvmStd, Turbulence2DEvmStdKL, &
  &         Turbulence3DEvmStd, Turbulence3DEvmStdKL
  PUBLIC :: Turbulence2DEvmLS, Turbulence3DEvmLS
  PUBLIC :: Turbulence3DEvmAKN
  PUBLIC :: Turbulence2DNlEvmkeCLS, Turbulence3DNlEvmkeCLS
  PUBLIC :: Turbulence2DRSMGL, Turbulence3DRSMGL
  PUBLIC :: Turbulence2DRSMSSG, Turbulence3DRSMSSG
  PUBLIC :: Turbulence2DRSMLS, Turbulence3DRSMLS
  PUBLIC :: Turbulence2DLESCSM, Turbulence3DLESCSM
  PUBLIC :: Limiter2DKEM, Limiter3DKEM
  PUBLIC :: Limiter2DRSM, Limiter3DRSM
  PUBLIC :: CalResidual2D, CalResidual3D
  PUBLIC :: YapCorrection2D, YapCorrection3D
  PUBLIC :: WallFunctionKEM1S, WallFunctionRSM1S, &
  &         WallFunctionKEM2S, WallFunctionRSM2S, CalUtauS
  PUBLIC :: WallFunctionKEM1R, WallFunctionRSM1R, &
  &         WallFunctionKEM2R, WallFunctionRSM2R, CalUtauR
  ! 共有サブルーチン (自作) +++++++++++++++++++++++++++++++++++++++++++
  PUBLIC :: SaveFlux3DKEM, SumDQH3D
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 物理量から保存量流束を計算                                    ****
!**** 計算対象 : 単相, 二次元, 圧縮性, 層流                         ****
!***********************************************************************
SUBROUTINE SetFlux2D( &
&            GAMMA, &
&            IS, IE, JS, JE, LS, LE, &
&            QH, AJA, RHO, U, V, P, &
&            YY &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)            :: GAMMA
  INTEGER, INTENT(IN)            :: IS, IE, JS, JE
  INTEGER, INTENT(IN)            :: LS, LE
  REAL,    INTENT(OUT)           :: QH (IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(IN)            :: AJA(IS:IE, JS:JE)
  REAL,    INTENT(IN)            :: RHO(IS:IE, JS:JE), &
  &                                 U  (IS:IE, JS:JE), &
  &                                 V  (IS:IE, JS:JE), &
  &                                 P  (IS:IE, JS:JE)
  REAL,    INTENT(IN),  OPTIONAL :: YY (IS:IE, JS:JE, 5: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 4) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 保存量流束から物理量を計算 +++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J) .GT. 0.0) THEN
    QH(I,J,1) = RHO(I,J) / AJA(I,J)
    QH(I,J,2) = QH(I,J,1) * U(I,J)
    QH(I,J,3) = QH(I,J,1) * V(I,J)
    QH(I,J,4) = P(I,J) / ((GAMMA - 1.0) * AJA(I,J)) &
    &         + QH(I,J,1) * 0.5 * (U(I,J)**2 + V(I,J)**2)
    IF(.NOT. PRESENT(YY)) CYCLE
    QH(I,J,5:LE) = QH(I,J,1) * YY(I,J,5:LE)
  ELSE
    QH(I,J,:) = 0.0
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SetFlux2D
!***********************************************************************
!**** 物理量から保存量流束を計算                                    ****
!**** 計算対象 : 単相, 三次元, 圧縮性, 層流                         ****
!***********************************************************************
SUBROUTINE SetFlux3D( &
&            GAMMA, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            QH, AJA, RHO, U, V, W, P, &
&            YY &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)            :: GAMMA
  INTEGER, INTENT(IN)            :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)            :: LS, LE
  REAL,    INTENT(OUT)           :: QH (IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(IN)            :: AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)            :: RHO(IS:IE, JS:JE, KS:KE), &
  &                                 U  (IS:IE, JS:JE, KS:KE), &
  &                                 V  (IS:IE, JS:JE, KS:KE), &
  &                                 W  (IS:IE, JS:JE, KS:KE), &
  &                                 P  (IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN),  OPTIONAL :: YY (IS:IE, JS:JE, KS:KE, 6: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 5) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 保存量流束から物理量を計算 +++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, K)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J,K) .GT. 0.0) THEN
    QH(I,J,K,1) = RHO(I,J,K) / AJA(I,J,K)
    QH(I,J,K,2) = QH(I,J,K,1) * U(I,J,K)
    QH(I,J,K,3) = QH(I,J,K,1) * V(I,J,K)
    QH(I,J,K,4) = QH(I,J,K,1) * W(I,J,K)
    QH(I,J,K,5) = P(I,J,K) / ((GAMMA - 1.0) * AJA(I,J,K)) &
    &           + QH(I,J,K,1) * 0.5 * ( &
    &               U(I,J,K)**2 + V(I,J,K)**2 + W(I,J,K)**2 &
    &           )
    IF(.NOT. PRESENT(YY)) CYCLE
    QH(I,J,K,6:LE) = QH(I,J,K,1) * YY(I,J,K,6:LE)
  ELSE
    QH(I,J,K,:) = 0.0
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SetFlux3D
!***********************************************************************
!**** 物理量から保存量流束を計算                                    ****
!**** 計算対象 : 単相, 二次元, 圧縮性, k-e                          ****
!***********************************************************************
SUBROUTINE SetFlux2DKEM( &
&            GAMMA, &
&            IS, IE, JS, JE, LS, LE, &
&            QH, AJA, RHO, U, V, P, AK, EPS, &
&            YY &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)            :: GAMMA
  INTEGER, INTENT(IN)            :: IS, IE, JS, JE
  INTEGER, INTENT(IN)            :: LS, LE
  REAL,    INTENT(OUT)           :: QH (IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(IN)            :: AJA(IS:IE, JS:JE)
  REAL,    INTENT(IN)            :: RHO(IS:IE, JS:JE), &
  &                                 U  (IS:IE, JS:JE), &
  &                                 V  (IS:IE, JS:JE), &
  &                                 P  (IS:IE, JS:JE), &
  &                                 AK (IS:IE, JS:JE), &
  &                                 EPS(IS:IE, JS:JE)
  REAL,    INTENT(IN),  OPTIONAL :: YY (IS:IE, JS:JE, 7: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 6) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 保存量流束から物理量を計算 +++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J) .GT. 0.0) THEN
    QH(I,J,1) = RHO(I,J) / AJA(I,J)
    QH(I,J,2) = QH(I,J,1) * U(I,J)
    QH(I,J,3) = QH(I,J,1) * V(I,J)
    QH(I,J,4) = P(I,J) / ((GAMMA - 1.0) * AJA(I,J)) &
    &         + QH(I,J,1) * 0.5 * (U(I,J)**2 + V(I,J)**2) &
    &         + QH(I,J,1) * AK(I,J)
    QH(I,J,5) = QH(I,J,1) * AK(I,J)
    QH(I,J,6) = QH(I,J,1) * EPS(I,J)
    IF(.NOT. PRESENT(YY)) CYCLE
    QH(I,J,7:LE) = QH(I,J,1) * YY(I,J,7:LE)
  ELSE
    QH(I,J,:) = 0.0
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SetFlux2DKEM
!***********************************************************************
!**** 物理量から保存量流束を計算                                    ****
!**** 計算対象 : 単相, 三次元, 圧縮性, k-e                          ****
!***********************************************************************
SUBROUTINE SetFlux3DKEM( &
&            GAMMA, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            QH, AJA, RHO, U, V, W, P, AK, EPS, &
&            YY &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)            :: GAMMA
  INTEGER, INTENT(IN)            :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)            :: LS, LE
  REAL,    INTENT(OUT)           :: QH (IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(IN)            :: AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)            :: RHO(IS:IE, JS:JE, KS:KE), &
  &                                 U  (IS:IE, JS:JE, KS:KE), &
  &                                 V  (IS:IE, JS:JE, KS:KE), &
  &                                 W  (IS:IE, JS:JE, KS:KE), &
  &                                 P  (IS:IE, JS:JE, KS:KE), &
  &                                 AK (IS:IE, JS:JE, KS:KE), &
  &                                 EPS(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN),  OPTIONAL :: YY (IS:IE, JS:JE, KS:KE, 8: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 7) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 保存量流束から物理量を計算 +++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, K)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J,K) .GT. 0.0) THEN
    QH(I,J,K,1) = RHO(I,J,K) / AJA(I,J,K)
    QH(I,J,K,2) = QH(I,J,K,1) * U(I,J,K)
    QH(I,J,K,3) = QH(I,J,K,1) * V(I,J,K)
    QH(I,J,K,4) = QH(I,J,K,1) * W(I,J,K)
    QH(I,J,K,5) = P(I,J,K) / ((GAMMA - 1.0) * AJA(I,J,K)) &
    &           + QH(I,J,K,1) * 0.5 * ( &
    &               U(I,J,K)**2 + V(I,J,K)**2 + W(I,J,K)**2 &
    &           ) &
    &           + QH(I,J,K,1) * AK(I,J,K)
    QH(I,J,K,6) = QH(I,J,K,1) * AK(I,J,K)
    QH(I,J,K,7) = QH(I,J,K,1) * EPS(I,J,K)
    IF(.NOT. PRESENT(YY)) CYCLE
    QH(I,J,K,8:LE) = QH(I,J,K,1) * YY(I,J,K,8:LE)
  ELSE
    QH(I,J,K,:) = 0.0
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SetFlux3DKEM
!***********************************************************************
!**** 物理量から保存量流束を計算                                    ****
!**** 計算対象 : 単相, 二次元, 圧縮性, RSM                          ****
!***********************************************************************
SUBROUTINE SetFlux2DRSM( &
&            GAMMA, &
&            IS, IE, JS, JE, LS, LE, &
&            QH, AJA, RHO, U, V, P, uu, vv, ww, uv, EPS, &
&            YY &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)            :: GAMMA
  INTEGER, INTENT(IN)            :: IS, IE, JS, JE
  INTEGER, INTENT(IN)            :: LS, LE
  REAL,    INTENT(OUT)           :: QH (IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(IN)            :: AJA(IS:IE, JS:JE)
  REAL,    INTENT(IN)            :: RHO(IS:IE, JS:JE), &
  &                                 U  (IS:IE, JS:JE), &
  &                                 V  (IS:IE, JS:JE), &
  &                                 P  (IS:IE, JS:JE), &
  &                                 uu (IS:IE, JS:JE), &
  &                                 vv (IS:IE, JS:JE), &
  &                                 ww (IS:IE, JS:JE), &
  &                                 uv (IS:IE, JS:JE), &
  &                                 EPS(IS:IE, JS:JE)
  REAL,    INTENT(IN),  OPTIONAL :: YY (IS:IE, JS:JE, 10: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 9) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 保存量流束から物理量を計算 +++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J) .GT. 0.0) THEN
    QH(I,J,1) = RHO(I,J) / AJA(I,J)
    QH(I,J,2) = QH(I,J,1) * U(I,J)
    QH(I,J,3) = QH(I,J,1) * V(I,J)
    QH(I,J,4) = P(I,J) / ((GAMMA - 1.0) * AJA(I,J)) &
    &         + QH(I,J,1) * 0.5 * ( &
    &             U(I,J)**2 + V(I,J)**2 &
    &           + uu(I,J) + vv(I,J) + ww(I,J) &
    &         )
    QH(I,J,5) = QH(I,J,1) * uu(I,J)
    QH(I,J,6) = QH(I,J,1) * vv(I,J)
    QH(I,J,7) = QH(I,J,1) * ww(I,J)
    QH(I,J,8) = QH(I,J,1) * uv(I,J)
    QH(I,J,9) = QH(I,J,1) * EPS(I,J)
    IF(.NOT. PRESENT(YY)) CYCLE
    QH(I,J,10:LE) = QH(I,J,1) * YY(I,J,10:LE)
  ELSE
    QH(I,J,:) = 0.0
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SetFlux2DRSM
!***********************************************************************
!**** 物理量から保存量流束を計算                                    ****
!**** 計算対象 : 単相, 三次元, 圧縮性, RSM                          ****
!***********************************************************************
SUBROUTINE SetFlux3DRSM( &
&            GAMMA, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            QH, AJA, RHO, U, V, W, P, uu, vv, ww, uv, vw, wu, EPS, &
&            YY &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)            :: GAMMA
  INTEGER, INTENT(IN)            :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)            :: LS, LE
  REAL,    INTENT(OUT)           :: QH (IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(IN)            :: AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)            :: RHO(IS:IE, JS:JE, KS:KE), &
  &                                 U  (IS:IE, JS:JE, KS:KE), &
  &                                 V  (IS:IE, JS:JE, KS:KE), &
  &                                 W  (IS:IE, JS:JE, KS:KE), &
  &                                 P  (IS:IE, JS:JE, KS:KE), &
  &                                 uu (IS:IE, JS:JE, KS:KE), &
  &                                 vv (IS:IE, JS:JE, KS:KE), &
  &                                 ww (IS:IE, JS:JE, KS:KE), &
  &                                 uv (IS:IE, JS:JE, KS:KE), &
  &                                 vw (IS:IE, JS:JE, KS:KE), &
  &                                 wu (IS:IE, JS:JE, KS:KE), &
  &                                 EPS(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN),  OPTIONAL :: YY (IS:IE, JS:JE, KS:KE, 13: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 12) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 保存量流束から物理量を計算 +++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, K)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J,K) .GT. 0.0) THEN
    QH(I,J,K, 1) = RHO(I,J,K) / AJA(I,J,K)
    QH(I,J,K, 2) = QH(I,J,K,1) * U(I,J,K)
    QH(I,J,K, 3) = QH(I,J,K,1) * V(I,J,K)
    QH(I,J,K, 4) = QH(I,J,K,1) * W(I,J,K)
    QH(I,J,K, 5) = P(I,J,K) / ((GAMMA - 1.0) * AJA(I,J,K)) &
    &            + QH(I,J,K,1) * 0.5 * ( &
    &                U(I,J,K)**2 + V(I,J,K)**2 + W(I,J,K)**2 &
    &              + uu(I,J,K) + vv(I,J,K) + ww(I,J,K) &
    &            )
    QH(I,J,K, 6) = QH(I,J,K,1) * uu(I,J,K)
    QH(I,J,K, 7) = QH(I,J,K,1) * vv(I,J,K)
    QH(I,J,K, 8) = QH(I,J,K,1) * ww(I,J,K)
    QH(I,J,K, 9) = QH(I,J,K,1) * uv(I,J,K)
    QH(I,J,K,10) = QH(I,J,K,1) * vw(I,J,K)
    QH(I,J,K,11) = QH(I,J,K,1) * wu(I,J,K)
    QH(I,J,K,12) = QH(I,J,K,1) * EPS(I,J,K)
    IF(.NOT. PRESENT(YY)) CYCLE
    QH(I,J,K,13:LE) = QH(I,J,K,1) * YY(I,J,K,13:LE)
  ELSE
    QH(I,J,K, :) = 0.0
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SetFlux3DRSM
!***********************************************************************
!**** 保存量流束から物理量を計算                                    ****
!**** 計算対象 : 単相, 二次元, 圧縮性, 層流                         ****
!***********************************************************************
SUBROUTINE SetPhysics2D( &
&            RG, GAMMA, &
&            IS, IE, JS, JE, LS, LE, &
&            QH, AJA, RHO, U, V, P, T, &
&            YY &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)            :: RG, GAMMA
  INTEGER, INTENT(IN)            :: IS, IE, JS, JE
  INTEGER, INTENT(IN)            :: LS, LE
  REAL,    INTENT(IN)            :: QH (IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(IN)            :: AJA(IS:IE, JS:JE)
  REAL,    INTENT(OUT)           :: RHO(IS:IE, JS:JE), &
  &                                 U  (IS:IE, JS:JE), &
  &                                 V  (IS:IE, JS:JE), &
  &                                 P  (IS:IE, JS:JE), &
  &                                 T  (IS:IE, JS:JE)
  REAL,    INTENT(OUT), OPTIONAL :: YY (IS:IE, JS:JE, 5: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 4) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 保存量流束から物理量を計算 +++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
  DO J = JS, JE
  DO I = IS, IE
  IF(QH(I,J,1) .GT. 0.0) THEN
    RHO(I,J) = AJA(I,J) * QH(I,J,1)
    U(I,J)   = QH(I,J,2) / QH(I,J,1)
    V(I,J)   = QH(I,J,3) / QH(I,J,1)
    P(I,J)   = ( QH(I,J,4) &
    &          - 0.5 * (QH(I,J,2)**2 + QH(I,J,3)**2) / QH(I,J,1) &
    &        ) * (GAMMA - 1.0) * AJA(I,J)
    T(I,J)   = P(I,J) / (RHO(I,J) * RG)
    IF(.NOT. PRESENT(YY)) CYCLE
    YY(I,J,5:LE) = YY(I,J,5:LE) / QH(I,J,1)
  ELSE
    RHO(I,J) = 0.0
    U(I,J)   = 0.0
    V(I,J)   = 0.0
    P(I,J)   = 0.0
    T(I,J)   = 0.0
    IF(.NOT. PRESENT(YY)) CYCLE
    YY(I,J,5:LE) = 0.0
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SetPhysics2D
!***********************************************************************
!**** 保存量流束から物理量を計算                                    ****
!**** 計算対象 : 単相, 三次元, 圧縮性, 層流                         ****
!***********************************************************************
SUBROUTINE SetPhysics3D( &
&            RG, GAMMA, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            QH, AJA, RHO, U, V, W, P, T, &
&            YY &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)            :: RG, GAMMA
  INTEGER, INTENT(IN)            :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)            :: LS, LE
  REAL,    INTENT(IN)            :: QH (IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(IN)            :: AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT)           :: RHO(IS:IE, JS:JE, KS:KE), &
  &                                 U  (IS:IE, JS:JE, KS:KE), &
  &                                 V  (IS:IE, JS:JE, KS:KE), &
  &                                 W  (IS:IE, JS:JE, KS:KE), &
  &                                 P  (IS:IE, JS:JE, KS:KE), &
  &                                 T  (IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT), OPTIONAL :: YY (IS:IE, JS:JE, KS:KE, 6: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 5) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 保存量流束から物理量を計算 +++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, K)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(QH(I,J,K,1) .GT. 0.0) THEN
    RHO(I,J,K) = AJA(I,J,K) * QH(I,J,K,1)
    U(I,J,K)   = QH(I,J,K,2) / QH(I,J,K,1)
    V(I,J,K)   = QH(I,J,K,3) / QH(I,J,K,1)
    W(I,J,K)   = QH(I,J,K,4) / QH(I,J,K,1)
    P(I,J,K)   = ( QH(I,J,K,5) &
    &            - 0.5 * ( QH(I,J,K,2)**2 &
    &                    + QH(I,J,K,3)**2 &
    &                    + QH(I,J,K,4)**2 ) / QH(I,J,K,1) &
    &          ) * (GAMMA - 1.0) * AJA(I,J,K)
    T(I,J,K)   = P(I,J,K) / (RHO(I,J,K) * RG)
    IF(.NOT. PRESENT(YY)) CYCLE
    YY(I,J,K,6:LE) = QH(I,J,K,6:LE) / QH(I,J,K,1)
  ELSE
    RHO(I,J,K) = 0.0
    U(I,J,K)   = 0.0
    V(I,J,K)   = 0.0
    W(I,J,K)   = 0.0
    P(I,J,K)   = 0.0
    T(I,J,K)   = 0.0
    IF(.NOT. PRESENT(YY)) CYCLE
    YY(I,J,K,6:LE) = 0.0
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SetPhysics3D
!***********************************************************************
!**** 保存量流束から物理量を計算                                    ****
!**** 計算対象 : 単相, 二次元, 圧縮性, k-e                          ****
!***********************************************************************
SUBROUTINE SetPhysics2DKEM( &
&            RG, GAMMA, &
&            IS, IE, JS, JE, LS, LE, &
&            QH, AJA, RHO, U, V, P, T, AK, EPS, &
&            YY &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)            :: RG, GAMMA
  INTEGER, INTENT(IN)            :: IS, IE, JS, JE
  INTEGER, INTENT(IN)            :: LS, LE
  REAL,    INTENT(IN)            :: QH (IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(IN)            :: AJA(IS:IE, JS:JE)
  REAL,    INTENT(OUT)           :: RHO(IS:IE, JS:JE), &
  &                                 U  (IS:IE, JS:JE), &
  &                                 V  (IS:IE, JS:JE), &
  &                                 P  (IS:IE, JS:JE), &
  &                                 T  (IS:IE, JS:JE), &
  &                                 AK (IS:IE, JS:JE), &
  &                                 EPS(IS:IE, JS:JE)
  REAL,    INTENT(OUT), OPTIONAL :: YY (IS:IE, JS:JE, 7: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 6) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 保存量流束から物理量を計算 +++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
  DO J = JS, JE
  DO I = IS, IE
  IF(QH(I,J,1) .GT. 0.0) THEN
    RHO(I,J) = AJA(I,J) * QH(I,J,1)
    U(I,J)   = QH(I,J,2) / QH(I,J,1)
    V(I,J)   = QH(I,J,3) / QH(I,J,1)
    P(I,J)   = ( QH(I,J,4) - QH(I,J,5) &
    &          - 0.5 * (QH(I,J,2)**2 + QH(I,J,3)**2) / QH(I,J,1) &
    &        ) * (GAMMA - 1.0) * AJA(I,J)
    T(I,J)   = P(I,J) / (RHO(I,J) * RG)
    AK(I,J)  = QH(I,J,5) / QH(I,J,1)
    EPS(I,J) = QH(I,J,6) / QH(I,J,1)
    IF(.NOT. PRESENT(YY)) CYCLE
    YY(I,J,7:LE) = YY(I,J,7:LE) / QH(I,J,1)
  ELSE
    RHO(I,J) = 0.0
    U(I,J)   = 0.0
    V(I,J)   = 0.0
    P(I,J)   = 0.0
    T(I,J)   = 0.0
    AK(I,J)  = 0.0
    EPS(I,J) = 0.0
    IF(.NOT. PRESENT(YY)) CYCLE
    YY(I,J,7:LE) = 0.0
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SetPhysics2DKEM
!***********************************************************************
!**** 保存量流束から物理量を計算                                    ****
!**** 計算対象 : 単相, 三次元, 圧縮性, k-e                          ****
!***********************************************************************
SUBROUTINE SetPhysics3DKEM( &
&            RG, GAMMA, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            QH, AJA, RHO, U, V, W, P, T, AK, EPS, &
&            YY &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)            :: RG, GAMMA
  INTEGER, INTENT(IN)            :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)            :: LS, LE
  REAL,    INTENT(IN)            :: QH (IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(IN)            :: AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT)           :: RHO(IS:IE, JS:JE, KS:KE), &
  &                                 U  (IS:IE, JS:JE, KS:KE), &
  &                                 V  (IS:IE, JS:JE, KS:KE), &
  &                                 W  (IS:IE, JS:JE, KS:KE), &
  &                                 P  (IS:IE, JS:JE, KS:KE), &
  &                                 T  (IS:IE, JS:JE, KS:KE), &
  &                                 AK (IS:IE, JS:JE, KS:KE), &
  &                                 EPS(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT), OPTIONAL :: YY (IS:IE, JS:JE, KS:KE, 8: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 7) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 保存量流束から物理量を計算 +++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, K)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF( QH(I,J,K,1) .GT. 0.0 ) THEN
    RHO(I,J,K) = AJA(I,J,K) * QH(I,J,K,1)
    U(I,J,K)   = QH(I,J,K,2) / QH(I,J,K,1)
    V(I,J,K)   = QH(I,J,K,3) / QH(I,J,K,1)
    W(I,J,K)   = QH(I,J,K,4) / QH(I,J,K,1)
    P(I,J,K)   = ( QH(I,J,K,5) &
    &            - QH(I,J,K,6) &
    &            - 0.5 * ( QH(I,J,K,2)**2 &
    &                    + QH(I,J,K,3)**2 &
    &                    + QH(I,J,K,4)**2 ) / QH(I,J,K,1) &
    &          ) * (GAMMA - 1.0) * AJA(I,J,K)
    T(I,J,K)   = P(I,J,K) / (RHO(I,J,K) * RG)
    AK(I,J,K)  = QH(I,J,K,6) / QH(I,J,K,1)
    EPS(I,J,K) = QH(I,J,K,7) / QH(I,J,K,1)
    IF(.NOT. PRESENT(YY)) CYCLE
    YY(I,J,K,8:LE) = QH(I,J,K,8:LE) / QH(I,J,K,1)
  ELSE
    RHO(I,J,K) = 0.0
    U(I,J,K)   = 0.0
    V(I,J,K)   = 0.0
    W(I,J,K)   = 0.0
    P(I,J,K)   = 0.0
    T(I,J,K)   = 0.0
    AK(I,J,K)  = 0.0
    EPS(I,J,K) = 0.0
    IF(.NOT. PRESENT(YY)) CYCLE
    YY(I,J,K,8:LE) = 0.0
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SetPhysics3DKEM
!***********************************************************************
!**** 保存量流束から物理量を計算                                    ****
!**** 計算対象 : 単相, 二次元, 圧縮性, RSM                          ****
!***********************************************************************
SUBROUTINE SetPhysics2DRSM( &
&            RG, GAMMA, &
&            IS, IE, JS, JE, LS, LE, &
&            QH, AJA, RHO, U, V, P, T, uu, vv, ww, uv, EPS, &
&            YY &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)            :: RG, GAMMA
  INTEGER, INTENT(IN)            :: IS, IE, JS, JE
  INTEGER, INTENT(IN)            :: LS, LE
  REAL,    INTENT(IN)            :: QH (IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(IN)            :: AJA(IS:IE, JS:JE)
  REAL,    INTENT(OUT)           :: RHO(IS:IE, JS:JE), &
  &                                 U  (IS:IE, JS:JE), &
  &                                 V  (IS:IE, JS:JE), &
  &                                 P  (IS:IE, JS:JE), &
  &                                 T  (IS:IE, JS:JE), &
  &                                 uu (IS:IE, JS:JE), &
  &                                 vv (IS:IE, JS:JE), &
  &                                 ww (IS:IE, JS:JE), &
  &                                 uv (IS:IE, JS:JE), &
  &                                 EPS(IS:IE, JS:JE)
  REAL,    INTENT(OUT), OPTIONAL :: YY (IS:IE, JS:JE, 10: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 9) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 保存量流束から物理量を計算 +++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
  DO J = JS, JE
  DO I = IS, IE
  IF(QH(I,J,1) .GT. 0.0) THEN
    RHO(I,J) = AJA(I,J) * QH(I,J,1)
    U(I,J)   = QH(I,J,2) / QH(I,J,1)
    V(I,J)   = QH(I,J,3) / QH(I,J,1)
    P(I,J)   = ( QH(I,J,4) &
    &          - 0.5 * (QH(I,J,2)**2 + QH(I,J,3)**2) / QH(I,J,1) &
    &          - 0.5 * (QH(I,J,5) + QH(I,J,6) + QH(I,J,7)) &
    &        ) * (GAMMA - 1.0) * AJA(I,J)
    T(I,J)   = P(I,J) / (RHO(I,J) * RG)
    uu(I,J)  = QH(I,J,5) / QH(I,J,1)
    vv(I,J)  = QH(I,J,6) / QH(I,J,1)
    ww(I,J)  = QH(I,J,7) / QH(I,J,1)
    uv(I,J)  = QH(I,J,8) / QH(I,J,1)
    EPS(I,J) = QH(I,J,9) / QH(I,J,1)
    IF(.NOT. PRESENT(YY)) CYCLE
    YY(I,J,10:LE) = YY(I,J,10:LE) / QH(I,J,1)
  ELSE
    RHO(I,J) = 0.0
    U(I,J)   = 0.0
    V(I,J)   = 0.0
    P(I,J)   = 0.0
    T(I,J)   = 0.0
    uu(I,J)  = 0.0
    vv(I,J)  = 0.0
    ww(I,J)  = 0.0
    uv(I,J)  = 0.0
    EPS(I,J) = 0.0
    IF(.NOT. PRESENT(YY)) CYCLE
    YY(I,J,10:LE) = 0.0
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SetPhysics2DRSM
!***********************************************************************
!**** 保存量流束から物理量を計算                                    ****
!**** 計算対象 : 単相, 三次元, 圧縮性, RSM                          ****
!***********************************************************************
SUBROUTINE SetPhysics3DRSM( &
&            RG, GAMMA, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            QH, AJA, RHO, U, V, W, P, T, uu, vv, ww, uv, vw, wu, EPS, &
&            YY &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)            :: RG, GAMMA
  INTEGER, INTENT(IN)            :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)            :: LS, LE
  REAL,    INTENT(IN)            :: QH (IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(IN)            :: AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT)           :: RHO(IS:IE, JS:JE, KS:KE), &
  &                                 U  (IS:IE, JS:JE, KS:KE), &
  &                                 V  (IS:IE, JS:JE, KS:KE), &
  &                                 W  (IS:IE, JS:JE, KS:KE), &
  &                                 P  (IS:IE, JS:JE, KS:KE), &
  &                                 T  (IS:IE, JS:JE, KS:KE), &
  &                                 uu (IS:IE, JS:JE, KS:KE), &
  &                                 vv (IS:IE, JS:JE, KS:KE), &
  &                                 ww (IS:IE, JS:JE, KS:KE), &
  &                                 uv (IS:IE, JS:JE, KS:KE), &
  &                                 vw (IS:IE, JS:JE, KS:KE), &
  &                                 wu (IS:IE, JS:JE, KS:KE), &
  &                                 EPS(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT), OPTIONAL :: YY (IS:IE, JS:JE, KS:KE, 13: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 12) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 保存量流束から物理量を計算 +++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, K)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(QH(I,J,K,1) .GT. 0.0) THEN
    RHO(I,J,K) = AJA(I,J,K) * QH(I,J,K,1)
    U(I,J,K)   = QH(I,J,K,2) / QH(I,J,K,1)
    V(I,J,K)   = QH(I,J,K,3) / QH(I,J,K,1)
    W(I,J,K)   = QH(I,J,K,4) / QH(I,J,K,1)
    P(I,J,K)   = ( QH(I,J,K,5) &
    &            - 0.5 * ( QH(I,J,K,2)**2 &
    &                    + QH(I,J,K,3)**2 &
    &                    + QH(I,J,K,4)**2 ) / QH(I,J,K,1) &
    &            - 0.5 * (QH(I,J,K,6) + QH(I,J,K,7) + QH(I,J,K,8)) &
    &          ) * (GAMMA - 1.0) * AJA(I,J,K)
    T(I,J,K)   = P(I,J,K) / (RHO(I,J,K) * RG)
    uu(I,J,K)  = QH(I,J,K, 6) / QH(I,J,K,1)
    vv(I,J,K)  = QH(I,J,K, 7) / QH(I,J,K,1)
    ww(I,J,K)  = QH(I,J,K, 8) / QH(I,J,K,1)
    uv(I,J,K)  = QH(I,J,K, 9) / QH(I,J,K,1)
    vw(I,J,K)  = QH(I,J,K,10) / QH(I,J,K,1)
    wu(I,J,K)  = QH(I,J,K,11) / QH(I,J,K,1)
    EPS(I,J,K) = QH(I,J,K,12) / QH(I,J,K,1)
    IF(.NOT. PRESENT(YY)) CYCLE
    YY(I,J,K,13:LE) = QH(I,J,K,13:LE) / QH(I,J,K,1)
  ELSE
    RHO(I,J,K) = 0.0
    U(I,J,K)   = 0.0
    V(I,J,K)   = 0.0
    W(I,J,K)   = 0.0
    P(I,J,K)   = 0.0
    T(I,J,K)   = 0.0
    uu(I,J,K)  = 0.0
    vv(I,J,K)  = 0.0
    ww(I,J,K)  = 0.0
    uv(I,J,K)  = 0.0
    vw(I,J,K)  = 0.0
    wu(I,J,K)  = 0.0
    EPS(I,J,K) = 0.0
    IF(.NOT. PRESENT(YY)) CYCLE
    YY(I,J,K,13:LE) = 0.0
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SetPhysics3DRSM
!***********************************************************************
!**** 粘性係数計算                                                  ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!***********************************************************************
SUBROUTINE ViscosityCoefficient2D( &
&            AMUINF, TINF, S1, &
&            IS, IE, JS, JE, &
&            T, AMU &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: AMUINF, TINF, S1
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  REAL,    INTENT(IN)  :: T(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: AMU(IS:IE, JS:JE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
  DO J = JS, JE
  DO I = IS, IE
  IF(T(I,J) .GT. 0.0) THEN
    AMU(I,J) = AMUINF * (T(I,J) / TINF)**1.5 &
    &        * (TINF + S1) / (T(I,J) + S1)
  ELSE
    AMU(I,J) = 0.0
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE ViscosityCoefficient2D
!***********************************************************************
!**** 粘性係数計算                                                  ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!***********************************************************************
SUBROUTINE ViscosityCoefficient3D( &
&            AMUINF, TINF, S1, &
&            IS, IE, JS, JE, KS, KE, &
&            T, AMU &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: AMUINF, TINF, S1
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  REAL,    INTENT(IN)  :: T(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: AMU(IS:IE, JS:JE, KS:KE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, K)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(T(I,J,K) .GT. 0.0) THEN
    AMU(I,J,K) = AMUINF * (T(I,J,K) / TINF)**1.5 &
    &          * (TINF + S1) / (T(I,J,K) + S1)
  ELSE
    AMU(I,J,K) = 0.0
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE ViscosityCoefficient3D
!***********************************************************************
!**** 時間刻み計算 - 局所時間刻み法(Local time step method)         ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!***********************************************************************
SUBROUTINE CalLocalDt2D( &
&            CN, RG, GAMMA, &
&            IS, IE, JS, JE, &
&            XIX, XIY, ETX, ETY, &
&            RHO, U, V, T, AMU, &
&            DTLOCL &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: CN, RG, GAMMA
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE), XIY(IS:IE, JS:JE), &
  &                       ETX(IS:IE, JS:JE), ETY(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: U(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: V(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: T(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: DTLOCL(IS:IE, JS:JE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: UU, VV, A, nueff
  REAL    :: ALCXI, ALCET
  REAL    :: ALDXI, ALDET
  REAL    :: ALXI, ALET
  ! 処理開始 ***********************************************************
  ! 各格子点の局所時間刻みを計算 ---------------------------------------
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, UU, VV, A, nueff, ALCXI, ALCET, ALDXI, ALDET, ALXI, ALET &
  !$)
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J) .GT. 0.0) THEN
    ! 対流項
    UU = XIX(I,J) * U(I,J) + XIY(I,J) * V(I,J)
    VV = ETX(I,J) * U(I,J) + ETY(I,J) * V(I,J)
    A  = SQRT(GAMMA * RG * T(I,J))
    ALCXI = ABS(UU) + A * SQRT(XIX(I,J)**2 + XIY(I,J)**2)
    ALCET = ABS(VV) + A * SQRT(ETX(I,J)**2 + ETY(I,J)**2)
    ! 拡散項
    nueff = 2.0 * AMU(I,J) / RHO(I,J)
    ALDXI = nueff * (XIX(I,J)**2 + XIY(I,J)**2)
    ALDET = nueff * (ETX(I,J)**2 + ETY(I,J)**2)
    ! 全ての項による局所時間刻み
    ALXI = ALCXI + ALDXI
    ALET = ALCET + ALDET
    DTLOCL(I,J) = CN / MAX(ZERO, ALXI, ALET)
  ELSE
    DTLOCL(I,J) = 0.0
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE CalLocalDt2D
!***********************************************************************
!**** 時間刻み計算 - 局所時間刻み法(Local time step method)         ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!***********************************************************************
SUBROUTINE CalLocalDt3D( &
&            CN, RG, GAMMA, &
&            IS, IE, JS, JE, KS, KE, &
&            XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, &
&            RHO, U, V, W, T, AMU, &
&            DTLOCL &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: CN, RG, GAMMA
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE, KS:KE), &
  &                       XIY(IS:IE, JS:JE, KS:KE), &
  &                       XIZ(IS:IE, JS:JE, KS:KE), &
  &                       ETX(IS:IE, JS:JE, KS:KE), &
  &                       ETY(IS:IE, JS:JE, KS:KE), &
  &                       ETZ(IS:IE, JS:JE, KS:KE), &
  &                       ZEX(IS:IE, JS:JE, KS:KE), &
  &                       ZEY(IS:IE, JS:JE, KS:KE), &
  &                       ZEZ(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: U(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: V(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: W(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: T(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: DTLOCL(IS:IE, JS:JE, KS:KE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: UU, VV, WW, A, nueff
  REAL    :: ALCXI, ALCET, ALCZE
  REAL    :: ALDXI, ALDET, ALDZE
  REAL    :: ALXI, ALET, ALZE
  ! 処理開始 ***********************************************************
  ! 各格子点の局所時間刻みを計算 ---------------------------------------
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, UU, VV, WW, A, nueff, &
  !$  ALCXI, ALCET, ALCZE, ALDXI, ALDET, ALDZE, ALXI, ALET, ALZE &
  !$)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J,K) .GT. 0.0) THEN
    ! 対流項
    UU = XIX(I,J,K)*U(I,J,K)+XIY(I,J,K)*V(I,J,K)+XIZ(I,J,K)*W(I,J,K)
    VV = ETX(I,J,K)*U(I,J,K)+ETY(I,J,K)*V(I,J,K)+ETZ(I,J,K)*W(I,J,K)
    WW = ZEX(I,J,K)*U(I,J,K)+ZEY(I,J,K)*V(I,J,K)+ZEZ(I,J,K)*W(I,J,K)
    A  = SQRT(GAMMA * RG * T(I,J,K))
    ALCXI = ABS(UU) &
    &     + A * SQRT(XIX(I,J,K)**2 + XIY(I,J,K)**2 + XIZ(I,J,K)**2)
    ALCET = ABS(VV) &
    &     + A * SQRT(ETX(I,J,K)**2 + ETY(I,J,K)**2 + ETZ(I,J,K)**2)
    ALCZE = ABS(WW) &
    &     + A * SQRT(ZEX(I,J,K)**2 + ZEY(I,J,K)**2 + ZEZ(I,J,K)**2)
    ! 拡散項
    nueff = 2.0 * AMU(I,J,K) / RHO(I,J,K)
    ALDXI = nueff * (XIX(I,J,K)**2 + XIY(I,J,K)**2 + XIZ(I,J,K)**2)
    ALDET = nueff * (ETX(I,J,K)**2 + ETY(I,J,K)**2 + ETZ(I,J,K)**2)
    ALDZE = nueff * (ZEX(I,J,K)**2 + ZEY(I,J,K)**2 + ZEZ(I,J,K)**2)
    ! 全ての項による局所時間刻み
    ALXI = ALCXI + ALDXI
    ALET = ALCET + ALDET
    ALZE = ALCZE + ALDZE
    DTLOCL(I,J,K) = CN / MAX(ZERO, ALXI, ALET, ALZE)
   ELSE
    DTLOCL(I,J,K) = 0.0
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE CalLocalDt3D
!***********************************************************************
!**** 時間刻み計算                                                  ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!***********************************************************************
SUBROUTINE CalDt2D( &
&            IS, IE, JS, JE, DTLOCL, DTMIN, DTMAX, DT &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  REAL,    INTENT(IN)  :: DTLOCL(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: DTMIN, DTMAX
  REAL,    INTENT(OUT) :: DT
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  ! 処理開始 ***********************************************************
  DT = DTMAX
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF(DTLOCL(I,J) .GT. 0.0) THEN
    DT = MIN(DT, DTLOCL(I,J))
  ENDIF
  ENDDO
  ENDDO
  DT = MAX(DT, DTMIN)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE CalDt2D
!***********************************************************************
!**** 時間刻み計算                                                  ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!***********************************************************************
SUBROUTINE CalDt3D( &
&            IS, IE, JS, JE, KS, KE, DTLOCL, DTMIN, DTMAX, DT &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  REAL,    INTENT(IN)  :: DTLOCL(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: DTMIN, DTMAX
  REAL,    INTENT(OUT) :: DT
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  ! 処理開始 ***********************************************************
  DT = DTMAX
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF(DTLOCL(I,J,K) .GT. 0.0) THEN
    DT = MIN(DT, DTLOCL(I,J,K))
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  DT = MAX(DT, DTMIN)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE CalDt3D
!***********************************************************************
!**** Runge-Kutta法(Jameson-Baker, 1983) : NMAX段階, N段目          ****
!**** 計算対象 : 二次元                                             ****
!***********************************************************************
SUBROUTINE RungeKutta2D( &
&            IS, IE, JS, JE, LS, LE, &
&            N, NMAX, DTLOCL, &
&            DQH, QH0, QH1 &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: LS, LE
  INTEGER, INTENT(IN)  :: N, NMAX
  REAL,    INTENT(IN)  :: DTLOCL(IS: IE, JS: JE)
  REAL,    INTENT(IN)  :: DQH(IS: IE, JS: JE, LS: LE)
  REAL,    INTENT(IN)  :: QH0(IS: IE, JS: JE, LS: LE)
  REAL,    INTENT(OUT) :: QH1(IS: IE, JS: JE, LS: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  REAL    :: A
  ! 処理開始 ***********************************************************
  A = 1.0 / REAL(NMAX + 1 - N)
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, L)
  DO L = LS, LE
  !$OMP DO
  DO J = JS, JE
  DO I = IS, IE
  IF(QH0(I,J,1) .GT. 0.0) THEN
    QH1(I,J,L) = QH0(I,J,L) + A * DTLOCL(I,J) * DQH(I,J,L)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE RungeKutta2D
!***********************************************************************
!**** Runge-Kutta法(Jameson-Baker, 1983) : NMAX段階, N段目          ****
!**** 計算対象 : 三次元                                             ****
!***********************************************************************
SUBROUTINE RungeKutta3D( &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            N, NMAX, DTLOCL, &
&            DQH, QH0, QH1 &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  INTEGER, INTENT(IN)  :: N, NMAX
  REAL,    INTENT(IN)  :: DTLOCL(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)  :: DQH(IS: IE, JS: JE, KS: KE, LS: LE)
  REAL,    INTENT(IN)  :: QH0(IS: IE, JS: JE, KS: KE, LS: LE)
  REAL,    INTENT(OUT) :: QH1(IS: IE, JS: JE, KS: KE, LS: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  REAL    :: A
  ! 処理開始 ***********************************************************
  A = 1.0 / REAL(NMAX + 1 - N)
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  !$OMP DO
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(QH0(I,J,K,1) .GT. 0.0) THEN
    QH1(I,J,K,L) = QH0(I,J,K,L) + A * DTLOCL(I,J,K) * DQH(I,J,K,L)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE RungeKutta3D
!***********************************************************************
!**** 陰解法                                                        ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!**** 計算スキーム : デルタ形式LU-ADI法                             ****
!***********************************************************************
SUBROUTINE LUADI2D( &
&            lambda, &
&            RG, GAMMA, &
&            IS, IE, JS, JE, LS, LE, &
&            XIX, XIY, ETX, ETY, AJA, &
&            RHO, U, V, T, AMU, &
&            DTLOCL, &
&            RHS0, RHS1, QH0, QH1 &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: lambda
  REAL,    INTENT(IN)    :: RG, GAMMA
  INTEGER, INTENT(IN)    :: IS, IE, JS, JE
  INTEGER, INTENT(IN)    :: LS, LE
  REAL,    INTENT(IN)    :: XIX   (IS: IE, JS: JE), &
  &                         XIY   (IS: IE, JS: JE), &
  &                         ETX   (IS: IE, JS: JE), &
  &                         ETY   (IS: IE, JS: JE), &
  &                         AJA   (IS: IE, JS: JE)
  REAL,    INTENT(IN)    :: RHO   (IS: IE, JS: JE)
  REAL,    INTENT(IN)    :: U     (IS: IE, JS: JE), &
  &                         V     (IS: IE, JS: JE)
  REAL,    INTENT(IN)    :: T     (IS: IE, JS: JE)
  REAL,    INTENT(IN)    :: AMU   (IS: IE, JS: JE)
  REAL,    INTENT(IN)    :: DTLOCL(IS: IE, JS: JE)
  REAL,    INTENT(IN)    :: RHS0  (IS: IE, JS: JE, LS: LE)
  REAL,    INTENT(IN)    :: RHS1  (IS: IE, JS: JE, LS: LE)
  REAL,    INTENT(IN)    :: QH0   (IS: IE, JS: JE, LS: LE)
  REAL,    INTENT(INOUT) :: QH1   (IS: IE, JS: JE, LS: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER              :: I, J, L
  REAL                 :: UU, VV, C
  REAL                 :: THI, PHI, BETA, TIL, AKXT, AKYT, THIT
  REAL                 :: CT, VIS
  REAL,    ALLOCATABLE :: DQH(:, :, :)
  REAL,    ALLOCATABLE :: D1(:, :), D2(:, :)
  REAL,    ALLOCATABLE :: EE(:, :), EA(:, :), EP(:, :), EM(:, :)
  REAL,    ALLOCATABLE :: RR(:, :, :), RI(:, :, :)
  ! 処理開始 ***********************************************************
  ALLOCATE(DQH(IS: IE, JS: JE, LS: LE))
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE( &
  !$  I, J, L, UU, VV, C, &
  !$  THI, PHI, BETA, TIL, AKXT, AKYT, THIT, CT, VIS, &
  !$  D1, D2, EE, EA, EP, EM, RR, RI &
  !$)
  ! 初期設定 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO L = LS, LE
  !$OMP DO
  DO J = JS, JE
  DO I = IS, IE
  IF(QH0(I,J,1) .GT. 0.0) THEN
    DQH(I,J,L) = QH0(I,J,L) - QH1(I,J,L) &
    &          + DTLOCL(I,J) * ( &
    &              (1.0 - lambda) * RHS0(I,J,L) &
    &            + lambda         * RHS1(I,J,L) &
    &          )
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  ! XI方向 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! メモリ確保 =========================================================
  ALLOCATE(D1(IS: IE, LS: LE))
  ALLOCATE(D2(IS: IE, LS: LE))
  ALLOCATE(EE(IS: IE, LS: LE))
  ALLOCATE(EA(IS: IE, LS: LE))
  ALLOCATE(EP(IS: IE, LS: LE))
  ALLOCATE(EM(IS: IE, LS: LE))
  ALLOCATE(RR(IS: IE, LS:  4, LS:  4))
  ALLOCATE(RI(IS: IE, LS:  4, LS:  4))
  ! 行列反転処理 =======================================================
  !$OMP DO
  DO J = JS + 1, JE - 1
    ! 代入 -------------------------------------------------------------
    DO L = LS, LE
    DO I = IS, IE
    IF(QH0(I,J,1) .GT. 0.0) THEN
      D1(I,L) = DQH(I,J,L)
    ENDIF
    ENDDO
    ENDDO
    ! 対角化 -----------------------------------------------------------
    DO I = IS, IE
    IF(QH0(I,J,1) .GT. 0.0) THEN
      ! 命名
      UU   = U(I,J)
      VV   = V(I,J)
      C    = SQRT(GAMMA * RG * T(I,J))
      THI  = XIX(I,J) * UU + XIY(I,J) * VV
      PHI  = 0.5 * (GAMMA - 1.0) * (UU**2 + VV**2)
      BETA = 1.0 / (2.0 * C**2)
      TIL  = SQRT(XIX(I,J)**2 + XIY(I,J)**2)
      AKXT = XIX(I,J) / TIL
      AKYT = XIY(I,J) / TIL
      THIT = AKXT * UU + AKYT * VV
      CT   = C * TIL
      VIS  = 2.0 * AMU(I,J) / RHO(I,J) * TIL**2
      ! RR
      RR(I,1,1) = 1.0
      RR(I,1,2) = 0.0
      RR(I,1,3) = 1.0
      RR(I,1,4) = 1.0
      RR(I,2,1) = UU
      RR(I,2,2) = AKYT
      RR(I,2,3) = UU + AKXT * C
      RR(I,2,4) = UU - AKXT * C
      RR(I,3,1) = VV
      RR(I,3,2) =-AKXT
      RR(I,3,3) = VV + AKYT * C
      RR(I,3,4) = VV - AKYT * C
      RR(I,4,1) = PHI / (GAMMA - 1.0)
      RR(I,4,2) = AKYT * UU - AKXT * VV
      RR(I,4,3) = (PHI + C**2) / (GAMMA - 1.0) + C * THIT
      RR(I,4,4) = (PHI + C**2) / (GAMMA - 1.0) - C * THIT
      ! RI
      RI(I,1,1) = 1.0 - PHI / C**2
      RI(I,1,2) = (GAMMA - 1.0) * UU / C**2
      RI(I,1,3) = (GAMMA - 1.0) * VV / C**2
      RI(I,1,4) =-(GAMMA - 1.0) / C**2
      RI(I,2,1) =-AKYT * UU + AKXT * VV
      RI(I,2,2) = AKYT
      RI(I,2,3) =-AKXT
      RI(I,2,4) = 0.0
      RI(I,3,1) = BETA * (PHI - C * THIT)
      RI(I,3,2) = BETA * (AKXT * C - (GAMMA - 1.0) * UU)
      RI(I,3,3) = BETA * (AKYT * C - (GAMMA - 1.0) * VV)
      RI(I,3,4) = BETA * (GAMMA - 1.0)
      RI(I,4,1) = BETA * (PHI + C * THIT)
      RI(I,4,2) =-BETA * (AKXT * C + (GAMMA - 1.0) * UU)
      RI(I,4,3) =-BETA * (AKYT * C + (GAMMA - 1.0) * VV)
      RI(I,4,4) = BETA * (GAMMA - 1.0)
      ! EE
      EE(I,1) = THI
      EE(I,2) = THI
      EE(I,3) = THI + CT
      EE(I,4) = THI - CT
      DO L = 5, LE
        EE(I,L) = THI
      ENDDO
      ! EA
      DO L = LS, LE
        EA(I,L) = ABS(EE(I,L)) + VIS ! + DEL * CT
      ENDDO
      ! EP
      DO L = LS, LE
        EP(I,L) = 0.5 * (EE(I,L) + EA(I,L))
      ENDDO
      ! EM
      DO L = LS, LE
        EM(I,L) = 0.5 * (EE(I,L) - EA(I,L))
      ENDDO
    ENDIF
    ENDDO
    ! RR反転 -----------------------------------------------------------
    DO I = IS, IE
    IF(QH0(I,J,1) .GT. 0.0) THEN
      DO L = LS, 4
        D2(I,L) = RI(I,L,1) * D1(I,1) + RI(I,L,2) * D1(I,2) &
        &       + RI(I,L,3) * D1(I,3) + RI(I,L,4) * D1(I,4)
      ENDDO
      DO L = 5, LE
        D2(I,L) = D1(I,L)
      ENDDO
    ENDIF
    ENDDO
    ! 前進スイープ -----------------------------------------------------
    D1(:,:) = 0.0
    DO I = IS + 1, IE - 1, + 1
    IF(QH0(I-1,J,1) .GT. 0.0 .AND. QH0(I,J,1) .GT. 0.0) THEN
      DO L = LS, LE
        D1(I,L) = ( D2(I,L) + DTLOCL(I,J) * lambda &
        &         * EP(I-1,L) * D1(I-1,L) ) &
        &       / (1.0 + DTLOCL(I,J) * lambda * EA(I,L))
      ENDDO
    ENDIF
    ENDDO
    ! スカラー演算 -----------------------------------------------------
    DO I = IS, IE
    IF(QH0(I,J,1) .GT. 0.0) THEN
      DO L = LS, LE
        D1(I,L) = D1(I,L) * (1.0 + DTLOCL(I,J) * lambda * EA(I,L))
      ENDDO
    ENDIF
    ENDDO
    ! 後退スイープ -----------------------------------------------------
    D2(:,:) = 0.0
    DO I = IE - 1, IS + 1, - 1
    IF(QH0(I,J,1) .GT. 0.0 .AND. QH0(I+1,J,1) .GT. 0.0) THEN
      DO L = LS, LE
        D2(I,L) = ( D1(I,L) - DTLOCL(I,J) * lambda &
        &         * EM(I+1,L) * D2(I+1,L) ) &
        &       / (1.0 + DTLOCL(I,J) * lambda * EA(I,L))
      ENDDO
    ENDIF
    ENDDO
    ! RI反転 -----------------------------------------------------------
    DO I = IS, IE
    IF(QH0(I,J,1) .GT. 0.0) THEN
      DO L = LS, 4
        D1(I,L) = RR(I,L,1) * D2(I,1) + RR(I,L,2) * D2(I,2) &
        &       + RR(I,L,3) * D2(I,3) + RR(I,L,4) * D2(I,4)
      ENDDO
      DO L = 5, LE
        D1(I,L) = D2(I,L)
      ENDDO
    ENDIF
    ENDDO
    ! 代入 -------------------------------------------------------------
    DO L = LS, LE
    DO I = IS, IE
    IF(QH0(I,J,1) .GT. 0.0) THEN
      DQH(I,J,L) = D1(I,L)
    ENDIF
    ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  ! メモリ解放 =========================================================
  DEALLOCATE(D1, D2, EE, EA, EP, EM, RR, RI)
  ! ET方向 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! メモリ確保 =========================================================
  ALLOCATE(D1(JS: JE, LS: LE))
  ALLOCATE(D2(JS: JE, LS: LE))
  ALLOCATE(EE(JS: JE, LS: LE))
  ALLOCATE(EA(JS: JE, LS: LE))
  ALLOCATE(EP(JS: JE, LS: LE))
  ALLOCATE(EM(JS: JE, LS: LE))
  ALLOCATE(RR(JS: JE, LS:  4, LS:  4))
  ALLOCATE(RI(JS: JE, LS:  4, LS:  4))
  ! 行列反転処理 =======================================================
  !$OMP DO
  DO I = IS + 1, IE - 1
    ! 代入 -------------------------------------------------------------
    DO L = LS, LE
    DO J = JS, JE
    IF(QH0(I,J,1) .GT. 0.0) THEN
      D1(J,L) = DQH(I,J,L)
    ENDIF
    ENDDO
    ENDDO
    ! 対角化 -----------------------------------------------------------
    DO J = JS, JE
    IF(QH0(I,J,1) .GT. 0.0) THEN
      ! 命名
      UU   = U(I,J)
      VV   = V(I,J)
      C    = SQRT(GAMMA * RG * T(I,J))
      THI  = ETX(I,J) * UU + ETY(I,J) * VV
      PHI  = 0.5 * (GAMMA - 1.0) * (UU**2 + VV**2)
      BETA = 1.0 / (2.0 * C**2)
      TIL  = SQRT(ETX(I,J)**2 + ETY(I,J)**2)
      AKXT = ETX(I,J) / TIL
      AKYT = ETY(I,J) / TIL
      THIT = AKXT * UU + AKYT * VV
      CT   = C * TIL
      VIS  = 2.0 * AMU(I,J) / RHO(I,J) * TIL**2
      ! RR
      RR(J,1,1) = 1.0
      RR(J,1,2) = 0.0
      RR(J,1,3) = 1.0
      RR(J,1,4) = 1.0
      RR(J,2,1) = UU
      RR(J,2,2) = AKYT
      RR(J,2,3) = UU + AKXT * C
      RR(J,2,4) = UU - AKXT * C
      RR(J,3,1) = VV
      RR(J,3,2) =-AKXT
      RR(J,3,3) = VV + AKYT * C
      RR(J,3,4) = VV - AKYT * C
      RR(J,4,1) = PHI / (GAMMA - 1.0)
      RR(J,4,2) = AKYT * UU - AKXT * VV
      RR(J,4,3) = (PHI + C**2) / (GAMMA - 1.0) + C * THIT
      RR(J,4,4) = (PHI + C**2) / (GAMMA - 1.0) - C * THIT
      ! RI
      RI(J,1,1) = 1.0 - PHI / C**2
      RI(J,1,2) = (GAMMA - 1.0) * UU / C**2
      RI(J,1,3) = (GAMMA - 1.0) * VV / C**2
      RI(J,1,4) =-(GAMMA - 1.0) / C**2
      RI(J,2,1) =-AKYT * UU + AKXT * VV
      RI(J,2,2) = AKYT
      RI(J,2,3) =-AKXT
      RI(J,2,4) = 0.0
      RI(J,3,1) = BETA * (PHI - C * THIT)
      RI(J,3,2) = BETA * (AKXT * C - (GAMMA - 1.0) * UU)
      RI(J,3,3) = BETA * (AKYT * C - (GAMMA - 1.0) * VV)
      RI(J,3,4) = BETA * (GAMMA - 1.0)
      RI(J,4,1) = BETA * (PHI + C * THIT)
      RI(J,4,2) =-BETA * (AKXT * C + (GAMMA - 1.0) * UU)
      RI(J,4,3) =-BETA * (AKYT * C + (GAMMA - 1.0) * VV)
      RI(J,4,4) = BETA * (GAMMA - 1.0)
      ! EE
      EE(J,1) = THI
      EE(J,2) = THI
      EE(J,3) = THI + CT
      EE(J,4) = THI - CT
      DO L = 5, LE
        EE(J,L) = THI
      ENDDO
      ! EA
      DO L = LS, LE
        EA(J,L) = ABS(EE(J,L)) + VIS ! + DEL * CT
      ENDDO
      ! EP
      DO L = LS, LE
        EP(J,L) = 0.5 * (EE(J,L) + EA(J,L))
      ENDDO
      ! EM
      DO L = LS, LE
        EM(J,L) = 0.5 * (EE(J,L) - EA(J,L))
      ENDDO
    ENDIF
    ENDDO
    ! RR反転 -----------------------------------------------------------
    DO J = JS, JE
    IF(QH0(I,J,1) .GT. 0.0) THEN
      DO L = LS, 4
        D2(J,L) = RI(J,L,1) * D1(J,1) + RI(J,L,2) * D1(J,2) &
        &       + RI(J,L,3) * D1(J,3) + RI(J,L,4) * D1(J,4)
      ENDDO
      DO L = 5, LE
        D2(J,L) = D1(J,L)
      ENDDO
    ENDIF
    ENDDO
    ! 前進スイープ -----------------------------------------------------
    D1(:,:) = 0.0
    DO J = JS + 1, JE - 1, + 1
    IF(QH0(I,J-1,1) .GT. 0.0 .AND. QH0(I,J,1) .GT. 0.0) THEN
      DO L = LS, LE
        D1(J,L) = ( D2(J,L) + DTLOCL(I,J) * lambda &
        &         * EP(J-1,L) * D1(J-1,L) ) &
        &       / (1.0 + DTLOCL(I,J) * lambda * EA(J,L))
      ENDDO
    ENDIF
    ENDDO
    ! スカラー演算 -----------------------------------------------------
    DO J = JS, JE
    IF(QH0(I,J,1) .GT. 0.0) THEN
      DO L = LS, LE
        D1(J,L) = D1(J,L) * (1.0 + DTLOCL(I,J) * lambda * EA(J,L))
      ENDDO
    ENDIF
    ENDDO
    ! 後退スイープ -----------------------------------------------------
    D2(:,:) = 0.0
    DO J = JE - 1, JS + 1, - 1
    IF(QH0(I,J,1) .GT. 0.0 .AND. QH0(I,J+1,1) .GT. 0.0) THEN
      DO L = LS, LE
        D2(J,L) = ( D1(J,L) - DTLOCL(I,J) * lambda &
        &         * EM(J+1,L) * D2(J+1,L) ) &
        &       / (1.0 + DTLOCL(I,J) * lambda * EA(J,L))
      ENDDO
    ENDIF
    ENDDO
    ! RI反転 -----------------------------------------------------------
    DO J = JS, JE
    IF(QH0(I,J,1) .GT. 0.0) THEN
      DO L = LS, 4
        D1(J,L) = RR(J,L,1) * D2(J,1) + RR(J,L,2) * D2(J,2) &
        &       + RR(J,L,3) * D2(J,3) + RR(J,L,4) * D2(J,4)
      ENDDO
      DO L = 5, LE
        D1(J,L) = D2(J,L)
      ENDDO
    ENDIF
    ENDDO
    ! 代入 -------------------------------------------------------------
    DO L = LS, LE
    DO J = JS, JE
    IF(QH0(I,J,1) .GT. 0.0) THEN
      DQH(I,J,L) = D1(J,L)
    ENDIF
    ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  ! メモリ解放 =========================================================
  DEALLOCATE(D1, D2, EE, EA, EP, EM, RR, RI)
  ! デルタ形式から流束関数の時間を進める +++++++++++++++++++++++++++++++
  DO L = LS, LE
  !$OMP DO
  DO J = JS, JE
  DO I = IS, IE
  IF(QH0(I,J,1) .GT. 0.0) THEN
    QH1(I,J,L) = QH1(I,J,L) + DQH(I,J,L)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  ! メモリ解放 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP END PARALLEL
  DEALLOCATE(DQH)
  ! 処理終了 ***********************************************************
  RETURN
! 定義終了 *************************************************************
END SUBROUTINE LUADI2D
!***********************************************************************
!**** 陰解法                                                        ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!**** 計算スキーム : デルタ形式LU-ADI法                             ****
!***********************************************************************
SUBROUTINE LUADI3D( &
&            lambda, &
&            RG, GAMMA, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, AJA, &
&            RHO, U, V, W, T, AMU, &
&            DTLOCL, &
&            RHS0, RHS1, QH0, QH1 &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: lambda
  REAL,    INTENT(IN)    :: RG, GAMMA
  INTEGER, INTENT(IN)    :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)    :: LS, LE
  REAL,    INTENT(IN)    :: XIX   (IS: IE, JS: JE, KS: KE), &
  &                         XIY   (IS: IE, JS: JE, KS: KE), &
  &                         XIZ   (IS: IE, JS: JE, KS: KE), &
  &                         ETX   (IS: IE, JS: JE, KS: KE), &
  &                         ETY   (IS: IE, JS: JE, KS: KE), &
  &                         ETZ   (IS: IE, JS: JE, KS: KE), &
  &                         ZEX   (IS: IE, JS: JE, KS: KE), &
  &                         ZEY   (IS: IE, JS: JE, KS: KE), &
  &                         ZEZ   (IS: IE, JS: JE, KS: KE), &
  &                         AJA   (IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)    :: RHO   (IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)    :: U     (IS: IE, JS: JE, KS: KE), &
  &                         V     (IS: IE, JS: JE, KS: KE), &
  &                         W     (IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)    :: T     (IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)    :: AMU   (IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)    :: DTLOCL(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)    :: RHS0  (IS: IE, JS: JE, KS: KE, LS: LE)
  REAL,    INTENT(IN)    :: RHS1  (IS: IE, JS: JE, KS: KE, LS: LE)
  REAL,    INTENT(IN)    :: QH0   (IS: IE, JS: JE, KS: KE, LS: LE)
  REAL,    INTENT(INOUT) :: QH1   (IS: IE, JS: JE, KS: KE, LS: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER              :: I, J, K, L
  REAL                 :: UU, VV, WW, C
  REAL                 :: THI, PHI, BETA, TIL, AKXT, AKYT, AKZT, THIT
  REAL                 :: CT, VIS
  REAL,    ALLOCATABLE :: DQH(:, :, :, :)
  REAL,    ALLOCATABLE :: D1(:, :), D2(:, :)
  REAL,    ALLOCATABLE :: EE(:, :), EA(:, :), EP(:, :), EM(:, :)
  REAL,    ALLOCATABLE :: RR(:, :, :), RI(:, :, :)
  ! 処理開始 ***********************************************************
  ALLOCATE(DQH(IS: IE, JS: JE, KS: KE, LS: LE))
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, L, UU, VV, WW, C, &
  !$  THI, PHI, BETA, TIL, AKXT, AKYT, AKZT, THIT, CT, VIS, &
  !$  D1, D2, EE, EA, EP, EM, RR, RI &
  !$)
  ! 初期設定 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO L = LS, LE
  !$OMP DO
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(QH0(I,J,K,1) .GT. 0.0) THEN
    DQH(I,J,K,L) = QH0(I,J,K,L) - QH1(I,J,K,L) &
    &            + DTLOCL(I,J,K) * ( &
    &                (1.0 - lambda) * RHS0(I,J,K,L) &
    &              + lambda         * RHS1(I,J,K,L) &
    &            )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  ! XI方向 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! メモリ確保 =========================================================
  ALLOCATE(D1(IS: IE, LS: LE))
  ALLOCATE(D2(IS: IE, LS: LE))
  ALLOCATE(EE(IS: IE, LS: LE))
  ALLOCATE(EA(IS: IE, LS: LE))
  ALLOCATE(EP(IS: IE, LS: LE))
  ALLOCATE(EM(IS: IE, LS: LE))
  ALLOCATE(RR(IS: IE, LS:  5, LS:  5))
  ALLOCATE(RI(IS: IE, LS:  5, LS:  5))
  ! 行列反転処理 =======================================================
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
    ! 代入 -------------------------------------------------------------
    DO L = LS, LE
    DO I = IS, IE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      D1(I,L) = DQH(I,J,K,L)
    ENDIF
    ENDDO
    ENDDO
    ! 対角化 -----------------------------------------------------------
    DO I = IS, IE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      ! 命名
      UU   = U(I,J,K)
      VV   = V(I,J,K)
      WW   = W(I,J,K)
      C    = SQRT(GAMMA * RG * T(I,J,K))
      THI  = XIX(I,J,K) * UU + XIY(I,J,K) * VV + XIZ(I,J,K) * WW
      PHI  = 0.5 * (GAMMA - 1.0) * (UU**2 + VV**2 + WW**2)
      BETA = 1.0 / (2.0 * C**2)
      TIL  = SQRT(XIX(I,J,K)**2 + XIY(I,J,K)**2 + XIZ(I,J,K)**2)
      AKXT = XIX(I,J,K) / TIL
      AKYT = XIY(I,J,K) / TIL
      AKZT = XIZ(I,J,K) / TIL
      THIT = AKXT * UU + AKYT * VV + AKZT * WW
      CT   = C * TIL
      VIS  = 2.0 * AMU(I,J,K) / RHO(I,J,K) * TIL**2
      ! RR
      RR(I,1,1) = AKXT
      RR(I,2,1) = AKXT * UU
      RR(I,3,1) = AKXT * VV + AKZT
      RR(I,4,1) = AKXT * WW - AKYT
      RR(I,5,1) = AKXT * PHI / (GAMMA - 1.0) + AKZT * VV - AKYT * WW
      RR(I,1,2) = AKYT
      RR(I,2,2) = AKYT * UU - AKZT
      RR(I,3,2) = AKYT * VV
      RR(I,4,2) = AKYT * WW + AKXT
      RR(I,5,2) = AKYT * PHI / (GAMMA - 1.0) + AKXT * WW - AKZT * UU
      RR(I,1,3) = AKZT
      RR(I,2,3) = AKZT * UU + AKYT
      RR(I,3,3) = AKZT * VV - AKXT
      RR(I,4,3) = AKZT * WW
      RR(I,5,3) = AKZT * PHI / (GAMMA - 1.0) + AKYT * UU - AKXT * VV
      RR(I,1,4) = 1.0
      RR(I,2,4) = UU + AKXT * C
      RR(I,3,4) = VV + AKYT * C
      RR(I,4,4) = WW + AKZT * C
      RR(I,5,4) = (PHI + C**2) / (GAMMA - 1.0) + C * THIT
      RR(I,1,5) = 1.0
      RR(I,2,5) = UU - AKXT * C
      RR(I,3,5) = VV - AKYT * C
      RR(I,4,5) = WW - AKZT * C
      RR(I,5,5) = (PHI + C**2) / (GAMMA - 1.0) - C * THIT
      ! RI
      RI(I,1,1) = AKXT - AKXT * PHI / C**2 + AKYT * WW - AKZT * VV
      RI(I,2,1) = AKYT - AKYT * PHI / C**2 + AKZT * UU - AKXT * WW
      RI(I,3,1) = AKZT - AKZT * PHI / C**2 + AKXT * VV - AKYT * UU
      RI(I,4,1) = BETA * (PHI - C * THIT)
      RI(I,5,1) = BETA * (PHI + C * THIT)
      RI(I,1,2) = (GAMMA - 1.0) * AKXT * UU / C**2
      RI(I,2,2) = (GAMMA - 1.0) * AKYT * UU / C**2 - AKZT
      RI(I,3,2) = (GAMMA - 1.0) * AKZT * UU / C**2 + AKYT
      RI(I,4,2) = BETA * (AKXT * C - (GAMMA - 1.0) * UU)
      RI(I,5,2) =-BETA * (AKXT * C + (GAMMA - 1.0) * UU)
      RI(I,1,3) = (GAMMA - 1.0) * AKXT * VV / C**2 + AKZT
      RI(I,2,3) = (GAMMA - 1.0) * AKYT * VV / C**2
      RI(I,3,3) = (GAMMA - 1.0) * AKZT * VV / C**2 - AKXT
      RI(I,4,3) = BETA * (AKYT * C - (GAMMA - 1.0) * VV)
      RI(I,5,3) =-BETA * (AKYT * C + (GAMMA - 1.0) * VV)
      RI(I,1,4) = (GAMMA - 1.0) * AKXT * WW / C**2 - AKYT
      RI(I,2,4) = (GAMMA - 1.0) * AKYT * WW / C**2 + AKXT
      RI(I,3,4) = (GAMMA - 1.0) * AKZT * WW / C**2
      RI(I,4,4) = BETA * (AKZT * C - (GAMMA - 1.0) * WW)
      RI(I,5,4) =-BETA * (AKZT * C + (GAMMA - 1.0) * WW)
      RI(I,1,5) =-(GAMMA - 1.0) * AKXT / C**2
      RI(I,2,5) =-(GAMMA - 1.0) * AKYT / C**2
      RI(I,3,5) =-(GAMMA - 1.0) * AKZT / C**2
      RI(I,4,5) = BETA * (GAMMA - 1.0)
      RI(I,5,5) = BETA * (GAMMA - 1.0)
      ! EE
      EE(I,1) = THI
      EE(I,2) = THI
      EE(I,3) = THI
      EE(I,4) = THI + CT
      EE(I,5) = THI - CT
      DO L = 6, LE
        EE(I,L) = THI
      ENDDO
      ! EA
      DO L = LS, LE
        EA(I,L) = ABS(EE(I,L)) + VIS ! + DEL * CT
      ENDDO
      ! EP
      DO L = LS, LE
        EP(I,L) = 0.5 * (EE(I,L) + EA(I,L))
      ENDDO
      ! EM
      DO L = LS, LE
        EM(I,L) = 0.5 * (EE(I,L) - EA(I,L))
      ENDDO
    ENDIF
    ENDDO
    ! RR反転 -----------------------------------------------------------
    DO I = IS, IE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      DO L = LS, 5
        D2(I,L) = RI(I,L,1) * D1(I,1) + RI(I,L,2) * D1(I,2) &
        &       + RI(I,L,3) * D1(I,3) + RI(I,L,4) * D1(I,4) &
        &       + RI(I,L,5) * D1(I,5)
      ENDDO
      DO L = 6, LE
        D2(I,L) = D1(I,L)
      ENDDO
    ENDIF
    ENDDO
    ! 前進スイープ -----------------------------------------------------
    D1(:,:) = 0.0
    DO I = IS + 1, IE - 1, + 1
    IF(QH0(I-1,J,K,1) .GT. 0.0 .AND. QH0(I,J,K,1) .GT. 0.0) THEN
      DO L = LS, LE
        D1(I,L) = ( D2(I,L) + DTLOCL(I,J,K) * lambda &
        &         * EP(I-1,L) * D1(I-1,L) ) &
        &       / (1.0 + DTLOCL(I,J,K) * lambda * EA(I,L))
      ENDDO
    ENDIF
    ENDDO
    ! スカラー演算 -----------------------------------------------------
    DO I = IS, IE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      DO L = LS, LE
        D1(I,L) = D1(I,L) * (1.0 + DTLOCL(I,J,K) * lambda * EA(I,L))
      ENDDO
    ENDIF
    ENDDO
    ! 後退スイープ -----------------------------------------------------
    D2(:,:) = 0.0
    DO I = IE - 1, IS + 1, - 1
    IF(QH0(I,J,K,1) .GT. 0.0 .AND. QH0(I+1,J,K,1) .GT. 0.0) THEN
      DO L = LS, LE
        D2(I,L) = ( D1(I,L) - DTLOCL(I,J,K) * lambda &
        &         * EM(I+1,L) * D2(I+1,L) ) &
        &       / (1.0 + DTLOCL(I,J,K) * lambda * EA(I,L))
      ENDDO
    ENDIF
    ENDDO
    ! RI反転 -----------------------------------------------------------
    DO I = IS, IE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      DO L = LS, 5
        D1(I,L) = RR(I,L,1) * D2(I,1) + RR(I,L,2) * D2(I,2) &
        &       + RR(I,L,3) * D2(I,3) + RR(I,L,4) * D2(I,4) &
        &       + RR(I,L,5) * D2(I,5)
      ENDDO
      DO L = 6, LE
        D1(I,L) = D2(I,L)
      ENDDO
    ENDIF
    ENDDO
    ! 代入 -------------------------------------------------------------
    DO L = LS, LE
    DO I = IS, IE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      DQH(I,J,K,L) = D1(I,L)
    ENDIF
    ENDDO
    ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ! メモリ解放 =========================================================
  DEALLOCATE(D1, D2, EE, EA, EP, EM, RR, RI)
  ! ET方向 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! メモリ確保 =========================================================
  ALLOCATE(D1(JS: JE, LS: LE))
  ALLOCATE(D2(JS: JE, LS: LE))
  ALLOCATE(EE(JS: JE, LS: LE))
  ALLOCATE(EA(JS: JE, LS: LE))
  ALLOCATE(EP(JS: JE, LS: LE))
  ALLOCATE(EM(JS: JE, LS: LE))
  ALLOCATE(RR(JS: JE, LS:  5, LS:  5))
  ALLOCATE(RI(JS: JE, LS:  5, LS:  5))
  ! 行列反転処理 =======================================================
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO I = IS + 1, IE - 1
    ! 代入 -------------------------------------------------------------
    DO L = LS, LE
    DO J = JS, JE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      D1(J,L) = DQH(I,J,K,L)
    ENDIF
    ENDDO
    ENDDO
    ! 対角化 -----------------------------------------------------------
    DO J = JS, JE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      ! 命名
      UU   = U(I,J,K)
      VV   = V(I,J,K)
      WW   = W(I,J,K)
      C    = SQRT(GAMMA * RG * T(I,J,K))
      THI  = ETX(I,J,K) * UU + ETY(I,J,K) * VV + ETZ(I,J,K) * WW
      PHI  = 0.5 * (GAMMA - 1.0) * (UU**2 + VV**2 + WW**2)
      BETA = 1.0 / (2.0 * C**2)
      TIL  = SQRT(ETX(I,J,K)**2 + ETY(I,J,K)**2 + ETZ(I,J,K)**2)
      AKXT = ETX(I,J,K) / TIL
      AKYT = ETY(I,J,K) / TIL
      AKZT = ETZ(I,J,K) / TIL
      THIT = AKXT * UU + AKYT * VV + AKZT * WW
      CT   = C * TIL
      VIS  = 2.0 * AMU(I,J,K) / RHO(I,J,K) * TIL**2
      ! RR
      RR(J,1,1) = AKXT
      RR(J,2,1) = AKXT * UU
      RR(J,3,1) = AKXT * VV + AKZT
      RR(J,4,1) = AKXT * WW - AKYT
      RR(J,5,1) = AKXT * PHI / (GAMMA - 1.0) + AKZT * VV - AKYT * WW
      RR(J,1,2) = AKYT
      RR(J,2,2) = AKYT * UU - AKZT
      RR(J,3,2) = AKYT * VV
      RR(J,4,2) = AKYT * WW + AKXT
      RR(J,5,2) = AKYT * PHI / (GAMMA - 1.0) + AKXT * WW - AKZT * UU
      RR(J,1,3) = AKZT
      RR(J,2,3) = AKZT * UU + AKYT
      RR(J,3,3) = AKZT * VV - AKXT
      RR(J,4,3) = AKZT * WW
      RR(J,5,3) = AKZT * PHI / (GAMMA - 1.0) + AKYT * UU - AKXT * VV
      RR(J,1,4) = 1.0
      RR(J,2,4) = UU + AKXT * C
      RR(J,3,4) = VV + AKYT * C
      RR(J,4,4) = WW + AKZT * C
      RR(J,5,4) = (PHI + C**2) / (GAMMA - 1.0) + C * THIT
      RR(J,1,5) = 1.0
      RR(J,2,5) = UU - AKXT * C
      RR(J,3,5) = VV - AKYT * C
      RR(J,4,5) = WW - AKZT * C
      RR(J,5,5) = (PHI + C**2) / (GAMMA - 1.0) - C * THIT
      ! RI
      RI(J,1,1) = AKXT - AKXT * PHI / C**2 + AKYT * WW - AKZT * VV
      RI(J,2,1) = AKYT - AKYT * PHI / C**2 + AKZT * UU - AKXT * WW
      RI(J,3,1) = AKZT - AKZT * PHI / C**2 + AKXT * VV - AKYT * UU
      RI(J,4,1) = BETA * (PHI - C * THIT)
      RI(J,5,1) = BETA * (PHI + C * THIT)
      RI(J,1,2) = (GAMMA - 1.0) * AKXT * UU / C**2
      RI(J,2,2) = (GAMMA - 1.0) * AKYT * UU / C**2 - AKZT
      RI(J,3,2) = (GAMMA - 1.0) * AKZT * UU / C**2 + AKYT
      RI(J,4,2) = BETA * (AKXT * C - (GAMMA - 1.0) * UU)
      RI(J,5,2) =-BETA * (AKXT * C + (GAMMA - 1.0) * UU)
      RI(J,1,3) = (GAMMA - 1.0) * AKXT * VV / C**2 + AKZT
      RI(J,2,3) = (GAMMA - 1.0) * AKYT * VV / C**2
      RI(J,3,3) = (GAMMA - 1.0) * AKZT * VV / C**2 - AKXT
      RI(J,4,3) = BETA * (AKYT * C - (GAMMA - 1.0) * VV)
      RI(J,5,3) =-BETA * (AKYT * C + (GAMMA - 1.0) * VV)
      RI(J,1,4) = (GAMMA - 1.0) * AKXT * WW / C**2 - AKYT
      RI(J,2,4) = (GAMMA - 1.0) * AKYT * WW / C**2 + AKXT
      RI(J,3,4) = (GAMMA - 1.0) * AKZT * WW / C**2
      RI(J,4,4) = BETA * (AKZT * C - (GAMMA - 1.0) * WW)
      RI(J,5,4) =-BETA * (AKZT * C + (GAMMA - 1.0) * WW)
      RI(J,1,5) =-(GAMMA - 1.0) * AKXT / C**2
      RI(J,2,5) =-(GAMMA - 1.0) * AKYT / C**2
      RI(J,3,5) =-(GAMMA - 1.0) * AKZT / C**2
      RI(J,4,5) = BETA * (GAMMA - 1.0)
      RI(J,5,5) = BETA * (GAMMA - 1.0)
      ! EE
      EE(J,1) = THI
      EE(J,2) = THI
      EE(J,3) = THI
      EE(J,4) = THI + CT
      EE(J,5) = THI - CT
      DO L = 6, LE
        EE(J,L) = THI
      ENDDO
      ! EA
      DO L = LS, LE
        EA(J,L) = ABS(EE(J,L)) + VIS ! + DEL * CT
      ENDDO
      ! EP
      DO L = LS, LE
        EP(J,L) = 0.5 * (EE(J,L) + EA(J,L))
      ENDDO
      ! EM
      DO L = LS, LE
        EM(J,L) = 0.5 * (EE(J,L) - EA(J,L))
      ENDDO
    ENDIF
    ENDDO
    ! RR反転 -----------------------------------------------------------
    DO J = JS, JE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      DO L = LS, 5
        D2(J,L) = RI(J,L,1) * D1(J,1) + RI(J,L,2) * D1(J,2) &
        &       + RI(J,L,3) * D1(J,3) + RI(J,L,4) * D1(J,4) &
        &       + RI(J,L,5) * D1(J,5)
      ENDDO
      DO L = 6, LE
        D2(J,L) = D1(J,L)
      ENDDO
    ENDIF
    ENDDO
    ! 前進スイープ -----------------------------------------------------
    D1(:,:) = 0.0
    DO J = JS + 1, JE - 1, + 1
    IF(QH0(I,J-1,K,1) .GT. 0.0 .AND. QH0(I,J,K,1) .GT. 0.0) THEN
      DO L = LS, LE
        D1(J,L) = ( D2(J,L) + DTLOCL(I,J,K) * lambda &
        &         * EP(J-1,L) * D1(J-1,L) ) &
        &       / (1.0 + DTLOCL(I,J,K) * lambda * EA(J,L))
      ENDDO
    ENDIF
    ENDDO
    ! スカラー演算 -----------------------------------------------------
    DO J = JS, JE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      DO L = LS, LE
        D1(J,L) = D1(J,L) * (1.0 + DTLOCL(I,J,K) * lambda * EA(J,L))
      ENDDO
    ENDIF
    ENDDO
    ! 後退スイープ -----------------------------------------------------
    D2(:,:) = 0.0
    DO J = JE - 1, JS + 1, - 1
    IF(QH0(I,J,K,1) .GT. 0.0 .AND. QH0(I,J+1,K,1) .GT. 0.0) THEN
      DO L = LS, LE
        D2(J,L) = ( D1(J,L) - DTLOCL(I,J,K) * lambda &
        &         * EM(J+1,L) * D2(J+1,L) ) &
        &       / (1.0 + DTLOCL(I,J,K) * lambda * EA(J,L))
      ENDDO
    ENDIF
    ENDDO
    ! RI反転 -----------------------------------------------------------
    DO J = JS, JE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      DO L = LS, 5
        D1(J,L) = RR(J,L,1) * D2(J,1) + RR(J,L,2) * D2(J,2) &
        &       + RR(J,L,3) * D2(J,3) + RR(J,L,4) * D2(J,4) &
        &       + RR(J,L,5) * D2(J,5)
      ENDDO
      DO L = 6, LE
        D1(J,L) = D2(J,L)
      ENDDO
    ENDIF
    ENDDO
    ! 代入 -------------------------------------------------------------
    DO L = LS, LE
    DO J = JS, JE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      DQH(I,J,K,L) = D1(J,L)
    ENDIF
    ENDDO
    ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ! メモリ解放 =========================================================
  DEALLOCATE(D1, D2, EE, EA, EP, EM, RR, RI)
  ! ZE方向 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! メモリ確保 =========================================================
  ALLOCATE(D1(KS: KE, LS: LE))
  ALLOCATE(D2(KS: KE, LS: LE))
  ALLOCATE(EE(KS: KE, LS: LE))
  ALLOCATE(EA(KS: KE, LS: LE))
  ALLOCATE(EP(KS: KE, LS: LE))
  ALLOCATE(EM(KS: KE, LS: LE))
  ALLOCATE(RR(KS: KE, LS:  5, LS:  5))
  ALLOCATE(RI(KS: KE, LS:  5, LS:  5))
  ! 行列反転処理 =======================================================
  !$OMP DO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    ! 代入 -------------------------------------------------------------
    DO L = LS, LE
    DO K = KS, KE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      D1(K,L) = DQH(I,J,K,L)
    ENDIF
    ENDDO
    ENDDO
    ! 対角化 -----------------------------------------------------------
    DO K = KS, KE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      ! 命名
      UU   = U(I,J,K)
      VV   = V(I,J,K)
      WW   = W(I,J,K)
      C    = SQRT(GAMMA * RG * T(I,J,K))
      THI  = ZEX(I,J,K) * UU + ZEY(I,J,K) * VV + ZEZ(I,J,K) * WW
      PHI  = 0.5 * (GAMMA - 1.0) * (UU**2 + VV**2 + WW**2)
      BETA = 1.0 / (2.0 * C**2)
      TIL  = SQRT(ZEX(I,J,K)**2 + ZEY(I,J,K)**2 + ZEZ(I,J,K)**2)
      AKXT = ZEX(I,J,K) / TIL
      AKYT = ZEY(I,J,K) / TIL
      AKZT = ZEZ(I,J,K) / TIL
      THIT = AKXT * UU + AKYT * VV + AKZT * WW
      CT   = C * TIL
      VIS  = 2.0 * AMU(I,J,K) / RHO(I,J,K) * TIL**2
      ! RR
      RR(K,1,1) = AKXT
      RR(K,2,1) = AKXT * UU
      RR(K,3,1) = AKXT * VV + AKZT
      RR(K,4,1) = AKXT * WW - AKYT
      RR(K,5,1) = AKXT * PHI / (GAMMA - 1.0) + AKZT * VV - AKYT * WW
      RR(K,1,2) = AKYT
      RR(K,2,2) = AKYT * UU - AKZT
      RR(K,3,2) = AKYT * VV
      RR(K,4,2) = AKYT * WW + AKXT
      RR(K,5,2) = AKYT * PHI / (GAMMA - 1.0) + AKXT * WW - AKZT * UU
      RR(K,1,3) = AKZT
      RR(K,2,3) = AKZT * UU + AKYT
      RR(K,3,3) = AKZT * VV - AKXT
      RR(K,4,3) = AKZT * WW
      RR(K,5,3) = AKZT * PHI / (GAMMA - 1.0) + AKYT * UU - AKXT * VV
      RR(K,1,4) = 1.0
      RR(K,2,4) = UU + AKXT * C
      RR(K,3,4) = VV + AKYT * C
      RR(K,4,4) = WW + AKZT * C
      RR(K,5,4) = (PHI + C**2) / (GAMMA - 1.0) + C * THIT
      RR(K,1,5) = 1.0
      RR(K,2,5) = UU - AKXT * C
      RR(K,3,5) = VV - AKYT * C
      RR(K,4,5) = WW - AKZT * C
      RR(K,5,5) = (PHI + C**2) / (GAMMA - 1.0) - C * THIT
      ! RI
      RI(K,1,1) = AKXT - AKXT * PHI / C**2 + AKYT * WW - AKZT * VV
      RI(K,2,1) = AKYT - AKYT * PHI / C**2 + AKZT * UU - AKXT * WW
      RI(K,3,1) = AKZT - AKZT * PHI / C**2 + AKXT * VV - AKYT * UU
      RI(K,4,1) = BETA * (PHI - C * THIT)
      RI(K,5,1) = BETA * (PHI + C * THIT)
      RI(K,1,2) = (GAMMA - 1.0) * AKXT * UU / C**2
      RI(K,2,2) = (GAMMA - 1.0) * AKYT * UU / C**2 - AKZT
      RI(K,3,2) = (GAMMA - 1.0) * AKZT * UU / C**2 + AKYT
      RI(K,4,2) = BETA * (AKXT * C - (GAMMA - 1.0) * UU)
      RI(K,5,2) =-BETA * (AKXT * C + (GAMMA - 1.0) * UU)
      RI(K,1,3) = (GAMMA - 1.0) * AKXT * VV / C**2 + AKZT
      RI(K,2,3) = (GAMMA - 1.0) * AKYT * VV / C**2
      RI(K,3,3) = (GAMMA - 1.0) * AKZT * VV / C**2 - AKXT
      RI(K,4,3) = BETA * (AKYT * C - (GAMMA - 1.0) * VV)
      RI(K,5,3) =-BETA * (AKYT * C + (GAMMA - 1.0) * VV)
      RI(K,1,4) = (GAMMA - 1.0) * AKXT * WW / C**2 - AKYT
      RI(K,2,4) = (GAMMA - 1.0) * AKYT * WW / C**2 + AKXT
      RI(K,3,4) = (GAMMA - 1.0) * AKZT * WW / C**2
      RI(K,4,4) = BETA * (AKZT * C - (GAMMA - 1.0) * WW)
      RI(K,5,4) =-BETA * (AKZT * C + (GAMMA - 1.0) * WW)
      RI(K,1,5) =-(GAMMA - 1.0) * AKXT / C**2
      RI(K,2,5) =-(GAMMA - 1.0) * AKYT / C**2
      RI(K,3,5) =-(GAMMA - 1.0) * AKZT / C**2
      RI(K,4,5) = BETA * (GAMMA - 1.0)
      RI(K,5,5) = BETA * (GAMMA - 1.0)
      ! EE
      EE(K,1) = THI
      EE(K,2) = THI
      EE(K,3) = THI
      EE(K,4) = THI + CT
      EE(K,5) = THI - CT
      DO L = 6, LE
        EE(K,L) = THI
      ENDDO
      ! EA
      DO L = LS, LE
        EA(K,L) = ABS(EE(K,L)) + VIS ! + DEL * CT
      ENDDO
      ! EP
      DO L = LS, LE
        EP(K,L) = 0.5 * (EE(K,L) + EA(K,L))
      ENDDO
      ! EM
      DO L = LS, LE
        EM(K,L) = 0.5 * (EE(K,L) - EA(K,L))
      ENDDO
    ENDIF
    ENDDO
    ! RR反転 -----------------------------------------------------------
    DO K = KS, KE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      DO L = LS, 5
        D2(K,L) = RI(K,L,1) * D1(K,1) + RI(K,L,2) * D1(K,2) &
        &       + RI(K,L,3) * D1(K,3) + RI(K,L,4) * D1(K,4) &
        &       + RI(K,L,5) * D1(K,5)
      ENDDO
      DO L = 6, LE
        D2(K,L) = D1(K,L)
      ENDDO
    ENDIF
    ENDDO
    ! 前進スイープ -----------------------------------------------------
    D1(:,:) = 0.0
    DO K = KS + 1, KE - 1, + 1
    IF(QH0(I,J,K-1,1) .GT. 0.0 .AND. QH0(I,J,K,1) .GT. 0.0) THEN
      DO L = LS, LE
        D1(K,L) = ( D2(K,L) + DTLOCL(I,J,K) * lambda &
        &         * EP(K-1,L) * D1(K-1,L) ) &
        &       / (1.0 + DTLOCL(I,J,K) * lambda * EA(K,L))
      ENDDO
    ENDIF
    ENDDO
    ! スカラー演算 -----------------------------------------------------
    DO K = KS, KE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      DO L = LS, LE
        D1(K,L) = D1(K,L) * (1.0 + DTLOCL(I,J,K) * lambda * EA(K,L))
      ENDDO
    ENDIF
    ENDDO
    ! 後退スイープ -----------------------------------------------------
    D2(:,:) = 0.0
    DO K = KE - 1, KS + 1, - 1
    IF(QH0(I,J,K,1) .GT. 0.0 .AND. QH0(I,J,K+1,1) .GT. 0.0) THEN
      DO L = LS, LE
        D2(K,L) = ( D1(K,L) - DTLOCL(I,J,K) * lambda &
        &         * EM(K+1,L) * D2(K+1,L) ) &
        &       / (1.0 + DTLOCL(I,J,K) * lambda * EA(K,L))
      ENDDO
    ENDIF
    ENDDO
    ! RI反転 -----------------------------------------------------------
    DO K = KS, KE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      DO L = LS, 5
        D1(K,L) = RR(K,L,1) * D2(K,1) + RR(K,L,2) * D2(K,2) &
        &       + RR(K,L,3) * D2(K,3) + RR(K,L,4) * D2(K,4) &
        &       + RR(K,L,5) * D2(K,5)
      ENDDO
      DO L = 6, LE
        D1(K,L) = D2(K,L)
      ENDDO
    ENDIF
    ENDDO
    ! 代入 -------------------------------------------------------------
    DO L = LS, LE
    DO K = KS, KE
    IF(QH0(I,J,K,1) .GT. 0.0) THEN
      DQH(I,J,K,L) = D1(K,L)
    ENDIF
    ENDDO
    ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ! メモリ解放 =========================================================
  DEALLOCATE(D1, D2, EE, EA, EP, EM, RR, RI)
  ! デルタ形式から流束関数の時間を進める +++++++++++++++++++++++++++++++
  DO L = LS, LE
  !$OMP DO
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(QH0(I,J,K,1) .GT. 0.0) THEN
    QH1(I,J,K,L) = QH1(I,J,K,L) + DQH(I,J,K,L)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  ! メモリ解放 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP END PARALLEL
  DEALLOCATE(DQH)
  ! 処理終了 ***********************************************************
  RETURN
! 定義終了 *************************************************************
END SUBROUTINE LUADI3D
!***********************************************************************
!**** 対流項計算                                                    ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分                               ****
!****                一次精度風上差分                               ****
!****                二次精度TVDスキーム (Yee-Harten)               ****
!****                三次精度TVDスキーム (Chakravarthy-Osher)       ****
!****                四次精度TVDスキーム (Yamamoto-Daiguji)         ****
!***********************************************************************
SUBROUTINE Convection2D( &
&            Order, DEL, RG, GAMMA, &
&            IS, IE, JS, JE, LS, LE, &
&            XIX, XIY, ETX, ETY, AJA, &
&            QH, U, V, P, T, &
&            DQC &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: b1 = 4.0, b2 = 2.0
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: Order
  REAL,    INTENT(IN)  :: DEL
  REAL,    INTENT(IN)  :: RG, GAMMA
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS: IE, JS: JE), XIY(IS: IE, JS: JE), &
  &                       ETX(IS: IE, JS: JE), ETY(IS: IE, JS: JE), &
  &                       AJA(IS: IE, JS: JE)
  REAL,    INTENT(IN)  :: QH(IS: IE, JS: JE, LS: LE)
  REAL,    INTENT(IN)  :: U(IS: IE, JS: JE), V(IS: IE, JS: JE), &
  &                       P(IS: IE, JS: JE), T(IS: IE, JS: JE)
  REAL,    INTENT(OUT) :: DQC(IS: IE, JS: JE, LS: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  REAL, ALLOCATABLE :: ET(:, :, :), FT(:, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 4) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(ET(IS: IE, JS: JE, LS: LE), FT(IS: IE, JS: JE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ET (:, :, :) = 0.0
  FT (:, :, :) = 0.0
  DQC(:, :, :) = 0.0
  ! 各方向の空間差分 +++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL TVDXI
  CALL TVDET
  ! 対流項ベクトルの空間差分 +++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, L)
  DO L = LS, LE
  !$OMP DO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (QH(I  ,J-1,1) .GT. 0.0) .AND. (QH(I-1,J  ,1) .GT. 0.0) .AND. &
  &   (QH(I  ,J  ,1) .GT. 0.0) .AND. &
  &   (QH(I+1,J  ,1) .GT. 0.0) .AND. (QH(I  ,J+1,1) .GT. 0.0) &
  & ) THEN
    DQC(I,J,L) = ET(I-1,J  ,L) - ET(I,J,L) &
    &          + FT(I  ,J-1,L) - FT(I,J,L)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** TVDスキーム(XI方向)                                           ****
!***********************************************************************
SUBROUTINE TVDXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  REAL    :: AKX, AKY
  REAL    :: RWL, RWR, RWD, UM, VM, HM, CM
  REAL    :: THI, PHI, BETA, TIL, AKXT, AKYT, THIT, PH
  REAL    :: FSIGMA, FGAMMA, ABSA, AP, AM
  REAL    :: D3FP1, D3FP2, D3FP3, D3FM1, D3FM2, D3FM3
  REAL, ALLOCATABLE :: F(:, :), EIG(:, :)
  REAL, ALLOCATABLE :: RR(:, :, :), RI(:, :, :), EIGM(:, :), DELTA(:, :)
  REAL, ALLOCATABLE :: DU(:, :), DUR(:, :)
  REAL, ALLOCATABLE :: GG(:, :), PHIM(:, :)
  REAL, ALLOCATABLE :: DFP(:, :), DFM(:, :)
  REAL, ALLOCATABLE :: D3FP(:, :), D3FM(:, :)
  REAL, ALLOCATABLE :: DFBP(:, :), DFBM(:, :)
  REAL, ALLOCATABLE :: DFTP(:, :), DFTM(:, :)
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE( &
  !$  I, J, L, &
  !$  AKX, AKY, RWL, RWR, RWD, UM, VM, HM, CM, &
  !$  THI, PHI, BETA, TIL, AKXT, AKYT, THIT, PH, &
  !$  FSIGMA, FGAMMA, ABSA, AP, AM, &
  !$  D3FP1, D3FP2, D3FP3, D3FM1, D3FM2, D3FM3, &
  !$  F, EIG, RR, RI, EIGM, DELTA, DU, DUR, &
  !$  GG, PHIM, DFP, DFM, D3FP, D3FM, &
  !$  DFBP, DFBM, DFTP, DFTM &
  !$)
  ALLOCATE(F(IS: IE, LS: LE), EIG(IS: IE, LS: LE))
  ALLOCATE(RR(IS: IE - 1, LS: 4, LS: 4))
  ALLOCATE(RI(IS: IE - 1, LS: 4, LS: 4))
  ALLOCATE(EIGM(IS: IE - 1, LS: LE))
  ALLOCATE(DELTA(IS: IE - 1, LS: LE))
  ALLOCATE(DU(IS: IE - 1, LS: LE), DUR(IS: IE - 1, LS: LE))
  ALLOCATE(GG(IS: IE, LS: LE), PHIM(IS: IE - 1, LS: LE))
  ALLOCATE(DFP (IS    : IE - 1, LS: LE), DFM (IS    : IE - 1, LS: LE))
  ALLOCATE(D3FP(IS + 1: IE - 2, LS: LE), D3FM(IS + 1: IE - 2, LS: LE))
  ALLOCATE(DFBP(IS    : IE - 3, LS: LE), DFBM(IS + 1: IE - 2, LS: LE))
  ALLOCATE(DFTP(IS + 1: IE - 2, LS: LE), DFTM(IS + 2: IE - 1, LS: LE))
  !$OMP DO
  DO J = JS + 1, JE - 1
    ! 二次精度中心差分 +++++++++++++++++++++++++++++++++++++++++++++++++
    DO I = IS, IE
    IF(QH(I,J,1) .GT. 0.0) THEN
      THI = XIX(I,J) * U(I,J) + XIY(I,J) * V(I,J)
      TIL = SQRT(XIX(I,J)**2 + XIY(I,J)**2)
      CM  = SQRT(GAMMA * RG * T(I,J))
      PH  = P(I,J) / AJA(I,J)
      F(I,1) = QH(I,J,1) * THI
      F(I,2) = QH(I,J,2) * THI + PH * XIX(I,J)
      F(I,3) = QH(I,J,3) * THI + PH * XIY(I,J)
      F(I,4) = QH(I,J,4) * THI + PH * THI
      DO L = 5, LE
        F(I,L) = QH(I,J,L) * THI
      ENDDO
      EIG(I,1) = THI
      EIG(I,2) = THI
      EIG(I,3) = THI + CM * TIL
      EIG(I,4) = THI - CM * TIL
      DO L = 5, LE
        EIG(I,L) = THI
      ENDDO
    ENDIF
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO I = IS, IE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I+1,J,1) .GT. 0.0) THEN
      ET(I,J,L) = 0.5 * (F(I,L) + F(I+1,L))
    ENDIF
    ENDDO
    ENDDO
    IF(Order .EQ. 0) CYCLE
    ! 対角化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO I = IS, IE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I+1,J,1) .GT. 0.0) THEN
      AKX  = 0.5 * (XIX(I,J) + XIX(I+1,J))
      AKY  = 0.5 * (XIY(I,J) + XIY(I+1,J))
      ! Roe 平均
      RWD  = 1.0 / ( SQRT(QH(I  ,J,1) * AJA(I  ,J)) &
      &            + SQRT(QH(I+1,J,1) * AJA(I+1,J)) )
      RWL  = SQRT(QH(I  ,J,1) * AJA(I  ,J)) * RWD
      RWR  = SQRT(QH(I+1,J,1) * AJA(I+1,J)) * RWD
      UM   = RWL * U(I,J) + RWR * U(I+1,J)
      VM   = RWL * V(I,J) + RWR * V(I+1,J)
      HM   = RWL * ( &
      &        0.5 * (U(I  ,J)**2 + V(I  ,J)**2) &
      &      + GAMMA * RG / (GAMMA - 1.0) * T(I  ,J) &
      &    ) &
      &    + RWR * ( &
      &        0.5 * (U(I+1,J)**2 + V(I+1,J)**2) &
      &      + GAMMA * RG / (GAMMA - 1.0) * T(I+1,J) &
      &    )
      CM   = SQRT((GAMMA - 1.0) * (HM - 0.5 * (UM**2 + VM**2)))
      ! 対角化
      THI  = AKX * UM + AKY * VM
      PHI  = 0.5 * (GAMMA - 1.0) * (UM**2 + VM**2)
      BETA = 1.0 / (2.0 * CM**2)
      TIL  = SQRT(AKX**2 + AKY**2)
      AKXT = AKX / TIL
      AKYT = AKY / TIL
      THIT = AKXT * UM + AKYT * VM
      RR(I,1,1) = 1.0
      RR(I,2,1) = UM
      RR(I,3,1) = VM
      RR(I,4,1) = PHI / (GAMMA - 1.0)
      RR(I,1,2) = 0.0
      RR(I,2,2) = AKYT
      RR(I,3,2) =-AKXT
      RR(I,4,2) = AKYT * UM - AKXT * VM
      RR(I,1,3) = 1.0
      RR(I,2,3) = UM + AKXT * CM
      RR(I,3,3) = VM + AKYT * CM
      RR(I,4,3) = (PHI + CM**2) / (GAMMA - 1.0) + CM * THIT
      RR(I,1,4) = 1.0
      RR(I,2,4) = UM - AKXT * CM
      RR(I,3,4) = VM - AKYT * CM
      RR(I,4,4) = (PHI + CM**2) / (GAMMA - 1.0) - CM * THIT
      RI(I,1,1) = 1.0 - PHI / CM**2
      RI(I,2,1) =-AKYT * UM + AKXT * VM
      RI(I,3,1) = BETA * (PHI - CM * THIT)
      RI(I,4,1) = BETA * (PHI + CM * THIT)
      RI(I,1,2) = (GAMMA - 1.0) * UM / CM**2
      RI(I,2,2) = AKYT
      RI(I,3,2) = BETA * (AKXT * CM - (GAMMA - 1.0) * UM)
      RI(I,4,2) =-BETA * (AKXT * CM + (GAMMA - 1.0) * UM)
      RI(I,1,3) = (GAMMA - 1.0) * VM / CM**2
      RI(I,2,3) =-AKXT
      RI(I,3,3) = BETA * (AKYT * CM - (GAMMA - 1.0) * VM)
      RI(I,4,3) =-BETA * (AKYT * CM + (GAMMA - 1.0) * VM)
      RI(I,1,4) =-(GAMMA - 1.0) / CM**2
      RI(I,2,4) = 0.0
      RI(I,3,4) = BETA * (GAMMA - 1.0)
      RI(I,4,4) = BETA * (GAMMA - 1.0)
      EIGM(I,1) = THI
      EIGM(I,2) = THI
      EIGM(I,3) = THI + CM * TIL
      EIGM(I,4) = THI - CM * TIL
      DO L = 5, LE
        EIGM(I,L) = THI
      ENDDO
      DO L = LS, LE
        DELTA(I,L) = MAX( &
        &              0.0, &
        &              EIGM(I,L) - EIG(I,L), EIG(I+1,L) - EIGM(I,L) &
        &            ) + DEL * CM * TIL
      ENDDO
    ENDIF
    ENDDO
    ! Delta u
    DO L = LS, LE
    DO I = IS, IE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I+1,J,1) .GT. 0.0) THEN
      DU(I,L) = - QH(I  ,J,L) * AJA(I  ,J) &
      &         + QH(I+1,J,L) * AJA(I+1,J)
      DU(I,L) = DU(I,L) * 2.0 / (AJA(I,J) + AJA(I+1,J))
    ENDIF
    ENDDO
    ENDDO
    ! R^-1 Delta u
    DO L = LS, LE
    DO I = IS, IE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I+1,J,1) .GT. 0.0) THEN
      IF(L .LE. 4) THEN
        DUR(I,L) = RI(I,L,1) * DU(I,1) + RI(I,L,2) * DU(I,2) &
        &        + RI(I,L,3) * DU(I,3) + RI(I,L,4) * DU(I,4)
      ELSE
        DUR(I,L) = DU(I,L)
      ENDIF
    ENDIF
    ENDDO
    ENDDO
    ! 二次精度TVDスキーム ++++++++++++++++++++++++++++++++++++++++++++++
    IF(Order .EQ. 2) THEN
      ! g
      DO L = LS, LE
      I = IS
      IF(QH(I,J,1) .GT. 0.0 .AND. QH(I+1,J,1) .GT. 0.0) THEN
        GG(I,L) = DUR(I,L)
      ENDIF
      DO I = IS + 1, IE - 1
      IF(QH(I,J,1) .LE. 0.0) CYCLE
      IF(QH(I-1,J,1) .GT. 0.0 .AND. QH(I+1,J,1) .GT. 0.0) THEN
        GG(I,L) = minmod2(DUR(I,L), DUR(I-1,L))
      ELSEIF(QH(I-1,J,1) .GT. 0.0) THEN
        GG(I,L) = DUR(I-1,L)
      ELSEIF(QH(I+1,J,1) .GT. 0.0) THEN
        GG(I,L) = DUR(I,L)
      ENDIF
      ENDDO
      I = IE
      IF(QH(I-1,J,1) .GT. 0.0 .AND. QH(I,J,1) .GT. 0.0) THEN
        GG(I,L) = DUR(I-1,L)
      ENDIF
      ENDDO
      ! phi
      DO L = LS, LE
      DO I = IS, IE - 1
      IF(QH(I,J,1) .GT. 0.0 .AND. QH(I+1,J,1) .GT. 0.0) THEN
        FSIGMA = Func_sigma(Func_psi(EIGM(I,L), DELTA(I,L)))
        FGAMMA = Func_gamma(FSIGMA, GG(I,L), GG(I+1,L), DUR(I,L))
        PHIM(I,L) = FSIGMA * (GG(I,L) + GG(I+1,L)) &
        &         - Func_psi(EIGM(I,L) + FGAMMA, DELTA(I,L)) * DUR(I,L)
      ENDIF
      ENDDO
      ENDDO
      ! 数値流束
      DO L = LS, LE
      DO I = IS, IE - 1
      IF(QH(I,J,1) .GT. 0.0 .AND. QH(I+1,J,1) .GT. 0.0) THEN
        IF(L .LE. 4) THEN
          ET(I,J,L) = ET(I,J,L) + 0.5 * ( RR(I,L,1) * PHIM(I,1) &
          &                             + RR(I,L,2) * PHIM(I,2) &
          &                             + RR(I,L,3) * PHIM(I,3) &
          &                             + RR(I,L,4) * PHIM(I,4) )
        ELSE
          ET(I,J,L) = ET(I,J,L) + 0.5 * PHIM(I,L)
        ENDIF
      ENDIF
      ENDDO
      ENDDO
      CYCLE
    ENDIF
    ! 一次精度風上差分 +++++++++++++++++++++++++++++++++++++++++++++++++
    ! Delta f^+-
    DO L = LS, LE
    DO I = IS, IE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I+1,J,1) .GT. 0.0) THEN
      ABSA  = Func_psi(EIGM(I,L), DELTA(I,L))
      AP    = 0.5 * (EIGM(I,L) + ABSA)
      AM    = 0.5 * (EIGM(I,L) - ABSA)
      DFP(I,L) = AP * DUR(I,L)
      DFM(I,L) = AM * DUR(I,L)
    ENDIF
    ENDDO
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO I = IS, IE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I+1,J,1) .GT. 0.0) THEN
      IF(L .LE. 4) THEN
        ET(I,J,L) = ET(I,J,L) - 0.5 * ( RR(I,L,1) * DFP(I,1) &
        &                             + RR(I,L,2) * DFP(I,2) &
        &                             + RR(I,L,3) * DFP(I,3) &
        &                             + RR(I,L,4) * DFP(I,4) &
        &                             - RR(I,L,1) * DFM(I,1) &
        &                             - RR(I,L,2) * DFM(I,2) &
        &                             - RR(I,L,3) * DFM(I,3) &
        &                             - RR(I,L,4) * DFM(I,4) )
      ELSE
        ET(I,J,L) = ET(I,J,L) - 0.5 * (DFP(I,L) - DFM(I,L))
      ENDIF
    ENDIF
    ENDDO
    ENDDO
    IF(Order .EQ. 1) CYCLE
    ! 三次&四次精度TVDスキーム +++++++++++++++++++++++++++++++++++++++++
    IF(Order .EQ. 4) THEN
      ! Delta^3 bar(f)^+-
      DO L = LS, LE
      DO I = IS + 1, IE - 2
      IF( QH(I-1,J,1) .GT. 0.0 .AND. QH(I  ,J,1) .GT. 0.0 .AND. &
      &   QH(I+1,J,1) .GT. 0.0 .AND. QH(I+2,J,1) .GT. 0.0 ) THEN
        D3FP1 = minmod3(DFP(I-1,L), b2 * DFP(I  ,L), b2 * DFP(I+1,L))
        D3FP2 = minmod3(DFP(I  ,L), b2 * DFP(I+1,L), b2 * DFP(I-1,L))
        D3FP3 = minmod3(DFP(I+1,L), b2 * DFP(I-1,L), b2 * DFP(I  ,L))
        D3FM1 = minmod3(DFM(I-1,L), b2 * DFM(I  ,L), b2 * DFM(I+1,L))
        D3FM2 = minmod3(DFM(I  ,L), b2 * DFM(I+1,L), b2 * DFM(I-1,L))
        D3FM3 = minmod3(DFM(I+1,L), b2 * DFM(I-1,L), b2 * DFM(I  ,L))
        D3FP(I,L) = D3FP1 - 2.0 * D3FP2 + D3FP3
        D3FM(I,L) = D3FM1 - 2.0 * D3FM2 + D3FM3
      ENDIF
      ENDDO
      ENDDO
      ! Delta^asterisk f^+-
      DO L = LS, LE
      DO I = IS + 1, IE - 2
      IF( QH(I-1,J,1) .GT. 0.0 .AND. QH(I  ,J,1) .GT. 0.0 .AND. &
      &   QH(I+1,J,1) .GT. 0.0 .AND. QH(I+2,J,1) .GT. 0.0 ) THEN
        DFP(I,L) = DFP(I,L) - D3FP(I,L) / 6.0
        DFM(I,L) = DFM(I,L) - D3FM(I,L) / 6.0
      ENDIF
      ENDDO
      ENDDO
    ENDIF
    ! Delta bar(f)^+-, Delta tilde(f)^+-
    DO L = LS, LE
    DO I = IS + 1, IE - 2
    IF( QH(I-1,J,1) .GT. 0.0 .AND. QH(I  ,J,1) .GT. 0.0 .AND. &
    &   QH(I+1,J,1) .GT. 0.0 .AND. QH(I+2,J,1) .GT. 0.0 ) THEN
      DFBP(I-1,L) = minmod2(DFP(I-1,L), b1 * DFP(I  ,L))
      DFBM(I  ,L) = minmod2(DFM(I  ,L), b1 * DFM(I+1,L))
      DFTP(I  ,L) = minmod2(DFP(I  ,L), b1 * DFP(I-1,L))
      DFTM(I+1,L) = minmod2(DFM(I+1,L), b1 * DFM(I  ,L))
    ENDIF
    ENDDO
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO I = IS + 1, IE - 2
    IF( QH(I-1,J,1) .GT. 0.0 .AND. QH(I  ,J,1) .GT. 0.0 .AND. &
    &   QH(I+1,J,1) .GT. 0.0 .AND. QH(I+2,J,1) .GT. 0.0 ) THEN
      IF(L .LE. 4) THEN
        ET(I,J,L) = ET(I,J,L) &
        &         + ( 1.0 * RR(I-1,L,1) * DFBP(I-1,1) &
        &           + 1.0 * RR(I-1,L,2) * DFBP(I-1,2) &
        &           + 1.0 * RR(I-1,L,3) * DFBP(I-1,3) &
        &           + 1.0 * RR(I-1,L,4) * DFBP(I-1,4) &
        &           + 2.0 * RR(I  ,L,1) * DFTP(I  ,1) &
        &           + 2.0 * RR(I  ,L,2) * DFTP(I  ,2) &
        &           + 2.0 * RR(I  ,L,3) * DFTP(I  ,3) &
        &           + 2.0 * RR(I  ,L,4) * DFTP(I  ,4) &
        &           - 2.0 * RR(I  ,L,1) * DFBM(I  ,1) &
        &           - 2.0 * RR(I  ,L,2) * DFBM(I  ,2) &
        &           - 2.0 * RR(I  ,L,3) * DFBM(I  ,3) &
        &           - 2.0 * RR(I  ,L,4) * DFBM(I  ,4) &
        &           - 1.0 * RR(I+1,L,1) * DFTM(I+1,1) &
        &           - 1.0 * RR(I+1,L,2) * DFTM(I+1,2) &
        &           - 1.0 * RR(I+1,L,3) * DFTM(I+1,3) &
        &           - 1.0 * RR(I+1,L,4) * DFTM(I+1,4) &
        &         ) / 6.0
      ELSE
        ET(I,J,L) = ET(I,J,L) &
        &         + ( 1.0 * DFBP(I-1,L) + 2.0 * DFTP(I  ,L) &
        &           - 2.0 * DFBM(I  ,L) - 1.0 * DFTM(I+1,L) ) / 6.0
      ENDIF
    ENDIF
    ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE TVDXI
!***********************************************************************
!**** TVDスキーム(ET方向)                                           ****
!***********************************************************************
SUBROUTINE TVDET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  REAL    :: AKX, AKY
  REAL    :: RWL, RWR, RWD, UM, VM, HM, CM
  REAL    :: THI, PHI, BETA, TIL, AKXT, AKYT, THIT, PH
  REAL    :: FSIGMA, FGAMMA, ABSA, AP, AM
  REAL    :: D3FP1, D3FP2, D3FP3, D3FM1, D3FM2, D3FM3
  REAL, ALLOCATABLE :: F(:, :), EIG(:, :)
  REAL, ALLOCATABLE :: RR(:, :, :), RI(:, :, :), EIGM(:, :), DELTA(:, :)
  REAL, ALLOCATABLE :: DU(:, :), DUR(:, :)
  REAL, ALLOCATABLE :: GG(:, :), PHIM(:, :)
  REAL, ALLOCATABLE :: DFP(:, :), DFM(:, :)
  REAL, ALLOCATABLE :: D3FP(:, :), D3FM(:, :)
  REAL, ALLOCATABLE :: DFBP(:, :), DFBM(:, :)
  REAL, ALLOCATABLE :: DFTP(:, :), DFTM(:, :)
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE( &
  !$  I, J, L, &
  !$  AKX, AKY, RWL, RWR, RWD, UM, VM, HM, CM, &
  !$  THI, PHI, BETA, TIL, AKXT, AKYT, THIT, PH, &
  !$  FSIGMA, FGAMMA, ABSA, AP, AM, &
  !$  D3FP1, D3FP2, D3FP3, D3FM1, D3FM2, D3FM3, &
  !$  F, EIG, RR, RI, EIGM, DELTA, DU, DUR, &
  !$  GG, PHIM, DFP, DFM, D3FP, D3FM, &
  !$  DFBP, DFBM, DFTP, DFTM &
  !$)
  ALLOCATE(F(JS: JE, LS: LE), EIG(JS: JE, LS: LE))
  ALLOCATE(RR(JS: JE - 1, LS: 4, LS: 4))
  ALLOCATE(RI(JS: JE - 1, LS: 4, LS: 4))
  ALLOCATE(EIGM(JS: JE - 1, LS: LE))
  ALLOCATE(DELTA(JS: JE - 1, LS: LE))
  ALLOCATE(DU(JS: JE - 1, LS: LE), DUR(JS: JE - 1, LS: LE))
  ALLOCATE(GG(JS: JE, LS: LE), PHIM(JS: JE - 1, LS: LE))
  ALLOCATE(DFP (JS    : JE - 1, LS: LE), DFM (JS    : JE - 1, LS: LE))
  ALLOCATE(D3FP(JS + 1: JE - 2, LS: LE), D3FM(JS + 1: JE - 2, LS: LE))
  ALLOCATE(DFBP(JS    : JE - 3, LS: LE), DFBM(JS + 1: JE - 2, LS: LE))
  ALLOCATE(DFTP(JS + 1: JE - 2, LS: LE), DFTM(JS + 2: JE - 1, LS: LE))
  !$OMP DO
  DO I = IS + 1, IE - 1
    ! 二次精度中心差分 +++++++++++++++++++++++++++++++++++++++++++++++++
    DO J = JS, JE
    IF(QH(I,J,1) .GT. 0.0) THEN
      THI = ETX(I,J) * U(I,J) + ETY(I,J) * V(I,J)
      TIL = SQRT(ETX(I,J)**2 + ETY(I,J)**2)
      CM  = SQRT(GAMMA * RG * T(I,J))
      PH  = P(I,J) / AJA(I,J)
      F(J,1) = QH(I,J,1) * THI
      F(J,2) = QH(I,J,2) * THI + PH * ETX(I,J)
      F(J,3) = QH(I,J,3) * THI + PH * ETY(I,J)
      F(J,4) = QH(I,J,4) * THI + PH * THI
      DO L = 5, LE
        F(J,L) = QH(I,J,L) * THI
      ENDDO
      EIG(J,1) = THI
      EIG(J,2) = THI
      EIG(J,3) = THI + CM * TIL
      EIG(J,4) = THI - CM * TIL
      DO L = 5, LE
        EIG(J,L) = THI
      ENDDO
    ENDIF
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO J = JS, JE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J+1,1) .GT. 0.0) THEN
      FT(I,J,L) = 0.5 * (F(J,L) + F(J+1,L))
    ENDIF
    ENDDO
    ENDDO
    IF(Order .EQ. 0) CYCLE
    ! 対角化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO J = JS, JE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J+1,1) .GT. 0.0) THEN
      AKX  = 0.5 * (ETX(I,J) + ETX(I,J+1))
      AKY  = 0.5 * (ETY(I,J) + ETY(I,J+1))
      ! Roe 平均
      RWD  = 1.0 / ( SQRT(QH(I,J  ,1) * AJA(I,J  )) &
      &            + SQRT(QH(I,J+1,1) * AJA(I,J+1)) )
      RWL  = SQRT(QH(I,J  ,1) * AJA(I,J  )) * RWD
      RWR  = SQRT(QH(I,J+1,1) * AJA(I,J+1)) * RWD
      UM   = RWL * U(I,J) + RWR * U(I,J+1)
      VM   = RWL * V(I,J) + RWR * V(I,J+1)
      HM   = RWL * ( &
      &        0.5 * (U(I,J  )**2 + V(I,J  )**2) &
      &      + GAMMA * RG / (GAMMA - 1.0) * T(I,J  ) &
      &    ) &
      &    + RWR * ( &
      &        0.5 * (U(I,J+1)**2 + V(I,J+1)**2) &
      &      + GAMMA * RG / (GAMMA - 1.0) * T(I,J+1) &
      &    )
      CM   = SQRT((GAMMA - 1.0) * (HM - 0.5 * (UM**2 + VM**2)))
      ! 対角化
      THI  = AKX * UM + AKY * VM
      PHI  = 0.5 * (GAMMA - 1.0) * (UM**2 + VM**2)
      BETA = 1.0 / (2.0 * CM**2)
      TIL  = SQRT(AKX**2 + AKY**2)
      AKXT = AKX / TIL
      AKYT = AKY / TIL
      THIT = AKXT * UM + AKYT * VM
      RR(J,1,1) = 1.0
      RR(J,2,1) = UM
      RR(J,3,1) = VM
      RR(J,4,1) = PHI / (GAMMA - 1.0)
      RR(J,1,2) = 0.0
      RR(J,2,2) = AKYT
      RR(J,3,2) =-AKXT
      RR(J,4,2) = AKYT * UM - AKXT * VM
      RR(J,1,3) = 1.0
      RR(J,2,3) = UM + AKXT * CM
      RR(J,3,3) = VM + AKYT * CM
      RR(J,4,3) = (PHI + CM**2) / (GAMMA - 1.0) + CM * THIT
      RR(J,1,4) = 1.0
      RR(J,2,4) = UM - AKXT * CM
      RR(J,3,4) = VM - AKYT * CM
      RR(J,4,4) = (PHI + CM**2) / (GAMMA - 1.0) - CM * THIT
      RI(J,1,1) = 1.0 - PHI / CM**2
      RI(J,2,1) =-AKYT * UM + AKXT * VM
      RI(J,3,1) = BETA * (PHI - CM * THIT)
      RI(J,4,1) = BETA * (PHI + CM * THIT)
      RI(J,1,2) = (GAMMA - 1.0) * UM / CM**2
      RI(J,2,2) = AKYT
      RI(J,3,2) = BETA * (AKXT * CM - (GAMMA - 1.0) * UM)
      RI(J,4,2) =-BETA * (AKXT * CM + (GAMMA - 1.0) * UM)
      RI(J,1,3) = (GAMMA - 1.0) * VM / CM**2
      RI(J,2,3) =-AKXT
      RI(J,3,3) = BETA * (AKYT * CM - (GAMMA - 1.0) * VM)
      RI(J,4,3) =-BETA * (AKYT * CM + (GAMMA - 1.0) * VM)
      RI(J,1,4) =-(GAMMA - 1.0) / CM**2
      RI(J,2,4) = 0.0
      RI(J,3,4) = BETA * (GAMMA - 1.0)
      RI(J,4,4) = BETA * (GAMMA - 1.0)
      EIGM(J,1) = THI
      EIGM(J,2) = THI
      EIGM(J,3) = THI + CM * TIL
      EIGM(J,4) = THI - CM * TIL
      DO L = 5, LE
        EIGM(J,L) = THI
      ENDDO
      DO L = LS, LE
        DELTA(J,L) = MAX( &
        &              0.0, &
        &              EIGM(J,L) - EIG(J,L), EIG(J+1,L) - EIGM(J,L) &
        &            ) + DEL * CM * TIL
      ENDDO
    ENDIF
    ENDDO
    ! Delta u
    DO L = LS, LE
    DO J = JS, JE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J+1,1) .GT. 0.0) THEN
      DU(J,L) = - QH(I,J  ,L) * AJA(I,J  ) &
      &         + QH(I,J+1,L) * AJA(I,J+1)
      DU(J,L) = DU(J,L) * 2.0 / (AJA(I,J) + AJA(I,J+1))
    ENDIF
    ENDDO
    ENDDO
    ! R^-1 Delta u
    DO L = LS, LE
    DO J = JS, JE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J+1,1) .GT. 0.0) THEN
      IF(L .LE. 4) THEN
        DUR(J,L) = RI(J,L,1) * DU(J,1) + RI(J,L,2) * DU(J,2) &
        &        + RI(J,L,3) * DU(J,3) + RI(J,L,4) * DU(J,4)
      ELSE
        DUR(J,L) = DU(J,L)
      ENDIF
    ENDIF
    ENDDO
    ENDDO
    ! 二次精度TVDスキーム ++++++++++++++++++++++++++++++++++++++++++++++
    IF(Order .EQ. 2) THEN
      ! g
      DO L = LS, LE
      J = JS
      IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J+1,1) .GT. 0.0) THEN
        GG(J,L) = DUR(J,L)
      ENDIF
      DO J = JS + 1, JE - 1
      IF(QH(I,J,1) .LE. 0.0) CYCLE
      IF(QH(I,J-1,1) .GT. 0.0 .AND. QH(I,J+1,1) .GT. 0.0) THEN
        GG(J,L) = minmod2(DUR(J,L), DUR(J-1,L))
      ELSEIF(QH(I,J-1,1) .GT. 0.0) THEN
        GG(J,L) = DUR(J-1,L)
      ELSEIF(QH(I,J+1,1) .GT. 0.0) THEN
        GG(J,L) = DUR(J,L)
      ENDIF
      ENDDO
      J = JE
      IF(QH(I,J-1,1) .GT. 0.0 .AND. QH(I,J,1) .GT. 0.0) THEN
        GG(J,L) = DUR(J-1,L)
      ENDIF
      ENDDO
      ! phi
      DO L = LS, LE
      DO J = JS, JE - 1
      IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J+1,1) .GT. 0.0) THEN
        FSIGMA = Func_sigma(Func_psi(EIGM(J,L), DELTA(J,L)))
        FGAMMA = Func_gamma(FSIGMA, GG(J,L), GG(J+1,L), DUR(J,L))
        PHIM(J,L) = FSIGMA * (GG(J,L) + GG(J+1,L)) &
        &         - Func_psi(EIGM(J,L) + FGAMMA, DELTA(J,L)) * DUR(J,L)
      ENDIF
      ENDDO
      ENDDO
      ! 数値流束
      DO L = LS, LE
      DO J = JS, JE - 1
      IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J+1,1) .GT. 0.0) THEN
        IF(L .LE. 4) THEN
          FT(I,J,L) = FT(I,J,L) + 0.5 * ( RR(J,L,1) * PHIM(J,1) &
          &                             + RR(J,L,2) * PHIM(J,2) &
          &                             + RR(J,L,3) * PHIM(J,3) &
          &                             + RR(J,L,4) * PHIM(J,4) )
        ELSE
          FT(I,J,L) = FT(I,J,L) + 0.5 * PHIM(J,L)
        ENDIF
      ENDIF
      ENDDO
      ENDDO
      CYCLE
    ENDIF
    ! 一次精度風上差分 +++++++++++++++++++++++++++++++++++++++++++++++++
    ! Delta f^+-
    DO L = LS, LE
    DO J = JS, JE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J+1,1) .GT. 0.0) THEN
      ABSA  = Func_psi(EIGM(J,L), DELTA(J,L))
      AP    = 0.5 * (EIGM(J,L) + ABSA)
      AM    = 0.5 * (EIGM(J,L) - ABSA)
      DFP(J,L) = AP * DUR(J,L)
      DFM(J,L) = AM * DUR(J,L)
    ENDIF
    ENDDO
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO J = JS, JE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J+1,1) .GT. 0.0) THEN
      IF(L .LE. 4) THEN
        FT(I,J,L) = FT(I,J,L) - 0.5 * ( RR(J,L,1) * DFP(J,1) &
        &                             + RR(J,L,2) * DFP(J,2) &
        &                             + RR(J,L,3) * DFP(J,3) &
        &                             + RR(J,L,4) * DFP(J,4) &
        &                             - RR(J,L,1) * DFM(J,1) &
        &                             - RR(J,L,2) * DFM(J,2) &
        &                             - RR(J,L,3) * DFM(J,3) &
        &                             - RR(J,L,4) * DFM(J,4) )
      ELSE
        FT(I,J,L) = FT(I,J,L) - 0.5 * (DFP(J,L) - DFM(J,L))
      ENDIF
    ENDIF
    ENDDO
    ENDDO
    IF(Order .EQ. 1) CYCLE
    ! 三次&四次精度TVDスキーム +++++++++++++++++++++++++++++++++++++++++
    IF(Order .EQ. 4) THEN
      ! Delta^3 bar(f)^+-
      DO L = LS, LE
      DO J = JS + 1, JE - 2
      IF( QH(I,J-1,1) .GT. 0.0 .AND. QH(I,J  ,1) .GT. 0.0 .AND. &
      &   QH(I,J+1,1) .GT. 0.0 .AND. QH(I,J+2,1) .GT. 0.0 ) THEN
        D3FP1 = minmod3(DFP(J-1,L), b2 * DFP(J  ,L), b2 * DFP(J+1,L))
        D3FP2 = minmod3(DFP(J  ,L), b2 * DFP(J+1,L), b2 * DFP(J-1,L))
        D3FP3 = minmod3(DFP(J+1,L), b2 * DFP(J-1,L), b2 * DFP(J  ,L))
        D3FM1 = minmod3(DFM(J-1,L), b2 * DFM(J  ,L), b2 * DFM(J+1,L))
        D3FM2 = minmod3(DFM(J  ,L), b2 * DFM(J+1,L), b2 * DFM(J-1,L))
        D3FM3 = minmod3(DFM(J+1,L), b2 * DFM(J-1,L), b2 * DFM(J  ,L))
        D3FP(J,L) = D3FP1 - 2.0 * D3FP2 + D3FP3
        D3FM(J,L) = D3FM1 - 2.0 * D3FM2 + D3FM3
      ENDIF
      ENDDO
      ENDDO
      ! Delta^asterisk f^+-
      DO L = LS, LE
      DO J = JS + 1, JE - 2
      IF( QH(I,J-1,1) .GT. 0.0 .AND. QH(I,J  ,1) .GT. 0.0 .AND. &
      &   QH(I,J+1,1) .GT. 0.0 .AND. QH(I,J+2,1) .GT. 0.0 ) THEN
        DFP(J,L) = DFP(J,L) - D3FP(J,L) / 6.0
        DFM(J,L) = DFM(J,L) - D3FM(J,L) / 6.0
      ENDIF
      ENDDO
      ENDDO
    ENDIF
    ! Delta bar(f)^+-, Delta tilde(f)^+-
    DO L = LS, LE
    DO J = JS + 1, JE - 2
    IF( QH(I,J-1,1) .GT. 0.0 .AND. QH(I,J  ,1) .GT. 0.0 .AND. &
    &   QH(I,J+1,1) .GT. 0.0 .AND. QH(I,J+2,1) .GT. 0.0 ) THEN
      DFBP(J-1,L) = minmod2(DFP(J-1,L), b1 * DFP(J  ,L))
      DFBM(J  ,L) = minmod2(DFM(J  ,L), b1 * DFM(J+1,L))
      DFTP(J  ,L) = minmod2(DFP(J  ,L), b1 * DFP(J-1,L))
      DFTM(J+1,L) = minmod2(DFM(J+1,L), b1 * DFM(J  ,L))
    ENDIF
    ENDDO
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO J = JS + 1, JE - 2
    IF( QH(I,J-1,1) .GT. 0.0 .AND. QH(I,J  ,1) .GT. 0.0 .AND. &
    &   QH(I,J+1,1) .GT. 0.0 .AND. QH(I,J+2,1) .GT. 0.0 ) THEN
      IF(L .LE. 4) THEN
        FT(I,J,L) = FT(I,J,L) &
        &         + ( 1.0 * RR(J-1,L,1) * DFBP(J-1,1) &
        &           + 1.0 * RR(J-1,L,2) * DFBP(J-1,2) &
        &           + 1.0 * RR(J-1,L,3) * DFBP(J-1,3) &
        &           + 1.0 * RR(J-1,L,4) * DFBP(J-1,4) &
        &           + 2.0 * RR(J  ,L,1) * DFTP(J  ,1) &
        &           + 2.0 * RR(J  ,L,2) * DFTP(J  ,2) &
        &           + 2.0 * RR(J  ,L,3) * DFTP(J  ,3) &
        &           + 2.0 * RR(J  ,L,4) * DFTP(J  ,4) &
        &           - 2.0 * RR(J  ,L,1) * DFBM(J  ,1) &
        &           - 2.0 * RR(J  ,L,2) * DFBM(J  ,2) &
        &           - 2.0 * RR(J  ,L,3) * DFBM(J  ,3) &
        &           - 2.0 * RR(J  ,L,4) * DFBM(J  ,4) &
        &           - 1.0 * RR(J+1,L,1) * DFTM(J+1,1) &
        &           - 1.0 * RR(J+1,L,2) * DFTM(J+1,2) &
        &           - 1.0 * RR(J+1,L,3) * DFTM(J+1,3) &
        &           - 1.0 * RR(J+1,L,4) * DFTM(J+1,4) &
        &         ) / 6.0
      ELSE
        FT(I,J,L) = FT(I,J,L) &
        &         + ( 1.0 * DFBP(J-1,L) + 2.0 * DFTP(J  ,L) &
        &           - 2.0 * DFBM(J  ,L) - 1.0 * DFTM(J+1,L) ) / 6.0
      ENDIF
    ENDIF
    ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE TVDET
! 定義終了 *************************************************************
END SUBROUTINE Convection2D
!***********************************************************************
!**** 対流項計算                                                    ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分                               ****
!****                一次精度風上差分                               ****
!****                二次精度TVDスキーム (Yee-Harten)               ****
!****                三次精度TVDスキーム (Chakravarthy-Osher)       ****
!****                四次精度TVDスキーム (Yamamoto-Daiguji)         ****
!***********************************************************************
SUBROUTINE Convection3D( &
&            Order, DEL, RG, GAMMA, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, AJA, &
&            QH, U, V, W, P, T, &
&            DQC &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: b1 = 4.0, b2 = 2.0
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: Order
  REAL,    INTENT(IN)  :: DEL
  REAL,    INTENT(IN)  :: RG, GAMMA
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS: IE, JS: JE, KS: KE), &
  &                       XIY(IS: IE, JS: JE, KS: KE), &
  &                       XIZ(IS: IE, JS: JE, KS: KE), &
  &                       ETX(IS: IE, JS: JE, KS: KE), &
  &                       ETY(IS: IE, JS: JE, KS: KE), &
  &                       ETZ(IS: IE, JS: JE, KS: KE), &
  &                       ZEX(IS: IE, JS: JE, KS: KE), &
  &                       ZEY(IS: IE, JS: JE, KS: KE), &
  &                       ZEZ(IS: IE, JS: JE, KS: KE), &
  &                       AJA(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)  :: QH(IS: IE, JS: JE, KS: KE, LS: LE)
  REAL,    INTENT(IN)  :: U(IS: IE, JS: JE, KS: KE), &
  &                       V(IS: IE, JS: JE, KS: KE), &
  &                       W(IS: IE, JS: JE, KS: KE), &
  &                       P(IS: IE, JS: JE, KS: KE), &
  &                       T(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(OUT) :: DQC(IS: IE, JS: JE, KS: KE, LS: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  REAL, ALLOCATABLE :: ET(:, :, :, :), FT(:, :, :, :), GT(:, :, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 5) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(ET(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(FT(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(GT(IS: IE, JS: JE, KS: KE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ET (:, :, :, :) = 0.0
  FT (:, :, :, :) = 0.0
  GT (:, :, :, :) = 0.0
  DQC(:, :, :, :) = 0.0
  ! 各方向の空間差分 +++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL TVDXI
  CALL TVDET
  CALL TVDZE
  ! 対流項ベクトルの空間差分 +++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( ( QH(I  ,J  ,K-1,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J-1,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I-1,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I+1,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J+1,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J  ,K+1,1) .GT. 0.0 ) &
  & ) THEN
    DQC(I,J,K,L) = ET(I-1,J  ,K  ,L) - ET(I,J,K,L) &
    &            + FT(I  ,J-1,K  ,L) - FT(I,J,K,L) &
    &            + GT(I  ,J  ,K-1,L) - GT(I,J,K,L)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! メモリ解放 =========================================================
  DEALLOCATE(ET,FT,GT)
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** TVDスキーム(XI方向)                                           ****
!***********************************************************************
SUBROUTINE TVDXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  REAL    :: AKX, AKY, AKZ
  REAL    :: RWL, RWR, RWD, UM, VM, WM, HM, CM
  REAL    :: THI, PHI, BETA, TIL, AKXT, AKYT, AKZT, THIT, PH
  REAL    :: FSIGMA, FGAMMA, ABSA, AP, AM
  REAL    :: D3FP1, D3FP2, D3FP3, D3FM1, D3FM2, D3FM3
  REAL, ALLOCATABLE :: F(:, :), EIG(:, :)
  REAL, ALLOCATABLE :: RR(:, :, :), RI(:, :, :), EIGM(:, :), DELTA(:, :)
  REAL, ALLOCATABLE :: DU(:, :), DUR(:, :)
  REAL, ALLOCATABLE :: GG(:, :), PHIM(:, :)
  REAL, ALLOCATABLE :: DFP(:, :), DFM(:, :)
  REAL, ALLOCATABLE :: D3FP(:, :), D3FM(:, :)
  REAL, ALLOCATABLE :: DFBP(:, :), DFBM(:, :)
  REAL, ALLOCATABLE :: DFTP(:, :), DFTM(:, :)
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, L, &
  !$  AKX, AKY, AKZ, RWL, RWR, RWD, UM, VM, WM, HM, CM, &
  !$  THI, PHI, BETA, TIL, AKXT, AKYT, AKZT, THIT, PH, &
  !$  FSIGMA, FGAMMA, ABSA, AP, AM, &
  !$  D3FP1, D3FP2, D3FP3, D3FM1, D3FM2, D3FM3, &
  !$  F, EIG, RR, RI, EIGM, DELTA, DU, DUR, &
  !$  GG, PHIM, DFP, DFM, D3FP, D3FM, &
  !$  DFBP, DFBM, DFTP, DFTM &
  !$)
  ALLOCATE(F(IS: IE, LS: LE), EIG(IS: IE, LS: LE))
  ALLOCATE(RR(IS: IE - 1, LS: 5, LS: 5))
  ALLOCATE(RI(IS: IE - 1, LS: 5, LS: 5))
  ALLOCATE(EIGM(IS: IE - 1, LS: LE))
  ALLOCATE(DELTA(IS: IE - 1, LS: LE))
  ALLOCATE(DU(IS: IE - 1, LS: LE), DUR(IS: IE - 1, LS: LE))
  ALLOCATE(GG(IS: IE, LS: LE), PHIM(IS: IE - 1, LS: LE))
  ALLOCATE(DFP (IS    : IE - 1, LS: LE), DFM (IS    : IE - 1, LS: LE))
  ALLOCATE(D3FP(IS + 1: IE - 2, LS: LE), D3FM(IS + 1: IE - 2, LS: LE))
  ALLOCATE(DFBP(IS    : IE - 3, LS: LE), DFBM(IS + 1: IE - 2, LS: LE))
  ALLOCATE(DFTP(IS + 1: IE - 2, LS: LE), DFTM(IS + 2: IE - 1, LS: LE))
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
    ! 二次精度中心差分 +++++++++++++++++++++++++++++++++++++++++++++++++
    DO I = IS, IE
    IF(QH(I,J,K,1) .GT. 0.0) THEN
      THI = XIX(I,J,K) * U(I,J,K) &
      &   + XIY(I,J,K) * V(I,J,K) &
      &   + XIZ(I,J,K) * W(I,J,K)
      TIL = SQRT(XIX(I,J,K)**2 + XIY(I,J,K)**2 + XIZ(I,J,K)**2)
      CM  = SQRT(GAMMA * RG * T(I,J,K))
      PH  = P(I,J,K) / AJA(I,J,K)
      F(I,1) = QH(I,J,K,1) * THI
      F(I,2) = QH(I,J,K,2) * THI + PH * XIX(I,J,K)
      F(I,3) = QH(I,J,K,3) * THI + PH * XIY(I,J,K)
      F(I,4) = QH(I,J,K,4) * THI + PH * XIZ(I,J,K)
      F(I,5) = QH(I,J,K,5) * THI + PH * THI
      DO L = 6, LE
        F(I,L) = QH(I,J,K,L) * THI
      ENDDO
      EIG(I,1) = THI
      EIG(I,2) = THI
      EIG(I,3) = THI
      EIG(I,4) = THI + CM * TIL
      EIG(I,5) = THI - CM * TIL
      DO L = 6, LE
        EIG(I,L) = THI
      ENDDO
    ENDIF
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO I = IS, IE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I+1,J,K,1) .GT. 0.0) THEN
      ET(I,J,K,L) = 0.5 * (F(I,L) + F(I+1,L))
    ENDIF
    ENDDO
    ENDDO
    IF(Order .EQ. 0) CYCLE
    ! 対角化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO I = IS, IE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I+1,J,K,1) .GT. 0.0) THEN
      AKX  = 0.5 * (XIX(I,J,K) + XIX(I+1,J,K))
      AKY  = 0.5 * (XIY(I,J,K) + XIY(I+1,J,K))
      AKZ  = 0.5 * (XIZ(I,J,K) + XIZ(I+1,J,K))
      ! Roe 平均
      RWD  = 1.0 / ( SQRT(QH(I  ,J,K,1) * AJA(I  ,J,K)) &
      &            + SQRT(QH(I+1,J,K,1) * AJA(I+1,J,K)) )
      RWL  = SQRT(QH(I  ,J,K,1) * AJA(I  ,J,K)) * RWD
      RWR  = SQRT(QH(I+1,J,K,1) * AJA(I+1,J,K)) * RWD
      UM   = RWL * U(I,J,K) + RWR * U(I+1,J,K)
      VM   = RWL * V(I,J,K) + RWR * V(I+1,J,K)
      WM   = RWL * W(I,J,K) + RWR * W(I+1,J,K)
      HM   = RWL * ( &
      &        0.5 * (U(I  ,J,K)**2 + V(I  ,J,K)**2 + W(I  ,J,K)**2) &
      &      + GAMMA * RG / (GAMMA - 1.0) * T(I  ,J,K) &
      &    ) &
      &    + RWR * ( &
      &        0.5 * (U(I+1,J,K)**2 + V(I+1,J,K)**2 + W(I+1,J,K)**2) &
      &      + GAMMA * RG / (GAMMA - 1.0) * T(I+1,J,K) &
      &    )
      CM   = SQRT((GAMMA - 1.0) * (HM - 0.5 * (UM**2 + VM**2 + WM**2)))
      ! 対角化
      THI  = AKX * UM + AKY * VM + AKZ * WM
      PHI  = 0.5 * (GAMMA - 1.0) * (UM**2 + VM**2 + WM**2)
      BETA = 1.0 / (2.0 * CM**2)
      TIL  = SQRT(AKX**2 + AKY**2 + AKZ**2)
      AKXT = AKX / TIL
      AKYT = AKY / TIL
      AKZT = AKZ / TIL
      THIT = AKXT * UM + AKYT * VM + AKZT * WM
      RR(I,1,1) = AKXT
      RR(I,2,1) = AKXT * UM
      RR(I,3,1) = AKXT * VM + AKZT
      RR(I,4,1) = AKXT * WM - AKYT
      RR(I,5,1) = AKXT * PHI / (GAMMA - 1.0) + AKZT * VM - AKYT * WM
      RR(I,1,2) = AKYT
      RR(I,2,2) = AKYT * UM - AKZT
      RR(I,3,2) = AKYT * VM
      RR(I,4,2) = AKYT * WM + AKXT
      RR(I,5,2) = AKYT * PHI / (GAMMA - 1.0) + AKXT * WM - AKZT * UM
      RR(I,1,3) = AKZT
      RR(I,2,3) = AKZT * UM + AKYT
      RR(I,3,3) = AKZT * VM - AKXT
      RR(I,4,3) = AKZT * WM
      RR(I,5,3) = AKZT * PHI / (GAMMA - 1.0) + AKYT * UM - AKXT * VM
      RR(I,1,4) = 1.0
      RR(I,2,4) = UM + AKXT * CM
      RR(I,3,4) = VM + AKYT * CM
      RR(I,4,4) = WM + AKZT * CM
      RR(I,5,4) = (PHI + CM**2) / (GAMMA - 1.0) + CM * THIT
      RR(I,1,5) = 1.0
      RR(I,2,5) = UM - AKXT * CM
      RR(I,3,5) = VM - AKYT * CM
      RR(I,4,5) = WM - AKZT * CM
      RR(I,5,5) = (PHI + CM**2) / (GAMMA - 1.0) - CM * THIT
      RI(I,1,1) = AKXT - AKXT * PHI / CM**2 + AKYT * WM - AKZT * VM
      RI(I,2,1) = AKYT - AKYT * PHI / CM**2 + AKZT * UM - AKXT * WM
      RI(I,3,1) = AKZT - AKZT * PHI / CM**2 + AKXT * VM - AKYT * UM
      RI(I,4,1) = BETA * (PHI - CM * THIT)
      RI(I,5,1) = BETA * (PHI + CM * THIT)
      RI(I,1,2) = (GAMMA - 1.0) * AKXT * UM / CM**2
      RI(I,2,2) = (GAMMA - 1.0) * AKYT * UM / CM**2 - AKZT
      RI(I,3,2) = (GAMMA - 1.0) * AKZT * UM / CM**2 + AKYT
      RI(I,4,2) = BETA * (AKXT * CM - (GAMMA - 1.0) * UM)
      RI(I,5,2) =-BETA * (AKXT * CM + (GAMMA - 1.0) * UM)
      RI(I,1,3) = (GAMMA - 1.0) * AKXT * VM / CM**2 + AKZT
      RI(I,2,3) = (GAMMA - 1.0) * AKYT * VM / CM**2
      RI(I,3,3) = (GAMMA - 1.0) * AKZT * VM / CM**2 - AKXT
      RI(I,4,3) = BETA * (AKYT * CM - (GAMMA - 1.0) * VM)
      RI(I,5,3) =-BETA * (AKYT * CM + (GAMMA - 1.0) * VM)
      RI(I,1,4) = (GAMMA - 1.0) * AKXT * WM / CM**2 - AKYT
      RI(I,2,4) = (GAMMA - 1.0) * AKYT * WM / CM**2 + AKXT
      RI(I,3,4) = (GAMMA - 1.0) * AKZT * WM / CM**2
      RI(I,4,4) = BETA * (AKZT * CM - (GAMMA - 1.0) * WM)
      RI(I,5,4) =-BETA * (AKZT * CM + (GAMMA - 1.0) * WM)
      RI(I,1,5) =-(GAMMA - 1.0) * AKXT / CM**2
      RI(I,2,5) =-(GAMMA - 1.0) * AKYT / CM**2
      RI(I,3,5) =-(GAMMA - 1.0) * AKZT / CM**2
      RI(I,4,5) = BETA * (GAMMA - 1.0)
      RI(I,5,5) = BETA * (GAMMA - 1.0)
      EIGM(I,1) = THI
      EIGM(I,2) = THI
      EIGM(I,3) = THI
      EIGM(I,4) = THI + CM * TIL
      EIGM(I,5) = THI - CM * TIL
      DO L = 6, LE
        EIGM(I,L) = THI
      ENDDO
      DO L = LS, LE
        DELTA(I,L) = MAX( &
        &              0.0, &
        &              EIGM(I,L) - EIG(I,L), EIG(I+1,L) - EIGM(I,L) &
        &            ) + DEL * CM * TIL
      ENDDO
    ENDIF
    ENDDO
    ! Delta u
    DO L = LS, LE
    DO I = IS, IE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I+1,J,K,1) .GT. 0.0) THEN
      DU(I,L) = - QH(I  ,J,K,L) * AJA(I  ,J,K) &
      &         + QH(I+1,J,K,L) * AJA(I+1,J,K)
      DU(I,L) = DU(I,L) * 2.0 / (AJA(I,J,K) + AJA(I+1,J,K))
    ENDIF
    ENDDO
    ENDDO
    ! R^-1 Delta u
    DO L = LS, LE
    DO I = IS, IE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I+1,J,K,1) .GT. 0.0) THEN
      IF(L .LE. 5) THEN
        DUR(I,L) = RI(I,L,1) * DU(I,1) + RI(I,L,2) * DU(I,2) &
        &        + RI(I,L,3) * DU(I,3) + RI(I,L,4) * DU(I,4) &
        &        + RI(I,L,5) * DU(I,5)
      ELSE
        DUR(I,L) = DU(I,L)
      ENDIF
    ENDIF
    ENDDO
    ENDDO
    ! 二次精度TVDスキーム ++++++++++++++++++++++++++++++++++++++++++++++
    IF(Order .EQ. 2) THEN
      ! g
      DO L = LS, LE
      I = IS
      IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I+1,J,K,1) .GT. 0.0) THEN
        GG(I,L) = DUR(I,L)
      ENDIF
      DO I = IS + 1, IE - 1
      IF(QH(I,J,K,1) .LE. 0.0) CYCLE
      IF(QH(I-1,J,K,1) .GT. 0.0 .AND. QH(I+1,J,K,1) .GT. 0.0) THEN
        GG(I,L) = minmod2(DUR(I,L), DUR(I-1,L))
      ELSEIF(QH(I-1,J,K,1) .GT. 0.0) THEN
        GG(I,L) = DUR(I-1,L)
      ELSEIF(QH(I+1,J,K,1) .GT. 0.0) THEN
        GG(I,L) = DUR(I,L)
      ENDIF
      ENDDO
      I = IE
      IF(QH(I-1,J,K,1) .GT. 0.0 .AND. QH(I,J,K,1) .GT. 0.0) THEN
        GG(I,L) = DUR(I-1,L)
      ENDIF
      ENDDO
      ! phi
      DO L = LS, LE
      DO I = IS, IE - 1
      IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I+1,J,K,1) .GT. 0.0) THEN
        FSIGMA = Func_sigma(Func_psi(EIGM(I,L), DELTA(I,L)))
        FGAMMA = Func_gamma(FSIGMA, GG(I,L), GG(I+1,L), DUR(I,L))
        PHIM(I,L) = FSIGMA * (GG(I,L) + GG(I+1,L)) &
        &         - Func_psi(EIGM(I,L) + FGAMMA, DELTA(I,L)) * DUR(I,L)
      ENDIF
      ENDDO
      ENDDO
      ! 数値流束
      DO L = LS, LE
      DO I = IS, IE - 1
      IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I+1,J,K,1) .GT. 0.0) THEN
        IF(L .LE. 5) THEN
          ET(I,J,K,L) = ET(I,J,K,L) + 0.5 * ( RR(I,L,1) * PHIM(I,1) &
          &                                 + RR(I,L,2) * PHIM(I,2) &
          &                                 + RR(I,L,3) * PHIM(I,3) &
          &                                 + RR(I,L,4) * PHIM(I,4) &
          &                                 + RR(I,L,5) * PHIM(I,5) )
        ELSE
          ET(I,J,K,L) = ET(I,J,K,L) + 0.5 * PHIM(I,L)
        ENDIF
      ENDIF
      ENDDO
      ENDDO
      CYCLE
    ENDIF
    ! 一次精度風上差分 +++++++++++++++++++++++++++++++++++++++++++++++++
    ! Delta f^+-
    DO L = LS, LE
    DO I = IS, IE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I+1,J,K,1) .GT. 0.0) THEN
      ABSA  = Func_psi(EIGM(I,L), DELTA(I,L))
      AP    = 0.5 * (EIGM(I,L) + ABSA)
      AM    = 0.5 * (EIGM(I,L) - ABSA)
      DFP(I,L) = AP * DUR(I,L)
      DFM(I,L) = AM * DUR(I,L)
    ENDIF
    ENDDO
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO I = IS, IE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I+1,J,K,1) .GT. 0.0) THEN
      IF(L .LE. 5) THEN
        ET(I,J,K,L) = ET(I,J,K,L) - 0.5 * ( RR(I,L,1) * DFP(I,1) &
        &                                 + RR(I,L,2) * DFP(I,2) &
        &                                 + RR(I,L,3) * DFP(I,3) &
        &                                 + RR(I,L,4) * DFP(I,4) &
        &                                 + RR(I,L,5) * DFP(I,5) &
        &                                 - RR(I,L,1) * DFM(I,1) &
        &                                 - RR(I,L,2) * DFM(I,2) &
        &                                 - RR(I,L,3) * DFM(I,3) &
        &                                 - RR(I,L,4) * DFM(I,4) &
        &                                 - RR(I,L,5) * DFM(I,5) )
      ELSE
        ET(I,J,K,L) = ET(I,J,K,L) - 0.5 * (DFP(I,L) - DFM(I,L))
      ENDIF
    ENDIF
    ENDDO
    ENDDO
    IF(Order .EQ. 1) CYCLE
    ! 三次&四次精度TVDスキーム +++++++++++++++++++++++++++++++++++++++++
    IF(Order .EQ. 4) THEN
      ! Delta^3 bar(f)^+-
      DO L = LS, LE
      DO I = IS + 1, IE - 2
      IF( QH(I-1,J,K,1) .GT. 0.0 .AND. QH(I  ,J,K,1) .GT. 0.0 .AND. &
      &   QH(I+1,J,K,1) .GT. 0.0 .AND. QH(I+2,J,K,1) .GT. 0.0 ) THEN
        D3FP1 = minmod3(DFP(I-1,L), b2 * DFP(I  ,L), b2 * DFP(I+1,L))
        D3FP2 = minmod3(DFP(I  ,L), b2 * DFP(I+1,L), b2 * DFP(I-1,L))
        D3FP3 = minmod3(DFP(I+1,L), b2 * DFP(I-1,L), b2 * DFP(I  ,L))
        D3FM1 = minmod3(DFM(I-1,L), b2 * DFM(I  ,L), b2 * DFM(I+1,L))
        D3FM2 = minmod3(DFM(I  ,L), b2 * DFM(I+1,L), b2 * DFM(I-1,L))
        D3FM3 = minmod3(DFM(I+1,L), b2 * DFM(I-1,L), b2 * DFM(I  ,L))
        D3FP(I,L) = D3FP1 - 2.0 * D3FP2 + D3FP3
        D3FM(I,L) = D3FM1 - 2.0 * D3FM2 + D3FM3
      ENDIF
      ENDDO
      ENDDO
      ! Delta^asterisk f^+-
      DO L = LS, LE
      DO I = IS + 1, IE - 2
      IF( QH(I-1,J,K,1) .GT. 0.0 .AND. QH(I  ,J,K,1) .GT. 0.0 .AND. &
      &   QH(I+1,J,K,1) .GT. 0.0 .AND. QH(I+2,J,K,1) .GT. 0.0 ) THEN
        DFP(I,L) = DFP(I,L) - D3FP(I,L) / 6.0
        DFM(I,L) = DFM(I,L) - D3FM(I,L) / 6.0
      ENDIF
      ENDDO
      ENDDO
    ENDIF
    ! Delta bar(f)^+-, Delta tilde(f)^+-
    DO L = LS, LE
    DO I = IS + 1, IE - 2
    IF( QH(I-1,J,K,1) .GT. 0.0 .AND. QH(I  ,J,K,1) .GT. 0.0 .AND. &
    &   QH(I+1,J,K,1) .GT. 0.0 .AND. QH(I+2,J,K,1) .GT. 0.0 ) THEN
      DFBP(I-1,L) = minmod2(DFP(I-1,L), b1 * DFP(I  ,L))
      DFBM(I  ,L) = minmod2(DFM(I  ,L), b1 * DFM(I+1,L))
      DFTP(I  ,L) = minmod2(DFP(I  ,L), b1 * DFP(I-1,L))
      DFTM(I+1,L) = minmod2(DFM(I+1,L), b1 * DFM(I  ,L))
    ENDIF
    ENDDO
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO I = IS + 1, IE - 2
    IF( QH(I-1,J,K,1) .GT. 0.0 .AND. QH(I  ,J,K,1) .GT. 0.0 .AND. &
    &   QH(I+1,J,K,1) .GT. 0.0 .AND. QH(I+2,J,K,1) .GT. 0.0 ) THEN
      IF(L .LE. 5) THEN
        ET(I,J,K,L) = ET(I,J,K,L) &
        &           + ( 1.0 * RR(I-1,L,1) * DFBP(I-1,1) &
        &             + 1.0 * RR(I-1,L,2) * DFBP(I-1,2) &
        &             + 1.0 * RR(I-1,L,3) * DFBP(I-1,3) &
        &             + 1.0 * RR(I-1,L,4) * DFBP(I-1,4) &
        &             + 1.0 * RR(I-1,L,5) * DFBP(I-1,5) &
        &             + 2.0 * RR(I  ,L,1) * DFTP(I  ,1) &
        &             + 2.0 * RR(I  ,L,2) * DFTP(I  ,2) &
        &             + 2.0 * RR(I  ,L,3) * DFTP(I  ,3) &
        &             + 2.0 * RR(I  ,L,4) * DFTP(I  ,4) &
        &             + 2.0 * RR(I  ,L,5) * DFTP(I  ,5) &
        &             - 2.0 * RR(I  ,L,1) * DFBM(I  ,1) &
        &             - 2.0 * RR(I  ,L,2) * DFBM(I  ,2) &
        &             - 2.0 * RR(I  ,L,3) * DFBM(I  ,3) &
        &             - 2.0 * RR(I  ,L,4) * DFBM(I  ,4) &
        &             - 2.0 * RR(I  ,L,5) * DFBM(I  ,5) &
        &             - 1.0 * RR(I+1,L,1) * DFTM(I+1,1) &
        &             - 1.0 * RR(I+1,L,2) * DFTM(I+1,2) &
        &             - 1.0 * RR(I+1,L,3) * DFTM(I+1,3) &
        &             - 1.0 * RR(I+1,L,4) * DFTM(I+1,4) &
        &             - 1.0 * RR(I+1,L,5) * DFTM(I+1,5) &
        &           ) / 6.0
      ELSE
        ET(I,J,K,L) = ET(I,J,K,L) &
        &           + ( 1.0 * DFBP(I-1,L) + 2.0 * DFTP(I  ,L) &
        &             - 2.0 * DFBM(I  ,L) - 1.0 * DFTM(I+1,L) ) / 6.0
      ENDIF
    ENDIF
    ENDDO
    ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ! メモリ解放 =========================================================
  DEALLOCATE(F,EIG,RR,RI,EIGM,DELTA,DU,DUR,GG,PHIM,DFP,DFM,D3FP,D3FM,DFBP,DFBM,DFTP,DFTM)
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE TVDXI
!***********************************************************************
!**** TVDスキーム(ET方向)                                           ****
!***********************************************************************
SUBROUTINE TVDET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  REAL    :: AKX, AKY, AKZ
  REAL    :: RWL, RWR, RWD, UM, VM, WM, HM, CM
  REAL    :: THI, PHI, BETA, TIL, AKXT, AKYT, AKZT, THIT, PH
  REAL    :: FSIGMA, FGAMMA, ABSA, AP, AM
  REAL    :: D3FP1, D3FP2, D3FP3, D3FM1, D3FM2, D3FM3
  REAL, ALLOCATABLE :: F(:, :), EIG(:, :)
  REAL, ALLOCATABLE :: RR(:, :, :), RI(:, :, :), EIGM(:, :), DELTA(:, :)
  REAL, ALLOCATABLE :: DU(:, :), DUR(:, :)
  REAL, ALLOCATABLE :: GG(:, :), PHIM(:, :)
  REAL, ALLOCATABLE :: DFP(:, :), DFM(:, :)
  REAL, ALLOCATABLE :: D3FP(:, :), D3FM(:, :)
  REAL, ALLOCATABLE :: DFBP(:, :), DFBM(:, :)
  REAL, ALLOCATABLE :: DFTP(:, :), DFTM(:, :)
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, L, &
  !$  AKX, AKY, AKZ, RWL, RWR, RWD, UM, VM, WM, HM, CM, &
  !$  THI, PHI, BETA, TIL, AKXT, AKYT, AKZT, THIT, PH, &
  !$  FSIGMA, FGAMMA, ABSA, AP, AM, &
  !$  D3FP1, D3FP2, D3FP3, D3FM1, D3FM2, D3FM3, &
  !$  F, EIG, RR, RI, EIGM, DELTA, DU, DUR, &
  !$  GG, PHIM, DFP, DFM, D3FP, D3FM, &
  !$  DFBP, DFBM, DFTP, DFTM &
  !$)
  ALLOCATE(F(JS: JE, LS: LE), EIG(JS: JE, LS: LE))
  ALLOCATE(RR(JS: JE - 1, LS: 5, LS: 5))
  ALLOCATE(RI(JS: JE - 1, LS: 5, LS: 5))
  ALLOCATE(EIGM(JS: JE - 1, LS: LE))
  ALLOCATE(DELTA(JS: JE - 1, LS: LE))
  ALLOCATE(DU(JS: JE - 1, LS: LE), DUR(JS: JE - 1, LS: LE))
  ALLOCATE(GG(JS: JE, LS: LE), PHIM(JS: JE - 1, LS: LE))
  ALLOCATE(DFP (JS    : JE - 1, LS: LE), DFM (JS    : JE - 1, LS: LE))
  ALLOCATE(D3FP(JS + 1: JE - 2, LS: LE), D3FM(JS + 1: JE - 2, LS: LE))
  ALLOCATE(DFBP(JS    : JE - 3, LS: LE), DFBM(JS + 1: JE - 2, LS: LE))
  ALLOCATE(DFTP(JS + 1: JE - 2, LS: LE), DFTM(JS + 2: JE - 1, LS: LE))
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO I = IS + 1, IE - 1
    ! 二次精度中心差分 +++++++++++++++++++++++++++++++++++++++++++++++++
    DO J = JS, JE
    IF(QH(I,J,K,1) .GT. 0.0) THEN
      THI = ETX(I,J,K) * U(I,J,K) &
      &   + ETY(I,J,K) * V(I,J,K) &
      &   + ETZ(I,J,K) * W(I,J,K)
      TIL = SQRT(ETX(I,J,K)**2 + ETY(I,J,K)**2 + ETZ(I,J,K)**2)
      CM  = SQRT(GAMMA * RG * T(I,J,K))
      PH  = P(I,J,K) / AJA(I,J,K)
      F(J,1) = QH(I,J,K,1) * THI
      F(J,2) = QH(I,J,K,2) * THI + PH * ETX(I,J,K)
      F(J,3) = QH(I,J,K,3) * THI + PH * ETY(I,J,K)
      F(J,4) = QH(I,J,K,4) * THI + PH * ETZ(I,J,K)
      F(J,5) = QH(I,J,K,5) * THI + PH * THI
      DO L = 6, LE
        F(J,L) = QH(I,J,K,L) * THI
      ENDDO
      EIG(J,1) = THI
      EIG(J,2) = THI
      EIG(J,3) = THI
      EIG(J,4) = THI + CM * TIL
      EIG(J,5) = THI - CM * TIL
      DO L = 6, LE
        EIG(J,L) = THI
      ENDDO
    ENDIF
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO J = JS, JE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J+1,K,1) .GT. 0.0) THEN
      FT(I,J,K,L) = 0.5 * (F(J,L) + F(J+1,L))
    ENDIF
    ENDDO
    ENDDO
    IF(Order .EQ. 0) CYCLE
    ! 対角化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO J = JS, JE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J+1,K,1) .GT. 0.0) THEN
      AKX  = 0.5 * (ETX(I,J,K) + ETX(I,J+1,K))
      AKY  = 0.5 * (ETY(I,J,K) + ETY(I,J+1,K))
      AKZ  = 0.5 * (ETZ(I,J,K) + ETZ(I,J+1,K))
      ! Roe 平均
      RWD  = 1.0 / ( SQRT(QH(I,J  ,K,1) * AJA(I,J  ,K)) &
      &            + SQRT(QH(I,J+1,K,1) * AJA(I,J+1,K)) )
      RWL  = SQRT(QH(I,J  ,K,1) * AJA(I,J  ,K)) * RWD
      RWR  = SQRT(QH(I,J+1,K,1) * AJA(I,J+1,K)) * RWD
      UM   = RWL * U(I,J,K) + RWR * U(I,J+1,K)
      VM   = RWL * V(I,J,K) + RWR * V(I,J+1,K)
      WM   = RWL * W(I,J,K) + RWR * W(I,J+1,K)
      HM   = RWL * ( &
      &        0.5 * (U(I,J  ,K)**2 + V(I,J  ,K)**2 + W(I,J  ,K)**2) &
      &      + GAMMA * RG / (GAMMA - 1.0) * T(I,J  ,K) &
      &    ) &
      &    + RWR * ( &
      &        0.5 * (U(I,J+1,K)**2 + V(I,J+1,K)**2 + W(I,J+1,K)**2) &
      &      + GAMMA * RG / (GAMMA - 1.0) * T(I,J+1,K) &
      &    )
      CM   = SQRT((GAMMA - 1.0) * (HM - 0.5 * (UM**2 + VM**2 + WM**2)))
      ! 対角化
      THI  = AKX * UM + AKY * VM + AKZ * WM
      PHI  = 0.5 * (GAMMA - 1.0) * (UM**2 + VM**2 + WM**2)
      BETA = 1.0 / (2.0 * CM**2)
      TIL  = SQRT(AKX**2 + AKY**2 + AKZ**2)
      AKXT = AKX / TIL
      AKYT = AKY / TIL
      AKZT = AKZ / TIL
      THIT = AKXT * UM + AKYT * VM + AKZT * WM
      RR(J,1,1) = AKXT
      RR(J,2,1) = AKXT * UM
      RR(J,3,1) = AKXT * VM + AKZT
      RR(J,4,1) = AKXT * WM - AKYT
      RR(J,5,1) = AKXT * PHI / (GAMMA - 1.0) + AKZT * VM - AKYT * WM
      RR(J,1,2) = AKYT
      RR(J,2,2) = AKYT * UM - AKZT
      RR(J,3,2) = AKYT * VM
      RR(J,4,2) = AKYT * WM + AKXT
      RR(J,5,2) = AKYT * PHI / (GAMMA - 1.0) + AKXT * WM - AKZT * UM
      RR(J,1,3) = AKZT
      RR(J,2,3) = AKZT * UM + AKYT
      RR(J,3,3) = AKZT * VM - AKXT
      RR(J,4,3) = AKZT * WM
      RR(J,5,3) = AKZT * PHI / (GAMMA - 1.0) + AKYT * UM - AKXT * VM
      RR(J,1,4) = 1.0
      RR(J,2,4) = UM + AKXT * CM
      RR(J,3,4) = VM + AKYT * CM
      RR(J,4,4) = WM + AKZT * CM
      RR(J,5,4) = (PHI + CM**2) / (GAMMA - 1.0) + CM * THIT
      RR(J,1,5) = 1.0
      RR(J,2,5) = UM - AKXT * CM
      RR(J,3,5) = VM - AKYT * CM
      RR(J,4,5) = WM - AKZT * CM
      RR(J,5,5) = (PHI + CM**2) / (GAMMA - 1.0) - CM * THIT
      RI(J,1,1) = AKXT - AKXT * PHI / CM**2 + AKYT * WM - AKZT * VM
      RI(J,2,1) = AKYT - AKYT * PHI / CM**2 + AKZT * UM - AKXT * WM
      RI(J,3,1) = AKZT - AKZT * PHI / CM**2 + AKXT * VM - AKYT * UM
      RI(J,4,1) = BETA * (PHI - CM * THIT)
      RI(J,5,1) = BETA * (PHI + CM * THIT)
      RI(J,1,2) = (GAMMA - 1.0) * AKXT * UM / CM**2
      RI(J,2,2) = (GAMMA - 1.0) * AKYT * UM / CM**2 - AKZT
      RI(J,3,2) = (GAMMA - 1.0) * AKZT * UM / CM**2 + AKYT
      RI(J,4,2) = BETA * (AKXT * CM - (GAMMA - 1.0) * UM)
      RI(J,5,2) =-BETA * (AKXT * CM + (GAMMA - 1.0) * UM)
      RI(J,1,3) = (GAMMA - 1.0) * AKXT * VM / CM**2 + AKZT
      RI(J,2,3) = (GAMMA - 1.0) * AKYT * VM / CM**2
      RI(J,3,3) = (GAMMA - 1.0) * AKZT * VM / CM**2 - AKXT
      RI(J,4,3) = BETA * (AKYT * CM - (GAMMA - 1.0) * VM)
      RI(J,5,3) =-BETA * (AKYT * CM + (GAMMA - 1.0) * VM)
      RI(J,1,4) = (GAMMA - 1.0) * AKXT * WM / CM**2 - AKYT
      RI(J,2,4) = (GAMMA - 1.0) * AKYT * WM / CM**2 + AKXT
      RI(J,3,4) = (GAMMA - 1.0) * AKZT * WM / CM**2
      RI(J,4,4) = BETA * (AKZT * CM - (GAMMA - 1.0) * WM)
      RI(J,5,4) =-BETA * (AKZT * CM + (GAMMA - 1.0) * WM)
      RI(J,1,5) =-(GAMMA - 1.0) * AKXT / CM**2
      RI(J,2,5) =-(GAMMA - 1.0) * AKYT / CM**2
      RI(J,3,5) =-(GAMMA - 1.0) * AKZT / CM**2
      RI(J,4,5) = BETA * (GAMMA - 1.0)
      RI(J,5,5) = BETA * (GAMMA - 1.0)
      EIGM(J,1) = THI
      EIGM(J,2) = THI
      EIGM(J,3) = THI
      EIGM(J,4) = THI + CM * TIL
      EIGM(J,5) = THI - CM * TIL
      DO L = 6, LE
        EIGM(J,L) = THI
      ENDDO
      DO L = LS, LE
        DELTA(J,L) = MAX( &
        &              0.0, &
        &              EIGM(J,L) - EIG(J,L), EIG(J+1,L) - EIGM(J,L) &
        &            ) + DEL * CM * TIL
      ENDDO
    ENDIF
    ENDDO
    ! Delta u
    DO L = LS, LE
    DO J = JS, JE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J+1,K,1) .GT. 0.0) THEN
      DU(J,L) = - QH(I,J  ,K,L) * AJA(I,J  ,K) &
      &         + QH(I,J+1,K,L) * AJA(I,J+1,K)
      DU(J,L) = DU(J,L) * 2.0 / (AJA(I,J,K) + AJA(I,J+1,K))
    ENDIF
    ENDDO
    ENDDO
    ! R^-1 Delta u
    DO L = LS, LE
    DO J = JS, JE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J+1,K,1) .GT. 0.0) THEN
      IF(L .LE. 5) THEN
        DUR(J,L) = RI(J,L,1) * DU(J,1) + RI(J,L,2) * DU(J,2) &
        &        + RI(J,L,3) * DU(J,3) + RI(J,L,4) * DU(J,4) &
        &        + RI(J,L,5) * DU(J,5)
      ELSE
        DUR(J,L) = DU(J,L)
      ENDIF
    ENDIF
    ENDDO
    ENDDO
    ! 二次精度TVDスキーム ++++++++++++++++++++++++++++++++++++++++++++++
    IF(Order .EQ. 2) THEN
      ! g
      DO L = LS, LE
      J = JS
      IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J+1,K,1) .GT. 0.0) THEN
        GG(J,L) = DUR(J,L)
      ENDIF
      DO J = JS + 1, JE - 1
      IF(QH(I,J,K,1) .LE. 0.0) CYCLE
      IF(QH(I,J-1,K,1) .GT. 0.0 .AND. QH(I,J+1,K,1) .GT. 0.0) THEN
        GG(J,L) = minmod2(DUR(J,L), DUR(J-1,L))
      ELSEIF(QH(I,J-1,K,1) .GT. 0.0) THEN
        GG(J,L) = DUR(J-1,L)
      ELSEIF(QH(I,J+1,K,1) .GT. 0.0) THEN
        GG(J,L) = DUR(J,L)
      ENDIF
      ENDDO
      J = JE
      IF(QH(I,J-1,K,1) .GT. 0.0 .AND. QH(I,J,K,1) .GT. 0.0) THEN
        GG(J,L) = DUR(J-1,L)
      ENDIF
      ENDDO
      ! phi
      DO L = LS, LE
      DO J = JS, JE - 1
      IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J+1,K,1) .GT. 0.0) THEN
        FSIGMA = Func_sigma(Func_psi(EIGM(J,L), DELTA(J,L)))
        FGAMMA = Func_gamma(FSIGMA, GG(J,L), GG(J+1,L), DUR(J,L))
        PHIM(J,L) = FSIGMA * (GG(J,L) + GG(J+1,L)) &
        &         - Func_psi(EIGM(J,L) + FGAMMA, DELTA(J,L)) * DUR(J,L)
      ENDIF
      ENDDO
      ENDDO
      ! 数値流束
      DO L = LS, LE
      DO J = JS, JE - 1
      IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J+1,K,1) .GT. 0.0) THEN
        IF(L .LE. 5) THEN
          FT(I,J,K,L) = FT(I,J,K,L) + 0.5 * ( RR(J,L,1) * PHIM(J,1) &
          &                                 + RR(J,L,2) * PHIM(J,2) &
          &                                 + RR(J,L,3) * PHIM(J,3) &
          &                                 + RR(J,L,4) * PHIM(J,4) &
          &                                 + RR(J,L,5) * PHIM(J,5) )
        ELSE
          FT(I,J,K,L) = FT(I,J,K,L) + 0.5 * PHIM(J,L)
        ENDIF
      ENDIF
      ENDDO
      ENDDO
      CYCLE
    ENDIF
    ! 一次精度風上差分 +++++++++++++++++++++++++++++++++++++++++++++++++
    ! Delta f^+-
    DO L = LS, LE
    DO J = JS, JE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J+1,K,1) .GT. 0.0) THEN
      ABSA  = Func_psi(EIGM(J,L), DELTA(J,L))
      AP    = 0.5 * (EIGM(J,L) + ABSA)
      AM    = 0.5 * (EIGM(J,L) - ABSA)
      DFP(J,L) = AP * DUR(J,L)
      DFM(J,L) = AM * DUR(J,L)
    ENDIF
    ENDDO
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO J = JS, JE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J+1,K,1) .GT. 0.0) THEN
      IF(L .LE. 5) THEN
        FT(I,J,K,L) = FT(I,J,K,L) - 0.5 * ( RR(J,L,1) * DFP(J,1) &
        &                                 + RR(J,L,2) * DFP(J,2) &
        &                                 + RR(J,L,3) * DFP(J,3) &
        &                                 + RR(J,L,4) * DFP(J,4) &
        &                                 + RR(J,L,5) * DFP(J,5) &
        &                                 - RR(J,L,1) * DFM(J,1) &
        &                                 - RR(J,L,2) * DFM(J,2) &
        &                                 - RR(J,L,3) * DFM(J,3) &
        &                                 - RR(J,L,4) * DFM(J,4) &
        &                                 - RR(J,L,5) * DFM(J,5) )
      ELSE
        FT(I,J,K,L) = FT(I,J,K,L) - 0.5 * (DFP(J,L) - DFM(J,L))
      ENDIF
    ENDIF
    ENDDO
    ENDDO
    IF(Order .EQ. 1) CYCLE
    ! 三次&四次精度TVDスキーム +++++++++++++++++++++++++++++++++++++++++
    IF(Order .EQ. 4) THEN
      ! Delta^3 bar(f)^+-
      DO L = LS, LE
      DO J = JS + 1, JE - 2
      IF( QH(I,J-1,K,1) .GT. 0.0 .AND. QH(I,J  ,K,1) .GT. 0.0 .AND. &
      &   QH(I,J+1,K,1) .GT. 0.0 .AND. QH(I,J+2,K,1) .GT. 0.0 ) THEN
        D3FP1 = minmod3(DFP(J-1,L), b2 * DFP(J  ,L), b2 * DFP(J+1,L))
        D3FP2 = minmod3(DFP(J  ,L), b2 * DFP(J+1,L), b2 * DFP(J-1,L))
        D3FP3 = minmod3(DFP(J+1,L), b2 * DFP(J-1,L), b2 * DFP(J  ,L))
        D3FM1 = minmod3(DFM(J-1,L), b2 * DFM(J  ,L), b2 * DFM(J+1,L))
        D3FM2 = minmod3(DFM(J  ,L), b2 * DFM(J+1,L), b2 * DFM(J-1,L))
        D3FM3 = minmod3(DFM(J+1,L), b2 * DFM(J-1,L), b2 * DFM(J  ,L))
        D3FP(J,L) = D3FP1 - 2.0 * D3FP2 + D3FP3
        D3FM(J,L) = D3FM1 - 2.0 * D3FM2 + D3FM3
      ENDIF
      ENDDO
      ENDDO
      ! Delta^asterisk f^+-
      DO L = LS, LE
      DO J = JS + 1, JE - 2
      IF( QH(I,J-1,K,1) .GT. 0.0 .AND. QH(I,J  ,K,1) .GT. 0.0 .AND. &
      &   QH(I,J+1,K,1) .GT. 0.0 .AND. QH(I,J+2,K,1) .GT. 0.0 ) THEN
        DFP(J,L) = DFP(J,L) - D3FP(J,L) / 6.0
        DFM(J,L) = DFM(J,L) - D3FM(J,L) / 6.0
      ENDIF
      ENDDO
      ENDDO
    ENDIF
    ! Delta bar(f)^+-, Delta tilde(f)^+-
    DO L = LS, LE
    DO J = JS + 1, JE - 2
    IF( QH(I,J-1,K,1) .GT. 0.0 .AND. QH(I,J  ,K,1) .GT. 0.0 .AND. &
    &   QH(I,J+1,K,1) .GT. 0.0 .AND. QH(I,J+2,K,1) .GT. 0.0 ) THEN
      DFBP(J-1,L) = minmod2(DFP(J-1,L), b1 * DFP(J  ,L))
      DFBM(J  ,L) = minmod2(DFM(J  ,L), b1 * DFM(J+1,L))
      DFTP(J  ,L) = minmod2(DFP(J  ,L), b1 * DFP(J-1,L))
      DFTM(J+1,L) = minmod2(DFM(J+1,L), b1 * DFM(J  ,L))
    ENDIF
    ENDDO
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO J = JS + 1, JE - 2
    IF( QH(I,J-1,K,1) .GT. 0.0 .AND. QH(I,J  ,K,1) .GT. 0.0 .AND. &
    &   QH(I,J+1,K,1) .GT. 0.0 .AND. QH(I,J+2,K,1) .GT. 0.0 ) THEN
      IF(L .LE. 5) THEN
        FT(I,J,K,L) = FT(I,J,K,L) &
        &           + ( 1.0 * RR(J-1,L,1) * DFBP(J-1,1) &
        &             + 1.0 * RR(J-1,L,2) * DFBP(J-1,2) &
        &             + 1.0 * RR(J-1,L,3) * DFBP(J-1,3) &
        &             + 1.0 * RR(J-1,L,4) * DFBP(J-1,4) &
        &             + 1.0 * RR(J-1,L,5) * DFBP(J-1,5) &
        &             + 2.0 * RR(J  ,L,1) * DFTP(J  ,1) &
        &             + 2.0 * RR(J  ,L,2) * DFTP(J  ,2) &
        &             + 2.0 * RR(J  ,L,3) * DFTP(J  ,3) &
        &             + 2.0 * RR(J  ,L,4) * DFTP(J  ,4) &
        &             + 2.0 * RR(J  ,L,5) * DFTP(J  ,5) &
        &             - 2.0 * RR(J  ,L,1) * DFBM(J  ,1) &
        &             - 2.0 * RR(J  ,L,2) * DFBM(J  ,2) &
        &             - 2.0 * RR(J  ,L,3) * DFBM(J  ,3) &
        &             - 2.0 * RR(J  ,L,4) * DFBM(J  ,4) &
        &             - 2.0 * RR(J  ,L,5) * DFBM(J  ,5) &
        &             - 1.0 * RR(J+1,L,1) * DFTM(J+1,1) &
        &             - 1.0 * RR(J+1,L,2) * DFTM(J+1,2) &
        &             - 1.0 * RR(J+1,L,3) * DFTM(J+1,3) &
        &             - 1.0 * RR(J+1,L,4) * DFTM(J+1,4) &
        &             - 1.0 * RR(J+1,L,5) * DFTM(J+1,5) &
        &           ) / 6.0
      ELSE
        FT(I,J,K,L) = FT(I,J,K,L) &
        &           + ( 1.0 * DFBP(J-1,L) + 2.0 * DFTP(J  ,L) &
        &             - 2.0 * DFBM(J  ,L) - 1.0 * DFTM(J+1,L) ) / 6.0
      ENDIF
    ENDIF
    ENDDO
    ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ! メモリ解放 =========================================================
  DEALLOCATE(F,EIG,RR,RI,EIGM,DELTA,DU,DUR,GG,PHIM,DFP,DFM,D3FP,D3FM,DFBP,DFBM,DFTP,DFTM)
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE TVDET
!***********************************************************************
!**** TVDスキーム(ZE方向)                                           ****
!***********************************************************************
SUBROUTINE TVDZE
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  REAL    :: AKX, AKY, AKZ
  REAL    :: RWL, RWR, RWD, UM, VM, WM, HM, CM
  REAL    :: THI, PHI, BETA, TIL, AKXT, AKYT, AKZT, THIT, PH
  REAL    :: FSIGMA, FGAMMA, ABSA, AP, AM
  REAL    :: D3FP1, D3FP2, D3FP3, D3FM1, D3FM2, D3FM3
  REAL, ALLOCATABLE :: F(:, :), EIG(:, :)
  REAL, ALLOCATABLE :: RR(:, :, :), RI(:, :, :), EIGM(:, :), DELTA(:, :)
  REAL, ALLOCATABLE :: DU(:, :), DUR(:, :)
  REAL, ALLOCATABLE :: GG(:, :), PHIM(:, :)
  REAL, ALLOCATABLE :: DFP(:, :), DFM(:, :)
  REAL, ALLOCATABLE :: D3FP(:, :), D3FM(:, :)
  REAL, ALLOCATABLE :: DFBP(:, :), DFBM(:, :)
  REAL, ALLOCATABLE :: DFTP(:, :), DFTM(:, :)
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, L, &
  !$  AKX, AKY, AKZ, RWL, RWR, RWD, UM, VM, WM, HM, CM, &
  !$  THI, PHI, BETA, TIL, AKXT, AKYT, AKZT, THIT, PH, &
  !$  FSIGMA, FGAMMA, ABSA, AP, AM, &
  !$  D3FP1, D3FP2, D3FP3, D3FM1, D3FM2, D3FM3, &
  !$  F, EIG, RR, RI, EIGM, DELTA, DU, DUR, &
  !$  GG, PHIM, DFP, DFM, D3FP, D3FM, &
  !$  DFBP, DFBM, DFTP, DFTM &
  !$)
  ALLOCATE(F(KS: KE, LS: LE), EIG(KS: KE, LS: LE))
  ALLOCATE(RR(KS: KE - 1, LS: 5, LS: 5))
  ALLOCATE(RI(KS: KE - 1, LS: 5, LS: 5))
  ALLOCATE(EIGM(KS: KE - 1, LS: LE))
  ALLOCATE(DELTA(KS: KE - 1, LS: LE))
  ALLOCATE(DU(KS: KE - 1, LS: LE), DUR(KS: KE - 1, LS: LE))
  ALLOCATE(GG(KS: KE, LS: LE), PHIM(KS: KE - 1, LS: LE))
  ALLOCATE(DFP (KS    : KE - 1, LS: LE), DFM (KS    : KE - 1, LS: LE))
  ALLOCATE(D3FP(KS + 1: KE - 2, LS: LE), D3FM(KS + 1: KE - 2, LS: LE))
  ALLOCATE(DFBP(KS    : KE - 3, LS: LE), DFBM(KS + 1: KE - 2, LS: LE))
  ALLOCATE(DFTP(KS + 1: KE - 2, LS: LE), DFTM(KS + 2: KE - 1, LS: LE))
  !$OMP DO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    ! 二次精度中心差分 +++++++++++++++++++++++++++++++++++++++++++++++++
    DO K = KS, KE
    IF(QH(I,J,K,1) .GT. 0.0) THEN
      THI = ZEX(I,J,K) * U(I,J,K) &
      &   + ZEY(I,J,K) * V(I,J,K) &
      &   + ZEZ(I,J,K) * W(I,J,K)
      TIL = SQRT(ZEX(I,J,K)**2 + ZEY(I,J,K)**2 + ZEZ(I,J,K)**2)
      CM  = SQRT(GAMMA * RG * T(I,J,K))
      PH  = P(I,J,K) / AJA(I,J,K)
      F(K,1) = QH(I,J,K,1) * THI
      F(K,2) = QH(I,J,K,2) * THI + PH * ZEX(I,J,K)
      F(K,3) = QH(I,J,K,3) * THI + PH * ZEY(I,J,K)
      F(K,4) = QH(I,J,K,4) * THI + PH * ZEZ(I,J,K)
      F(K,5) = QH(I,J,K,5) * THI + PH * THI
      DO L = 6, LE
        F(K,L) = QH(I,J,K,L) * THI
      ENDDO
      EIG(K,1) = THI
      EIG(K,2) = THI
      EIG(K,3) = THI
      EIG(K,4) = THI + CM * TIL
      EIG(K,5) = THI - CM * TIL
      DO L = 6, LE
        EIG(K,L) = THI
      ENDDO
    ENDIF
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO K = KS, KE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K+1,1) .GT. 0.0) THEN
      GT(I,J,K,L) = 0.5 * (F(K,L) + F(K+1,L))
    ENDIF
    ENDDO
    ENDDO
    IF(Order .EQ. 0) CYCLE
    ! 対角化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO K = KS, KE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K+1,1) .GT. 0.0) THEN
      AKX  = 0.5 * (ZEX(I,J,K) + ZEX(I,J,K+1))
      AKY  = 0.5 * (ZEY(I,J,K) + ZEY(I,J,K+1))
      AKZ  = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J,K+1))
      ! Roe 平均
      RWD  = 1.0 / ( SQRT(QH(I,J,K  ,1) * AJA(I,J,K  )) &
      &            + SQRT(QH(I,J,K+1,1) * AJA(I,J,K+1)) )
      RWL  = SQRT(QH(I,J,K  ,1) * AJA(I,J,K  )) * RWD
      RWR  = SQRT(QH(I,J,K+1,1) * AJA(I,J,K+1)) * RWD
      UM   = RWL * U(I,J,K) + RWR * U(I,J,K+1)
      VM   = RWL * V(I,J,K) + RWR * V(I,J,K+1)
      WM   = RWL * W(I,J,K) + RWR * W(I,J,K+1)
      HM   = RWL * ( &
      &        0.5 * (U(I,J,K  )**2 + V(I,J,K  )**2 + W(I,J,K  )**2) &
      &      + GAMMA * RG / (GAMMA - 1.0) * T(I,J,K  ) &
      &    ) &
      &    + RWR * ( &
      &        0.5 * (U(I,J,K+1)**2 + V(I,J,K+1)**2 + W(I,J,K+1)**2) &
      &      + GAMMA * RG / (GAMMA - 1.0) * T(I,J,K+1) &
      &    )
      CM   = SQRT((GAMMA - 1.0) * (HM - 0.5 * (UM**2 + VM**2 + WM**2)))
      ! 対角化
      THI  = AKX * UM + AKY * VM + AKZ * WM
      PHI  = 0.5 * (GAMMA - 1.0) * (UM**2 + VM**2 + WM**2)
      BETA = 1.0 / (2.0 * CM**2)
      TIL  = SQRT(AKX**2 + AKY**2 + AKZ**2)
      AKXT = AKX / TIL
      AKYT = AKY / TIL
      AKZT = AKZ / TIL
      THIT = AKXT * UM + AKYT * VM + AKZT * WM
      RR(K,1,1) = AKXT
      RR(K,2,1) = AKXT * UM
      RR(K,3,1) = AKXT * VM + AKZT
      RR(K,4,1) = AKXT * WM - AKYT
      RR(K,5,1) = AKXT * PHI / (GAMMA - 1.0) + AKZT * VM - AKYT * WM
      RR(K,1,2) = AKYT
      RR(K,2,2) = AKYT * UM - AKZT
      RR(K,3,2) = AKYT * VM
      RR(K,4,2) = AKYT * WM + AKXT
      RR(K,5,2) = AKYT * PHI / (GAMMA - 1.0) + AKXT * WM - AKZT * UM
      RR(K,1,3) = AKZT
      RR(K,2,3) = AKZT * UM + AKYT
      RR(K,3,3) = AKZT * VM - AKXT
      RR(K,4,3) = AKZT * WM
      RR(K,5,3) = AKZT * PHI / (GAMMA - 1.0) + AKYT * UM - AKXT * VM
      RR(K,1,4) = 1.0
      RR(K,2,4) = UM + AKXT * CM
      RR(K,3,4) = VM + AKYT * CM
      RR(K,4,4) = WM + AKZT * CM
      RR(K,5,4) = (PHI + CM**2) / (GAMMA - 1.0) + CM * THIT
      RR(K,1,5) = 1.0
      RR(K,2,5) = UM - AKXT * CM
      RR(K,3,5) = VM - AKYT * CM
      RR(K,4,5) = WM - AKZT * CM
      RR(K,5,5) = (PHI + CM**2) / (GAMMA - 1.0) - CM * THIT
      RI(K,1,1) = AKXT - AKXT * PHI / CM**2 + AKYT * WM - AKZT * VM
      RI(K,2,1) = AKYT - AKYT * PHI / CM**2 + AKZT * UM - AKXT * WM
      RI(K,3,1) = AKZT - AKZT * PHI / CM**2 + AKXT * VM - AKYT * UM
      RI(K,4,1) = BETA * (PHI - CM * THIT)
      RI(K,5,1) = BETA * (PHI + CM * THIT)
      RI(K,1,2) = (GAMMA - 1.0) * AKXT * UM / CM**2
      RI(K,2,2) = (GAMMA - 1.0) * AKYT * UM / CM**2 - AKZT
      RI(K,3,2) = (GAMMA - 1.0) * AKZT * UM / CM**2 + AKYT
      RI(K,4,2) = BETA * (AKXT * CM - (GAMMA - 1.0) * UM)
      RI(K,5,2) =-BETA * (AKXT * CM + (GAMMA - 1.0) * UM)
      RI(K,1,3) = (GAMMA - 1.0) * AKXT * VM / CM**2 + AKZT
      RI(K,2,3) = (GAMMA - 1.0) * AKYT * VM / CM**2
      RI(K,3,3) = (GAMMA - 1.0) * AKZT * VM / CM**2 - AKXT
      RI(K,4,3) = BETA * (AKYT * CM - (GAMMA - 1.0) * VM)
      RI(K,5,3) =-BETA * (AKYT * CM + (GAMMA - 1.0) * VM)
      RI(K,1,4) = (GAMMA - 1.0) * AKXT * WM / CM**2 - AKYT
      RI(K,2,4) = (GAMMA - 1.0) * AKYT * WM / CM**2 + AKXT
      RI(K,3,4) = (GAMMA - 1.0) * AKZT * WM / CM**2
      RI(K,4,4) = BETA * (AKZT * CM - (GAMMA - 1.0) * WM)
      RI(K,5,4) =-BETA * (AKZT * CM + (GAMMA - 1.0) * WM)
      RI(K,1,5) =-(GAMMA - 1.0) * AKXT / CM**2
      RI(K,2,5) =-(GAMMA - 1.0) * AKYT / CM**2
      RI(K,3,5) =-(GAMMA - 1.0) * AKZT / CM**2
      RI(K,4,5) = BETA * (GAMMA - 1.0)
      RI(K,5,5) = BETA * (GAMMA - 1.0)
      EIGM(K,1) = THI
      EIGM(K,2) = THI
      EIGM(K,3) = THI
      EIGM(K,4) = THI + CM * TIL
      EIGM(K,5) = THI - CM * TIL
      DO L = 6, LE
        EIGM(K,L) = THI
      ENDDO
      DO L = LS, LE
        DELTA(K,L) = MAX( &
        &              0.0, &
        &              EIGM(K,L) - EIG(K,L), EIG(K+1,L) - EIGM(K,L) &
        &            ) + DEL * CM * TIL
      ENDDO
    ENDIF
    ENDDO
    ! Delta u
    DO L = LS, LE
    DO K = KS, KE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K+1,1) .GT. 0.0) THEN
      DU(K,L) = - QH(I,J,K  ,L) * AJA(I,J,K  ) &
      &         + QH(I,J,K+1,L) * AJA(I,J,K+1)
      DU(K,L) = DU(K,L) * 2.0 / (AJA(I,J,K) + AJA(I,J,K+1))
    ENDIF
    ENDDO
    ENDDO
    ! R^-1 Delta u
    DO L = LS, LE
    DO K = KS, KE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K+1,1) .GT. 0.0) THEN
      IF(L .LE. 5) THEN
        DUR(K,L) = RI(K,L,1) * DU(K,1) + RI(K,L,2) * DU(K,2) &
        &        + RI(K,L,3) * DU(K,3) + RI(K,L,4) * DU(K,4) &
        &        + RI(K,L,5) * DU(K,5)
      ELSE
        DUR(K,L) = DU(K,L)
      ENDIF
    ENDIF
    ENDDO
    ENDDO
    ! 二次精度TVDスキーム ++++++++++++++++++++++++++++++++++++++++++++++
    IF(Order .EQ. 2) THEN
      ! g
      DO L = LS, LE
      K = KS
      IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K+1,1) .GT. 0.0) THEN
        GG(K,L) = DUR(K,L)
      ENDIF
      DO K = KS + 1, KE - 1
      IF(QH(I,J,K,1) .LE. 0.0) CYCLE
      IF(QH(I,J,K-1,1) .GT. 0.0 .AND. QH(I,J,K+1,1) .GT. 0.0) THEN
        GG(K,L) = minmod2(DUR(K,L), DUR(K-1,L))
      ELSEIF(QH(I,J,K-1,1) .GT. 0.0) THEN
        GG(K,L) = DUR(K-1,L)
      ELSEIF(QH(I,J,K+1,1) .GT. 0.0) THEN
        GG(K,L) = DUR(K,L)
      ENDIF
      ENDDO
      K = KE
      IF(QH(I,J,K-1,1) .GT. 0.0 .AND. QH(I,J,K,1) .GT. 0.0) THEN
        GG(K,L) = DUR(K-1,L)
      ENDIF
      ENDDO
      ! phi
      DO L = LS, LE
      DO K = KS, KE - 1
      IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K+1,1) .GT. 0.0) THEN
        FSIGMA = Func_sigma(Func_psi(EIGM(K,L), DELTA(K,L)))
        FGAMMA = Func_gamma(FSIGMA, GG(K,L), GG(K+1,L), DUR(K,L))
        PHIM(K,L) = FSIGMA * (GG(K,L) + GG(K+1,L)) &
        &         - Func_psi(EIGM(K,L) + FGAMMA, DELTA(K,L)) * DUR(K,L)
      ENDIF
      ENDDO
      ENDDO
      ! 数値流束
      DO L = LS, LE
      DO K = KS, KE - 1
      IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K+1,1) .GT. 0.0) THEN
        IF(L .LE. 5) THEN
          GT(I,J,K,L) = GT(I,J,K,L) + 0.5 * ( RR(K,L,1) * PHIM(K,1) &
          &                                 + RR(K,L,2) * PHIM(K,2) &
          &                                 + RR(K,L,3) * PHIM(K,3) &
          &                                 + RR(K,L,4) * PHIM(K,4) &
          &                                 + RR(K,L,5) * PHIM(K,5) )
        ELSE
          GT(I,J,K,L) = GT(I,J,K,L) + 0.5 * PHIM(K,L)
        ENDIF
      ENDIF
      ENDDO
      ENDDO
      CYCLE
    ENDIF
    ! 一次精度風上差分 +++++++++++++++++++++++++++++++++++++++++++++++++
    ! Delta f^+-
    DO L = LS, LE
    DO K = KS, KE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K+1,1) .GT. 0.0) THEN
      ABSA  = Func_psi(EIGM(K,L), DELTA(K,L))
      AP    = 0.5 * (EIGM(K,L) + ABSA)
      AM    = 0.5 * (EIGM(K,L) - ABSA)
      DFP(K,L) = AP * DUR(K,L)
      DFM(K,L) = AM * DUR(K,L)
    ENDIF
    ENDDO
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO K = KS, KE - 1
    IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K+1,1) .GT. 0.0) THEN
      IF(L .LE. 5) THEN
        GT(I,J,K,L) = GT(I,J,K,L) - 0.5 * ( RR(K,L,1) * DFP(K,1) &
        &                                 + RR(K,L,2) * DFP(K,2) &
        &                                 + RR(K,L,3) * DFP(K,3) &
        &                                 + RR(K,L,4) * DFP(K,4) &
        &                                 + RR(K,L,5) * DFP(K,5) &
        &                                 - RR(K,L,1) * DFM(K,1) &
        &                                 - RR(K,L,2) * DFM(K,2) &
        &                                 - RR(K,L,3) * DFM(K,3) &
        &                                 - RR(K,L,4) * DFM(K,4) &
        &                                 - RR(K,L,5) * DFM(K,5) )
      ELSE
        GT(I,J,K,L) = GT(I,J,K,L) - 0.5 * (DFP(K,L) - DFM(K,L))
      ENDIF
    ENDIF
    ENDDO
    ENDDO
    IF(Order .EQ. 1) CYCLE
    ! 三次&四次精度TVDスキーム +++++++++++++++++++++++++++++++++++++++++
    IF(Order .EQ. 4) THEN
      ! Delta^3 bar(f)^+-
      DO L = LS, LE
      DO K = KS + 1, KE - 2
      IF( QH(I,J,K-1,1) .GT. 0.0 .AND. QH(I,J,K  ,1) .GT. 0.0 .AND. &
      &   QH(I,J,K+1,1) .GT. 0.0 .AND. QH(I,J,K+2,1) .GT. 0.0 ) THEN
        D3FP1 = minmod3(DFP(K-1,L), b2 * DFP(K  ,L), b2 * DFP(K+1,L))
        D3FP2 = minmod3(DFP(K  ,L), b2 * DFP(K+1,L), b2 * DFP(K-1,L))
        D3FP3 = minmod3(DFP(K+1,L), b2 * DFP(K-1,L), b2 * DFP(K  ,L))
        D3FM1 = minmod3(DFM(K-1,L), b2 * DFM(K  ,L), b2 * DFM(K+1,L))
        D3FM2 = minmod3(DFM(K  ,L), b2 * DFM(K+1,L), b2 * DFM(K-1,L))
        D3FM3 = minmod3(DFM(K+1,L), b2 * DFM(K-1,L), b2 * DFM(K  ,L))
        D3FP(K,L) = D3FP1 - 2.0 * D3FP2 + D3FP3
        D3FM(K,L) = D3FM1 - 2.0 * D3FM2 + D3FM3
      ENDIF
      ENDDO
      ENDDO
      ! Delta^asterisk f^+-
      DO L = LS, LE
      DO K = KS + 1, KE - 2
      IF( QH(I,J,K-1,1) .GT. 0.0 .AND. QH(I,J,K  ,1) .GT. 0.0 .AND. &
      &   QH(I,J,K+1,1) .GT. 0.0 .AND. QH(I,J,K+2,1) .GT. 0.0 ) THEN
        DFP(K,L) = DFP(K,L) - D3FP(K,L) / 6.0
        DFM(K,L) = DFM(K,L) - D3FM(K,L) / 6.0
      ENDIF
      ENDDO
      ENDDO
    ENDIF
    ! Delta bar(f)^+-, Delta tilde(f)^+-
    DO L = LS, LE
    DO K = KS + 1, KE - 2
    IF( QH(I,J,K-1,1) .GT. 0.0 .AND. QH(I,J,K  ,1) .GT. 0.0 .AND. &
    &   QH(I,J,K+1,1) .GT. 0.0 .AND. QH(I,J,K+2,1) .GT. 0.0 ) THEN
      DFBP(K-1,L) = minmod2(DFP(K-1,L), b1 * DFP(K  ,L))
      DFBM(K  ,L) = minmod2(DFM(K  ,L), b1 * DFM(K+1,L))
      DFTP(K  ,L) = minmod2(DFP(K  ,L), b1 * DFP(K-1,L))
      DFTM(K+1,L) = minmod2(DFM(K+1,L), b1 * DFM(K  ,L))
    ENDIF
    ENDDO
    ENDDO
    ! 数値流束
    DO L = LS, LE
    DO K = KS + 1, KE - 2
    IF( QH(I,J,K-1,1) .GT. 0.0 .AND. QH(I,J,K  ,1) .GT. 0.0 .AND. &
    &   QH(I,J,K+1,1) .GT. 0.0 .AND. QH(I,J,K+2,1) .GT. 0.0 ) THEN
      IF(L .LE. 5) THEN
        GT(I,J,K,L) = GT(I,J,K,L) &
        &           + ( 1.0 * RR(K-1,L,1) * DFBP(K-1,1) &
        &             + 1.0 * RR(K-1,L,2) * DFBP(K-1,2) &
        &             + 1.0 * RR(K-1,L,3) * DFBP(K-1,3) &
        &             + 1.0 * RR(K-1,L,4) * DFBP(K-1,4) &
        &             + 1.0 * RR(K-1,L,5) * DFBP(K-1,5) &
        &             + 2.0 * RR(K  ,L,1) * DFTP(K  ,1) &
        &             + 2.0 * RR(K  ,L,2) * DFTP(K  ,2) &
        &             + 2.0 * RR(K  ,L,3) * DFTP(K  ,3) &
        &             + 2.0 * RR(K  ,L,4) * DFTP(K  ,4) &
        &             + 2.0 * RR(K  ,L,5) * DFTP(K  ,5) &
        &             - 2.0 * RR(K  ,L,1) * DFBM(K  ,1) &
        &             - 2.0 * RR(K  ,L,2) * DFBM(K  ,2) &
        &             - 2.0 * RR(K  ,L,3) * DFBM(K  ,3) &
        &             - 2.0 * RR(K  ,L,4) * DFBM(K  ,4) &
        &             - 2.0 * RR(K  ,L,5) * DFBM(K  ,5) &
        &             - 1.0 * RR(K+1,L,1) * DFTM(K+1,1) &
        &             - 1.0 * RR(K+1,L,2) * DFTM(K+1,2) &
        &             - 1.0 * RR(K+1,L,3) * DFTM(K+1,3) &
        &             - 1.0 * RR(K+1,L,4) * DFTM(K+1,4) &
        &             - 1.0 * RR(K+1,L,5) * DFTM(K+1,5) &
        &           ) / 6.0
      ELSE
        GT(I,J,K,L) = GT(I,J,K,L) &
        &           + ( 1.0 * DFBP(K-1,L) + 2.0 * DFTP(K  ,L) &
        &             - 2.0 * DFBM(K  ,L) - 1.0 * DFTM(K+1,L) ) / 6.0
      ENDIF
    ENDIF
    ENDDO
    ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE TVDZE
! 定義終了 *************************************************************
END SUBROUTINE Convection3D
!***********************************************************************
!**** 関数psi(エントロピー修正)                                     ****
!***********************************************************************
REAL FUNCTION Func_psi(x, delta)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN) :: x, delta
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL :: absx
  ! 処理開始 ***********************************************************
  absx = ABS(x)
  IF(absx .GE. delta) THEN
    Func_psi = absx
  ELSE
    Func_psi = (x**2 + delta**2) / (2.0 * delta)
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
! 定義終了 *************************************************************
END FUNCTION Func_psi
!***********************************************************************
!**** 関数sigma                                                     ****
!***********************************************************************
REAL FUNCTION Func_sigma(c)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN) :: c
  ! 処理開始 ***********************************************************
  Func_sigma = 0.5 * c
  ! 処理終了 ***********************************************************
  RETURN
! 定義終了 *************************************************************
END FUNCTION Func_sigma
!***********************************************************************
!**** 関数gamma                                                     ****
!***********************************************************************
REAL FUNCTION Func_gamma(sigma, gl, gr, Deltau)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN) :: sigma, gl, gr, Deltau
  ! 処理開始 ***********************************************************
  IF(Deltau .NE. 0.0) THEN
    Func_gamma = sigma * (gr - gl) / Deltau
  ELSE
    Func_gamma = 0.0
  ENDIF
  ! 処理終了 ***********************************************************
  RETURN
! 定義終了 *************************************************************
END FUNCTION Func_gamma
!***********************************************************************
!**** 二変数minmod関数の定義                                        ****
!***********************************************************************
REAL FUNCTION minmod2(x, y)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN) :: x, y
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL :: s
  ! 処理開始 ***********************************************************
  s = SIGN(1.0, x)
  minmod2 = s * MAX(0.0, MIN(ABS(x), s * y))
  ! 処理終了 ***********************************************************
  RETURN
END FUNCTION minmod2
!***********************************************************************
!**** 三変数minmod関数の定義                                        ****
!***********************************************************************
REAL FUNCTION minmod3(x, y, z)
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN) :: x, y, z
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL :: s
  ! 処理開始 ***********************************************************
  s = SIGN(1.0, x)
  minmod3 = s * MAX(0.0, MIN(ABS(x), s * y, s * z))
  ! 処理終了 ***********************************************************
  RETURN
END FUNCTION minmod3
!***********************************************************************
!**** 回転座標系の見かけの力を計算(二次元)                          ****
!***********************************************************************
SUBROUTINE RotationForce2D( &
&            OmegaZ, &
&            IS, IE, JS, JE, LS, LE, &
&            X, Y, AJA, RHO, U, V, &
&            DQR &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: OmegaZ
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: X(IS: IE, JS: JE), Y(IS: IE, JS: JE)
  REAL,    INTENT(IN)  :: AJA(IS: IE, JS: JE)
  REAL,    INTENT(IN)  :: RHO(IS: IE, JS: JE)
  REAL,    INTENT(IN)  :: U(IS: IE, JS: JE), V(IS: IE, JS: JE)
  REAL,    INTENT(OUT) :: DQR(IS: IE, JS: JE, LS: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: QH1
  REAL    :: CoriolisX, CoriolisY
  REAL    :: CentrifugalX, CentrifugalY
  ! 処理開始 ***********************************************************
  DQR = 0.0
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, QH1, CoriolisX, CoriolisY, CentrifugalX, CentrifugalY &
  !$)
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J) .GT. 0.0) THEN
    QH1          = RHO(I,J) / AJA(I,J)
    CoriolisX    =-2.0 * OmegaZ * V(I,J)
    CoriolisY    = 2.0 * OmegaZ * U(I,J)
    CentrifugalX =-OmegaZ**2 * X(I,J)
    CentrifugalY =-OmegaZ**2 * Y(I,J)
    DQR(I,J,2) =-QH1 * (CoriolisX + CentrifugalX)
    DQR(I,J,3) =-QH1 * (CoriolisY + CentrifugalY)
    DQR(I,J,4) =-QH1 * (CentrifugalX * U(I,J) + CentrifugalY * V(I,J))
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE RotationForce2D
!***********************************************************************
!**** 回転座標系の見かけの力を計算(三次元)                          ****
!***********************************************************************
SUBROUTINE RotationForce3D( &
&            OmegaX, OmegaY, OmegaZ, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            X, Y, Z, AJA, RHO, U, V, W, &
&            DQR &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: OmegaX, OmegaY, OmegaZ
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: X(IS: IE, JS: JE, KS: KE), &
  &                       Y(IS: IE, JS: JE, KS: KE), &
  &                       Z(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)  :: AJA(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)  :: RHO(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)  :: U(IS: IE, JS: JE, KS: KE), &
  &                       V(IS: IE, JS: JE, KS: KE), &
  &                       W(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(OUT) :: DQR(IS: IE, JS: JE, KS: KE, LS: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: QH1
  REAL    :: CoriolisX, CoriolisY, CoriolisZ
  REAL    :: CentrifugalX, CentrifugalY, CentrifugalZ
  ! 処理開始 ***********************************************************
  DQR = 0.0
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, QH1, &
  !$  CoriolisX, CoriolisY, CoriolisZ, &
  !$  CentrifugalX, CentrifugalY, CentrifugalZ &
  !$)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J,K) .GT. 0.0) THEN
    QH1          = RHO(I,J,K) / AJA(I,J,K)
    CoriolisX    = 2.0 * (OmegaY * W(I,J,K) - OmegaZ * V(I,J,K))
    CoriolisY    = 2.0 * (OmegaZ * U(I,J,K) - OmegaX * W(I,J,K))
    CoriolisZ    = 2.0 * (OmegaX * V(I,J,K) - OmegaY * U(I,J,K))
    CentrifugalX = OmegaY * (OmegaX * Y(I,J,K) - OmegaY * X(I,J,K)) &
    &            - OmegaZ * (OmegaZ * X(I,J,K) - OmegaX * Z(I,J,K))
    CentrifugalY = OmegaZ * (OmegaY * Z(I,J,K) - OmegaZ * Y(I,J,K)) &
    &            - OmegaX * (OmegaX * Y(I,J,K) - OmegaY * X(I,J,K))
    CentrifugalZ = OmegaX * (OmegaZ * X(I,J,K) - OmegaX * Z(I,J,K)) &
    &            - OmegaY * (OmegaY * Z(I,J,K) - OmegaZ * Y(I,J,K))
    DQR(I,J,K,2) =-QH1 * (CoriolisX + CentrifugalX)
    DQR(I,J,K,3) =-QH1 * (CoriolisY + CentrifugalY)
    DQR(I,J,K,4) =-QH1 * (CoriolisZ + CentrifugalZ)
    DQR(I,J,K,5) =-QH1 * ( CentrifugalX * U(I,J,K) &
    &                    + CentrifugalY * V(I,J,K) &
    &                    + CentrifugalZ * W(I,J,K) )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE RotationForce3D
!***********************************************************************
!**** 拡散項(粘性項)の計算(スタガード格子上の点を使用)              ****
!**** 計算対象 : 単相, 二次元, 圧縮性, 層流                         ****
!**** 計算スキーム : 二次精度中心差分                               ****
!***********************************************************************
SUBROUTINE Viscosity2D( &
&            RG, GAMMA, PR, &
&            IS, IE, JS, JE, LS, LE, &
&            XIX, XIY, ETX, ETY, AJA, &
&            QH, U, V, T, AMU, &
&            DQD &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: RG, GAMMA, PR
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE), XIY(IS:IE, JS:JE), &
  &                       ETX(IS:IE, JS:JE), ETY(IS:IE, JS:JE), &
  &                       AJA(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: QH(IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(IN)  :: U(IS:IE, JS:JE), &
  &                       V(IS:IE, JS:JE), &
  &                       T(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  REAL, ALLOCATABLE :: RH(:, :, :), SH(:, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 4) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS:IE, JS:JE, LS:LE))
  ALLOCATE(SH(IS:IE, JS:JE, LS:LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH (:,:,:) = 0.0
  SH (:,:,:) = 0.0
  DQD(:,:,:) = 0.0
  ! 各方向の空間差分 +++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL VISXI
  CALL VISET
  ! 拡散項ベクトルの空間差分 +++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, L)
  DO L = LS, LE
  !$OMP DO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (QH(I  ,J  ,1) .GT. 0.0) .AND. &
  &   (QH(I-1,J  ,1) .GT. 0.0) .AND. (QH(I+1,J  ,1) .GT. 0.0) .AND. &
  &   (QH(I  ,J-1,1) .GT. 0.0) .AND. (QH(I  ,J+1,1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,L) = (-RH(I-1,J  ,L) + RH(I,J,L)) &
    &          + (-SH(I  ,J-1,L) + SH(I,J,L))
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 拡散項(粘性項)の計算(XI方向)                                  ****
!***********************************************************************
SUBROUTINE VISXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL :: UM, VM, TM
  REAL :: UXIM, VXIM, TXIM, &
  &       UETM, VETM, TETM
  REAL :: UXM, VXM, TXM, &
  &       UYM, VYM, TYM
  REAL :: AMUM
  REAL :: DELVM, TAUXX, TAUYY, TAUXY
  REAL :: QX, QY
  REAL :: R4, S4
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, UM, VM, TM, &
  !$  UXIM, VXIM, TXIM, UETM, VETM, TETM, &
  !$  UXM, VXM, TXM, UYM, VYM, TYM, &
  !$  AMUM, DELVM, TAUXX, TAUYY, TAUXY, QX, QY, R4, S4 &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I+1,J,1) .GT. 0.0) THEN
      ! XI方向の平均量の計算 -------------------------------------------
      XIXM  = 0.5 * (XIX(I,J) + XIX(I+1,J))
      XIYM  = 0.5 * (XIY(I,J) + XIY(I+1,J))
      ETXM  = 0.5 * (ETX(I,J) + ETX(I+1,J))
      ETYM  = 0.5 * (ETY(I,J) + ETY(I+1,J))
      AJAM  = 0.5 * (AJA(I,J) + AJA(I+1,J))
      UM    = 0.5 * (  U(I,J) +   U(I+1,J))
      VM    = 0.5 * (  V(I,J) +   V(I+1,J))
      TM    = 0.5 * (  T(I,J) +   T(I+1,J))
      ! x,y微分の計算(XI方向の平均量) ----------------------------------
      UXIM = (-U(I,J) + U(I+1,J))
      VXIM = (-V(I,J) + V(I+1,J))
      TXIM = (-T(I,J) + T(I+1,J))
      UETM = 0.5 * ( - 0.5 * (U(I,J-1) + U(I+1,J-1)) &
      &              + 0.5 * (U(I,J+1) + U(I+1,J+1)) )
      VETM = 0.5 * ( - 0.5 * (V(I,J-1) + V(I+1,J-1)) &
      &              + 0.5 * (V(I,J+1) + V(I+1,J+1)) )
      TETM = 0.5 * ( - 0.5 * (T(I,J-1) + T(I+1,J-1)) &
      &              + 0.5 * (T(I,J+1) + T(I+1,J+1)) )
      UXM  = XIXM * UXIM + ETXM * UETM
      UYM  = XIYM * UXIM + ETYM * UETM
      VXM  = XIXM * VXIM + ETXM * VETM
      VYM  = XIYM * VXIM + ETYM * VETM
      TXM  = XIXM * TXIM + ETXM * TETM
      TYM  = XIYM * TXIM + ETYM * TETM
      ! 粘性係数のXI方向の平均量 ---------------------------------------
      AMUM = 0.5 * (AMU(I,J) + AMU(I+1,J))
      ! 拡散項ベクトルの計算 -------------------------------------------
      DELVM = UXM + VYM
      TAUXX = AMUM * (2.0 * UXM - 2.0 / 3.0 * DELVM)
      TAUYY = AMUM * (2.0 * VYM - 2.0 / 3.0 * DELVM)
      TAUXY = AMUM * (UYM + VXM)
      QX    =-(AMUM / PR) / (GAMMA - 1.0) * GAMMA * RG * TXM
      QY    =-(AMUM / PR) / (GAMMA - 1.0) * GAMMA * RG * TYM
      R4    = TAUXX * UM + TAUXY * VM - QX
      S4    = TAUXY * UM + TAUYY * VM - QY
      RH(I,J,1) = 0.0
      RH(I,J,2) = (XIXM * TAUXX + XIYM * TAUXY) / AJAM
      RH(I,J,3) = (XIXM * TAUXY + XIYM * TAUYY) / AJAM
      RH(I,J,4) = (XIXM * R4    + XIYM * S4   ) / AJAM
    ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISXI
!***********************************************************************
!**** 拡散項(粘性項)の計算(ET方向)                                  ****
!***********************************************************************
SUBROUTINE VISET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL :: UM, VM, TM
  REAL :: UXIM, VXIM, TXIM, &
  &       UETM, VETM, TETM
  REAL :: UXM, VXM, TXM, &
  &       UYM, VYM, TYM
  REAL :: AMUM
  REAL :: DELVM, TAUXX, TAUYY, TAUXY
  REAL :: QX, QY
  REAL :: R4, S4
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, UM, VM, TM, &
  !$  UXIM, VXIM, TXIM, UETM, VETM, TETM, &
  !$  UXM, VXM, TXM, UYM, VYM, TYM, &
  !$  AMUM, DELVM, TAUXX, TAUYY, TAUXY, QX, QY, R4, S4 &
  !$)
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
    IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J+1,1) .GT. 0.0) THEN
      ! ET方向の平均量の計算 -------------------------------------------
      XIXM  = 0.5 * (XIX(I,J) + XIX(I,J+1))
      XIYM  = 0.5 * (XIY(I,J) + XIY(I,J+1))
      ETXM  = 0.5 * (ETX(I,J) + ETX(I,J+1))
      ETYM  = 0.5 * (ETY(I,J) + ETY(I,J+1))
      AJAM  = 0.5 * (AJA(I,J) + AJA(I,J+1))
      UM    = 0.5 * (  U(I,J) +   U(I,J+1))
      VM    = 0.5 * (  V(I,J) +   V(I,J+1))
      TM    = 0.5 * (  T(I,J) +   T(I,J+1))
      ! x,y微分の計算(ET方向の平均量) ----------------------------------
      UXIM = 0.5 * ( - 0.5 * (U(I-1,J) + U(I-1,J+1)) &
      &              + 0.5 * (U(I+1,J) + U(I+1,J+1)) )
      VXIM = 0.5 * ( - 0.5 * (V(I-1,J) + V(I-1,J+1)) &
      &              + 0.5 * (V(I+1,J) + V(I+1,J+1)) )
      TXIM = 0.5 * ( - 0.5 * (T(I-1,J) + T(I-1,J+1)) &
      &              + 0.5 * (T(I+1,J) + T(I+1,J+1)) )
      UETM = (-U(I,J) + U(I,J+1))
      VETM = (-V(I,J) + V(I,J+1))
      TETM = (-T(I,J) + T(I,J+1))
      UXM  = XIXM * UXIM + ETXM * UETM
      UYM  = XIYM * UXIM + ETYM * UETM
      VXM  = XIXM * VXIM + ETXM * VETM
      VYM  = XIYM * VXIM + ETYM * VETM
      TXM  = XIXM * TXIM + ETXM * TETM
      TYM  = XIYM * TXIM + ETYM * TETM
      ! 粘性係数のET方向の平均量 ---------------------------------------
      AMUM = 0.5 * (AMU(I,J) + AMU(I,J+1))
      ! 拡散項ベクトルの計算 -------------------------------------------
      DELVM = UXM + VYM
      TAUXX = AMUM * (2.0 * UXM - 2.0 / 3.0 * DELVM)
      TAUYY = AMUM * (2.0 * VYM - 2.0 / 3.0 * DELVM)
      TAUXY = AMUM * (UYM + VXM)
      QX    =-(AMUM / PR) / (GAMMA-1.0) * GAMMA * RG * TXM
      QY    =-(AMUM / PR) / (GAMMA-1.0) * GAMMA * RG * TYM
      R4    = TAUXX * UM + TAUXY * VM - QX
      S4    = TAUXY * UM + TAUYY * VM - QY
      SH(I,J,1) = 0.0
      SH(I,J,2) = (ETXM * TAUXX + ETYM * TAUXY) / AJAM
      SH(I,J,3) = (ETXM * TAUXY + ETYM * TAUYY) / AJAM
      SH(I,J,4) = (ETXM * R4    + ETYM * S4   ) / AJAM
    ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISET
! 定義終了 *************************************************************
END SUBROUTINE Viscosity2D
!***********************************************************************
!**** 拡散項(粘性項)の計算(スタガード格子上の点を使用)              ****
!**** 計算対象 : 単相, 三次元, 圧縮性, 層流                         ****
!**** 計算スキーム : 二次精度中心差分                               ****
!***********************************************************************
SUBROUTINE Viscosity3D( &
&            RG, GAMMA, PR, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, AJA, &
&            QH, U, V, W, T, AMU, &
&            DQD &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: RG, GAMMA, PR
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE, KS:KE), &
  &                       XIY(IS:IE, JS:JE, KS:KE), &
  &                       XIZ(IS:IE, JS:JE, KS:KE), &
  &                       ETX(IS:IE, JS:JE, KS:KE), &
  &                       ETY(IS:IE, JS:JE, KS:KE), &
  &                       ETZ(IS:IE, JS:JE, KS:KE), &
  &                       ZEX(IS:IE, JS:JE, KS:KE), &
  &                       ZEY(IS:IE, JS:JE, KS:KE), &
  &                       ZEZ(IS:IE, JS:JE, KS:KE), &
  &                       AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: QH(IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(IN)  :: U(IS:IE, JS:JE, KS:KE), &
  &                       V(IS:IE, JS:JE, KS:KE), &
  &                       W(IS:IE, JS:JE, KS:KE), &
  &                       T(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, KS:KE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  REAL, ALLOCATABLE :: RH(:, :, :, :), SH(:, :, :, :), TH(:, :, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 5) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS:IE, JS:JE, KS:KE, LS:LE))
  ALLOCATE(SH(IS:IE, JS:JE, KS:KE, LS:LE))
  ALLOCATE(TH(IS:IE, JS:JE, KS:KE, LS:LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH (:,:,:,:) = 0.0
  SH (:,:,:,:) = 0.0
  TH (:,:,:,:) = 0.0
  DQD(:,:,:,:) = 0.0
  ! 各方向の空間差分 +++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL VISXI
  CALL VISET
  CALL VISZE
  ! 拡散項ベクトルの空間差分 +++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( ( QH(I  ,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I-1,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I+1,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J-1,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J+1,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J  ,K-1,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J  ,K+1,1) .GT. 0.0 ) &
  & ) THEN
    DQD(I,J,K,L) = ( -RH(I-1,J  ,K  ,L) + RH(I,J,K,L) ) &
    &            + ( -SH(I  ,J-1,K  ,L) + SH(I,J,K,L) ) &
    &            + ( -TH(I  ,J  ,K-1,L) + TH(I,J,K,L) )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 拡散項(粘性項)の計算(XI方向)                                  ****
!***********************************************************************
SUBROUTINE VISXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL :: UM, VM, WM, TM
  REAL :: UXIM, VXIM, WXIM, TXIM, &
  &       UETM, VETM, WETM, TETM, &
  &       UZEM, VZEM, WZEM, TZEM
  REAL :: UXM, VXM, WXM, TXM, &
  &       UYM, VYM, WYM, TYM, &
  &       UZM, VZM, WZM, TZM
  REAL :: AMUM
  REAL :: DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX
  REAL :: QX, QY, QZ
  REAL :: R5, S5, T5
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  UM, VM, WM, TM, &
  !$  UXIM, VXIM, WXIM, TXIM, &
  !$  UETM, VETM, WETM, TETM, &
  !$  UZEM, VZEM, WZEM, TZEM, &
  !$  UXM, VXM, WXM, TXM, &
  !$  UYM, VYM, WYM, TYM, &
  !$  UZM, VZM, WZM, TZM, &
  !$  AMUM, DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, &
  !$  QX, QY, QZ, R5, S5, T5 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
    IF( QH(I,J,K,1) .GT. 0.0 .AND. QH(I+1,J,K,1) .GT. 0.0 ) THEN
      ! XI方向の平均量の計算 -------------------------------------------
      XIXM  = 0.5 * ( XIX(I,J,K) + XIX(I+1,J,K) )
      XIYM  = 0.5 * ( XIY(I,J,K) + XIY(I+1,J,K) )
      XIZM  = 0.5 * ( XIZ(I,J,K) + XIZ(I+1,J,K) )
      ETXM  = 0.5 * ( ETX(I,J,K) + ETX(I+1,J,K) )
      ETYM  = 0.5 * ( ETY(I,J,K) + ETY(I+1,J,K) )
      ETZM  = 0.5 * ( ETZ(I,J,K) + ETZ(I+1,J,K) )
      ZEXM  = 0.5 * ( ZEX(I,J,K) + ZEX(I+1,J,K) )
      ZEYM  = 0.5 * ( ZEY(I,J,K) + ZEY(I+1,J,K) )
      ZEZM  = 0.5 * ( ZEZ(I,J,K) + ZEZ(I+1,J,K) )
      AJAM  = 0.5 * ( AJA(I,J,K) + AJA(I+1,J,K) )
      UM    = 0.5 * (   U(I,J,K) +   U(I+1,J,K) )
      VM    = 0.5 * (   V(I,J,K) +   V(I+1,J,K) )
      WM    = 0.5 * (   W(I,J,K) +   W(I+1,J,K) )
      TM    = 0.5 * (   T(I,J,K) +   T(I+1,J,K) )
      ! x,y,z微分の計算(XI方向の平均量) --------------------------------
      UXIM   = ( -  U(I,J,K) +  U(I+1,J,K) )
      VXIM   = ( -  V(I,J,K) +  V(I+1,J,K) )
      WXIM   = ( -  W(I,J,K) +  W(I+1,J,K) )
      TXIM   = ( -  T(I,J,K) +  T(I+1,J,K) )
      UETM   = 0.5 * ( - 0.5*(   U(I,J-1,K) +   U(I+1,J-1,K) ) &
      &                + 0.5*(   U(I,J+1,K) +   U(I+1,J+1,K) )  )
      VETM   = 0.5 * ( - 0.5*(   V(I,J-1,K) +   V(I+1,J-1,K) ) &
      &                + 0.5*(   V(I,J+1,K) +   V(I+1,J+1,K) )  )
      WETM   = 0.5 * ( - 0.5*(   W(I,J-1,K) +   W(I+1,J-1,K) ) &
      &                + 0.5*(   W(I,J+1,K) +   W(I+1,J+1,K) )  )
      TETM   = 0.5 * ( - 0.5*(   T(I,J-1,K) +   T(I+1,J-1,K) ) &
      &                + 0.5*(   T(I,J+1,K) +   T(I+1,J+1,K) )  )
      UZEM   = 0.5 * ( - 0.5*(   U(I,J,K-1) +   U(I+1,J,K-1) ) &
      &                + 0.5*(   U(I,J,K+1) +   U(I+1,J,K+1) )  )
      VZEM   = 0.5 * ( - 0.5*(   V(I,J,K-1) +   V(I+1,J,K-1) ) &
      &                + 0.5*(   V(I,J,K+1) +   V(I+1,J,K+1) )  )
      WZEM   = 0.5 * ( - 0.5*(   W(I,J,K-1) +   W(I+1,J,K-1) ) &
      &                + 0.5*(   W(I,J,K+1) +   W(I+1,J,K+1) )  )
      TZEM   = 0.5 * ( - 0.5*(   T(I,J,K-1) +   T(I+1,J,K-1) ) &
      &                + 0.5*(   T(I,J,K+1) +   T(I+1,J,K+1) )  )
      UXM   = XIXM *  UXIM + ETXM *  UETM + ZEXM *  UZEM
      UYM   = XIYM *  UXIM + ETYM *  UETM + ZEYM *  UZEM
      UZM   = XIZM *  UXIM + ETZM *  UETM + ZEZM *  UZEM
      VXM   = XIXM *  VXIM + ETXM *  VETM + ZEXM *  VZEM
      VYM   = XIYM *  VXIM + ETYM *  VETM + ZEYM *  VZEM
      VZM   = XIZM *  VXIM + ETZM *  VETM + ZEZM *  VZEM
      WXM   = XIXM *  WXIM + ETXM *  WETM + ZEXM *  WZEM
      WYM   = XIYM *  WXIM + ETYM *  WETM + ZEYM *  WZEM
      WZM   = XIZM *  WXIM + ETZM *  WETM + ZEZM *  WZEM
      TXM   = XIXM *  TXIM + ETXM *  TETM + ZEXM *  TZEM
      TYM   = XIYM *  TXIM + ETYM *  TETM + ZEYM *  TZEM
      TZM   = XIZM *  TXIM + ETZM *  TETM + ZEZM *  TZEM
      ! 粘性係数のXI方向の平均量 ---------------------------------------
      AMUM  = 0.5 * (AMU(I,J,K) + AMU(I+1,J,K))
      ! 拡散項ベクトルの計算 -------------------------------------------
      DELVM = UXM+VYM+WZM
      TAUXX = AMUM*(2.0*UXM-2.0/3.0*DELVM)
      TAUYY = AMUM*(2.0*VYM-2.0/3.0*DELVM)
      TAUZZ = AMUM*(2.0*WZM-2.0/3.0*DELVM)
      TAUXY = AMUM*(UYM+VXM)
      TAUYZ = AMUM*(VZM+WYM)
      TAUZX = AMUM*(WXM+UZM)
      QX    =-(AMUM/PR)/(GAMMA-1.0)*GAMMA*RG*TXM
      QY    =-(AMUM/PR)/(GAMMA-1.0)*GAMMA*RG*TYM
      QZ    =-(AMUM/PR)/(GAMMA-1.0)*GAMMA*RG*TZM
      R5    = TAUXX*UM + TAUXY*VM + TAUZX*WM-QX
      S5    = TAUXY*UM + TAUYY*VM + TAUYZ*WM-QY
      T5    = TAUZX*UM + TAUYZ*VM + TAUZZ*WM-QZ
      RH(I,J,K,1) = 0.0
      RH(I,J,K,2) = (XIXM*TAUXX + XIYM*TAUXY + XIZM*TAUZX)/AJAM
      RH(I,J,K,3) = (XIXM*TAUXY + XIYM*TAUYY + XIZM*TAUYZ)/AJAM
      RH(I,J,K,4) = (XIXM*TAUZX + XIYM*TAUYZ + XIZM*TAUZZ)/AJAM
      RH(I,J,K,5) = (XIXM*R5    + XIYM*S5    + XIZM*T5   )/AJAM
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISXI
!***********************************************************************
!**** 拡散項(粘性項)の計算(ET方向)                                  ****
!***********************************************************************
SUBROUTINE VISET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL :: UM, VM, WM, TM
  REAL :: UXIM, VXIM, WXIM, TXIM, &
  &       UETM, VETM, WETM, TETM, &
  &       UZEM, VZEM, WZEM, TZEM
  REAL :: UXM, VXM, WXM, TXM, &
  &       UYM, VYM, WYM, TYM, &
  &       UZM, VZM, WZM, TZM
  REAL :: AMUM
  REAL :: DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX
  REAL :: QX, QY, QZ
  REAL :: R5, S5, T5
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  UM, VM, WM, TM, &
  !$  UXIM, VXIM, WXIM, TXIM, &
  !$  UETM, VETM, WETM, TETM, &
  !$  UZEM, VZEM, WZEM, TZEM, &
  !$  UXM, VXM, WXM, TXM, &
  !$  UYM, VYM, WYM, TYM, &
  !$  UZM, VZM, WZM, TZM, &
  !$  AMUM, DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, &
  !$  QX, QY, QZ, R5, S5, T5 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
    IF( QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J+1,K,1) .GT. 0.0 ) THEN
      ! ET方向の平均量の計算 -------------------------------------------
      XIXM  = 0.5 * ( XIX(I,J,K) + XIX(I,J+1,K) )
      XIYM  = 0.5 * ( XIY(I,J,K) + XIY(I,J+1,K) )
      XIZM  = 0.5 * ( XIZ(I,J,K) + XIZ(I,J+1,K) )
      ETXM  = 0.5 * ( ETX(I,J,K) + ETX(I,J+1,K) )
      ETYM  = 0.5 * ( ETY(I,J,K) + ETY(I,J+1,K) )
      ETZM  = 0.5 * ( ETZ(I,J,K) + ETZ(I,J+1,K) )
      ZEXM  = 0.5 * ( ZEX(I,J,K) + ZEX(I,J+1,K) )
      ZEYM  = 0.5 * ( ZEY(I,J,K) + ZEY(I,J+1,K) )
      ZEZM  = 0.5 * ( ZEZ(I,J,K) + ZEZ(I,J+1,K) )
      AJAM  = 0.5 * ( AJA(I,J,K) + AJA(I,J+1,K) )
      UM    = 0.5 * (   U(I,J,K) +   U(I,J+1,K) )
      VM    = 0.5 * (   V(I,J,K) +   V(I,J+1,K) )
      WM    = 0.5 * (   W(I,J,K) +   W(I,J+1,K) )
      TM    = 0.5 * (   T(I,J,K) +   T(I,J+1,K) )
      ! x,y,z微分の計算(ET方向の平均量) --------------------------------
      UXIM   = 0.5 * ( - 0.5*(   U(I-1,J,K) +   U(I-1,J+1,K) ) &
      &                + 0.5*(   U(I+1,J,K) +   U(I+1,J+1,K) )  )
      VXIM   = 0.5 * ( - 0.5*(   V(I-1,J,K) +   V(I-1,J+1,K) ) &
      &                + 0.5*(   V(I+1,J,K) +   V(I+1,J+1,K) )  )
      WXIM   = 0.5 * ( - 0.5*(   W(I-1,J,K) +   W(I-1,J+1,K) ) &
      &                + 0.5*(   W(I+1,J,K) +   W(I+1,J+1,K) )  )
      TXIM   = 0.5 * ( - 0.5*(   T(I-1,J,K) +   T(I-1,J+1,K) ) &
      &                + 0.5*(   T(I+1,J,K) +   T(I+1,J+1,K) )  )
      UETM   = ( -  U(I,J,K) +  U(I,J+1,K) )
      VETM   = ( -  V(I,J,K) +  V(I,J+1,K) )
      WETM   = ( -  W(I,J,K) +  W(I,J+1,K) )
      TETM   = ( -  T(I,J,K) +  T(I,J+1,K) )
      UZEM   = 0.5 * ( - 0.5*(   U(I,J,K-1) +   U(I,J+1,K-1) ) &
      &                + 0.5*(   U(I,J,K+1) +   U(I,J+1,K+1) )  )
      VZEM   = 0.5 * ( - 0.5*(   V(I,J,K-1) +   V(I,J+1,K-1) ) &
      &                + 0.5*(   V(I,J,K+1) +   V(I,J+1,K+1) )  )
      WZEM   = 0.5 * ( - 0.5*(   W(I,J,K-1) +   W(I,J+1,K-1) ) &
      &                + 0.5*(   W(I,J,K+1) +   W(I,J+1,K+1) )  )
      TZEM   = 0.5 * ( - 0.5*(   T(I,J,K-1) +   T(I,J+1,K-1) ) &
      &                + 0.5*(   T(I,J,K+1) +   T(I,J+1,K+1) )  )
      UXM   = XIXM *  UXIM + ETXM *  UETM + ZEXM *  UZEM
      UYM   = XIYM *  UXIM + ETYM *  UETM + ZEYM *  UZEM
      UZM   = XIZM *  UXIM + ETZM *  UETM + ZEZM *  UZEM
      VXM   = XIXM *  VXIM + ETXM *  VETM + ZEXM *  VZEM
      VYM   = XIYM *  VXIM + ETYM *  VETM + ZEYM *  VZEM
      VZM   = XIZM *  VXIM + ETZM *  VETM + ZEZM *  VZEM
      WXM   = XIXM *  WXIM + ETXM *  WETM + ZEXM *  WZEM
      WYM   = XIYM *  WXIM + ETYM *  WETM + ZEYM *  WZEM
      WZM   = XIZM *  WXIM + ETZM *  WETM + ZEZM *  WZEM
      TXM   = XIXM *  TXIM + ETXM *  TETM + ZEXM *  TZEM
      TYM   = XIYM *  TXIM + ETYM *  TETM + ZEYM *  TZEM
      TZM   = XIZM *  TXIM + ETZM *  TETM + ZEZM *  TZEM
      ! 粘性係数のET方向の平均量 ---------------------------------------
      AMUM  = 0.5 * (AMU(I,J,K) + AMU(I,J+1,K))
      ! 拡散項ベクトルの計算 -------------------------------------------
      DELVM = UXM+VYM+WZM
      TAUXX = AMUM*(2.0*UXM-2.0/3.0*DELVM)
      TAUYY = AMUM*(2.0*VYM-2.0/3.0*DELVM)
      TAUZZ = AMUM*(2.0*WZM-2.0/3.0*DELVM)
      TAUXY = AMUM*(UYM+VXM)
      TAUYZ = AMUM*(VZM+WYM)
      TAUZX = AMUM*(WXM+UZM)
      QX    =-(AMUM/PR)/(GAMMA-1.0)*GAMMA*RG*TXM
      QY    =-(AMUM/PR)/(GAMMA-1.0)*GAMMA*RG*TYM
      QZ    =-(AMUM/PR)/(GAMMA-1.0)*GAMMA*RG*TZM
      R5    = TAUXX*UM + TAUXY*VM + TAUZX*WM-QX
      S5    = TAUXY*UM + TAUYY*VM + TAUYZ*WM-QY
      T5    = TAUZX*UM + TAUYZ*VM + TAUZZ*WM-QZ
      SH(I,J,K,1) = 0.0
      SH(I,J,K,2) = (ETXM*TAUXX + ETYM*TAUXY + ETZM*TAUZX)/AJAM
      SH(I,J,K,3) = (ETXM*TAUXY + ETYM*TAUYY + ETZM*TAUYZ)/AJAM
      SH(I,J,K,4) = (ETXM*TAUZX + ETYM*TAUYZ + ETZM*TAUZZ)/AJAM
      SH(I,J,K,5) = (ETXM*R5    + ETYM*S5    + ETZM*T5   )/AJAM
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISET
!***********************************************************************
!**** 拡散項(粘性項)の計算(ZE方向)                                  ****
!***********************************************************************
SUBROUTINE VISZE
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL :: UM, VM, WM, TM
  REAL :: UXIM, VXIM, WXIM, TXIM, &
  &       UETM, VETM, WETM, TETM, &
  &       UZEM, VZEM, WZEM, TZEM
  REAL :: UXM, VXM, WXM, TXM, &
  &       UYM, VYM, WYM, TYM, &
  &       UZM, VZM, WZM, TZM
  REAL :: AMUM
  REAL :: DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX
  REAL :: QX, QY, QZ
  REAL :: R5, S5, T5
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  UM, VM, WM, TM, &
  !$  UXIM, VXIM, WXIM, TXIM, &
  !$  UETM, VETM, WETM, TETM, &
  !$  UZEM, VZEM, WZEM, TZEM, &
  !$  UXM, VXM, WXM, TXM, &
  !$  UYM, VYM, WYM, TYM, &
  !$  UZM, VZM, WZM, TZM, &
  !$  AMUM, DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, &
  !$  QX, QY, QZ, R5, S5, T5 &
  !$)
  DO K = KS    , KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    IF( QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K+1,1) .GT. 0.0 ) THEN
      ! ZE方向の平均量の計算 -------------------------------------------
      XIXM  = 0.5 * ( XIX(I,J,K) + XIX(I,J,K+1) )
      XIYM  = 0.5 * ( XIY(I,J,K) + XIY(I,J,K+1) )
      XIZM  = 0.5 * ( XIZ(I,J,K) + XIZ(I,J,K+1) )
      ETXM  = 0.5 * ( ETX(I,J,K) + ETX(I,J,K+1) )
      ETYM  = 0.5 * ( ETY(I,J,K) + ETY(I,J,K+1) )
      ETZM  = 0.5 * ( ETZ(I,J,K) + ETZ(I,J,K+1) )
      ZEXM  = 0.5 * ( ZEX(I,J,K) + ZEX(I,J,K+1) )
      ZEYM  = 0.5 * ( ZEY(I,J,K) + ZEY(I,J,K+1) )
      ZEZM  = 0.5 * ( ZEZ(I,J,K) + ZEZ(I,J,K+1) )
      AJAM  = 0.5 * ( AJA(I,J,K) + AJA(I,J,K+1) )
      UM    = 0.5 * (   U(I,J,K) +   U(I,J,K+1) )
      VM    = 0.5 * (   V(I,J,K) +   V(I,J,K+1) )
      WM    = 0.5 * (   W(I,J,K) +   W(I,J,K+1) )
      TM    = 0.5 * (   T(I,J,K) +   T(I,J,K+1) )
      ! x,y,z微分の計算(ZE方向の平均量) --------------------------------
      UXIM   = 0.5 * ( - 0.5*(   U(I-1,J,K) +   U(I-1,J,K+1) ) &
      &                + 0.5*(   U(I+1,J,K) +   U(I+1,J,K+1) )  )
      VXIM   = 0.5 * ( - 0.5*(   V(I-1,J,K) +   V(I-1,J,K+1) ) &
      &                + 0.5*(   V(I+1,J,K) +   V(I+1,J,K+1) )  )
      WXIM   = 0.5 * ( - 0.5*(   W(I-1,J,K) +   W(I-1,J,K+1) ) &
      &                + 0.5*(   W(I+1,J,K) +   W(I+1,J,K+1) )  )
      TXIM   = 0.5 * ( - 0.5*(   T(I-1,J,K) +   T(I-1,J,K+1) ) &
      &                + 0.5*(   T(I+1,J,K) +   T(I+1,J,K+1) )  )
      UETM   = 0.5 * ( - 0.5*(   U(I,J-1,K) +   U(I,J-1,K+1) ) &
      &                + 0.5*(   U(I,J+1,K) +   U(I,J+1,K+1) )  )
      VETM   = 0.5 * ( - 0.5*(   V(I,J-1,K) +   V(I,J-1,K+1) ) &
      &                + 0.5*(   V(I,J+1,K) +   V(I,J+1,K+1) )  )
      WETM   = 0.5 * ( - 0.5*(   W(I,J-1,K) +   W(I,J-1,K+1) ) &
      &                + 0.5*(   W(I,J+1,K) +   W(I,J+1,K+1) )  )
      TETM   = 0.5 * ( - 0.5*(   T(I,J-1,K) +   T(I,J-1,K+1) ) &
      &                + 0.5*(   T(I,J+1,K) +   T(I,J+1,K+1) )  )
      UZEM   = ( -  U(I,J,K) +  U(I,J,K+1) )
      VZEM   = ( -  V(I,J,K) +  V(I,J,K+1) )
      WZEM   = ( -  W(I,J,K) +  W(I,J,K+1) )
      TZEM   = ( -  T(I,J,K) +  T(I,J,K+1) )
      UXM   = XIXM *  UXIM + ETXM *  UETM + ZEXM *  UZEM
      UYM   = XIYM *  UXIM + ETYM *  UETM + ZEYM *  UZEM
      UZM   = XIZM *  UXIM + ETZM *  UETM + ZEZM *  UZEM
      VXM   = XIXM *  VXIM + ETXM *  VETM + ZEXM *  VZEM
      VYM   = XIYM *  VXIM + ETYM *  VETM + ZEYM *  VZEM
      VZM   = XIZM *  VXIM + ETZM *  VETM + ZEZM *  VZEM
      WXM   = XIXM *  WXIM + ETXM *  WETM + ZEXM *  WZEM
      WYM   = XIYM *  WXIM + ETYM *  WETM + ZEYM *  WZEM
      WZM   = XIZM *  WXIM + ETZM *  WETM + ZEZM *  WZEM
      TXM   = XIXM *  TXIM + ETXM *  TETM + ZEXM *  TZEM
      TYM   = XIYM *  TXIM + ETYM *  TETM + ZEYM *  TZEM
      TZM   = XIZM *  TXIM + ETZM *  TETM + ZEZM *  TZEM
      ! 粘性係数のZE方向の平均量 ---------------------------------------
      AMUM  = 0.5 * (AMU(I,J,K) + AMU(I,J,K+1))
      ! 拡散項ベクトルの計算 -------------------------------------------
      DELVM = UXM+VYM+WZM
      TAUXX = AMUM*(2.0*UXM-2.0/3.0*DELVM)
      TAUYY = AMUM*(2.0*VYM-2.0/3.0*DELVM)
      TAUZZ = AMUM*(2.0*WZM-2.0/3.0*DELVM)
      TAUXY = AMUM*(UYM+VXM)
      TAUYZ = AMUM*(VZM+WYM)
      TAUZX = AMUM*(WXM+UZM)
      QX    =-(AMUM/PR)/(GAMMA-1.0)*GAMMA*RG*TXM
      QY    =-(AMUM/PR)/(GAMMA-1.0)*GAMMA*RG*TYM
      QZ    =-(AMUM/PR)/(GAMMA-1.0)*GAMMA*RG*TZM
      R5    = TAUXX*UM + TAUXY*VM + TAUZX*WM-QX
      S5    = TAUXY*UM + TAUYY*VM + TAUYZ*WM-QY
      T5    = TAUZX*UM + TAUYZ*VM + TAUZZ*WM-QZ
      TH(I,J,K,1) = 0.0
      TH(I,J,K,2) = (ZEXM*TAUXX + ZEYM*TAUXY + ZEZM*TAUZX)/AJAM
      TH(I,J,K,3) = (ZEXM*TAUXY + ZEYM*TAUYY + ZEZM*TAUYZ)/AJAM
      TH(I,J,K,4) = (ZEXM*TAUZX + ZEYM*TAUYZ + ZEZM*TAUZZ)/AJAM
      TH(I,J,K,5) = (ZEXM*R5    + ZEYM*S5    + ZEZM*T5   )/AJAM
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISZE
! 定義終了 *************************************************************
END SUBROUTINE Viscosity3D
!***********************************************************************
!**** 乱流モデル : RANS, EVM, Standard k-e Model (1974)             ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence2DEvmStd( &
&            PELIM, RG, GAMMA, PR, PRT, &
&            IS, IE, JS, JE, LS, LE, &
&            XIX, XIY, ETX, ETY, AJA, &
&            QH, U, V, T, AK, EPS, AMU, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: CMU = 0.09
  REAL, PARAMETER :: CE1 = 1.44, CE2 = 1.92
  REAL, PARAMETER :: SIGK = 1.0, SIGE = 1.3
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR, PRT
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS: IE, JS: JE), XIY(IS: IE, JS: JE), &
  &                       ETX(IS: IE, JS: JE), ETY(IS: IE, JS: JE), &
  &                       AJA(IS: IE, JS: JE)
  REAL,    INTENT(IN)  :: QH(IS: IE, JS: JE, LS: LE)
  REAL,    INTENT(IN)  :: U(IS: IE, JS: JE), V(IS: IE, JS: JE), &
  &                       T(IS: IE, JS: JE), &
  &                       AK(IS: IE, JS: JE), EPS(IS: IE, JS: JE)
  REAL,    INTENT(IN)  :: AMU(IS: IE, JS: JE)
  REAL,    INTENT(OUT) :: AMUT(IS: IE, JS: JE)
  REAL,    INTENT(OUT) :: DQD(IS: IE, JS: JE, LS: LE)
  REAL,    INTENT(OUT) :: DQP(IS: IE, JS: JE, LS: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :), SH(:, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 6) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, LS: LE), SH(IS: IE, JS: JE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: UXI, VXI, UET, VET
  REAL    :: UX, VX, UY, VY
  REAL    :: DELV, SS2, PRO, EPSPK
  ! 処理開始 ***********************************************************
  ! 渦粘性係数の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
  DO J = JS, JE
  DO I = IS, IE
  IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J,6) .GT. ZERO) THEN
    AMUT(I,J) = CMU * QH(I,J,5)**2 / QH(I,J,6) * AJA(I,J)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 生成項と散逸項の計算 +++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, UXI, VXI, UET, VET, UX, VX, UY, VY, DELV, SS2, PRO, EPSPK &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (QH(I,J-1,1) .GT. 0.0) .AND. (QH(I-1,J,1) .GT. 0.0) .AND. &
  &   (QH(I,J,1)   .GT. 0.0) .AND. &
  &   (QH(I+1,J,1) .GT. 0.0) .AND. (QH(I,J+1,1) .GT. 0.0) &
  & ) THEN
    ! ひずみ速度テンソルの計算 -----------------------------------------
    UXI = (-U(I-1,J  ) + U(I+1,J  )) * 0.5
    VXI = (-V(I-1,J  ) + V(I+1,J  )) * 0.5
    UET = (-U(I  ,J-1) + U(I  ,J+1)) * 0.5
    VET = (-V(I  ,J-1) + V(I  ,J+1)) * 0.5
    UX  = XIX(I,J) * UXI + ETX(I,J) * UET
    UY  = XIY(I,J) * UXI + ETY(I,J) * UET
    VX  = XIX(I,J) * VXI + ETX(I,J) * VET
    VY  = XIY(I,J) * VXI + ETY(I,J) * VET
    ! 生成項の計算 -----------------------------------------------------
    DELV = UX + VY
    SS2  = 2.0 * (UX**2 + VY**2) + (UY + VX)**2 - 2.0 / 3.0 * DELV**2
    PRO  = AMUT(I,J) * SS2 &
    &    - 2.0 / 3.0 * QH(I,J,5) * AJA(I,J) * DELV
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      PRO = MAX(0.0, MIN(QH(I,J,6) * AJA(I,J) / PELIM, PRO))
    ENDIF
    ! 生産と散逸の和 ---------------------------------------------------
    EPSPK = MIN( &
    &       MAX( &
    &         SQRT(SS2), SQRT(QH(I,J,6) * AJA(I,J) / AMU(I,J)) &
    &       ), &
    &       QH(I,J,6) / MAX(ZERO, QH(I,J,5)) &
    &     )
    DQP(I,J,5) = PRO / AJA(I,J) - QH(I,J,6)
    DQP(I,J,6) = (CE1 * PRO / AJA(I,J) - CE2 * QH(I,J,6)) * EPSPK
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  ! 処理開始 ***********************************************************
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL VISXI
  CALL VISET
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, L)
  DO L = LS, LE
  !$OMP DO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (QH(I,J-1,1) .GT. 0.0) .AND. (QH(I-1,J,1) .GT. 0.0) .AND. &
  &   (QH(I,J,1)   .GT. 0.0) .AND. &
  &   (QH(I+1,J,1) .GT. 0.0) .AND. (QH(I,J+1,1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,L) = (-RH(I-1,J  ,L) + RH(I,J,L)) &
    &          + (-SH(I  ,J-1,L) + SH(I,J,L))
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項(粘性項)の計算(XI方向)                                  ****
!***********************************************************************
SUBROUTINE VISXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: UM, VM, TM, RHOKM
  REAL    :: UXIM, VXIM, TXIM, AKXIM, EPSXIM, &
  &          UETM, VETM, TETM, AKETM, EPSETM
  REAL    :: UXM, VXM, TXM, AKXM, EPSXM, &
  &          UYM, VYM, TYM, AKYM, EPSYM
  REAL    :: AMUM, AMUTM, AMUKM, AMUEM
  REAL    :: DELVM, TAUXX, TAUYY, TAUXY
  REAL    :: QX, QY
  REAL    :: R4, S4
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, XIXM, XIYM, ETXM, ETYM, AJAM, UM, VM, TM, RHOKM, &
  !$  UXIM, VXIM, TXIM, AKXIM, EPSXIM, &
  !$  UETM, VETM, TETM, AKETM, EPSETM, &
  !$  UXM, VXM, TXM, AKXM, EPSXM, &
  !$  UYM, VYM, TYM, AKYM, EPSYM, &
  !$  AMUM, AMUTM, AMUKM, AMUEM, &
  !$  DELVM, TAUXX, TAUYY, TAUXY, QX, QY, R4, S4 &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF(QH(I,J,1) .GT. 0.0 .AND. QH(I+1,J,1) .GT. 0.0) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    XIXM  = 0.5 * (XIX(I,J) + XIX(I+1,J))
    XIYM  = 0.5 * (XIY(I,J) + XIY(I+1,J))
    ETXM  = 0.5 * (ETX(I,J) + ETX(I+1,J))
    ETYM  = 0.5 * (ETY(I,J) + ETY(I+1,J))
    AJAM  = 0.5 * (AJA(I,J) + AJA(I+1,J))
    UM    = 0.5 * (  U(I,J) +   U(I+1,J))
    VM    = 0.5 * (  V(I,J) +   V(I+1,J))
    TM    = 0.5 * (  T(I,J) +   T(I+1,J))
    RHOKM = 0.5 * (QH(I,J,5) * AJA(I,J) + QH(I+1,J,5) * AJA(I+1,J))
    ! x,y微分の計算(XI方向の平均量) ------------------------------------
    UXIM   = (-  U(I,J) +   U(I+1,J))
    VXIM   = (-  V(I,J) +   V(I+1,J))
    TXIM   = (-  T(I,J) +   T(I+1,J))
    AKXIM  = (- AK(I,J) +  AK(I+1,J))
    EPSXIM = (-EPS(I,J) + EPS(I+1,J))
    UETM   = 0.5 * ( - 0.5 * (  U(I,J-1) +   U(I+1,J-1)) &
    &                + 0.5 * (  U(I,J+1) +   U(I+1,J+1)) )
    VETM   = 0.5 * ( - 0.5 * (  V(I,J-1) +   V(I+1,J-1)) &
    &                + 0.5 * (  V(I,J+1) +   V(I+1,J+1)) )
    TETM   = 0.5 * ( - 0.5 * (  T(I,J-1) +   T(I+1,J-1)) &
    &                + 0.5 * (  T(I,J+1) +   T(I+1,J+1)) )
    AKETM  = 0.5 * ( - 0.5 * ( AK(I,J-1) +  AK(I+1,J-1)) &
    &                + 0.5 * ( AK(I,J+1) +  AK(I+1,J+1)) )
    EPSETM = 0.5 * ( - 0.5 * (EPS(I,J-1) + EPS(I+1,J-1)) &
    &                + 0.5 * (EPS(I,J+1) + EPS(I+1,J+1)) )
    UXM   = XIXM *   UXIM + ETXM *   UETM
    UYM   = XIYM *   UXIM + ETYM *   UETM
    VXM   = XIXM *   VXIM + ETXM *   VETM
    VYM   = XIYM *   VXIM + ETYM *   VETM
    TXM   = XIXM *   TXIM + ETXM *   TETM
    TYM   = XIYM *   TXIM + ETYM *   TETM
    AKXM  = XIXM *  AKXIM + ETXM *  AKETM
    AKYM  = XIYM *  AKXIM + ETYM *  AKETM
    EPSXM = XIXM * EPSXIM + ETXM * EPSETM
    EPSYM = XIYM * EPSXIM + ETYM * EPSETM
    ! 粘性係数のXI方向の平均量 -----------------------------------------
    AMUM  = 0.5 * ( AMU(I,J) +  AMU(I+1,J))
    AMUTM = 0.5 * (AMUT(I,J) + AMUT(I+1,J))
    AMUKM = AMUM + AMUTM / SIGK
    AMUEM = AMUM + AMUTM / SIGE
    ! 拡散項ベクトルの計算 ---------------------------------------------
    DELVM = UXM + VYM
    TAUXX = (AMUM + AMUTM) * (2.0 * UXM - 2.0 / 3.0 * DELVM) &
    &     - 2.0 / 3.0 * RHOKM
    TAUYY = (AMUM + AMUTM) * (2.0 * VYM - 2.0 / 3.0 * DELVM) &
    &     - 2.0 / 3.0 * RHOKM
    TAUXY = (AMUM + AMUTM) * (UYM + VXM)
    QX    =-(AMUM / PR + AMUTM / PRT) &
    &     / (GAMMA - 1.0) * GAMMA * RG * TXM
    QY    =-(AMUM / PR + AMUTM / PRT) &
    &     / (GAMMA - 1.0) * GAMMA * RG * TYM
    R4    = TAUXX * UM + TAUXY * VM - QX
    S4    = TAUXY * UM + TAUYY * VM - QY
    RH(I,J,1) = 0.0
    RH(I,J,2) = (XIXM * TAUXX + XIYM * TAUXY) / AJAM
    RH(I,J,3) = (XIXM * TAUXY + XIYM * TAUYY) / AJAM
    RH(I,J,4) = (XIXM * R4    + XIYM * S4   ) / AJAM
    RH(I,J,5) = AMUKM * (XIXM *  AKXM + XIYM *  AKYM) / AJAM
    RH(I,J,6) = AMUEM * (XIXM * EPSXM + XIYM * EPSYM) / AJAM
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISXI
!***********************************************************************
!**** 拡散項(粘性項)の計算(ET方向)                                  ****
!***********************************************************************
SUBROUTINE VISET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: UM, VM, TM, RHOKM
  REAL    :: UXIM, VXIM, TXIM, AKXIM, EPSXIM, &
  &          UETM, VETM, TETM, AKETM, EPSETM
  REAL    :: UXM, VXM, TXM, AKXM, EPSXM, &
  &          UYM, VYM, TYM, AKYM, EPSYM
  REAL    :: AMUM, AMUTM, AMUKM, AMUEM
  REAL    :: DELVM, TAUXX, TAUYY, TAUXY
  REAL    :: QX, QY
  REAL    :: R4, S4
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, XIXM, XIYM, ETXM, ETYM, AJAM, UM, VM, TM, RHOKM, &
  !$  UXIM, VXIM, TXIM, AKXIM, EPSXIM, &
  !$  UETM, VETM, TETM, AKETM, EPSETM, &
  !$  UXM, VXM, TXM, AKXM, EPSXM, &
  !$  UYM, VYM, TYM, AKYM, EPSYM, &
  !$  AMUM, AMUTM, AMUKM, AMUEM, &
  !$  DELVM, TAUXX, TAUYY, TAUXY, QX, QY, R4, S4 &
  !$)
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J+1,1) .GT. 0.0) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    XIXM  = 0.5 * (XIX(I,J) + XIX(I,J+1))
    XIYM  = 0.5 * (XIY(I,J) + XIY(I,J+1))
    ETXM  = 0.5 * (ETX(I,J) + ETX(I,J+1))
    ETYM  = 0.5 * (ETY(I,J) + ETY(I,J+1))
    AJAM  = 0.5 * (AJA(I,J) + AJA(I,J+1))
    UM    = 0.5 * (  U(I,J) +   U(I,J+1))
    VM    = 0.5 * (  V(I,J) +   V(I,J+1))
    TM    = 0.5 * (  T(I,J) +   T(I,J+1))
    RHOKM = 0.5 * (QH(I,J,5) * AJA(I,J) + QH(I,J+1,5) * AJA(I,J+1))
    ! x,y微分の計算(ET方向の平均量) ------------------------------------
    UXIM   = 0.5 * ( - 0.5*(  U(I-1,J) +   U(I-1,J+1)) &
    &                + 0.5*(  U(I+1,J) +   U(I+1,J+1)) )
    VXIM   = 0.5 * ( - 0.5*(  V(I-1,J) +   V(I-1,J+1)) &
    &                + 0.5*(  V(I+1,J) +   V(I+1,J+1)) )
    TXIM   = 0.5 * ( - 0.5*(  T(I-1,J) +   T(I-1,J+1)) &
    &                + 0.5*(  T(I+1,J) +   T(I+1,J+1)) )
    AKXIM  = 0.5 * ( - 0.5*( AK(I-1,J) +  AK(I-1,J+1)) &
    &                + 0.5*( AK(I+1,J) +  AK(I+1,J+1)) )
    EPSXIM = 0.5 * ( - 0.5*(EPS(I-1,J) + EPS(I-1,J+1)) &
    &                + 0.5*(EPS(I+1,J) + EPS(I+1,J+1)) )
    UETM   = (-  U(I,J) +   U(I,J+1))
    VETM   = (-  V(I,J) +   V(I,J+1))
    TETM   = (-  T(I,J) +   T(I,J+1))
    AKETM  = (- AK(I,J) +  AK(I,J+1))
    EPSETM = (-EPS(I,J) + EPS(I,J+1))
    UXM   = XIXM *   UXIM + ETXM *   UETM
    UYM   = XIYM *   UXIM + ETYM *   UETM
    VXM   = XIXM *   VXIM + ETXM *   VETM
    VYM   = XIYM *   VXIM + ETYM *   VETM
    TXM   = XIXM *   TXIM + ETXM *   TETM
    TYM   = XIYM *   TXIM + ETYM *   TETM
    AKXM  = XIXM *  AKXIM + ETXM *  AKETM
    AKYM  = XIYM *  AKXIM + ETYM *  AKETM
    EPSXM = XIXM * EPSXIM + ETXM * EPSETM
    EPSYM = XIYM * EPSXIM + ETYM * EPSETM
    ! 粘性係数のET方向の平均量 -----------------------------------------
    AMUM  = 0.5 * ( AMU(I,J) +  AMU(I,J+1))
    AMUTM = 0.5 * (AMUT(I,J) + AMUT(I,J+1))
    AMUKM = AMUM + AMUTM / SIGK
    AMUEM = AMUM + AMUTM / SIGE
    ! 拡散項ベクトルの計算 ---------------------------------------------
    DELVM = UXM + VYM
    TAUXX = (AMUM + AMUTM) * (2.0 * UXM - 2.0 / 3.0 * DELVM) &
    &     - 2.0 / 3.0 * RHOKM
    TAUYY = (AMUM + AMUTM) * (2.0 * VYM - 2.0 / 3.0 * DELVM) &
    &     - 2.0 / 3.0 * RHOKM
    TAUXY = (AMUM + AMUTM) * (UYM + VXM)
    QX    =-(AMUM / PR + AMUTM / PRT) &
    &     / (GAMMA - 1.0) * GAMMA * RG * TXM
    QY    =-(AMUM / PR + AMUTM / PRT) &
    &     / (GAMMA - 1.0) * GAMMA * RG * TYM
    R4    = TAUXX * UM + TAUXY * VM - QX
    S4    = TAUXY * UM + TAUYY * VM - QY
    SH(I,J,1) = 0.0
    SH(I,J,2) = (ETXM * TAUXX + ETYM * TAUXY) / AJAM
    SH(I,J,3) = (ETXM * TAUXY + ETYM * TAUYY) / AJAM
    SH(I,J,4) = (ETXM * R4    + ETYM * S4   ) / AJAM
    SH(I,J,5) = AMUKM * (ETXM *  AKXM + ETYM *  AKYM) / AJAM
    SH(I,J,6) = AMUEM * (ETXM * EPSXM + ETYM * EPSYM) / AJAM
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISET
! 定義終了 *************************************************************
END SUBROUTINE Turbulence2DEvmStd
!***********************************************************************
!**** 乱流モデル : RANS, EVM, Standard k-e Model (1974)             ****
!****                         with Kato-Launder Modification (1993) ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence2DEvmStdKL( &
&            PELIM, RG, GAMMA, PR, PRT, &
&            IS, IE, JS, JE, LS, LE, &
&            XIX, XIY, ETX, ETY, AJA, &
&            QH, U, V, T, AK, EPS, AMU, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: CMU = 0.09
  REAL, PARAMETER :: CE1 = 1.44, CE2 = 1.92
  REAL, PARAMETER :: SIGK = 1.0, SIGE = 1.3
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR, PRT
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS: IE, JS: JE), XIY(IS: IE, JS: JE), &
  &                       ETX(IS: IE, JS: JE), ETY(IS: IE, JS: JE), &
  &                       AJA(IS: IE, JS: JE)
  REAL,    INTENT(IN)  :: QH(IS: IE, JS: JE, LS: LE)
  REAL,    INTENT(IN)  :: U(IS: IE, JS: JE), V(IS: IE, JS: JE), &
  &                       T(IS: IE, JS: JE), &
  &                       AK(IS: IE, JS: JE), EPS(IS: IE, JS: JE)
  REAL,    INTENT(IN)  :: AMU(IS: IE, JS: JE)
  REAL,    INTENT(OUT) :: AMUT(IS: IE, JS: JE)
  REAL,    INTENT(OUT) :: DQD(IS: IE, JS: JE, LS: LE)
  REAL,    INTENT(OUT) :: DQP(IS: IE, JS: JE, LS: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :), SH(:, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 6) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, LS: LE), SH(IS: IE, JS: JE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: UXI, VXI, UET, VET
  REAL    :: UX, VX, UY, VY
  REAL    :: DELV, SS2, SC, OMG, OS, PRO, EPSPK
  ! 処理開始 ***********************************************************
  ! 渦粘性係数の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
  DO J = JS, JE
  DO I = IS, IE
  IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J,6) .GT. ZERO) THEN
    AMUT(I,J) = CMU * QH(I,J,5)**2 / QH(I,J,6) * AJA(I,J)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 生成項と散逸項の計算 +++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  UXI, VXI, UET, VET, UX, UY, VX, VY, &
  !$  DELV, SS2, SC, OMG, OS, PRO, EPSPK &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (QH(I,J-1,1) .GT. 0.0) .AND. (QH(I-1,J,1) .GT. 0.0) .AND. &
  &   (QH(I,J,1)   .GT. 0.0) .AND. &
  &   (QH(I+1,J,1) .GT. 0.0) .AND. (QH(I,J+1,1) .GT. 0.0) &
  & ) THEN
    ! ひずみ速度テンソルの計算 -----------------------------------------
    UXI = (-U(I-1,J  ) + U(I+1,J  )) * 0.5
    VXI = (-V(I-1,J  ) + V(I+1,J  )) * 0.5
    UET = (-U(I  ,J-1) + U(I  ,J+1)) * 0.5
    VET = (-V(I  ,J-1) + V(I  ,J+1)) * 0.5
    UX  = XIX(I,J) * UXI + ETX(I,J) * UET
    UY  = XIY(I,J) * UXI + ETY(I,J) * UET
    VX  = XIX(I,J) * VXI + ETX(I,J) * VET
    VY  = XIY(I,J) * VXI + ETY(I,J) * VET
    ! 生成項の計算 -----------------------------------------------------
    DELV = UX + VY
    SS2  = 2.0 * (UX**2 + VY**2) + (UY + VX)**2 - 2.0 / 3.0 * DELV**2
    SC  = SQRT(SS2)
    OMG = SQRT((UY - VX)**2)
    ! OMG / SC > 1 の場合、Kato-Launderの修正を適用しない
    IF(SC .GT. ZERO) THEN
      OS = MIN(OMG / SC, 1.0)
    ELSE
      OS = 1.0
    ENDIF
    PRO = OS * ( &
    &     AMUT(I,J) * SS2 &
    &   - 2.0 / 3.0 * DELV * QH(I,J,5) * AJA(I,J) &
    &   )
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      PRO = MAX(0.0, MIN(QH(I,J,6) * AJA(I,J) / PELIM, PRO))
    ENDIF
    ! 生産と散逸の和 ---------------------------------------------------
    EPSPK = MIN( &
    &       MAX(SC, SQRT(QH(I,J,6) * AJA(I,J) / AMU(I,J))), &
    &       QH(I,J,6) / MAX(ZERO, QH(I,J,5)) &
    &     )
    DQP(I,J,5) = PRO / AJA(I,J) - QH(I,J,6)
    DQP(I,J,6) = (CE1 * PRO / AJA(I,J) - CE2 * QH(I,J,6)) * EPSPK
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  ! 処理開始 ***********************************************************
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL VISXI
  CALL VISET
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, L)
  DO L = LS, LE
  !$OMP DO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (QH(I,J-1,1) .GT. 0.0) .AND. (QH(I-1,J,1) .GT. 0.0) .AND. &
  &   (QH(I,J,1)   .GT. 0.0) .AND. &
  &   (QH(I+1,J,1) .GT. 0.0) .AND. (QH(I,J+1,1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,L) = (-RH(I-1,J  ,L) + RH(I,J,L)) &
    &          + (-SH(I  ,J-1,L) + SH(I,J,L))
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項(粘性項)の計算(XI方向)                                  ****
!***********************************************************************
SUBROUTINE VISXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: UM, VM, TM, RHOKM
  REAL    :: UXIM, VXIM, TXIM, AKXIM, EPSXIM, &
  &          UETM, VETM, TETM, AKETM, EPSETM
  REAL    :: UXM, VXM, TXM, AKXM, EPSXM, &
  &          UYM, VYM, TYM, AKYM, EPSYM
  REAL    :: AMUM, AMUTM, AMUKM, AMUEM
  REAL    :: DELVM, TAUXX, TAUYY, TAUXY
  REAL    :: QX, QY
  REAL    :: R4, S4
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, XIXM, XIYM, ETXM, ETYM, AJAM, UM, VM, TM, RHOKM, &
  !$  UXIM, VXIM, TXIM, AKXIM, EPSXIM, &
  !$  UETM, VETM, TETM, AKETM, EPSETM, &
  !$  UXM, VXM, TXM, AKXM, EPSXM, &
  !$  UYM, VYM, TYM, AKYM, EPSYM, &
  !$  AMUM, AMUTM, AMUKM, AMUEM, &
  !$  DELVM, TAUXX, TAUYY, TAUXY, QX, QY, R4, S4 &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF(QH(I,J,1) .GT. 0.0 .AND. QH(I+1,J,1) .GT. 0.0) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    XIXM  = 0.5 * (XIX(I,J) + XIX(I+1,J))
    XIYM  = 0.5 * (XIY(I,J) + XIY(I+1,J))
    ETXM  = 0.5 * (ETX(I,J) + ETX(I+1,J))
    ETYM  = 0.5 * (ETY(I,J) + ETY(I+1,J))
    AJAM  = 0.5 * (AJA(I,J) + AJA(I+1,J))
    UM    = 0.5 * (  U(I,J) +   U(I+1,J))
    VM    = 0.5 * (  V(I,J) +   V(I+1,J))
    TM    = 0.5 * (  T(I,J) +   T(I+1,J))
    RHOKM = 0.5 * (QH(I,J,5) * AJA(I,J) + QH(I+1,J,5) * AJA(I+1,J))
    ! x,y微分の計算(XI方向の平均量) ------------------------------------
    UXIM   = (-  U(I,J) +   U(I+1,J))
    VXIM   = (-  V(I,J) +   V(I+1,J))
    TXIM   = (-  T(I,J) +   T(I+1,J))
    AKXIM  = (- AK(I,J) +  AK(I+1,J))
    EPSXIM = (-EPS(I,J) + EPS(I+1,J))
    UETM   = 0.5 * ( - 0.5 * (  U(I,J-1) +   U(I+1,J-1)) &
    &                + 0.5 * (  U(I,J+1) +   U(I+1,J+1)) )
    VETM   = 0.5 * ( - 0.5 * (  V(I,J-1) +   V(I+1,J-1)) &
    &                + 0.5 * (  V(I,J+1) +   V(I+1,J+1)) )
    TETM   = 0.5 * ( - 0.5 * (  T(I,J-1) +   T(I+1,J-1)) &
    &                + 0.5 * (  T(I,J+1) +   T(I+1,J+1)) )
    AKETM  = 0.5 * ( - 0.5 * ( AK(I,J-1) +  AK(I+1,J-1)) &
    &                + 0.5 * ( AK(I,J+1) +  AK(I+1,J+1)) )
    EPSETM = 0.5 * ( - 0.5 * (EPS(I,J-1) + EPS(I+1,J-1)) &
    &                + 0.5 * (EPS(I,J+1) + EPS(I+1,J+1)) )
    UXM   = XIXM *   UXIM + ETXM *   UETM
    UYM   = XIYM *   UXIM + ETYM *   UETM
    VXM   = XIXM *   VXIM + ETXM *   VETM
    VYM   = XIYM *   VXIM + ETYM *   VETM
    TXM   = XIXM *   TXIM + ETXM *   TETM
    TYM   = XIYM *   TXIM + ETYM *   TETM
    AKXM  = XIXM *  AKXIM + ETXM *  AKETM
    AKYM  = XIYM *  AKXIM + ETYM *  AKETM
    EPSXM = XIXM * EPSXIM + ETXM * EPSETM
    EPSYM = XIYM * EPSXIM + ETYM * EPSETM
    ! 粘性係数のXI方向の平均量 -----------------------------------------
    AMUM  = 0.5 * ( AMU(I,J) +  AMU(I+1,J))
    AMUTM = 0.5 * (AMUT(I,J) + AMUT(I+1,J))
    AMUKM = AMUM + AMUTM / SIGK
    AMUEM = AMUM + AMUTM / SIGE
    ! 拡散項ベクトルの計算 ---------------------------------------------
    DELVM = UXM + VYM
    TAUXX = (AMUM + AMUTM) * (2.0 * UXM - 2.0 / 3.0 * DELVM) &
    &     - 2.0 / 3.0 * RHOKM
    TAUYY = (AMUM + AMUTM) * (2.0 * VYM - 2.0 / 3.0 * DELVM) &
    &     - 2.0 / 3.0 * RHOKM
    TAUXY = (AMUM + AMUTM) * (UYM + VXM)
    QX    =-(AMUM / PR + AMUTM / PRT) &
    &     / (GAMMA - 1.0) * GAMMA * RG * TXM
    QY    =-(AMUM / PR + AMUTM / PRT) &
    &     / (GAMMA - 1.0) * GAMMA * RG * TYM
    R4    = TAUXX * UM + TAUXY * VM - QX
    S4    = TAUXY * UM + TAUYY * VM - QY
    RH(I,J,1) = 0.0
    RH(I,J,2) = (XIXM * TAUXX + XIYM * TAUXY) / AJAM
    RH(I,J,3) = (XIXM * TAUXY + XIYM * TAUYY) / AJAM
    RH(I,J,4) = (XIXM * R4    + XIYM * S4   ) / AJAM
    RH(I,J,5) = AMUKM * (XIXM *  AKXM + XIYM *  AKYM) / AJAM
    RH(I,J,6) = AMUEM * (XIXM * EPSXM + XIYM * EPSYM) / AJAM
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISXI
!***********************************************************************
!**** 拡散項(粘性項)の計算(ET方向)                                  ****
!***********************************************************************
SUBROUTINE VISET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: UM, VM, TM, RHOKM
  REAL    :: UXIM, VXIM, TXIM, AKXIM, EPSXIM, &
  &          UETM, VETM, TETM, AKETM, EPSETM
  REAL    :: UXM, VXM, TXM, AKXM, EPSXM, &
  &          UYM, VYM, TYM, AKYM, EPSYM
  REAL    :: AMUM, AMUTM, AMUKM, AMUEM
  REAL    :: DELVM, TAUXX, TAUYY, TAUXY
  REAL    :: QX, QY
  REAL    :: R4, S4
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, XIXM, XIYM, ETXM, ETYM, AJAM, UM, VM, TM, RHOKM, &
  !$  UXIM, VXIM, TXIM, AKXIM, EPSXIM, &
  !$  UETM, VETM, TETM, AKETM, EPSETM, &
  !$  UXM, VXM, TXM, AKXM, EPSXM, &
  !$  UYM, VYM, TYM, AKYM, EPSYM, &
  !$  AMUM, AMUTM, AMUKM, AMUEM, &
  !$  DELVM, TAUXX, TAUYY, TAUXY, QX, QY, R4, S4 &
  !$)
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF(QH(I,J,1) .GT. 0.0 .AND. QH(I,J+1,1) .GT. 0.0) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    XIXM  = 0.5 * (XIX(I,J) + XIX(I,J+1))
    XIYM  = 0.5 * (XIY(I,J) + XIY(I,J+1))
    ETXM  = 0.5 * (ETX(I,J) + ETX(I,J+1))
    ETYM  = 0.5 * (ETY(I,J) + ETY(I,J+1))
    AJAM  = 0.5 * (AJA(I,J) + AJA(I,J+1))
    UM    = 0.5 * (  U(I,J) +   U(I,J+1))
    VM    = 0.5 * (  V(I,J) +   V(I,J+1))
    TM    = 0.5 * (  T(I,J) +   T(I,J+1))
    RHOKM = 0.5 * (QH(I,J,5) * AJA(I,J) + QH(I,J+1,5) * AJA(I,J+1))
    ! x,y微分の計算(ET方向の平均量) ------------------------------------
    UXIM   = 0.5 * ( - 0.5*(  U(I-1,J) +   U(I-1,J+1)) &
    &                + 0.5*(  U(I+1,J) +   U(I+1,J+1)) )
    VXIM   = 0.5 * ( - 0.5*(  V(I-1,J) +   V(I-1,J+1)) &
    &                + 0.5*(  V(I+1,J) +   V(I+1,J+1)) )
    TXIM   = 0.5 * ( - 0.5*(  T(I-1,J) +   T(I-1,J+1)) &
    &                + 0.5*(  T(I+1,J) +   T(I+1,J+1)) )
    AKXIM  = 0.5 * ( - 0.5*( AK(I-1,J) +  AK(I-1,J+1)) &
    &                + 0.5*( AK(I+1,J) +  AK(I+1,J+1)) )
    EPSXIM = 0.5 * ( - 0.5*(EPS(I-1,J) + EPS(I-1,J+1)) &
    &                + 0.5*(EPS(I+1,J) + EPS(I+1,J+1)) )
    UETM   = (-  U(I,J) +   U(I,J+1))
    VETM   = (-  V(I,J) +   V(I,J+1))
    TETM   = (-  T(I,J) +   T(I,J+1))
    AKETM  = (- AK(I,J) +  AK(I,J+1))
    EPSETM = (-EPS(I,J) + EPS(I,J+1))
    UXM   = XIXM *   UXIM + ETXM *   UETM
    UYM   = XIYM *   UXIM + ETYM *   UETM
    VXM   = XIXM *   VXIM + ETXM *   VETM
    VYM   = XIYM *   VXIM + ETYM *   VETM
    TXM   = XIXM *   TXIM + ETXM *   TETM
    TYM   = XIYM *   TXIM + ETYM *   TETM
    AKXM  = XIXM *  AKXIM + ETXM *  AKETM
    AKYM  = XIYM *  AKXIM + ETYM *  AKETM
    EPSXM = XIXM * EPSXIM + ETXM * EPSETM
    EPSYM = XIYM * EPSXIM + ETYM * EPSETM
    ! 粘性係数のET方向の平均量 -----------------------------------------
    AMUM  = 0.5 * ( AMU(I,J) +  AMU(I,J+1))
    AMUTM = 0.5 * (AMUT(I,J) + AMUT(I,J+1))
    AMUKM = AMUM + AMUTM / SIGK
    AMUEM = AMUM + AMUTM / SIGE
    ! 拡散項ベクトルの計算 ---------------------------------------------
    DELVM = UXM + VYM
    TAUXX = (AMUM + AMUTM) * (2.0 * UXM - 2.0 / 3.0 * DELVM) &
    &     - 2.0 / 3.0 * RHOKM
    TAUYY = (AMUM + AMUTM) * (2.0 * VYM - 2.0 / 3.0 * DELVM) &
    &     - 2.0 / 3.0 * RHOKM
    TAUXY = (AMUM + AMUTM) * (UYM + VXM)
    QX    =-(AMUM / PR + AMUTM / PRT) &
    &     / (GAMMA - 1.0) * GAMMA * RG * TXM
    QY    =-(AMUM / PR + AMUTM / PRT) &
    &     / (GAMMA - 1.0) * GAMMA * RG * TYM
    R4    = TAUXX * UM + TAUXY * VM - QX
    S4    = TAUXY * UM + TAUYY * VM - QY
    SH(I,J,1) = 0.0
    SH(I,J,2) = (ETXM * TAUXX + ETYM * TAUXY) / AJAM
    SH(I,J,3) = (ETXM * TAUXY + ETYM * TAUYY) / AJAM
    SH(I,J,4) = (ETXM * R4    + ETYM * S4   ) / AJAM
    SH(I,J,5) = AMUKM * (ETXM *  AKXM + ETYM *  AKYM) / AJAM
    SH(I,J,6) = AMUEM * (ETXM * EPSXM + ETYM * EPSYM) / AJAM
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISET
! 定義終了 *************************************************************
END SUBROUTINE Turbulence2DEvmStdKL
!***********************************************************************
!**** 乱流モデル : RANS, EVM, Standard k-e Model (1974)             ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence3DEvmStd( &
&            PELIM, RG, GAMMA, PR, PRT, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, AJA, &
&            QH, U, V, W, T, AK, EPS, AMU, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: CMU = 0.09
  REAL, PARAMETER :: CE1 = 1.44, CE2 = 1.92
  REAL, PARAMETER :: SIGK = 1.0, SIGE = 1.3
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR, PRT
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS: IE, JS: JE, KS: KE), &
  &                       XIY(IS: IE, JS: JE, KS: KE), &
  &                       XIZ(IS: IE, JS: JE, KS: KE), &
  &                       ETX(IS: IE, JS: JE, KS: KE), &
  &                       ETY(IS: IE, JS: JE, KS: KE), &
  &                       ETZ(IS: IE, JS: JE, KS: KE), &
  &                       ZEX(IS: IE, JS: JE, KS: KE), &
  &                       ZEY(IS: IE, JS: JE, KS: KE), &
  &                       ZEZ(IS: IE, JS: JE, KS: KE), &
  &                       AJA(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)  :: QH(IS: IE, JS: JE, KS: KE, LS: LE)
  REAL,    INTENT(IN)  :: U(IS: IE, JS: JE, KS: KE), &
  &                       V(IS: IE, JS: JE, KS: KE), &
  &                       W(IS: IE, JS: JE, KS: KE), &
  &                       T(IS: IE, JS: JE, KS: KE), &
  &                       AK(IS: IE, JS: JE, KS: KE), &
  &                       EPS(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)  :: AMU(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(OUT) :: AMUT(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(OUT) :: DQD(IS: IE, JS: JE, KS: KE, LS: LE)
  REAL,    INTENT(OUT) :: DQP(IS: IE, JS: JE, KS: KE, LS: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :, :), SH(:, :, :, :), TH(:, :, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 7) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(SH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(TH(IS: IE, JS: JE, KS: KE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  REAL    :: UX, VX, WX, UY, VY, WY, UZ, VZ, WZ
  REAL    :: DELV, SS2, PRO, EPSPK
  ! 処理開始 ***********************************************************
  ! 渦粘性係数の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, K)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K,7) .GT. ZERO) THEN
    AMUT(I,J,K) = CMU * QH(I,J,K,6)**2 / QH(I,J,K,7) * AJA(I,J,K)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 生成項と散逸項の計算 +++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE, &
  !$  UX, VX, WX, UY, VY, WY, UZ, VZ, WZ, &
  !$  DELV, PRO, EPSPK &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( ( QH(I  ,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I-1,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I+1,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J-1,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J+1,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J  ,K-1,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J  ,K+1,1) .GT. 0.0 ) &
  & ) THEN
    ! ひずみ速度テンソルの計算 -----------------------------------------
    UXI = ( -  U(I-1,J  ,K  ) +  U(I+1,J  ,K  ) )*0.5
    VXI = ( -  V(I-1,J  ,K  ) +  V(I+1,J  ,K  ) )*0.5
    WXI = ( -  W(I-1,J  ,K  ) +  W(I+1,J  ,K  ) )*0.5
    UET = ( -  U(I  ,J-1,K  ) +  U(I  ,J+1,K  ) )*0.5
    VET = ( -  V(I  ,J-1,K  ) +  V(I  ,J+1,K  ) )*0.5
    WET = ( -  W(I  ,J-1,K  ) +  W(I  ,J+1,K  ) )*0.5
    UZE = ( -  U(I  ,J  ,K-1) +  U(I  ,J  ,K+1) )*0.5
    VZE = ( -  V(I  ,J  ,K-1) +  V(I  ,J  ,K+1) )*0.5
    WZE = ( -  W(I  ,J  ,K-1) +  W(I  ,J  ,K+1) )*0.5
    UX  = XIX(I,J,K)*UXI + ETX(I,J,K)*UET + ZEX(I,J,K)*UZE
    UY  = XIY(I,J,K)*UXI + ETY(I,J,K)*UET + ZEY(I,J,K)*UZE
    UZ  = XIZ(I,J,K)*UXI + ETZ(I,J,K)*UET + ZEZ(I,J,K)*UZE
    VX  = XIX(I,J,K)*VXI + ETX(I,J,K)*VET + ZEX(I,J,K)*VZE
    VY  = XIY(I,J,K)*VXI + ETY(I,J,K)*VET + ZEY(I,J,K)*VZE
    VZ  = XIZ(I,J,K)*VXI + ETZ(I,J,K)*VET + ZEZ(I,J,K)*VZE
    WX  = XIX(I,J,K)*WXI + ETX(I,J,K)*WET + ZEX(I,J,K)*WZE
    WY  = XIY(I,J,K)*WXI + ETY(I,J,K)*WET + ZEY(I,J,K)*WZE
    WZ  = XIZ(I,J,K)*WXI + ETZ(I,J,K)*WET + ZEZ(I,J,K)*WZE
    ! 生成項の計算 -----------------------------------------------------
    DELV = UX+VY+WZ
    SS2  = 2.0*(UX**2 + VY**2 + WZ**2) &
    &    + (UY+VX)**2 + (VZ+WY)**2 + (WX+UZ)**2 &
    &    - 2.0/3.0 * DELV**2
    PRO  = AMUT(I,J,K)*SS2 &
    &    - 2.0/3.0 * QH(I,J,K,6)*AJA(I,J,K) * DELV
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      PRO = MAX(0.0, MIN(QH(I,J,K,7) * AJA(I,J,K) / PELIM, PRO))
    ENDIF
    ! 生成と散逸の和 ---------------------------------------------------
    EPSPK = MIN( &
    &       MAX( &
    &         SQRT(SS2), &
    &         SQRT(QH(I,J,K,7) * AJA(I,J,K) / AMU(I,J,K)) &
    &       ), &
    &       QH(I,J,K,7) / MAX(ZERO, QH(I,J,K,6)) &
    &     )
    DQP(I,J,K,6) = PRO/AJA(I,J,K) - QH(I,J,K,7)
    DQP(I,J,K,7) = (CE1*PRO/AJA(I,J,K)-CE2*QH(I,J,K,7))*EPSPK
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  ! 処理開始 ***********************************************************
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL VISXI
  CALL VISET
  CALL VISZE
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (QH(I,J,K-1,1) .GT. 0.0) .AND. (QH(I,J-1,K,1) .GT. 0.0) .AND. &
  &   (QH(I-1,J,K,1) .GT. 0.0) .AND. (QH(I,J,K,1)   .GT. 0.0) .AND. &
  &   (QH(I+1,J,K,1) .GT. 0.0) .AND. (QH(I,J+1,K,1) .GT. 0.0) .AND. &
  &   (QH(I,J,K+1,1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,K,L) = (-RH(I-1,J  ,K  ,L) + RH(I,J,K,L) ) &
    &            + (-SH(I  ,J-1,K  ,L) + SH(I,J,K,L) ) &
    &            + (-TH(I  ,J  ,K-1,L) + TH(I,J,K,L) )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項(粘性項)の計算(XI方向)                                  ****
!***********************************************************************
SUBROUTINE VISXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: UM, VM, WM, TM, RHOKM
  REAL    :: UXIM, VXIM, WXIM, TXIM, AKXIM, EPSXIM, &
  &          UETM, VETM, WETM, TETM, AKETM, EPSETM, &
  &          UZEM, VZEM, WZEM, TZEM, AKZEM, EPSZEM
  REAL    :: UXM, VXM, WXM, TXM, AKXM, EPSXM, &
  &          UYM, VYM, WYM, TYM, AKYM, EPSYM, &
  &          UZM, VZM, WZM, TZM, AKZM, EPSZM
  REAL    :: AMUM, AMUTM, AMUKM, AMUEM
  REAL    :: DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX
  REAL    :: QX, QY, QZ
  REAL    :: R5, S5, T5
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  UM, VM, WM, TM, RHOKM, &
  !$  UXIM, VXIM, WXIM, TXIM, AKXIM, EPSXIM, &
  !$  UETM, VETM, WETM, TETM, AKETM, EPSETM, &
  !$  UZEM, VZEM, WZEM, TZEM, AKZEM, EPSZEM, &
  !$  UXM, VXM, WXM, TXM, AKXM, EPSXM, &
  !$  UYM, VYM, WYM, TYM, AKYM, EPSYM, &
  !$  UZM, VZM, WZM, TZM, AKZM, EPSZM, &
  !$  AMUM, AMUTM, AMUKM, AMUEM, &
  !$  DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, &
  !$  QX, QY, QZ, R5, S5, T5 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF( QH(I,J,K,1) .GT. 0.0 .AND. QH(I+1,J,K,1) .GT. 0.0 ) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    XIXM  = 0.5 * ( XIX(I,J,K) + XIX(I+1,J,K) )
    XIYM  = 0.5 * ( XIY(I,J,K) + XIY(I+1,J,K) )
    XIZM  = 0.5 * ( XIZ(I,J,K) + XIZ(I+1,J,K) )
    ETXM  = 0.5 * ( ETX(I,J,K) + ETX(I+1,J,K) )
    ETYM  = 0.5 * ( ETY(I,J,K) + ETY(I+1,J,K) )
    ETZM  = 0.5 * ( ETZ(I,J,K) + ETZ(I+1,J,K) )
    ZEXM  = 0.5 * ( ZEX(I,J,K) + ZEX(I+1,J,K) )
    ZEYM  = 0.5 * ( ZEY(I,J,K) + ZEY(I+1,J,K) )
    ZEZM  = 0.5 * ( ZEZ(I,J,K) + ZEZ(I+1,J,K) )
    AJAM  = 0.5 * ( AJA(I,J,K) + AJA(I+1,J,K) )
    UM    = 0.5 * (   U(I,J,K) +   U(I+1,J,K) )
    VM    = 0.5 * (   V(I,J,K) +   V(I+1,J,K) )
    WM    = 0.5 * (   W(I,J,K) +   W(I+1,J,K) )
    TM    = 0.5 * (   T(I,J,K) +   T(I+1,J,K) )
    RHOKM = 0.5 * (  QH(I  ,J,K,6)*AJA(I  ,J,K) &
    &             +  QH(I+1,J,K,6)*AJA(I+1,J,K) )
    ! x,y,z微分の計算(XI方向の平均量) ----------------------------------
    UXIM   = ( -  U(I,J,K) +  U(I+1,J,K) )
    VXIM   = ( -  V(I,J,K) +  V(I+1,J,K) )
    WXIM   = ( -  W(I,J,K) +  W(I+1,J,K) )
    TXIM   = ( -  T(I,J,K) +  T(I+1,J,K) )
    AKXIM  = ( - AK(I,J,K) + AK(I+1,J,K) )
    EPSXIM = ( -EPS(I,J,K) +EPS(I+1,J,K) )
    UETM   = 0.5 * ( - 0.5*(   U(I,J-1,K) +   U(I+1,J-1,K) ) &
    &                + 0.5*(   U(I,J+1,K) +   U(I+1,J+1,K) )  )
    VETM   = 0.5 * ( - 0.5*(   V(I,J-1,K) +   V(I+1,J-1,K) ) &
    &                + 0.5*(   V(I,J+1,K) +   V(I+1,J+1,K) )  )
    WETM   = 0.5 * ( - 0.5*(   W(I,J-1,K) +   W(I+1,J-1,K) ) &
    &                + 0.5*(   W(I,J+1,K) +   W(I+1,J+1,K) )  )
    TETM   = 0.5 * ( - 0.5*(   T(I,J-1,K) +   T(I+1,J-1,K) ) &
    &                + 0.5*(   T(I,J+1,K) +   T(I+1,J+1,K) )  )
    AKETM  = 0.5 * ( - 0.5*(  AK(I,J-1,K) +  AK(I+1,J-1,K) ) &
    &                + 0.5*(  AK(I,J+1,K) +  AK(I+1,J+1,K) )  )
    EPSETM = 0.5 * ( - 0.5*( EPS(I,J-1,K) + EPS(I+1,J-1,K) ) &
    &                + 0.5*( EPS(I,J+1,K) + EPS(I+1,J+1,K) )  )
    UZEM   = 0.5 * ( - 0.5*(   U(I,J,K-1) +   U(I+1,J,K-1) ) &
    &                + 0.5*(   U(I,J,K+1) +   U(I+1,J,K+1) )  )
    VZEM   = 0.5 * ( - 0.5*(   V(I,J,K-1) +   V(I+1,J,K-1) ) &
    &                + 0.5*(   V(I,J,K+1) +   V(I+1,J,K+1) )  )
    WZEM   = 0.5 * ( - 0.5*(   W(I,J,K-1) +   W(I+1,J,K-1) ) &
    &                + 0.5*(   W(I,J,K+1) +   W(I+1,J,K+1) )  )
    TZEM   = 0.5 * ( - 0.5*(   T(I,J,K-1) +   T(I+1,J,K-1) ) &
    &                + 0.5*(   T(I,J,K+1) +   T(I+1,J,K+1) )  )
    AKZEM  = 0.5 * ( - 0.5*(  AK(I,J,K-1) +  AK(I+1,J,K-1) ) &
    &                + 0.5*(  AK(I,J,K+1) +  AK(I+1,J,K+1) )  )
    EPSZEM = 0.5 * ( - 0.5*( EPS(I,J,K-1) + EPS(I+1,J,K-1) ) &
    &                + 0.5*( EPS(I,J,K+1) + EPS(I+1,J,K+1) )  )
    UXM   = XIXM *  UXIM + ETXM *  UETM + ZEXM *  UZEM
    UYM   = XIYM *  UXIM + ETYM *  UETM + ZEYM *  UZEM
    UZM   = XIZM *  UXIM + ETZM *  UETM + ZEZM *  UZEM
    VXM   = XIXM *  VXIM + ETXM *  VETM + ZEXM *  VZEM
    VYM   = XIYM *  VXIM + ETYM *  VETM + ZEYM *  VZEM
    VZM   = XIZM *  VXIM + ETZM *  VETM + ZEZM *  VZEM
    WXM   = XIXM *  WXIM + ETXM *  WETM + ZEXM *  WZEM
    WYM   = XIYM *  WXIM + ETYM *  WETM + ZEYM *  WZEM
    WZM   = XIZM *  WXIM + ETZM *  WETM + ZEZM *  WZEM
    TXM   = XIXM *  TXIM + ETXM *  TETM + ZEXM *  TZEM
    TYM   = XIYM *  TXIM + ETYM *  TETM + ZEYM *  TZEM
    TZM   = XIZM *  TXIM + ETZM *  TETM + ZEZM *  TZEM
    AKXM  = XIXM * AKXIM + ETXM * AKETM + ZEXM * AKZEM
    AKYM  = XIYM * AKXIM + ETYM * AKETM + ZEYM * AKZEM
    AKZM  = XIZM * AKXIM + ETZM * AKETM + ZEZM * AKZEM
    EPSXM = XIXM *EPSXIM + ETXM *EPSETM + ZEXM *EPSZEM
    EPSYM = XIYM *EPSXIM + ETYM *EPSETM + ZEYM *EPSZEM
    EPSZM = XIZM *EPSXIM + ETZM *EPSETM + ZEZM *EPSZEM
    ! 粘性係数のXI方向の平均量 -----------------------------------------
    AMUM  = 0.5 * (   AMU(I,J,K) +   AMU(I+1,J,K) )
    AMUTM = 0.5 * (  AMUT(I,J,K) +  AMUT(I+1,J,K) )
    AMUKM = AMUM+AMUTM/SIGK
    AMUEM = AMUM+AMUTM/SIGE
    ! 拡散項ベクトルの計算 ---------------------------------------------
    DELVM = UXM+VYM+WZM
    TAUXX = (AMUM+AMUTM)*(2.0*UXM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUYY = (AMUM+AMUTM)*(2.0*VYM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUZZ = (AMUM+AMUTM)*(2.0*WZM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUXY = (AMUM+AMUTM)*(UYM+VXM)
    TAUYZ = (AMUM+AMUTM)*(VZM+WYM)
    TAUZX = (AMUM+AMUTM)*(WXM+UZM)
    QX    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TXM
    QY    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TYM
    QZ    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TZM
    R5    = TAUXX*UM + TAUXY*VM + TAUZX*WM-QX
    S5    = TAUXY*UM + TAUYY*VM + TAUYZ*WM-QY
    T5    = TAUZX*UM + TAUYZ*VM + TAUZZ*WM-QZ
    RH(I,J,K,1) = 0.0
    RH(I,J,K,2) = (XIXM*TAUXX + XIYM*TAUXY + XIZM*TAUZX)/AJAM
    RH(I,J,K,3) = (XIXM*TAUXY + XIYM*TAUYY + XIZM*TAUYZ)/AJAM
    RH(I,J,K,4) = (XIXM*TAUZX + XIYM*TAUYZ + XIZM*TAUZZ)/AJAM
    RH(I,J,K,5) = (XIXM*R5    + XIYM*S5    + XIZM*T5   )/AJAM
    RH(I,J,K,6) = AMUKM*(XIXM* AKXM +XIYM* AKYM +XIZM* AKZM)/AJAM
    RH(I,J,K,7) = AMUEM*(XIXM*EPSXM +XIYM*EPSYM +XIZM*EPSZM)/AJAM
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISXI
!***********************************************************************
!**** 拡散項(粘性項)の計算(ET方向)                                  ****
!***********************************************************************
SUBROUTINE VISET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: UM, VM, WM, TM, RHOKM
  REAL    :: UXIM, VXIM, WXIM, TXIM, AKXIM, EPSXIM, &
  &          UETM, VETM, WETM, TETM, AKETM, EPSETM, &
  &          UZEM, VZEM, WZEM, TZEM, AKZEM, EPSZEM
  REAL    :: UXM, VXM, WXM, TXM, AKXM, EPSXM, &
  &          UYM, VYM, WYM, TYM, AKYM, EPSYM, &
  &          UZM, VZM, WZM, TZM, AKZM, EPSZM
  REAL    :: AMUM, AMUTM, AMUKM, AMUEM
  REAL    :: DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX
  REAL    :: QX, QY, QZ
  REAL    :: R5, S5, T5
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  UM, VM, WM, TM, RHOKM, &
  !$  UXIM, VXIM, WXIM, TXIM, AKXIM, EPSXIM, &
  !$  UETM, VETM, WETM, TETM, AKETM, EPSETM, &
  !$  UZEM, VZEM, WZEM, TZEM, AKZEM, EPSZEM, &
  !$  UXM, VXM, WXM, TXM, AKXM, EPSXM, &
  !$  UYM, VYM, WYM, TYM, AKYM, EPSYM, &
  !$  UZM, VZM, WZM, TZM, AKZM, EPSZM, &
  !$  AMUM, AMUTM, AMUKM, AMUEM, &
  !$  DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, &
  !$  QX, QY, QZ, R5, S5, T5 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF( QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J+1,K,1) .GT. 0.0 ) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    XIXM  = 0.5 * ( XIX(I,J,K) + XIX(I,J+1,K) )
    XIYM  = 0.5 * ( XIY(I,J,K) + XIY(I,J+1,K) )
    XIZM  = 0.5 * ( XIZ(I,J,K) + XIZ(I,J+1,K) )
    ETXM  = 0.5 * ( ETX(I,J,K) + ETX(I,J+1,K) )
    ETYM  = 0.5 * ( ETY(I,J,K) + ETY(I,J+1,K) )
    ETZM  = 0.5 * ( ETZ(I,J,K) + ETZ(I,J+1,K) )
    ZEXM  = 0.5 * ( ZEX(I,J,K) + ZEX(I,J+1,K) )
    ZEYM  = 0.5 * ( ZEY(I,J,K) + ZEY(I,J+1,K) )
    ZEZM  = 0.5 * ( ZEZ(I,J,K) + ZEZ(I,J+1,K) )
    AJAM  = 0.5 * ( AJA(I,J,K) + AJA(I,J+1,K) )
    UM    = 0.5 * (   U(I,J,K) +   U(I,J+1,K) )
    VM    = 0.5 * (   V(I,J,K) +   V(I,J+1,K) )
    WM    = 0.5 * (   W(I,J,K) +   W(I,J+1,K) )
    TM    = 0.5 * (   T(I,J,K) +   T(I,J+1,K) )
    RHOKM = 0.5 * (  QH(I,J  ,K,6)*AJA(I,J  ,K) &
    &             +  QH(I,J+1,K,6)*AJA(I,J+1,K) )
    ! x,y,z微分の計算(ET方向の平均量) ----------------------------------
    UXIM   = 0.5 * ( - 0.5*(   U(I-1,J,K) +   U(I-1,J+1,K) ) &
    &                + 0.5*(   U(I+1,J,K) +   U(I+1,J+1,K) )  )
    VXIM   = 0.5 * ( - 0.5*(   V(I-1,J,K) +   V(I-1,J+1,K) ) &
    &                + 0.5*(   V(I+1,J,K) +   V(I+1,J+1,K) )  )
    WXIM   = 0.5 * ( - 0.5*(   W(I-1,J,K) +   W(I-1,J+1,K) ) &
    &                + 0.5*(   W(I+1,J,K) +   W(I+1,J+1,K) )  )
    TXIM   = 0.5 * ( - 0.5*(   T(I-1,J,K) +   T(I-1,J+1,K) ) &
    &                + 0.5*(   T(I+1,J,K) +   T(I+1,J+1,K) )  )
    AKXIM  = 0.5 * ( - 0.5*(  AK(I-1,J,K) +  AK(I-1,J+1,K) ) &
    &                + 0.5*(  AK(I+1,J,K) +  AK(I+1,J+1,K) )  )
    EPSXIM = 0.5 * ( - 0.5*( EPS(I-1,J,K) + EPS(I-1,J+1,K) ) &
    &                + 0.5*( EPS(I+1,J,K) + EPS(I+1,J+1,K) )  )
    UETM   = ( -  U(I,J,K) +  U(I,J+1,K) )
    VETM   = ( -  V(I,J,K) +  V(I,J+1,K) )
    WETM   = ( -  W(I,J,K) +  W(I,J+1,K) )
    TETM   = ( -  T(I,J,K) +  T(I,J+1,K) )
    AKETM  = ( - AK(I,J,K) + AK(I,J+1,K) )
    EPSETM = ( -EPS(I,J,K) +EPS(I,J+1,K) )
    UZEM   = 0.5 * ( - 0.5*(   U(I,J,K-1) +   U(I,J+1,K-1) ) &
    &                + 0.5*(   U(I,J,K+1) +   U(I,J+1,K+1) )  )
    VZEM   = 0.5 * ( - 0.5*(   V(I,J,K-1) +   V(I,J+1,K-1) ) &
    &                + 0.5*(   V(I,J,K+1) +   V(I,J+1,K+1) )  )
    WZEM   = 0.5 * ( - 0.5*(   W(I,J,K-1) +   W(I,J+1,K-1) ) &
    &                + 0.5*(   W(I,J,K+1) +   W(I,J+1,K+1) )  )
    TZEM   = 0.5 * ( - 0.5*(   T(I,J,K-1) +   T(I,J+1,K-1) ) &
    &                + 0.5*(   T(I,J,K+1) +   T(I,J+1,K+1) )  )
    AKZEM  = 0.5 * ( - 0.5*(  AK(I,J,K-1) +  AK(I,J+1,K-1) ) &
    &                + 0.5*(  AK(I,J,K+1) +  AK(I,J+1,K+1) )  )
    EPSZEM = 0.5 * ( - 0.5*( EPS(I,J,K-1) + EPS(I,J+1,K-1) ) &
    &                + 0.5*( EPS(I,J,K+1) + EPS(I,J+1,K+1) )  )
    UXM   = XIXM *  UXIM + ETXM *  UETM + ZEXM *  UZEM
    UYM   = XIYM *  UXIM + ETYM *  UETM + ZEYM *  UZEM
    UZM   = XIZM *  UXIM + ETZM *  UETM + ZEZM *  UZEM
    VXM   = XIXM *  VXIM + ETXM *  VETM + ZEXM *  VZEM
    VYM   = XIYM *  VXIM + ETYM *  VETM + ZEYM *  VZEM
    VZM   = XIZM *  VXIM + ETZM *  VETM + ZEZM *  VZEM
    WXM   = XIXM *  WXIM + ETXM *  WETM + ZEXM *  WZEM
    WYM   = XIYM *  WXIM + ETYM *  WETM + ZEYM *  WZEM
    WZM   = XIZM *  WXIM + ETZM *  WETM + ZEZM *  WZEM
    TXM   = XIXM *  TXIM + ETXM *  TETM + ZEXM *  TZEM
    TYM   = XIYM *  TXIM + ETYM *  TETM + ZEYM *  TZEM
    TZM   = XIZM *  TXIM + ETZM *  TETM + ZEZM *  TZEM
    AKXM  = XIXM * AKXIM + ETXM * AKETM + ZEXM * AKZEM
    AKYM  = XIYM * AKXIM + ETYM * AKETM + ZEYM * AKZEM
    AKZM  = XIZM * AKXIM + ETZM * AKETM + ZEZM * AKZEM
    EPSXM = XIXM *EPSXIM + ETXM *EPSETM + ZEXM *EPSZEM
    EPSYM = XIYM *EPSXIM + ETYM *EPSETM + ZEYM *EPSZEM
    EPSZM = XIZM *EPSXIM + ETZM *EPSETM + ZEZM *EPSZEM
    ! 粘性係数のET方向の平均量 -----------------------------------------
    AMUM  = 0.5 * (   AMU(I,J,K) +   AMU(I,J+1,K) )
    AMUTM = 0.5 * (  AMUT(I,J,K) +  AMUT(I,J+1,K) )
    AMUKM = AMUM+AMUTM/SIGK
    AMUEM = AMUM+AMUTM/SIGE
    ! 拡散項ベクトルの計算 ---------------------------------------------
    DELVM = UXM+VYM+WZM
    TAUXX = (AMUM+AMUTM)*(2.0*UXM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUYY = (AMUM+AMUTM)*(2.0*VYM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUZZ = (AMUM+AMUTM)*(2.0*WZM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUXY = (AMUM+AMUTM)*(UYM+VXM)
    TAUYZ = (AMUM+AMUTM)*(VZM+WYM)
    TAUZX = (AMUM+AMUTM)*(WXM+UZM)
    QX    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TXM
    QY    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TYM
    QZ    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TZM
    R5    = TAUXX*UM + TAUXY*VM + TAUZX*WM-QX
    S5    = TAUXY*UM + TAUYY*VM + TAUYZ*WM-QY
    T5    = TAUZX*UM + TAUYZ*VM + TAUZZ*WM-QZ
    SH(I,J,K,1) = 0.0
    SH(I,J,K,2) = (ETXM*TAUXX + ETYM*TAUXY + ETZM*TAUZX)/AJAM
    SH(I,J,K,3) = (ETXM*TAUXY + ETYM*TAUYY + ETZM*TAUYZ)/AJAM
    SH(I,J,K,4) = (ETXM*TAUZX + ETYM*TAUYZ + ETZM*TAUZZ)/AJAM
    SH(I,J,K,5) = (ETXM*R5    + ETYM*S5    + ETZM*T5   )/AJAM
    SH(I,J,K,6) = AMUKM*(ETXM* AKXM +ETYM* AKYM +ETZM* AKZM)/AJAM
    SH(I,J,K,7) = AMUEM*(ETXM*EPSXM +ETYM*EPSYM +ETZM*EPSZM)/AJAM
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISET
!***********************************************************************
!**** 拡散項(粘性項)の計算(ZE方向)                                  ****
!***********************************************************************
SUBROUTINE VISZE
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: UM, VM, WM, TM, RHOKM
  REAL    :: UXIM, VXIM, WXIM, TXIM, AKXIM, EPSXIM, &
  &          UETM, VETM, WETM, TETM, AKETM, EPSETM, &
  &          UZEM, VZEM, WZEM, TZEM, AKZEM, EPSZEM
  REAL    :: UXM, VXM, WXM, TXM, AKXM, EPSXM, &
  &          UYM, VYM, WYM, TYM, AKYM, EPSYM, &
  &          UZM, VZM, WZM, TZM, AKZM, EPSZM
  REAL    :: AMUM, AMUTM, AMUKM, AMUEM
  REAL    :: DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX
  REAL    :: QX, QY, QZ
  REAL    :: R5, S5, T5
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  UM, VM, WM, TM, RHOKM, &
  !$  UXIM, VXIM, WXIM, TXIM, AKXIM, EPSXIM, &
  !$  UETM, VETM, WETM, TETM, AKETM, EPSETM, &
  !$  UZEM, VZEM, WZEM, TZEM, AKZEM, EPSZEM, &
  !$  UXM, VXM, WXM, TXM, AKXM, EPSXM, &
  !$  UYM, VYM, WYM, TYM, AKYM, EPSYM, &
  !$  UZM, VZM, WZM, TZM, AKZM, EPSZM, &
  !$  AMUM, AMUTM, AMUKM, AMUEM, &
  !$  DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, &
  !$  QX, QY, QZ, R5, S5, T5 &
  !$)
  DO K = KS    , KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K+1,1) .GT. 0.0 ) THEN
    ! ZE方向の平均量の計算 ---------------------------------------------
    XIXM  = 0.5 * ( XIX(I,J,K) + XIX(I,J,K+1) )
    XIYM  = 0.5 * ( XIY(I,J,K) + XIY(I,J,K+1) )
    XIZM  = 0.5 * ( XIZ(I,J,K) + XIZ(I,J,K+1) )
    ETXM  = 0.5 * ( ETX(I,J,K) + ETX(I,J,K+1) )
    ETYM  = 0.5 * ( ETY(I,J,K) + ETY(I,J,K+1) )
    ETZM  = 0.5 * ( ETZ(I,J,K) + ETZ(I,J,K+1) )
    ZEXM  = 0.5 * ( ZEX(I,J,K) + ZEX(I,J,K+1) )
    ZEYM  = 0.5 * ( ZEY(I,J,K) + ZEY(I,J,K+1) )
    ZEZM  = 0.5 * ( ZEZ(I,J,K) + ZEZ(I,J,K+1) )
    AJAM  = 0.5 * ( AJA(I,J,K) + AJA(I,J,K+1) )
    UM    = 0.5 * (   U(I,J,K) +   U(I,J,K+1) )
    VM    = 0.5 * (   V(I,J,K) +   V(I,J,K+1) )
    WM    = 0.5 * (   W(I,J,K) +   W(I,J,K+1) )
    TM    = 0.5 * (   T(I,J,K) +   T(I,J,K+1) )
    RHOKM = 0.5 * (  QH(I,J,K  ,6)*AJA(I,J,K  ) &
    &             +  QH(I,J,K+1,6)*AJA(I,J,K+1) )
    ! x,y,z微分の計算(ZE方向の平均量) ----------------------------------
    UXIM   = 0.5 * ( - 0.5*(   U(I-1,J,K) +   U(I-1,J,K+1) ) &
    &                + 0.5*(   U(I+1,J,K) +   U(I+1,J,K+1) )  )
    VXIM   = 0.5 * ( - 0.5*(   V(I-1,J,K) +   V(I-1,J,K+1) ) &
    &                + 0.5*(   V(I+1,J,K) +   V(I+1,J,K+1) )  )
    WXIM   = 0.5 * ( - 0.5*(   W(I-1,J,K) +   W(I-1,J,K+1) ) &
    &                + 0.5*(   W(I+1,J,K) +   W(I+1,J,K+1) )  )
    TXIM   = 0.5 * ( - 0.5*(   T(I-1,J,K) +   T(I-1,J,K+1) ) &
    &                + 0.5*(   T(I+1,J,K) +   T(I+1,J,K+1) )  )
    AKXIM  = 0.5 * ( - 0.5*(  AK(I-1,J,K) +  AK(I-1,J,K+1) ) &
    &                + 0.5*(  AK(I+1,J,K) +  AK(I+1,J,K+1) )  )
    EPSXIM = 0.5 * ( - 0.5*( EPS(I-1,J,K) + EPS(I-1,J,K+1) ) &
    &                + 0.5*( EPS(I+1,J,K) + EPS(I+1,J,K+1) )  )
    UETM   = 0.5 * ( - 0.5*(   U(I,J-1,K) +   U(I,J-1,K+1) ) &
    &                + 0.5*(   U(I,J+1,K) +   U(I,J+1,K+1) )  )
    VETM   = 0.5 * ( - 0.5*(   V(I,J-1,K) +   V(I,J-1,K+1) ) &
    &                + 0.5*(   V(I,J+1,K) +   V(I,J+1,K+1) )  )
    WETM   = 0.5 * ( - 0.5*(   W(I,J-1,K) +   W(I,J-1,K+1) ) &
    &                + 0.5*(   W(I,J+1,K) +   W(I,J+1,K+1) )  )
    TETM   = 0.5 * ( - 0.5*(   T(I,J-1,K) +   T(I,J-1,K+1) ) &
    &                + 0.5*(   T(I,J+1,K) +   T(I,J+1,K+1) )  )
    AKETM  = 0.5 * ( - 0.5*(  AK(I,J-1,K) +  AK(I,J-1,K+1) ) &
    &                + 0.5*(  AK(I,J+1,K) +  AK(I,J+1,K+1) )  )
    EPSETM = 0.5 * ( - 0.5*( EPS(I,J-1,K) + EPS(I,J-1,K+1) ) &
    &                + 0.5*( EPS(I,J+1,K) + EPS(I,J+1,K+1) )  )
    UZEM   = ( -  U(I,J,K) +  U(I,J,K+1) )
    VZEM   = ( -  V(I,J,K) +  V(I,J,K+1) )
    WZEM   = ( -  W(I,J,K) +  W(I,J,K+1) )
    TZEM   = ( -  T(I,J,K) +  T(I,J,K+1) )
    AKZEM  = ( - AK(I,J,K) + AK(I,J,K+1) )
    EPSZEM = ( -EPS(I,J,K) +EPS(I,J,K+1) )
    UXM   = XIXM *  UXIM + ETXM *  UETM + ZEXM *  UZEM
    UYM   = XIYM *  UXIM + ETYM *  UETM + ZEYM *  UZEM
    UZM   = XIZM *  UXIM + ETZM *  UETM + ZEZM *  UZEM
    VXM   = XIXM *  VXIM + ETXM *  VETM + ZEXM *  VZEM
    VYM   = XIYM *  VXIM + ETYM *  VETM + ZEYM *  VZEM
    VZM   = XIZM *  VXIM + ETZM *  VETM + ZEZM *  VZEM
    WXM   = XIXM *  WXIM + ETXM *  WETM + ZEXM *  WZEM
    WYM   = XIYM *  WXIM + ETYM *  WETM + ZEYM *  WZEM
    WZM   = XIZM *  WXIM + ETZM *  WETM + ZEZM *  WZEM
    TXM   = XIXM *  TXIM + ETXM *  TETM + ZEXM *  TZEM
    TYM   = XIYM *  TXIM + ETYM *  TETM + ZEYM *  TZEM
    TZM   = XIZM *  TXIM + ETZM *  TETM + ZEZM *  TZEM
    AKXM  = XIXM * AKXIM + ETXM * AKETM + ZEXM * AKZEM
    AKYM  = XIYM * AKXIM + ETYM * AKETM + ZEYM * AKZEM
    AKZM  = XIZM * AKXIM + ETZM * AKETM + ZEZM * AKZEM
    EPSXM = XIXM *EPSXIM + ETXM *EPSETM + ZEXM *EPSZEM
    EPSYM = XIYM *EPSXIM + ETYM *EPSETM + ZEYM *EPSZEM
    EPSZM = XIZM *EPSXIM + ETZM *EPSETM + ZEZM *EPSZEM
    ! 粘性係数のZE方向の平均量 -----------------------------------------
    AMUM  = 0.5 * (   AMU(I,J,K) +   AMU(I,J,K+1) )
    AMUTM = 0.5 * (  AMUT(I,J,K) +  AMUT(I,J,K+1) )
    AMUKM = AMUM+AMUTM/SIGK
    AMUEM = AMUM+AMUTM/SIGE
    ! 拡散項ベクトルの計算 ---------------------------------------------
    DELVM = UXM+VYM+WZM
    TAUXX = (AMUM+AMUTM)*(2.0*UXM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUYY = (AMUM+AMUTM)*(2.0*VYM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUZZ = (AMUM+AMUTM)*(2.0*WZM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUXY = (AMUM+AMUTM)*(UYM+VXM)
    TAUYZ = (AMUM+AMUTM)*(VZM+WYM)
    TAUZX = (AMUM+AMUTM)*(WXM+UZM)
    QX    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TXM
    QY    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TYM
    QZ    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TZM
    R5    = TAUXX*UM + TAUXY*VM + TAUZX*WM-QX
    S5    = TAUXY*UM + TAUYY*VM + TAUYZ*WM-QY
    T5    = TAUZX*UM + TAUYZ*VM + TAUZZ*WM-QZ
    TH(I,J,K,1) = 0.0
    TH(I,J,K,2) = (ZEXM*TAUXX + ZEYM*TAUXY + ZEZM*TAUZX)/AJAM
    TH(I,J,K,3) = (ZEXM*TAUXY + ZEYM*TAUYY + ZEZM*TAUYZ)/AJAM
    TH(I,J,K,4) = (ZEXM*TAUZX + ZEYM*TAUYZ + ZEZM*TAUZZ)/AJAM
    TH(I,J,K,5) = (ZEXM*R5    + ZEYM*S5    + ZEZM*T5   )/AJAM
    TH(I,J,K,6) = AMUKM*(ZEXM* AKXM +ZEYM* AKYM +ZEZM* AKZM)/AJAM
    TH(I,J,K,7) = AMUEM*(ZEXM*EPSXM +ZEYM*EPSYM +ZEZM*EPSZM)/AJAM
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISZE
! 定義終了 *************************************************************
END SUBROUTINE Turbulence3DEvmStd
!***********************************************************************
!**** 乱流モデル : RANS, EVM, Standard k-e Model (1974)             ****
!****                         with Kato-Launder Modification (1993) ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence3DEvmStdKL( &
&            PELIM, RG, GAMMA, PR, PRT, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, AJA, &
&            QH, U, V, W, T, AK, EPS, AMU, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: CMU = 0.09
  REAL, PARAMETER :: CE1 = 1.44, CE2 = 1.92
  REAL, PARAMETER :: SIGK = 1.0, SIGE = 1.3
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR, PRT
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS: IE, JS: JE, KS: KE), &
  &                       XIY(IS: IE, JS: JE, KS: KE), &
  &                       XIZ(IS: IE, JS: JE, KS: KE), &
  &                       ETX(IS: IE, JS: JE, KS: KE), &
  &                       ETY(IS: IE, JS: JE, KS: KE), &
  &                       ETZ(IS: IE, JS: JE, KS: KE), &
  &                       ZEX(IS: IE, JS: JE, KS: KE), &
  &                       ZEY(IS: IE, JS: JE, KS: KE), &
  &                       ZEZ(IS: IE, JS: JE, KS: KE), &
  &                       AJA(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)  :: QH(IS: IE, JS: JE, KS: KE, LS: LE)
  REAL,    INTENT(IN)  :: U(IS: IE, JS: JE, KS: KE), &
  &                       V(IS: IE, JS: JE, KS: KE), &
  &                       W(IS: IE, JS: JE, KS: KE), &
  &                       T(IS: IE, JS: JE, KS: KE), &
  &                       AK(IS: IE, JS: JE, KS: KE), &
  &                       EPS(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)  :: AMU(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(OUT) :: AMUT(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(OUT) :: DQD(IS: IE, JS: JE, KS: KE, LS: LE)
  REAL,    INTENT(OUT) :: DQP(IS: IE, JS: JE, KS: KE, LS: LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :, :), SH(:, :, :, :), TH(:, :, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 7) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(SH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(TH(IS: IE, JS: JE, KS: KE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL NonDiffTerm
  CALL DiffTerm
  ! メモリ解放 =========================================================
  DEALLOCATE(RH,SH,TH)
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  REAL    :: UX, VX, WX, UY, VY, WY, UZ, VZ, WZ
  REAL    :: DELV, SS2, SC, OMG, OS, PRO, EPSPK
  ! 処理開始 ***********************************************************
  ! 渦粘性係数の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, K)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K,7) .GT. ZERO) THEN
    AMUT(I,J,K) = CMU * QH(I,J,K,6)**2 / QH(I,J,K,7) * AJA(I,J,K)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 生成項と散逸項の計算 +++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE, &
  !$  UX, VX, WX, UY, VY, WY, UZ, VZ, WZ, &
  !$  DELV, SS2, SC, OMG, OS, PRO, EPSPK &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( ( QH(I  ,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I-1,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I+1,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J-1,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J+1,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J  ,K-1,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J  ,K+1,1) .GT. 0.0 ) &
  & ) THEN
    ! ひずみ速度テンソルの計算 -----------------------------------------
    UXI = ( -  U(I-1,J  ,K  ) +  U(I+1,J  ,K  ) )*0.5
    VXI = ( -  V(I-1,J  ,K  ) +  V(I+1,J  ,K  ) )*0.5
    WXI = ( -  W(I-1,J  ,K  ) +  W(I+1,J  ,K  ) )*0.5
    UET = ( -  U(I  ,J-1,K  ) +  U(I  ,J+1,K  ) )*0.5
    VET = ( -  V(I  ,J-1,K  ) +  V(I  ,J+1,K  ) )*0.5
    WET = ( -  W(I  ,J-1,K  ) +  W(I  ,J+1,K  ) )*0.5
    UZE = ( -  U(I  ,J  ,K-1) +  U(I  ,J  ,K+1) )*0.5
    VZE = ( -  V(I  ,J  ,K-1) +  V(I  ,J  ,K+1) )*0.5
    WZE = ( -  W(I  ,J  ,K-1) +  W(I  ,J  ,K+1) )*0.5
    UX  = XIX(I,J,K)*UXI + ETX(I,J,K)*UET + ZEX(I,J,K)*UZE
    UY  = XIY(I,J,K)*UXI + ETY(I,J,K)*UET + ZEY(I,J,K)*UZE
    UZ  = XIZ(I,J,K)*UXI + ETZ(I,J,K)*UET + ZEZ(I,J,K)*UZE
    VX  = XIX(I,J,K)*VXI + ETX(I,J,K)*VET + ZEX(I,J,K)*VZE
    VY  = XIY(I,J,K)*VXI + ETY(I,J,K)*VET + ZEY(I,J,K)*VZE
    VZ  = XIZ(I,J,K)*VXI + ETZ(I,J,K)*VET + ZEZ(I,J,K)*VZE
    WX  = XIX(I,J,K)*WXI + ETX(I,J,K)*WET + ZEX(I,J,K)*WZE
    WY  = XIY(I,J,K)*WXI + ETY(I,J,K)*WET + ZEY(I,J,K)*WZE
    WZ  = XIZ(I,J,K)*WXI + ETZ(I,J,K)*WET + ZEZ(I,J,K)*WZE
    ! 生成項の計算 -----------------------------------------------------
    DELV = UX + VY + WZ
    SS2 = (UY + VX)**2 + (VZ + WY)**2 + (WX + UZ)**2 &
    &   + 2.0 * (UX**2 + VY**2 + WZ**2) &
    &   - 2.0 / 3.0 * DELV**2
    SC  = SQRT(SS2)
    OMG = SQRT((UY - VX)**2 + (VZ - WY)**2 + (WX - UZ)**2)
    ! OMG / SC > 1 の場合、Kato-Launderの修正を適用しない
    IF(SC .GT. ZERO) THEN
      OS = MIN(OMG / SC, 1.0)
    ELSE
      OS = 1.0
    ENDIF
    PRO = OS * ( &
    &     AMUT(I, J, K) * SS2 &
    &   - 2.0 / 3.0 * DELV * QH(I, J, K, 6)*AJA(I, J, K) &
    &   )
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      PRO = MAX(0.0, MIN(QH(I,J,K,7) * AJA(I,J,K) / PELIM, PRO))
    ENDIF
    ! 生成と散逸の和 ---------------------------------------------------
    EPSPK = MIN( &
    &       MAX( &
    &         SC, SQRT(QH(I,J,K,7) * AJA(I,J,K) / AMU(I,J,K)) &
    &       ), &
    &       QH(I,J,K,7) / MAX(ZERO, QH(I,J,K,6)) &
    &     )
    DQP(I,J,K,6) = PRO/AJA(I,J,K) - QH(I,J,K,7)
    DQP(I,J,K,7) = (CE1*PRO/AJA(I,J,K)-CE2*QH(I,J,K,7))*EPSPK
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  ! 処理開始 ***********************************************************
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL VISXI
  CALL VISET
  CALL VISZE
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (QH(I,J,K-1,1) .GT. 0.0) .AND. (QH(I,J-1,K,1) .GT. 0.0) .AND. &
  &   (QH(I-1,J,K,1) .GT. 0.0) .AND. (QH(I,J,K,1)   .GT. 0.0) .AND. &
  &   (QH(I+1,J,K,1) .GT. 0.0) .AND. (QH(I,J+1,K,1) .GT. 0.0) .AND. &
  &   (QH(I,J,K+1,1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,K,L) = (-RH(I-1,J  ,K  ,L) + RH(I,J,K,L) ) &
    &            + (-SH(I  ,J-1,K  ,L) + SH(I,J,K,L) ) &
    &            + (-TH(I  ,J  ,K-1,L) + TH(I,J,K,L) )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項(粘性項)の計算(XI方向)                                  ****
!***********************************************************************
SUBROUTINE VISXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: UM, VM, WM, TM, RHOKM
  REAL    :: UXIM, VXIM, WXIM, TXIM, AKXIM, EPSXIM, &
  &          UETM, VETM, WETM, TETM, AKETM, EPSETM, &
  &          UZEM, VZEM, WZEM, TZEM, AKZEM, EPSZEM
  REAL    :: UXM, VXM, WXM, TXM, AKXM, EPSXM, &
  &          UYM, VYM, WYM, TYM, AKYM, EPSYM, &
  &          UZM, VZM, WZM, TZM, AKZM, EPSZM
  REAL    :: AMUM, AMUTM, AMUKM, AMUEM
  REAL    :: DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX
  REAL    :: QX, QY, QZ
  REAL    :: R5, S5, T5
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  UM, VM, WM, TM, RHOKM, &
  !$  UXIM, VXIM, WXIM, TXIM, AKXIM, EPSXIM, &
  !$  UETM, VETM, WETM, TETM, AKETM, EPSETM, &
  !$  UZEM, VZEM, WZEM, TZEM, AKZEM, EPSZEM, &
  !$  UXM, VXM, WXM, TXM, AKXM, EPSXM, &
  !$  UYM, VYM, WYM, TYM, AKYM, EPSYM, &
  !$  UZM, VZM, WZM, TZM, AKZM, EPSZM, &
  !$  AMUM, AMUTM, AMUKM, AMUEM, &
  !$  DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, &
  !$  QX, QY, QZ, R5, S5, T5 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF( QH(I,J,K,1) .GT. 0.0 .AND. QH(I+1,J,K,1) .GT. 0.0 ) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    XIXM  = 0.5 * ( XIX(I,J,K) + XIX(I+1,J,K) )
    XIYM  = 0.5 * ( XIY(I,J,K) + XIY(I+1,J,K) )
    XIZM  = 0.5 * ( XIZ(I,J,K) + XIZ(I+1,J,K) )
    ETXM  = 0.5 * ( ETX(I,J,K) + ETX(I+1,J,K) )
    ETYM  = 0.5 * ( ETY(I,J,K) + ETY(I+1,J,K) )
    ETZM  = 0.5 * ( ETZ(I,J,K) + ETZ(I+1,J,K) )
    ZEXM  = 0.5 * ( ZEX(I,J,K) + ZEX(I+1,J,K) )
    ZEYM  = 0.5 * ( ZEY(I,J,K) + ZEY(I+1,J,K) )
    ZEZM  = 0.5 * ( ZEZ(I,J,K) + ZEZ(I+1,J,K) )
    AJAM  = 0.5 * ( AJA(I,J,K) + AJA(I+1,J,K) )
    UM    = 0.5 * (   U(I,J,K) +   U(I+1,J,K) )
    VM    = 0.5 * (   V(I,J,K) +   V(I+1,J,K) )
    WM    = 0.5 * (   W(I,J,K) +   W(I+1,J,K) )
    TM    = 0.5 * (   T(I,J,K) +   T(I+1,J,K) )
    RHOKM = 0.5 * (  QH(I  ,J,K,6)*AJA(I  ,J,K) &
    &             +  QH(I+1,J,K,6)*AJA(I+1,J,K) )
    ! x,y,z微分の計算(XI方向の平均量) ----------------------------------
    UXIM   = ( -  U(I,J,K) +  U(I+1,J,K) )
    VXIM   = ( -  V(I,J,K) +  V(I+1,J,K) )
    WXIM   = ( -  W(I,J,K) +  W(I+1,J,K) )
    TXIM   = ( -  T(I,J,K) +  T(I+1,J,K) )
    AKXIM  = ( - AK(I,J,K) + AK(I+1,J,K) )
    EPSXIM = ( -EPS(I,J,K) +EPS(I+1,J,K) )
    UETM   = 0.5 * ( - 0.5*(   U(I,J-1,K) +   U(I+1,J-1,K) ) &
    &                + 0.5*(   U(I,J+1,K) +   U(I+1,J+1,K) )  )
    VETM   = 0.5 * ( - 0.5*(   V(I,J-1,K) +   V(I+1,J-1,K) ) &
    &                + 0.5*(   V(I,J+1,K) +   V(I+1,J+1,K) )  )
    WETM   = 0.5 * ( - 0.5*(   W(I,J-1,K) +   W(I+1,J-1,K) ) &
    &                + 0.5*(   W(I,J+1,K) +   W(I+1,J+1,K) )  )
    TETM   = 0.5 * ( - 0.5*(   T(I,J-1,K) +   T(I+1,J-1,K) ) &
    &                + 0.5*(   T(I,J+1,K) +   T(I+1,J+1,K) )  )
    AKETM  = 0.5 * ( - 0.5*(  AK(I,J-1,K) +  AK(I+1,J-1,K) ) &
    &                + 0.5*(  AK(I,J+1,K) +  AK(I+1,J+1,K) )  )
    EPSETM = 0.5 * ( - 0.5*( EPS(I,J-1,K) + EPS(I+1,J-1,K) ) &
    &                + 0.5*( EPS(I,J+1,K) + EPS(I+1,J+1,K) )  )
    UZEM   = 0.5 * ( - 0.5*(   U(I,J,K-1) +   U(I+1,J,K-1) ) &
    &                + 0.5*(   U(I,J,K+1) +   U(I+1,J,K+1) )  )
    VZEM   = 0.5 * ( - 0.5*(   V(I,J,K-1) +   V(I+1,J,K-1) ) &
    &                + 0.5*(   V(I,J,K+1) +   V(I+1,J,K+1) )  )
    WZEM   = 0.5 * ( - 0.5*(   W(I,J,K-1) +   W(I+1,J,K-1) ) &
    &                + 0.5*(   W(I,J,K+1) +   W(I+1,J,K+1) )  )
    TZEM   = 0.5 * ( - 0.5*(   T(I,J,K-1) +   T(I+1,J,K-1) ) &
    &                + 0.5*(   T(I,J,K+1) +   T(I+1,J,K+1) )  )
    AKZEM  = 0.5 * ( - 0.5*(  AK(I,J,K-1) +  AK(I+1,J,K-1) ) &
    &                + 0.5*(  AK(I,J,K+1) +  AK(I+1,J,K+1) )  )
    EPSZEM = 0.5 * ( - 0.5*( EPS(I,J,K-1) + EPS(I+1,J,K-1) ) &
    &                + 0.5*( EPS(I,J,K+1) + EPS(I+1,J,K+1) )  )
    UXM   = XIXM *  UXIM + ETXM *  UETM + ZEXM *  UZEM
    UYM   = XIYM *  UXIM + ETYM *  UETM + ZEYM *  UZEM
    UZM   = XIZM *  UXIM + ETZM *  UETM + ZEZM *  UZEM
    VXM   = XIXM *  VXIM + ETXM *  VETM + ZEXM *  VZEM
    VYM   = XIYM *  VXIM + ETYM *  VETM + ZEYM *  VZEM
    VZM   = XIZM *  VXIM + ETZM *  VETM + ZEZM *  VZEM
    WXM   = XIXM *  WXIM + ETXM *  WETM + ZEXM *  WZEM
    WYM   = XIYM *  WXIM + ETYM *  WETM + ZEYM *  WZEM
    WZM   = XIZM *  WXIM + ETZM *  WETM + ZEZM *  WZEM
    TXM   = XIXM *  TXIM + ETXM *  TETM + ZEXM *  TZEM
    TYM   = XIYM *  TXIM + ETYM *  TETM + ZEYM *  TZEM
    TZM   = XIZM *  TXIM + ETZM *  TETM + ZEZM *  TZEM
    AKXM  = XIXM * AKXIM + ETXM * AKETM + ZEXM * AKZEM
    AKYM  = XIYM * AKXIM + ETYM * AKETM + ZEYM * AKZEM
    AKZM  = XIZM * AKXIM + ETZM * AKETM + ZEZM * AKZEM
    EPSXM = XIXM *EPSXIM + ETXM *EPSETM + ZEXM *EPSZEM
    EPSYM = XIYM *EPSXIM + ETYM *EPSETM + ZEYM *EPSZEM
    EPSZM = XIZM *EPSXIM + ETZM *EPSETM + ZEZM *EPSZEM
    ! 粘性係数のXI方向の平均量 -----------------------------------------
    AMUM  = 0.5 * (   AMU(I,J,K) +   AMU(I+1,J,K) )
    AMUTM = 0.5 * (  AMUT(I,J,K) +  AMUT(I+1,J,K) )
    AMUKM = AMUM+AMUTM/SIGK
    AMUEM = AMUM+AMUTM/SIGE
    ! 拡散項ベクトルの計算 ---------------------------------------------
    DELVM = UXM+VYM+WZM
    TAUXX = (AMUM+AMUTM)*(2.0*UXM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUYY = (AMUM+AMUTM)*(2.0*VYM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUZZ = (AMUM+AMUTM)*(2.0*WZM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUXY = (AMUM+AMUTM)*(UYM+VXM)
    TAUYZ = (AMUM+AMUTM)*(VZM+WYM)
    TAUZX = (AMUM+AMUTM)*(WXM+UZM)
    QX    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TXM
    QY    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TYM
    QZ    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TZM
    R5    = TAUXX*UM + TAUXY*VM + TAUZX*WM-QX
    S5    = TAUXY*UM + TAUYY*VM + TAUYZ*WM-QY
    T5    = TAUZX*UM + TAUYZ*VM + TAUZZ*WM-QZ
    RH(I,J,K,1) = 0.0
    RH(I,J,K,2) = (XIXM*TAUXX + XIYM*TAUXY + XIZM*TAUZX)/AJAM
    RH(I,J,K,3) = (XIXM*TAUXY + XIYM*TAUYY + XIZM*TAUYZ)/AJAM
    RH(I,J,K,4) = (XIXM*TAUZX + XIYM*TAUYZ + XIZM*TAUZZ)/AJAM
    RH(I,J,K,5) = (XIXM*R5    + XIYM*S5    + XIZM*T5   )/AJAM
    RH(I,J,K,6) = AMUKM*(XIXM* AKXM +XIYM* AKYM +XIZM* AKZM)/AJAM
    RH(I,J,K,7) = AMUEM*(XIXM*EPSXM +XIYM*EPSYM +XIZM*EPSZM)/AJAM
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISXI
!***********************************************************************
!**** 拡散項(粘性項)の計算(ET方向)                                  ****
!***********************************************************************
SUBROUTINE VISET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: UM, VM, WM, TM, RHOKM
  REAL    :: UXIM, VXIM, WXIM, TXIM, AKXIM, EPSXIM, &
  &          UETM, VETM, WETM, TETM, AKETM, EPSETM, &
  &          UZEM, VZEM, WZEM, TZEM, AKZEM, EPSZEM
  REAL    :: UXM, VXM, WXM, TXM, AKXM, EPSXM, &
  &          UYM, VYM, WYM, TYM, AKYM, EPSYM, &
  &          UZM, VZM, WZM, TZM, AKZM, EPSZM
  REAL    :: AMUM, AMUTM, AMUKM, AMUEM
  REAL    :: DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX
  REAL    :: QX, QY, QZ
  REAL    :: R5, S5, T5
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  UM, VM, WM, TM, RHOKM, &
  !$  UXIM, VXIM, WXIM, TXIM, AKXIM, EPSXIM, &
  !$  UETM, VETM, WETM, TETM, AKETM, EPSETM, &
  !$  UZEM, VZEM, WZEM, TZEM, AKZEM, EPSZEM, &
  !$  UXM, VXM, WXM, TXM, AKXM, EPSXM, &
  !$  UYM, VYM, WYM, TYM, AKYM, EPSYM, &
  !$  UZM, VZM, WZM, TZM, AKZM, EPSZM, &
  !$  AMUM, AMUTM, AMUKM, AMUEM, &
  !$  DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, &
  !$  QX, QY, QZ, R5, S5, T5 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF( QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J+1,K,1) .GT. 0.0 ) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    XIXM  = 0.5 * ( XIX(I,J,K) + XIX(I,J+1,K) )
    XIYM  = 0.5 * ( XIY(I,J,K) + XIY(I,J+1,K) )
    XIZM  = 0.5 * ( XIZ(I,J,K) + XIZ(I,J+1,K) )
    ETXM  = 0.5 * ( ETX(I,J,K) + ETX(I,J+1,K) )
    ETYM  = 0.5 * ( ETY(I,J,K) + ETY(I,J+1,K) )
    ETZM  = 0.5 * ( ETZ(I,J,K) + ETZ(I,J+1,K) )
    ZEXM  = 0.5 * ( ZEX(I,J,K) + ZEX(I,J+1,K) )
    ZEYM  = 0.5 * ( ZEY(I,J,K) + ZEY(I,J+1,K) )
    ZEZM  = 0.5 * ( ZEZ(I,J,K) + ZEZ(I,J+1,K) )
    AJAM  = 0.5 * ( AJA(I,J,K) + AJA(I,J+1,K) )
    UM    = 0.5 * (   U(I,J,K) +   U(I,J+1,K) )
    VM    = 0.5 * (   V(I,J,K) +   V(I,J+1,K) )
    WM    = 0.5 * (   W(I,J,K) +   W(I,J+1,K) )
    TM    = 0.5 * (   T(I,J,K) +   T(I,J+1,K) )
    RHOKM = 0.5 * (  QH(I,J  ,K,6)*AJA(I,J  ,K) &
    &             +  QH(I,J+1,K,6)*AJA(I,J+1,K) )
    ! x,y,z微分の計算(ET方向の平均量) ----------------------------------
    UXIM   = 0.5 * ( - 0.5*(   U(I-1,J,K) +   U(I-1,J+1,K) ) &
    &                + 0.5*(   U(I+1,J,K) +   U(I+1,J+1,K) )  )
    VXIM   = 0.5 * ( - 0.5*(   V(I-1,J,K) +   V(I-1,J+1,K) ) &
    &                + 0.5*(   V(I+1,J,K) +   V(I+1,J+1,K) )  )
    WXIM   = 0.5 * ( - 0.5*(   W(I-1,J,K) +   W(I-1,J+1,K) ) &
    &                + 0.5*(   W(I+1,J,K) +   W(I+1,J+1,K) )  )
    TXIM   = 0.5 * ( - 0.5*(   T(I-1,J,K) +   T(I-1,J+1,K) ) &
    &                + 0.5*(   T(I+1,J,K) +   T(I+1,J+1,K) )  )
    AKXIM  = 0.5 * ( - 0.5*(  AK(I-1,J,K) +  AK(I-1,J+1,K) ) &
    &                + 0.5*(  AK(I+1,J,K) +  AK(I+1,J+1,K) )  )
    EPSXIM = 0.5 * ( - 0.5*( EPS(I-1,J,K) + EPS(I-1,J+1,K) ) &
    &                + 0.5*( EPS(I+1,J,K) + EPS(I+1,J+1,K) )  )
    UETM   = ( -  U(I,J,K) +  U(I,J+1,K) )
    VETM   = ( -  V(I,J,K) +  V(I,J+1,K) )
    WETM   = ( -  W(I,J,K) +  W(I,J+1,K) )
    TETM   = ( -  T(I,J,K) +  T(I,J+1,K) )
    AKETM  = ( - AK(I,J,K) + AK(I,J+1,K) )
    EPSETM = ( -EPS(I,J,K) +EPS(I,J+1,K) )
    UZEM   = 0.5 * ( - 0.5*(   U(I,J,K-1) +   U(I,J+1,K-1) ) &
    &                + 0.5*(   U(I,J,K+1) +   U(I,J+1,K+1) )  )
    VZEM   = 0.5 * ( - 0.5*(   V(I,J,K-1) +   V(I,J+1,K-1) ) &
    &                + 0.5*(   V(I,J,K+1) +   V(I,J+1,K+1) )  )
    WZEM   = 0.5 * ( - 0.5*(   W(I,J,K-1) +   W(I,J+1,K-1) ) &
    &                + 0.5*(   W(I,J,K+1) +   W(I,J+1,K+1) )  )
    TZEM   = 0.5 * ( - 0.5*(   T(I,J,K-1) +   T(I,J+1,K-1) ) &
    &                + 0.5*(   T(I,J,K+1) +   T(I,J+1,K+1) )  )
    AKZEM  = 0.5 * ( - 0.5*(  AK(I,J,K-1) +  AK(I,J+1,K-1) ) &
    &                + 0.5*(  AK(I,J,K+1) +  AK(I,J+1,K+1) )  )
    EPSZEM = 0.5 * ( - 0.5*( EPS(I,J,K-1) + EPS(I,J+1,K-1) ) &
    &                + 0.5*( EPS(I,J,K+1) + EPS(I,J+1,K+1) )  )
    UXM   = XIXM *  UXIM + ETXM *  UETM + ZEXM *  UZEM
    UYM   = XIYM *  UXIM + ETYM *  UETM + ZEYM *  UZEM
    UZM   = XIZM *  UXIM + ETZM *  UETM + ZEZM *  UZEM
    VXM   = XIXM *  VXIM + ETXM *  VETM + ZEXM *  VZEM
    VYM   = XIYM *  VXIM + ETYM *  VETM + ZEYM *  VZEM
    VZM   = XIZM *  VXIM + ETZM *  VETM + ZEZM *  VZEM
    WXM   = XIXM *  WXIM + ETXM *  WETM + ZEXM *  WZEM
    WYM   = XIYM *  WXIM + ETYM *  WETM + ZEYM *  WZEM
    WZM   = XIZM *  WXIM + ETZM *  WETM + ZEZM *  WZEM
    TXM   = XIXM *  TXIM + ETXM *  TETM + ZEXM *  TZEM
    TYM   = XIYM *  TXIM + ETYM *  TETM + ZEYM *  TZEM
    TZM   = XIZM *  TXIM + ETZM *  TETM + ZEZM *  TZEM
    AKXM  = XIXM * AKXIM + ETXM * AKETM + ZEXM * AKZEM
    AKYM  = XIYM * AKXIM + ETYM * AKETM + ZEYM * AKZEM
    AKZM  = XIZM * AKXIM + ETZM * AKETM + ZEZM * AKZEM
    EPSXM = XIXM *EPSXIM + ETXM *EPSETM + ZEXM *EPSZEM
    EPSYM = XIYM *EPSXIM + ETYM *EPSETM + ZEYM *EPSZEM
    EPSZM = XIZM *EPSXIM + ETZM *EPSETM + ZEZM *EPSZEM
    ! 粘性係数のET方向の平均量 -----------------------------------------
    AMUM  = 0.5 * (   AMU(I,J,K) +   AMU(I,J+1,K) )
    AMUTM = 0.5 * (  AMUT(I,J,K) +  AMUT(I,J+1,K) )
    AMUKM = AMUM+AMUTM/SIGK
    AMUEM = AMUM+AMUTM/SIGE
    ! 拡散項ベクトルの計算 ---------------------------------------------
    DELVM = UXM+VYM+WZM
    TAUXX = (AMUM+AMUTM)*(2.0*UXM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUYY = (AMUM+AMUTM)*(2.0*VYM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUZZ = (AMUM+AMUTM)*(2.0*WZM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUXY = (AMUM+AMUTM)*(UYM+VXM)
    TAUYZ = (AMUM+AMUTM)*(VZM+WYM)
    TAUZX = (AMUM+AMUTM)*(WXM+UZM)
    QX    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TXM
    QY    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TYM
    QZ    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TZM
    R5    = TAUXX*UM + TAUXY*VM + TAUZX*WM-QX
    S5    = TAUXY*UM + TAUYY*VM + TAUYZ*WM-QY
    T5    = TAUZX*UM + TAUYZ*VM + TAUZZ*WM-QZ
    SH(I,J,K,1) = 0.0
    SH(I,J,K,2) = (ETXM*TAUXX + ETYM*TAUXY + ETZM*TAUZX)/AJAM
    SH(I,J,K,3) = (ETXM*TAUXY + ETYM*TAUYY + ETZM*TAUYZ)/AJAM
    SH(I,J,K,4) = (ETXM*TAUZX + ETYM*TAUYZ + ETZM*TAUZZ)/AJAM
    SH(I,J,K,5) = (ETXM*R5    + ETYM*S5    + ETZM*T5   )/AJAM
    SH(I,J,K,6) = AMUKM*(ETXM* AKXM +ETYM* AKYM +ETZM* AKZM)/AJAM
    SH(I,J,K,7) = AMUEM*(ETXM*EPSXM +ETYM*EPSYM +ETZM*EPSZM)/AJAM
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISET
!***********************************************************************
!**** 拡散項(粘性項)の計算(ZE方向)                                  ****
!***********************************************************************
SUBROUTINE VISZE
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: UM, VM, WM, TM, RHOKM
  REAL    :: UXIM, VXIM, WXIM, TXIM, AKXIM, EPSXIM, &
  &          UETM, VETM, WETM, TETM, AKETM, EPSETM, &
  &          UZEM, VZEM, WZEM, TZEM, AKZEM, EPSZEM
  REAL    :: UXM, VXM, WXM, TXM, AKXM, EPSXM, &
  &          UYM, VYM, WYM, TYM, AKYM, EPSYM, &
  &          UZM, VZM, WZM, TZM, AKZM, EPSZM
  REAL    :: AMUM, AMUTM, AMUKM, AMUEM
  REAL    :: DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX
  REAL    :: QX, QY, QZ
  REAL    :: R5, S5, T5
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  UM, VM, WM, TM, RHOKM, &
  !$  UXIM, VXIM, WXIM, TXIM, AKXIM, EPSXIM, &
  !$  UETM, VETM, WETM, TETM, AKETM, EPSETM, &
  !$  UZEM, VZEM, WZEM, TZEM, AKZEM, EPSZEM, &
  !$  UXM, VXM, WXM, TXM, AKXM, EPSXM, &
  !$  UYM, VYM, WYM, TYM, AKYM, EPSYM, &
  !$  UZM, VZM, WZM, TZM, AKZM, EPSZM, &
  !$  AMUM, AMUTM, AMUKM, AMUEM, &
  !$  DELVM, TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, &
  !$  QX, QY, QZ, R5, S5, T5 &
  !$)
  DO K = KS    , KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( QH(I,J,K,1) .GT. 0.0 .AND. QH(I,J,K+1,1) .GT. 0.0 ) THEN
    ! ZE方向の平均量の計算 ---------------------------------------------
    XIXM  = 0.5 * ( XIX(I,J,K) + XIX(I,J,K+1) )
    XIYM  = 0.5 * ( XIY(I,J,K) + XIY(I,J,K+1) )
    XIZM  = 0.5 * ( XIZ(I,J,K) + XIZ(I,J,K+1) )
    ETXM  = 0.5 * ( ETX(I,J,K) + ETX(I,J,K+1) )
    ETYM  = 0.5 * ( ETY(I,J,K) + ETY(I,J,K+1) )
    ETZM  = 0.5 * ( ETZ(I,J,K) + ETZ(I,J,K+1) )
    ZEXM  = 0.5 * ( ZEX(I,J,K) + ZEX(I,J,K+1) )
    ZEYM  = 0.5 * ( ZEY(I,J,K) + ZEY(I,J,K+1) )
    ZEZM  = 0.5 * ( ZEZ(I,J,K) + ZEZ(I,J,K+1) )
    AJAM  = 0.5 * ( AJA(I,J,K) + AJA(I,J,K+1) )
    UM    = 0.5 * (   U(I,J,K) +   U(I,J,K+1) )
    VM    = 0.5 * (   V(I,J,K) +   V(I,J,K+1) )
    WM    = 0.5 * (   W(I,J,K) +   W(I,J,K+1) )
    TM    = 0.5 * (   T(I,J,K) +   T(I,J,K+1) )
    RHOKM = 0.5 * (  QH(I,J,K  ,6)*AJA(I,J,K  ) &
    &             +  QH(I,J,K+1,6)*AJA(I,J,K+1) )
    ! x,y,z微分の計算(ZE方向の平均量) ----------------------------------
    UXIM   = 0.5 * ( - 0.5*(   U(I-1,J,K) +   U(I-1,J,K+1) ) &
    &                + 0.5*(   U(I+1,J,K) +   U(I+1,J,K+1) )  )
    VXIM   = 0.5 * ( - 0.5*(   V(I-1,J,K) +   V(I-1,J,K+1) ) &
    &                + 0.5*(   V(I+1,J,K) +   V(I+1,J,K+1) )  )
    WXIM   = 0.5 * ( - 0.5*(   W(I-1,J,K) +   W(I-1,J,K+1) ) &
    &                + 0.5*(   W(I+1,J,K) +   W(I+1,J,K+1) )  )
    TXIM   = 0.5 * ( - 0.5*(   T(I-1,J,K) +   T(I-1,J,K+1) ) &
    &                + 0.5*(   T(I+1,J,K) +   T(I+1,J,K+1) )  )
    AKXIM  = 0.5 * ( - 0.5*(  AK(I-1,J,K) +  AK(I-1,J,K+1) ) &
    &                + 0.5*(  AK(I+1,J,K) +  AK(I+1,J,K+1) )  )
    EPSXIM = 0.5 * ( - 0.5*( EPS(I-1,J,K) + EPS(I-1,J,K+1) ) &
    &                + 0.5*( EPS(I+1,J,K) + EPS(I+1,J,K+1) )  )
    UETM   = 0.5 * ( - 0.5*(   U(I,J-1,K) +   U(I,J-1,K+1) ) &
    &                + 0.5*(   U(I,J+1,K) +   U(I,J+1,K+1) )  )
    VETM   = 0.5 * ( - 0.5*(   V(I,J-1,K) +   V(I,J-1,K+1) ) &
    &                + 0.5*(   V(I,J+1,K) +   V(I,J+1,K+1) )  )
    WETM   = 0.5 * ( - 0.5*(   W(I,J-1,K) +   W(I,J-1,K+1) ) &
    &                + 0.5*(   W(I,J+1,K) +   W(I,J+1,K+1) )  )
    TETM   = 0.5 * ( - 0.5*(   T(I,J-1,K) +   T(I,J-1,K+1) ) &
    &                + 0.5*(   T(I,J+1,K) +   T(I,J+1,K+1) )  )
    AKETM  = 0.5 * ( - 0.5*(  AK(I,J-1,K) +  AK(I,J-1,K+1) ) &
    &                + 0.5*(  AK(I,J+1,K) +  AK(I,J+1,K+1) )  )
    EPSETM = 0.5 * ( - 0.5*( EPS(I,J-1,K) + EPS(I,J-1,K+1) ) &
    &                + 0.5*( EPS(I,J+1,K) + EPS(I,J+1,K+1) )  )
    UZEM   = ( -  U(I,J,K) +  U(I,J,K+1) )
    VZEM   = ( -  V(I,J,K) +  V(I,J,K+1) )
    WZEM   = ( -  W(I,J,K) +  W(I,J,K+1) )
    TZEM   = ( -  T(I,J,K) +  T(I,J,K+1) )
    AKZEM  = ( - AK(I,J,K) + AK(I,J,K+1) )
    EPSZEM = ( -EPS(I,J,K) +EPS(I,J,K+1) )
    UXM   = XIXM *  UXIM + ETXM *  UETM + ZEXM *  UZEM
    UYM   = XIYM *  UXIM + ETYM *  UETM + ZEYM *  UZEM
    UZM   = XIZM *  UXIM + ETZM *  UETM + ZEZM *  UZEM
    VXM   = XIXM *  VXIM + ETXM *  VETM + ZEXM *  VZEM
    VYM   = XIYM *  VXIM + ETYM *  VETM + ZEYM *  VZEM
    VZM   = XIZM *  VXIM + ETZM *  VETM + ZEZM *  VZEM
    WXM   = XIXM *  WXIM + ETXM *  WETM + ZEXM *  WZEM
    WYM   = XIYM *  WXIM + ETYM *  WETM + ZEYM *  WZEM
    WZM   = XIZM *  WXIM + ETZM *  WETM + ZEZM *  WZEM
    TXM   = XIXM *  TXIM + ETXM *  TETM + ZEXM *  TZEM
    TYM   = XIYM *  TXIM + ETYM *  TETM + ZEYM *  TZEM
    TZM   = XIZM *  TXIM + ETZM *  TETM + ZEZM *  TZEM
    AKXM  = XIXM * AKXIM + ETXM * AKETM + ZEXM * AKZEM
    AKYM  = XIYM * AKXIM + ETYM * AKETM + ZEYM * AKZEM
    AKZM  = XIZM * AKXIM + ETZM * AKETM + ZEZM * AKZEM
    EPSXM = XIXM *EPSXIM + ETXM *EPSETM + ZEXM *EPSZEM
    EPSYM = XIYM *EPSXIM + ETYM *EPSETM + ZEYM *EPSZEM
    EPSZM = XIZM *EPSXIM + ETZM *EPSETM + ZEZM *EPSZEM
    ! 粘性係数のZE方向の平均量 -----------------------------------------
    AMUM  = 0.5 * (   AMU(I,J,K) +   AMU(I,J,K+1) )
    AMUTM = 0.5 * (  AMUT(I,J,K) +  AMUT(I,J,K+1) )
    AMUKM = AMUM+AMUTM/SIGK
    AMUEM = AMUM+AMUTM/SIGE
    ! 拡散項ベクトルの計算 ---------------------------------------------
    DELVM = UXM+VYM+WZM
    TAUXX = (AMUM+AMUTM)*(2.0*UXM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUYY = (AMUM+AMUTM)*(2.0*VYM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUZZ = (AMUM+AMUTM)*(2.0*WZM-2.0/3.0*DELVM) &
    &     - 2.0/3.0*RHOKM
    TAUXY = (AMUM+AMUTM)*(UYM+VXM)
    TAUYZ = (AMUM+AMUTM)*(VZM+WYM)
    TAUZX = (AMUM+AMUTM)*(WXM+UZM)
    QX    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TXM
    QY    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TYM
    QZ    =-(AMUM/PR+AMUTM/PRT)/(GAMMA-1.0)*GAMMA*RG*TZM
    R5    = TAUXX*UM + TAUXY*VM + TAUZX*WM-QX
    S5    = TAUXY*UM + TAUYY*VM + TAUYZ*WM-QY
    T5    = TAUZX*UM + TAUYZ*VM + TAUZZ*WM-QZ
    TH(I,J,K,1) = 0.0
    TH(I,J,K,2) = (ZEXM*TAUXX + ZEYM*TAUXY + ZEZM*TAUZX)/AJAM
    TH(I,J,K,3) = (ZEXM*TAUXY + ZEYM*TAUYY + ZEZM*TAUYZ)/AJAM
    TH(I,J,K,4) = (ZEXM*TAUZX + ZEYM*TAUYZ + ZEZM*TAUZZ)/AJAM
    TH(I,J,K,5) = (ZEXM*R5    + ZEYM*S5    + ZEZM*T5   )/AJAM
    TH(I,J,K,6) = AMUKM*(ZEXM* AKXM +ZEYM* AKYM +ZEZM* AKZM)/AJAM
    TH(I,J,K,7) = AMUEM*(ZEXM*EPSXM +ZEYM*EPSYM +ZEZM*EPSZM)/AJAM
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE VISZE
! 定義終了 *************************************************************
END SUBROUTINE Turbulence3DEvmStdKL
!***********************************************************************
!**** 乱流モデル : RANS, EVM, Launder-Sharma Model (1974)           ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence2DEvmLS( &
&            PELIM, RG, GAMMA, PR, PRT, &
&            IS, IE, JS, JE, LS, LE, &
&            XIX, XIY, ETX, ETY, AJA, &
&            RHO, U, V, T, AK, EPST, AMU, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: CMU = 0.09
  REAL, PARAMETER :: CE1 = 1.44, CE2 = 1.92
  REAL, PARAMETER :: SIGK = 1.0, SIGE = 1.3
  REAL, PARAMETER :: ZERO = 1.0E-20
  REAL, PARAMETER :: TwoThird = 0.66666667
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR, PRT
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE), XIY(IS:IE, JS:JE), &
  &                       ETX(IS:IE, JS:JE), ETY(IS:IE, JS:JE), &
  &                       AJA(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE), &
  &                       U(IS:IE, JS:JE), V(IS:IE, JS:JE), &
  &                       T(IS:IE, JS:JE), &
  &                       AK(IS:IE, JS:JE), EPST(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: AMUT(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(OUT) :: DQP(IS:IE, JS:JE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :), SH(:, :, :)
  REAL, ALLOCATABLE :: UXXI(:, :), UYXI(:, :), VXXI(:, :), VYXI(:, :)
  REAL, ALLOCATABLE :: UXET(:, :), UYET(:, :), VXET(:, :), VYET(:, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 6) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, LS: LE), SH(IS: IE, JS: JE, LS: LE))
  ALLOCATE(UXXI(IS: IE, JS: JE), UYXI(IS: IE, JS: JE))
  ALLOCATE(VXXI(IS: IE, JS: JE), VYXI(IS: IE, JS: JE))
  ALLOCATE(UXET(IS: IE, JS: JE), UYET(IS: IE, JS: JE))
  ALLOCATE(VXET(IS: IE, JS: JE), VYET(IS: IE, JS: JE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL PreVel2ndDiff
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 二階微分計算用の前準備                                        ****
!***********************************************************************
SUBROUTINE PreVel2ndDiff
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  UXXI = 0.0
  UYXI = 0.0
  VXXI = 0.0
  VYXI = 0.0
  UXET = 0.0
  UYET = 0.0
  VXET = 0.0
  VYET = 0.0
  ! 速度の二階微分の前準備 +++++++++++++++++++++++++++++++++++++++++++++
  CALL PreVel2ndDiffXI
  CALL PreVel2ndDiffET
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiff
!***********************************************************************
!**** 二階微分計算用の前準備(XI方向)                                ****
!***********************************************************************
SUBROUTINE PreVel2ndDiffXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM
  REAL    :: UXI, VXI, UET, VET
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, XIXM, XIYM, ETXM, ETYM, UXI, VXI, UET, VET &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF(RHO(I,J) .GT. 0.0 .AND. RHO(I+1,J) .GT. 0.0) THEN
    ! XI方向の中間値 ---------------------------------------------------
    XIXM = 0.5 * (XIX(I,J) + XIX(I+1,J))
    XIYM = 0.5 * (XIY(I,J) + XIY(I+1,J))
    ETXM = 0.5 * (ETX(I,J) + ETX(I+1,J))
    ETYM = 0.5 * (ETY(I,J) + ETY(I+1,J))
    ! 計算空間一階微分 -------------------------------------------------
    UXI = - U(I,J) + U(I+1,J)
    UET = - 0.25 * (U(I,J-1) + U(I+1,J-1) - U(I,J+1) - U(I+1,J+1))
    VXI = - V(I,J) + V(I+1,J)
    VET = - 0.25 * (V(I,J-1) + V(I+1,J-1) - V(I,J+1) - V(I+1,J+1))
    ! 物理空間一階微分 -------------------------------------------------
    UXXI(I,J) = UXI * XIXM + UET * ETXM
    UYXI(I,J) = UXI * XIYM + UET * ETYM
    VXXI(I,J) = VXI * XIXM + VET * ETXM
    VYXI(I,J) = VXI * XIYM + VET * ETYM
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiffXI
!***********************************************************************
!**** 二階微分計算用の前準備(ET方向)                                ****
!***********************************************************************
SUBROUTINE PreVel2ndDiffET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM
  REAL    :: UXI, VXI, UET, VET
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, XIXM, XIYM, ETXM, ETYM, UXI, VXI, UET, VET &
  !$)
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF(RHO(I,J) .GT. 0.0 .AND. RHO(I,J+1) .GT. 0.0) THEN
    ! ET方向の中間値 ---------------------------------------------------
    XIXM = 0.5 * (XIX(I,J) + XIX(I,J+1))
    XIYM = 0.5 * (XIY(I,J) + XIY(I,J+1))
    ETXM = 0.5 * (ETX(I,J) + ETX(I,J+1))
    ETYM = 0.5 * (ETY(I,J) + ETY(I,J+1))
    ! 計算空間一階微分 -------------------------------------------------
    UXI = - 0.25 * ( U(I-1,J) - U(I+1,J) + U(I-1,J+1) - U(I+1,J+1) )
    UET = - U(I,J) + U(I,J+1)
    VXI = - 0.25 * ( V(I-1,J) - V(I+1,J) + V(I-1,J+1) - V(I+1,J+1) )
    VET = - V(I,J) + V(I,J+1)
    ! 物理空間一階微分 -------------------------------------------------
    UXET(I,J) = UXI * XIXM + UET * ETXM
    UYET(I,J) = UXI * XIYM + UET * ETYM
    VXET(I,J) = VXI * XIXM + VET * ETXM
    VYET(I,J) = VXI * XIYM + VET * ETYM
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiffET
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: UXI, VXI, UET, VET
  REAL    :: UX, VX, UY, VY, DELV
  REAL    :: KRXI, KRET, KRX, KRY
  REAL    :: UXX, UXY, UYX, UYY, VXX, VXY, VYX, VYY
  REAL    :: Rtt, f_mu, f_2
  REAL    :: nu, nu_t
  REAL    :: SS2, PP_k, D, E
  REAL    :: QH1, AEAK
  ! 処理開始 ***********************************************************
  ! 渦粘性係数の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, nu, Rtt, f_mu, nu_t &
  !$)
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J) .GT. 0.0) THEN
    IF(EPST(I,J) .GT. ZERO) THEN
      nu        = AMU(I,J) / RHO(I,J)
      Rtt       = AK(I,J)**2 / (nu * EPST(I,J))
      f_mu      = EXP(-3.4 / (1.0 + Rtt / 50.0)**2)
      nu_t      = CMU * f_mu * AK(I,J)**2 / EPST(I,J)
      AMUT(I,J) = RHO(I,J) * nu_t
    ELSE
      AMUT(I,J) = 0.0
    ENDIF
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 非拡散項の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  UXI, VXI, UET, VET, &
  !$  UX, VX, UY, VY, DELV, &
  !$  KRXI, KRET, KRX, KRY, &
  !$  UXX, UXY, UYX, UYY, VXX, VXY, VYX, VYY, &
  !$  Rtt, f_2, nu, nu_t, SS2, PP_k, D, E, QH1, AEAK &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J-1) .GT. 0.0) .AND. (RHO(I-1,J) .GT. 0.0) .AND. &
  &   (RHO(I,J)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0) &
  & ) THEN
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.5 * (U(I+1,J) - U(I-1,J))
    VXI = 0.5 * (V(I+1,J) - V(I-1,J))
    UET = 0.5 * (U(I,J+1) - U(I,J-1))
    VET = 0.5 * (V(I,J+1) - V(I,J-1))
    ! 物理空間方向一階微分
    UX = UXI * XIX(I,J) + UET * ETX(I,J)
    UY = UXI * XIY(I,J) + UET * ETY(I,J)
    VX = VXI * XIX(I,J) + VET * ETX(I,J)
    VY = VXI * XIY(I,J) + VET * ETY(I,J)
    ! 乱れの速度スケール sqrt(k) の一階微分の計算 ----------------------
    ! 計算空間方向一階微分
    KRXI = 0.5 * (SQRT(AK(I+1,J)) - SQRT(AK(I-1,J)))
    KRET = 0.5 * (SQRT(AK(I,J+1)) - SQRT(AK(I,J-1)))
    ! 物理空間方向一階微分
    KRX = KRXI * XIX(I,J) + KRET * ETX(I,J)
    KRY = KRXI * XIY(I,J) + KRET * ETY(I,J)
    ! 速度の二階微分の計算 ---------------------------------------------
    UXX = - (UXXI(I-1,J) - UXXI(I,J)) * XIX(I,J) &
    &     - (UXET(I,J-1) - UXET(I,J)) * ETX(I,J)
    UXY = - (UXXI(I-1,J) - UXXI(I,J)) * XIY(I,J) &
    &     - (UXET(I,J-1) - UXET(I,J)) * ETY(I,J)
    UYX = - (UYXI(I-1,J) - UYXI(I,J)) * XIX(I,J) &
    &     - (UYET(I,J-1) - UYET(I,J)) * ETX(I,J)
    UYY = - (UYXI(I-1,J) - UYXI(I,J)) * XIY(I,J) &
    &     - (UYET(I,J-1) - UYET(I,J)) * ETY(I,J)
    VXX = - (VXXI(I-1,J) - VXXI(I,J)) * XIX(I,J) &
    &     - (VXET(I,J-1) - VXET(I,J)) * ETX(I,J)
    VXY = - (VXXI(I-1,J) - VXXI(I,J)) * XIY(I,J) &
    &     - (VXET(I,J-1) - VXET(I,J)) * ETY(I,J)
    VYX = - (VYXI(I-1,J) - VYXI(I,J)) * XIX(I,J) &
    &     - (VYET(I,J-1) - VYET(I,J)) * ETX(I,J)
    VYY = - (VYXI(I-1,J) - VYXI(I,J)) * XIY(I,J) &
    &     - (VYET(I,J-1) - VYET(I,J)) * ETY(I,J)
    ! モデル関数の計算 -------------------------------------------------
    nu   =  AMU(I,J) / RHO(I,J)
    nu_t = AMUT(I,J) / RHO(I,J)
    IF(EPST(I,J) .GT. ZERO) THEN
      Rtt = AK(I,J)**2 / (nu * EPST(I,J))
    ELSE
      Rtt = 0.0
    ENDIF
    f_2  = 1.0 - 0.3 * EXP(-Rtt**2)
    ! 補正項 D, E の計算 -----------------------------------------------
    D = 2.0 * nu * (KRX**2 + KRY**2)
    E = 2.0 * nu * nu_t * ( &
    &     UXX**2 + VXX**2 + UYY**2 + VYY**2 &
    &   + UXY**2 + VXY**2 + UYX**2 + VYX**2 &
    & )
    ! レイノルズ応力の生成(平均歪み増加率) -----------------------------
    DELV = UX + VY
    SS2  = 2.0 * (UX**2 + VY**2) + (UY + VX)**2 - TwoThird * DELV**2
    PP_k = nu_t * SS2 - TwoThird * AK(I,J) * DELV
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      PP_k = MAX(0.0, MIN((EPST(I,J) + D) / PELIM, PP_k))
    ENDIF
    ! 生成、散逸の和 ---------------------------------------------------
    QH1 = RHO(I,J) / AJA(I,J)
    AEAK = MIN( &
    &      MAX(SQRT(SS2), SQRT((EPST(I,J) + D) / nu)), &
    &      EPST(I,J) / MAX(ZERO, AK(I,J)) &
    &    )
    DQP(I,J,5) = QH1 * (PP_k - (EPST(I,J) + D))
    DQP(I,J,6) = QH1 * (AEAK * (CE1 * PP_k - CE2 * f_2 * EPST(I,J)) + E)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH = 0.0
  SH = 0.0
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL DIFFXI
  CALL DIFFET
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, L)
  DO L = LS, LE
  !$OMP DO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J-1) .GT. 0.0) .AND. (RHO(I-1,J) .GT. 0.0) .AND. &
  &   (RHO(I,J)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,L) = (-RH(I-1,J  ,L) + RH(I,J,L)) &
    &          + (-SH(I  ,J-1,L) + SH(I,J,L))
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項の計算(xi方向差分用)                                    ****
!***********************************************************************
SUBROUTINE DIFFXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: RHOM, UM, VM, TM, AKM, EPSTM, AMUM, AMUTM
  REAL    :: UXI, UET, UX, UY
  REAL    :: VXI, VET, VX, VY
  REAL    :: TXI, TET, TX, TY
  REAL    :: AKXI, AKET, AKX, AKY
  REAL    :: AEXI, AEET, AEX, AEY
  REAL    :: nu, nu_t
  REAL    :: c_u, c_t, c_k, c_e
  REAL    :: TAUXX, TAUYY, TAUXY, DELV
  REAL    :: R4, S4
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, &
  !$  RHOM, UM, VM, TM, AKM, EPSTM, AMUM, AMUTM, &
  !$  UXI, UET, UX, UY, &
  !$  VXI, VET, VX, VY, &
  !$  TXI, TET, TX, TY, &
  !$  AKXI, AKET, AKX, AKY, &
  !$  AEXI, AEET, AEX, AEY, &
  !$  nu, nu_t, c_u, c_t, c_k, c_e, &
  !$  TAUXX, TAUYY, TAUXY, DELV, R4, S4, QH1 &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF( (RHO(I,J-1) .GT. 0.0) .AND. (RHO(I,J)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0)  &
  & ) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J) + XIX(I+1,J))
    XIYM = 0.5 * (XIY(I,J) + XIY(I+1,J))
    ETXM = 0.5 * (ETX(I,J) + ETX(I+1,J))
    ETYM = 0.5 * (ETY(I,J) + ETY(I+1,J))
    AJAM = 0.5 * (AJA(I,J) + AJA(I+1,J))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J) +  RHO(I+1,J))
    UM    = 0.5 * (   U(I,J) +    U(I+1,J))
    VM    = 0.5 * (   V(I,J) +    V(I+1,J))
    TM    = 0.5 * (   T(I,J) +    T(I+1,J))
    AKM   = 0.5 * (  AK(I,J) +   AK(I+1,J))
    EPSTM = 0.5 * (EPST(I,J) + EPST(I+1,J))
    AMUM  = 0.5 * ( AMU(I,J) +  AMU(I+1,J))
    AMUTM = 0.5 * (AMUT(I,J) + AMUT(I+1,J))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = - U(I,J) + U(I+1,J)
    VXI = - V(I,J) + V(I+1,J)
    UET = 0.25 * (- U(I,J-1) - U(I+1,J-1) + U(I,J+1) + U(I+1,J+1))
    VET = 0.25 * (- V(I,J-1) - V(I+1,J-1) + V(I,J+1) + V(I+1,J+1))
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM
    UY = UXI * XIYM + UET * ETYM
    VX = VXI * XIXM + VET * ETXM
    VY = VXI * XIYM + VET * ETYM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = - T(I,J) + T(I+1,J)
    TET = 0.25 * (- T(I,J-1) - T(I+1,J-1) + T(I,J+1) + T(I+1,J+1))
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM
    TY = TXI * XIYM + TET * ETYM
    ! k の一階微分の計算 -----------------------------------------------
    ! 計算空間方向一階微分
    AKXI = - AK(I,J) + AK(I+1,J)
    AKET = 0.25 * (- AK(I,J-1) - AK(I+1,J-1) + AK(I,J+1) + AK(I+1,J+1))
    ! 物理空間方向一階微分
    AKX = AKXI * XIXM + AKET * ETXM
    AKY = AKXI * XIYM + AKET * ETYM
    ! epsilon tilde の一階微分の計算 -----------------------------------
    ! 計算空間方向一階微分
    AEXI = - EPST(I,J) + EPST(I+1,J)
    AEET = 0.25 * (-EPST(I,J-1)-EPST(I+1,J-1)+EPST(I,J+1)+EPST(I+1,J+1))
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM
    AEY = AEXI * XIYM + AEET * ETYM
    ! 係数計算 ---------------------------------------------------------
    nu   = AMUM  / RHOM
    nu_t = AMUTM / RHOM
    c_u  = nu + nu_t
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    c_k  = nu + nu_t / SIGK
    c_e  = nu + nu_t / SIGE
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    DELV  = UX + VY
    TAUXX = c_u * (2.0 * UX - TwoThird * DELV) - TwoThird * AKM
    TAUYY = c_u * (2.0 * VY - TwoThird * DELV) - TwoThird * AKM
    TAUXY = c_u * (UY + VX)
    ! エネルギーの拡散 -------------------------------------------------
    R4 = TAUXX * UM + TAUXY * VM + c_t * TX
    S4 = TAUXY * UM + TAUYY * VM + c_t * TY
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    RH(I,J,1) = 0.0
    RH(I,J,2) = QH1 * (TAUXX * XIXM + TAUXY * XIYM)
    RH(I,J,3) = QH1 * (TAUXY * XIXM + TAUYY * XIYM)
    RH(I,J,4) = QH1 * (R4 * XIXM + S4 * XIYM)
    RH(I,J,5) = QH1 * c_k * (AKX * XIXM + AKY * XIYM)
    RH(I,J,6) = QH1 * c_e * (AEX * XIXM + AEY * XIYM)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFXI
!***********************************************************************
!**** 拡散項の計算(eta方向差分用)                                   ****
!***********************************************************************
SUBROUTINE DIFFET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: RHOM, UM, VM, TM, AKM, EPSTM, AMUM, AMUTM
  REAL    :: UXI, UET, UX, UY
  REAL    :: VXI, VET, VX, VY
  REAL    :: TXI, TET, TX, TY
  REAL    :: AKXI, AKET, AKX, AKY
  REAL    :: AEXI, AEET, AEX, AEY
  REAL    :: nu, nu_t
  REAL    :: c_u, c_t, c_k, c_e
  REAL    :: TAUXX, TAUYY, TAUXY, DELV
  REAL    :: R4, S4
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, &
  !$  RHOM, UM, VM, TM, AKM, EPSTM, AMUM, AMUTM, &
  !$  UXI, UET, UX, UY, &
  !$  VXI, VET, VX, VY, &
  !$  TXI, TET, TX, TY, &
  !$  AKXI, AKET, AKX, AKY, &
  !$  AEXI, AEET, AEX, AEY, &
  !$  nu, nu_t, c_u, c_t, c_k, c_e, &
  !$  TAUXX, TAUYY, TAUXY, DELV, R4, S4, QH1 &
  !$)
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I-1,J) .GT. 0.0) .AND. (RHO(I,J)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0) &
  & ) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J) + XIX(I,J+1))
    XIYM = 0.5 * (XIY(I,J) + XIY(I,J+1))
    ETXM = 0.5 * (ETX(I,J) + ETX(I,J+1))
    ETYM = 0.5 * (ETY(I,J) + ETY(I,J+1))
    AJAM = 0.5 * (AJA(I,J) + AJA(I,J+1))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J) +  RHO(I,J+1))
    UM    = 0.5 * (   U(I,J) +    U(I,J+1))
    VM    = 0.5 * (   V(I,J) +    V(I,J+1))
    TM    = 0.5 * (   T(I,J) +    T(I,J+1))
    AKM   = 0.5 * (  AK(I,J) +   AK(I,J+1))
    EPSTM = 0.5 * (EPST(I,J) + EPST(I,J+1))
    AMUM  = 0.5 * ( AMU(I,J) +  AMU(I,J+1))
    AMUTM = 0.5 * (AMUT(I,J) + AMUT(I,J+1))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * (- U(I-1,J) - U(I-1,J+1) + U(I+1,J) + U(I+1,J+1))
    VXI = 0.25 * (- V(I-1,J) - V(I-1,J+1) + V(I+1,J) + V(I+1,J+1))
    UET = - U(I,J) + U(I,J+1)
    VET = - V(I,J) + V(I,J+1)
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM
    UY = UXI * XIYM + UET * ETYM
    VX = VXI * XIXM + VET * ETXM
    VY = VXI * XIYM + VET * ETYM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * (- T(I-1,J) - T(I-1,J+1) + T(I+1,J) + T(I+1,J+1))
    TET = - T(I,J) + T(I,J+1)
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM
    TY = TXI * XIYM + TET * ETYM
    ! k の一階微分の計算 -----------------------------------------------
    ! 計算空間方向一階微分
    AKXI = 0.25 * (- AK(I-1,J) - AK(I-1,J+1) + AK(I+1,J) + AK(I+1,J+1))
    AKET = - AK(I,J) + AK(I,J+1)
    ! 物理空間方向一階微分
    AKX = AKXI * XIXM + AKET * ETXM
    AKY = AKXI * XIYM + AKET * ETYM
    ! epsilon tilde の一階微分の計算 -----------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * (-EPST(I-1,J)-EPST(I-1,J+1)+EPST(I+1,J)+EPST(I+1,J+1))
    AEET = - EPST(I,J) + EPST(I,J+1)
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM
    AEY = AEXI * XIYM + AEET * ETYM
    ! 係数計算 ---------------------------------------------------------
    nu   = AMUM  / RHOM
    nu_t = AMUTM / RHOM
    c_u  = nu + nu_t
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    c_k  = nu + nu_t / SIGK
    c_e  = nu + nu_t / SIGE
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    DELV  = UX + VY
    TAUXX = c_u * (2.0 * UX - TwoThird * DELV) - TwoThird * AKM
    TAUYY = c_u * (2.0 * VY - TwoThird * DELV) - TwoThird * AKM
    TAUXY = c_u * (UY + VX)
    ! エネルギーの拡散 -------------------------------------------------
    R4 = TAUXX * UM + TAUXY * VM + c_t * TX
    S4 = TAUXY * UM + TAUYY * VM + c_t * TY
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    SH(I,J,1) = 0.0
    SH(I,J,2) = QH1 * (TAUXX * ETXM + TAUXY * ETYM)
    SH(I,J,3) = QH1 * (TAUXY * ETXM + TAUYY * ETYM)
    SH(I,J,4) = QH1 * (R4 * ETXM + S4 * ETYM)
    SH(I,J,5) = QH1 * c_k * (AKX * ETXM + AKY * ETYM)
    SH(I,J,6) = QH1 * c_e * (AEX * ETXM + AEY * ETYM)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFET
! 定義終了 *************************************************************
END SUBROUTINE Turbulence2DEvmLS
!***********************************************************************
!**** 乱流モデル : RANS, EVM, Launder-Sharma Model (1974)           ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence3DEvmLS( &
&            PELIM, RG, GAMMA, PR, PRT, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, AJA, &
&            RHO, U, V, W, T, AK, EPST, AMU, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: CMU = 0.09
  REAL, PARAMETER :: CE1 = 1.44, CE2 = 1.92
  REAL, PARAMETER :: SIGK = 1.0, SIGE = 1.3
  REAL, PARAMETER :: ZERO = 1.0E-20
  REAL, PARAMETER :: TwoThird = 0.66666667
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR, PRT
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE, KS:KE), &
  &                       XIY(IS:IE, JS:JE, KS:KE), &
  &                       XIZ(IS:IE, JS:JE, KS:KE), &
  &                       ETX(IS:IE, JS:JE, KS:KE), &
  &                       ETY(IS:IE, JS:JE, KS:KE), &
  &                       ETZ(IS:IE, JS:JE, KS:KE), &
  &                       ZEX(IS:IE, JS:JE, KS:KE), &
  &                       ZEY(IS:IE, JS:JE, KS:KE), &
  &                       ZEZ(IS:IE, JS:JE, KS:KE), &
  &                       AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE, KS:KE), &
  &                       U(IS:IE, JS:JE, KS:KE), &
  &                       V(IS:IE, JS:JE, KS:KE), &
  &                       W(IS:IE, JS:JE, KS:KE), &
  &                       T(IS:IE, JS:JE, KS:KE), &
  &                       AK(IS:IE, JS:JE, KS:KE), &
  &                       EPST(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: AMUT(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(OUT) :: DQP(IS:IE, JS:JE, KS:KE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :, :), SH(:, :, :, :), TH(:, :, :, :)
  REAL, ALLOCATABLE :: UXXI(:, :, :), UYXI(:, :, :), UZXI(:, :, :), &
  &                    VXXI(:, :, :), VYXI(:, :, :), VZXI(:, :, :), &
  &                    WXXI(:, :, :), WYXI(:, :, :), WZXI(:, :, :)
  REAL, ALLOCATABLE :: UXET(:, :, :), UYET(:, :, :), UZET(:, :, :), &
  &                    VXET(:, :, :), VYET(:, :, :), VZET(:, :, :), &
  &                    WXET(:, :, :), WYET(:, :, :), WZET(:, :, :)
  REAL, ALLOCATABLE :: UXZE(:, :, :), UYZE(:, :, :), UZZE(:, :, :), &
  &                    VXZE(:, :, :), VYZE(:, :, :), VZZE(:, :, :), &
  &                    WXZE(:, :, :), WYZE(:, :, :), WZZE(:, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 7) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(SH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(TH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(UXXI(IS: IE, JS: JE, KS: KE), UYXI(IS: IE, JS: JE, KS: KE))
  ALLOCATE(UZXI(IS: IE, JS: JE, KS: KE), VXXI(IS: IE, JS: JE, KS: KE))
  ALLOCATE(VYXI(IS: IE, JS: JE, KS: KE), VZXI(IS: IE, JS: JE, KS: KE))
  ALLOCATE(WXXI(IS: IE, JS: JE, KS: KE), WYXI(IS: IE, JS: JE, KS: KE))
  ALLOCATE(WZXI(IS: IE, JS: JE, KS: KE))
  ALLOCATE(UXET(IS: IE, JS: JE, KS: KE), UYET(IS: IE, JS: JE, KS: KE))
  ALLOCATE(UZET(IS: IE, JS: JE, KS: KE), VXET(IS: IE, JS: JE, KS: KE))
  ALLOCATE(VYET(IS: IE, JS: JE, KS: KE), VZET(IS: IE, JS: JE, KS: KE))
  ALLOCATE(WXET(IS: IE, JS: JE, KS: KE), WYET(IS: IE, JS: JE, KS: KE))
  ALLOCATE(WZET(IS: IE, JS: JE, KS: KE))
  ALLOCATE(UXZE(IS: IE, JS: JE, KS: KE), UYZE(IS: IE, JS: JE, KS: KE))
  ALLOCATE(UZZE(IS: IE, JS: JE, KS: KE), VXZE(IS: IE, JS: JE, KS: KE))
  ALLOCATE(VYZE(IS: IE, JS: JE, KS: KE), VZZE(IS: IE, JS: JE, KS: KE))
  ALLOCATE(WXZE(IS: IE, JS: JE, KS: KE), WYZE(IS: IE, JS: JE, KS: KE))
  ALLOCATE(WZZE(IS: IE, JS: JE, KS: KE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL PreVel2ndDiff
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 二階微分計算用の前準備                                        ****
!***********************************************************************
SUBROUTINE PreVel2ndDiff
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  UXXI = 0.0
  UYXI = 0.0
  UZXI = 0.0
  VXXI = 0.0
  VYXI = 0.0
  VZXI = 0.0
  WXXI = 0.0
  WYXI = 0.0
  WZXI = 0.0
  UXET = 0.0
  UYET = 0.0
  UZET = 0.0
  VXET = 0.0
  VYET = 0.0
  VZET = 0.0
  WXET = 0.0
  WYET = 0.0
  WZET = 0.0
  UXZE = 0.0
  UYZE = 0.0
  UZZE = 0.0
  VXZE = 0.0
  VYZE = 0.0
  VZZE = 0.0
  WXZE = 0.0
  WYZE = 0.0
  WZZE = 0.0
  ! 速度の二階微分の前準備 +++++++++++++++++++++++++++++++++++++++++++++
  CALL PreVel2ndDiffXI
  CALL PreVel2ndDiffET
  CALL PreVel2ndDiffZE
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiff
!***********************************************************************
!**** 二階微分計算用の前準備(XI方向)                                ****
!***********************************************************************
SUBROUTINE PreVel2ndDiffXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF(RHO(I,J,K) .GT. 0.0 .AND. RHO(I+1,J,K) .GT. 0.0) THEN
    ! XI方向の中間値 ---------------------------------------------------
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I+1,J,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I+1,J,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I+1,J,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I+1,J,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I+1,J,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I+1,J,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I+1,J,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I+1,J,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I+1,J,K))
    ! 計算空間一階微分 -------------------------------------------------
    UXI = - U(I,J,K) + U(I+1,J,K)
    UET = - 0.25 * ( U(I,J-1,K) + U(I+1,J-1,K) &
    &              - U(I,J+1,K) - U(I+1,J+1,K) )
    UZE = - 0.25 * ( U(I,J,K-1) + U(I+1,J,K-1) &
    &              - U(I,J,K+1) - U(I+1,J,K+1) )
    VXI = - V(I,J,K) + V(I+1,J,K)
    VET = - 0.25 * ( V(I,J-1,K) + V(I+1,J-1,K) &
    &              - V(I,J+1,K) - V(I+1,J+1,K) )
    VZE = - 0.25 * ( V(I,J,K-1) + V(I+1,J,K-1) &
    &              - V(I,J,K+1) - V(I+1,J,K+1) )
    WXI = - W(I,J,K) + W(I+1,J,K)
    WET = - 0.25 * ( W(I,J-1,K) + W(I+1,J-1,K) &
    &              - W(I,J+1,K) - W(I+1,J+1,K) )
    WZE = - 0.25 * ( W(I,J,K-1) + W(I+1,J,K-1) &
    &              - W(I,J,K+1) - W(I+1,J,K+1) )
    ! 物理空間一階微分 -------------------------------------------------
    UXXI(I,J,K) = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UYXI(I,J,K) = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZXI(I,J,K) = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VXXI(I,J,K) = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VYXI(I,J,K) = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZXI(I,J,K) = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WXXI(I,J,K) = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WYXI(I,J,K) = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZXI(I,J,K) = WXI * XIZM + WET * ETZM + WZE * ZEZM
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiffXI
!***********************************************************************
!**** 二階微分計算用の前準備(ET方向)                                ****
!***********************************************************************
SUBROUTINE PreVel2ndDiffET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF(RHO(I,J,K) .GT. 0.0 .AND. RHO(I,J+1,K) .GT. 0.0) THEN
    ! ET方向の中間値 ---------------------------------------------------
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J+1,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J+1,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J+1,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J+1,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J+1,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J+1,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J+1,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J+1,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J+1,K))
    ! 計算空間一階微分 -------------------------------------------------
    UXI = - 0.25 * ( U(I-1,J  ,K) - U(I+1,J  ,K) &
    &              + U(I-1,J+1,K) - U(I+1,J+1,K) )
    UET = - U(I,J,K) + U(I,J+1,K)
    UZE = - 0.25 * ( U(I,J,K-1) + U(I,J+1,K-1) &
    &              - U(I,J,K+1) - U(I,J+1,K+1) )
    VXI = - 0.25 * ( V(I-1,J  ,K) - V(I+1,J  ,K) &
    &              + V(I-1,J+1,K) - V(I+1,J+1,K) )
    VET = - V(I,J,K) + V(I,J+1,K)
    VZE = - 0.25 * ( V(I,J,K-1) + V(I,J+1,K-1) &
    &              - V(I,J,K+1) - V(I,J+1,K+1) )
    WXI = - 0.25 * ( W(I-1,J  ,K) - W(I+1,J  ,K) &
    &              + W(I-1,J+1,K) - W(I+1,J+1,K) )
    WET = - W(I,J,K) + W(I,J+1,K)
    WZE = - 0.25 * ( W(I,J,K-1) + W(I,J+1,K-1) &
    &              - W(I,J,K+1) - W(I,J+1,K+1) )
    ! 物理空間一階微分 -------------------------------------------------
    UXET(I,J,K) = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UYET(I,J,K) = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZET(I,J,K) = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VXET(I,J,K) = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VYET(I,J,K) = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZET(I,J,K) = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WXET(I,J,K) = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WYET(I,J,K) = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZET(I,J,K) = WXI * XIZM + WET * ETZM + WZE * ZEZM
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiffET
!***********************************************************************
!**** 二階微分計算用の前準備(ZE方向)                                ****
!***********************************************************************
SUBROUTINE PreVel2ndDiffZE
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE &
  !$)
  DO K = KS    , KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF(RHO(I,J,K) .GT. 0.0 .AND. RHO(I,J,K+1) .GT. 0.0) THEN
    ! ZE方向の中間値 ---------------------------------------------------
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J,K+1))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J,K+1))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J,K+1))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J,K+1))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J,K+1))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J,K+1))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J,K+1))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J,K+1))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J,K+1))
    ! 計算空間一階微分 -------------------------------------------------
    UXI = - 0.25 * ( U(I-1,J,K  ) - U(I+1,J,K  ) &
    &              + U(I-1,J,K+1) - U(I+1,J,K+1) )
    UET = - 0.25 * ( U(I,J-1,K  ) - U(I,J+1,K  ) &
    &              + U(I,J-1,K+1) - U(I,J+1,K+1) )
    UZE = - U(I,J,K) + U(I,J,K+1)
    VXI = - 0.25 * ( V(I-1,J,K  ) - V(I+1,J,K  ) &
    &              + V(I-1,J,K+1) - V(I+1,J,K+1) )
    VET = - 0.25 * ( V(I,J-1,K  ) - V(I,J+1,K  ) &
    &              + V(I,J-1,K+1) - V(I,J+1,K+1) )
    VZE = - V(I,J,K) + V(I,J,K+1)
    WXI = - 0.25 * ( W(I-1,J,K  ) - W(I+1,J,K  ) &
    &              + W(I-1,J,K+1) - W(I+1,J,K+1) )
    WET = - 0.25 * ( W(I,J-1,K  ) - W(I,J+1,K  ) &
    &              + W(I,J-1,K+1) - W(I,J+1,K+1) )
    WZE = - W(I,J,K) + W(I,J,K+1)
    ! 物理空間一階微分 -------------------------------------------------
    UXZE(I,J,K) = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UYZE(I,J,K) = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZZE(I,J,K) = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VXZE(I,J,K) = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VYZE(I,J,K) = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZZE(I,J,K) = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WXZE(I,J,K) = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WYZE(I,J,K) = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZZE(I,J,K) = WXI * XIZM + WET * ETZM + WZE * ZEZM
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiffZE
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  REAL    :: UX, VX, WX, UY, VY, WY, UZ, VZ, WZ, DELV
  REAL    :: KRXI, KRET, KRZE, KRX, KRY, KRZ
  REAL    :: UXX, UXY, UXZ, UYX, UYY, UYZ, UZX, UZY, UZZ, &
  &          VXX, VXY, VXZ, VYX, VYY, VYZ, VZX, VZY, VZZ, &
  &          WXX, WXY, WXZ, WYX, WYY, WYZ, WZX, WZY, WZZ
  REAL    :: Rtt, f_mu, f_2
  REAL    :: nu, nu_t
  REAL    :: SS2, PP_k, D, E
  REAL    :: QH1, AEAK
  ! 処理開始 ***********************************************************
  ! 渦粘性係数の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, nu, Rtt, f_mu, nu_t &
  !$)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J,K) .GT. 0.0) THEN
    IF(EPST(I,J,K) .GT. ZERO) THEN
      nu          = AMU(I,J,K) / RHO(I,J,K)
      Rtt         = AK(I,J,K)**2 / (nu * EPST(I,J,K))
      f_mu        = EXP(-3.4 / (1.0 + Rtt / 50.0)**2)
      nu_t        = CMU * f_mu * AK(I,J,K)**2 / EPST(I,J,K)
      AMUT(I,J,K) = RHO(I,J,K) * nu_t
    ELSE
      AMUT(I,J,K) = 0.0
    ENDIF
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 非拡散項の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE, &
  !$  UX, VX, WX, UY, VY, WY, UZ, VZ, WZ, DELV, &
  !$  KRXI, KRET, KRZE, KRX, KRY, KRZ, &
  !$  UXX, UXY, UXZ, UYX, UYY, UYZ, UZX, UZY, UZZ, &
  !$  VXX, VXY, VXZ, VYX, VYY, VYZ, VZX, VZY, VZZ, &
  !$  WXX, WXY, WXZ, WYX, WYY, WYZ, WZX, WZY, WZZ, &
  !$  Rtt, f_2, nu, nu_t, SS2, PP_k, D, E, QH1, AEAK &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I,J-1,K) .GT. 0.0) .AND. &
  &   (RHO(I-1,J,K) .GT. 0.0) .AND. (RHO(I,J,K)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J,K) .GT. 0.0) .AND. (RHO(I,J+1,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.5 * (U(I+1,J,K) - U(I-1,J,K))
    VXI = 0.5 * (V(I+1,J,K) - V(I-1,J,K))
    WXI = 0.5 * (W(I+1,J,K) - W(I-1,J,K))
    UET = 0.5 * (U(I,J+1,K) - U(I,J-1,K))
    VET = 0.5 * (V(I,J+1,K) - V(I,J-1,K))
    WET = 0.5 * (W(I,J+1,K) - W(I,J-1,K))
    UZE = 0.5 * (U(I,J,K+1) - U(I,J,K-1))
    VZE = 0.5 * (V(I,J,K+1) - V(I,J,K-1))
    WZE = 0.5 * (W(I,J,K+1) - W(I,J,K-1))
    ! 物理空間方向一階微分
    UX = UXI * XIX(I,J,K) + UET * ETX(I,J,K) + UZE * ZEX(I,J,K)
    UY = UXI * XIY(I,J,K) + UET * ETY(I,J,K) + UZE * ZEY(I,J,K)
    UZ = UXI * XIZ(I,J,K) + UET * ETZ(I,J,K) + UZE * ZEZ(I,J,K)
    VX = VXI * XIX(I,J,K) + VET * ETX(I,J,K) + VZE * ZEX(I,J,K)
    VY = VXI * XIY(I,J,K) + VET * ETY(I,J,K) + VZE * ZEY(I,J,K)
    VZ = VXI * XIZ(I,J,K) + VET * ETZ(I,J,K) + VZE * ZEZ(I,J,K)
    WX = WXI * XIX(I,J,K) + WET * ETX(I,J,K) + WZE * ZEX(I,J,K)
    WY = WXI * XIY(I,J,K) + WET * ETY(I,J,K) + WZE * ZEY(I,J,K)
    WZ = WXI * XIZ(I,J,K) + WET * ETZ(I,J,K) + WZE * ZEZ(I,J,K)
    ! 乱れの速度スケール sqrt(k) の一階微分の計算 ----------------------
    ! 計算空間方向一階微分
    KRXI = 0.5 * (SQRT(AK(I+1,J,K)) - SQRT(AK(I-1,J,K)))
    KRET = 0.5 * (SQRT(AK(I,J+1,K)) - SQRT(AK(I,J-1,K)))
    KRZE = 0.5 * (SQRT(AK(I,J,K+1)) - SQRT(AK(I,J,K-1)))
    ! 物理空間方向一階微分
    KRX = KRXI * XIX(I,J,K) + KRET * ETX(I,J,K) + KRZE * ZEX(I,J,K)
    KRY = KRXI * XIY(I,J,K) + KRET * ETY(I,J,K) + KRZE * ZEY(I,J,K)
    KRZ = KRXI * XIZ(I,J,K) + KRET * ETZ(I,J,K) + KRZE * ZEZ(I,J,K)
    ! 速度の二階微分の計算 ---------------------------------------------
    UXX = - (UXXI(I-1,J,K) - UXXI(I,J,K)) * XIX(I,J,K) &
    &     - (UXET(I,J-1,K) - UXET(I,J,K)) * ETX(I,J,K) &
    &     - (UXZE(I,J,K-1) - UXZE(I,J,K)) * ZEX(I,J,K)
    UXY = - (UXXI(I-1,J,K) - UXXI(I,J,K)) * XIY(I,J,K) &
    &     - (UXET(I,J-1,K) - UXET(I,J,K)) * ETY(I,J,K) &
    &     - (UXZE(I,J,K-1) - UXZE(I,J,K)) * ZEY(I,J,K)
    UXZ = - (UXXI(I-1,J,K) - UXXI(I,J,K)) * XIZ(I,J,K) &
    &     - (UXET(I,J-1,K) - UXET(I,J,K)) * ETZ(I,J,K) &
    &     - (UXZE(I,J,K-1) - UXZE(I,J,K)) * ZEZ(I,J,K)
    UYX = - (UYXI(I-1,J,K) - UYXI(I,J,K)) * XIX(I,J,K) &
    &     - (UYET(I,J-1,K) - UYET(I,J,K)) * ETX(I,J,K) &
    &     - (UYZE(I,J,K-1) - UYZE(I,J,K)) * ZEX(I,J,K)
    UYY = - (UYXI(I-1,J,K) - UYXI(I,J,K)) * XIY(I,J,K) &
    &     - (UYET(I,J-1,K) - UYET(I,J,K)) * ETY(I,J,K) &
    &     - (UYZE(I,J,K-1) - UYZE(I,J,K)) * ZEY(I,J,K)
    UYZ = - (UYXI(I-1,J,K) - UYXI(I,J,K)) * XIZ(I,J,K) &
    &     - (UYET(I,J-1,K) - UYET(I,J,K)) * ETZ(I,J,K) &
    &     - (UYZE(I,J,K-1) - UYZE(I,J,K)) * ZEZ(I,J,K)
    UZX = - (UZXI(I-1,J,K) - UZXI(I,J,K)) * XIX(I,J,K) &
    &     - (UZET(I,J-1,K) - UZET(I,J,K)) * ETX(I,J,K) &
    &     - (UZZE(I,J,K-1) - UZZE(I,J,K)) * ZEX(I,J,K)
    UZY = - (UZXI(I-1,J,K) - UZXI(I,J,K)) * XIY(I,J,K) &
    &     - (UZET(I,J-1,K) - UZET(I,J,K)) * ETY(I,J,K) &
    &     - (UZZE(I,J,K-1) - UZZE(I,J,K)) * ZEY(I,J,K)
    UZZ = - (UZXI(I-1,J,K) - UZXI(I,J,K)) * XIZ(I,J,K) &
    &     - (UZET(I,J-1,K) - UZET(I,J,K)) * ETZ(I,J,K) &
    &     - (UZZE(I,J,K-1) - UZZE(I,J,K)) * ZEZ(I,J,K)
    VXX = - (VXXI(I-1,J,K) - VXXI(I,J,K)) * XIX(I,J,K) &
    &     - (VXET(I,J-1,K) - VXET(I,J,K)) * ETX(I,J,K) &
    &     - (VXZE(I,J,K-1) - VXZE(I,J,K)) * ZEX(I,J,K)
    VXY = - (VXXI(I-1,J,K) - VXXI(I,J,K)) * XIY(I,J,K) &
    &     - (VXET(I,J-1,K) - VXET(I,J,K)) * ETY(I,J,K) &
    &     - (VXZE(I,J,K-1) - VXZE(I,J,K)) * ZEY(I,J,K)
    VXZ = - (VXXI(I-1,J,K) - VXXI(I,J,K)) * XIZ(I,J,K) &
    &     - (VXET(I,J-1,K) - VXET(I,J,K)) * ETZ(I,J,K) &
    &     - (VXZE(I,J,K-1) - VXZE(I,J,K)) * ZEZ(I,J,K)
    VYX = - (VYXI(I-1,J,K) - VYXI(I,J,K)) * XIX(I,J,K) &
    &     - (VYET(I,J-1,K) - VYET(I,J,K)) * ETX(I,J,K) &
    &     - (VYZE(I,J,K-1) - VYZE(I,J,K)) * ZEX(I,J,K)
    VYY = - (VYXI(I-1,J,K) - VYXI(I,J,K)) * XIY(I,J,K) &
    &     - (VYET(I,J-1,K) - VYET(I,J,K)) * ETY(I,J,K) &
    &     - (VYZE(I,J,K-1) - VYZE(I,J,K)) * ZEY(I,J,K)
    VYZ = - (VYXI(I-1,J,K) - VYXI(I,J,K)) * XIZ(I,J,K) &
    &     - (VYET(I,J-1,K) - VYET(I,J,K)) * ETZ(I,J,K) &
    &     - (VYZE(I,J,K-1) - VYZE(I,J,K)) * ZEZ(I,J,K)
    VZX = - (VZXI(I-1,J,K) - VZXI(I,J,K)) * XIX(I,J,K) &
    &     - (VZET(I,J-1,K) - VZET(I,J,K)) * ETX(I,J,K) &
    &     - (VZZE(I,J,K-1) - VZZE(I,J,K)) * ZEX(I,J,K)
    VZY = - (VZXI(I-1,J,K) - VZXI(I,J,K)) * XIY(I,J,K) &
    &     - (VZET(I,J-1,K) - VZET(I,J,K)) * ETY(I,J,K) &
    &     - (VZZE(I,J,K-1) - VZZE(I,J,K)) * ZEY(I,J,K)
    VZZ = - (VZXI(I-1,J,K) - VZXI(I,J,K)) * XIZ(I,J,K) &
    &     - (VZET(I,J-1,K) - VZET(I,J,K)) * ETZ(I,J,K) &
    &     - (VZZE(I,J,K-1) - VZZE(I,J,K)) * ZEZ(I,J,K)
    WXX = - (WXXI(I-1,J,K) - WXXI(I,J,K)) * XIX(I,J,K) &
    &     - (WXET(I,J-1,K) - WXET(I,J,K)) * ETX(I,J,K) &
    &     - (WXZE(I,J,K-1) - WXZE(I,J,K)) * ZEX(I,J,K)
    WXY = - (WXXI(I-1,J,K) - WXXI(I,J,K)) * XIY(I,J,K) &
    &     - (WXET(I,J-1,K) - WXET(I,J,K)) * ETY(I,J,K) &
    &     - (WXZE(I,J,K-1) - WXZE(I,J,K)) * ZEY(I,J,K)
    WXZ = - (WXXI(I-1,J,K) - WXXI(I,J,K)) * XIZ(I,J,K) &
    &     - (WXET(I,J-1,K) - WXET(I,J,K)) * ETZ(I,J,K) &
    &     - (WXZE(I,J,K-1) - WXZE(I,J,K)) * ZEZ(I,J,K)
    WYX = - (WYXI(I-1,J,K) - WYXI(I,J,K)) * XIX(I,J,K) &
    &     - (WYET(I,J-1,K) - WYET(I,J,K)) * ETX(I,J,K) &
    &     - (WYZE(I,J,K-1) - WYZE(I,J,K)) * ZEX(I,J,K)
    WYY = - (WYXI(I-1,J,K) - WYXI(I,J,K)) * XIY(I,J,K) &
    &     - (WYET(I,J-1,K) - WYET(I,J,K)) * ETY(I,J,K) &
    &     - (WYZE(I,J,K-1) - WYZE(I,J,K)) * ZEY(I,J,K)
    WYZ = - (WYXI(I-1,J,K) - WYXI(I,J,K)) * XIZ(I,J,K) &
    &     - (WYET(I,J-1,K) - WYET(I,J,K)) * ETZ(I,J,K) &
    &     - (WYZE(I,J,K-1) - WYZE(I,J,K)) * ZEZ(I,J,K)
    WZX = - (WZXI(I-1,J,K) - WZXI(I,J,K)) * XIX(I,J,K) &
    &     - (WZET(I,J-1,K) - WZET(I,J,K)) * ETX(I,J,K) &
    &     - (WZZE(I,J,K-1) - WZZE(I,J,K)) * ZEX(I,J,K)
    WZY = - (WZXI(I-1,J,K) - WZXI(I,J,K)) * XIY(I,J,K) &
    &     - (WZET(I,J-1,K) - WZET(I,J,K)) * ETY(I,J,K) &
    &     - (WZZE(I,J,K-1) - WZZE(I,J,K)) * ZEY(I,J,K)
    WZZ = - (WZXI(I-1,J,K) - WZXI(I,J,K)) * XIZ(I,J,K) &
    &     - (WZET(I,J-1,K) - WZET(I,J,K)) * ETZ(I,J,K) &
    &     - (WZZE(I,J,K-1) - WZZE(I,J,K)) * ZEZ(I,J,K)
    ! モデル関数の計算 -------------------------------------------------
    nu   = AMU(I,J,K)  / RHO(I,J,K)
    nu_t = AMUT(I,J,K) / RHO(I,J,K)
    IF(EPST(I,J,K) .GT. ZERO) THEN
      Rtt = AK(I,J,K)**2 / (nu * EPST(I,J,K))
    ELSE
      Rtt = 0.0
    ENDIF
    f_2  = 1.0 - 0.3 * EXP(-Rtt**2)
    ! 補正項 D, E の計算 -----------------------------------------------
    D = 2.0 * nu * (KRX**2 + KRY**2 + KRZ**2)
    E = 2.0 * nu * nu_t * ( &
    &     UXX**2 + VXX**2 + WXX**2 &
    &   + UYY**2 + VYY**2 + WYY**2 &
    &   + UZZ**2 + VZZ**2 + WZZ**2 &
    &   + UXY**2 + VXY**2 + WXY**2 &
    &   + UYZ**2 + VYZ**2 + WYZ**2 &
    &   + UZX**2 + VZX**2 + WZX**2 &
    &   + UXZ**2 + VXZ**2 + WXZ**2 &
    &   + UYX**2 + VYX**2 + WYX**2 &
    &   + UZY**2 + VZY**2 + WZY**2 &
    & )
    ! レイノルズ応力の生成(平均歪み増加率) -----------------------------
    DELV = UX + VY + WZ
    SS2  = 2.0 * (UX**2 + VY**2 + WZ**2) &
    &    + (UY + VX)**2 + (VZ + WY)**2 + (WX + UZ)**2 &
    &    - TwoThird * DELV**2
    PP_k = nu_t * SS2 - TwoThird * AK(I,J,K) * DELV
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      PP_k = MAX(0.0, MIN((EPST(I,J,K) + D) / PELIM, PP_k))
    ENDIF
    ! 生成、散逸の和 ---------------------------------------------------
    QH1 = RHO(I,J,K) / AJA(I,J,K)
    AEAK = MIN( &
    &      MAX(SQRT(SS2), SQRT((EPST(I,J,K) + D) / nu)), &
    &      EPST(I,J,K) / MAX(ZERO, AK(I,J,K)) &
    &    )
    DQP(I,J,K,6) = QH1 * (PP_k - (EPST(I,J,K) + D))
    DQP(I,J,K,7) = QH1 * ( &
    &                AEAK * (CE1 * PP_k - CE2 * f_2 * EPST(I,J,K)) + E &
    &            )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH = 0.0
  SH = 0.0
  TH = 0.0
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL DIFFXI
  CALL DIFFET
  CALL DIFFZE
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I,J-1,K) .GT. 0.0) .AND. &
  &   (RHO(I-1,J,K) .GT. 0.0) .AND. (RHO(I,J,K)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J,K) .GT. 0.0) .AND. (RHO(I,J+1,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,K,L) = (-RH(I-1,J  ,K  ,L) + RH(I,J,K,L) ) &
    &            + (-SH(I  ,J-1,K  ,L) + SH(I,J,K,L) ) &
    &            + (-TH(I  ,J  ,K-1,L) + TH(I,J,K,L) )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項の計算(xi方向差分用)                                    ****
!***********************************************************************
SUBROUTINE DIFFXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, TM, AKM, EPSTM, AMUM, AMUTM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: AKXI, AKET, AKZE, AKX, AKY, AKZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu, nu_t
  REAL    :: c_u, c_t, c_k, c_e
  REAL    :: TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV
  REAL    :: R5, S5, T5
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, TM, AKM, EPSTM, AMUM, AMUTM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  AKXI, AKET, AKZE, AKX, AKY, AKZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, nu_t, c_u, c_t, c_k, c_e, &
  !$  TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV, R5, S5, T5, QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I,J-1,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K)   .GT. 0.0) .AND. (RHO(I+1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J+1,K) .GT. 0.0) .AND. (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I+1,J,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I+1,J,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I+1,J,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I+1,J,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I+1,J,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I+1,J,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I+1,J,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I+1,J,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I+1,J,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I+1,J,K))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J,K) +  RHO(I+1,J,K))
    UM    = 0.5 * (   U(I,J,K) +    U(I+1,J,K))
    VM    = 0.5 * (   V(I,J,K) +    V(I+1,J,K))
    WM    = 0.5 * (   W(I,J,K) +    W(I+1,J,K))
    TM    = 0.5 * (   T(I,J,K) +    T(I+1,J,K))
    AKM   = 0.5 * (  AK(I,J,K) +   AK(I+1,J,K))
    EPSTM = 0.5 * (EPST(I,J,K) + EPST(I+1,J,K))
    AMUM  = 0.5 * ( AMU(I,J,K) +  AMU(I+1,J,K))
    AMUTM = 0.5 * (AMUT(I,J,K) + AMUT(I+1,J,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = - U(I,J,K) + U(I+1,J,K)
    VXI = - V(I,J,K) + V(I+1,J,K)
    WXI = - W(I,J,K) + W(I+1,J,K)
    UET = 0.25 * ( - U(I,J-1,K) - U(I+1,J-1,K) &
    &              + U(I,J+1,K) + U(I+1,J+1,K) )
    VET = 0.25 * ( - V(I,J-1,K) - V(I+1,J-1,K) &
    &              + V(I,J+1,K) + V(I+1,J+1,K) )
    WET = 0.25 * ( - W(I,J-1,K) - W(I+1,J-1,K) &
    &              + W(I,J+1,K) + W(I+1,J+1,K) )
    UZE = 0.25 * ( - U(I,J,K-1) - U(I+1,J,K-1) &
    &              + U(I,J,K+1) + U(I+1,J,K+1) )
    VZE = 0.25 * ( - V(I,J,K-1) - V(I+1,J,K-1) &
    &              + V(I,J,K+1) + V(I+1,J,K+1) )
    WZE = 0.25 * ( - W(I,J,K-1) - W(I+1,J,K-1) &
    &              + W(I,J,K+1) + W(I+1,J,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = - T(I,J,K) + T(I+1,J,K)
    TET = 0.25 * ( - T(I,J-1,K) - T(I+1,J-1,K) &
    &              + T(I,J+1,K) + T(I+1,J+1,K) )
    TZE = 0.25 * ( - T(I,J,K-1) - T(I+1,J,K-1) &
    &              + T(I,J,K+1) + T(I+1,J,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! k の一階微分の計算 -----------------------------------------------
    ! 計算空間方向一階微分
    AKXI = - AK(I,J,K) + AK(I+1,J,K)
    AKET = 0.25 * ( - AK(I,J-1,K) - AK(I+1,J-1,K) &
    &               + AK(I,J+1,K) + AK(I+1,J+1,K) )
    AKZE = 0.25 * ( - AK(I,J,K-1) - AK(I+1,J,K-1) &
    &               + AK(I,J,K+1) + AK(I+1,J,K+1) )
    ! 物理空間方向一階微分
    AKX = AKXI * XIXM + AKET * ETXM + AKZE * ZEXM
    AKY = AKXI * XIYM + AKET * ETYM + AKZE * ZEYM
    AKZ = AKXI * XIZM + AKET * ETZM + AKZE * ZEZM
    ! epsilon tilde の一階微分の計算 -----------------------------------
    ! 計算空間方向一階微分
    AEXI = - EPST(I,J,K) + EPST(I+1,J,K)
    AEET = 0.25 * ( - EPST(I,J-1,K) - EPST(I+1,J-1,K) &
    &               + EPST(I,J+1,K) + EPST(I+1,J+1,K) )
    AEZE = 0.25 * ( - EPST(I,J,K-1) - EPST(I+1,J,K-1) &
    &               + EPST(I,J,K+1) + EPST(I+1,J,K+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 係数計算 ---------------------------------------------------------
    nu   = AMUM  / RHOM
    nu_t = AMUTM / RHOM
    c_u  = nu + nu_t
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    c_k  = nu + nu_t / SIGK
    c_e  = nu + nu_t / SIGE
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    DELV  = UX + VY + WZ
    TAUXX = c_u * (2.0 * UX - TwoThird * DELV) - TwoThird * AKM
    TAUYY = c_u * (2.0 * VY - TwoThird * DELV) - TwoThird * AKM
    TAUZZ = c_u * (2.0 * WZ - TwoThird * DELV) - TwoThird * AKM
    TAUXY = c_u * (UY + VX)
    TAUYZ = c_u * (VZ + WY)
    TAUZX = c_u * (WX + UZ)
    ! エネルギーの拡散 -------------------------------------------------
    R5 = TAUXX * UM + TAUXY * VM + TAUZX * WM + c_t * TX
    S5 = TAUXY * UM + TAUYY * VM + TAUYZ * WM + c_t * TY
    T5 = TAUZX * UM + TAUYZ * VM + TAUZZ * WM + c_t * TZ
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    RH(I,J,K,1) = 0.0
    RH(I,J,K,2) = QH1 * (TAUXX * XIXM + TAUXY * XIYM + TAUZX * XIZM)
    RH(I,J,K,3) = QH1 * (TAUXY * XIXM + TAUYY * XIYM + TAUYZ * XIZM)
    RH(I,J,K,4) = QH1 * (TAUZX * XIXM + TAUYZ * XIYM + TAUZZ * XIZM)
    RH(I,J,K,5) = QH1 * (R5 * XIXM + S5 * XIYM + T5 * XIZM)
    RH(I,J,K,6) = QH1 * c_k * (AKX * XIXM + AKY * XIYM + AKZ * XIZM)
    RH(I,J,K,7) = QH1 * c_e * (AEX * XIXM + AEY * XIYM + AEZ * XIZM)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFXI
!***********************************************************************
!**** 拡散項の計算(eta方向差分用)                                   ****
!***********************************************************************
SUBROUTINE DIFFET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, TM, AKM, EPSTM, AMUM, AMUTM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: AKXI, AKET, AKZE, AKX, AKY, AKZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu, nu_t
  REAL    :: c_u, c_t, c_k, c_e
  REAL    :: TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV
  REAL    :: R5, S5, T5
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, TM, AKM, EPSTM, AMUM, AMUTM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  AKXI, AKET, AKZE, AKX, AKY, AKZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, nu_t, c_u, c_t, c_k, c_e, &
  !$  TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV, R5, S5, T5, QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I-1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K)   .GT. 0.0) .AND. (RHO(I+1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J+1,K) .GT. 0.0) .AND. (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J+1,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J+1,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J+1,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J+1,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J+1,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J+1,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J+1,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J+1,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J+1,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J+1,K))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J,K) +  RHO(I,J+1,K))
    UM    = 0.5 * (   U(I,J,K) +    U(I,J+1,K))
    VM    = 0.5 * (   V(I,J,K) +    V(I,J+1,K))
    WM    = 0.5 * (   W(I,J,K) +    W(I,J+1,K))
    TM    = 0.5 * (   T(I,J,K) +    T(I,J+1,K))
    AKM   = 0.5 * (  AK(I,J,K) +   AK(I,J+1,K))
    EPSTM = 0.5 * (EPST(I,J,K) + EPST(I,J+1,K))
    AMUM  = 0.5 * ( AMU(I,J,K) +  AMU(I,J+1,K))
    AMUTM = 0.5 * (AMUT(I,J,K) + AMUT(I,J+1,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J+1,K) &
    &            + U(I+1,J,K) + U(I+1,J+1,K) )
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J+1,K) &
    &            + V(I+1,J,K) + V(I+1,J+1,K) )
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J+1,K) &
    &            + W(I+1,J,K) + W(I+1,J+1,K) )
    UET = - U(I,J,K) + U(I,J+1,K)
    VET = - V(I,J,K) + V(I,J+1,K)
    WET = - W(I,J,K) + W(I,J+1,K)
    UZE = 0.25 * ( - U(I,J,K-1) - U(I,J+1,K-1) &
    &              + U(I,J,K+1) + U(I,J+1,K+1) )
    VZE = 0.25 * ( - V(I,J,K-1) - V(I,J+1,K-1) &
    &              + V(I,J,K+1) + V(I,J+1,K+1) )
    WZE = 0.25 * ( - W(I,J,K-1) - W(I,J+1,K-1) &
    &              + W(I,J,K+1) + W(I,J+1,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J+1,K) &
    &              + T(I+1,J,K) + T(I+1,J+1,K) )
    TET = - T(I,J,K) + T(I,J+1,K)
    TZE = 0.25 * ( - T(I,J,K-1) - T(I,J+1,K-1) &
    &              + T(I,J,K+1) + T(I,J+1,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! k の一階微分の計算 -----------------------------------------------
    ! 計算空間方向一階微分
    AKXI = 0.25 * ( - AK(I-1,J,K) - AK(I-1,J+1,K) &
    &               + AK(I+1,J,K) + AK(I+1,J+1,K) )
    AKET = - AK(I,J,K) + AK(I,J+1,K)
    AKZE = 0.25 * ( - AK(I,J,K-1) - AK(I,J+1,K-1) &
    &               + AK(I,J,K+1) + AK(I,J+1,K+1) )
    ! 物理空間方向一階微分
    AKX = AKXI * XIXM + AKET * ETXM + AKZE * ZEXM
    AKY = AKXI * XIYM + AKET * ETYM + AKZE * ZEYM
    AKZ = AKXI * XIZM + AKET * ETZM + AKZE * ZEZM
    ! epsilon tilde の一階微分の計算 -----------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPST(I-1,J,K) - EPST(I-1,J+1,K) &
    &               + EPST(I+1,J,K) + EPST(I+1,J+1,K) )
    AEET = - EPST(I,J,K) + EPST(I,J+1,K)
    AEZE = 0.25 * ( - EPST(I,J,K-1) - EPST(I,J+1,K-1) &
    &               + EPST(I,J,K+1) + EPST(I,J+1,K+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 係数計算 ---------------------------------------------------------
    nu   = AMUM  / RHOM
    nu_t = AMUTM / RHOM
    c_u  = nu + nu_t
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    c_k  = nu + nu_t / SIGK
    c_e  = nu + nu_t / SIGE
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    DELV  = UX + VY + WZ
    TAUXX = c_u * (2.0 * UX - TwoThird * DELV) - TwoThird * AKM
    TAUYY = c_u * (2.0 * VY - TwoThird * DELV) - TwoThird * AKM
    TAUZZ = c_u * (2.0 * WZ - TwoThird * DELV) - TwoThird * AKM
    TAUXY = c_u * (UY + VX)
    TAUYZ = c_u * (VZ + WY)
    TAUZX = c_u * (WX + UZ)
    ! エネルギーの拡散 -------------------------------------------------
    R5 = TAUXX * UM + TAUXY * VM + TAUZX * WM + c_t * TX
    S5 = TAUXY * UM + TAUYY * VM + TAUYZ * WM + c_t * TY
    T5 = TAUZX * UM + TAUYZ * VM + TAUZZ * WM + c_t * TZ
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    SH(I,J,K,1) = 0.0
    SH(I,J,K,2) = QH1 * (TAUXX * ETXM + TAUXY * ETYM + TAUZX * ETZM)
    SH(I,J,K,3) = QH1 * (TAUXY * ETXM + TAUYY * ETYM + TAUYZ * ETZM)
    SH(I,J,K,4) = QH1 * (TAUZX * ETXM + TAUYZ * ETYM + TAUZZ * ETZM)
    SH(I,J,K,5) = QH1 * (R5 * ETXM + S5 * ETYM + T5 * ETZM)
    SH(I,J,K,6) = QH1 * c_k * (AKX * ETXM + AKY * ETYM + AKZ * ETZM)
    SH(I,J,K,7) = QH1 * c_e * (AEX * ETXM + AEY * ETYM + AEZ * ETZM)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFET
!***********************************************************************
!**** 拡散項の計算(zeta方向差分用)                                  ****
!***********************************************************************
SUBROUTINE DIFFZE
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, TM, AKM, EPSTM, AMUM, AMUTM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: AKXI, AKET, AKZE, AKX, AKY, AKZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu, nu_t
  REAL    :: c_u, c_t, c_k, c_e
  REAL    :: TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV
  REAL    :: R5, S5, T5
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, TM, AKM, EPSTM, AMUM, AMUTM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  AKXI, AKET, AKZE, AKX, AKY, AKZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, nu_t, c_u, c_t, c_k, c_e, &
  !$  TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV, R5, S5, T5, QH1 &
  !$)
  DO K = KS    , KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J-1,K) .GT. 0.0) .AND. (RHO(I-1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K)   .GT. 0.0) .AND. (RHO(I+1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J+1,K) .GT. 0.0) .AND. (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    ! ZE方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J,K+1))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J,K+1))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J,K+1))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J,K+1))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J,K+1))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J,K+1))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J,K+1))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J,K+1))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J,K+1))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J,K+1))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J,K) +  RHO(I,J,K+1))
    UM    = 0.5 * (   U(I,J,K) +    U(I,J,K+1))
    VM    = 0.5 * (   V(I,J,K) +    V(I,J,K+1))
    WM    = 0.5 * (   W(I,J,K) +    W(I,J,K+1))
    TM    = 0.5 * (   T(I,J,K) +    T(I,J,K+1))
    AKM   = 0.5 * (  AK(I,J,K) +   AK(I,J,K+1))
    EPSTM = 0.5 * (EPST(I,J,K) + EPST(I,J,K+1))
    AMUM  = 0.5 * ( AMU(I,J,K) +  AMU(I,J,K+1))
    AMUTM = 0.5 * (AMUT(I,J,K) + AMUT(I,J,K+1))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J,K+1) &
    &              + U(I+1,J,K) + U(I+1,J,K+1) )
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J,K+1) &
    &              + V(I+1,J,K) + V(I+1,J,K+1) )
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J,K+1) &
    &              + W(I+1,J,K) + W(I+1,J,K+1) )
    UET = 0.25 * ( - U(I,J-1,K) - U(I,J-1,K+1) &
    &              + U(I,J+1,K) + U(I,J+1,K+1) )
    VET = 0.25 * ( - V(I,J-1,K) - V(I,J-1,K+1) &
    &              + V(I,J+1,K) + V(I,J+1,K+1) )
    WET = 0.25 * ( - W(I,J-1,K) - W(I,J-1,K+1) &
    &              + W(I,J+1,K) + W(I,J+1,K+1) )
    UZE = - U(I,J,K) + U(I,J,K+1)
    VZE = - V(I,J,K) + V(I,J,K+1)
    WZE = - W(I,J,K) + W(I,J,K+1)
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J,K+1) &
    &              + T(I+1,J,K) + T(I+1,J,K+1) )
    TET = 0.25 * ( - T(I,J-1,K) - T(I,J-1,K+1) &
    &              + T(I,J+1,K) + T(I,J+1,K+1) )
    TZE = - T(I,J,K) + T(I,J,K+1)
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! k の一階微分の計算 -----------------------------------------------
    ! 計算空間方向一階微分
    AKXI = 0.25 * ( - AK(I-1,J,K) - AK(I-1,J,K+1) &
    &               + AK(I+1,J,K) + AK(I+1,J,K+1) )
    AKET = 0.25 * ( - AK(I,J-1,K) - AK(I,J-1,K+1) &
    &               + AK(I,J+1,K) + AK(I,J+1,K+1) )
    AKZE = - AK(I,J,K) + AK(I,J,K+1)
    ! 物理空間方向一階微分
    AKX = AKXI * XIXM + AKET * ETXM + AKZE * ZEXM
    AKY = AKXI * XIYM + AKET * ETYM + AKZE * ZEYM
    AKZ = AKXI * XIZM + AKET * ETZM + AKZE * ZEZM
    ! epsilon tilde の一階微分の計算 -----------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPST(I-1,J,K) - EPST(I-1,J,K+1) &
    &               + EPST(I+1,J,K) + EPST(I+1,J,K+1) )
    AEET = 0.25 * ( - EPST(I,J-1,K) - EPST(I,J-1,K+1) &
    &               + EPST(I,J+1,K) + EPST(I,J+1,K+1) )
    AEZE = - EPST(I,J,K) + EPST(I,J,K+1)
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 係数計算 ---------------------------------------------------------
    nu   = AMUM  / RHOM
    nu_t = AMUTM / RHOM
    c_u  = nu + nu_t
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    c_k  = nu + nu_t / SIGK
    c_e  = nu + nu_t / SIGE
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    DELV  = UX + VY + WZ
    TAUXX = c_u * (2.0 * UX - TwoThird * DELV) - TwoThird * AKM
    TAUYY = c_u * (2.0 * VY - TwoThird * DELV) - TwoThird * AKM
    TAUZZ = c_u * (2.0 * WZ - TwoThird * DELV) - TwoThird * AKM
    TAUXY = c_u * (UY + VX)
    TAUYZ = c_u * (VZ + WY)
    TAUZX = c_u * (WX + UZ)
    ! エネルギーの拡散 -------------------------------------------------
    R5 = TAUXX * UM + TAUXY * VM + TAUZX * WM + c_t * TX
    S5 = TAUXY * UM + TAUYY * VM + TAUYZ * WM + c_t * TY
    T5 = TAUZX * UM + TAUYZ * VM + TAUZZ * WM + c_t * TZ
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    TH(I,J,K,1) = 0.0
    TH(I,J,K,2) = QH1 * (TAUXX * ZEXM + TAUXY * ZEYM + TAUZX * ZEZM)
    TH(I,J,K,3) = QH1 * (TAUXY * ZEXM + TAUYY * ZEYM + TAUYZ * ZEZM)
    TH(I,J,K,4) = QH1 * (TAUZX * ZEXM + TAUYZ * ZEYM + TAUZZ * ZEZM)
    TH(I,J,K,5) = QH1 * (R5 * ZEXM + S5 * ZEYM + T5 * ZEZM)
    TH(I,J,K,6) = QH1 * c_k * (AKX * ZEXM + AKY * ZEYM + AKZ * ZEZM)
    TH(I,J,K,7) = QH1 * c_e * (AEX * ZEXM + AEY * ZEYM + AEZ * ZEZM)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFZE
! 定義終了 *************************************************************
END SUBROUTINE Turbulence3DEvmLS
!***********************************************************************
!**** 乱流モデル : RANS, EVM, Abe-Kondo-Nagano Model (1994)         ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence3DEvmAKN( &
&            PELIM, RG, GAMMA, PR, PRT, &
&            IS, IE, I1, I3, JS, JE, KS, KE, LS, LE, &
&            X, Y, Z, &
&            XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, AJA, &
&            RHO, U, V, W, T, AK, EPST, AMU, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: CMU = 0.09
  REAL, PARAMETER :: CE1 = 1.5, CE2 = 1.9
  REAL, PARAMETER :: SIGK = 1.4, SIGE = 1.4
  REAL, PARAMETER :: ZERO = 1.0E-20
  REAL, PARAMETER :: TwoThird = 0.66666667
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR, PRT
  INTEGER, INTENT(IN)  :: IS, IE, I1, I3, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: X(IS:IE, JS:JE, KS:KE), &
  &                       Y(IS:IE, JS:JE, KS:KE), &
  &                       Z(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE, KS:KE), &
  &                       XIY(IS:IE, JS:JE, KS:KE), &
  &                       XIZ(IS:IE, JS:JE, KS:KE), &
  &                       ETX(IS:IE, JS:JE, KS:KE), &
  &                       ETY(IS:IE, JS:JE, KS:KE), &
  &                       ETZ(IS:IE, JS:JE, KS:KE), &
  &                       ZEX(IS:IE, JS:JE, KS:KE), &
  &                       ZEY(IS:IE, JS:JE, KS:KE), &
  &                       ZEZ(IS:IE, JS:JE, KS:KE), &
  &                       AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE, KS:KE), &
  &                       U(IS:IE, JS:JE, KS:KE), &
  &                       V(IS:IE, JS:JE, KS:KE), &
  &                       W(IS:IE, JS:JE, KS:KE), &
  &                       T(IS:IE, JS:JE, KS:KE), &
  &                       AK(IS:IE, JS:JE, KS:KE), &
  &                       EPST(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: AMUT(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(OUT) :: DQP(IS:IE, JS:JE, KS:KE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :, :), SH(:, :, :, :), TH(:, :, :, :)
  REAL, ALLOCATABLE :: UXXI(:, :, :), UYXI(:, :, :), UZXI(:, :, :), &
  &                    VXXI(:, :, :), VYXI(:, :, :), VZXI(:, :, :), &
  &                    WXXI(:, :, :), WYXI(:, :, :), WZXI(:, :, :)
  REAL, ALLOCATABLE :: UXET(:, :, :), UYET(:, :, :), UZET(:, :, :), &
  &                    VXET(:, :, :), VYET(:, :, :), VZET(:, :, :), &
  &                    WXET(:, :, :), WYET(:, :, :), WZET(:, :, :)
  REAL, ALLOCATABLE :: UXZE(:, :, :), UYZE(:, :, :), UZZE(:, :, :), &
  &                    VXZE(:, :, :), VYZE(:, :, :), VZZE(:, :, :), &
  &                    WXZE(:, :, :), WYZE(:, :, :), WZZE(:, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 7) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(SH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(TH(IS: IE, JS: JE, KS: KE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL NonDiffTerm
  CALL DiffTerm
  ! メモリ解放 =========================================================
  DEALLOCATE(RH,SH,TH)
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  REAL    :: UX, VX, WX, UY, VY, WY, UZ, VZ, WZ, DELV
  DOUBLE PRECISION :: ETA, YST
  REAL    :: f_mu, f_1
  REAL,DIMENSION(IS:IE,JS:JE,KS:KE) :: Rtt, f_2
  REAL    :: nu, nu_t
  REAL    :: SS2, PP_k, D, E
  REAL    :: QH1, AEAK
  ! 処理開始 ***********************************************************
  ! 渦粘性係数,補正項の計算 ++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, nu, f_mu, nu_t &
  !$)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J,K) .GT. 0.0) THEN
    IF(EPST(I,J,K) .GT. ZERO .AND. AK(I,J,K) .GT. ZERO) THEN
      ETA = (DBLE(AMU(I,J,K) / RHO(I,J,K))**3.0 / &
      &     DBLE(EPST(I,J,K)))**0.25
      IF(I .GE. I1 .AND. I .LT. I3) THEN
        YST = DSQRT(DBLE(X(I,J,K) - X(I,JS,K))**2.0 + &
        &           DBLE(Y(I,J,K) - Y(I,JS,K))**2.0 + &
        &           DBLE(Z(I,J,K) - Z(I,JS,K))**2.0)
      ELSE
        YST = DSQRT(DBLE(X(I,J,K) - X(I1,JS,K))**2.0 + &
        &           DBLE(Y(I,J,K) - Y(I1,JS,K))**2.0 + &
        &           DBLE(Z(I,J,K) - Z(I1,JS,K))**2.0)
      END IF
      YST = YST / ETA
      nu          = AMU(I,J,K) / RHO(I,J,K)
      Rtt(I,J,K)  = REAL(DBLE(AK(I,J,K))**2.0 / DBLE(nu * EPST(I,J,K)))
      f_mu        = REAL(1.0D0 - DEXP(-YST / 14.0D0))**2.0 * &
      &             (1.0 + 5.0 / Rtt(I,J,K)**0.75 * EXP(-(Rtt(I,J,K) / 200.0)**2.0))
      f_2(I,J,K)  = REAL(1.0D0 - DEXP(-YST / 3.1D0))**2.0 * &
      &             (1.0 - 0.3 * EXP(-(Rtt(I,J,K) / 6.5)**2.0))
      nu_t        = CMU * f_mu * AK(I,J,K)**2 / EPST(I,J,K)
      AMUT(I,J,K) = RHO(I,J,K) * nu_t
    ELSE
      Rtt(I,J,K)  = 0.0
      f_2(I,J,K)  = 1.0
      AMUT(I,J,K) = 0.0
    ENDIF
    IF(Rtt(I,J,K) .LE. ZERO) then
      Rtt(I,J,K)  = 0.0
      f_2(I,J,K)  = 1.0
      AMUT(I,J,K) = 0.0
    ENDIF
!if(k .eq. ks .and. j .eq. js) write(*,*) i,j,'fmu,f2',f_mu,f_2(I,J,K),Rtt(I,J,K),nu,AMU(I,J,K),RHO(I,J,K)
!if(k .eq. ks .and. j .eq. js+1) write(*,*) i,j,'fmu,f2',f_mu,f_2(I,J,K),Rtt(I,J,K),nu,AMU(I,J,K),RHO(I,J,K)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 非拡散項の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE, &
  !$  UX, VX, WX, UY, VY, WY, UZ, VZ, WZ, DELV, &
  !$  SS2, PP_k, D, E, f_1, QH1, AEAK &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I,J-1,K) .GT. 0.0) .AND. &
  &   (RHO(I-1,J,K) .GT. 0.0) .AND. (RHO(I,J,K)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J,K) .GT. 0.0) .AND. (RHO(I,J+1,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.5 * (U(I+1,J,K) - U(I-1,J,K))
    VXI = 0.5 * (V(I+1,J,K) - V(I-1,J,K))
    WXI = 0.5 * (W(I+1,J,K) - W(I-1,J,K))
    UET = 0.5 * (U(I,J+1,K) - U(I,J-1,K))
    VET = 0.5 * (V(I,J+1,K) - V(I,J-1,K))
    WET = 0.5 * (W(I,J+1,K) - W(I,J-1,K))
    UZE = 0.5 * (U(I,J,K+1) - U(I,J,K-1))
    VZE = 0.5 * (V(I,J,K+1) - V(I,J,K-1))
    WZE = 0.5 * (W(I,J,K+1) - W(I,J,K-1))
    ! 物理空間方向一階微分
    UX = UXI * XIX(I,J,K) + UET * ETX(I,J,K) + UZE * ZEX(I,J,K)
    UY = UXI * XIY(I,J,K) + UET * ETY(I,J,K) + UZE * ZEY(I,J,K)
    UZ = UXI * XIZ(I,J,K) + UET * ETZ(I,J,K) + UZE * ZEZ(I,J,K)
    VX = VXI * XIX(I,J,K) + VET * ETX(I,J,K) + VZE * ZEX(I,J,K)
    VY = VXI * XIY(I,J,K) + VET * ETY(I,J,K) + VZE * ZEY(I,J,K)
    VZ = VXI * XIZ(I,J,K) + VET * ETZ(I,J,K) + VZE * ZEZ(I,J,K)
    WX = WXI * XIX(I,J,K) + WET * ETX(I,J,K) + WZE * ZEX(I,J,K)
    WY = WXI * XIY(I,J,K) + WET * ETY(I,J,K) + WZE * ZEY(I,J,K)
    WZ = WXI * XIZ(I,J,K) + WET * ETZ(I,J,K) + WZE * ZEZ(I,J,K)
    ! 補正項 D, E の計算 -----------------------------------------------
    D = 0.0
    E = 0.0
    ! 補正項 f_1の計算 -------------------------------------------------
    f_1 = 1.0
    ! レイノルズ応力の生成(平均歪み増加率) -----------------------------
    DELV = UX + VY + WZ
    SS2  = 2.0 * (UX**2 + VY**2 + WZ**2) &
    &    + (UY + VX)**2 + (VZ + WY)**2 + (WX + UZ)**2 &
    &    - TwoThird * DELV**2
    PP_k = AMUT(I,J,K) / RHO(I,J,K) * SS2 - TwoThird * AK(I,J,K) * DELV
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      PP_k = MAX(0.0, MIN((EPST(I,J,K) + D) / PELIM, PP_k))
    ENDIF
    ! 生成、散逸の和 ---------------------------------------------------
    QH1 = RHO(I,J,K) / AJA(I,J,K)
if(QH1 .ne. QH1) then
 write(*,*) 'QH1 is diverged'
 write(*,*) i,j,k,RHO(I,J,K),AJA(I,J,K)
 stop
end if
    AEAK = MIN( &
    &      MAX(SQRT(SS2), SQRT((EPST(I,J,K) + D) / (AMU(I,J,K) / RHO(I,J,K)))), &
    &      EPST(I,J,K) / MAX(ZERO, AK(I,J,K)) &
    &    )
    DQP(I,J,K,6) = QH1 * (PP_k - (EPST(I,J,K) + D))
    DQP(I,J,K,7) = QH1 * ( &
    &                AEAK * (CE1 * f_1 * PP_k - CE2 * f_2(I,J,K) * EPST(I,J,K)) + E &
    &            )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH = 0.0
  SH = 0.0
  TH = 0.0
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL DIFFXI
  CALL DIFFET
  CALL DIFFZE
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I,J-1,K) .GT. 0.0) .AND. &
  &   (RHO(I-1,J,K) .GT. 0.0) .AND. (RHO(I,J,K)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J,K) .GT. 0.0) .AND. (RHO(I,J+1,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,K,L) = (-RH(I-1,J  ,K  ,L) + RH(I,J,K,L) ) &
    &            + (-SH(I  ,J-1,K  ,L) + SH(I,J,K,L) ) &
    &            + (-TH(I  ,J  ,K-1,L) + TH(I,J,K,L) )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項の計算(xi方向差分用)                                    ****
!***********************************************************************
SUBROUTINE DIFFXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, TM, AKM, EPSTM, AMUM, AMUTM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: AKXI, AKET, AKZE, AKX, AKY, AKZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu, nu_t
  REAL    :: c_u, c_t, c_k, c_e
  REAL    :: TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV
  REAL    :: R5, S5, T5
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, TM, AKM, EPSTM, AMUM, AMUTM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  AKXI, AKET, AKZE, AKX, AKY, AKZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, nu_t, c_u, c_t, c_k, c_e, &
  !$  TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV, R5, S5, T5, QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I,J-1,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K)   .GT. 0.0) .AND. (RHO(I+1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J+1,K) .GT. 0.0) .AND. (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I+1,J,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I+1,J,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I+1,J,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I+1,J,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I+1,J,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I+1,J,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I+1,J,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I+1,J,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I+1,J,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I+1,J,K))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J,K) +  RHO(I+1,J,K))
    UM    = 0.5 * (   U(I,J,K) +    U(I+1,J,K))
    VM    = 0.5 * (   V(I,J,K) +    V(I+1,J,K))
    WM    = 0.5 * (   W(I,J,K) +    W(I+1,J,K))
    TM    = 0.5 * (   T(I,J,K) +    T(I+1,J,K))
    AKM   = 0.5 * (  AK(I,J,K) +   AK(I+1,J,K))
    EPSTM = 0.5 * (EPST(I,J,K) + EPST(I+1,J,K))
    AMUM  = 0.5 * ( AMU(I,J,K) +  AMU(I+1,J,K))
    AMUTM = 0.5 * (AMUT(I,J,K) + AMUT(I+1,J,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = - U(I,J,K) + U(I+1,J,K)
    VXI = - V(I,J,K) + V(I+1,J,K)
    WXI = - W(I,J,K) + W(I+1,J,K)
    UET = 0.25 * ( - U(I,J-1,K) - U(I+1,J-1,K) &
    &              + U(I,J+1,K) + U(I+1,J+1,K) )
    VET = 0.25 * ( - V(I,J-1,K) - V(I+1,J-1,K) &
    &              + V(I,J+1,K) + V(I+1,J+1,K) )
    WET = 0.25 * ( - W(I,J-1,K) - W(I+1,J-1,K) &
    &              + W(I,J+1,K) + W(I+1,J+1,K) )
    UZE = 0.25 * ( - U(I,J,K-1) - U(I+1,J,K-1) &
    &              + U(I,J,K+1) + U(I+1,J,K+1) )
    VZE = 0.25 * ( - V(I,J,K-1) - V(I+1,J,K-1) &
    &              + V(I,J,K+1) + V(I+1,J,K+1) )
    WZE = 0.25 * ( - W(I,J,K-1) - W(I+1,J,K-1) &
    &              + W(I,J,K+1) + W(I+1,J,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = - T(I,J,K) + T(I+1,J,K)
    TET = 0.25 * ( - T(I,J-1,K) - T(I+1,J-1,K) &
    &              + T(I,J+1,K) + T(I+1,J+1,K) )
    TZE = 0.25 * ( - T(I,J,K-1) - T(I+1,J,K-1) &
    &              + T(I,J,K+1) + T(I+1,J,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! k の一階微分の計算 -----------------------------------------------
    ! 計算空間方向一階微分
    AKXI = - AK(I,J,K) + AK(I+1,J,K)
    AKET = 0.25 * ( - AK(I,J-1,K) - AK(I+1,J-1,K) &
    &               + AK(I,J+1,K) + AK(I+1,J+1,K) )
    AKZE = 0.25 * ( - AK(I,J,K-1) - AK(I+1,J,K-1) &
    &               + AK(I,J,K+1) + AK(I+1,J,K+1) )
    ! 物理空間方向一階微分
    AKX = AKXI * XIXM + AKET * ETXM + AKZE * ZEXM
    AKY = AKXI * XIYM + AKET * ETYM + AKZE * ZEYM
    AKZ = AKXI * XIZM + AKET * ETZM + AKZE * ZEZM
    ! epsilon tilde の一階微分の計算 -----------------------------------
    ! 計算空間方向一階微分
    AEXI = - EPST(I,J,K) + EPST(I+1,J,K)
    AEET = 0.25 * ( - EPST(I,J-1,K) - EPST(I+1,J-1,K) &
    &               + EPST(I,J+1,K) + EPST(I+1,J+1,K) )
    AEZE = 0.25 * ( - EPST(I,J,K-1) - EPST(I+1,J,K-1) &
    &               + EPST(I,J,K+1) + EPST(I+1,J,K+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 係数計算 ---------------------------------------------------------
    nu   = AMUM  / RHOM
    nu_t = AMUTM / RHOM
    c_u  = nu + nu_t
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    c_k  = nu + nu_t / SIGK
    c_e  = nu + nu_t / SIGE
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    DELV  = UX + VY + WZ
    TAUXX = c_u * (2.0 * UX - TwoThird * DELV) - TwoThird * AKM
    TAUYY = c_u * (2.0 * VY - TwoThird * DELV) - TwoThird * AKM
    TAUZZ = c_u * (2.0 * WZ - TwoThird * DELV) - TwoThird * AKM
    TAUXY = c_u * (UY + VX)
    TAUYZ = c_u * (VZ + WY)
    TAUZX = c_u * (WX + UZ)
    ! エネルギーの拡散 -------------------------------------------------
    R5 = TAUXX * UM + TAUXY * VM + TAUZX * WM + c_t * TX
    S5 = TAUXY * UM + TAUYY * VM + TAUYZ * WM + c_t * TY
    T5 = TAUZX * UM + TAUYZ * VM + TAUZZ * WM + c_t * TZ
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    RH(I,J,K,1) = 0.0
    RH(I,J,K,2) = QH1 * (TAUXX * XIXM + TAUXY * XIYM + TAUZX * XIZM)
    RH(I,J,K,3) = QH1 * (TAUXY * XIXM + TAUYY * XIYM + TAUYZ * XIZM)
    RH(I,J,K,4) = QH1 * (TAUZX * XIXM + TAUYZ * XIYM + TAUZZ * XIZM)
    RH(I,J,K,5) = QH1 * (R5 * XIXM + S5 * XIYM + T5 * XIZM)
    RH(I,J,K,6) = QH1 * c_k * (AKX * XIXM + AKY * XIYM + AKZ * XIZM)
    RH(I,J,K,7) = QH1 * c_e * (AEX * XIXM + AEY * XIYM + AEZ * XIZM)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFXI
!***********************************************************************
!**** 拡散項の計算(eta方向差分用)                                   ****
!***********************************************************************
SUBROUTINE DIFFET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, TM, AKM, EPSTM, AMUM, AMUTM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: AKXI, AKET, AKZE, AKX, AKY, AKZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu, nu_t
  REAL    :: c_u, c_t, c_k, c_e
  REAL    :: TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV
  REAL    :: R5, S5, T5
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, TM, AKM, EPSTM, AMUM, AMUTM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  AKXI, AKET, AKZE, AKX, AKY, AKZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, nu_t, c_u, c_t, c_k, c_e, &
  !$  TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV, R5, S5, T5, QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I-1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K)   .GT. 0.0) .AND. (RHO(I+1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J+1,K) .GT. 0.0) .AND. (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J+1,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J+1,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J+1,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J+1,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J+1,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J+1,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J+1,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J+1,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J+1,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J+1,K))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J,K) +  RHO(I,J+1,K))
    UM    = 0.5 * (   U(I,J,K) +    U(I,J+1,K))
    VM    = 0.5 * (   V(I,J,K) +    V(I,J+1,K))
    WM    = 0.5 * (   W(I,J,K) +    W(I,J+1,K))
    TM    = 0.5 * (   T(I,J,K) +    T(I,J+1,K))
    AKM   = 0.5 * (  AK(I,J,K) +   AK(I,J+1,K))
    EPSTM = 0.5 * (EPST(I,J,K) + EPST(I,J+1,K))
    AMUM  = 0.5 * ( AMU(I,J,K) +  AMU(I,J+1,K))
    AMUTM = 0.5 * (AMUT(I,J,K) + AMUT(I,J+1,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J+1,K) &
    &            + U(I+1,J,K) + U(I+1,J+1,K) )
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J+1,K) &
    &            + V(I+1,J,K) + V(I+1,J+1,K) )
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J+1,K) &
    &            + W(I+1,J,K) + W(I+1,J+1,K) )
    UET = - U(I,J,K) + U(I,J+1,K)
    VET = - V(I,J,K) + V(I,J+1,K)
    WET = - W(I,J,K) + W(I,J+1,K)
    UZE = 0.25 * ( - U(I,J,K-1) - U(I,J+1,K-1) &
    &              + U(I,J,K+1) + U(I,J+1,K+1) )
    VZE = 0.25 * ( - V(I,J,K-1) - V(I,J+1,K-1) &
    &              + V(I,J,K+1) + V(I,J+1,K+1) )
    WZE = 0.25 * ( - W(I,J,K-1) - W(I,J+1,K-1) &
    &              + W(I,J,K+1) + W(I,J+1,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J+1,K) &
    &              + T(I+1,J,K) + T(I+1,J+1,K) )
    TET = - T(I,J,K) + T(I,J+1,K)
    TZE = 0.25 * ( - T(I,J,K-1) - T(I,J+1,K-1) &
    &              + T(I,J,K+1) + T(I,J+1,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! k の一階微分の計算 -----------------------------------------------
    ! 計算空間方向一階微分
    AKXI = 0.25 * ( - AK(I-1,J,K) - AK(I-1,J+1,K) &
    &               + AK(I+1,J,K) + AK(I+1,J+1,K) )
    AKET = - AK(I,J,K) + AK(I,J+1,K)
    AKZE = 0.25 * ( - AK(I,J,K-1) - AK(I,J+1,K-1) &
    &               + AK(I,J,K+1) + AK(I,J+1,K+1) )
    ! 物理空間方向一階微分
    AKX = AKXI * XIXM + AKET * ETXM + AKZE * ZEXM
    AKY = AKXI * XIYM + AKET * ETYM + AKZE * ZEYM
    AKZ = AKXI * XIZM + AKET * ETZM + AKZE * ZEZM
    ! epsilon tilde の一階微分の計算 -----------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPST(I-1,J,K) - EPST(I-1,J+1,K) &
    &               + EPST(I+1,J,K) + EPST(I+1,J+1,K) )
    AEET = - EPST(I,J,K) + EPST(I,J+1,K)
    AEZE = 0.25 * ( - EPST(I,J,K-1) - EPST(I,J+1,K-1) &
    &               + EPST(I,J,K+1) + EPST(I,J+1,K+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 係数計算 ---------------------------------------------------------
    nu   = AMUM  / RHOM
    nu_t = AMUTM / RHOM
    c_u  = nu + nu_t
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    c_k  = nu + nu_t / SIGK
    c_e  = nu + nu_t / SIGE
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    DELV  = UX + VY + WZ
    TAUXX = c_u * (2.0 * UX - TwoThird * DELV) - TwoThird * AKM
    TAUYY = c_u * (2.0 * VY - TwoThird * DELV) - TwoThird * AKM
    TAUZZ = c_u * (2.0 * WZ - TwoThird * DELV) - TwoThird * AKM
    TAUXY = c_u * (UY + VX)
    TAUYZ = c_u * (VZ + WY)
    TAUZX = c_u * (WX + UZ)
    ! エネルギーの拡散 -------------------------------------------------
    R5 = TAUXX * UM + TAUXY * VM + TAUZX * WM + c_t * TX
    S5 = TAUXY * UM + TAUYY * VM + TAUYZ * WM + c_t * TY
    T5 = TAUZX * UM + TAUYZ * VM + TAUZZ * WM + c_t * TZ
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    SH(I,J,K,1) = 0.0
    SH(I,J,K,2) = QH1 * (TAUXX * ETXM + TAUXY * ETYM + TAUZX * ETZM)
    SH(I,J,K,3) = QH1 * (TAUXY * ETXM + TAUYY * ETYM + TAUYZ * ETZM)
    SH(I,J,K,4) = QH1 * (TAUZX * ETXM + TAUYZ * ETYM + TAUZZ * ETZM)
    SH(I,J,K,5) = QH1 * (R5 * ETXM + S5 * ETYM + T5 * ETZM)
    SH(I,J,K,6) = QH1 * c_k * (AKX * ETXM + AKY * ETYM + AKZ * ETZM)
    SH(I,J,K,7) = QH1 * c_e * (AEX * ETXM + AEY * ETYM + AEZ * ETZM)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFET
!***********************************************************************
!**** 拡散項の計算(zeta方向差分用)                                  ****
!***********************************************************************
SUBROUTINE DIFFZE
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, TM, AKM, EPSTM, AMUM, AMUTM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: AKXI, AKET, AKZE, AKX, AKY, AKZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu, nu_t
  REAL    :: c_u, c_t, c_k, c_e
  REAL    :: TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV
  REAL    :: R5, S5, T5
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, TM, AKM, EPSTM, AMUM, AMUTM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  AKXI, AKET, AKZE, AKX, AKY, AKZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, nu_t, c_u, c_t, c_k, c_e, &
  !$  TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV, R5, S5, T5, QH1 &
  !$)
  DO K = KS    , KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J-1,K) .GT. 0.0) .AND. (RHO(I-1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K)   .GT. 0.0) .AND. (RHO(I+1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J+1,K) .GT. 0.0) .AND. (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    ! ZE方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J,K+1))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J,K+1))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J,K+1))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J,K+1))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J,K+1))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J,K+1))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J,K+1))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J,K+1))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J,K+1))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J,K+1))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J,K) +  RHO(I,J,K+1))
    UM    = 0.5 * (   U(I,J,K) +    U(I,J,K+1))
    VM    = 0.5 * (   V(I,J,K) +    V(I,J,K+1))
    WM    = 0.5 * (   W(I,J,K) +    W(I,J,K+1))
    TM    = 0.5 * (   T(I,J,K) +    T(I,J,K+1))
    AKM   = 0.5 * (  AK(I,J,K) +   AK(I,J,K+1))
    EPSTM = 0.5 * (EPST(I,J,K) + EPST(I,J,K+1))
    AMUM  = 0.5 * ( AMU(I,J,K) +  AMU(I,J,K+1))
    AMUTM = 0.5 * (AMUT(I,J,K) + AMUT(I,J,K+1))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J,K+1) &
    &              + U(I+1,J,K) + U(I+1,J,K+1) )
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J,K+1) &
    &              + V(I+1,J,K) + V(I+1,J,K+1) )
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J,K+1) &
    &              + W(I+1,J,K) + W(I+1,J,K+1) )
    UET = 0.25 * ( - U(I,J-1,K) - U(I,J-1,K+1) &
    &              + U(I,J+1,K) + U(I,J+1,K+1) )
    VET = 0.25 * ( - V(I,J-1,K) - V(I,J-1,K+1) &
    &              + V(I,J+1,K) + V(I,J+1,K+1) )
    WET = 0.25 * ( - W(I,J-1,K) - W(I,J-1,K+1) &
    &              + W(I,J+1,K) + W(I,J+1,K+1) )
    UZE = - U(I,J,K) + U(I,J,K+1)
    VZE = - V(I,J,K) + V(I,J,K+1)
    WZE = - W(I,J,K) + W(I,J,K+1)
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J,K+1) &
    &              + T(I+1,J,K) + T(I+1,J,K+1) )
    TET = 0.25 * ( - T(I,J-1,K) - T(I,J-1,K+1) &
    &              + T(I,J+1,K) + T(I,J+1,K+1) )
    TZE = - T(I,J,K) + T(I,J,K+1)
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! k の一階微分の計算 -----------------------------------------------
    ! 計算空間方向一階微分
    AKXI = 0.25 * ( - AK(I-1,J,K) - AK(I-1,J,K+1) &
    &               + AK(I+1,J,K) + AK(I+1,J,K+1) )
    AKET = 0.25 * ( - AK(I,J-1,K) - AK(I,J-1,K+1) &
    &               + AK(I,J+1,K) + AK(I,J+1,K+1) )
    AKZE = - AK(I,J,K) + AK(I,J,K+1)
    ! 物理空間方向一階微分
    AKX = AKXI * XIXM + AKET * ETXM + AKZE * ZEXM
    AKY = AKXI * XIYM + AKET * ETYM + AKZE * ZEYM
    AKZ = AKXI * XIZM + AKET * ETZM + AKZE * ZEZM
    ! epsilon tilde の一階微分の計算 -----------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPST(I-1,J,K) - EPST(I-1,J,K+1) &
    &               + EPST(I+1,J,K) + EPST(I+1,J,K+1) )
    AEET = 0.25 * ( - EPST(I,J-1,K) - EPST(I,J-1,K+1) &
    &               + EPST(I,J+1,K) + EPST(I,J+1,K+1) )
    AEZE = - EPST(I,J,K) + EPST(I,J,K+1)
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 係数計算 ---------------------------------------------------------
    nu   = AMUM  / RHOM
    nu_t = AMUTM / RHOM
    c_u  = nu + nu_t
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    c_k  = nu + nu_t / SIGK
    c_e  = nu + nu_t / SIGE
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    DELV  = UX + VY + WZ
    TAUXX = c_u * (2.0 * UX - TwoThird * DELV) - TwoThird * AKM
    TAUYY = c_u * (2.0 * VY - TwoThird * DELV) - TwoThird * AKM
    TAUZZ = c_u * (2.0 * WZ - TwoThird * DELV) - TwoThird * AKM
    TAUXY = c_u * (UY + VX)
    TAUYZ = c_u * (VZ + WY)
    TAUZX = c_u * (WX + UZ)
    ! エネルギーの拡散 -------------------------------------------------
    R5 = TAUXX * UM + TAUXY * VM + TAUZX * WM + c_t * TX
    S5 = TAUXY * UM + TAUYY * VM + TAUYZ * WM + c_t * TY
    T5 = TAUZX * UM + TAUYZ * VM + TAUZZ * WM + c_t * TZ
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    TH(I,J,K,1) = 0.0
    TH(I,J,K,2) = QH1 * (TAUXX * ZEXM + TAUXY * ZEYM + TAUZX * ZEZM)
    TH(I,J,K,3) = QH1 * (TAUXY * ZEXM + TAUYY * ZEYM + TAUYZ * ZEZM)
    TH(I,J,K,4) = QH1 * (TAUZX * ZEXM + TAUYZ * ZEYM + TAUZZ * ZEZM)
    TH(I,J,K,5) = QH1 * (R5 * ZEXM + S5 * ZEYM + T5 * ZEZM)
    TH(I,J,K,6) = QH1 * c_k * (AKX * ZEXM + AKY * ZEYM + AKZ * ZEZM)
    TH(I,J,K,7) = QH1 * c_e * (AEX * ZEXM + AEY * ZEYM + AEZ * ZEZM)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFZE
! 定義終了 *************************************************************
END SUBROUTINE Turbulence3DEvmAKN
!***********************************************************************
!**** 乱流モデル : RANS, NL-EVM(ke), Craft-Launder-Suga Model(1996) ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence2DNlEvmkeCLS( &
&            PELIM, RG, GAMMA, PR, PRT, &
&            IS, IE, JS, JE, LS, LE, &
&            XIX, XIY, ETX, ETY, AJA, &
&            RHO, U, V, T, AK, EPST, AMU, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: c1 = - 0.1, c2 = 0.1, c3 = 0.26
  REAL, PARAMETER :: c_e1 = 1.44
  REAL, PARAMETER :: SIGK = 1.0, SIGE = 1.3
  REAL, PARAMETER :: aiimin = - 0.66666667, aiimax = 1.3333333
  REAL, PARAMETER :: TwoThird = 0.66666667
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR, PRT
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE), XIY(IS:IE, JS:JE), &
  &                       ETX(IS:IE, JS:JE), ETY(IS:IE, JS:JE), &
  &                       AJA(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE), &
  &                       U(IS:IE, JS:JE), V(IS:IE, JS:JE), &
  &                       T(IS:IE, JS:JE), &
  &                       AK(IS:IE, JS:JE), EPST(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: AMUT(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(OUT) :: DQP(IS:IE, JS:JE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :), SH(:, :, :)
  REAL, ALLOCATABLE :: UXXI(:, :), UYXI(:, :), VXXI(:, :), VYXI(:, :)
  REAL, ALLOCATABLE :: UXET(:, :), UYET(:, :), VXET(:, :), VYET(:, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 6) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, LS: LE), SH(IS: IE, JS: JE, LS: LE))
  ALLOCATE(UXXI(IS: IE, JS: JE), UYXI(IS: IE, JS: JE))
  ALLOCATE(VXXI(IS: IE, JS: JE), VYXI(IS: IE, JS: JE))
  ALLOCATE(UXET(IS: IE, JS: JE), UYET(IS: IE, JS: JE))
  ALLOCATE(VXET(IS: IE, JS: JE), VYET(IS: IE, JS: JE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL PreVel2ndDiff
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 二階微分計算用の前準備                                        ****
!***********************************************************************
SUBROUTINE PreVel2ndDiff
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  UXXI = 0.0
  UYXI = 0.0
  VXXI = 0.0
  VYXI = 0.0
  UXET = 0.0
  UYET = 0.0
  VXET = 0.0
  VYET = 0.0
  ! 速度の二階微分の前準備 +++++++++++++++++++++++++++++++++++++++++++++
  CALL PreVel2ndDiffXI
  CALL PreVel2ndDiffET
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiff
!***********************************************************************
!**** 二階微分計算用の前準備(XI方向)                                ****
!***********************************************************************
SUBROUTINE PreVel2ndDiffXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM
  REAL    :: UXI, VXI, UET, VET
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, XIXM, XIYM, ETXM, ETYM, UXI, VXI, UET, VET &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF(RHO(I,J) .GT. 0.0 .AND. RHO(I+1,J) .GT. 0.0) THEN
    ! XI方向の中間値 ---------------------------------------------------
    XIXM = 0.5 * (XIX(I,J) + XIX(I+1,J))
    XIYM = 0.5 * (XIY(I,J) + XIY(I+1,J))
    ETXM = 0.5 * (ETX(I,J) + ETX(I+1,J))
    ETYM = 0.5 * (ETY(I,J) + ETY(I+1,J))
    ! 計算空間一階微分 -------------------------------------------------
    UXI = - U(I,J) + U(I+1,J)
    UET = - 0.25 * (U(I,J-1) + U(I+1,J-1) - U(I,J+1) - U(I+1,J+1))
    VXI = - V(I,J) + V(I+1,J)
    VET = - 0.25 * (V(I,J-1) + V(I+1,J-1) - V(I,J+1) - V(I+1,J+1))
    ! 物理空間一階微分 -------------------------------------------------
    UXXI(I,J) = UXI * XIXM + UET * ETXM
    UYXI(I,J) = UXI * XIYM + UET * ETYM
    VXXI(I,J) = VXI * XIXM + VET * ETXM
    VYXI(I,J) = VXI * XIYM + VET * ETYM
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiffXI
!***********************************************************************
!**** 二階微分計算用の前準備(ET方向)                                ****
!***********************************************************************
SUBROUTINE PreVel2ndDiffET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM
  REAL    :: UXI, VXI, UET, VET
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, XIXM, XIYM, ETXM, ETYM, UXI, VXI, UET, VET &
  !$)
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF(RHO(I,J) .GT. 0.0 .AND. RHO(I,J+1) .GT. 0.0) THEN
    ! ET方向の中間値 ---------------------------------------------------
    XIXM = 0.5 * (XIX(I,J) + XIX(I,J+1))
    XIYM = 0.5 * (XIY(I,J) + XIY(I,J+1))
    ETXM = 0.5 * (ETX(I,J) + ETX(I,J+1))
    ETYM = 0.5 * (ETY(I,J) + ETY(I,J+1))
    ! 計算空間一階微分 -------------------------------------------------
    UXI = - 0.25 * ( U(I-1,J) - U(I+1,J) + U(I-1,J+1) - U(I+1,J+1) )
    UET = - U(I,J) + U(I,J+1)
    VXI = - 0.25 * ( V(I-1,J) - V(I+1,J) + V(I-1,J+1) - V(I+1,J+1) )
    VET = - V(I,J) + V(I,J+1)
    ! 物理空間一階微分 -------------------------------------------------
    UXET(I,J) = UXI * XIXM + UET * ETXM
    UYET(I,J) = UXI * XIYM + UET * ETYM
    VXET(I,J) = VXI * XIXM + VET * ETXM
    VYET(I,J) = VXI * XIYM + VET * ETYM
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiffET
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: UXI, UET, UX, UY
  REAL    :: VXI, VET, VX, VY
  REAL    :: UXIXI, UETET, UXIET, VXIXI, VETET, VXIET
  REAL    :: UXX, UXY, UYX, UYY, VXX, VXY, VYX, VYY
  REAL    :: KRXI, KRET, KRX, KRY
  REAL    :: S11, S22, S33, S12, SS2, ST, SS
  REAL    :: O12, OO2, OT
  REAL    :: eta
  REAL    :: nu, nu_t
  REAL    :: Rtt
  REAL    :: c_mu, f_mu, c_e2
  REAL    :: c4, c6, c7
  REAL    :: nu_t1, nu_t2, nu_t3, AKAE, AEAK
  REAL    :: a11, a22, a33, a12
  REAL    :: uu, vv, ww, uv
  REAL    :: PP_k
  REAL    :: D, E
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  ! 非拡散項の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, UXI, UET, UX, UY, VXI, VET, VX, VY, &
  !$  UXIXI, UETET, UXIET, VXIXI, VETET, VXIET, &
  !$  UXX, UXY, UYX, UYY, VXX, VXY, VYX, VYY, &
  !$  KRXI, KRET, KRX, KRY, &
  !$  S11, S22, S33, S12, SS2, ST, SS, O12, OO2, OT, &
  !$  eta, nu, nu_t, Rtt, c_mu, f_mu, c_e2, c4, c6, c7, &
  !$  nu_t1, nu_t2, nu_t3, AKAE, AEAK, &
  !$  a11, a22, a33, a12, uu, vv, ww, uv, PP_k, D, E, QH1 &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF(RHO(I,J) .GT. 0.0) THEN
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.5 * (- U(I-1,J) + U(I+1,J))
    VXI = 0.5 * (- V(I-1,J) + V(I+1,J))
    UET = 0.5 * (- U(I,J-1) + U(I,J+1))
    VET = 0.5 * (- V(I,J-1) + V(I,J+1))
    ! 物理空間方向一階微分
    UX  = XIX(I,J) * UXI + ETX(I,J) * UET
    VX  = XIX(I,J) * VXI + ETX(I,J) * VET
    UY  = XIY(I,J) * UXI + ETY(I,J) * UET
    VY  = XIY(I,J) * VXI + ETY(I,J) * VET
    ! 速度の二階微分の計算 ---------------------------------------------
    UXX = - (UXXI(I-1,J) - UXXI(I,J)) * XIX(I,J) &
    &     - (UXET(I,J-1) - UXET(I,J)) * ETX(I,J)
    UXY = - (UXXI(I-1,J) - UXXI(I,J)) * XIY(I,J) &
    &     - (UXET(I,J-1) - UXET(I,J)) * ETY(I,J)
    UYX = - (UYXI(I-1,J) - UYXI(I,J)) * XIX(I,J) &
    &     - (UYET(I,J-1) - UYET(I,J)) * ETX(I,J)
    UYY = - (UYXI(I-1,J) - UYXI(I,J)) * XIY(I,J) &
    &     - (UYET(I,J-1) - UYET(I,J)) * ETY(I,J)
    VXX = - (VXXI(I-1,J) - VXXI(I,J)) * XIX(I,J) &
    &     - (VXET(I,J-1) - VXET(I,J)) * ETX(I,J)
    VXY = - (VXXI(I-1,J) - VXXI(I,J)) * XIY(I,J) &
    &     - (VXET(I,J-1) - VXET(I,J)) * ETY(I,J)
    VYX = - (VYXI(I-1,J) - VYXI(I,J)) * XIX(I,J) &
    &     - (VYET(I,J-1) - VYET(I,J)) * ETX(I,J)
    VYY = - (VYXI(I-1,J) - VYXI(I,J)) * XIY(I,J) &
    &     - (VYET(I,J-1) - VYET(I,J)) * ETY(I,J)
    ! 乱れの速度スケール sqrt(k) の一階微分の計算 ----------------------
    ! 計算空間方向一階微分
    KRXI = 0.5 * (- SQRT(AK(I-1,J)) + SQRT(AK(I+1,J)))
    KRET = 0.5 * (- SQRT(AK(I,J-1)) + SQRT(AK(I,J+1)))
    ! 物理空間方向一階微分
    KRX = KRXI * XIX(I,J) + KRET * ETX(I,J)
    KRY = KRXI * XIY(I,J) + KRET * ETY(I,J)
    ! 歪み速度と渦度 ---------------------------------------------------
    SS  = (UX + VY) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 =          - SS
    S12 = UY + VX
    SS2 = S11 * S11 + S22 * S22 + S33 * S33 + 2.0 * S12 * S12
    O12 = (UY - VX)
    OO2 = 2.0 * O12 * O12
    IF(EPST(I,J) .GT. ZERO) THEN
      AKAE = AK(I,J) / EPST(I,J)
    ELSE
      AKAE = 0.0
    ENDIF
    ST = AKAE * SQRT(0.5 * SS2)
    OT = AKAE * SQRT(0.5 * OO2)
    ! 各モデル係数の計算 -----------------------------------------------
    nu = AMU(I,J) / RHO(I,J)
    IF(EPST(I,J) .GT. ZERO) THEN
      Rtt = AK(I,J)**2 / (nu * EPST(I,J))
    ELSE
      Rtt = 0.0
    ENDIF
    eta  = MAX(ST, OT)
    c_mu = 0.3 / (1.0 + 0.35 * eta**1.5) &
    &    * (1.0 - EXP(-0.36 / MAX(ZERO, EXP(-0.75 * eta))))
    f_mu = 1.0 - EXP(-SQRT(Rtt / 90.0) - (Rtt / 400.0)**2)
    c_e2 = 1.92 * (1.0 - 0.3 * EXP(-Rtt**2))
    c4   =-10.0 * c_mu**2
    c6   =- 5.0 * c_mu**2
    c7   =  5.0 * c_mu**2
    IF(EPST(I,J) .GT. ZERO) THEN
      nu_t = c_mu * f_mu * AK(I,J)**2 / EPST(I,J)
    ELSE
      nu_t = 0.0
    ENDIF
    AMUT(I,J) = RHO(I,J) * nu_t
    ! 非等方レイノルズ応力テンソル a_ij --------------------------------
    IF(AK(I,J) .GT. ZERO .OR. EPST(I,J)**2 .GT. ZERO) THEN
      nu_t1 =-nu_t / AK(I,J)
      nu_t2 = nu_t / EPST(I,J)
      nu_t3 = nu_t * AK(I,J) / EPST(I,J)**2
    ELSE
      nu_t1 = 0.0
      nu_t2 = 0.0
      nu_t3 = 0.0
    ENDIF
    a11 = nu_t1 * S11 &
    &   + nu_t2 * ( &
    &       c1 * (S11 * S11 + S12 * S12 - SS2 / 3.0) &
    &     + c2 * O12 * S12 * 2.0 &
    &     + c3 * (O12 * O12 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &     - c4 * (S11 + S22) * S12 * O12 * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S11 &
    &   )
    a22 = nu_t1 * S22 &
    &   + nu_t2 * ( &
    &       c1 * (S22 * S22 + S12 * S12 - SS2 / 3.0) &
    &     - c2 * O12 * S12 * 2.0 &
    &     + c3 * (O12 * O12 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * (S11 + S22) * S12 * O12 * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S22 &
    &   )
    a33 = nu_t1 * S33 &
    &   + nu_t2 * ( &
    &       c1 * (S33 * S33 - SS2 / 3.0) - c3 * OO2 / 3.0 &
    &   ) &
    &   + nu_t3 * (c6 * SS2 + c7 * OO2) * S33
    a12 = nu_t1 * S12 &
    &   + nu_t2 * ( &
    &       c1 * (S11 + S22) * S12 - c2 * O12 * (S11 - S22) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * (S11 * S11 - S22 * S22) * O12 &
    &     + (c6 * SS2 + c7 * OO2) * S12 &
    &   )
    a11 = MIN(aiimax, MAX(aiimin, a11))
    a22 = MIN(aiimax, MAX(aiimin, a22))
    a33 = MIN(aiimax, MAX(aiimin, a33))
    a12 = MIN(1.0, MAX(-1.0, a12))
    ! レイノルズ応力 ---------------------------------------------------
    uu = AK(I,J) * (a11 + TwoThird)
    vv = AK(I,J) * (a22 + TwoThird)
    ww = AK(I,J) * (a33 + TwoThird)
    uv = AK(I,J) * a12
    ! 低レイノルズ補正項 D ---------------------------------------------
    D = 2.0 * nu * (KRX**2 + KRY**2)
    ! 低レイノルズ補正項 E ---------------------------------------------
    IF(Rtt .LE. 250.0) THEN
      E = 0.0022 * ST * nu_t * AK(I,J) * AKAE * ( &
      &     UXX**2 + VXX**2 + UYY**2 + VYY**2 &
      &   + UXY**2 + VXY**2 + UYX**2 + VYX**2 &
      & )
    ELSE
      E = 0.0
    ENDIF
    ! レイノルズ応力の生成(平均歪み増加率) -----------------------------
    PP_k = - uu * UX - vv * VY - uv * (VX + UY)
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      PP_k = MAX(0.0, MIN((EPST(I,J) + D) / PELIM, PP_k))
    ENDIF
    ! ソース項の和 -----------------------------------------------------
    QH1 = RHO(I,J) / AJA(I,J)
    AEAK = MIN( &
    &      MAX(SQRT(0.5 * SS2), SQRT((EPST(I,J) + D) / nu)), &
    &      EPST(I,J) / MAX(ZERO, AK(I,J)) &
    &    )
    DQP(I,J,5) = QH1 * (PP_k - (EPST(I,J) + D))
    DQP(I,J,6) = QH1 * ( &
    &            AEAK * (c_e1 * PP_k - c_e2 * EPST(I,J)) &
    &          + E &
    &          )
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH = 0.0
  SH = 0.0
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL DIFFXI
  CALL DIFFET
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, L)
  DO L = LS, LE
  !$OMP DO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J-1) .GT. 0.0) .AND. (RHO(I-1,J) .GT. 0.0) .AND. &
  &   (RHO(I,J)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,L) = (-RH(I-1,J  ,L) + RH(I,J,L)) &
    &          + (-SH(I  ,J-1,L) + SH(I,J,L))
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項の計算(xi方向差分用)                                    ****
!***********************************************************************
SUBROUTINE DIFFXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: RHOM, UM, VM, AKM, EPSTM, AMUM
  REAL    :: UXI, UET, UX, UY
  REAL    :: VXI, VET, VX, VY
  REAL    :: TXI, TET, TX, TY
  REAL    :: AKXI, AKET, AKX, AKY
  REAL    :: AEXI, AEET, AEX, AEY
  REAL    :: S11, S22, S33, S12, SS2, ST, SS
  REAL    :: O12, OO2, OT
  REAL    :: eta
  REAL    :: nu, nu_t
  REAL    :: Rtt, c_mu, f_mu
  REAL    :: c4, c6, c7
  REAL    :: nu_t1, nu_t2, nu_t3, AKAE
  REAL    :: a11, a22, a33, a12
  REAL    :: uu, vv, ww, uv
  REAL    :: c_t, c_k, c_e
  REAL    :: TAUXX, TAUYY, TAUXY
  REAL    :: R4, S4
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, &
  !$  RHOM, UM, VM, AKM, EPSTM, AMUM, &
  !$  UXI, UET, UX, UY, VXI, VET, VX, VY, TXI, TET, TX, TY, &
  !$  AKXI, AKET, AKX, AKY, AEXI, AEET, AEX, AEY, &
  !$  S11, S22, S33, S12, SS2, ST, SS, O12, OO2, OT, &
  !$  eta, nu, nu_t, Rtt, c_mu, f_mu, c4, c6, c7, &
  !$  nu_t1, nu_t2, nu_t3, AKAE, &
  !$  a11, a22, a33, a12, uu, vv, ww, uv, &
  !$  c_t, c_k, c_e, TAUXX, TAUYY, TAUXY, R4, S4, QH1 &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF((RHO(I,J) .GT. 0.0) .AND. (RHO(I+1,J) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J) + XIX(I+1,J))
    XIYM = 0.5 * (XIY(I,J) + XIY(I+1,J))
    ETXM = 0.5 * (ETX(I,J) + ETX(I+1,J))
    ETYM = 0.5 * (ETY(I,J) + ETY(I+1,J))
    AJAM = 0.5 * (AJA(I,J) + AJA(I+1,J))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J) +  RHO(I+1,J))
    UM    = 0.5 * (   U(I,J) +    U(I+1,J))
    VM    = 0.5 * (   V(I,J) +    V(I+1,J))
    AKM   = 0.5 * (  AK(I,J) +   AK(I+1,J))
    EPSTM = 0.5 * (EPST(I,J) + EPST(I+1,J))
    AMUM  = 0.5 * ( AMU(I,J) +  AMU(I+1,J))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = - U(I,J) + U(I+1,J)
    UET = 0.25 * (- U(I,J-1) - U(I+1,J-1) + U(I,J+1) + U(I+1,J+1))
    VXI = - V(I,J) + V(I+1,J)
    VET = 0.25 * (- V(I,J-1) - V(I+1,J-1) + V(I,J+1) + V(I+1,J+1))
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM
    UY = UXI * XIYM + UET * ETYM
    VX = VXI * XIXM + VET * ETXM
    VY = VXI * XIYM + VET * ETYM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = - T(I,J) + T(I+1,J)
    TET = 0.25 * (- T(I,J-1) - T(I+1,J-1) + T(I,J+1) + T(I+1,J+1))
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM
    TY = TXI * XIYM + TET * ETYM
    ! k の一階微分の計算 -----------------------------------------------
    ! 計算空間方向一階微分
    AKXI = - AK(I,J) + AK(I+1,J)
    AKET = 0.25 * (- AK(I,J-1) - AK(I+1,J-1) + AK(I,J+1) + AK(I+1,J+1))
    ! 物理空間方向一階微分
    AKX = AKXI * XIXM + AKET * ETXM
    AKY = AKXI * XIYM + AKET * ETYM
    ! epsilon tilde の一階微分の計算 -----------------------------------
    ! 計算空間方向一階微分
    AEXI = - EPST(I,J) + EPST(I+1,J)
    AEET = 0.25 * ( - EPST(I,J-1) - EPST(I+1,J-1) &
    &               + EPST(I,J+1) + EPST(I+1,J+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM
    AEY = AEXI * XIYM + AEET * ETYM
    ! 歪み速度と渦度 ---------------------------------------------------
    SS  = (UX + VY) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 =          - SS
    S12 = UY + VX
    SS2 = S11 * S11 + S22 * S22 + S33 * S33 + 2.0 * S12 * S12
    O12 = (UY - VX)
    OO2 = 2.0 * O12 * O12
    IF(EPSTM .GT. ZERO) THEN
      AKAE = AKM / EPSTM
    ELSE
      AKAE = 0.0
    ENDIF
    ST = AKAE * SQRT(0.5 * SS2)
    OT = AKAE * SQRT(0.5 * OO2)
    ! 非等方レイノルズ応力テンソル a_ij の係数の計算 -------------------
    nu = AMUM / RHOM
    IF(EPSTM .GT. ZERO) THEN
      Rtt = AKM**2 / (nu * EPSTM)
    ELSE
      Rtt = 0.0
    ENDIF
    eta  = MAX(ST, OT)
    c_mu = 0.3 / (1.0 + 0.35 * eta**1.5) &
    &    * (1.0 - EXP(-0.36 / MAX(ZERO, EXP(-0.75 * eta))))
    f_mu = 1.0 - EXP(-SQRT(Rtt / 90.0) - (Rtt / 400.0)**2)
    c4   =-10.0 * c_mu**2
    c6   =- 5.0 * c_mu**2
    c7   =  5.0 * c_mu**2
    IF(EPSTM .GT. ZERO) THEN
      nu_t = c_mu * f_mu * AKM**2 / EPSTM
    ELSE
      nu_t = 0.0
    ENDIF
    ! 非等方レイノルズ応力テンソル a_ij --------------------------------
    IF(AKM .GT. ZERO .OR. EPSTM**2 .GT. ZERO) THEN
      nu_t1 =-nu_t / AKM
      nu_t2 = nu_t / EPSTM
      nu_t3 = nu_t * AKM / EPSTM**2
    ELSE
      nu_t1 = 0.0
      nu_t2 = 0.0
      nu_t3 = 0.0
    ENDIF
    a11 = nu_t1 * S11 &
    &   + nu_t2 * ( &
    &       c1 * (S11 * S11 + S12 * S12 - SS2 / 3.0) &
    &     + c2 * O12 * S12 * 2.0 &
    &     + c3 * (O12 * O12 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &     - c4 * (S11 + S22) * S12 * O12 * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S11 &
    &   )
    a22 = nu_t1 * S22 &
    &   + nu_t2 * ( &
    &       c1 * (S22 * S22 + S12 * S12 - SS2 / 3.0) &
    &     - c2 * O12 * S12 * 2.0 &
    &     + c3 * (O12 * O12 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * (S11 + S22) * S12 * O12 * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S22 &
    &   )
    a33 = nu_t1 * S33 &
    &   + nu_t2 * ( &
    &       c1 * (S33 * S33 - SS2 / 3.0) - c3 * OO2 / 3.0 &
    &   ) &
    &   + nu_t3 * (c6 * SS2 + c7 * OO2) * S33
    a12 = nu_t1 * S12 &
    &   + nu_t2 * ( &
    &       c1 * (S11 + S22) * S12 - c2 * O12 * (S11 - S22) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * (S11 * S11 - S22 * S22) * O12 &
    &     + (c6 * SS2 + c7 * OO2) * S12 &
    &   )
    a11 = MIN(aiimax, MAX(aiimin, a11))
    a22 = MIN(aiimax, MAX(aiimin, a22))
    a33 = MIN(aiimax, MAX(aiimin, a33))
    a12 = MIN(1.0, MAX(-1.0, a12))
    ! レイノルズ応力 ---------------------------------------------------
    uu = AKM * (a11 + TwoThird)
    vv = AKM * (a22 + TwoThird)
    ww = AKM * (a33 + TwoThird)
    uv = AKM * a12
    ! 拡散係数 ---------------------------------------------------------
    c_t = (nu / PR + nu_t / PRT) * GAMMA * RG / (GAMMA - 1.0)
    c_k = nu + nu_t / SIGK
    c_e = nu + nu_t / SIGE
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    TAUXX = nu * S11 - uu
    TAUYY = nu * S22 - vv
    TAUXY = nu * S12 - uv
    ! エネルギーの拡散 -------------------------------------------------
    R4 = TAUXX * UM + TAUXY * VM + c_t * TX
    S4 = TAUXY * UM + TAUYY * VM + c_t * TY
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    RH(I,J,1) = 0.0
    RH(I,J,2) = QH1 * (TAUXX * XIXM + TAUXY * XIYM)
    RH(I,J,3) = QH1 * (TAUXY * XIXM + TAUYY * XIYM)
    RH(I,J,4) = QH1 * (R4 * XIXM + S4 * XIYM)
    RH(I,J,5) = QH1 * c_k * (AKX * XIXM + AKY * XIYM)
    RH(I,J,6) = QH1 * c_e * (AEX * XIXM + AEY * XIYM)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFXI
!***********************************************************************
!**** 拡散項の計算(eta方向差分用)                                   ****
!***********************************************************************
SUBROUTINE DIFFET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: RHOM, UM, VM, AKM, EPSTM, AMUM
  REAL    :: UXI, UET, UX, UY
  REAL    :: VXI, VET, VX, VY
  REAL    :: TXI, TET, TX, TY
  REAL    :: AKXI, AKET, AKX, AKY
  REAL    :: AEXI, AEET, AEX, AEY
  REAL    :: S11, S22, S33, S12, SS2, ST, SS
  REAL    :: O12, OO2, OT
  REAL    :: eta
  REAL    :: nu, nu_t
  REAL    :: Rtt, c_mu, f_mu
  REAL    :: c4, c6, c7
  REAL    :: nu_t1, nu_t2, nu_t3, AKAE
  REAL    :: a11, a22, a33, a12
  REAL    :: uu, vv, ww, uv
  REAL    :: c_t, c_k, c_e
  REAL    :: TAUXX, TAUYY, TAUXY
  REAL    :: R4, S4
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, &
  !$  RHOM, UM, VM, AKM, EPSTM, AMUM, &
  !$  UXI, UET, UX, UY, VXI, VET, VX, VY, TXI, TET, TX, TY, &
  !$  AKXI, AKET, AKX, AKY, AEXI, AEET, AEX, AEY, &
  !$  S11, S22, S33, S12, SS2, ST, SS, O12, OO2, OT, &
  !$  eta, nu, nu_t, Rtt, c_mu, f_mu, c4, c6, c7, &
  !$  nu_t1, nu_t2, nu_t3, AKAE, &
  !$  a11, a22, a33, a12, uu, vv, ww, uv, &
  !$  c_t, c_k, c_e, TAUXX, TAUYY, TAUXY, R4, S4, QH1 &
  !$)
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF((RHO(I,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0)) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J) + XIX(I,J+1))
    XIYM = 0.5 * (XIY(I,J) + XIY(I,J+1))
    ETXM = 0.5 * (ETX(I,J) + ETX(I,J+1))
    ETYM = 0.5 * (ETY(I,J) + ETY(I,J+1))
    AJAM = 0.5 * (AJA(I,J) + AJA(I,J+1))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J) +  RHO(I,J+1))
    UM    = 0.5 * (   U(I,J) +    U(I,J+1))
    VM    = 0.5 * (   V(I,J) +    V(I,J+1))
    AKM   = 0.5 * (  AK(I,J) +   AK(I,J+1))
    EPSTM = 0.5 * (EPST(I,J) + EPST(I,J+1))
    AMUM  = 0.5 * ( AMU(I,J) +  AMU(I,J+1))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * (- U(I-1,J) - U(I-1,J+1) + U(I+1,J) + U(I+1,J+1))
    UET = - U(I,J) + U(I,J+1)
    VXI = 0.25 * (- V(I-1,J) - V(I-1,J+1) + V(I+1,J) + V(I+1,J+1))
    VET = - V(I,J) + V(I,J+1)
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM
    UY = UXI * XIYM + UET * ETYM
    VX = VXI * XIXM + VET * ETXM
    VY = VXI * XIYM + VET * ETYM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * (- T(I-1,J) - T(I-1,J+1) + T(I+1,J) + T(I+1,J+1))
    TET = - T(I,J) + T(I,J+1)
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM
    TY = TXI * XIYM + TET * ETYM
    ! k の一階微分の計算 -----------------------------------------------
    ! 計算空間方向一階微分
    AKXI = 0.25 * (- AK(I-1,J) - AK(I-1,J+1) + AK(I+1,J) + AK(I+1,J+1))
    AKET = - AK(I,J) + AK(I,J+1)
    ! 物理空間方向一階微分
    AKX = AKXI * XIXM + AKET * ETXM
    AKY = AKXI * XIYM + AKET * ETYM
    ! epsilon tilde の一階微分の計算 -----------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPST(I-1,J) - EPST(I-1,J+1) &
    &               + EPST(I+1,J) + EPST(I+1,J+1) )
    AEET = - EPST(I,J) + EPST(I,J+1)
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM
    AEY = AEXI * XIYM + AEET * ETYM
    ! 歪み速度と渦度 ---------------------------------------------------
    SS  = (UX + VY) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 =          - SS
    S12 = UY + VX
    SS2 = S11 * S11 + S22 * S22 + S33 * S33 + 2.0 * S12 * S12
    O12 = (UY - VX)
    OO2 = 2.0 * O12 * O12
    IF(EPSTM .GT. ZERO) THEN
      AKAE = AKM / EPSTM
    ELSE
      AKAE = 0.0
    ENDIF
    ST = AKAE * SQRT(0.5 * SS2)
    OT = AKAE * SQRT(0.5 * OO2)
    ! 非等方レイノルズ応力テンソル a_ij の係数の計算 -------------------
    nu = AMUM / RHOM
    IF(EPSTM .GT. ZERO) THEN
      Rtt = AKM**2 / (nu * EPSTM)
    ELSE
      Rtt = 0.0
    ENDIF
    eta  = MAX(ST, OT)
    c_mu = 0.3 / (1.0 + 0.35 * eta**1.5) &
    &    * (1.0 - EXP(-0.36 / MAX(ZERO, EXP(-0.75 * eta))))
    f_mu = 1.0 - EXP(-SQRT(Rtt / 90.0) - (Rtt / 400.0)**2)
    c4   =-10.0 * c_mu**2
    c6   =- 5.0 * c_mu**2
    c7   =  5.0 * c_mu**2
    IF(EPSTM .GT. ZERO) THEN
      nu_t = c_mu * f_mu * AKM**2 / EPSTM
    ELSE
      nu_t = 0.0
    ENDIF
    ! 非等方レイノルズ応力テンソル a_ij --------------------------------
    IF(AKM .GT. ZERO .OR. EPSTM**2 .GT. ZERO) THEN
      nu_t1 =-nu_t / AKM
      nu_t2 = nu_t / EPSTM
      nu_t3 = nu_t * AKM / EPSTM**2
    ELSE
      nu_t1 = 0.0
      nu_t2 = 0.0
      nu_t3 = 0.0
    ENDIF
    a11 = nu_t1 * S11 &
    &   + nu_t2 * ( &
    &       c1 * (S11 * S11 + S12 * S12 - SS2 / 3.0) &
    &     + c2 * O12 * S12 * 2.0 &
    &     + c3 * (O12 * O12 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &     - c4 * (S11 + S22) * S12 * O12 * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S11 &
    &   )
    a22 = nu_t1 * S22 &
    &   + nu_t2 * ( &
    &       c1 * (S22 * S22 + S12 * S12 - SS2 / 3.0) &
    &     - c2 * O12 * S12 * 2.0 &
    &     + c3 * (O12 * O12 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * (S11 + S22) * S12 * O12 * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S22 &
    &   )
    a33 = nu_t1 * S33 &
    &   + nu_t2 * ( &
    &       c1 * (S33 * S33 - SS2 / 3.0) - c3 * OO2 / 3.0 &
    &   ) &
    &   + nu_t3 * (c6 * SS2 + c7 * OO2) * S33
    a12 = nu_t1 * S12 &
    &   + nu_t2 * ( &
    &       c1 * (S11 + S22) * S12 - c2 * O12 * (S11 - S22) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * (S11 * S11 - S22 * S22) * O12 &
    &     + (c6 * SS2 + c7 * OO2) * S12 &
    &   )
    a11 = MIN(aiimax, MAX(aiimin, a11))
    a22 = MIN(aiimax, MAX(aiimin, a22))
    a33 = MIN(aiimax, MAX(aiimin, a33))
    a12 = MIN(1.0, MAX(-1.0, a12))
    ! レイノルズ応力 ---------------------------------------------------
    uu = AKM * (a11 + TwoThird)
    vv = AKM * (a22 + TwoThird)
    ww = AKM * (a33 + TwoThird)
    uv = AKM * a12
    ! 拡散係数 ---------------------------------------------------------
    c_t = (nu / PR + nu_t / PRT) * GAMMA * RG / (GAMMA - 1.0)
    c_k = nu + nu_t / SIGK
    c_e = nu + nu_t / SIGE
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    TAUXX = nu * S11 - uu
    TAUYY = nu * S22 - vv
    TAUXY = nu * S12 - uv
    ! エネルギーの拡散 -------------------------------------------------
    R4 = TAUXX * UM + TAUXY * VM + c_t * TX
    S4 = TAUXY * UM + TAUYY * VM + c_t * TY
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    SH(I,J,1) = 0.0
    SH(I,J,2) = QH1 * (TAUXX * ETXM + TAUXY * ETYM)
    SH(I,J,3) = QH1 * (TAUXY * ETXM + TAUYY * ETYM)
    SH(I,J,4) = QH1 * (R4 * ETXM + S4 * ETYM)
    SH(I,J,5) = QH1 * c_k * (AKX * ETXM + AKY * ETYM)
    SH(I,J,6) = QH1 * c_e * (AEX * ETXM + AEY * ETYM)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFET
! 定義終了 *************************************************************
END SUBROUTINE Turbulence2DNlEvmkeCLS
!***********************************************************************
!**** 乱流モデル : RANS, NL-EVM(ke), Craft-Launder-Suga Model(1996) ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence3DNlEvmkeCLS( &
&            PELIM, RG, GAMMA, PR, PRT, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, AJA, &
&            RHO, U, V, W, T, AK, EPST, AMU, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: c1 = - 0.1, c2 = 0.1, c3 = 0.26
  REAL, PARAMETER :: c_e1 = 1.44
  REAL, PARAMETER :: SIGK = 1.0, SIGE = 1.3
  REAL, PARAMETER :: aiimin = - 0.66666667, aiimax = 1.3333333
  REAL, PARAMETER :: TwoThird = 0.66666667
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR, PRT
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE, KS:KE), &
  &                       XIY(IS:IE, JS:JE, KS:KE), &
  &                       XIZ(IS:IE, JS:JE, KS:KE), &
  &                       ETX(IS:IE, JS:JE, KS:KE), &
  &                       ETY(IS:IE, JS:JE, KS:KE), &
  &                       ETZ(IS:IE, JS:JE, KS:KE), &
  &                       ZEX(IS:IE, JS:JE, KS:KE), &
  &                       ZEY(IS:IE, JS:JE, KS:KE), &
  &                       ZEZ(IS:IE, JS:JE, KS:KE), &
  &                       AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE, KS:KE), &
  &                       U(IS:IE, JS:JE, KS:KE), &
  &                       V(IS:IE, JS:JE, KS:KE), &
  &                       W(IS:IE, JS:JE, KS:KE), &
  &                       T(IS:IE, JS:JE, KS:KE), &
  &                       AK(IS:IE, JS:JE, KS:KE), &
  &                       EPST(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: AMUT(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(OUT) :: DQP(IS:IE, JS:JE, KS:KE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :, :), SH(:, :, :, :), TH(:, :, :, :)
  REAL, ALLOCATABLE :: UXXI(:, :, :), UYXI(:, :, :), UZXI(:, :, :), &
  &                    VXXI(:, :, :), VYXI(:, :, :), VZXI(:, :, :), &
  &                    WXXI(:, :, :), WYXI(:, :, :), WZXI(:, :, :)
  REAL, ALLOCATABLE :: UXET(:, :, :), UYET(:, :, :), UZET(:, :, :), &
  &                    VXET(:, :, :), VYET(:, :, :), VZET(:, :, :), &
  &                    WXET(:, :, :), WYET(:, :, :), WZET(:, :, :)
  REAL, ALLOCATABLE :: UXZE(:, :, :), UYZE(:, :, :), UZZE(:, :, :), &
  &                    VXZE(:, :, :), VYZE(:, :, :), VZZE(:, :, :), &
  &                    WXZE(:, :, :), WYZE(:, :, :), WZZE(:, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 7) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(SH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(TH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(UXXI(IS: IE, JS: JE, KS: KE), UYXI(IS: IE, JS: JE, KS: KE))
  ALLOCATE(UZXI(IS: IE, JS: JE, KS: KE), VXXI(IS: IE, JS: JE, KS: KE))
  ALLOCATE(VYXI(IS: IE, JS: JE, KS: KE), VZXI(IS: IE, JS: JE, KS: KE))
  ALLOCATE(WXXI(IS: IE, JS: JE, KS: KE), WYXI(IS: IE, JS: JE, KS: KE))
  ALLOCATE(WZXI(IS: IE, JS: JE, KS: KE))
  ALLOCATE(UXET(IS: IE, JS: JE, KS: KE), UYET(IS: IE, JS: JE, KS: KE))
  ALLOCATE(UZET(IS: IE, JS: JE, KS: KE), VXET(IS: IE, JS: JE, KS: KE))
  ALLOCATE(VYET(IS: IE, JS: JE, KS: KE), VZET(IS: IE, JS: JE, KS: KE))
  ALLOCATE(WXET(IS: IE, JS: JE, KS: KE), WYET(IS: IE, JS: JE, KS: KE))
  ALLOCATE(WZET(IS: IE, JS: JE, KS: KE))
  ALLOCATE(UXZE(IS: IE, JS: JE, KS: KE), UYZE(IS: IE, JS: JE, KS: KE))
  ALLOCATE(UZZE(IS: IE, JS: JE, KS: KE), VXZE(IS: IE, JS: JE, KS: KE))
  ALLOCATE(VYZE(IS: IE, JS: JE, KS: KE), VZZE(IS: IE, JS: JE, KS: KE))
  ALLOCATE(WXZE(IS: IE, JS: JE, KS: KE), WYZE(IS: IE, JS: JE, KS: KE))
  ALLOCATE(WZZE(IS: IE, JS: JE, KS: KE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL PreVel2ndDiff
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 二階微分計算用の前準備                                        ****
!***********************************************************************
SUBROUTINE PreVel2ndDiff
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  UXXI = 0.0
  UYXI = 0.0
  UZXI = 0.0
  VXXI = 0.0
  VYXI = 0.0
  VZXI = 0.0
  WXXI = 0.0
  WYXI = 0.0
  WZXI = 0.0
  UXET = 0.0
  UYET = 0.0
  UZET = 0.0
  VXET = 0.0
  VYET = 0.0
  VZET = 0.0
  WXET = 0.0
  WYET = 0.0
  WZET = 0.0
  UXZE = 0.0
  UYZE = 0.0
  UZZE = 0.0
  VXZE = 0.0
  VYZE = 0.0
  VZZE = 0.0
  WXZE = 0.0
  WYZE = 0.0
  WZZE = 0.0
  ! 速度の二階微分の前準備 +++++++++++++++++++++++++++++++++++++++++++++
  CALL PreVel2ndDiffXI
  CALL PreVel2ndDiffET
  CALL PreVel2ndDiffZE
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiff
!***********************************************************************
!**** 二階微分計算用の前準備(XI方向)                                ****
!***********************************************************************
SUBROUTINE PreVel2ndDiffXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF(RHO(I,J,K) .GT. 0.0 .AND. RHO(I+1,J,K) .GT. 0.0) THEN
    ! XI方向の中間値 ---------------------------------------------------
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I+1,J,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I+1,J,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I+1,J,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I+1,J,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I+1,J,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I+1,J,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I+1,J,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I+1,J,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I+1,J,K))
    ! 計算空間一階微分 -------------------------------------------------
    UXI = - U(I,J,K) + U(I+1,J,K)
    UET = - 0.25 * ( U(I,J-1,K) + U(I+1,J-1,K) &
    &              - U(I,J+1,K) - U(I+1,J+1,K) )
    UZE = - 0.25 * ( U(I,J,K-1) + U(I+1,J,K-1) &
    &              - U(I,J,K+1) - U(I+1,J,K+1) )
    VXI = - V(I,J,K) + V(I+1,J,K)
    VET = - 0.25 * ( V(I,J-1,K) + V(I+1,J-1,K) &
    &              - V(I,J+1,K) - V(I+1,J+1,K) )
    VZE = - 0.25 * ( V(I,J,K-1) + V(I+1,J,K-1) &
    &              - V(I,J,K+1) - V(I+1,J,K+1) )
    WXI = - W(I,J,K) + W(I+1,J,K)
    WET = - 0.25 * ( W(I,J-1,K) + W(I+1,J-1,K) &
    &              - W(I,J+1,K) - W(I+1,J+1,K) )
    WZE = - 0.25 * ( W(I,J,K-1) + W(I+1,J,K-1) &
    &              - W(I,J,K+1) - W(I+1,J,K+1) )
    ! 物理空間一階微分 -------------------------------------------------
    UXXI(I,J,K) = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UYXI(I,J,K) = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZXI(I,J,K) = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VXXI(I,J,K) = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VYXI(I,J,K) = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZXI(I,J,K) = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WXXI(I,J,K) = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WYXI(I,J,K) = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZXI(I,J,K) = WXI * XIZM + WET * ETZM + WZE * ZEZM
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiffXI
!***********************************************************************
!**** 二階微分計算用の前準備(ET方向)                                ****
!***********************************************************************
SUBROUTINE PreVel2ndDiffET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF(RHO(I,J,K) .GT. 0.0 .AND. RHO(I,J+1,K) .GT. 0.0) THEN
    ! ET方向の中間値 ---------------------------------------------------
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J+1,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J+1,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J+1,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J+1,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J+1,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J+1,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J+1,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J+1,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J+1,K))
    ! 計算空間一階微分 -------------------------------------------------
    UXI = - 0.25 * ( U(I-1,J  ,K) - U(I+1,J  ,K) &
    &              + U(I-1,J+1,K) - U(I+1,J+1,K) )
    UET = - U(I,J,K) + U(I,J+1,K)
    UZE = - 0.25 * ( U(I,J,K-1) + U(I,J+1,K-1) &
    &              - U(I,J,K+1) - U(I,J+1,K+1) )
    VXI = - 0.25 * ( V(I-1,J  ,K) - V(I+1,J  ,K) &
    &              + V(I-1,J+1,K) - V(I+1,J+1,K) )
    VET = - V(I,J,K) + V(I,J+1,K)
    VZE = - 0.25 * ( V(I,J,K-1) + V(I,J+1,K-1) &
    &              - V(I,J,K+1) - V(I,J+1,K+1) )
    WXI = - 0.25 * ( W(I-1,J  ,K) - W(I+1,J  ,K) &
    &              + W(I-1,J+1,K) - W(I+1,J+1,K) )
    WET = - W(I,J,K) + W(I,J+1,K)
    WZE = - 0.25 * ( W(I,J,K-1) + W(I,J+1,K-1) &
    &              - W(I,J,K+1) - W(I,J+1,K+1) )
    ! 物理空間一階微分 -------------------------------------------------
    UXET(I,J,K) = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UYET(I,J,K) = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZET(I,J,K) = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VXET(I,J,K) = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VYET(I,J,K) = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZET(I,J,K) = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WXET(I,J,K) = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WYET(I,J,K) = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZET(I,J,K) = WXI * XIZM + WET * ETZM + WZE * ZEZM
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiffET
!***********************************************************************
!**** 二階微分計算用の前準備(ZE方向)                                ****
!***********************************************************************
SUBROUTINE PreVel2ndDiffZE
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE &
  !$)
  DO K = KS    , KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF(RHO(I,J,K) .GT. 0.0 .AND. RHO(I,J,K+1) .GT. 0.0) THEN
    ! ZE方向の中間値 ---------------------------------------------------
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J,K+1))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J,K+1))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J,K+1))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J,K+1))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J,K+1))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J,K+1))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J,K+1))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J,K+1))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J,K+1))
    ! 計算空間一階微分 -------------------------------------------------
    UXI = - 0.25 * ( U(I-1,J,K  ) - U(I+1,J,K  ) &
    &              + U(I-1,J,K+1) - U(I+1,J,K+1) )
    UET = - 0.25 * ( U(I,J-1,K  ) - U(I,J+1,K  ) &
    &              + U(I,J-1,K+1) - U(I,J+1,K+1) )
    UZE = - U(I,J,K) + U(I,J,K+1)
    VXI = - 0.25 * ( V(I-1,J,K  ) - V(I+1,J,K  ) &
    &              + V(I-1,J,K+1) - V(I+1,J,K+1) )
    VET = - 0.25 * ( V(I,J-1,K  ) - V(I,J+1,K  ) &
    &              + V(I,J-1,K+1) - V(I,J+1,K+1) )
    VZE = - V(I,J,K) + V(I,J,K+1)
    WXI = - 0.25 * ( W(I-1,J,K  ) - W(I+1,J,K  ) &
    &              + W(I-1,J,K+1) - W(I+1,J,K+1) )
    WET = - 0.25 * ( W(I,J-1,K  ) - W(I,J+1,K  ) &
    &              + W(I,J-1,K+1) - W(I,J+1,K+1) )
    WZE = - W(I,J,K) + W(I,J,K+1)
    ! 物理空間一階微分 -------------------------------------------------
    UXZE(I,J,K) = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UYZE(I,J,K) = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZZE(I,J,K) = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VXZE(I,J,K) = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VYZE(I,J,K) = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZZE(I,J,K) = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WXZE(I,J,K) = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WYZE(I,J,K) = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZZE(I,J,K) = WXI * XIZM + WET * ETZM + WZE * ZEZM
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE PreVel2ndDiffZE
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  REAL    :: UX, VX, WX, UY, VY, WY, UZ, VZ, WZ
  REAL    :: UXIXI, UETET, UZEZE, UXIET, UETZE, UZEXI, &
  &          VXIXI, VETET, VZEZE, VXIET, VETZE, VZEXI, &
  &          WXIXI, WETET, WZEZE, WXIET, WETZE, WZEXI
  REAL    :: UXX, UXY, UXZ, UYX, UYY, UYZ, UZX, UZY, UZZ, &
  &          VXX, VXY, VXZ, VYX, VYY, VYZ, VZX, VZY, VZZ, &
  &          WXX, WXY, WXZ, WYX, WYY, WYZ, WZX, WZY, WZZ
  REAL    :: KRXI, KRET, KRZE, KRX, KRY, KRZ
  REAL    :: S11, S22, S33, S12, S23, S31, SS2, ST, SS
  REAL    :: O12, O23, O31, OO2, OT
  REAL    :: eta
  REAL    :: nu, nu_t
  REAL    :: Rtt
  REAL    :: c_mu, f_mu, c_e2
  REAL    :: c4, c6, c7
  REAL    :: nu_t1, nu_t2, nu_t3, AKAE, AEAK
  REAL    :: a11, a22, a33, a12, a23, a31
  REAL    :: uu, vv, ww, uv, vw, wu
  REAL    :: PP_k
  REAL    :: D, E
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  ! 非拡散項の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE, &
  !$  UX, VX, WX, UY, VY, WY, UZ, VZ, WZ, &
  !$  UXIXI, UETET, UZEZE, UXIET, UETZE, UZEXI, &
  !$  VXIXI, VETET, VZEZE, VXIET, VETZE, VZEXI, &
  !$  WXIXI, WETET, WZEZE, WXIET, WETZE, WZEXI, &
  !$  UXX, UXY, UXZ, UYX, UYY, UYZ, UZX, UZY, UZZ, &
  !$  VXX, VXY, VXZ, VYX, VYY, VYZ, VZX, VZY, VZZ, &
  !$  WXX, WXY, WXZ, WYX, WYY, WYZ, WZX, WZY, WZZ, &
  !$  KRXI, KRET, KRZE, KRX, KRY, KRZ, &
  !$  S11, S22, S33, S12, S23, S31, SS2, ST, SS, &
  !$  O12, O23, O31, OO2, OT, &
  !$  eta, nu, nu_t, Rtt, c_mu, f_mu, c_e2, c4, c6, c7, &
  !$  nu_t1, nu_t2, nu_t3, AKAE, AEAK, &
  !$  a11, a22, a33, a12, a23, a31, uu, vv, ww, uv, vw, wu, &
  !$  PP_k, D, E, QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
    IF(RHO(I,J,K) .GT. 0.0) THEN
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = (-U(I-1,J,K) + U(I+1,J,K)) * 0.5
    UET = (-U(I,J-1,K) + U(I,J+1,K)) * 0.5
    UZE = (-U(I,J,K-1) + U(I,J,K+1)) * 0.5
    VXI = (-V(I-1,J,K) + V(I+1,J,K)) * 0.5
    VET = (-V(I,J-1,K) + V(I,J+1,K)) * 0.5
    VZE = (-V(I,J,K-1) + V(I,J,K+1)) * 0.5
    WXI = (-W(I-1,J,K) + W(I+1,J,K)) * 0.5
    WET = (-W(I,J-1,K) + W(I,J+1,K)) * 0.5
    WZE = (-W(I,J,K-1) + W(I,J,K+1)) * 0.5
    ! 物理空間方向一階微分
    UX  = XIX(I,J,K) * UXI + ETX(I,J,K) * UET + ZEX(I,J,K) * UZE
    VX  = XIX(I,J,K) * VXI + ETX(I,J,K) * VET + ZEX(I,J,K) * VZE
    WX  = XIX(I,J,K) * WXI + ETX(I,J,K) * WET + ZEX(I,J,K) * WZE
    UY  = XIY(I,J,K) * UXI + ETY(I,J,K) * UET + ZEY(I,J,K) * UZE
    VY  = XIY(I,J,K) * VXI + ETY(I,J,K) * VET + ZEY(I,J,K) * VZE
    WY  = XIY(I,J,K) * WXI + ETY(I,J,K) * WET + ZEY(I,J,K) * WZE
    UZ  = XIZ(I,J,K) * UXI + ETZ(I,J,K) * UET + ZEZ(I,J,K) * UZE
    VZ  = XIZ(I,J,K) * VXI + ETZ(I,J,K) * VET + ZEZ(I,J,K) * VZE
    WZ  = XIZ(I,J,K) * WXI + ETZ(I,J,K) * WET + ZEZ(I,J,K) * WZE
    ! 速度の二階微分の計算 ---------------------------------------------
    UXX = - (UXXI(I-1,J,K) - UXXI(I,J,K)) * XIX(I,J,K) &
    &     - (UXET(I,J-1,K) - UXET(I,J,K)) * ETX(I,J,K) &
    &     - (UXZE(I,J,K-1) - UXZE(I,J,K)) * ZEX(I,J,K)
    UXY = - (UXXI(I-1,J,K) - UXXI(I,J,K)) * XIY(I,J,K) &
    &     - (UXET(I,J-1,K) - UXET(I,J,K)) * ETY(I,J,K) &
    &     - (UXZE(I,J,K-1) - UXZE(I,J,K)) * ZEY(I,J,K)
    UXZ = - (UXXI(I-1,J,K) - UXXI(I,J,K)) * XIZ(I,J,K) &
    &     - (UXET(I,J-1,K) - UXET(I,J,K)) * ETZ(I,J,K) &
    &     - (UXZE(I,J,K-1) - UXZE(I,J,K)) * ZEZ(I,J,K)
    UYX = - (UYXI(I-1,J,K) - UYXI(I,J,K)) * XIX(I,J,K) &
    &     - (UYET(I,J-1,K) - UYET(I,J,K)) * ETX(I,J,K) &
    &     - (UYZE(I,J,K-1) - UYZE(I,J,K)) * ZEX(I,J,K)
    UYY = - (UYXI(I-1,J,K) - UYXI(I,J,K)) * XIY(I,J,K) &
    &     - (UYET(I,J-1,K) - UYET(I,J,K)) * ETY(I,J,K) &
    &     - (UYZE(I,J,K-1) - UYZE(I,J,K)) * ZEY(I,J,K)
    UYZ = - (UYXI(I-1,J,K) - UYXI(I,J,K)) * XIZ(I,J,K) &
    &     - (UYET(I,J-1,K) - UYET(I,J,K)) * ETZ(I,J,K) &
    &     - (UYZE(I,J,K-1) - UYZE(I,J,K)) * ZEZ(I,J,K)
    UZX = - (UZXI(I-1,J,K) - UZXI(I,J,K)) * XIX(I,J,K) &
    &     - (UZET(I,J-1,K) - UZET(I,J,K)) * ETX(I,J,K) &
    &     - (UZZE(I,J,K-1) - UZZE(I,J,K)) * ZEX(I,J,K)
    UZY = - (UZXI(I-1,J,K) - UZXI(I,J,K)) * XIY(I,J,K) &
    &     - (UZET(I,J-1,K) - UZET(I,J,K)) * ETY(I,J,K) &
    &     - (UZZE(I,J,K-1) - UZZE(I,J,K)) * ZEY(I,J,K)
    UZZ = - (UZXI(I-1,J,K) - UZXI(I,J,K)) * XIZ(I,J,K) &
    &     - (UZET(I,J-1,K) - UZET(I,J,K)) * ETZ(I,J,K) &
    &     - (UZZE(I,J,K-1) - UZZE(I,J,K)) * ZEZ(I,J,K)
    VXX = - (VXXI(I-1,J,K) - VXXI(I,J,K)) * XIX(I,J,K) &
    &     - (VXET(I,J-1,K) - VXET(I,J,K)) * ETX(I,J,K) &
    &     - (VXZE(I,J,K-1) - VXZE(I,J,K)) * ZEX(I,J,K)
    VXY = - (VXXI(I-1,J,K) - VXXI(I,J,K)) * XIY(I,J,K) &
    &     - (VXET(I,J-1,K) - VXET(I,J,K)) * ETY(I,J,K) &
    &     - (VXZE(I,J,K-1) - VXZE(I,J,K)) * ZEY(I,J,K)
    VXZ = - (VXXI(I-1,J,K) - VXXI(I,J,K)) * XIZ(I,J,K) &
    &     - (VXET(I,J-1,K) - VXET(I,J,K)) * ETZ(I,J,K) &
    &     - (VXZE(I,J,K-1) - VXZE(I,J,K)) * ZEZ(I,J,K)
    VYX = - (VYXI(I-1,J,K) - VYXI(I,J,K)) * XIX(I,J,K) &
    &     - (VYET(I,J-1,K) - VYET(I,J,K)) * ETX(I,J,K) &
    &     - (VYZE(I,J,K-1) - VYZE(I,J,K)) * ZEX(I,J,K)
    VYY = - (VYXI(I-1,J,K) - VYXI(I,J,K)) * XIY(I,J,K) &
    &     - (VYET(I,J-1,K) - VYET(I,J,K)) * ETY(I,J,K) &
    &     - (VYZE(I,J,K-1) - VYZE(I,J,K)) * ZEY(I,J,K)
    VYZ = - (VYXI(I-1,J,K) - VYXI(I,J,K)) * XIZ(I,J,K) &
    &     - (VYET(I,J-1,K) - VYET(I,J,K)) * ETZ(I,J,K) &
    &     - (VYZE(I,J,K-1) - VYZE(I,J,K)) * ZEZ(I,J,K)
    VZX = - (VZXI(I-1,J,K) - VZXI(I,J,K)) * XIX(I,J,K) &
    &     - (VZET(I,J-1,K) - VZET(I,J,K)) * ETX(I,J,K) &
    &     - (VZZE(I,J,K-1) - VZZE(I,J,K)) * ZEX(I,J,K)
    VZY = - (VZXI(I-1,J,K) - VZXI(I,J,K)) * XIY(I,J,K) &
    &     - (VZET(I,J-1,K) - VZET(I,J,K)) * ETY(I,J,K) &
    &     - (VZZE(I,J,K-1) - VZZE(I,J,K)) * ZEY(I,J,K)
    VZZ = - (VZXI(I-1,J,K) - VZXI(I,J,K)) * XIZ(I,J,K) &
    &     - (VZET(I,J-1,K) - VZET(I,J,K)) * ETZ(I,J,K) &
    &     - (VZZE(I,J,K-1) - VZZE(I,J,K)) * ZEZ(I,J,K)
    WXX = - (WXXI(I-1,J,K) - WXXI(I,J,K)) * XIX(I,J,K) &
    &     - (WXET(I,J-1,K) - WXET(I,J,K)) * ETX(I,J,K) &
    &     - (WXZE(I,J,K-1) - WXZE(I,J,K)) * ZEX(I,J,K)
    WXY = - (WXXI(I-1,J,K) - WXXI(I,J,K)) * XIY(I,J,K) &
    &     - (WXET(I,J-1,K) - WXET(I,J,K)) * ETY(I,J,K) &
    &     - (WXZE(I,J,K-1) - WXZE(I,J,K)) * ZEY(I,J,K)
    WXZ = - (WXXI(I-1,J,K) - WXXI(I,J,K)) * XIZ(I,J,K) &
    &     - (WXET(I,J-1,K) - WXET(I,J,K)) * ETZ(I,J,K) &
    &     - (WXZE(I,J,K-1) - WXZE(I,J,K)) * ZEZ(I,J,K)
    WYX = - (WYXI(I-1,J,K) - WYXI(I,J,K)) * XIX(I,J,K) &
    &     - (WYET(I,J-1,K) - WYET(I,J,K)) * ETX(I,J,K) &
    &     - (WYZE(I,J,K-1) - WYZE(I,J,K)) * ZEX(I,J,K)
    WYY = - (WYXI(I-1,J,K) - WYXI(I,J,K)) * XIY(I,J,K) &
    &     - (WYET(I,J-1,K) - WYET(I,J,K)) * ETY(I,J,K) &
    &     - (WYZE(I,J,K-1) - WYZE(I,J,K)) * ZEY(I,J,K)
    WYZ = - (WYXI(I-1,J,K) - WYXI(I,J,K)) * XIZ(I,J,K) &
    &     - (WYET(I,J-1,K) - WYET(I,J,K)) * ETZ(I,J,K) &
    &     - (WYZE(I,J,K-1) - WYZE(I,J,K)) * ZEZ(I,J,K)
    WZX = - (WZXI(I-1,J,K) - WZXI(I,J,K)) * XIX(I,J,K) &
    &     - (WZET(I,J-1,K) - WZET(I,J,K)) * ETX(I,J,K) &
    &     - (WZZE(I,J,K-1) - WZZE(I,J,K)) * ZEX(I,J,K)
    WZY = - (WZXI(I-1,J,K) - WZXI(I,J,K)) * XIY(I,J,K) &
    &     - (WZET(I,J-1,K) - WZET(I,J,K)) * ETY(I,J,K) &
    &     - (WZZE(I,J,K-1) - WZZE(I,J,K)) * ZEY(I,J,K)
    WZZ = - (WZXI(I-1,J,K) - WZXI(I,J,K)) * XIZ(I,J,K) &
    &     - (WZET(I,J-1,K) - WZET(I,J,K)) * ETZ(I,J,K) &
    &     - (WZZE(I,J,K-1) - WZZE(I,J,K)) * ZEZ(I,J,K)
    ! 乱れの速度スケール sqrt(k) の一階微分の計算 ----------------------
    ! 計算空間方向一階微分
    KRXI = 0.5 * (- SQRT(AK(I-1,J,K)) + SQRT(AK(I+1,J,K)))
    KRET = 0.5 * (- SQRT(AK(I,J-1,K)) + SQRT(AK(I,J+1,K)))
    KRZE = 0.5 * (- SQRT(AK(I,J,K-1)) + SQRT(AK(I,J,K+1)))
    ! 物理空間方向一階微分
    KRX = KRXI * XIX(I,J,K) + KRET * ETX(I,J,K) + KRZE * ZEX(I,J,K)
    KRY = KRXI * XIY(I,J,K) + KRET * ETY(I,J,K) + KRZE * ZEY(I,J,K)
    KRZ = KRXI * XIZ(I,J,K) + KRET * ETZ(I,J,K) + KRZE * ZEZ(I,J,K)
    ! 歪み速度と渦度 ---------------------------------------------------
    SS  = (UX + VY + WZ) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 = 2.0 * WZ - SS
    S12 = UY + VX
    S23 = VZ + WY
    S31 = WX + UZ
    SS2 = S11 * S11 + S22 * S22 + S33 * S33 &
    &   + 2.0 * (S12 * S12 + S23 * S23 + S31 * S31)
    O12 = UY - VX
    O23 = VZ - WY
    O31 = WX - UZ
    OO2 = 2.0 * (O12 * O12 + O23 * O23 + O31 * O31)
    IF(EPST(I,J,K) .GT. ZERO) THEN
      AKAE = AK(I,J,K) / EPST(I,J,K)
    ELSE
      AKAE = 0.0
    ENDIF
    ST = AKAE * SQRT(0.5 * SS2)
    OT = AKAE * SQRT(0.5 * OO2)
    ! 各モデル係数の計算 -----------------------------------------------
    nu = AMU(I,J,K) / RHO(I,J,K)
    IF(EPST(I,J,K) .GT. ZERO) THEN
      Rtt = AK(I,J,K)**2 / (nu * EPST(I,J,K))
    ELSE
      Rtt = 0.0
    ENDIF
    eta  = MAX(ST, OT)
    c_mu = 0.3 / (1.0 + 0.35 * eta**1.5) &
    &    * (1.0 - EXP(-0.36 / MAX(ZERO, EXP(-0.75 * eta))))
    f_mu = 1.0 - EXP(-SQRT(Rtt / 90.0) - (Rtt / 400.0)**2)
    c_e2 = 1.92 * (1.0 - 0.3 * EXP(-Rtt**2))
    c4   =-10.0 * c_mu**2
    c6   =- 5.0 * c_mu**2
    c7   =  5.0 * c_mu**2
    IF(EPST(I,J,K) .GT. ZERO) THEN
      nu_t = c_mu * f_mu * AK(I,J,K)**2 / EPST(I,J,K)
    ELSE
      nu_t = 0.0
    ENDIF
    AMUT(I,J,K) = RHO(I,J,K) * nu_t
    ! 非等方レイノルズ応力テンソル a_ij --------------------------------
    IF(AK(I,J,K) .GT. ZERO .OR. EPST(I,J,K)**2 .GT. ZERO) THEN
      nu_t1 =-nu_t / AK(I,J,K)
      nu_t2 = nu_t / EPST(I,J,K)
      nu_t3 = nu_t * AK(I,J,K) / EPST(I,J,K)**2
    ELSE
      nu_t1 = 0.0
      nu_t2 = 0.0
      nu_t3 = 0.0
    ENDIF
    a11 = nu_t1 * S11 &
    &   + nu_t2 * ( &
    &       c1 * (S11 * S11 + S31 * S31 + S12 * S12 - SS2 / 3.0) &
    &     + c2 * (O12 * S12 - O31 * S31) * 2.0 &
    &     + c3 * (O12 * O12 + O31 * O31 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( (S12 * S23 + (S33 + S11) * S31) * O31 &
    &            - (S23 * S31 + (S11 + S22) * S12) * O12 ) * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S11 &
    &   )
    a22 = nu_t1 * S22 &
    &   + nu_t2 * ( &
    &       c1 * (S22 * S22 + S12 * S12 + S23 * S23 - SS2 / 3.0) &
    &     + c2 * (O23 * S23 - O12 * S12) * 2.0 &
    &     + c3 * (O23 * O23 + O12 * O12 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( (S23 * S31 + (S11 + S22) * S12) * O12 &
    &            - (S31 * S12 + (S22 + S33) * S23) * O23 ) * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S22 &
    &   )
    a33 = nu_t1 * S33 &
    &   + nu_t2 * ( &
    &       c1 * (S33 * S33 + S23 * S23 + S31 * S31 - SS2 / 3.0) &
    &     + c2 * (O31 * S31 - O23 * S23) * 2.0 &
    &     + c3 * (O31 * O31 + O23 * O23 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( (S31 * S12 + (S22 + S33) * S23) * O23 &
    &            - (S12 * S23 + (S33 + S11) * S31) * O31) * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S33 &
    &   )
    a12 = nu_t1 * S12 &
    &   + nu_t2 * ( &
    &       c1 * ((S11 + S22) * S12 + S23 * S31) &
    &     + c2 * (O23 * S31 - O31 * S23 - O12 * (S11 - S22)) &
    &     - c3 * (O23 * O31) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( &
    &         (S11 * S11 - S22 * S22 - S23 * S23 + S31 * S31) * O12 &
    &       + ((S22 + S33) * S23 + S31 * S12) * O31 &
    &       - ((S33 + S11) * S31 + S12 * S23) * O23 &
    &     ) &
    &     + (c6 * SS2 + c7 * OO2) * S12 &
    &   )
    a23 = nu_t1 * S23 &
    &   + nu_t2 * ( &
    &       c1 * ((S22 + S33) * S23 + S31 * S12) &
    &     + c2 * (O31 * S12 - O12 * S31 - O23 * (S22 - S33)) &
    &     - c3 * (O31 * O12) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( &
    &         (S22 * S22 - S33 * S33 - S31 * S31 + S12 * S12) * O23 &
    &       + ((S11 + S33) * S31 + S12 * S23) * O12 &
    &       - ((S11 + S22) * S12 + S23 * S31) * O31 &
    &     ) &
    &     + (c6 * SS2 + c7 * OO2) * S23 &
    &   )
    a31 = nu_t1 * S31 &
    &   + nu_t2 * ( &
    &       c1 * ((S11 + S33) * S31 + S12 * S23) &
    &     + c2 * (O12 * S23 - O23 * S12 - O31 * (S33 - S11)) &
    &     - c3 * (O12 * O23) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( &
    &         (S33 * S33 - S11 * S11 - S12 * S12 + S23 * S23) * O31 &
    &       + ((S11 + S22) * S12 + S31 * S23) * O23 &
    &       - ((S22 + S33) * S23 + S12 * S31) * O12 &
    &     ) &
    &     + (c6 * SS2 + c7 * OO2) * S31 &
    &   )
    a11 = MIN(aiimax, MAX(aiimin, a11))
    a22 = MIN(aiimax, MAX(aiimin, a22))
    a33 = MIN(aiimax, MAX(aiimin, a33))
    a12 = MIN(1.0, MAX(-1.0, a12))
    a23 = MIN(1.0, MAX(-1.0, a23))
    a31 = MIN(1.0, MAX(-1.0, a31))
    ! レイノルズ応力 ---------------------------------------------------
    uu = AK(I,J,K) * (a11 + TwoThird)
    vv = AK(I,J,K) * (a22 + TwoThird)
    ww = AK(I,J,K) * (a33 + TwoThird)
    uv = AK(I,J,K) * a12
    vw = AK(I,J,K) * a23
    wu = AK(I,J,K) * a31
    ! 低レイノルズ補正項 D ---------------------------------------------
    D = 2.0 * nu * (KRX**2 + KRY**2 + KRZ**2)
    ! 低レイノルズ補正項 E ---------------------------------------------
    IF(Rtt .LE. 250.0) THEN
      E = 0.0022 * ST * nu_t * AK(I,J,K) * AKAE * ( &
      &     UXX**2 + VXX**2 + WXX**2 &
      &   + UYY**2 + VYY**2 + WYY**2 &
      &   + UZZ**2 + VZZ**2 + WZZ**2 &
      &   + UXY**2 + VXY**2 + WXY**2 &
      &   + UYZ**2 + VYZ**2 + WYZ**2 &
      &   + UZX**2 + VZX**2 + WZX**2 &
      &   + UXZ**2 + VXZ**2 + WXZ**2 &
      &   + UYX**2 + VYX**2 + WYX**2 &
      &   + UZY**2 + VZY**2 + WZY**2 &
      & )
    ELSE
      E = 0.0
    ENDIF
    ! レイノルズ応力の生成(平均歪み増加率) -----------------------------
    PP_k = - uu * UX - vv * VY - ww * WZ &
    &      - uv * (VX + UY) - vw * (WY + VZ) - wu * (WX + UZ)
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      PP_k = MAX(0.0, MIN((EPST(I,J,K) + D) / PELIM, PP_k))
    ENDIF
    ! ソース項の和 -----------------------------------------------------
    QH1 = RHO(I,J,K) / AJA(I,J,K)
    AEAK = MIN( &
    &      MAX(SQRT(0.5 * SS2), SQRT((EPST(I,J,K) + D) / nu)), &
    &      EPST(I,J,K) / MAX(ZERO, AK(I,J,K)) &
    &    )
    DQP(I,J,K,6) = QH1 * (PP_k - (EPST(I,J,K) + D))
    DQP(I,J,K,7) = QH1 * (AEAK * (c_e1 * PP_k - c_e2 * EPST(I,J,K)) + E)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH = 0.0
  SH = 0.0
  TH = 0.0
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL DIFFXI
  CALL DIFFET
  CALL DIFFZE
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I,J-1,K) .GT. 0.0) .AND. &
  &   (RHO(I-1,J,K) .GT. 0.0) .AND. (RHO(I,J,K)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J,K) .GT. 0.0) .AND. (RHO(I,J+1,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,K,L) = (-RH(I-1,J  ,K  ,L) + RH(I,J,K,L) ) &
    &            + (-SH(I  ,J-1,K  ,L) + SH(I,J,K,L) ) &
    &            + (-TH(I  ,J  ,K-1,L) + TH(I,J,K,L) )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項の計算(xi方向差分用)                                    ****
!***********************************************************************
SUBROUTINE DIFFXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, AKM, EPSTM, AMUM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: AKXI, AKET, AKZE, AKX, AKY, AKZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: S11, S22, S33, S12, S23, S31, SS2, ST, SS
  REAL    :: O12, O23, O31, OO2, OT
  REAL    :: eta
  REAL    :: nu, nu_t
  REAL    :: Rtt, c_mu, f_mu
  REAL    :: c4, c6, c7
  REAL    :: nu_t1, nu_t2, nu_t3, AKAE
  REAL    :: a11, a22, a33, a12, a23, a31
  REAL    :: uu, vv, ww, uv, vw, wu
  REAL    :: c_t, c_k, c_e
  REAL    :: TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX
  REAL    :: R5, S5, T5
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, AKM, EPSTM, AMUM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  AKXI, AKET, AKZE, AKX, AKY, AKZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  S11, S22, S33, S12, S23, S31, SS2, ST, SS, &
  !$  O12, O23, O31, OO2, OT, &
  !$  eta, nu, nu_t, Rtt, c_mu, f_mu, c4, c6, c7, &
  !$  nu_t1, nu_t2, nu_t3, AKAE, &
  !$  a11, a22, a33, a12, a23, a31, uu, vv, ww, uv, vw, wu, &
  !$  c_t, c_k, c_e, &
  !$  TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, R5, S5, T5, QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF((RHO(I,J,K) .GT. 0.0) .AND. (RHO(I+1,J,K) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I+1,J,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I+1,J,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I+1,J,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I+1,J,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I+1,J,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I+1,J,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I+1,J,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I+1,J,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I+1,J,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I+1,J,K))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J,K) +  RHO(I+1,J,K))
    UM    = 0.5 * (   U(I,J,K) +    U(I+1,J,K))
    VM    = 0.5 * (   V(I,J,K) +    V(I+1,J,K))
    WM    = 0.5 * (   W(I,J,K) +    W(I+1,J,K))
    AKM   = 0.5 * (  AK(I,J,K) +   AK(I+1,J,K))
    EPSTM = 0.5 * (EPST(I,J,K) + EPST(I+1,J,K))
    AMUM  = 0.5 * ( AMU(I,J,K) +  AMU(I+1,J,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = - U(I,J,K) + U(I+1,J,K)
    UET = 0.25 * ( - U(I,J-1,K) - U(I+1,J-1,K) &
    &              + U(I,J+1,K) + U(I+1,J+1,K) )
    UZE = 0.25 * ( - U(I,J,K-1) - U(I+1,J,K-1) &
    &              + U(I,J,K+1) + U(I+1,J,K+1) )
    VXI = - V(I,J,K) + V(I+1,J,K)
    VET = 0.25 * ( - V(I,J-1,K) - V(I+1,J-1,K) &
    &              + V(I,J+1,K) + V(I+1,J+1,K) )
    VZE = 0.25 * ( - V(I,J,K-1) - V(I+1,J,K-1) &
    &              + V(I,J,K+1) + V(I+1,J,K+1) )
    WXI = - W(I,J,K) + W(I+1,J,K)
    WET = 0.25 * ( - W(I,J-1,K) - W(I+1,J-1,K) &
    &              + W(I,J+1,K) + W(I+1,J+1,K) )
    WZE = 0.25 * ( - W(I,J,K-1) - W(I+1,J,K-1) &
    &              + W(I,J,K+1) + W(I+1,J,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = - T(I,J,K) + T(I+1,J,K)
    TET = 0.25 * ( - T(I,J-1,K) - T(I+1,J-1,K) &
    &              + T(I,J+1,K) + T(I+1,J+1,K) )
    TZE = 0.25 * ( - T(I,J,K-1) - T(I+1,J,K-1) &
    &              + T(I,J,K+1) + T(I+1,J,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! k の一階微分の計算 -----------------------------------------------
    ! 計算空間方向一階微分
    AKXI = - AK(I,J,K) + AK(I+1,J,K)
    AKET = 0.25 * ( - AK(I,J-1,K) - AK(I+1,J-1,K) &
    &               + AK(I,J+1,K) + AK(I+1,J+1,K) )
    AKZE = 0.25 * ( - AK(I,J,K-1) - AK(I+1,J,K-1) &
    &               + AK(I,J,K+1) + AK(I+1,J,K+1) )
    ! 物理空間方向一階微分
    AKX = AKXI * XIXM + AKET * ETXM + AKZE * ZEXM
    AKY = AKXI * XIYM + AKET * ETYM + AKZE * ZEYM
    AKZ = AKXI * XIZM + AKET * ETZM + AKZE * ZEZM
    ! epsilon tilde の一階微分の計算 -----------------------------------
    ! 計算空間方向一階微分
    AEXI = - EPST(I,J,K) + EPST(I+1,J,K)
    AEET = 0.25 * ( - EPST(I,J-1,K) - EPST(I+1,J-1,K) &
    &               + EPST(I,J+1,K) + EPST(I+1,J+1,K) )
    AEZE = 0.25 * ( - EPST(I,J,K-1) - EPST(I+1,J,K-1) &
    &               + EPST(I,J,K+1) + EPST(I+1,J,K+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 歪み速度と渦度 ---------------------------------------------------
    SS  = (UX + VY + WZ) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 = 2.0 * WZ - SS
    S12 = UY + VX
    S23 = VZ + WY
    S31 = WX + UZ
    SS2 = S11 * S11 + S22 * S22 + S33 * S33 &
    &   + 2.0 * (S12 * S12 + S23 * S23 + S31 * S31)
    O12 = UY - VX
    O23 = VZ - WY
    O31 = WX - UZ
    OO2 = 2.0 * (O12 * O12 + O23 * O23 + O31 * O31)
    IF(EPSTM .GT. ZERO) THEN
      AKAE = AKM / EPSTM
    ELSE
      AKAE = 0.0
    ENDIF
    ST = AKAE * SQRT(0.5 * SS2)
    OT = AKAE * SQRT(0.5 * OO2)
    ! 非等方レイノルズ応力テンソル a_ij の係数の計算 -------------------
    nu = AMUM / RHOM
    IF(EPSTM .GT. ZERO) THEN
      Rtt = AKM**2 / (nu * EPSTM)
    ELSE
      Rtt = 0.0
    ENDIF
    eta  = MAX(ST, OT)
    c_mu = 0.3 / (1.0 + 0.35 * eta**1.5) &
    &    * (1.0 - EXP(-0.36 / MAX(ZERO, EXP(-0.75 * eta))))
    f_mu = 1.0 - EXP(-SQRT(Rtt / 90.0) - (Rtt / 400.0)**2)
    c4   =-10.0 * c_mu**2
    c6   =- 5.0 * c_mu**2
    c7   =  5.0 * c_mu**2
    IF(EPSTM .GT. ZERO) THEN
      nu_t = c_mu * f_mu * AKM**2 / EPSTM
    ELSE
      nu_t = 0.0
    ENDIF
    ! 非等方レイノルズ応力テンソル a_ij --------------------------------
    IF(AKM .GT. ZERO .OR. EPSTM**2 .GT. ZERO) THEN
      nu_t1 =-nu_t / AKM
      nu_t2 = nu_t / EPSTM
      nu_t3 = nu_t * AKM / EPSTM**2
    ELSE
      nu_t1 = 0.0
      nu_t2 = 0.0
      nu_t3 = 0.0
    ENDIF
    a11 = nu_t1 * S11 &
    &   + nu_t2 * ( &
    &       c1 * (S11 * S11 + S31 * S31 + S12 * S12 - SS2 / 3.0) &
    &     + c2 * (O12 * S12 - O31 * S31) * 2.0 &
    &     + c3 * (O12 * O12 + O31 * O31 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( (S12 * S23 + (S33 + S11) * S31) * O31 &
    &            - (S23 * S31 + (S11 + S22) * S12) * O12 ) * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S11 &
    &   )
    a22 = nu_t1 * S22 &
    &   + nu_t2 * ( &
    &       c1 * (S22 * S22 + S12 * S12 + S23 * S23 - SS2 / 3.0) &
    &     + c2 * (O23 * S23 - O12 * S12) * 2.0 &
    &     + c3 * (O23 * O23 + O12 * O12 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( (S23 * S31 + (S11 + S22) * S12) * O12 &
    &            - (S31 * S12 + (S22 + S33) * S23) * O23 ) * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S22 &
    &   )
    a33 = nu_t1 * S33 &
    &   + nu_t2 * ( &
    &       c1 * (S33 * S33 + S23 * S23 + S31 * S31 - SS2 / 3.0) &
    &     + c2 * (O31 * S31 - O23 * S23) * 2.0 &
    &     + c3 * (O31 * O31 + O23 * O23 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( (S31 * S12 + (S22 + S33) * S23) * O23 &
    &            - (S12 * S23 + (S33 + S11) * S31) * O31) * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S33 &
    &   )
    a12 = nu_t1 * S12 &
    &   + nu_t2 * ( &
    &       c1 * ((S11 + S22) * S12 + S23 * S31) &
    &     + c2 * (O23 * S31 - O31 * S23 - O12 * (S11 - S22)) &
    &     - c3 * (O23 * O31) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( &
    &         (S11 * S11 - S22 * S22 - S23 * S23 + S31 * S31) * O12 &
    &       + ((S22 + S33) * S23 + S31 * S12) * O31 &
    &       - ((S33 + S11) * S31 + S12 * S23) * O23 &
    &     ) &
    &     + (c6 * SS2 + c7 * OO2) * S12 &
    &   )
    a23 = nu_t1 * S23 &
    &   + nu_t2 * ( &
    &       c1 * ((S22 + S33) * S23 + S31 * S12) &
    &     + c2 * (O31 * S12 - O12 * S31 - O23 * (S22 - S33)) &
    &     - c3 * (O31 * O12) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( &
    &         (S22 * S22 - S33 * S33 - S31 * S31 + S12 * S12) * O23 &
    &       + ((S11 + S33) * S31 + S12 * S23) * O12 &
    &       - ((S11 + S22) * S12 + S23 * S31) * O31 &
    &     ) &
    &     + (c6 * SS2 + c7 * OO2) * S23 &
    &   )
    a31 = nu_t1 * S31 &
    &   + nu_t2 * ( &
    &       c1 * ((S11 + S33) * S31 + S12 * S23) &
    &     + c2 * (O12 * S23 - O23 * S12 - O31 * (S33 - S11)) &
    &     - c3 * (O12 * O23) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( &
    &         (S33 * S33 - S11 * S11 - S12 * S12 + S23 * S23) * O31 &
    &       + ((S11 + S22) * S12 + S31 * S23) * O23 &
    &       - ((S22 + S33) * S23 + S12 * S31) * O12 &
    &     ) &
    &     + (c6 * SS2 + c7 * OO2) * S31 &
    &   )
    a11 = MIN(aiimax, MAX(aiimin, a11))
    a22 = MIN(aiimax, MAX(aiimin, a22))
    a33 = MIN(aiimax, MAX(aiimin, a33))
    a12 = MIN(1.0, MAX(-1.0, a12))
    a23 = MIN(1.0, MAX(-1.0, a23))
    a31 = MIN(1.0, MAX(-1.0, a31))
    ! レイノルズ応力 ---------------------------------------------------
    uu = AKM * (a11 + TwoThird)
    vv = AKM * (a22 + TwoThird)
    ww = AKM * (a33 + TwoThird)
    uv = AKM * a12
    vw = AKM * a23
    wu = AKM * a31
    ! 拡散係数 ---------------------------------------------------------
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    c_k  = nu + nu_t / SIGK
    c_e  = nu + nu_t / SIGE
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    TAUXX = nu * S11 - uu
    TAUYY = nu * S22 - vv
    TAUZZ = nu * S33 - ww
    TAUXY = nu * S12 - uv
    TAUYZ = nu * S23 - vw
    TAUZX = nu * S31 - wu
    ! エネルギーの拡散 -------------------------------------------------
    R5 = TAUXX * UM + TAUXY * VM + TAUZX * WM + c_t * TX
    S5 = TAUXY * UM + TAUYY * VM + TAUYZ * WM + c_t * TY
    T5 = TAUZX * UM + TAUYZ * VM + TAUZZ * WM + c_t * TZ
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    RH(I,J,K,1) = 0.0
    RH(I,J,K,2) = QH1 * (TAUXX * XIXM + TAUXY * XIYM + TAUZX * XIZM)
    RH(I,J,K,3) = QH1 * (TAUXY * XIXM + TAUYY * XIYM + TAUYZ * XIZM)
    RH(I,J,K,4) = QH1 * (TAUZX * XIXM + TAUYZ * XIYM + TAUZZ * XIZM)
    RH(I,J,K,5) = QH1 * (R5 * XIXM + S5 * XIYM + T5 * XIZM)
    RH(I,J,K,6) = QH1 * c_k * (AKX * XIXM + AKY * XIYM + AKZ * XIZM)
    RH(I,J,K,7) = QH1 * c_e * (AEX * XIXM + AEY * XIYM + AEZ * XIZM)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFXI
!***********************************************************************
!**** 拡散項の計算(eta方向差分用)                                   ****
!***********************************************************************
SUBROUTINE DIFFET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, AKM, EPSTM, AMUM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: AKXI, AKET, AKZE, AKX, AKY, AKZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: S11, S22, S33, S12, S23, S31, SS2, ST, SS
  REAL    :: O12, O23, O31, OO2, OT
  REAL    :: eta
  REAL    :: nu, nu_t
  REAL    :: Rtt, c_mu, f_mu
  REAL    :: c4, c6, c7
  REAL    :: nu_t1, nu_t2, nu_t3, AKAE
  REAL    :: a11, a22, a33, a12, a23, a31
  REAL    :: uu, vv, ww, uv, vw, wu
  REAL    :: c_t, c_k, c_e
  REAL    :: TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX
  REAL    :: R5, S5, T5
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, AKM, EPSTM, AMUM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  AKXI, AKET, AKZE, AKX, AKY, AKZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  S11, S22, S33, S12, S23, S31, SS2, ST, SS, &
  !$  O12, O23, O31, OO2, OT, &
  !$  eta, nu, nu_t, Rtt, c_mu, f_mu, c4, c6, c7, &
  !$  nu_t1, nu_t2, nu_t3, AKAE, &
  !$  a11, a22, a33, a12, a23, a31, uu, vv, ww, uv, vw, wu, &
  !$  c_t, c_k, c_e, &
  !$  TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, R5, S5, T5, QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF((RHO(I,J,K) .GT. 0.0) .AND. (RHO(I,J+1,K) .GT. 0.0)) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J+1,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J+1,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J+1,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J+1,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J+1,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J+1,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J+1,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J+1,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J+1,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J+1,K))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J,K) +  RHO(I,J+1,K))
    UM    = 0.5 * (   U(I,J,K) +    U(I,J+1,K))
    VM    = 0.5 * (   V(I,J,K) +    V(I,J+1,K))
    WM    = 0.5 * (   W(I,J,K) +    W(I,J+1,K))
    AKM   = 0.5 * (  AK(I,J,K) +   AK(I,J+1,K))
    EPSTM = 0.5 * (EPST(I,J,K) + EPST(I,J+1,K))
    AMUM  = 0.5 * ( AMU(I,J,K) +  AMU(I,J+1,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J+1,K) &
    &              + U(I+1,J,K) + U(I+1,J+1,K) )
    UET = - U(I,J,K) + U(I,J+1,K)
    UZE = 0.25 * ( - U(I,J,K-1) - U(I,J+1,K-1) &
    &              + U(I,J,K+1) + U(I,J+1,K+1) )
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J+1,K) &
    &              + V(I+1,J,K) + V(I+1,J+1,K) )
    VET = - V(I,J,K) + V(I,J+1,K)
    VZE = 0.25 * ( - V(I,J,K-1) - V(I,J+1,K-1) &
    &              + V(I,J,K+1) + V(I,J+1,K+1) )
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J+1,K) &
    &              + W(I+1,J,K) + W(I+1,J+1,K) )
    WET = - W(I,J,K) + W(I,J+1,K)
    WZE = 0.25 * ( - W(I,J,K-1) - W(I,J+1,K-1) &
    &              + W(I,J,K+1) + W(I,J+1,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J+1,K) &
    &              + T(I+1,J,K) + T(I+1,J+1,K) )
    TET = - T(I,J,K) + T(I,J+1,K)
    TZE = 0.25 * ( - T(I,J,K-1) - T(I,J+1,K-1) &
    &              + T(I,J,K+1) + T(I,J+1,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! k の一階微分の計算 -----------------------------------------------
    ! 計算空間方向一階微分
    AKXI = 0.25 * ( - AK(I-1,J,K) - AK(I-1,J+1,K) &
    &               + AK(I+1,J,K) + AK(I+1,J+1,K) )
    AKET = - AK(I,J,K) + AK(I,J+1,K)
    AKZE = 0.25 * ( - AK(I,J,K-1) - AK(I,J+1,K-1) &
    &               + AK(I,J,K+1) + AK(I,J+1,K+1) )
    ! 物理空間方向一階微分
    AKX = AKXI * XIXM + AKET * ETXM + AKZE * ZEXM
    AKY = AKXI * XIYM + AKET * ETYM + AKZE * ZEYM
    AKZ = AKXI * XIZM + AKET * ETZM + AKZE * ZEZM
    ! epsilon tilde の一階微分の計算 -----------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPST(I-1,J,K) - EPST(I-1,J+1,K) &
    &               + EPST(I+1,J,K) + EPST(I+1,J+1,K) )
    AEET = - EPST(I,J,K) + EPST(I,J+1,K)
    AEZE = 0.25 * ( - EPST(I,J,K-1) - EPST(I,J+1,K-1) &
    &               + EPST(I,J,K+1) + EPST(I,J+1,K+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 歪み速度と渦度 ---------------------------------------------------
    SS  = (UX + VY + WZ) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 = 2.0 * WZ - SS
    S12 = UY + VX
    S23 = VZ + WY
    S31 = WX + UZ
    SS2 = S11 * S11 + S22 * S22 + S33 * S33 &
    &   + 2.0 * (S12 * S12 + S23 * S23 + S31 * S31)
    O12 = UY - VX
    O23 = VZ - WY
    O31 = WX - UZ
    OO2 = 2.0 * (O12 * O12 + O23 * O23 + O31 * O31)
    IF(EPSTM .GT. ZERO) THEN
      AKAE = AKM / EPSTM
    ELSE
      AKAE = 0.0
    ENDIF
    ST = AKAE * SQRT(0.5 * SS2)
    OT = AKAE * SQRT(0.5 * OO2)
    ! 非等方レイノルズ応力テンソル a_ij の係数の計算 -------------------
    nu = AMUM / RHOM
    IF(EPSTM .GT. ZERO) THEN
      Rtt = AKM**2 / (nu * EPSTM)
    ELSE
      Rtt = 0.0
    ENDIF
    eta  = MAX(ST, OT)
    c_mu = 0.3 / (1.0 + 0.35 * eta**1.5) &
    &    * (1.0 - EXP(-0.36 / MAX(ZERO, EXP(-0.75 * eta))))
    f_mu = 1.0 - EXP(-SQRT(Rtt / 90.0) - (Rtt / 400.0)**2)
    c4   =-10.0 * c_mu**2
    c6   =- 5.0 * c_mu**2
    c7   =  5.0 * c_mu**2
    IF(EPSTM .GT. ZERO) THEN
      nu_t = c_mu * f_mu * AKM**2 / EPSTM
    ELSE
      nu_t = 0.0
    ENDIF
    ! 非等方レイノルズ応力テンソル a_ij --------------------------------
    IF(AKM .GT. ZERO .OR. EPSTM**2 .GT. ZERO) THEN
      nu_t1 =-nu_t / AKM
      nu_t2 = nu_t / EPSTM
      nu_t3 = nu_t * AKM / EPSTM**2
    ELSE
      nu_t1 = 0.0
      nu_t2 = 0.0
      nu_t3 = 0.0
    ENDIF
    a11 = nu_t1 * S11 &
    &   + nu_t2 * ( &
    &       c1 * (S11 * S11 + S31 * S31 + S12 * S12 - SS2 / 3.0) &
    &     + c2 * (O12 * S12 - O31 * S31) * 2.0 &
    &     + c3 * (O12 * O12 + O31 * O31 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( (S12 * S23 + (S33 + S11) * S31) * O31 &
    &            - (S23 * S31 + (S11 + S22) * S12) * O12 ) * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S11 &
    &   )
    a22 = nu_t1 * S22 &
    &   + nu_t2 * ( &
    &       c1 * (S22 * S22 + S12 * S12 + S23 * S23 - SS2 / 3.0) &
    &     + c2 * (O23 * S23 - O12 * S12) * 2.0 &
    &     + c3 * (O23 * O23 + O12 * O12 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( (S23 * S31 + (S11 + S22) * S12) * O12 &
    &            - (S31 * S12 + (S22 + S33) * S23) * O23 ) * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S22 &
    &   )
    a33 = nu_t1 * S33 &
    &   + nu_t2 * ( &
    &       c1 * (S33 * S33 + S23 * S23 + S31 * S31 - SS2 / 3.0) &
    &     + c2 * (O31 * S31 - O23 * S23) * 2.0 &
    &     + c3 * (O31 * O31 + O23 * O23 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( (S31 * S12 + (S22 + S33) * S23) * O23 &
    &            - (S12 * S23 + (S33 + S11) * S31) * O31) * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S33 &
    &   )
    a12 = nu_t1 * S12 &
    &   + nu_t2 * ( &
    &       c1 * ((S11 + S22) * S12 + S23 * S31) &
    &     + c2 * (O23 * S31 - O31 * S23 - O12 * (S11 - S22)) &
    &     - c3 * (O23 * O31) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( &
    &         (S11 * S11 - S22 * S22 - S23 * S23 + S31 * S31) * O12 &
    &       + ((S22 + S33) * S23 + S31 * S12) * O31 &
    &       - ((S33 + S11) * S31 + S12 * S23) * O23 &
    &     ) &
    &     + (c6 * SS2 + c7 * OO2) * S12 &
    &   )
    a23 = nu_t1 * S23 &
    &   + nu_t2 * ( &
    &       c1 * ((S22 + S33) * S23 + S31 * S12) &
    &     + c2 * (O31 * S12 - O12 * S31 - O23 * (S22 - S33)) &
    &     - c3 * (O31 * O12) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( &
    &         (S22 * S22 - S33 * S33 - S31 * S31 + S12 * S12) * O23 &
    &       + ((S11 + S33) * S31 + S12 * S23) * O12 &
    &       - ((S11 + S22) * S12 + S23 * S31) * O31 &
    &     ) &
    &     + (c6 * SS2 + c7 * OO2) * S23 &
    &   )
    a31 = nu_t1 * S31 &
    &   + nu_t2 * ( &
    &       c1 * ((S11 + S33) * S31 + S12 * S23) &
    &     + c2 * (O12 * S23 - O23 * S12 - O31 * (S33 - S11)) &
    &     - c3 * (O12 * O23) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( &
    &         (S33 * S33 - S11 * S11 - S12 * S12 + S23 * S23) * O31 &
    &       + ((S11 + S22) * S12 + S31 * S23) * O23 &
    &       - ((S22 + S33) * S23 + S12 * S31) * O12 &
    &     ) &
    &     + (c6 * SS2 + c7 * OO2) * S31 &
    &   )
    a11 = MIN(aiimax, MAX(aiimin, a11))
    a22 = MIN(aiimax, MAX(aiimin, a22))
    a33 = MIN(aiimax, MAX(aiimin, a33))
    a12 = MIN(1.0, MAX(-1.0, a12))
    a23 = MIN(1.0, MAX(-1.0, a23))
    a31 = MIN(1.0, MAX(-1.0, a31))
    ! レイノルズ応力 ---------------------------------------------------
    uu = AKM * (a11 + TwoThird)
    vv = AKM * (a22 + TwoThird)
    ww = AKM * (a33 + TwoThird)
    uv = AKM * a12
    vw = AKM * a23
    wu = AKM * a31
    ! 拡散係数 ---------------------------------------------------------
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    c_k  = nu + nu_t / SIGK
    c_e  = nu + nu_t / SIGE
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    TAUXX = nu * S11 - uu
    TAUYY = nu * S22 - vv
    TAUZZ = nu * S33 - ww
    TAUXY = nu * S12 - uv
    TAUYZ = nu * S23 - vw
    TAUZX = nu * S31 - wu
    ! エネルギーの拡散 -------------------------------------------------
    R5 = TAUXX * UM + TAUXY * VM + TAUZX * WM + c_t * TX
    S5 = TAUXY * UM + TAUYY * VM + TAUYZ * WM + c_t * TY
    T5 = TAUZX * UM + TAUYZ * VM + TAUZZ * WM + c_t * TZ
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    SH(I,J,K,1) = 0.0
    SH(I,J,K,2) = QH1 * (TAUXX * ETXM + TAUXY * ETYM + TAUZX * ETZM)
    SH(I,J,K,3) = QH1 * (TAUXY * ETXM + TAUYY * ETYM + TAUYZ * ETZM)
    SH(I,J,K,4) = QH1 * (TAUZX * ETXM + TAUYZ * ETYM + TAUZZ * ETZM)
    SH(I,J,K,5) = QH1 * (R5 * ETXM + S5 * ETYM + T5 * ETZM)
    SH(I,J,K,6) = QH1 * c_k * (AKX * ETXM + AKY * ETYM + AKZ * ETZM)
    SH(I,J,K,7) = QH1 * c_e * (AEX * ETXM + AEY * ETYM + AEZ * ETZM)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFET
!***********************************************************************
!**** 拡散項の計算(zeta方向差分用)                                  ****
!***********************************************************************
SUBROUTINE DIFFZE
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, AKM, EPSTM, AMUM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: AKXI, AKET, AKZE, AKX, AKY, AKZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: S11, S22, S33, S12, S23, S31, SS2, ST, SS
  REAL    :: O12, O23, O31, OO2, OT
  REAL    :: eta
  REAL    :: nu, nu_t
  REAL    :: Rtt, c_mu, f_mu
  REAL    :: c4, c6, c7
  REAL    :: nu_t1, nu_t2, nu_t3, AKAE
  REAL    :: a11, a22, a33, a12, a23, a31
  REAL    :: uu, vv, ww, uv, vw, wu
  REAL    :: c_t, c_k, c_e
  REAL    :: TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX
  REAL    :: R5, S5, T5
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, AKM, EPSTM, AMUM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  AKXI, AKET, AKZE, AKX, AKY, AKZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  S11, S22, S33, S12, S23, S31, SS2, ST, SS, &
  !$  O12, O23, O31, OO2, OT, &
  !$  eta, nu, nu_t, Rtt, c_mu, f_mu, c4, c6, c7, &
  !$  nu_t1, nu_t2, nu_t3, AKAE, &
  !$  a11, a22, a33, a12, a23, a31, uu, vv, ww, uv, vw, wu, &
  !$  c_t, c_k, c_e, &
  !$  TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, R5, S5, T5, QH1 &
  !$)
  DO K = KS    , KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF((RHO(I,J,K) .GT. 0.0) .AND. (RHO(I,J,K+1) .GT. 0.0)) THEN
    ! ZE方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J,K+1))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J,K+1))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J,K+1))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J,K+1))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J,K+1))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J,K+1))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J,K+1))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J,K+1))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J,K+1))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J,K+1))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J,K) +  RHO(I,J,K+1))
    UM    = 0.5 * (   U(I,J,K) +    U(I,J,K+1))
    VM    = 0.5 * (   V(I,J,K) +    V(I,J,K+1))
    WM    = 0.5 * (   W(I,J,K) +    W(I,J,K+1))
    AKM   = 0.5 * (  AK(I,J,K) +   AK(I,J,K+1))
    EPSTM = 0.5 * (EPST(I,J,K) + EPST(I,J,K+1))
    AMUM  = 0.5 * ( AMU(I,J,K) +  AMU(I,J,K+1))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J,K+1) &
    &              + U(I+1,J,K) + U(I+1,J,K+1) )
    UET = 0.25 * ( - U(I,J-1,K) - U(I,J-1,K+1) &
    &              + U(I,J+1,K) + U(I,J+1,K+1) )
    UZE = - U(I,J,K) + U(I,J,K+1)
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J,K+1) &
    &              + V(I+1,J,K) + V(I+1,J,K+1) )
    VET = 0.25 * ( - V(I,J-1,K) - V(I,J-1,K+1) &
    &              + V(I,J+1,K) + V(I,J+1,K+1) )
    VZE = - V(I,J,K) + V(I,J,K+1)
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J,K+1) &
    &              + W(I+1,J,K) + W(I+1,J,K+1) )
    WET = 0.25 * ( - W(I,J-1,K) - W(I,J-1,K+1) &
    &              + W(I,J+1,K) + W(I,J+1,K+1) )
    WZE = - W(I,J,K) + W(I,J,K+1)
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J,K+1) &
    &              + T(I+1,J,K) + T(I+1,J,K+1) )
    TET = 0.25 * ( - T(I,J-1,K) - T(I,J-1,K+1) &
    &              + T(I,J+1,K) + T(I,J+1,K+1) )
    TZE = - T(I,J,K) + T(I,J,K+1)
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! k の一階微分の計算 -----------------------------------------------
    ! 計算空間方向一階微分
    AKXI = 0.25 * ( - AK(I-1,J,K) - AK(I-1,J,K+1) &
    &               + AK(I+1,J,K) + AK(I+1,J,K+1) )
    AKET = 0.25 * ( - AK(I,J-1,K) - AK(I,J-1,K+1) &
    &               + AK(I,J+1,K) + AK(I,J+1,K+1) )
    AKZE = - AK(I,J,K) + AK(I,J,K+1)
    ! 物理空間方向一階微分
    AKX = AKXI * XIXM + AKET * ETXM + AKZE * ZEXM
    AKY = AKXI * XIYM + AKET * ETYM + AKZE * ZEYM
    AKZ = AKXI * XIZM + AKET * ETZM + AKZE * ZEZM
    ! epsilon tilde の一階微分の計算 -----------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPST(I-1,J,K) - EPST(I-1,J,K+1) &
    &               + EPST(I+1,J,K) + EPST(I+1,J,K+1) )
    AEET = 0.25 * ( - EPST(I,J-1,K) - EPST(I,J-1,K+1) &
    &               + EPST(I,J+1,K) + EPST(I,J+1,K+1) )
    AEZE = - EPST(I,J,K) + EPST(I,J,K+1)
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 歪み速度と渦度 ---------------------------------------------------
    SS  = (UX + VY + WZ) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 = 2.0 * WZ - SS
    S12 = UY + VX
    S23 = VZ + WY
    S31 = WX + UZ
    SS2 = S11 * S11 + S22 * S22 + S33 * S33 &
    &   + 2.0 * (S12 * S12 + S23 * S23 + S31 * S31)
    O12 = UY - VX
    O23 = VZ - WY
    O31 = WX - UZ
    OO2 = 2.0 * (O12 * O12 + O23 * O23 + O31 * O31)
    IF(EPSTM .GT. ZERO) THEN
      AKAE = AKM / EPSTM
    ELSE
      AKAE = 0.0
    ENDIF
    ST = AKAE * SQRT(0.5 * SS2)
    OT = AKAE * SQRT(0.5 * OO2)
    ! 非等方レイノルズ応力テンソル a_ij の係数の計算 -------------------
    nu = AMUM / RHOM
    IF(EPSTM .GT. ZERO) THEN
      Rtt = AKM**2 / (nu * EPSTM)
    ELSE
      Rtt = 0.0
    ENDIF
    eta  = MAX(ST, OT)
    c_mu = 0.3 / (1.0 + 0.35 * eta**1.5) &
    &    * (1.0 - EXP(-0.36 / MAX(ZERO, EXP(-0.75 * eta))))
    f_mu = 1.0 - EXP(-SQRT(Rtt / 90.0) - (Rtt / 400.0)**2)
    c4   =-10.0 * c_mu**2
    c6   =- 5.0 * c_mu**2
    c7   =  5.0 * c_mu**2
    IF(EPSTM .GT. ZERO) THEN
      nu_t = c_mu * f_mu * AKM**2 / EPSTM
    ELSE
      nu_t = 0.0
    ENDIF
    ! 非等方レイノルズ応力テンソル a_ij --------------------------------
    IF(AKM .GT. ZERO .OR. EPSTM**2 .GT. ZERO) THEN
      nu_t1 =-nu_t / AKM
      nu_t2 = nu_t / EPSTM
      nu_t3 = nu_t * AKM / EPSTM**2
    ELSE
      nu_t1 = 0.0
      nu_t2 = 0.0
      nu_t3 = 0.0
    ENDIF
    a11 = nu_t1 * S11 &
    &   + nu_t2 * ( &
    &       c1 * (S11 * S11 + S31 * S31 + S12 * S12 - SS2 / 3.0) &
    &     + c2 * (O12 * S12 - O31 * S31) * 2.0 &
    &     + c3 * (O12 * O12 + O31 * O31 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( (S12 * S23 + (S33 + S11) * S31) * O31 &
    &            - (S23 * S31 + (S11 + S22) * S12) * O12 ) * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S11 &
    &   )
    a22 = nu_t1 * S22 &
    &   + nu_t2 * ( &
    &       c1 * (S22 * S22 + S12 * S12 + S23 * S23 - SS2 / 3.0) &
    &     + c2 * (O23 * S23 - O12 * S12) * 2.0 &
    &     + c3 * (O23 * O23 + O12 * O12 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( (S23 * S31 + (S11 + S22) * S12) * O12 &
    &            - (S31 * S12 + (S22 + S33) * S23) * O23 ) * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S22 &
    &   )
    a33 = nu_t1 * S33 &
    &   + nu_t2 * ( &
    &       c1 * (S33 * S33 + S23 * S23 + S31 * S31 - SS2 / 3.0) &
    &     + c2 * (O31 * S31 - O23 * S23) * 2.0 &
    &     + c3 * (O31 * O31 + O23 * O23 - OO2 / 3.0) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( (S31 * S12 + (S22 + S33) * S23) * O23 &
    &            - (S12 * S23 + (S33 + S11) * S31) * O31) * 2.0 &
    &     + (c6 * SS2 + c7 * OO2) * S33 &
    &   )
    a12 = nu_t1 * S12 &
    &   + nu_t2 * ( &
    &       c1 * ((S11 + S22) * S12 + S23 * S31) &
    &     + c2 * (O23 * S31 - O31 * S23 - O12 * (S11 - S22)) &
    &     - c3 * (O23 * O31) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( &
    &         (S11 * S11 - S22 * S22 - S23 * S23 + S31 * S31) * O12 &
    &       + ((S22 + S33) * S23 + S31 * S12) * O31 &
    &       - ((S33 + S11) * S31 + S12 * S23) * O23 &
    &     ) &
    &     + (c6 * SS2 + c7 * OO2) * S12 &
    &   )
    a23 = nu_t1 * S23 &
    &   + nu_t2 * ( &
    &       c1 * ((S22 + S33) * S23 + S31 * S12) &
    &     + c2 * (O31 * S12 - O12 * S31 - O23 * (S22 - S33)) &
    &     - c3 * (O31 * O12) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( &
    &         (S22 * S22 - S33 * S33 - S31 * S31 + S12 * S12) * O23 &
    &       + ((S11 + S33) * S31 + S12 * S23) * O12 &
    &       - ((S11 + S22) * S12 + S23 * S31) * O31 &
    &     ) &
    &     + (c6 * SS2 + c7 * OO2) * S23 &
    &   )
    a31 = nu_t1 * S31 &
    &   + nu_t2 * ( &
    &       c1 * ((S11 + S33) * S31 + S12 * S23) &
    &     + c2 * (O12 * S23 - O23 * S12 - O31 * (S33 - S11)) &
    &     - c3 * (O12 * O23) &
    &   ) &
    &   + nu_t3 * ( &
    &       c4 * ( &
    &         (S33 * S33 - S11 * S11 - S12 * S12 + S23 * S23) * O31 &
    &       + ((S11 + S22) * S12 + S31 * S23) * O23 &
    &       - ((S22 + S33) * S23 + S12 * S31) * O12 &
    &     ) &
    &     + (c6 * SS2 + c7 * OO2) * S31 &
    &   )
    a11 = MIN(aiimax, MAX(aiimin, a11))
    a22 = MIN(aiimax, MAX(aiimin, a22))
    a33 = MIN(aiimax, MAX(aiimin, a33))
    a12 = MIN(1.0, MAX(-1.0, a12))
    a23 = MIN(1.0, MAX(-1.0, a23))
    a31 = MIN(1.0, MAX(-1.0, a31))
    ! レイノルズ応力 ---------------------------------------------------
    uu = AKM * (a11 + TwoThird)
    vv = AKM * (a22 + TwoThird)
    ww = AKM * (a33 + TwoThird)
    uv = AKM * a12
    vw = AKM * a23
    wu = AKM * a31
    ! 拡散係数 ---------------------------------------------------------
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    c_k  = nu + nu_t / SIGK
    c_e  = nu + nu_t / SIGE
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    TAUXX = nu * S11 - uu
    TAUYY = nu * S22 - vv
    TAUZZ = nu * S33 - ww
    TAUXY = nu * S12 - uv
    TAUYZ = nu * S23 - vw
    TAUZX = nu * S31 - wu
    ! エネルギーの拡散 -------------------------------------------------
    R5 = TAUXX * UM + TAUXY * VM + TAUZX * WM + c_t * TX
    S5 = TAUXY * UM + TAUYY * VM + TAUYZ * WM + c_t * TY
    T5 = TAUZX * UM + TAUYZ * VM + TAUZZ * WM + c_t * TZ
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    TH(I,J,K,1) = 0.0
    TH(I,J,K,2) = QH1 * (TAUXX * ZEXM + TAUXY * ZEYM + TAUZX * ZEZM)
    TH(I,J,K,3) = QH1 * (TAUXY * ZEXM + TAUYY * ZEYM + TAUYZ * ZEZM)
    TH(I,J,K,4) = QH1 * (TAUZX * ZEXM + TAUYZ * ZEYM + TAUZZ * ZEZM)
    TH(I,J,K,5) = QH1 * (R5 * ZEXM + S5 * ZEYM + T5 * ZEZM)
    TH(I,J,K,6) = QH1 * c_k * (AKX * ZEXM + AKY * ZEYM + AKZ * ZEZM)
    TH(I,J,K,7) = QH1 * c_e * (AEX * ZEXM + AEY * ZEYM + AEZ * ZEZM)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFZE
! 定義終了 *************************************************************
END SUBROUTINE Turbulence3DNlEvmkeCLS
!***********************************************************************
!**** 乱流モデル : RANS, RSM Gibson-Launder Model (1978)            ****
!****              乱流熱伝達にはVandrommeのモデルを使用            ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence2DRSMGL( &
&            PELIM, RG, GAMMA, PR, &
&            IS, IE, JS, JE, LS, LE, &
&            XIX, XIY, ETX, ETY, AJA, &
&            RHO, U, V, T, uu, vv, ww, uv, EPS, AMU, &
&            n1, n2, d, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: C_se = 0.313
  REAL, PARAMETER :: C_s = 0.22, C_e = 0.18
  REAL, PARAMETER :: C_1 = 1.8, C_2 = 0.6
  REAL, PARAMETER :: C_l = 2.5
  REAL, PARAMETER :: C_w1 = 0.5, C_w2 = 0.3
  REAL, PARAMETER :: C_e1 = 1.44, C_e2 = 1.92
  REAL, PARAMETER :: TwoThird = 0.66666667
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE), XIY(IS:IE, JS:JE), &
  &                       ETX(IS:IE, JS:JE), ETY(IS:IE, JS:JE), &
  &                       AJA(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE), &
  &                       U(IS:IE, JS:JE), &
  &                       V(IS:IE, JS:JE), &
  &                       T(IS:IE, JS:JE), &
  &                       uu(IS:IE, JS:JE), &
  &                       vv(IS:IE, JS:JE), &
  &                       ww(IS:IE, JS:JE), &
  &                       uv(IS:IE, JS:JE), &
  &                       EPS(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: n1(IS:IE, JS:JE), n2(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: d(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: AMUT(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(OUT) :: DQP(IS:IE, JS:JE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :), SH(:, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 9) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, LS: LE), SH(IS: IE, JS: JE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: TMP, TMP1, TMPw, TMPw1, TMPw2, TMPw1km, TMPw2km
  REAL    :: UXI, VXI, UET, VET, UX, VX, UY, VY
  REAL    :: P11, P22, P33, P12, P_k
  REAL    :: eps11, eps22, eps33, eps12, e_k
  REAL    :: phi1_11, phi1_22, phi1_33, phi1_12
  REAL    :: phi2_11, phi2_22, phi2_33, phi2_12
  REAL    :: phiw1_11, phiw1_22, phiw1_33, phiw1_12
  REAL    :: phiw2_11, phiw2_22, phiw2_33, phiw2_12
  REAL    :: phi11, phi22, phi33, phi12
  REAL    :: AK, SS2, nu, AEAK
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  TMP, TMP1, TMPw, TMPw1, TMPw2, TMPw1km, TMPw2km, &
  !$  UXI, VXI, UET, VET, UX, VX, UY, VY, &
  !$  P11, P22, P33, P12, P_k, &
  !$  eps11, eps22, eps33, eps12, e_k, &
  !$  phi1_11, phi1_22, phi1_33, phi1_12, &
  !$  phi2_11, phi2_22, phi2_33, phi2_12, &
  !$  phiw1_11, phiw1_22, phiw1_33, phiw1_12, &
  !$  phiw2_11, phiw2_22, phiw2_33, phiw2_12, &
  !$  phi11, phi22, phi33, phi12, &
  !$  AK, SS2, nu, AEAK, QH1 &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF(RHO(I,J) .GT. 0.0) THEN
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.5 * (U(I+1,J) - U(I-1,J))
    VXI = 0.5 * (V(I+1,J) - V(I-1,J))
    UET = 0.5 * (U(I,J+1) - U(I,J-1))
    VET = 0.5 * (V(I,J+1) - V(I,J-1))
    ! 物理空間方向一階微分
    UX = UXI * XIX(I,J) + UET * ETX(I,J)
    UY = UXI * XIY(I,J) + UET * ETY(I,J)
    VX = VXI * XIX(I,J) + VET * ETX(I,J)
    VY = VXI * XIY(I,J) + VET * ETY(I,J)
    ! 生成項と散逸項 ---------------------------------------------------
    P11 = - 2.0 * (uu(I,J) * UX + uv(I,J) * UY)
    P22 = - 2.0 * (uv(I,J) * VX + vv(I,J) * VY)
    P33 = 0.0
    P12 = - (uu(I,J) * VX + uv(I,J) * (UX + VY) + vv(I,J) * UY)
    P_k = 0.5 * (P11 + P22 + P33)
    e_k = EPS(I,J)
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      TMP = TwoThird * e_k / PELIM
      P11 = MIN(TMP, MAX(-TMP, P11))
      P22 = MIN(TMP, MAX(-TMP, P22))
      P33 = MIN(TMP, MAX(-TMP, P33))
      P_k = 0.5 * (P11 + P22 + P33)
    ENDIF
    TMP   = TwoThird * e_k
    eps11 = TMP
    eps22 = TMP
    eps33 = TMP
    eps12 = 0.0
    ! e / k (時間スケールリミッター) -----------------------------------
    AK   = 0.5 * (uu(I,J) + vv(I,J) + ww(I,J))
    SS2  = 2.0 * (UX**2 + VY**2) + (UY + VX)**2 &
    &    - TwoThird * (UX + VY)**2
    nu   = AMU(I,J) / RHO(I,J)
    AEAK = MIN(MAX(SQRT(SS2), SQRT(e_k / nu)), e_k / MAX(ZERO, AK))
    ! 再分配項 ---------------------------------------------------------
    TMP1    = - C_1 * AEAK
    phi1_11 = TMP1 * (uu(I,J) - TwoThird * AK)
    phi1_22 = TMP1 * (vv(I,J) - TwoThird * AK)
    phi1_33 = TMP1 * (ww(I,J) - TwoThird * AK)
    phi1_12 = TMP1 * uv(I,J)
    phi2_11 = - C_2 * (P11 - TwoThird * P_k)
    phi2_22 = - C_2 * (P22 - TwoThird * P_k)
    phi2_33 = - C_2 * (P33 - TwoThird * P_k)
    phi2_12 = - C_2 * P12
    IF(C_l * e_k * d(I,J) .GT. ZERO) THEN
      TMPw = AK**1.5 / MAX(ZERO, (C_l * e_k * d(I,J)))
      ! 壁面反射項のダンピング関数の制限
      ! TMPw = MIN(1.0, MAX(0.0, TMPw))
    ELSE
      TMPw = 0.0
    ENDIF
    TMPw1   = C_w1 * TMPw * AEAK
    TMPw2   = C_w2 * TMPw
    TMPw1km = uu(I,J) * n1(I,J) * n1(I,J) &
    &       + vv(I,J) * n2(I,J) * n2(I,J) &
    &       + 2.0 * uv(I,J) * n1(I,J) * n2(I,J)
    TMPw2km = phi2_11 * n1(I,J) * n1(I,J) &
    &       + phi2_22 * n2(I,J) * n2(I,J) &
    &       + 2.0 * phi2_12 * n1(I,J) * n2(I,J)
    phiw1_11 = TMPw1 * (TMPw1km - 3.0 * n1(I,J) &
    &          * (uu(I,J) * n1(I,J) + uv(I,J) * n2(I,J)) &
    &        )
    phiw1_22 = TMPw1 * (TMPw1km - 3.0 * n2(I,J) &
    &          * (uv(I,J) * n1(I,J) + vv(I,J) * n2(I,J)) &
    &        )
    phiw1_33 = TMPw1 * TMPw1km
    phiw1_12 = - 1.5 * TMPw1 * ( &
    &          n1(I,J) * (uv(I,J) * n1(I,J) + vv(I,J) * n2(I,J)) &
    &        + n2(I,J) * (uu(I,J) * n1(I,J) + uv(I,J) * n2(I,J)) &
    &        )
    phiw2_11 = TMPw2 * (TMPw2km - 3.0 * n1(I,J) &
    &          * (phi2_11 * n1(I,J) + phi2_12 * n2(I,J)) &
    &        )
    phiw2_22 = TMPw2 * (TMPw2km - 3.0 * n2(I,J) &
    &          * (phi2_12 * n1(I,J) + phi2_22 * n2(I,J)) &
    &        )
    phiw2_33 = TMPw2 * TMPw2km
    phiw2_12 = - 1.5 * TMPw2 * ( &
    &          n1(I,J) * (phi2_12 * n1(I,J) + phi2_22 * n2(I,J)) &
    &        + n2(I,J) * (phi2_11 * n1(I,J) + phi2_12 * n2(I,J)) &
    &        )
    phi11 = phi1_11 + phi2_11 + phiw1_11 + phiw2_11
    phi22 = phi1_22 + phi2_22 + phiw1_22 + phiw2_22
    phi33 = phi1_33 + phi2_33 + phiw1_33 + phiw2_33
    phi12 = phi1_12 + phi2_12 + phiw1_12 + phiw2_12
    ! 渦粘性係数の計算 -------------------------------------------------
    IF(e_k .GT. ZERO) THEN
      AMUT(I,J) = 0.09 * AK**2 / e_k * RHO(I,J)
    ENDIF
    ! 生成項、圧力再分配項、散逸項の和 ---------------------------------
    QH1 = RHO(I,J) / AJA(I,J)
    DQP(I,J,5) = QH1 * (P11 - eps11 + phi11)
    DQP(I,J,6) = QH1 * (P22 - eps22 + phi22)
    DQP(I,J,7) = QH1 * (P33 - eps33 + phi33)
    DQP(I,J,8) = QH1 * (P12 - eps12 + phi12)
    DQP(I,J,9) = QH1 * AEAK * (C_e1 * P_k - C_e2 * e_k)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH = 0.0
  SH = 0.0
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL DIFFXI
  CALL DIFFET
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, L)
  DO L = LS, LE
  !$OMP DO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J-1) .GT. 0.0) .AND. (RHO(I-1,J) .GT. 0.0) .AND. &
  &   (RHO(I,J)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,L) = (-RH(I-1,J  ,L) + RH(I,J,L)) &
    &          + (-SH(I  ,J-1,L) + SH(I,J,L))
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項の計算(xi方向差分用)                                    ****
!***********************************************************************
SUBROUTINE DIFFXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: RHOM, UM, VM, uuM, vvM, wwM, uvM, EPSM, AMUM
  REAL    :: UXI, UET, UX, UY
  REAL    :: VXI, VET, VX, VY
  REAL    :: TXI, TET, TX, TY
  REAL    :: uuXI, uuET, uuX, uuY
  REAL    :: vvXI, vvET, vvX, vvY
  REAL    :: wwXI, wwET, wwX, wwY
  REAL    :: uvXI, uvET, uvX, uvY
  REAL    :: AEXI, AEET, AEX, AEY
  REAL    :: nu
  REAL    :: S11, S22, S12, SS
  REAL    :: C_p, AKAE
  REAL    :: R2M, S2M, R2T, S2T
  REAL    :: R3M, S3M, R3T, S3T
  REAL    :: R4M, S4M, R4T, S4T
  REAL    :: R5M, S5M, R5T, S5T
  REAL    :: R6M, S6M, R6T, S6T
  REAL    :: R7M, S7M, R7T, S7T
  REAL    :: R8M, S8M, R8T, S8T
  REAL    :: R9M, S9M, R9T, S9T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, &
  !$  RHOM, UM, VM, uuM, vvM, wwM, uvM, EPSM, AMUM, &
  !$  UXI, UET, UX, UY, VXI, VET, VX, VY, TXI, TET, TX, TY, &
  !$  uuXI, uuET, uuX, uuY, vvXI, vvET, vvX, vvY, &
  !$  wwXI, wwET, wwX, wwY, uvXI, uvET, uvX, uvY, &
  !$  AEXI, AEET, AEX, AEY, &
  !$  nu, S11, S22, S12, SS, C_p, AKAE, &
  !$  R2M, S2M, R2T, S2T, R3M, S3M, R3T, S3T, R4M, S4M, R4T, S4T, &
  !$  R5M, S5M, R5T, S5T, R6M, S6M, R6T, S6T, R7M, S7M, R7T, S7T, &
  !$  R8M, S8M, R8T, S8T, R9M, S9M, R9T, S9T, QH1 &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF((RHO(I,J) .GT. 0.0) .AND. (RHO(I+1,J) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J) + XIX(I+1,J))
    XIYM = 0.5 * (XIY(I,J) + XIY(I+1,J))
    ETXM = 0.5 * (ETX(I,J) + ETX(I+1,J))
    ETYM = 0.5 * (ETY(I,J) + ETY(I+1,J))
    AJAM = 0.5 * (AJA(I,J) + AJA(I+1,J))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J) + RHO(I+1,J))
    UM   = 0.5 * (  U(I,J) +   U(I+1,J))
    VM   = 0.5 * (  V(I,J) +   V(I+1,J))
    uuM  = 0.5 * ( uu(I,J) +  uu(I+1,J))
    vvM  = 0.5 * ( vv(I,J) +  vv(I+1,J))
    wwM  = 0.5 * ( ww(I,J) +  ww(I+1,J))
    uvM  = 0.5 * ( uv(I,J) +  uv(I+1,J))
    EPSM = 0.5 * (EPS(I,J) + EPS(I+1,J))
    AMUM = 0.5 * (AMU(I,J) + AMU(I+1,J))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = - U(I,J) + U(I+1,J)
    UET = 0.25 * (- U(I,J-1) - U(I+1,J-1) + U(I,J+1) + U(I+1,J+1))
    VXI = - V(I,J) + V(I+1,J)
    VET = 0.25 * (- V(I,J-1) - V(I+1,J-1) + V(I,J+1) + V(I+1,J+1))
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM
    UY = UXI * XIYM + UET * ETYM
    VX = VXI * XIXM + VET * ETXM
    VY = VXI * XIYM + VET * ETYM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = - T(I,J) + T(I+1,J)
    TET = 0.25 * (- T(I,J-1) - T(I+1,J-1) + T(I,J+1) + T(I+1,J+1))
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM
    TY = TXI * XIYM + TET * ETYM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = - uu(I,J) + uu(I+1,J)
    uuET = 0.25 * ( - uu(I,J-1) - uu(I+1,J-1) &
    &               + uu(I,J+1) + uu(I+1,J+1) )
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM
    uuY = uuXI * XIYM + uuET * ETYM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = - vv(I,J) + vv(I+1,J)
    vvET = 0.25 * ( - vv(I,J-1) - vv(I+1,J-1) &
    &               + vv(I,J+1) + vv(I+1,J+1) )
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM
    vvY = vvXI * XIYM + vvET * ETYM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = - ww(I,J) + ww(I+1,J)
    wwET = 0.25 * ( - ww(I,J-1) - ww(I+1,J-1) &
    &               + ww(I,J+1) + ww(I+1,J+1) )
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM
    wwY = wwXI * XIYM + wwET * ETYM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = - uv(I,J) + uv(I+1,J)
    uvET = 0.25 * ( - uv(I,J-1) - uv(I+1,J-1) &
    &               + uv(I,J+1) + uv(I+1,J+1) )
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM
    uvY = uvXI * XIYM + uvET * ETYM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = - EPS(I,J) + EPS(I+1,J)
    AEET = 0.25 * ( - EPS(I,J-1) - EPS(I+1,J-1) &
    &               + EPS(I,J+1) + EPS(I+1,J+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM
    AEY = AEXI * XIYM + AEET * ETYM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S12 = UY + VX
    ! U の拡散 ---------------------------------------------------------
    R2M = nu * S11
    S2M = nu * S12
    R2T =-uuM
    S2T =-uvM
    ! V の拡散 ---------------------------------------------------------
    R3M = nu * S12
    S3M = nu * S22
    R3T =-uvM
    S3T =-vvM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R4M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM)
    S4M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM)
    R4T = AKAE * C_se * C_p * (uuM * TX + uvM * TY) &
    &   + AKAE * C_s  * ( &
    &     uuM * (uuX + vvX + wwX) + uvM * (uuY + vvY + wwY) &
    &   ) &
    &   - (uuM * UM + uvM * VM)
    S4T = AKAE * C_se * C_p * (uvM * TX + vvM * TY) &
    &   + AKAE * C_s  * ( &
    &     uvM * (uuX + vvX + wwX) + vvM * (uuY + vvY + wwY) &
    &   ) &
    &   - (uvM * UM + vvM * VM)
    ! uu の拡散 --------------------------------------------------------
    R5M = nu * uuX
    S5M = nu * uuY
    R5T = C_s * AKAE * (uuM * uuX + uvM * uuY)
    S5T = C_s * AKAE * (uvM * uuX + vvM * uuY)
    ! vv の拡散 --------------------------------------------------------
    R6M = nu * vvX
    S6M = nu * vvY
    R6T = C_s * AKAE * (uuM * vvX + uvM * vvY)
    S6T = C_s * AKAE * (uvM * vvX + vvM * vvY)
    ! ww の拡散 --------------------------------------------------------
    R7M = nu * wwX
    S7M = nu * wwY
    R7T = C_s * AKAE * (uuM * wwX + uvM * wwY)
    S7T = C_s * AKAE * (uvM * wwX + vvM * wwY)
    ! uv の拡散 --------------------------------------------------------
    R8M = nu * uvX
    S8M = nu * uvY
    R8T = C_s * AKAE * (uuM * uvX + uvM * uvY)
    S8T = C_s * AKAE * (uvM * uvX + vvM * uvY)
    ! epsilon の拡散 ---------------------------------------------------
    R9M = nu * AEX
    S9M = nu * AEY
    R9T = C_e * AKAE * (uuM * AEX + uvM * AEY)
    S9T = C_e * AKAE * (uvM * AEX + vvM * AEY)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    RH(I,J,1) = 0.0
    RH(I,J,2) = QH1 * ((R2M + R2T) * XIXM + (S2M + S2T) * XIYM)
    RH(I,J,3) = QH1 * ((R3M + R3T) * XIXM + (S3M + S3T) * XIYM)
    RH(I,J,4) = QH1 * ((R4M + R4T) * XIXM + (S4M + S4T) * XIYM)
    RH(I,J,5) = QH1 * ((R5M + R5T) * XIXM + (S5M + S5T) * XIYM)
    RH(I,J,6) = QH1 * ((R6M + R6T) * XIXM + (S6M + S6T) * XIYM)
    RH(I,J,7) = QH1 * ((R7M + R7T) * XIXM + (S7M + S7T) * XIYM)
    RH(I,J,8) = QH1 * ((R8M + R8T) * XIXM + (S8M + S8T) * XIYM)
    RH(I,J,9) = QH1 * ((R9M + R9T) * XIXM + (S9M + S9T) * XIYM)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFXI
!***********************************************************************
!**** 拡散項の計算(eta方向差分用)                                   ****
!***********************************************************************
SUBROUTINE DIFFET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: RHOM, UM, VM, uuM, vvM, wwM, uvM, EPSM, AMUM
  REAL    :: UXI, UET, UX, UY
  REAL    :: VXI, VET, VX, VY
  REAL    :: TXI, TET, TX, TY
  REAL    :: uuXI, uuET, uuX, uuY
  REAL    :: vvXI, vvET, vvX, vvY
  REAL    :: wwXI, wwET, wwX, wwY
  REAL    :: uvXI, uvET, uvX, uvY
  REAL    :: AEXI, AEET, AEX, AEY
  REAL    :: nu
  REAL    :: S11, S22, S12, SS
  REAL    :: C_p, AKAE
  REAL    :: R2M, S2M, R2T, S2T
  REAL    :: R3M, S3M, R3T, S3T
  REAL    :: R4M, S4M, R4T, S4T
  REAL    :: R5M, S5M, R5T, S5T
  REAL    :: R6M, S6M, R6T, S6T
  REAL    :: R7M, S7M, R7T, S7T
  REAL    :: R8M, S8M, R8T, S8T
  REAL    :: R9M, S9M, R9T, S9T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, &
  !$  RHOM, UM, VM, uuM, vvM, wwM, uvM, EPSM, AMUM, &
  !$  UXI, UET, UX, UY, VXI, VET, VX, VY, TXI, TET, TX, TY, &
  !$  uuXI, uuET, uuX, uuY, vvXI, vvET, vvX, vvY, &
  !$  wwXI, wwET, wwX, wwY, uvXI, uvET, uvX, uvY, &
  !$  AEXI, AEET, AEX, AEY, &
  !$  nu, S11, S22, S12, SS, C_p, AKAE, &
  !$  R2M, S2M, R2T, S2T, R3M, S3M, R3T, S3T, R4M, S4M, R4T, S4T, &
  !$  R5M, S5M, R5T, S5T, R6M, S6M, R6T, S6T, R7M, S7M, R7T, S7T, &
  !$  R8M, S8M, R8T, S8T, R9M, S9M, R9T, S9T, QH1 &
  !$)
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF((RHO(I,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0)) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J) + XIX(I,J+1))
    XIYM = 0.5 * (XIY(I,J) + XIY(I,J+1))
    ETXM = 0.5 * (ETX(I,J) + ETX(I,J+1))
    ETYM = 0.5 * (ETY(I,J) + ETY(I,J+1))
    AJAM = 0.5 * (AJA(I,J) + AJA(I,J+1))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J) + RHO(I,J+1))
    UM   = 0.5 * (  U(I,J) +   U(I,J+1))
    VM   = 0.5 * (  V(I,J) +   V(I,J+1))
    uuM  = 0.5 * ( uu(I,J) +  uu(I,J+1))
    vvM  = 0.5 * ( vv(I,J) +  vv(I,J+1))
    wwM  = 0.5 * ( ww(I,J) +  ww(I,J+1))
    uvM  = 0.5 * ( uv(I,J) +  uv(I,J+1))
    EPSM = 0.5 * (EPS(I,J) + EPS(I,J+1))
    AMUM = 0.5 * (AMU(I,J) + AMU(I,J+1))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * (- U(I-1,J) - U(I-1,J+1) + U(I+1,J) + U(I+1,J+1))
    UET = - U(I,J) + U(I,J+1)
    VXI = 0.25 * (- V(I-1,J) - V(I-1,J+1) + V(I+1,J) + V(I+1,J+1))
    VET = - V(I,J) + V(I,J+1)
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM
    UY = UXI * XIYM + UET * ETYM
    VX = VXI * XIXM + VET * ETXM
    VY = VXI * XIYM + VET * ETYM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * (- T(I-1,J) - T(I-1,J+1) + T(I+1,J) + T(I+1,J+1))
    TET = - T(I,J) + T(I,J+1)
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM
    TY = TXI * XIYM + TET * ETYM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = 0.25 * (- uu(I-1,J) - uu(I-1,J+1) + uu(I+1,J) + uu(I+1,J+1))
    uuET = - uu(I,J) + uu(I,J+1)
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM
    uuY = uuXI * XIYM + uuET * ETYM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = 0.25 * (- vv(I-1,J) - vv(I-1,J+1) + vv(I+1,J) + vv(I+1,J+1))
    vvET = - vv(I,J) + vv(I,J+1)
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM
    vvY = vvXI * XIYM + vvET * ETYM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = 0.25 * (- ww(I-1,J) - ww(I-1,J+1) + ww(I+1,J) + ww(I+1,J+1))
    wwET = - ww(I,J) + ww(I,J+1)
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM
    wwY = wwXI * XIYM + wwET * ETYM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = 0.25 * (- uv(I-1,J) - uv(I-1,J+1) + uv(I+1,J) + uv(I+1,J+1))
    uvET = - uv(I,J) + uv(I,J+1)
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM
    uvY = uvXI * XIYM + uvET * ETYM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPS(I-1,J) - EPS(I-1,J+1) &
    &               + EPS(I+1,J) + EPS(I+1,J+1) )
    AEET = - EPS(I,J) + EPS(I,J+1)
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM
    AEY = AEXI * XIYM + AEET * ETYM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S12 = UY + VX
    ! U の拡散 ---------------------------------------------------------
    R2M = nu * S11
    S2M = nu * S12
    R2T =-uuM
    S2T =-uvM
    ! V の拡散 ---------------------------------------------------------
    R3M = nu * S12
    S3M = nu * S22
    R3T =-uvM
    S3T =-vvM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R4M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM)
    S4M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM)
    R4T = AKAE * C_se * C_p * (uuM * TX + uvM * TY) &
    &   + AKAE * C_s  * ( &
    &     uuM * (uuX + vvX + wwX) + uvM * (uuY + vvY + wwY) &
    &   ) &
    &   - (uuM * UM + uvM * VM)
    S4T = AKAE * C_se * C_p * (uvM * TX + vvM * TY) &
    &   + AKAE * C_s  * ( &
    &     uvM * (uuX + vvX + wwX) + vvM * (uuY + vvY + wwY) &
    &   ) &
    &   - (uvM * UM + vvM * VM)
    ! uu の拡散 --------------------------------------------------------
    R5M = nu * uuX
    S5M = nu * uuY
    R5T = C_s * AKAE * (uuM * uuX + uvM * uuY)
    S5T = C_s * AKAE * (uvM * uuX + vvM * uuY)
    ! vv の拡散 --------------------------------------------------------
    R6M = nu * vvX
    S6M = nu * vvY
    R6T = C_s * AKAE * (uuM * vvX + uvM * vvY)
    S6T = C_s * AKAE * (uvM * vvX + vvM * vvY)
    ! ww の拡散 --------------------------------------------------------
    R7M = nu * wwX
    S7M = nu * wwY
    R7T = C_s * AKAE * (uuM * wwX + uvM * wwY)
    S7T = C_s * AKAE * (uvM * wwX + vvM * wwY)
    ! uv の拡散 --------------------------------------------------------
    R8M = nu * uvX
    S8M = nu * uvY
    R8T = C_s * AKAE * (uuM * uvX + uvM * uvY)
    S8T = C_s * AKAE * (uvM * uvX + vvM * uvY)
    ! epsilon の拡散 ---------------------------------------------------
    R9M = nu * AEX
    S9M = nu * AEY
    R9T = C_e * AKAE * (uuM * AEX + uvM * AEY)
    S9T = C_e * AKAE * (uvM * AEX + vvM * AEY)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    SH(I,J,1) = 0.0
    SH(I,J,2) = QH1 * ((R2M + R2T) * ETXM + (S2M + S2T) * ETYM)
    SH(I,J,3) = QH1 * ((R3M + R3T) * ETXM + (S3M + S3T) * ETYM)
    SH(I,J,4) = QH1 * ((R4M + R4T) * ETXM + (S4M + S4T) * ETYM)
    SH(I,J,5) = QH1 * ((R5M + R5T) * ETXM + (S5M + S5T) * ETYM)
    SH(I,J,6) = QH1 * ((R6M + R6T) * ETXM + (S6M + S6T) * ETYM)
    SH(I,J,7) = QH1 * ((R7M + R7T) * ETXM + (S7M + S7T) * ETYM)
    SH(I,J,8) = QH1 * ((R8M + R8T) * ETXM + (S8M + S8T) * ETYM)
    SH(I,J,9) = QH1 * ((R9M + R9T) * ETXM + (S9M + S9T) * ETYM)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFET
! 定義終了 *************************************************************
END SUBROUTINE Turbulence2DRSMGL
!***********************************************************************
!**** 乱流モデル : RANS, RSM Gibson-Launder Model (1978)            ****
!****              乱流熱伝達にはVandrommeのモデルを使用            ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence3DRSMGL( &
&            PELIM, RG, GAMMA, PR, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, AJA, &
&            RHO, U, V, W, T, uu, vv, ww, uv, vw, wu, EPS, AMU, &
&            n1, n2, n3, d, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: C_se = 0.313
  REAL, PARAMETER :: C_s = 0.22, C_e = 0.18
  REAL, PARAMETER :: C_1 = 1.8, C_2 = 0.6
  REAL, PARAMETER :: C_l = 2.5
  REAL, PARAMETER :: C_w1 = 0.5, C_w2 = 0.3
  REAL, PARAMETER :: C_e1 = 1.44, C_e2 = 1.92
  REAL, PARAMETER :: TwoThird = 0.66666667
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE, KS:KE), &
  &                       XIY(IS:IE, JS:JE, KS:KE), &
  &                       XIZ(IS:IE, JS:JE, KS:KE), &
  &                       ETX(IS:IE, JS:JE, KS:KE), &
  &                       ETY(IS:IE, JS:JE, KS:KE), &
  &                       ETZ(IS:IE, JS:JE, KS:KE), &
  &                       ZEX(IS:IE, JS:JE, KS:KE), &
  &                       ZEY(IS:IE, JS:JE, KS:KE), &
  &                       ZEZ(IS:IE, JS:JE, KS:KE), &
  &                       AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE, KS:KE), &
  &                       U(IS:IE, JS:JE, KS:KE), &
  &                       V(IS:IE, JS:JE, KS:KE), &
  &                       W(IS:IE, JS:JE, KS:KE), &
  &                       T(IS:IE, JS:JE, KS:KE), &
  &                       uu(IS:IE, JS:JE, KS:KE), &
  &                       vv(IS:IE, JS:JE, KS:KE), &
  &                       ww(IS:IE, JS:JE, KS:KE), &
  &                       uv(IS:IE, JS:JE, KS:KE), &
  &                       vw(IS:IE, JS:JE, KS:KE), &
  &                       wu(IS:IE, JS:JE, KS:KE), &
  &                       EPS(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: n1(IS:IE, JS:JE, KS:KE), &
  &                       n2(IS:IE, JS:JE, KS:KE), &
  &                       n3(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: d(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: AMUT(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(OUT) :: DQP(IS:IE, JS:JE, KS:KE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :, :), SH(:, :, :, :), TH(:, :, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 12) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(SH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(TH(IS: IE, JS: JE, KS: KE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: TMP, TMP1, TMPw, TMPw1, TMPw2, TMPw1km, TMPw2km
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  REAL    :: UX, VX, WX, UY, VY, WY, UZ, VZ, WZ
  REAL    :: P11, P22, P33, P12, P23, P31, P_k
  REAL    :: eps11, eps22, eps33, eps12, eps23, eps31, e_k
  REAL    :: phi1_11, phi1_22, phi1_33, phi1_12, phi1_23, phi1_31
  REAL    :: phi2_11, phi2_22, phi2_33, phi2_12, phi2_23, phi2_31
  REAL    :: phiw1_11, phiw1_22, phiw1_33, phiw1_12, phiw1_23, phiw1_31
  REAL    :: phiw2_11, phiw2_22, phiw2_33, phiw2_12, phiw2_23, phiw2_31
  REAL    :: phi11, phi22, phi33, phi12, phi23, phi31
  REAL    :: AK, SS2, nu, AEAK
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  TMP, TMP1, TMPw, TMPw1, TMPw2, TMPw1km, TMPw2km, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE, &
  !$  UX, VX, WX, UY, VY, WY, UZ, VZ, WZ, &
  !$  P11, P22, P33, P12, P23, P31, P_k, &
  !$  eps11, eps22, eps33, eps12, eps23, eps31, e_k, &
  !$  phi1_11, phi1_22, phi1_33, phi1_12, phi1_23, phi1_31, &
  !$  phi2_11, phi2_22, phi2_33, phi2_12, phi2_23, phi2_31, &
  !$  phiw1_11, phiw1_22, phiw1_33, phiw1_12, phiw1_23, phiw1_31, &
  !$  phiw2_11, phiw2_22, phiw2_33, phiw2_12, phiw2_23, phiw2_31, &
  !$  phi11, phi22, phi33, phi12, phi23, phi31, AK, SS2, nu, AEAK, QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF(RHO(I,J,K) .GT. 0.0) THEN
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.5 * (U(I+1,J,K) - U(I-1,J,K))
    UET = 0.5 * (U(I,J+1,K) - U(I,J-1,K))
    UZE = 0.5 * (U(I,J,K+1) - U(I,J,K-1))
    VXI = 0.5 * (V(I+1,J,K) - V(I-1,J,K))
    VET = 0.5 * (V(I,J+1,K) - V(I,J-1,K))
    VZE = 0.5 * (V(I,J,K+1) - V(I,J,K-1))
    WXI = 0.5 * (W(I+1,J,K) - W(I-1,J,K))
    WET = 0.5 * (W(I,J+1,K) - W(I,J-1,K))
    WZE = 0.5 * (W(I,J,K+1) - W(I,J,K-1))
    ! 物理空間方向一階微分
    UX = UXI * XIX(I,J,K) + UET * ETX(I,J,K) + UZE * ZEX(I,J,K)
    UY = UXI * XIY(I,J,K) + UET * ETY(I,J,K) + UZE * ZEY(I,J,K)
    UZ = UXI * XIZ(I,J,K) + UET * ETZ(I,J,K) + UZE * ZEZ(I,J,K)
    VX = VXI * XIX(I,J,K) + VET * ETX(I,J,K) + VZE * ZEX(I,J,K)
    VY = VXI * XIY(I,J,K) + VET * ETY(I,J,K) + VZE * ZEY(I,J,K)
    VZ = VXI * XIZ(I,J,K) + VET * ETZ(I,J,K) + VZE * ZEZ(I,J,K)
    WX = WXI * XIX(I,J,K) + WET * ETX(I,J,K) + WZE * ZEX(I,J,K)
    WY = WXI * XIY(I,J,K) + WET * ETY(I,J,K) + WZE * ZEY(I,J,K)
    WZ = WXI * XIZ(I,J,K) + WET * ETZ(I,J,K) + WZE * ZEZ(I,J,K)
    ! 生成項と散逸項 ---------------------------------------------------
    P11 = - 2.0 * (uu(I,J,K) * UX + uv(I,J,K) * UY + wu(I,J,K) * UZ)
    P22 = - 2.0 * (uv(I,J,K) * VX + vv(I,J,K) * VY + vw(I,J,K) * VZ)
    P33 = - 2.0 * (wu(I,J,K) * WX + vw(I,J,K) * WY + ww(I,J,K) * WZ)
    P12 = - (uv(I,J,K) * UX + uu(I,J,K) * VX) &
    &     - (vv(I,J,K) * UY + uv(I,J,K) * VY) &
    &     - (vw(I,J,K) * UZ + wu(I,J,K) * VZ)
    P23 = - (wu(I,J,K) * VX + uv(I,J,K) * WX) &
    &     - (vw(I,J,K) * VY + vv(I,J,K) * WY) &
    &     - (ww(I,J,K) * VZ + vw(I,J,K) * WZ)
    P31 = - (uu(I,J,K) * WX + wu(I,J,K) * UX) &
    &     - (uv(I,J,K) * WY + vw(I,J,K) * UY) &
    &     - (wu(I,J,K) * WZ + ww(I,J,K) * UZ)
    P_k = 0.5 * (P11 + P22 + P33)
    e_k = EPS(I,J,K)
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      TMP = TwoThird * e_k / PELIM
      P11 = MIN(TMP, MAX(-TMP, P11))
      P22 = MIN(TMP, MAX(-TMP, P22))
      P33 = MIN(TMP, MAX(-TMP, P33))
      P_k = 0.5 * (P11 + P22 + P33)
    ENDIF
    TMP   = TwoThird * e_k
    eps11 = TMP
    eps22 = TMP
    eps33 = TMP
    eps12 = 0.0
    eps23 = 0.0
    eps31 = 0.0
    ! e / k (時間スケールリミッター) -----------------------------------
    AK   = 0.5 * (uu(I,J,K) + vv(I,J,K) + ww(I,J,K))
    SS2  = 2.0 * (UX**2 + VY**2 + WZ**2) &
    &    + (UY + VX)**2 + (VZ + WY)**2 + (WX + UZ)**2 &
    &    - TwoThird * (UX + VY + WZ)**2
    nu   = AMU(I,J,K) / RHO(I,J,K)
    AEAK = MIN(MAX(SQRT(SS2), SQRT(e_k / nu)), e_k / MAX(ZERO, AK))
    ! 再分配項 ---------------------------------------------------------
    TMP1    = - C_1 * AEAK
    phi1_11 = TMP1 * (uu(I,J,K) - TwoThird * AK)
    phi1_22 = TMP1 * (vv(I,J,K) - TwoThird * AK)
    phi1_33 = TMP1 * (ww(I,J,K) - TwoThird * AK)
    phi1_12 = TMP1 * uv(I,J,K)
    phi1_23 = TMP1 * vw(I,J,K)
    phi1_31 = TMP1 * wu(I,J,K)
    phi2_11 = - C_2 * (P11 - TwoThird * P_k)
    phi2_22 = - C_2 * (P22 - TwoThird * P_k)
    phi2_33 = - C_2 * (P33 - TwoThird * P_k)
    phi2_12 = - C_2 * P12
    phi2_23 = - C_2 * P23
    phi2_31 = - C_2 * P31
    IF(C_l * e_k * d(I,J,K) .GT. ZERO) THEN
      TMPw = AK**1.5 / MAX(ZERO, (C_l * e_k * d(I,J,K)))
      ! 壁面反射項のダンピング関数の制限
      ! TMPw = MIN(1.0, MAX(0.0, TMPw))
    ELSE
      TMPw = 0.0
    ENDIF
    TMPw1   = C_w1 * TMPw * AEAK
    TMPw2   = C_w2 * TMPw
    TMPw1km = uu(I,J,K) * n1(I,J,K) * n1(I,J,K) &
    &       + vv(I,J,K) * n2(I,J,K) * n2(I,J,K) &
    &       + ww(I,J,K) * n3(I,J,K) * n3(I,J,K) &
    &       + 2.0 * ( &
    &         uv(I,J,K) * n1(I,J,K) * n2(I,J,K) &
    &       + vw(I,J,K) * n2(I,J,K) * n3(I,J,K) &
    &       + wu(I,J,K) * n3(I,J,K) * n1(I,J,K) &
    &       )
    TMPw2km = phi2_11 * n1(I,J,K) * n1(I,J,K) &
    &       + phi2_22 * n2(I,J,K) * n2(I,J,K) &
    &       + phi2_33 * n3(I,J,K) * n3(I,J,K) &
    &       + 2.0 * ( &
    &         phi2_12 * n1(I,J,K) * n2(I,J,K) &
    &       + phi2_23 * n2(I,J,K) * n3(I,J,K) &
    &       + phi2_31 * n3(I,J,K) * n1(I,J,K) &
    &       )
    phiw1_11 = TMPw1 * (TMPw1km - 3.0 * n1(I,J,K) &
    &          * ( uu(I,J,K) * n1(I,J,K) &
    &            + uv(I,J,K) * n2(I,J,K) &
    &            + wu(I,J,K) * n3(I,J,K) ) &
    &        )
    phiw1_22 = TMPw1 * (TMPw1km - 3.0 * n2(I,J,K) &
    &          * ( uv(I,J,K) * n1(I,J,K) &
    &            + vv(I,J,K) * n2(I,J,K) &
    &            + vw(I,J,K) * n3(I,J,K) ) &
    &        )
    phiw1_33 = TMPw1 * (TMPw1km - 3.0 * n3(I,J,K) &
    &          * ( wu(I,J,K) * n1(I,J,K) &
    &            + vw(I,J,K) * n2(I,J,K) &
    &            + ww(I,J,K) * n3(I,J,K) ) &
    &        )
    phiw1_12 = - 1.5 * TMPw1 * ( &
    &          n1(I,J,K) * ( uv(I,J,K) * n1(I,J,K) &
    &                      + vv(I,J,K) * n2(I,J,K) &
    &                      + vw(I,J,K) * n3(I,J,K) ) &
    &        + n2(I,J,K) * ( uu(I,J,K) * n1(I,J,K) &
    &                      + uv(I,J,K) * n2(I,J,K) &
    &                      + wu(I,J,K) * n3(I,J,K) ) &
    &        )
    phiw1_23 = - 1.5 * TMPw1 * ( &
    &          n2(I,J,K) * ( wu(I,J,K) * n1(I,J,K) &
    &                      + vw(I,J,K) * n2(I,J,K) &
    &                      + ww(I,J,K) * n3(I,J,K) ) &
    &        + n3(I,J,K) * ( uv(I,J,K) * n1(I,J,K) &
    &                      + vv(I,J,K) * n2(I,J,K) &
    &                      + vw(I,J,K) * n3(I,J,K) ) &
    &        )
    phiw1_31 = - 1.5 * TMPw1 * ( &
    &          n3(I,J,K) * ( uu(I,J,K) * n1(I,J,K) &
    &                      + uv(I,J,K) * n2(I,J,K) &
    &                      + wu(I,J,K) * n3(I,J,K) ) &
    &        + n1(I,J,K) * ( wu(I,J,K) * n1(I,J,K) &
    &                      + vw(I,J,K) * n2(I,J,K) &
    &                      + ww(I,J,K) * n3(I,J,K) ) &
    &        )
    phiw2_11 = TMPw2 * (TMPw2km - 3.0 * n1(I,J,K) &
    &          * ( phi2_11 * n1(I,J,K) &
    &            + phi2_12 * n2(I,J,K) &
    &            + phi2_31 * n3(I,J,K) ) &
    &        )
    phiw2_22 = TMPw2 * (TMPw2km - 3.0 * n2(I,J,K) &
    &          * ( phi2_12 * n1(I,J,K) &
    &            + phi2_22 * n2(I,J,K) &
    &            + phi2_23 * n3(I,J,K) ) &
    &        )
    phiw2_33 = TMPw2 * (TMPw2km - 3.0 * n3(I,J,K) &
    &          * ( phi2_31 * n1(I,J,K) &
    &            + phi2_23 * n2(I,J,K) &
    &            + phi2_33 * n3(I,J,K) ) &
    &        )
    phiw2_12 = - 1.5 * TMPw2 * ( &
    &          n1(I,J,K) * ( phi2_12 * n1(I,J,K) &
    &                      + phi2_22 * n2(I,J,K) &
    &                      + phi2_23 * n3(I,J,K) ) &
    &        + n2(I,J,K) * ( phi2_11 * n1(I,J,K) &
    &                      + phi2_12 * n2(I,J,K) &
    &                      + phi2_31 * n3(I,J,K) ) &
    &        )
    phiw2_23 = - 1.5 * TMPw2 * ( &
    &          n2(I,J,K) * ( phi2_31 * n1(I,J,K) &
    &                      + phi2_23 * n2(I,J,K) &
    &                      + phi2_33 * n3(I,J,K) ) &
    &        + n3(I,J,K) * ( phi2_12 * n1(I,J,K) &
    &                      + phi2_22 * n2(I,J,K) &
    &                      + phi2_23 * n3(I,J,K) ) &
    &        )
    phiw2_31 = - 1.5 * TMPw2 * ( &
    &          n3(I,J,K) * ( phi2_11 * n1(I,J,K) &
    &                      + phi2_12 * n2(I,J,K) &
    &                      + phi2_31 * n3(I,J,K) ) &
    &        + n1(I,J,K) * ( phi2_31 * n1(I,J,K) &
    &                      + phi2_23 * n2(I,J,K) &
    &                      + phi2_33 * n3(I,J,K) ) &
    &        )
    phi11 = phi1_11 + phi2_11 + phiw1_11 + phiw2_11
    phi22 = phi1_22 + phi2_22 + phiw1_22 + phiw2_22
    phi33 = phi1_33 + phi2_33 + phiw1_33 + phiw2_33
    phi12 = phi1_12 + phi2_12 + phiw1_12 + phiw2_12
    phi23 = phi1_23 + phi2_23 + phiw1_23 + phiw2_23
    phi31 = phi1_31 + phi2_31 + phiw1_31 + phiw2_31
    ! 渦粘性係数の計算 -------------------------------------------------
    IF(e_k .GT. ZERO) THEN
      AMUT(I,J,K) = 0.09 * AK**2 / e_k * RHO(I,J,K)
    ENDIF
    ! 生成項、圧力再分配項、散逸項の和 ---------------------------------
    QH1 = RHO(I,J,K) / AJA(I,J,K)
    DQP(I,J,K, 6) = QH1 * (P11 - eps11 + phi11)
    DQP(I,J,K, 7) = QH1 * (P22 - eps22 + phi22)
    DQP(I,J,K, 8) = QH1 * (P33 - eps33 + phi33)
    DQP(I,J,K, 9) = QH1 * (P12 - eps12 + phi12)
    DQP(I,J,K,10) = QH1 * (P23 - eps23 + phi23)
    DQP(I,J,K,11) = QH1 * (P31 - eps31 + phi31)
    DQP(I,J,K,12) = QH1 * AEAK * (C_e1 * P_k - C_e2 * e_k)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH = 0.0
  SH = 0.0
  TH = 0.0
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL DIFFXI
  CALL DIFFET
  CALL DIFFZE
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I,J-1,K) .GT. 0.0) .AND. &
  &   (RHO(I-1,J,K) .GT. 0.0) .AND. (RHO(I,J,K)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J,K) .GT. 0.0) .AND. (RHO(I,J+1,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,K,L) = (-RH(I-1,J  ,K  ,L) + RH(I,J,K,L) ) &
    &            + (-SH(I  ,J-1,K  ,L) + SH(I,J,K,L) ) &
    &            + (-TH(I  ,J  ,K-1,L) + TH(I,J,K,L) )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項の計算(xi方向差分用)                                    ****
!***********************************************************************
SUBROUTINE DIFFXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: uuXI, uuET, uuZE, uuX, uuY, uuZ
  REAL    :: vvXI, vvET, vvZE, vvX, vvY, vvZ
  REAL    :: wwXI, wwET, wwZE, wwX, wwY, wwZ
  REAL    :: uvXI, uvET, uvZE, uvX, uvY, uvZ
  REAL    :: vwXI, vwET, vwZE, vwX, vwY, vwZ
  REAL    :: wuXI, wuET, wuZE, wuX, wuY, wuZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu
  REAL    :: S11, S22, S33, S12, S23, S31, SS
  REAL    :: C_p, AKAE
  REAL    :: R02M, S02M, T02M, R02T, S02T, T02T
  REAL    :: R03M, S03M, T03M, R03T, S03T, T03T
  REAL    :: R04M, S04M, T04M, R04T, S04T, T04T
  REAL    :: R05M, S05M, T05M, R05T, S05T, T05T
  REAL    :: R06M, S06M, T06M, R06T, S06T, T06T
  REAL    :: R07M, S07M, T07M, R07T, S07T, T07T
  REAL    :: R08M, S08M, T08M, R08T, S08T, T08T
  REAL    :: R09M, S09M, T09M, R09T, S09T, T09T
  REAL    :: R10M, S10M, T10M, R10T, S10T, T10T
  REAL    :: R11M, S11M, T11M, R11T, S11T, T11T
  REAL    :: R12M, S12M, T12M, R12T, S12T, T12T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  uuXI, uuET, uuZE, uuX, uuY, uuZ, &
  !$  vvXI, vvET, vvZE, vvX, vvY, vvZ, &
  !$  wwXI, wwET, wwZE, wwX, wwY, wwZ, &
  !$  uvXI, uvET, uvZE, uvX, uvY, uvZ, &
  !$  vwXI, vwET, vwZE, vwX, vwY, vwZ, &
  !$  wuXI, wuET, wuZE, wuX, wuY, wuZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, S11, S22, S33, S12, S23, S31, SS, C_p, AKAE, &
  !$  R02M, S02M, T02M, R02T, S02T, T02T, &
  !$  R03M, S03M, T03M, R03T, S03T, T03T, &
  !$  R04M, S04M, T04M, R04T, S04T, T04T, &
  !$  R05M, S05M, T05M, R05T, S05T, T05T, &
  !$  R06M, S06M, T06M, R06T, S06T, T06T, &
  !$  R07M, S07M, T07M, R07T, S07T, T07T, &
  !$  R08M, S08M, T08M, R08T, S08T, T08T, &
  !$  R09M, S09M, T09M, R09T, S09T, T09T, &
  !$  R10M, S10M, T10M, R10T, S10T, T10T, &
  !$  R11M, S11M, T11M, R11T, S11T, T11T, &
  !$  R12M, S12M, T12M, R12T, S12T, T12T, &
  !$  QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF((RHO(I,J,K) .GT. 0.0) .AND. (RHO(I+1,J,K) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I+1,J,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I+1,J,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I+1,J,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I+1,J,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I+1,J,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I+1,J,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I+1,J,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I+1,J,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I+1,J,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I+1,J,K))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J,K) + RHO(I+1,J,K))
    UM   = 0.5 * (  U(I,J,K) +   U(I+1,J,K))
    VM   = 0.5 * (  V(I,J,K) +   V(I+1,J,K))
    WM   = 0.5 * (  W(I,J,K) +   W(I+1,J,K))
    uuM  = 0.5 * ( uu(I,J,K) +  uu(I+1,J,K))
    vvM  = 0.5 * ( vv(I,J,K) +  vv(I+1,J,K))
    wwM  = 0.5 * ( ww(I,J,K) +  ww(I+1,J,K))
    uvM  = 0.5 * ( uv(I,J,K) +  uv(I+1,J,K))
    vwM  = 0.5 * ( vw(I,J,K) +  vw(I+1,J,K))
    wuM  = 0.5 * ( wu(I,J,K) +  wu(I+1,J,K))
    EPSM = 0.5 * (EPS(I,J,K) + EPS(I+1,J,K))
    AMUM = 0.5 * (AMU(I,J,K) + AMU(I+1,J,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = - U(I,J,K) + U(I+1,J,K)
    UET = 0.25 * ( - U(I,J-1,K) - U(I+1,J-1,K) &
    &              + U(I,J+1,K) + U(I+1,J+1,K) )
    UZE = 0.25 * ( - U(I,J,K-1) - U(I+1,J,K-1) &
    &              + U(I,J,K+1) + U(I+1,J,K+1) )
    VXI = - V(I,J,K) + V(I+1,J,K)
    VET = 0.25 * ( - V(I,J-1,K) - V(I+1,J-1,K) &
    &              + V(I,J+1,K) + V(I+1,J+1,K) )
    VZE = 0.25 * ( - V(I,J,K-1) - V(I+1,J,K-1) &
    &              + V(I,J,K+1) + V(I+1,J,K+1) )
    WXI = - W(I,J,K) + W(I+1,J,K)
    WET = 0.25 * ( - W(I,J-1,K) - W(I+1,J-1,K) &
    &              + W(I,J+1,K) + W(I+1,J+1,K) )
    WZE = 0.25 * ( - W(I,J,K-1) - W(I+1,J,K-1) &
    &              + W(I,J,K+1) + W(I+1,J,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = - T(I,J,K) + T(I+1,J,K)
    TET = 0.25 * ( - T(I,J-1,K) - T(I+1,J-1,K) &
    &              + T(I,J+1,K) + T(I+1,J+1,K) )
    TZE = 0.25 * ( - T(I,J,K-1) - T(I+1,J,K-1) &
    &              + T(I,J,K+1) + T(I+1,J,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = - uu(I,J,K) + uu(I+1,J,K)
    uuET = 0.25 * ( - uu(I,J-1,K) - uu(I+1,J-1,K) &
    &               + uu(I,J+1,K) + uu(I+1,J+1,K) )
    uuZE = 0.25 * ( - uu(I,J,K-1) - uu(I+1,J,K-1) &
    &               + uu(I,J,K+1) + uu(I+1,J,K+1) )
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM + uuZE * ZEXM
    uuY = uuXI * XIYM + uuET * ETYM + uuZE * ZEYM
    uuZ = uuXI * XIZM + uuET * ETZM + uuZE * ZEZM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = - vv(I,J,K) + vv(I+1,J,K)
    vvET = 0.25 * ( - vv(I,J-1,K) - vv(I+1,J-1,K) &
    &               + vv(I,J+1,K) + vv(I+1,J+1,K) )
    vvZE = 0.25 * ( - vv(I,J,K-1) - vv(I+1,J,K-1) &
    &               + vv(I,J,K+1) + vv(I+1,J,K+1) )
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM + vvZE * ZEXM
    vvY = vvXI * XIYM + vvET * ETYM + vvZE * ZEYM
    vvZ = vvXI * XIZM + vvET * ETZM + vvZE * ZEZM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = - ww(I,J,K) + ww(I+1,J,K)
    wwET = 0.25 * ( - ww(I,J-1,K) - ww(I+1,J-1,K) &
    &               + ww(I,J+1,K) + ww(I+1,J+1,K) )
    wwZE = 0.25 * ( - ww(I,J,K-1) - ww(I+1,J,K-1) &
    &               + ww(I,J,K+1) + ww(I+1,J,K+1) )
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM + wwZE * ZEXM
    wwY = wwXI * XIYM + wwET * ETYM + wwZE * ZEYM
    wwZ = wwXI * XIZM + wwET * ETZM + wwZE * ZEZM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = - uv(I,J,K) + uv(I+1,J,K)
    uvET = 0.25 * ( - uv(I,J-1,K) - uv(I+1,J-1,K) &
    &               + uv(I,J+1,K) + uv(I+1,J+1,K) )
    uvZE = 0.25 * ( - uv(I,J,K-1) - uv(I+1,J,K-1) &
    &               + uv(I,J,K+1) + uv(I+1,J,K+1) )
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM + uvZE * ZEXM
    uvY = uvXI * XIYM + uvET * ETYM + uvZE * ZEYM
    uvZ = uvXI * XIZM + uvET * ETZM + uvZE * ZEZM
    ! vw の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vwXI = - vw(I,J,K) + vw(I+1,J,K)
    vwET = 0.25 * ( - vw(I,J-1,K) - vw(I+1,J-1,K) &
    &               + vw(I,J+1,K) + vw(I+1,J+1,K) )
    vwZE = 0.25 * ( - vw(I,J,K-1) - vw(I+1,J,K-1) &
    &               + vw(I,J,K+1) + vw(I+1,J,K+1) )
    ! 物理空間方向一階微分
    vwX = vwXI * XIXM + vwET * ETXM + vwZE * ZEXM
    vwY = vwXI * XIYM + vwET * ETYM + vwZE * ZEYM
    vwZ = vwXI * XIZM + vwET * ETZM + vwZE * ZEZM
    ! wu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wuXI = - wu(I,J,K) + wu(I+1,J,K)
    wuET = 0.25 * ( - wu(I,J-1,K) - wu(I+1,J-1,K) &
    &               + wu(I,J+1,K) + wu(I+1,J+1,K) )
    wuZE = 0.25 * ( - wu(I,J,K-1) - wu(I+1,J,K-1) &
    &               + wu(I,J,K+1) + wu(I+1,J,K+1) )
    ! 物理空間方向一階微分
    wuX = wuXI * XIXM + wuET * ETXM + wuZE * ZEXM
    wuY = wuXI * XIYM + wuET * ETYM + wuZE * ZEYM
    wuZ = wuXI * XIZM + wuET * ETZM + wuZE * ZEZM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = - EPS(I,J,K) + EPS(I+1,J,K)
    AEET = 0.25 * ( - EPS(I,J-1,K) - EPS(I+1,J-1,K) &
    &               + EPS(I,J+1,K) + EPS(I+1,J+1,K) )
    AEZE = 0.25 * ( - EPS(I,J,K-1) - EPS(I+1,J,K-1) &
    &               + EPS(I,J,K+1) + EPS(I+1,J,K+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY + WZ) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 = 2.0 * WZ - SS
    S12 = UY + VX
    S23 = VZ + WY
    S31 = WX + UZ
    ! U の拡散 ---------------------------------------------------------
    R02M = nu * S11
    S02M = nu * S12
    T02M = nu * S31
    R02T =-uuM
    S02T =-uvM
    T02T =-wuM
    ! V の拡散 ---------------------------------------------------------
    R03M = nu * S12
    S03M = nu * S22
    T03M = nu * S23
    R03T =-uvM
    S03T =-vvM
    T03T =-vwM
    ! W の拡散 ---------------------------------------------------------
    R04M = nu * S31
    S04M = nu * S23
    T04M = nu * S33
    R04T =-wuM
    S04T =-vwM
    T04T =-wwM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R05M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM + S31 * WM)
    S05M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM + S23 * WM)
    T05M = C_p * nu / PR * TZ + nu * (S31 * UM + S23 * VM + S33 * WM)
    R05T = AKAE * C_se * C_p * (uuM * TX + uvM * TY + wuM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uuM * (uuX + vvX + wwX) &
    &    + uvM * (uuY + vvY + wwY) &
    &    + wuM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uuM * UM + uvM * VM + wuM * WM)
    S05T = AKAE * C_se * C_p * (uvM * TX + vvM * TY + vwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uvM * (uuX + vvX + wwX) &
    &    + vvM * (uuY + vvY + wwY) &
    &    + vwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uvM * UM + vvM * VM + vwM * WM)
    T05T = AKAE * C_se * C_p * (wuM * TX + vwM * TY + wwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      wuM * (uuX + vvX + wwX) &
    &    + vwM * (uuY + vvY + wwY) &
    &    + wwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (wuM * UM + vwM * VM + wwM * WM)
    ! uu の拡散 --------------------------------------------------------
    R06M = nu * uuX
    S06M = nu * uuY
    T06M = nu * uuZ
    R06T = C_s * AKAE * (uuM * uuX + uvM * uuY + wuM * uuZ)
    S06T = C_s * AKAE * (uvM * uuX + vvM * uuY + vwM * uuZ)
    T06T = C_s * AKAE * (wuM * uuX + vwM * uuY + wwM * uuZ)
    ! vv の拡散 --------------------------------------------------------
    R07M = nu * vvX
    S07M = nu * vvY
    T07M = nu * vvZ
    R07T = C_s * AKAE * (uuM * vvX + uvM * vvY + wuM * vvZ)
    S07T = C_s * AKAE * (uvM * vvX + vvM * vvY + vwM * vvZ)
    T07T = C_s * AKAE * (wuM * vvX + vwM * vvY + wwM * vvZ)
    ! ww の拡散 --------------------------------------------------------
    R08M = nu * wwX
    S08M = nu * wwY
    T08M = nu * wwZ
    R08T = C_s * AKAE * (uuM * wwX + uvM * wwY + wuM * wwZ)
    S08T = C_s * AKAE * (uvM * wwX + vvM * wwY + vwM * wwZ)
    T08T = C_s * AKAE * (wuM * wwX + vwM * wwY + wwM * wwZ)
    ! uv の拡散 --------------------------------------------------------
    R09M = nu * uvX
    S09M = nu * uvY
    T09M = nu * uvZ
    R09T = C_s * AKAE * (uuM * uvX + uvM * uvY + wuM * uvZ)
    S09T = C_s * AKAE * (uvM * uvX + vvM * uvY + vwM * uvZ)
    T09T = C_s * AKAE * (wuM * uvX + vwM * uvY + wwM * uvZ)
    ! vw の拡散 --------------------------------------------------------
    R10M = nu * vwX
    S10M = nu * vwY
    T10M = nu * vwZ
    R10T = C_s * AKAE * (uuM * vwX + uvM * vwY + wuM * vwZ)
    S10T = C_s * AKAE * (uvM * vwX + vvM * vwY + vwM * vwZ)
    T10T = C_s * AKAE * (wuM * vwX + vwM * vwY + wwM * vwZ)
    ! wu の拡散 --------------------------------------------------------
    R11M = nu * wuX
    S11M = nu * wuY
    T11M = nu * wuZ
    R11T = C_s * AKAE * (uuM * wuX + uvM * wuY + wuM * wuZ)
    S11T = C_s * AKAE * (uvM * wuX + vvM * wuY + vwM * wuZ)
    T11T = C_s * AKAE * (wuM * wuX + vwM * wuY + wwM * wuZ)
    ! epsilon の拡散 ---------------------------------------------------
    R12M = nu * AEX
    S12M = nu * AEY
    T12M = nu * AEZ
    R12T = C_e * AKAE * (uuM * AEX + uvM * AEY + wuM * AEZ)
    S12T = C_e * AKAE * (uvM * AEX + vvM * AEY + vwM * AEZ)
    T12T = C_e * AKAE * (wuM * AEX + vwM * AEY + wwM * AEZ)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    RH(I,J,K, 1) = 0.0
    RH(I,J,K, 2) = QH1 * ( (R02M + R02T) * XIXM &
    &                    + (S02M + S02T) * XIYM &
    &                    + (T02M + T02T) * XIZM )
    RH(I,J,K, 3) = QH1 * ( (R03M + R03T) * XIXM &
    &                    + (S03M + S03T) * XIYM &
    &                    + (T03M + T03T) * XIZM )
    RH(I,J,K, 4) = QH1 * ( (R04M + R04T) * XIXM &
    &                    + (S04M + S04T) * XIYM &
    &                    + (T04M + T04T) * XIZM )
    RH(I,J,K, 5) = QH1 * ( (R05M + R05T) * XIXM &
    &                    + (S05M + S05T) * XIYM &
    &                    + (T05M + T05T) * XIZM )
    RH(I,J,K, 6) = QH1 * ( (R06M + R06T) * XIXM &
    &                    + (S06M + S06T) * XIYM &
    &                    + (T06M + T06T) * XIZM )
    RH(I,J,K, 7) = QH1 * ( (R07M + R07T) * XIXM &
    &                    + (S07M + S07T) * XIYM &
    &                    + (T07M + T07T) * XIZM )
    RH(I,J,K, 8) = QH1 * ( (R08M + R08T) * XIXM &
    &                    + (S08M + S08T) * XIYM &
    &                    + (T08M + T08T) * XIZM )
    RH(I,J,K, 9) = QH1 * ( (R09M + R09T) * XIXM &
    &                    + (S09M + S09T) * XIYM &
    &                    + (T09M + T09T) * XIZM )
    RH(I,J,K,10) = QH1 * ( (R10M + R10T) * XIXM &
    &                    + (S10M + S10T) * XIYM &
    &                    + (T10M + T10T) * XIZM )
    RH(I,J,K,11) = QH1 * ( (R11M + R11T) * XIXM &
    &                    + (S11M + S11T) * XIYM &
    &                    + (T11M + T11T) * XIZM )
    RH(I,J,K,12) = QH1 * ( (R12M + R12T) * XIXM &
    &                    + (S12M + S12T) * XIYM &
    &                    + (T12M + T12T) * XIZM )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFXI
!***********************************************************************
!**** 拡散項の計算(eta方向差分用)                                   ****
!***********************************************************************
SUBROUTINE DIFFET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: uuXI, uuET, uuZE, uuX, uuY, uuZ
  REAL    :: vvXI, vvET, vvZE, vvX, vvY, vvZ
  REAL    :: wwXI, wwET, wwZE, wwX, wwY, wwZ
  REAL    :: uvXI, uvET, uvZE, uvX, uvY, uvZ
  REAL    :: vwXI, vwET, vwZE, vwX, vwY, vwZ
  REAL    :: wuXI, wuET, wuZE, wuX, wuY, wuZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu
  REAL    :: S11, S22, S33, S12, S23, S31, SS
  REAL    :: C_p, AKAE
  REAL    :: R02M, S02M, T02M, R02T, S02T, T02T
  REAL    :: R03M, S03M, T03M, R03T, S03T, T03T
  REAL    :: R04M, S04M, T04M, R04T, S04T, T04T
  REAL    :: R05M, S05M, T05M, R05T, S05T, T05T
  REAL    :: R06M, S06M, T06M, R06T, S06T, T06T
  REAL    :: R07M, S07M, T07M, R07T, S07T, T07T
  REAL    :: R08M, S08M, T08M, R08T, S08T, T08T
  REAL    :: R09M, S09M, T09M, R09T, S09T, T09T
  REAL    :: R10M, S10M, T10M, R10T, S10T, T10T
  REAL    :: R11M, S11M, T11M, R11T, S11T, T11T
  REAL    :: R12M, S12M, T12M, R12T, S12T, T12T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  uuXI, uuET, uuZE, uuX, uuY, uuZ, &
  !$  vvXI, vvET, vvZE, vvX, vvY, vvZ, &
  !$  wwXI, wwET, wwZE, wwX, wwY, wwZ, &
  !$  uvXI, uvET, uvZE, uvX, uvY, uvZ, &
  !$  vwXI, vwET, vwZE, vwX, vwY, vwZ, &
  !$  wuXI, wuET, wuZE, wuX, wuY, wuZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, S11, S22, S33, S12, S23, S31, SS, C_p, AKAE, &
  !$  R02M, S02M, T02M, R02T, S02T, T02T, &
  !$  R03M, S03M, T03M, R03T, S03T, T03T, &
  !$  R04M, S04M, T04M, R04T, S04T, T04T, &
  !$  R05M, S05M, T05M, R05T, S05T, T05T, &
  !$  R06M, S06M, T06M, R06T, S06T, T06T, &
  !$  R07M, S07M, T07M, R07T, S07T, T07T, &
  !$  R08M, S08M, T08M, R08T, S08T, T08T, &
  !$  R09M, S09M, T09M, R09T, S09T, T09T, &
  !$  R10M, S10M, T10M, R10T, S10T, T10T, &
  !$  R11M, S11M, T11M, R11T, S11T, T11T, &
  !$  R12M, S12M, T12M, R12T, S12T, T12T, &
  !$  QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF((RHO(I,J,K) .GT. 0.0) .AND. (RHO(I,J+1,K) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J+1,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J+1,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J+1,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J+1,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J+1,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J+1,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J+1,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J+1,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J+1,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J+1,K))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J,K) + RHO(I,J+1,K))
    UM   = 0.5 * (  U(I,J,K) +   U(I,J+1,K))
    VM   = 0.5 * (  V(I,J,K) +   V(I,J+1,K))
    WM   = 0.5 * (  W(I,J,K) +   W(I,J+1,K))
    uuM  = 0.5 * ( uu(I,J,K) +  uu(I,J+1,K))
    vvM  = 0.5 * ( vv(I,J,K) +  vv(I,J+1,K))
    wwM  = 0.5 * ( ww(I,J,K) +  ww(I,J+1,K))
    uvM  = 0.5 * ( uv(I,J,K) +  uv(I,J+1,K))
    vwM  = 0.5 * ( vw(I,J,K) +  vw(I,J+1,K))
    wuM  = 0.5 * ( wu(I,J,K) +  wu(I,J+1,K))
    EPSM = 0.5 * (EPS(I,J,K) + EPS(I,J+1,K))
    AMUM = 0.5 * (AMU(I,J,K) + AMU(I,J+1,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J+1,K) &
    &              + U(I+1,J,K) + U(I+1,J+1,K) )
    UET = - U(I,J,K) + U(I,J+1,K)
    UZE = 0.25 * ( - U(I,J,K-1) - U(I,J+1,K-1) &
    &              + U(I,J,K+1) + U(I,J+1,K+1) )
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J+1,K) &
    &              + V(I+1,J,K) + V(I+1,J+1,K) )
    VET = - V(I,J,K) + V(I,J+1,K)
    VZE = 0.25 * ( - V(I,J,K-1) - V(I,J+1,K-1) &
    &              + V(I,J,K+1) + V(I,J+1,K+1) )
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J+1,K) &
    &              + W(I+1,J,K) + W(I+1,J+1,K) )
    WET = - W(I,J,K) + W(I,J+1,K)
    WZE = 0.25 * ( - W(I,J,K-1) - W(I,J+1,K-1) &
    &              + W(I,J,K+1) + W(I,J+1,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J+1,K) &
    &              + T(I+1,J,K) + T(I+1,J+1,K) )
    TET = - T(I,J,K) + T(I,J+1,K)
    TZE = 0.25 * ( - T(I,J,K-1) - T(I,J+1,K-1) &
    &              + T(I,J,K+1) + T(I,J+1,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = 0.25 * ( - uu(I-1,J,K) - uu(I-1,J+1,K) &
    &               + uu(I+1,J,K) + uu(I+1,J+1,K) )
    uuET = - uu(I,J,K) + uu(I,J+1,K)
    uuZE = 0.25 * ( - uu(I,J,K-1) - uu(I,J+1,K-1) &
    &               + uu(I,J,K+1) + uu(I,J+1,K+1) )
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM + uuZE * ZEXM
    uuY = uuXI * XIYM + uuET * ETYM + uuZE * ZEYM
    uuZ = uuXI * XIZM + uuET * ETZM + uuZE * ZEZM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = 0.25 * ( - vv(I-1,J,K) - vv(I-1,J+1,K) &
    &               + vv(I+1,J,K) + vv(I+1,J+1,K) )
    vvET = - vv(I,J,K) + vv(I,J+1,K)
    vvZE = 0.25 * ( - vv(I,J,K-1) - vv(I,J+1,K-1) &
    &               + vv(I,J,K+1) + vv(I,J+1,K+1) )
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM + vvZE * ZEXM
    vvY = vvXI * XIYM + vvET * ETYM + vvZE * ZEYM
    vvZ = vvXI * XIZM + vvET * ETZM + vvZE * ZEZM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = 0.25 * ( - ww(I-1,J,K) - ww(I-1,J+1,K) &
    &               + ww(I+1,J,K) + ww(I+1,J+1,K) )
    wwET = - ww(I,J,K) + ww(I,J+1,K)
    wwZE = 0.25 * ( - ww(I,J,K-1) - ww(I,J+1,K-1) &
    &               + ww(I,J,K+1) + ww(I,J+1,K+1) )
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM + wwZE * ZEXM
    wwY = wwXI * XIYM + wwET * ETYM + wwZE * ZEYM
    wwZ = wwXI * XIZM + wwET * ETZM + wwZE * ZEZM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = 0.25 * ( - uv(I-1,J,K) - uv(I-1,J+1,K) &
    &               + uv(I+1,J,K) + uv(I+1,J+1,K) )
    uvET = - uv(I,J,K) + uv(I,J+1,K)
    uvZE = 0.25 * ( - uv(I,J,K-1) - uv(I,J+1,K-1) &
    &               + uv(I,J,K+1) + uv(I,J+1,K+1) )
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM + uvZE * ZEXM
    uvY = uvXI * XIYM + uvET * ETYM + uvZE * ZEYM
    uvZ = uvXI * XIZM + uvET * ETZM + uvZE * ZEZM
    ! vw の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vwXI = 0.25 * ( - vw(I-1,J,K) - vw(I-1,J+1,K) &
    &               + vw(I+1,J,K) + vw(I+1,J+1,K) )
    vwET = - vw(I,J,K) + vw(I,J+1,K)
    vwZE = 0.25 * ( - vw(I,J,K-1) - vw(I,J+1,K-1) &
    &               + vw(I,J,K+1) + vw(I,J+1,K+1) )
    ! 物理空間方向一階微分
    vwX = vwXI * XIXM + vwET * ETXM + vwZE * ZEXM
    vwY = vwXI * XIYM + vwET * ETYM + vwZE * ZEYM
    vwZ = vwXI * XIZM + vwET * ETZM + vwZE * ZEZM
    ! wu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wuXI = 0.25 * ( - wu(I-1,J,K) - wu(I-1,J+1,K) &
    &               + wu(I+1,J,K) + wu(I+1,J+1,K) )
    wuET = - wu(I,J,K) + wu(I,J+1,K)
    wuZE = 0.25 * ( - wu(I,J,K-1) - wu(I,J+1,K-1) &
    &               + wu(I,J,K+1) + wu(I,J+1,K+1) )
    ! 物理空間方向一階微分
    wuX = wuXI * XIXM + wuET * ETXM + wuZE * ZEXM
    wuY = wuXI * XIYM + wuET * ETYM + wuZE * ZEYM
    wuZ = wuXI * XIZM + wuET * ETZM + wuZE * ZEZM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPS(I-1,J,K) - EPS(I-1,J+1,K) &
    &               + EPS(I+1,J,K) + EPS(I+1,J+1,K) )
    AEET = - EPS(I,J,K) + EPS(I,J+1,K)
    AEZE = 0.25 * ( - EPS(I,J,K-1) - EPS(I,J+1,K-1) &
    &               + EPS(I,J,K+1) + EPS(I,J+1,K+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY + WZ) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 = 2.0 * WZ - SS
    S12 = UY + VX
    S23 = VZ + WY
    S31 = WX + UZ
    ! U の拡散 ---------------------------------------------------------
    R02M = nu * S11
    S02M = nu * S12
    T02M = nu * S31
    R02T =-uuM
    S02T =-uvM
    T02T =-wuM
    ! V の拡散 ---------------------------------------------------------
    R03M = nu * S12
    S03M = nu * S22
    T03M = nu * S23
    R03T =-uvM
    S03T =-vvM
    T03T =-vwM
    ! W の拡散 ---------------------------------------------------------
    R04M = nu * S31
    S04M = nu * S23
    T04M = nu * S33
    R04T =-wuM
    S04T =-vwM
    T04T =-wwM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R05M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM + S31 * WM)
    S05M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM + S23 * WM)
    T05M = C_p * nu / PR * TZ + nu * (S31 * UM + S23 * VM + S33 * WM)
    R05T = AKAE * C_se * C_p * (uuM * TX + uvM * TY + wuM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uuM * (uuX + vvX + wwX) &
    &    + uvM * (uuY + vvY + wwY) &
    &    + wuM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uuM * UM + uvM * VM + wuM * WM)
    S05T = AKAE * C_se * C_p * (uvM * TX + vvM * TY + vwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uvM * (uuX + vvX + wwX) &
    &    + vvM * (uuY + vvY + wwY) &
    &    + vwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uvM * UM + vvM * VM + vwM * WM)
    T05T = AKAE * C_se * C_p * (wuM * TX + vwM * TY + wwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      wuM * (uuX + vvX + wwX) &
    &    + vwM * (uuY + vvY + wwY) &
    &    + wwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (wuM * UM + vwM * VM + wwM * WM)
    ! uu の拡散 --------------------------------------------------------
    R06M = nu * uuX
    S06M = nu * uuY
    T06M = nu * uuZ
    R06T = C_s * AKAE * (uuM * uuX + uvM * uuY + wuM * uuZ)
    S06T = C_s * AKAE * (uvM * uuX + vvM * uuY + vwM * uuZ)
    T06T = C_s * AKAE * (wuM * uuX + vwM * uuY + wwM * uuZ)
    ! vv の拡散 --------------------------------------------------------
    R07M = nu * vvX
    S07M = nu * vvY
    T07M = nu * vvZ
    R07T = C_s * AKAE * (uuM * vvX + uvM * vvY + wuM * vvZ)
    S07T = C_s * AKAE * (uvM * vvX + vvM * vvY + vwM * vvZ)
    T07T = C_s * AKAE * (wuM * vvX + vwM * vvY + wwM * vvZ)
    ! ww の拡散 --------------------------------------------------------
    R08M = nu * wwX
    S08M = nu * wwY
    T08M = nu * wwZ
    R08T = C_s * AKAE * (uuM * wwX + uvM * wwY + wuM * wwZ)
    S08T = C_s * AKAE * (uvM * wwX + vvM * wwY + vwM * wwZ)
    T08T = C_s * AKAE * (wuM * wwX + vwM * wwY + wwM * wwZ)
    ! uv の拡散 --------------------------------------------------------
    R09M = nu * uvX
    S09M = nu * uvY
    T09M = nu * uvZ
    R09T = C_s * AKAE * (uuM * uvX + uvM * uvY + wuM * uvZ)
    S09T = C_s * AKAE * (uvM * uvX + vvM * uvY + vwM * uvZ)
    T09T = C_s * AKAE * (wuM * uvX + vwM * uvY + wwM * uvZ)
    ! vw の拡散 --------------------------------------------------------
    R10M = nu * vwX
    S10M = nu * vwY
    T10M = nu * vwZ
    R10T = C_s * AKAE * (uuM * vwX + uvM * vwY + wuM * vwZ)
    S10T = C_s * AKAE * (uvM * vwX + vvM * vwY + vwM * vwZ)
    T10T = C_s * AKAE * (wuM * vwX + vwM * vwY + wwM * vwZ)
    ! wu の拡散 --------------------------------------------------------
    R11M = nu * wuX
    S11M = nu * wuY
    T11M = nu * wuZ
    R11T = C_s * AKAE * (uuM * wuX + uvM * wuY + wuM * wuZ)
    S11T = C_s * AKAE * (uvM * wuX + vvM * wuY + vwM * wuZ)
    T11T = C_s * AKAE * (wuM * wuX + vwM * wuY + wwM * wuZ)
    ! epsilon の拡散 ---------------------------------------------------
    R12M = nu * AEX
    S12M = nu * AEY
    T12M = nu * AEZ
    R12T = C_e * AKAE * (uuM * AEX + uvM * AEY + wuM * AEZ)
    S12T = C_e * AKAE * (uvM * AEX + vvM * AEY + vwM * AEZ)
    T12T = C_e * AKAE * (wuM * AEX + vwM * AEY + wwM * AEZ)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    SH(I,J,K, 1) = 0.0
    SH(I,J,K, 2) = QH1 * ( (R02M + R02T) * ETXM &
    &                    + (S02M + S02T) * ETYM &
    &                    + (T02M + T02T) * ETZM )
    SH(I,J,K, 3) = QH1 * ( (R03M + R03T) * ETXM &
    &                    + (S03M + S03T) * ETYM &
    &                    + (T03M + T03T) * ETZM )
    SH(I,J,K, 4) = QH1 * ( (R04M + R04T) * ETXM &
    &                    + (S04M + S04T) * ETYM &
    &                    + (T04M + T04T) * ETZM )
    SH(I,J,K, 5) = QH1 * ( (R05M + R05T) * ETXM &
    &                    + (S05M + S05T) * ETYM &
    &                    + (T05M + T05T) * ETZM )
    SH(I,J,K, 6) = QH1 * ( (R06M + R06T) * ETXM &
    &                    + (S06M + S06T) * ETYM &
    &                    + (T06M + T06T) * ETZM )
    SH(I,J,K, 7) = QH1 * ( (R07M + R07T) * ETXM &
    &                    + (S07M + S07T) * ETYM &
    &                    + (T07M + T07T) * ETZM )
    SH(I,J,K, 8) = QH1 * ( (R08M + R08T) * ETXM &
    &                    + (S08M + S08T) * ETYM &
    &                    + (T08M + T08T) * ETZM )
    SH(I,J,K, 9) = QH1 * ( (R09M + R09T) * ETXM &
    &                    + (S09M + S09T) * ETYM &
    &                    + (T09M + T09T) * ETZM )
    SH(I,J,K,10) = QH1 * ( (R10M + R10T) * ETXM &
    &                    + (S10M + S10T) * ETYM &
    &                    + (T10M + T10T) * ETZM )
    SH(I,J,K,11) = QH1 * ( (R11M + R11T) * ETXM &
    &                    + (S11M + S11T) * ETYM &
    &                    + (T11M + T11T) * ETZM )
    SH(I,J,K,12) = QH1 * ( (R12M + R12T) * ETXM &
    &                    + (S12M + S12T) * ETYM &
    &                    + (T12M + T12T) * ETZM )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFET
!***********************************************************************
!**** 拡散項の計算(zeta方向差分用)                                  ****
!***********************************************************************
SUBROUTINE DIFFZE
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: uuXI, uuET, uuZE, uuX, uuY, uuZ
  REAL    :: vvXI, vvET, vvZE, vvX, vvY, vvZ
  REAL    :: wwXI, wwET, wwZE, wwX, wwY, wwZ
  REAL    :: uvXI, uvET, uvZE, uvX, uvY, uvZ
  REAL    :: vwXI, vwET, vwZE, vwX, vwY, vwZ
  REAL    :: wuXI, wuET, wuZE, wuX, wuY, wuZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu
  REAL    :: S11, S22, S33, S12, S23, S31, SS
  REAL    :: C_p, AKAE
  REAL    :: R02M, S02M, T02M, R02T, S02T, T02T
  REAL    :: R03M, S03M, T03M, R03T, S03T, T03T
  REAL    :: R04M, S04M, T04M, R04T, S04T, T04T
  REAL    :: R05M, S05M, T05M, R05T, S05T, T05T
  REAL    :: R06M, S06M, T06M, R06T, S06T, T06T
  REAL    :: R07M, S07M, T07M, R07T, S07T, T07T
  REAL    :: R08M, S08M, T08M, R08T, S08T, T08T
  REAL    :: R09M, S09M, T09M, R09T, S09T, T09T
  REAL    :: R10M, S10M, T10M, R10T, S10T, T10T
  REAL    :: R11M, S11M, T11M, R11T, S11T, T11T
  REAL    :: R12M, S12M, T12M, R12T, S12T, T12T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  uuXI, uuET, uuZE, uuX, uuY, uuZ, &
  !$  vvXI, vvET, vvZE, vvX, vvY, vvZ, &
  !$  wwXI, wwET, wwZE, wwX, wwY, wwZ, &
  !$  uvXI, uvET, uvZE, uvX, uvY, uvZ, &
  !$  vwXI, vwET, vwZE, vwX, vwY, vwZ, &
  !$  wuXI, wuET, wuZE, wuX, wuY, wuZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, S11, S22, S33, S12, S23, S31, SS, C_p, AKAE, &
  !$  R02M, S02M, T02M, R02T, S02T, T02T, &
  !$  R03M, S03M, T03M, R03T, S03T, T03T, &
  !$  R04M, S04M, T04M, R04T, S04T, T04T, &
  !$  R05M, S05M, T05M, R05T, S05T, T05T, &
  !$  R06M, S06M, T06M, R06T, S06T, T06T, &
  !$  R07M, S07M, T07M, R07T, S07T, T07T, &
  !$  R08M, S08M, T08M, R08T, S08T, T08T, &
  !$  R09M, S09M, T09M, R09T, S09T, T09T, &
  !$  R10M, S10M, T10M, R10T, S10T, T10T, &
  !$  R11M, S11M, T11M, R11T, S11T, T11T, &
  !$  R12M, S12M, T12M, R12T, S12T, T12T, &
  !$  QH1 &
  !$)
  DO K = KS    , KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF((RHO(I,J,K) .GT. 0.0) .AND. (RHO(I,J,K+1) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J,K+1))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J,K+1))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J,K+1))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J,K+1))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J,K+1))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J,K+1))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J,K+1))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J,K+1))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J,K+1))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J,K+1))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J,K) + RHO(I,J,K+1))
    UM   = 0.5 * (  U(I,J,K) +   U(I,J,K+1))
    VM   = 0.5 * (  V(I,J,K) +   V(I,J,K+1))
    WM   = 0.5 * (  W(I,J,K) +   W(I,J,K+1))
    uuM  = 0.5 * ( uu(I,J,K) +  uu(I,J,K+1))
    vvM  = 0.5 * ( vv(I,J,K) +  vv(I,J,K+1))
    wwM  = 0.5 * ( ww(I,J,K) +  ww(I,J,K+1))
    uvM  = 0.5 * ( uv(I,J,K) +  uv(I,J,K+1))
    vwM  = 0.5 * ( vw(I,J,K) +  vw(I,J,K+1))
    wuM  = 0.5 * ( wu(I,J,K) +  wu(I,J,K+1))
    EPSM = 0.5 * (EPS(I,J,K) + EPS(I,J,K+1))
    AMUM = 0.5 * (AMU(I,J,K) + AMU(I,J,K+1))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J,K+1) &
    &              + U(I+1,J,K) + U(I+1,J,K+1) )
    UET = 0.25 * ( - U(I,J-1,K) - U(I,J-1,K+1) &
    &              + U(I,J+1,K) + U(I,J+1,K+1) )
    UZE = - U(I,J,K) + U(I,J,K+1)
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J,K+1) &
    &              + V(I+1,J,K) + V(I+1,J,K+1) )
    VET = 0.25 * ( - V(I,J-1,K) - V(I,J-1,K+1) &
    &              + V(I,J+1,K) + V(I,J+1,K+1) )
    VZE = - V(I,J,K) + V(I,J,K+1)
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J,K+1) &
    &              + W(I+1,J,K) + W(I+1,J,K+1) )
    WET = 0.25 * ( - W(I,J-1,K) - W(I,J-1,K+1) &
    &              + W(I,J+1,K) + W(I,J+1,K+1) )
    WZE = - W(I,J,K) + W(I,J,K+1)
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J,K+1) &
    &              + T(I+1,J,K) + T(I+1,J,K+1) )
    TET = 0.25 * ( - T(I,J-1,K) - T(I,J-1,K+1) &
    &              + T(I,J+1,K) + T(I,J+1,K+1) )
    TZE = - T(I,J,K) + T(I,J,K+1)
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = 0.25 * ( - uu(I-1,J,K) - uu(I-1,J,K+1) &
    &               + uu(I+1,J,K) + uu(I+1,J,K+1) )
    uuET = 0.25 * ( - uu(I,J-1,K) - uu(I,J-1,K+1) &
    &               + uu(I,J+1,K) + uu(I,J+1,K+1) )
    uuZE = - uu(I,J,K) + uu(I,J,K+1)
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM + uuZE * ZEXM
    uuY = uuXI * XIYM + uuET * ETYM + uuZE * ZEYM
    uuZ = uuXI * XIZM + uuET * ETZM + uuZE * ZEZM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = 0.25 * ( - vv(I-1,J,K) - vv(I-1,J,K+1) &
    &               + vv(I+1,J,K) + vv(I+1,J,K+1) )
    vvET = 0.25 * ( - vv(I,J-1,K) - vv(I,J-1,K+1) &
    &               + vv(I,J+1,K) + vv(I,J+1,K+1) )
    vvZE = - vv(I,J,K) + vv(I,J,K+1)
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM + vvZE * ZEXM
    vvY = vvXI * XIYM + vvET * ETYM + vvZE * ZEYM
    vvZ = vvXI * XIZM + vvET * ETZM + vvZE * ZEZM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = 0.25 * ( - ww(I-1,J,K) - ww(I-1,J,K+1) &
    &               + ww(I+1,J,K) + ww(I+1,J,K+1) )
    wwET = 0.25 * ( - ww(I,J-1,K) - ww(I,J-1,K+1) &
    &               + ww(I,J+1,K) + ww(I,J+1,K+1) )
    wwZE = - ww(I,J,K) + ww(I,J,K+1)
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM + wwZE * ZEXM
    wwY = wwXI * XIYM + wwET * ETYM + wwZE * ZEYM
    wwZ = wwXI * XIZM + wwET * ETZM + wwZE * ZEZM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = 0.25 * ( - uv(I-1,J,K) - uv(I-1,J,K+1) &
    &               + uv(I+1,J,K) + uv(I+1,J,K+1) )
    uvET = 0.25 * ( - uv(I,J-1,K) - uv(I,J-1,K+1) &
    &               + uv(I,J+1,K) + uv(I,J+1,K+1) )
    uvZE = - uv(I,J,K) + uv(I,J,K+1)
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM + uvZE * ZEXM
    uvY = uvXI * XIYM + uvET * ETYM + uvZE * ZEYM
    uvZ = uvXI * XIZM + uvET * ETZM + uvZE * ZEZM
    ! vw の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vwXI = 0.25 * ( - vw(I-1,J,K) - vw(I-1,J,K+1) &
    &               + vw(I+1,J,K) + vw(I+1,J,K+1) )
    vwET = 0.25 * ( - vw(I,J-1,K) - vw(I,J-1,K+1) &
    &               + vw(I,J+1,K) + vw(I,J+1,K+1) )
    vwZE = - vw(I,J,K) + vw(I,J,K+1)
    ! 物理空間方向一階微分
    vwX = vwXI * XIXM + vwET * ETXM + vwZE * ZEXM
    vwY = vwXI * XIYM + vwET * ETYM + vwZE * ZEYM
    vwZ = vwXI * XIZM + vwET * ETZM + vwZE * ZEZM
    ! wu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wuXI = 0.25 * ( - wu(I-1,J,K) - wu(I-1,J,K+1) &
    &               + wu(I+1,J,K) + wu(I+1,J,K+1) )
    wuET = 0.25 * ( - wu(I,J-1,K) - wu(I,J-1,K+1) &
    &               + wu(I,J+1,K) + wu(I,J+1,K+1) )
    wuZE = - wu(I,J,K) + wu(I,J,K+1)
    ! 物理空間方向一階微分
    wuX = wuXI * XIXM + wuET * ETXM + wuZE * ZEXM
    wuY = wuXI * XIYM + wuET * ETYM + wuZE * ZEYM
    wuZ = wuXI * XIZM + wuET * ETZM + wuZE * ZEZM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPS(I-1,J,K) - EPS(I-1,J,K+1) &
    &               + EPS(I+1,J,K) + EPS(I+1,J,K+1) )
    AEET = 0.25 * ( - EPS(I,J-1,K) - EPS(I,J-1,K+1) &
    &               + EPS(I,J+1,K) + EPS(I,J+1,K+1) )
    AEZE = - EPS(I,J,K) + EPS(I,J,K+1)
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY + WZ) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 = 2.0 * WZ - SS
    S12 = UY + VX
    S23 = VZ + WY
    S31 = WX + UZ
    ! U の拡散 ---------------------------------------------------------
    R02M = nu * S11
    S02M = nu * S12
    T02M = nu * S31
    R02T =-uuM
    S02T =-uvM
    T02T =-wuM
    ! V の拡散 ---------------------------------------------------------
    R03M = nu * S12
    S03M = nu * S22
    T03M = nu * S23
    R03T =-uvM
    S03T =-vvM
    T03T =-vwM
    ! W の拡散 ---------------------------------------------------------
    R04M = nu * S31
    S04M = nu * S23
    T04M = nu * S33
    R04T =-wuM
    S04T =-vwM
    T04T =-wwM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R05M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM + S31 * WM)
    S05M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM + S23 * WM)
    T05M = C_p * nu / PR * TZ + nu * (S31 * UM + S23 * VM + S33 * WM)
    R05T = AKAE * C_se * C_p * (uuM * TX + uvM * TY + wuM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uuM * (uuX + vvX + wwX) &
    &    + uvM * (uuY + vvY + wwY) &
    &    + wuM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uuM * UM + uvM * VM + wuM * WM)
    S05T = AKAE * C_se * C_p * (uvM * TX + vvM * TY + vwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uvM * (uuX + vvX + wwX) &
    &    + vvM * (uuY + vvY + wwY) &
    &    + vwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uvM * UM + vvM * VM + vwM * WM)
    T05T = AKAE * C_se * C_p * (wuM * TX + vwM * TY + wwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      wuM * (uuX + vvX + wwX) &
    &    + vwM * (uuY + vvY + wwY) &
    &    + wwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (wuM * UM + vwM * VM + wwM * WM)
    ! uu の拡散 --------------------------------------------------------
    R06M = nu * uuX
    S06M = nu * uuY
    T06M = nu * uuZ
    R06T = C_s * AKAE * (uuM * uuX + uvM * uuY + wuM * uuZ)
    S06T = C_s * AKAE * (uvM * uuX + vvM * uuY + vwM * uuZ)
    T06T = C_s * AKAE * (wuM * uuX + vwM * uuY + wwM * uuZ)
    ! vv の拡散 --------------------------------------------------------
    R07M = nu * vvX
    S07M = nu * vvY
    T07M = nu * vvZ
    R07T = C_s * AKAE * (uuM * vvX + uvM * vvY + wuM * vvZ)
    S07T = C_s * AKAE * (uvM * vvX + vvM * vvY + vwM * vvZ)
    T07T = C_s * AKAE * (wuM * vvX + vwM * vvY + wwM * vvZ)
    ! ww の拡散 --------------------------------------------------------
    R08M = nu * wwX
    S08M = nu * wwY
    T08M = nu * wwZ
    R08T = C_s * AKAE * (uuM * wwX + uvM * wwY + wuM * wwZ)
    S08T = C_s * AKAE * (uvM * wwX + vvM * wwY + vwM * wwZ)
    T08T = C_s * AKAE * (wuM * wwX + vwM * wwY + wwM * wwZ)
    ! uv の拡散 --------------------------------------------------------
    R09M = nu * uvX
    S09M = nu * uvY
    T09M = nu * uvZ
    R09T = C_s * AKAE * (uuM * uvX + uvM * uvY + wuM * uvZ)
    S09T = C_s * AKAE * (uvM * uvX + vvM * uvY + vwM * uvZ)
    T09T = C_s * AKAE * (wuM * uvX + vwM * uvY + wwM * uvZ)
    ! vw の拡散 --------------------------------------------------------
    R10M = nu * vwX
    S10M = nu * vwY
    T10M = nu * vwZ
    R10T = C_s * AKAE * (uuM * vwX + uvM * vwY + wuM * vwZ)
    S10T = C_s * AKAE * (uvM * vwX + vvM * vwY + vwM * vwZ)
    T10T = C_s * AKAE * (wuM * vwX + vwM * vwY + wwM * vwZ)
    ! wu の拡散 --------------------------------------------------------
    R11M = nu * wuX
    S11M = nu * wuY
    T11M = nu * wuZ
    R11T = C_s * AKAE * (uuM * wuX + uvM * wuY + wuM * wuZ)
    S11T = C_s * AKAE * (uvM * wuX + vvM * wuY + vwM * wuZ)
    T11T = C_s * AKAE * (wuM * wuX + vwM * wuY + wwM * wuZ)
    ! epsilon の拡散 ---------------------------------------------------
    R12M = nu * AEX
    S12M = nu * AEY
    T12M = nu * AEZ
    R12T = C_e * AKAE * (uuM * AEX + uvM * AEY + wuM * AEZ)
    S12T = C_e * AKAE * (uvM * AEX + vvM * AEY + vwM * AEZ)
    T12T = C_e * AKAE * (wuM * AEX + vwM * AEY + wwM * AEZ)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    TH(I,J,K, 1) = 0.0
    TH(I,J,K, 2) = QH1 * ( (R02M + R02T) * ZEXM &
    &                    + (S02M + S02T) * ZEYM &
    &                    + (T02M + T02T) * ZEZM )
    TH(I,J,K, 3) = QH1 * ( (R03M + R03T) * ZEXM &
    &                    + (S03M + S03T) * ZEYM &
    &                    + (T03M + T03T) * ZEZM )
    TH(I,J,K, 4) = QH1 * ( (R04M + R04T) * ZEXM &
    &                    + (S04M + S04T) * ZEYM &
    &                    + (T04M + T04T) * ZEZM )
    TH(I,J,K, 5) = QH1 * ( (R05M + R05T) * ZEXM &
    &                    + (S05M + S05T) * ZEYM &
    &                    + (T05M + T05T) * ZEZM )
    TH(I,J,K, 6) = QH1 * ( (R06M + R06T) * ZEXM &
    &                    + (S06M + S06T) * ZEYM &
    &                    + (T06M + T06T) * ZEZM )
    TH(I,J,K, 7) = QH1 * ( (R07M + R07T) * ZEXM &
    &                    + (S07M + S07T) * ZEYM &
    &                    + (T07M + T07T) * ZEZM )
    TH(I,J,K, 8) = QH1 * ( (R08M + R08T) * ZEXM &
    &                    + (S08M + S08T) * ZEYM &
    &                    + (T08M + T08T) * ZEZM )
    TH(I,J,K, 9) = QH1 * ( (R09M + R09T) * ZEXM &
    &                    + (S09M + S09T) * ZEYM &
    &                    + (T09M + T09T) * ZEZM )
    TH(I,J,K,10) = QH1 * ( (R10M + R10T) * ZEXM &
    &                    + (S10M + S10T) * ZEYM &
    &                    + (T10M + T10T) * ZEZM )
    TH(I,J,K,11) = QH1 * ( (R11M + R11T) * ZEXM &
    &                    + (S11M + S11T) * ZEYM &
    &                    + (T11M + T11T) * ZEZM )
    TH(I,J,K,12) = QH1 * ( (R12M + R12T) * ZEXM &
    &                    + (S12M + S12T) * ZEYM &
    &                    + (T12M + T12T) * ZEZM )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFZE
! 定義終了 *************************************************************
END SUBROUTINE Turbulence3DRSMGL
!***********************************************************************
!**** 乱流モデル : RANS, RSM Speziale-Sarkar-Gatski Model (1991)    ****
!****              epsilon方程式の定数をGibson-Launderの値に変更    ****
!****              乱流熱伝達にはVandrommeのモデルを使用            ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence2DRSMSSG( &
&            PELIM, RG, GAMMA, PR, &
&            IS, IE, JS, JE, LS, LE, &
&            XIX, XIY, ETX, ETY, AJA, &
&            RHO, U, V, T, uu, vv, ww, uv, EPS, AMU, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: C_se = 0.313
  ! REAL, PARAMETER :: C_s = 0.22, C_e = 0.183
  REAL, PARAMETER :: C_s = 0.22, C_e = 0.18
  REAL, PARAMETER :: C_1 = 3.4, C_2 = 4.2, C_3 = 0.8, &
  &                  C_4 = 1.25, C_5 = 0.4
  REAL, PARAMETER :: C_1a = 1.8, C_3a = 1.3
  ! REAL, PARAMETER :: C_e1 = 1.44, C_e2 = 1.83
  REAL, PARAMETER :: C_e1 = 1.44, C_e2 = 1.92
  REAL, PARAMETER :: OneThird = 0.33333333, TwoThird = 0.66666667
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE), XIY(IS:IE, JS:JE), &
  &                       ETX(IS:IE, JS:JE), ETY(IS:IE, JS:JE), &
  &                       AJA(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE), &
  &                       U(IS:IE, JS:JE), &
  &                       V(IS:IE, JS:JE), &
  &                       T(IS:IE, JS:JE), &
  &                       uu(IS:IE, JS:JE), &
  &                       vv(IS:IE, JS:JE), &
  &                       ww(IS:IE, JS:JE), &
  &                       uv(IS:IE, JS:JE), &
  &                       EPS(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: AMUT(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(OUT) :: DQP(IS:IE, JS:JE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :), SH(:, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 9) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, LS: LE), SH(IS: IE, JS: JE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: TMP, TMP1, TMP2
  REAL    :: UXI, VXI, UET, VET, UX, VX, UY, VY
  REAL    :: P11, P22, P33, P12, P_k
  REAL    :: eps11, eps22, eps33, eps12, e_k
  REAL    :: S11, S22, S33, S12, SS2
  REAL    :: O12
  REAL    :: b11, b22, b33, b12
  REAL    :: phi11, phi22, phi33, phi12
  REAL    :: AK, nu, AEAK
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, TMP, TMP1, TMP2, &
  !$  UXI, VXI, UET, VET, UX, VX, UY, VY, &
  !$  P11, P22, P33, P12, P_k, eps11, eps22, eps33, eps12, e_k, &
  !$  S11, S22, S33, S12, SS2, O12, b11, b22, b33, b12, &
  !$  phi11, phi22, phi33, phi12, AK, nu, AEAK, QH1 &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF(RHO(I,J) .GT. 0.0) THEN
    ! 乱流エネルギー ---------------------------------------------------
    AK = 0.5 * (uu(I,J) + vv(I,J) + ww(I,J))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.5 * (U(I+1,J) - U(I-1,J))
    VXI = 0.5 * (V(I+1,J) - V(I-1,J))
    UET = 0.5 * (U(I,J+1) - U(I,J-1))
    VET = 0.5 * (V(I,J+1) - V(I,J-1))
    ! 物理空間方向一階微分
    UX = UXI * XIX(I,J) + UET * ETX(I,J)
    UY = UXI * XIY(I,J) + UET * ETY(I,J)
    VX = VXI * XIX(I,J) + VET * ETX(I,J)
    VY = VXI * XIY(I,J) + VET * ETY(I,J)
    ! 生成項と散逸項 ---------------------------------------------------
    P11 = - 2.0 * (uu(I,J) * UX + uv(I,J) * UY)
    P22 = - 2.0 * (uv(I,J) * VX + vv(I,J) * VY)
    P33 = 0.0
    P12 = - (uu(I,J) * VX + uv(I,J) * (UX + VY) + vv(I,J) * UY)
    P_k = 0.5 * (P11 + P22 + P33)
    e_k = EPS(I,J)
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      TMP = TwoThird * e_k / PELIM
      P11 = MIN(TMP, MAX(-TMP, P11))
      P22 = MIN(TMP, MAX(-TMP, P22))
      P33 = MIN(TMP, MAX(-TMP, P33))
      P_k = 0.5 * (P11 + P22 + P33)
    ENDIF
    TMP   = TwoThird * e_k
    eps11 = TMP
    eps22 = TMP
    eps33 = TMP
    eps12 = 0.0
    ! 歪み速度 ---------------------------------------------------------
    TMP = (UX + VY) * OneThird
    S11 = UX - TMP
    S22 = VY - TMP
    S33 =    - TMP
    S12 = 0.5 * (UY + VX)
    SS2 = S11 * S11 + S22 * S22 + S33 * S33 + 2.0 * S12 * S12
    ! 渦度 -------------------------------------------------------------
    O12 = 0.5 * (UY - VX)
    ! 非等方レイノルズ応力テンソル b_ij --------------------------------
    IF(AK .GT. ZERO) THEN
      b11 = 0.5 * uu(I,J) / AK - OneThird
      b22 = 0.5 * vv(I,J) / AK - OneThird
      b33 = 0.5 * ww(I,J) / AK - OneThird
      b12 = 0.5 * uv(I,J) / AK
      b11 = MIN(TwoThird, MAX(- OneThird, b11))
      b22 = MIN(TwoThird, MAX(- OneThird, b22))
      b33 = MIN(TwoThird, MAX(- OneThird, b33))
      b12 = MIN(0.5, MAX(-0.5, b12))
    ELSE
      b11 =-OneThird
      b22 =-OneThird
      b33 =-OneThird
      b12 = 0.0
    ENDIF
    ! 再分配項 ---------------------------------------------------------
    TMP1 = b11 * b11 + b22 * b22 + b33 * b33 + 2.0 * b12 * b12
    TMP2 = b11 * S11 + b22 * S22 + b33 * S33 + 2.0 * b12 * S12
    TMP2 = TMP2 * TwoThird
    phi11 =-(C_1 * e_k + C_1a * P_k) * b11 &
    &     + C_2 * e_k * (b11**2 + b12**2 - TMP1 * OneThird) &
    &     + (C_3 - C_3a * SQRT(TMP1)) * AK * S11 &
    &     + C_4 * AK * (2.0 * (b11 * S11 + b12 * S12) - TMP2) &
    &     + C_5 * AK * (2.0 * b12 * O12)
    phi22 =-(C_1 * e_k + C_1a * P_k) * b22 &
    &     + C_2 * e_k * (b12**2 + b22**2 - TMP1 * OneThird) &
    &     + (C_3 - C_3a * SQRT(TMP1)) * AK * S22 &
    &     + C_4 * AK * (2.0 * (b12 * S12 + b22 * S22) - TMP2) &
    &     - C_5 * AK * (2.0 * b12 * O12)
    phi33 =-(C_1 * e_k + C_1a * P_k) * b33 &
    &     + C_2 * e_k * (b33 * b33 - TMP1 * OneThird) &
    &     + (C_3 - C_3a * SQRT(TMP1)) * AK * S33 &
    &     + C_4 * AK * (2.0 * b33 * S33 - TMP2)
    phi12 =-(C_1 * e_k + C_1a * P_k) * b12 &
    &     + C_2 * e_k * (b11 + b22) * b12 &
    &     + (C_3 - C_3a * SQRT(TMP1)) * AK * S12 &
    &     + C_4 * AK * ((b11 + b22) * S12 + b12 * (S11 + S22)) &
    &     + C_5 * AK * (b22 - b11) * O12
    ! 渦粘性係数の計算 -------------------------------------------------
    IF(e_k .GT. ZERO) THEN
      AMUT(I,J) = 0.09 * AK**2 / e_k * RHO(I,J)
    ENDIF
    ! 生成項、圧力再分配項、散逸項の和 ---------------------------------
    QH1  = RHO(I,J) / AJA(I,J)
    nu   = AMU(I,J) / RHO(I,J)
    AEAK = MIN( &
    &        MAX(SQRT(2.0 * SS2), SQRT(e_k / nu)), &
    &        e_k / MAX(ZERO, AK) &
    &    )
    DQP(I,J,5) = QH1 * (P11 - eps11 + phi11)
    DQP(I,J,6) = QH1 * (P22 - eps22 + phi22)
    DQP(I,J,7) = QH1 * (P33 - eps33 + phi33)
    DQP(I,J,8) = QH1 * (P12 - eps12 + phi12)
    DQP(I,J,9) = QH1 * AEAK * (C_e1 * P_k - C_e2 * e_k)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH = 0.0
  SH = 0.0
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL DIFFXI
  CALL DIFFET
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, L)
  DO L = LS, LE
  !$OMP DO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J-1) .GT. 0.0) .AND. (RHO(I-1,J) .GT. 0.0) .AND. &
  &   (RHO(I,J)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,L) = (-RH(I-1,J  ,L) + RH(I,J,L)) &
    &          + (-SH(I  ,J-1,L) + SH(I,J,L))
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項の計算(xi方向差分用)                                    ****
!***********************************************************************
SUBROUTINE DIFFXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: RHOM, UM, VM, uuM, vvM, wwM, uvM, EPSM, AMUM
  REAL    :: UXI, UET, UX, UY
  REAL    :: VXI, VET, VX, VY
  REAL    :: TXI, TET, TX, TY
  REAL    :: uuXI, uuET, uuX, uuY
  REAL    :: vvXI, vvET, vvX, vvY
  REAL    :: wwXI, wwET, wwX, wwY
  REAL    :: uvXI, uvET, uvX, uvY
  REAL    :: AEXI, AEET, AEX, AEY
  REAL    :: nu
  REAL    :: S11, S22, S12, SS
  REAL    :: C_p, AKAE
  REAL    :: R2M, S2M, R2T, S2T
  REAL    :: R3M, S3M, R3T, S3T
  REAL    :: R4M, S4M, R4T, S4T
  REAL    :: R5M, S5M, R5T, S5T
  REAL    :: R6M, S6M, R6T, S6T
  REAL    :: R7M, S7M, R7T, S7T
  REAL    :: R8M, S8M, R8T, S8T
  REAL    :: R9M, S9M, R9T, S9T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, &
  !$  RHOM, UM, VM, uuM, vvM, wwM, uvM, EPSM, AMUM, &
  !$  UXI, UET, UX, UY, VXI, VET, VX, VY, TXI, TET, TX, TY, &
  !$  uuXI, uuET, uuX, uuY, vvXI, vvET, vvX, vvY, &
  !$  wwXI, wwET, wwX, wwY, uvXI, uvET, uvX, uvY, &
  !$  AEXI, AEET, AEX, AEY, &
  !$  nu, S11, S22, S12, SS, C_p, AKAE, &
  !$  R2M, S2M, R2T, S2T, R3M, S3M, R3T, S3T, R4M, S4M, R4T, S4T, &
  !$  R5M, S5M, R5T, S5T, R6M, S6M, R6T, S6T, R7M, S7M, R7T, S7T, &
  !$  R8M, S8M, R8T, S8T, R9M, S9M, R9T, S9T, QH1 &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF((RHO(I,J) .GT. 0.0) .AND. (RHO(I+1,J) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J) + XIX(I+1,J))
    XIYM = 0.5 * (XIY(I,J) + XIY(I+1,J))
    ETXM = 0.5 * (ETX(I,J) + ETX(I+1,J))
    ETYM = 0.5 * (ETY(I,J) + ETY(I+1,J))
    AJAM = 0.5 * (AJA(I,J) + AJA(I+1,J))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J) + RHO(I+1,J))
    UM   = 0.5 * (  U(I,J) +   U(I+1,J))
    VM   = 0.5 * (  V(I,J) +   V(I+1,J))
    uuM  = 0.5 * ( uu(I,J) +  uu(I+1,J))
    vvM  = 0.5 * ( vv(I,J) +  vv(I+1,J))
    wwM  = 0.5 * ( ww(I,J) +  ww(I+1,J))
    uvM  = 0.5 * ( uv(I,J) +  uv(I+1,J))
    EPSM = 0.5 * (EPS(I,J) + EPS(I+1,J))
    AMUM = 0.5 * (AMU(I,J) + AMU(I+1,J))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = - U(I,J) + U(I+1,J)
    UET = 0.25 * (- U(I,J-1) - U(I+1,J-1) + U(I,J+1) + U(I+1,J+1))
    VXI = - V(I,J) + V(I+1,J)
    VET = 0.25 * (- V(I,J-1) - V(I+1,J-1) + V(I,J+1) + V(I+1,J+1))
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM
    UY = UXI * XIYM + UET * ETYM
    VX = VXI * XIXM + VET * ETXM
    VY = VXI * XIYM + VET * ETYM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = - T(I,J) + T(I+1,J)
    TET = 0.25 * (- T(I,J-1) - T(I+1,J-1) + T(I,J+1) + T(I+1,J+1))
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM
    TY = TXI * XIYM + TET * ETYM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = - uu(I,J) + uu(I+1,J)
    uuET = 0.25 * ( - uu(I,J-1) - uu(I+1,J-1) &
    &               + uu(I,J+1) + uu(I+1,J+1) )
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM
    uuY = uuXI * XIYM + uuET * ETYM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = - vv(I,J) + vv(I+1,J)
    vvET = 0.25 * ( - vv(I,J-1) - vv(I+1,J-1) &
    &               + vv(I,J+1) + vv(I+1,J+1) )
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM
    vvY = vvXI * XIYM + vvET * ETYM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = - ww(I,J) + ww(I+1,J)
    wwET = 0.25 * ( - ww(I,J-1) - ww(I+1,J-1) &
    &               + ww(I,J+1) + ww(I+1,J+1) )
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM
    wwY = wwXI * XIYM + wwET * ETYM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = - uv(I,J) + uv(I+1,J)
    uvET = 0.25 * ( - uv(I,J-1) - uv(I+1,J-1) &
    &               + uv(I,J+1) + uv(I+1,J+1) )
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM
    uvY = uvXI * XIYM + uvET * ETYM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = - EPS(I,J) + EPS(I+1,J)
    AEET = 0.25 * ( - EPS(I,J-1) - EPS(I+1,J-1) &
    &               + EPS(I,J+1) + EPS(I+1,J+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM
    AEY = AEXI * XIYM + AEET * ETYM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S12 = UY + VX
    ! U の拡散 ---------------------------------------------------------
    R2M = nu * S11
    S2M = nu * S12
    R2T =-uuM
    S2T =-uvM
    ! V の拡散 ---------------------------------------------------------
    R3M = nu * S12
    S3M = nu * S22
    R3T =-uvM
    S3T =-vvM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R4M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM)
    S4M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM)
    R4T = AKAE * C_se * C_p * (uuM * TX + uvM * TY) &
    &   + AKAE * C_s  * ( &
    &     uuM * (uuX + vvX + wwX) + uvM * (uuY + vvY + wwY) &
    &   ) &
    &   - (uuM * UM + uvM * VM)
    S4T = AKAE * C_se * C_p * (uvM * TX + vvM * TY) &
    &   + AKAE * C_s  * ( &
    &     uvM * (uuX + vvX + wwX) + vvM * (uuY + vvY + wwY) &
    &   ) &
    &   - (uvM * UM + vvM * VM)
    ! uu の拡散 --------------------------------------------------------
    R5M = nu * uuX
    S5M = nu * uuY
    R5T = C_s * AKAE * (uuM * uuX + uvM * uuY)
    S5T = C_s * AKAE * (uvM * uuX + vvM * uuY)
    ! vv の拡散 --------------------------------------------------------
    R6M = nu * vvX
    S6M = nu * vvY
    R6T = C_s * AKAE * (uuM * vvX + uvM * vvY)
    S6T = C_s * AKAE * (uvM * vvX + vvM * vvY)
    ! ww の拡散 --------------------------------------------------------
    R7M = nu * wwX
    S7M = nu * wwY
    R7T = C_s * AKAE * (uuM * wwX + uvM * wwY)
    S7T = C_s * AKAE * (uvM * wwX + vvM * wwY)
    ! uv の拡散 --------------------------------------------------------
    R8M = nu * uvX
    S8M = nu * uvY
    R8T = C_s * AKAE * (uuM * uvX + uvM * uvY)
    S8T = C_s * AKAE * (uvM * uvX + vvM * uvY)
    ! epsilon の拡散 ---------------------------------------------------
    R9M = nu * AEX
    S9M = nu * AEY
    R9T = C_e * AKAE * (uuM * AEX + uvM * AEY)
    S9T = C_e * AKAE * (uvM * AEX + vvM * AEY)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    RH(I,J,1) = 0.0
    RH(I,J,2) = QH1 * ((R2M + R2T) * XIXM + (S2M + S2T) * XIYM)
    RH(I,J,3) = QH1 * ((R3M + R3T) * XIXM + (S3M + S3T) * XIYM)
    RH(I,J,4) = QH1 * ((R4M + R4T) * XIXM + (S4M + S4T) * XIYM)
    RH(I,J,5) = QH1 * ((R5M + R5T) * XIXM + (S5M + S5T) * XIYM)
    RH(I,J,6) = QH1 * ((R6M + R6T) * XIXM + (S6M + S6T) * XIYM)
    RH(I,J,7) = QH1 * ((R7M + R7T) * XIXM + (S7M + S7T) * XIYM)
    RH(I,J,8) = QH1 * ((R8M + R8T) * XIXM + (S8M + S8T) * XIYM)
    RH(I,J,9) = QH1 * ((R9M + R9T) * XIXM + (S9M + S9T) * XIYM)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFXI
!***********************************************************************
!**** 拡散項の計算(eta方向差分用)                                   ****
!***********************************************************************
SUBROUTINE DIFFET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: RHOM, UM, VM, uuM, vvM, wwM, uvM, EPSM, AMUM
  REAL    :: UXI, UET, UX, UY
  REAL    :: VXI, VET, VX, VY
  REAL    :: TXI, TET, TX, TY
  REAL    :: uuXI, uuET, uuX, uuY
  REAL    :: vvXI, vvET, vvX, vvY
  REAL    :: wwXI, wwET, wwX, wwY
  REAL    :: uvXI, uvET, uvX, uvY
  REAL    :: AEXI, AEET, AEX, AEY
  REAL    :: nu
  REAL    :: S11, S22, S12, SS
  REAL    :: C_p, AKAE
  REAL    :: R2M, S2M, R2T, S2T
  REAL    :: R3M, S3M, R3T, S3T
  REAL    :: R4M, S4M, R4T, S4T
  REAL    :: R5M, S5M, R5T, S5T
  REAL    :: R6M, S6M, R6T, S6T
  REAL    :: R7M, S7M, R7T, S7T
  REAL    :: R8M, S8M, R8T, S8T
  REAL    :: R9M, S9M, R9T, S9T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, &
  !$  RHOM, UM, VM, uuM, vvM, wwM, uvM, EPSM, AMUM, &
  !$  UXI, UET, UX, UY, VXI, VET, VX, VY, TXI, TET, TX, TY, &
  !$  uuXI, uuET, uuX, uuY, vvXI, vvET, vvX, vvY, &
  !$  wwXI, wwET, wwX, wwY, uvXI, uvET, uvX, uvY, &
  !$  AEXI, AEET, AEX, AEY, &
  !$  nu, S11, S22, S12, SS, C_p, AKAE, &
  !$  R2M, S2M, R2T, S2T, R3M, S3M, R3T, S3T, R4M, S4M, R4T, S4T, &
  !$  R5M, S5M, R5T, S5T, R6M, S6M, R6T, S6T, R7M, S7M, R7T, S7T, &
  !$  R8M, S8M, R8T, S8T, R9M, S9M, R9T, S9T, QH1 &
  !$)
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF((RHO(I,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0)) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J) + XIX(I,J+1))
    XIYM = 0.5 * (XIY(I,J) + XIY(I,J+1))
    ETXM = 0.5 * (ETX(I,J) + ETX(I,J+1))
    ETYM = 0.5 * (ETY(I,J) + ETY(I,J+1))
    AJAM = 0.5 * (AJA(I,J) + AJA(I,J+1))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J) + RHO(I,J+1))
    UM   = 0.5 * (  U(I,J) +   U(I,J+1))
    VM   = 0.5 * (  V(I,J) +   V(I,J+1))
    uuM  = 0.5 * ( uu(I,J) +  uu(I,J+1))
    vvM  = 0.5 * ( vv(I,J) +  vv(I,J+1))
    wwM  = 0.5 * ( ww(I,J) +  ww(I,J+1))
    uvM  = 0.5 * ( uv(I,J) +  uv(I,J+1))
    EPSM = 0.5 * (EPS(I,J) + EPS(I,J+1))
    AMUM = 0.5 * (AMU(I,J) + AMU(I,J+1))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * (- U(I-1,J) - U(I-1,J+1) + U(I+1,J) + U(I+1,J+1))
    UET = - U(I,J) + U(I,J+1)
    VXI = 0.25 * (- V(I-1,J) - V(I-1,J+1) + V(I+1,J) + V(I+1,J+1))
    VET = - V(I,J) + V(I,J+1)
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM
    UY = UXI * XIYM + UET * ETYM
    VX = VXI * XIXM + VET * ETXM
    VY = VXI * XIYM + VET * ETYM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * (- T(I-1,J) - T(I-1,J+1) + T(I+1,J) + T(I+1,J+1))
    TET = - T(I,J) + T(I,J+1)
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM
    TY = TXI * XIYM + TET * ETYM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = 0.25 * (- uu(I-1,J) - uu(I-1,J+1) + uu(I+1,J) + uu(I+1,J+1))
    uuET = - uu(I,J) + uu(I,J+1)
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM
    uuY = uuXI * XIYM + uuET * ETYM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = 0.25 * (- vv(I-1,J) - vv(I-1,J+1) + vv(I+1,J) + vv(I+1,J+1))
    vvET = - vv(I,J) + vv(I,J+1)
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM
    vvY = vvXI * XIYM + vvET * ETYM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = 0.25 * (- ww(I-1,J) - ww(I-1,J+1) + ww(I+1,J) + ww(I+1,J+1))
    wwET = - ww(I,J) + ww(I,J+1)
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM
    wwY = wwXI * XIYM + wwET * ETYM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = 0.25 * (- uv(I-1,J) - uv(I-1,J+1) + uv(I+1,J) + uv(I+1,J+1))
    uvET = - uv(I,J) + uv(I,J+1)
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM
    uvY = uvXI * XIYM + uvET * ETYM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPS(I-1,J) - EPS(I-1,J+1) &
    &               + EPS(I+1,J) + EPS(I+1,J+1) )
    AEET = - EPS(I,J) + EPS(I,J+1)
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM
    AEY = AEXI * XIYM + AEET * ETYM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S12 = UY + VX
    ! U の拡散 ---------------------------------------------------------
    R2M = nu * S11
    S2M = nu * S12
    R2T =-uuM
    S2T =-uvM
    ! V の拡散 ---------------------------------------------------------
    R3M = nu * S12
    S3M = nu * S22
    R3T =-uvM
    S3T =-vvM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R4M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM)
    S4M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM)
    R4T = AKAE * C_se * C_p * (uuM * TX + uvM * TY) &
    &   + AKAE * C_s  * ( &
    &     uuM * (uuX + vvX + wwX) + uvM * (uuY + vvY + wwY) &
    &   ) &
    &   - (uuM * UM + uvM * VM)
    S4T = AKAE * C_se * C_p * (uvM * TX + vvM * TY) &
    &   + AKAE * C_s  * ( &
    &     uvM * (uuX + vvX + wwX) + vvM * (uuY + vvY + wwY) &
    &   ) &
    &   - (uvM * UM + vvM * VM)
    ! uu の拡散 --------------------------------------------------------
    R5M = nu * uuX
    S5M = nu * uuY
    R5T = C_s * AKAE * (uuM * uuX + uvM * uuY)
    S5T = C_s * AKAE * (uvM * uuX + vvM * uuY)
    ! vv の拡散 --------------------------------------------------------
    R6M = nu * vvX
    S6M = nu * vvY
    R6T = C_s * AKAE * (uuM * vvX + uvM * vvY)
    S6T = C_s * AKAE * (uvM * vvX + vvM * vvY)
    ! ww の拡散 --------------------------------------------------------
    R7M = nu * wwX
    S7M = nu * wwY
    R7T = C_s * AKAE * (uuM * wwX + uvM * wwY)
    S7T = C_s * AKAE * (uvM * wwX + vvM * wwY)
    ! uv の拡散 --------------------------------------------------------
    R8M = nu * uvX
    S8M = nu * uvY
    R8T = C_s * AKAE * (uuM * uvX + uvM * uvY)
    S8T = C_s * AKAE * (uvM * uvX + vvM * uvY)
    ! epsilon の拡散 ---------------------------------------------------
    R9M = nu * AEX
    S9M = nu * AEY
    R9T = C_e * AKAE * (uuM * AEX + uvM * AEY)
    S9T = C_e * AKAE * (uvM * AEX + vvM * AEY)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    SH(I,J,1) = 0.0
    SH(I,J,2) = QH1 * ((R2M + R2T) * ETXM + (S2M + S2T) * ETYM)
    SH(I,J,3) = QH1 * ((R3M + R3T) * ETXM + (S3M + S3T) * ETYM)
    SH(I,J,4) = QH1 * ((R4M + R4T) * ETXM + (S4M + S4T) * ETYM)
    SH(I,J,5) = QH1 * ((R5M + R5T) * ETXM + (S5M + S5T) * ETYM)
    SH(I,J,6) = QH1 * ((R6M + R6T) * ETXM + (S6M + S6T) * ETYM)
    SH(I,J,7) = QH1 * ((R7M + R7T) * ETXM + (S7M + S7T) * ETYM)
    SH(I,J,8) = QH1 * ((R8M + R8T) * ETXM + (S8M + S8T) * ETYM)
    SH(I,J,9) = QH1 * ((R9M + R9T) * ETXM + (S9M + S9T) * ETYM)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFET
! 定義終了 *************************************************************
END SUBROUTINE Turbulence2DRSMSSG
!***********************************************************************
!**** 乱流モデル : RANS, RSM Speziale-Sarkar-Gatski Model (1991)    ****
!****              epsilon方程式の定数をGibson-Launderの値に変更    ****
!****              乱流熱伝達にはVandrommeのモデルを使用            ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence3DRSMSSG( &
&            PELIM, RG, GAMMA, PR, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, AJA, &
&            RHO, U, V, W, T, uu, vv, ww, uv, vw, wu, EPS, AMU, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: C_se = 0.313
  ! REAL, PARAMETER :: C_s = 0.22, C_e = 0.183
  REAL, PARAMETER :: C_s = 0.22, C_e = 0.18
  REAL, PARAMETER :: C_1 = 3.4, C_2 = 4.2, C_3 = 0.8, &
  &                  C_4 = 1.25, C_5 = 0.4
  REAL, PARAMETER :: C_1a = 1.8, C_3a = 1.3
  ! REAL, PARAMETER :: C_e1 = 1.44, C_e2 = 1.83
  REAL, PARAMETER :: C_e1 = 1.44, C_e2 = 1.92
  REAL, PARAMETER :: OneThird = 0.33333333, TwoThird = 0.66666667
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE, KS:KE), &
  &                       XIY(IS:IE, JS:JE, KS:KE), &
  &                       XIZ(IS:IE, JS:JE, KS:KE), &
  &                       ETX(IS:IE, JS:JE, KS:KE), &
  &                       ETY(IS:IE, JS:JE, KS:KE), &
  &                       ETZ(IS:IE, JS:JE, KS:KE), &
  &                       ZEX(IS:IE, JS:JE, KS:KE), &
  &                       ZEY(IS:IE, JS:JE, KS:KE), &
  &                       ZEZ(IS:IE, JS:JE, KS:KE), &
  &                       AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE, KS:KE), &
  &                       U(IS:IE, JS:JE, KS:KE), &
  &                       V(IS:IE, JS:JE, KS:KE), &
  &                       W(IS:IE, JS:JE, KS:KE), &
  &                       T(IS:IE, JS:JE, KS:KE), &
  &                       uu(IS:IE, JS:JE, KS:KE), &
  &                       vv(IS:IE, JS:JE, KS:KE), &
  &                       ww(IS:IE, JS:JE, KS:KE), &
  &                       uv(IS:IE, JS:JE, KS:KE), &
  &                       vw(IS:IE, JS:JE, KS:KE), &
  &                       wu(IS:IE, JS:JE, KS:KE), &
  &                       EPS(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: AMUT(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(OUT) :: DQP(IS:IE, JS:JE, KS:KE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :, :), SH(:, :, :, :), TH(:, :, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 12) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(SH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(TH(IS: IE, JS: JE, KS: KE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: TMP, TMPBB, TMPBS, TMP1, TMP3
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  REAL    :: UX, VX, WX, UY, VY, WY, UZ, VZ, WZ
  REAL    :: P11, P22, P33, P12, P23, P31, P_k
  REAL    :: eps11, eps22, eps33, eps12, eps23, eps31, e_k
  REAL    :: S11, S22, S33, S12, S23, S31, SS2
  REAL    :: O12, O23, O31
  REAL    :: b11, b22, b33, b12, b23, b31
  REAL    :: phi11, phi22, phi33, phi12, phi23, phi31
  REAL    :: AK, nu, AEAK
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  TMP, TMPBB, TMPBS, TMP1, TMP3, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE, &
  !$  UX, VX, WX, UY, VY, WY, UZ, VZ, WZ, &
  !$  P11, P22, P33, P12, P23, P31, P_k, &
  !$  eps11, eps22, eps33, eps12, eps23, eps31, e_k, &
  !$  S11, S22, S33, S12, S23, S31, SS2, O12, O23, O31, &
  !$  b11, b22, b33, b12, b23, b31, &
  !$  phi11, phi22, phi33, phi12, phi23, phi31, AK, nu, AEAK, QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF(RHO(I,J,K) .GT. 0.0) THEN
    ! 乱流エネルギー ---------------------------------------------------
    AK = 0.5 * (uu(I,J,K) + vv(I,J,K) + ww(I,J,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.5 * (U(I+1,J,K) - U(I-1,J,K))
    UET = 0.5 * (U(I,J+1,K) - U(I,J-1,K))
    UZE = 0.5 * (U(I,J,K+1) - U(I,J,K-1))
    VXI = 0.5 * (V(I+1,J,K) - V(I-1,J,K))
    VET = 0.5 * (V(I,J+1,K) - V(I,J-1,K))
    VZE = 0.5 * (V(I,J,K+1) - V(I,J,K-1))
    WXI = 0.5 * (W(I+1,J,K) - W(I-1,J,K))
    WET = 0.5 * (W(I,J+1,K) - W(I,J-1,K))
    WZE = 0.5 * (W(I,J,K+1) - W(I,J,K-1))
    ! 物理空間方向一階微分
    UX = UXI * XIX(I,J,K) + UET * ETX(I,J,K) + UZE * ZEX(I,J,K)
    UY = UXI * XIY(I,J,K) + UET * ETY(I,J,K) + UZE * ZEY(I,J,K)
    UZ = UXI * XIZ(I,J,K) + UET * ETZ(I,J,K) + UZE * ZEZ(I,J,K)
    VX = VXI * XIX(I,J,K) + VET * ETX(I,J,K) + VZE * ZEX(I,J,K)
    VY = VXI * XIY(I,J,K) + VET * ETY(I,J,K) + VZE * ZEY(I,J,K)
    VZ = VXI * XIZ(I,J,K) + VET * ETZ(I,J,K) + VZE * ZEZ(I,J,K)
    WX = WXI * XIX(I,J,K) + WET * ETX(I,J,K) + WZE * ZEX(I,J,K)
    WY = WXI * XIY(I,J,K) + WET * ETY(I,J,K) + WZE * ZEY(I,J,K)
    WZ = WXI * XIZ(I,J,K) + WET * ETZ(I,J,K) + WZE * ZEZ(I,J,K)
    ! 生成項と散逸項 ---------------------------------------------------
    P11 = - 2.0 * (uu(I,J,K) * UX + uv(I,J,K) * UY + wu(I,J,K) * UZ)
    P22 = - 2.0 * (uv(I,J,K) * VX + vv(I,J,K) * VY + vw(I,J,K) * VZ)
    P33 = - 2.0 * (wu(I,J,K) * WX + vw(I,J,K) * WY + ww(I,J,K) * WZ)
    P12 = - (uv(I,J,K) * UX + uu(I,J,K) * VX) &
    &     - (vv(I,J,K) * UY + uv(I,J,K) * VY) &
    &     - (vw(I,J,K) * UZ + wu(I,J,K) * VZ)
    P23 = - (wu(I,J,K) * VX + uv(I,J,K) * WX) &
    &     - (vw(I,J,K) * VY + vv(I,J,K) * WY) &
    &     - (ww(I,J,K) * VZ + vw(I,J,K) * WZ)
    P31 = - (uu(I,J,K) * WX + wu(I,J,K) * UX) &
    &     - (uv(I,J,K) * WY + vw(I,J,K) * UY) &
    &     - (wu(I,J,K) * WZ + ww(I,J,K) * UZ)
    P_k = 0.5 * (P11 + P22 + P33)
    e_k = EPS(I,J,K)
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      TMP = TwoThird * e_k / PELIM
      P11 = MIN(TMP, MAX(-TMP, P11))
      P22 = MIN(TMP, MAX(-TMP, P22))
      P33 = MIN(TMP, MAX(-TMP, P33))
      P_k = 0.5 * (P11 + P22 + P33)
    ENDIF
    TMP   = TwoThird * e_k
    eps11 = TMP
    eps22 = TMP
    eps33 = TMP
    eps12 = 0.0
    eps23 = 0.0
    eps31 = 0.0
    ! 歪み速度 ---------------------------------------------------------
    TMP = (UX + VY + WZ) * OneThird
    S11 = UX - TMP
    S22 = VY - TMP
    S33 = WZ - TMP
    S12 = 0.5 * (UY + VX)
    S23 = 0.5 * (VZ + WY)
    S31 = 0.5 * (WX + UZ)
    SS2 = S11 * S11 + S22 * S22 + S33 * S33 &
    &   + 2.0 * (S12 * S12 + S23 * S23 + S31 * S31)
    ! 渦度 -------------------------------------------------------------
    O12 = 0.5 * (UY - VX)
    O23 = 0.5 * (VZ - WY)
    O31 = 0.5 * (WX - UZ)
    ! 非等方レイノルズ応力テンソル b_ij --------------------------------
    IF(AK .GT. ZERO) THEN
      b11 = 0.5 * uu(I,J,K) / AK - OneThird
      b22 = 0.5 * vv(I,J,K) / AK - OneThird
      b33 = 0.5 * ww(I,J,K) / AK - OneThird
      b12 = 0.5 * uv(I,J,K) / AK
      b23 = 0.5 * vw(I,J,K) / AK
      b31 = 0.5 * wu(I,J,K) / AK
      b11 = MIN(TwoThird, MAX(- OneThird, b11))
      b22 = MIN(TwoThird, MAX(- OneThird, b22))
      b33 = MIN(TwoThird, MAX(- OneThird, b33))
      b12 = MIN(0.5, MAX(-0.5, b12))
      b23 = MIN(0.5, MAX(-0.5, b23))
      b31 = MIN(0.5, MAX(-0.5, b31))
    ELSE
      b11 =-OneThird
      b22 =-OneThird
      b33 =-OneThird
      b12 = 0.0
      b23 = 0.0
      b31 = 0.0
    ENDIF
    ! 再分配項 ---------------------------------------------------------
    TMPBB = b11 * b11 + b22 * b22 + b33 * b33 &
    &     + 2.0 * (b12 * b12 + b23 * b23 + b31 * b31)
    TMPBS = b11 * S11 + b22 * S22 + b33 * S33 &
    &     + 2.0 * (b12 * S12 + b23 * S23 + b31 * S31)
    TMP1  = - (C_1 * e_k + C_1a * P_k)
    TMP3  = (C_3 - C_3a * SQRT(TMPBB)) * AK
    phi11 = TMP1 * b11 &
    &     + C_2 * e_k * ( &
    &       b11 * b11 + b31 * b31 + b12 * b12 - OneThird * TMPBB &
    &     ) + TMP3 * S11 &
    &     + C_4 * AK * ( &
    &       2.0 * (b11 * S11 + b31 * S31 + b12 * S12) &
    &     - TwoThird * TMPBS &
    &     ) &
    &     + C_5 * AK * 2.0 * (b12 * O12 - b31 * O31)
    phi22 = TMP1 * b22 &
    &     + C_2 * e_k * ( &
    &       b22 * b22 + b12 * b12 + b23 * b23 - OneThird * TMPBB &
    &     ) + TMP3 * S22 &
    &     + C_4 * AK * ( &
    &       2.0 * (b22 * S22 + b12 * S12 + b23 * S23) &
    &     - TwoThird * TMPBS &
    &     ) &
    &     + C_5 * AK * 2.0 * (b23 * O23 - b12 * O12)
    phi33 = TMP1 * b33 &
    &     + C_2 * e_k * ( &
    &       b33 * b33 + b23 * b23 + b31 * b31 - OneThird * TMPBB &
    &     ) + TMP3 * S33 &
    &     + C_4 * AK * ( &
    &       2.0 * (b33 * S33 + b23 * S23 + b31 * S31) &
    &     - TwoThird * TMPBS &
    &     ) &
    &     + C_5 * AK * 2.0 * (b31 * O31 - b23 * O23)
    phi12 = TMP1 * b12 &
    &     + C_2 * e_k * (b12 * (b11 + b22) + b23 * b31) &
    &     + TMP3 * S12 &
    &     + C_4 * AK * ( &
    &       (b11 + b22) * S12 + b12 * (S11 + S22) &
    &     + b23 * S31 + b31 * S23 &
    &     ) &
    &     - C_5 * AK * ((b11 - b22) * O12 + b23 * O31 - b31 * O23)
    phi23 = TMP1 * b23 &
    &     + C_2 * e_k * (b23 * (b22 + b33) + b31 * b12) &
    &     + TMP3 * S23 &
    &     + C_4 * AK * ( &
    &       (b22 + b33) * S23 + b23 * (S22 + S33) &
    &     + b31 * S12 + b12 * S31 &
    &     ) &
    &     - C_5 * AK * ((b22 - b33) * O23 + b31 * O12 - b12 * O31)
    phi31 = TMP1 * b31 &
    &     + C_2 * e_k * (b31 * (b33 + b11) + b12 * b23) &
    &     + TMP3 * S31 &
    &     + C_4 * AK * ( &
    &       (b33 + b11) * S31 + b31 * (S33 + S11) &
    &     + b12 * S23 + b23 * S12 &
    &     ) &
    &     - C_5 * AK * ((b33 - b11) * O31 + b12 * O23 - b23 * O12)
    ! 渦粘性係数の計算 -------------------------------------------------
    IF(e_k .GT. ZERO) THEN
      AMUT(I,J,K) = 0.09 * AK**2 / e_k * RHO(I,J,K)
    ENDIF
    ! 生成項、圧力再分配項、散逸項の和 ---------------------------------
    QH1  = RHO(I,J,K) / AJA(I,J,K)
    nu   = AMU(I,J,K) / RHO(I,J,K)
    AEAK = MIN( &
    &        MAX(SQRT(2.0 * SS2), SQRT(e_k / nu)), &
    &        e_k / MAX(ZERO, AK) &
    &    )
    DQP(I,J,K, 6) = QH1 * (P11 - eps11 + phi11)
    DQP(I,J,K, 7) = QH1 * (P22 - eps22 + phi22)
    DQP(I,J,K, 8) = QH1 * (P33 - eps33 + phi33)
    DQP(I,J,K, 9) = QH1 * (P12 - eps12 + phi12)
    DQP(I,J,K,10) = QH1 * (P23 - eps23 + phi23)
    DQP(I,J,K,11) = QH1 * (P31 - eps31 + phi31)
    DQP(I,J,K,12) = QH1 * AEAK * (C_e1 * P_k - C_e2 * e_k)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH = 0.0
  SH = 0.0
  TH = 0.0
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL DIFFXI
  CALL DIFFET
  CALL DIFFZE
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I,J-1,K) .GT. 0.0) .AND. &
  &   (RHO(I-1,J,K) .GT. 0.0) .AND. (RHO(I,J,K)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J,K) .GT. 0.0) .AND. (RHO(I,J+1,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,K,L) = (-RH(I-1,J  ,K  ,L) + RH(I,J,K,L) ) &
    &            + (-SH(I  ,J-1,K  ,L) + SH(I,J,K,L) ) &
    &            + (-TH(I  ,J  ,K-1,L) + TH(I,J,K,L) )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項の計算(xi方向差分用)                                    ****
!***********************************************************************
SUBROUTINE DIFFXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: uuXI, uuET, uuZE, uuX, uuY, uuZ
  REAL    :: vvXI, vvET, vvZE, vvX, vvY, vvZ
  REAL    :: wwXI, wwET, wwZE, wwX, wwY, wwZ
  REAL    :: uvXI, uvET, uvZE, uvX, uvY, uvZ
  REAL    :: vwXI, vwET, vwZE, vwX, vwY, vwZ
  REAL    :: wuXI, wuET, wuZE, wuX, wuY, wuZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu
  REAL    :: S11, S22, S33, S12, S23, S31, SS
  REAL    :: C_p, AKAE
  REAL    :: R02M, S02M, T02M, R02T, S02T, T02T
  REAL    :: R03M, S03M, T03M, R03T, S03T, T03T
  REAL    :: R04M, S04M, T04M, R04T, S04T, T04T
  REAL    :: R05M, S05M, T05M, R05T, S05T, T05T
  REAL    :: R06M, S06M, T06M, R06T, S06T, T06T
  REAL    :: R07M, S07M, T07M, R07T, S07T, T07T
  REAL    :: R08M, S08M, T08M, R08T, S08T, T08T
  REAL    :: R09M, S09M, T09M, R09T, S09T, T09T
  REAL    :: R10M, S10M, T10M, R10T, S10T, T10T
  REAL    :: R11M, S11M, T11M, R11T, S11T, T11T
  REAL    :: R12M, S12M, T12M, R12T, S12T, T12T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  uuXI, uuET, uuZE, uuX, uuY, uuZ, &
  !$  vvXI, vvET, vvZE, vvX, vvY, vvZ, &
  !$  wwXI, wwET, wwZE, wwX, wwY, wwZ, &
  !$  uvXI, uvET, uvZE, uvX, uvY, uvZ, &
  !$  vwXI, vwET, vwZE, vwX, vwY, vwZ, &
  !$  wuXI, wuET, wuZE, wuX, wuY, wuZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, S11, S22, S33, S12, S23, S31, SS, C_p, AKAE, &
  !$  R02M, S02M, T02M, R02T, S02T, T02T, &
  !$  R03M, S03M, T03M, R03T, S03T, T03T, &
  !$  R04M, S04M, T04M, R04T, S04T, T04T, &
  !$  R05M, S05M, T05M, R05T, S05T, T05T, &
  !$  R06M, S06M, T06M, R06T, S06T, T06T, &
  !$  R07M, S07M, T07M, R07T, S07T, T07T, &
  !$  R08M, S08M, T08M, R08T, S08T, T08T, &
  !$  R09M, S09M, T09M, R09T, S09T, T09T, &
  !$  R10M, S10M, T10M, R10T, S10T, T10T, &
  !$  R11M, S11M, T11M, R11T, S11T, T11T, &
  !$  R12M, S12M, T12M, R12T, S12T, T12T, &
  !$  QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF((RHO(I,J,K) .GT. 0.0) .AND. (RHO(I+1,J,K) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I+1,J,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I+1,J,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I+1,J,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I+1,J,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I+1,J,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I+1,J,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I+1,J,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I+1,J,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I+1,J,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I+1,J,K))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J,K) + RHO(I+1,J,K))
    UM   = 0.5 * (  U(I,J,K) +   U(I+1,J,K))
    VM   = 0.5 * (  V(I,J,K) +   V(I+1,J,K))
    WM   = 0.5 * (  W(I,J,K) +   W(I+1,J,K))
    uuM  = 0.5 * ( uu(I,J,K) +  uu(I+1,J,K))
    vvM  = 0.5 * ( vv(I,J,K) +  vv(I+1,J,K))
    wwM  = 0.5 * ( ww(I,J,K) +  ww(I+1,J,K))
    uvM  = 0.5 * ( uv(I,J,K) +  uv(I+1,J,K))
    vwM  = 0.5 * ( vw(I,J,K) +  vw(I+1,J,K))
    wuM  = 0.5 * ( wu(I,J,K) +  wu(I+1,J,K))
    EPSM = 0.5 * (EPS(I,J,K) + EPS(I+1,J,K))
    AMUM = 0.5 * (AMU(I,J,K) + AMU(I+1,J,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = - U(I,J,K) + U(I+1,J,K)
    UET = 0.25 * ( - U(I,J-1,K) - U(I+1,J-1,K) &
    &              + U(I,J+1,K) + U(I+1,J+1,K) )
    UZE = 0.25 * ( - U(I,J,K-1) - U(I+1,J,K-1) &
    &              + U(I,J,K+1) + U(I+1,J,K+1) )
    VXI = - V(I,J,K) + V(I+1,J,K)
    VET = 0.25 * ( - V(I,J-1,K) - V(I+1,J-1,K) &
    &              + V(I,J+1,K) + V(I+1,J+1,K) )
    VZE = 0.25 * ( - V(I,J,K-1) - V(I+1,J,K-1) &
    &              + V(I,J,K+1) + V(I+1,J,K+1) )
    WXI = - W(I,J,K) + W(I+1,J,K)
    WET = 0.25 * ( - W(I,J-1,K) - W(I+1,J-1,K) &
    &              + W(I,J+1,K) + W(I+1,J+1,K) )
    WZE = 0.25 * ( - W(I,J,K-1) - W(I+1,J,K-1) &
    &              + W(I,J,K+1) + W(I+1,J,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = - T(I,J,K) + T(I+1,J,K)
    TET = 0.25 * ( - T(I,J-1,K) - T(I+1,J-1,K) &
    &              + T(I,J+1,K) + T(I+1,J+1,K) )
    TZE = 0.25 * ( - T(I,J,K-1) - T(I+1,J,K-1) &
    &              + T(I,J,K+1) + T(I+1,J,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = - uu(I,J,K) + uu(I+1,J,K)
    uuET = 0.25 * ( - uu(I,J-1,K) - uu(I+1,J-1,K) &
    &               + uu(I,J+1,K) + uu(I+1,J+1,K) )
    uuZE = 0.25 * ( - uu(I,J,K-1) - uu(I+1,J,K-1) &
    &               + uu(I,J,K+1) + uu(I+1,J,K+1) )
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM + uuZE * ZEXM
    uuY = uuXI * XIYM + uuET * ETYM + uuZE * ZEYM
    uuZ = uuXI * XIZM + uuET * ETZM + uuZE * ZEZM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = - vv(I,J,K) + vv(I+1,J,K)
    vvET = 0.25 * ( - vv(I,J-1,K) - vv(I+1,J-1,K) &
    &               + vv(I,J+1,K) + vv(I+1,J+1,K) )
    vvZE = 0.25 * ( - vv(I,J,K-1) - vv(I+1,J,K-1) &
    &               + vv(I,J,K+1) + vv(I+1,J,K+1) )
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM + vvZE * ZEXM
    vvY = vvXI * XIYM + vvET * ETYM + vvZE * ZEYM
    vvZ = vvXI * XIZM + vvET * ETZM + vvZE * ZEZM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = - ww(I,J,K) + ww(I+1,J,K)
    wwET = 0.25 * ( - ww(I,J-1,K) - ww(I+1,J-1,K) &
    &               + ww(I,J+1,K) + ww(I+1,J+1,K) )
    wwZE = 0.25 * ( - ww(I,J,K-1) - ww(I+1,J,K-1) &
    &               + ww(I,J,K+1) + ww(I+1,J,K+1) )
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM + wwZE * ZEXM
    wwY = wwXI * XIYM + wwET * ETYM + wwZE * ZEYM
    wwZ = wwXI * XIZM + wwET * ETZM + wwZE * ZEZM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = - uv(I,J,K) + uv(I+1,J,K)
    uvET = 0.25 * ( - uv(I,J-1,K) - uv(I+1,J-1,K) &
    &               + uv(I,J+1,K) + uv(I+1,J+1,K) )
    uvZE = 0.25 * ( - uv(I,J,K-1) - uv(I+1,J,K-1) &
    &               + uv(I,J,K+1) + uv(I+1,J,K+1) )
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM + uvZE * ZEXM
    uvY = uvXI * XIYM + uvET * ETYM + uvZE * ZEYM
    uvZ = uvXI * XIZM + uvET * ETZM + uvZE * ZEZM
    ! vw の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vwXI = - vw(I,J,K) + vw(I+1,J,K)
    vwET = 0.25 * ( - vw(I,J-1,K) - vw(I+1,J-1,K) &
    &               + vw(I,J+1,K) + vw(I+1,J+1,K) )
    vwZE = 0.25 * ( - vw(I,J,K-1) - vw(I+1,J,K-1) &
    &               + vw(I,J,K+1) + vw(I+1,J,K+1) )
    ! 物理空間方向一階微分
    vwX = vwXI * XIXM + vwET * ETXM + vwZE * ZEXM
    vwY = vwXI * XIYM + vwET * ETYM + vwZE * ZEYM
    vwZ = vwXI * XIZM + vwET * ETZM + vwZE * ZEZM
    ! wu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wuXI = - wu(I,J,K) + wu(I+1,J,K)
    wuET = 0.25 * ( - wu(I,J-1,K) - wu(I+1,J-1,K) &
    &               + wu(I,J+1,K) + wu(I+1,J+1,K) )
    wuZE = 0.25 * ( - wu(I,J,K-1) - wu(I+1,J,K-1) &
    &               + wu(I,J,K+1) + wu(I+1,J,K+1) )
    ! 物理空間方向一階微分
    wuX = wuXI * XIXM + wuET * ETXM + wuZE * ZEXM
    wuY = wuXI * XIYM + wuET * ETYM + wuZE * ZEYM
    wuZ = wuXI * XIZM + wuET * ETZM + wuZE * ZEZM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = - EPS(I,J,K) + EPS(I+1,J,K)
    AEET = 0.25 * ( - EPS(I,J-1,K) - EPS(I+1,J-1,K) &
    &               + EPS(I,J+1,K) + EPS(I+1,J+1,K) )
    AEZE = 0.25 * ( - EPS(I,J,K-1) - EPS(I+1,J,K-1) &
    &               + EPS(I,J,K+1) + EPS(I+1,J,K+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY + WZ) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 = 2.0 * WZ - SS
    S12 = UY + VX
    S23 = VZ + WY
    S31 = WX + UZ
    ! U の拡散 ---------------------------------------------------------
    R02M = nu * S11
    S02M = nu * S12
    T02M = nu * S31
    R02T =-uuM
    S02T =-uvM
    T02T =-wuM
    ! V の拡散 ---------------------------------------------------------
    R03M = nu * S12
    S03M = nu * S22
    T03M = nu * S23
    R03T =-uvM
    S03T =-vvM
    T03T =-vwM
    ! W の拡散 ---------------------------------------------------------
    R04M = nu * S31
    S04M = nu * S23
    T04M = nu * S33
    R04T =-wuM
    S04T =-vwM
    T04T =-wwM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R05M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM + S31 * WM)
    S05M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM + S23 * WM)
    T05M = C_p * nu / PR * TZ + nu * (S31 * UM + S23 * VM + S33 * WM)
    R05T = AKAE * C_se * C_p * (uuM * TX + uvM * TY + wuM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uuM * (uuX + vvX + wwX) &
    &    + uvM * (uuY + vvY + wwY) &
    &    + wuM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uuM * UM + uvM * VM + wuM * WM)
    S05T = AKAE * C_se * C_p * (uvM * TX + vvM * TY + vwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uvM * (uuX + vvX + wwX) &
    &    + vvM * (uuY + vvY + wwY) &
    &    + vwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uvM * UM + vvM * VM + vwM * WM)
    T05T = AKAE * C_se * C_p * (wuM * TX + vwM * TY + wwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      wuM * (uuX + vvX + wwX) &
    &    + vwM * (uuY + vvY + wwY) &
    &    + wwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (wuM * UM + vwM * VM + wwM * WM)
    ! uu の拡散 --------------------------------------------------------
    R06M = nu * uuX
    S06M = nu * uuY
    T06M = nu * uuZ
    R06T = C_s * AKAE * (uuM * uuX + uvM * uuY + wuM * uuZ)
    S06T = C_s * AKAE * (uvM * uuX + vvM * uuY + vwM * uuZ)
    T06T = C_s * AKAE * (wuM * uuX + vwM * uuY + wwM * uuZ)
    ! vv の拡散 --------------------------------------------------------
    R07M = nu * vvX
    S07M = nu * vvY
    T07M = nu * vvZ
    R07T = C_s * AKAE * (uuM * vvX + uvM * vvY + wuM * vvZ)
    S07T = C_s * AKAE * (uvM * vvX + vvM * vvY + vwM * vvZ)
    T07T = C_s * AKAE * (wuM * vvX + vwM * vvY + wwM * vvZ)
    ! ww の拡散 --------------------------------------------------------
    R08M = nu * wwX
    S08M = nu * wwY
    T08M = nu * wwZ
    R08T = C_s * AKAE * (uuM * wwX + uvM * wwY + wuM * wwZ)
    S08T = C_s * AKAE * (uvM * wwX + vvM * wwY + vwM * wwZ)
    T08T = C_s * AKAE * (wuM * wwX + vwM * wwY + wwM * wwZ)
    ! uv の拡散 --------------------------------------------------------
    R09M = nu * uvX
    S09M = nu * uvY
    T09M = nu * uvZ
    R09T = C_s * AKAE * (uuM * uvX + uvM * uvY + wuM * uvZ)
    S09T = C_s * AKAE * (uvM * uvX + vvM * uvY + vwM * uvZ)
    T09T = C_s * AKAE * (wuM * uvX + vwM * uvY + wwM * uvZ)
    ! vw の拡散 --------------------------------------------------------
    R10M = nu * vwX
    S10M = nu * vwY
    T10M = nu * vwZ
    R10T = C_s * AKAE * (uuM * vwX + uvM * vwY + wuM * vwZ)
    S10T = C_s * AKAE * (uvM * vwX + vvM * vwY + vwM * vwZ)
    T10T = C_s * AKAE * (wuM * vwX + vwM * vwY + wwM * vwZ)
    ! wu の拡散 --------------------------------------------------------
    R11M = nu * wuX
    S11M = nu * wuY
    T11M = nu * wuZ
    R11T = C_s * AKAE * (uuM * wuX + uvM * wuY + wuM * wuZ)
    S11T = C_s * AKAE * (uvM * wuX + vvM * wuY + vwM * wuZ)
    T11T = C_s * AKAE * (wuM * wuX + vwM * wuY + wwM * wuZ)
    ! epsilon の拡散 ---------------------------------------------------
    R12M = nu * AEX
    S12M = nu * AEY
    T12M = nu * AEZ
    R12T = C_e * AKAE * (uuM * AEX + uvM * AEY + wuM * AEZ)
    S12T = C_e * AKAE * (uvM * AEX + vvM * AEY + vwM * AEZ)
    T12T = C_e * AKAE * (wuM * AEX + vwM * AEY + wwM * AEZ)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    RH(I,J,K, 1) = 0.0
    RH(I,J,K, 2) = QH1 * ( (R02M + R02T) * XIXM &
    &                    + (S02M + S02T) * XIYM &
    &                    + (T02M + T02T) * XIZM )
    RH(I,J,K, 3) = QH1 * ( (R03M + R03T) * XIXM &
    &                    + (S03M + S03T) * XIYM &
    &                    + (T03M + T03T) * XIZM )
    RH(I,J,K, 4) = QH1 * ( (R04M + R04T) * XIXM &
    &                    + (S04M + S04T) * XIYM &
    &                    + (T04M + T04T) * XIZM )
    RH(I,J,K, 5) = QH1 * ( (R05M + R05T) * XIXM &
    &                    + (S05M + S05T) * XIYM &
    &                    + (T05M + T05T) * XIZM )
    RH(I,J,K, 6) = QH1 * ( (R06M + R06T) * XIXM &
    &                    + (S06M + S06T) * XIYM &
    &                    + (T06M + T06T) * XIZM )
    RH(I,J,K, 7) = QH1 * ( (R07M + R07T) * XIXM &
    &                    + (S07M + S07T) * XIYM &
    &                    + (T07M + T07T) * XIZM )
    RH(I,J,K, 8) = QH1 * ( (R08M + R08T) * XIXM &
    &                    + (S08M + S08T) * XIYM &
    &                    + (T08M + T08T) * XIZM )
    RH(I,J,K, 9) = QH1 * ( (R09M + R09T) * XIXM &
    &                    + (S09M + S09T) * XIYM &
    &                    + (T09M + T09T) * XIZM )
    RH(I,J,K,10) = QH1 * ( (R10M + R10T) * XIXM &
    &                    + (S10M + S10T) * XIYM &
    &                    + (T10M + T10T) * XIZM )
    RH(I,J,K,11) = QH1 * ( (R11M + R11T) * XIXM &
    &                    + (S11M + S11T) * XIYM &
    &                    + (T11M + T11T) * XIZM )
    RH(I,J,K,12) = QH1 * ( (R12M + R12T) * XIXM &
    &                    + (S12M + S12T) * XIYM &
    &                    + (T12M + T12T) * XIZM )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFXI
!***********************************************************************
!**** 拡散項の計算(eta方向差分用)                                   ****
!***********************************************************************
SUBROUTINE DIFFET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: uuXI, uuET, uuZE, uuX, uuY, uuZ
  REAL    :: vvXI, vvET, vvZE, vvX, vvY, vvZ
  REAL    :: wwXI, wwET, wwZE, wwX, wwY, wwZ
  REAL    :: uvXI, uvET, uvZE, uvX, uvY, uvZ
  REAL    :: vwXI, vwET, vwZE, vwX, vwY, vwZ
  REAL    :: wuXI, wuET, wuZE, wuX, wuY, wuZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu
  REAL    :: S11, S22, S33, S12, S23, S31, SS
  REAL    :: C_p, AKAE
  REAL    :: R02M, S02M, T02M, R02T, S02T, T02T
  REAL    :: R03M, S03M, T03M, R03T, S03T, T03T
  REAL    :: R04M, S04M, T04M, R04T, S04T, T04T
  REAL    :: R05M, S05M, T05M, R05T, S05T, T05T
  REAL    :: R06M, S06M, T06M, R06T, S06T, T06T
  REAL    :: R07M, S07M, T07M, R07T, S07T, T07T
  REAL    :: R08M, S08M, T08M, R08T, S08T, T08T
  REAL    :: R09M, S09M, T09M, R09T, S09T, T09T
  REAL    :: R10M, S10M, T10M, R10T, S10T, T10T
  REAL    :: R11M, S11M, T11M, R11T, S11T, T11T
  REAL    :: R12M, S12M, T12M, R12T, S12T, T12T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  uuXI, uuET, uuZE, uuX, uuY, uuZ, &
  !$  vvXI, vvET, vvZE, vvX, vvY, vvZ, &
  !$  wwXI, wwET, wwZE, wwX, wwY, wwZ, &
  !$  uvXI, uvET, uvZE, uvX, uvY, uvZ, &
  !$  vwXI, vwET, vwZE, vwX, vwY, vwZ, &
  !$  wuXI, wuET, wuZE, wuX, wuY, wuZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, S11, S22, S33, S12, S23, S31, SS, C_p, AKAE, &
  !$  R02M, S02M, T02M, R02T, S02T, T02T, &
  !$  R03M, S03M, T03M, R03T, S03T, T03T, &
  !$  R04M, S04M, T04M, R04T, S04T, T04T, &
  !$  R05M, S05M, T05M, R05T, S05T, T05T, &
  !$  R06M, S06M, T06M, R06T, S06T, T06T, &
  !$  R07M, S07M, T07M, R07T, S07T, T07T, &
  !$  R08M, S08M, T08M, R08T, S08T, T08T, &
  !$  R09M, S09M, T09M, R09T, S09T, T09T, &
  !$  R10M, S10M, T10M, R10T, S10T, T10T, &
  !$  R11M, S11M, T11M, R11T, S11T, T11T, &
  !$  R12M, S12M, T12M, R12T, S12T, T12T, &
  !$  QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF((RHO(I,J,K) .GT. 0.0) .AND. (RHO(I,J+1,K) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J+1,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J+1,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J+1,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J+1,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J+1,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J+1,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J+1,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J+1,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J+1,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J+1,K))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J,K) + RHO(I,J+1,K))
    UM   = 0.5 * (  U(I,J,K) +   U(I,J+1,K))
    VM   = 0.5 * (  V(I,J,K) +   V(I,J+1,K))
    WM   = 0.5 * (  W(I,J,K) +   W(I,J+1,K))
    uuM  = 0.5 * ( uu(I,J,K) +  uu(I,J+1,K))
    vvM  = 0.5 * ( vv(I,J,K) +  vv(I,J+1,K))
    wwM  = 0.5 * ( ww(I,J,K) +  ww(I,J+1,K))
    uvM  = 0.5 * ( uv(I,J,K) +  uv(I,J+1,K))
    vwM  = 0.5 * ( vw(I,J,K) +  vw(I,J+1,K))
    wuM  = 0.5 * ( wu(I,J,K) +  wu(I,J+1,K))
    EPSM = 0.5 * (EPS(I,J,K) + EPS(I,J+1,K))
    AMUM = 0.5 * (AMU(I,J,K) + AMU(I,J+1,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J+1,K) &
    &              + U(I+1,J,K) + U(I+1,J+1,K) )
    UET = - U(I,J,K) + U(I,J+1,K)
    UZE = 0.25 * ( - U(I,J,K-1) - U(I,J+1,K-1) &
    &              + U(I,J,K+1) + U(I,J+1,K+1) )
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J+1,K) &
    &              + V(I+1,J,K) + V(I+1,J+1,K) )
    VET = - V(I,J,K) + V(I,J+1,K)
    VZE = 0.25 * ( - V(I,J,K-1) - V(I,J+1,K-1) &
    &              + V(I,J,K+1) + V(I,J+1,K+1) )
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J+1,K) &
    &              + W(I+1,J,K) + W(I+1,J+1,K) )
    WET = - W(I,J,K) + W(I,J+1,K)
    WZE = 0.25 * ( - W(I,J,K-1) - W(I,J+1,K-1) &
    &              + W(I,J,K+1) + W(I,J+1,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J+1,K) &
    &              + T(I+1,J,K) + T(I+1,J+1,K) )
    TET = - T(I,J,K) + T(I,J+1,K)
    TZE = 0.25 * ( - T(I,J,K-1) - T(I,J+1,K-1) &
    &              + T(I,J,K+1) + T(I,J+1,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = 0.25 * ( - uu(I-1,J,K) - uu(I-1,J+1,K) &
    &               + uu(I+1,J,K) + uu(I+1,J+1,K) )
    uuET = - uu(I,J,K) + uu(I,J+1,K)
    uuZE = 0.25 * ( - uu(I,J,K-1) - uu(I,J+1,K-1) &
    &               + uu(I,J,K+1) + uu(I,J+1,K+1) )
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM + uuZE * ZEXM
    uuY = uuXI * XIYM + uuET * ETYM + uuZE * ZEYM
    uuZ = uuXI * XIZM + uuET * ETZM + uuZE * ZEZM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = 0.25 * ( - vv(I-1,J,K) - vv(I-1,J+1,K) &
    &               + vv(I+1,J,K) + vv(I+1,J+1,K) )
    vvET = - vv(I,J,K) + vv(I,J+1,K)
    vvZE = 0.25 * ( - vv(I,J,K-1) - vv(I,J+1,K-1) &
    &               + vv(I,J,K+1) + vv(I,J+1,K+1) )
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM + vvZE * ZEXM
    vvY = vvXI * XIYM + vvET * ETYM + vvZE * ZEYM
    vvZ = vvXI * XIZM + vvET * ETZM + vvZE * ZEZM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = 0.25 * ( - ww(I-1,J,K) - ww(I-1,J+1,K) &
    &               + ww(I+1,J,K) + ww(I+1,J+1,K) )
    wwET = - ww(I,J,K) + ww(I,J+1,K)
    wwZE = 0.25 * ( - ww(I,J,K-1) - ww(I,J+1,K-1) &
    &               + ww(I,J,K+1) + ww(I,J+1,K+1) )
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM + wwZE * ZEXM
    wwY = wwXI * XIYM + wwET * ETYM + wwZE * ZEYM
    wwZ = wwXI * XIZM + wwET * ETZM + wwZE * ZEZM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = 0.25 * ( - uv(I-1,J,K) - uv(I-1,J+1,K) &
    &               + uv(I+1,J,K) + uv(I+1,J+1,K) )
    uvET = - uv(I,J,K) + uv(I,J+1,K)
    uvZE = 0.25 * ( - uv(I,J,K-1) - uv(I,J+1,K-1) &
    &               + uv(I,J,K+1) + uv(I,J+1,K+1) )
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM + uvZE * ZEXM
    uvY = uvXI * XIYM + uvET * ETYM + uvZE * ZEYM
    uvZ = uvXI * XIZM + uvET * ETZM + uvZE * ZEZM
    ! vw の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vwXI = 0.25 * ( - vw(I-1,J,K) - vw(I-1,J+1,K) &
    &               + vw(I+1,J,K) + vw(I+1,J+1,K) )
    vwET = - vw(I,J,K) + vw(I,J+1,K)
    vwZE = 0.25 * ( - vw(I,J,K-1) - vw(I,J+1,K-1) &
    &               + vw(I,J,K+1) + vw(I,J+1,K+1) )
    ! 物理空間方向一階微分
    vwX = vwXI * XIXM + vwET * ETXM + vwZE * ZEXM
    vwY = vwXI * XIYM + vwET * ETYM + vwZE * ZEYM
    vwZ = vwXI * XIZM + vwET * ETZM + vwZE * ZEZM
    ! wu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wuXI = 0.25 * ( - wu(I-1,J,K) - wu(I-1,J+1,K) &
    &               + wu(I+1,J,K) + wu(I+1,J+1,K) )
    wuET = - wu(I,J,K) + wu(I,J+1,K)
    wuZE = 0.25 * ( - wu(I,J,K-1) - wu(I,J+1,K-1) &
    &               + wu(I,J,K+1) + wu(I,J+1,K+1) )
    ! 物理空間方向一階微分
    wuX = wuXI * XIXM + wuET * ETXM + wuZE * ZEXM
    wuY = wuXI * XIYM + wuET * ETYM + wuZE * ZEYM
    wuZ = wuXI * XIZM + wuET * ETZM + wuZE * ZEZM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPS(I-1,J,K) - EPS(I-1,J+1,K) &
    &               + EPS(I+1,J,K) + EPS(I+1,J+1,K) )
    AEET = - EPS(I,J,K) + EPS(I,J+1,K)
    AEZE = 0.25 * ( - EPS(I,J,K-1) - EPS(I,J+1,K-1) &
    &               + EPS(I,J,K+1) + EPS(I,J+1,K+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY + WZ) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 = 2.0 * WZ - SS
    S12 = UY + VX
    S23 = VZ + WY
    S31 = WX + UZ
    ! U の拡散 ---------------------------------------------------------
    R02M = nu * S11
    S02M = nu * S12
    T02M = nu * S31
    R02T =-uuM
    S02T =-uvM
    T02T =-wuM
    ! V の拡散 ---------------------------------------------------------
    R03M = nu * S12
    S03M = nu * S22
    T03M = nu * S23
    R03T =-uvM
    S03T =-vvM
    T03T =-vwM
    ! W の拡散 ---------------------------------------------------------
    R04M = nu * S31
    S04M = nu * S23
    T04M = nu * S33
    R04T =-wuM
    S04T =-vwM
    T04T =-wwM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R05M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM + S31 * WM)
    S05M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM + S23 * WM)
    T05M = C_p * nu / PR * TZ + nu * (S31 * UM + S23 * VM + S33 * WM)
    R05T = AKAE * C_se * C_p * (uuM * TX + uvM * TY + wuM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uuM * (uuX + vvX + wwX) &
    &    + uvM * (uuY + vvY + wwY) &
    &    + wuM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uuM * UM + uvM * VM + wuM * WM)
    S05T = AKAE * C_se * C_p * (uvM * TX + vvM * TY + vwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uvM * (uuX + vvX + wwX) &
    &    + vvM * (uuY + vvY + wwY) &
    &    + vwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uvM * UM + vvM * VM + vwM * WM)
    T05T = AKAE * C_se * C_p * (wuM * TX + vwM * TY + wwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      wuM * (uuX + vvX + wwX) &
    &    + vwM * (uuY + vvY + wwY) &
    &    + wwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (wuM * UM + vwM * VM + wwM * WM)
    ! uu の拡散 --------------------------------------------------------
    R06M = nu * uuX
    S06M = nu * uuY
    T06M = nu * uuZ
    R06T = C_s * AKAE * (uuM * uuX + uvM * uuY + wuM * uuZ)
    S06T = C_s * AKAE * (uvM * uuX + vvM * uuY + vwM * uuZ)
    T06T = C_s * AKAE * (wuM * uuX + vwM * uuY + wwM * uuZ)
    ! vv の拡散 --------------------------------------------------------
    R07M = nu * vvX
    S07M = nu * vvY
    T07M = nu * vvZ
    R07T = C_s * AKAE * (uuM * vvX + uvM * vvY + wuM * vvZ)
    S07T = C_s * AKAE * (uvM * vvX + vvM * vvY + vwM * vvZ)
    T07T = C_s * AKAE * (wuM * vvX + vwM * vvY + wwM * vvZ)
    ! ww の拡散 --------------------------------------------------------
    R08M = nu * wwX
    S08M = nu * wwY
    T08M = nu * wwZ
    R08T = C_s * AKAE * (uuM * wwX + uvM * wwY + wuM * wwZ)
    S08T = C_s * AKAE * (uvM * wwX + vvM * wwY + vwM * wwZ)
    T08T = C_s * AKAE * (wuM * wwX + vwM * wwY + wwM * wwZ)
    ! uv の拡散 --------------------------------------------------------
    R09M = nu * uvX
    S09M = nu * uvY
    T09M = nu * uvZ
    R09T = C_s * AKAE * (uuM * uvX + uvM * uvY + wuM * uvZ)
    S09T = C_s * AKAE * (uvM * uvX + vvM * uvY + vwM * uvZ)
    T09T = C_s * AKAE * (wuM * uvX + vwM * uvY + wwM * uvZ)
    ! vw の拡散 --------------------------------------------------------
    R10M = nu * vwX
    S10M = nu * vwY
    T10M = nu * vwZ
    R10T = C_s * AKAE * (uuM * vwX + uvM * vwY + wuM * vwZ)
    S10T = C_s * AKAE * (uvM * vwX + vvM * vwY + vwM * vwZ)
    T10T = C_s * AKAE * (wuM * vwX + vwM * vwY + wwM * vwZ)
    ! wu の拡散 --------------------------------------------------------
    R11M = nu * wuX
    S11M = nu * wuY
    T11M = nu * wuZ
    R11T = C_s * AKAE * (uuM * wuX + uvM * wuY + wuM * wuZ)
    S11T = C_s * AKAE * (uvM * wuX + vvM * wuY + vwM * wuZ)
    T11T = C_s * AKAE * (wuM * wuX + vwM * wuY + wwM * wuZ)
    ! epsilon の拡散 ---------------------------------------------------
    R12M = nu * AEX
    S12M = nu * AEY
    T12M = nu * AEZ
    R12T = C_e * AKAE * (uuM * AEX + uvM * AEY + wuM * AEZ)
    S12T = C_e * AKAE * (uvM * AEX + vvM * AEY + vwM * AEZ)
    T12T = C_e * AKAE * (wuM * AEX + vwM * AEY + wwM * AEZ)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    SH(I,J,K, 1) = 0.0
    SH(I,J,K, 2) = QH1 * ( (R02M + R02T) * ETXM &
    &                    + (S02M + S02T) * ETYM &
    &                    + (T02M + T02T) * ETZM )
    SH(I,J,K, 3) = QH1 * ( (R03M + R03T) * ETXM &
    &                    + (S03M + S03T) * ETYM &
    &                    + (T03M + T03T) * ETZM )
    SH(I,J,K, 4) = QH1 * ( (R04M + R04T) * ETXM &
    &                    + (S04M + S04T) * ETYM &
    &                    + (T04M + T04T) * ETZM )
    SH(I,J,K, 5) = QH1 * ( (R05M + R05T) * ETXM &
    &                    + (S05M + S05T) * ETYM &
    &                    + (T05M + T05T) * ETZM )
    SH(I,J,K, 6) = QH1 * ( (R06M + R06T) * ETXM &
    &                    + (S06M + S06T) * ETYM &
    &                    + (T06M + T06T) * ETZM )
    SH(I,J,K, 7) = QH1 * ( (R07M + R07T) * ETXM &
    &                    + (S07M + S07T) * ETYM &
    &                    + (T07M + T07T) * ETZM )
    SH(I,J,K, 8) = QH1 * ( (R08M + R08T) * ETXM &
    &                    + (S08M + S08T) * ETYM &
    &                    + (T08M + T08T) * ETZM )
    SH(I,J,K, 9) = QH1 * ( (R09M + R09T) * ETXM &
    &                    + (S09M + S09T) * ETYM &
    &                    + (T09M + T09T) * ETZM )
    SH(I,J,K,10) = QH1 * ( (R10M + R10T) * ETXM &
    &                    + (S10M + S10T) * ETYM &
    &                    + (T10M + T10T) * ETZM )
    SH(I,J,K,11) = QH1 * ( (R11M + R11T) * ETXM &
    &                    + (S11M + S11T) * ETYM &
    &                    + (T11M + T11T) * ETZM )
    SH(I,J,K,12) = QH1 * ( (R12M + R12T) * ETXM &
    &                    + (S12M + S12T) * ETYM &
    &                    + (T12M + T12T) * ETZM )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFET
!***********************************************************************
!**** 拡散項の計算(zeta方向差分用)                                  ****
!***********************************************************************
SUBROUTINE DIFFZE
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: uuXI, uuET, uuZE, uuX, uuY, uuZ
  REAL    :: vvXI, vvET, vvZE, vvX, vvY, vvZ
  REAL    :: wwXI, wwET, wwZE, wwX, wwY, wwZ
  REAL    :: uvXI, uvET, uvZE, uvX, uvY, uvZ
  REAL    :: vwXI, vwET, vwZE, vwX, vwY, vwZ
  REAL    :: wuXI, wuET, wuZE, wuX, wuY, wuZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu
  REAL    :: S11, S22, S33, S12, S23, S31, SS
  REAL    :: C_p, AKAE
  REAL    :: R02M, S02M, T02M, R02T, S02T, T02T
  REAL    :: R03M, S03M, T03M, R03T, S03T, T03T
  REAL    :: R04M, S04M, T04M, R04T, S04T, T04T
  REAL    :: R05M, S05M, T05M, R05T, S05T, T05T
  REAL    :: R06M, S06M, T06M, R06T, S06T, T06T
  REAL    :: R07M, S07M, T07M, R07T, S07T, T07T
  REAL    :: R08M, S08M, T08M, R08T, S08T, T08T
  REAL    :: R09M, S09M, T09M, R09T, S09T, T09T
  REAL    :: R10M, S10M, T10M, R10T, S10T, T10T
  REAL    :: R11M, S11M, T11M, R11T, S11T, T11T
  REAL    :: R12M, S12M, T12M, R12T, S12T, T12T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  uuXI, uuET, uuZE, uuX, uuY, uuZ, &
  !$  vvXI, vvET, vvZE, vvX, vvY, vvZ, &
  !$  wwXI, wwET, wwZE, wwX, wwY, wwZ, &
  !$  uvXI, uvET, uvZE, uvX, uvY, uvZ, &
  !$  vwXI, vwET, vwZE, vwX, vwY, vwZ, &
  !$  wuXI, wuET, wuZE, wuX, wuY, wuZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, S11, S22, S33, S12, S23, S31, SS, C_p, AKAE, &
  !$  R02M, S02M, T02M, R02T, S02T, T02T, &
  !$  R03M, S03M, T03M, R03T, S03T, T03T, &
  !$  R04M, S04M, T04M, R04T, S04T, T04T, &
  !$  R05M, S05M, T05M, R05T, S05T, T05T, &
  !$  R06M, S06M, T06M, R06T, S06T, T06T, &
  !$  R07M, S07M, T07M, R07T, S07T, T07T, &
  !$  R08M, S08M, T08M, R08T, S08T, T08T, &
  !$  R09M, S09M, T09M, R09T, S09T, T09T, &
  !$  R10M, S10M, T10M, R10T, S10T, T10T, &
  !$  R11M, S11M, T11M, R11T, S11T, T11T, &
  !$  R12M, S12M, T12M, R12T, S12T, T12T, &
  !$  QH1 &
  !$)
  DO K = KS    , KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF((RHO(I,J,K) .GT. 0.0) .AND. (RHO(I,J,K+1) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J,K+1))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J,K+1))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J,K+1))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J,K+1))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J,K+1))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J,K+1))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J,K+1))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J,K+1))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J,K+1))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J,K+1))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J,K) + RHO(I,J,K+1))
    UM   = 0.5 * (  U(I,J,K) +   U(I,J,K+1))
    VM   = 0.5 * (  V(I,J,K) +   V(I,J,K+1))
    WM   = 0.5 * (  W(I,J,K) +   W(I,J,K+1))
    uuM  = 0.5 * ( uu(I,J,K) +  uu(I,J,K+1))
    vvM  = 0.5 * ( vv(I,J,K) +  vv(I,J,K+1))
    wwM  = 0.5 * ( ww(I,J,K) +  ww(I,J,K+1))
    uvM  = 0.5 * ( uv(I,J,K) +  uv(I,J,K+1))
    vwM  = 0.5 * ( vw(I,J,K) +  vw(I,J,K+1))
    wuM  = 0.5 * ( wu(I,J,K) +  wu(I,J,K+1))
    EPSM = 0.5 * (EPS(I,J,K) + EPS(I,J,K+1))
    AMUM = 0.5 * (AMU(I,J,K) + AMU(I,J,K+1))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J,K+1) &
    &              + U(I+1,J,K) + U(I+1,J,K+1) )
    UET = 0.25 * ( - U(I,J-1,K) - U(I,J-1,K+1) &
    &              + U(I,J+1,K) + U(I,J+1,K+1) )
    UZE = - U(I,J,K) + U(I,J,K+1)
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J,K+1) &
    &              + V(I+1,J,K) + V(I+1,J,K+1) )
    VET = 0.25 * ( - V(I,J-1,K) - V(I,J-1,K+1) &
    &              + V(I,J+1,K) + V(I,J+1,K+1) )
    VZE = - V(I,J,K) + V(I,J,K+1)
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J,K+1) &
    &              + W(I+1,J,K) + W(I+1,J,K+1) )
    WET = 0.25 * ( - W(I,J-1,K) - W(I,J-1,K+1) &
    &              + W(I,J+1,K) + W(I,J+1,K+1) )
    WZE = - W(I,J,K) + W(I,J,K+1)
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J,K+1) &
    &              + T(I+1,J,K) + T(I+1,J,K+1) )
    TET = 0.25 * ( - T(I,J-1,K) - T(I,J-1,K+1) &
    &              + T(I,J+1,K) + T(I,J+1,K+1) )
    TZE = - T(I,J,K) + T(I,J,K+1)
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = 0.25 * ( - uu(I-1,J,K) - uu(I-1,J,K+1) &
    &               + uu(I+1,J,K) + uu(I+1,J,K+1) )
    uuET = 0.25 * ( - uu(I,J-1,K) - uu(I,J-1,K+1) &
    &               + uu(I,J+1,K) + uu(I,J+1,K+1) )
    uuZE = - uu(I,J,K) + uu(I,J,K+1)
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM + uuZE * ZEXM
    uuY = uuXI * XIYM + uuET * ETYM + uuZE * ZEYM
    uuZ = uuXI * XIZM + uuET * ETZM + uuZE * ZEZM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = 0.25 * ( - vv(I-1,J,K) - vv(I-1,J,K+1) &
    &               + vv(I+1,J,K) + vv(I+1,J,K+1) )
    vvET = 0.25 * ( - vv(I,J-1,K) - vv(I,J-1,K+1) &
    &               + vv(I,J+1,K) + vv(I,J+1,K+1) )
    vvZE = - vv(I,J,K) + vv(I,J,K+1)
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM + vvZE * ZEXM
    vvY = vvXI * XIYM + vvET * ETYM + vvZE * ZEYM
    vvZ = vvXI * XIZM + vvET * ETZM + vvZE * ZEZM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = 0.25 * ( - ww(I-1,J,K) - ww(I-1,J,K+1) &
    &               + ww(I+1,J,K) + ww(I+1,J,K+1) )
    wwET = 0.25 * ( - ww(I,J-1,K) - ww(I,J-1,K+1) &
    &               + ww(I,J+1,K) + ww(I,J+1,K+1) )
    wwZE = - ww(I,J,K) + ww(I,J,K+1)
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM + wwZE * ZEXM
    wwY = wwXI * XIYM + wwET * ETYM + wwZE * ZEYM
    wwZ = wwXI * XIZM + wwET * ETZM + wwZE * ZEZM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = 0.25 * ( - uv(I-1,J,K) - uv(I-1,J,K+1) &
    &               + uv(I+1,J,K) + uv(I+1,J,K+1) )
    uvET = 0.25 * ( - uv(I,J-1,K) - uv(I,J-1,K+1) &
    &               + uv(I,J+1,K) + uv(I,J+1,K+1) )
    uvZE = - uv(I,J,K) + uv(I,J,K+1)
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM + uvZE * ZEXM
    uvY = uvXI * XIYM + uvET * ETYM + uvZE * ZEYM
    uvZ = uvXI * XIZM + uvET * ETZM + uvZE * ZEZM
    ! vw の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vwXI = 0.25 * ( - vw(I-1,J,K) - vw(I-1,J,K+1) &
    &               + vw(I+1,J,K) + vw(I+1,J,K+1) )
    vwET = 0.25 * ( - vw(I,J-1,K) - vw(I,J-1,K+1) &
    &               + vw(I,J+1,K) + vw(I,J+1,K+1) )
    vwZE = - vw(I,J,K) + vw(I,J,K+1)
    ! 物理空間方向一階微分
    vwX = vwXI * XIXM + vwET * ETXM + vwZE * ZEXM
    vwY = vwXI * XIYM + vwET * ETYM + vwZE * ZEYM
    vwZ = vwXI * XIZM + vwET * ETZM + vwZE * ZEZM
    ! wu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wuXI = 0.25 * ( - wu(I-1,J,K) - wu(I-1,J,K+1) &
    &               + wu(I+1,J,K) + wu(I+1,J,K+1) )
    wuET = 0.25 * ( - wu(I,J-1,K) - wu(I,J-1,K+1) &
    &               + wu(I,J+1,K) + wu(I,J+1,K+1) )
    wuZE = - wu(I,J,K) + wu(I,J,K+1)
    ! 物理空間方向一階微分
    wuX = wuXI * XIXM + wuET * ETXM + wuZE * ZEXM
    wuY = wuXI * XIYM + wuET * ETYM + wuZE * ZEYM
    wuZ = wuXI * XIZM + wuET * ETZM + wuZE * ZEZM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPS(I-1,J,K) - EPS(I-1,J,K+1) &
    &               + EPS(I+1,J,K) + EPS(I+1,J,K+1) )
    AEET = 0.25 * ( - EPS(I,J-1,K) - EPS(I,J-1,K+1) &
    &               + EPS(I,J+1,K) + EPS(I,J+1,K+1) )
    AEZE = - EPS(I,J,K) + EPS(I,J,K+1)
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY + WZ) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 = 2.0 * WZ - SS
    S12 = UY + VX
    S23 = VZ + WY
    S31 = WX + UZ
    ! U の拡散 ---------------------------------------------------------
    R02M = nu * S11
    S02M = nu * S12
    T02M = nu * S31
    R02T =-uuM
    S02T =-uvM
    T02T =-wuM
    ! V の拡散 ---------------------------------------------------------
    R03M = nu * S12
    S03M = nu * S22
    T03M = nu * S23
    R03T =-uvM
    S03T =-vvM
    T03T =-vwM
    ! W の拡散 ---------------------------------------------------------
    R04M = nu * S31
    S04M = nu * S23
    T04M = nu * S33
    R04T =-wuM
    S04T =-vwM
    T04T =-wwM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R05M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM + S31 * WM)
    S05M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM + S23 * WM)
    T05M = C_p * nu / PR * TZ + nu * (S31 * UM + S23 * VM + S33 * WM)
    R05T = AKAE * C_se * C_p * (uuM * TX + uvM * TY + wuM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uuM * (uuX + vvX + wwX) &
    &    + uvM * (uuY + vvY + wwY) &
    &    + wuM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uuM * UM + uvM * VM + wuM * WM)
    S05T = AKAE * C_se * C_p * (uvM * TX + vvM * TY + vwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uvM * (uuX + vvX + wwX) &
    &    + vvM * (uuY + vvY + wwY) &
    &    + vwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uvM * UM + vvM * VM + vwM * WM)
    T05T = AKAE * C_se * C_p * (wuM * TX + vwM * TY + wwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      wuM * (uuX + vvX + wwX) &
    &    + vwM * (uuY + vvY + wwY) &
    &    + wwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (wuM * UM + vwM * VM + wwM * WM)
    ! uu の拡散 --------------------------------------------------------
    R06M = nu * uuX
    S06M = nu * uuY
    T06M = nu * uuZ
    R06T = C_s * AKAE * (uuM * uuX + uvM * uuY + wuM * uuZ)
    S06T = C_s * AKAE * (uvM * uuX + vvM * uuY + vwM * uuZ)
    T06T = C_s * AKAE * (wuM * uuX + vwM * uuY + wwM * uuZ)
    ! vv の拡散 --------------------------------------------------------
    R07M = nu * vvX
    S07M = nu * vvY
    T07M = nu * vvZ
    R07T = C_s * AKAE * (uuM * vvX + uvM * vvY + wuM * vvZ)
    S07T = C_s * AKAE * (uvM * vvX + vvM * vvY + vwM * vvZ)
    T07T = C_s * AKAE * (wuM * vvX + vwM * vvY + wwM * vvZ)
    ! ww の拡散 --------------------------------------------------------
    R08M = nu * wwX
    S08M = nu * wwY
    T08M = nu * wwZ
    R08T = C_s * AKAE * (uuM * wwX + uvM * wwY + wuM * wwZ)
    S08T = C_s * AKAE * (uvM * wwX + vvM * wwY + vwM * wwZ)
    T08T = C_s * AKAE * (wuM * wwX + vwM * wwY + wwM * wwZ)
    ! uv の拡散 --------------------------------------------------------
    R09M = nu * uvX
    S09M = nu * uvY
    T09M = nu * uvZ
    R09T = C_s * AKAE * (uuM * uvX + uvM * uvY + wuM * uvZ)
    S09T = C_s * AKAE * (uvM * uvX + vvM * uvY + vwM * uvZ)
    T09T = C_s * AKAE * (wuM * uvX + vwM * uvY + wwM * uvZ)
    ! vw の拡散 --------------------------------------------------------
    R10M = nu * vwX
    S10M = nu * vwY
    T10M = nu * vwZ
    R10T = C_s * AKAE * (uuM * vwX + uvM * vwY + wuM * vwZ)
    S10T = C_s * AKAE * (uvM * vwX + vvM * vwY + vwM * vwZ)
    T10T = C_s * AKAE * (wuM * vwX + vwM * vwY + wwM * vwZ)
    ! wu の拡散 --------------------------------------------------------
    R11M = nu * wuX
    S11M = nu * wuY
    T11M = nu * wuZ
    R11T = C_s * AKAE * (uuM * wuX + uvM * wuY + wuM * wuZ)
    S11T = C_s * AKAE * (uvM * wuX + vvM * wuY + vwM * wuZ)
    T11T = C_s * AKAE * (wuM * wuX + vwM * wuY + wwM * wuZ)
    ! epsilon の拡散 ---------------------------------------------------
    R12M = nu * AEX
    S12M = nu * AEY
    T12M = nu * AEZ
    R12T = C_e * AKAE * (uuM * AEX + uvM * AEY + wuM * AEZ)
    S12T = C_e * AKAE * (uvM * AEX + vvM * AEY + vwM * AEZ)
    T12T = C_e * AKAE * (wuM * AEX + vwM * AEY + wwM * AEZ)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    TH(I,J,K, 1) = 0.0
    TH(I,J,K, 2) = QH1 * ( (R02M + R02T) * ZEXM &
    &                    + (S02M + S02T) * ZEYM &
    &                    + (T02M + T02T) * ZEZM )
    TH(I,J,K, 3) = QH1 * ( (R03M + R03T) * ZEXM &
    &                    + (S03M + S03T) * ZEYM &
    &                    + (T03M + T03T) * ZEZM )
    TH(I,J,K, 4) = QH1 * ( (R04M + R04T) * ZEXM &
    &                    + (S04M + S04T) * ZEYM &
    &                    + (T04M + T04T) * ZEZM )
    TH(I,J,K, 5) = QH1 * ( (R05M + R05T) * ZEXM &
    &                    + (S05M + S05T) * ZEYM &
    &                    + (T05M + T05T) * ZEZM )
    TH(I,J,K, 6) = QH1 * ( (R06M + R06T) * ZEXM &
    &                    + (S06M + S06T) * ZEYM &
    &                    + (T06M + T06T) * ZEZM )
    TH(I,J,K, 7) = QH1 * ( (R07M + R07T) * ZEXM &
    &                    + (S07M + S07T) * ZEYM &
    &                    + (T07M + T07T) * ZEZM )
    TH(I,J,K, 8) = QH1 * ( (R08M + R08T) * ZEXM &
    &                    + (S08M + S08T) * ZEYM &
    &                    + (T08M + T08T) * ZEZM )
    TH(I,J,K, 9) = QH1 * ( (R09M + R09T) * ZEXM &
    &                    + (S09M + S09T) * ZEYM &
    &                    + (T09M + T09T) * ZEZM )
    TH(I,J,K,10) = QH1 * ( (R10M + R10T) * ZEXM &
    &                    + (S10M + S10T) * ZEYM &
    &                    + (T10M + T10T) * ZEZM )
    TH(I,J,K,11) = QH1 * ( (R11M + R11T) * ZEXM &
    &                    + (S11M + S11T) * ZEYM &
    &                    + (T11M + T11T) * ZEZM )
    TH(I,J,K,12) = QH1 * ( (R12M + R12T) * ZEXM &
    &                    + (S12M + S12T) * ZEYM &
    &                    + (T12M + T12T) * ZEZM )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFZE
! 定義終了 *************************************************************
END SUBROUTINE Turbulence3DRSMSSG
!***********************************************************************
!**** 乱流モデル : RANS, RSM Launder-Shima Model (1989, 1993)       ****
!****              乱流熱伝達にはVandrommeのモデルを使用            ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence2DRSMLS( &
&            PELIM, RG, GAMMA, PR, &
&            IS, IE, JS, JE, LS, LE, &
&            XIX, XIY, ETX, ETY, AJA, &
&            RHO, U, V, T, uu, vv, ww, uv, EPS, AMU, &
&            n1, n2, d, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: C_se = 0.313
  REAL, PARAMETER :: C_s = 0.22, C_e = 0.18
  REAL, PARAMETER :: C_l = 2.5
  REAL, PARAMETER :: C_e1 = 1.45, C_e2 = 1.9
  REAL, PARAMETER :: aiimin = - 0.66666667, aiimax = 1.3333333
  REAL, PARAMETER :: TwoThird = 0.66666667
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE), XIY(IS:IE, JS:JE), &
  &                       ETX(IS:IE, JS:JE), ETY(IS:IE, JS:JE), &
  &                       AJA(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE), &
  &                       U(IS:IE, JS:JE), &
  &                       V(IS:IE, JS:JE), &
  &                       T(IS:IE, JS:JE), &
  &                       uu(IS:IE, JS:JE), &
  &                       vv(IS:IE, JS:JE), &
  &                       ww(IS:IE, JS:JE), &
  &                       uv(IS:IE, JS:JE), &
  &                       EPS(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: n1(IS:IE, JS:JE), n2(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: d(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: AMUT(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(OUT) :: DQP(IS:IE, JS:JE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :), SH(:, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 9) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, LS: LE), SH(IS: IE, JS: JE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: TMP, TMP1, TMPw, TMPw1, TMPw2, TMPw1km, TMPw2km
  REAL    :: UXI, VXI, UET, VET, UX, VX, UY, VY
  REAL    :: KRXI, KRET, KRX, KRY
  REAL    :: P11, P22, P33, P12, P_k
  REAL    :: eps11, eps22, eps33, eps12, e_k
  REAL    :: C_1, C_2, C_w1, C_w2
  REAL    :: psi1, psi2
  REAL    :: A, A2, A3
  REAL    :: a11, a22, a33, a12
  REAL    :: phi1_11, phi1_22, phi1_33, phi1_12
  REAL    :: phi2_11, phi2_22, phi2_33, phi2_12
  REAL    :: phiw1_11, phiw1_22, phiw1_33, phiw1_12
  REAL    :: phiw2_11, phiw2_22, phiw2_33, phiw2_12
  REAL    :: phi11, phi22, phi33, phi12
  REAL    :: AK, EPST, SS2, AEAK, nu, Rt
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  TMP, TMP1, TMPw, TMPw1, TMPw2, TMPw1km, TMPw2km, &
  !$  UXI, VXI, UET, VET, UX, VX, UY, VY, KRXI, KRET, KRX, KRY, &
  !$  P11, P22, P33, P12, P_k, eps11, eps22, eps33, eps12, e_k, &
  !$  C_1, C_2, C_w1, C_w2, psi1, psi2, A, A2, A3, &
  !$  a11, a22, a33, a12, &
  !$  phi1_11, phi1_22, phi1_33, phi1_12, &
  !$  phi2_11, phi2_22, phi2_33, phi2_12, &
  !$  phiw1_11, phiw1_22, phiw1_33, phiw1_12, &
  !$  phiw2_11, phiw2_22, phiw2_33, phiw2_12, &
  !$  phi11, phi22, phi33, phi12, &
  !$  AK, EPST, SS2, AEAK, nu, Rt, QH1 &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF(RHO(I,J) .GT. 0.0) THEN
    ! 乱流エネルギー ---------------------------------------------------
    AK = 0.5 * (uu(I,J) + vv(I,J) + ww(I,J))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.5 * (U(I+1,J) - U(I-1,J))
    VXI = 0.5 * (V(I+1,J) - V(I-1,J))
    UET = 0.5 * (U(I,J+1) - U(I,J-1))
    VET = 0.5 * (V(I,J+1) - V(I,J-1))
    ! 物理空間方向一階微分
    UX = UXI * XIX(I,J) + UET * ETX(I,J)
    UY = UXI * XIY(I,J) + UET * ETY(I,J)
    VX = VXI * XIX(I,J) + VET * ETX(I,J)
    VY = VXI * XIY(I,J) + VET * ETY(I,J)
    ! 生成項と散逸項 ---------------------------------------------------
    P11 = - 2.0 * (uu(I,J) * UX + uv(I,J) * UY)
    P22 = - 2.0 * (uv(I,J) * VX + vv(I,J) * VY)
    P33 = 0.0
    P12 = - (uu(I,J) * VX + uv(I,J) * (UX + VY) + vv(I,J) * UY)
    P_k = 0.5 * (P11 + P22 + P33)
    e_k = EPS(I,J)
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      TMP = TwoThird * e_k / PELIM
      P11 = MIN(TMP, MAX(-TMP, P11))
      P22 = MIN(TMP, MAX(-TMP, P22))
      P33 = MIN(TMP, MAX(-TMP, P33))
      P_k = 0.5 * (P11 + P22 + P33)
    ENDIF
    TMP   = TwoThird * e_k
    eps11 = TMP
    eps22 = TMP
    eps33 = TMP
    eps12 = 0.0
    ! 非等方レイノルズ応力テンソル a_ij --------------------------------
    IF(AK .GT. ZERO) THEN
      a11 = uu(I,J) / AK - TwoThird
      a22 = vv(I,J) / AK - TwoThird
      a33 = ww(I,J) / AK - TwoThird
      a12 = uv(I,J) / AK
      a11 = MIN(aiimax, MAX(aiimin, a11))
      a22 = MIN(aiimax, MAX(aiimin, a22))
      a33 = MIN(aiimax, MAX(aiimin, a33))
      a12 = MIN(1.0, MAX(-1.0, a12))
    ELSE
      a11 = aiimin
      a22 = aiimin
      a33 = aiimin
      a12 = 0.0
    ENDIF
    ! A : Lumley's stress flatness factor ------------------------------
    A2 = a11 * a11 + a22 * a22 + a33 * a33 + 2.0 * a12 * a12
    A2 = MIN(8.0 / 3.0, MAX(0.0, A2))
    A3 = a11 * a11 * a11 + a22 * a22 * a22 + a33 * a33 * a33 &
    &  + 3.0 * a12 * a12 * (a11 + a22)
    A3 = MIN(A2, MAX(A2 - 8.0 / 9.0, A3))
    A  = 1.0 - 9.0 / 8.0 * (A2 - A3)
    A  = MIN(1.0, MAX(0.0, A))
    ! 乱れの散逸率など -------------------------------------------------
    KRXI   = 0.5 * ( &
    &        - SQRT(0.5 * (uu(I-1,J) + vv(I-1,J) + ww(I-1,J))) &
    &        + SQRT(0.5 * (uu(I+1,J) + vv(I+1,J) + ww(I+1,J))) &
    &      )
    KRET   = 0.5 * ( &
    &        - SQRT(0.5 * (uu(I,J-1) + vv(I,J-1) + ww(I,J-1))) &
    &        + SQRT(0.5 * (uu(I,J+1) + vv(I,J+1) + ww(I,J+1))) &
    &      )
    KRX    = KRXI * XIX(I,J) + KRET * ETX(I,J)
    KRY    = KRXI * XIY(I,J) + KRET * ETY(I,J)
    nu     = AMU(I,J) / RHO(I,J)
    EPST   = e_k - 2.0 * nu * (KRX**2 + KRY**2)
    EPST   = MAX(0.0, EPST)
    IF(e_k .GT. ZERO) THEN
      Rt = AK**2 / (nu * e_k)
    ELSE
      Rt = 0.0
    ENDIF
    ! モデル関数 -------------------------------------------------------
    C_1  = 1.0 + 2.58 * A * A2**0.25 * (1.0 - EXP(-(0.0067 * Rt)**2))
    C_2  = 0.75 * SQRT(A)
    C_w1 =-TwoThird * C_1 + 1.67
    IF(C_2 .GT. 0.25) THEN
      C_w2 = (TwoThird * (C_2 - 1.0) + 0.5) / C_2
    ELSE
      C_w2 = 0.0
    ENDIF
    IF(e_k .GT. ZERO) THEN
      psi1 = 1.5 * A * (P_k / e_k - 1.0)
    ELSE
      psi1 = 0.0
    ENDIF
    psi2 = 0.35 * (1.0 - 0.3 * A2) * EXP(-SQRT(0.002 * Rt))
    SS2  = 2.0 * (UX**2 + VY**2) + (UY + VX)**2 &
    &    - TwoThird * (UX + VY)**2
    AEAK = MIN(MAX(SQRT(SS2), SQRT(e_k / nu)), e_k / MAX(ZERO, AK))
    ! 再分配項 ---------------------------------------------------------
    TMP1    = - C_1 * AEAK
    phi1_11 = TMP1 * (uu(I,J) - TwoThird * AK)
    phi1_22 = TMP1 * (vv(I,J) - TwoThird * AK)
    phi1_33 = TMP1 * (ww(I,J) - TwoThird * AK)
    phi1_12 = TMP1 * uv(I,J)
    phi2_11 = - C_2 * (P11 - TwoThird * P_k)
    phi2_22 = - C_2 * (P22 - TwoThird * P_k)
    phi2_33 = - C_2 * (P33 - TwoThird * P_k)
    phi2_12 = - C_2 * P12
    IF(C_l * e_k * d(I,J) .GT. ZERO) THEN
      TMPw = AK**1.5 / (C_l * e_k * d(I,J))
      ! 壁面反射項のダンピング関数の制限
      ! TMPw = MIN(1.0, MAX(0.0, TMPw))
    ELSE
      TMPw = 0.0
    ENDIF
    TMPw1   = C_w1 * TMPw * AEAK
    TMPw2   = C_w2 * TMPw
    TMPw1km = uu(I,J) * n1(I,J) * n1(I,J) &
    &       + vv(I,J) * n2(I,J) * n2(I,J) &
    &       + 2.0 * uv(I,J) * n1(I,J) * n2(I,J)
    TMPw2km = phi2_11 * n1(I,J) * n1(I,J) &
    &       + phi2_22 * n2(I,J) * n2(I,J) &
    &       + 2.0 * phi2_12 * n1(I,J) * n2(I,J)
    phiw1_11 = TMPw1 * (TMPw1km - 3.0 * n1(I,J) &
    &          * (uu(I,J) * n1(I,J) + uv(I,J) * n2(I,J)) &
    &        )
    phiw1_22 = TMPw1 * (TMPw1km - 3.0 * n2(I,J) &
    &          * (uv(I,J) * n1(I,J) + vv(I,J) * n2(I,J)) &
    &        )
    phiw1_33 = TMPw1 * TMPw1km
    phiw1_12 = - 1.5 * TMPw1 * ( &
    &          n1(I,J) * (uv(I,J) * n1(I,J) + vv(I,J) * n2(I,J)) &
    &        + n2(I,J) * (uu(I,J) * n1(I,J) + uv(I,J) * n2(I,J)) &
    &        )
    phiw2_11 = TMPw2 * (TMPw2km - 3.0 * n1(I,J) &
    &          * (phi2_11 * n1(I,J) + phi2_12 * n2(I,J)) &
    &        )
    phiw2_22 = TMPw2 * (TMPw2km - 3.0 * n2(I,J) &
    &          * (phi2_12 * n1(I,J) + phi2_22 * n2(I,J)) &
    &        )
    phiw2_33 = TMPw2 * TMPw2km
    phiw2_12 = - 1.5 * TMPw2 * ( &
    &          n1(I,J) * (phi2_12 * n1(I,J) + phi2_22 * n2(I,J)) &
    &        + n2(I,J) * (phi2_11 * n1(I,J) + phi2_12 * n2(I,J)) &
    &        )
    phi11 = phi1_11 + phi2_11 + phiw1_11 + phiw2_11
    phi22 = phi1_22 + phi2_22 + phiw1_22 + phiw2_22
    phi33 = phi1_33 + phi2_33 + phiw1_33 + phiw2_33
    phi12 = phi1_12 + phi2_12 + phiw1_12 + phiw2_12
    ! 渦粘性係数の計算 -------------------------------------------------
    IF(e_k .GT. ZERO) THEN
      AMUT(I,J) = 0.09 * AK**2 / e_k * RHO(I,J)
    ENDIF
    ! 生成項、圧力再分配項、散逸項の和 ---------------------------------
    QH1 = RHO(I,J) / AJA(I,J)
    DQP(I,J,5) = QH1 * (P11 - eps11 + phi11)
    DQP(I,J,6) = QH1 * (P22 - eps22 + phi22)
    DQP(I,J,7) = QH1 * (P33 - eps33 + phi33)
    DQP(I,J,8) = QH1 * (P12 - eps12 + phi12)
    DQP(I,J,9) = QH1 * AEAK * ( &
    &            (C_e1 + psi1 + psi2) * P_k - C_e2 * EPST &
    &          )
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH = 0.0
  SH = 0.0
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL DIFFXI
  CALL DIFFET
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, L)
  DO L = LS, LE
  !$OMP DO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J-1) .GT. 0.0) .AND. (RHO(I-1,J) .GT. 0.0) .AND. &
  &   (RHO(I,J)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,L) = (-RH(I-1,J  ,L) + RH(I,J,L)) &
    &          + (-SH(I  ,J-1,L) + SH(I,J,L))
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項の計算(xi方向差分用)                                    ****
!***********************************************************************
SUBROUTINE DIFFXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: RHOM, UM, VM, uuM, vvM, wwM, uvM, EPSM, AMUM
  REAL    :: UXI, UET, UX, UY
  REAL    :: VXI, VET, VX, VY
  REAL    :: TXI, TET, TX, TY
  REAL    :: uuXI, uuET, uuX, uuY
  REAL    :: vvXI, vvET, vvX, vvY
  REAL    :: wwXI, wwET, wwX, wwY
  REAL    :: uvXI, uvET, uvX, uvY
  REAL    :: AEXI, AEET, AEX, AEY
  REAL    :: nu
  REAL    :: S11, S22, S12, SS
  REAL    :: C_p, AKAE
  REAL    :: R2M, S2M, R2T, S2T
  REAL    :: R3M, S3M, R3T, S3T
  REAL    :: R4M, S4M, R4T, S4T
  REAL    :: R5M, S5M, R5T, S5T
  REAL    :: R6M, S6M, R6T, S6T
  REAL    :: R7M, S7M, R7T, S7T
  REAL    :: R8M, S8M, R8T, S8T
  REAL    :: R9M, S9M, R9T, S9T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, &
  !$  RHOM, UM, VM, uuM, vvM, wwM, uvM, EPSM, AMUM, &
  !$  UXI, UET, UX, UY, VXI, VET, VX, VY, TXI, TET, TX, TY, &
  !$  uuXI, uuET, uuX, uuY, vvXI, vvET, vvX, vvY, &
  !$  wwXI, wwET, wwX, wwY, uvXI, uvET, uvX, uvY, &
  !$  AEXI, AEET, AEX, AEY, &
  !$  nu, S11, S22, S12, SS, C_p, AKAE, &
  !$  R2M, S2M, R2T, S2T, R3M, S3M, R3T, S3T, R4M, S4M, R4T, S4T, &
  !$  R5M, S5M, R5T, S5T, R6M, S6M, R6T, S6T, R7M, S7M, R7T, S7T, &
  !$  R8M, S8M, R8T, S8T, R9M, S9M, R9T, S9T, QH1 &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF((RHO(I,J) .GT. 0.0) .AND. (RHO(I+1,J) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J) + XIX(I+1,J))
    XIYM = 0.5 * (XIY(I,J) + XIY(I+1,J))
    ETXM = 0.5 * (ETX(I,J) + ETX(I+1,J))
    ETYM = 0.5 * (ETY(I,J) + ETY(I+1,J))
    AJAM = 0.5 * (AJA(I,J) + AJA(I+1,J))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J) + RHO(I+1,J))
    UM   = 0.5 * (  U(I,J) +   U(I+1,J))
    VM   = 0.5 * (  V(I,J) +   V(I+1,J))
    uuM  = 0.5 * ( uu(I,J) +  uu(I+1,J))
    vvM  = 0.5 * ( vv(I,J) +  vv(I+1,J))
    wwM  = 0.5 * ( ww(I,J) +  ww(I+1,J))
    uvM  = 0.5 * ( uv(I,J) +  uv(I+1,J))
    EPSM = 0.5 * (EPS(I,J) + EPS(I+1,J))
    AMUM = 0.5 * (AMU(I,J) + AMU(I+1,J))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = - U(I,J) + U(I+1,J)
    UET = 0.25 * (- U(I,J-1) - U(I+1,J-1) + U(I,J+1) + U(I+1,J+1))
    VXI = - V(I,J) + V(I+1,J)
    VET = 0.25 * (- V(I,J-1) - V(I+1,J-1) + V(I,J+1) + V(I+1,J+1))
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM
    UY = UXI * XIYM + UET * ETYM
    VX = VXI * XIXM + VET * ETXM
    VY = VXI * XIYM + VET * ETYM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = - T(I,J) + T(I+1,J)
    TET = 0.25 * (- T(I,J-1) - T(I+1,J-1) + T(I,J+1) + T(I+1,J+1))
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM
    TY = TXI * XIYM + TET * ETYM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = - uu(I,J) + uu(I+1,J)
    uuET = 0.25 * ( - uu(I,J-1) - uu(I+1,J-1) &
    &               + uu(I,J+1) + uu(I+1,J+1) )
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM
    uuY = uuXI * XIYM + uuET * ETYM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = - vv(I,J) + vv(I+1,J)
    vvET = 0.25 * ( - vv(I,J-1) - vv(I+1,J-1) &
    &               + vv(I,J+1) + vv(I+1,J+1) )
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM
    vvY = vvXI * XIYM + vvET * ETYM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = - ww(I,J) + ww(I+1,J)
    wwET = 0.25 * ( - ww(I,J-1) - ww(I+1,J-1) &
    &               + ww(I,J+1) + ww(I+1,J+1) )
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM
    wwY = wwXI * XIYM + wwET * ETYM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = - uv(I,J) + uv(I+1,J)
    uvET = 0.25 * ( - uv(I,J-1) - uv(I+1,J-1) &
    &               + uv(I,J+1) + uv(I+1,J+1) )
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM
    uvY = uvXI * XIYM + uvET * ETYM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = - EPS(I,J) + EPS(I+1,J)
    AEET = 0.25 * ( - EPS(I,J-1) - EPS(I+1,J-1) &
    &               + EPS(I,J+1) + EPS(I+1,J+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM
    AEY = AEXI * XIYM + AEET * ETYM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S12 = UY + VX
    ! U の拡散 ---------------------------------------------------------
    R2M = nu * S11
    S2M = nu * S12
    R2T =-uuM
    S2T =-uvM
    ! V の拡散 ---------------------------------------------------------
    R3M = nu * S12
    S3M = nu * S22
    R3T =-uvM
    S3T =-vvM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R4M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM)
    S4M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM)
    R4T = AKAE * C_se * C_p * (uuM * TX + uvM * TY) &
    &   + AKAE * C_s  * ( &
    &     uuM * (uuX + vvX + wwX) + uvM * (uuY + vvY + wwY) &
    &   ) &
    &   - (uuM * UM + uvM * VM)
    S4T = AKAE * C_se * C_p * (uvM * TX + vvM * TY) &
    &   + AKAE * C_s  * ( &
    &     uvM * (uuX + vvX + wwX) + vvM * (uuY + vvY + wwY) &
    &   ) &
    &   - (uvM * UM + vvM * VM)
    ! uu の拡散 --------------------------------------------------------
    R5M = nu * uuX
    S5M = nu * uuY
    R5T = C_s * AKAE * (uuM * uuX + uvM * uuY)
    S5T = C_s * AKAE * (uvM * uuX + vvM * uuY)
    ! vv の拡散 --------------------------------------------------------
    R6M = nu * vvX
    S6M = nu * vvY
    R6T = C_s * AKAE * (uuM * vvX + uvM * vvY)
    S6T = C_s * AKAE * (uvM * vvX + vvM * vvY)
    ! ww の拡散 --------------------------------------------------------
    R7M = nu * wwX
    S7M = nu * wwY
    R7T = C_s * AKAE * (uuM * wwX + uvM * wwY)
    S7T = C_s * AKAE * (uvM * wwX + vvM * wwY)
    ! uv の拡散 --------------------------------------------------------
    R8M = nu * uvX
    S8M = nu * uvY
    R8T = C_s * AKAE * (uuM * uvX + uvM * uvY)
    S8T = C_s * AKAE * (uvM * uvX + vvM * uvY)
    ! epsilon の拡散 ---------------------------------------------------
    R9M = nu * AEX
    S9M = nu * AEY
    R9T = C_e * AKAE * (uuM * AEX + uvM * AEY)
    S9T = C_e * AKAE * (uvM * AEX + vvM * AEY)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    RH(I,J,1) = 0.0
    RH(I,J,2) = QH1 * ((R2M + R2T) * XIXM + (S2M + S2T) * XIYM)
    RH(I,J,3) = QH1 * ((R3M + R3T) * XIXM + (S3M + S3T) * XIYM)
    RH(I,J,4) = QH1 * ((R4M + R4T) * XIXM + (S4M + S4T) * XIYM)
    RH(I,J,5) = QH1 * ((R5M + R5T) * XIXM + (S5M + S5T) * XIYM)
    RH(I,J,6) = QH1 * ((R6M + R6T) * XIXM + (S6M + S6T) * XIYM)
    RH(I,J,7) = QH1 * ((R7M + R7T) * XIXM + (S7M + S7T) * XIYM)
    RH(I,J,8) = QH1 * ((R8M + R8T) * XIXM + (S8M + S8T) * XIYM)
    RH(I,J,9) = QH1 * ((R9M + R9T) * XIXM + (S9M + S9T) * XIYM)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFXI
!***********************************************************************
!**** 拡散項の計算(eta方向差分用)                                   ****
!***********************************************************************
SUBROUTINE DIFFET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: RHOM, UM, VM, uuM, vvM, wwM, uvM, EPSM, AMUM
  REAL    :: UXI, UET, UX, UY
  REAL    :: VXI, VET, VX, VY
  REAL    :: TXI, TET, TX, TY
  REAL    :: uuXI, uuET, uuX, uuY
  REAL    :: vvXI, vvET, vvX, vvY
  REAL    :: wwXI, wwET, wwX, wwY
  REAL    :: uvXI, uvET, uvX, uvY
  REAL    :: AEXI, AEET, AEX, AEY
  REAL    :: nu
  REAL    :: S11, S22, S12, SS
  REAL    :: C_p, AKAE
  REAL    :: R2M, S2M, R2T, S2T
  REAL    :: R3M, S3M, R3T, S3T
  REAL    :: R4M, S4M, R4T, S4T
  REAL    :: R5M, S5M, R5T, S5T
  REAL    :: R6M, S6M, R6T, S6T
  REAL    :: R7M, S7M, R7T, S7T
  REAL    :: R8M, S8M, R8T, S8T
  REAL    :: R9M, S9M, R9T, S9T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, &
  !$  RHOM, UM, VM, uuM, vvM, wwM, uvM, EPSM, AMUM, &
  !$  UXI, UET, UX, UY, VXI, VET, VX, VY, TXI, TET, TX, TY, &
  !$  uuXI, uuET, uuX, uuY, vvXI, vvET, vvX, vvY, &
  !$  wwXI, wwET, wwX, wwY, uvXI, uvET, uvX, uvY, &
  !$  AEXI, AEET, AEX, AEY, &
  !$  nu, S11, S22, S12, SS, C_p, AKAE, &
  !$  R2M, S2M, R2T, S2T, R3M, S3M, R3T, S3T, R4M, S4M, R4T, S4T, &
  !$  R5M, S5M, R5T, S5T, R6M, S6M, R6T, S6T, R7M, S7M, R7T, S7T, &
  !$  R8M, S8M, R8T, S8T, R9M, S9M, R9T, S9T, QH1 &
  !$)
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF((RHO(I,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0)) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J) + XIX(I,J+1))
    XIYM = 0.5 * (XIY(I,J) + XIY(I,J+1))
    ETXM = 0.5 * (ETX(I,J) + ETX(I,J+1))
    ETYM = 0.5 * (ETY(I,J) + ETY(I,J+1))
    AJAM = 0.5 * (AJA(I,J) + AJA(I,J+1))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J) + RHO(I,J+1))
    UM   = 0.5 * (  U(I,J) +   U(I,J+1))
    VM   = 0.5 * (  V(I,J) +   V(I,J+1))
    uuM  = 0.5 * ( uu(I,J) +  uu(I,J+1))
    vvM  = 0.5 * ( vv(I,J) +  vv(I,J+1))
    wwM  = 0.5 * ( ww(I,J) +  ww(I,J+1))
    uvM  = 0.5 * ( uv(I,J) +  uv(I,J+1))
    EPSM = 0.5 * (EPS(I,J) + EPS(I,J+1))
    AMUM = 0.5 * (AMU(I,J) + AMU(I,J+1))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * (- U(I-1,J) - U(I-1,J+1) + U(I+1,J) + U(I+1,J+1))
    UET = - U(I,J) + U(I,J+1)
    VXI = 0.25 * (- V(I-1,J) - V(I-1,J+1) + V(I+1,J) + V(I+1,J+1))
    VET = - V(I,J) + V(I,J+1)
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM
    UY = UXI * XIYM + UET * ETYM
    VX = VXI * XIXM + VET * ETXM
    VY = VXI * XIYM + VET * ETYM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * (- T(I-1,J) - T(I-1,J+1) + T(I+1,J) + T(I+1,J+1))
    TET = - T(I,J) + T(I,J+1)
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM
    TY = TXI * XIYM + TET * ETYM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = 0.25 * (- uu(I-1,J) - uu(I-1,J+1) + uu(I+1,J) + uu(I+1,J+1))
    uuET = - uu(I,J) + uu(I,J+1)
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM
    uuY = uuXI * XIYM + uuET * ETYM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = 0.25 * (- vv(I-1,J) - vv(I-1,J+1) + vv(I+1,J) + vv(I+1,J+1))
    vvET = - vv(I,J) + vv(I,J+1)
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM
    vvY = vvXI * XIYM + vvET * ETYM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = 0.25 * (- ww(I-1,J) - ww(I-1,J+1) + ww(I+1,J) + ww(I+1,J+1))
    wwET = - ww(I,J) + ww(I,J+1)
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM
    wwY = wwXI * XIYM + wwET * ETYM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = 0.25 * (- uv(I-1,J) - uv(I-1,J+1) + uv(I+1,J) + uv(I+1,J+1))
    uvET = - uv(I,J) + uv(I,J+1)
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM
    uvY = uvXI * XIYM + uvET * ETYM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPS(I-1,J) - EPS(I-1,J+1) &
    &               + EPS(I+1,J) + EPS(I+1,J+1) )
    AEET = - EPS(I,J) + EPS(I,J+1)
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM
    AEY = AEXI * XIYM + AEET * ETYM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S12 = UY + VX
    ! U の拡散 ---------------------------------------------------------
    R2M = nu * S11
    S2M = nu * S12
    R2T =-uuM
    S2T =-uvM
    ! V の拡散 ---------------------------------------------------------
    R3M = nu * S12
    S3M = nu * S22
    R3T =-uvM
    S3T =-vvM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R4M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM)
    S4M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM)
    R4T = AKAE * C_se * C_p * (uuM * TX + uvM * TY) &
    &   + AKAE * C_s  * ( &
    &     uuM * (uuX + vvX + wwX) + uvM * (uuY + vvY + wwY) &
    &   ) &
    &   - (uuM * UM + uvM * VM)
    S4T = AKAE * C_se * C_p * (uvM * TX + vvM * TY) &
    &   + AKAE * C_s  * ( &
    &     uvM * (uuX + vvX + wwX) + vvM * (uuY + vvY + wwY) &
    &   ) &
    &   - (uvM * UM + vvM * VM)
    ! uu の拡散 --------------------------------------------------------
    R5M = nu * uuX
    S5M = nu * uuY
    R5T = C_s * AKAE * (uuM * uuX + uvM * uuY)
    S5T = C_s * AKAE * (uvM * uuX + vvM * uuY)
    ! vv の拡散 --------------------------------------------------------
    R6M = nu * vvX
    S6M = nu * vvY
    R6T = C_s * AKAE * (uuM * vvX + uvM * vvY)
    S6T = C_s * AKAE * (uvM * vvX + vvM * vvY)
    ! ww の拡散 --------------------------------------------------------
    R7M = nu * wwX
    S7M = nu * wwY
    R7T = C_s * AKAE * (uuM * wwX + uvM * wwY)
    S7T = C_s * AKAE * (uvM * wwX + vvM * wwY)
    ! uv の拡散 --------------------------------------------------------
    R8M = nu * uvX
    S8M = nu * uvY
    R8T = C_s * AKAE * (uuM * uvX + uvM * uvY)
    S8T = C_s * AKAE * (uvM * uvX + vvM * uvY)
    ! epsilon の拡散 ---------------------------------------------------
    R9M = nu * AEX
    S9M = nu * AEY
    R9T = C_e * AKAE * (uuM * AEX + uvM * AEY)
    S9T = C_e * AKAE * (uvM * AEX + vvM * AEY)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    SH(I,J,1) = 0.0
    SH(I,J,2) = QH1 * ((R2M + R2T) * ETXM + (S2M + S2T) * ETYM)
    SH(I,J,3) = QH1 * ((R3M + R3T) * ETXM + (S3M + S3T) * ETYM)
    SH(I,J,4) = QH1 * ((R4M + R4T) * ETXM + (S4M + S4T) * ETYM)
    SH(I,J,5) = QH1 * ((R5M + R5T) * ETXM + (S5M + S5T) * ETYM)
    SH(I,J,6) = QH1 * ((R6M + R6T) * ETXM + (S6M + S6T) * ETYM)
    SH(I,J,7) = QH1 * ((R7M + R7T) * ETXM + (S7M + S7T) * ETYM)
    SH(I,J,8) = QH1 * ((R8M + R8T) * ETXM + (S8M + S8T) * ETYM)
    SH(I,J,9) = QH1 * ((R9M + R9T) * ETXM + (S9M + S9T) * ETYM)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFET
! 定義終了 *************************************************************
END SUBROUTINE Turbulence2DRSMLS
!***********************************************************************
!**** 乱流モデル : RANS, RSM Launder-Shima Model (1989, 1993)       ****
!****              乱流熱伝達にはVandrommeのモデルを使用            ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence3DRSMLS( &
&            PELIM, RG, GAMMA, PR, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, AJA, &
&            RHO, U, V, W, T, uu, vv, ww, uv, vw, wu, EPS, AMU, &
&            n1, n2, n3, d, &
&            AMUT, DQD, DQP &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: C_se = 0.313
  REAL, PARAMETER :: C_s = 0.22, C_e = 0.18
  REAL, PARAMETER :: C_l = 2.5
  REAL, PARAMETER :: C_e1 = 1.45, C_e2 = 1.9
  REAL, PARAMETER :: aiimin = - 0.66666667, aiimax = 1.3333333
  REAL, PARAMETER :: TwoThird = 0.66666667
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: PELIM
  REAL,    INTENT(IN)  :: RG, GAMMA, PR
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE, KS:KE), &
  &                       XIY(IS:IE, JS:JE, KS:KE), &
  &                       XIZ(IS:IE, JS:JE, KS:KE), &
  &                       ETX(IS:IE, JS:JE, KS:KE), &
  &                       ETY(IS:IE, JS:JE, KS:KE), &
  &                       ETZ(IS:IE, JS:JE, KS:KE), &
  &                       ZEX(IS:IE, JS:JE, KS:KE), &
  &                       ZEY(IS:IE, JS:JE, KS:KE), &
  &                       ZEZ(IS:IE, JS:JE, KS:KE), &
  &                       AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE, KS:KE), &
  &                       U(IS:IE, JS:JE, KS:KE), &
  &                       V(IS:IE, JS:JE, KS:KE), &
  &                       W(IS:IE, JS:JE, KS:KE), &
  &                       T(IS:IE, JS:JE, KS:KE), &
  &                       uu(IS:IE, JS:JE, KS:KE), &
  &                       vv(IS:IE, JS:JE, KS:KE), &
  &                       ww(IS:IE, JS:JE, KS:KE), &
  &                       uv(IS:IE, JS:JE, KS:KE), &
  &                       vw(IS:IE, JS:JE, KS:KE), &
  &                       wu(IS:IE, JS:JE, KS:KE), &
  &                       EPS(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: n1(IS:IE, JS:JE, KS:KE), &
  &                       n2(IS:IE, JS:JE, KS:KE), &
  &                       n3(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: d(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: AMUT(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(OUT) :: DQP(IS:IE, JS:JE, KS:KE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :, :), SH(:, :, :, :), TH(:, :, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 12) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(SH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(TH(IS: IE, JS: JE, KS: KE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  DQP  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: TMP, TMP1, TMPw, TMPw1, TMPw2, TMPw1km, TMPw2km
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  REAL    :: UX, VX, WX, UY, VY, WY, UZ, VZ, WZ
  REAL    :: KRXI, KRET, KRZE, KRX, KRY, KRZ
  REAL    :: P11, P22, P33, P12, P23, P31, P_k
  REAL    :: eps11, eps22, eps33, eps12, eps23, eps31, e_k
  REAL    :: C_1, C_2, C_w1, C_w2
  REAL    :: psi1, psi2
  REAL    :: A, A2, A3
  REAL    :: a11, a22, a33, a12, a23, a31
  REAL    :: phi1_11, phi1_22, phi1_33, phi1_12, phi1_23, phi1_31
  REAL    :: phi2_11, phi2_22, phi2_33, phi2_12, phi2_23, phi2_31
  REAL    :: phiw1_11, phiw1_22, phiw1_33, phiw1_12, phiw1_23, phiw1_31
  REAL    :: phiw2_11, phiw2_22, phiw2_33, phiw2_12, phiw2_23, phiw2_31
  REAL    :: phi11, phi22, phi33, phi12, phi23, phi31
  REAL    :: AK, EPST, SS2, AEAK, nu, Rt
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  TMP, TMP1, TMPw, TMPw1, TMPw2, TMPw1km, TMPw2km, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE, &
  !$  UX, VX, WX, UY, VY, WY, UZ, VZ, WZ, &
  !$  KRXI, KRET, KRZE, KRX, KRY, KRZ, &
  !$  P11, P22, P33, P12, P23, P31, P_k, &
  !$  eps11, eps22, eps33, eps12, eps23, eps31, e_k, &
  !$  C_1, C_2, C_w1, C_w2, psi1, psi2, A, A2, A3, &
  !$  a11, a22, a33, a12, a23, a31, &
  !$  phi1_11, phi1_22, phi1_33, phi1_12, phi1_23, phi1_31, &
  !$  phi2_11, phi2_22, phi2_33, phi2_12, phi2_23, phi2_31, &
  !$  phiw1_11, phiw1_22, phiw1_33, phiw1_12, phiw1_23, phiw1_31, &
  !$  phiw2_11, phiw2_22, phiw2_33, phiw2_12, phiw2_23, phiw2_31, &
  !$  phi11, phi22, phi33, phi12, phi23, phi31, &
  !$  AK, EPST, SS2, AEAK, nu, Rt, QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF(RHO(I,J,K) .GT. 0.0) THEN
    ! 乱流エネルギー ---------------------------------------------------
    AK = 0.5 * (uu(I,J,K) + vv(I,J,K) + ww(I,J,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.5 * (U(I+1,J,K) - U(I-1,J,K))
    UET = 0.5 * (U(I,J+1,K) - U(I,J-1,K))
    UZE = 0.5 * (U(I,J,K+1) - U(I,J,K-1))
    VXI = 0.5 * (V(I+1,J,K) - V(I-1,J,K))
    VET = 0.5 * (V(I,J+1,K) - V(I,J-1,K))
    VZE = 0.5 * (V(I,J,K+1) - V(I,J,K-1))
    WXI = 0.5 * (W(I+1,J,K) - W(I-1,J,K))
    WET = 0.5 * (W(I,J+1,K) - W(I,J-1,K))
    WZE = 0.5 * (W(I,J,K+1) - W(I,J,K-1))
    ! 物理空間方向一階微分
    UX = UXI * XIX(I,J,K) + UET * ETX(I,J,K) + UZE * ZEX(I,J,K)
    UY = UXI * XIY(I,J,K) + UET * ETY(I,J,K) + UZE * ZEY(I,J,K)
    UZ = UXI * XIZ(I,J,K) + UET * ETZ(I,J,K) + UZE * ZEZ(I,J,K)
    VX = VXI * XIX(I,J,K) + VET * ETX(I,J,K) + VZE * ZEX(I,J,K)
    VY = VXI * XIY(I,J,K) + VET * ETY(I,J,K) + VZE * ZEY(I,J,K)
    VZ = VXI * XIZ(I,J,K) + VET * ETZ(I,J,K) + VZE * ZEZ(I,J,K)
    WX = WXI * XIX(I,J,K) + WET * ETX(I,J,K) + WZE * ZEX(I,J,K)
    WY = WXI * XIY(I,J,K) + WET * ETY(I,J,K) + WZE * ZEY(I,J,K)
    WZ = WXI * XIZ(I,J,K) + WET * ETZ(I,J,K) + WZE * ZEZ(I,J,K)
    ! 生成項と散逸項 ---------------------------------------------------
    P11 = - 2.0 * (uu(I,J,K) * UX + uv(I,J,K) * UY + wu(I,J,K) * UZ)
    P22 = - 2.0 * (uv(I,J,K) * VX + vv(I,J,K) * VY + vw(I,J,K) * VZ)
    P33 = - 2.0 * (wu(I,J,K) * WX + vw(I,J,K) * WY + ww(I,J,K) * WZ)
    P12 = - (uv(I,J,K) * UX + uu(I,J,K) * VX) &
    &     - (vv(I,J,K) * UY + uv(I,J,K) * VY) &
    &     - (vw(I,J,K) * UZ + wu(I,J,K) * VZ)
    P23 = - (wu(I,J,K) * VX + uv(I,J,K) * WX) &
    &     - (vw(I,J,K) * VY + vv(I,J,K) * WY) &
    &     - (ww(I,J,K) * VZ + vw(I,J,K) * WZ)
    P31 = - (uu(I,J,K) * WX + wu(I,J,K) * UX) &
    &     - (uv(I,J,K) * WY + vw(I,J,K) * UY) &
    &     - (wu(I,J,K) * WZ + ww(I,J,K) * UZ)
    P_k = 0.5 * (P11 + P22 + P33)
    e_k = EPS(I,J,K)
    ! 局所平衡仮定(生成と散逸のバランス)によるリミッター
    IF(PELIM .GT. 0.0) THEN
      TMP = TwoThird * e_k / PELIM
      P11 = MIN(TMP, MAX(-TMP, P11))
      P22 = MIN(TMP, MAX(-TMP, P22))
      P33 = MIN(TMP, MAX(-TMP, P33))
      P_k = 0.5 * (P11 + P22 + P33)
    ENDIF
    TMP   = TwoThird * e_k
    eps11 = TMP
    eps22 = TMP
    eps33 = TMP
    eps12 = 0.0
    eps23 = 0.0
    eps31 = 0.0
    ! 非等方レイノルズ応力テンソル a_ij --------------------------------
    IF(AK .GT. ZERO) THEN
      a11 = uu(I,J,K) / AK - TwoThird
      a22 = vv(I,J,K) / AK - TwoThird
      a33 = ww(I,J,K) / AK - TwoThird
      a12 = uv(I,J,K) / AK
      a23 = vw(I,J,K) / AK
      a31 = wu(I,J,K) / AK
      a11 = MIN(aiimax, MAX(aiimin, a11))
      a22 = MIN(aiimax, MAX(aiimin, a22))
      a33 = MIN(aiimax, MAX(aiimin, a33))
      a12 = MIN(1.0, MAX(-1.0, a12))
      a23 = MIN(1.0, MAX(-1.0, a23))
      a31 = MIN(1.0, MAX(-1.0, a31))
    ELSE
      a11 = aiimin
      a22 = aiimin
      a33 = aiimin
      a12 = 0.0
      a23 = 0.0
      a31 = 0.0
    ENDIF
    ! A : Lumley's stress flatness factor ------------------------------
    A2 = a11 * a11 + a22 * a22 + a33 * a33 &
    &  + 2.0 * (a12 * a12 + a23 * a23 + a31 * a31)
    A2 = MIN(8.0 / 3.0, MAX(0.0, A2))
    A3 = a11 * (a11 * a11 + a12 * a12 + a31 * a31) &
    &  + a22 * (a12 * a12 + a22 * a22 + a23 * a23) &
    &  + a33 * (a31 * a31 + a23 * a23 + a33 * a33) &
    &  + 2.0 * a12 * (a11 * a12 + a12 * a22 + a31 * a23) &
    &  + 2.0 * a23 * (a12 * a31 + a22 * a23 + a23 * a33) &
    &  + 2.0 * a31 * (a31 * a11 + a23 * a12 + a33 * a31)
    A3 = MIN(A2, MAX(A2 - 8.0 / 9.0, A3))
    A  = 1.0 - 9.0 / 8.0 * (A2 - A3)
    A  = MIN(1.0, MAX(0.0, A))
    ! 乱れの散逸率など -------------------------------------------------
    KRXI   = 0.5 * ( &
    &        - SQRT(0.5 * (uu(I-1,J,K) + vv(I-1,J,K) + ww(I-1,J,K))) &
    &        + SQRT(0.5 * (uu(I+1,J,K) + vv(I+1,J,K) + ww(I+1,J,K))) &
    &      )
    KRET   = 0.5 * ( &
    &        - SQRT(0.5 * (uu(I,J-1,K) + vv(I,J-1,K) + ww(I,J-1,K))) &
    &        + SQRT(0.5 * (uu(I,J+1,K) + vv(I,J+1,K) + ww(I,J+1,K))) &
    &      )
    KRZE   = 0.5 * ( &
    &        - SQRT(0.5 * (uu(I,J,K-1) + vv(I,J,K-1) + ww(I,J,K-1))) &
    &        + SQRT(0.5 * (uu(I,J,K+1) + vv(I,J,K+1) + ww(I,J,K+1))) &
    &      )
    KRX    = KRXI * XIX(I,J,K) + KRET * ETX(I,J,K) + KRZE * ZEX(I,J,K)
    KRY    = KRXI * XIY(I,J,K) + KRET * ETY(I,J,K) + KRZE * ZEY(I,J,K)
    KRZ    = KRXI * XIZ(I,J,K) + KRET * ETZ(I,J,K) + KRZE * ZEZ(I,J,K)
    nu     = AMU(I,J,K) / RHO(I,J,K)
    EPST   = e_k - 2.0 * nu * (KRX**2 + KRY**2 + KRZ**2)
    EPST   = MAX(0.0, EPST)
    IF(e_k .GT. ZERO) THEN
      Rt = AK**2 / (nu * e_k)
    ELSE
      Rt = 0.0
    ENDIF
    ! モデル関数 -------------------------------------------------------
    C_1  = 1.0 + 2.58 * A * A2**0.25 * (1.0 - EXP(-(0.0067 * Rt)**2))
    C_2  = 0.75 * SQRT(A)
    C_w1 =-TwoThird * C_1 + 1.67
    IF(C_2 .GT. 0.25) THEN
      C_w2 = (TwoThird * (C_2 - 1.0) + 0.5) / C_2
    ELSE
      C_w2 = 0.0
    ENDIF
    IF(e_k .GT. ZERO) THEN
      psi1 = 1.5 * A * (P_k / e_k - 1.0)
    ELSE
      psi1 = 0.0
    ENDIF
    psi2 = 0.35 * (1.0 - 0.3 * A2) * EXP(-SQRT(0.002 * Rt))
    SS2  = 2.0 * (UX**2 + VY**2 + WZ**2) &
    &    + (UY + VX)**2 + (VZ + WY)**2 + (WX + UZ)**2 &
    &    - TwoThird * (UX + VY + WZ)**2
    AEAK = MIN(MAX(SQRT(SS2), SQRT(e_k / nu)), e_k / MAX(ZERO, AK))
    ! 再分配項 ---------------------------------------------------------
    TMP1    = - C_1 * AEAK
    phi1_11 = TMP1 * (uu(I,J,K) - TwoThird * AK)
    phi1_22 = TMP1 * (vv(I,J,K) - TwoThird * AK)
    phi1_33 = TMP1 * (ww(I,J,K) - TwoThird * AK)
    phi1_12 = TMP1 * uv(I,J,K)
    phi1_23 = TMP1 * vw(I,J,K)
    phi1_31 = TMP1 * wu(I,J,K)
    phi2_11 = - C_2 * (P11 - TwoThird * P_k)
    phi2_22 = - C_2 * (P22 - TwoThird * P_k)
    phi2_33 = - C_2 * (P33 - TwoThird * P_k)
    phi2_12 = - C_2 * P12
    phi2_23 = - C_2 * P23
    phi2_31 = - C_2 * P31
    IF(C_l * e_k * d(I,J,K) .GT. ZERO) THEN
      TMPw = AK**1.5 / (C_l * e_k * d(I,J,K))
      ! 壁面反射項のダンピング関数の制限
      ! TMPw = MIN(1.0, MAX(0.0, TMPw))
    ELSE
      TMPw = 0.0
    ENDIF
    TMPw1   = C_w1 * TMPw * AEAK
    TMPw2   = C_w2 * TMPw
    TMPw1km = uu(I,J,K) * n1(I,J,K) * n1(I,J,K) &
    &       + vv(I,J,K) * n2(I,J,K) * n2(I,J,K) &
    &       + ww(I,J,K) * n3(I,J,K) * n3(I,J,K) &
    &       + 2.0 * ( &
    &         uv(I,J,K) * n1(I,J,K) * n2(I,J,K) &
    &       + vw(I,J,K) * n2(I,J,K) * n3(I,J,K) &
    &       + wu(I,J,K) * n3(I,J,K) * n1(I,J,K) &
    &       )
    TMPw2km = phi2_11 * n1(I,J,K) * n1(I,J,K) &
    &       + phi2_22 * n2(I,J,K) * n2(I,J,K) &
    &       + phi2_33 * n3(I,J,K) * n3(I,J,K) &
    &       + 2.0 * ( &
    &         phi2_12 * n1(I,J,K) * n2(I,J,K) &
    &       + phi2_23 * n2(I,J,K) * n3(I,J,K) &
    &       + phi2_31 * n3(I,J,K) * n1(I,J,K) &
    &       )
    phiw1_11 = TMPw1 * (TMPw1km - 3.0 * n1(I,J,K) &
    &          * ( uu(I,J,K) * n1(I,J,K) &
    &            + uv(I,J,K) * n2(I,J,K) &
    &            + wu(I,J,K) * n3(I,J,K) ) &
    &        )
    phiw1_22 = TMPw1 * (TMPw1km - 3.0 * n2(I,J,K) &
    &          * ( uv(I,J,K) * n1(I,J,K) &
    &            + vv(I,J,K) * n2(I,J,K) &
    &            + vw(I,J,K) * n3(I,J,K) ) &
    &        )
    phiw1_33 = TMPw1 * (TMPw1km - 3.0 * n3(I,J,K) &
    &          * ( wu(I,J,K) * n1(I,J,K) &
    &            + vw(I,J,K) * n2(I,J,K) &
    &            + ww(I,J,K) * n3(I,J,K) ) &
    &        )
    phiw1_12 = - 1.5 * TMPw1 * ( &
    &          n1(I,J,K) * ( uv(I,J,K) * n1(I,J,K) &
    &                      + vv(I,J,K) * n2(I,J,K) &
    &                      + vw(I,J,K) * n3(I,J,K) ) &
    &        + n2(I,J,K) * ( uu(I,J,K) * n1(I,J,K) &
    &                      + uv(I,J,K) * n2(I,J,K) &
    &                      + wu(I,J,K) * n3(I,J,K) ) &
    &        )
    phiw1_23 = - 1.5 * TMPw1 * ( &
    &          n2(I,J,K) * ( wu(I,J,K) * n1(I,J,K) &
    &                      + vw(I,J,K) * n2(I,J,K) &
    &                      + ww(I,J,K) * n3(I,J,K) ) &
    &        + n3(I,J,K) * ( uv(I,J,K) * n1(I,J,K) &
    &                      + vv(I,J,K) * n2(I,J,K) &
    &                      + vw(I,J,K) * n3(I,J,K) ) &
    &        )
    phiw1_31 = - 1.5 * TMPw1 * ( &
    &          n3(I,J,K) * ( uu(I,J,K) * n1(I,J,K) &
    &                      + uv(I,J,K) * n2(I,J,K) &
    &                      + wu(I,J,K) * n3(I,J,K) ) &
    &        + n1(I,J,K) * ( wu(I,J,K) * n1(I,J,K) &
    &                      + vw(I,J,K) * n2(I,J,K) &
    &                      + ww(I,J,K) * n3(I,J,K) ) &
    &        )
    phiw2_11 = TMPw2 * (TMPw2km - 3.0 * n1(I,J,K) &
    &          * ( phi2_11 * n1(I,J,K) &
    &            + phi2_12 * n2(I,J,K) &
    &            + phi2_31 * n3(I,J,K) ) &
    &        )
    phiw2_22 = TMPw2 * (TMPw2km - 3.0 * n2(I,J,K) &
    &          * ( phi2_12 * n1(I,J,K) &
    &            + phi2_22 * n2(I,J,K) &
    &            + phi2_23 * n3(I,J,K) ) &
    &        )
    phiw2_33 = TMPw2 * (TMPw2km - 3.0 * n3(I,J,K) &
    &          * ( phi2_31 * n1(I,J,K) &
    &            + phi2_23 * n2(I,J,K) &
    &            + phi2_33 * n3(I,J,K) ) &
    &        )
    phiw2_12 = - 1.5 * TMPw2 * ( &
    &          n1(I,J,K) * ( phi2_12 * n1(I,J,K) &
    &                      + phi2_22 * n2(I,J,K) &
    &                      + phi2_23 * n3(I,J,K) ) &
    &        + n2(I,J,K) * ( phi2_11 * n1(I,J,K) &
    &                      + phi2_12 * n2(I,J,K) &
    &                      + phi2_31 * n3(I,J,K) ) &
    &        )
    phiw2_23 = - 1.5 * TMPw2 * ( &
    &          n2(I,J,K) * ( phi2_31 * n1(I,J,K) &
    &                      + phi2_23 * n2(I,J,K) &
    &                      + phi2_33 * n3(I,J,K) ) &
    &        + n3(I,J,K) * ( phi2_12 * n1(I,J,K) &
    &                      + phi2_22 * n2(I,J,K) &
    &                      + phi2_23 * n3(I,J,K) ) &
    &        )
    phiw2_31 = - 1.5 * TMPw2 * ( &
    &          n3(I,J,K) * ( phi2_11 * n1(I,J,K) &
    &                      + phi2_12 * n2(I,J,K) &
    &                      + phi2_31 * n3(I,J,K) ) &
    &        + n1(I,J,K) * ( phi2_31 * n1(I,J,K) &
    &                      + phi2_23 * n2(I,J,K) &
    &                      + phi2_33 * n3(I,J,K) ) &
    &        )
    phi11 = phi1_11 + phi2_11 + phiw1_11 + phiw2_11
    phi22 = phi1_22 + phi2_22 + phiw1_22 + phiw2_22
    phi33 = phi1_33 + phi2_33 + phiw1_33 + phiw2_33
    phi12 = phi1_12 + phi2_12 + phiw1_12 + phiw2_12
    phi23 = phi1_23 + phi2_23 + phiw1_23 + phiw2_23
    phi31 = phi1_31 + phi2_31 + phiw1_31 + phiw2_31
    ! 渦粘性係数の計算 -------------------------------------------------
    IF(e_k .GT. ZERO) THEN
      AMUT(I,J,K) = 0.09 * AK**2 / e_k * RHO(I,J,K)
    ENDIF
    ! 生成項、圧力再分配項、散逸項の和 ---------------------------------
    QH1 = RHO(I,J,K) / AJA(I,J,K)
    DQP(I,J,K, 6) = QH1 * (P11 - eps11 + phi11)
    DQP(I,J,K, 7) = QH1 * (P22 - eps22 + phi22)
    DQP(I,J,K, 8) = QH1 * (P33 - eps33 + phi33)
    DQP(I,J,K, 9) = QH1 * (P12 - eps12 + phi12)
    DQP(I,J,K,10) = QH1 * (P23 - eps23 + phi23)
    DQP(I,J,K,11) = QH1 * (P31 - eps31 + phi31)
    DQP(I,J,K,12) = QH1 * AEAK * ( &
    &               (C_e1 + psi1 + psi2) * P_k - C_e2 * EPST &
    &             )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH = 0.0
  SH = 0.0
  TH = 0.0
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL DIFFXI
  CALL DIFFET
  CALL DIFFZE
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I,J-1,K) .GT. 0.0) .AND. &
  &   (RHO(I-1,J,K) .GT. 0.0) .AND. (RHO(I,J,K)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J,K) .GT. 0.0) .AND. (RHO(I,J+1,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,K,L) = (-RH(I-1,J  ,K  ,L) + RH(I,J,K,L) ) &
    &            + (-SH(I  ,J-1,K  ,L) + SH(I,J,K,L) ) &
    &            + (-TH(I  ,J  ,K-1,L) + TH(I,J,K,L) )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項の計算(xi方向差分用)                                    ****
!***********************************************************************
SUBROUTINE DIFFXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: uuXI, uuET, uuZE, uuX, uuY, uuZ
  REAL    :: vvXI, vvET, vvZE, vvX, vvY, vvZ
  REAL    :: wwXI, wwET, wwZE, wwX, wwY, wwZ
  REAL    :: uvXI, uvET, uvZE, uvX, uvY, uvZ
  REAL    :: vwXI, vwET, vwZE, vwX, vwY, vwZ
  REAL    :: wuXI, wuET, wuZE, wuX, wuY, wuZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu
  REAL    :: S11, S22, S33, S12, S23, S31, SS
  REAL    :: C_p, AKAE
  REAL    :: R02M, S02M, T02M, R02T, S02T, T02T
  REAL    :: R03M, S03M, T03M, R03T, S03T, T03T
  REAL    :: R04M, S04M, T04M, R04T, S04T, T04T
  REAL    :: R05M, S05M, T05M, R05T, S05T, T05T
  REAL    :: R06M, S06M, T06M, R06T, S06T, T06T
  REAL    :: R07M, S07M, T07M, R07T, S07T, T07T
  REAL    :: R08M, S08M, T08M, R08T, S08T, T08T
  REAL    :: R09M, S09M, T09M, R09T, S09T, T09T
  REAL    :: R10M, S10M, T10M, R10T, S10T, T10T
  REAL    :: R11M, S11M, T11M, R11T, S11T, T11T
  REAL    :: R12M, S12M, T12M, R12T, S12T, T12T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  uuXI, uuET, uuZE, uuX, uuY, uuZ, &
  !$  vvXI, vvET, vvZE, vvX, vvY, vvZ, &
  !$  wwXI, wwET, wwZE, wwX, wwY, wwZ, &
  !$  uvXI, uvET, uvZE, uvX, uvY, uvZ, &
  !$  vwXI, vwET, vwZE, vwX, vwY, vwZ, &
  !$  wuXI, wuET, wuZE, wuX, wuY, wuZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, S11, S22, S33, S12, S23, S31, SS, C_p, AKAE, &
  !$  R02M, S02M, T02M, R02T, S02T, T02T, &
  !$  R03M, S03M, T03M, R03T, S03T, T03T, &
  !$  R04M, S04M, T04M, R04T, S04T, T04T, &
  !$  R05M, S05M, T05M, R05T, S05T, T05T, &
  !$  R06M, S06M, T06M, R06T, S06T, T06T, &
  !$  R07M, S07M, T07M, R07T, S07T, T07T, &
  !$  R08M, S08M, T08M, R08T, S08T, T08T, &
  !$  R09M, S09M, T09M, R09T, S09T, T09T, &
  !$  R10M, S10M, T10M, R10T, S10T, T10T, &
  !$  R11M, S11M, T11M, R11T, S11T, T11T, &
  !$  R12M, S12M, T12M, R12T, S12T, T12T, &
  !$  QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF((RHO(I,J,K) .GT. 0.0) .AND. (RHO(I+1,J,K) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I+1,J,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I+1,J,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I+1,J,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I+1,J,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I+1,J,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I+1,J,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I+1,J,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I+1,J,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I+1,J,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I+1,J,K))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J,K) + RHO(I+1,J,K))
    UM   = 0.5 * (  U(I,J,K) +   U(I+1,J,K))
    VM   = 0.5 * (  V(I,J,K) +   V(I+1,J,K))
    WM   = 0.5 * (  W(I,J,K) +   W(I+1,J,K))
    uuM  = 0.5 * ( uu(I,J,K) +  uu(I+1,J,K))
    vvM  = 0.5 * ( vv(I,J,K) +  vv(I+1,J,K))
    wwM  = 0.5 * ( ww(I,J,K) +  ww(I+1,J,K))
    uvM  = 0.5 * ( uv(I,J,K) +  uv(I+1,J,K))
    vwM  = 0.5 * ( vw(I,J,K) +  vw(I+1,J,K))
    wuM  = 0.5 * ( wu(I,J,K) +  wu(I+1,J,K))
    EPSM = 0.5 * (EPS(I,J,K) + EPS(I+1,J,K))
    AMUM = 0.5 * (AMU(I,J,K) + AMU(I+1,J,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = - U(I,J,K) + U(I+1,J,K)
    UET = 0.25 * ( - U(I,J-1,K) - U(I+1,J-1,K) &
    &              + U(I,J+1,K) + U(I+1,J+1,K) )
    UZE = 0.25 * ( - U(I,J,K-1) - U(I+1,J,K-1) &
    &              + U(I,J,K+1) + U(I+1,J,K+1) )
    VXI = - V(I,J,K) + V(I+1,J,K)
    VET = 0.25 * ( - V(I,J-1,K) - V(I+1,J-1,K) &
    &              + V(I,J+1,K) + V(I+1,J+1,K) )
    VZE = 0.25 * ( - V(I,J,K-1) - V(I+1,J,K-1) &
    &              + V(I,J,K+1) + V(I+1,J,K+1) )
    WXI = - W(I,J,K) + W(I+1,J,K)
    WET = 0.25 * ( - W(I,J-1,K) - W(I+1,J-1,K) &
    &              + W(I,J+1,K) + W(I+1,J+1,K) )
    WZE = 0.25 * ( - W(I,J,K-1) - W(I+1,J,K-1) &
    &              + W(I,J,K+1) + W(I+1,J,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = - T(I,J,K) + T(I+1,J,K)
    TET = 0.25 * ( - T(I,J-1,K) - T(I+1,J-1,K) &
    &              + T(I,J+1,K) + T(I+1,J+1,K) )
    TZE = 0.25 * ( - T(I,J,K-1) - T(I+1,J,K-1) &
    &              + T(I,J,K+1) + T(I+1,J,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = - uu(I,J,K) + uu(I+1,J,K)
    uuET = 0.25 * ( - uu(I,J-1,K) - uu(I+1,J-1,K) &
    &               + uu(I,J+1,K) + uu(I+1,J+1,K) )
    uuZE = 0.25 * ( - uu(I,J,K-1) - uu(I+1,J,K-1) &
    &               + uu(I,J,K+1) + uu(I+1,J,K+1) )
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM + uuZE * ZEXM
    uuY = uuXI * XIYM + uuET * ETYM + uuZE * ZEYM
    uuZ = uuXI * XIZM + uuET * ETZM + uuZE * ZEZM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = - vv(I,J,K) + vv(I+1,J,K)
    vvET = 0.25 * ( - vv(I,J-1,K) - vv(I+1,J-1,K) &
    &               + vv(I,J+1,K) + vv(I+1,J+1,K) )
    vvZE = 0.25 * ( - vv(I,J,K-1) - vv(I+1,J,K-1) &
    &               + vv(I,J,K+1) + vv(I+1,J,K+1) )
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM + vvZE * ZEXM
    vvY = vvXI * XIYM + vvET * ETYM + vvZE * ZEYM
    vvZ = vvXI * XIZM + vvET * ETZM + vvZE * ZEZM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = - ww(I,J,K) + ww(I+1,J,K)
    wwET = 0.25 * ( - ww(I,J-1,K) - ww(I+1,J-1,K) &
    &               + ww(I,J+1,K) + ww(I+1,J+1,K) )
    wwZE = 0.25 * ( - ww(I,J,K-1) - ww(I+1,J,K-1) &
    &               + ww(I,J,K+1) + ww(I+1,J,K+1) )
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM + wwZE * ZEXM
    wwY = wwXI * XIYM + wwET * ETYM + wwZE * ZEYM
    wwZ = wwXI * XIZM + wwET * ETZM + wwZE * ZEZM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = - uv(I,J,K) + uv(I+1,J,K)
    uvET = 0.25 * ( - uv(I,J-1,K) - uv(I+1,J-1,K) &
    &               + uv(I,J+1,K) + uv(I+1,J+1,K) )
    uvZE = 0.25 * ( - uv(I,J,K-1) - uv(I+1,J,K-1) &
    &               + uv(I,J,K+1) + uv(I+1,J,K+1) )
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM + uvZE * ZEXM
    uvY = uvXI * XIYM + uvET * ETYM + uvZE * ZEYM
    uvZ = uvXI * XIZM + uvET * ETZM + uvZE * ZEZM
    ! vw の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vwXI = - vw(I,J,K) + vw(I+1,J,K)
    vwET = 0.25 * ( - vw(I,J-1,K) - vw(I+1,J-1,K) &
    &               + vw(I,J+1,K) + vw(I+1,J+1,K) )
    vwZE = 0.25 * ( - vw(I,J,K-1) - vw(I+1,J,K-1) &
    &               + vw(I,J,K+1) + vw(I+1,J,K+1) )
    ! 物理空間方向一階微分
    vwX = vwXI * XIXM + vwET * ETXM + vwZE * ZEXM
    vwY = vwXI * XIYM + vwET * ETYM + vwZE * ZEYM
    vwZ = vwXI * XIZM + vwET * ETZM + vwZE * ZEZM
    ! wu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wuXI = - wu(I,J,K) + wu(I+1,J,K)
    wuET = 0.25 * ( - wu(I,J-1,K) - wu(I+1,J-1,K) &
    &               + wu(I,J+1,K) + wu(I+1,J+1,K) )
    wuZE = 0.25 * ( - wu(I,J,K-1) - wu(I+1,J,K-1) &
    &               + wu(I,J,K+1) + wu(I+1,J,K+1) )
    ! 物理空間方向一階微分
    wuX = wuXI * XIXM + wuET * ETXM + wuZE * ZEXM
    wuY = wuXI * XIYM + wuET * ETYM + wuZE * ZEYM
    wuZ = wuXI * XIZM + wuET * ETZM + wuZE * ZEZM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = - EPS(I,J,K) + EPS(I+1,J,K)
    AEET = 0.25 * ( - EPS(I,J-1,K) - EPS(I+1,J-1,K) &
    &               + EPS(I,J+1,K) + EPS(I+1,J+1,K) )
    AEZE = 0.25 * ( - EPS(I,J,K-1) - EPS(I+1,J,K-1) &
    &               + EPS(I,J,K+1) + EPS(I+1,J,K+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY + WZ) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 = 2.0 * WZ - SS
    S12 = UY + VX
    S23 = VZ + WY
    S31 = WX + UZ
    ! U の拡散 ---------------------------------------------------------
    R02M = nu * S11
    S02M = nu * S12
    T02M = nu * S31
    R02T =-uuM
    S02T =-uvM
    T02T =-wuM
    ! V の拡散 ---------------------------------------------------------
    R03M = nu * S12
    S03M = nu * S22
    T03M = nu * S23
    R03T =-uvM
    S03T =-vvM
    T03T =-vwM
    ! W の拡散 ---------------------------------------------------------
    R04M = nu * S31
    S04M = nu * S23
    T04M = nu * S33
    R04T =-wuM
    S04T =-vwM
    T04T =-wwM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R05M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM + S31 * WM)
    S05M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM + S23 * WM)
    T05M = C_p * nu / PR * TZ + nu * (S31 * UM + S23 * VM + S33 * WM)
    R05T = AKAE * C_se * C_p * (uuM * TX + uvM * TY + wuM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uuM * (uuX + vvX + wwX) &
    &    + uvM * (uuY + vvY + wwY) &
    &    + wuM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uuM * UM + uvM * VM + wuM * WM)
    S05T = AKAE * C_se * C_p * (uvM * TX + vvM * TY + vwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uvM * (uuX + vvX + wwX) &
    &    + vvM * (uuY + vvY + wwY) &
    &    + vwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uvM * UM + vvM * VM + vwM * WM)
    T05T = AKAE * C_se * C_p * (wuM * TX + vwM * TY + wwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      wuM * (uuX + vvX + wwX) &
    &    + vwM * (uuY + vvY + wwY) &
    &    + wwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (wuM * UM + vwM * VM + wwM * WM)
    ! uu の拡散 --------------------------------------------------------
    R06M = nu * uuX
    S06M = nu * uuY
    T06M = nu * uuZ
    R06T = C_s * AKAE * (uuM * uuX + uvM * uuY + wuM * uuZ)
    S06T = C_s * AKAE * (uvM * uuX + vvM * uuY + vwM * uuZ)
    T06T = C_s * AKAE * (wuM * uuX + vwM * uuY + wwM * uuZ)
    ! vv の拡散 --------------------------------------------------------
    R07M = nu * vvX
    S07M = nu * vvY
    T07M = nu * vvZ
    R07T = C_s * AKAE * (uuM * vvX + uvM * vvY + wuM * vvZ)
    S07T = C_s * AKAE * (uvM * vvX + vvM * vvY + vwM * vvZ)
    T07T = C_s * AKAE * (wuM * vvX + vwM * vvY + wwM * vvZ)
    ! ww の拡散 --------------------------------------------------------
    R08M = nu * wwX
    S08M = nu * wwY
    T08M = nu * wwZ
    R08T = C_s * AKAE * (uuM * wwX + uvM * wwY + wuM * wwZ)
    S08T = C_s * AKAE * (uvM * wwX + vvM * wwY + vwM * wwZ)
    T08T = C_s * AKAE * (wuM * wwX + vwM * wwY + wwM * wwZ)
    ! uv の拡散 --------------------------------------------------------
    R09M = nu * uvX
    S09M = nu * uvY
    T09M = nu * uvZ
    R09T = C_s * AKAE * (uuM * uvX + uvM * uvY + wuM * uvZ)
    S09T = C_s * AKAE * (uvM * uvX + vvM * uvY + vwM * uvZ)
    T09T = C_s * AKAE * (wuM * uvX + vwM * uvY + wwM * uvZ)
    ! vw の拡散 --------------------------------------------------------
    R10M = nu * vwX
    S10M = nu * vwY
    T10M = nu * vwZ
    R10T = C_s * AKAE * (uuM * vwX + uvM * vwY + wuM * vwZ)
    S10T = C_s * AKAE * (uvM * vwX + vvM * vwY + vwM * vwZ)
    T10T = C_s * AKAE * (wuM * vwX + vwM * vwY + wwM * vwZ)
    ! wu の拡散 --------------------------------------------------------
    R11M = nu * wuX
    S11M = nu * wuY
    T11M = nu * wuZ
    R11T = C_s * AKAE * (uuM * wuX + uvM * wuY + wuM * wuZ)
    S11T = C_s * AKAE * (uvM * wuX + vvM * wuY + vwM * wuZ)
    T11T = C_s * AKAE * (wuM * wuX + vwM * wuY + wwM * wuZ)
    ! epsilon の拡散 ---------------------------------------------------
    R12M = nu * AEX
    S12M = nu * AEY
    T12M = nu * AEZ
    R12T = C_e * AKAE * (uuM * AEX + uvM * AEY + wuM * AEZ)
    S12T = C_e * AKAE * (uvM * AEX + vvM * AEY + vwM * AEZ)
    T12T = C_e * AKAE * (wuM * AEX + vwM * AEY + wwM * AEZ)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    RH(I,J,K, 1) = 0.0
    RH(I,J,K, 2) = QH1 * ( (R02M + R02T) * XIXM &
    &                    + (S02M + S02T) * XIYM &
    &                    + (T02M + T02T) * XIZM )
    RH(I,J,K, 3) = QH1 * ( (R03M + R03T) * XIXM &
    &                    + (S03M + S03T) * XIYM &
    &                    + (T03M + T03T) * XIZM )
    RH(I,J,K, 4) = QH1 * ( (R04M + R04T) * XIXM &
    &                    + (S04M + S04T) * XIYM &
    &                    + (T04M + T04T) * XIZM )
    RH(I,J,K, 5) = QH1 * ( (R05M + R05T) * XIXM &
    &                    + (S05M + S05T) * XIYM &
    &                    + (T05M + T05T) * XIZM )
    RH(I,J,K, 6) = QH1 * ( (R06M + R06T) * XIXM &
    &                    + (S06M + S06T) * XIYM &
    &                    + (T06M + T06T) * XIZM )
    RH(I,J,K, 7) = QH1 * ( (R07M + R07T) * XIXM &
    &                    + (S07M + S07T) * XIYM &
    &                    + (T07M + T07T) * XIZM )
    RH(I,J,K, 8) = QH1 * ( (R08M + R08T) * XIXM &
    &                    + (S08M + S08T) * XIYM &
    &                    + (T08M + T08T) * XIZM )
    RH(I,J,K, 9) = QH1 * ( (R09M + R09T) * XIXM &
    &                    + (S09M + S09T) * XIYM &
    &                    + (T09M + T09T) * XIZM )
    RH(I,J,K,10) = QH1 * ( (R10M + R10T) * XIXM &
    &                    + (S10M + S10T) * XIYM &
    &                    + (T10M + T10T) * XIZM )
    RH(I,J,K,11) = QH1 * ( (R11M + R11T) * XIXM &
    &                    + (S11M + S11T) * XIYM &
    &                    + (T11M + T11T) * XIZM )
    RH(I,J,K,12) = QH1 * ( (R12M + R12T) * XIXM &
    &                    + (S12M + S12T) * XIYM &
    &                    + (T12M + T12T) * XIZM )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFXI
!***********************************************************************
!**** 拡散項の計算(eta方向差分用)                                   ****
!***********************************************************************
SUBROUTINE DIFFET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: uuXI, uuET, uuZE, uuX, uuY, uuZ
  REAL    :: vvXI, vvET, vvZE, vvX, vvY, vvZ
  REAL    :: wwXI, wwET, wwZE, wwX, wwY, wwZ
  REAL    :: uvXI, uvET, uvZE, uvX, uvY, uvZ
  REAL    :: vwXI, vwET, vwZE, vwX, vwY, vwZ
  REAL    :: wuXI, wuET, wuZE, wuX, wuY, wuZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu
  REAL    :: S11, S22, S33, S12, S23, S31, SS
  REAL    :: C_p, AKAE
  REAL    :: R02M, S02M, T02M, R02T, S02T, T02T
  REAL    :: R03M, S03M, T03M, R03T, S03T, T03T
  REAL    :: R04M, S04M, T04M, R04T, S04T, T04T
  REAL    :: R05M, S05M, T05M, R05T, S05T, T05T
  REAL    :: R06M, S06M, T06M, R06T, S06T, T06T
  REAL    :: R07M, S07M, T07M, R07T, S07T, T07T
  REAL    :: R08M, S08M, T08M, R08T, S08T, T08T
  REAL    :: R09M, S09M, T09M, R09T, S09T, T09T
  REAL    :: R10M, S10M, T10M, R10T, S10T, T10T
  REAL    :: R11M, S11M, T11M, R11T, S11T, T11T
  REAL    :: R12M, S12M, T12M, R12T, S12T, T12T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  uuXI, uuET, uuZE, uuX, uuY, uuZ, &
  !$  vvXI, vvET, vvZE, vvX, vvY, vvZ, &
  !$  wwXI, wwET, wwZE, wwX, wwY, wwZ, &
  !$  uvXI, uvET, uvZE, uvX, uvY, uvZ, &
  !$  vwXI, vwET, vwZE, vwX, vwY, vwZ, &
  !$  wuXI, wuET, wuZE, wuX, wuY, wuZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, S11, S22, S33, S12, S23, S31, SS, C_p, AKAE, &
  !$  R02M, S02M, T02M, R02T, S02T, T02T, &
  !$  R03M, S03M, T03M, R03T, S03T, T03T, &
  !$  R04M, S04M, T04M, R04T, S04T, T04T, &
  !$  R05M, S05M, T05M, R05T, S05T, T05T, &
  !$  R06M, S06M, T06M, R06T, S06T, T06T, &
  !$  R07M, S07M, T07M, R07T, S07T, T07T, &
  !$  R08M, S08M, T08M, R08T, S08T, T08T, &
  !$  R09M, S09M, T09M, R09T, S09T, T09T, &
  !$  R10M, S10M, T10M, R10T, S10T, T10T, &
  !$  R11M, S11M, T11M, R11T, S11T, T11T, &
  !$  R12M, S12M, T12M, R12T, S12T, T12T, &
  !$  QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF((RHO(I,J,K) .GT. 0.0) .AND. (RHO(I,J+1,K) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J+1,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J+1,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J+1,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J+1,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J+1,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J+1,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J+1,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J+1,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J+1,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J+1,K))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J,K) + RHO(I,J+1,K))
    UM   = 0.5 * (  U(I,J,K) +   U(I,J+1,K))
    VM   = 0.5 * (  V(I,J,K) +   V(I,J+1,K))
    WM   = 0.5 * (  W(I,J,K) +   W(I,J+1,K))
    uuM  = 0.5 * ( uu(I,J,K) +  uu(I,J+1,K))
    vvM  = 0.5 * ( vv(I,J,K) +  vv(I,J+1,K))
    wwM  = 0.5 * ( ww(I,J,K) +  ww(I,J+1,K))
    uvM  = 0.5 * ( uv(I,J,K) +  uv(I,J+1,K))
    vwM  = 0.5 * ( vw(I,J,K) +  vw(I,J+1,K))
    wuM  = 0.5 * ( wu(I,J,K) +  wu(I,J+1,K))
    EPSM = 0.5 * (EPS(I,J,K) + EPS(I,J+1,K))
    AMUM = 0.5 * (AMU(I,J,K) + AMU(I,J+1,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J+1,K) &
    &              + U(I+1,J,K) + U(I+1,J+1,K) )
    UET = - U(I,J,K) + U(I,J+1,K)
    UZE = 0.25 * ( - U(I,J,K-1) - U(I,J+1,K-1) &
    &              + U(I,J,K+1) + U(I,J+1,K+1) )
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J+1,K) &
    &              + V(I+1,J,K) + V(I+1,J+1,K) )
    VET = - V(I,J,K) + V(I,J+1,K)
    VZE = 0.25 * ( - V(I,J,K-1) - V(I,J+1,K-1) &
    &              + V(I,J,K+1) + V(I,J+1,K+1) )
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J+1,K) &
    &              + W(I+1,J,K) + W(I+1,J+1,K) )
    WET = - W(I,J,K) + W(I,J+1,K)
    WZE = 0.25 * ( - W(I,J,K-1) - W(I,J+1,K-1) &
    &              + W(I,J,K+1) + W(I,J+1,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J+1,K) &
    &              + T(I+1,J,K) + T(I+1,J+1,K) )
    TET = - T(I,J,K) + T(I,J+1,K)
    TZE = 0.25 * ( - T(I,J,K-1) - T(I,J+1,K-1) &
    &              + T(I,J,K+1) + T(I,J+1,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = 0.25 * ( - uu(I-1,J,K) - uu(I-1,J+1,K) &
    &               + uu(I+1,J,K) + uu(I+1,J+1,K) )
    uuET = - uu(I,J,K) + uu(I,J+1,K)
    uuZE = 0.25 * ( - uu(I,J,K-1) - uu(I,J+1,K-1) &
    &               + uu(I,J,K+1) + uu(I,J+1,K+1) )
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM + uuZE * ZEXM
    uuY = uuXI * XIYM + uuET * ETYM + uuZE * ZEYM
    uuZ = uuXI * XIZM + uuET * ETZM + uuZE * ZEZM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = 0.25 * ( - vv(I-1,J,K) - vv(I-1,J+1,K) &
    &               + vv(I+1,J,K) + vv(I+1,J+1,K) )
    vvET = - vv(I,J,K) + vv(I,J+1,K)
    vvZE = 0.25 * ( - vv(I,J,K-1) - vv(I,J+1,K-1) &
    &               + vv(I,J,K+1) + vv(I,J+1,K+1) )
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM + vvZE * ZEXM
    vvY = vvXI * XIYM + vvET * ETYM + vvZE * ZEYM
    vvZ = vvXI * XIZM + vvET * ETZM + vvZE * ZEZM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = 0.25 * ( - ww(I-1,J,K) - ww(I-1,J+1,K) &
    &               + ww(I+1,J,K) + ww(I+1,J+1,K) )
    wwET = - ww(I,J,K) + ww(I,J+1,K)
    wwZE = 0.25 * ( - ww(I,J,K-1) - ww(I,J+1,K-1) &
    &               + ww(I,J,K+1) + ww(I,J+1,K+1) )
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM + wwZE * ZEXM
    wwY = wwXI * XIYM + wwET * ETYM + wwZE * ZEYM
    wwZ = wwXI * XIZM + wwET * ETZM + wwZE * ZEZM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = 0.25 * ( - uv(I-1,J,K) - uv(I-1,J+1,K) &
    &               + uv(I+1,J,K) + uv(I+1,J+1,K) )
    uvET = - uv(I,J,K) + uv(I,J+1,K)
    uvZE = 0.25 * ( - uv(I,J,K-1) - uv(I,J+1,K-1) &
    &               + uv(I,J,K+1) + uv(I,J+1,K+1) )
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM + uvZE * ZEXM
    uvY = uvXI * XIYM + uvET * ETYM + uvZE * ZEYM
    uvZ = uvXI * XIZM + uvET * ETZM + uvZE * ZEZM
    ! vw の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vwXI = 0.25 * ( - vw(I-1,J,K) - vw(I-1,J+1,K) &
    &               + vw(I+1,J,K) + vw(I+1,J+1,K) )
    vwET = - vw(I,J,K) + vw(I,J+1,K)
    vwZE = 0.25 * ( - vw(I,J,K-1) - vw(I,J+1,K-1) &
    &               + vw(I,J,K+1) + vw(I,J+1,K+1) )
    ! 物理空間方向一階微分
    vwX = vwXI * XIXM + vwET * ETXM + vwZE * ZEXM
    vwY = vwXI * XIYM + vwET * ETYM + vwZE * ZEYM
    vwZ = vwXI * XIZM + vwET * ETZM + vwZE * ZEZM
    ! wu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wuXI = 0.25 * ( - wu(I-1,J,K) - wu(I-1,J+1,K) &
    &               + wu(I+1,J,K) + wu(I+1,J+1,K) )
    wuET = - wu(I,J,K) + wu(I,J+1,K)
    wuZE = 0.25 * ( - wu(I,J,K-1) - wu(I,J+1,K-1) &
    &               + wu(I,J,K+1) + wu(I,J+1,K+1) )
    ! 物理空間方向一階微分
    wuX = wuXI * XIXM + wuET * ETXM + wuZE * ZEXM
    wuY = wuXI * XIYM + wuET * ETYM + wuZE * ZEYM
    wuZ = wuXI * XIZM + wuET * ETZM + wuZE * ZEZM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPS(I-1,J,K) - EPS(I-1,J+1,K) &
    &               + EPS(I+1,J,K) + EPS(I+1,J+1,K) )
    AEET = - EPS(I,J,K) + EPS(I,J+1,K)
    AEZE = 0.25 * ( - EPS(I,J,K-1) - EPS(I,J+1,K-1) &
    &               + EPS(I,J,K+1) + EPS(I,J+1,K+1) )
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY + WZ) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 = 2.0 * WZ - SS
    S12 = UY + VX
    S23 = VZ + WY
    S31 = WX + UZ
    ! U の拡散 ---------------------------------------------------------
    R02M = nu * S11
    S02M = nu * S12
    T02M = nu * S31
    R02T =-uuM
    S02T =-uvM
    T02T =-wuM
    ! V の拡散 ---------------------------------------------------------
    R03M = nu * S12
    S03M = nu * S22
    T03M = nu * S23
    R03T =-uvM
    S03T =-vvM
    T03T =-vwM
    ! W の拡散 ---------------------------------------------------------
    R04M = nu * S31
    S04M = nu * S23
    T04M = nu * S33
    R04T =-wuM
    S04T =-vwM
    T04T =-wwM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R05M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM + S31 * WM)
    S05M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM + S23 * WM)
    T05M = C_p * nu / PR * TZ + nu * (S31 * UM + S23 * VM + S33 * WM)
    R05T = AKAE * C_se * C_p * (uuM * TX + uvM * TY + wuM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uuM * (uuX + vvX + wwX) &
    &    + uvM * (uuY + vvY + wwY) &
    &    + wuM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uuM * UM + uvM * VM + wuM * WM)
    S05T = AKAE * C_se * C_p * (uvM * TX + vvM * TY + vwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uvM * (uuX + vvX + wwX) &
    &    + vvM * (uuY + vvY + wwY) &
    &    + vwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uvM * UM + vvM * VM + vwM * WM)
    T05T = AKAE * C_se * C_p * (wuM * TX + vwM * TY + wwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      wuM * (uuX + vvX + wwX) &
    &    + vwM * (uuY + vvY + wwY) &
    &    + wwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (wuM * UM + vwM * VM + wwM * WM)
    ! uu の拡散 --------------------------------------------------------
    R06M = nu * uuX
    S06M = nu * uuY
    T06M = nu * uuZ
    R06T = C_s * AKAE * (uuM * uuX + uvM * uuY + wuM * uuZ)
    S06T = C_s * AKAE * (uvM * uuX + vvM * uuY + vwM * uuZ)
    T06T = C_s * AKAE * (wuM * uuX + vwM * uuY + wwM * uuZ)
    ! vv の拡散 --------------------------------------------------------
    R07M = nu * vvX
    S07M = nu * vvY
    T07M = nu * vvZ
    R07T = C_s * AKAE * (uuM * vvX + uvM * vvY + wuM * vvZ)
    S07T = C_s * AKAE * (uvM * vvX + vvM * vvY + vwM * vvZ)
    T07T = C_s * AKAE * (wuM * vvX + vwM * vvY + wwM * vvZ)
    ! ww の拡散 --------------------------------------------------------
    R08M = nu * wwX
    S08M = nu * wwY
    T08M = nu * wwZ
    R08T = C_s * AKAE * (uuM * wwX + uvM * wwY + wuM * wwZ)
    S08T = C_s * AKAE * (uvM * wwX + vvM * wwY + vwM * wwZ)
    T08T = C_s * AKAE * (wuM * wwX + vwM * wwY + wwM * wwZ)
    ! uv の拡散 --------------------------------------------------------
    R09M = nu * uvX
    S09M = nu * uvY
    T09M = nu * uvZ
    R09T = C_s * AKAE * (uuM * uvX + uvM * uvY + wuM * uvZ)
    S09T = C_s * AKAE * (uvM * uvX + vvM * uvY + vwM * uvZ)
    T09T = C_s * AKAE * (wuM * uvX + vwM * uvY + wwM * uvZ)
    ! vw の拡散 --------------------------------------------------------
    R10M = nu * vwX
    S10M = nu * vwY
    T10M = nu * vwZ
    R10T = C_s * AKAE * (uuM * vwX + uvM * vwY + wuM * vwZ)
    S10T = C_s * AKAE * (uvM * vwX + vvM * vwY + vwM * vwZ)
    T10T = C_s * AKAE * (wuM * vwX + vwM * vwY + wwM * vwZ)
    ! wu の拡散 --------------------------------------------------------
    R11M = nu * wuX
    S11M = nu * wuY
    T11M = nu * wuZ
    R11T = C_s * AKAE * (uuM * wuX + uvM * wuY + wuM * wuZ)
    S11T = C_s * AKAE * (uvM * wuX + vvM * wuY + vwM * wuZ)
    T11T = C_s * AKAE * (wuM * wuX + vwM * wuY + wwM * wuZ)
    ! epsilon の拡散 ---------------------------------------------------
    R12M = nu * AEX
    S12M = nu * AEY
    T12M = nu * AEZ
    R12T = C_e * AKAE * (uuM * AEX + uvM * AEY + wuM * AEZ)
    S12T = C_e * AKAE * (uvM * AEX + vvM * AEY + vwM * AEZ)
    T12T = C_e * AKAE * (wuM * AEX + vwM * AEY + wwM * AEZ)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    SH(I,J,K, 1) = 0.0
    SH(I,J,K, 2) = QH1 * ( (R02M + R02T) * ETXM &
    &                    + (S02M + S02T) * ETYM &
    &                    + (T02M + T02T) * ETZM )
    SH(I,J,K, 3) = QH1 * ( (R03M + R03T) * ETXM &
    &                    + (S03M + S03T) * ETYM &
    &                    + (T03M + T03T) * ETZM )
    SH(I,J,K, 4) = QH1 * ( (R04M + R04T) * ETXM &
    &                    + (S04M + S04T) * ETYM &
    &                    + (T04M + T04T) * ETZM )
    SH(I,J,K, 5) = QH1 * ( (R05M + R05T) * ETXM &
    &                    + (S05M + S05T) * ETYM &
    &                    + (T05M + T05T) * ETZM )
    SH(I,J,K, 6) = QH1 * ( (R06M + R06T) * ETXM &
    &                    + (S06M + S06T) * ETYM &
    &                    + (T06M + T06T) * ETZM )
    SH(I,J,K, 7) = QH1 * ( (R07M + R07T) * ETXM &
    &                    + (S07M + S07T) * ETYM &
    &                    + (T07M + T07T) * ETZM )
    SH(I,J,K, 8) = QH1 * ( (R08M + R08T) * ETXM &
    &                    + (S08M + S08T) * ETYM &
    &                    + (T08M + T08T) * ETZM )
    SH(I,J,K, 9) = QH1 * ( (R09M + R09T) * ETXM &
    &                    + (S09M + S09T) * ETYM &
    &                    + (T09M + T09T) * ETZM )
    SH(I,J,K,10) = QH1 * ( (R10M + R10T) * ETXM &
    &                    + (S10M + S10T) * ETYM &
    &                    + (T10M + T10T) * ETZM )
    SH(I,J,K,11) = QH1 * ( (R11M + R11T) * ETXM &
    &                    + (S11M + S11T) * ETYM &
    &                    + (T11M + T11T) * ETZM )
    SH(I,J,K,12) = QH1 * ( (R12M + R12T) * ETXM &
    &                    + (S12M + S12T) * ETYM &
    &                    + (T12M + T12T) * ETZM )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFET
!***********************************************************************
!**** 拡散項の計算(zeta方向差分用)                                  ****
!***********************************************************************
SUBROUTINE DIFFZE
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: uuXI, uuET, uuZE, uuX, uuY, uuZ
  REAL    :: vvXI, vvET, vvZE, vvX, vvY, vvZ
  REAL    :: wwXI, wwET, wwZE, wwX, wwY, wwZ
  REAL    :: uvXI, uvET, uvZE, uvX, uvY, uvZ
  REAL    :: vwXI, vwET, vwZE, vwX, vwY, vwZ
  REAL    :: wuXI, wuET, wuZE, wuX, wuY, wuZ
  REAL    :: AEXI, AEET, AEZE, AEX, AEY, AEZ
  REAL    :: nu
  REAL    :: S11, S22, S33, S12, S23, S31, SS
  REAL    :: C_p, AKAE
  REAL    :: R02M, S02M, T02M, R02T, S02T, T02T
  REAL    :: R03M, S03M, T03M, R03T, S03T, T03T
  REAL    :: R04M, S04M, T04M, R04T, S04T, T04T
  REAL    :: R05M, S05M, T05M, R05T, S05T, T05T
  REAL    :: R06M, S06M, T06M, R06T, S06T, T06T
  REAL    :: R07M, S07M, T07M, R07T, S07T, T07T
  REAL    :: R08M, S08M, T08M, R08T, S08T, T08T
  REAL    :: R09M, S09M, T09M, R09T, S09T, T09T
  REAL    :: R10M, S10M, T10M, R10T, S10T, T10T
  REAL    :: R11M, S11M, T11M, R11T, S11T, T11T
  REAL    :: R12M, S12M, T12M, R12T, S12T, T12T
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, uuM, vvM, wwM, uvM, vwM, wuM, EPSM, AMUM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  uuXI, uuET, uuZE, uuX, uuY, uuZ, &
  !$  vvXI, vvET, vvZE, vvX, vvY, vvZ, &
  !$  wwXI, wwET, wwZE, wwX, wwY, wwZ, &
  !$  uvXI, uvET, uvZE, uvX, uvY, uvZ, &
  !$  vwXI, vwET, vwZE, vwX, vwY, vwZ, &
  !$  wuXI, wuET, wuZE, wuX, wuY, wuZ, &
  !$  AEXI, AEET, AEZE, AEX, AEY, AEZ, &
  !$  nu, S11, S22, S33, S12, S23, S31, SS, C_p, AKAE, &
  !$  R02M, S02M, T02M, R02T, S02T, T02T, &
  !$  R03M, S03M, T03M, R03T, S03T, T03T, &
  !$  R04M, S04M, T04M, R04T, S04T, T04T, &
  !$  R05M, S05M, T05M, R05T, S05T, T05T, &
  !$  R06M, S06M, T06M, R06T, S06T, T06T, &
  !$  R07M, S07M, T07M, R07T, S07T, T07T, &
  !$  R08M, S08M, T08M, R08T, S08T, T08T, &
  !$  R09M, S09M, T09M, R09T, S09T, T09T, &
  !$  R10M, S10M, T10M, R10T, S10T, T10T, &
  !$  R11M, S11M, T11M, R11T, S11T, T11T, &
  !$  R12M, S12M, T12M, R12T, S12T, T12T, &
  !$  QH1 &
  !$)
  DO K = KS    , KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF((RHO(I,J,K) .GT. 0.0) .AND. (RHO(I,J,K+1) .GT. 0.0)) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J,K+1))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J,K+1))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J,K+1))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J,K+1))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J,K+1))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J,K+1))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J,K+1))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J,K+1))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J,K+1))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J,K+1))
    ! 物理量
    RHOM = 0.5 * (RHO(I,J,K) + RHO(I,J,K+1))
    UM   = 0.5 * (  U(I,J,K) +   U(I,J,K+1))
    VM   = 0.5 * (  V(I,J,K) +   V(I,J,K+1))
    WM   = 0.5 * (  W(I,J,K) +   W(I,J,K+1))
    uuM  = 0.5 * ( uu(I,J,K) +  uu(I,J,K+1))
    vvM  = 0.5 * ( vv(I,J,K) +  vv(I,J,K+1))
    wwM  = 0.5 * ( ww(I,J,K) +  ww(I,J,K+1))
    uvM  = 0.5 * ( uv(I,J,K) +  uv(I,J,K+1))
    vwM  = 0.5 * ( vw(I,J,K) +  vw(I,J,K+1))
    wuM  = 0.5 * ( wu(I,J,K) +  wu(I,J,K+1))
    EPSM = 0.5 * (EPS(I,J,K) + EPS(I,J,K+1))
    AMUM = 0.5 * (AMU(I,J,K) + AMU(I,J,K+1))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J,K+1) &
    &              + U(I+1,J,K) + U(I+1,J,K+1) )
    UET = 0.25 * ( - U(I,J-1,K) - U(I,J-1,K+1) &
    &              + U(I,J+1,K) + U(I,J+1,K+1) )
    UZE = - U(I,J,K) + U(I,J,K+1)
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J,K+1) &
    &              + V(I+1,J,K) + V(I+1,J,K+1) )
    VET = 0.25 * ( - V(I,J-1,K) - V(I,J-1,K+1) &
    &              + V(I,J+1,K) + V(I,J+1,K+1) )
    VZE = - V(I,J,K) + V(I,J,K+1)
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J,K+1) &
    &              + W(I+1,J,K) + W(I+1,J,K+1) )
    WET = 0.25 * ( - W(I,J-1,K) - W(I,J-1,K+1) &
    &              + W(I,J+1,K) + W(I,J+1,K+1) )
    WZE = - W(I,J,K) + W(I,J,K+1)
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J,K+1) &
    &              + T(I+1,J,K) + T(I+1,J,K+1) )
    TET = 0.25 * ( - T(I,J-1,K) - T(I,J-1,K+1) &
    &              + T(I,J+1,K) + T(I,J+1,K+1) )
    TZE = - T(I,J,K) + T(I,J,K+1)
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! uu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uuXI = 0.25 * ( - uu(I-1,J,K) - uu(I-1,J,K+1) &
    &               + uu(I+1,J,K) + uu(I+1,J,K+1) )
    uuET = 0.25 * ( - uu(I,J-1,K) - uu(I,J-1,K+1) &
    &               + uu(I,J+1,K) + uu(I,J+1,K+1) )
    uuZE = - uu(I,J,K) + uu(I,J,K+1)
    ! 物理空間方向一階微分
    uuX = uuXI * XIXM + uuET * ETXM + uuZE * ZEXM
    uuY = uuXI * XIYM + uuET * ETYM + uuZE * ZEYM
    uuZ = uuXI * XIZM + uuET * ETZM + uuZE * ZEZM
    ! vv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vvXI = 0.25 * ( - vv(I-1,J,K) - vv(I-1,J,K+1) &
    &               + vv(I+1,J,K) + vv(I+1,J,K+1) )
    vvET = 0.25 * ( - vv(I,J-1,K) - vv(I,J-1,K+1) &
    &               + vv(I,J+1,K) + vv(I,J+1,K+1) )
    vvZE = - vv(I,J,K) + vv(I,J,K+1)
    ! 物理空間方向一階微分
    vvX = vvXI * XIXM + vvET * ETXM + vvZE * ZEXM
    vvY = vvXI * XIYM + vvET * ETYM + vvZE * ZEYM
    vvZ = vvXI * XIZM + vvET * ETZM + vvZE * ZEZM
    ! ww の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wwXI = 0.25 * ( - ww(I-1,J,K) - ww(I-1,J,K+1) &
    &               + ww(I+1,J,K) + ww(I+1,J,K+1) )
    wwET = 0.25 * ( - ww(I,J-1,K) - ww(I,J-1,K+1) &
    &               + ww(I,J+1,K) + ww(I,J+1,K+1) )
    wwZE = - ww(I,J,K) + ww(I,J,K+1)
    ! 物理空間方向一階微分
    wwX = wwXI * XIXM + wwET * ETXM + wwZE * ZEXM
    wwY = wwXI * XIYM + wwET * ETYM + wwZE * ZEYM
    wwZ = wwXI * XIZM + wwET * ETZM + wwZE * ZEZM
    ! uv の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    uvXI = 0.25 * ( - uv(I-1,J,K) - uv(I-1,J,K+1) &
    &               + uv(I+1,J,K) + uv(I+1,J,K+1) )
    uvET = 0.25 * ( - uv(I,J-1,K) - uv(I,J-1,K+1) &
    &               + uv(I,J+1,K) + uv(I,J+1,K+1) )
    uvZE = - uv(I,J,K) + uv(I,J,K+1)
    ! 物理空間方向一階微分
    uvX = uvXI * XIXM + uvET * ETXM + uvZE * ZEXM
    uvY = uvXI * XIYM + uvET * ETYM + uvZE * ZEYM
    uvZ = uvXI * XIZM + uvET * ETZM + uvZE * ZEZM
    ! vw の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    vwXI = 0.25 * ( - vw(I-1,J,K) - vw(I-1,J,K+1) &
    &               + vw(I+1,J,K) + vw(I+1,J,K+1) )
    vwET = 0.25 * ( - vw(I,J-1,K) - vw(I,J-1,K+1) &
    &               + vw(I,J+1,K) + vw(I,J+1,K+1) )
    vwZE = - vw(I,J,K) + vw(I,J,K+1)
    ! 物理空間方向一階微分
    vwX = vwXI * XIXM + vwET * ETXM + vwZE * ZEXM
    vwY = vwXI * XIYM + vwET * ETYM + vwZE * ZEYM
    vwZ = vwXI * XIZM + vwET * ETZM + vwZE * ZEZM
    ! wu の一階微分の計算 ----------------------------------------------
    ! 計算空間方向一階微分
    wuXI = 0.25 * ( - wu(I-1,J,K) - wu(I-1,J,K+1) &
    &               + wu(I+1,J,K) + wu(I+1,J,K+1) )
    wuET = 0.25 * ( - wu(I,J-1,K) - wu(I,J-1,K+1) &
    &               + wu(I,J+1,K) + wu(I,J+1,K+1) )
    wuZE = - wu(I,J,K) + wu(I,J,K+1)
    ! 物理空間方向一階微分
    wuX = wuXI * XIXM + wuET * ETXM + wuZE * ZEXM
    wuY = wuXI * XIYM + wuET * ETYM + wuZE * ZEYM
    wuZ = wuXI * XIZM + wuET * ETZM + wuZE * ZEZM
    ! epsilon の一階微分の計算 -----------------------------------------
    ! 計算空間方向一階微分
    AEXI = 0.25 * ( - EPS(I-1,J,K) - EPS(I-1,J,K+1) &
    &               + EPS(I+1,J,K) + EPS(I+1,J,K+1) )
    AEET = 0.25 * ( - EPS(I,J-1,K) - EPS(I,J-1,K+1) &
    &               + EPS(I,J+1,K) + EPS(I,J+1,K+1) )
    AEZE = - EPS(I,J,K) + EPS(I,J,K+1)
    ! 物理空間方向一階微分
    AEX = AEXI * XIXM + AEET * ETXM + AEZE * ZEXM
    AEY = AEXI * XIYM + AEET * ETYM + AEZE * ZEYM
    AEZ = AEXI * XIZM + AEET * ETZM + AEZE * ZEZM
    ! 動粘度 -----------------------------------------------------------
    nu = AMUM / RHOM
    ! 歪み速度 ---------------------------------------------------------
    SS  = (UX + VY + WZ) * TwoThird
    S11 = 2.0 * UX - SS
    S22 = 2.0 * VY - SS
    S33 = 2.0 * WZ - SS
    S12 = UY + VX
    S23 = VZ + WY
    S31 = WX + UZ
    ! U の拡散 ---------------------------------------------------------
    R02M = nu * S11
    S02M = nu * S12
    T02M = nu * S31
    R02T =-uuM
    S02T =-uvM
    T02T =-wuM
    ! V の拡散 ---------------------------------------------------------
    R03M = nu * S12
    S03M = nu * S22
    T03M = nu * S23
    R03T =-uvM
    S03T =-vvM
    T03T =-vwM
    ! W の拡散 ---------------------------------------------------------
    R04M = nu * S31
    S04M = nu * S23
    T04M = nu * S33
    R04T =-wuM
    S04T =-vwM
    T04T =-wwM
    ! エネルギーの拡散 -------------------------------------------------
    C_p  = GAMMA * RG / (GAMMA - 1.0)
    IF(EPSM .GT. ZERO) THEN
      AKAE = 0.5 * (uuM + vvM + wwM) / EPSM
    ELSE
      AKAE = 0.0
    ENDIF
    R05M = C_p * nu / PR * TX + nu * (S11 * UM + S12 * VM + S31 * WM)
    S05M = C_p * nu / PR * TY + nu * (S12 * UM + S22 * VM + S23 * WM)
    T05M = C_p * nu / PR * TZ + nu * (S31 * UM + S23 * VM + S33 * WM)
    R05T = AKAE * C_se * C_p * (uuM * TX + uvM * TY + wuM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uuM * (uuX + vvX + wwX) &
    &    + uvM * (uuY + vvY + wwY) &
    &    + wuM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uuM * UM + uvM * VM + wuM * WM)
    S05T = AKAE * C_se * C_p * (uvM * TX + vvM * TY + vwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      uvM * (uuX + vvX + wwX) &
    &    + vvM * (uuY + vvY + wwY) &
    &    + vwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (uvM * UM + vvM * VM + vwM * WM)
    T05T = AKAE * C_se * C_p * (wuM * TX + vwM * TY + wwM * TZ) &
    &    + AKAE * C_s  * ( &
    &      wuM * (uuX + vvX + wwX) &
    &    + vwM * (uuY + vvY + wwY) &
    &    + wwM * (uuZ + vvZ + wwZ) &
    &    ) &
    &    - (wuM * UM + vwM * VM + wwM * WM)
    ! uu の拡散 --------------------------------------------------------
    R06M = nu * uuX
    S06M = nu * uuY
    T06M = nu * uuZ
    R06T = C_s * AKAE * (uuM * uuX + uvM * uuY + wuM * uuZ)
    S06T = C_s * AKAE * (uvM * uuX + vvM * uuY + vwM * uuZ)
    T06T = C_s * AKAE * (wuM * uuX + vwM * uuY + wwM * uuZ)
    ! vv の拡散 --------------------------------------------------------
    R07M = nu * vvX
    S07M = nu * vvY
    T07M = nu * vvZ
    R07T = C_s * AKAE * (uuM * vvX + uvM * vvY + wuM * vvZ)
    S07T = C_s * AKAE * (uvM * vvX + vvM * vvY + vwM * vvZ)
    T07T = C_s * AKAE * (wuM * vvX + vwM * vvY + wwM * vvZ)
    ! ww の拡散 --------------------------------------------------------
    R08M = nu * wwX
    S08M = nu * wwY
    T08M = nu * wwZ
    R08T = C_s * AKAE * (uuM * wwX + uvM * wwY + wuM * wwZ)
    S08T = C_s * AKAE * (uvM * wwX + vvM * wwY + vwM * wwZ)
    T08T = C_s * AKAE * (wuM * wwX + vwM * wwY + wwM * wwZ)
    ! uv の拡散 --------------------------------------------------------
    R09M = nu * uvX
    S09M = nu * uvY
    T09M = nu * uvZ
    R09T = C_s * AKAE * (uuM * uvX + uvM * uvY + wuM * uvZ)
    S09T = C_s * AKAE * (uvM * uvX + vvM * uvY + vwM * uvZ)
    T09T = C_s * AKAE * (wuM * uvX + vwM * uvY + wwM * uvZ)
    ! vw の拡散 --------------------------------------------------------
    R10M = nu * vwX
    S10M = nu * vwY
    T10M = nu * vwZ
    R10T = C_s * AKAE * (uuM * vwX + uvM * vwY + wuM * vwZ)
    S10T = C_s * AKAE * (uvM * vwX + vvM * vwY + vwM * vwZ)
    T10T = C_s * AKAE * (wuM * vwX + vwM * vwY + wwM * vwZ)
    ! wu の拡散 --------------------------------------------------------
    R11M = nu * wuX
    S11M = nu * wuY
    T11M = nu * wuZ
    R11T = C_s * AKAE * (uuM * wuX + uvM * wuY + wuM * wuZ)
    S11T = C_s * AKAE * (uvM * wuX + vvM * wuY + vwM * wuZ)
    T11T = C_s * AKAE * (wuM * wuX + vwM * wuY + wwM * wuZ)
    ! epsilon の拡散 ---------------------------------------------------
    R12M = nu * AEX
    S12M = nu * AEY
    T12M = nu * AEZ
    R12T = C_e * AKAE * (uuM * AEX + uvM * AEY + wuM * AEZ)
    S12T = C_e * AKAE * (uvM * AEX + vvM * AEY + vwM * AEZ)
    T12T = C_e * AKAE * (wuM * AEX + vwM * AEY + wwM * AEZ)
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    TH(I,J,K, 1) = 0.0
    TH(I,J,K, 2) = QH1 * ( (R02M + R02T) * ZEXM &
    &                    + (S02M + S02T) * ZEYM &
    &                    + (T02M + T02T) * ZEZM )
    TH(I,J,K, 3) = QH1 * ( (R03M + R03T) * ZEXM &
    &                    + (S03M + S03T) * ZEYM &
    &                    + (T03M + T03T) * ZEZM )
    TH(I,J,K, 4) = QH1 * ( (R04M + R04T) * ZEXM &
    &                    + (S04M + S04T) * ZEYM &
    &                    + (T04M + T04T) * ZEZM )
    TH(I,J,K, 5) = QH1 * ( (R05M + R05T) * ZEXM &
    &                    + (S05M + S05T) * ZEYM &
    &                    + (T05M + T05T) * ZEZM )
    TH(I,J,K, 6) = QH1 * ( (R06M + R06T) * ZEXM &
    &                    + (S06M + S06T) * ZEYM &
    &                    + (T06M + T06T) * ZEZM )
    TH(I,J,K, 7) = QH1 * ( (R07M + R07T) * ZEXM &
    &                    + (S07M + S07T) * ZEYM &
    &                    + (T07M + T07T) * ZEZM )
    TH(I,J,K, 8) = QH1 * ( (R08M + R08T) * ZEXM &
    &                    + (S08M + S08T) * ZEYM &
    &                    + (T08M + T08T) * ZEZM )
    TH(I,J,K, 9) = QH1 * ( (R09M + R09T) * ZEXM &
    &                    + (S09M + S09T) * ZEYM &
    &                    + (T09M + T09T) * ZEZM )
    TH(I,J,K,10) = QH1 * ( (R10M + R10T) * ZEXM &
    &                    + (S10M + S10T) * ZEYM &
    &                    + (T10M + T10T) * ZEZM )
    TH(I,J,K,11) = QH1 * ( (R11M + R11T) * ZEXM &
    &                    + (S11M + S11T) * ZEYM &
    &                    + (T11M + T11T) * ZEZM )
    TH(I,J,K,12) = QH1 * ( (R12M + R12T) * ZEXM &
    &                    + (S12M + S12T) * ZEYM &
    &                    + (T12M + T12T) * ZEZM )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFZE
! 定義終了 *************************************************************
END SUBROUTINE Turbulence3DRSMLS
!***********************************************************************
!**** 乱流モデル : LES, 0Eq, Coherent-Structure Smagorinsky Model   ****
!****              (Kobayashi, 2005, Phys. Fluids, 17, 045104)      ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence2DLESCSM( &
&            RG, GAMMA, PR, PRT, &
&            OmegaZ, &
&            IS, IE, JS, JE, LS, LE, &
&            XIX, XIY, ETX, ETY, AJA, &
&            RHO, U, V, T, AMU, &
&            AMUT, DQD &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: C1 = 0.05, C2 = 0.045
  REAL, PARAMETER :: ZERO = 1.0E-20
  REAL, PARAMETER :: TwoThird = 0.66666667
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: RG, GAMMA, PR, PRT
  REAL,    INTENT(IN)  :: OmegaZ
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE), XIY(IS:IE, JS:JE), &
  &                       ETX(IS:IE, JS:JE), ETY(IS:IE, JS:JE), &
  &                       AJA(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE), &
  &                       U(IS:IE, JS:JE), V(IS:IE, JS:JE), &
  &                       T(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: AMUT(IS:IE, JS:JE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :), SH(:, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 4) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, LS: LE))
  ALLOCATE(SH(IS: IE, JS: JE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: UXI, VXI, UET, VET
  REAL    :: UX, VX, UY, VY
  REAL    :: S11, S22, S12, SS
  REAL    :: O12
  REAL    :: Q, E, FCS, FO, C, SS2
  ! 処理開始 ***********************************************************
  ! 渦粘性係数の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  UXI, VXI, UET, VET, UX, VX, UY, VY, &
  !$  S11, S22, S12, SS, O12, &
  !$  Q, E, FCS, FO, C, SS2 &
  !$)
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J) .GT. 0.0) THEN
    ! 速度勾配
    IF(I .EQ. IS)THEN
      IF(RHO(I+1,J) .GT. 0.0 .AND. RHO(I+2,J) .GT. 0.0) THEN
        UXI =-1.5 * U(I,J) + 2.0 * U(I+1,J) - 0.5 * U(I+2,J)
        VXI =-1.5 * V(I,J) + 2.0 * V(I+1,J) - 0.5 * V(I+2,J)
      ELSE
        UXI = 0.0
        VXI = 0.0
      ENDIF
    ELSE IF(I .EQ. IE)THEN
      IF(RHO(I-2,J) .GT. 0.0 .AND. RHO(I-1,J) .GT. 0.0) THEN
        UXI = 0.5 * U(I-2,J) - 2.0 * U(I-1,J) + 1.5 * U(I,J)
        VXI = 0.5 * V(I-2,J) - 2.0 * V(I-1,J) + 1.5 * V(I,J)
      ELSE
        UXI = 0.0
        VXI = 0.0
      ENDIF
    ELSE
      IF(RHO(I-1,J) .GT. 0.0 .AND. RHO(I+1,J) .GT. 0.0) THEN
        UXI =-0.5 * (U(I-1,J) - U(I+1,J))
        VXI =-0.5 * (V(I-1,J) - V(I+1,J))
      ELSE
        UXI = 0.0
        VXI = 0.0
      ENDIF
    ENDIF
    IF(J .EQ. JS)THEN
      IF(RHO(I,J+1) .GT. 0.0 .AND. RHO(I,J+2) .GT. 0.0) THEN
        UET =-1.5 * U(I,J) + 2.0 * U(I,J+1) - 0.5 * U(I,J+2)
        VET =-1.5 * V(I,J) + 2.0 * V(I,J+1) - 0.5 * V(I,J+2)
      ELSE
        UET = 0.0
        VET = 0.0
      ENDIF
    ELSE IF(J .EQ. JE)THEN
      IF(RHO(I,J-2) .GT. 0.0 .AND. RHO(I,J-1) .GT. 0.0) THEN
        UET = 0.5 * U(I,J-2) - 2.0 * U(I,J-1) + 1.5 * U(I,J)
        VET = 0.5 * V(I,J-2) - 2.0 * V(I,J-1) + 1.5 * V(I,J)
      ELSE
        UET = 0.0
        VET = 0.0
      ENDIF
    ELSE
      IF(RHO(I,J-1) .GT. 0.0 .AND. RHO(I,J+1) .GT. 0.0) THEN
        UET =-0.5 * (U(I,J-1) - U(I,J+1))
        VET =-0.5 * (V(I,J-1) - V(I,J+1))
      ELSE
        UET = 0.0
        VET = 0.0
      ENDIF
    ENDIF
    UX = XIX(I,J) * UXI + ETX(I,J) * UET
    UY = XIY(I,J) * UXI + ETY(I,J) * UET
    VX = XIX(I,J) * VXI + ETX(I,J) * VET
    VY = XIY(I,J) * VXI + ETY(I,J) * VET
    ! 歪み速度テンソル, 角速度テンソル
    SS  = (UX + VY) / 3.0
    S11 = UX - SS
    S22 = VY - SS
    S12 = 0.5 * (UY + VX)
    O12 = 0.5 * (UY - VX) - OmegaZ
    ! 速度勾配の第二不変量
    Q   = O12 * O12 - S12 * S12 &
    &   - 0.5 * (S11 * S11 + S22 * S22)
    ! 速度勾配の大きさ
    E   = O12 * O12 + S12 * S12 &
    &   + 0.5 * (S11 * S11 + S22 * S22)
    ! 係数
    IF(E .GT. ZERO) THEN
      FCS = Q / E
    ELSE
      FCS = 1.0
    ENDIF
    FO  = 1.0 - FCS
    C   = C2 * (ABS(FCS))**1.5 * FO
    SS2 = S11 * S11 + S22 * S22 + 2.0 * S12 * S12
    SS2 = SQRT(2.0 * SS2)
    AMUT(I,J) = RHO(I,J) * 2.0 * C * SS2 / AJA(I,J)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH = 0.0
  SH = 0.0
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL DIFFXI
  CALL DIFFET
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, L)
  DO L = LS, LE
  !$OMP DO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J-1) .GT. 0.0) .AND. (RHO(I-1,J) .GT. 0.0) .AND. &
  &   (RHO(I,J)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,L) = (-RH(I-1,J  ,L) + RH(I,J,L)) &
    &          + (-SH(I  ,J-1,L) + SH(I,J,L))
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項の計算(xi方向差分用)                                    ****
!***********************************************************************
SUBROUTINE DIFFXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: RHOM, UM, VM, TM, AMUM, AMUTM
  REAL    :: UXI, UET, UX, UY
  REAL    :: VXI, VET, VX, VY
  REAL    :: TXI, TET, TX, TY
  REAL    :: nu, nu_t
  REAL    :: c_u, c_t
  REAL    :: TAUXX, TAUYY, TAUXY, DELV
  REAL    :: R4, S4
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, &
  !$  RHOM, UM, VM, TM, AMUM, AMUTM, &
  !$  UXI, UET, UX, UY, &
  !$  VXI, VET, VX, VY, &
  !$  TXI, TET, TX, TY, &
  !$  nu, nu_t, c_u, c_t, &
  !$  TAUXX, TAUYY, TAUXY, DELV, R4, S4, QH1 &
  !$)
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF( (RHO(I,J-1) .GT. 0.0) .AND. (RHO(I,J)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0)  &
  & ) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J) + XIX(I+1,J))
    XIYM = 0.5 * (XIY(I,J) + XIY(I+1,J))
    ETXM = 0.5 * (ETX(I,J) + ETX(I+1,J))
    ETYM = 0.5 * (ETY(I,J) + ETY(I+1,J))
    AJAM = 0.5 * (AJA(I,J) + AJA(I+1,J))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J) +  RHO(I+1,J))
    UM    = 0.5 * (   U(I,J) +    U(I+1,J))
    VM    = 0.5 * (   V(I,J) +    V(I+1,J))
    TM    = 0.5 * (   T(I,J) +    T(I+1,J))
    AMUM  = 0.5 * ( AMU(I,J) +  AMU(I+1,J))
    AMUTM = 0.5 * (AMUT(I,J) + AMUT(I+1,J))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = - U(I,J) + U(I+1,J)
    VXI = - V(I,J) + V(I+1,J)
    UET = 0.25 * (- U(I,J-1) - U(I+1,J-1) + U(I,J+1) + U(I+1,J+1))
    VET = 0.25 * (- V(I,J-1) - V(I+1,J-1) + V(I,J+1) + V(I+1,J+1))
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM
    UY = UXI * XIYM + UET * ETYM
    VX = VXI * XIXM + VET * ETXM
    VY = VXI * XIYM + VET * ETYM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = - T(I,J) + T(I+1,J)
    TET = 0.25 * (- T(I,J-1) - T(I+1,J-1) + T(I,J+1) + T(I+1,J+1))
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM
    TY = TXI * XIYM + TET * ETYM
    ! 係数計算 ---------------------------------------------------------
    nu   = AMUM  / RHOM
    nu_t = AMUTM / RHOM
    c_u  = nu + nu_t
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    DELV  = UX + VY
    TAUXX = c_u * (2.0 * UX - TwoThird * DELV) ! - TwoThird * AKM
    TAUYY = c_u * (2.0 * VY - TwoThird * DELV) ! - TwoThird * AKM
    TAUXY = c_u * (UY + VX)
    ! エネルギーの拡散 -------------------------------------------------
    R4 = TAUXX * UM + TAUXY * VM + c_t * TX
    S4 = TAUXY * UM + TAUYY * VM + c_t * TY
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    RH(I,J,1) = 0.0
    RH(I,J,2) = QH1 * (TAUXX * XIXM + TAUXY * XIYM)
    RH(I,J,3) = QH1 * (TAUXY * XIXM + TAUYY * XIYM)
    RH(I,J,4) = QH1 * (R4 * XIXM + S4 * XIYM)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFXI
!***********************************************************************
!**** 拡散項の計算(eta方向差分用)                                   ****
!***********************************************************************
SUBROUTINE DIFFET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: XIXM, XIYM, ETXM, ETYM, AJAM
  REAL    :: RHOM, UM, VM, TM, AMUM, AMUTM
  REAL    :: UXI, UET, UX, UY
  REAL    :: VXI, VET, VX, VY
  REAL    :: TXI, TET, TX, TY
  REAL    :: nu, nu_t
  REAL    :: c_u, c_t
  REAL    :: TAUXX, TAUYY, TAUXY, DELV
  REAL    :: R4, S4
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, &
  !$  XIXM, XIYM, ETXM, ETYM, AJAM, &
  !$  RHOM, UM, VM, TM, AMUM, AMUTM, &
  !$  UXI, UET, UX, UY, &
  !$  VXI, VET, VX, VY, &
  !$  TXI, TET, TX, TY, &
  !$  nu, nu_t, c_u, c_t, &
  !$  TAUXX, TAUYY, TAUXY, DELV, R4, S4, QH1 &
  !$)
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I-1,J) .GT. 0.0) .AND. (RHO(I,J)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J) .GT. 0.0) .AND. (RHO(I,J+1) .GT. 0.0) &
  & ) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J) + XIX(I,J+1))
    XIYM = 0.5 * (XIY(I,J) + XIY(I,J+1))
    ETXM = 0.5 * (ETX(I,J) + ETX(I,J+1))
    ETYM = 0.5 * (ETY(I,J) + ETY(I,J+1))
    AJAM = 0.5 * (AJA(I,J) + AJA(I,J+1))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J) +  RHO(I,J+1))
    UM    = 0.5 * (   U(I,J) +    U(I,J+1))
    VM    = 0.5 * (   V(I,J) +    V(I,J+1))
    TM    = 0.5 * (   T(I,J) +    T(I,J+1))
    AMUM  = 0.5 * ( AMU(I,J) +  AMU(I,J+1))
    AMUTM = 0.5 * (AMUT(I,J) + AMUT(I,J+1))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * (- U(I-1,J) - U(I-1,J+1) + U(I+1,J) + U(I+1,J+1))
    VXI = 0.25 * (- V(I-1,J) - V(I-1,J+1) + V(I+1,J) + V(I+1,J+1))
    UET = - U(I,J) + U(I,J+1)
    VET = - V(I,J) + V(I,J+1)
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM
    UY = UXI * XIYM + UET * ETYM
    VX = VXI * XIXM + VET * ETXM
    VY = VXI * XIYM + VET * ETYM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * (- T(I-1,J) - T(I-1,J+1) + T(I+1,J) + T(I+1,J+1))
    TET = - T(I,J) + T(I,J+1)
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM
    TY = TXI * XIYM + TET * ETYM
    ! 係数計算 ---------------------------------------------------------
    nu   = AMUM  / RHOM
    nu_t = AMUTM / RHOM
    c_u  = nu + nu_t
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    DELV  = UX + VY
    TAUXX = c_u * (2.0 * UX - TwoThird * DELV) ! - TwoThird * AKM
    TAUYY = c_u * (2.0 * VY - TwoThird * DELV) ! - TwoThird * AKM
    TAUXY = c_u * (UY + VX)
    ! エネルギーの拡散 -------------------------------------------------
    R4 = TAUXX * UM + TAUXY * VM + c_t * TX
    S4 = TAUXY * UM + TAUYY * VM + c_t * TY
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    SH(I,J,1) = 0.0
    SH(I,J,2) = QH1 * (TAUXX * ETXM + TAUXY * ETYM)
    SH(I,J,3) = QH1 * (TAUXY * ETXM + TAUYY * ETYM)
    SH(I,J,4) = QH1 * (R4 * ETXM + S4 * ETYM)
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFET
! 定義終了 *************************************************************
END SUBROUTINE Turbulence2DLESCSM
!***********************************************************************
!**** 乱流モデル : LES, 0Eq, Coherent-Structure Smagorinsky Model   ****
!****              (Kobayashi, 2005, Phys. Fluids, 17, 045104)      ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!**** 計算スキーム : 二次精度中心差分(スタガード格子上の点を使用)   ****
!***********************************************************************
SUBROUTINE Turbulence3DLESCSM( &
&            RG, GAMMA, PR, PRT, &
&            OmegaX, OmegaY, OmegaZ, &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            XIX, XIY, XIZ, ETX, ETY, ETZ, ZEX, ZEY, ZEZ, AJA, &
&            RHO, U, V, W, T, AMU, &
&            AMUT, DQD &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: C1 = 0.05, C2 = 0.045
  REAL, PARAMETER :: ZERO = 1.0E-20
  REAL, PARAMETER :: TwoThird = 0.66666667
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: RG, GAMMA, PR, PRT
  REAL,    INTENT(IN)  :: OmegaX, OmegaY, OmegaZ
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: XIX(IS:IE, JS:JE, KS:KE), &
  &                       XIY(IS:IE, JS:JE, KS:KE), &
  &                       XIZ(IS:IE, JS:JE, KS:KE), &
  &                       ETX(IS:IE, JS:JE, KS:KE), &
  &                       ETY(IS:IE, JS:JE, KS:KE), &
  &                       ETZ(IS:IE, JS:JE, KS:KE), &
  &                       ZEX(IS:IE, JS:JE, KS:KE), &
  &                       ZEY(IS:IE, JS:JE, KS:KE), &
  &                       ZEZ(IS:IE, JS:JE, KS:KE), &
  &                       AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: RHO(IS:IE, JS:JE, KS:KE), &
  &                       U(IS:IE, JS:JE, KS:KE), &
  &                       V(IS:IE, JS:JE, KS:KE), &
  &                       W(IS:IE, JS:JE, KS:KE), &
  &                       T(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: AMU(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: AMUT(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(OUT) :: DQD(IS:IE, JS:JE, KS:KE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, ALLOCATABLE :: RH(:, :, :, :), SH(:, :, :, :), TH(:, :, :, :)
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 5) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! メモリ確保 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ALLOCATE(RH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(SH(IS: IE, JS: JE, KS: KE, LS: LE))
  ALLOCATE(TH(IS: IE, JS: JE, KS: KE, LS: LE))
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AMUT = 0.0
  DQD  = 0.0
  ! 内部手続きの呼び出し +++++++++++++++++++++++++++++++++++++++++++++++
  CALL NonDiffTerm
  CALL DiffTerm
  ! 処理終了 ***********************************************************
  RETURN
! 内部手続き ***********************************************************
CONTAINS
!***********************************************************************
!**** 非拡散項の計算                                                ****
!***********************************************************************
SUBROUTINE NonDiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE
  REAL    :: UX, VX, WX, UY, VY, WY, UZ, VZ, WZ
  REAL    :: S11, S22, S33, S12, S23, S31, SS
  REAL    :: O12, O23, O31
  REAL    :: Q, E, FCS, FO, C, SS2
  ! 処理開始 ***********************************************************
  ! 渦粘性係数の計算 +++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  UXI, VXI, WXI, UET, VET, WET, UZE, VZE, WZE, &
  !$  UX, VX, WX, UY, VY, WY, UZ, VZ, WZ, &
  !$  S11, S22, S33, S12, S23, S31, SS, &
  !$  O12, O23, O31, &
  !$  Q, E, FCS, FO, C, SS2 &
  !$)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J,K) .GT. 0.0) THEN
    ! 速度勾配
    IF(I .EQ. IS)THEN
      IF(RHO(I+1,J,K) .GT. 0.0 .AND. RHO(I+2,J,K) .GT. 0.0) THEN
        UXI =-1.5 * U(I,J,K) + 2.0 * U(I+1,J,K) - 0.5 * U(I+2,J,K)
        VXI =-1.5 * V(I,J,K) + 2.0 * V(I+1,J,K) - 0.5 * V(I+2,J,K)
        WXI =-1.5 * W(I,J,K) + 2.0 * W(I+1,J,K) - 0.5 * W(I+2,J,K)
      ELSE
        UXI = 0.0
        VXI = 0.0
        WXI = 0.0
      ENDIF
    ELSE IF(I .EQ. IE)THEN
      IF(RHO(I-2,J,K) .GT. 0.0 .AND. RHO(I-1,J,K) .GT. 0.0) THEN
        UXI = 0.5 * U(I-2,J,K) - 2.0 * U(I-1,J,K) + 1.5 * U(I,J,K)
        VXI = 0.5 * V(I-2,J,K) - 2.0 * V(I-1,J,K) + 1.5 * V(I,J,K)
        WXI = 0.5 * W(I-2,J,K) - 2.0 * W(I-1,J,K) + 1.5 * W(I,J,K)
      ELSE
        UXI = 0.0
        VXI = 0.0
        WXI = 0.0
      ENDIF
    ELSE
      IF(RHO(I-1,J,K) .GT. 0.0 .AND. RHO(I+1,J,K) .GT. 0.0) THEN
        UXI =-0.5 * (U(I-1,J,K) - U(I+1,J,K))
        VXI =-0.5 * (V(I-1,J,K) - V(I+1,J,K))
        WXI =-0.5 * (W(I-1,J,K) - W(I+1,J,K))
      ELSE
        UXI = 0.0
        VXI = 0.0
        WXI = 0.0
      ENDIF
    ENDIF
    IF(J .EQ. JS)THEN
      IF(RHO(I,J+1,K) .GT. 0.0 .AND. RHO(I,J+2,K) .GT. 0.0) THEN
        UET =-1.5 * U(I,J,K) + 2.0 * U(I,J+1,K) - 0.5 * U(I,J+2,K)
        VET =-1.5 * V(I,J,K) + 2.0 * V(I,J+1,K) - 0.5 * V(I,J+2,K)
        WET =-1.5 * W(I,J,K) + 2.0 * W(I,J+1,K) - 0.5 * W(I,J+2,K)
      ELSE
        UET = 0.0
        VET = 0.0
        WET = 0.0
      ENDIF
    ELSE IF(J .EQ. JE)THEN
      IF(RHO(I,J-2,K) .GT. 0.0 .AND. RHO(I,J-1,K) .GT. 0.0) THEN
        UET = 0.5 * U(I,J-2,K) - 2.0 * U(I,J-1,K) + 1.5 * U(I,J,K)
        VET = 0.5 * V(I,J-2,K) - 2.0 * V(I,J-1,K) + 1.5 * V(I,J,K)
        WET = 0.5 * W(I,J-2,K) - 2.0 * W(I,J-1,K) + 1.5 * W(I,J,K)
      ELSE
        UET = 0.0
        VET = 0.0
        WET = 0.0
      ENDIF
    ELSE
      IF(RHO(I,J-1,K) .GT. 0.0 .AND. RHO(I,J+1,K) .GT. 0.0) THEN
        UET =-0.5 * (U(I,J-1,K) - U(I,J+1,K))
        VET =-0.5 * (V(I,J-1,K) - V(I,J+1,K))
        WET =-0.5 * (W(I,J-1,K) - W(I,J+1,K))
      ELSE
        UET = 0.0
        VET = 0.0
        WET = 0.0
      ENDIF
    ENDIF
    IF(K .EQ. KS)THEN
      IF(RHO(I,J,K+1) .GT. 0.0 .AND. RHO(I,J,K+2) .GT. 0.0) THEN
        UZE =-1.5 * U(I,J,K) + 2.0 * U(I,J,K+1) - 0.5 * U(I,J,K+2)
        VZE =-1.5 * V(I,J,K) + 2.0 * V(I,J,K+1) - 0.5 * V(I,J,K+2)
        WZE =-1.5 * W(I,J,K) + 2.0 * W(I,J,K+1) - 0.5 * W(I,J,K+2)
      ELSE
        UZE = 0.0
        VZE = 0.0
        WZE = 0.0
      ENDIF
    ELSE IF(K .EQ. KE)THEN
      IF(RHO(I,J,K-2) .GT. 0.0 .AND. RHO(I,J,K-1) .GT. 0.0) THEN
        UZE = 0.5 * U(I,J,K-2) - 2.0 * U(I,J,K-1) + 1.5 * U(I,J,K)
        VZE = 0.5 * V(I,J,K-2) - 2.0 * V(I,J,K-1) + 1.5 * V(I,J,K)
        WZE = 0.5 * W(I,J,K-2) - 2.0 * W(I,J,K-1) + 1.5 * W(I,J,K)
      ELSE
        UZE = 0.0
        VZE = 0.0
        WZE = 0.0
      ENDIF
    ELSE
      IF(RHO(I,J,K-1) .GT. 0.0 .AND. RHO(I,J,K+1) .GT. 0.0) THEN
        UZE =-0.5 * (U(I,J,K-1) - U(I,J,K+1))
        VZE =-0.5 * (V(I,J,K-1) - V(I,J,K+1))
        WZE =-0.5 * (W(I,J,K-1) - W(I,J,K+1))
      ELSE
        UZE = 0.0
        VZE = 0.0
        WZE = 0.0
      ENDIF
    ENDIF
    UX = XIX(I,J,K) * UXI + ETX(I,J,K) * UET + ZEX(I,J,K) * UZE
    UY = XIY(I,J,K) * UXI + ETY(I,J,K) * UET + ZEY(I,J,K) * UZE
    UZ = XIZ(I,J,K) * UXI + ETZ(I,J,K) * UET + ZEZ(I,J,K) * UZE
    VX = XIX(I,J,K) * VXI + ETX(I,J,K) * VET + ZEX(I,J,K) * VZE
    VY = XIY(I,J,K) * VXI + ETY(I,J,K) * VET + ZEY(I,J,K) * VZE
    VZ = XIZ(I,J,K) * VXI + ETZ(I,J,K) * VET + ZEZ(I,J,K) * VZE
    WX = XIX(I,J,K) * WXI + ETX(I,J,K) * WET + ZEX(I,J,K) * WZE
    WY = XIY(I,J,K) * WXI + ETY(I,J,K) * WET + ZEY(I,J,K) * WZE
    WZ = XIZ(I,J,K) * WXI + ETZ(I,J,K) * WET + ZEZ(I,J,K) * WZE
    ! 歪み速度テンソル, 角速度テンソル
    SS  = (UX + VY + WZ) / 3.0
    S11 = UX - SS
    S22 = VY - SS
    S33 = WZ - SS
    S12 = 0.5 * (UY + VX)
    S23 = 0.5 * (VZ + WY)
    S31 = 0.5 * (WX + UZ)
    O12 = 0.5 * (UY - VX) - OmegaZ
    O23 = 0.5 * (VZ - WY) - OmegaX
    O31 = 0.5 * (WX - UZ) - OmegaY
    ! 速度勾配の第二不変量
    Q   = O12 * O12 - S12 * S12 &
    &   + O23 * O23 - S23 * S23 &
    &   + O31 * O31 - S31 * S31 &
    &   - 0.5 * (S11 * S11 + S22 * S22 + S33 * S33)
    ! 速度勾配の大きさ
    E   = O12 * O12 + S12 * S12 &
    &   + O23 * O23 + S23 * S23 &
    &   + O31 * O31 + S31 * S31 &
    &   + 0.5 * (S11 * S11 + S22 * S22 + S33 * S33)
    ! 係数
    IF(E .GT. ZERO) THEN
      FCS = Q / E
    ELSE
      FCS = 1.0
    ENDIF
    FO  = 1.0 - FCS
    C   = C2 * (ABS(FCS))**1.5 * FO
    SS2 = S11 * S11 + S22 * S22 + S33 * S33 &
    &   + 2.0 * (S12 * S12 + S23 * S23 + S31 * S31)
    SS2 = SQRT(2.0 * SS2)
    AMUT(I,J,K) = RHO(I,J,K) * 2.0 * C * SS2 / AJA(I,J,K)**(2.0 / 3.0)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE NonDiffTerm
!***********************************************************************
!**** 拡散項の計算                                                  ****
!***********************************************************************
SUBROUTINE DiffTerm
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  ! 処理開始 ***********************************************************
  ! 初期化 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RH = 0.0
  SH = 0.0
  TH = 0.0
  ! 拡散項の要素を計算 +++++++++++++++++++++++++++++++++++++++++++++++++
  CALL DIFFXI
  CALL DIFFET
  CALL DIFFZE
  ! 拡散項の差分 +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I,J-1,K) .GT. 0.0) .AND. &
  &   (RHO(I-1,J,K) .GT. 0.0) .AND. (RHO(I,J,K)   .GT. 0.0) .AND. &
  &   (RHO(I+1,J,K) .GT. 0.0) .AND. (RHO(I,J+1,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    DQD(I,J,K,L) = (-RH(I-1,J  ,K  ,L) + RH(I,J,K,L) ) &
    &            + (-SH(I  ,J-1,K  ,L) + SH(I,J,K,L) ) &
    &            + (-TH(I  ,J  ,K-1,L) + TH(I,J,K,L) )
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DiffTerm
!***********************************************************************
!**** 拡散項の計算(xi方向差分用)                                    ****
!***********************************************************************
SUBROUTINE DIFFXI
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, TM, AMUM, AMUTM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: nu, nu_t
  REAL    :: c_u, c_t
  REAL    :: TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV
  REAL    :: R5, S5, T5
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, TM, AMUM, AMUTM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  nu, nu_t, c_u, c_t, &
  !$  TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV, R5, S5, T5, QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS    , IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I,J-1,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K)   .GT. 0.0) .AND. (RHO(I+1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J+1,K) .GT. 0.0) .AND. (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    ! XI方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I+1,J,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I+1,J,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I+1,J,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I+1,J,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I+1,J,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I+1,J,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I+1,J,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I+1,J,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I+1,J,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I+1,J,K))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J,K) +  RHO(I+1,J,K))
    UM    = 0.5 * (   U(I,J,K) +    U(I+1,J,K))
    VM    = 0.5 * (   V(I,J,K) +    V(I+1,J,K))
    WM    = 0.5 * (   W(I,J,K) +    W(I+1,J,K))
    TM    = 0.5 * (   T(I,J,K) +    T(I+1,J,K))
    AMUM  = 0.5 * ( AMU(I,J,K) +  AMU(I+1,J,K))
    AMUTM = 0.5 * (AMUT(I,J,K) + AMUT(I+1,J,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = - U(I,J,K) + U(I+1,J,K)
    VXI = - V(I,J,K) + V(I+1,J,K)
    WXI = - W(I,J,K) + W(I+1,J,K)
    UET = 0.25 * ( - U(I,J-1,K) - U(I+1,J-1,K) &
    &              + U(I,J+1,K) + U(I+1,J+1,K) )
    VET = 0.25 * ( - V(I,J-1,K) - V(I+1,J-1,K) &
    &              + V(I,J+1,K) + V(I+1,J+1,K) )
    WET = 0.25 * ( - W(I,J-1,K) - W(I+1,J-1,K) &
    &              + W(I,J+1,K) + W(I+1,J+1,K) )
    UZE = 0.25 * ( - U(I,J,K-1) - U(I+1,J,K-1) &
    &              + U(I,J,K+1) + U(I+1,J,K+1) )
    VZE = 0.25 * ( - V(I,J,K-1) - V(I+1,J,K-1) &
    &              + V(I,J,K+1) + V(I+1,J,K+1) )
    WZE = 0.25 * ( - W(I,J,K-1) - W(I+1,J,K-1) &
    &              + W(I,J,K+1) + W(I+1,J,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = - T(I,J,K) + T(I+1,J,K)
    TET = 0.25 * ( - T(I,J-1,K) - T(I+1,J-1,K) &
    &              + T(I,J+1,K) + T(I+1,J+1,K) )
    TZE = 0.25 * ( - T(I,J,K-1) - T(I+1,J,K-1) &
    &              + T(I,J,K+1) + T(I+1,J,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! 係数計算 ---------------------------------------------------------
    nu   = AMUM  / RHOM
    nu_t = AMUTM / RHOM
    c_u  = nu + nu_t
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    DELV  = UX + VY + WZ
    TAUXX = c_u * (2.0 * UX - TwoThird * DELV) ! - TwoThird * AKM
    TAUYY = c_u * (2.0 * VY - TwoThird * DELV) ! - TwoThird * AKM
    TAUZZ = c_u * (2.0 * WZ - TwoThird * DELV) ! - TwoThird * AKM
    TAUXY = c_u * (UY + VX)
    TAUYZ = c_u * (VZ + WY)
    TAUZX = c_u * (WX + UZ)
    ! エネルギーの拡散 -------------------------------------------------
    R5 = TAUXX * UM + TAUXY * VM + TAUZX * WM + c_t * TX
    S5 = TAUXY * UM + TAUYY * VM + TAUYZ * WM + c_t * TY
    T5 = TAUZX * UM + TAUYZ * VM + TAUZZ * WM + c_t * TZ
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    RH(I,J,K,1) = 0.0
    RH(I,J,K,2) = QH1 * (TAUXX * XIXM + TAUXY * XIYM + TAUZX * XIZM)
    RH(I,J,K,3) = QH1 * (TAUXY * XIXM + TAUYY * XIYM + TAUYZ * XIZM)
    RH(I,J,K,4) = QH1 * (TAUZX * XIXM + TAUYZ * XIYM + TAUZZ * XIZM)
    RH(I,J,K,5) = QH1 * (R5 * XIXM + S5 * XIYM + T5 * XIZM)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFXI
!***********************************************************************
!**** 拡散項の計算(eta方向差分用)                                   ****
!***********************************************************************
SUBROUTINE DIFFET
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, TM, AMUM, AMUTM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: nu, nu_t
  REAL    :: c_u, c_t
  REAL    :: TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV
  REAL    :: R5, S5, T5
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, TM, AMUM, AMUTM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  nu, nu_t, c_u, c_t, &
  !$  TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV, R5, S5, T5, QH1 &
  !$)
  DO K = KS + 1, KE - 1
  DO J = JS    , JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J,K-1) .GT. 0.0) .AND. (RHO(I-1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K)   .GT. 0.0) .AND. (RHO(I+1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J+1,K) .GT. 0.0) .AND. (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    ! ET方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J+1,K))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J+1,K))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J+1,K))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J+1,K))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J+1,K))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J+1,K))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J+1,K))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J+1,K))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J+1,K))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J+1,K))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J,K) +  RHO(I,J+1,K))
    UM    = 0.5 * (   U(I,J,K) +    U(I,J+1,K))
    VM    = 0.5 * (   V(I,J,K) +    V(I,J+1,K))
    WM    = 0.5 * (   W(I,J,K) +    W(I,J+1,K))
    TM    = 0.5 * (   T(I,J,K) +    T(I,J+1,K))
    AMUM  = 0.5 * ( AMU(I,J,K) +  AMU(I,J+1,K))
    AMUTM = 0.5 * (AMUT(I,J,K) + AMUT(I,J+1,K))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J+1,K) &
    &            + U(I+1,J,K) + U(I+1,J+1,K) )
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J+1,K) &
    &            + V(I+1,J,K) + V(I+1,J+1,K) )
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J+1,K) &
    &            + W(I+1,J,K) + W(I+1,J+1,K) )
    UET = - U(I,J,K) + U(I,J+1,K)
    VET = - V(I,J,K) + V(I,J+1,K)
    WET = - W(I,J,K) + W(I,J+1,K)
    UZE = 0.25 * ( - U(I,J,K-1) - U(I,J+1,K-1) &
    &              + U(I,J,K+1) + U(I,J+1,K+1) )
    VZE = 0.25 * ( - V(I,J,K-1) - V(I,J+1,K-1) &
    &              + V(I,J,K+1) + V(I,J+1,K+1) )
    WZE = 0.25 * ( - W(I,J,K-1) - W(I,J+1,K-1) &
    &              + W(I,J,K+1) + W(I,J+1,K+1) )
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J+1,K) &
    &              + T(I+1,J,K) + T(I+1,J+1,K) )
    TET = - T(I,J,K) + T(I,J+1,K)
    TZE = 0.25 * ( - T(I,J,K-1) - T(I,J+1,K-1) &
    &              + T(I,J,K+1) + T(I,J+1,K+1) )
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! 係数計算 ---------------------------------------------------------
    nu   = AMUM  / RHOM
    nu_t = AMUTM / RHOM
    c_u  = nu + nu_t
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    DELV  = UX + VY + WZ
    TAUXX = c_u * (2.0 * UX - TwoThird * DELV) ! - TwoThird * AKM
    TAUYY = c_u * (2.0 * VY - TwoThird * DELV) ! - TwoThird * AKM
    TAUZZ = c_u * (2.0 * WZ - TwoThird * DELV) ! - TwoThird * AKM
    TAUXY = c_u * (UY + VX)
    TAUYZ = c_u * (VZ + WY)
    TAUZX = c_u * (WX + UZ)
    ! エネルギーの拡散 -------------------------------------------------
    R5 = TAUXX * UM + TAUXY * VM + TAUZX * WM + c_t * TX
    S5 = TAUXY * UM + TAUYY * VM + TAUYZ * WM + c_t * TY
    T5 = TAUZX * UM + TAUYZ * VM + TAUZZ * WM + c_t * TZ
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    SH(I,J,K,1) = 0.0
    SH(I,J,K,2) = QH1 * (TAUXX * ETXM + TAUXY * ETYM + TAUZX * ETZM)
    SH(I,J,K,3) = QH1 * (TAUXY * ETXM + TAUYY * ETYM + TAUYZ * ETZM)
    SH(I,J,K,4) = QH1 * (TAUZX * ETXM + TAUYZ * ETYM + TAUZZ * ETZM)
    SH(I,J,K,5) = QH1 * (R5 * ETXM + S5 * ETYM + T5 * ETZM)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFET
!***********************************************************************
!**** 拡散項の計算(zeta方向差分用)                                  ****
!***********************************************************************
SUBROUTINE DIFFZE
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM
  REAL    :: RHOM, UM, VM, WM, TM, AMUM, AMUTM
  REAL    :: UXI, UET, UZE, UX, UY, UZ
  REAL    :: VXI, VET, VZE, VX, VY, VZ
  REAL    :: WXI, WET, WZE, WX, WY, WZ
  REAL    :: TXI, TET, TZE, TX, TY, TZ
  REAL    :: nu, nu_t
  REAL    :: c_u, c_t
  REAL    :: TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV
  REAL    :: R5, S5, T5
  REAL    :: QH1
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
  !$  I, J, K, &
  !$  XIXM, XIYM, XIZM, ETXM, ETYM, ETZM, ZEXM, ZEYM, ZEZM, AJAM, &
  !$  RHOM, UM, VM, WM, TM, AMUM, AMUTM, &
  !$  UXI, UET, UZE, UX, UY, UZ, &
  !$  VXI, VET, VZE, VX, VY, VZ, &
  !$  WXI, WET, WZE, WX, WY, WZ, &
  !$  TXI, TET, TZE, TX, TY, TZ, &
  !$  nu, nu_t, c_u, c_t, &
  !$  TAUXX, TAUYY, TAUZZ, TAUXY, TAUYZ, TAUZX, DELV, R5, S5, T5, QH1 &
  !$)
  DO K = KS    , KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (RHO(I,J-1,K) .GT. 0.0) .AND. (RHO(I-1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J,K)   .GT. 0.0) .AND. (RHO(I+1,J,K) .GT. 0.0) .AND. &
  &   (RHO(I,J+1,K) .GT. 0.0) .AND. (RHO(I,J,K+1) .GT. 0.0) &
  & ) THEN
    ! ZE方向の平均量の計算 ---------------------------------------------
    ! 座標変換パラメータ
    XIXM = 0.5 * (XIX(I,J,K) + XIX(I,J,K+1))
    XIYM = 0.5 * (XIY(I,J,K) + XIY(I,J,K+1))
    XIZM = 0.5 * (XIZ(I,J,K) + XIZ(I,J,K+1))
    ETXM = 0.5 * (ETX(I,J,K) + ETX(I,J,K+1))
    ETYM = 0.5 * (ETY(I,J,K) + ETY(I,J,K+1))
    ETZM = 0.5 * (ETZ(I,J,K) + ETZ(I,J,K+1))
    ZEXM = 0.5 * (ZEX(I,J,K) + ZEX(I,J,K+1))
    ZEYM = 0.5 * (ZEY(I,J,K) + ZEY(I,J,K+1))
    ZEZM = 0.5 * (ZEZ(I,J,K) + ZEZ(I,J,K+1))
    AJAM = 0.5 * (AJA(I,J,K) + AJA(I,J,K+1))
    ! 物理量
    RHOM  = 0.5 * ( RHO(I,J,K) +  RHO(I,J,K+1))
    UM    = 0.5 * (   U(I,J,K) +    U(I,J,K+1))
    VM    = 0.5 * (   V(I,J,K) +    V(I,J,K+1))
    WM    = 0.5 * (   W(I,J,K) +    W(I,J,K+1))
    TM    = 0.5 * (   T(I,J,K) +    T(I,J,K+1))
    AMUM  = 0.5 * ( AMU(I,J,K) +  AMU(I,J,K+1))
    AMUTM = 0.5 * (AMUT(I,J,K) + AMUT(I,J,K+1))
    ! 速度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    UXI = 0.25 * ( - U(I-1,J,K) - U(I-1,J,K+1) &
    &              + U(I+1,J,K) + U(I+1,J,K+1) )
    VXI = 0.25 * ( - V(I-1,J,K) - V(I-1,J,K+1) &
    &              + V(I+1,J,K) + V(I+1,J,K+1) )
    WXI = 0.25 * ( - W(I-1,J,K) - W(I-1,J,K+1) &
    &              + W(I+1,J,K) + W(I+1,J,K+1) )
    UET = 0.25 * ( - U(I,J-1,K) - U(I,J-1,K+1) &
    &              + U(I,J+1,K) + U(I,J+1,K+1) )
    VET = 0.25 * ( - V(I,J-1,K) - V(I,J-1,K+1) &
    &              + V(I,J+1,K) + V(I,J+1,K+1) )
    WET = 0.25 * ( - W(I,J-1,K) - W(I,J-1,K+1) &
    &              + W(I,J+1,K) + W(I,J+1,K+1) )
    UZE = - U(I,J,K) + U(I,J,K+1)
    VZE = - V(I,J,K) + V(I,J,K+1)
    WZE = - W(I,J,K) + W(I,J,K+1)
    ! 物理空間方向一階微分
    UX = UXI * XIXM + UET * ETXM + UZE * ZEXM
    UY = UXI * XIYM + UET * ETYM + UZE * ZEYM
    UZ = UXI * XIZM + UET * ETZM + UZE * ZEZM
    VX = VXI * XIXM + VET * ETXM + VZE * ZEXM
    VY = VXI * XIYM + VET * ETYM + VZE * ZEYM
    VZ = VXI * XIZM + VET * ETZM + VZE * ZEZM
    WX = WXI * XIXM + WET * ETXM + WZE * ZEXM
    WY = WXI * XIYM + WET * ETYM + WZE * ZEYM
    WZ = WXI * XIZM + WET * ETZM + WZE * ZEZM
    ! 温度の一階微分の計算 ---------------------------------------------
    ! 計算空間方向一階微分
    TXI = 0.25 * ( - T(I-1,J,K) - T(I-1,J,K+1) &
    &              + T(I+1,J,K) + T(I+1,J,K+1) )
    TET = 0.25 * ( - T(I,J-1,K) - T(I,J-1,K+1) &
    &              + T(I,J+1,K) + T(I,J+1,K+1) )
    TZE = - T(I,J,K) + T(I,J,K+1)
    ! 物理空間方向一階微分
    TX = TXI * XIXM + TET * ETXM + TZE * ZEXM
    TY = TXI * XIYM + TET * ETYM + TZE * ZEYM
    TZ = TXI * XIZM + TET * ETZM + TZE * ZEZM
    ! 係数計算 ---------------------------------------------------------
    nu   = AMUM  / RHOM
    nu_t = AMUTM / RHOM
    c_u  = nu + nu_t
    c_t  = (nu / PR + nu_t / PRT) / (GAMMA - 1.0) * GAMMA * RG
    ! 粘性応力とレイノルズ応力の和 -------------------------------------
    DELV  = UX + VY + WZ
    TAUXX = c_u * (2.0 * UX - TwoThird * DELV) ! - TwoThird * AKM
    TAUYY = c_u * (2.0 * VY - TwoThird * DELV) ! - TwoThird * AKM
    TAUZZ = c_u * (2.0 * WZ - TwoThird * DELV) ! - TwoThird * AKM
    TAUXY = c_u * (UY + VX)
    TAUYZ = c_u * (VZ + WY)
    TAUZX = c_u * (WX + UZ)
    ! エネルギーの拡散 -------------------------------------------------
    R5 = TAUXX * UM + TAUXY * VM + TAUZX * WM + c_t * TX
    S5 = TAUXY * UM + TAUYY * VM + TAUYZ * WM + c_t * TY
    T5 = TAUZX * UM + TAUYZ * VM + TAUZZ * WM + c_t * TZ
    ! 各方程式の拡散項 -------------------------------------------------
    QH1 = RHOM / AJAM
    TH(I,J,K,1) = 0.0
    TH(I,J,K,2) = QH1 * (TAUXX * ZEXM + TAUXY * ZEYM + TAUZX * ZEZM)
    TH(I,J,K,3) = QH1 * (TAUXY * ZEXM + TAUYY * ZEYM + TAUYZ * ZEZM)
    TH(I,J,K,4) = QH1 * (TAUZX * ZEXM + TAUYZ * ZEYM + TAUZZ * ZEZM)
    TH(I,J,K,5) = QH1 * (R5 * ZEXM + S5 * ZEYM + T5 * ZEZM)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE DIFFZE
! 定義終了 *************************************************************
END SUBROUTINE Turbulence3DLESCSM
!***********************************************************************
!**** 流束関数用リミッター                                          ****
!**** 計算対象 : 単相, 二次元, 圧縮性, k-e                          ****
!***********************************************************************
SUBROUTINE Limiter2DKEM( &
&            AVELIM, &
&            IS, IE, JS, JE, LS, LE, AJA, &
&            QH &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: AVELIM
  INTEGER, INTENT(IN)    :: IS, IE, JS, JE
  INTEGER, INTENT(IN)    :: LS, LE
  REAL,    INTENT(IN)    :: AJA(IS:IE, JS:JE)
  REAL,    INTENT(INOUT) :: QH(IS:IE, JS:JE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  REAL    :: AKMAX
  REAL    :: QHAVE
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 6) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 上限と下限に制限をかける +++++++++++++++++++++++++++++++++++++++++++
  DO J = JS, JE
  DO I = IS, IE
  IF(QH(I,J,1) .GT. 0.0) THEN
    AKMAX     = 0.5 * (QH(I,J,2)**2 + QH(I,J,3)**2) / QH(I,J,1)
    QH(I,J,4) = QH(I,J,4) - QH(I,J,5)
    QH(I,J,5) = MIN(AKMAX, MAX(0.0, QH(I,J,5)))
    QH(I,J,6) = MAX(0.0, QH(I,J,6))
    QH(I,J,4) = QH(I,J,4) + QH(I,J,5)
  ENDIF
  ENDDO
  ENDDO
  ! 平均値から明らかに逸脱する点に制限をかける +++++++++++++++++++++++++
  IF(AVELIM .LE. 0.0) RETURN
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (QH(I  ,J-1,1) .GT. 0.0) .AND. (QH(I-1,J  ,1) .GT. 0.0) .AND. &
  &   (QH(I  ,J  ,1) .GT. 0.0) .AND. &
  &   (QH(I+1,J  ,1) .GT. 0.0) .AND. (QH(I  ,J+1,1) .GT. 0.0) &
  & ) THEN
  DO L = 5, 6
    QHAVE = ( QH(I  ,J-1,L) * AJA(I  ,J-1) &
    &       + QH(I-1,J  ,L) * AJA(I-1,J  ) &
    &       + QH(I+1,J  ,L) * AJA(I+1,J  ) &
    &       + QH(I  ,J+1,L) * AJA(I  ,J+1) &
    &       ) / (4.0 * AJA(I,J))
    QH(I,J,L) = MAX(QHAVE / AVELIM, MIN(QHAVE * AVELIM, QH(I,J,L)))
  ENDDO
  ENDIF
  ENDDO
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE Limiter2DKEM
!***********************************************************************
!**** 流束関数用リミッター                                          ****
!**** 計算対象 : 単相, 三次元, 圧縮性, k-e                          ****
!***********************************************************************
SUBROUTINE Limiter3DKEM( &
&            AVELIM, &
&            IS, IE, JS, JE, KS, KE, LS, LE, AJA, &
&            QH &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: AVELIM
  INTEGER, INTENT(IN)    :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)    :: LS, LE
  REAL,    INTENT(IN)    :: AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(INOUT) :: QH(IS:IE, JS:JE, KS:KE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  REAL    :: AKMAX
  REAL    :: QHAVE
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 7) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 上限と下限に制限をかける +++++++++++++++++++++++++++++++++++++++++++
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(QH(I,J,K,1) .GT. 0.0) THEN
    AKMAX       = 0.5 * ( QH(I,J,K,2)**2 &
    &                   + QH(I,J,K,3)**2 &
    &                   + QH(I,J,K,4)**2) / QH(I,J,K,1)
    QH(I,J,K,5) = QH(I,J,K,5) - QH(I,J,K,6)
    QH(I,J,K,6) = MIN(AKMAX, MAX(0.0, QH(I,J,K,6)))
    QH(I,J,K,7) = MAX(0.0, QH(I,J,K,7))
    QH(I,J,K,5) = QH(I,J,K,5) + QH(I,J,K,6)
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  ! 平均値から明らかに逸脱する点に制限をかける +++++++++++++++++++++++++
  IF(AVELIM .LE. 0.0) RETURN
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( ( QH(I  ,J  ,K-1,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J-1,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I-1,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I+1,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J+1,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J  ,K+1,1) .GT. 0.0 ) &
  & ) THEN
  DO L = 6, 7
    QHAVE = ( QH(I  ,J  ,K-1,L)*AJA(I  ,J  ,K-1) &
    &       + QH(I  ,J-1,K  ,L)*AJA(I  ,J-1,K  ) &
    &       + QH(I-1,J  ,K  ,L)*AJA(I-1,J  ,K  ) &
    &       + QH(I+1,J  ,K  ,L)*AJA(I+1,J  ,K  ) &
    &       + QH(I  ,J+1,K  ,L)*AJA(I  ,J+1,K  ) &
    &       + QH(I  ,J  ,K+1,L)*AJA(I  ,J  ,K+1) &
    &       ) / (6.0 * AJA(I,J,K))
    QH(I,J,K,L) = MAX(QHAVE / AVELIM, MIN(QHAVE * AVELIM, QH(I,J,K,L)))
  ENDDO
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE Limiter3DKEM
!***********************************************************************
!**** 流束関数用リミッター                                          ****
!**** 計算対象 : 単相, 二次元, 圧縮性, RSM                          ****
!***********************************************************************
SUBROUTINE Limiter2DRSM( &
&            AVELIM, &
&            IS, IE, JS, JE, LS, LE, AJA, &
&            QH &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: AVELIM
  INTEGER, INTENT(IN)    :: IS, IE, JS, JE
  INTEGER, INTENT(IN)    :: LS, LE
  REAL,    INTENT(IN)    :: AJA(IS:IE, JS:JE)
  REAL,    INTENT(INOUT) :: QH(IS:IE, JS:JE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  REAL    :: uiuimax, uiujmax
  REAL    :: QHAVE
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 9) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 上限と下限に制限をかける +++++++++++++++++++++++++++++++++++++++++++
  DO J = JS, JE
  DO I = IS, IE
  IF(QH(I,J,1) .GT. 0.0) THEN
    uiuimax   = (QH(I,J,2)**2 + QH(I,J,3)**2) / QH(I,J,1)
    QH(I,J,4) = QH(I,J,4) - 0.5 * (QH(I,J,5) + QH(I,J,6) + QH(I,J,7))
    QH(I,J,5) = MIN(uiuimax, MAX(0.0, QH(I,J,5)))
    QH(I,J,6) = MIN(uiuimax, MAX(0.0, QH(I,J,6)))
    QH(I,J,7) = MIN(uiuimax, MAX(0.0, QH(I,J,7)))
    uiujmax   = 0.5 * (QH(I,J,5) + QH(I,J,6) + QH(I,J,7))
    QH(I,J,8) = MIN(uiujmax, MAX(- uiujmax, QH(I,J,8)))
    QH(I,J,9) = MAX(0.0, QH(I,J,9))
    QH(I,J,4) = QH(I,J,4) + 0.5 * (QH(I,J,5) + QH(I,J,6) + QH(I,J,7))
  ENDIF
  ENDDO
  ENDDO
  ! 平均値から明らかに逸脱する点に制限をかける +++++++++++++++++++++++++
  IF(AVELIM .LE. 0.0) RETURN
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( (QH(I  ,J-1,1) .GT. 0.0) .AND. (QH(I-1,J  ,1) .GT. 0.0) .AND. &
  &   (QH(I  ,J  ,1) .GT. 0.0) .AND. &
  &   (QH(I+1,J  ,1) .GT. 0.0) .AND. (QH(I  ,J+1,1) .GT. 0.0) &
  & ) THEN
  L = 9
    QHAVE = ( QH(I  ,J-1,L) * AJA(I  ,J-1) &
    &       + QH(I-1,J  ,L) * AJA(I-1,J  ) &
    &       + QH(I+1,J  ,L) * AJA(I+1,J  ) &
    &       + QH(I  ,J+1,L) * AJA(I  ,J+1) &
    &       ) / (4.0 * AJA(I,J))
    QH(I,J,L) = MAX(QHAVE / AVELIM, MIN(QHAVE * AVELIM, QH(I,J,L)))
  ENDIF
  ENDDO
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE Limiter2DRSM
!***********************************************************************
!**** 流束関数用リミッター                                          ****
!**** 計算対象 : 単相, 三次元, 圧縮性, RSM                          ****
!***********************************************************************
SUBROUTINE Limiter3DRSM( &
&            AVELIM, &
&            IS, IE, JS, JE, KS, KE, LS, LE, AJA, &
&            QH &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)    :: AVELIM
  INTEGER, INTENT(IN)    :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)    :: LS, LE
  REAL,    INTENT(IN)    :: AJA(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(INOUT) :: QH(IS:IE, JS:JE, KS:KE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  REAL    :: uiuimax, uiujmax
  REAL    :: QHAVE
  ! 処理開始 ***********************************************************
  ! 例外処理(計算対象が適合していない場合は終了) +++++++++++++++++++++++
  IF(LS .NE. 1 .OR. LE .LT. 12) THEN
    WRITE(*,'(A)') 'Error : LS or LE Irregal'
    STOP
  ENDIF
  ! 上限と下限に制限をかける +++++++++++++++++++++++++++++++++++++++++++
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(QH(I,J,K,1) .GT. 0.0) THEN
    uiuimax = ( QH(I,J,K,2)**2 &
    &         + QH(I,J,K,3)**2 &
    &         + QH(I,J,K,4)**2 ) / QH(I,J,K,1)
    QH(I,J,K, 5) = QH(I,J,K,5) &
    &            - 0.5 * (QH(I,J,K,6) + QH(I,J,K,7) + QH(I,J,K,8))
    QH(I,J,K, 6) = MIN(uiuimax, MAX(0.0, QH(I,J,K, 6)))
    QH(I,J,K, 7) = MIN(uiuimax, MAX(0.0, QH(I,J,K, 7)))
    QH(I,J,K, 8) = MIN(uiuimax, MAX(0.0, QH(I,J,K, 8)))
    uiujmax = 0.5 * (QH(I,J,K,6) + QH(I,J,K,7) + QH(I,J,K,8))
    QH(I,J,K, 9) = MIN(uiujmax, MAX(- uiujmax, QH(I,J,K, 9)))
    QH(I,J,K,10) = MIN(uiujmax, MAX(- uiujmax, QH(I,J,K,10)))
    QH(I,J,K,11) = MIN(uiujmax, MAX(- uiujmax, QH(I,J,K,11)))
    QH(I,J,K,12) = MAX(0.0, QH(I,J,K,12))
    QH(I,J,K, 5) = QH(I,J,K,5) &
    &            + 0.5 * (QH(I,J,K,6) + QH(I,J,K,7) + QH(I,J,K,8))
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  ! 平均値から明らかに逸脱する点に制限をかける +++++++++++++++++++++++++
  IF(AVELIM .LE. 0.0) RETURN
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF( ( QH(I  ,J  ,K-1,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J-1,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I-1,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I+1,J  ,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J+1,K  ,1) .GT. 0.0 ) .AND. &
  &   ( QH(I  ,J  ,K+1,1) .GT. 0.0 ) &
  & ) THEN
  L = 12
    QHAVE = ( QH(I  ,J  ,K-1,L)*AJA(I  ,J  ,K-1) &
    &       + QH(I  ,J-1,K  ,L)*AJA(I  ,J-1,K  ) &
    &       + QH(I-1,J  ,K  ,L)*AJA(I-1,J  ,K  ) &
    &       + QH(I+1,J  ,K  ,L)*AJA(I+1,J  ,K  ) &
    &       + QH(I  ,J+1,K  ,L)*AJA(I  ,J+1,K  ) &
    &       + QH(I  ,J  ,K+1,L)*AJA(I  ,J  ,K+1) &
    &       ) / (6.0 * AJA(I,J,K))
    QH(I,J,K,L) = MAX(QHAVE / AVELIM, MIN(QHAVE * AVELIM, QH(I,J,K,L)))
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE Limiter3DRSM
!***********************************************************************
!**** 残差計算                                                      ****
!**** 計算対象 : 単相, 二次元, 圧縮性                               ****
!***********************************************************************
SUBROUTINE CalResidual2D( &
&            IS, IE, JS, JE, LS, LE, &
&            DTLOCL, QH0, QH, RES &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: DTLOCL(IS:IE, JS:JE)
  REAL,    INTENT(IN)  :: QH0(IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(IN)  :: QH(IS:IE, JS:JE, LS:LE)
  REAL,    INTENT(OUT) :: RES(IS:IE, JS:JE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, L
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, L)
  DO L = LS, LE
  !$OMP DO
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF(QH0(I,J,1) .GT. 0.0 .AND. DTLOCL(I,J) .GT. 0.0) THEN
    RES(I,J,L) = ((QH(I,J,L) - QH0(I,J,L)) / QH(I,J,1))**2 &
    &          / DTLOCL(I,J)
  ELSE
    RES(I,J,L) = 0.0
  ENDIF
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE CalResidual2D
!***********************************************************************
!**** 残差計算                                                      ****
!**** 計算対象 : 単相, 三次元, 圧縮性                               ****
!***********************************************************************
SUBROUTINE CalResidual3D( &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            DTLOCL, QH0, QH, RES &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  INTEGER, INTENT(IN)  :: LS, LE
  REAL,    INTENT(IN)  :: DTLOCL(IS:IE, JS:JE, KS:KE)
  REAL,    INTENT(IN)  :: QH0(IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(IN)  :: QH(IS:IE, JS:JE, KS:KE, LS:LE)
  REAL,    INTENT(OUT) :: RES(IS:IE, JS:JE, KS:KE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  !$OMP DO
  DO K = KS + 1, KE - 1
  DO J = JS + 1, JE - 1
  DO I = IS + 1, IE - 1
  IF(QH0(I,J,K,1) .GT. 0.0 .AND. DTLOCL(I,J,K) .GT. 0.0) THEN
    RES(I,J,K,L) = ((QH(I,J,K,L) - QH0(I,J,K,L)) / QH(I,J,K,1))**2 &
    &            / DTLOCL(I,J,K)
  ELSE
    RES(I,J,K,L) = 0.0
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END DO
  ENDDO
  !$OMP END PARALLEL
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE CalResidual3D
!***********************************************************************
!**** Yapの長さスケール補正(1987)                                   ****
!****   KEMの剥離流れにおける壁近傍の長さスケールの過大評価、       ****
!****   すなわちepsilonの過小評価を補正するepsilonの付加的な生成項  ****
!***********************************************************************
SUBROUTINE YapCorrection2D( &
&            IS, IE, JS, JE, Y, AJA, RHO, AK, EPS, Yc &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE
  REAL,    INTENT(IN)  :: Y(IS: IE, JS: JE)
  REAL,    INTENT(IN)  :: AJA(IS: IE, JS: JE)
  REAL,    INTENT(IN)  :: RHO(IS: IE, JS: JE)
  REAL,    INTENT(IN)  :: AK(IS: IE, JS: JE)
  REAL,    INTENT(IN)  :: EPS(IS: IE, JS: JE)
  REAL,    INTENT(OUT) :: Yc(IS: IE, JS: JE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J
  REAL    :: Le, k32PeLe
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, Le, k32PeLe)
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J) .GT. 0.0) THEN
    Le      = 2.5 * Y(I,J)
    k32PeLe = AK(I,J)**1.5 / MAX(ZERO, Le * EPS(I,J))
    Yc(I,J) = 0.83 * EPS(I,J)**2 / MAX(ZERO, AK(I,J)) &
    &       * (k32PeLe - 1.0) * k32PeLe**2
    Yc(I,J) = RHO(I,J) / AJA(I,J) * MAX(Yc(I,J), 0.0)
  ELSE
    Yc(I,J) = 0.0
  ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE YapCorrection2D
!***********************************************************************
!**** Yapの長さスケール補正(1987)                                   ****
!****   KEMの剥離流れにおける壁近傍の長さスケールの過大評価、       ****
!****   すなわちepsilonの過小評価を補正するepsilonの付加的な生成項  ****
!***********************************************************************
SUBROUTINE YapCorrection3D( &
&            IS, IE, JS, JE, KS, KE, Y, AJA, RHO, AK, EPS, Yc &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: ZERO = 1.0E-20
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE
  REAL,    INTENT(IN)  :: Y(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)  :: AJA(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)  :: RHO(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)  :: AK(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(IN)  :: EPS(IS: IE, JS: JE, KS: KE)
  REAL,    INTENT(OUT) :: Yc(IS: IE, JS: JE, KS: KE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K
  REAL    :: Le, k32PeLe
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, K, Le, k32PeLe)
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
  IF(RHO(I,J,K) .GT. 0.0) THEN
    Le        = 2.5 * Y(I,J,K)
    k32PeLe   = AK(I,J,K)**1.5 / MAX(ZERO, Le * EPS(I,J,K))
    Yc(I,J,K) = 0.83 * EPS(I,J,K)**2 / MAX(ZERO, AK(I,J,K)) &
    &         * (k32PeLe - 1.0) * k32PeLe**2
    Yc(I,J,K) = RHO(I,J,K) / AJA(I,J,K) * MAX(Yc(I,J,K), 0.0)
  ELSE
    Yc(I,J,K) = 0.0
  ENDIF
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE YapCorrection3D
!***********************************************************************
!**** KEM用の壁関数(滑面)                                           ****
!***********************************************************************
SUBROUTINE WallFunctionKEM1S( &
&            yp, up, nup, kp0, epsp0, &
&            utau, dudy1, dudy2, kp, epsp &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: cmu = 0.09
  ! REAL, PARAMETER :: kappa = 0.41, E = 9.0
  REAL, PARAMETER :: kappa = 0.40, B = 5.5
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN)  :: yp                 ! 壁からの距離   (yp > 0)
  REAL, INTENT(IN)  :: up                 ! 壁水平方向速度 (up > 0)
  REAL, INTENT(IN)  :: nup                ! 動粘度         (nup > 0)
  REAL, INTENT(IN)  :: kp0, epsp0         ! 前ステップの乱流量
  REAL, INTENT(OUT) :: utau               ! 摩擦速度       (utau > 0)
  REAL, INTENT(OUT) :: dudy1, dudy2       ! 速度の一階微分と二階微分
  REAL, INTENT(OUT) :: kp, epsp           ! 乱流量
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL :: utau0, yplusp
  ! 処理開始 ***********************************************************
  utau0  = cmu / (kappa * yp) * kp0**2 / epsp0
  yplusp = yp * utau0 / nup
  IF(yplusp .LE. 11.635)THEN
    utau = SQRT(nup * up / yp)
  ELSE
    ! utau = kappa * up / LOG(yplusp * E)
    utau = up / (LOG(yplusp) / kappa + B)
  ENDIF
  dudy1 = utau / (kappa * yp)
  dudy2 =-dudy1 / yp
  kp    = utau**2 / SQRT(cmu)
  epsp  = utau**3 / (kappa * yp)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE WallFunctionKEM1S
!***********************************************************************
!**** RSM用の壁関数(滑面)                                           ****
!***********************************************************************
SUBROUTINE WallFunctionRSM1S( &
&            yp, up, nup, uv0, epsp0, &
&            utau, dudy1, dudy2, uup, vvp, wwp, uvp, epsp &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! REAL, PARAMETER :: kappa = 0.41, E = 9.0
  REAL, PARAMETER :: kappa = 0.40, B = 5.5
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN)  :: yp                   ! 壁からの距離   (yp > 0)
  REAL, INTENT(IN)  :: up                   ! 壁水平方向速度 (up > 0)
  REAL, INTENT(IN)  :: nup                  ! 動粘度         (nup > 0)
  REAL, INTENT(IN)  :: uv0                  ! 前ステップのレイノルズ応力
  REAL, INTENT(IN)  :: epsp0                ! 前ステップの乱れの散逸率
  REAL, INTENT(OUT) :: utau                 ! 摩擦速度       (utau > 0)
  REAL, INTENT(OUT) :: dudy1, dudy2         ! 速度の一階微分と二階微分
  REAL, INTENT(OUT) :: uup, vvp, wwp, uvp   ! レイノルズ応力
  REAL, INTENT(OUT) :: epsp                 ! 乱れの散逸率
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL :: utau0, yplusp
  ! 処理開始 ***********************************************************
  utau0  = SQRT(ABS(uv0))
  yplusp = yp * utau0 / nup
  IF(yplusp .LE. 11.635)THEN
    utau = SQRT(nup * up / yp)
  ELSE
    ! utau = kappa * up / LOG(yplusp * E)
    utau = up / (LOG(yplusp) / kappa + B)
  ENDIF
  dudy1 = utau / (kappa * yp)
  dudy2 =-dudy1 / yp
  uvp   =-utau**2
  uup   =-4.9 * uvp
  vvp   =-1.0 * uvp
  wwp   =-2.4 * uvp
  epsp  = utau**3 / (kappa * yp)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE WallFunctionRSM1S
!***********************************************************************
!**** KEM用の壁関数(滑面)                                           ****
!***********************************************************************
SUBROUTINE WallFunctionKEM2S( &
&            yp, up, nup, &
&            utau, dudy1, dudy2, kp, epsp &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: cmu = 0.09
  REAL, PARAMETER :: kappa = 0.40, B = 5.5
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: yp                 ! 壁からの距離   (yp > 0)
  REAL,    INTENT(IN)  :: up                 ! 壁水平方向速度 (up > 0)
  REAL,    INTENT(IN)  :: nup                ! 動粘度         (nup > 0)
  REAL,    INTENT(OUT) :: utau               ! 摩擦速度       (utau > 0)
  REAL,    INTENT(OUT) :: dudy1, dudy2       ! 速度の一階微分と二階微分
  REAL,    INTENT(OUT) :: kp, epsp           ! 乱流量
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! 処理開始 ***********************************************************
  CALL CalUtauS(kappa, B, yp, up, nup, utau, dudy1, dudy2)
  kp   = utau**2 / SQRT(cmu)
  epsp = utau**3 / (kappa * yp)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE WallFunctionKEM2S
!***********************************************************************
!**** RSM用の壁関数(滑面)                                           ****
!***********************************************************************
SUBROUTINE WallFunctionRSM2S( &
&            yp, up, nup, &
&            utau, dudy1, dudy2, uup, vvp, wwp, uvp, epsp &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: kappa = 0.40, B = 5.5
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: yp                 ! 壁からの距離   (yp > 0)
  REAL,    INTENT(IN)  :: up                 ! 壁水平方向速度 (up > 0)
  REAL,    INTENT(IN)  :: nup                ! 動粘度         (nup > 0)
  REAL,    INTENT(OUT) :: utau               ! 摩擦速度       (utau > 0)
  REAL,    INTENT(OUT) :: dudy1, dudy2       ! 速度の一階微分と二階微分
  REAL,    INTENT(OUT) :: uup, vvp, wwp, uvp ! レイノルズ応力
  REAL,    INTENT(OUT) :: epsp               ! 乱れの散逸率
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! 処理開始 ***********************************************************
  CALL CalUtauS(kappa, B, yp, up, nup, utau, dudy1, dudy2)
  uvp  = - utau**2
  uup  = - 4.9 * uvp
  vvp  = - 1.0 * uvp
  wwp  = - 2.4 * uvp
  epsp = utau**3 / (kappa * yp)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE WallFunctionRSM2S
!***********************************************************************
!**** 摩擦速度の計算(対数則)                                        ****
!***********************************************************************
SUBROUTINE CalUtauS( &
&            kappa, B, yp, up, nup, &
&            utau, dudy1, dudy2 &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, PARAMETER :: NMAX = 100
  REAL,    PARAMETER :: ZERO = 1.0E-20, RESMIN = 1.0E-8
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: kappa, B           ! 対数則のモデル定数
  REAL,    INTENT(IN)  :: yp                 ! 壁からの距離   (yp > 0)
  REAL,    INTENT(IN)  :: up                 ! 壁水平方向速度 (up > 0)
  REAL,    INTENT(IN)  :: nup                ! 動粘度         (nup > 0)
  REAL,    INTENT(OUT) :: utau               ! 摩擦速度       (utau > 0)
  REAL,    INTENT(OUT) :: dudy1, dudy2       ! 速度の一階微分と二階微分
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: n
  REAL    :: f, futau, dutau
  ! 処理開始 ***********************************************************
  utau = SQRT(nup * up / yp)
  DO n = 1, NMAX
    f     = up / utau - LOG(yp * utau / nup) / kappa - B
    futau = - (up / utau + 1 / kappa) / utau
    IF(ABS(futau) .LE. ZERO) futau = SIGN(ZERO, futau)
    dutau = - f / futau
    utau  = MAX(ZERO, utau + dutau)
    IF(ABS(dutau) .LT. utau * RESMIN) EXIT
  ENDDO
  IF(yp * utau / nup .LE. 11.635) THEN
    utau = SQRT(nup * up / yp)
  ENDIF
  dudy1 = utau / (kappa * yp)
  dudy2 =-dudy1 / yp
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE CalUtauS
!***********************************************************************
!**** KEM用の壁関数(粗面)                                           ****
!***********************************************************************
SUBROUTINE WallFunctionKEM1R( &
&            ks, yp, up, nup, kp0, epsp0, &
&            utau, dudy1, dudy2, kp, epsp &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: cmu = 0.09
  REAL, PARAMETER :: kappa = 0.40, B = 5.5
  REAL, PARAMETER :: ksp1 = 4.0, ksp2 = 15.0, ksp3 = 55.0
  REAL, PARAMETER :: At = 9.5, Ar = 8.5
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN)  :: ks                    ! 粗度要素高さ   (ks > 0)
  REAL, INTENT(IN)  :: yp                    ! 壁からの距離   (yp > 0)
  REAL, INTENT(IN)  :: up                    ! 壁水平方向速度 (up > 0)
  REAL, INTENT(IN)  :: nup                   ! 動粘度         (nup > 0)
  REAL, INTENT(IN)  :: kp0, epsp0            ! 前ステップの乱流量
  REAL, INTENT(OUT) :: utau                  ! 摩擦速度       (utau > 0)
  REAL, INTENT(OUT) :: dudy1, dudy2          ! 速度の一階微分と二階微分
  REAL, INTENT(OUT) :: kp, epsp              ! 乱流量
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL :: utau0, yplusp
  REAL :: ksp
  REAL :: lnksp1, lnksp2, lnksp3
  REAL :: alpha, beta, gamma
  REAL :: DB
  ! 処理開始 ***********************************************************
  utau0  = cmu / (kappa * yp) * kp0**2 / epsp0
  yplusp = yp * utau0 / nup
  ksp    = ks * utau0 / nup
  IF(yplusp .LE. 11.635)THEN
    utau = SQRT(nup * up / yp)
  ELSE
    lnksp1 = LOG(ksp1)
    lnksp2 = LOG(ksp2)
    lnksp3 = LOG(ksp3)
    beta   = 0.5 &
    &      * ( (B - Ar + lnksp3 / kappa) * (lnksp2**2 - lnksp1**2) &
    &        - (B - At + lnksp2 / kappa) * (lnksp3**2 - lnksp1**2) ) &
    &      / ( (B - At + lnksp2 / kappa) * (lnksp1 - lnksp3) &
    &        - (B - Ar + lnksp3 / kappa) * (lnksp1 - lnksp2) )
    alpha  = (B - At + lnksp2 / kappa) &
    &      / ((lnksp2 - beta)**2 - (lnksp1 - beta)**2)
    gamma  = - alpha * (lnksp1 - beta)**2
    IF(ksp .LT. ksp1) THEN
      DB = 0.0
    ELSEIF(ksp .GT. ksp3) THEN
      DB = B - Ar + LOG(ksp) / kappa
    ELSE
      DB = alpha * (LOG(ksp) - beta)**2 + gamma
    ENDIF
    utau = up / (LOG(yplusp) / kappa + B - DB)
  ENDIF
  dudy1 = utau / (kappa * yp)
  dudy2 =-dudy1 / yp
  kp    = utau**2 / SQRT(cmu)
  epsp  = utau**3 / (kappa * yp)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE WallFunctionKEM1R
!***********************************************************************
!**** RSM用の壁関数(粗面)                                           ****
!***********************************************************************
SUBROUTINE WallFunctionRSM1R( &
&            ks, yp, up, nup, uv0, epsp0, &
&            utau, dudy1, dudy2, uup, vvp, wwp, uvp, epsp &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: kappa = 0.40, B = 5.5
  REAL, PARAMETER :: ksp1 = 4.0, ksp2 = 15.0, ksp3 = 55.0
  REAL, PARAMETER :: At = 9.5, Ar = 8.5
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN)  :: ks                   ! 粗度要素高さ   (ks > 0)
  REAL, INTENT(IN)  :: yp                   ! 壁からの距離   (yp > 0)
  REAL, INTENT(IN)  :: up                   ! 壁水平方向速度 (up > 0)
  REAL, INTENT(IN)  :: nup                  ! 動粘度         (nup > 0)
  REAL, INTENT(IN)  :: uv0                  ! 前ステップのレイノルズ応力
  REAL, INTENT(IN)  :: epsp0                ! 前ステップの乱れの散逸率
  REAL, INTENT(OUT) :: utau                 ! 摩擦速度       (utau > 0)
  REAL, INTENT(OUT) :: dudy1, dudy2         ! 速度の一階微分と二階微分
  REAL, INTENT(OUT) :: uup, vvp, wwp, uvp   ! レイノルズ応力
  REAL, INTENT(OUT) :: epsp                 ! 乱れの散逸率
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL :: utau0, yplusp
  REAL :: ksp
  REAL :: lnksp1, lnksp2, lnksp3
  REAL :: alpha, beta, gamma
  REAL :: DB
  ! 処理開始 ***********************************************************
  utau0  = SQRT(ABS(uv0))
  yplusp = yp * utau0 / nup
  ksp    = ks * utau0 / nup
  IF(yplusp .LE. 11.635)THEN
    utau = SQRT(nup * up / yp)
  ELSE
    lnksp1 = LOG(ksp1)
    lnksp2 = LOG(ksp2)
    lnksp3 = LOG(ksp3)
    beta   = 0.5 &
    &      * ( (B - Ar + lnksp3 / kappa) * (lnksp2**2 - lnksp1**2) &
    &        - (B - At + lnksp2 / kappa) * (lnksp3**2 - lnksp1**2) ) &
    &      / ( (B - At + lnksp2 / kappa) * (lnksp1 - lnksp3) &
    &        - (B - Ar + lnksp3 / kappa) * (lnksp1 - lnksp2) )
    alpha  = (B - At + lnksp2 / kappa) &
    &      / ((lnksp2 - beta)**2 - (lnksp1 - beta)**2)
    gamma  = - alpha * (lnksp1 - beta)**2
    IF(ksp .LT. ksp1) THEN
      DB = 0.0
    ELSEIF(ksp .GT. ksp3) THEN
      DB = B - Ar + LOG(ksp) / kappa
    ELSE
      DB = alpha * (LOG(ksp) - beta)**2 + gamma
    ENDIF
    utau = up / (LOG(yplusp) / kappa + B - DB)
  ENDIF
  dudy1 = utau / (kappa * yp)
  dudy2 =-dudy1 / yp
  uvp  = - utau**2
  uup  = - 4.9 * uvp
  vvp  = - 1.0 * uvp
  wwp  = - 2.4 * uvp
  epsp = utau**3 / (kappa * yp)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE WallFunctionRSM1R
!***********************************************************************
!**** KEM用の壁関数(粗面)                                           ****
!***********************************************************************
SUBROUTINE WallFunctionKEM2R( &
&            ks, yp, up, nup, &
&            utau, dudy1, dudy2, kp, epsp &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: cmu = 0.09
  REAL, PARAMETER :: kappa = 0.40, B = 5.5
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN)  :: ks                    ! 粗度要素高さ   (ks > 0)
  REAL, INTENT(IN)  :: yp                    ! 壁からの距離   (yp > 0)
  REAL, INTENT(IN)  :: up                    ! 壁水平方向速度 (up > 0)
  REAL, INTENT(IN)  :: nup                   ! 動粘度         (nup > 0)
  REAL, INTENT(OUT) :: utau                  ! 摩擦速度       (utau > 0)
  REAL, INTENT(OUT) :: dudy1, dudy2          ! 速度の一階微分と二階微分
  REAL, INTENT(OUT) :: kp, epsp              ! 乱流量
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! 処理開始 ***********************************************************
  CALL CalUtauR(kappa, B, ks, yp, up, nup, utau, dudy1, dudy2)
  kp   = utau**2 / SQRT(cmu)
  epsp = utau**3 / (kappa * yp)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE WallFunctionKEM2R
!***********************************************************************
!**** RSM用の壁関数(粗面)                                           ****
!***********************************************************************
SUBROUTINE WallFunctionRSM2R( &
&            ks, yp, up, nup, &
&            utau, dudy1, dudy2, uup, vvp, wwp, uvp, epsp &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, PARAMETER :: kappa = 0.40, B = 5.5
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL, INTENT(IN)  :: ks                    ! 粗度要素高さ   (ks > 0)
  REAL, INTENT(IN)  :: yp                    ! 壁からの距離   (yp > 0)
  REAL, INTENT(IN)  :: up                    ! 壁水平方向速度 (up > 0)
  REAL, INTENT(IN)  :: nup                   ! 動粘度         (nup > 0)
  REAL, INTENT(OUT) :: utau                  ! 摩擦速度       (utau > 0)
  REAL, INTENT(OUT) :: dudy1, dudy2          ! 速度の一階微分と二階微分
  REAL, INTENT(OUT) :: uup, vvp, wwp, uvp    ! レイノルズ応力
  REAL, INTENT(OUT) :: epsp                  ! 乱れの散逸率
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! 処理開始 ***********************************************************
  CALL CalUtauR(kappa, B, ks, yp, up, nup, utau, dudy1, dudy2)
  uvp  = - utau**2
  uup  = - 4.9 * uvp
  vvp  = - 1.0 * uvp
  wwp  = - 2.4 * uvp
  epsp = utau**3 / (kappa * yp)
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE WallFunctionRSM2R
!***********************************************************************
!**** 摩擦速度の計算(粗面の対数則)                                  ****
!***********************************************************************
SUBROUTINE CalUtauR( &
&            kappa, B, ks, yp, up, nup, &
&            utau, dudy1, dudy2 &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 局所定数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, PARAMETER :: NMAX = 100
  REAL,    PARAMETER :: ZERO = 1.0E-20, RESMIN = 1.0E-8
  REAL,    PARAMETER :: ksp1 = 4.0, ksp2 = 15.0, ksp3 = 55.0
  REAL,    PARAMETER :: At = 9.5, Ar = 8.5
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL,    INTENT(IN)  :: kappa, B           ! 対数則のモデル定数
  REAL,    INTENT(IN)  :: ks                 ! 粗度要素高さ   (ks > 0)
  REAL,    INTENT(IN)  :: yp                 ! 壁からの距離   (yp > 0)
  REAL,    INTENT(IN)  :: up                 ! 壁水平方向速度 (up > 0)
  REAL,    INTENT(IN)  :: nup                ! 動粘度         (nup > 0)
  REAL,    INTENT(OUT) :: utau               ! 摩擦速度       (utau > 0)
  REAL,    INTENT(OUT) :: dudy1, dudy2       ! 速度の一階微分と二階微分
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: n
  REAL    :: ksp
  REAL    :: lnksp1, lnksp2, lnksp3
  REAL    :: alpha, beta, gamma
  REAL    :: DB, DButau, f, futau, dutau
  ! 処理開始 ***********************************************************
  lnksp1 = LOG(ksp1)
  lnksp2 = LOG(ksp2)
  lnksp3 = LOG(ksp3)
  beta   = 0.5 &
  &      * ( (B - Ar + lnksp3 / kappa) * (lnksp2**2 - lnksp1**2) &
  &        - (B - At + lnksp2 / kappa) * (lnksp3**2 - lnksp1**2) ) &
  &      / ( (B - At + lnksp2 / kappa) * (lnksp1 - lnksp3) &
  &        - (B - Ar + lnksp3 / kappa) * (lnksp1 - lnksp2) )
  alpha  = (B - At + lnksp2 / kappa) &
  &      / ((lnksp2 - beta)**2 - (lnksp1 - beta)**2)
  gamma  = - alpha * (lnksp1 - beta)**2
  utau = SQRT(nup * up / yp)
  DO n = 1, NMAX
    ksp  = ks * utau / nup
    IF(ksp .LT. ksp1) THEN
      DB     = 0.0
      DButau = 0.0
    ELSEIF(ksp .GT. ksp3) THEN
      DB     = B - Ar + LOG(ksp) / kappa
      DButau = 1.0 / (kappa * ksp)
    ELSE
      DB     = alpha * (LOG(ksp) - beta)**2 + gamma
      DButau = 2.0 * alpha * (LOG(ksp) - beta) / ksp
    ENDIF
    f     = up / utau - LOG(yp * utau / nup) / kappa - B + DB
    futau = - (up / utau + 1 / kappa) / utau + DButau
    IF(ABS(futau) .LE. ZERO) futau = SIGN(ZERO, futau)
    dutau = - f / futau
    utau  = MAX(ZERO, utau + dutau)
    IF(ABS(dutau) .LT. utau * RESMIN) EXIT
  ENDDO
  IF(yp * utau / nup .LE. 11.635) THEN
    utau = SQRT(nup * up / yp)
  ENDIF
  dudy1 = utau / (kappa * yp)
  dudy2 =-dudy1 / yp
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE CalUtauR
!***********************************************************************
!**** 一ステップ前の流束をを保存                                    ****
!**** 計算対象 : 単相, 三次元, 圧縮性, k-e                          ****
!***********************************************************************
SUBROUTINE SaveFlux3DKEM( &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            QH, DQH, QH0, DQH0  &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE, LS, LE
  REAL   , INTENT(IN)  :: QH  (IS:IE, JS:JE, KS:KE, LS:LE)
  REAL   , INTENT(IN)  :: DQH (IS:IE, JS:JE, KS:KE, LS:LE)
  REAL   , INTENT(OUT) :: QH0 (IS:IE, JS:JE, KS:KE, LS:LE)
  REAL   , INTENT(OUT) :: DQH0(IS:IE, JS:JE, KS:KE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
    QH0(I,J,K,L)  = QH(I,J,K,L)
    DQH0(I,J,K,L) = DQH(I,J,K,L)
  ENDDO
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SaveFlux3DKEM
!***********************************************************************
!**** 各項の総和		                                    ****
!**** 計算対象 : 単相, 三次元, 圧縮性, k-e                          ****
!***********************************************************************
SUBROUTINE SumDQH3D( &
&            IS, IE, JS, JE, KS, KE, LS, LE, &
&            DQC, DQD, DQP, DQR, DQH &
&          )
  ! 変数宣言 ***********************************************************
  IMPLICIT NONE
  ! 引数変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER, INTENT(IN)  :: IS, IE, JS, JE, KS, KE, LS, LE
  REAL   , INTENT(IN)  :: DQC(IS:IE, JS:JE, KS:KE, LS:LE)
  REAL   , INTENT(IN)  :: DQD(IS:IE, JS:JE, KS:KE, LS:LE)
  REAL   , INTENT(IN)  :: DQP(IS:IE, JS:JE, KS:KE, LS:LE)
  REAL   , INTENT(IN)  :: DQR(IS:IE, JS:JE, KS:KE, LS:LE)
  REAL   , INTENT(OUT) :: DQH(IS:IE, JS:JE, KS:KE, LS:LE)
  ! 局所変数 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  INTEGER :: I, J, K, L
  ! 処理開始 ***********************************************************
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, K, L)
  DO L = LS, LE
  DO K = KS, KE
  DO J = JS, JE
  DO I = IS, IE
    DQH(I,J,K,L) = DQC(I,J,K,L) + DQD(I,J,K,L) &
  &              + DQP(I,J,K,L) + DQR(I,J,K,L)
  ENDDO
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ! 処理終了 ***********************************************************
  RETURN
END SUBROUTINE SumDQH3D
! 定義終了 *************************************************************
END MODULE Package_Flow
