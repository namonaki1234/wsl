!*********************************************************************
!****** MODULE_FOR_INITIAL_FILE_MAKER(2D)                          ***
!****** VER. 2014.05.26                                            ***
!*********************************************************************
MODULE MOD_IS_DIM
  IMPLICIT NONE
  ! SUBROUTINE_DEFINITION ********************************************
  PUBLIC :: GRID_SETTING_2D
  PUBLIC :: BOX_POSITIONING_2D
  PUBLIC :: CIRCLE_POSITIONING_2D
  PUBLIC :: WALL_SETTING

  ! FUNCTION_DEFINITION ***********************************************
  PUBLIC :: DIST_2D
CONTAINS

!*********************************************************************
!*** WALL_SETTING_FOR_2D                                           ***
!*********************************************************************
SUBROUTINE WALL_SETTING                         &
&          (NUM,X1,Y1,X2,Y2,IDIM,WALL_N,WALL_ANC,WALL_OUT,WALL_NOR,DIS)
  IMPLICIT NONE
  ! MAINROUTINE_DEFINITION *******************************************
  INTEGER, INTENT(IN)    :: WALL_N
  REAL   , INTENT(IN)    :: DIS
  INTEGER, INTENT(IN)    :: IDIM
  REAL   , INTENT(IN)    :: X1,Y1,X2,Y2
  INTEGER, INTENT(IN)    :: NUM
  REAL   , INTENT(INOUT) :: WALL_OUT(0:IDIM-1,0:WALL_N-1)
  REAL   , INTENT(INOUT) :: WALL_NOR(0:IDIM-1,0:WALL_N-1)
  REAL   , INTENT(INOUT) :: WALL_ANC(0:IDIM-1,0:WALL_N-1,0:1)

  ! SUBROUTINE_DEFINITION ********************************************
  REAL :: SIN_TH,COS_TH
  REAL :: CENTER1,CENTER2
  REAL :: ZERO = 1.0E-10
  
  ! WALL_ANCHOR_SETTING **********************************************
  WALL_ANC(0,NUM,0) = X1
  WALL_ANC(1,NUM,0) = Y1
  WALL_ANC(0,NUM,1) = X2
  WALL_ANC(1,NUM,1) = Y2

  ! WALL_OUT_CALCULATION ********************************************
  CENTER1 = (X1+X2)*0.5
  CENTER2 = (Y1+Y2)*0.5
  SIN_TH = (Y2-Y1)/MAX(ZERO,DIST_2D(X1,Y1,X2,Y2))
  COS_TH = (X2-X1)/MAX(ZERO,DIST_2D(X1,Y1,X2,Y2))

  WALL_OUT(0,NUM) = CENTER1+SIN_TH*DIS
  WALL_OUT(1,NUM) = CENTER2-COS_TH*DIS

  ! WALL_NOR_SETTING *************************************************
  WALL_NOR(0,NUM) = SIN_TH
  WALL_NOR(1,NUM) = -COS_TH

  RETURN
END SUBROUTINE WALL_SETTING

!*********************************************************************
!*** GRID_AND_PARTICLE_INTERFERENCE                                ***
!*********************************************************************
SUBROUTINE GRID_AND_PARTICLE                 &
&          (IDIM,NN,IFLUID,WALL_N,WALL_ANC,  &
&           NUMP,DIS,TYPE,X,ITYPEP,WALL_OUT, &
&           WALL_NOR,DIS_MIN)
  IMPLICIT NONE
  ! MAINROUTINE_VARIABLE *********************************************
  INTEGER, INTENT(IN)    :: IDIM
  INTEGER, INTENT(IN)    :: NN
  INTEGER, INTENT(IN)    :: WALL_N
  REAL   , INTENT(IN)    :: WALL_ANC(0:IDIM-1,0:WALL_N-1,0:IDIM-1)
  INTEGER, INTENT(IN)    :: IFLUID
  INTEGER, INTENT(IN)    :: NUMP
  REAL   , INTENT(IN)    :: DIS
  INTEGER, INTENT(IN)    :: TYPE(0:IFLUID)
  REAL   , INTENT(IN)    :: WALL_OUT(0:IDIM-1,0:WALL_N-1)
  REAL   , INTENT(INOUT) :: X(0:IDIM-1,1:NN)
  INTEGER, INTENT(INOUT) :: ITYPEP(1:NN)
  REAL   , INTENT(IN)    :: WALL_NOR(0:IDIM-1,0:WALL_N-1) 
  REAL   , INTENT(IN)    :: DIS_MIN

  ! SUBROUTINE_VARIABLE **********************************************
  INTEGER :: I,J,K,L
  REAL    :: DISTANCE
  INTEGER :: SWITCH

  ! GRID_AND_PARTICLE_INTERFERENCE_CHECK *****************************
  DO I=1,NUMP
    DO J=0,WALL_N-1
      DISTANCE = PERP_2D(X(0,I),X(1,I),                   &
      &                  WALL_NOR(0,J),WALL_NOR(1,J),     &
      &                  WALL_ANC(0,J,0),WALL_ANC(1,J,0), &
      &                  WALL_ANC(0,J,1),WALL_ANC(1,J,1))

      SWITCH = FRONT_OR_BACK                     &
      &        (X(0,I),X(1,I),                   &
      &         WALL_NOR(0,K),WALL_NOR(1,K),     &
      &         WALL_ANC(0,K,0),WALL_ANC(1,K,0), &
      &         WALL_ANC(0,K,1),WALL_ANC(1,K,1))

      ! BEHIND_THE_WALL **********************************************
      IF(SWITCH .EQ. 0)THEN
        ITYPEP(I) = TYPE(0)
        DO L=0,IDIM-1
          X(L,I) = 0.0
        END DO
      END IF

      ! DISTANCE_TO_NEAREST__WALL ************************************
      IF(DISTANCE .LT. DIS*DIS_MIN*0.5)THEN
        ITYPEP(I) = TYPE(0)
        DO L=0,IDIM-1
          X(L,I) = 0.0
        END DO
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE GRID_AND_PARTICLE

!*********************************************************************
!*** CALCLATION_OF_GRID_INFORMATION                                ***
!*********************************************************************
SUBROUTINE GRID_INFORMATION                             &
&          (IDIM,WALL_N,WALL_ANC,WALL_OUT,WALL_NOR,DIS, &
&           GRID_DIST,NEAR_WALL,GRID,GRID_NUM,KER_C)
  IMPLICIT NONE
  ! MAINROUTINE_VARIABLE *********************************************
  INTEGER, INTENT(IN)    :: WALL_N
  INTEGER, INTENT(IN)    :: IDIM
  INTEGER, INTENT(IN)    :: GRID_NUM
  REAL   , INTENT(IN)    :: DIS
  REAL   , INTENT(IN)    :: WALL_ANC(0:IDIM-1,0:WALL_N-1,0:IDIM-1)
  REAL   , INTENT(IN)    :: GRID(1:GRID_NUM,0:IDIM-1,0:1)
  REAL   , INTENT(IN)    :: WALL_OUT(0:IDIM-1,0:WALL_N-1)
  REAL   , INTENT(IN)    :: WALL_NOR(0:IDIM-1,0:WALL_N-1) 
  REAL   , INTENT(IN) :: KER_C(0:3)
  REAL   , ALLOCATABLE, INTENT(OUT) :: GRID_DIST(:,:)
  INTEGER, ALLOCATABLE, INTENT(OUT) :: NEAR_WALL(:,:)

  ! SUBROUTINE_VARIABLE **********************************************
  INTEGER :: I,J,K,L
  INTEGER :: SWITCH
  INTEGER :: K_MIN
  REAL    :: R_MIN
  REAL    :: GRID_MAX
  REAL, ALLOCATABLE :: G_POINT(:)

  ! ARRAY_ALLOCATION *************************************************
  ALLOCATE(GRID_DIST(1:GRID_NUM,0:4*(IDIM-1)-1))
  ALLOCATE(NEAR_WALL(1:GRID_NUM,0:4*(IDIM-1)-1))
  ALLOCATE(G_POINT(0:IDIM-1))

  ! INITIAL_SETTING **************************************************
  GRID_MAX = DIS*(MAXVAL(KER_C(0:3))+3.0)

  ! ARRAY_INITIALIZE *************************************************
  GRID_DIST(1:GRID_NUM,0:4*(IDIM-1)-1) = GRID_MAX
  NEAR_WALL(1:GRID_NUM,0:4*(IDIM-1)-1) = -1

  ! GRID_DIST_SETTING ************************************************
  DO I=1,GRID_NUM
    DO J=0,4*(IDIM-1)-1
      IF(J .EQ. 0)THEN
        G_POINT(0) = GRID(I,0,0)
        G_POINT(1) = GRID(I,1,0)
      END IF

      IF(J .EQ. 1)THEN
        G_POINT(0) = GRID(I,0,1)
        G_POINT(1) = GRID(I,1,0)
      END IF

      IF(J .EQ. 2)THEN
        G_POINT(0) = GRID(I,0,1)
        G_POINT(1) = GRID(I,1,1)
      END IF

      IF(J .EQ. 3)THEN
        G_POINT(0) = GRID(I,0,0)
        G_POINT(1) = GRID(I,1,1)
      END IF

      ! SEARCH_NEAREST_WALL ******************************************
      NEAR_WALL(I,J) = -1
      R_MIN = GRID_MAX 

      DO K=0,WALL_N-1
        R_MIN = FOOT_OF_PERP_LINE_2D              &
        &       (G_POINT(0),G_POINT(1),           &
        &        WALL_NOR(0,K),WALL_NOR(1,K),     &
        &        WALL_ANC(0,K,0),WALL_ANC(1,K,0), &
        &        WALL_ANC(0,K,1),WALL_ANC(1,K,1), &
        &        GRID_MAX)

        ! NEAREST_WALL_NUMBER ****************************************
        IF(R_MIN .LT. ABS(GRID_DIST(I,J)))THEN
          NEAR_WALL(I,J) = K
          GRID_DIST(I,J) = R_MIN
        END IF
      END DO

      ! GRID_DIST **************************************************
      GRID_DIST(I,J) = GRID_MAX

      DO K=0,WALL_N-1
        R_MIN = PERP_2D                           &
        &       (G_POINT(0),G_POINT(1),           &
        &        WALL_NOR(0,K),WALL_NOR(1,K),     &
        &        WALL_ANC(0,K,0),WALL_ANC(1,K,0), &
        &        WALL_ANC(0,K,1),WALL_ANC(1,K,1))

        SWITCH = FRONT_OR_BACK                     &
        &        (G_POINT(0),G_POINT(1),           &
        &         WALL_NOR(0,K),WALL_NOR(1,K),     &
        &         WALL_ANC(0,K,0),WALL_ANC(1,K,0), &
        &         WALL_ANC(0,K,1),WALL_ANC(1,K,1))

        ! IF_GRID_POINT_IS_BEHIND_THE_WAL ****************************
        IF(SWITCH .EQ. 0)THEN
          R_MIN = -R_MIN
        END IF

        GRID_DIST(I,J) = MIN(GRID_DIST(I,J),R_MIN)
      END DO

      GRID_DIST(I,J) = MIN(GRID_MAX,GRID_DIST(I,J))

! DELETE *******************************************
      OPEN(888,FILE='./data/MPSini/CHECK.DAT')
      WRITE(888,*)G_POINT(0),G_POINT(1),GRID_DIST(I,J),NEAR_WALL(I,J)
! DELETE *******************************************

    END DO
  END DO

  ! ARRAY_DEALLOCATION ***********************************************
  DEALLOCATE(G_POINT)

  RETURN
END SUBROUTINE GRID_INFORMATION

!*********************************************************************
!*** FOOT_OF_PERPENDICULAR_LINE_FOR_2D                             ***
!*********************************************************************
REAL FUNCTION FOOT_OF_PERP_LINE_2D &
&             (X0,Y0,A,B,X1,Y1,X2,Y2,MAX)
  IMPLICIT NONE
  ! MAINROUTINE_VARIABLE *********************************************
  REAL :: X0,Y0
  REAL :: X1,X2
  REAL :: Y1,Y2
  REAL :: A,B
  REAL :: MAX

  ! FUNCTION_VARIABLE ************************************************
  REAL    :: C
  REAL    :: R_REF
  REAL    :: X_SUB,Y_SUB
  REAL    :: X_SUB1,Y_SUB1
  REAL    :: X_SUB2,Y_SUB2
  INTEGER :: I,J,K
  REAL    :: CP_ALL
  REAL    :: CP1,CP2,CP3
  INTEGER :: SWITCH
  REAL   , ALLOCATABLE :: D0(:,:)
  REAL    :: R_REF0,R_REF1,R_REF2
  REAL    :: VEC0,VEC1,VEC2

  ! ARRAY_ALLOCATION *************************************************
  ALLOCATE(D0(0:1,0:2))

  ! PARPENDICULAR_LENGTH *********************************************
  C = -(A*X1+B*Y1)

  IF(A .EQ. 0.0)THEN
    IF(B .EQ. 0.0)THEN
      WRITE(900,*)'ERROR IN FUNCTION PERP'
      WRITE(900,*)'A ,B AND C ARE ZERO'
    END IF
  END IF

  R_REF = ABS(A*X0+B*Y0+C)/SQRT(A**2+B**2)

  ! VARIABLE_INITIALIZING ********************************************
  SWITCH  = 0

  ! FOOT_OF_PERPENDICULAR_LINE ***************************************
  ! POSITION_OF_FOOT
  X_SUB1 = X0+R_REF*A/SQRT(A**2+B**2)
  Y_SUB1 = Y0+R_REF*B/SQRT(A**2+B**2)
  X_SUB2 = X0-R_REF*A/SQRT(A**2+B**2)
  Y_SUB2 = Y0-R_REF*B/SQRT(A**2+B**2)

  R_REF1 = ABS(A*X_SUB1+B*Y_SUB1+C)/SQRT(A**2+B**2)
  R_REF2 = ABS(A*X_SUB2+B*Y_SUB2+C)/SQRT(A**2+B**2)

  IF(R_REF1 .LE. R_REF2)THEN
    X_SUB = X_SUB1
    Y_SUB = Y_SUB1
    R_REF0 = R_REF1
  ELSE
    X_SUB = X_SUB2
    Y_SUB = Y_SUB2
    R_REF0 = R_REF2
  END IF

  ! VECTORS_BETWEEN_ACHORS
  D0(0,0) = X2-X1
  D0(1,0) = Y2-Y1

  D0(0,1) = X_SUB-X1
  D0(1,1) = Y_SUB-Y1

  D0(0,2) = X_SUB-X2
  D0(1,2) = Y_SUB-Y2

  ! MAGNITUDE_OF_VECTORS
  VEC0 = SQRT(D0(0,0)**2+D0(1,0)**2)
  VEC1 = SQRT(D0(0,1)**2+D0(1,1)**2)
  VEC1 = SQRT(D0(0,2)**2+D0(1,2)**2)

  FOOT_OF_PERP_LINE_2D = MAX

  IF(VEC0 .GE. VEC1)THEN
    IF(VEC0 .GE. VEC2)THEN
      FOOT_OF_PERP_LINE_2D = MIN(R_REF,MAX)
    END IF
  END IF

  ! ARRAY_DEALLOCATION ***********************************************
  DEALLOCATE(D0)

END FUNCTION FOOT_OF_PERP_LINE_2D

!*********************************************************************
!*** INSIDE_OF_OUTSIDE                                             ***
!*********************************************************************
INTEGER FUNCTION FRONT_OR_BACK         &
&                (X0,Y0,A,B,X1,Y1,X2,Y2)
  IMPLICIT NONE
  ! MAINROUTINE_VARIABLE *********************************************
  REAL :: X0,Y0
  REAL :: A,B
  REAL :: X1,Y1
  REAL :: X2,Y2

  ! FUNCTION_VARIABLE ************************************************
  REAL :: C
  REAL :: R_REF
  REAL :: X_SUB ,Y_SUB
  REAL :: X_SUB1,Y_SUB1
  REAL :: X_SUB2,Y_SUB2
  REAL :: R_REF0,R_REF1,R_REF2
  REAL :: IP 
  REAL, ALLOCATABLE :: D0(:)

  ! ARRAY_ALLOCATION *************************************************
  ALLOCATE(D0(0:2))

  ! INSIDE_OR_OUTSIDE ************************************************
  C = -(A*X1+B*Y1)

  ! PARPENDICULAR_LENGTH
  R_REF = ABS(A*X0+B*Y0+C)/SQRT(A**2+B**2)

  ! FOOT_OF_PERPENDICULAR_LINE
  X_SUB1 = X0+R_REF*A/SQRT(A**2+B**2)
  Y_SUB1 = Y0+R_REF*B/SQRT(A**2+B**2)
  X_SUB2 = X0-R_REF*A/SQRT(A**2+B**2)
  Y_SUB2 = Y0-R_REF*B/SQRT(A**2+B**2)

  R_REF1 = ABS(A*X_SUB1+B*Y_SUB1+C)/SQRT(A**2+B**2)
  R_REF2 = ABS(A*X_SUB2+B*Y_SUB2+C)/SQRT(A**2+B**2)

  IF(R_REF1 .LE. R_REF2)THEN
    X_SUB = X_SUB1
    Y_SUB = Y_SUB1
    R_REF0 = R_REF1
  ELSE
    X_SUB = X_SUB2
    Y_SUB = Y_SUB2
    R_REF0 = R_REF2
  END IF
    
  ! VECTORS(FOOT_TO_GRID_POINT)
  D0(0) = X0-X_SUB
  D0(1) = Y0-Y_SUB
  
  ! INNER_PRODUCT
  IP = D0(0)*A+D0(1)*B  

  ! IF_GRID_POINT_IS_ON_THE_FRONT_OR_BACK_OF_WALL
  FRONT_OR_BACK = 0

  IF(IP .LT. 0.0)THEN
    FRONT_OR_BACK = 1
  END IF

  ! ARRAY_DEALLOCATION ***********************************************
  DEALLOCATE(D0)

  RETURN
END FUNCTION FRONT_OR_BACK

!*********************************************************************
!*** PERPENDICULAR_LENGTH_FOR_2D                                   ***
!*********************************************************************
REAL FUNCTION PERP_2D(X0,Y0,A,B,X1,Y1,X2,Y2)
  IMPLICIT NONE
  ! MAINROUTINE_VARIABLE *********************************************
  REAL :: X0,Y0
  REAL :: X1,X2
  REAL :: Y1,Y2
  REAL :: A,B

  ! FUNCTION_VARIABLE ************************************************
  REAL :: C

  ! PARPENDICULAR_LENGTH *********************************************
  C = -(A*X1+B*Y1)

  IF(A .EQ. 0.0)THEN
    IF(B .EQ. 0.0)THEN
      WRITE(900,*)'ERROR IN FUNCTION PERP'
      WRITE(900,*)'A ,B AND C ARE ZERO'
    END IF
  END IF

  PERP_2D = ABS(A*X0+B*Y0+C)/SQRT(A**2+B**2)

END FUNCTION PERP_2D

!*********************************************************************
!*** CHECK_IF_THE_GRID_POINT_IS_INSIDE_CALCULATION_DOMAIN          ***
!*********************************************************************
INTEGER FUNCTION INTERSECTION_2D(X1,Y1,X2,Y2,X3,Y3,X4,Y4)
  IMPLICIT NONE
  ! MAINROUTINE_VARIABLE *********************************************
  REAL :: X1,X2,X3,X4
  REAL :: Y1,Y2,Y3,Y4
  REAL :: R,S

  ! SUBROUTINE_VARIABLE **********************************************
  IF((X2-X1)*(Y4-Y3)-(Y2-Y1)*(X4-X3) .EQ. 0.0)THEN
    INTERSECTION_2D = 0
    RETURN
  END IF

  ! CALCULATION_COEFFICIENT ******************************************
  R = ((Y4-Y3)*(X3-X1)-(X4-X3)*(Y3-Y1))/((X2-X1)*(Y4-Y3)-(Y2-Y1)*(X4-X3))
  S = ((Y2-Y1)*(X3-X1)-(X2-X1)*(Y3-Y1))/((X2-X1)*(Y4-Y3)-(Y2-Y1)*(X4-X3))

  ! IF_GRID_POINT_IS_INSIDE_CALCULATION_DOMAIN ***********************
  INTERSECTION_2D = 0

  IF(0.0 .LT. R)THEN
    IF(R .LE. 1.0)THEN
      IF(0.0 .LT. S)THEN
        IF(S .LE. 1.0)THEN
          INTERSECTION_2D = 1
        END IF
      END IF
    END IF
  END IF

END FUNCTION INTERSECTION_2D

!*********************************************************************
!*** CALCULATION_DISTANCE                                          ***
!*********************************************************************
REAL FUNCTION DIST_2D(X1,Y1,X2,Y2)
  IMPLICIT NONE
  ! MAINTOURINE_VARIABLE *********************************************
  REAL :: X1,X2
  REAL :: Y1,Y2

  ! CALCULATION_DISTANCE *********************************************
  DIST_2D = (X2-X1)**2+(Y2-Y1)**2
  DIST_2D = SQRT(DIST_2D)
END FUNCTION DIST_2D

!*********************************************************************
!*** GRID_SETTING_FOR_2D                                           ***
!*********************************************************************
SUBROUTINE GRID_SETTING_2D         &
&          (NN,IDIM,NUMP,DIS,X,    &
&           GRID_N,GRID_NUM,KER_C, &
&           GRID,DOM_SIZE)
  IMPLICIT NONE
  ! MAINROUTINE_VARIABLE *********************************************
  INTEGER, INTENT(IN)              :: NN,IDIM,NUMP
  REAL   , INTENT(IN)              :: DIS
  REAL   , INTENT(IN)              :: X(0:IDIM-1,1:NN)
  REAL   , INTENT(IN)              :: KER_C(0:3)
  INTEGER, INTENT(INOUT)           :: GRID_NUM
  REAL   , INTENT(OUT), ALLOCATABLE:: GRID(:,:,:)
  INTEGER, INTENT(INOUT)           :: GRID_N(0:IDIM-1)
  REAL   , INTENT(INOUT)           :: DOM_SIZE(0:IDIM-1,0:1)

  ! SUBROUTINE_VARIBALE **********************************************
  INTEGER              :: I,J,L
  INTEGER              :: Z_NUM
  INTEGER              :: Z_INC
  REAL                 :: KERNEL_MAX
  REAL                 :: MARGE
  INTEGER, ALLOCATABLE :: SW(:)
  REAL   , ALLOCATABLE :: DIV(:)
  REAL   , ALLOCATABLE :: X_RANGE(:)
  INTEGER, ALLOCATABLE :: MOD_X(:)
  real   , allocatable :: grid_x(:,:)
  real   , allocatable :: grid_y(:,:)

  ! CALCULATION_START ************************************************
  WRITE(*,'(A)')"********** SUBROUTINE GRID_SETTING_2D START"

  ! ARRAY_ALLOCATION *************************************************
  ALLOCATE(SW(0:IDIM-1))
  ALLOCATE(DIV(0:IDIM-1))
  ALLOCATE(X_RANGE(0:IDIM-1))
  ALLOCATE(MOD_X(0:2))

  ! ARRAY_INITIALIZING ***********************************************
  DO I=0,IDIM-1
    GRID_N(I)  = 1
    SW(I)      = 0
    DIV(I)     = 0.0
    X_RANGE(I) = 0.0
    MOD_X(I)   = 0
  END DO

  ! MARGIN_SETTING ***************************************************
!  MARGE = 2.0*DIS
  MARGE = 0.0
  ! MAXIMUM_KERNEL_RADIUS ********************************************
  KERNEL_MAX = MAXVAL(KER_C(0:3))*DIS

  ! CALCULATION_GRID_NUMBER ******************************************
  DO I=0,IDIM-1
    ! CALCULATION_DOMAIN_SIZE
    DOM_SIZE(I,0) = MIN(DOM_SIZE(I,0),MINVAL(X(I,1:NN)))-MARGE
    DOM_SIZE(I,1) = MAX(DOM_SIZE(I,1),MAXVAL(X(I,1:NN)))+MARGE
    X_RANGE(I) = ABS(DOM_SIZE(I,1)-DOM_SIZE(I,0))

    ! GRID_NUMBER_CHECK
    DO WHILE(SW(I) == 0)
      DIV(I) = X_RANGE(I)/REAL(GRID_N(I))
      IF(DIV(I) > KERNEL_MAX)THEN
        GRID_N(I) = GRID_N(I)+1
      ELSE
        GRID_N(I) = GRID_N(I)-1
        DIV(I) = X_RANGE(I)/REAL(GRID_N(I))
        SW(I) = 1
      END IF
    END DO

    ! CHECK_GRID_DISTANCE
    IF((DIV(I) < KERNEL_MAX))THEN
      WRITE(900,'(A   )')"ERROR IN SUBROUTINE GRID_NUMBER"
      WRITE(900,'(A,I2)')"GRID_N : TOO LARGE",I
      STOP
    END IF
  END DO

  ! CALCULATION_GRID_NUMBER ******************************************
  GRID_NUM = 1

  DO I=0,IDIM-1
    GRID_NUM = GRID_NUM*GRID_N(I)
  END DO

  ! OUTPUT_GRID_NUMBER ***********************************************
  DO I=0,IDIM-1
    WRITE(900,'(A,I1,A,I16)')"GRID_N(",I,")       =",GRID_N(I)
  END DO

  WRITE(900,'(A,I16)')"GRID_NUM        =",GRID_NUM
  WRITE(900,'(A)')"---------------------------------"

  ! ARRAY_ALLOCATION *************************************************
  ALLOCATE(GRID(1:GRID_NUM,0:IDIM-1,0:1))
  allocate(grid_x(0:grid_n(0),0:grid_n(1)))
  allocate(grid_y(0:grid_n(0),0:grid_n(1)))

  ! GRID_INITIALIZING ************************************************
  DO I=1,GRID_NUM
    DO J=0,IDIM-1
      DO L=0,1
        GRID(I,J,L) = 0.0
      END DO
    END DO
  END DO

  ! DOMAIN_SIZE_SETTING **********************************************
  DO I=0,IDIM-1
    DIV(I) = ABS(X_RANGE(I))/REAL(GRID_N(I))
  END DO

  ! GRID_POSITION ****************************************************
  DO I=1,GRID_NUM
    MOD_X(0) = MOD(I,GRID_N(0))
    MOD_X(1) = (I-MOD(I,GRID_N(0)))/GRID_N(0)

    IF(MOD_X(0) .EQ. 0)THEN
     MOD_X(0) = GRID_N(0)
     MOD_X(1) = MOD_X(1)-1
    END IF

    MOD_X(0) = MOD_X(0)-1

    ! CALCULATION_GRID_POSITION
    DO L=0,IDIM-1
      GRID(I,L,0) = DOM_SIZE(L,0)+REAL(MOD_X(L)  )*DIV(L)
      GRID(I,L,1) = DOM_SIZE(L,0)+REAL(MOD_X(L)+1)*DIV(L)
    END DO
  END DO

  ! GRID_POSITION_CHECK **********************************************
  DO I=0,IDIM-1
    ! MAXVAL
    IF(MAXVAL(X(I,1:NUMP)) > MAXVAL(GRID(1:GRID_NUM,I,0:1)))THEN
      WRITE(900,'(A   )')"ERROR IN SUBROUTINE GRID_SETTING_2D"
      WRITE(900,'(A,I2)')"PARTILCLE : OUT OF GRID(MAXVAL)",I
      WRITE(900,'(F10.5,F10.5)')MAXVAL(X(I,1:NUMP)), &
      &                       MAXVAL(GRID(1:GRID_NUM,I,0:1))
    END IF

    ! MINVAL
    IF(MINVAL(X(I,1:NUMP)) < MINVAL(GRID(1:GRID_NUM,I,0:1)))THEN
      WRITE(900,'(A   )')"ERROR IN SUBROUTINE GRID_SETTING_2D"
      WRITE(900,'(A,I2)')"PARTILCLE : OUT OF GRID(MINVAL)",I
      WRITE(900,'(F10.5,F10.5)')MINVAL(X(I,1:NUMP)), &
      &                       MINVAL(GRID(1:GRID_NUM,I,0:1))
    END IF
  END DO
  
  do i = 0,grid_n(0)
    do j = 0,grid_n(1)
        grid_x(i,j) = dom_size(0,0) + div(0)*real(i)
        grid_y(i,j) = dom_size(1,0) + div(1)*real(j)
    end do
  end do
  
  open(1, file = './data/MPSini/grid_visualize.g' , form = 'unformatted', status = 'replace')
  
    write(1)grid_n(0)+1, grid_n(1)+1
    
    write(1) ((grid_x(i,j), i = 0, grid_n(0)), j = 0, grid_n(1)), &
    &        ((grid_y(i,j), i = 0, grid_n(0)), j = 0, grid_n(1))
    
  close(1)
  
  ! ARRAY_DEALLOCATION ***********************************************
  DEALLOCATE(SW)
  DEALLOCATE(DIV)
  DEALLOCATE(X_RANGE)
  DEALLOCATE(MOD_X)

  ! CALCULATION_END **************************************************
  WRITE(*,'(A)')"********** SUBROUTINE GRID_SETTING_2D END"
  WRITE(*,'(A)')""

  RETURN
END SUBROUTINE GRID_SETTING_2D

!*********************************************************************
!*** BOX_POSITIONING_FOR_2D                                        ***
!*********************************************************************
SUBROUTINE BOX_POSITIONING_2D   &
&          (IDIM,NN,NUMP,X,DIS, &
&           OUT_TYPE,P_TYPE,    &
&           WALL_FIX,OUT_FIX,   &
&           HOR0,HOR1,          &
&           VER0,VER1,          &
&           INTERVAL,ITYPEP,DIS_MIN,TEMP,T)
  IMPLICIT NONE
  ! MAINROUTINE_VARIABLE *********************************************
  INTEGER, INTENT(IN)    :: IDIM
  INTEGER, INTENT(IN)    :: NN
  REAL   , INTENT(IN)    :: INTERVAL
  REAL   , INTENT(IN)    :: DIS
  REAL   , INTENT(IN)    :: HOR0
  REAL   , INTENT(IN)    :: HOR1
  REAL   , INTENT(IN)    :: VER0
  REAL   , INTENT(IN)    :: VER1
  INTEGER, INTENT(IN)    :: WALL_FIX
  INTEGER, INTENT(IN)    :: OUT_FIX
  INTEGER, INTENT(IN)    :: P_TYPE
  REAL   , INTENT(IN)    :: DIS_MIN
  REAL   , INTENT(IN)    :: TEMP
  INTEGER, INTENT(INOUT) :: ITYPEP(1:NN)
  INTEGER, INTENT(INOUT) :: NUMP
  REAL   , INTENT(INOUT) :: X(0:IDIM-1,1:NN)
  INTEGER, INTENT(INOUT) :: OUT_TYPE(1:NN)
  REAL   , INTENT(INOUT) :: T(1:NN)

  ! SUBROUTINE_VARIABLE **********************************************
  INTEGER         :: I,J,K,N
  INTEGER         :: NUM(0:1)
  INTEGER         :: POS_SW(0:1)
  REAL            :: X_POS(0:1)
  REAL            :: RR,RR_MIN
  REAL            :: X_MIN(0:1)
  REAL            :: X_MAX(0:1)

  ! CALCULATION_START ************************************************
  WRITE(*,'(A)')"********** SUBROUTINE BOX_POSITIONING_2D START"

  ! DATA_COPY ********************************************************
  X_MIN(0) = MIN(HOR0,HOR1)
  X_MAX(0) = MAX(HOR0,HOR1)
  X_MIN(1) = MIN(VER0,VER1)
  X_MAX(1) = MAX(VER0,VER1)

  ! REFERENCE_POSITION ***********************************************
  ! X_DIRECTION
  IF(X_MIN(0) == HOR0)THEN
    POS_SW(0) = 0
  ELSE
    POS_SW(0) = 1
  END IF

  ! Y_DIRECTION
  IF(X_MIN(1) == VER0)THEN
    POS_SW(1) = 0
  ELSE
    POS_SW(1) = 1
  END IF

  ! VARIABLE_CHECK ***************************************************
  DO K=0,1
    IF(X_MAX(K) <= X_MIN(K))THEN
      WRITE(900,'(A   )')"ERROR IN SUBROUTINE BOX_POSITIONING_2D"
      WRITE(900,'(A,I2)')"X_MIN IS GREATER THAN X_MAX",K
      STOP
    END IF
  END DO

  IF(INTERVAL <= 0.0)THEN
    WRITE(900,'(A)')"ERROR IN SUBROUTINE BOX_POSITIONING_2D"
    WRITE(900,'(A)')"INTERVAL IS NEGATIVE OR ZERO"
    STOP
  END IF

  ! NUMBER_OF_SEARCHING_POINT ****************************************
  DO K=0,1
    NUM(K) = INT(ABS(X_MAX(K)-X_MIN(K))/INTERVAL)
  END DO

  ! PARTICLE_POSITIONING *********************************************
  DO I=0,NUM(0)
    DO J=0,NUM(1)
      ! SEARCHING_POSITION
      IF(POS_SW(0) == 0)THEN
        X_POS(0) = X_MIN(0)+REAL(I)*INTERVAL
      ELSE
        X_POS(0) = X_MAX(0)-REAL(I)*INTERVAL
      END IF

      IF(POS_SW(1) == 0)THEN
        X_POS(1) = X_MIN(1)+REAL(J)*INTERVAL
      ELSE
        X_POS(1) = X_MAX(1)-REAL(J)*INTERVAL
      END IF

      ! INTERPARTICLE_DISTANCE
      RR_MIN = 1.0E10
      DO N=1,NUMP
        RR     = DIST_2D(X(0,N),X(1,N),X_POS(0),X_POS(1))
        RR_MIN = MIN(RR,RR_MIN)
      END DO

      ! PARTICLE_POSITIONING
      IF(RR_MIN > DIS_MIN*DIS)THEN
        NUMP = NUMP+1
        DO K=0,1
          X(K,NUMP) = X_POS(K)
          ITYPEP(NUMP) = P_TYPE
          T(NUMP)      = TEMP
        END DO
        OUT_TYPE(NUMP)  = OUT_FIX
      END IF
    END DO
  END DO

  ! CALCULATION_END **************************************************
  WRITE(*,'(A)')"********** SUBROUTINE BOX_POSITIONING_2D END"
  WRITE(*,'(A)')""

  RETURN
END SUBROUTINE BOX_POSITIONING_2D

!*********************************************************************
!*** CIRCLE_POSITIONING_FOR_2D                                     ***
!*********************************************************************
SUBROUTINE CIRCLE_POSITIONING_2D &
&          (IDIM,NN,NUMP,X,DIS,  &
&           OUT_TYPE,P_TYPE,     &
&           WALL_FIX,OUT_FIX,    &
&           CENT0,CENT1,RADIUS,  &
&           INTERVAL,ITYPEP,DIS_MIN,TEMP,T)
  IMPLICIT NONE
  ! MAINROUTINE_VARIABLE *********************************************
  INTEGER, INTENT(IN)    :: IDIM
  INTEGER, INTENT(IN)    :: NN
  REAL   , INTENT(IN)    :: INTERVAL
  REAL   , INTENT(IN)    :: DIS
  REAL   , INTENT(IN)    :: CENT0
  REAL   , INTENT(IN)    :: CENT1
  REAL   , INTENT(IN)    :: RADIUS
  INTEGER, INTENT(IN)    :: WALL_FIX
  INTEGER, INTENT(IN)    :: OUT_FIX
  INTEGER, INTENT(IN)    :: P_TYPE
  REAL   , INTENT(IN)    :: DIS_MIN
  REAL   , INTENT(IN)    :: TEMP
  INTEGER, INTENT(INOUT) :: ITYPEP(1:NN)
  INTEGER, INTENT(INOUT) :: NUMP
  REAL   , INTENT(INOUT) :: X(0:IDIM-1,1:NN)
  INTEGER, INTENT(INOUT) :: OUT_TYPE(1:NN)
  REAL   , INTENT(INOUT) :: T(1:NN)

  ! SUBROUTINE_VARIABLE **********************************************
  INTEGER         :: I,J,K,N
  INTEGER         :: NUM(0:1)
  REAL            :: X_POS(0:1)
  REAL            :: RR,RR_MIN
  REAL            :: INTER_R
  REAL            :: R_POS
  REAL            :: T_POS
  REAL, PARAMETER :: PI = ACOS(0.0)*2.0

  ! CALCULATION_START ************************************************
  WRITE(*,'(A)')"********** SUBROUTINE CIRCLE_POSITIONING_2D START"

  ! NUMBER_OF_SEARCHING_POINT ****************************************
  NUM(0) = 360
  NUM(1) = INT(RADIUS/INTERVAL)

  INTER_R = 360.0/REAL(NUM(0))*2.0*PI/360.0

  ! PARTICLE_POSITIONING *********************************************
  DO J=NUM(1)-1,0,-1
    DO I=0,NUM(0)-1
      ! SEARCHING_POSITION
      T_POS = REAL(I)*INTER_R
      R_POS = REAL(J)*INTERVAL
      X_POS(0) = CENT0+R_POS*COS(T_POS)
      X_POS(1) = CENT1+R_POS*SIN(T_POS)

      ! INTERPARTICLE_DISTANCE
      RR_MIN = 1.0E10
      DO N=1,NUMP
        RR     = DIST_2D(X(0,N),X(1,N),X_POS(0),X_POS(1))
        RR_MIN = MIN(RR,RR_MIN)
      END DO

      ! PARTICLE_POSITIONING
      IF(RR_MIN > DIS_MIN*DIS)THEN
        NUMP = NUMP+1
        DO K=0,1
          X(K,NUMP) = X_POS(K)
          ITYPEP(NUMP) = P_TYPE
          T(NUMP)      = TEMP
        END DO
        OUT_TYPE(NUMP)  = OUT_FIX
      END IF
    END DO
  END DO

  ! CALCULATION_END **************************************************
  WRITE(*,'(A)')"********** SUBROUTINE CIRCLE_POSITIONING_2D END"
  WRITE(*,'(A)')""

  RETURN
END SUBROUTINE CIRCLE_POSITIONING_2D


END MODULE MOD_IS_DIM
