      PROGRAM GOODING
!-----MAIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /GIVEN/ OX1, OY1, OZ1, OX2, OY2, OZ2, OX3, OY3,
     1          OZ3, EL1, EM1, EN1, EL2, EM2, EN2, EL3, EM3, EN3, T12,
     1          T13, GM
      COMMON /USEFUL/ AL, Q, EI, BOM, R1, R2, R3
      COMMON /CNTROL/ LOBDEP, HN, MAXIT, PDINC, CRIVAL, CRITGM,
     1          NFAIL, ITNUM, CR, CRIT, AXRTIO
      COMMON /KNWNSL/ IKN, RO1KN(9), RO3KN(9), R1KNSQ(9),
     1          R3KNSQ(9)

      LOGICAL LOBDEP

!-----OPEN THE INpUT FILE AND PREPARE THE OUTPUT FILE
      OPEN(UNIT=5, FILE='in.data', STATUS='old')
      OPEN(UNIT=6, FILE='out.txt', STATUS='unknown')

      !READ IN ALL OF OUR DATA
      LOBDEP = .TRUE.
      PRINT *, "LOBDEP = ", LOBDEP

      READ(5,*) OX1, OY1, OZ1
      PRINT *, "R1: ", OX1, OY1, OZ1

      READ(5,*) OX2, OY2, OZ2
      PRINT *, "R2: ", OX2, OY2, OZ2

      READ(5,*) OX3, OY3, OZ3
      PRINT *, "R3: ", OX3, OY3, OZ3

      READ(5,*) EL1, EM1, EN1
      PRINT *, "LAMBDA1: ", EL1, EM1, EN1

      READ(5,*) EL2, EM2, EN2
      PRINT *, "LAMBDA2: ", EL2, EM2, EN2

      READ(5,*) EL3, EM3, EN3
      PRINT *, "LAMBDA3: ", EL3, EM3, EN3

      READ(5,*) T12
      PRINT *, "T12: ", T12

      READ(5,*) T13
      PRINT *, "T13: ", T13

      READ(5,*) GM
      PRINT *, "GM: ", GM

      READ(5,*) HN
      PRINT *, "HN: ", HN

      READ(5,*) MAXIT
      PRINT *, "MAXIT: ", MAXIT

      READ(5,*) PDINC
      PRINT *, "PDINC: ", PDINC

      READ(5,*) CRIVAL
      PRINT *, "CRIVAL: ", CRIVAL
      
      READ(5,*) CRITGM
      PRINT *, "CRITGM: ", CRITGM

      READ(5,*) NHREV
      PRINT *, "NHREV = ", NHREV

      READ(5,*) RO1
      PRINT *, "INITIAL RO1: ", RO1

      READ(5,*) RO3
      PRINT *, "INITIAL RO3: ", RO3

      READ(5,*) IND
      PRINT *, "IND: ", IND

      CALL OBS3LS(NHREV, IND, RO1, RO3)

      PRINT *, "FINAL RO1 = ", RO1
      PRINT *, "FINAL RO3 = ", RO3
      PRINT *, "NFAIL = ", NFAIL
      PRINT *, "EOF"

      RO1 = RO1 - 2
      RO3 = RO3 - 2

      CALL OBS3LS(NHREV, IND, RO1, RO3)

      PRINT *, "FINAL RO1 = ", RO1
      PRINT *, "FINAL RO3 = ", RO3
      PRINT *, "NFAIL = ", NFAIL
      PRINT *, "EOF"

      RO1 = RO1 - 2
      RO3 = RO3 - 2

      CALL OBS3LS(NHREV, IND, RO1, RO3)

      PRINT *, "FINAL RO1 = ", RO1
      PRINT *, "FINAL RO3 = ", RO3
      PRINT *, "NFAIL = ", NFAIL
      PRINT *, "EOF"

      RO1 = RO1 - 2
      RO3 = RO3 - 2
      
      CALL OBS3LS(NHREV, IND, RO1, RO3)

      PRINT *, "FINAL RO1 = ", RO1
      PRINT *, "FINAL RO3 = ", RO3
      PRINT *, "NFAIL = ", NFAIL
      PRINT *, "EOF"

      STOP
      END

      SUBROUTINE TLAMB (M, Q, QSQFM1, X, N, T, DT, D2T, D3T)

              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              LOGICAL LM1, L1, L2, L3

              PI = 4.0D0*DATAN(1.0D0)
              SW = 0.4D0

              LM1 = N.EQ.-1
              L1 = N.GE.1
              L2 = N.GE.2
              L3 = N.EQ.3
              QSQ = Q*Q
              XSQ = X*X
              U = (1D0 - X)*(1D0 + X)

              IF (.NOT.LM1) THEN
!---------------------NEEDED IF SERIES AND OTHERWISE USEFUL WHEN Z=?
                      DT = 0D0
                      D2T = 0D0
                      D3T = 0D0
              END IF

              IF (LM1 .OR. M.GT.0 .OR. X.LT.0D0 .OR. DABS(U).GT.SW) THEN

!---------------------DIRECT COMPUTATION (NOT SERIES)
                      Y = DSQRT(DABS(U))
                      Z = DSQRT(QSQFM1 + QSQ*XSQ)
                      QX = Q*X
                      IF (QX.LE.0D0) THEN
                              A = Z - QX
                              B = Q*Z - X
                      END IF

                      IF (QX.LT.0D0 .AND. LM1) THEN
                              AA = QSQFM1/A
                              BB = QSQFM1*(QSQ*U - XSQ)/B
                      END IF

                      IF (QX.EQ.0D0.AND.LM1 .OR. QX.GT.0D0) THEN

                              AA = Z + QX
                              BB = Q*Z + X
                      END IF

                      IF (QX.GT.0D0) THEN
                              A = QSQFM1/AA
                              B = QSQFM1*(QSQ*U - XSQ)/BB
                      END IF

                      IF (.NOT.LM1) THEN
                              IF (QX*U.GE.0D0) THEN
                                      G = X*Z + Q*U
                              ELSE
                                      G = (XSQ - QSQ*U)/(X*Z - Q*U)
                              END IF

                              F = A*Y
                              
                              IF (X.LE.1D0) THEN
                                      T = M*PI + DATAN2(F, G)
                              ELSE
                                      IF (F.GT.SW) THEN
                                              T = DLOG(F + G)
                                      ELSE
                                              FG1 = F/(G + 1D0)
                                              TERM = 2D0*FG1
                                              FG1SQ = FG1*FG1
                                              T = TERM
                                              TWOI1 = 1D0
    1                                         TWOI1 = TWOI1 + 2D0
                                              TERM = TERM*FG1SQ
                                              TOLD = T
                                              T = T + TERM/TWOI1
                                              IF (T.NE.TOLD) GO TO 1
!--------------------------------(CONTINUE LOOPING FOR INVERSE TANH)
                                      END IF
                              END IF

                              T = 2D0*(T/Y + B)/U
                              IF (L1 .AND. Z.NE.0D0) THEN
                                      QZ = Q/Z
                                      QZ2 = QZ*QZ
                                      QZ = QZ*QZ2
                                      DT = (3D0*X*T - 4D0*(A +
     1                                      QX*QSQFM1)/Z)/U
                                      IF (L2) D2T = (3D0*T + 5D0*X*DT +
     1                                               4D0*QZ*QSQFM1)/U
                                      IF (L3) D3T = (8D0*DT + 7D0*X*D2T
     1                                   - 12D0*QZ*QZ2*X*QSQFM1)/U
                              END IF
                      ELSE
                              DT = B
                              D2T = BB
                              D3T = AA
                      END IF
              ELSE
!---------------------COMPUTE BY SERIES
                      U0I = 1D0
                      IF (L1) U1I = 1D0
                      IF (L2) U2I = 1D0
                      IF (L3) U3I = 1D0
                      TERM = 4D0
                      TQ = Q*QSQFM1
                      I = 0
                      IF (Q.LT.5D-1) TQSUM = 1D0 - Q*QSQ
                      IF (Q.GE.5D-1) TQSUM = (1D0/(1D0 + Q) + Q)*QSQFM1
                      TTMOLD = TERM/3D0
                      T = TTMOLD*TQSUM
                      
!---------------------START OF LOOP
    2                 I = I + 1
                      P = I
                      U0I = U0I*U
                      IF (L1 .AND. I.GT.1) U1I = U1I*U
                      IF (L2 .AND. I.GT.2) U2I = U2I*U
                      IF (L3 .AND. I.GT.3) U3I = U3I*U
                      TERM = TERM*(P - 0.5D0)/P
                      TQ = TQ*QSQ
                      TQSUM = TQSUM + TQ
                      TOLD = T
                      TTERM = TERM/(2D0*P + 3D0)
                      TQTERM = TTERM*TQSUM
                      T = T - U0I*((1.5D0*P + 0.25D0)*TQTERM/(P*P -
     1                    0.25D0) - TTMOLD*TQ)
                      TTMOLD = TTERM
                      TQTERM = TQTERM*P
                      IF (L1) DT = DT + TQTERM*U1I
                      IF (L2) D2T = D2T + TQTERM*U2I*(P - 1D0)
                      IF (L3) D3T = D3T + TQTERM*U3I*(P - 1D0)*(P - 2D0)
                      IF (I.LT.N .OR. T.NE.TOLD) GO TO 2
                      
!---------------------(END OF LOOP)

                      IF (L3) D3T = 8D0*X*(1.5D0*D2T - XSQ*D3T)
                      IF (L2) D2T = 2D0*(2D0*XSQ*D2T - DT)
                      IF (L1) DT = -2D0*X*DT
                      T = T/XSQ
              END IF
              RETURN
              END SUBROUTINE TLAMB

      SUBROUTINE XLAMB (M, Q, QSQFM1, TIN, N, X, XPL)
              
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              PARAMETER (PI=4.0D0*DATAN(1.0D0), TOL=3D-7, C0=1.7D0)
              PARAMETER (C1=0.5D0, C2=0.03D0, C3=0.15D0, C41=1D0)
              PARAMETER (C42=0.24D0)
              
              D8RT(X) = DSQRT(DSQRT(DSQRT(X)))
              THR2 = DATAN2(QSQFM1, 2D0*Q)/PI
              
              IF (M.EQ.0) THEN
!-------------SINGLE-REV STARTER FROM T (AT X=0) & BILINEAR USUALLY
                      N = 1
                      
                      CALL TLAMB(M,Q,QSQFM1,0D0,0,T0,DT,D2T,D3T)

                      TDIFF = TIN - T0
                      IF (TDIFF.LE.0D0) THEN
                              X = T0*TDIFF/(-4D0*TIN)
!-----------------------------( -4 IS THE VALUE OF DT, FOR X=0)
                      ELSE
                              X = -TDIFF/(TDIFF + 4D0)
                              W = X + C0*DSQRT(2D0*(1D0 - THR2))
                              IF (W.LT.0D0)
     1                             X = X - DSQRT(D8RT(-W))*(X +
     1                                 DSQRT(TDIFF/(TDIFF + 1.5D0*T0)))
                              W = 4D0/(4D0 + TDIFF)
                              X = X*(1D0 + X*(C1*W - C2*X*DSQRT(W)))
                      END IF
              ELSE
!-------------WITH MULTIREVS, FIRST GET T(MIN) AS BASIS FOR STARTER
                      XM = 1D0/(1.5D0*(M + 5D-1)*PI)
                      IF (THR2.LT.5D-1) XM = D8RT(2D0*THR2)*XM
                      IF (THR2.GT.5D-1) XM = (2D0 - D8RT(2D0 -
     1                                       2D0*THR2))*XM
!---------------------STARTER FOR TMIN
                      DO 1 I=1,12

                      CALL TLAMB (M,Q,QSQFM1,XM,3,TMIN,DT,D2T,D3T)

                      IF (D2T.EQ.0D0) GO TO 2

                      XMOLD = XM
                      XM = XM - DT*D2T/(D2T*D2T - DT*D3T/2D0)
                      XTEST = DABS(XMOLD/XM - 1D0)
                      IF (XTEST.LE.TOL) GO TO 2
    1                 CONTINUE
                      N = -1
                      RETURN

!-------(BREAK OFF AND EXIT IF TMIN NOT LOCATED - SHOULD NEVER HAPPEN)
!-------NOW PROCEED FROM T(MIN) TO FULL STARTER

    2                 CONTINUE
                      TDIFFM = TIN - TMIN
                      IF (TDIFFM.LT.0D0) THEN
                              N = 0
                              RETURN

!-----------------------------(EXIT IF NO SOLUTION WITH THIS M)

                      ELSE IF (TDIFFM.EQ.0D0) THEN
                              X = XM
                              N = 1
                              RETURN

!----------------(EXIT IF UNIQUE SOLUTION ALREADY FROM X(TMIN))

                      ELSE
                              N = 3
                              IF (D2T.EQ.0D0) D2T = 6D0*M*PI
                              X = DSQRT(TDIFFM/(D2T/2D0 + TDIFFM/(1D0 -
     1                            XM)**2))
                              W = XM + X
                              W = W*4D0/(4D0 + TDIFFM) + (1D0 - W)**2
                              X = X*(1D0 - (1D0 + M + C41*(THR2 -
     1                            0.5D0))/(1D0 + C3*M)*X*(C1*W +
     1                            C2*X*DSQRT(W))) + XM
                              D2T2 = D2T/2D0
                              
                              IF (X.GE.1D0) THEN
                                      N = 1
                                      GO TO 3
                              END IF

!-----------------------------(NO FINITE SOLUTION WITH X > XM)

                      END IF

              END IF
              
!-------------(NOW HAVE A STARTER, SO PRECEED BY HALLEY)

    5         CONTINUE
              
              DO 4 I=1,3

              CALL TLAMB (M, Q, QSQFM1, X, 2, T, DT, D2T, D3T)
              
              T = TIN - T
              IF (DT.NE.0D0) X = X + T*DT/(DT*DT + T*D2T/2D0)
    4         CONTINUE

              IF (N.NE.3) RETURN

!-------------(EXIT IF ONLY ONE SOLUTION, NORMALLY WHEN M=0)
              
              N = 2
              XPL = X

!-------------(SECOND MULTI-REV STARTER)

    3         CALL TLAMB (M, Q, QSQFM1, 0D0, 0, T0, DT, D2T, D3T)

              TDIFF0 = T0 - TMIN
              TDIFF = TIN - T0
              IF (TDIFF.LE.0) THEN
                      X = XM - DSQRT(TDIFFM/(D2T2 - TDIFFM*(D2T2/TDIFF0
     1                    - 1D0/XM**2)))
              ELSE
                      X = -TDIFF/(TDIFF + 4D0)
                      !IJ = 200 IN ORIGINAL PAPER BUT NOT 1990 PAPER
                      W = X + C0*DSQRT(2D0*(1D0 - THR2))

                      IF (W.LT.0D0) X = X - DSQRT(D8RT(-W))*(X +
     1                                  DSQRT(TDIFF/(TDIFF + 1.5D0*T0)))
                      W = 4D0/(4D0 + TDIFF)
                      X = X*(1D0 + (1D0 + M + C42*(THR2 - 0.5D0))/(1D0 +
     1                    C3*M)*X*(C1*W - C2*X*DSQRT(W)))
                      IF (X.LE.-1D0) THEN
                              N = N - 1
!-----------------------------(NO FINITE SOLUTION WITH X < XM)
                              IF (N.EQ.1) X = XPL
                      END IF
              END IF

              GO TO 5

              END SUBROUTINE XLAMB

      SUBROUTINE VALAMB (GM, R1, R2, TH, TDELT, N, VR11, VT11, VR12,
     1           VT12, VR21, VT21, VR22, VT22)
              
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              PARAMETER (PI=4.0D0*DATAN(1.0D0), TWOPI=2D0*PI)
              
              M = TH/TWOPI
              THR2 = TH/2D0 - M*PI

              DR = R1 - R2
              R1R2 = R1*R2
              R1R2TH = 4D0*R1R2*DSIN(THR2)**2
              CSQ = DR**2 + R1R2TH
              C = DSQRT(CSQ)
              S = (R1 + R2 + C)/2D0
              GMS = DSQRT(GM*S/2D0)
              QSQFM1 = C/S
              Q = DSQRT(R1R2)*DCOS(THR2)/S
              
              IF (C.NE.0D0) THEN
                      RHO = DR/C
                      SIG = R1R2TH/CSQ
              ELSE
                      RHO = 0D0
                      SIG = 1D0
              END IF
              
              T = 4D0*GMS*TDELT/S**2
              
              CALL XLAMB (M, Q, QSQFM1, T, N, X1, X2)

!-------------PROCEED FOR SINGLE SOLUTION, OR A PAIR

              DO 1 I=1,N
             
              IF (I.EQ.1) THEN
                      X = X1
              ELSE
                      X = X2
              END IF

              CALL TLAMB (M,Q,QSQFM1,X,-1,UNUSED,QZMINX,QZPLX,ZPLQX)

              VT2 = GMS*ZPLQX*DSQRT(SIG)
              VR1 = GMS*(QZMINX - QZPLX*RHO)/R1
              VT1 = VT2/R1
              VR2 = -GMS*(QZMINX + QZPLX*RHO)/R2
              VT2 = VT2/R2

              IF (I.EQ.1) THEN
                      VR11 = VR1
                      VT11 = VT1
                      VR12 = VR2
                      VT12 = VT2
              ELSE
                      VR21 = VR1
                      VT21 = VT1
                      VR22 = VR2
                      VT22 = VT2
              END IF

    1         CONTINUE

              RETURN
              END SUBROUTINE VALAMB

      SUBROUTINE CALCPS (NHREV, RO1, RO3, IND, NUM, CX, CY, CZ)
!-------------(CALCULATE POSITION ALONG SIGHT-LINE TO COMPARE WITH
!-------------EXTERNALLY GIVEN LINE OF SIGHT)
              
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              COMMON /GIVEN/ OX1, OY1, OZ1, OX2, OY2, OZ2, OX3, OY3,
     1         OZ3, EL1, EM1, EN1, EL2, EM2, EN2, EL3, EM3, EN3, T12,
     1         T13, GM
              COMMON /USEFUL/ AL, Q, EI, BOM, R1, R2, R3
              COMMON /CNTROL/ LOBDEP, HN, MAXIT, PDINC, CRIVAL, CRITGM,
     1          NFAIL, ITNUM, CR, CRIT, AXRTIO
              COMMON /KNWNSL/ IKN, RO1KN(9), RO3KN(9), R1KNSQ(9),
     1          R3KNSQ(9)

              PARAMETER (PI=4.0D0*DATAN(1.0D0))
              
              RNORM(X,Y,Z) = DSQRT(X*X + Y*Y + Z*Z)
              X1 = OX1 + RO1*EL1
              Y1 = OY1 + RO1*EM1
              Z1 = OZ1 + RO1*EN1
              R1 = RNORM(X1, Y1, Z1)
              X3 = OX3 + RO3*EL3
              Y3 = OY3 + RO3*EM3
              Z3 = OZ3 + RO3*EN3
              R3 = RNORM(X3, Y3, Z3)
              
              CALL VECMUL (X1, Y1, Z1, X3, Y3, Z3, X13, Y13, Z13)
              
              TH = DATAN2(RNORM(X13, Y13, Z13), X1*X3 + Y1*Y3 + Z1*Z3)

!-------------(FAILS ONLY IF EITHER R1 OR R2 IS ZERO

              M = MOD(NHREV,2)
              IF (M.EQ.1) TH = PI - TH
              TH = TH + NHREV*PI

              CALL VALAMB (GM, R1, R3, TH, T13, NUM, VR1, VT1, VR3, VT3,
     1                     WR1, WT1, WR3, WT3)

              IF (NUM.GT.0) THEN

                      IF (IND.EQ.1) THEN
                              VR1 = WR1
                              VT1 = WT1
                              ALV = ALW
                      END IF

                      CALL VECMUL (X1, Y1, Z1, X3, Y3, Z3, XN, YN, ZN)

                      CALL VECMUL (XN, YN, ZN, X1, Y1, Z1, XT, YT, ZT)

                      RT = RNORM(XT, YT, ZT)
                      IF (M.EQ.1) RT = -RT
                      IF (RT.NE.0D0) RT = 1D0/RT
                      
                      V1X1 = VR1*X1/R1 + VT1*XT*RT
                      V1Y = VR1*Y1/R1 + VT1*YT*RT
                      V1Z = VR1*Z1/R1 + VT1*ZT*RT

                      CALL PV3ELS (GM, X1, Y1, Z1, V1X1, V1Y, V1Z, AL,
     1                             Q, EI, BOM, OM, TAU)
                      CALL ELS3PV (GM, AL, Q, EI, BOM, OM, TAU+T12, X2,
     1                             Y2, Z2, W1, W2, W3)
!----------(ALV GIVES BETTER ACCURACY THAN AL FOR NEAR-PARABOLIC ORBITS)
                      
                      R2 = RNORM(X2, Y2, Z2)
                      
!------------(ONLY FOR THE CONVERGENCE CRITERION IN THE CALLING ROUTINE)

                      CX = X2 - OX2
                      CY = Y2 - OY2
                      CZ = Z2 - OZ2
              END IF
              RETURN
              END SUBROUTINE CALCPS

      SUBROUTINE OBS3LS (NHREV, IND, RO1, RO3)
!-----ORBIT FROM THREE OBSERVED LINES OF SIGHT (ANGLES ONLY)
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              LOGICAL LOBDEP
              DIMENSION NLDF(5), W(9), W1(9), W3(9), UW(9), RO1SQ(9),
     1         RO3SQ(9), RO1KNO(9), RO3KNO(9), RO1QU(9), RO3QU(9),
     1         RO13SQ(9), RO13QU(9)
              COMMON /GIVEN/ OX1, OY1, OZ1, OX2, OY2, OZ2, OX3, OY3,
     1         OZ3, EL1, EM1, EN1, EL2, EM2, EN2, EL3, EM3, EN3, T12,
     1         T13, GM
              COMMON /USEFUL/ AL, Q, EI, BOM, R1, R2, R3
              COMMON /CNTROL/ LOBDEP, HN, MAXIT, PDINC, CRIVAL, CRITGM,
     1         NFAIL, ITNUM, CR, CRIT, AXRTIO
              COMMON /KNWNSL/ IKN, RO1KN(9), RO3KN(9), R1KNSQ(9),
     1         R3KNSQ(9)

              PARAMETER (PI=4.0D0*DATAN(1.0D0), RADIAN=180D0/PI)

              RNORM(X,Y,Z) = DSQRT(X*X + Y*Y + Z*Z)

              AXRTIO = 0D0
              FCOLD = 0D0
              NFAIL = 0

    1         PDINCM = PDINC
              IF (LOBDEP .AND. IKN.EQ.0 .AND. RO1.EQ.RO3) THEN


                      CALL CALCPS(NHREV,RO1,RO3,IND,NUM,CX,CY,CZ)

                      IF (NUM.GT.0) THEN
                              CR = EL2*CX + EM2*CY + EN2*CZ

                              IF (CR.LT.0D0) THEN
                                      DR = PDINC*(R1 + R3)

                                      CALL CALCPS(NHREV,RO1+DR,RO3+DR,
     1                                            IND,NUM,CX,CY,CZ)

                                      RO1 = RO1 + DR*CR/(CR - EL2*CX -
     1                                      EM2*CY - EN2*CZ)

                                      IF (RO1.LT.0D0) RO1 = 0D0
                                      RO3 = RO1

                              END IF
                      END IF
              END IF

              DO 8 ITNUM=1,MAXIT
              PRINT *, "ITERATION # = ", ITNUM
              NPDF = 0
    2         NMODS = 0
              NFLMOD = 0

    3         CALL CALCPS(NHREV,RO1,RO3,IND,NUM,CX,CY,CZ)

              R10 = R1
              R30 = R3
              IF (NUM.EQ.0) THEN

                      IF (ITNUM.GT.1 .AND. NFLMOD.LT.3) THEN

                              NFLMOD = NFLMOD + 1
                              D3 = D3/3D0
                              D1 = D1/3D0
                              RO1 = RO1OLD + D3
                              RO3 = RO3OLD + D1
                              GO TO 3
                      ELSE
                              NFAIL = NFAIL + 1
                              IF (NFAIL.EQ.1) THEN

                                      R01 = -(EL1*OX1+EM1*OY1+EN1*OZ1)
                                      RO3 = -(EL3*OX3+EM3*OY3+EN3*OZ3)

                                      IF (LOBDEP .AND. RO1.LT.0D0) 
     1                                    RO1 = 0D0
                                      IF (LOBDEP .AND. RO3.LT.0D0) 
     1                                    RO3 = 0D0

                                      GO TO 1

                              ELSE IF (NFAIL.EQ.2) THEN

                                      X13 = OX3 - OX1
                                      Y13 = OY3 - OY1
                                      Z13 = OZ3 - OZ1
                                      D1 = X13*EL1+Y13*EM1+Z13*EN1
                                      D3 = X13*EL3+Y13*EM3+Z13*EN3
                                      D2 = EL1*EL3+EM1*EM3+EN1*EN3
                                      D4 = 1D0 - D2**2
                                      RO1 = (D1 - D3*D2)/D4
                                      RO3 = (D1*D2 - D3)/D4

                                      IF (LOBDEP .AND. RO1.LT.0D0)
     1                                    RO1 = 0D0
                                      IF (LOBDEP .AND. RO3.LT.0D0)
     1                                    RO3 = 0D0

                                      GO TO 1

                              ELSE IF (NFAIL.EQ.3) THEN

                                      RO1 = 0
                                      RO3 = 0
                                      GO TO 1

                              ELSE
                                      RETURN
                              END IF
                      END IF
              END IF

              CR = EL2*CX + EM2*CY + EN2*CZ

              CALL VECMUL(EL2,EM2,EN2,CX,CY,CZ,ENX,ENY,ENZ)

              CALL VECMUL(ENX,ENY,ENZ,EL2,EM2,EN2,PX,PY,PZ)

              PR = RNORM(PX, PY, PZ)

              CALL VECMUL(EL2,EM2,EN2,PX,PY,PZ,ENX,ENY,ENZ)

              ENR = RNORM(ENX, ENY, ENZ)

              IF (ENR.EQ.0D0) THEN
                      CRIT = 0D0
                      GO TO 7

              END IF

              F = (PX*CX + PY*CY + PZ*CZ)/PR
              FC = F

              DO 4 I=1,IKN
              RO1KNO(I) = RO1 - RO1KN(I)
              RO3KNO(I) = RO3 - RO3KN(I)
              RO1SQ(I) = RO1KNO(I)**2
              RO3SQ(I) = RO3KNO(I)**2
              RO1QU(I) = RO1SQ(I) + R1KNSQ(I)
              RO3QU(I) = RO3SQ(I) + R3KNSQ(I)
              RO13SQ(I) = RO1SQ(I) + RO3SQ(I)
              RO13QU(I) = RO1QU(I) + RO3QU(I)
              FC = FC*DSQRT(RO13QU(I)/RO13SQ(I))
    4         CONTINUE

              IF (ITNUM.GT.1 .AND. FC.GT.2D0*FCOLD .AND. 
     1             NMODS.LT.2) THEN
                      FSUM = FC + FCOLD
                      RO1 = (FC*RO1OLD + FCOLD*RO1)/FSUM
                      RO3 = (FC*RO3OLD + FCOLD*RO3)/FSUM
                      NMODS = NMODS + 1

                      GO TO 3

              END IF

              DRO1 = PDINCM*R10
              DRO3 = PDINCM*R30
              D2RO1 = 2D0*DRO1
              D2RO3 = 2D0*DRO3
              DRO1SQ = DRO1**2
              DRO3SQ = DRO3**2
              
              CALL CALCPS(NHREV,RO1-DRO1,RO3,IND,NLDF(1),CX,CY,CZ)

              FM1 = (PX*CX + PY*CY + PZ*CZ)/PR - F
              GM1 = (ENX*CX + ENY*CY + ENZ*CZ)/ENR

              CALL CALCPS(NHREV,RO1+DRO1,RO3,IND,NLDF(2),CX,CY,CZ)

              FP1 = (PX*CX + PY*CY + PZ*CZ)/PR - F
              GP1 = (ENX*CX + ENY*CY + ENZ*CZ)/ENR
              FD1 = (FP1 - FM1)/D2RO1
              FDD1 = (FP1 + FM1)/DRO1SQ
              GD1 = (GP1 - GM1)/D2RO1
              GDD1 = (GP1 + GM1)/DRO1SQ

              CALL CALCPS(NHREV,RO1,RO3-DRO3,IND,NLDF(3),CX,CY,CZ)

              FM3 = (PX*CX + PY*CY + PZ*CZ)/PR - F
              GM3 = (ENX*CX + ENY*CY + ENZ*CZ)/ENR

              CALL CALCPS(NHREV,RO1,RO3+DRO3,IND,NLDF(4),CX,CY,CZ)

              FP3 = (PX*CX + PY*CY + PZ*CZ)/PR - F
              GP3 = (ENX*CX + ENY*CY + ENZ*CZ)/ENR
              FD3 = (FP3 - FM3)/D2RO3
              FDD3 = (FP3 + FM3)/DRO3SQ
              GD3 = (GP3 - GM3)/D2RO3
              GDD3 = (GP3 + GM3)/DRO3SQ

              CALL CALCPS(NHREV,RO1+DRO1,RO3+DRO3,IND,NLDF(5),CX,CY,CZ)

              F13 = (PX*CX + PY*CY + PZ*CZ)/PR - F
              G13 = (ENX*CX + ENY*CY + ENZ*CZ)/ENR
              ROFAC = DRO1/DRO3
              FD13 = F13/(DRO1*DRO3) - (FD1/DRO3 + FD3/DRO1) -
     1               0.5D0*(FDD1*ROFAC + FDD3/ROFAC)
              GD13 = G13/(DRO1*DRO3) - (GD1/DRO3 + GD3/DRO1) -
     1               0.5D0*(GDD1*ROFAC + GDD3/ROFAC)

              DO 5 I=1,5
              IF (NLDF(I).EQ.0) THEN
                      NPDF = NPDF + 1

                      IF (NPDF.LE.3) THEN
                              PDINCM = PDINCM/10D0
                              GO TO 2

                      ELSE

                              NFAIL = 1

                              RETURN
                      END IF
              END IF

    5         CONTINUE

              DO 6 I=1,IKN
              W(I) = 1D0/RO13SQ(I) - 1D0/RO13QU(I)
              UW(I) = W(I) - (2D0/RO13SQ(I) + 2D0/RO13QU(I))
              W1(I) = W(I)*RO1KNO(I)
              W3(I) = W(I)*RO3KNO(I)
              FD1 = FD1 - F*W1(I)
              FD3 = FD3 - F*W3(I)
              FDD1 = FDD1 - (2D0*FD1*W1(I) + W(I)*F*(1D0 +
     1               RO1SQ(I)*UW(I)))
              GDD1 = GDD1 - 2D0*W1(I)*GD1
              FDD3 = FDD3 - (2D0*FD3*W3(I) + W(I)*F*(1D0 +
     1               RO3SQ(I)*UW(I)))
              GDD3 = GDD3 - 2D0*W3(I)*GD3
              FD13 = FD13 - (FD3*W1(I) + FD1*W3(I) + W(I)*F*
     1               RO1KNO(I)*RO3KNO(I)*UW(I))
              GD13 = GD13 - (W1(I)*GD3 + W3(I)*GD1)
    6         CONTINUE
              
              DEL = FD1*GD3 - FD3*GD1
              D3NR = -GD3*F/DEL
              D1NR = GD1*F/DEL
              FD1H = FD1 + HN*(FDD1*D3NR + FD13*D1NR)
              FD3H = FD3 + HN*(FD13*D3NR + FDD3*D1NR)
              GD1H = GD1 + HN*(GDD1*D3NR + GD13*D1NR)
              GD3H = GD3 + HN*(GD13*D3NR + GDD3*D1NR)
              DELH = FD1H*GD3H - FD3H*GD1H
              D3 = -GD3H*F/DELH
              D1 = GD1H*F/DELH
              FGXY = FD1**2 + GD1**2 + FD3**2 + GD3**2
              AXRTIO = 2D0*DABS(DEL)/(FGXY+DSQRT(FGXY**2 - 4D0*DEL**2))

              IF (DABS(AXRTIO).LT.CRITGM) THEN
                      D3 = DSIGN(DSQRT(DABS(D3NR*D3)), D3NR)
                      D1 = DSIGN(DSQRT(DABS(D1NR*D1)), D1NR)

              END IF

              RO1OLD = RO1
              RO3OLD = RO3
              RO1 = RO1 + D3
              RO3 = RO3 + D1
              DEN = MAX1(CR,R2)
              IF (AL.GT.0D0) DEN = MAX1(DEN, T12*DSQRT(AL*(2D0*GM/AL -
     1                             R2)/R2))
              CRIT = F/DEN
              CRITSQ = CRIT**2
    7         IF (CRITSQ.LT.CRIVAL) GO TO 9
              FCOLD = FC
    8         CONTINUE

              NFAIL = -1
              RETURN
    9         NFAIL = 0
              
              CALL CALCPS(NHREV,RO1,RO3,IND,NUM,CX,CY,CZ)
              
              CR = EL2*CX + EM2*CY + EN2*CZ

              RETURN
              END SUBROUTINE OBS3LS


      SUBROUTINE PV3ELS(GM,X,Y,Z,XDOT,YDOT,ZDOT,AL,Q,EI,BOM,
     1                       OM,TAU)
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              PARAMETER (PI=4.0D0*DATAN(1.0D0))

              XSQYSQ = X*X + Y*Y
              RSQ = XSQYSQ + Z*Z
              R = DSQRT(RSQ)
              VR = (X*XDOT + Y*YDOT + Z*ZDOT)/R
              HX = Y*ZDOT - Z*YDOT
              HY = Z*XDOT - X*ZDOT
              HZ = X*YDOT - Y*XDOT
              HSQ = HX*HX + HY*HY + HZ*HZ

              IF (HSQ.EQ.0D0) THEN
                      !(RECTILINEAR ORBIT)
                      EI = PI/2D0
                      IF (XSQYSQ.EQ.0D0) THEN
                              !(AXIAL ORBIT)
                              BOM = 0D0
                      ELSE
                              !(GENERAL RECTILINEAR ORBIT)
                              BOM = DATAN2(Y, X)
                      END IF

                      U = DATAN2(Z, DSQRT(XSQYSQ))
                      VT = 0D0
              ELSE
                      !(NON-DEGENERATE ORBIT)
                      BX = HY*Z - HZ*Y
                      BY = HZ*X - HX*Z
                      BZ = HX*Y - HY*X
                      HX = Y*BZ - Z*BY
                      HY = Z*BX - X*BZ
                      HZ = X*BY - Y*BX
                      W = HX*HX + HY*HY
                      H = DSQRT(W + HZ*HZ)
                      EI = DATAN2(DSQRT(W), HZ)
                      IF (W.EQ.0D0) THEN
                              !(ORBIT IN REFERENCE PLANE)
                              BOM = 0D0
                              U = DATAN2(Y*DSIGN(1D0, HZ), X)
                      ELSE
                              !(GENERAL ORBIT)
                              BOM = DATAN2(HX, -HY)
                              U = DATAN2(H*Z, RSQ*BZ)
                      END IF

                      VT = H/(R*RSQ)
              END IF

              CALL PV2ELS(GM,R,U,VR,VT,AL,Q,OM,TAU)

              RETURN

              END SUBROUTINE PV3ELS

      SUBROUTINE PV2ELS(GM,R,U,VR,VT,AL,Q,OM,TAU)
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              LOGICAL,PARAMETER :: L = .FALSE.

              PARAMETER (PI=4.0D0*DATAN(1.0D0), SW=0.25D0)

              !(ALL ORBITS)
              VSQ = VR*VR + VT*VT

              AL = 2D0*GM/R - VSQ
              ALP = DABS(AL)
              RTAL = DSQRT(ALP)
              D = R*VR
              H = R*VT
              P = H*H
              ESQ1 = P*AL
              ES = D*RTAL
              ESES = ES*ES
              EC = R*VSQ - GM
              ECEC = EC*EC
              IF (AL.GT.0D0) THEN
                      !(ONE ESQ FORMULA SUPERIOR FOR THE ELLIPSE)
                      ESQ = ECEC + ESES
              ELSE
                      !(DIFFERENT FORMULA SUPERIOR FOR THE HYPERBOLA)
                      ESQ = GM*GM - ESQ1
              END IF

              E = DSQRT(ESQ)
              Q = P/(GM + E)

              IF (AL.EQ.0D0) THEN
                      !(PARABOLA)
                      TAU = D*(2D0*Q + R)/(3D0*GM)
                      V = 2D0*DATAN2(VR, VT)
              ELSE IF (E.EQ.0D0) THEN
                      !(CIRCLE)
                      TAU = 0D0
                      V = 0D0
              ELSE
                      !(ELLIPSE OR HYPERBOLA)
                      E1 = AL*Q
                      IF (AL.GT.0D0) THEN
                              !(ELLIPSE)
                              EH = DATAN2(ES, EC)
                              IF (GM*EH*EH/6D0+E1 .GE. GM*SW) THEN
                                      !(GENERAL CASE)
                                      EM = GM*EH - ES
                                      ECESQ = GM*EC - ESQ
                              ELSE
                                      !(FOR E1 AND EH BOTH NEAR 0)
                                      EM = GM*EMKEP(E1/GM, EH)
                                      ECESQ = (ESQ1*ECEC -
     1                                        ESQ*ESES)/(ESQ + GM*EC)
                              END IF
                      ELSE
                              !(HYRPERBOLA)
                              EH = ASINH(ES/E)
                              IF (GM*EH*EH/6D0-E1 .GE. GM*SW) THEN
                                      !(GENERAL CASE)
                                      EM = ES - GM*EH
                                      ECESQ = ESQ - GM*EC
                              ELSE
                                      !(FOR E1 AND EH BOTH NEAR 0)
                                      EM = E*SHMKEP(-E1/E, ES/E)
                                      ECESQ = -(ESQ1*ECEC +
     1                                        ESQ*ESES)/(ESQ + GM*EC)
                              END IF
                      END IF
                      !(ELLIPSE OR HYPERBOLA STILL)
                      EN = ALP*RTAL
                      TAU = EM/EN
                      V = DATAN2(ES*H*RTAL, ECESQ)
              END IF

              !(ALL ORBITS)
              OM = U - V

              !(THE FOLLOWING IS ONLY FOR TESTING IF L=TRUE)

              IF (L .AND. AL.GT.0D0) THEN
                      !(FOR ELLIPSE, ADJUST REVOLUTIONS IF REQUIRED
                      !USING L)
                      ADJ = 2D0*PI*DSIGN(DINT(DABS(OM/(2D0*PI)) + 
     1                      0.5D0),OM)
                      OM = OM - ADJ
                      TAU = TAU + ADJ/EN
              END IF

              RETURN
              
              END SUBROUTINE PV2ELS

      SUBROUTINE ELS3PV(GM,AL,Q,EI,BOM,OM,TAU,X,Y,Z,XDOT,YDOT,ZDOT)
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              PARAMETER (PI=4.0D0*DATAN(1.0D0))

              CALL ELS2PV(GM,AL,Q,OM,TAU,R,U,VR,VT)

              C = DCOS(U)
              S = DSIN(U)
              X1 = R*C
              Y1 = R*S
              X2 = VR*C - VT*S
              Y2 = VR*S + VT*C
              C = DCOS(EI)
              S = DSIN(EI)
              Z = Y1*S
              Y1 = Y1*C
              ZDOT = Y2*S
              Y2 = Y2*C
              C = DCOS(BOM)
              S = DSIN(BOM)
              X = X1*C - Y1*S
              Y = X1*S + Y1*C
              XDOT = X2*C - Y2*S
              YDOT = X2*S + Y2*C

              RETURN

              END SUBROUTINE ELS3PV

      SUBROUTINE ELS2PV(GM,AL,Q,OM,TAU,R,U,VR,VT)

              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              PARAMETER (PI=4.0D0*DATAN(1.0D0))

              IF (AL.EQ.0D0) THEN
                      !(PARABOLA - GM CANNOT BE 0)

                      D = DCBSOL(0.5D0/GM, Q, 1.5D0*GM*TAU)

                      R = Q + 0.5D0*D*D/GM
                      H = DSQRT(2D0*GM*Q)
                      V = 2D0*DATAN2(D, H)
              ELSE
                      !(ELLIPSE OR HYPERBOLA)
                      E1 = AL*Q
                      E = GM - E1
                      EP1 = GM + E
                      H = DSQRT(Q*EP1)
                      ALP = DABS(AL)
                      RTAL = DSQRT(ALP)
                      !(LAST 6 ITEMS COULD BE SAVED IF REPEATING
                      !AM,ALV,Q)
                      EM = TAU*ALP*RTAL
                      IF (AL.GT.0D0) THEN
                              !(ELLIPSE - GM CANNOT BE 0)
                              !(MAKE SURE E1 ARGUMENT TO EKEPL IS
                              !BETWEEN [0,1])
                              EE2 = 0.5D0*EKEPL(EM/GM, E1/GM)

                              S2 = DSIN(EE2)
                              C2 = DCOS(EE2)
                              R = Q + 2D0*E*S2*S2/AL
                              D = 2D0*E*S2*C2/RTAL
                              V = 2D0*DATAN2(EP1*S2, H*RTAL*C2)
                              EMV = EM/GM - V
                              V = V + 4D0*PI*DSIGN(DINT(DABS(EMV/
     1                            (4D0*PI)) + 0.5D0), EMV)
                      ELSE
                              !(HYPERBOLA)
                              S = SHKEPL(EM/E, -E1/E)

                              S2 = S*S
                              C = DSQRT(1D0 + S2)
                              S2 = S2/(C + 1D0)
                              R = Q - E*S2/AL
                              D = E*S/RTAL
                              V = DATAN2(S*H*RTAL, -GM*S2 - E1)
                      END IF
              END IF

              !(ALL ORBITS)
              U = OM + V
              VR = D/R
              VT = H/R
              
              RETURN
              
              END SUBROUTINE ELS2PV

      FUNCTION SHKEPL(EL,G1)
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              PARAMETER (PI=4.0D0*DATAN(1.0D0), SW=0.5D0)

              S = EL

              IF (EL.EQ.0D0) GO TO 2

              !START BASED ON LAGRANGE'S THEOREM
              G = 1D0 - G1
              CL = DSQRT(1D0 + EL**2)
              AL = ASINH(EL)
              W = G**2*AL/CL**3
              S = 1D0 - G/CL
              S = EL + G*AL/DCUBRT(S**3 + W*EL*(1.5D0 -
     1            G/0.75D0))

              !TWO ITERATIONS (AT MOST) OF HALLEY-THEN-NEWTON
              DO 1 ITER=1,2
              S0 = S*S
              S1 = S0 + 1D0
              S2 = DSQRT(S1)
              S3 = S1*S2
              FDD = G*S/S3
              FDDD = G*(1D0 - 2D0*S0)/(S1*S3)
              IF ((1D0/6D0)*S0+G1 .GE. SW) THEN
                      F = (S - G*ASINH(S)) - EL
                      FD = 1D0 - G/S2
              ELSE
                      F = SHMKEP(G1, S) - EL

                      FD = (S0/(S2 + 1D0) + G1)/S2
              END IF

              DS = F*FD/(0.5D0*F*FDD - FD*FD)
              STEMP = S + DS

              IF (STEMP.EQ.S) GO TO 2
              F = F + DS*(FD + 0.5D0*DS*(FDD + 
     1            (1D0/3D0)*DS*FDDD))
              FD = FD + DS*(FDD + 0.5D0*DS*FDDD)
              S = STEMP - F/FD
    1         CONTINUE

    2         SHKEPL = S

              RETURN

              END FUNCTION SHKEPL

      FUNCTION DCUBRT(X) RESULT(C)
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)

              IF (X.EQ.0) THEN
                      C = 0D0
              ELSE
                      Y = DABS(X)
                      C = Y**(1D0/3D0)
                      C = C - (1D0/3D0)*(C - Y/C**2)
                      C = DSIGN(C, X)
              END IF

              RETURN

              END FUNCTION DCUBRT

      FUNCTION DCBSOL(A,B,C) RESULT(X)
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)

              IF (A.EQ.0 .AND. B.EQ.0 .OR. C.EQ.0) THEN
                      X = 0D0
              ELSE
                      BSQ = B*B
                      D = DSQRT(A)*DABS(C)
                      D = DCUBRT(D+DSQRT(B*BSQ+D*D))**2

                      X = 2D0*C/(D + B + BSQ/D)
              END IF

              RETURN

              END FUNCTION DCBSOL

      FUNCTION EKEPL(EM,E1)
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)

              INTEGER :: ITER
              PARAMETER (PI=4.0D0*ATAN(1.0D0), SW=0.25D0)

              !RANGE-REDUCE EM TO LIE IN REANGE -PI TO PI
              EMR = MOD(EM, 2D0*PI)
              IF (EMR.LT.-PI) EMR = EMR + 2D0*PI
              IF (EMR.GT.PI) EMR = EMR - 2D0*PI

              EE = EMR

              IF (EE.NE.0D0) THEN
                      IF (EE.LE.0D0) EE = -EE
                      !(EMR IS RANGE-REDUCED EM AND EE IS ABSOLUTE VALUE
                      !OF EMR)
                      !STARTER BY FIRST SOLVING CUBIC EQUATION

                      E = 1D0 - E1
                      W = DCBSOL(E,2D0*E1,3D0*EE)

                      !EFFECTIVELY INTERPOLATE IN EMR (ABS)
                      EE = (EE*EE + (PI - EE)*W)/PI
                      IF (EMR.LT.0D0) EE = -EE

                      !DO TWO ITERATIONS OF HALLEY EACH FOLLOWED BY
                      !NEWTON
                      DO ITER=1,2
                      FDD = E*DSIN(EE)
                      FDDD = E*DCOS(EE)
                      IF (EE*EE/6D0+E1 .GE. SW) THEN
                              F = (EE - FDD) - EMR
                              FD = 1D0 - FDDD
                      ELSE
                              F = EMKEP(E1,EE) - EMR

                              FD = 2D0*E*DSIN(0.5D0*EE)**2 + E1
                      END IF

                      DEE = F*FD/(0.5D0*F*FDD - FD*FD)
                      F = F + DEE*(FD + 0.5D0*DEE*(FDD + 
     1                    (1D0/3D0)*DEE*FDDD))
                      !TO REDUCE THE DANGER OF THE UNDERFLOW REPLACE THE
                      !LAST LINE BY;
                      !W = FD + 0.5*DEE*(FDD + (1/3)*DEE*FDDD)
                      FD = FD + DEE*(FDD + 0.5D0*DEE*FDDD)
                      EE = EE + DEE - F/FD
                      !IF REPLACING AS ABOVE, THEN ALSO REPLACE THE LAST
                      !LINE BY;
                      !EE = EE - (F - DEE*(FD - W))/FD
                      END DO

              END IF

              !RANGE-EXPAND
              EKEPL = EE + (EM - EMR)

              RETURN

              END FUNCTION EKEPL

      FUNCTION EMKEP(E1,EE)
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              PARAMETER (PI=4.0D0*DATAN(1.0D0))

              X = E1*DSIN(EE)
              EE2 = -EE*EE
              TERM = EE
              D = 0D0

              DO
    1         D = D + 2D0
              TERM = TERM*EE2/(D*(D + 1D0))
              X0 = X
              X = X - TERM
              IF(X.NE.X0) GO TO 1
              IF(X.EQ.X0) EXIT
              END DO

              EMKEP = X

              RETURN

              END FUNCTION EMKEP

      FUNCTION SHMKEP(G1, S)
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              PARAMETER (PI=4.0D0*DATAN(1.0D0))

              G = 1D0 - G1
              T = S/(1D0 + DSQRT(1D0 + S*S))
              TSQ = T*T
              X = S*(G1 + G*TSQ)
              TERM = 2D0*G*T
              TWOI1 = 1D0
 
    1         TWOI1 = TWOI1 + 2D0
              TERM = TERM*TSQ
              X0 = X
              X = X - TERM/TWOI1
              IF(X.NE.X0) GO TO 1

              SHMKEP = X

              RETURN

              END FUNCTION SHMKEP

      SUBROUTINE VECMUL(X1,Y1,Z1,X3,Y3,Z3,X13,Y13,Z13)
              !CROSS-PRODUCT
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              INTENT(IN) :: X1,Y1,Z1
              INTENT(IN) :: X3,Y3,Z3
              INTENT(OUT) :: X13,Y13,Z13

              X13 = Y1*Z3 - Y3*Z1
              Y13 = Z1*X3 - Z3*X1
              Z13 = X1*Y3 - X3*Y1

              RETURN

              END SUBROUTINE VECMUL
