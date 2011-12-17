
!   PURPOSE             - COMPUTE THE FAST FOURIER TRANSFORM OF A
!                           COMPLEX VALUED SEQUENCE
!
!   USAGE               - CALL FFTCC (A,N,IWK,WK)
!
!   ARGUMENTS    A      - COMPLEX VECTOR OF LENGTH N. ON INPUT A
!                           CONTAINS THE COMPLEX VALUED SEQUENCE TO BE
!                           TRANSFORMED. ON OUTPUT A IS REPLACED BY THE
!                           FOURIER TRANSFORM.
!                N      - INPUT NUMBER OF DATA POINTS TO BE
!                           TRANSFORMED. N MAY BE ANY POSITIVE
!                           INTEGER.
!                IWK    - INTEGER WORK VECTOR OF LENGTH 6*N+150.
!                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS)
!                WK     - REAL WORK VECTOR OF LENGTH 6*N+150.
!                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS)
!
!   REMARKS  1.  FFTCC COMPUTES THE FOURIER TRANSFORM, X, ACCORDING
!                TO THE FOLLOWING FORMULA;
!
!                  X(K+1) = SUM FROM J = 0 TO N-1 OF
!                           A(J+1)*CEXP((0.0,(2.0*PI*J*K)/N))
!                  FOR K=0,1,...,N-1 AND PI=3.1415...
!
!                NOTE THAT X OVERWRITES A ON OUTPUT.
!            2.  FFTCC CAN BE USED TO COMPUTE
!
!                  X(K+1) = (1/N)*SUM FROM J = 0 TO N-1 OF
!                           A(J+1)*CEXP((0.0,(-2.0*PI*J*K)/N))
!                  FOR K=0,1,...,N-1 AND PI=3.1415...
!
!                BY PERFORMING THE FOLLOWING STEPS;
!
!                     DO 10 I=1,N
!                        A(I) = CONJG(A(I))
!                  10 CONTINUE
!                     CALL FFTCC (A,N,IWK,WK)
!                     DO 20 I=1,N
!                        A(I) = CONJG(A(I))/N
!                  20 CONTINUE
!

!-----------------------------------------------------------------------
!
      SUBROUTINE FFTCC (A,N,IWK,WK)
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IWK(6*N+150)
      REAL               WK(6*N+150)
      COMPLEX            A(N)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IAM,IAP,IBM,IBP,IC,ICC,ICF,ICK,ID,IDM1,II, &
                         IJA,IKB,IKT,ILL,IM,IRD,ISF,ISK,ISP,ISS,ITA,ITB, &
                         J,JA,JF,JJ,JK,K,K0,K1,K2,K3,KA,KB,KD2,KF,KH,KN, &
                         KT,KTP,L,L1,M,MM,MM1,MP
      REAL               CM,SM,C1,C2,C3,S1,S2,S3,C30,RAD,A0,A1,A4,B4, &
                         A2,A3,B0,B1,B2,B3,ZERO,HALF,ONE,TWO,Z0(2), &
                         Z1(2),Z2(2),Z3(2),Z4(2)
      COMPLEX            ZA0,ZA1,ZA2,ZA3,ZA4,AK2
      EQUIVALENCE        (ZA0,Z0(1)),(ZA1,Z1(1)),(ZA2,Z2(1)), &
                         (ZA3,Z3(1)),(A0,Z0(1)),(B0,Z0(2)),(A1,Z1(1)), &
                         (B1,Z1(2)),(A2,Z2(1)),(B2,Z2(2)),(A3,Z3(1)), &
                         (B3,Z3(2)),(ZA4,Z4(1)),(Z4(1),A4),(Z4(2),B4)
      DATA               RAD/6.283185/, &
                         C30/.8660254/
      DATA               ZERO,HALF,ONE,TWO/0.0,0.5,1.0,2.0/
!                                  FIRST EXECUTABLE STATEMENT
      IF (N .EQ. 1) GO TO 9005
      K = N
      M = 0
      J = 2
      JJ = 4
      JF = 0
!                                  DETERMINE THE SQUARE FACTORS OF N
      IWK(1) = 1
    5 I = K/JJ
      IF (I*JJ .NE. K) GO TO 10
      M = M+1
      IWK(M+1) = J
      K = I
      GO TO 5
   10 J = J + 2
      IF (J .EQ. 4) J = 3
      JJ = J * J
      IF (JJ .LE. K) GO TO 5
      KT = M
!                                  DETERMINE THE REMAINING FACTORS OF N
      J = 2
   15 I = K / J
      IF (I*J .NE. K) GO TO 20
      M = M + 1
      IWK(M+1) = J
      K = I
      GO TO 15
   20 J = J + 1
      IF (J .EQ. 3) GO TO 15
      J = J + 1
      IF(J.LE.K) GO TO 15
      K = IWK(M+1)
      IF (IWK(KT+1) .GT. IWK(M+1)) K = IWK(KT+1)
      IF(KT.LE.0) GO TO 30
      KTP = KT + 2
      DO 25  I = 1,KT
         J = KTP - I
         M = M+1
         IWK(M+1) = IWK(J)
   25 CONTINUE
   30 MP = M+1
      IC = MP+1
      ID = IC+MP
      ILL = ID+MP
      IRD = ILL+MP+1
      ICC = IRD+MP
      ISS = ICC+MP
      ICK = ISS+MP
      ISK = ICK+K
      ICF = ISK+K
      ISF = ICF+K
      IAP = ISF+K
      KD2 = (K-1) / 2 + 1
      IBP = IAP + KD2
      IAM = IBP + KD2
      IBM = IAM + KD2
      MM1 = M-1
      I = 1
   35 L = MP - I
      J = IC - I
      IWK(ILL+L) = 0
      IF ((IWK(J-1) + IWK(J)) .EQ. 4) IWK(ILL+L) = 1
      IF (IWK(ILL+L) .EQ. 0) GO TO 40
      I = I + 1
      L = L - 1
      IWK(ILL+L) = 0
   40 I = I + 1
      IF(I.LE.MM1) GO TO 35
      IWK(ILL+1) = 0
      IWK(ILL+MP) = 0
      IWK(IC) = 1
      IWK(ID) = N
      DO 45  J = 1,M
         K = IWK(J+1)
         IWK(IC+J) = IWK(IC+J-1) * K
         IWK(ID+J) = IWK(ID+J-1) / K
         WK(IRD+J) = RAD/IWK(IC+J)
         C1 = RAD/K
         IF (K .LE. 2) GO TO 45
         WK(ICC+J) = COS(C1)
         WK(ISS+J) = SIN(C1)
   45 CONTINUE
      MM = M
      IF (IWK(ILL+M) .EQ. 1) MM = M - 1
      IF (MM .LE. 1) GO TO 50
      SM = IWK(IC+MM-2) * WK(IRD+M)
      CM = COS(SM)
      SM = SIN(SM)
   50 KB = 0
      KN = N
      JJ = 0
      I = 1
      C1 = ONE
      S1 = ZERO
      L1 = 1
   55 IF (IWK(ILL+I+1) .EQ. 1) GO TO 60
      KF = IWK(I+1)
      GO TO 65
   60 KF = 4
      I = I+1
   65 ISP = IWK(ID+I)
      IF (L1 .EQ. 1) GO TO 70
      S1 = JJ * WK(IRD+I)
      C1 = COS(S1)
      S1 = SIN(S1)
!                                  FACTORS OF 2, 3, AND 4 ARE
!                                  HANDLED SEPARATELY.
   70 IF (KF .GT. 4) GO TO 140
      GO TO (75,75,90,115), KF
   75 K0 = KB + ISP
      K2 = K0 + ISP
      IF (L1 .EQ. 1) GO TO 85
   80 K0 = K0 - 1
      IF (K0 .LT. KB) GO TO 190
      K2 = K2 - 1
      ZA4 = A(K2+1)
      A0 = A4*C1-B4*S1
      B0 = A4*S1+B4*C1
      A(K2+1) = A(K0+1)-ZA0
      A(K0+1) = A(K0+1)+ZA0
      GO TO 80
   85 K0 = K0 - 1
      IF (K0 .LT. KB) GO TO 190
      K2 = K2 - 1
      AK2 = A(K2+1)
      A(K2+1) = A(K0+1)-AK2
      A(K0+1) = A(K0+1)+AK2
      GO TO 85
   90 IF (L1 .EQ. 1) GO TO 95
      C2 = C1 * C1 - S1 * S1
      S2 = TWO * C1 * S1
   95 JA = KB + ISP - 1
      KA = JA + KB
      IKB = KB+1
      IJA = JA+1
      DO 110 II = IKB,IJA
         K0 = KA - II + 1
         K1 = K0 + ISP
         K2 = K1 + ISP
         ZA0 = A(K0+1)
         IF (L1 .EQ. 1) GO TO 100
         ZA4 = A(K1+1)
         A1 = A4*C1-B4*S1
         B1 = A4*S1+B4*C1
         ZA4 = A(K2+1)
         A2 = A4*C2-B4*S2
         B2 = A4*S2+B4*C2
         GO TO 105
  100    ZA1 = A(K1+1)
         ZA2 = A(K2+1)
  105    A(K0+1) = CMPLX(A0+A1+A2,B0+B1+B2)
         A0 = -HALF * (A1+A2) + A0
         A1 = (A1-A2) * C30
         B0 = -HALF * (B1+B2) + B0
         B1 = (B1-B2) * C30
         A(K1+1) = CMPLX(A0-B1,B0+A1)
         A(K2+1) = CMPLX(A0+B1,B0-A1)
  110 CONTINUE
      GO TO 190
  115 IF (L1 .EQ. 1) GO TO 120
      C2 = C1 * C1 - S1 * S1
      S2 = TWO * C1 * S1
      C3 = C1 * C2 - S1 * S2
      S3 = S1 * C2 + C1 * S2
  120 JA = KB + ISP - 1
      KA = JA + KB
      IKB = KB+1
      IJA = JA+1
      DO 135 II = IKB,IJA
         K0 = KA - II + 1
         K1 = K0 + ISP
         K2 = K1 + ISP
         K3 = K2 + ISP
         ZA0 = A(K0+1)
         IF (L1 .EQ. 1) GO TO 125
         ZA4 = A(K1+1)
         A1 = A4*C1-B4*S1
         B1 = A4*S1+B4*C1
         ZA4 = A(K2+1)
         A2 = A4*C2-B4*S2
         B2 = A4*S2+B4*C2
         ZA4 = A(K3+1)
         A3 = A4*C3-B4*S3
         B3 = A4*S3+B4*C3
         GO TO 130
  125    ZA1 = A(K1+1)
         ZA2 = A(K2+1)
         ZA3 = A(K3+1)
  130    A(K0+1) = CMPLX(A0+A2+A1+A3,B0+B2+B1+B3)
         A(K1+1) = CMPLX(A0+A2-A1-A3,B0+B2-B1-B3)
         A(K2+1) = CMPLX(A0-A2-B1+B3,B0-B2+A1-A3)
         A(K3+1) = CMPLX(A0-A2+B1-B3,B0-B2-A1+A3)
  135 CONTINUE
      GO TO 190
  140 JK = KF - 1
      KH = JK/2
      K3 = IWK(ID+I-1)
      K0 = KB + ISP
      IF (L1 .EQ. 1) GO TO 150
      K = JK - 1
      WK(ICF+1) = C1
      WK(ISF+1) = S1
      DO 145 J = 1,K
         WK(ICF+J+1) = WK(ICF+J) * C1 - WK(ISF+J) * S1
         WK(ISF+J+1) = WK(ICF+J) * S1 + WK(ISF+J) * C1
  145 CONTINUE
  150 IF (KF .EQ. JF) GO TO 160
      C2 = WK(ICC+I)
      WK(ICK+1) = C2
      WK(ICK+JK) = C2
      S2 = WK(ISS+I)
      WK(ISK+1) = S2
      WK(ISK+JK) = -S2
      DO 155 J = 1,KH
         K = JK - J
         WK(ICK+K) = WK(ICK+J) * C2 - WK(ISK+J) * S2
         WK(ICK+J+1) = WK(ICK+K)
         WK(ISK+J+1) = WK(ICK+J) * S2 + WK(ISK+J) * C2
         WK(ISK+K) = -WK(ISK+J+1)
  155 CONTINUE
  160 K0 = K0 - 1
      K1 = K0
      K2 = K0 + K3
      ZA0 = A(K0+1)
      A3 = A0
      B3 = B0
      DO 175 J = 1,KH
         K1 = K1 + ISP
         K2 = K2 - ISP
         IF (L1 .EQ. 1) GO TO 165
         K = KF - J
         ZA4 = A(K1+1)
         A1 = A4*WK(ICF+J)-B4*WK(ISF+J)
         B1 = A4*WK(ISF+J)+B4*WK(ICF+J)
         ZA4 = A(K2+1)
         A2 = A4*WK(ICF+K)-B4*WK(ISF+K)
         B2 = A4*WK(ISF+K)+B4*WK(ICF+K)
         GO TO 170
  165    ZA1 = A(K1+1)
         ZA2 = A(K2+1)
  170    WK(IAP+J) = A1 + A2
         WK(IAM+J) = A1 - A2
         WK(IBP+J) = B1 + B2
         WK(IBM+J) = B1 - B2
         A3 = A1 + A2 + A3
         B3 = B1 + B2 + B3
  175 CONTINUE
      A(K0+1) = CMPLX(A3,B3)
      K1 = K0
      K2 = K0 + K3
      DO 185 J = 1,KH
         K1 = K1 + ISP
         K2 = K2 - ISP
         JK = J
         A1 = A0
         B1 = B0
         A2 = ZERO
         B2 = ZERO
         DO 180  K = 1,KH
            A1 = A1 + WK(IAP+K) * WK(ICK+JK)
            A2 = A2 + WK(IAM+K) * WK(ISK+JK)
            B1 = B1 + WK(IBP+K) * WK(ICK+JK)
            B2 = B2 + WK(IBM+K) * WK(ISK+JK)
            JK = JK + J
            IF (JK .GE. KF) JK = JK - KF
  180    CONTINUE
         A(K1+1) = CMPLX(A1-B2,B1+A2)
         A(K2+1) = CMPLX(A1+B2,B1-A2)
  185 CONTINUE
      IF (K0 .GT. KB) GO TO 160
      JF = KF
  190 IF ( I .GE. MM ) GO TO 195
      I = I + 1
      GO TO 55
  195 I = MM
      L1 = 0
      KB = IWK(ID+I-1) + KB
      IF (KB .GE. KN) GO TO 215
  200 JJ = IWK(IC+I-2) + JJ
      IF (JJ .LT. IWK(IC+I-1)) GO TO 205
      I = I - 1
      JJ = JJ - IWK(IC+I)
      GO TO 200
  205 IF (I .NE. MM) GO TO 210
      C2 = C1
      C1 = CM * C1 - SM * S1
      S1 = SM * C2 + CM * S1
      GO TO 70
  210 IF (IWK(ILL+I) .EQ. 1) I = I + 1
      GO TO 55
  215 I = 1
      JA = KT - 1
      KA = JA + 1
      IF(JA.LT.1) GO TO 225
      DO 220  II = 1,JA
         J = KA - II
         IWK(J+1) = IWK(J+1) - 1
         I = IWK(J+1) + I
  220 CONTINUE
!                                  THE RESULT IS NOW PERMUTED TO
!                                  NORMAL ORDER.
  225 IF (KT .LE. 0) GO TO 270
      J = 1
      I = 0
      KB = 0
  230 K2 = IWK(ID+J) + KB
      K3 = K2
      JJ = IWK(IC+J-1)
      JK = JJ
      K0 = KB + JJ
      ISP = IWK(IC+J) - JJ
  235 K = K0 + JJ
  240 ZA4 = A(K0+1)
      A(K0+1) = A(K2+1)
      A(K2+1) = ZA4
      K0 = K0 + 1
      K2 = K2 + 1
      IF (K0 .LT. K) GO TO 240
      K0 = K0 + ISP
      K2 = K2 + ISP
      IF (K0 .LT. K3) GO TO 235
      IF (K0 .GE. K3 + ISP) GO TO 245
      K0 = K0 - IWK(ID+J) + JJ
      GO TO 235
  245 K3 = IWK(ID+J) + K3
      IF (K3 - KB .GE. IWK(ID+J-1)) GO TO 250
      K2 = K3 + JK
      JK = JK + JJ
      K0 = K3 - IWK(ID+J) + JK
      GO TO 235
  250 IF (J .GE. KT) GO TO 260
      K = IWK(J+1) + I
      J = J + 1
  255 I = I + 1
      IWK(ILL+I) = J
      IF (I .LT. K) GO TO 255
      GO TO 230
  260 KB = K3
      IF (I .LE. 0) GO TO 265
      J = IWK(ILL+I)
      I = I - 1
      GO TO 230
  265 IF (KB .GE. N) GO TO 270
      J = 1
      GO TO 230
  270 JK = IWK(IC+KT)
      ISP = IWK(ID+KT)
      M = M - KT
      KB = ISP/JK-2
      IF (KT .GE. M-1 ) GO TO 9005
      ITA = ILL+KB+1
      ITB = ITA+JK
      IDM1 = ID-1
      IKT = KT+1
      IM = M+1
      DO 275 J = IKT,IM
         IWK(IDM1+J) = IWK(IDM1+J)/JK
  275 CONTINUE
      JJ = 0
      DO 290 J = 1,KB
         K = KT
  280    JJ = IWK(ID+K+1) + JJ
         IF (JJ .LT. IWK(ID+K)) GO TO 285
         JJ = JJ - IWK(ID+K)
         K = K + 1
         GO TO 280
  285    IWK(ILL+J) = JJ
         IF (JJ .EQ. J) IWK(ILL+J) = -J
  290 CONTINUE
!                                  DETERMINE THE PERMUTATION CYCLES
!                                  OF LENGTH GREATER THAN OR EQUAL
!                                  TO TWO.
      DO 300  J = 1,KB
         IF (IWK(ILL+J) .LE. 0) GO TO 300
         K2 = J
  295    K2 = IABS(IWK(ILL+K2))
         IF (K2 .EQ. J) GO TO 300
         IWK(ILL+K2) = -IWK(ILL+K2)
         GO TO 295
  300 CONTINUE
!                                  REORDER A FOLLOWING THE
!                                  PERMUTATION CYCLES
      I = 0
      J = 0
      KB = 0
      KN = N
  305 J = J + 1
      IF (IWK(ILL+J) .LT. 0) GO TO 305
      K = IWK(ILL+J)
      K0 = JK * K + KB
  310 ZA4 = A(K0+I+1)
      WK(ITA+I) = A4
      WK(ITB+I) = B4
      I = I + 1
      IF (I .LT. JK) GO TO 310
      I = 0
  315 K = -IWK(ILL+K)
      JJ = K0
      K0 = JK * K + KB
  320 A(JJ+I+1) = A(K0+I+1)
      I = I + 1
      IF (I .LT. JK) GO TO 320
      I = 0
      IF (K .NE. J) GO TO 315
  325 A(K0+I+1) = CMPLX(WK(ITA+I),WK(ITB+I))
      I = I + 1
      IF (I .LT. JK) GO TO 325
      I = 0
      IF (J .LT. K2) GO TO 305
      J = 0
      KB = KB + ISP
      IF (KB .LT. KN) GO TO 305
 9005 RETURN
      END
!
!-----------------------------------------------------------------------

!   PURPOSE             - COMPUTE THE FAST FOURIER TRANSFORM OF A
!                           REAL VALUED SEQUENCE
!
!   USAGE               - CALL FFTRC (A,N,X,IWK,WK)
!
!   ARGUMENTS    A      - INPUT REAL VECTOR OF LENGTH N WHICH
!                           CONTAINS THE DATA TO BE TRANSFORMED.
!                N      - INPUT NUMBER OF DATA POINTS TO BE TRANSFORMED.
!                            N MUST BE A POSITIVE EVEN INTEGER.
!                X      - OUTPUT COMPLEX VECTOR OF LENGTH N/2+1
!                           CONTAINING THE FIRST N/2+1 COEFFICIENTS OF
!                           THE FOURIER TRANSFORM. THE REMAINING
!                           COEFFICIENTS MAY BE DETERMINED BY
!                             X(N+2-I) = CONJG(X(I)), FOR I=2,...,N/2.
!                IWK    - INTEGER WORK VECTOR.
!                           IF N IS A POWER OF 2, THEN IWK SHOULD BE OF
!                           LENGTH M WHERE N=2**M.
!                           OTHERWISE, IWK SHOULD BE OF LENGTH
!                           6*(N/2)+150.
!                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS)
!                WK     - REAL WORK VECTOR OF LENGTH 6*(N/2)+150.
!                           WK IS NOT USED IF N IS A POWER OF 2.
!                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS)

!

!   REMARKS  1.  FFTRC COMPUTES THE FOURIER TRANSFORM, X, ACCORDING
!                TO THE FOLLOWING FORMULA;
!
!                  X(K+1) = SUM FROM J = 0 TO N-1 OF
!                           A(J+1)*CEXP((0.0,(2.0*PI*J*K)/N))
!                  FOR K=0,1,...,N/2 AND PI=3.1415...
!
!                THE USER CAN COMPUTE THE REMAINING X VALUES BY
!                PERFORMING THE FOLLOWING STEPS;
!
!                     ND2 = N/2
!                     DO 10 I=2,ND2
!                        X(N+2-I) = CONJG(X(I))
!                  10 CONTINUE
!            2.  FFTRC CAN BE USED TO COMPUTE
!
!                  X(K+1) = (1/N)*SUM FROM J = 0 TO N-1 OF
!                           A(J+1)*CEXP((0.0,(-2.0*PI*J*K)/N))
!                  FOR K=0,1,...,N/2 AND PI=3.1415...
!
!                BY PERFORMING THE FOLLOWING STEPS;
!
!                     CALL FFTRC (A,N,X,IWK,WK)
!                     ND2P1 = N/2+1
!                     DO 10 I=1,ND2P1
!                        X(I) = CONJG(X(I))/N
!                  10 CONTINUE
!

!-----------------------------------------------------------------------
!
      SUBROUTINE FFTRC  (A,N,X,IWK,WK)
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IWK(6*N+150)
      REAL               A(N),WK(6*N+150)
      COMPLEX            X(N)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ND2P1,ND2,I,MTWO,M,IMAX,ND4,NP2,K,NMK,J
      REAL               RPI,ZERO,ONE,HALF,THETA,TP,G(2),B(2),Z(2),AI, &
                         AR
      COMPLEX            XIMAG,ALPH,BETA,GAM,S1,ZD
      EQUIVALENCE        (GAM,G(1)),(ALPH,B(1)),(Z(1),AR),(Z(2),AI), &
                         (ZD,Z(1))
      DATA               ZERO/0.0/,HALF/0.5/,ONE/1.0/,IMAX/24/
      DATA               RPI/3.141593/
!                                  FIRST EXECUTABLE STATEMENT
      IF (N .NE. 2) GO TO 5
!                                  N EQUAL TO 2
      ZD = CMPLX(A(1),A(2))
      THETA = AR
      TP = AI
      X(2) = CMPLX(THETA-TP,ZERO)
      X(1) = CMPLX(THETA+TP,ZERO)
      GO TO 9005
    5 CONTINUE
!                                  N GREATER THAN 2
      ND2 = N/2
      ND2P1 = ND2+1
!                                  MOVE A TO X
      J = 1
      DO 6 I=1,ND2
         X(I) = CMPLX(A(J),A(J+1))
         J = J+2
    6 CONTINUE
!                                  COMPUTE THE CENTER COEFFICIENT
      GAM = CMPLX(ZERO,ZERO)
      DO 10 I=1,ND2
         GAM = GAM + X(I)
   10 CONTINUE
      TP = G(1)-G(2)
      GAM = CMPLX(TP,ZERO)
!                                  DETERMINE THE SMALLEST M SUCH THAT
!                                  N IS LESS THAN OR EQUAL TO 2**M
      MTWO = 2
      M = 1
      DO 15 I=1,IMAX
         IF (ND2 .LE. MTWO) GO TO 20
         MTWO = MTWO+MTWO
         M = M+1
   15 CONTINUE
   20 IF (ND2 .EQ. MTWO) GO TO 25
!                                  N IS NOT A POWER OF TWO, CALL FFTCC
      CALL FFTCC (X,ND2,IWK,WK)
      GO TO 30
!                                  N IS A POWER OF TWO, CALL FFT2C
   25 CALL FFT2C (X,M,IWK)
   30 ALPH = X(1)
      X(1) = B(1) + B(2)
      ND4 = (ND2+1)/2
      IF (ND4 .LT. 2) GO TO 40
      NP2 = ND2 + 2
      THETA = RPI/ND2
      TP = THETA
      XIMAG = CMPLX(ZERO,ONE)
!                                  DECOMPOSE THE COMPLEX VECTOR X
!                                  INTO THE COMPONENTS OF THE TRANSFORM
!                                  OF THE INPUT DATA.
      DO 35 K = 2,ND4
         NMK = NP2 - K
         S1 = CONJG(X(NMK))
         ALPH = X(K) + S1
         BETA = XIMAG*(S1-X(K))
         S1 = CMPLX(COS(THETA),SIN(THETA))
         X(K) = (ALPH+BETA*S1)*HALF
         X(NMK) = CONJG(ALPH-BETA*S1)*HALF
         THETA = THETA + TP
   35 CONTINUE
   40 CONTINUE
      X(ND2P1) = GAM
 9005 RETURN
      END
!
!-----------------------------------------------------------------------

!   PURPOSE             - COMPUTE THE FAST FOURIER TRANSFORM OF A
!                           COMPLEX VALUED SEQUENCE OF LENGTH EQUAL TO
!                           A POWER OF TWO
!
!   USAGE               - CALL FFT2C (A,M,IWK)
!
!   ARGUMENTS    A      - COMPLEX VECTOR OF LENGTH N, WHERE N=2**M.
!                           ON INPUT A CONTAINS THE COMPLEX VALUED
!                           SEQUENCE TO BE TRANSFORMED.
!                           ON OUTPUT A IS REPLACED BY THE
!                           FOURIER TRANSFORM.
!                M      - INPUT EXPONENT TO WHICH 2 IS RAISED TO
!                           PRODUCE THE NUMBER OF DATA POINTS, N
!                           (I.E. N = 2**M).
!                IWK    - WORK VECTOR OF LENGTH M+1.
!


!   REMARKS  1.  FFT2C COMPUTES THE FOURIER TRANSFORM, X, ACCORDING
!                TO THE FOLLOWING FORMULA;
!
!                  X(K+1) = SUM FROM J = 0 TO N-1 OF
!                           A(J+1)*CEXP((0.0,(2.0*PI*J*K)/N))
!                  FOR K=0,1,...,N-1 AND PI=3.1415...
!
!                NOTE THAT X OVERWRITES A ON OUTPUT.
!            2.  FFT2C CAN BE USED TO COMPUTE THE INVERSE FOURIER
!                TRANSFORM, X, ACCORDING TO THE FOLLOWING FORMULA;
!
!                  X(K+1) = (1/N)*SUM FROM J = 0 TO N-1 OF
!                           A(J+1)*CEXP((0.0,(-2.0*PI*J*K)/N))
!                  FOR K=0,1,...,N-1 AND PI=3.1415...
!
!                BY PERFORMING THE FOLLOWING STEPS;
!
!                     DO 10 I=1,N
!                        A(I) = CONJG(A(I))
!                  10 CONTINUE
!                     CALL FFT2C (A,M,IWK)
!                     DO 20 I=1,N
!                        A(I) = CONJG(A(I))/N
!                  20 CONTINUE
!

!-----------------------------------------------------------------------
!
      SUBROUTINE  FFT2C (A,M,IWK)
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M,IWK(1)
      COMPLEX            A(1)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,ISP,J,JJ,JSP,K,K0,K1,K2,K3,KB,KN,MK,MM,MP,N, &
                         N4,N8,N2,LM,NN,JK
      REAL               RAD,C1,C2,C3,S1,S2,S3,CK,SK,SQ,A0,A1,A2,A3, &
                         B0,B1,B2,B3,TWOPI,TEMP, &
                         ZERO,ONE,Z0(2),Z1(2),Z2(2),Z3(2)
      COMPLEX            ZA0,ZA1,ZA2,ZA3,AK2
      EQUIVALENCE        (ZA0,Z0(1)),(ZA1,Z1(1)),(ZA2,Z2(1)), &
                         (ZA3,Z3(1)),(A0,Z0(1)),(B0,Z0(2)),(A1,Z1(1)), &
                         (B1,Z1(2)),(A2,Z2(1)),(B2,Z2(2)),(A3,Z3(1)), &
                         (B3,Z3(2))
      DATA               SQ/.7071068/, &
                         SK/.3826834/, &
                         CK/.9238795/, &
                         TWOPI/6.283185/
      DATA               ZERO/0.0/,ONE/1.0/
!                                  SQ=SQRT2/2,SK=SIN(PI/8),CK=COS(PI/8)
!                                  TWOPI=2*PI
!                                  FIRST EXECUTABLE STATEMENT
      MP = M+1
      N = 2**M
      IWK(1) = 1
      MM = (M/2)*2
      KN = N+1
!                                  INITIALIZE WORK VECTOR
      DO 5  I=2,MP
         IWK(I) = IWK(I-1)+IWK(I-1)
    5 CONTINUE
      RAD = TWOPI/N
      MK = M - 4
      KB = 1
      IF (MM .EQ. M) GO TO 15
      K2 = KN
      K0 = IWK(MM+1) + KB
   10 K2 = K2 - 1
      K0 = K0 - 1
      AK2 = A(K2)
      A(K2) = A(K0) - AK2
      A(K0) = A(K0) + AK2
      IF (K0 .GT. KB) GO TO 10
   15 C1 = ONE
      S1 = ZERO
      JJ = 0
      K = MM - 1
      J = 4
      IF (K .GE. 1) GO TO 30
      GO TO 70
   20 IF (IWK(J) .GT. JJ) GO TO 25
      JJ = JJ - IWK(J)
      J = J-1
      IF (IWK(J) .GT. JJ) GO TO 25
      JJ = JJ - IWK(J)
      J = J - 1
      K = K + 2
      GO TO 20
   25 JJ = IWK(J) + JJ
      J = 4
   30 ISP = IWK(K)
      IF (JJ .EQ. 0) GO TO 40
!                                  RESET TRIGONOMETRIC PARAMETERS
      C2 = JJ * ISP * RAD
      C1 = COS(C2)
      S1 = SIN(C2)
   35 C2 = C1 * C1 - S1 * S1
      S2 = C1 * (S1 + S1)
      C3 = C2 * C1 - S2 * S1
      S3 = C2 * S1 + S2 * C1
   40 JSP = ISP + KB
!                                  DETERMINE FOURIER COEFFICIENTS
!                                    IN GROUPS OF 4
      DO 50 I=1,ISP
         K0 = JSP - I
         K1 = K0 + ISP
         K2 = K1 + ISP
         K3 = K2 + ISP
         ZA0 = A(K0)
         ZA1 = A(K1)
         ZA2 = A(K2)
         ZA3 = A(K3)
         IF (S1 .EQ. ZERO) GO TO 45
         TEMP = A1
         A1 = A1 * C1 - B1 * S1
         B1 = TEMP * S1 + B1 * C1
         TEMP = A2
         A2 = A2 * C2 - B2 * S2
         B2 = TEMP * S2 + B2 * C2
         TEMP = A3
         A3 = A3 * C3 - B3 * S3
         B3 = TEMP * S3 + B3 * C3
   45    TEMP = A0 + A2
         A2 = A0 - A2
         A0 = TEMP
         TEMP = A1 + A3
         A3 = A1 - A3
         A1 = TEMP
         TEMP = B0 + B2
         B2 = B0 - B2
         B0 = TEMP
         TEMP = B1 + B3
         B3 = B1 - B3
         B1 = TEMP
         A(K0) = CMPLX(A0+A1,B0+B1)
         A(K1) = CMPLX(A0-A1,B0-B1)
         A(K2) = CMPLX(A2-B3,B2+A3)
         A(K3) = CMPLX(A2+B3,B2-A3)
   50 CONTINUE
      IF (K .LE. 1) GO TO 55
      K = K - 2
      GO TO 30
   55 KB = K3 + ISP
!                                  CHECK FOR COMPLETION OF FINAL
!                                    ITERATION
      IF (KN .LE. KB) GO TO 70
      IF (J .NE. 1) GO TO 60
      K = 3
      J = MK
      GO TO 20
   60 J = J - 1
      C2 = C1
      IF (J .NE. 2) GO TO 65
      C1 = C1 * CK + S1 * SK
      S1 = S1 * CK - C2 * SK
      GO TO 35
   65 C1 = (C1 - S1) * SQ
      S1 = (C2 + S1) * SQ
      GO TO 35
   70 CONTINUE
!                                  PERMUTE THE COMPLEX VECTOR IN
!                                    REVERSE BINARY ORDER TO NORMAL
!                                    ORDER
      IF(M .LE. 1) GO TO 9005
      MP = M+1
      JJ = 1
!                                  INITIALIZE WORK VECTOR
      IWK(1) = 1
      DO 75  I = 2,MP
         IWK(I) = IWK(I-1) * 2
   75 CONTINUE
      N4 = IWK(MP-2)
      IF (M .GT. 2) N8 = IWK(MP-3)
      N2 = IWK(MP-1)
      LM = N2
      NN = IWK(MP)+1
      MP = MP-4
!                                  DETERMINE INDICES AND SWITCH A
      J = 2
   80 JK = JJ + N2
      AK2 = A(J)
      A(J) = A(JK)
      A(JK) = AK2
      J = J+1
      IF (JJ .GT. N4) GO TO 85
      JJ = JJ + N4
      GO TO 105
   85 JJ = JJ - N4
      IF (JJ .GT. N8) GO TO 90
      JJ = JJ + N8
      GO TO 105
   90 JJ = JJ - N8
      K = MP
   95 IF (IWK(K) .GE. JJ) GO TO 100
      JJ = JJ - IWK(K)
      K = K - 1
      GO TO 95
  100 JJ = IWK(K) + JJ
  105 IF (JJ .LE. J) GO TO 110
      K = NN - J
      JK = NN - JJ
      AK2 = A(J)
      A(J) = A(JJ)
      A(JJ) = AK2
      AK2 = A(K)
      A(K) = A(JK)
      A(JK) = AK2
  110 J = J + 1
!                                  CYCLE REPEATED UNTIL LIMITING NUMBER
!                                    OF CHANGES IS ACHIEVED
      IF (J .LE. LM) GO TO 80
!
 9005 RETURN
      END
