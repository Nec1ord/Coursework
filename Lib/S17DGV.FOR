      SUBROUTINE S17DGV(Y,NZ,ASCLE,TOL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
C      EXP(-ALIM)=ASCLE=1.0E+3*X02AKF()/TOL. THE TEST IS MADE TO SEE
C      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
C      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
C      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
C      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
C      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
C
C     .. Scalar Arguments ..
      COMPLEX*16        Y
      DOUBLE PRECISION  ASCLE, TOL
      INTEGER           NZ
C     .. Local Scalars ..
      DOUBLE PRECISION  SS, ST, YI, YR
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, MAX, MIN, DBLE
C     .. Executable Statements ..
C
      NZ = 0
      YR = DBLE(Y)
      YI = DIMAG(Y)
      YR = ABS(YR)
      YI = ABS(YI)
      ST = MIN(YR,YI)
      IF (ST.LE.ASCLE) THEN
         SS = MAX(YR,YI)
         ST = ST/TOL
         IF (SS.LT.ST) NZ = 1
      END IF
      RETURN
      END
