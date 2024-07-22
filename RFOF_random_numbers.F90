module RFOF_random_numbers

#include "config.h"

contains

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  !--------------------------------------------------------------------------------
  function rand_uniform_var0mean1() result(y)
    real(8) :: y
    real(8) :: r(1)

    call RM48(r,1)
    y=sqrt(3d0)*2d0*(r(1)-0.5d0)

  end function rand_uniform_var0mean1


      SUBROUTINE RM48(RVEC,LENV)
! for 32-bit machines, use IMPLICIT DOUBLE PRECISION
!      INCLUDE 'dblprc'
!      INCLUDE 'dimpar'
!      INCLUDE 'iounit'
!
      IMPLICIT NONE
!
      INTEGER I,I24,I97,IDUM,II,IJ,IJKL,IJKLIN,IJKLUT,IOSEED,IVEC
      INTEGER J,J97,JJ,K,KALLED,KL,L,LENV,LOOP2,M,NOW
      INTEGER MODCNS,NTOT,NTOT2,NTOT2N,NTOTIN,NTOTUT,NTOT2T
      DOUBLE PRECISION C,CD,CM,HALF,ONE,RVEC,S,T,TWOM24,U,UNI,ZERO

      DIMENSION RVEC(*)
      COMMON/R48ST1/U(97),C,I97,J97
      PARAMETER (MODCNS=1000000000)
      SAVE CD, CM, TWOM24,  ZERO, ONE, NTOT, NTOT2, IJKL
      DATA NTOT,NTOT2,IJKL/-1,0,0/
!
      IF (NTOT .GE. 0)  GO TO 50
!
!        Default initialization. User has called RM48 without RM48IN.
      IJKL = 54217137
      NTOT = 0
      NTOT2 = 0
      KALLED = 0
      GO TO 1
!
      ENTRY      RM48IN(IJKLIN, NTOTIN,NTOT2N)
!         Initializing routine for RM48, may be called before
!         generating pseudorandom numbers with RM48.   The input
      IJKL = IJKLIN
      NTOT = MAX(NTOTIN,0)
      NTOT2= MAX(NTOT2N,0)
      KALLED = 1
!          always come here to initialize
    1 CONTINUE
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ, 177)     + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL, 169)
!     & ' RM48 INITIALIZED:',IJKL,NTOT,NTOT2
      ONE = 1.D+00
      HALF = 0.5D+00
      ZERO = 0.D+00
      DO 2 II= 1, 97
      S = 0.D+00
      T = HALF
      DO 3 JJ= 1, 48
         M = MOD(MOD(I*J,179)*K, 179)
         I = J
         J = K
         K = M
         L = MOD(53*L+1, 169)
         IF (MOD(L*M,64) .GE. 32)  S = S+T
    3    T = HALF*T
    2 U(II) = S
      TWOM24 = ONE
      DO 4 I24= 1, 24
    4 TWOM24 = HALF*TWOM24
      C  =   362436.D+00*TWOM24
      CD =  7654321.D+00*TWOM24
      CM = 16777213.D+00*TWOM24
      I97 = 97
      J97 = 33
!       Complete initialization by skipping
      DO 45 LOOP2= 1, NTOT2+1
      NOW = MODCNS
      IF (LOOP2 .EQ. NTOT2+1)  NOW=NTOT
      IF (NOW .GT. 0)  THEN
          DO 40 IDUM = 1, NTOT
          UNI = U(I97)-U(J97)
          IF (UNI .LT. ZERO)  UNI=UNI+ONE
          U(I97) = UNI
          I97 = I97-1
          IF (I97 .EQ. 0)  I97=97
          J97 = J97-1
          IF (J97 .EQ. 0)  J97=97
          C = C - CD
          IF (C .LT. ZERO)  C=C+CM
   40     CONTINUE
      ENDIF
   45 CONTINUE
      IF (KALLED .EQ. 1)  RETURN
!
!          Normal entry to generate LENV random numbers
   50 CONTINUE
      DO 100 IVEC= 1, LENV
      UNI = U(I97)-U(J97)
      IF (UNI .LT. ZERO)  UNI=UNI+ONE
      U(I97) = UNI
      I97 = I97-1
      IF (I97 .EQ. 0)  I97=97
      J97 = J97-1
      IF (J97 .EQ. 0)  J97=97
      C = C - CD
      IF (C .LT. ZERO)  C=C+CM
      UNI = UNI-C
      IF (UNI .LT. ZERO) UNI=UNI+ONE
      RVEC(IVEC) = UNI
!C             An exact zero here is very unlikely, but let's be safe.
  100 CONTINUE
      NTOT = NTOT + LENV
         IF (NTOT .GE. MODCNS)  THEN
         NTOT2 = NTOT2 + 1
         NTOT = NTOT - MODCNS
         ENDIF
      RETURN
!           Entry to output current status
      ENTRY RM48UT(IJKLUT,NTOTUT,NTOT2T)
      IJKLUT = IJKL
      NTOTUT = NTOT
      NTOT2T = NTOT2
      RETURN
!
      ENTRY      RM48WR(IOSEED)
!         Output routine for RM48, without skipping numbers
      WRITE (IOSEED,'(2Z8)') NTOT,NTOT2
      WRITE (IOSEED,'(2Z8,Z16)') I97,J97,C
      WRITE (IOSEED,'(24(4Z16,/),Z16)') U
      RETURN
!
      ENTRY      RM48RD(IOSEED)
!         Initializing routine for RM48, without skipping numbers
      READ (IOSEED,'(2Z8)') NTOT,NTOT2
      READ (IOSEED,'(2Z8,Z16)') I97,J97,C
      READ (IOSEED,'(24(4Z16,/),Z16)') U
      CLOSE (UNIT=IOSEED)
      IJKL = 54217137
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ, 177)     + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL, 169)
!     &  ' RM48 INITIALIZED:',IJKL,NTOT,NTOT2
      ONE  = 1.D+00
      HALF = 0.5D+00
      ZERO = 0.D+00
      TWOM24 = ONE
      DO 400 I24= 1, 24
  400 TWOM24 = HALF*TWOM24
      CD =  7654321.D+00*TWOM24
      CM = 16777213.D+00*TWOM24
      RETURN
      END SUBROUTINE


  !--------------------------------------------------------------------------------
  function gaussian_random_number() result(y)

    ! Output
    real(8) :: y

    ! Local
    real(8) :: RVEC(1)

    CALL RM48(RVEC,1)

    y = sqrt(3d0) * 2d0 * ( RVEC(1) - 0.5d0 )

  end function gaussian_random_number

end module RFOF_random_numbers
