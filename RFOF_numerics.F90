#include "config.h"

module RFOF_numerics

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  implicit none

contains


  !--------------------------------------------------------------------------------
   subroutine Ridders_derivative(FCN, x, h, rel_tol, err, fprim)

    ! Input
    real(8), intent(in) :: x        ! Value where the derivative should be evaluated
    real(8), intent(in) :: h        ! Initial step size
    real(8), intent(in) :: rel_tol  ! Relative tolerance

    ! Output
    real(8), intent(out) :: err
    real(8), intent(out) :: fprim

    ! Local
    real(8), allocatable :: a(:,:)
    real(8) :: hh, errt, fac
    integer :: i, j

    real(8), parameter :: BIG = 1d50
    real(8), parameter :: CON = 1.4d0
    real(8), parameter :: CON2 = CON*CON
    real(8), parameter :: SAFE = 2d0
    integer, parameter :: NTAB = 10

    interface 
       function FCN(x) result(y)
         real(8) :: x
         real(8) :: y
       end function FCN
    end interface

    allocate(a(NTAB,NTAB))
    hh=h
    a(1,1)=(FCN(x+hh)-FCN(x-hh))/(2d0*hh)
    err=BIG

    do i=2, NTAB

       !print *, i, err

       !Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of extrapolation.
       hh = hh / CON
       a(1,i)=(FCN(x+hh)-FCN(x-hh))/(2.0*hh)

       !Try new, smaller step-size. 
       fac=CON2;

       do j=2, i
          !Compute extrapolations of various orders, requiring
          a(j,i)=(a(j-1,i)*fac - a(j-1,i-1))/(fac-1.0)

          !no new function evaluations
          fac=CON2*fac
          errt=MAX( ABS(a(j,i)-a(j-1,i)), ABS(a(j,i)-a(j-1,i-1)) )

          !The error strategy is to compare each new extrapolation to one order lower, 
          !both at the present stepsize and the previous one. if (errt <= *err) { 
          !If error is decreased, save the improved answer.
          err=errt 
          fprim=a(j,i)
       enddo
       
       !print *, err/fprim

       !If higher order is worse by a significant factor SAFE, then quit early. 
       if ( (ABS(a(i,i)-a(i-1,i-1)) .ge. SAFE*err) .or. &
            (err .lt. rel_tol * fprim) ) then
          return
       endif
    enddo

    deallocate(a)

   end subroutine Ridders_derivative



  !--------------------------------------------------------------------------------
  !  subroutine polycoef_3points
  !--------------------------------------------------------------------------------

  subroutine polycoef_3points(t,x,coef,error)

    ! Input
    real(8), intent(in) :: t(3), x(3)

    ! Output
    real(8), intent(out) :: coef(3)
    integer, intent(out) :: error

    ! Local
    real(8) :: xbiss, xprim, x0
    real(8) :: xp1, xp2, y1, y3

    if ( abs(t(1)-t(2)) .eq. 0d0 .or. &
         abs(t(2)-t(3)) .eq. 0d0 .or. &
         abs(t(1)-t(3)) .eq. 0d0 ) then
       error = 1
       return
    endif

    error = 0

    xp1   = ( x(2) - x(1) ) / ( t(2) - t(1) )
    xp2   = ( x(3) - x(2) ) / ( t(3) - t(2) )
    xbiss = ( xp2  - xp1) * 2. / ( t(3) - t(1) )

    y1 = x(1) - 0.5*t(1)**2 * xbiss
    y3 = x(3) - 0.5*t(3)**2 * xbiss

    xprim  = (y3-y1) / ( t(3) - t(1) )

    x0 = y1 - xprim*t(1)

    coef(1) = x0
    coef(2) = xprim
    coef(3) = 0.5 * xbiss

    !    print *, 'polycoef', xp1,xp2,t

  end subroutine polycoef_3points


  !--------------------------------------------------------------------------------
  !  subroutine solve_quadratic_polynomial_real_roots
  !--------------------------------------------------------------------------------

  subroutine solve_quadratic_polynomial_real_roots(coef,root1,root2,Nsolutions)

    ! Input
    real(8), intent(in) :: coef(3)

    ! Output
    real(8), intent(out) :: root1, root2
    integer, intent(out) :: Nsolutions

    ! Local
    real(8) :: D, sqrtD

    ! In case "coef(3)==0" then the equation is linear
    if ( abs(coef(3)) .lt. 1e-30) then
       if ( abs(coef(2)) .lt. 1e-30) then
          Nsolutions = 0
       else

          !print *, '2nd order solve; linear:',coef
          root1 = - coef(1)/coef(2)
          Nsolutions = 1
       endif
       return
    endif

    D = coef(2)*coef(2) - 4.0*coef(3)*coef(1)
    if (D .lt. 0.) then
       Nsolutions = 0
       return
    endif
    Nsolutions = 2

    sqrtD = sqrt(D)

    root1 = 0.5 * (-coef(2) + sqrt(D)) / coef(3)
    root2 = 0.5 * (-coef(2) - sqrt(D)) / coef(3)

  end subroutine solve_quadratic_polynomial_real_roots


  !--------------------------------------------------------------------------------
  !  subroutine eval_quadratic_polynomial
  !--------------------------------------------------------------------------------

  function eval_quadratic_polynomial(coef,t) result(x)

    ! Input
    real(8), intent(in) :: coef(3), t

    ! Output
    real(8) :: x

    x = coef(1) + coef(2) * t + coef(3) * t*t

  end function eval_quadratic_polynomial



  !--------------------------------------------------------------------------------
  !  subroutine eval_derivative_of_quadratic_polynomial
  !--------------------------------------------------------------------------------

  function eval_derivative_of_quadratic_polynomial(coef,t) result(x)

    ! Input
    real(8), intent(in) :: coef(3), t

    ! Output
    real(8) :: x

    x = coef(2) + 2d0  * coef(3) * t

  end function eval_derivative_of_quadratic_polynomial


  !******************************************************************************
  !
  !  function AIRY_SPECIAL_ARG(PX) result(AAIRY)
  !
  !     Approximate absolute value of Airy function
  !     for the special type of complex arguments:
  !        (1/2 - sqrt(3)/2*i)*abs(PX).
  !
  !******************************************************************************
  function AIRY_SPECIAL_ARG(PX) result(AAIRY)

    ! Input
    real(8), intent(in) :: PX

    ! Output
    real(8) :: AAIRY

    ! Local
    real(8) :: ZX,ZZ,ZAIRY

    ZX=ABS(PX)

    IF (ZX.LE.0.6D0)THEN
       AAIRY=0.355028052626D+00 + ZX*( &
            -0.12940957613D+00 + ZX*( &
            0.0707529220372D+00 + ZX*( &
            -0.0333528064715D+00 + ZX*( &
            0.0130038946939D+00 + ZX*( &
            -0.00395567028802D+00 + ZX*( &
            0.000724240169477D+00))))))
       RETURN
    ENDIF

    IF (ZX.LE.1.3D0)THEN
       AAIRY=0.354988229372D+00 + ZX*( &
            -0.129085408213D+00 + ZX*( &
            0.0696258262121D+00 + ZX*( &
            -0.0311826451966D+00 + ZX*( &
            0.0105258764479D+00 + ZX*( &
            -0.00234196202216D+00 + ZX*( &
            0.000252467277846D+00))))))
       RETURN
    ENDIF

    IF (ZX.LE.2.0D0)THEN
       AAIRY=0.352168834148D+00 + ZX*( &
            -0.118022458348D+00 + ZX*( &
            0.0513455951605D+00 + ZX*( &
            -0.0148348531668D+00 + ZX*( &
            0.00214258956552D+00 + ZX*( &
            6.38563856304D-06 + ZX*( &
            -2.94389942827D-05))))))
       RETURN
    ENDIF

    IF (ZX.LE.5.0D0)THEN
       ZZ=1.D0+2.D0*ZX**0.75
       ZAIRY=0.3550280D0/ZZ**0.333333D0

       AAIRY=ZAIRY + &
            0.0434619905813D+00 + ZX*( &
            -0.0224540300024D+00 + ZX*( &
            0.00717044704144D+00 + ZX*( &
            -0.00144369232055D+00 + ZX*( &
            0.00017825427158D+00 + ZX*( &
            -1.22716747708D-05 + ZX*( &
            3.57995227599D-07))))))
       RETURN
    ENDIF


    IF (ZX.LE.15.0D0)THEN
       ZZ=1.D0+2.D0*ZX**0.75
       ZAIRY=0.3550280D0/ZZ**0.333333D0

       AAIRY=ZAIRY + 1.D0/( &
            14.6761321409D+00 + ZX*( &
            19.2153006129D+00 + ZX*( &
            0.383545352445D+00 + ZX*( &
            -0.0518094348427D+00 + ZX*( &
            0.00345465900274D+00  + ZX*( &
            -0.000121904306268D+00 + ZX*( &
            1.7875509282D-06)))))))
       RETURN
    ENDIF

    IF (ZX.LE.80.0D0)THEN
       ZZ=1.D0+2.D0*ZX**0.75D0
       ZAIRY=0.3550280D0/ZZ**0.333333D0

       AAIRY=ZAIRY + 1.D0/( &
            10.8684533775D+00 + ZX*( &
            21.2543350522D+00 + ZX*( &
            -0.0643027495275D+00 + ZX*( &
            0.000622205098646D+00 + ZX*( &
            -6.17628680505D-06 + ZX*( &
            3.86095322274D-08 + ZX*( &
            -1.06389328913D-10)))))))
       RETURN
    ENDIF


    IF (ZX.LE.500.0D0)THEN
       ZZ=1.D0+2.D0*ZX**0.75
       ZAIRY=0.3550280D0/ZZ**0.333333D0

       AAIRY=ZAIRY + 1.D0/( &
            36.1627244264D+00 + ZX*( &
            19.9454634004D+00 + ZX*( &
            -0.030717930085D+00 + ZX*( &
            6.41826065885D-05 + ZX*( &
            -1.07923135637D-07 + ZX*( &
            1.10910364144D-10 + ZX*( &
            -5.00429643327D-14)))))))
       RETURN
    ENDIF

    IF (ZX.LE.2000.0D0)THEN
       ZZ=1.D0+2.D0*ZX**0.75
       ZAIRY=0.3550280D0/ZZ**0.333333D0

       AAIRY=ZAIRY + 1.D0/( &
            451.117628726D+00 + ZX*( &
            15.8143302253D+00 + ZX*( &
            -0.0111008684172D+00 + ZX*( &
            6.98608370923D-06 + ZX*( &
            -2.9894045035D-09 + ZX*( &
            7.40926179563D-13 + ZX*( &
            -7.94638712766D-17)))))))
       RETURN
    ENDIF

    IF (ZX.LE.10000.0D0)THEN
       ZZ=1.D0+2.D0*ZX**0.75
       ZAIRY=0.3550280D0/ZZ**0.333333D0

       AAIRY=ZAIRY + 1.D0/( &
            3697.19614782D+00 + ZX*( &
            7.80398364788D+00 + ZX*( &
            -0.00167518652892D+00 + ZX*( &
            2.67434934504D-07 + ZX*( &
            -2.66237264056D-11 + ZX*( &
            1.46810677768D-15 + ZX*( &
            -3.41255618292D-20)))))))
       RETURN
    ENDIF

    ZZ=1.D0+2.D0*ZX**0.75D0
    AAIRY=0.3550280D0/ZZ**0.333333D0+7.D-6


    !      IF (PX .GT. 1.5448) THEN
    !!     Descending approximation "Abramowitz and Stegun, eq 10.4.59, page 448"
    !         ZXI=0.6666667E0*PX*SQRT(PX)
    !         ZSUM=1.E0-0.06944444E0/ZXI
    !         AAIRY=0.2829497E0*EXP(-ZXI)*ZSUM*PX**-0.25E0
    !      ELSE
    !!     Ascending approximation "Abramowitz and Stegun, eq 10.4.2, page 446"
    !         ZF=1.E0+0.1666666E0*PX**3 + 5.55555556E-3*PX**6
    !         ZG=PX  +0.0833333E0*PX**4 + 1.98412698E-3*PX**7
    !         AAIRY=0.3550280E0*ZF-0.2588194E0*ZG
    !      ENDIF

  end function AIRY_SPECIAL_ARG

  !******************************************************************************
  ! End of *AAIRY*
  !******************************************************************************


  !************************************************************************
  !*                                                                      *
  !*    Program to calculate the first kind Bessel function of integer    *
  !*    order N, for any REAL X, using the function BESSELJ(N,X).         *
  !*                                                                      *
  !* -------------------------------------------------------------------- *
  !*   Reference: From Numath Library By Tuan Dang Trong in Fortran 77.   *
  !*                                                                      *
  !*                               F90 Release 1.0 By J-P Moreau, Paris.  *
  !************************************************************************
  FUNCTION BESSELJ(N,X) result(BESSJ)

    !     This subroutine calculates the first kind modified Bessel function
    !     of integer order N, for any REAL X. We use here the classical
    !     recursion formula, when X > N. For X < N, the Miller's algorithm
    !     is used to avoid overflows. 
    !     REFERENCE:
    !     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
    !     MATHEMATICAL TABLES, VOL.5, 1962.

    ! Input
    integer, intent(in) :: N
    real(8), intent(in) :: X

    ! Output
    real(8) :: BESSJ

    ! Local
    INTEGER :: IACC, J, JSUM, M
    REAL(8) :: BIGNO, BIGNI
    PARAMETER (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
    REAL(8) BESSJ0,BESSJ1,TOX,BJM,BJ,BJP,SUM

    IF (N.EQ.0) THEN
       BESSJ = BESSELJ0(X)
       RETURN
    ENDIF
    IF (N.EQ.1) THEN
       BESSJ = BESSELJ1(X)
       RETURN
    ENDIF
    IF (X.EQ.0.) THEN
       BESSJ = 0.
       RETURN
    ENDIF
    TOX = 2./X
    IF (X.GT.FLOAT(N)) THEN
       BJM = BESSELJ0(X)
       BJ  = BESSELJ1(X)
       DO J = 1,N-1
          BJP = J*TOX*BJ-BJM
          BJM = BJ
          BJ  = BJP
       ENDDO
       BESSJ = BJ
    ELSE
       M = 2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
       BESSJ = 0.
       JSUM = 0
       SUM = 0.
       BJP = 0.
       BJ  = 1.
       DO J = M,1,-1
          BJM = J*TOX*BJ-BJP
          BJP = BJ
          BJ  = BJM
          IF (ABS(BJ).GT.BIGNO) THEN
             BJ  = BJ*BIGNI
             BJP = BJP*BIGNI
             BESSJ = BESSJ*BIGNI
             SUM = SUM*BIGNI
          ENDIF
          IF (JSUM.NE.0) SUM = SUM+BJ
          JSUM = 1-JSUM
          IF (J.EQ.N) BESSJ = BJP
       ENDDO
       SUM = 2.*SUM-BJ
       BESSJ = BESSJ/SUM
    ENDIF
    RETURN
  END FUNCTION BESSELJ
   
  ! ---------------------------------------------------------------------------
  FUNCTION BESSELJ0(X) result(BESSJ0)

    !     This subroutine calculates the First Kind Bessel Function of
    !     order 0, for any real number X. The polynomial approximation by
    !     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
    !     REFERENCES:
    !     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
    !     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
    !     VOL.5, 1962.

    ! Input
    REAL(8) X

    ! Output
    REAL(8) BESSJ0

    ! Local
    REAL(8) AX,FR,FS,Z,FP,FQ,XX
    REAL(8) Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  
    REAL(8) Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6

    DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
         -.2073370639D-5,.2093887211D-6 /
    DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
         -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
    DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
         651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
    DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
         9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /

    IF(X.EQ.0.D0) GO TO 1
    AX = ABS (X)
    IF (AX.LT.8.) THEN
       Y = X*X
       FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
       FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
       BESSJ0 = FR/FS
    ELSE
       Z = 8./AX
       Y = Z*Z
       XX = AX-.785398164
       FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
       FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
       BESSJ0 = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
    ENDIF
    RETURN
1   BESSJ0 = 1.D0
    RETURN
  END FUNCTION BESSELJ0

  ! ---------------------------------------------------------------------------
  FUNCTION BESSELJ1(X) result(BESSJ1)

    !     This subroutine calculates the First Kind Bessel Function of
    !     order 1, for any real number X. The polynomial approximation by
    !     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
    !     REFERENCES:
    !     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
    !     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
    !     VOL.5, 1962.

    ! Input
    REAL(8) X

    ! Output
    REAL(8) BESSJ1

    ! Local
    REAL(8) AX,FR,FS,Z,FP,FQ,XX
    REAL(8) Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  
    REAL(8) Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6

    DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
         .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
    DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
         .8449199096D-5,-.88228987D-6,.105787412D-6 /
    DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, & 
         242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
    DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
         18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /

    AX = ABS(X)
    IF (AX.LT.8.) THEN
       Y = X*X
       FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
       FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
       BESSJ1 = X*(FR/FS)
    ELSE
       Z = 8./AX
       Y = Z*Z
       XX = AX-2.35619491
       FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
       FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
       BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
    ENDIF
    RETURN
  END FUNCTION BESSELJ1
  !End of file bessel functions
  !************************************************************************



  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine add_to_sorted_list_3columns(list, nelements, t, x, y)

    ! Input
    real(8), intent(in) :: t,x,y

    ! Input/Output
    real(8), pointer, intent(inout) :: list(:,:)
    integer, intent(inout) :: nelements

    ! Local
    integer :: j, jcnt

    !   print *, "in add_to_list...", nelements, t, x, y
    !   print *, list

    do jcnt=1,nelements
       j = nelements-jcnt+1

       if (t .gt. list(j,1)) then
          list(j+1,1) = t
          list(j+1,2) = x
          list(j+1,3) = y
          nelements = nelements + 1
          return
       else
          list(j+1,1) = list(j,1)
          list(j+1,2) = list(j,2)
          list(j+1,3) = list(j,3)
       endif
    enddo
    list(1,1) = t
    list(1,2) = x
    list(1,3) = y
    nelements = nelements + 1

  end subroutine add_to_sorted_list_3columns

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine test_list()

    integer :: nelements, j, k
    real(8), pointer :: list(:,:)
    real(8) :: t, x, y

    allocate(list(10,3))
    nelements=0

    do j=1,10

       t=(real(j)-3.3)**2
       x=1.
       y=1.

       call add_to_sorted_list_3columns(list, nelements, t, x, y)
       !print *, nelements, t
       !print *, (list(k,1), k=1,nelements)
       !                do k=1,nelements
       !                print *, list(k,1), list(k,2), list(k,3)
       !                enddo
    enddo
    deallocate(list)

  end subroutine test_list


end module RFOF_numerics
