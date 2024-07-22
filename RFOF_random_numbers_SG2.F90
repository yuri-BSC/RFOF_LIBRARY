module RFOF_random_numbers_SG2
  use RFOF_random_numbers
contains


  subroutine init_RANDOM_NUMBER( init_parameter) 
    implicit none
    intrinsic random_seed, random_number, system_clock
    integer(4), intent(in):: init_parameter
    INTEGER, ALLOCATABLE :: seed(:), gseed(:) 
    integer(4) ::ic4, crate4, cmax4
    integer::isize, size
    real:: harvest(7)

    call random_seed(SIZE=size)  
    ALLOCATE (seed(size)); ALLOCATE (gseed(size))

    IF(init_parameter.EQ.0) THEN  
       call system_clock( ic4,  crate4,  cmax4)

    ElSE
       ic4=init_parameter
    ENDIF
    do isize=1, size
       seed(isize)= ic4 
    enddo
    call random_seed(PUT=seed(1:size)) 
    call random_seed(GET=gseed(1:size)) 

    !Print *, '  getseed=', gseed
  end subroutine init_RANDOM_NUMBER
  !----------------------------------------------
  subroutine init1_RM48( init_parameter) 
    implicit none
    intrinsic   system_clock
    integer, intent(in):: init_parameter

    integer:: IJKLIN, NTOTIN, NTOT2N
    integer(4) ::ic4, crate4, cmax4 
    IF(init_parameter.EQ.0) THEN  
       call system_clock( ic4,  crate4,  cmax4)
       cmax4=32000; crate4=1003
       IJKLIN=ic4; NTOTIN=2; NTOT2N=3
    ElSE
       IJKLIN=init_parameter
    ENDIF
    !print *,'RND_INIT!! count-clock, count  rate, maxcount  ', ic4, crate4, cmax4
    call RM48IN(IJKLIN, NTOTIN,NTOT2N)
  end subroutine init1_RM48


  !____________________________________

  !--------------------------------------------------------------------------------
  subroutine gaussian_random_vectorRM48(grw)

    implicit none
    real(8), intent(out)::grw(2)
    real(8):: r(2),v1,v2,w !, v(2)
1   w=2.d0
    do  while (w>1.d0)
       call RM48(r,2)
       v1=2.d0*r(1)-1.d0;
       v2=2.d0*r(2)-1.d0;  
       w=v1*v1+v2*v2
    end do
    if(w.EQ.0) GO TO 1
    w=sqrt(-2.d0*log(w)/w)
    !grw(1)=v(1)*w; grw(2)=v(2)*w
    grw(1)=v1*w; grw(2)=v2*w
  end subroutine gaussian_random_vectorRM48

  !----------------------------------------
  !--------------------------------------------------------------------------------
  subroutine gaussian_random_vector(grw)

    implicit none
    real(8), intent(out)::grw(2)
    !real(8):: r(2)
    real(8)::v1,v2,w 
    real::r(2) 
1   w=2.d0
    do  while (w>1.d0)
       call random_number(HARVEST=r)
       ! call RM48(r,2)
       v1=2.d0*r(1)-1.d0;
       v2=2.d0*r(2)-1.d0;  
       w=v1*v1+v2*v2
    end do
    if(w.EQ.0) GO TO 1
    w=sqrt(-2.d0*log(w)/w)
    !grw(1)=v(1)*w; grw(2)=v(2)*w
    grw(1)=v1*w; grw(2)=v2*w
  end subroutine gaussian_random_vector



  !-----------------------------------------


  !-----------------------------------------------------------

  subroutine Marsaglia_Bray_RM48( X )

    implicit none
    real(8), intent(out)::X 
    real(8)::   sum,UU(1),U , V, VV(1), W, XU(2),VW(2)
    real(8)::f0p8638, f0p9745,f0p9973002,f0p9986501
    real(8)::f1, f1p5,f2, f2p3153508, f3,f4p5, f6
    real(8)::f6p0432809, f6p6313339,f9p0334237,f13p2626678, f49p0024445

    data f0p8638/0.8638d0/, f0p9745/0.9745d0/
    data   f0p9973002/0.9973002d0/, f0p9986501/0.9986501d0/
    data f1/1.d0/,f1p5/1.5d0/, f2/2.d0/,f2p3153508 /2.3153508d0/,f3/3.0d0/,f4p5/4.5d0/ 
    data f6/6.d0/, f6p0432809/6.0432809d0/, f6p6313339/6.6313339d0/,f9p0334237/9.0334237d0/
    data  f13p2626678/13.2626678d0/, f49p0024445/49.0024445d0/

    call RM48(UU,1)
    U=UU(1);  
    if(U.LE.(f0p8638)) THEN
       call RM48(VW,2)
       V=f2*VW(1)-f1 ; W=f2*VW(2)-f1
       X=f2p3153508*U-f1+V+W; return
    elseif(U.LE.f0p9745) THEN
       call RM48(VV,1)
       X=f1p5*(VV(1)-f1+f9p0334237*(U-f0p8638)); return
    elseif(U.LE.f0p9973002) THEN
       do 
          call RM48(XU,2)
          X=f6*XU(1)-f3; U=XU(2)
          V=dabs(x); W=f6p6313339*(f3-V)**2
          if(v<f1p5) then
             sum=f6p0432809*(f1p5-v) 
             if (v<f1) then
                sum=sum +f13p2626678*(f3-v**2)-W 
             endif
          else 
             sum=0  
          endif
          if(u<f49p0024445*exp(-v**2/f2)-sum-w) exit
       enddo
       return
    else       ! U>9973002d0 , very low probability
       do
          call RM48(VW,2)  
          V=VW(1); W=VW(2)
          if(W.EQ.0) cycle
          X=f4p5 -log(w)
          if(X*V**2<f4p5) exit
       enddo
       if(U>f0p9986501) then
          X=sqrt(f2*X) 
       else
          X=-sqrt(f2*X)  
       endif
    endif
    return
  end subroutine Marsaglia_Bray_RM48


  subroutine MarsagliaBray( X )

    implicit none
    real(8), intent(out)::X 
    real(8)::   sum,UU(1),U , V, VV(1), W, XU(2),VW(2)
    real(8)::f0p8638, f0p9745,f0p9973002,f0p9986501
    real(8)::f1, f1p5,f2, f2p3153508, f3,f4p5, f6
    real(8)::f6p0432809, f6p6313339,f9p0334237,f13p2626678, f49p0024445

    data f0p8638/0.8638d0/, f0p9745/0.9745d0/
    data   f0p9973002/0.9973002d0/, f0p9986501/0.9986501d0/
    data f1/1.d0/,f1p5/1.5d0/, f2/2.d0/,f2p3153508 /2.3153508d0/,f3/3.0d0/,f4p5/4.5d0/ 
    data f6/6.d0/, f6p0432809/6.0432809d0/, f6p6313339/6.6313339d0/,f9p0334237/9.0334237d0/
    data  f13p2626678/13.2626678d0/, f49p0024445/49.0024445d0/

    !call RM48(UU,1)
    call random_number(HARVEST=UU)
    U=UU(1);  
    if(U.LE.(f0p8638)) THEN
       call RM48(VW,2)
       V=f2*VW(1)-f1 ; W=f2*VW(2)-f1
       X=f2p3153508*U-f1+V+W; return
    elseif(U.LE.f0p9745) THEN
       call RM48(VV,1)
       X=f1p5*(VV(1)-f1+f9p0334237*(U-f0p8638)); return
    elseif(U.LE.f0p9973002) THEN
       do 
          !call RM48(XU,2)
          call random_number(HARVEST=XU)
          X=f6*XU(1)-f3; U=XU(2)
          V=dabs(x); W=f6p6313339*(f3-V)**2
          if(v<f1p5) then
             sum=f6p0432809*(f1p5-v) 
             if (v<f1) then
                sum=sum +f13p2626678*(f3-v**2)-W 
             endif
          else 
             sum=0  
          endif
          if(u<f49p0024445*exp(-v**2/f2)-sum-w) exit
       enddo
       return
    else       ! U>9973002d0 , very low probability
       do
          !call RM48(VW,2) 
          call random_number(HARVEST=VW) 
          V=VW(1); W=VW(2)
          if(W.EQ.0) cycle
          X=f4p5 -log(w)
          if(X*V**2<f4p5) exit
       enddo
       if(U>f0p9986501) then
          X=sqrt(f2*X) 
       else
          X=-sqrt(f2*X)  
       endif
    endif
    return
  end subroutine MarsagliaBray
  !-----------------------------------------------------------------------

  subroutine Correlated_Gauss_Gen( delta_t,normal_vect, delta_w, delta_Z   )
    implicit none
    real(8),intent(in):: delta_t !>   the time step
    real(8),intent(in)::normal_vect(2)!> N(0,1) independent random numbers, denoted by U(2) in Platen,Kloeden

    real(8),intent(out)::delta_w, delta_Z 
    real(8) ::sqrt_delta_t !> $\sqrt(\Delta t)$, the square- root of the time step
    real(8)::  sqrt_delta_t3 !> $(\Delta t)^{3/2}$ 
    real(8)::half, inv_sqrt3 !$1/2; \ / 1/ \sqrt(3)$
    data half /0.5d0/, inv_sqrt3 /0.57735026918962576451D00/ 
    sqrt_delta_t=sqrt(delta_t)
    sqrt_delta_t3=sqrt_delta_t**3  
    delta_w=normal_vect(1)*sqrt_delta_t;
    delta_Z=half*sqrt_delta_t3*(normal_vect(1)+inv_sqrt3*normal_vect(2))
    return
  end subroutine Correlated_Gauss_Gen



  subroutine PolarGaussAccuracyTest1( nit)
    implicit none
    real(8), intent(in) ::nit 

    real(8)  :: ni, m1x, m1y, m2x, m2y
    real(8):: mxy, m4x, m4y
    real(8):: x, y, x2, y2,xy
    real(8) x4,y4 
    real(8)::grw(2)
    call init1_RM48(0)
    ni=0; m1x=0; m1y=0;m2x=0; m2y=0;mxy=0
    m4x=0; m4y=0;
    do while (ni<nit)
       ni=ni+1;
       call gaussian_random_vectorRM48(grw)
       x=grw(1);y=grw(2)
       x2=x*x;y2=y*y ;xy=x*y; 
       x4=x2*x2; y4=y2*y2
       m1x=m1x+x;m1y=m1y+y; m2x=m2x+x2; m2y=m2y+y2
       mxy=mxy+xy; 
       m4x=m4x+x4; m4y=m4y+y4
    end do
    m1x=m1x/nit;m1y=m1y/nit;m2x=m2x/nit;m2y=m2y/nit;
    mxy=mxy/nit ;
    m4x=m4x/nit;m4y=m4y/nit;
    print *,' m1x=',m1x,' m1y=',m1y
    print *,' m2x=',m2x,' m2y=',m2y
    print *,' m4x=',m4x,' m4y=', m4y
    print *,' mxy=',mxy
  end subroutine PolarGaussAccuracyTest1

  subroutine MarsagliaBrayAccuracyTest1( nit)
    implicit none
    real(8), intent(in) ::nit 

    real(8)  :: ni, m1x, m1y, m2x, m2y
    real(8):: mxy, m4x, m4y
    real(8):: x, y, x2, y2,xy
    real(8) x4,y4 
    real(8)::grw1,grw2
    call init1_RM48(0)
    ni=0; m1x=0; m1y=0;m2x=0; m2y=0;mxy=0
    m4x=0; m4y=0;
    do while (ni<nit)
       ni=ni+1;
       call Marsaglia_Bray_RM48(grw1)
       call Marsaglia_Bray_RM48(grw2)         
       x=grw1;y=grw2 
       x2=x*x;y2=y*y ;xy=x*y; 
       x4=x2*x2; y4=y2*y2
       m1x=m1x+x;m1y=m1y+y; m2x=m2x+x2; m2y=m2y+y2
       mxy=mxy+xy; 
       m4x=m4x+x4; m4y=m4y+y4
    end do
    m1x=m1x/nit;m1y=m1y/nit;m2x=m2x/nit;m2y=m2y/nit;
    mxy=mxy/nit ;
    m4x=m4x/nit;m4y=m4y/nit;
    print *,' m1x=',m1x,' m1y=',m1y
    print *,' m2x=',m2x,' m2y=',m2y
    print *,' m4x=',m4x,' m4y=', m4y
    print *,' mxy=',mxy
  end subroutine MarsagliaBrayAccuracyTest1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine  PolarGaussSpeedTest1( nit)

    implicit none
    real(8), intent(in) ::nit ! number of iterations, 
    real(8)  :: ni, m1x, m1y  
    real(8)::grw(2)
    call init1_RM48(0)
    ni=0; m1x=0; m1y=0; 
    do while (ni<nit)
       ni=ni+1;
       call gaussian_random_vectorRM48(grw)
       m1x=m1x+grw(1);m1y=m1y+grw(2);
    end do
    m1x=m1x/nit;m1y=m1y/nit;
    print *,' m1x=',m1x,' m1y=',m1y
  end subroutine  PolarGaussSpeedTest1

  subroutine MarsagliaBraySpeedTest1( nit)

    implicit none
    real(8), intent(in) ::nit ! number of iterations,
    real(8)  :: ni, m1x, m1y
    real(8)::grw1,grw2
    call init1_RM48(0)
    ni=0; m1x=0; m1y=0; 
    do while (ni<nit)
       ni=ni+1;
       call Marsaglia_Bray_RM48(grw1); call Marsaglia_Bray_RM48(grw2)
       m1x=m1x+grw1;m1y=m1y+grw2;
    end do
    m1x=m1x/nit;m1y=m1y/nit;
    print *,' m1x=',m1x,' m1y=',m1y
  end subroutine MarsagliaBraySpeedTest1


  !----------------------------------------
  subroutine PolarGaussAccuracyTest2( nit)

    implicit none
    real(8), intent(in) ::nit

    real(8)  :: ni, m1x, m1y, m2x, m2y
    real(8):: mxy, m4x, m4y
    real(8):: x, y, x2, y2,xy
    real(8) x4,y4
    real(8)::grw(2)
    call init_RANDOM_NUMBER(0)
    ni=0; m1x=0; m1y=0;m2x=0; m2y=0;mxy=0
    m4x=0; m4y=0;
    do while (ni<nit)
       ni=ni+1;
       call gaussian_random_vector(grw)
       x=grw(1);y=grw(2)
       x2=x*x;y2=y*y ;xy=x*y; 
       x4=x2*x2; y4=y2*y2
       m1x=m1x+x;m1y=m1y+y; m2x=m2x+x2; m2y=m2y+y2
       mxy=mxy+xy; 
       m4x=m4x+x4; m4y=m4y+y4
    end do
    m1x=m1x/nit;m1y=m1y/nit;m2x=m2x/nit;m2y=m2y/nit;
    mxy=mxy/nit ;
    m4x=m4x/nit;m4y=m4y/nit;
    print *,' m1x=',m1x,' m1y=',m1y
    print *,' m2x=',m2x,' m2y=',m2y
    print *,' m4x=',m4x,' m4y=', m4y
    print *,' mxy=',mxy
  end subroutine PolarGaussAccuracyTest2

  subroutine MarsagliaBrayAccuracyTest2(nit)

    implicit none
    real(8), intent(in) ::nit 

    real(8)  :: ni, m1x, m1y, m2x, m2y
    real(8):: mxy, m4x, m4y
    real(8):: x, y, x2, y2,xy
    real(8) x4,y4 
    real(8)::grw1,grw2
    call init_RANDOM_NUMBER(0)
    ni=0; m1x=0; m1y=0;m2x=0; m2y=0;mxy=0
    m4x=0; m4y=0;
    do while (ni<nit)
       ni=ni+1;
       call MarsagliaBray(grw1)
       call MarsagliaBray(grw2)
       x=grw1;y=grw2 
       x2=x*x;y2=y*y ;xy=x*y; 
       x4=x2*x2; y4=y2*y2
       m1x=m1x+x;m1y=m1y+y; m2x=m2x+x2; m2y=m2y+y2
       mxy=mxy+xy; 
       m4x=m4x+x4; m4y=m4y+y4
    end do
    m1x=m1x/nit;m1y=m1y/nit;m2x=m2x/nit;m2y=m2y/nit;
    mxy=mxy/nit ;
    m4x=m4x/nit;m4y=m4y/nit;
    print *,' m1x=',m1x,' m1y=',m1y
    print *,' m2x=',m2x,' m2y=',m2y
    print *,' m4x=',m4x,' m4y=', m4y
    print *,' mxy=',mxy
  end subroutine MarsagliaBrayAccuracyTest2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine  PolarGaussSpeedTest2( nit)

    implicit none
    real(8), intent(in) ::nit ! number of iterations, intent(in)  
    real(8)  :: ni, m1x, m1y  
    real(8)::grw(2)
    call init_RANDOM_NUMBER(0)
    ni=0; m1x=0; m1y=0; 
    do while (ni<nit)
       ni=ni+1;
       call gaussian_random_vector(grw)
       m1x=m1x+grw(1);m1y=m1y+grw(2);  
    end do
    m1x=m1x/nit;m1y=m1y/nit; 
    print *,' m1x=',m1x,' m1y=',m1y
  end subroutine  PolarGaussSpeedTest2

  subroutine MarsagliaBraySpeedTest2( nit)

    implicit none
    real(8), intent(in) ::nit ! number of iterations, intent(in)  
    real(8)  :: ni, m1x, m1y  
    real(8)::grw1,grw2
    call init_RANDOM_NUMBER(0)
    ni=0; m1x=0; m1y=0; 
    do while (ni<nit)
       ni=ni+1;
       call MarsagliaBray(grw1); call MarsagliaBray(grw2) 
       m1x=m1x+grw1;m1y=m1y+grw2;  
    end do
    m1x=m1x/nit;m1y=m1y/nit; 
    print *,' m1x=',m1x,' m1y=',m1y
  end subroutine MarsagliaBraySpeedTest2

  subroutine correlated_GaussTest(nit)
    implicit none
    real(8), intent(in) ::nit ! number of iterations, intent(in)
    real(8)  :: ni, m1x, m1y, m2x,m2y, mxy, m4x, m4y, delta_t
    real(8)  ::x,y, x2,y2,xy,x4,y4, delta_w, delta_Z
    real(8)::normal_vector(2)
    m1x=0; m1y=0; m2x=0; m2y=0; mxy=0; m4x=0; m4y=0;
    call init_RANDOM_NUMBER(0)
    ni=0
    do while (ni<nit)
       ni=ni+1;
       delta_t=4.d0
       call gaussian_random_vectorRM48(normal_vector)
       call Correlated_Gauss_Gen(delta_t, normal_vector, delta_w, delta_Z  )

       x=delta_w;y=delta_Z
       x2=x*x;y2=y*y ;xy=x*y; 
       x4=x2*x2; y4=y2*y2
       m1x=m1x+x;m1y=m1y+y; m2x=m2x+x2; m2y=m2y+y2
       mxy=mxy+xy; 
       m4x=m4x+x4; m4y=m4y+y4
    end do
    m1x=m1x/nit;m1y=m1y/nit;m2x=m2x/nit;m2y=m2y/nit;
    mxy=mxy/nit ;
    m4x=m4x/nit;m4y=m4y/nit;
    print *,' m1x=',m1x,' m1y=',m1y
    print *,' m2x=  ',m2x,' exact=',delta_t, ' m2y=',m2y, ' exact=',delta_t**3/(3.d0) 
    print *,' m4x=',m4x,' exact=',3.d0*delta_t**2,  ' m4y=', m4y,' exact=',delta_t**6/(3.d0)
    print *,' mxy=',mxy, ' exact=',delta_t**2/(2.d0)

    return
  end subroutine correlated_GaussTest


end module RFOF_random_numbers_SG2
