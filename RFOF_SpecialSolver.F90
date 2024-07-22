module RFOF_SpecialSolver

  use diffusion_coef
  use RFOF_random_numbers_SG2

  common /liniar_aprox/AA, BB, lambda

contains

  !------------------------------------------
  function aprox_equal(x,y) result(yn)
    implicit none
    real(8), intent(in):: x,y
    real(8):: internal_error !>The internal relative error
    real(8):: a !>internal
    logical  ::yn 
    data internal_error/1.0d-8/ ! For Monte-Carlo it is sufficient
    a= dabs(x-y)-(dabs(x)+dabs(y) )*internal_error
    if(a.LE.0) then
       yn=.TRUE.
    else
       yn=.FALSE.
    endif
    return
  end function aprox_equal


  !------------------------------------------------
  subroutine RFOF_sigma(x, sigma, Error_flag) 
    implicit none
    real(8), intent(in):: x
    real(8), intent(out)  ::sigma

    integer, intent(out):: Error_Flag
    real(8)::diff_c !< Diffusion coefficient
    !Error_Flag=0
    call diffusion_coefficient(x,diff_c,Error_flag)
    if (error_Flag.NE.0) return
    if (diff_c.LE.0) then
       Error_flag=9
       return
    endif
    sigma=sqrt(2.d0*diff_c) 
    return
  end subroutine RFOF_sigma

  !--------------------------------

  subroutine linear_fit( x, y,A,B,  ser)
    implicit none
    real(8), intent(in)::y(2),x(2)
    real(8), intent(out)::A, B 
    integer, intent(out)::ser  !< In the degenerate case, when the order is different from \f$x1<x2\f$ , then ser=8. In the normal case ser=0

    A=0;B=0   
    if(  x(1).GE.x(2) )  then
       ser=8; return
    endif
    ser=0

    A=(y(2)-y(1))/(x(2)-x(1))
    B=y(1)-A*x(1) 
    return
  end subroutine linear_fit



  !--------------------------------------------------------------------
  subroutine linear_aprox_minmax( x, y,A,B,minmaxerr, ser)
    implicit none
    real(8), intent(in)::y(3),x(3)
    real(8), intent(out)::A, B, minmaxerr
    integer, intent(out)::ser
    real(8):: norm  !internal 


    A=0;B=0;minmaxerr=1.0d10  
    if(  (x(1).GE.x(2)) .OR.(x(2).GE.x(3)))  then
       ser=8; return
    endif
    ser=0
    norm=max(dabs(y(1)-y(3)),dabs(y(1)-y(2))) 
    norm=max(norm,dabs(y(3)-y(2)))
    if(norm.EQ.0) then
       A=0; minmaxerr=0; B=y(1);
       return
    else 
       A=(y(3)-y(1))/(x(3)-x(1))
       B=0.5d0*(y(1)+y(2)-A*(x(1)+x(2)))
       minmaxerr=dabs(0.5d0*(A*(x(2)-x(1))+y(1)-y(2)))
       minmaxerr=minmaxerr/norm
    endif
    return
  end subroutine linear_aprox_minmax


  !--------------------------------------------------------------------------------
  subroutine linear_aprox_minmax2( xc,delta, y,A,B,minmaxerr, ser)
    implicit none
    real(8), intent(in)::y(3),xc, delta
    real(8), intent(out)::A, B, minmaxerr
    integer, intent(out)::ser
    real(8):: norm  !internal 

    A=0;B=0;minmaxerr=1.0d10  

    if(  delta<=0)  then

       print*, ' delta=', delta, ' xc=',xc,' y=',y 
       ser=8; return
    endif
    ser=0
    norm=max(dabs(y(1)-y(3)),dabs(y(1)-y(2))) 
    norm=max(norm,dabs(y(3)-y(2)))
    if(norm.EQ.0) then
       A=0; minmaxerr=0; B=y(1);
       return
    else 
       A=0.5d0*(y(3)-y(1))/delta
       B=0.5d0*(y(1)+y(2)-A*(2.d0*xc-delta)  )
       minmaxerr=dabs(0.5d0*(A*delta+y(1)-y(2)))
       minmaxerr=minmaxerr/norm
    endif
    return
  end subroutine linear_aprox_minmax2


  !----------------------------------------------------------------------

  subroutine quadratic_interpolation( xc,delta_x, y,a,b,c, ser)
    implicit none
    real(8), intent(in)::y(3),xc, delta_x
    real(8), intent(out)::a !< The returned coefficients of the interpolating polynomial
    real(8), intent(out)::b !< The returned coefficients of the interpolating polynomial
    real(8), intent(out)::c !< The returned coefficients of the interpolating polynomial
    integer, intent(out)::ser !> In the degenerate case, when the order is different from \f$x1<x2<x3\f$, then ser=8. In the normal case ser=0
    real(8)::  d2,dd, dd2,y13,y123 !internal 

    a=0;b=0; c=0; ser=0 
    if(  delta_x<=0)  then
       ser=8; return
    endif
    dd =delta_x*delta_x 
    dd2=dd*2.d0
    d2=2.0d0*delta_x
    y13=y(1)-y(3)
    y123=y(1)-2.d0*y(2)+y(3)
    a=y123/dd2
    b=-y13/d2-xc*y123/dd
    c=y(2)+xc*y13/d2+xc*xc*y123/dd2
    return
  end subroutine quadratic_interpolation



  !--------------------------------------------
  subroutine RFOF_exact_solution(lambda,x_in, delta_t,W,A,B, x_fin,ser)
    implicit none
    real(8), intent(in)::lambda,x_in,delta_t, W,A,B
    real(8), intent(out):: x_fin
    integer, intent(out)::ser !>Error signal: ser=7 when the exponent is too large, ser=0 when OK.
    real(8)::BA, expp
    ser=0
    if(A.EQ.0) then
       x_fin=x_in+B*W; return
    endif
    BA=B/A
    expp=delta_t*A**2*(lambda-0.5d0)+A*w
    if (dabs(expp)>3) then 
       ser=7; return
    endif
    x_fin=(x_in+BA)*exp(expp)-BA

    return
  end subroutine RFOF_exact_solution

  !------------------------------------------------------------------
  subroutine reflect(miror, object_image) 
    real(8), intent(in):: miror !Mirror position 
    real(8), intent(inout)::object_image !>

    object_image=2.d0*miror-object_image
  end subroutine reflect
  !----------------------------------------------------------



  subroutine RFOF_one_step_generalcase(method, x_min, delta_t_in,x_in, tolerance,&
       delta_t_out,x_fin, output_error,ser)
    implicit none
    integer, intent(in)::method

    real(8), intent(in)::  x_min !>Esential cut-off: allowed lower limit for x(t). When x(t)<x_min we reset

    real(8), intent(in):: delta_t_in, x_in !>initial time step, initial value of variable
    real(8), intent(in)::tolerance

    real(8), intent(out)::x_fin !> The final position, robust, error controlled result, 0.5 strong order 
    real(8), intent(out)::  delta_t_out,output_error  

    integer:: maxdivision, ncount !>maximal number of reducing the time step to attain the requested accuracy; the counter
    real(8)::linearity_approx_limit,log_lin_apprx_lim

    real(8):: delta_t_curent, nonlin_error_curent
    real(8)::A,B, W, dW, lambda

    integer, intent(out)::ser !>Error signal
    real(8):: image !>Image of reflection at the "miror" x=x_min  called at end points 

    real(8)::sigma(3), sgm 

    real(8)::sqrtdt, delta_x, xf !internal
    DATA maxdivision /10/
    COMMON /RFOF_boundary/linearity_approx_limit
    ser=0
    !print *,' ONE step start x_in=', x_in
    if(x_in.LE.0) then
       print*, ' bug in the input x_in in the subroutine RFOF_one_step_gen'
       ser=2
       return
    endif
    call RFOF_sigma(x_in, sgm, ser) 
    if(ser.NE.0) return
    sigma(2)=sgm 
    if(sigma(2)<=0) then
       x_fin=x_in
       delta_t_out=delta_t_in
       output_error=0
       ser=9
       return
    endif
    delta_t_curent=delta_t_in
    sqrtdt=sqrt(delta_t_curent)
    delta_x=sigma(2)*sqrtdt

    ! Diagnostic
    if (delta_x==0) then
       print *, ' DIAGN: delta_t_in', delta_t_in, ' sigma2=',sigma(2)
    endif

    ! end diagnostic

    !Previous version, this produce sticking
    if(x_in-delta_x<0) then
       delta_x=dabs(x_in)*0.5d0
       delta_t_curent=(delta_x/sigma(2))**2
    endif

    ncount=0 
    !DIAGNOSTIC
    !PRINT*,' '
    !PRINT*,' '
    !PRINT*,' START one step ', ' xin=',x_in, ' delta_x=', delta_x, ' delta_t=', delta_t_curent
    !PRINT*,' Start iteration time step'
    DO  
       ncount=ncount+1
       if(ncount.EQ.maxdivision) then
          ser=1
          exit
       endif
       call RFOF_sigma(x_in-delta_x, sgm, ser) 
       if(ser.NE.0) exit
       sigma(1)=sgm 
       call RFOF_sigma(x_in+delta_x, sgm, ser) 
       sigma(3)=sgm 
       if(ser.NE.0)exit
       !print *, ' TIME step ITERATION, nr=', ncount
       !print *, 'sigma =', sigma
       ! print *, 'delta_x=', delta_x

       call linear_aprox_minmax2( x_in,delta_x, sigma,A,B,nonlin_error_curent, ser)
       !print *, 'error=', nonlin_error_curent

       output_error=nonlin_error_curent
       if(ser>0) exit
       output_error=nonlin_error_curent

       if(nonlin_error_curent.LT.tolerance) then 
          exit
       endif
       delta_x=delta_x/2 
       !print *, 'delta-x new', delta_x
       delta_t_curent= delta_t_curent/4
    enddo !end of the time step decreasing loop
    ! PRINT*,' END iteration time step'
    !print *, ' delta-x final=', delta_x
    if(ser>0) return 
    lambda=1.d0
    !delta_t_curent=(delta_x/sigma(2))**2
    sqrtdt=sqrt(delta_t_curent) 

    IF(method==1)  THEN
       ! We use the strong order 0.5 method, robust 

       call Marsaglia_Bray_RM48( W )
       dW=W*sqrtdt
       Call RFOF_exact_solution(lambda,x_in, delta_t_curent,dW,A,B, xf,ser)
       !Print *, 'ONE STEP Ex SOL, x_in=', x_in, ' xf=',xf
       !if(ser>0) exit
       if(xf<x_min) then
          !reflect(miror, object_image) 
          call reflect(x_min, xf )
          !Print *, 'REFLECT', xf
       endif
       !if(ser>0) exit
    ELSE

       Call RFOF_SDE_integrator_corection(x_min, x_in, delta_x, sigma, delta_t_curent,xf,ser) 
       !if(xf>x_min) then 
       !x_fin=xf
       !else
       !x_fin=x_min
       !endif 
    ENDIF
    delta_t_out=delta_t_curent
    x_fin=xf
    !print*, ' Iteration Mthod=',method, ' XF=', xf,  ' increment=',xf-x_in
    !PRINT*,' END one step'
    ! PRINT*,' ' 
    !PRINT*,' '  
    return
  end  subroutine RFOF_one_step_generalcase


  !-----------------------------------------------------------

  subroutine RFOF_SDE_integrator_corection(x_min,x_in, delta_x, sigma, dt,  x_final, ser)
    implicit none
    real(8), intent(in)::x_min ! The minimal admisible value of x. Cutoff parameter
    real(8), intent(in)::x_in,delta_x, dt, sigma(3) 
    real(8), intent(out)::x_final !> The improved value
    real(8):: m,n,p !>Coefficients of the quadratic interpolation of the Sqrt(diffusion_coeff(x))
    real(8)::dw, dz 

    real(8)::normal_vector(2)

    real(8):: s1,s, sms, dy1,dy2,dy3,dy4 ! internal variables
    integer::ser !> Error signal,  equal 0 when OK .

    call quadratic_interpolation( x_in,delta_x, sigma,m,n,p, ser)
    if (ser>0) return 

    call gaussian_random_vectorRM48(normal_vector)  
    call Correlated_Gauss_Gen( dt ,normal_vector,dw, dz  )
    s=m*x_in**2+n*x_in+p
    s1=2.d0*m*x_in+n 
    sms=(s1**2+2.d0*m*s)
    dy1=dw*s+dt*s1*s+ 0.5d0*(-dt+dw**2)*s*s1
    dy2=dz*s*sms +&
         0.5d0*dt**2* (3*m*s1*s**2+&
         s1*s*sms )
    dy3=(dt*dw-dz)*(s1**2*s+m*s**2)
    dy3=(dt*dw-dz)* s*(s1**2+m*s)
    dy4=0.5*dw*(-dt+dw**2/3)*s*sms
    x_final=x_in+dy1+dy2+dy3+dy4

    if(x_final<x_min)  then 
       !reflect(miror, object_image) 
       call reflect(x_min, x_final )
    endif

    return
  end subroutine RFOF_SDE_integrator_corection


  !----------------------------------------------------------------------------------

  subroutine RFOF_multi_step_generalcase( method, x_min, delta_t_in, delta_t_out, x_in,x_fin, &
       time_in, time_fin, T_max, tolerance,ser)

    implicit none
    integer, intent(in)::method

    real(8), intent(in)::  x_min !>Esential cut-off: allowed lower limit for x(t). When x(t)<x_min we reset

    real(8), intent(in):: delta_t_in, x_in !>initial time step, initial value of variable
    real(8), intent(in)::tolerance

    real(8), intent(out)::   x_fin 
    real(8)::   x_final_extrapol  !> INternal variable, the extrapolated return value of one -step integrator

    real(8), intent(in):: time_in, T_max 

    real(8), intent(out):: time_fin

    real(8), intent(out):: delta_t_out !> The readjusted time step, according to the output "accuracy"
    real(8) :: output_error !> Internal, the actual error from the one-step subroutine. It



's used     ! to readjust the time step.    integer:: ncount, maxdivision !>maximal subdivision allowed for halving the time step    real(8)::linearity_approx_limit !,log_lin_apprx_lim    !>This subroutine is designed for x above 'linearity_approx_limit























'.Below this value we wrote the    !> subroutine RFOF_multi_step_boundarycase    real(8):: delta_t_curent !> the adjusted ccurent time step    real(8)::time_curent    real(8)::x_curent,x_new !>internal variables    real(8)::relative_error !>internal variable returned by one step subroutine.The estimated relative error    real(8)::new_time !>internal variable, the new time  after adjusting according to given input rell. error    integer, intent(out)::ser    !>The error signal ser=0 when OK.       !>ser=1 when by succesive reducing the time step the accuracy is sstill greater then tolerance    !>ser=2 when the input initial position is wrong    !>ser=3 when the input time step is wrong    !>ser=4 when the input initial time is wrong    !>ser=7: Too large exponent in the subroutine RFOF_Exact_Solution    !>ser=8 when the error signal from called subroutine  linear_aprox_minmax2  is nonzero    !>ser=9 Improbable value of the diffusion coefficient    !>ser=10 when the error signal from called subroutine  for diffusion coefficient  is nonzero    DATA maxdivision /10/ !> The time step delta_t si decreased at most by a factor 2**maxdivision    ser=0 !>Normal return    !> This maybe desactivated and replaced by automatic correctin to cut-off x_min    if(x_in.LT.x_min) then       print*, ' bug in the input y_in in the subroutine RFOF_multi_step_boundarycase




'        ser=2       return    endif    if(delta_t_in.LE.0) then       print*, ' bug in the input delta_t in the subroutine RFOF_multi_step_boundarycase





'       ser=3       return    endif    if(time_in.LT.0) then       print*, ' bug in the input time_in the subroutine RFOF_multi_step_boundarycase




























'       ser=4       return    endif    !!ATTENTION, automatic correction, error catching eliminated    ! if (x_in<x_min) then    !x_curent=x_min     !else    !x_curent=x_in    ! endif    !!___________________    time_curent=time_in    if(aprox_equal(time_curent,T_max)) then       time_fin=T_max       x_fin=x_in       delta_t_out=delta_t_in       return    endif    delta_t_curent=delta_t_in    x_curent=x_in    !> Main time loop    do       !>Adjusting the time step to reach exactly the final value T_max, when time_curent is close to T_max       new_time=time_curent+delta_t_curent       if(new_time.GT.T_max) then          delta_t_curent=T_max-time_curent            ! DIAGNOSTIC          !  if(delta_t_curent.LT.(1.0d-5)) then      !print *, ' !!MULTISTEP!! time=', time_curent,' T_max=',T_max, 'delta_t=', delta_t_curent
          ! ser=100
          ! exit
          ! endif
          ! END DIAGNOSTIC
       endif

       !RFOF_one_step_generalcase(method, x_min, delta_t_in,   x_in, tolerance,&
       !delta_t_out,x_fin, output_error,ser)
       call RFOF_one_step_generalcase(method,x_min, delta_t_curent,x_curent, tolerance,delta_t_out,x_new,&
            output_error,ser)
       ! print *, ' '
       ! print *, ' '
       !  Print *, 'multistep x_curent=',x_curent, ' x_new=', x_new, ' time=',time_curent

       if(ser>0) exit

       delta_t_curent=delta_t_out
       ! DIAGNOSTIC
       !if(delta_t_curent.LT.(1.0d-5)) then 
       !print *, ' !!MULTISTEP, after call one step!!time=', time_curent,  'delta_t_curent=', delta_t_curent
       ! ser=100
       ! exit
       !endif
       ! END DIAGNOSTIC
       time_curent=time_curent+delta_t_curent
       if(x_new<x_min) then

          !print *,' before reflect in multi-step, xnew=',x_new
          call reflect(x_min,x_new)
       endif
       x_curent=x_new
       if(aprox_equal(time_curent,T_max)) then
          time_Curent=T_max
          exit
       endif

       if( output_error.LT.tolerance*1.0d-1)  then
          delta_t_curent=delta_t_curent*2.0d0
       endif
       !print *, ' end multistep time loop'
    enddo

    if(ser.NE.0) then
       print *, 'Anomalous return'
       if(ser==1) print *,' Too much accuracy required'
       return
    endif
    delta_t_out=delta_t_curent 
    x_fin=x_curent
    time_fin=time_curent
    ! DIAGNOSTIC
    !print *, ' end multistep x_fin=',x_fin
    !print *, ' '
    !print *, ' '
    ! END DIAGNOSTIC
    return
  end subroutine RFOF_multi_step_generalcase

  !----------------------------------------------------
  subroutine RFOF_boundary_ymin(d_prime,y_min )
    implicit none
    real(8), intent(in)::d_prime

    real(8), intent(out)::y_min !> minimal admisible value of y 
    real(8)::linearity_approx_limit
    COMMON /RFOF_boundary/linearity_approx_limit
    y_min=log(d_prime**2/linearity_approx_limit)/2
    return
  end subroutine RFOF_boundary_ymin



  !----------------------------------------------



  subroutine RFOF_multi_step_boundarycase(y_min,y_max,delta_t_in,y_in,time_in,T_max, &
       tolerance,delta_t_out,y_fin, time_fin,ser)
    implicit none
    real(8), intent(in):: y_min

    real(8), intent(in):: y_max

    real(8), intent(in):: delta_t_in, y_in


    real(8), intent(in):: time_in, T_max

    real(8), intent(out):: time_fin

    real(8), intent(in)::tolerance

    real(8), intent(out)::  delta_t_out, y_fin  

    integer:: ncount, maxdivision !>maximal subdivision allowed for halving the time step
    real(8)::linearity_approx_limit !,log_lin_apprx_lim

    real(8):: delta_t_curent !> the adjusted ccurent time step
    real(8)::normal_vector(2)!> internal, N(0,1) independent random nmbers, denoted by U(2)

    real(8)::time_curent
    real(8)::y_new !>internal variable
    ! real(8)::y_min !>
    real(8)::relative_error !>internal variable returned by one step subroutine.The estimated relative error
    real(8)::new_time !>internal variable, the new time  after adjusting according to given input rell. error
    integer, intent(out)::ser


    real(8)::sqrtdt, delta_y, yf, y_curent !internal
    logical :: endof_main_loop
    DATA maxdivision /10/


    ser=0
    if(y_in.LT.y_min) then
       print*, ' bug in the input y_in in the subroutine RFOF_multi_step_boundarycase'
       ser=2
       return
    endif

    if(delta_t_in.LE.0) then
       print*, ' bug in the input delta_t in the subroutine RFOF_multi_step_boundarycase'
       ser=3
       return
    endif

    if(time_in.LT.0) then
       print*, ' bug in the input time_in the subroutine RFOF_multi_step_boundarycase'
       ser=4
       return
    endif
    time_curent=time_in
    y_curent=y_in
    delta_t_curent=delta_t_in

    do
       new_time=time_curent+delta_t_curent
       if(new_time.GT.T_max) then
          delta_t_curent=T_max-time_curent   
       endif

       call  Gaussian_random_vectorRM48(normal_vector)
       ncount=0

       do

          call RFOF_one_step_boundarycase(y_max, delta_t_curent,normal_vector,&
               y_curent, y_new, relative_error) 
          if(relative_error.LT.tolerance) exit
          delta_t_curent=delta_t_curent/2

          ncount=ncount+1
          if(ncount.GT.maxdivision) then
             ser=1;
             print *, ' TOO MUCH ACCURACY REQUESTED ser=', ser
             exit
          endif
       enddo
       if(ser>0) exit

       time_curent=time_curent+delta_t_curent
       y_curent=y_new

       if(y_curent.LT.y_min ) then
          time_fin=time_curent
          exit
       endif

       if(aprox_equal(time_curent,T_max)) then
          time_curent=T_max
          exit
       endif


       if(relative_error.LT.tolerance*1.d-2)  then
          delta_t_curent=delta_t_curent*2.d0
          ! print *, ' DOUBLING, time step=', delta_t_curent
       endif
    enddo

    if(ser>0) then
       print *, 'Anomalous return'
       return
    endif
    delta_t_out=delta_t_curent   
    y_fin=y_curent
    time_fin=time_curent
    return
  end subroutine RFOF_multi_step_boundarycase

  !-------------------------------------------------------------------


  subroutine RFOF_one_step_boundarycase(y_max,delta_t,normal_vector, y_in, y_fin,  error)

    implicit none

    real(8), intent(in):: y_max

    real(8), intent(in)::normal_vector(2) !>N(0,1) independent random nmbers, denoted by U(2)
    real(8), intent(in):: delta_t, y_in !> time step, initial value
    real(8), intent(out)::   y_fin,  error  !> final value, estimated  error (absolute error for y, relative error for x)
    real(8)::y_new !>Internal, output to y_fin when y_new<y_max
    real(8):: delta_w, delta_z !> internal correlated Gaussian random variables
    real(8)::w2, de, de2, de3,  dy1,dy3, dy4, mhalf, third, one_8,mone_16,dw2
    data mhalf /-0.5d0/, third /0.33333333333333333d00/, one_8/0.125d0/, mone_16 /-0.0625d0/  
    ! print *,' subr one step y_in=', y_in 
    de= exp(y_in) 
    de2=de*de; de3=de2*de  !;de4=de3*de
    !call gaussian_random_vectorRM48(normal_vector) previously called
    call Correlated_Gauss_Gen( delta_t,normal_vector, delta_w, delta_z  )
    dw2=delta_w*delta_w

    dy1=mhalf*de*delta_w  + one_8* (-delta_t+dw2)*de2 
    dy3=mone_16*(delta_t*delta_w-delta_z)*de3
    dy4= -one_8*delta_w*( -delta_t+third*dw2 )*de3
    y_new=y_in+dy1+dy3+dy4
    error= dabs(dy3+dy4)

    if( y_new<y_max) then
       y_fin=y_new
    else
       y_fin=y_max   
    endif
    return
  end subroutine RFOF_one_step_boundarycase

  !---------------------------------------------------------------------------------
  ! END OF WORKING SUBROUTIONS  

  !-----------------------------------------------------------------------------------------
  ! Test program for the subroutines  RFOF_multi_step_boundarycase and subroutine RFOF_one_step_boundarycase 
  !Author:Gyorgy Steinbrecher, EURATOM-MEdC

  !------------------------------------------------------------------------------
  ! MAIN :CALL OF THE TEST PROGRAM FOR BOUNDARY APPROXIMATIONS
  !program specialSolvertest
  !use RFOF_SpecialSolver 
  !implicit none
  !intrinsic CPU_TIME 
  !real(8):: ntraject !> Number of MC trajectory sampling
  !ntraject=1.0d4
  !call multi_step_boundarytest1(ntraject)
  !end program specialSolvertest

  !END OF THE MAIN :CALL OF THE TEST PROGRAM 
  !------------------------------------------------------------------

  !------------------------------------------------
  function exact_moments(x) result(y)
    real(8), intent(in)::x
    real(8)::y
    y=(2.d0+x)**(-2)
    return
  end function exact_moments
  !---------------------------------------------------------------------------


  subroutine multi_step_boundarytest1(ntraject)

    use RFOF_random_numbers
    use RFOF_random_numbers_SG2
    implicit none
    intrinsic CPU_TIME 
    real(8), intent(in)::ntraject !>Number of trajectories in MC sampling
    real(8) ::count_traject !>counter of trajectories in MC sampling
    real(8) :: y_min

    real(8):: y_max

    real(8) :: delta_t_in, y_in

    real(8) :: time_in, T_max

    real(8):: moments(4) !> first four exp- moments of the random variable \f$y(T_{max})\f$, computed by new method
    real(8):: moments2(4) !>first four exp- moments of the random variable \f$y(T_{max})\f$, computed by Euler method

    real(8) :: time_fin

    real(8) ::tolerance

    real(8) ::  delta_t_out, y_fin  

    integer ::ser
    !integer::error_signal_linaprox  

    integer:: k !> internal
    real(8)::exactval !> internal, the exact mean value in stationary state
    real(8), parameter:: alpha=0.1d0 !> internal
    real(8)::exact_moments
    logical::escape
    real(8)::ndivision ! The number of time discretization points in the Euler method, used for test
    real(8)::time_begin,time_end, cputime_new, cputime_Euler
    ! Fix the input values 
    y_min=0  
    y_in=1.d0 
    y_max=2.d0
    delta_t_in=0.01d0 ! initial time step for new method
    time_in=0
    T_max=0.1d0
    !print * , ' T_max=', T_max
    tolerance=1.0d-3 
    ndivision=400
    ! Begin test 
    ! Initialise the seed of random number generators
    call init_RANDOM_NUMBER( 0) 
    call init1_RM48(0) 

    do k = 1, 4
       moments(k)=0
       moments2(k)=0
    enddo
    CALL CPU_TIME ( time_begin )  
    count_traject=0
    ! Start of Monte-Carlo sampling 
    do ! Main loop
       count_traject=count_traject+1.d0 
       call RFOF_multi_step_boundarycase(y_min,y_max,delta_t_in,y_in,time_in,T_max, &
            tolerance,delta_t_out,y_fin, time_fin,ser) 
       if(ser.NE.0) exit ! particle killed
       if((y_fin.LT.y_min) ) goto 100
       ! Compute the moments
       do  k=1,4
          moments(k)=moments(k)+exp(-k*alpha*y_fin) 
       enddo
100    continue
       if (count_traject.GE.ntraject)  exit
    enddo ! End of main loop 
    if(ser.NE.0) then
       ! Anomalous ending
       Print *, 'ANOMALOUS RETURN of subroutine RFOF_multi_step_boundarycase '
       print *,' ERROR SIGNAL=',ser
       return
    endif
    do k=1,4
       moments(k)=moments(k)/ntraject 
    enddo

    Print *, 'delta_t_out=',delta_t_out,' ser=',ser
    CALL CPU_TIME ( time_end)
    cputime_new=time_end-time_begin  
    Print *,' Moments computed by new methods' 
    do k=1,4
       !exactval=exact_moments(k*alpha
       print *, 'moments(', k, ')=',moments(k) 
    enddo
    ! The simplest EULER method
    do k=1,4
       moments2(k)=0 ! Will be computed by a Euler method 
    enddo
    count_traject=0
    CALL CPU_TIME ( time_begin )  
    do ! Main loop
       escape=.false.
       count_traject=count_traject+1.d0 
       call boundarytest1_Euler(y_min, y_max, y_in, y_fin, time_in, T_max, ndivision, escape) 
       if(escape) goto 200 
       if((y_fin.LT.y_min) ) goto 200
       ! Compute the moments
       do  k=1,4
          moments2(k)=moments2(k)+exp(-k*alpha*y_fin)
       enddo

200    continue
       if (count_traject.GE.ntraject)  exit
    enddo

    CALL CPU_TIME ( time_end)
    cputime_Euler=time_end-time_begin  
    do k=1,4
       moments2(k)=moments2(k)/ntraject 
    enddo
    Print *,' Moments computed by both methods'
    do k=1,4
       !exactval=exact_moments(k*alpha)
       exactval=moments2(k)
       print *, 'moments(', k, ')=',moments(k),' exact=',exactval 
    enddo
    print *,'cputime new=',cputime_new,'cputime_Euler=',cputime_Euler
  end subroutine multi_step_boundarytest1

  !---------------------------------------------------------------------------------------
  subroutine boundarytest1_Euler(y_min,y_max,y_in, y_fin, t_in, t_fin, ndivision, escape)
    use RFOF_random_numbers_SG2
    implicit none
    real(8), intent(in)::y_min !< The lowest limit of the allowed domain
    !< If \f$ y(t)<y_min \f$ the subroutine stops
    real(8), intent(in):: y_max
    !<If the updated value has y(t)>y_max, it is reset to \f$ y(t)=y_max \f$, in order to avoid overflow
    real(8), intent(in)::y_in !< initial position
    real(8), intent(in)::t_in !< start time
    real(8), intent(in)::t_fin !< final time
    real(8), intent(in):: ndivision !< Number of discretisation of the time interval
    real(8), intent(out):: y_fin !< final position
    logical, intent(out)::  escape !< True if \f$ y(t)<y_{min} \f$ for some \f$ t_{min}<t<t_{fin} \f$
    real(8):: dt, sqrtdt, dy, y_curent,y_curent_new,n !< internal
    real(8):: w !< \f$ N(0,1) \f$ variable

    dt=(t_fin-t_in)/ndivision; sqrtdt=sqrt(dt) 
    y_curent=y_in
    escape=.FALSE. 
    n=0
    do ! Main loop 
       if(y_curent.LE.y_min) then
          escape=.TRUE.
          exit
       endif
       call Marsaglia_Bray_RM48( w ) ! w is N(0,1) random number
       dy=0.5d0*exp(y_curent)*w*sqrtdt  ! The increment
       y_curent_new=y_curent+dy
       if(y_curent_new<y_max) then
          y_curent=y_curent_new
       else
          y_curent=y_max
       endif
       n=n+1
       if (n.GE.ndivision) exit
    enddo
    y_fin=y_curent
    return
  end subroutine boundarytest1_Euler
  ! END OF TEST OF BOUNDARY EFFECTS

  !---------------------------------------------------------------------------------------
  !TEST PROGRAMS FOR GENERAL POSITION

  !--------------------------------------------------------------------------------------
  !MAIN CALL PROGRAM 

  !-----------------------------------------------------------------------------------------
  ! Test program for the subroutines  RFOF_multi_step_generalcase and subroutine RFOF_one_step_generalcase 
  !Author:Gyorgy Steinbrecher, EURATOM-MEdC 

  !!> Test of the step solution of the SDE
  !!> $dx(t)=\sigma(x)* d\sigma(x)/dx* dt=\sigma(x)*dw$
  !!>where $\sigma(x)^{2}/2=D(x)=Diffusion_Coefficient(x)$
  !program specialSolvertest3
  !use RFOF_SpecialSolver 
  !implicit none
  !intrinsic CPU_TIME 
  !real(8):: ntraject !> Number of MC trajectory sampling
  !ntraject=100.0d0
  !call multi_step_generaltest3(ntraject)
  !end program specialSolvertest3


  !---------------------------------------------------------------------------
  ! TEST PROGRAMS 
  !---------------------------------------------------------------------------

  !------------------------------------------------
  subroutine multi_step_generaltest3( ntraject)
    use RFOF_random_numbers
    use RFOF_random_numbers_SG2
    implicit none
    intrinsic CPU_TIME 
    real(8), intent(in)::ntraject !< Number of trajectories in MC sampling
    integer ::method !<The method used in one step integration method=1 => use of 0.5 strong order method
    !<method=2 => use of 1.5 strong order method 
    real(8) ::count_traject !< counter of trajectories in MC sampling

    real(8)::  x_min !< Esential cut-off: allowed lower limit for \f$ x(t) \f$. When \f$ x(t)<x_{min} \f$ we reset

    real(8) :: delta_t_in, x_in !> initial time step, initial value of variable
    real(8) ::tolerance

    real(8) ::   x_fin 

    real(8) :: time_in, T_max 

    real(8) :: time_fin


    !integer::error_signal_linaprox
    real(8) :: delta_t_out !> The readjusted time step, according to the output "accuracy"
    real(8) :: output_error !> Internal, the actual error from the one-step subroutine. 

    integer:: ncount, maxdivision !> maximal subdivision allowed for halving the time step
    real(8)::linearity_approx_limit !> log_lin_apprx_lim

    real(8):: delta_t_curent !> the adjusted ccurent time step
    real(8)::time_curent
    real(8)::x_curent,x_new,mean_moments2(3,4) !> internal variables
    real(8)::relative_error !>internal variable returned by one step subroutine.The estimated relative error
    real(8)::new_time !>internal variable, the new time  after adjusting according to given input rell. error
    integer::ser


    real(8):: moments(100,3,4) !> first four exp- moments of the random variable \f$ y(T_max) \f$, computed by 2 of new methods

    real(8):: mean_moments(3,4) !> Mean values of previous moments(10,3,4)
    real(8)::dispersion(3,4) !>  Mean square error in previous moments(10,3,4)
    integer:: k,m,batch !> internal
    real(8)::exactval !>internal, the exact mean value in stationary state
    real(8), parameter:: alpha=0.1d0 !>internal
    real(8)::exact_moments
    real(8)::ndivision !> The number of time discretization points in the Euler method, used for test
    real(8)::time_begin,time_end, cputime(3) 
    integer:: n_batch !>Number of batch in Monte-Carlo  error estimation
    real(8):: disp !internal

    ! Fix the input values 
    n_batch=10
    x_min=  2.d0 !Lower cut-off 
    x_in=   3.d0 
    delta_t_in=1.0d-2 ! initial time step for new method
    time_in=0
    !T_max=0.05d0
    T_max= 5.0d-1
    tolerance=1.0d-2 
    ndivision=10000
    ! Begin test 
    ! Initialise the seed of random number generators
    call init_RANDOM_NUMBER(0) 
    call init1_RM48(0) 
    DO 100, method=1,3
       CALL CPU_TIME ( time_begin )  
       print*, 'START METHOD =',method 
       DO 90 batch=1,n_batch
          print*, 'START METHOD =',method , ' batch=', batch
          do k = 1, 4
             moments(batch, method,k)=0;  
          enddo
          count_traject=0
          ! Start of Monte-Carlo sampling 
          !if (method.NE.2) exit
          DO ! Main loop, traject
             count_traject=count_traject+1.d0 
             !Print*, ''
             !Print*, ''
             !Print*, 'NEW TRAJECT=',count_traject, ' batch=',batch, 'method=',method
             !RFOF_multi_step_generalcase( method, x_min, delta_t_in, delta_t_out, x_in,x_fin, &
             !time_in, time_fin, T_max, tolerance,ser)

             If(method<3) then
                call RFOF_multi_step_generalcase(method, x_min, delta_t_in, delta_t_out, x_in,x_fin,&
                     time_in, time_fin, T_max, tolerance,ser) 
             else
                call generaltest3_Euler(x_min,x_in, x_fin, time_in, T_max, ndivision,ser) 
             endif
             if(ser.NE.0) then
                print *, ' bug in method nr',method, 'ser=',ser
                exit ! Signal of bug, stop Mone-Carlo
             endif

             do  k=1,4
                moments(batch,method,k)=moments(batch,method,k)+  exp(-k*alpha*(x_fin-x_min) ) 
             enddo

             !print *, 'xf=',x_fin, 'ncount=',count_traject
             if (count_traject.GE.ntraject)  exit

          enddo ! End traject conting loop 
          if (ser.NE.0) exit
          Print *, 'method=',method,'batch=',batch,  'ser=',ser
90        CONTINUE ! End batch loop
          Print *, 'method=',method,'batch=',batch,  'ser=',ser
          CALL CPU_TIME ( time_end)
          cputime(method )=time_end-time_begin  
          if (ser.NE.0) exit
100       CONTINUE !End of main loop method
          if(ser.NE.0) then
             ! Anomalous ending
             Print *, 'ANOMALOUS RETURN of subroutine RFOF_multi_step_boundarycase '
             print *,' ERROR SIGNAL=',ser
             return
          endif
          do batch=1,n_batch
             do method=1,3 
                do k=1,4
                   moments(batch, method,k)=moments(batch,method,k)/ntraject 
                enddo
             enddo
          enddo

          do method=1,3 
             do k=1,4
                mean_moments(method,k)=0
                mean_moments2(method,k)=0
                do batch=1,n_batch
                   mean_moments(method,k)=mean_moments(method,k)+moments(batch, method,k)
                   mean_moments2(method,k)=mean_moments2(method,k)+moments(batch, method,k)**2
                enddo
                mean_moments(method,k)=mean_moments(method,k)/n_batch 
                mean_moments2(method,k)=mean_moments2(method,k)/n_batch 
                disp=mean_moments2(method,k)-mean_moments(method,k)**2
                if (disp>=0) then
                   dispersion(method,k)= sqrt(disp  )
                else
                   print*, ' ERROR in specialSolvrtest3 disp=',disp, 'method=',method,' order=',k
                endif
             enddo; enddo
             Print *,' Moments computed by new methods'
             Print *,' Moments computed by 3 methods'
110          format('method',  I2, 'moments(',I2 , ')=',D11.4,' sigma=',D11.4 )
             do k=1,4
                do method=1,3
                   write(*,110) method, k, mean_moments(method,k), dispersion(method,k)
                enddo
                print *,''
                print *,''   
             enddo
             Print *,'time1=',cputime(1),'time2=',cputime(2), 'time Euler=',cputime(3)
           end subroutine multi_step_generaltest3


           !---------------------------------------------------------------------------------------

           subroutine generaltest3_Euler(x_min,x_in, x_fin, t_in, t_fin, ndivision, ser)
             use RFOF_random_numbers_SG2
             use diffusion_coef
             implicit none
             real(8), intent(in)::x_min !> The lowest limit of the allowed domain

             real(8), intent(in)::x_in !>initial position
             real(8), intent(in)::t_in !> start time
             real(8), intent(in)::t_fin !> final time
             real(8), intent(in):: ndivision !> Number of discretisation of the time interval
             real(8), intent(out):: x_fin !> final position 
             real(8):: dt, sqrtdt, dx, x_curent,x_curent_new,n !>internal
             real(8):: w !> \f$ N(0,1) \f$ variable
             real(8)::diff_c, D_diff_c !>The diffusion coefficient and its  derivative. 

             integer::ser !error flag=0 in normal return

             dt=(t_fin-t_in)/ndivision; sqrtdt=sqrt(dt) 
             x_curent=x_in
             n=0
             do ! Main loop, n is the counter 
                if(x_curent<=x_min) then
                   x_curent=x_min
                endif
                call Marsaglia_Bray_RM48( w ) ! w is N(0,1) random number
                call diffusion_coefficient(x_curent,diff_c,ser) 
                if(ser.NE.0) exit
                call D_diffusion_coefficient(x_curent,D_diff_c)
                dx=D_diff_c*dt+sqrt(2.d0*diff_c)*sqrtdt*w ! The increment
                x_curent_new=x_curent+dx
                if(x_curent_new>=x_min) then
                   x_curent=x_curent_new
                else
                   x_curent=x_min
                endif
                n=n+1
                if (n>=ndivision) exit
             enddo
             x_fin=x_curent
             return
           end subroutine generaltest3_Euler
           ! END TEST PROGRAMES

         end module RFOF_SpecialSolver
