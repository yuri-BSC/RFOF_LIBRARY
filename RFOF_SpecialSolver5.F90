!--------------------------------------------------------------------------------------------
module RFOF_SpecialSolver5

use diffusion_coef5 ! Diff_coef(x)=O(x) for x->0, the realistic model
!use diffusion_coef4b ! Diff_coef(x)=O(x**2) for x->0, exact solution exists
use RFOF_random_numbers_SG2
use RFOF_random_numbers
! use RFOF_numerics  .Some problems, but the 3 point subroutine is OK and is
! copied at the end of this file

real(8), parameter::mean_diff_coef=1.d0

real(8), parameter::small_adimensional_diff_coef_fact=1.0d-3
!Recomended value::Sqrt(Monte_Carlo_Relative_Tolerance)
real(8), parameter:: reg_param_diff_coef=mean_diff_coef*small_adimensional_diff_coef_fact
! regularization parameter: 0<reg_param_diff_coef<<typical order of magnitude of the diffusion coefficient  
real(8), parameter:: reg_param_diff_coef2=reg_param_diff_coef*reg_param_diff_coef
real(8), parameter:: X_MIN=1.0d-3

real(8), parameter::X_MAX=1.0d4 * X_MIN

!  The alowed domain of the trajectory x(t) is: X_MAX>x(t)> X_MIN ; When the new value x(t+delta_t)
! is outside, then a) Either discard setting x(t+delta_t)=x(t); b) either use the reflected value: 
! x(t+delta_t)=2*X_MAX-x(t+delta_t) when x(t+delta_t)>X_MAX, similarly when x(t+delta_t)<X_MIN
! The reflection is described in the subroutine  "correction_domain_boundary "
 
contains

!------------------------------------------


 subroutine  correction_domain_boundary(X_inout )  
 implicit none
 real(8), intent(inout):: X_inout  
 if (X_inout<X_MAX) then
    if(X_inout>X_MIN ) then
      return
 else
   X_inout=2.d0*X_MIN-X_inout
   return
 endif 
 else      
    X_inout=2.d0*X_MAX-X_inout   
 endif
 return
 end subroutine correction_domain_boundary

!----------------------------------------------------

function aprox_equal(x,y) result(yes_no)
implicit none
 real(8), intent(in):: x,y
 real(8), parameter:: internal_relative_error=1.0d-6
  ! For Monte-Carlo this is more than sufficient  
 real(8):: a !internal
 logical  ::yes_no ! 
a= dabs(x-y)-(dabs(x)+dabs(y) )*internal_relative_error
 if(a.LE.0) then
   yes_no=.TRUE.
   else
  yes_no=.FALSE.
   endif
return
end function aprox_equal 

!----------------------------------------------------------------------------------------


subroutine RFOF_regularized_diffusion_coef(diff_coef, regularization_method,&
                               regularized_diff_coef ) 
implicit none 
real(8), intent(in):: diff_coef !the exact diffusion coef
integer, intent(in)::regularization_method
real(8), intent(out)::regularized_diff_coef
real(8)::diff2,sqrtdiff 

  if ((regularization_method==0).OR.(regularization_method==1)) then
      regularized_diff_coef=diff_coef
  elseif (regularization_method==2) then
      diff2=diff_coef*diff_coef
      regularized_diff_coef=sqrt(diff2+reg_param_diff_coef2)
  else
       print*, 'BUG IN subroutine RFOF_regularized_diffusion_coef'
  endif

  
return
end subroutine RFOF_regularized_diffusion_coef



subroutine RFOF_regularized_diffusion_coef_Euler(x,  regularization_method,&
                               regularized_diff_coef, derivative_reg_diff_coef,Error_Flag ) 
implicit none 
real(8), intent(in)::x ! magnetic momentum
integer, intent(in)::regularization_method
real(8), intent(out)::regularized_diff_coef
real(8), intent(out)::derivative_reg_diff_coef
integer, intent(out)::Error_Flag
real(8):: diff_coef !the exact diffusion coef
real(8)::diff2,sqrtdiff,derivative_diff_coef

call diffusion_coefficient5(x,diff_coef,Error_Flag)
if(Error_Flag.NE.0)  then
print*, 'BUG in subroutine RFOF_regularized_diffusion_coef call diff_coef'
return
endif
Call  D_diffusion_coefficient5(x,derivative_diff_coef)
if( (regularization_method==0).OR.(regularization_method==1)) then
  regularized_diff_coef=diff_coef 
derivative_reg_diff_coef= derivative_diff_coef
  return
  endif

 diff2=diff_coef*diff_coef
 sqrtdiff=sqrt(diff2+reg_param_diff_coef2)
 
 if (regularization_method==2) then
  regularized_diff_coef=sqrtdiff
  derivative_reg_diff_coef=diff_coef/sqrtdiff 
  ! Derivative of the regularization function, method 2
       else
       print*, 'BUG IN subroutine RFOF_regularized_diffusion_coef'
   endif
!and multiply with the derivative of diffusion coefficient
derivative_reg_diff_coef=derivative_reg_diff_coef*derivative_diff_coef
return
end subroutine RFOF_regularized_diffusion_coef_Euler


!------------------------------------------------
subroutine RFOF_sigma(regularization_method,x, sigma, Error_flag) 
 implicit none
 integer, intent(in)::regularization_method
 real(8), intent(in):: x
 real(8), intent(out)  ::sigma
 integer, intent(out):: Error_Flag ! denoted also by "Error_Flag" in the folowing
 !Error_Flag=0 normal case. Error_Flag=9 for wrong value(negative) of the diffusion coefficinet
 !Error_Flag=10 for anomalous return from the diffusion coefficient routine
 !Error_Flag=14 when the integer "regularization_method" from subroutine RFOF_sigma is not 1 or2
 real(8)::diff_coeficient !Diffusion coefficient
 real(8):: regularized_diff_coef
Error_Flag=0
if (x.LE.0.d0) then
    if (regularization_method==2) then
    sigma=sqrt(2.d0*reg_param_diff_coef)
    return
    endif
endif
    

call diffusion_coefficient5(x,diff_coeficient,Error_flag)
if (error_Flag.NE.0) return
! For diagnostic, at test, later can be removed
if (diff_coeficient.LE.0) then
Error_flag=9
 return
endif 
! For diagnostic, at test, later can be removed
if ( (regularization_method<0).OR.(regularization_method> 2) ) then
Error_flag=14
 return
endif 
! End removable diagnostic
call RFOF_regularized_diffusion_coef(diff_coeficient, regularization_method,&
                               regularized_diff_coef) 
 
!sigma=sqrt(2.d0*diff_c) is the exact form, now it is regularized
sigma=sqrt(2.d0*regularized_diff_coef) 
return
 end subroutine RFOF_sigma 

!--------------------------------
 
subroutine linear_fit( x, y,A,B,  Error_Flag)
implicit none
real(8), intent(in)::y(2),x(2)
real(8), intent(out)::A, B 
integer, intent(out)::Error_Flag   
 
! In the degenerate case, when the order is different from x1<x2 , then Error_Flag=8 
!In the normal case Error_Flag=0

A=0;B=0   
if(  x(1).GE.x(2) )  then
Error_Flag=8; return
endif
Error_Flag=0

A=(y(2)-y(1))/(x(2)-x(1))
B=y(1)-A*x(1) 
 return
end subroutine linear_fit 


!--------------------------------------------------------------------
subroutine linear_aprox_minmax( x, y,A,B,minmaxerr, Error_Flag)
  implicit none
  real(8), intent(in)::y(3),x(3)
  real(8), intent(out)::A, B, minmaxerr
  integer, intent(out)::Error_Flag
  real(8):: norm  !internal 
  ! In the degenerate case, when the order is different from x1<x2<x3, then Error_Flag=8 
  ! In the normal case Error_Flag=0

  A=0;B=0;minmaxerr=1.0d10  
  if(  (x(1).GE.x(2)) .OR.(x(2).GE.x(3)))  then
     Error_Flag=8; return
  endif
  Error_Flag=0
  norm=dmax1(dabs(y(1)),dabs( y(2)),dabs( y(3))) 
  if(norm.EQ.0) then
     A=0; minmaxerr=0; B=0;
     Error_Flag=12
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

subroutine linear_aprox_minmax2( xc,delta, y,A,B,minmaxerr, Error_Flag)
  implicit none
  real(8), intent(in)::y(3),xc, delta
  real(8), intent(out)::A, B, minmaxerr
  integer, intent(out)::Error_Flag
  real(8):: norm  !internal 
  ! In the degenerate case, when the order is different from x1<x2<x3, then Error_Flag=8 
  !In the normal case Error_Flag=0
  A=0;B=0;minmaxerr=1.0d10  
  if(  delta<=0)  then
     print*, ' delta=', delta, ' xc=',xc,' y=',y 
     Error_Flag=8; return
  endif
  Error_Flag=0
  norm=dmax1(dabs(y(1) ),dabs(y(2) ),dabs(y(3) )  ) 
  if(norm.EQ.0) then
     A=0; minmaxerr=0; B=0;
     Error_Flag=12
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

subroutine quadratic_interpolation( xc,delta_x, y,a,b,c, Error_Flag)
  implicit none
  real(8), intent(in)::y(3),xc, delta_x
  real(8), intent(out)::a,b,c !>The returned coeficients of interpolating polinomial
  integer, intent(out)::Error_Flag
  ! In the degenerate case, when the order is different from x1<x2<x3, then Error_Flag=8 
  !In the normal case Error_Flag=0
  real(8)::  d2,dd, dd2,y13,y123 !internal 
  ! Try to detect anomaly
  a=0;b=0; c=0; Error_Flag=0 
  if(  delta_x<=0)  then
     Error_Flag=8; return
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
subroutine RFOF_exact_solution(lambda,x_in, delta_t,W,A,B, x_fin,Error_Flag)
  implicit none
  real(8), intent(in)::lambda,x_in,delta_t, W,A,B
  real(8), intent(out):: x_fin
  integer, intent(out)::Error_Flag !Error signal: Error_Flag=7 when the exponent is too large, Error_Flag=0 when OK.
  real(8)::BA, expp
  Error_Flag=0
  if(A.EQ.0) then
     x_fin=x_in+B*W; return
  endif
  BA=B/A
  expp=delta_t*A**2*(lambda-0.5d0)+A*w
  if (dabs(expp)>1.0d0) then ! The exponent must be infinitesimal
     print* ,' Exponent in exact solution too large=', expp
     print*,'Increase one of global parameters: X_MIN or Small_adimensional_diff_coef_fact'
     print*, 'Or decrease the tolerance'
     Error_Flag=7; return
  endif

  x_fin=(x_in+BA)*dexp(expp)-BA

  return
end subroutine RFOF_exact_solution

!------------------------------------------------------------------
subroutine reflect(mirror, object_image) 
  implicit none
  real(8), intent(in):: mirror !Mirror position 
  real(8), intent(inout)::object_image 
  !Object/image position
  object_image=2.d0*mirror-object_image
end subroutine reflect

!----------------------------------------------------------
 
  subroutine RFOF_one_step_generalcase(integration_method,regularization_method,delta_t_in,&
  x_in, tolerance,delta_t_out,x_fin, output_error,nr_subdivision ,Error_Flag)
  implicit none
   integer, intent(in)::integration_method
   ! When method=1 the robust, strong order 0.5 method is used, recomended near extreme values of x
   ! When method=2 , then the strong order 1.5 method is used.
  integer, intent(in)::regularization_method
  !Specifies the call of the subroutine RFOF_sigma()
  real(8), intent(in):: delta_t_in, x_in !>initial time step, initial value of variable
  real(8), intent(in)::tolerance  !> The allowed relative error
   real(8), intent(out)::x_fin !> The final position, robust, error controlled result, 0.5 strong order 
  real(8), intent(out)::  delta_t_out,output_error  
  !The new, possibly reduced, adapted, value of the time step, new position and its  relative error  
  !real(8)::linearity_approx_limit,log_lin_apprx_lim
  !This subroutine is designed for x above 'linearity_approx_limit'.Below this value we wrote the
  ! subroutine RFOF_multi_step_boundarycase
  real(8):: nonlin_error_curent
  real(8):: delta_t_curent  
  real(8)::A,B 
  !Coefficients of the linear approximation
  real(8)::W, dW ! internal,W, dW are N(0,1) random numbers
  real(8), parameter:: lambda=1 ! parameter in the exact solution, to include
!  more general cases
    integer, intent(out)::Error_Flag !>Error signal
  integer, intent(out)::nr_subdivision !> the number of multiplying the time step with a time step reducing factor
  real(8)::sigma(3), sgma 
  !The square root of the diffusion coefficient at the selected points
  real(8)::xx(3) ! internal
  !The points were the diffusion coefficient is calculated 
  real(8)::sqrtdt, delta_x, xf, x_curent !internal
  integer, parameter:: Nr_x_subdivision=2
  real(8),parameter::division_factor_dx=1.d0/(Nr_x_subdivision+0.0d0)
  ! delta_t_new is delta_t_old/ Nr_time_division, to achieve accuracy
   real(8),parameter::division_factor_dt=division_factor_dx*division_factor_dx
 ! Factors for reducing the space step or time step, in order to achieve the requested linear aproximation accuracys 
   integer, parameter::maxdivision=50 !The allowed maximal subdivision , in order to
   !avoid looping
   
   x_curent=dabs(x_in)
 If(regularization_method==1)  call correction_domain_boundary(x_curent)

  Error_Flag=0
 ! print *,' ONE step start x_in=', x_in
 ! REMOVABLE IF
 ! if(x_curent.LE.0) then
 ! print*, ' bug in the input x_in in the subroutine RFOF_one_step_gen, x_in=',x_in
!  Error_Flag=2
 ! return
  !endif

  ! NORMAL CASE
  call RFOF_sigma(regularization_method, x_curent, sgma, Error_Flag) 
  sigma(2)=sgma
  if(Error_Flag.NE.0)  then
  print*, 'Signal from subroutine one_generalcase: Error flag=',Error_Flag
  return
  endif


   xx(2)=x_curent
  if(sigma(2)<=0) then
  !Not probable, but be sure
   x_fin=x_curent
   delta_t_out=delta_t_in
      Error_Flag=9
   return
  endif
  delta_t_curent=delta_t_in
  sqrtdt=sqrt(delta_t_curent)
  delta_x=sigma(2)*sqrtdt

   if (delta_x.LE.0) then
  print *, ' DIAGN: delta_t_in', delta_t_in, ' sigma2=',sigma(2)
  Error_flag=3
     return     
   endif

     nr_subdivision=0
    DO  
   ! In this main loop the time step is decreased such that the relative 
   ! error is less then the input value of "tolerance"
 
      if(nr_subdivision.EQ.maxdivision) then
            Error_Flag=1
          exit
      endif
  !sigma**2/2=Diffusion coefficient
   
       xx(1)=x_curent-delta_x ;xx(3)=x_curent+delta_x
  If(regularization_method==1) then 
        if(xx(1)<X_MIN) xx(1)=X_MIN
    else ! Regularization method==2
        if(xx(1)<0)  xx(1)=0
    endif

  call RFOF_sigma(regularization_method,xx(1), sgma, Error_Flag) 
  sigma(1)=sgma 
   
  if(Error_Flag.NE.0) exit
  call RFOF_sigma(regularization_method,xx(3), sgma, Error_Flag) 
  if(Error_Flag.NE.0) exit
  sigma(3)=sgma 
    
   call linear_aprox_minmax( xx, sigma,A,B,nonlin_error_curent, Error_Flag)
    if(Error_Flag>0) exit
     output_error=nonlin_error_curent
     if(nonlin_error_curent.LT.tolerance) then ! the time step was reduced sufficiently
     exit
  else 
        delta_t_curent= delta_t_curent*division_factor_dt
      delta_x=delta_x*division_factor_dx
      nr_subdivision=nr_subdivision+1
  endif
   
  enddo !end of the time step decreasing loop

  ! May be desactivated, for orientation only
  if (nr_subdivision>4)Then
      print*, '                TIME STEP Red factor=',division_factor_dt**nr_subdivision
      print*,' NsubdivionSTEP=', nr_subdivision
 print *, ' delta-x final=', delta_x
  endif
  ! PRINT*,' END  time step reducing'
    if(Error_Flag>0) return 
  ! The correct time step is established. Now we perform 1step , using linear fit and exact solution
  ! of the linear SDE, that approximates locally the non-linear SDE
 
    sqrtdt=sqrt(delta_t_curent) 

    IF(integration_method==1)  THEN
 ! We use the strong order 0.5 method       
         call Marsaglia_Bray_RM48( W )
          dW=W*sqrtdt
      Call RFOF_exact_solution(lambda,x_curent, delta_t_curent,dW,A,B, xf,Error_Flag)
       ELSE
          !Now the corection that uses the full information on diffusion coefficient, is added
             ! We use instead of linear interpolation a quadratic one 
           Call RFOF_SDE_integrator_corection( xx, sigma, delta_t_curent,xf,Error_Flag) 
      ENDIF

         delta_t_out=delta_t_curent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    xf=dabs(xf)
 If(regularization_method==1)  call correction_domain_boundary(xf)
     x_fin=xf
  
  return
end  subroutine RFOF_one_step_generalcase 


!-----------------------------------------------------------
 
subroutine RFOF_SDE_integrator_corection(xx, sigma, dt,  x_final, Error_Flag)
implicit none
real(8), intent(in)::xx(3), sigma(3) ! input for quadratic interpolation of sigma(x)
real(8), intent(in)::dt ! The time step
real(8):: x_in ! identical with xx(2), input value 
real(8), intent(out)::x_final ! The improved value
real(8)::coef(3), m,n,p !Coefficients  of the quadratic interpolation of sigma(x)== Sqrt(diffusion_coeff(x))
!coef(1)=p; coef(2)=n; coef(3)=m 
 real(8)::dw, dz 
 !Correlated Gaussian variables
 real(8)::normal_vector(2)
 !Uncorrelated  Gaussian variables
real(8):: s1,s, sms, dy1,dy2,dy3,dy4 , xf! internal variables
integer::Error_Flag ! Error signal,  equal 0 when OK 
integer::error ! Error signal from the quadratic interpolation subroutine

!call quadratic_interpolation( x_in,delta_x, sigma,m,n,p, Error_Flag)
!call quadratic_interpolation( xx, sigma,m,n,p, Error_Flag)
call polycoef_3points(xx,sigma,coef,error)
! We use the approximation $\sigma(x)=\sqrt(DiffCoef(x)=mx^2+nx+p)$
 !or $\sigma(x)=\sqrt(DiffCoef(x)=coef(3)x^2+coef(2)x+coef(1)$ 
 ! Link 
 x_in=xx(2)
 m=coef(3); n=coef(2) ; p=coef(1)
 ! End link
 if(error.NE.0)  Error_Flag=15
if (Error_Flag>0) return 
  call gaussian_random_vectorRM48(normal_vector)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!DIAGN
  if(dt<0)  then
  print* ,'PROBLEMS IN CORRELATED GAUSS' 
  Error_Flag=1111
  return
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
xf=x_in+dy1+dy2+dy3+dy4
 x_final=xf
             
return
end subroutine RFOF_SDE_integrator_corection


!----------------------------------------------------------------------------------
subroutine RFOF_multi_step_generalcase( integration_method, regularization_method, delta_t_in,&
 delta_t_out, x_in,x_fin,time_in, time_fin, T_max, tolerance,Error_Flag)
  implicit none
  integer, intent(in)::integration_method
  !   The method used in one step integration method=1 => use of 0.5 strong order method
  !method=2 => use of 1.5 strong order method 
integer, intent(in)::regularization_method 
! Regularization used to deal with the behavior Diff_coef(x=0)=0 in the subroutine sigma()
   real(8), intent(in):: delta_t_in, x_in !>initial time step, initial value of variable
  real(8), intent(in)::tolerance
  ! The allowed, imposed relative error
  real(8), intent(out)::   x_fin 
  real(8)::   x_final_extrapol  ! Internal variable, the extrapolated return value of one -step integrator
  !The  new position  
  real(8), intent(in):: time_in, T_max 
   !The initial time, T_max is the allowend maximal time,
  real(8), intent(out):: time_fin
  ! We have: time_in<time_fin<T_max; The time when the particle leave the boundary domain.
  ! t_fin is the time when the movement is stopped due to some event, like too 
  ! close to singularity ( in this version is not used)   
  !integer::error_signal_linaprox
  real(8), intent(out):: delta_t_out ! The readjusted time step, according to the output "accuracy"
    real(8) :: output_error ! Internal, the actual error from the one-step subroutine. It's used 
! to readjust the time step.
 ! integer:: ncount, maxdivision !maximal subdivision allowed for halving the time step
 ! real(8)::linearity_approx_limit !,log_lin_apprx_lim
  !This subroutine is designed for x above 'linearity_approx_limit'.Below this value we wrote the
  ! subroutine RFOF_multi_step_boundarycase
  real(8):: delta_t_curent ! the adjusted ccurent time step
   real(8)::time_curent
   real(8)::x_curent,x_new !internal variables
  ! real(8)::x_gpar,x_par ! For method 1, to store the previous value of x_curent
   real(8)::relative_error !internal variable returned by one step subroutine.The estimated relative error
  integer::nr_subdivision !internal variable
 ! logical::end_signal !True if the time is close to final time, up to one time step
  integer, intent(out)::Error_Flag
  real(8):: nstep ! internal, number of steps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! real(8):: ntot_subdiv, N_tot_call !For test runs

 ! logical::stopped ! common with one-step integrator
  
  !The error signal Error_Flag=0 when OK.   
  !Error_Flag=1 when by succesive reducing the time step the accuracy is sstill greater then tolerance
  !Error_Flag=2 when the input initial position is wrong
  !Error_Flag=3 when the input time step is wrong
  !Error_Flag=4 when the input initial time is wrong
  !Error_Flag=7: Too large exponent in the subroutine RFOF_Exact_Solution
  !Error_Flag=8 when the error signal from called subroutine  linear_aprox_minmax2  is nonzero
  !Error_Flag=9 Improbable value of the diffusion coefficient
  !Error_Flag=10 when the error signal from called subroutine  for diffusion coefficient  is nonzero
  !Error_Flag=11, when the normalized magnetic momentum, "x", is negative or zero (x=0 is a sticking point )
  !Error_Flag=12, when the values of the diffusion coefficient are zero
  !Error_Flag=14 when the integer "smoothing_method" from subroutine RFOF_sigma is not 1 or2
  !Error_Flag=15  anomalous return to subroutine RFOF_SDE_integrator_corection from the quadratic interpolatin
  ! subroutine polycoef_3points, stored in the file RFOF_numerics.f90 (module RFOF_numerics)
  !Error_Flag=16 The initial point is outside the allowed domain, imposed by regularization_method 1
  ! The signal is from subroutine RFOF_multi_step_generalcase  
!DATA maxdivision /10/ !> The time step delta_t si decreased at most by a factor 2**maxdivision

Error_Flag=0

  if(x_in.LT.0) then
  print*, ' bug in the input x_in in the subroutine RFOF_multi_step_generalcase' 
  Error_Flag=2
  return
  endif
  if(delta_t_in.LE.0) then
  print*, ' bug in the input delta_t in the subroutine RFOF_multi_step_generalcase'
  Error_Flag=3
  return
  endif
  
  if ((time_in.LT.0).OR.(time_in.GT.T_max)) then
  print*, ' bug in the input time_in the subroutine RFOF_multi_step_generalcase'
  Error_Flag=4
  return
  endif 

   time_curent=time_in
    x_curent= x_in
    delta_t_curent=delta_t_in

      if(aprox_equal(time_curent,T_max)) then
           time_fin=T_max
           x_fin=x_curent
           delta_t_out=delta_t_curent
           return
       endif
  delta_t_curent= dmin1(delta_t_in, T_max-time_curent)

  DO   ! main time loop
!Adjusting the time step to reach exactly the final value T_max, when time_curent is close to T_max
     delta_t_curent= dmin1(delta_t_curent, T_max-time_curent)
      call RFOF_one_step_generalcase(integration_method, regularization_method, delta_t_curent,&
      x_curent, tolerance,delta_t_out,x_new, output_error, nr_subdivision, Error_Flag)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   delta_t_curent=delta_t_out
   x_curent=x_new
   time_curent=time_curent+delta_t_curent
   if(aprox_equal(time_curent,T_max)) then
   time_curent=T_max
     exit
   endif
   nstep=nstep+1.0d0
          if(Error_Flag>0) exit
! Increase the time step when it is too small, for a given accuracy
 
        if( output_error.LT.tolerance*1.0d-1)  then
           delta_t_curent=delta_t_curent*2.0d0
       endif
 !print *, ' end multistep time loop'
 ENDDO ! end of main time loop 

     if(Error_Flag.NE.0) then
      print *, 'Anomalous return'
      if(Error_Flag==1) print *,' Too much accuracy required'
        return
     endif

     delta_t_out=delta_t_curent 
     x_fin=x_curent
    time_fin=time_curent
    return
  end subroutine RFOF_multi_step_generalcase

!----------------------------------------------------
 
end module RFOF_SpecialSolver5



  !--------------------------------------------------------------------------------
  !  subroutine polycoef_3points
  !--------------------------------------------------------------------------------


 ! Extracted from module RFOF_numerics

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

    y1 = x(1) - 0.5d0*t(1)**2 * xbiss
    y3 = x(3) - 0.5d0*t(3)**2 * xbiss

    xprim  = (y3-y1) / ( t(3) - t(1) )

    x0 = y1 - xprim*t(1)

    coef(1) = x0
    coef(2) = xprim
    coef(3) = 0.5d0 * xbiss

    !    print *, 'polycoef', xp1,xp2,t

  end subroutine polycoef_3points

