Module RFOF_BesselProcess0 
use RFOF_stochastics
contains

!----------------------------------------------------------------------------
! Subroutine Bessel0(x_in,delta_t,C_dif,x_fin)
Subroutine Bessel0(x_in,delta_t,C_dif,x_fin)
use RFOF_stochastics
implicit none
!Input
real(8), intent(in)::x_in !initial position
real(8), intent(in)::delta_t !the time interval
real(8),intent(in)::C_dif !From local aproximation: Diff_coef=2*C_dif*x+O(x**2)
!Output
real(8), intent(out)::x_fin ! the final position
!Local
real(8):: grw(2) !Random normal vector, uncorrelated components
real(8)::x1,x2,csqrtdt ! local variables
csqrtdt=sqrt(delta_t*C_dif)
call gaussian_random_vectorRM48(grw)
x1=grw(1)*csqrtdt+sqrt(x_in)
x2=grw(2)*csqrtdt
x_fin=x1*x1+x2*x2
return
end Subroutine Bessel0
!------------------------------------------------------------------------


!-----------------------------------------------------------------
! Subroutine Bessel0_General(x_in,delta_t,C_dif,D_0,x_fin)
Subroutine Bessel0_General(x_in,delta_t,C_dif,D_0,x_fin)
use RFOF_stochastics
implicit none
!Input
real(8), intent(in)::x_in !initial position
real(8), intent(in)::delta_t !the time interval
real(8),intent(in)::C_dif, D_0 !From local aproximation: Diff_coef=2*C_dif*x +D_0+error

!Output
real(8), intent(out)::x_fin ! the final position
!Local
real(8):: grw(2) !Random normal vector, uncorrelated components
real(8)::normal_random ! random normal numbers
real(8)::x1,x2,csqrtdt, sqrtxin,d2c ! local variables
if(C_dif.NE.0) Then 
d2c=D_0/(2.0d0*C_dif)
csqrtdt=sqrt(delta_t*dabs(C_dif) ) 
   if(C_dif>0) Then
             sqrtxin=sqrt(d2c+x_in)
             call gaussian_random_vectorRM48(grw)
             x1=grw(1)*csqrtdt+sqrtxin
             x2=grw(2)*csqrtdt
              x_fin=x1*x1+x2*x2-d2c
     else    ! C_dif<0        
             sqrtxin=sqrt(-d2c-x_in)
             call gaussian_random_vectorRM48(grw)
             x1=grw(1)*csqrtdt+sqrtxin
             x2=grw(2)*csqrtdt
              x_fin=-x1*x1-x2*x2-d2c
          endif          
else
! C_dif=0
call Marsaglia_Bray_RM48(normal_random)
x_fin=x_in+sqrt(2.0d0*delta_t*D_0)*normal_random
endif
return
end Subroutine Bessel0_General
!


!------------------------------------------------------------------------
!  subroutine exact_moments(x_in, delta_t,C_dif, moments)
subroutine exact_moments(x_in, delta_t,C_dif, moments)
implicit none
!Input
real(8), intent(in)::x_in !starting point of  0 order bessel process
real(8), intent(in)::delta_t ! duration of Bessel procees
real(8), intent(in)::C_dif!From local aproximation: Diff_coef=2*C_dif*x 
!Fxed parameter, adjustable
integer, parameter::max_moment_order=2!Returns all moments from 1 to this value
!Output
real(8), intent(out)::moments(max_moment_order)
! The first 2 exact moments of the 0 order Bessel Procees, starting at , used for test 
!local
real(8) :: cdt ! local
cdt=C_dif*delta_t
moments(1)=x_in+2.d0*cdt
moments(2)=8.0d0*cdt*cdt+8.0d0*cdt*x_in+x_in*x_in
return
end subroutine exact_moments

!------------------------------------------------------------------------
! subroutine exact_momentsG(x_in, delta_t,C_dif,D_0 , moments)
subroutine exact_momentsG(x_in, delta_t,C_dif,D_0 , moments)
implicit none
!Input
real(8), intent(in)::x_in !starting point of  0 order Bessel process
real(8), intent(in)::delta_t ! duration of Bessel procees
real(8), intent(in)::C_dif,D_0 !From local aproximation: Diff_coef=2*C_dif*x +D_0
!Fixed parameter, adjustable
integer, parameter::max_moment_order=2!Returns all moments from 1 to this value

!Output
real(8), intent(out)::moments(max_moment_order)
! The first 2 exact moments of the 0 order Bessel Procees, starting at , used for test 
!Local
real(8) :: cdt ! local
cdt=C_dif*delta_t
moments(1)=x_in+2.d0*cdt
moments(2)=8.0d0*cdt*cdt+8.0d0*cdt*x_in+x_in*x_in+2.d0*D_0*delta_t
return
end subroutine exact_momentsG
!-------------------------------------------------------------------------
!------------------------------------------------------------------------
!Subroutine RFOF_MeanExitTime1(a,x_0,b,relative_tolerance,C_dif, &
!             &mean_first_exit_time, Error_Flag)


Subroutine RFOF_MeanExitTime0(x_0,x_min,C_dif,Time_exit,ErrorFlag)
Implicit none
!Input
real(8), intent(in)::x_0!starting point of  0 order Bessel process
real(8), intent(in)::C_dif !From local aproximation: Diff_coef=2*C_dif*x O(x**2)
real(8), intent(in)::x_min ! Domain when the previos linear approximation holds
!Output
Real(8), intent(out)::Time_exit ! The mean first exit time
Integer, intent(out)::ErrorFlag ! ErrorFlag=0 when 0<=x_0<=x_min;ErrorFlag=101 otherwise
if( (0.LE.x_0).AND.(x_0.LE.x_min) ) then
    ErrorFlag=0
Else
    ErrorFlag=101
endif
Time_exit=(x_min-x_0)/(2.d0*C_dif)
return
End Subroutine RFOF_MeanExitTime0
!---------------------------------------------------------------------------------




!------------------------------------------------------------------------
!Subroutine RFOF_MeanExitTime1(a,x_0,b,relative_tolerance,C_dif, D_0,&
!             &mean_first_exit_time, Error_Flag)


Subroutine RFOF_MeanExitTime1(a,x_0,b,relative_tolerance,C_dif, D_0,&
             &mean_first_exit_time, Error_Flag)
Implicit none
!Input
real(8), intent(in)::a,b ! Endpoints of the interval
real(8), intent(in)::x_0!starting point of  0 order Bessel process
real(8), intent(in)::relative_tolerance ! Admissible relative error
real(8), intent(in)::C_dif,D_0 !From local aproximat:ion: Diff_coef=2*C_dif*x +D_0
!Output
Real(8), intent(out):: mean_first_exit_time ! The mean first exit time from (a,b) starting
!at x_0
integer, intent(out)::Error_Flag !zero when OK, 102 if not a<=x_0<=b, wrong ordering
! 103 if the diffusion coefficient is negative on interval (a,b)
real(8), parameter::tolerance_factor=1.d0 ! Adjustable parameter
!Local
real(8)::aa, bb,xx,cc,dif_coef_a,dif_coef_b
if((a<0).OR.(x_0.LE.a).OR.(x_0.GE.b) ) then
    Error_Flag=102
    return
endif
dif_coef_a=2.d0*C_dif*a+D_0
if(dif_coef_a<0)  then
    Error_Flag=103
    return
endif

dif_coef_b=2.d0*C_dif*b+D_0
if(dif_coef_b<0)  then
    Error_Flag=103
    return
endif
Error_Flag=0
IF(C_dif.NE.0) THEN

    IF((a.NE.0).OR.(D_0.NE.0)) THEN
                 if((b-a)/a<relative_tolerance*tolerance_factor) then
                    ! Small interval, tipical , better to use approximation
                     mean_first_exit_time=(b-x_0)*(x_0-a)/(2.d0*(2.0d0*C_dif*x_0+D_0))
                  ! The simple, quick approximate formula is used
                   else ! No aproximation, large interval (a,b) 
            ! General  case exact formula, more time consuming
                   aa=dlog(dif_coef_a)
                   bb=dlog(dif_coef_b)
                   xx=dlog(2.d0*x_0*C_dif+D_0 )
                   cc=2.d0*C_dif*(aa-bb)
                   mean_first_exit_time=-x_0/(2.d0*C_dif)+(b*aa-a*bb+(a-b)*xx)/cc
                 endif
     ELSE ! a=D_0=0   
                    mean_first_exit_time=(b-x_0)/(2.d0*C_dif)
                  ! The local approximation near x=0 is made by enforcing D_0==0
     ENDIF   ! a==0 or not


 ELSE  ! C_dif=0
     if(D_0.NE.0) then
 mean_first_exit_time= (x_0-a)*(b-x_0)/(2.0d0*D_0)
      else
          mean_first_exit_time=1.0d99
          Error_Flag=103
      endif
ENDIF  ! C_dif = or not 0    
return
End Subroutine RFOF_MeanExitTime1

!----------------------------------------------------------------------------
!------------------------------------------------------------------------
!Subroutine RFOF_MeanExitTime2(a,x_0,b,relative_tolerance,C_dif, D_0,&
!             &mean_first_exit_time, Error_Flag)


Subroutine RFOF_MeanExitTime2(a,x_0,b,relative_tolerance,C_dif, D_0, mean_first_exit_time)
Implicit none
!Input
real(8), intent(in)::a,b ! Endpoints of the interval
real(8), intent(in)::x_0!starting point of  0 order Bessel process
real(8), intent(in)::relative_tolerance ! Admissible relative error
real(8), intent(in)::C_dif,D_0 !From local aproximat:ion: Diff_coef=2*C_dif*x +D_0
!Output
Real(8), intent(out):: mean_first_exit_time ! The mean first exit time from (a,b) starting
!at x_0
real(8), parameter::tolerance_factor=1.d0 ! Adjustable parameter
!Local
real(8)::aa, bb,xx,cc

if(C_dif.NE.0) THEN

     if((a.NE.0).OR.(D_0.NE.0)) then
              if((b-a)/a<relative_tolerance*tolerance_factor) then
                    ! Small interval, tipical , better to use approximation
                 mean_first_exit_time=(b-x_0)*(x_0-a)/(2.d0*(2.0d0*C_dif*x_0+D_0))
                  ! The simple, quick approximate formula is used
                   else ! No aproximation, large interval (a,b) 
            ! General  case exact formula, more time consuming
                   aa=dlog(2.d0*a*C_dif+D_0)
                   bb=dlog(2.d0*b*C_dif+D_0)
                   xx=dlog(2.d0*x_0*C_dif+D_0 )
                   cc=2.d0*C_dif*(aa-bb)
                   mean_first_exit_time=-x_0/(2.d0*C_dif)+(b*aa-a*bb+(a-b)*xx)/cc
                 endif
       else ! a=D_0=0   
                    mean_first_exit_time=(b-x_0)/(2.d0*C_dif)
                  ! The local approximation near x=0 is made by enforcing D_0==0
       endif   ! a==0 or not


else  ! C_dif=0
 mean_first_exit_time= (x_0-a)*(b-x_0)/(2.0d0*D_0)  
ENDIF  ! C_dif = or not 0    
return
End Subroutine RFOF_MeanExitTime2





!-----------------------------------------------------------------------------
!SUBROUTINE RFOF_Bessel_0_move(a,x_init,b,max_duration,relative_tolerance,C_dif,&
!        D_0,duration,x_final, Error_Flag)



SUBROUTINE RFOF_Bessel_0_move(a,x_init,b,max_duration,relative_tolerance,C_dif,&
        D_0,duration,x_final, Error_Flag)
Implicit none
!Input
real(8), intent(in)::a,b ! Endpoints of the interval
real(8), intent(in)::x_init!starting point of  0 order Bessel process
real(8), intent(in)::max_duration ! The maximal time allowed
real(8), intent(in)::relative_tolerance ! Admissible relative error
real(8), intent(in)::C_dif,D_0 !From local linear  aproximation: Diff_coef=2*C_dif*x +D_0
!Output
Real(8), intent(out)::duration!The first exit time from (a,b) starting
                              !at x_0
real(8),intent(out)::x_final ! final position
Integer, intent(out)::Error_Flag !ErrorFlag=0 when a<=x_0<=b;ErrorFlag=102
!When the Diffusionn coefficient is negative on (a,b), then Error_Flag=103
!otherwise Error_Flag=0
!Local
integer::error_flag_exittime
real(8), parameter:: max_iterations=1.0d6 !5 ! maximal number of iterations allowed
real(8), parameter::time_fraction=5.0d-5  ! 2.0d-4 !  time step compared do mean exit time
! see definition of delta_t

!For test, by first exit time comparation, these 2 parameters will be set very large
real(8):: mean_exit_time1,time_step,time_step_prop ,  time_curent, time_remained
real(8)::  x_curent1,x_curent2
real(8)::count_iterations 
        
              x_curent1=x_init
              time_curent=0
              count_iterations=0

     call RFOF_MeanExitTime1(a,(a+b)*0.5d0,b,relative_tolerance,C_dif, D_0,mean_exit_time1,error_flag_exittime)
     Error_Flag=error_flag_exittime
      time_step_prop=mean_exit_time1*time_fraction ! proposed time step, but maybe toolong so we compare
DO ! count_iteration=1, max_iterations
             count_iterations=count_iterations+1.0d0
            time_remained=max_duration-time_curent
            time_step=min(time_step_prop,time_remained) ! Actual time step
             call Bessel0_General(x_curent1,time_step,C_dif,D_0,x_curent2)
            time_curent=time_curent+time_step 
       if((x_curent2<a).OR.(x_curent2>b).OR.(time_curent.GE.max_duration)) exit
       x_curent1=x_curent2 
      IF(count_iterations.GE.max_iterations) exit
ENDDO
duration=time_curent
x_final=x_curent2
! REMOVE
  ! if(count_iterations>(max_iterations/100)*99) THEN
  ! print*, 'iterat=', count_iterations,' duration=',duration,' t-rmaind=',time_remained 
  ! endif
      ! End remove

return
END SUBROUTINE RFOF_Bessel_0_move

!________________________________________________________________

!SUBROUTINE RFOF_Bessel_0_move_test(a,x_0,b,max_duration,relative_tolerance,C_dif, D_0,duration,x_final, ErrorFlag)
!Simplified version, for test




SUBROUTINE RFOF_Bessel_0_move_test(a,x_0,b,max_duration,relative_tolerance,C_dif, D_0,duration,x_final, ErrorFlag)
Implicit none
!Input
real(8), intent(in)::a,b ! Endpoints of the interval
real(8), intent(in)::x_0!starting point of  0 order Bessel process
real(8), intent(in)::max_duration ! The maximal time allowed
real(8), intent(in)::relative_tolerance ! Admissible relative error
real(8), intent(in)::C_dif,D_0 !From local aproximation: Diff_coef=2*C_dif*x +D_0
!Output
Real(8), intent(out)::duration!The first exit time from (a,b) starting
                              !at x_0
real(8),intent(out)::x_final ! final position
Integer, intent(out)::ErrorFlag !ErrorFlag=0 when a<=x_0<=b;ErrorFlag=10 otherwise
                                 !)
!Local
integer, parameter:: max_iterations=50! maximal number of iterations allowed
integer, parameter::time_fraction=10!time step compared do mean exit time
real(8):: time1,time2, time3, delta_t ,  time_curent, x_curent1,x_curent2
integer::count_iteration
x_curent1=x_0
time_curent=0
DO count_iteration=1, max_iterations
!call RFOF_MeanExitTime(a,x,b,relative_tolerance,C_dif, D_0,time1,ErrorFlag)
!time2=time1/time_fraction ! proposed time step, but maybe toolong so we compare
time2=max_duration/time_fraction
time3=max_duration-time_curent
delta_t=min(time2, time3)
call Bessel0_General(x_curent1,delta_t,C_dif,D_0,x_curent2)
time_curent=time_curent+delta_t

if((x_curent2<a).OR.(x_curent2>b).OR.(time_curent.GE.max_duration)) exit
x_curent1=x_curent2
ENDDO
duration=time_curent
x_final=x_curent2
return
END SUBROUTINE RFOF_Bessel_0_move_test


End Module RFOF_BesselProcess0 
