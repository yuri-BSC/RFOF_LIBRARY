!-----------------------------------------------------------------------------------------
! Test program for the subroutine Bessel0 from the file RFOF_BesselProcess0.f90 
!Author:Gyorgy Steinbrecher, EURATOM-MEdC 

program RFOF_BesselProcessTestMain
use RFOF_BesselProcess0
use RFOF_Random_Numbers_SG2
implicit none
intrinsic CPU_TIME 
real(8):: ntraject ! Number of MC trajectory sampling
integer::n_group 
!Number of  subgroups of trajectories, used to estimate the Monte-Carlo error
!Total number of trajectories=n_group*ntraject
real(8)::x_in
! The initial point of the trajectories,
real(8)::delta_t !The time step in the tested adaptive methods
real(8)::C_dif   !Input parameter in subroutine Bessel0
!From local aproximation: Diff_coef=C_dif*x+O(x**2)

integer, parameter::nr_moments=2 
!Number of mean values (moments of exp random variable) computed, to be compared for test
real(8):: moments_ex(nr_moments)
! exact values, to be compared with the computed
real(8):: mean_moment(nr_moments),dispersion(nr_moments) ! Computed by MonteCarlo
real(8):: loc_mean_moments( nr_moments),loc_dispersion( nr_moments)
integer::Error_Flag 
integer:: k_moments 
real(8):: time_in,time_fin, cputime
integer::ios !internal , used for input-output
character(len=10)::date , time, zone ! For timing, text output
integer::values(8)  ! For timing, integer output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character (len=*),parameter::Out_filename='B0Test1.txt'

real(8):: sum1, sum2,ww(2), ntr,xx
! For testing random number generators, at the end
 !!!!!!!!!!
!INPUT DATA
 !!!!!!!!!
 ntraject=1.0d5
n_group=100
x_in=1.234d1  ! 3.23456789d-8! Initial value of trajectory
delta_t=5.200d2
C_dif=4.000d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!END INPUT DATA
! More input,X_MIN, mean_diff_coef, and small_adim_diff_coef_fact
!in this version, is at the heading of RFOF_SpecialSolver5

open(UNIT=3,file=Out_filename,iostat=ios,err=202,form='formatted',status='new',action='readwrite')
202 Print*, 'file opening status: IOS=',ios
601 format(8I5)
602 format(a30)
call date_and_time(date, time, zone, values)
Write(*, fmt=602)' START TIME-file default unit' 
Write(*, fmt=601)values
Write(UNIT=3, fmt=602)' START TIME-file unit=3' 
Write(UNIT=3, fmt=601)values
write(UNIT=3, fmt=200)'INPUT DATA' 
write(*, fmt=200)'INPUT DATA' 
200 format(a50)
write(UNIT=3, fmt=201) 'Ntraject=', ntraject,' n_group=',n_group, ' start x=',x_in,' delta_t=', delta_t 
write(*, fmt=201) 'Ntraject=', ntraject,' n_group=',n_group,'  start x=',x_in, ' delta_t=', delta_t 
201 format(a12, d12.4,a12, I4,a12,d12.4,a12, d12.4)
300 format(a17,d12.4)
write(UNIT=3,fmt=300)'C_dif const=',C_dif  
write(*,fmt=300)'C_dif=',C_dif  
write(UNIT=3, fmt=204)'RESULTS: ' 
write(*, fmt=204)'RESULTS' 
204 format(a20)

 ! Computation
  !______________________________________________________________________
do k_moments=1,nr_moments
mean_moment(k_moments)=0.0d0 
dispersion(k_moments)=0.0d0
enddo
CALL CPU_TIME ( time_in)
 Error_Flag=0
call Bessel0_Monte_Carlo(ntraject, n_group, x_in,delta_t,C_dif, loc_mean_moments, loc_dispersion,Error_Flag)
 CALL CPU_TIME ( time_fin)
 cputime=time_fin-time_in
 if(Error_Flag==0 ) then
      do k_moments=1,4
        mean_moment(k_moments)=loc_mean_moments(k_moments)
        dispersion(k_moments)=loc_dispersion(k_moments)
      enddo
 else 
      9 format(' err flag=',i2)
      write(*,9)Error_Flag 
      write(UNIT=3,fmt=9)Error_Flag
 endif ! Error flag if
 
! END COMPUTATION
! BEGIN WRITE MOMENTs

10 format('ORDER MOMENTS=', I3)
11 format('Exact moment=',D11.4, ' computed moment=',D11.4, ' Disp=',D11.4)
12 format(' CPUtime=',D11.4 )
! Write moments
Call exact_moments(x_in, delta_t,C_dif, moments_ex) 
 do  k_moments=1,nr_moments
     write(UNIT=3, fmt=204)'  ' 
     write(*, fmt=204)' ' 
     write(UNIT=3, fmt=204)'  ' 
     write(*, fmt=204)' ' 
     write(*,10) k_moments
     write(UNIT=3,fmt=10) k_moments
write(*, 11) moments_ex(k_moments), mean_moment(k_moments), dispersion(k_moments)
write(UNIT=3, fmt=11) moments_ex(k_moments),mean_moment(k_moments), dispersion(k_moments)

 enddo
 !-------------------------------------------------
  !END WRITE MOMENTS
  ! BEGIN WRITE CPU-TIME
 write(UNIT=3, fmt=204)'  ' 
write(*, fmt=204)' ' 
write(UNIT=3,fmt=206) ' CPU Time '
write(*,fmt=206) ' CPU Time '
206 format(A20)
write(UNIT=3, fmt=204)'  ' 
write(*, fmt=204)' ' 
write(*, 12)  cputime
write(UNIT=3,fmt= 12)cputime
write(UNIT=3, fmt=205)'ERROR_FLAG=', Error_Flag
write(*, fmt=205)'ERROR_FLAG=', Error_Flag
205 format(a20,I4)
 call date_and_time(date, time, zone, values)
Write(*, fmt=602)'FINAL TIME ' 
Write(*, fmt=601)values
Write(UNIT=3, fmt=602)'FINAL TIME' 
Write(UNIT=3, fmt=601)values
call CPU_TIME(time_fin)
close(UNIT=3,iostat=ios,err=1001,  status='keep' )
1001 print *, 'File closed, iostat=',ios

! Test again, to be sure, of the Gaussian generator
sum1=0 ;sum2=0;
!Exact values to be returned:sum1=2 and sum2=8,also the independence of the components is tested
ntr=0
do 
ntr=ntr+1
call gaussian_random_vectorRM48(ww)
xx=ww(1)**2+ww(2)**2
sum1=sum1+xx
xx=xx*xx
sum2=sum2+xx
if(ntr>Ntraject) exit
enddo
 sum1=sum1/Ntraject
 sum2=sum2/Ntraject
 print*,'sum1=',sum1,'  sum2=',sum2
 stop  ' End of test'
end program  RFOF_BesselProcessTestMain
!---------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!------------------------------------------------
! TEST PROGRAMS 
!---------------------------------------------------------------------------

subroutine  Bessel0_Monte_Carlo( ntraject, n_group,x_in,delta_t,C_dif, loc_mean_moments, loc_dispersion,Err_Flag)
use RFOF_BesselProcess0
use RFOF_random_numbers
use RFOF_random_numbers_SG2
!Test the solver for the SDE , by computing mean values in stationary state
! and comparing with the exact values given by simpest Euler method
! Test of subroutine RFOF_multi_step_boundarycase for y_min=0,  
implicit none
intrinsic CPU_TIME !<used for speed test, and seed initialization in the random number generator
real(8), intent(in)::ntraject !>Number of trajectories in MC sampling
real(8) ::count_traject !>counter of trajectories in MC sampling
real(8), intent(in) :: delta_t, x_in !> time step, initial value of variable
  ! `in the tested Bessel procees, input for subroutine bessel0
real(8), intent(in)::  C_dif! input in the tested Bessel procees, input for subroutine bessel0
real(8) :: x_fin   !The  new position, output of subroutine Bessel0 
integer,intent(out)::Err_Flag
  !The error signal Error_Flag=0 when OK.
integer, parameter::nr_moments=2
integer, parameter::nr_max_group=1000
real(8):: loc_moments(nr_max_group,nr_moments) !local, moments for each group real(8):: loc_mean_moments2(nr_moments) 
! Local variable used to compute the dispersion
real(8), intent(out):: loc_mean_moments(nr_moments) ! Mean values of previous moments(1000,4)
real(8):: loc_mean_moments2(nr_moments) !Internal, used for compute dispersion
real(8), intent(out)::loc_dispersion(nr_moments) !  Mean square error in previous moments(10,4)
integer:: k,m,group ! internal
real(8)::time_begin,time_end, cputime 
integer, intent(in):: n_group !Number of group in Monte-Carlo  error estimation
real(8):: disp , m1 !internal, for computing the dispersion
! Fix the input values 
!print*, 'Semnal 3'
Err_Flag=0
!print*, 'after error flag'
if( n_group>nr_max_group) then
print *,' ERROR in main: n_group too large '
return
endif
!Begin test 
! Initialise the seed of random number generators
call init_RANDOM_NUMBER(0) 
  CALL CPU_TIME ( time_begin )  
      print*, 'START' 
      ! Error_Flag=0
  DO 90 group=1,n_group
   do k = 1, nr_moments
     loc_moments(group, k)=0; ! loc_mean_moments2(k)=0
     enddo
count_traject=0
!Start of Monte-Carlo sampling 
DO ! Main loop, traject
count_traject=count_traject+1.d0 
call Bessel0(x_in, delta_t, C_dif, x_fin)  
m1=x_fin 
loc_moments(group,1)=loc_moments(group,1)+x_fin 
loc_moments(group,2)=loc_moments(group,2)+x_fin*x_fin 

if (count_traject.GE.ntraject)  exit

enddo !!!!!!!!!!!!!!!!!!!!!!!!!!!! INNER MAIN LOOP 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 End MAIN traject counting loop 
if (Err_Flag.NE.0) then
Print *, 'group=',group,  'Error_Flag=',Err_Flag
 endif
 90 CONTINUE !!!!!!!!!!!!!!!!!!!!!!!!End group loop
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Outer main loop
if(Err_Flag.NE.0) then
!Anomalous ending
print *,' ERROR SIGNAL=',Err_Flag
return
endif 

  do group=1,n_group
   do k=1,nr_moments
loc_moments(group, k)=loc_moments(group,k)/ntraject 
 enddo
 enddo 

 
 do k=1,nr_moments
loc_mean_moments(k)=0
loc_mean_moments2(k)=0
do group=1,n_group
loc_mean_moments(k)=loc_mean_moments(k)+ loc_moments(group, k)
loc_mean_moments2(k)=loc_mean_moments2(k)+ loc_moments(group,k)**2
enddo
loc_mean_moments(k)=loc_mean_moments(k)/n_group 
loc_mean_moments2(k)=loc_mean_moments2(k)/n_group 
 disp=loc_mean_moments2(k)-loc_mean_moments(k)**2
if (disp>=0) then
loc_dispersion(k)= sqrt(disp/(n_group-1)  )
else
print*, ' ERROR in at moment order=',k
endif
 enddo

end subroutine Bessel0_Monte_Carlo
