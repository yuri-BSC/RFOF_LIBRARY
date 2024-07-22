!-----------------------------------------------------------------------------------------
! Test program for the subroutine Bessel0_General from the file RFOF_BesselProcess0.f90 
!Author:Gyorgy Steinbrecher, EURATOM-MEdC 

program RFOF_BesselProcessTest2Main
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
real(8)::C_dif, D_0 !Input parameter in subroutine Bessel0
!From local aproximation: Diff_coef=2*C_dif*x+D_0
real(8)::diff_coef ! diffusion coefficient
integer, parameter::nr_moments=2 
!Number of mean values (moments of exp random variable) computed, to be compared for test
real(8):: moments_ex(nr_moments)
! exact values, to be compared with the computed
real(8):: mean_moment(nr_moments),dispersion(nr_moments) ! Computed by MonteCarlo
real(8):: loc_mean_moments( nr_moments),loc_dispersion( nr_moments)
integer::Error_Flag 
integer:: k_moments 
real(8)::time_in,time_fin, cputime
integer::ios !internal , used for input-output
integer:: k_test, Nr_test ! number of tests
character(len=10)::date , time, zone ! For timing, text output
integer::values(8)  ! For timing, integer output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character (len=*),parameter::Out_filename='TestBessel0General_4.txt'

real(8):: sum1, sum2,ww(2), ntr,xx, z1,z2,zz, sumz1, sumz2
! For testing random number generators, at the end
 !!!!!!!!!!
!INPUT DATA
 !!!!!!!!! 
 Nr_test=20
 ntraject=1.0d5
n_group=100
x_in=0.134d0  ! 3.23456789d-8! Initial value of trajectory
delta_t=0.8200d0
C_dif=0.3d0    ; D_0=1.010987d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!END INPUT DATA

open(UNIT=3,file=Out_filename,iostat=ios,err=202,form='formatted',status='new',action='readwrite')
202 Print*, 'file opening status: IOS=',ios
601 format(8I5)
602 format(a30)
write(UNIT=3, fmt=200)'TEST: Subroutine Bessel0_General from module RFOF_BesselProcess0' 
write(*, fmt=200)'TEST: Subroutine Bessel0_General from module RFOF_BesselProcess0' 
 
call date_and_time(date, time, zone, values)
Write(*, fmt=602)' START TIME-file default unit' 
Write(*, fmt=601)values
Write(UNIT=3, fmt=602)' START TIME-file unit=3' 
Write(UNIT=3, fmt=601)values
write(UNIT=3, fmt=200)'INPUT DATA' 
write(*, fmt=200)'INPUT DATA' 
200 format(a80)
write(UNIT=3, fmt=201) 'Ntraject=', ntraject,' n_group=',n_group, ' start x=',x_in,' delta_t=', delta_t 
write(*, fmt=201) 'Ntraject=', ntraject,' n_group=',n_group,'  start x=',x_in, ' delta_t=', delta_t 
201 format(a12, d12.4,a12, I4,a12,d12.4,a12, d12.4)
300 format(2(a11,d14.7) )
write(UNIT=3, fmt=204)'RESULTS: ' 
write(*, fmt=204)'RESULTS' 
204 format(a20)
CALL CPU_TIME ( time_in)
C_dif=2.0d0; D_0=3.2d0

!!!!!!!!!!!!!!!!!!!!Loop to modify C_diff and D_0
DO k_test=1, Nr_test
write(UNIT=3, fmt=204)'  ' 
write(*, fmt=204)' '
C_dif=C_dif-0.5d0; D_0=D_0+1.2d0
diff_coef=2.d0*C_dif*x_in+D_0
If (diff_coef.LE.0) then
print*,' Wrong input, diffusion coefficient not positive:', diff_coef
goto 1022
endif
write(UNIT=3,fmt=300)'C_dif=',C_dif, ' D_0=',D_0  
write(*,fmt=300)'C_dif=',C_dif,' D_0=',D_0   

 ! Computation
do k_moments=1,nr_moments
mean_moment(k_moments)=0.0d0 
dispersion(k_moments)=0.0d0
enddo
 Error_Flag=0
call Bessel0_Monte_Carlo(ntraject, n_group, x_in,delta_t,C_dif,D_0,&
                                       &loc_mean_moments, loc_dispersion,Error_Flag)
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

111 format('Order moment=',I3,' Exact moment=',D11.4, ' computed=',&
& D11.4, ' Error=',D11.4)

12 format(' CPUtime=',D11.4 )
! Write moments
Call exact_momentsG(x_in, delta_t,C_dif,D_0, moments_ex) 
 do  k_moments=1,nr_moments
    !    ! write(UNIT=3, fmt=204)'  ' 
     !write(*, fmt=204)' ' 
     !write(*,10) k_moments
     !write(UNIT=3,fmt=10) k_moments
write(*, 111)k_moments, moments_ex(k_moments), mean_moment(k_moments), dispersion(k_moments)
write(UNIT=3, fmt=111)k_moments, moments_ex(k_moments),mean_moment(k_moments), dispersion(k_moments)

 enddo
 !-------------------------------------------------
1022  Continue
 END DO ! End loop for k_test
  !END WRITE MOMENTS
  ! BEGIN WRITE CPU-TIME
 CALL CPU_TIME ( time_fin)
 cputime=time_fin-time_in

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
sumz1=0; sumz2=0
!Exact values to be returned:sum1=2 and sum2=8,also the independence of the components is tested
ntr=0
do 
ntr=ntr+1
call gaussian_random_vectorRM48(ww)
Call Marsaglia_Bray_RM48(z1)
Call Marsaglia_Bray_RM48(z2)
xx=ww(1)**2+ww(2)**2
zz=z1*z1+z2*z2
sum1=sum1+xx
sumz1=sumz1+zz
xx=xx*xx
zz=zz*zz
sum2=sum2+xx
sumz2=sumz2+zz
if(ntr>Ntraject) exit
enddo
 sum1=sum1/Ntraject
 sum2=sum2/Ntraject
 print*,'sum1=',sum1,'  sum2=',sum2
sumz1=sumz1/Ntraject
 sumz2=sumz2/Ntraject
 print*,'sumz1=',sumz1,'  sumz2=',sum2

1111 stop  ' End of test'
end program  RFOF_BesselProcessTest2Main
!---------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!------------------------------------------------
! TEST PROGRAMS 
!---------------------------------------------------------------------------

subroutine  Bessel0_Monte_Carlo( ntraject, n_group,x_in,delta_t,C_dif,D_0, loc_mean_moments, loc_dispersion,Err_Flag)
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
real(8), intent(in)::  C_dif,D_0! input in the tested Bessel procees, input for subroutine bessel0
!From local aproximation: Diff_coef=2*C_dif*x +D_0+error
real(8) ::   x_fin   !The  new position, output of subroutine Bessel0 
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
      ! Error_Flag=0
  DO 90 group=1,n_group
   do k = 1, nr_moments
     loc_moments(group, k)=0; ! loc_mean_moments2(k)=0
     enddo
count_traject=0
!Start of Monte-Carlo sampling 
DO ! Main loop, traject
count_traject=count_traject+1.d0 
call Bessel0_General(x_in, delta_t, C_dif,D_0, x_fin)  
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
