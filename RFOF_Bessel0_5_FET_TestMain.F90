program RFOF_BesselProcessTest4main
  use RFOF_BesselProcess0
  use RFOF_Stochastics
  implicit none
  intrinsic CPU_TIME 
  real(8):: ntraject ! Number of MC trajectory sampling
  integer::n_group 
  !Number of  subgroups of trajectories, used to estimate the Monte-Carlo error
  !Total number of trajectories=n_group*ntraject
  real(8)::x_in 
  ! The initial point of the trajectories,
  real(8):: a,b ! The interval (a,b), to exit from x_in 
  real(8)::max_duration !The time step and the maximal alowed test time, mach
  !max_duration much larger cmpared to the  mean exit time
  real(8)::factor_max_duration ! factor_ max_duration*exact_FET=max_duration
  real(8)::C_dif, D_0 !Input parameter in subroutine Bessel0
  !From local aproximation: Diff_coef=2*C_dif*x+D_0
  real(8)::delta_C_Dif, delta_D_0 ! The step of modifcation of C_dif, D_0 in the
  !test loop
  real(8):: exact_FET, computed_FET, disp_FET! First mean exit time values: exact,
  ! computed, dispersion
  real(8)::diff_coef
  integer::Error_Flag, Error_Flag_FET ! zero, when OK, 102 when the intervals and internal points
  ! are not ordered, 103 when the diffusion coefficient is negative on the interval
  ! (a,b)
  real(8)::time_in,time_fin, cputime ! To measure speed of algorithm
  integer::ios !internal , used for input-output
  integer:: k_test, Nr_test !counter for test with different paramaters; number of tests
  character(len=10)::date , time, zone ! For timing, text output
  integer::values(8)  ! For timing, integer output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  character (len=*),parameter::Out_filename='Bessel0_5_FET_TESTresult.txt'

  real(8):: sum1, sum2,ww(2), ntr,xx, z1,z2,zz, sumz1, sumz2,relative_tolerance
  ! For testing random number generators, at the end
!!!!!!!!!!
  !INPUT DATA
!!!!!!!!! 
  Nr_test=2 ! NUmber of different values of parameters C_dif and  D_0
  ntraject=  1.0d3
  n_group= 10
  a=1.d0  
  ! a=a+1.0d-17
  !b=2.d0
  !b=a*(1.d0+1.0d-5)
  b=a*2.d0
  x_in=(a+b)/2.d0  ! 3.23456789d-8! Initial value of trajectory
  ! delta_t=1.000d-1 ;! max_duration=50;
  ! Initial test:max_duration>>first mean exit time
  C_dif = 1.1000d0; D_0=0.12d0 !   1.321d0 ! Initial values, changed in the following loop

  delta_C_dif=C_dif/3.0d0 ; delta_D_0=D_0/3.01d0 
  ! The increement of C_dif, D_0 in the N_test tests

  relative_tolerance=1.0d-3
  factor_max_duration=1.0d5 ! In this test is set very large, compared to meean
  !exit time. In such case the MC routine must give the correct value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !END INPUT DATA
  open(UNIT=3,file=Out_filename,iostat=ios,err=202,form='formatted',status='new',action='readwrite')
202 Print*, 'file opening status: IOS=',ios
601 format(8I5)
602 format(a30)
  write(UNIT=3, fmt=200)'TESTING Subroutines: RFOF_MeanExitTime, RFOF_Bessel_0_Move from module RFOF_BesselProcess0' 
  write(*, fmt=200)'TEST: Subs: Subs.RFOF_MeanExitTime, RFOF_Bessel_0_Move from module RFOF_BesselProcess0' 

  call date_and_time(date, time, zone, values)
  Write(*, fmt=602)' START TIME-file default unit' 
  Write(*, fmt=601)values
  Write(UNIT=3, fmt=602)' START TIME-file unit=3' 
  Write(UNIT=3, fmt=601)values
  write(UNIT=3, fmt=200)'INPUT DATA' 
  write(*, fmt=200)'INPUT DATA' 
200 format(a50)
  write(UNIT=3, fmt=201) 'Ntraject=', ntraject,' n_group=',n_group, ' start x=',x_in,&
       &' a=',a, '  b=',b 
  write(*, fmt=201) 'Ntraject=', ntraject,' n_group=',n_group,'  start x=',x_in,&
       &' a=',a, '  b=',b 
201 format(a12, d12.4,a12, I4,a12,d12.4,a6, d12.4,a6, d12.4)
300 format(2(a11,d14.7) )
  write(UNIT=3, fmt=204)'RESULTS: ' 
  write(*, fmt=204)'RESULTS' 
204 format(a20)
  CALL CPU_TIME ( time_in)

  !Loop to modify C_diff and D_0
  DO 1,  k_test=1, Nr_test
     write(UNIT=3, fmt=204)'  ' 
     write(*, fmt=204)' '
     ! Here we modify the input parameters
     C_dif=C_dif+delta_C_dif
     D_0=D_0+delta_D_0; ! Test negative and zero values too
     ! Test first the monotone dependencies 
     ! D_0=D_0 ! +1.2d0
     !diff_coef=2.d0*C_dif*x_in+D_0
     !If (diff_coef.LE.0) then
     !print*,' Wrong input, negative diffusion coefficient:', diff_coef
     !goto 1022
     ! endif
     write(UNIT=3,fmt=300)'C_dif=',C_dif, ' D_0=',D_0  
     write(*,fmt=300)'C_dif=',C_dif,' D_0=',D_0   
     ! Computation
     Error_Flag=0
     ! exact_FET, aprox_FET, disp_FET!
     ! Call the subroutine that give the exact value
     call RFOF_MeanExitTime1(a,x_in,b,relative_tolerance,C_dif, D_0,&
          &exact_FET, error_flag_FET)

     if(error_flag_FET.NE.0) Then
        print*,'Wrong input for RFOF_MeanExitTime2 subroutine'
        goto  1  ! No calculations, jump to next iteration
     endif
     max_duration=exact_FET*factor_max_duration
     print*, ' MAX DURATION=', max_duration
     ! Call   the local test subroutine, see below
     call Bessel0_FirstExitTime( ntraject, n_group,a,x_in,b,relative_tolerance,&
          max_duration,C_dif,D_0, computed_FET,disp_FET,Error_Flag)
     if(Error_Flag.NE.0 ) then
9       format(' err flag=',i4)
        write(*,9)Error_Flag 
        write(UNIT=3,fmt=9)Error_Flag
     endif ! Error flag if

     ! END COMPUTATION
     ! BEGIN WRITE RESULTS
111  format( 'Exact FET=',D11.4, ' computed=', D11.4, ' Stat. Error=',D11.4)
12   format(' CPUtime=',D11.4,'sec' )
     ! Write 

     if(Error_Flag.NE.0) then
        print*,'error in Bessel0_FirstExitTime() Time subroutine. ErrorFlag=',Error_Flag

     else
        write(*, 111)exact_FET, computed_FET, disp_FET
        write(UNIT=3, fmt=111)exact_FET, computed_FET, disp_FET
     endif
     !-------------------------------------------------
1022 Continue
1    CONTINUE ! End loop for k_test
     !END WRITE 
     ! BEGIN WRITE CPU-TIME
     CALL CPU_TIME ( time_fin)
     cputime=time_fin-time_in

     write(UNIT=3, fmt=204)'  ' 
     write(*, fmt=204)' ' 
     write(UNIT=3,fmt=206) ' CPU Time '
     write(*,fmt=206) ' CPU Time '
206  format(A20)
     write(UNIT=3, fmt=204)'  ' 
     write(*, fmt=204)' ' 
     write(*, 12)  cputime
     write(UNIT=3,fmt= 12)cputime
     !write(UNIT=3, fmt=205)'ERROR_FLAG=', Err_Flag
     !write(*, fmt=205)'ERROR_FLAG=', Err_Flag
205  format(a20,I4)
     call date_and_time(date, time, zone, values)
     Write(*, fmt=602)'FINAL TIME ' 
     Write(*, fmt=601)values
     Write(UNIT=3, fmt=602)'FINAL TIME' 
     Write(UNIT=3, fmt=601)values
     close(UNIT=3,iostat=ios,err=1001,  status='keep' )
1001 print *, 'File closed, iostat=',ios

1111 stop  ' End of test'
   end program  RFOF_BesselProcessTest4main
   !---------------------------------------------------------------------------------------
   !--------------------------------------------------------------------------------------
   !------------------------------------------------
   ! TEST PROGRAMS 
   !---------------------------------------------------------------------------

   subroutine Bessel0_FirstExitTime( ntraject, n_group,a,x_in,b,rel_tolerance,&
        &max_duration,C_dif,D_0, computed_FET,mean_square_err,Err_Flag)
     use RFOF_BesselProcess0
     use RFOF_random_numbers
     use RFOF_stochastics
     !Test the solver for the SDE, from module RFOF_BesselProcess0_4 , by computing mean values of the first exit time
     ! and comparing with the exact values
     ! Test of subroutine RFOF_multi_step_boundarycase for y_min=0,  
     implicit none
     intrinsic CPU_TIME !<used for speed test, and seed initialization in the random number generator
     !INPUT
     real(8), intent(in)::ntraject !>Number of trajectories in MC sampling
     integer, intent(in)::n_group ! Number of groups in MC sampling
     real(8), intent(in):: a, b  ! Interval end points
     real(8), intent(in) ::  x_in ! initial value of the variable x
     real(8), intent(in) ::rel_tolerance ! relative tolerance
     real(8), intent(in) :: max_duration ! Allowed maximal duration of the process
     ! inside (a,b) 
     real(8) ::count_traject !>counter of trajectories in MC sampling
     ! `in the tested Bessel procees, input for subroutine bessel0
     real(8), intent(in)::  C_dif,D_0! input in the tested Bessel procees, input for subroutine bessel0
     !From local aproximation: Diff_coef=2*C_dif*x +D_0+error
     !OUTPUT
     real(8), intent(out) :: computed_FET ! Computed value of the mean first exit time 
     real(8), intent(out) :: mean_square_err ! Computed value of the dispersion of the first exit time 
     integer,intent(out)::Err_Flag
     !The error signal Error_Flag=0 when OK.
     !LOCAL
     integer, parameter::nr_max_group=1000
     real(8):: loc_meanFET(nr_max_group) !local, Mean first exit time for each group real(8):: loc_mean_moments2(nr_moments) 
     ! Local variable used to compute the dispersion
     integer:: m,group , err_flag1! internal
     real(8)::time_begin,time_end, cputime 
     real(8)::   mean_square_err2 , m1,duration,x_fin !internal, for computing mean and the dispersion
     ! Fix the input values 
     integer::init_mode
     real(8)::init_parameter1, init_parameter2
     if( n_group>nr_max_group) then
        print *,' ERROR in main: n_group too large '
        return
     endif
     !Begin test 
     ! Initialise the seed of random number generators: 
     !call init_RANDOM_NUMBER(0) 
     init_mode=0
     init_parameter1=1.234d0
     init_parameter2=2.51d9
     CALL init_RM48(init_mode, init_parameter1, init_parameter2) 
     ! Initialization by using the intrinsic system_clock() subroutine
     CALL CPU_TIME ( time_begin )  
     ! Error_Flag=0
     DO 1001  group=1,n_group

        loc_meanFET(group)=0;   
        count_traject=0
        ! print*, 'GROUP=',group
        !Start of Monte-Carlo sampling 
        DO  ! Main loop, traject
           count_traject=count_traject+1.d0 
           Call RFOF_Bessel_0_move(a,x_in,b, max_duration,&
                &rel_tolerance,C_dif,D_0,duration,x_fin, err_flag1)
           Err_Flag=err_flag1

           loc_meanFET(group)=loc_meanFET(group)+duration 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! TO REMOVE
           !if(group.EQ.3) then
           !print*, 'duration=', duration, '  xfin=',x_fin
           !endif
           if((duration.GT.100.d0).OR.(duration.LE.0)) THEN
              print*, 'INNER LOOP, duration=', duration,' xfin=', x_fin
           endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           if(err_flag1.NE.0) then ! Exit from traject counting and group counting loop
              exit ! From traject counting loop
           endif

           if (count_traject.GE.ntraject)  exit
        ENDDO !!!!!!!!!!!!!!!!!!!!!!!!!!!! INNER MAIN LOOP 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! End MAIN traject counting loop 
        if (err_flag1.NE.0) then
           Print *, 'group=',group,  'Error_Flag=',err_flag1
           exit ! From group counting loop
        endif

1001    CONTINUE !!!!!!!!!!!!!!!!!!!!!!!!End group loop
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Outer main loop
        if(Err_Flag.NE.0) then
           !Anomalous ending
           print *,' ERROR SIGNAL=',Err_Flag
           return
        endif
        ! Monte-Carlo summation completed
        !Follows the processing of the results
        do group=1,n_group
           loc_meanFET(group)=loc_meanFET(group)/ntraject 
        enddo

        !computed_FET,disp_FET
        computed_FET=0
        do group=1,n_group
           computed_FET=computed_FET+ loc_meanFET(group)
        enddo
        computed_FET=computed_FET/n_group

        mean_square_err2=0
        do group=1,n_group
           mean_square_err2=mean_square_err2+ (loc_meanFET(group)-computed_FET)**2
        enddo

        mean_square_err2=mean_square_err2/n_group
        mean_square_err=sqrt(mean_square_err2)
        return
      end subroutine Bessel0_FirstExitTime
