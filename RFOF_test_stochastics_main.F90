!TEST PROGRAM


Program RFOF_test_stochastics_main
    use RFOF_stochastics
    implicit none
    integer(4)::init_mode 
    real(8):: nr_iterations, init_parameter1, init_parameter2

    nr_iterations=1.0d7
    ! Iteration number must be larger than 1.0d5
    init_parameter2=2.11d0
    ! init_parameter1=0.1d0

    print*,' TEST  INITIALIZATION OF RNADOM NUMBER GENERATOR'
    !call init_RM48(0, init_parameter1,init_parameter2)

    do init_mode=0,1
        init_parameter1=0.1001111d0
        do  while (init_parameter1<1.d0)
            print*,''
            print*,''
            print*,''
            print*, 'init_mode=', init_mode,'init_parameter1= ', init_parameter1,&
            & ' init_paramter2=', init_parameter2

            !------------------------------------------------------------
            ! For more visibility the subsequent tests are better to perform by 
            ! de- comment and commnent of the corresponding groups of instructions
            !-----------------------------------------------------------
            !------------------------------------------------
            ! Group 1 for test
            ! Remove comment from this group, for separate tests
            ! This version is used in BesselProcess generators
            print*,'  PolarGaussAccuracyTest1, testing the correct initialization of RND generators' 
            Call  PolarGaussAccuracyTest1(nr_iterations, init_mode,&
            &init_parameter1, init_parameter2)
            if(init_mode==0) then
                print*, ' seed shifted by clock; results must be different from ::: '
            else
                print*,'seed fixed by init_parameter only, prev results must be equal to '
            endif

            Call  PolarGaussAccuracyTest1(nr_iterations, init_mode,&
            &init_parameter1, init_parameter2)
            ! End group 1 test
            !----------------------------------------------------------------
            
            !-----------------------------------------------------
            ! Group 2 for test
            ! Remove comment from this group, for separate tests
            !  print*,'  PolarGaussAccuracyTest2, testing the correct initialization of RND generators' 
            ! Call  PolarGaussAccuracyTest2(nr_iterations, init_mode,&
            ! &init_parameter1,init_parameter2)
            ! if(init_mode==0) then
            !     print*, ' seed shifted by clock; results must be different from ::: '
            ! else
            !    print*,'seed fixed by init_parameter only, prev results must be equal to '
            ! endif
            ! Call  PolarGaussAccuracyTest2(nr_iterations, init_mode,&
            ! &init_parameter1,init_parameter2)
            ! End test group 2
            !-----------------------------------------------------
            
            !-----------------------------------------------------
            ! Group 3 for test
            ! Remove comment from this group, for separate tests
            ! print*,'  MarsagliaBrayAccuracyTest1, testing the correct initialization of RND generators' 
            ! Call  MarsagliaBrayAccuracyTest1(nr_iterations, init_mode,&
            ! &init_parameter1,init_parameter2)
            ! if(init_mode==0) then
            !     print*, ' seed shifted by clock; results must be different from ::: '
            ! else
            !   print*,'seed fixed by init_parameter only, prev results must be equal to '
            ! endif
            ! Call    MarsagliaBrayAccuracyTest1(nr_iterations, init_mode,&
            ! &init_parameter1,init_parameter2)
            ! End test group 3
            !-----------------------------------------------------


            !-----------------------------------------------------
            ! Group 4 test
            ! Remove comment from this group, for separate tests
            ! print*,'  MarsagliaBrayAccuracyTest2, testing the correct initialization of RND generators' 
            ! Call  MarsagliaBrayAccuracyTest2(nr_iterations, init_mode,&
            ! &init_parameter1,init_parameter2)
            ! Call  MarsagliaBrayAccuracyTest2(nr_iterations, init_mode,&
            ! &init_parameter1,init_parameter2)

            ! End test group 4
            !-----------------------------------------------------

            !-----------------------------------------------------
            ! Group 5 test
            ! Remove comment from this group, for separate tests
            ! print*,'  CorrelatedGaussTest, testing the correct initialization of RND generators' 
            ! Call  Correlated_GaussTest(nr_iterations, init_mode,&

            ! &init_parameter1,init_parameter2)
            ! if(init_mode==0) then
            !     print*, ' seed shifted by clock; results must be different from ::: '
            ! else
            !    print*,'seed fixed by init_parameter only, prev results must be equal to '
            ! endif
            ! Call   Correlated_GaussTest(nr_iterations, init_mode,&
            ! &init_parameter1,init_parameter2)
            ! End test group 5
            !-----------------------------------------------------



            init_parameter1=init_parameter1+0.5d0
        enddo

        print*,''
        print*,''
    enddo

    Stop 'End Test stochastics'
END Program RFOF_test_stochastics_main
