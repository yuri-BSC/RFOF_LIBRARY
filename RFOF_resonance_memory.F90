#include "config.h"

module RFOF_resonance_memory

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  use RFOF_parameters
  use RFOF_numerics

  implicit none

  type, public :: resonance_memory
     integer :: Number_points_in_memory              !< Counter to keep track on the number of point 
                                                     !< there are in the resonance memory.
     real(8), dimension(:), allocatable :: time      !< Simulation time; (not the orbit time; for suggesting time step in orbit code)
                                                     !< (in which the wish to take derivatives)
     real(8), dimension(:), allocatable :: NACC      !< Rate of time acceleration
     real(8), dimension(:), allocatable :: r         !< Major radius
     real(8), dimension(:), allocatable :: phi       !< Toroidal angle
     real(8), dimension(:), allocatable :: z         !< Vertical position
     real(8), dimension(:), allocatable :: omega_res !< Resonance function \f$ \omega - n\Omega_c - \mathbf{k}\cdot\mathbf{v}\f$
     real(8), dimension(:), allocatable :: omega_c   !< Angular cyclotron frequency \f$ \Omega_c \f$
     integer, dimension(:), allocatable :: nharm     !< Number of the nearest harmonic resonance.
     logical :: previous_resonance_exists            !< TRUE is previous resonance has been found.
     real(8) :: time_of_last_resonance_crossing      !<  The simulation time of the previous resonance.
     integer :: sign_omega_dot_at_last_resonance_crossing = 0  !< The sign of the time derivative of the 
                                                               !< resonance function at the previous resonance
  end type resonance_memory

  type(resonance_memory), dimension(:,:), allocatable, save :: resonanceMemoryAllWaves

contains


  !--------------------------------------------------------------------------------
  ! subroutine constructor_rf_resonance_memory(mem,nStoreTimes)
  !
  subroutine constructor_rf_resonance_memory(mem,nStoreTimes)

    ! Input/Output
    type(resonance_memory), intent(inout) :: mem

    ! Input
    integer, intent(inout) :: nStoreTimes

    allocate(mem%time(nStoreTimes))
    allocate(mem%NACC(nStoreTimes))
    allocate(mem%R(nStoreTimes))
    allocate(mem%phi(nStoreTimes))
    allocate(mem%z(nStoreTimes))
    allocate(mem%omega_res(nStoreTimes))
    allocate(mem%omega_c(nStoreTimes))
    allocate(mem%nharm(nStoreTimes))

    mem%Number_points_in_memory = 0
    mem%previous_resonance_exists = .FALSE.

  end subroutine constructor_rf_resonance_memory


  !--------------------------------------------------------------------------------
  ! subroutine destructor_rf_resonance_memory(mem)
  !
  subroutine destructor_rf_resonance_memory(mem)

    ! Input/Output
    type(resonance_memory), intent(inout) :: mem

    integer :: DeAllocateStatus

    if ( allocated(mem%time) ) then 
       deallocate(mem%time, STAT = DeAllocateStatus)
       if (DeAllocateStatus /= 0) STOP "*** Trouble deallocating in subroutine destructor_rf_resonance_memory ***"
    endif
    if ( allocated(mem%NACC) ) then
       deallocate(mem%NACC, STAT = DeAllocateStatus)
       if (DeAllocateStatus /= 0) STOP "*** Trouble deallocating in subroutine destructor_rf_resonance_memory ***"
    endif

    if ( allocated(mem%R) ) then
       deallocate(mem%R, STAT = DeAllocateStatus)
       if (DeAllocateStatus /= 0) STOP "*** Trouble deallocating in subroutine destructor_rf_resonance_memory ***"
    endif
    if ( allocated(mem%phi) ) then
       deallocate(mem%phi, STAT = DeAllocateStatus)
       if (DeAllocateStatus /= 0) STOP "*** Trouble deallocating in subroutine destructor_rf_resonance_memory ***"
    endif
    if ( allocated(mem%z) ) then
       deallocate(mem%z, STAT = DeAllocateStatus)
       if (DeAllocateStatus /= 0) STOP "*** Trouble deallocating in subroutine destructor_rf_resonance_memory ***"
    endif

    if ( allocated(mem%omega_res) ) then
       deallocate(mem%omega_res, STAT = DeAllocateStatus)
       if (DeAllocateStatus /= 0) STOP "*** Trouble deallocating in subroutine destructor_rf_resonance_memory ***"
    endif
    if ( allocated(mem%omega_c) ) then
       deallocate(mem%omega_c, STAT = DeAllocateStatus)
       if (DeAllocateStatus /= 0) STOP "*** Trouble deallocating in subroutine destructor_rf_resonance_memory ***"
    endif
    if ( allocated(mem%nharm) ) then
       deallocate(mem%nharm, STAT = DeAllocateStatus)
       if (DeAllocateStatus /= 0) STOP "*** Trouble deallocating in subroutine destructor_rf_resonance_memory ***"
    endif

  end subroutine destructor_rf_resonance_memory


  !--------------------------------------------------------------------------------
  ! subroutine reset_rf_resonance_memory(mem)
  !
  subroutine reset_rf_resonance_memory(mem)

    ! Output/Input
    type(resonance_memory), intent(inout) :: mem

    mem%Number_points_in_memory = 0

  end subroutine reset_rf_resonance_memory


  !--------------------------------------------------------------------------------
  ! subroutine reset_rf_resonance_memory_matrix(mem)
  !
  subroutine reset_rf_resonance_memory_matrix(mem)

    ! Output/Input
    type(resonance_memory), pointer, intent(inout) :: mem(:,:)

    ! Local
    integer :: j1,j2

    do j1=1,size(mem,1)
       do j2=1,size(mem,2)
          call reset_rf_resonance_memory(mem(j1,j2))
       enddo
    enddo

  end subroutine reset_rf_resonance_memory_matrix


  !--------------------------------------------------------------------------------
  ! subroutine constructor_rf_resonance_memory_matrix(mem,nWaves1,nWaves2)
  !
  subroutine constructor_rf_resonance_memory_matrix(mem,nWaves1,nWaves2)

    ! Input
    integer, intent(in) :: nWaves1
    integer, intent(in) :: nWaves2

    ! Output
    type(resonance_memory), pointer, intent(out) :: mem(:,:)

    integer :: nStoreTimes
    integer :: jWave1,jWave2

    NAMELIST/input_resonance_memory/nStoreTimes

    open( io_channel_3872, FILE='input.rfof')
    read( io_channel_3872, input_resonance_memory)
    close(io_channel_3872)

    allocate(mem(nWaves1,nWaves2))

    do jWave1=1,size(mem,1)
       do jWave2=1,size(mem,2)
          call constructor_rf_resonance_memory(mem(jWave1,jWave2) , nStoreTimes)
       end do
    end do

  end subroutine constructor_rf_resonance_memory_matrix

  !--------------------------------------------------------------------------------
  ! subroutine destructor_rf_resonance_memory_matrix(mem)
  !
  subroutine destructor_rf_resonance_memory_matrix(mem)

    type(resonance_memory), pointer, intent(inout) :: mem(:,:)

    integer :: DeAllocateStatus
    integer :: jWave1,jWave2

    if ( .not. associated(mem) ) then
       return
    endif

    do jWave1=1,size(mem,1)
       do jWave2=1,size(mem,2)
          call destructor_rf_resonance_memory(mem(jWave1,jWave2))
       end do
    end do

    deallocate(mem, STAT = DeAllocateStatus)
    if (DeAllocateStatus /= 0) STOP "*** Trouble deallocating in subroutine destructor_rf_resonance_memory_matrix ***"

  end subroutine destructor_rf_resonance_memory_matrix

  !--------------------------------------------------------------------------------
  ! subroutine storeNewPointResonanceMemorySort
  subroutine storeNewPointResonanceMemorySort(mem,time_at_resonance,time,nharm,NACC,R,phi,z, &
       omega_res,omega_c,was_stored)

    ! Input
    real(8)   , intent(in) :: time_at_resonance
    real(8)   , intent(in) :: time
    integer   , intent(in) :: nharm
    real(8)   , intent(in) :: NACC
    real(8)   , intent(in) :: R
    real(8)   , intent(in) :: z
    real(8)   , intent(in) :: phi
    real(8)   , intent(in) :: omega_res
    real(8)   , intent(in) :: omega_c

    ! Output
    logical, intent(out) :: was_stored

    ! Input/Output
    type(resonance_memory), intent(inout) :: mem

    ! Local
    integer Nmem, Nstored, j, jcnt

    was_stored=.TRUE.

    Nmem=size(mem%time)
    Nstored=mem%Number_points_in_memory ! Note that mem%Number_points_in_memory changes later in the routine,
    ! thus the loop below is written in Nstored, which remains unchanged.

    do jcnt=1,Nstored

       j=Nstored-jcnt+1

       if ( time .lt. mem%time(j) ) then
          if (j .lt. Nmem) then
             call storeNewPointResonanceMemoryLocationJ(mem,j+1,time,R,phi,z,omega_res,omega_c,nharm,NACC)
          else
             was_stored=.FALSE.
          endif
          return
       else
          if (j .lt. Nmem) then
             call storeNewPointResonanceMemoryLocationJ(mem,j+1,&
                  mem%time(j),mem%R(j),mem%phi(j),mem%z(j),&
                  mem%omega_res(j),mem%omega_c(j),mem%nharm(j),mem%NACC(j))
          endif
       endif
    enddo

    ! In case the new data is the latest in time, or mem%Number_points_in_memory=0:
    j=1
    call storeNewPointResonanceMemoryLocationJ(mem,j,time,R,phi,z,omega_res,omega_c,nharm,NACC)

  end subroutine storeNewPointResonanceMemorySort


  !--------------------------------------------------------------------------------
  ! subroutine storeNewPointResonanceMemory
  subroutine storeNewPointResonanceMemory(mem,time_at_resonance,time,nharm,NACC,R,phi,z, &
       omega_res,omega_c,was_stored)

    ! Input
    real(8)   , intent(in) :: time_at_resonance
    real(8)   , intent(in) :: time
    integer   , intent(in) :: nharm
    real(8)   , intent(in) :: NACC
    real(8)   , intent(in) :: R
    real(8)   , intent(in) :: z
    real(8)   , intent(in) :: phi
    real(8)   , intent(in) :: omega_res
    real(8)   , intent(in) :: omega_c

    ! Output
    logical, intent(out) :: was_stored

    ! Input/Output
    type(resonance_memory), intent(inout) :: mem

    ! Local
    integer Nmem, Nstored, j, jcnt

    was_stored=.TRUE.

    Nmem=size(mem%time)
    Nstored=mem%Number_points_in_memory ! Note that mem%Number_points_in_memory changes later in the routine,
    ! thus the loop below is written in Nstored, which remains unchanged.

    ! Move all points all ready in the list one position to leave the first position empty; to be fill in with the new point
    do jcnt=1,Nstored

       j=Nstored-jcnt+1

       if (j .lt. Nmem) then
          call storeNewPointResonanceMemoryLocationJ(mem,j+1,&
               mem%time(j),mem%R(j),mem%phi(j),mem%z(j),&
               mem%omega_res(j),mem%omega_c(j),mem%nharm(j),mem%NACC(j))
       endif
    enddo

    ! In case the new data is the latest in time, or mem%Number_points_in_memory=0:
    j=1
    call storeNewPointResonanceMemoryLocationJ(mem,j,time,R,phi,z,omega_res,omega_c,nharm,NACC)

  end subroutine storeNewPointResonanceMemory


  !--------------------------------------------------------------------------------
  subroutine storeNewPointResonanceMemoryLocationJ(mem,j,time,R,phi,z, &
       omega_res,omega_c,nharm,NACC)

    ! Input
    integer, intent(in) :: j
    real(8), intent(in) :: time
    real(8), intent(in) :: R
    real(8), intent(in) :: phi
    real(8), intent(in) :: z
    real(8), intent(in) :: omega_res
    real(8), intent(in) :: omega_c
    integer, intent(in) :: nharm
    real(8), intent(in) :: NACC

    ! Input/Output
    type(resonance_memory), intent(inout) :: mem

    mem%time(j)      = time
    mem%R(j)         = R
    mem%phi(j)       = phi
    mem%z(j)         = z
    mem%omega_res(j) = omega_res
    mem%omega_c(j)   = omega_c
    mem%nharm(j)     = nharm
    mem%NACC(j)      = NACC

    if (j .gt. mem%Number_points_in_memory) then
       mem%Number_points_in_memory = j
    endif

  end subroutine storeNewPointResonanceMemoryLocationJ


  !--------------------------------------------------------------------------------
  ! subroutine estimate_resonance_location
  !
  subroutine estimate_resonance_location(mem,time,R,phi,z,omega_res,nharm,NACC, &
       time_previous_step, &
       time_at_resonance,R_at_resonance,phi_at_resonance,z_at_resonance, &
       Dot_omega_res,DotDot_omega_res,MPI_node_Id,errorFlag)

    !    use RFOF_diagnostics

    ! Input
    type(resonance_memory), intent(in) :: mem
    real(8), intent(in) :: time
    real(8), intent(in) :: R
    real(8), intent(in) :: phi
    real(8), intent(in) :: z
    real(8), intent(in) :: omega_res
    integer, intent(in) :: nharm
    real(8), intent(in) :: NACC
    real(8), intent(in) :: time_previous_step
    integer, intent(in) :: MPI_node_Id

    ! Output
    real(8), intent(out) :: time_at_resonance
    real(8), intent(out) :: R_at_resonance
    real(8), intent(out) :: phi_at_resonance
    real(8), intent(out) :: z_at_resonance
    real(8), intent(out) :: Dot_omega_res
    real(8), intent(out) :: DotDot_omega_res
    integer, intent(out) :: errorFlag

    ! Local
    real(8) :: tACC(3),tORB(3)
    real(8) :: NACCvec(3)
    real(8) :: x(3), f(3), coef(3), tnorm
    real(8) :: root1,root2,xroot,time_root1,time_root2
    integer :: Nsolutions

    interface
       subroutine wrap_save2file_resonace_predictions(time,R,phi,z,omega_res,nharm, &
            time_at_resonance,R_at_resonance,phi_at_resonance,z_at_resonance, &
            Dot_omega_res,DotDot_omega_res,MPI_node_Id)

         ! Input
         real(8),    intent(in) :: time
         real(8),    intent(in) :: R
         real(8),    intent(in) :: phi
         real(8),    intent(in) :: z
         real(8),    intent(in) :: omega_res
         integer, intent(in) :: nharm
         integer, intent(in) :: MPI_node_Id

         ! Output
         real(8), intent(out) :: time_at_resonance
         real(8), intent(out) :: R_at_resonance
         real(8), intent(out) :: phi_at_resonance
         real(8), intent(out) :: z_at_resonance
         real(8), intent(out) :: Dot_omega_res
         real(8), intent(out) :: DotDot_omega_res

       end subroutine wrap_save2file_resonace_predictions
    end interface

    if ( mem%Number_points_in_memory .lt. 2 ) then
       errorFlag = 1
       return
    endif
    errorFlag=0

    ! Generate time vector without acceleration
    tACC(1)      = time
    tACC(2:3)    = mem%time(1:2)
    NACCvec(1)   = NACC
    NACCvec(2:3) = mem%NACC(1:2)
    call calc_non_accelerated_orbit_time_vec(tACC,NACCvec,tORB)

    !write(0,*)"acc-test",tORB(1)-tORB(2),tORB(2)-tORB(3),NACC

    ! Before matching to polynomial; fill the vectors x and f. 
    ! Note that the time x is normalised by tnorm for better convergence.
    tnorm = abs(tORB(1)) + abs(tORB(2)) + abs(tORB(3)) + 1e-12

    x(1)=0.
    x(2)=(tORB(2) - tORB(1)) / tnorm
    x(3)=(tORB(3) - tORB(1)) / tnorm

    f(1)=omega_res
    f(2)=mem%omega_res(1) + real(mem%nharm(1) - nharm)*mem%omega_c(1)
    f(3)=mem%omega_res(2) + real(mem%nharm(2) - nharm)*mem%omega_c(2)

    !    if (R.gt.3d0 .and. R.lt.3.3d0 .and. time.gt.6.9e-4) then
    !       write(0,'(A,200E15.5)')'estres:W',f,mem%omega_c(1:2)
    !    endif

    ! Find the polynomial coefficients and solve the quadratic equation
    call polycoef_3points(x,f,coef,errorFlag)
    if (errorFlag .ne. 0) return
    call solve_quadratic_polynomial_real_roots(coef,root1,root2,Nsolutions)

    if (Nsolutions .eq. 0) then
       errorFlag = 1
       return
    endif

    time_root1 = root1*tnorm*NACC + time
    time_root2 = root2*tnorm*NACC + time

    ! Remove root1 if it occures before time_previous_step
    if ( time_root1 .lt. time_previous_step ) then

       ! If also root2 occure before the time_previous_step, then 
       ! no valid resonance found
       if ( time_root2 .lt. time_previous_step  ) then
          errorFlag = 1
          return
       endif
       ! Root1 is removed by overwriting it with root2 data 
       ! (makes the logic simple; as if root1 was valid)
       root1 = root2
       time_root1 = time_root2
    endif

    ! Remove root2 if it occures before time_previous_step
    if ( time_root2 .lt. time_previous_step ) then
       ! Root2 is removed by overwriting it with root1 data 
       ! (makes the logic simple; as if root2 was valid)
       root2 = root1
       time_root2 = time_root1
    endif

    ! Choose the one of the two roots which is closed to the present time.
    ! Note that x was defined such that present time is x=0.
    if (abs(root1) .gt. abs(root2)) then
       xroot = root2
       time_at_resonance = time_root2
    else
       xroot = root1
       time_at_resonance = time_root1
    endif

    ! time_at_resonance - estimated times of the resonance. 
    ! Note that this time is an accelerated time (not orbit time), 
    ! thus the multiplication with NACC.
    Dot_omega_res = coef(2)
    DotDot_omega_res = 2. * coef(3)

    ! Find R_at_resonance
    f(1)=R
    f(2)=mem%R(1)
    f(3)=mem%R(2)
    call polycoef_3points(x,f,coef,errorFlag)
    if (errorFlag .ne. 0) return
    R_at_resonance = eval_quadratic_polynomial(coef,xroot)


    !if (R.gt.3d0 .and. R.lt.3.3d0 .and. time.gt.6.9e-4) then
    !   write(0,'(A,200E15.5)')'estres:t', x
    !write(0,'(A,200E15.5)')'estres:R',R,R_at_resonance,time,time_previous_step,time_at_resonance
    !   write(0,'(A,200E15.5)')'estres:Rc', coef
    !endif
    !   print *, "t=", x, time_at_resonance
    ! print *, "RES-MEM: R=", f, R_at_resonance,root1,root2
    !   print *, "coef", coef
    !   R_at_resonance = eval_quadratic_polynomial(coef,root2)
    !   print *, "2. R=", f, R_at_resonance
    !write(0,'(A,200E16.6)') 'estRres',mem%R,R_at_resonance,time_at_resonance

    ! Find phi_at_resonance
    f(1)=phi
    f(2)=mem%phi(1)
    f(3)=mem%phi(2)
    call polycoef_3points(x,f,coef,errorFlag)
    if (errorFlag .ne. 0) return
    phi_at_resonance = eval_quadratic_polynomial(coef,xroot)

    ! Find z_at_resonance
    f(1)=z
    f(2)=mem%z(1)
    f(3)=mem%z(2)
    call polycoef_3points(x,f,coef,errorFlag)
    if (errorFlag .ne. 0) return
    z_at_resonance = eval_quadratic_polynomial(coef,xroot)

    if (output__resonace_predictions) then
       call wrap_save2file_resonace_predictions(time,R,phi,z,omega_res,nharm, &
            time_at_resonance,R_at_resonance,phi_at_resonance,z_at_resonance, &
            Dot_omega_res,DotDot_omega_res,MPI_node_Id)
    endif

  end subroutine estimate_resonance_location



  !--------------------------------------------------------------------------------
  ! subroutine estimate_resonance_location
  !
  subroutine calc_non_accelerated_orbit_time_vec(timeACC,NACC,timeORBIT)

    ! Input
    real(8), intent(in) :: timeACC(:)
    real(8), intent(in) :: NACC(:)

    ! Output
    real(8), intent(inout) :: timeORBIT(:)

    ! Local
    integer :: Ntimes, j

    Ntimes = size(NACC)

    timeORBIT(1) = timeACC(1)
    do j = 2, Ntimes
       timeORBIT(j) = timeORBIT(j-1) + (timeACC(j)-timeACC(j-1)) / NACC(j-1)
    enddo

  end subroutine calc_non_accelerated_orbit_time_vec


end module RFOF_resonance_memory
