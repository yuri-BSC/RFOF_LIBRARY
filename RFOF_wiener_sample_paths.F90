#include "config.h"

module RFOF_wiener_sample_paths

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  implicit none

  type :: Wiener_process_sample_path
     real(8), pointer :: time(:)    !< Times [s]
     real(8), pointer :: Wiener(:)  !< Value of Wiener process [\f$s^{1/2}\f$]
     integer :: Nr_stored_points    !< Number of values stored in the sample path
  end type Wiener_process_sample_path

contains

  !--------------------------------------------------------------------------------
  subroutine constructor_Wiener_process_sample_path(samplePath,Nelements, &
       time1,time2,rand_normal_distribution)

    ! Input
    type(Wiener_process_sample_path), intent(inout) :: samplePath
    integer, intent(in) :: Nelements
    real(8), intent(in) :: time1
    real(8), intent(in) :: time2
    real(8), intent(in) :: rand_normal_distribution

    if (Nelements .lt. 2) then
       write(0,*)"Error in constructor_Wiener_process_sample_path:"
       write(0,*)"Attempting to allocate a sample path with less than 2 point"
       write(0,*)"ABORT"
       stop
    endif

    if (time1 .gt. time2) then
       write(0,*)"Error in constructor_Wiener_process_sample_path:"
       write(0,*)"Invalid order of time points."
       write(0,*)"ABORT"
       stop
    endif

    allocate(samplePath%time(Nelements))
    allocate(samplePath%Wiener(Nelements))

    samplePath%time(2)   = time1
    samplePath%Wiener(2) = 0
    samplePath%time(1)   = time2
    samplePath%Wiener(1) = sqrt(time2-time1) * rand_normal_distribution
    samplePath%Nr_stored_points = 2

  end subroutine constructor_Wiener_process_sample_path

  !--------------------------------------------------------------------------------
  subroutine destructor_Wiener_process_sample_path(samplePath)

    ! Input
    type(Wiener_process_sample_path), intent(inout) :: samplePath
    deallocate(samplePath%time)
    deallocate(samplePath%Wiener)

  end subroutine destructor_Wiener_process_sample_path

  !--------------------------------------------------------------------------------
  subroutine add_point_on_sample_path(samplePath, new_time, rand_normal_distribution)

    ! Input
    type(Wiener_process_sample_path), intent(inout) :: samplePath
    real(8), intent(in) :: new_time
    real(8), intent(in) :: rand_normal_distribution

    ! Local
    integer :: j,k

    if (new_time .lt. samplePath%time(samplePath%Nr_stored_points)) then
       write(0,*)"Error in add_point_on_sample_path:"
       write(0,*)"Attempting to add point before the start of the Wiener process"
       write(0,*)"ABORT"
       stop
    endif

    ! Make sure there's memory allocated to add one extra point on sample path
    if (samplePath%Nr_stored_points .eq. size(samplePath%time)) then
       call extendVectorAllocation(samplePath%time  ,samplePath%Nr_stored_points)
       call extendVectorAllocation(samplePath%Wiener,samplePath%Nr_stored_points)
    endif

    ! Find where in the list to add new point
    do j=1,samplePath%Nr_stored_points-1

       ! Reversed counter; counting backwards from samplePath%Nr_stored_points to 2
       k = samplePath%Nr_stored_points - j + 1

       ! Check if new point sort in between k and k+1 
       ! (note that samplePath%time(k) < samplePath%time(k-1))
       if ((new_time .ge. samplePath%time(k)) .and. &
            ( new_time .lt. samplePath%time(k-1)) ) then
          call add_point_on_sample_path_index_k(k, samplePath, new_time, &
               rand_normal_distribution)
          return
       endif
    enddo

    ! Add point at end of list
    call add_point_on_sample_path_index_k(samplePath%Nr_stored_points+1, &
         samplePath, new_time, rand_normal_distribution)

  end subroutine add_point_on_sample_path

  !--------------------------------------------------------------------------------
  subroutine add_point_on_sample_path_index_k(k, samplePath, new_time, &
       rand_normal_distribution)

    ! Input
    integer, intent(in) :: k
    type(Wiener_process_sample_path), intent(inout) :: samplePath
    real(8), intent(in) :: new_time
    real(8), intent(in) :: rand_normal_distribution

    ! Local
    real(8) :: newWiener

    if (k .le. samplePath%Nr_stored_points) then

       newWiener = Brownian_bridge( &
            samplePath%time(k),samplePath%time(k-1), &
            samplePath%Wiener(k),samplePath%Wiener(k-1), &
            new_time,rand_normal_distribution)

       ! Move points to make room to insert the new point
       samplePath%time(k+1:samplePath%Nr_stored_points+1) = &
            samplePath%time(k:samplePath%Nr_stored_points)
       samplePath%Wiener(k+1:samplePath%Nr_stored_points+1) = &
            samplePath%Wiener(k:samplePath%Nr_stored_points)

       ! Insert new point
       samplePath%time(k) = new_time
       samplePath%Wiener(k) = newWiener

    else

       samplePath%time(k) = new_time
       samplePath%Wiener(k) = sqrt(new_time - samplePath%time(k-1)) * &
            rand_normal_distribution

    endif

    samplePath%Nr_stored_points = samplePath%Nr_stored_points + 1

    write(0,*)"new samplePath%Nr_stored_points:",samplePath%Nr_stored_points

  end subroutine add_point_on_sample_path_index_k

  !--------------------------------------------------------------------------------
  subroutine remove_N_points_from_sample_path(samplePath,Nr_points_to_remove)

    type(Wiener_process_sample_path), intent(inout) :: samplePath
    integer, intent(in) :: Nr_points_to_remove

    samplePath%Nr_stored_points = samplePath%Nr_stored_points - Nr_points_to_remove

  end subroutine remove_N_points_from_sample_path

  !--------------------------------------------------------------------------------
  function Brownian_bridge(t0,t1,W0,W1,tn,rand_normal_distribution) result(Wn)

    ! Input
    real(8), intent(in) :: t0
    real(8), intent(in) :: t1
    real(8), intent(in) :: W0
    real(8), intent(in) :: W1
    real(8), intent(in) :: tn
    real(8), intent(in) :: rand_normal_distribution

    ! Output
    real(8) :: Wn

!    Wn=W0 + ( (W1-W0) * tn &
!         +   rand_normal_distribution * (tn-t0) * (t1-tn) ) / (t1-t0)
    Wn=W0 + (W1-W0) * (tn-t0) / (t1-t0) &
         + rand_normal_distribution * (tn-t0) * (t1-tn) / (t1-t0)

    !write(0,*)"Brownian_bridge",t0,tn,t1,w0,wn,w1

  end function Brownian_bridge

  !--------------------------------------------------------------------------------
  subroutine extendVectorAllocation(v,Nextra)

    ! Input
    real(8), pointer, intent(inout) :: v(:)
    integer, intent(in) :: Nextra

    ! Local
    real(8), allocatable, target :: w(:)
    integer :: N, j
    integer :: DeAllocateStatus

    N=size(v)

    allocate(w(size(v)+Nextra))
    do j=1,N
       w(j)=v(j)
    end do

    deallocate(v, STAT = DeAllocateStatus)
    if (DeAllocateStatus /= 0) STOP "*** Trouble deallocating in subroutine extendVector ***"

    allocate(v(size(w)))
    v = w

  end subroutine extendVectorAllocation

  !--------------------------------------------------------------------------------
  subroutine test_sample_path

    ! Local
    type(Wiener_process_sample_path) :: samplePath
    integer :: Nelements
    real(8) :: tstart,tend,t3,t4,t5,t6,rend,r3,r4,r5,r6
    real(8) :: dt,dW

    tstart = 1d0
    tend = 2d0
    t3 = 1.4d0
    t4 = 1.8d0
    t5 = 1.2d0
    t6 = 1.7d0

    rend = 1d0
    r3 = -1.1d0
    r4 = 0.74d0
    r5 = -0.24d0
    r6 = -0.24d0

    Nelements = 3

    call constructor_Wiener_process_sample_path(samplePath,Nelements,tstart,tend,rend)

    write(0,*)"path after constructor",samplePath%Nr_stored_points,size(samplePath%time),size(samplePath%Wiener)
    write(0,*)samplePath%time(1:samplePath%Nr_stored_points)
    write(0,*)samplePath%Wiener(1:samplePath%Nr_stored_points)

    ! Differentials
    dt = samplePath%time(  samplePath%Nr_stored_points-1) - samplePath%time(  samplePath%Nr_stored_points)
    dW = samplePath%Wiener(samplePath%Nr_stored_points-1) - samplePath%Wiener(samplePath%Nr_stored_points)

    call add_point_on_sample_path(samplePath, t3, r3)
    write(0,*)"1 point added",samplePath%Nr_stored_points,size(samplePath%time),size(samplePath%Wiener)
    write(0,*)samplePath%time(1:samplePath%Nr_stored_points)
    write(0,*)samplePath%Wiener(1:samplePath%Nr_stored_points)

    call add_point_on_sample_path(samplePath, t4, r4)
    write(0,*)"2 point added",samplePath%Nr_stored_points,size(samplePath%time),size(samplePath%Wiener)
    write(0,*)samplePath%time(1:samplePath%Nr_stored_points)
    write(0,*)samplePath%Wiener(1:samplePath%Nr_stored_points)

    call add_point_on_sample_path(samplePath, t4, r4)
    write(0,*)"3 point added",samplePath%Nr_stored_points,size(samplePath%time),size(samplePath%Wiener)
    write(0,*)samplePath%time(1:samplePath%Nr_stored_points)
    write(0,*)samplePath%Wiener(1:samplePath%Nr_stored_points)

    ! Once a part of the stepping has been completed; run this routine remove_N_points_from_sample_path
    ! and the then "start" of the sample path is move forward
    call remove_N_points_from_sample_path(samplePath,1)

    call add_point_on_sample_path(samplePath, t4, r4)
    write(0,*)"3 point added",samplePath%Nr_stored_points,size(samplePath%time),size(samplePath%Wiener)
    write(0,*)samplePath%time(1:samplePath%Nr_stored_points)
    write(0,*)samplePath%Wiener(1:samplePath%Nr_stored_points)

    call add_point_on_sample_path(samplePath, t5, r5) ! This should make the program stop since t5 is before t3, which is the "start time" (after we've removed the initial start point)
    write(0,*)"4 point added",samplePath%Nr_stored_points,size(samplePath%time),size(samplePath%Wiener)
    write(0,*)samplePath%time(1:samplePath%Nr_stored_points)
    write(0,*)samplePath%Wiener(1:samplePath%Nr_stored_points)

    call destructor_Wiener_process_sample_path(samplePath)

  end subroutine test_sample_path

end module RFOF_wiener_sample_paths
