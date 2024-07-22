!!! A module that allows to use the same source for both serial and parallel 
!!! versions. Serial code can be compiled without linking to mpi or having 
!!! a mpi compiler. 

#include "config.h"

module RFOF_mpi_module

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  IMPLICIT NONE
  !INCLUDE 'mpif.h' !not used in serial version!
  INTEGER, PARAMETER ::  MPI_SUCCESS=0

!!! RFOF variables:
  INTEGER, PARAMETER :: R4 = SELECTED_REAL_KIND(6,37) ! Real*4
  INTEGER, PARAMETER :: R8 = SELECTED_REAL_KIND(15,300) ! Real*8
  INTEGER, PARAMETER :: params_wp = R8
  INTEGER, save :: mpivar_id
  INTEGER, save :: mpivar_numproc
  LOGICAL, save :: mpivar_independentParallelJob
!!! RFOF variables.


!!!*******************************************************************
!!!
!!! Dummy routine for mpi init
!!!
!!!*******************************************************************
contains
  
  subroutine my_mpi_init(ierr)
    !use mpivar
    !use parseArgs
    implicit none
    integer :: ierr
    character(len=30) :: argString,numstr1,numstr2

!!$    ! Check if we are doing independent parallel processing.
!!$    ! Many CPUs, but no MPI to take care of parallelising.
!!$    if (parseArguments(argString, mpivar_numproc, mpivar_id)) then
!!$       IF (INDEX(argString,'--only-process')==1) THEN
!!$          WRITE(numstr1,*) mpivar_id+1
!!$       ELSE
!!$          WRITE(numstr1,*) mpivar_id
!!$       END IF
!!$       WRITE(numstr2,*) mpivar_numproc
!!$       WRITE(*,*)
!!$       WRITE(*,*) 'NOTE: single-process command argument '// &
!!$            TRIM(argString)//'...'
!!$       WRITE(*,*) '      given. Using process id and number of '// &
!!$            'processors'
!!$       WRITE(*,*) '      to simulate specific parallel process '// &
!!$       TRIM(ADJUSTL(numstr1))//'/'//TRIM(ADJUSTL(numstr2))//'.'
!!$       mpivar_independentParallelJob = .true.
!!$    else

       ! Well we are not doing anything in parallel!
       mpivar_id=0
       mpivar_numproc=1
       mpivar_independentParallelJob = .false.

!!$    end if

    ierr=MPI_SUCCESS
  end subroutine my_mpi_init

!!!*******************************************************************
!!!
!!! Dummy routine for mpi finalize
!!!
!!!*******************************************************************

  subroutine my_mpi_finalize(ierr)
    implicit none
    integer, intent(out) :: ierr
    ierr=MPI_SUCCESS
    return
  end subroutine my_mpi_finalize

!!!*******************************************************************
!!!
!!! Dummy routine for mpi barrier
!!!
!!!*******************************************************************

  SUBroutine my_mpi_barrier(ierr)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: ierr
    ierr=MPI_SUCCESS
    return
  end SUBroutine my_mpi_barrier

!!!*******************************************************************
!!!
!!! Dummy routine for 0D mpi bcast (INTEGER)
!!!
!!!*******************************************************************

  SUBroutine my_mpi_bcast0d_int(data,nsize,bcaster,ierr)
    implicit none
    integer, intent(in)  :: data
    integer, intent(in) :: nsize, bcaster
    integer, intent(out) :: ierr
    ierr=MPI_SUCCESS

    ! Shut up the compiler 
    if (nsize==0 .or. data==0 .or. bcaster==0) then 
    end if

    return
  end SUBroutine my_mpi_bcast0d_int

!!!*******************************************************************
!!!
!!! Dummy routine for 1D mpi bcast (INTEGER)
!!!
!!!*******************************************************************

  SUBroutine my_mpi_bcast1d_int(data,nsize,bcaster,ierr)
    implicit none
    integer, intent(in), dimension(:) :: data
    integer, intent(in) :: nsize, bcaster
    integer, intent(out) :: ierr
    ierr=MPI_SUCCESS
    ! Shut up the compiler 
    if (nsize==0 .or. data(1)==0 .or. bcaster==0) then 
    end if
    return
  end SUBroutine my_mpi_bcast1d_int

!!!*******************************************************************
!!!
!!! Dummy routine for 2D mpi bcast (INTEGER)
!!!
!!!*******************************************************************

  SUBroutine my_mpi_bcast2d_int(data,nsize,bcaster,ierr)
    implicit none
    integer, intent(in), dimension(:,:) :: data
    integer, intent(in) :: nsize, bcaster
    integer, intent(out) :: ierr
    ierr=MPI_SUCCESS
    ! Shut up the compiler 
    if (nsize==0 .or. data(1,1)==0 .or. bcaster==0) then 
    end if
    return
  end SUBroutine my_mpi_bcast2d_int
  

!!!*******************************************************************
!!!
!!! Dummy routine for mpi bcast (working precision)
!!!
!!!*******************************************************************

  SUBroutine my_mpi_bcast_wp(data,nsize,bcaster,ierr)
    !use params, only : params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:) :: data
    integer, intent(in) :: nsize, bcaster
    integer, intent(out) :: ierr
    ! Shut up the compiler 
    if (nsize==0 .or. data(1)==0 .or. bcaster==0) then 
    end if
    ierr=MPI_SUCCESS
    return
  end SUBroutine my_mpi_bcast_wp

!!!*******************************************************************
!!!
!!! Dummy routine for mpi bcast (working precision)
!!!
!!!*******************************************************************

  SUBroutine my_mpi_bcast2d_wp(data,nsize,bcaster,ierr)
    !use params, only : params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:,:) :: data
    integer, intent(in) :: nsize, bcaster
    integer, intent(out) :: ierr
    ! Shut up the compiler 
    if (nsize==0 .or. data(1,1)==0 .or. bcaster==0) then 
    end if
    ierr=MPI_SUCCESS
    return
  end SUBroutine my_mpi_bcast2d_wp

!!!*******************************************************************
!!!
!!! Dummy routine for mpi scatter (working precision)
!!!
!!!*******************************************************************

  SUBroutine my_mpi_scatter_wp(sendbuf,sendsize,recvbuf,recvsize,sender,ierr)
    !use params, only : params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:) :: recvbuf
    integer, intent(in) :: sendsize, recvsize, sender
    integer, intent(out) :: ierr


    ! Shut up the compiler 
    if (sendsize==0 .or. sender==0 .or. sendbuf(1)==0 .or. recvsize==0 .or. recvbuf(1)==0) then 
    end if

    ierr=MPI_SUCCESS
    return
  end SUBroutine my_mpi_scatter_wp

!!!*******************************************************************
!!!
!!! Dummy routine for mpi gather (working precision)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_gather_wp(sendbuf,sendsize,recvbuf,recvsize,gatherer,ierr)
    !use params, only : params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:) :: recvbuf
    integer, intent(in) :: sendsize, recvsize, gatherer
    integer, intent(out) :: ierr
    ! Shut up the compiler 
    if (sendsize==0 .or. gatherer==0 .or. sendbuf(1)==0 .or. recvsize==0 .or. recvbuf(1)==0) then 
    end if

    ierr=MPI_SUCCESS
    return
  end SUBroutine my_mpi_gather_wp

!!!*******************************************************************
!!!
!!! Routine for 1D mpi reduce (SUM + integer)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_reduce1d_sum_int(sendbuf,recvbuf,sendsize,root,ierr)
    implicit none
    integer, intent(in), dimension(:) :: sendbuf
    integer, intent(inout), dimension(:) :: recvbuf
    integer, intent(in) :: sendsize, root
    integer, intent(out) :: ierr

    !shut up the compiler 
    if (root==0) then
    end if
    ierr=MPI_SUCCESS
    recvbuf(1:sendsize)=sendbuf

    return
  end SUBroutine my_mpi_reduce1d_sum_int

!!!*******************************************************************
!!!
!!! Routine for 2D mpi reduce (SUM + integer)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_reduce2d_sum_int(sendbuf,recvbuf,sendsize,root,ierr)
    implicit none
    integer, intent(in), dimension(:,:) :: sendbuf
    integer, intent(inout), dimension(:,:) :: recvbuf
    integer, intent(in) :: sendsize, root
    integer, intent(out) :: ierr

    !shut up the compiler 
    if (root==0) then
    end if
    ierr=MPI_SUCCESS
    recvbuf=sendbuf

    return
  end SUBroutine my_mpi_reduce2d_sum_int

!!!*******************************************************************
!!!
!!! Dummy routine for 1D mpi reduce (SUM + working precision)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_reduce1d_sum_wp(sendbuf,recvbuf,sendsize,root,ierr)
    !use params, only : params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:) :: recvbuf
    integer, intent(in) :: sendsize, root
    integer, intent(out) :: ierr
    
    !shut up the compiler 
    if (sendsize==0) then
    end if
    if (root==0) then
    end if
    ierr=MPI_SUCCESS
    recvbuf(1:sendsize)=sendbuf
    return
  end SUBroutine my_mpi_reduce1d_sum_wp


!!!*******************************************************************
!!!
!!! Dummy routine for 2D mpi reduce (SUM + working precision)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_reduce2d_sum_wp(sendbuf,recvbuf,sendsize,root,ierr)
    !use params, only : params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:,:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:,:) :: recvbuf
    integer, intent(in) :: sendsize, root
    integer, intent(out) :: ierr

    !shut up the compiler 
    if (sendsize==0) then
    end if
    if (root==0) then
    end if
    ierr=MPI_SUCCESS
    recvbuf=sendbuf
    return
  end SUBroutine my_mpi_reduce2d_sum_wp
  

!!!*******************************************************************
!!!
!!! Dummy routine for 3D mpi reduce (SUM + working precision)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_reduce3d_sum_wp(sendbuf,recvbuf,sendsize,root,ierr)
    !use params, only : params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:,:,:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:,:,:) :: recvbuf
    integer, intent(in) :: sendsize, root
    integer, intent(out) :: ierr
    
    !shut up the compiler 
    if (sendsize==0) then
    end if
    if (root==0) then
    end if
    ierr=MPI_SUCCESS
    recvbuf=sendbuf
    return
  end SUBroutine my_mpi_reduce3d_sum_wp


!!!*******************************************************************
!!!
!!! Dummy routine for 4D mpi reduce (SUM + working precision)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_reduce4d_sum_wp(sendbuf,recvbuf,sendsize,root,ierr)
    !use params, only : params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:,:,:,:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:,:,:,:) :: recvbuf
    integer, intent(in) :: sendsize, root
    integer, intent(out) :: ierr
    
    !shut up the compiler 
    if (sendsize==0) then
    end if
    if (root==0) then
    end if
    ierr=MPI_SUCCESS
    recvbuf=sendbuf
    return
  end SUBroutine my_mpi_reduce4d_sum_wp

!!!*******************************************************************
!!!
!!! Dummy routine for 5D mpi reduce (SUM + working precision)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_reduce5d_sum_wp(sendbuf,recvbuf,sendsize,root,ierr)
    !use params, only : params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:,:,:,:,:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:,:,:,:,:) :: recvbuf
    integer, intent(in) :: sendsize, root
    integer, intent(out) :: ierr
   
    !shut up the compiler 
    if (sendsize==0) then
    end if
    if (root==0) then
    end if
    ierr=MPI_SUCCESS
    recvbuf=sendbuf
    return
  end SUBroutine my_mpi_reduce5d_sum_wp


END module RFOF_mpi_module
