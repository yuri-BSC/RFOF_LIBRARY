!!! A module that allows to use the same source for both serial and parallel 
!!! versions. Serial code can be compiled without linking to mpi or having 
!!! a mpi compiler. Note in Makefile use mpi_serial.f90 for serial and this 
!!! for parallel

#include "config.h"

module RFOF_mpi_module

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  IMPLICIT NONE
#include "mpif.h"  

!   INCLUDE 'mpif.h'

!!! RFOF variables:
  INTEGER, PARAMETER :: R4 = SELECTED_REAL_KIND(6,37) ! Real*4
  INTEGER, PARAMETER :: R8 = SELECTED_REAL_KIND(15,300) ! Real*8
  INTEGER, PARAMETER :: params_wp = R8
  INTEGER, PARAMETER :: params_dp = R8
  INTEGER, save :: mpivar_id
  INTEGER, save :: mpivar_numproc
  LOGICAL, save :: mpivar_independentParallelJob
!!! RFOF variables.


contains

!!!*******************************************************************
!!!
!!! Routine for mpi init
!!!
!!!*******************************************************************

  subroutine my_mpi_init(ierr)
    !use mpivar
    implicit none
    integer, intent(out) :: ierr
    
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,mpivar_id,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,mpivar_numproc,ierr)
    IF (mpivar_numproc.EQ.0) mpivar_numproc=1

    mpivar_independentParallelJob = .false.
    
    return
  end subroutine my_mpi_init

!!!*******************************************************************
!!!
!!! Routine for mpi finalize
!!!
!!!*******************************************************************

  subroutine my_mpi_finalize(ierr)
    implicit none
    integer, intent(out) :: ierr

    CALL MPI_FINALIZE(ierr)

    return
  end subroutine my_mpi_finalize

!!!*******************************************************************
!!!
!!! Routine for mpi barrier
!!!
!!!*******************************************************************

  SUBroutine my_mpi_barrier(ierr)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: ierr

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    return
  end SUBroutine my_mpi_barrier

!!!*******************************************************************
!!!
!!! Routine for 0D mpi bcast (INTEGER)
!!!
!!!*******************************************************************

  SUBroutine my_mpi_bcast0d_int(data,nsize,bcaster,ierr)
    implicit none
    integer, intent(inout) :: data
    integer, intent(in) :: nsize, bcaster
    integer, intent(out) :: ierr

    CALL MPI_BCAST(data,nsize,MPI_INTEGER,bcaster,MPI_COMM_WORLD,ierr)

    return
  end SUBroutine my_mpi_bcast0d_int

!!!*******************************************************************
!!!
!!! Routine for 1D mpi bcast (INTEGER)
!!!
!!!*******************************************************************

  SUBroutine my_mpi_bcast1d_int(data,nsize,bcaster,ierr)
    implicit none
    integer, intent(inout), dimension(:) :: data
    integer, intent(in) :: nsize, bcaster
    integer, intent(out) :: ierr

    CALL MPI_BCAST(data,nsize,MPI_INTEGER,bcaster,MPI_COMM_WORLD,ierr)

    return
  end SUBroutine my_mpi_bcast1d_int

!!!*******************************************************************
!!!
!!! Routine for 2D mpi bcast (INTEGER)
!!!
!!!*******************************************************************

  SUBroutine my_mpi_bcast2d_int(data,nsize,bcaster,ierr)
    implicit none
    integer, intent(inout), dimension(:,:) :: data
    integer, intent(in) :: nsize, bcaster
    integer, intent(out) :: ierr

    CALL MPI_BCAST(data,nsize,MPI_INTEGER,bcaster,MPI_COMM_WORLD,ierr)

    return
  end SUBroutine my_mpi_bcast2d_int


!!!*******************************************************************
!!!
!!! Routine for mpi bcast (working precision)
!!!
!!!*******************************************************************

  SUBroutine my_mpi_bcast_wp(data,nsize,bcaster,ierr)
    !use params, only : params_dp, params_sp, params_wp
    implicit none
    real(kind=params_wp), intent(inout), dimension(:) :: data
    integer, intent(in) :: nsize, bcaster
    integer, intent(out) :: ierr
    IF (params_wp==params_dp) THEN
       CALL MPI_BCAST(data,nsize,MPI_DOUBLE_PRECISION, &
                      bcaster,MPI_COMM_WORLD,ierr)
    ELSE
       CALL MPI_BCAST(data,nsize,MPI_REAL, &
                      bcaster,MPI_COMM_WORLD,ierr)
    END IF
    return
  end SUBroutine my_mpi_bcast_wp


!!!*******************************************************************
!!!
!!! Routine for mpi bcast (working precision)
!!!
!!!*******************************************************************

  SUBroutine my_mpi_bcast2d_wp(data,nsize,bcaster,ierr)
    !use params, only : params_dp, params_sp, params_wp
    implicit none
    real(kind=params_wp), intent(inout), dimension(:,:) :: data
    integer, intent(in) :: nsize, bcaster
    integer, intent(out) :: ierr
    IF (params_wp==params_dp) THEN
       CALL MPI_BCAST(data,nsize,MPI_DOUBLE_PRECISION, &
                      bcaster,MPI_COMM_WORLD,ierr)
    ELSE
       CALL MPI_BCAST(data,nsize,MPI_REAL, &
                      bcaster,MPI_COMM_WORLD,ierr)
    END IF
    return
  end SUBroutine my_mpi_bcast2d_wp

!!!*******************************************************************
!!!
!!! Routine for mpi scatter (working precision)
!!!
!!!*******************************************************************

  SUBroutine my_mpi_scatter_wp(sendbuf,sendsize,recvbuf,recvsize,sender,ierr)
    !use params, only : params_dp, params_sp, params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:) :: recvbuf
    integer, intent(in) :: sendsize, recvsize, sender
    integer, intent(out) :: ierr

    IF (params_wp==params_dp) THEN
       CALL MPI_SCATTER(sendbuf,sendsize,MPI_DOUBLE_PRECISION,&
            recvbuf,recvsize, MPI_DOUBLE_PRECISION,&
            sender,MPI_COMM_WORLD,ierr)
    ELSE
       CALL MPI_SCATTER(sendbuf,sendsize,MPI_REAL,&
            recvbuf,recvsize, MPI_DOUBLE_PRECISION,&
            sender,MPI_COMM_WORLD,ierr)
    END IF
    
    return
  end SUBroutine my_mpi_scatter_wp

!!!*******************************************************************
!!!
!!! Routine for mpi gather (working precision)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_gather_wp(sendbuf,sendsize,recvbuf,recvsize,gatherer,ierr)
    !use params, only : params_dp, params_sp, params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:) :: recvbuf
    integer, intent(in) :: sendsize, recvsize, gatherer
    integer, intent(out) :: ierr
    
    IF (params_wp==params_dp) THEN
       CALL MPI_GATHER(sendbuf,sendsize,MPI_DOUBLE_PRECISION,&
            recvbuf,recvsize,MPI_DOUBLE_PRECISION,&
            gatherer,MPI_COMM_WORLD,ierr)
    ELSE
       CALL MPI_GATHER(sendbuf,sendsize,MPI_REAL,&
            recvbuf,recvsize,MPI_DOUBLE_PRECISION,&
            gatherer,MPI_COMM_WORLD,ierr)
    END IF

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
    
    ! shut up the compiler
    if (root==0) then
    end if
    recvbuf=0
    CALL MPI_Reduce(sendbuf, recvbuf, sendsize, MPI_INTEGER, &
         MPI_SUM, 0, MPI_COMM_WORLD, ierr)

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
    
    ! shut up the compiler
    if (root==0) then
    end if
    recvbuf=0
    CALL MPI_Reduce(sendbuf, recvbuf, sendsize, MPI_INTEGER, &
         MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    return
  end SUBroutine my_mpi_reduce2d_sum_int

!!!*******************************************************************
!!!
!!! Routine for 1D mpi reduce (SUM + working precision)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_reduce1d_sum_wp(sendbuf,recvbuf,sendsize,root,ierr)
    !use params, only : params_dp, params_sp, params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:) :: recvbuf
    integer, intent(in) :: sendsize, root
    integer, intent(out) :: ierr
    
    ! shut up the compiler
    if (root==0) then
    end if
    recvbuf=0.0_params_wp
    IF (params_wp==params_dp) THEN
       CALL MPI_Reduce(sendbuf, recvbuf, sendsize, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ELSE
       CALL MPI_Reduce(sendbuf, recvbuf, sendsize, MPI_REAL, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    END IF

    return
  end SUBroutine my_mpi_reduce1d_sum_wp

!!$!!! Routine for 2D mpi reduce (SUM + working precision)
!!$!!!
!!$!!!*******************************************************************
!!$  
!!$  SUBroutine my_mpi_reduce2d_sum_wp(sendbuf,recvbuf,sendsize1,sendsize2,root,ierr)
!!$    !use params, only : params_dp, params_sp, params_wp
!!$    implicit none
!!$    real(kind=params_wp), intent(in), dimension(:) :: sendbuf
!!$    real(kind=params_wp), intent(inout), dimension(:) :: recvbuf
!!$    integer, intent(in) :: sendsize1, sendsize2, root
!!$    integer, intent(out) :: ierr
!!$    
!!$    ! Local
!!$    integer :: j1, j2
!!$
!!$    ! Local - for MPI call
!!$    real(8), allocatable :: svec(:)
!!$    real(8), allocatable :: rvec(:)
!!$    integer :: slen
!!$
!!$    n1 = size(sendbuf,1)
!!$    n2 = size(sendbuf,2)
!!$    slen = sendsize1 * sendsize2
!!$    allocate( svec(sendlen) , rvec(sendlen) )
!!$
!!$    do j1=1,sendsize1
!!$       do j2=1,sendsize2
!!$          svec(j2+(sendsize2-1)*j1)=sendbuf(j1,j2)
!!$       enddo
!!$    enddo
!!$    rvec(:)=0
!!$
!!$    call my_mpi_reduce1d_sum_wp(svec,rvec,slen,root,ierr)
!!$
!!$    do jf=1,sendsize1
!!$       do jn=1,sendsize2
!!$          recvbuf(j1,j2) = sendbuf(j2+(sendsize2-1)*j1)
!!$       enddo
!!$    enddo
!!$
!!$    deallocate( svec , rvec )
!!$
!!$    return
!!$  end SUBroutine my_mpi_reduce2d_sum_wp


!!!*******************************************************************
!!!
!!! Routine for 2D mpi reduce (SUM + working precision)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_reduce2d_sum_wp(sendbuf,recvbuf,sendsize,root,ierr)
    !use params, only : params_dp, params_sp, params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:,:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:,:) :: recvbuf
    integer, intent(in) :: sendsize, root
    integer, intent(out) :: ierr
    
    ! shut up the compiler
    if (root==0) then
    end if
    recvbuf=0.0_params_wp
    IF (params_wp==params_dp) THEN
       CALL MPI_Reduce(sendbuf, recvbuf, sendsize, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ELSE
       CALL MPI_Reduce(sendbuf, recvbuf, sendsize, MPI_REAL, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    END IF
    return
  end SUBroutine my_mpi_reduce2d_sum_wp


!!!*******************************************************************
!!!
!!! Routine for 3D mpi reduce (SUM + working precision)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_reduce3d_sum_wp(sendbuf,recvbuf,sendsize,root,ierr)
    !use params, only : params_dp, params_sp, params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:,:,:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:,:,:) :: recvbuf
    integer, intent(in) :: sendsize, root
    integer, intent(out) :: ierr
    
    ! shut up the compiler
    if (root==0) then
    end if
    recvbuf=0.0_params_wp
    IF (params_wp==params_dp) THEN
       CALL MPI_Reduce(sendbuf, recvbuf, sendsize, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ELSE
       CALL MPI_Reduce(sendbuf, recvbuf, sendsize, MPI_REAL, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    END IF

    return
  end SUBroutine my_mpi_reduce3d_sum_wp


!!!*******************************************************************
!!!
!!! Routine for 4D mpi reduce (SUM + working precision)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_reduce4d_sum_wp(sendbuf,recvbuf,sendsize,root,ierr)
    !use params, only : params_dp, params_sp, params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:,:,:,:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:,:,:,:) :: recvbuf
    integer, intent(in) :: sendsize, root
    integer, intent(out) :: ierr
    
    ! shut up the compiler
    if (root==0) then
    end if
    recvbuf=0.0_params_wp
    IF (params_wp==params_dp) THEN
       CALL MPI_Reduce(sendbuf, recvbuf, sendsize, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ELSE
       CALL MPI_Reduce(sendbuf, recvbuf, sendsize, MPI_REAL, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    END IF
    
    return
  end SUBroutine my_mpi_reduce4d_sum_wp

!!!*******************************************************************
!!!
!!! Routine for 5D mpi reduce (SUM + working precision)
!!!
!!!*******************************************************************
  
  SUBroutine my_mpi_reduce5d_sum_wp(sendbuf,recvbuf,sendsize,root,ierr)
    !use params, only : params_dp, params_sp, params_wp
    implicit none
    real(kind=params_wp), intent(in), dimension(:,:,:,:,:) :: sendbuf
    real(kind=params_wp), intent(inout), dimension(:,:,:,:,:) :: recvbuf
    integer, intent(in) :: sendsize, root
    integer, intent(out) :: ierr
    
    ! shut up the compiler
    if (root==0) then
    end if
    recvbuf=0.0_params_wp
    IF (params_wp==params_dp) THEN
       CALL MPI_Reduce(sendbuf, recvbuf, sendsize, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ELSE
       CALL MPI_Reduce(sendbuf, recvbuf, sendsize, MPI_REAL, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    END IF

    return
  end SUBroutine my_mpi_reduce5d_sum_wp


END module RFOF_mpi_module
