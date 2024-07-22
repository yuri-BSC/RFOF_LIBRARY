#include "config.h"

module RFOF_Efield_update

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  implicit none

  real(8), allocatable :: recvbuf(:,:)
  integer, allocatable :: recvbuf_int2d(:,:)

contains

  !--------------------------------------------------------------------------------
  subroutine destructor_RFOF_Efield_update

    if (allocated(recvbuf)) then
       deallocate(recvbuf)
    endif
    if (allocated(recvbuf_int2d)) then
       deallocate(recvbuf_int2d)
    endif

  end subroutine destructor_RFOF_Efield_update

  !--------------------------------------------------------------------------------
  ! subroutine update_efield_normalisation(calcRFpower,wave)
  subroutine update_efield_normalisation(dt, wave,MPI_node_Id,ierr,RFpower,diagno)

    use RFOF_waves
    use RFOF_diagnostics

    ! Input and Output
    real(8), intent(in) :: dt
    type(rf_wave_global), intent(inout) :: wave
    integer, intent(in) :: MPI_node_Id
    integer, intent(out) :: ierr

    ! OPTIONAL Input and Output
    real(8), intent(in), optional :: RFpower(:,:)
    type(RFOF_cumlative_diagnostics), intent(inout), optional :: diagno

    ! Local
    real(8) :: P, dP, dPmin, sigmaP, c, Pcmp
    integer :: jf, jn
    logical :: do_update

    real(8), allocatable :: pabs(:,:), pabs2(:,:)
    integer, allocatable :: nkicks(:,:)

    do_update = .true.

    allocate( pabs(   size(wave%RFpower,1) , size(wave%RFpower,2) ) )
    allocate( pabs2(  size(wave%RFpower,1) , size(wave%RFpower,2) ) )
    allocate( nkicks( size(wave%RFpower,1) , size(wave%RFpower,2) ) )

    !write(0,*)'Calling update_efield_normalisation',present(RFpower),present(diagno),diagno%time_end_sum_diagnostics,&
    !     diagno%time_start_sum_diagnostics, allocated(diagno%energy_to_markers_from_mode), &
    !     size(diagno%energy_to_markers_from_mode)

    ierr=0
    if ( present(RFpower) ) then

       do jf = 1, wave%nfreq
          do jn = 1, wave%nnphi(jf)
             if (RFpower(jf,jn) .gt. 0.) THEN
                ! Update E-field normalisation
                wave%EnormalisationFactor(jf,jn) = &
                     wave%EnormalisationFactor(jf,jn) * &
                     sqrt( wave%RFpower(jf,jn) / RFpower(jf,jn) )
             endif
          enddo
       enddo

    elseif ( present(diagno) ) then

       ! Check that the time step since previous update of the E-field absorption is finite
       if ( dt < 1d-40 ) then
          write(0,*)'WARNING in update_efield_normalisation: Diagnostic integration time = ',dt
          write(0,*)'        diagno%time_end_sum_diagnostics   =',diagno%time_end_sum_diagnostics
          write(0,*)'        diagno%time_start_sum_diagnostics =',diagno%time_start_sum_diagnostics
          pabs(:,:) = 0d0          
          pabs2(:,:) = 0d0          
       else
          pabs(:,:)  = diagno%energy_to_markers_from_mode(:,:) / dt
          pabs2(:,:) = diagno%energy_square_to_markers_from_mode(:,:) / dt
       endif
       nkicks = diagno%kick_counter

       call RFOF_merge_power_absorption_diagno_from_mpi_nodes(pabs,pabs2,nkicks,ierr)

       if (output__efield_normalization) then
          write(0,*)
          write(0,*)'Calling save2file_RFOF_efield_normalization'
          write(0,*)'  time:',diagno%time_end_sum_diagnostics
          write(0,*)'  nkickc:', nkicks
          write(0,*)'  Pabs:',Pabs, ', Prequested:',wave%RFpower
          write(0,*)'  E-normalisation:',wave%EnormalisationFactor
          write(0,*)
          call save2file_RFOF_efield_normalization(diagno%time_end_sum_diagnostics, Pabs, wave, 0)
       endif

       do jf = 1, wave%nfreq
          do jn = 1, wave%nnphi(jf)
             if ( (wave%RFpower(jf,jn) .gt. 0.) .and. &
                  (pabs(jf,jn) .le. 0.) ) then
                write(0,*)'Warning1: skipping update_Efield_normalisation',jf,jn,wave%RFpower(jf,jn),pabs(jf,jn)
                do_update = .false.
             endif
          enddo
       enddo

       if ( do_update ) then
          do jf = 1, wave%nfreq
             do jn = 1, wave%nnphi(jf)
                if (wave%RFpower(jf,jn) .gt. 0.) then

                   ! The power absorption in RFOF / the difference cmp to the prescibed RF power
                   P = pabs(jf,jn)
                   dP = P - wave%RFpower(jf,jn)

                   ! The absolute error-bar should be compared to either P or dP ; use the smallest value
                   dPmin = min( abs(P) , abs(dP) )

                   ! Statistically evaluate standard deviation (uncertainty) of the power absorption
                   sigmaP = sqrt(abs( pabs2(jf,jn) - pabs(jf,jn)**2 ) / nkicks(jf,jn) )

                   !write(0,*)"time since last normalisation:",dt, ' , sigma(P)=',sigmaP, ' , dP=', dPmin

                   ! Since there is a statistical error, how trust-worthy is this change in the absorption?
                   ! A function "c" is here generated to estimate the fraction of the power change that can be concidered trustworthy:
                   if ( dPmin > sigmaP) then
                      c=1d0-0.5d0*sigmaP/dPmin
                   else
                      if ( sigmaP .eq. 0d0) then
                         c=0.5d0
                      else
                         c=0.5d0*dPmin/sigmaP
                      endif
                   endif

                   ! The approximate power absorption, using the absorption from previous time steps 
                   !to make up for the uncertainty of the present time steps
                   Pcmp = wave%RFpower(jf,jn) + c * dP

                   ! Update E-field normalisation
                   wave%EnormalisationFactor(jf,jn) = &
                        wave%EnormalisationFactor(jf,jn) * &
                        sqrt( wave%RFpower(jf,jn) / Pcmp )

                   write(0,*)
                   write(0,*)'Update E-field normalisation jfreq,jn:',jf,jn
                   write(0,*)'         Pabs[MW]=', pabs(jf,jn)*1d-6
                   write(0,*)'   Prequested[MW]=', wave%RFpower(jf,jn)*1d-6
                   write(0,*)'             dP/P=', dP     / wave%RFpower(jf,jn)
                   write(0,*)'          dPmin/P=', dPmin  / wave%RFpower(jf,jn)
                   write(0,*)'         sigmaP/P=', sigmaP / wave%RFpower(jf,jn)
                   write(0,*)'          dPeff/P=', Pcmp / wave%RFpower(jf,jn)
                   write(0,*)'          nkicks =', diagno%kick_counter(jf,jn)
                   write(0,*)' New: wave%EnormalisationFactor   =',wave%EnormalisationFactor(jf,jn)
                   write(0,*)


                else
                   write(0,*)'Warning2: skipping update_Efield_normalisation',jf,jn,wave%RFpower(jf,jn),pabs(jf,jn)

                endif
             enddo  !! do jn = 1, wave%nnphi(jf)
          enddo   !! do jf = 1, wave%nfreq

          call reset_to_zero_RFOF_cumlative_diagnostics(diagno)

       endif
    endif

    deallocate( pabs  )
    deallocate( pabs2 )

  end subroutine update_efield_normalisation


  !--------------------------------------------------------------------------------
  ! subroutine update_efield_normalisation(calcRFpower,wave)
  subroutine RFOF_merge_absorption_diagno_from_mpi_nodes(diagno,ierr)

    use RFOF_diagnostics

    !Input/Output
    type(RFOF_cumlative_diagnostics), INTENT(INOUT) :: diagno
    integer, intent(out) :: ierr

    ! Local
    integer, parameter :: rfof_mpi_root = 0
    integer :: ierr1=0, ierr2=0, ierr3=0, ierr4=0, ierr5=0, ierr6=0
    
    call RFOF_mpi_reduce_and_bcast_2d_wp(diagno%energy_to_markers_from_mode,        rfof_mpi_root, ierr1)

    call RFOF_mpi_reduce_and_bcast_2d_wp(diagno%energy_square_to_markers_from_mode, rfof_mpi_root, ierr2)
    call RFOF_mpi_reduce_and_bcast_2d_wp(diagno%sum_weight_at_resonance_with_mode,  rfof_mpi_root, ierr3)


    !call RFOF_mpi_reduce_and_bcast_0d_wp(diagno%time_start_sum_diagnostics, rfof_mpi_root, ierr4)
    !call RFOF_mpi_reduce_and_bcast_0d_wp(diagno%time_end_sum_diagnostics  , rfof_mpi_root, ierr5)

    call RFOF_mpi_reduce_and_bcast_2d_int(diagno%kick_counter, rfof_mpi_root, ierr6)

    ierr = min(ierr1 , min(ierr2 , min(ierr3 , min(ierr4 , min(ierr5 , ierr6)))))
    if ( ierr > -1) then
       ierr = max(ierr1 , max( ierr2 , max(ierr3 , max(ierr4 , max(ierr5 , ierr6)))))
    endif

  end subroutine RFOF_merge_absorption_diagno_from_mpi_nodes


  !--------------------------------------------------------------------------------
  ! subroutine RFOF_merge_power_absorption_diagno_from_mpi_nodes(calcRFpower,wave)
  subroutine RFOF_merge_power_absorption_diagno_from_mpi_nodes(pabs,pabs2,count,ierr)

    use RFOF_diagnostics

    !Input/Output
    real(8), intent(inout) :: pabs(:,:), pabs2(:,:)
    integer, intent(inout) :: count(:,:)
    integer, intent(out) :: ierr

    ! Local
    integer, parameter :: rfof_mpi_root = 0
    integer :: ierr1=0, ierr2=0, ierr3=0

    call RFOF_mpi_reduce_and_bcast_2d_wp( pabs , rfof_mpi_root, ierr1)
    call RFOF_mpi_reduce_and_bcast_2d_wp( pabs2, rfof_mpi_root, ierr2)
    call RFOF_mpi_reduce_and_bcast_2d_int( count, rfof_mpi_root, ierr3)

    ierr = min(ierr1 , min(ierr2 , ierr3))
    if ( ierr > -1) then
       ierr = max(ierr1 , max(ierr2 , ierr3))
    endif

  end subroutine RFOF_merge_power_absorption_diagno_from_mpi_nodes



  !--------------------------------------------------------------------------------
  subroutine RFOF_mpi_reduce_and_bcast_2d_wp(data,rfof_mpi_root,ier)

    USE RFOF_mpi_module ! wrapper module for 'include mpif.h'

    real(8), intent(inout) :: data(:,:)
    integer, intent(in) :: rfof_mpi_root
    integer, intent(out) :: ier

    integer :: sendsize

    if (.not. ALLOCATED(recvbuf)) then
       allocate(recvbuf( size(data, 1), size(data, 2) ))
    endif

    sendsize = size(data)

    call my_mpi_reduce2d_sum_wp(data, recvbuf, sendsize, rfof_mpi_root, ier)
    if ( ier < 0 ) then
       write(0,*)"ERROR in RFOF/RFOF_mpi_reduce_and_bcast_2d_wp: my_mpi_reduce2d_sum_wp failed; ier=",ier
       return
    endif

    call my_mpi_bcast2d_wp(recvbuf, sendsize, rfof_mpi_root, ier)
    if ( ier < 0 ) then
       write(0,*)"ERROR in RFOF/RFOF_mpi_reduce_and_bcast_2d_wp: my_mpi_bcast2d_wp failed; ier=",ier
       return
    endif

    data = recvbuf

  end subroutine RFOF_mpi_reduce_and_bcast_2d_wp


  !--------------------------------------------------------------------------------
  subroutine RFOF_mpi_reduce_and_bcast_0d_wp(data,rfof_mpi_root,ier)

    USE RFOF_mpi_module ! wrapper module for 'include mpif.h'

    real(8), intent(inout) :: data
    integer, intent(in) :: rfof_mpi_root
    integer, intent(out) :: ier

    real(8) :: send(1)
    real(8) :: recv(1)

    send(1)=data

    call my_mpi_reduce1d_sum_wp(send, recv, 1, rfof_mpi_root, ier)
    if ( ier < 0 ) then
       write(0,*)"ERROR in RFOF/RFOF_mpi_reduce_and_bcast_0d_wp: my_mpi_reduce1d_sum_wp failed; ier=",ier
       return
    endif

    call my_mpi_bcast_wp(recv, 1, rfof_mpi_root, ier)
    if ( ier < 0 ) then
       write(0,*)"ERROR in RFOF/RFOF_mpi_reduce_and_bcast_0d_wp: my_mpi_bcast_wp failed; ier=",ier
       return
    endif

    data = recv(1)

  end subroutine RFOF_mpi_reduce_and_bcast_0d_wp


  !--------------------------------------------------------------------------------
  subroutine RFOF_mpi_reduce_and_bcast_0d_int(data,rfof_mpi_root,ier)

    USE RFOF_mpi_module ! wrapper module for 'include mpif.h'

    integer, intent(inout) :: data
    integer, intent(in) :: rfof_mpi_root
    integer, intent(out) :: ier

    integer :: send(1)
    integer :: recv(1)

    send(1)=data

    call my_mpi_reduce1d_sum_int(send, recv, 1, rfof_mpi_root, ier)
    if ( ier < 0 ) then
       write(0,*)"ERROR in RFOF/RFOF_mpi_reduce_and_bcast_0d_int: my_mpi_reduce1d_sum_int failed; ier=",ier
       return
    endif

    call my_mpi_bcast1d_int(recv, 1, rfof_mpi_root, ier)
    if ( ier < 0 ) then
       write(0,*)"ERROR in RFOF/RFOF_mpi_reduce_and_bcast_0d_int: my_mpi_bcast1d_int failed; ier=",ier
       return
    endif

    data = recv(1)

  end subroutine RFOF_mpi_reduce_and_bcast_0d_int

  !--------------------------------------------------------------------------------
  subroutine RFOF_mpi_reduce_and_bcast_2d_int(data,rfof_mpi_root,ier)

    USE RFOF_mpi_module ! wrapper module for 'include mpif.h'

    integer, intent(inout) :: data(:,:)
    integer, intent(in) :: rfof_mpi_root
    integer, intent(out) :: ier

    integer :: sendsize

    if (.not. ALLOCATED(recvbuf_int2d)) then
       allocate(recvbuf_int2d( size(data, 1), size(data, 2) ))
    endif

    sendsize = size(data)

    call my_mpi_reduce2d_sum_int(data, recvbuf_int2d, sendsize, rfof_mpi_root, ier)
    if ( ier < 0 ) then
       write(0,*)"ERROR in RFOF/RFOF_mpi_reduce_and_bcast_0d_int: my_mpi_reduce1d_sum_int failed; ier=",ier
       return
    endif

    call my_mpi_bcast2d_int(recvbuf_int2d, sendsize, rfof_mpi_root, ier)
    if ( ier < 0 ) then
       write(0,*)"ERROR in RFOF/RFOF_mpi_reduce_and_bcast_0d_int: my_mpi_bcast1d_int failed; ier=",ier
       return
    endif

    data = recvbuf_int2d

  end subroutine RFOF_mpi_reduce_and_bcast_2d_int


  !--------------------------------------------------------------------------------
  ! enorm_stats_test - test routine
  !--------------------------------------------------------------------------------
  subroutine enorm_stats_test

    use RFOF_waves
    use RFOF_diagnostics

    ! Input and Output
    type(rf_wave_global) :: wave
    type(RFOF_cumlative_diagnostics) :: diagno
    integer :: ierr

    integer :: nn = 1
    integer :: nf = 1
    integer :: jn, jf, j

    integer :: nsamples = 4
    real(8) :: x, w

    jn = nn
    jf = nf

    ! SET WAVE:
    allocate( wave%nnphi( nf ) , &
         wave%RFpower(nf,nn), &
         wave%EnormalisationFactor(nf,nn) )

    wave%nfreq = nf
    wave%nnphi(jf) = 1

    wave%RFpower(jf,jn) = 1.0
    wave%EnormalisationFactor(jf,jn) = 1.0


    ! SET DIAGNO
    diagno%time_start_sum_diagnostics = 0.0
    diagno%time_end_sum_diagnostics   = 1.0

    allocate( &
         diagno%energy_to_markers_from_mode(nn,nf), &
         diagno%energy_square_to_markers_from_mode(nn,nf), &
         diagno%sum_weight_at_resonance_with_mode(nn,nf), &
         diagno%kick_counter(nn,nf), &
         diagno%toroidal_momentum_to_markers_from_mode(nf,nn) )

    diagno%energy_to_markers_from_mode(jf,jn) = 0d0
    diagno%sum_weight_at_resonance_with_mode(jf,jn) = 0d0
    diagno%energy_square_to_markers_from_mode (jf,jn) = 0d0
    diagno%kick_counter(jf,jn) = 0

    

    do j=1,nsamples
       x = 4.0d0 + 4.*(-0.5 + 1.0*real(j-1)/real(nsamples-1))
       w = 1d0/real(nsamples)

       diagno%sum_weight_at_resonance_with_mode(jf,jn) = diagno%sum_weight_at_resonance_with_mode(jf,jn) + w
       diagno%energy_to_markers_from_mode(jf,jn) = diagno%energy_to_markers_from_mode(jf,jn) + w*x
       diagno%energy_square_to_markers_from_mode (jf,jn) = diagno%energy_square_to_markers_from_mode (jf,jn) + w*x**2
       diagno%kick_counter(jf,jn) = diagno%kick_counter(jf,jn) + 1
       
    enddo

    call update_efield_normalisation(1d0,wave,0,ierr,diagno=diagno)

    deallocate( wave%nnphi, &
         wave%RFpower, &
         wave%EnormalisationFactor)

    deallocate( diagno%energy_to_markers_from_mode, &
         diagno%energy_square_to_markers_from_mode, &
         diagno%toroidal_momentum_to_markers_from_mode, &
         diagno%sum_weight_at_resonance_with_mode, &
         diagno%kick_counter)

  end subroutine enorm_stats_test




end module RFOF_Efield_update
