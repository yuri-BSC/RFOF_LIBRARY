#include "config.h"

module RFOF_resonance_condition

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  use RFOF_parameters
  use RFOF_waves
  use RFOF_markers
  use RFOF_resonance_memory
  use RFOF_numerics

  implicit none

contains


  !--------------------------------------------------------------------------------
  ! subroutine rf_wave_resonance
  !
  subroutine rf_wave_particles_resonance(time,time_previous_step, &
       marker,resonanceMemory,RFglobal, &
       nr_resonant_waves, list_resonant_waves, &
       time_at_resonance, &
       particle_overshot_resonance,MPI_nod_Id)

    ! Input
    real(8), intent(in) :: time
    real(8), intent(in) :: time_previous_step
    type(particle), intent(in) :: marker
    type(rf_wave_global), intent(inout) :: RFglobal
    integer, intent(in) :: MPI_nod_Id

    ! Input/Output
    type(resonance_memory), intent(inout) :: resonanceMemory(:,:)

    ! Output
    logical, intent(out) :: particle_overshot_resonance
    integer, intent(out) :: nr_resonant_waves
    real(8), pointer, intent(out) :: list_resonant_waves(:,:)
    real(8),    intent(out) :: time_at_resonance

    ! Local
    type(rf_wave_local) :: RFlocal
    real(8) :: omega_res
    integer :: nharm, jfreq, jnphi
    logical :: closeEnoughToResonance
    logical :: hasCrossedResonance
    logical :: resonanceAlreadyTreatedInPreviousTimeStep
    real(8) :: t_at_resonance
    real(8) :: R_at_resonance
    real(8) :: phi_at_resonance
    real(8) :: z_at_resonance
    real(8) :: Dot_omega_res, DotDot_omega_res
    real(8) :: timeAtFirstResonance
    logical :: was_stored
    integer :: errorFlag
    !real(8) :: delta_test

    ! Preliminary value
    particle_overshot_resonance = .false.
    time_at_resonance = time + very_large_time_in_seconds    ! Dummy time - must be larger than any possible resonance time
    timeAtFirstResonance = time + very_large_time_in_seconds ! Dummy time - must be larger than any possible resonance time

    if ( abs(marker%omega_gyro) .lt. 1e-10) then
       write(0,*)
       write(0,*)"ERROR in rf_wave_particles_resonance: marker%omega_gyro",marker%omega_gyro, " is too small"
       write(0,*)
       nr_resonant_waves = 0
       return
    endif

    do jfreq=1,RFglobal%nfreq
       do jnphi=1,RFglobal%nnphi(jfreq)

          ! Local wave field needed for local wave vector
          RFglobal%j_freq_present_wave = jfreq
          RFglobal%j_nphi_present_wave = jnphi
          RFlocal = get_rf_wave_local(marker%R,marker%phi,marker%z,RFglobal)

          call wave_particles_resonance_function(marker, RFlocal, omega_res, nharm)

          call estimate_resonance_location(resonanceMemory(jfreq,jnphi), &
               time,marker%R,marker%phi,marker%z,omega_res, nharm, &
               marker%time_acceleration, time_previous_step, &
               t_at_resonance,R_at_resonance,phi_at_resonance,z_at_resonance, &
               Dot_omega_res,DotDot_omega_res, MPI_nod_Id,errorFlag)

          if ( (errorFlag .gt. 0) .or. (t_at_resonance.lt.time_previous_step) ) then
             t_at_resonance = time + very_large_time_in_seconds
          endif

          ! Store in memory
          call storeNewPointResonanceMemory(resonanceMemory(jfreq,jnphi), &
               t_at_resonance,time,nharm,marker%time_acceleration,marker%R,marker%phi,marker%z, &
               omega_res,marker%omega_gyro,was_stored)

          ! Add to list of resonances that the particle is at, but only if:
          !  1. Resonance was correctly calculate; no erroFlag
          !  2. Are we close enough to the resonance? If not we have to see if particle has overshot the resonance.
          !  3. Has resonance been treated in previous timesteps?
          if ( (errorFlag .eq. 0) .and. (was_stored .eqv. .TRUE.) ) then
             if (close_enough_to_resonance(t_at_resonance,R_at_resonance,z_at_resonance,time,time_previous_step,marker)) then
                if (resonance_already_treated_in_previous_timestep(resonanceMemory(jfreq,jnphi), &
                     time,marker%tauBounce,Dot_omega_res) .eqv. .false.) then

                   call add_to_sorted_list_3columns(list_resonant_waves, nr_resonant_waves, &
                        t_at_resonance, dble(jfreq), dble(jnphi))

                   resonanceMemory(jfreq,jnphi)%previous_resonance_exists = .TRUE.
                   resonanceMemory(jfreq,jnphi)%time_of_last_resonance_crossing = t_at_resonance
                   resonanceMemory(jfreq,jnphi)%sign_omega_dot_at_last_resonance_crossing = sign(1d0, Dot_omega_res)

                endif

             else

                !  Has particle overshot resonance?
                if ( t_at_resonance .lt. time ) then
                   ! Either the marker has crossed the resonance (refine time step!), or we're looking at a resonance that already's been accounted for (do nothing)
                   if ( resonance_already_treated_in_previous_timestep(resonanceMemory(jfreq,jnphi), &
                        time,marker%tauBounce,Dot_omega_res) .eqv. .false. ) then
                      nr_resonant_waves = 0
                      time_at_resonance = t_at_resonance
                      particle_overshot_resonance = .TRUE.
                      !write(0,'(A,20E14.5)')'OVERSHOT!',time,t_at_resonance,time_previous_step
                      return
                   endif
                else
                   ! Update the prediction for the upcoming resonance position: timeAtFirstResonance
                   timeAtFirstResonance = min( timeAtFirstResonance , t_at_resonance )
                   !!if ( abs(marker%R-3.0) .lt.0.07) &
                   !write(0,*)'not overshot - not res',time,t_at_resonance,timeAtFirstResonance, marker%R
                endif
             endif
          endif
       end do
    end do

    ! Store "time_at_resonance" as the time at the first resonance which we still haven't reached ()
    time_at_resonance = timeAtFirstResonance

    !write(0,'(A,2E13.5)') "rescon:estimate next resonance at time=", time_at_resonance, time_previous_step

    if (particle_overshot_resonance) then
       write(0,*)'Need to refine resonance',marker%R-R_at_resonance
    !elseif (nr_resonant_waves.gt.0) then
    !   write(0,*)'nr_resonant_waves',nr_resonant_waves
    endif

  end subroutine rf_wave_particles_resonance


  !--------------------------------------------------------------------------------
  ! subroutine wave_particles_resonance_function
  subroutine wave_particles_resonance_function(marker, RFlocal, omega_res, nharm)

    ! Input
    type(particle), intent(in) :: marker
    type(rf_wave_local), intent(in) :: RFlocal

    ! Output
    real(8), intent(out) :: omega_res
    integer, intent(out) :: nharm

    ! Local
    real(8) :: doppler

    ! Doppler shift, $k \cdot v$
    if (simplify__parallel_velocity_does_not_affect_resonance_condition) then
       doppler = 0d0
    else
       doppler = marker%vpar * RFlocal%kpar
    endif
    if ( .not. simplify__drift_velocity_does_not_affect_resonance_condition ) then
       doppler = doppler + &
            marker%vDriftRho * RFlocal%krho + &
            marker%vDriftDia * RFlocal%kdia
    endif

    ! Harmonic number closest to resonance
    nharm = int( 0.5 + abs( (RFlocal%omega - doppler) / marker%omega_gyro ) )

    ! How far from resonance?
    ! Difference between the Doppler shifted wave frequency and the harmonic cyclotron frequency
    omega_res = RFlocal%omega - doppler - real(nharm) * marker%omega_gyro

  end subroutine wave_particles_resonance_function


  !--------------------------------------------------------------------------------
  ! function close_enough_to_resonance
  !
  function close_enough_to_resonance(t_res,R_res,z_res,time,time_previous_step,marker) result(close_enough)

    ! Input
    real(8), intent(in) :: t_res
    real(8), intent(in) :: R_res
    real(8), intent(in) :: z_res
    real(8), intent(in) :: time
    real(8), intent(in) :: time_previous_step
    type(particle), intent(in) :: marker

    ! Output
    logical close_enough

    if ( (marker%R - R_res)**2 + (marker%z - z_res)**2 .lt. width_of_resonance_in_R(marker)**2 ) then
       close_enough=.true.
    else
       if ( abs(t_res-time_previous_step) .lt. 1e-3*(time-time_previous_step) ) then
          if (ignore_non_serious_error_messages .eqv. .FALSE.) then
             write(0,*)'----------------------------------------------'
             write(0,*)'WARNING in RFOF::close_enough_to_resonance' 
             write(0,*)'    Bad convergence; Particle far from resonance '
             write(0,*)'    in space, but close to in time'
             write(0,*)'      dt=',t_res-time_previous_step,' (compare to time step=',&
                  time-time_previous_step,')'
             write(0,*)'      dR=',R_res-marker%R,' (compare to the width of the resonance-layer=',&
                  width_of_resonance_in_R(marker),')'
             write(0,*)'----------------------------------------------'
          endif
          close_enough=.true.
       else
          close_enough=.false.
       end if
    end if

    !!if (t_res.gt.6.2d-7) &
    !!  print *, "close_enough_to_resonance", close_enough, marker%R, R_res, marker%R-R_res, width_of_resonance_in_R(marker)

  end function close_enough_to_resonance

  !--------------------------------------------------------------------------------
  ! function width_of_resonance_in_R
  !
  function width_of_resonance_in_R(marker) result(width)

    ! Input
    type(particle), intent(in) :: marker

    ! Output
    real(8) :: width

    width = width_of_rf_resonance_layer * marker%R

  end function width_of_resonance_in_R


  !--------------------------------------------------------------------------------
  ! function resonance_already_treated_in_previous_timestep
  !
  function resonance_already_treated_in_previous_timestep(mem,time,Torbit,Dot_omega_res) result(alreadyTreated)

    type(resonance_memory), intent(in) :: mem
    real(8), intent(in) :: time, Torbit, Dot_omega_res

    logical :: alreadyTreated

    alreadyTreated = .FALSE.

    if (mem%previous_resonance_exists) then

       !!if (time.gt.6.2d-7) &
       !!     print *, 'already_res... 1',time , mem%time_of_last_resonance_crossing, &
       !!     abs(time - mem%time_of_last_resonance_crossing) ,  0.3 * Torbit, &
       !!     (abs(time - mem%time_of_last_resonance_crossing) .lt. 0.3 * Torbit )

       if ( abs(time - mem%time_of_last_resonance_crossing) .lt. 0.5 * Torbit ) then

          !!if (time.gt.6.2d-7) &
          !!print *, 'already_res... 2',int(sign(1d0, Dot_omega_res)) , mem%sign_omega_dot_at_last_resonance_crossing

          if ( int(sign(1d0, Dot_omega_res)) .eq. mem%sign_omega_dot_at_last_resonance_crossing) then
             !!print *, 'true'
             alreadyTreated = .TRUE.
          endif
       endif
    endif

  end function resonance_already_treated_in_previous_timestep


!!$  !--------------------------------------------------------------------------------
!!$  !> function has_particle_overshot_resonance
!!$  !--------------------------------------------------------------------------------
!!$  function has_particle_overshot_resonance(time,time_at_resonance,mem) result(overshot)
!!$
!!$    ! Input
!!$    real(8), intent(in) :: time,time_at_resonance
!!$    type(resonance_memory), intent(in) :: mem
!!$
!!$    ! Output
!!$    logical :: overshot
!!$
!!$    if (time.gt.6.2d-7) &
!!$         print *, 'in overshoot1:',time_at_resonance,time_at_resonance, &
!!$         time,mem%previous_resonance_exists
!!$
!!$    !> Overshoot condition 1: time_at_resonance ash to be in the past
!!$    if (time .gt. time_at_resonance) then
!!$
!!$       !> Overshoot condition 2.1: there has to be a previous resoance or...
!!$       if (mem%previous_resonance_exists .eqv. .FALSE.) then
!!$          if (time.gt.6.2d-7) &
!!$               print * , 'OVERSHOT!!! (no res)'
!!$          overshot = .TRUE.
!!$          return
!!$       else
!!$          if (time.gt.6.2d-7) &
!!$               print * ,  'overshoot?', &
!!$               abs(time_at_resonance - mem%time_of_last_resonance_crossing), &
!!$               ' > ',&
!!$               abs(time_at_resonance - time)
!!$
!!$          !> Overshoot condition 2.2: ...or ???
!!$          if ( abs(time_at_resonance - mem%time_of_last_resonance_crossing) .gt. &
!!$               abs(time_at_resonance - time) ) then
!!$             if (time.gt.6.2d-7) &
!!$                  print * , 'OVERSHOT!!!'
!!$             ! overshot = .TRUE.
!!$             return
!!$          endif
!!$       endif
!!$    endif
!!$    overshot = .FALSE.
!!$
!!$  end function has_particle_overshot_resonance

end module RFOF_resonance_condition
