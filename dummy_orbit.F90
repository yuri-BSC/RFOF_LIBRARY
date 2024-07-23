#include "../src/config.h"

module dummy_orbit

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  use dum_magnetic_field
  use RFOF_magnetic_field
  use RFOF_local_magnetic_field, only: get_local_magnetic_field
  use RFOF_waves
  use RFOF_resonance_memory
  use RFOF_markers
  use RFOF_main
  use euitm_schemas
  use euitm_waves_interface
  use RFOF_diagnostics
  use RFOF_constants
  use RFOF_parameters
  use RFOF_wiener_sample_paths

  implicit none

  logical :: theta_step_constant = .true.
  real(8) :: d_theta = 1e-2
contains

  !--------------------------------------------------------------------------------
  subroutine run_dumorb

    use RFOF_Efield_update

    ! Local
    type(rf_wave_global) :: RFglobal
    type(particle) :: marker
    type(particle_static), target :: marker_static
    type(resonance_memory), pointer :: mem(:,:)
    type(RFOF_cumlative_diagnostics) :: diagno
    logical :: interaction_failed_particle_overshot_resonance
    real(8) :: time_at_resonance
    integer :: jtime, NtimeSteps
    real(8) :: time, dt, dtACC, time_at_last_Erenorm
    integer :: MPI_nod_Id = 0
    real(8), target :: NACC = 1
    integer :: Nstep_for_updating_Efield = 20000
    integer :: ierr ! Error flag

    !call test_sample_path

    print *, "init_dumorb"
    call init_dumorb(RFglobal,marker,marker_static,mem,diagno,NtimeSteps,dt,NACC)
    
    time=0.
    time_at_last_Erenorm = time
    do jtime=1,NtimeSteps
       !print *, "  --  step_orbit_dumorb -- ", time, dt

       !---------------------------
       ! START LOOP over markers

       ! Time stepping
       NACC = real(mod(jtime , 3)+1)  ! For testing acceleration; can we find resonances also with acceleration?
       dtACC=dt*NACC
       !write(0,*)"nacc=",nacc,dtACC,time
       call step_orbit_dumorb(marker,dtACC,NACC)
       time = time + dtACC

       ! Add RF kicks
       call RFOF_master(time,time-1e-6, &
            RFglobal,marker,mem,diagno,MPI_nod_Id, &
            interaction_failed_particle_overshot_resonance, &
            time_at_resonance)

       ! END LOOP over markers
       !---------------------------


       if ( (mod(jtime , Nstep_for_updating_Efield) == 0) .or. (jtime.eq.NtimeSteps) ) then
          call output_heating_to_screen(RFglobal, diagno, time)

          call store2file_RFOF_cumlative_diagnostics(diagno,RFglobal)
          call update_efield_normalisation(time - time_at_last_Erenorm,RFglobal,MPI_nod_Id,ierr,diagno=diagno)
          time_at_last_Erenorm = time
       endif

       if (interaction_failed_particle_overshot_resonance) then
          write(0,*)'...'
          print *, 'kick failed: t_res, t',time_at_resonance, time
       endif
    enddo
    write(0,*)'loop done'

    !call store2file_RFOF_cumlative_diagnostics(diagno,RFglobal)

    print *, "exit_dumorb"
    call exit_dumorb(RFglobal,mem,diagno)

    !call test_ridder

  end subroutine run_dumorb


  !--------------------------------------------------------------------------------
  ! subroutine init_dumorb
  !--------------------------------------------------------------------------------
  subroutine init_dumorb(RFglobal,marker,marker_static,mem,diagno,NtimeSteps,dt,time_acceleration)

    ! Output
    type(rf_wave_global), intent(out) :: RFglobal
    type(particle), intent(inout) :: marker
    type(particle_static), target, intent(inout) :: marker_static
    type(resonance_memory), pointer, intent(out) :: mem(:,:)
    type(RFOF_cumlative_diagnostics), intent(inout) :: diagno
    integer, intent(out) ::NtimeSteps
    real(8), intent(out) :: dt
    real(8), target, intent(in) :: time_acceleration

    type(type_waves) :: itm_waves

    ! Local
    real(8) :: R0, B0, q, aminor
    real(8) :: mass = 2.0
    real(8) :: weight, charge, R,phi,z,E,xi,tauBounce
    type(magnetic_field_local) :: Blocal
    real(8) :: Rmin
    real(8) :: Rmax
    real(8) :: Zmin
    real(8) :: Zmax

    NAMELIST/input_control/dt, NtimeSteps
    NAMELIST/input_magnetifield/R0, aminor, B0, q
    NAMELIST/input_marker/weight,R,phi,z,mass,charge,E,xi
    NAMELIST/simplify_rfof/ &
         simplify__static_resonance_position_during_RF_kick , &
         simplify__drift_velocity_does_not_affect_resonance_condition , &
         simplify__parallel_velocity_does_not_affect_resonance_condition , &
         simplify__assume_zero_larmor_radius_in_KPERPxRHO , &
         simplify__kpar_is_nphi_over_R
    NAMELIST/IO_control/&
         start_time_event_output , &
         output__2D_RZ_out , &
         output__Orbit , &
         MAX_number_of_points_stored_in_the_Orbit , &
         output__rf_kicks , &
         MAX_number_of_points_stored_in_the_rf_kick , &
         output__resonace_predictions , &
         MAX_number_of_points_stored_in_the_resonance_memory
    NAMELIST/rz_boundingbox/Rmin, Rmax, Zmin, Zmax

    ! Default:
    dt = 1e-4
    NtimeSteps = 10
!    time_acceleration = NACC

    open( io_channel_3872, FILE='input.rfof')
    read( io_channel_3872, input_control)
    read( io_channel_3872, input_magnetifield)
    read( io_channel_3872, input_marker)
    !read( io_channel_3872, simplify_rfof)
    !read( io_channel_3872, IO_control)
    !read( io_channel_3872, rz_boundingbox)
    close(io_channel_3872)

    ! Set plasma bounding box in R,Z
    call initialise_RFOF_parameters()

    ! Get B-field
    print *, " ---Init B field---"
    call magnetic_field_constructor_simple_field(R0,B0,q,aminor,Bglobal)
    !print *, " B0=", Bglobal%B0
    !print *, " R0=", Bglobal%R0
    ! Get wave field
    print *, 'mass:', mass 
    print *, " ---Init Wave field---"
    call dummy_rf_wave_field(RFglobal)

    call construct_dummy_euitm_waves(itm_waves)


    ! Get resonance memory
    print *, " ---Init resonance memory---"
    call constructor_rf_resonance_memory_matrix(mem,RFglobal%nfreq,RFglobal%max_nnphi)

    ! Get marker
    print *, " ---Init marker---"
    E = E * 1.6022e-19
    print *, "la massa es " , mass
    Blocal = get_local_magnetic_field(R,phi,z, marker)
    print *, "test-1"
    tauBounce = R0 / sqrt(2*E/1.66e-27)

    call make_marker(marker_static, weight,charge,mass,E,xi,tauBounce,Blocal)
    !call set_marker_pointers_from_marker(marker_static , marker)
    call set_marker_pointers(marker,   marker_static%Id,              marker_static%weight, &
         marker_static%R,              marker_static%phi,             marker_static%z,  &
         marker_static%charge,         marker_static%mass,            marker_static%energy, &
         marker_static%energy_kinetic, marker_static%velocity,        marker_static%magneticMoment, &
         marker_static%Pphi,           marker_static%vpar,            marker_static%vperp, &
         marker_static%omega_gyro,     marker_static%tauBounce,       marker_static%vDrift, &
         marker_static%vDriftRho,      marker_static%vDriftDia,       marker_static%d_vpar_d_rho, &
         marker_static%d_vpar_d_dia,   marker_static%d_vperp_d_rho,   marker_static%d_vperp_d_dia, &
         marker_static%d_vDriftRho_d_rho, marker_static%d_vDriftRho_d_dia, &
         marker_static%d_vDriftDia_d_rho, marker_static%d_vDriftDia_d_dia, &
         time_acceleration)

    ! Get diagnostics
    print *, " ---Init diagnostics---"
    call contructor_RFOF_cumlative_diagnostics(diagno,RFglobal%nfreq,RFglobal%max_nnphi)

  end subroutine init_dumorb

  !--------------------------------------------------------------------------------
  ! subroutine exit_dumorb
  !--------------------------------------------------------------------------------
  subroutine exit_dumorb(RFglobal,mem,diagno)

    ! Input/Output
    type(rf_wave_global), intent(inout) :: RFglobal
    type(resonance_memory), pointer, intent(out) :: mem(:,:)
    type(RFOF_cumlative_diagnostics), intent(inout) :: diagno

    call RFOF_destructor(RFglobal,mem,diagno)

  end subroutine exit_dumorb

  !--------------------------------------------------------------------------------
  ! subroutine step_orbit_dumorb
  !--------------------------------------------------------------------------------
  subroutine step_orbit_dumorb(marker,dt,NACC)

    ! Input
    real(8), intent(in) :: dt
    real(8), intent(in) :: NACC

    ! Input/Output
    type(particle), intent(inout) :: marker

    ! Local
    real(8) :: Bphi, Bpol, Bmod
    real(8) :: vpol, vphi
    real(8) :: rho, theta0, theta, R, phi, z
    real(8) :: dtORB
    type(magnetic_field_local) :: Blocal

    dtORB = (dt/NACC) ! Non-accelerated time step - orbit time step

    rho = sqrt( (marker%R-Bglobal%R0)**2 + marker%z**2 + 1e-5 )

!    rho = 0.30166
    theta0 = sign(1d0,marker%z)*acos(( marker%R-Bglobal%R0)/rho)

    ! Local B-field before kick
    Blocal = get_local_magnetic_field(marker%R,marker%phi,marker%z)

!    print *, "-"

    Bphi = Blocal%F / marker%R
    Bpol = Bglobal%B0 * rho / ( Bglobal%q * marker%R )
    Bmod = Blocal%Bmod

    vpol = marker%vpar * Bpol / Blocal%Bmod
    vphi = marker%vpar * Bphi / Blocal%Bmod

    if ( theta_step_constant ) then
       theta = theta0 + d_theta
       marker%tauBounce = dtORB * 2.*3.141592/d_theta
    else
       theta = theta0 + abs(vpol * dtORB / rho) + 1e-5
    endif
    phi = Blocal%phi + vphi * dtORB / Blocal%R
    R = Bglobal%R0 + rho*cos(theta)
    z =              rho*sin(theta)

    ! Local B-field after kick
    Blocal = get_local_magnetic_field(R,phi,z)
    call update_marker(marker,Blocal)

  end subroutine step_orbit_dumorb


  !--------------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------------
  subroutine output_heating_to_screen(RFglobal, diagno, time)

    type(rf_wave_global) :: RFglobal
    type(RFOF_cumlative_diagnostics) :: diagno
    real(8) :: time

    real(8) :: dt_diag

    if (sum(diagno%kick_counter) == 0) then
       write(0,*)" "
       write(0,*)" NO KICK IN THE CUMULATIVE DIAGNO, TIME=",time
       write(0,*)" "
       return
    endif

    dt_diag =  max(1d-40 , (diagno%time_end_sum_diagnostics - diagno%time_start_sum_diagnostics))

    write(0,*)" "
    print *, "Power  transfered to markers [W]    =", diagno%energy_to_markers / dt_diag
    print *, "Torque transfered to markers [Nm/s] =", diagno%toroidal_momentum_to_markers / dt_diag
    write(0,*)" "
    print *, "Energy transfered to markers [J]  =", diagno%energy_to_markers
    print *, "Pphi   transfered to markers [Nm] =", diagno%toroidal_momentum_to_markers

    if (RFglobal%nfreq .eq. 1 .and. RFglobal%max_nnphi .eq. 1) then
       print *, "               nphi/omega dE [Nm] =", &
            diagno%energy_to_markers * RFglobal%waves%coherentwave(1)%global_param%ntor(1) / &
            (2d0*rfof_pi*RFglobal%waves%coherentwave(1)%global_param%frequency)
    endif
    write(0,*)" "

  end subroutine output_heating_to_screen


  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  function anyfunc(x) result(y)
    real(8) :: x
    real(8) :: y    
    y=x**5
  end function anyfunc

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine test_ridder

    real(8) :: x, h, rel_tol, err, fprim

    !interface
    !   function anyfunc(x) result(y)
    !     real(8) :: x
    !     real(8) :: y
    !   end function anyfunc
    !end interface

    x=1.d2
    h=0.5
    rel_tol = 1e-6

    call Ridders_derivative(anyfunc, x, h, rel_tol, err, fprim)

    print *,  x, err, fprim, err/fprim

  end subroutine test_ridder

end module dummy_orbit
