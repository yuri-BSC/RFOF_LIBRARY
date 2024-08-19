#include "config.h"

module RFOF_diagnostics

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  use RFOF_parameters

  implicit none

  type :: RFOF_cumlative_diagnostics
     logical :: diagnostic_active = .false.

     integer, allocatable :: kick_counter(:,:)

     real(8) :: time_start_sum_diagnostics

     real(8) :: time_end_sum_diagnostics

     real(8) :: energy_to_markers = 0

     real(8) :: toroidal_momentum_to_markers = 0

     real(8), allocatable :: energy_to_markers_from_mode(:,:)

     real(8), allocatable :: energy_square_to_markers_from_mode(:,:)

     real(8), allocatable :: toroidal_momentum_to_markers_from_mode(:,:)

     real(8), allocatable :: sum_weight_at_resonance_with_mode(:,:)

     real(8), allocatable :: R_2DgridRZ(:,:)

     real(8), allocatable :: Z_2DgridRZ(:,:)

     real(8), allocatable :: energy_2DgridRZ(:,:)

     real(8), allocatable :: momentum_2DgridRZ(:,:)

  end type RFOF_cumlative_diagnostics

contains

  !--------------------------------------------------------------------------------
  !  subroutine contructor_RFOF_cumlative_diagnostics(diagno,nfreq,nnphi)
  subroutine contructor_RFOF_cumlative_diagnostics(diagno,nfreq,nnphi)

    ! Input
    integer, intent(in) :: nfreq, nnphi

    ! Input/Output
    type(RFOF_cumlative_diagnostics), intent(inout) :: diagno

    ! Local
    integer :: jf, jn
    integer :: jr, jz

    call read_diagnostic_parameters_from_namelist_file



    if (.not. ALLOCATED(diagno%kick_counter) ) then
       allocate(diagno%kick_counter(nfreq,nnphi))
    endif
    allocate(diagno%energy_to_markers_from_mode(nfreq,nnphi))
    allocate(diagno%energy_square_to_markers_from_mode(nfreq,nnphi))
    allocate(diagno%toroidal_momentum_to_markers_from_mode(nfreq,nnphi))
    allocate(diagno%sum_weight_at_resonance_with_mode(nfreq,nnphi))

    diagno%kick_counter(:,:) = 0
    diagno%energy_to_markers_from_mode(:,:) = 0d0
    diagno%energy_square_to_markers_from_mode(:,:) = 0d0
    diagno%toroidal_momentum_to_markers_from_mode(:,:) = 0d0
    diagno%sum_weight_at_resonance_with_mode(:,:) = 0d0

    if (output__2D_RZ_out) then

       print *, "Allocating memory for 2D-RZ grid. Dimension:",&
            NRedges_2DgridRZ , NZedges_2DgridRZ

       allocate( diagno%R_2DgridRZ(        NRedges_2DgridRZ  , NZedges_2DgridRZ   ) )
       allocate( diagno%Z_2DgridRZ(        NRedges_2DgridRZ  , NZedges_2DgridRZ   ) )
       allocate( diagno%energy_2DgridRZ(   NRedges_2DgridRZ-1, NZedges_2DgridRZ-1 ) )
       allocate( diagno%momentum_2DgridRZ( NRedges_2DgridRZ-1, NZedges_2DgridRZ-1 ) )
       
       diagno%energy_2DgridRZ(:,:)   = 0d0
       diagno%momentum_2DgridRZ(:,:) = 0d0

       do jr = 1, NRedges_2DgridRZ
          do jz = 1, NZedges_2DgridRZ
             diagno%R_2DgridRZ(jr,jz) = plasma_boundingbox%Rmin + &
                  (plasma_boundingbox%Rmax - plasma_boundingbox%Rmin) * &
                  dble(jr-1)/dble(NRedges_2DgridRZ-1)
             diagno%Z_2DgridRZ(jr,jz) = plasma_boundingbox%Zmin + &
                  (plasma_boundingbox%Zmax - plasma_boundingbox%Zmin) * &
                  dble(jz-1)/dble(NZedges_2DgridRZ-1)
          enddo
       enddo
    endif

  end subroutine contructor_RFOF_cumlative_diagnostics

  !--------------------------------------------------------------------------------
  ! subroutine destructor_RFOF_cumlative_diagnostics(diagno)
  subroutine destructor_RFOF_cumlative_diagnostics(diagno)

    ! Input/Output
    type(RFOF_cumlative_diagnostics), intent(inout) :: diagno

    if ( allocated(diagno%kick_counter) ) then
       deallocate(diagno%kick_counter)
    endif
    if ( allocated(diagno%energy_to_markers_from_mode) ) then
       deallocate(diagno%energy_to_markers_from_mode)
    endif
    if ( allocated(diagno%energy_square_to_markers_from_mode) ) then
       deallocate(diagno%energy_square_to_markers_from_mode)
    endif
    if ( allocated(diagno%toroidal_momentum_to_markers_from_mode) ) then
       deallocate(diagno%toroidal_momentum_to_markers_from_mode)
    endif
    if ( allocated(diagno%sum_weight_at_resonance_with_mode) ) then
       deallocate(diagno%sum_weight_at_resonance_with_mode)
    endif

    if (output__2D_RZ_out) then
       if ( allocated(diagno%R_2DgridRZ) ) then
          deallocate( diagno%R_2DgridRZ )
       endif
       if ( allocated(diagno%Z_2DgridRZ) ) then
          deallocate( diagno%Z_2DgridRZ )
       endif
       if ( allocated(diagno%energy_2DgridRZ) ) then
          deallocate( diagno%energy_2DgridRZ )
       endif
       if ( allocated(diagno%momentum_2DgridRZ) ) then
          deallocate( diagno%momentum_2DgridRZ )
       endif
    endif

  end subroutine destructor_RFOF_cumlative_diagnostics


  !--------------------------------------------------------------------------------
  ! subroutine read_diagnostic_parameters_from_namelist_file()
  subroutine read_diagnostic_parameters_from_namelist_file()

    NAMELIST/IO_control/&
         start_time_event_output , &
         output__2D_RZ_out , &
         output__Orbit , &
         MAX_number_of_points_stored_in_the_Orbit , &
         output__rf_kicks , &
         MAX_number_of_points_stored_in_the_rf_kick , &
         output__resonace_predictions , &
         MAX_number_of_points_stored_in_the_resonance_memory

    open( io_channel_3872, FILE='input.rfof')
    read( io_channel_3872, IO_control)
    close(io_channel_3872)

  end subroutine read_diagnostic_parameters_from_namelist_file


  !--------------------------------------------------------------------------------
  ! subroutine reset_to_zero_RFOF_cumlative_diagnostics(diagno)
  subroutine reset_to_zero_RFOF_cumlative_diagnostics(diagno)

    ! Input/Output
    type(RFOF_cumlative_diagnostics), intent(inout) :: diagno

    diagno%energy_to_markers_from_mode(:,:) = 0d0
    diagno%energy_square_to_markers_from_mode(:,:) = 0d0
    diagno%toroidal_momentum_to_markers_from_mode(:,:) = 0d0
    diagno%sum_weight_at_resonance_with_mode(:,:) = 0d0

    if (output__2D_RZ_out) then
       diagno%energy_2DgridRZ(:,:)   = 0d0
       diagno%momentum_2DgridRZ(:,:) = 0d0
    endif

    diagno%kick_counter(:,:) = 0
    diagno%energy_to_markers = 0d0
    diagno%toroidal_momentum_to_markers = 0d0

    diagno%time_start_sum_diagnostics = -very_large_time_in_seconds
    diagno%time_end_sum_diagnostics   = -very_large_time_in_seconds

    diagno%diagnostic_active = .false.

  end subroutine reset_to_zero_RFOF_cumlative_diagnostics

  !--------------------------------------------------------------------------------
  ! subroutine add_kick_info_to_diagnostics
  subroutine add_kick_info_to_diagnostics(diagno, marker, mem, &
       jfreq, jnphi, dE, dPphi, RFop1, RFop2, RFopType, MPI_node_Id)

    use RFOF_markers, only: particle
    use RFOF_resonance_memory, only: resonance_memory

    ! Input
    type(particle), intent(inout) :: marker
    type(resonance_memory), intent(in) :: mem
    integer, intent(in) :: jfreq, jnphi
    real(8), intent(in) :: dE, dPphi
    real(8), intent(in) :: RFop1, RFop2
    character(128), intent(in) :: RFopType
    integer, intent(in) :: MPI_node_Id

    ! Input/Output
    type(RFOF_cumlative_diagnostics), intent(inout) :: diagno

    ! Local
    integer :: nR, nZ, iR, iZ
    real(8) :: Rnorm, Znorm

    diagno%kick_counter(jfreq,jnphi) = diagno%kick_counter(jfreq,jnphi) + 1

    ! Global diagnostics
    diagno%sum_weight_at_resonance_with_mode = diagno%sum_weight_at_resonance_with_mode + &
         marker%weight

    diagno%energy_to_markers = diagno%energy_to_markers + &
         marker%weight * dE

    diagno%toroidal_momentum_to_markers = &
         diagno%toroidal_momentum_to_markers + &
         marker%weight * dPphi

    diagno%energy_to_markers_from_mode(jfreq,jnphi) = &
         diagno%energy_to_markers_from_mode(jfreq,jnphi) + &
         marker%weight * dE

    diagno%energy_square_to_markers_from_mode(jfreq,jnphi) = &
         diagno%energy_square_to_markers_from_mode(jfreq,jnphi) + &
         marker%weight * dE**2

    diagno%toroidal_momentum_to_markers_from_mode(jfreq,jnphi) = &
         diagno%toroidal_momentum_to_markers_from_mode(jfreq,jnphi) + &
         marker%weight * dPphi

    if (output__2D_RZ_out) then

       ! Local 2D diagnostics
       nR=size(diagno%R_2DgridRZ,1)
       nZ=size(diagno%R_2DgridRZ,2)

       Rnorm = (marker%R - diagno%R_2DgridRZ(1,1)) / &
            (diagno%R_2DgridRZ(nR,1) - diagno%R_2DgridRZ(1,1))

       Znorm = (marker%z - diagno%Z_2DgridRZ(1,1)) / &
            (diagno%Z_2DgridRZ(1,nZ) - diagno%Z_2DgridRZ(1,1))

       iR = int( Rnorm * dble(nR) + 1d0)
       iZ = int( Znorm * dble(nR) + 1d0)

       if (iR.gt.0 .and. iR.lt.nR .and. iZ.gt.0 .and. iZ.lt.nZ) then

          diagno%energy_2DgridRZ(ir,iz)   = &
               diagno%energy_2DgridRZ(ir,iz)   + dE * marker%weight
          diagno%momentum_2DgridRZ(ir,iz) = &
               diagno%momentum_2DgridRZ(ir,iz) + dPphi * marker%weight

       endif
    endif

    !if ( (output__rf_kicks) .and. &
    !     (time.gt.start_time_event_output) ) then
    if (output__rf_kicks) then
       ! Output kick data to file
       call save2file_rf_kicks(marker, mem, &
            jfreq, jnphi, dE, dPphi, RFop1, RFop2, RFopType,MPI_node_Id)
    endif

  end subroutine add_kick_info_to_diagnostics


  !--------------------------------------------------------------------------------
  ! subroutine store2file_RFOF_cumlative_diagnostics(diagno)
  subroutine store2file_RFOF_cumlative_diagnostics(diagno,RFglobal)

    use RFOF_waves, only: rf_wave_global

    ! Input
    type(RFOF_cumlative_diagnostics), intent(in) :: diagno
    type(rf_wave_global), intent(in) :: RFglobal

    ! Local
    integer :: jf, jn, ir, iz
    real(8) :: dt
    integer :: io_channel


  !--------------------------------------------------------------------------------
    write(0,*) "Store RFOF output to file"
    io_channel = io_channel_3872
    dt = diagno%time_end_sum_diagnostics - diagno%time_start_sum_diagnostics
    if ( dt < 1d-40 ) then
       write(0,*)'WARNING in store2file_RFOF_cumlative_diagnostics: Diagnostic integration time = ',dt
       write(0,*)'        diagno%time_end_sum_diagnostics   =',diagno%time_end_sum_diagnostics
       write(0,*)'        diagno%time_start_sum_diagnostics =',diagno%time_start_sum_diagnostics
       dt=1d-40
    endif
    open(file = "out.RFOF_integrate_1D", unit = io_channel)

    write(io_channel,*) "RFOF output" 
    write(io_channel,*) "",RFOF_version
    write(io_channel,*) "Time integrated 1D diagnostics. No. kicks=" , diagno%kick_counter
    write(io_channel,*) " "

    write(io_channel,*) "Number_of_frequencies = ", RFglobal%nfreq
    write(io_channel,'(A,1000E15.5)') " Frequencies = ", &
         (RFglobal%waves%coherentwave(jf)%global_param%frequency, jf=1,RFglobal%nfreq)

    write(io_channel,*) " " 
    write(io_channel,*) "% Number of toroidal modes for each frequency; (jfreq , nphi(jfreq)) : "
    do jf=1,RFglobal%nfreq
       write(io_channel,*) jf, size(RFglobal%waves%coherentwave(jf)%global_param%ntor)
    enddo
    write(io_channel,*) "% Toroidal modes numbers for each frequency; (jfreq , jnphi, nphi(jfreq,jnphi)) : "
    do jf=1,RFglobal%nfreq
       do jn=1,RFglobal%nnphi(jf)
          write(io_channel,*) jf , jn, RFglobal%waves%coherentwave(jf)%global_param%ntor(jn)
       enddo
    enddo

    write(io_channel,*) " "
    write(io_channel,*) "Length of diagnostic time integration interval = ",dt

    write(io_channel,*) " "
    write(io_channel,'(A,E15.8)') "Total power absorbed  [W]    = ", diagno%energy_to_markers
    write(io_channel,'(A,E15.8)') "Total torque absorbed [Nm/s] = ", diagno%toroidal_momentum_to_markers

    write(io_channel,*) " "
    write(io_channel,*) "% Total absorbed power and torque, for each frequency and toroidal mode; "
    write(io_channel,*) "%    (jfreq , jnphi, POWER(jfreq,jnphi), TORQUE(jfreq,jnphi)) : "
    do jf=1,RFglobal%nfreq
       do jn=1,RFglobal%nnphi(jf)
          write(io_channel,'(2I7,2E15.8)') jf , jn, &
               diagno%energy_to_markers_from_mode(jf,jn) / dt, &
               diagno%toroidal_momentum_to_markers_from_mode(jf,jn) / dt
       enddo
    enddo

    close(unit = io_channel)

  !--------------------------------------------------------------------------------

    if (output__2D_RZ_out) then
       write(0,*) "Store 2D cumulative diagnostics output to file"
       dt = max(diagno%time_end_sum_diagnostics - diagno%time_start_sum_diagnostics , 1d-40)
       
       io_channel = io_channel_3872
       open(file = "out.RFOF_integrate_2DRZ", unit = io_channel)
       
       write(io_channel,*) "RFOF output" 
       write(io_channel,*) "",RFOF_version
       write(io_channel,*) "Time integrated 2D-RZ diagnostics" 
       write(io_channel,*) " "
       
       write(io_channel,*) "Number_of_R_grid_cells = ", size(diagno%R_2DgridRZ,1)-1
       write(io_channel,*) "Number_of_Z_grid_cells = ", size(diagno%R_2DgridRZ,2)-1
       
       do ir = 1, size(diagno%R_2DgridRZ,1)-1
          do iz = 1, size(diagno%Z_2DgridRZ,2)-1
             write(io_channel,*)ir,iz, &
                  0.5d0 * ( diagno%R_2DgridRZ(ir,iz) + diagno%R_2DgridRZ(ir+1,iz+1) ), &
                  0.5d0 * ( diagno%Z_2DgridRZ(ir,iz) + diagno%Z_2DgridRZ(ir+1,iz+1) ), &
                  diagno%energy_2DgridRZ(ir,iz), &
                  diagno%momentum_2DgridRZ(ir,iz)
          enddo
       enddo

       close(unit = io_channel)
    endif

  end subroutine store2file_RFOF_cumlative_diagnostics


  !--------------------------------------------------------------------------------
  ! subroutine store2file_RFOF_orbit_info
  subroutine store2file_RFOF_orbit_info(marker,time,timestep,MPI_node_Id)

    use RFOF_markers, only: particle
    use RFOF_magnetic_field, only: magnetic_field_local 
    use RFOF_local_magnetic_field, only: get_local_magnetic_field

    ! Input
    type(particle), intent(in) :: marker
    real(8), intent(in) :: time
    real(8), intent(in) :: timestep
    integer, intent(in) :: MPI_node_Id

    ! Local
    real(8) :: data(17)
    character(128) :: filename, calling_routine_name, data_description, list_variable_names
    type(magnetic_field_local), target :: Blocal

    ! Get the magnetic field values at the location of the test particle
    Blocal = get_local_magnetic_field(marker%R,marker%phi,marker%z, marker%mass)

    data( 1) = time
    data( 2) = marker%weight
    data( 3) = marker%mass
    data( 4) = marker%charge
    data( 5) = marker%R
    data( 6) = marker%phi
    data( 7) = marker%z
    data( 8) = Blocal%psi
    data( 9) = Blocal%theta
    data(10) = marker%energy
    data(11) = marker%vpar
    data(12) = marker%vperp
    data(13) = marker%magneticMoment
    data(14) = marker%Pphi
    data(15) = marker%omega_gyro
    data(16) = marker%tauBounce
    data(17) = timestep

    if (MPI_node_Id.lt.10) then
       write(filename,'(A,I1)') 'out.orbits_',MPI_node_Id
    elseif (MPI_node_Id.lt.100) then
       write(filename,'(A,I2)') 'out.orbits_',MPI_node_Id
    elseif (MPI_node_Id.lt.1000) then
       write(filename,'(A,I3)') 'out.orbits_',MPI_node_Id
    elseif (MPI_node_Id.lt.10000) then
       write(filename,'(A,I4)') 'out.orbits_',MPI_node_Id
    elseif (MPI_node_Id.lt.100000) then
       write(filename,'(A,I5)') 'out.orbits_',MPI_node_Id
    else
       write(filename,*) 'out.orbits_',MPI_node_Id
    endif
    calling_routine_name = 'store2file_RFOF_orbit_info'
    data_description = 'Orbits'
    list_variable_names = 'time,weight,mass,charge,R,phi,z,psi,theta,energy,vpar,vperp,mu,Pphi,omega_c,tauBounce,timestep'
    !                      1    2      3    4      5 6   7 8   9     10     11   12    13 14   15      16        17

    call store2file_orbit_event( &
         number_of_points_stored_in_the_Orbit, &
         MAX_number_of_points_stored_in_the_Orbit, &
         filename, &
         calling_routine_name, &
         data_description, &
         list_variable_names, &
         data, &
         io_channel_3875)

  end subroutine store2file_RFOF_orbit_info


  !--------------------------------------------------------------------------------
  ! subroutine store2file_RFOF_efield_normalization
  subroutine save2file_RFOF_efield_normalization(time,Pabs,wave,MPI_node_Id)

    use RFOF_waves

    ! Input
    type(rf_wave_global), intent(in) :: wave
    real(8), intent(in) :: time
    real(8), intent(in) :: Pabs(:,:)
    !real(8), intent(in) :: Enorm
    integer, intent(in) :: MPI_node_Id

    ! Local
    real(8), allocatable :: data(:)
    character(128) :: filename
    character(128) :: calling_routine_name = 'save2file_RFOF_efield_normalization'
    character(128) :: data_description = 'E-field normalization'
    character(128) :: list_variable_names = 'time,Power_absorption,E_normalisation_factor'
    integer :: data_size, jcount
    integer :: jf, jn

    data_size=1
    do jf = 1, wave%nfreq
       data_size = data_size + 2*wave%nnphi(jf)
    enddo
    allocate( data(data_size) )

    data( 1 ) = time
    jcount = 1

    list_variable_names = 'time'
    do jf = 1, wave%nfreq
       do jn = 1, wave%nnphi(jf)
          write(list_variable_names,'(2A,I1,A,I1)') trim(list_variable_names),',Power_',jf,'_',jn
          write(list_variable_names,'(2A,I1,A,I1)') trim(list_variable_names),',Enorm_',jf,'_',jn

          jcount = jcount + 1
          data( jcount ) = Pabs(jf,jn)
          jcount = jcount + 1
          data( jcount ) = wave%EnormalisationFactor(jf,jn)
       enddo
    enddo

    if (MPI_node_Id.lt.10) then
       write(filename,'(A,I1)') 'out.enorm_',MPI_node_Id
    elseif (MPI_node_Id.lt.100) then
       write(filename,'(A,I2)') 'out.enorm_',MPI_node_Id
    elseif (MPI_node_Id.lt.1000) then
       write(filename,'(A,I3)') 'out.enorm_',MPI_node_Id
    elseif (MPI_node_Id.lt.10000) then
       write(filename,'(A,I4)') 'out.enorm_',MPI_node_Id
    elseif (MPI_node_Id.lt.100000) then
       write(filename,'(A,I5)') 'out.enorm_',MPI_node_Id
    else
       write(filename,*) 'out.enorm_',MPI_node_Id
    endif

    call store2file_orbit_event( &
         number_of_points_stored_in_the_efield_normalization, &
         MAX_number_of_points_stored_in_the_efield_normalization, &
         filename, &
         calling_routine_name, &
         data_description, &
         list_variable_names, &
         data, &
         io_channel_3879)

    deallocate( data )

  end subroutine save2file_RFOF_efield_normalization


  !--------------------------------------------------------------------------------
  subroutine save2file_resonace_predictions(time,R,phi,z,omega_res,nharm, &
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

    real(8), intent(in) :: time_at_resonance
    real(8), intent(in) :: R_at_resonance
    real(8), intent(in) :: phi_at_resonance
    real(8), intent(in) :: z_at_resonance
    real(8), intent(in) :: Dot_omega_res
    real(8), intent(in) :: DotDot_omega_res

    ! Local
    real(8) :: data(8)
    character(128) :: filename, calling_routine_name, data_description, list_variable_names
!    integer :: io_channel
!    character(128) filename

    data( 1) = time
    data( 2) = R
    data( 3) = phi
    data( 4) = z
    data( 5) = time_at_resonance
    data( 6) = R_at_resonance
    data( 7) = phi_at_resonance
    data( 8) = z_at_resonance

    if (MPI_node_Id.lt.10) then
       write(filename,'(A,I1)') 'out.resonance_',MPI_node_Id
    elseif (MPI_node_Id.lt.100) then
       write(filename,'(A,I2)') 'out.resonance_',MPI_node_Id
    elseif (MPI_node_Id.lt.1000) then
       write(filename,'(A,I3)') 'out.resonance_',MPI_node_Id
    elseif (MPI_node_Id.lt.10000) then
       write(filename,'(A,I4)') 'out.resonance_',MPI_node_Id
    elseif (MPI_node_Id.lt.100000) then
       write(filename,'(A,I5)') 'out.resonance_',MPI_node_Id
    else
       write(filename,*) 'out.resonance_',MPI_node_Id
    endif
    calling_routine_name = 'save2file_resonace_predictions'
    data_description = 'Resonance predictions'
    list_variable_names = 'time, R, phi, z, time_res, R_res, phi_res, z_res'

    call store2file_orbit_event( &
         number_of_points_stored_in_the_resonance_memory, &
         MAX_number_of_points_stored_in_the_resonance_memory, &
         filename, &
         calling_routine_name, &
         data_description, &
         list_variable_names, &
         data, &
         io_channel_3874)

  end subroutine save2file_resonace_predictions


  !--------------------------------------------------------------------------------
  subroutine save2file_rf_kicks(marker, mem, &
       jfreq, jnphi, dE, dPphi, RFop1, RFop2, RFopType,MPI_node_Id)

    use RFOF_markers, only: particle
    use RFOF_resonance_memory, only: resonance_memory

    ! Input
    type(particle), intent(inout) :: marker
    type(resonance_memory), intent(in) :: mem
    integer, intent(in) :: jfreq, jnphi
    real(8), intent(in) :: dE, dPphi
    real(8), intent(in) :: RFop1, RFop2
    character(128), intent(in) :: RFopType
    integer, intent(in) :: MPI_node_Id

    ! Local
    real(8) :: tres, fres
    integer :: j
    real(8) :: data(14)
    character(128) :: filename, calling_routine_name, data_description, list_variable_names

    fres=mem%omega_res(1)
    tres=mem%time(1)
    do j=1,mem%Number_points_in_memory
       if (abs(mem%omega_res(j)) .lt. fres) then
          fres=mem%omega_res(j)
          tres=mem%time(j)
       endif
    enddo

    data( 1) = tres
    data( 2) = marker%weight
    data( 3) = marker%mass
    data( 4) = marker%charge

    data( 5) = marker%R
    data( 6) = marker%z
    data( 7) = marker%vpar
    data( 8) = marker%vperp
    data( 9) = dble(jfreq)
    data(10) = dble(jnphi)
    data(11) = dE
    data(12) = dPphi
    data(13) = RFop1
    data(14) = RFop2

    if (MPI_node_Id.lt.10) then
       write(filename,'(A,I1)') 'out.rfkicks_',MPI_node_Id
    elseif (MPI_node_Id.lt.100) then
       write(filename,'(A,I2)') 'out.rfkicks_',MPI_node_Id
    elseif (MPI_node_Id.lt.1000) then
       write(filename,'(A,I3)') 'out.rfkicks_',MPI_node_Id
    elseif (MPI_node_Id.lt.10000) then
       write(filename,'(A,I4)') 'out.rfkicks_',MPI_node_Id
    elseif (MPI_node_Id.lt.100000) then
       write(filename,'(A,I5)') 'out.rfkicks_',MPI_node_Id
    else
       write(filename,*) 'out.orbits_',MPI_node_Id
    endif
    !write(calling_routine_name,'(2A)') 'save2file_rf_kicks/', RFopType(1:109)
    calling_routine_name = 'save2file_rf_kicks'
    !write(data_description,'(2A)')     'RF kicks from ', RFopType(1:114)
    !data_description = RFopType
    list_variable_names = 'time,R,z,vpar,vperp,jfreq,jnphi,dE,dPphi,RF1,RF2'

    calling_routine_name = 'save2file_rf_kicks'
    data_description = 'rf_kicks'
    list_variable_names = 'time,weight,mass,charge,R,z,vpar,vperp,jfreq,jnphi,dE,dPphi,RFop1,RFop2'

    call store2file_orbit_event( &
         number_of_points_stored_in_the_rf_kick, &
         MAX_number_of_points_stored_in_the_rf_kick, &
         filename, &
         calling_routine_name, &
         RFopType, &      !data_description, &
         list_variable_names, &
         data, &
         io_channel_3873)

  end subroutine save2file_rf_kicks


  !--------------------------------------------------------------------------------
  subroutine store2file_orbit_event( &
       no_points, &
       MAX_no_points, &
       filename, &
       calling_routine_name, &
       data_description, &
       list_variable_names, &
       data, &
       io_channel)

    ! Input/Output
    integer  , intent(inout) :: no_points

    ! Input
    integer  , intent(in) :: MAX_no_points
    character(128), intent(in) :: filename
    character(128), intent(in) :: calling_routine_name
    character(128), intent(in) :: data_description
    character(128), intent(in) :: list_variable_names
    real(8)  , intent(in) :: data(:)
    integer  , intent(in) :: io_channel

    ! Local
    integer :: j

    if (no_points .gt. MAX_no_points) then
       if (no_points .eq. MAX_no_points+1) then
          if (ignore_non_serious_error_messages .eqv. .FALSE.) then
             write(0,*) "------------------------------------------------------"
             write(0,*) "WARNING in ", calling_routine_name
             write(0,*) "  Number of points saved to file"
             write(0,*) "  has exceeded the limit of ", &
                  MAX_no_points
             write(0,*) "  No more points will be stored to file."
             write(0,*) "------------------------------------------------------"
             no_points = no_points + 1
          endif
       endif
       return
    endif

    if (no_points .eq. 0) then
       write(0,*)'Opening file:', filename
       open(file = filename, unit = io_channel)
       write(io_channel,*) 'RFOF ', RFOF_version
       write(io_channel,*) data_description
       write(io_channel,*) ' '
       write(io_channel,*) 'Nr_columns=',size(data)
       write(io_channel,*) ' '
       write(io_channel,*) 'Variables (1 per column):'
       write(io_channel,*) list_variable_names
       write(io_channel,*) ' '
    endif

    write(io_channel,'(200E15.7)')( data(j) , j=1,size(data) )

    no_points = no_points + 1

  end subroutine store2file_orbit_event


end module RFOF_diagnostics



!--------------------------------------------------------------------------------
subroutine wrap_save2file_resonace_predictions(time,R,phi,z,omega_res,nharm, &
       time_at_resonance,R_at_resonance,phi_at_resonance,z_at_resonance, &
       Dot_omega_res,DotDot_omega_res,MPI_node_Id)

  use RFOF_diagnostics

    ! Input
    real(8),    intent(in) :: time
    real(8),    intent(in) :: R
    real(8),    intent(in) :: phi
    real(8),    intent(in) :: z
    real(8),    intent(in) :: omega_res
    integer, intent(in) :: nharm
    integer, intent(in) :: MPI_node_Id

    real(8), intent(in) :: time_at_resonance
    real(8), intent(in) :: R_at_resonance
    real(8), intent(in) :: phi_at_resonance
    real(8), intent(in) :: z_at_resonance
    real(8), intent(in) :: Dot_omega_res
    real(8), intent(in) :: DotDot_omega_res

  call save2file_resonace_predictions(time,R,phi,z,omega_res,nharm, &
       time_at_resonance,R_at_resonance,phi_at_resonance,z_at_resonance, &
       Dot_omega_res,DotDot_omega_res,MPI_node_Id)

end subroutine wrap_save2file_resonace_predictions
