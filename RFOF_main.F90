#include "config.h"

module RFOF_main

#ifdef USE_ISO_C_BINDING
   use, intrinsic :: ISO_C_BINDING
#endif

   use RFOF_parameters
   use RFOF_waves
   use RFOF_magnetic_field
   use RFOF_markers
   use RFOF_resonance_memory
   use RFOF_resonance_condition
   use RFOF_numerics
   use RFOF_kick

   implicit none

   contains

   !--------------------------------------------------------------------------------
   ! subroutine RFOF_master(time,time_previous_step, &
   !     RFglobal,marker,resonanceMemory,diagno,MPI_nod_Id, &
   !     interaction_failed_particle_overshot_resonance, time_at_resonance)
   !
   subroutine RFOF_master(time,time_previous_step, &
         RFglobal,marker,resonanceMemory,diagno,MPI_nod_Id, &
         interaction_failed_particle_overshot_resonance, time_at_resonance)

      use RFOF_random_numbers

      ! Input
      real(8), intent(in) :: time
      real(8), intent(in) :: time_previous_step
      integer, intent(in) :: MPI_nod_Id

      ! Input/Output
      type(particle), intent(inout) :: marker
      type(resonance_memory), intent(inout) :: resonanceMemory(:,:)
      type(rf_wave_global), intent(inout) :: RFglobal
      type(RFOF_cumlative_diagnostics), intent(inout) :: diagno

      ! Output
      logical, intent(out) :: interaction_failed_particle_overshot_resonance
      real(8), intent(out) :: time_at_resonance

      ! Local variables
      integer :: nr_resonant_waves
      real(8), pointer :: list_resonant_waves(:,:)

      allocate(list_resonant_waves(RFglobal%max_nnphi*RFglobal%nfreq,3))
      nr_resonant_waves = 0
      time_at_resonance = time + very_large_time_in_seconds

      !--------------------------------------------------------------------------------
      ! Check resonance condition; is the marker at/has crossed the resonance?
      !--------------------------------------------------------------------------------
      call rf_wave_particles_resonance(time,time_previous_step, &
            marker,resonanceMemory,RFglobal, &
            nr_resonant_waves, list_resonant_waves, &
            time_at_resonance, &
            interaction_failed_particle_overshot_resonance,MPI_nod_Id)

      !--------------------------------------------------------------------------------
      ! If marker is in resonance (not overshot resonance), then give RF-Monte Carlo kick
      !--------------------------------------------------------------------------------
      if ( interaction_failed_particle_overshot_resonance .eqv. .FALSE.) then
         if ( nr_resonant_waves .gt. 0 ) then
            call split_RFstep_in_single_mode(marker,resonanceMemory,RFglobal, &
                  nr_resonant_waves,list_resonant_waves,diagno,MPI_nod_Id)
         endif
      endif
      deallocate(list_resonant_waves)

      !--------------------------------------------------------------------------------
      ! DIAGNOSTICS:
      !--------------------------------------------------------------------------------
      if ( .not. interaction_failed_particle_overshot_resonance ) then
         if ( .not. diagno%diagnostic_active ) then
            diagno%diagnostic_active = .true.
            diagno%time_start_sum_diagnostics = time_previous_step
         endif
         !diagno%time_end_sum_diagnostics = time
         diagno%time_end_sum_diagnostics = max( diagno%time_end_sum_diagnostics , time )
      endif

      if (output__Orbit) then
         if ( ( .not. interaction_failed_particle_overshot_resonance ) .and. &
               (time.gt.start_time_event_output) ) then
            !if (abs(rand_uniform_var0mean1()) .lt. 4d-2) then  ! Only to reduce the amount of out put for long simulations
            call store2file_RFOF_orbit_info(marker,time,time-time_previous_step,MPI_nod_Id)
            !endif
         endif
      endif

   end subroutine RFOF_master


   !--------------------------------------------------------------------------------
   ! subroutine split_RFstep_in_single_mode()
   !

   subroutine split_RFstep_in_single_mode(marker,resonanceMemory,RFglobal, &
         nr_resonant_waves,list_resonant_waves,diagno,MPI_nod_Id)

      use RFOF_local_magnetic_field, only: get_local_magnetic_field

      ! Input
      type(resonance_memory), intent(in) :: resonanceMemory(:,:)
      integer, intent(in) :: nr_resonant_waves
      real(8), pointer, intent(in) :: list_resonant_waves(:,:)
      integer, intent(in) :: MPI_nod_Id

      ! Input/Output
      type(particle), intent(inout) :: marker
      type(rf_wave_global), intent(inout) :: RFglobal
      type(RFOF_cumlative_diagnostics), intent(inout) :: diagno

      ! Local
      type(magnetic_field_local), target :: Blocal    ! Local properties of the magnetic field
      type(rf_wave_local) :: RFlocal
      integer :: jWave, jfreq, jnphi

      ! Get the magnetic field values at the location of the test particle
      Blocal = get_local_magnetic_field(marker%R,marker%phi,marker%z)


      do jWave=1,nr_resonant_waves
         jfreq = int(list_resonant_waves(jWave,2)+0.5)
         jnphi = int(list_resonant_waves(jWave,3)+0.5)

         !write(0,'(A, A,I2, A,I2, A,E11.2, A,I3)') &
         !     " -- Resonant wave field:", &
         !     "  id1=",jfreq, &
         !     ", id2=",jnphi, &
         !     ", time=", real(list_resonant_waves(jWave,1)) , &
         !     ", nphi=", int(RFglobal%waves%global_param%ntor(jfreq,jnphi)+0.01)

         RFglobal%j_freq_present_wave = jfreq
         RFglobal%j_nphi_present_wave = jnphi
         RFlocal = get_rf_wave_local(marker%R,marker%phi,marker%z,RFglobal)

         call quasilinear_RF_kick_steinbrecher_integrator(marker,&
               resonanceMemory(jfreq,jnphi),Blocal,RFlocal,RFglobal,diagno,MPI_nod_Id)
      enddo

   end subroutine split_RFstep_in_single_mode

   !--------------------------------------------------------------------------------
   subroutine RFOF_destructor(wave, mem, diagno)

      use RFOF_Efield_update

      type(rf_wave_global), intent(inout) :: wave
      type(resonance_memory), pointer, intent(out) :: mem(:,:)
      type(RFOF_cumlative_diagnostics), intent(inout) :: diagno

      call destructor_RFOF_Efield_update
      call rf_wave_destructor(wave)
      call destructor_rf_resonance_memory_matrix(mem)
      call destructor_RFOF_cumlative_diagnostics(diagno)

   end subroutine RFOF_destructor


end module RFOF_main
