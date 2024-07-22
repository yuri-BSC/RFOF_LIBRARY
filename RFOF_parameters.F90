#include "config.h"

module RFOF_parameters

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  implicit none

  !--------------------------------------------------------------------------------
  type :: rz_boundingbox
     real(8) :: Rmin = 0.1d0
     real(8) :: Rmax = 20.0d0
     real(8) :: Zmin = -20.0d0
     real(8) :: Zmax = 20.0d0
  end type rz_boundingbox

  !--------------------------------------------------------------------------------

  character(128), PARAMETER :: RFOF_version = "version 0"

  real(8), PARAMETER :: very_large_time_in_seconds = 1d10

  !--------------------------------------------------------------------------------
  ! I/O channels
  !--------------------------------------------------------------------------------

  integer, parameter :: io_channel_3871 = 3871  ! Kept closed
  integer, parameter :: io_channel_3872 = 3872  ! Kept closed. Used by store2file_RFOF_cumlative_diagnostics, files: out.RFOF_integrate_1D, out.RFOF_integrate_2DRZ
  integer, parameter :: io_channel_3873 = 3873  ! Kept open by save2file_rf_kicks
  integer, parameter :: io_channel_3874 = 3874  ! Kept open by save2file_resonace_predictions
  integer, parameter :: io_channel_3875 = 3875  ! Kept open by store2file_RFOF_orbit_info
  integer, parameter :: io_channel_3876 = 3876  ! Kept open by 
  integer, parameter :: io_channel_3877 = 3877  ! Kept open by 
  integer, parameter :: io_channel_3878 = 3878  ! Kept open by 
  integer, parameter :: io_channel_3879 = 3879  ! Kept open by store2file_RFOF_efield_normalization

  !--------------------------------------------------------------------------------
  ! Quasilinear model parameters
  !--------------------------------------------------------------------------------

  real(8), save :: width_of_rf_resonance_layer = 1.0d-2


  !--------------------------------------------------------------------------------
  ! Simplifying/ assumptions - default setting
  !--------------------------------------------------------------------------------

  logical, save :: simplify__static_resonance_position_during_RF_kick = .TRUE.

  logical, save :: simplify__drift_velocity_does_not_affect_resonance_condition = .TRUE.

  logical, save :: simplify__parallel_velocity_does_not_affect_resonance_condition = .TRUE.

  logical, save :: simplify__assume_zero_larmor_radius_in_KPERPxRHO = .FALSE.

  logical, save :: simplify__kpar_is_nphi_over_R = .TRUE.


  !--------------------------------------------------------------------------------
  ! I/O control
  !--------------------------------------------------------------------------------

  real(8), save :: start_time_event_output = 0e-4 !- 1e-4

  logical, save :: output__2D_RZ_out = .FALSE.

  logical, save :: output__Orbit = .FALSE.
  integer, save :: number_of_points_stored_in_the_Orbit = 0  ! Counter; updated each time a new point is added
  integer, save :: MAX_number_of_points_stored_in_the_Orbit = 10000

  logical, save :: output__rf_kicks = .TRUE.
  integer, save :: number_of_points_stored_in_the_rf_kick = 0  ! Counter; updated each time a new point is added
  integer, save :: MAX_number_of_points_stored_in_the_rf_kick = 10000

  logical, save :: output__resonace_predictions = .FALSE.
  integer, save :: number_of_points_stored_in_the_resonance_memory = 0  ! Counter; updated each time a new point is added
  integer, save :: MAX_number_of_points_stored_in_the_resonance_memory = 10000

  logical, save :: output__efield_normalization = .TRUE.
  integer, save :: number_of_points_stored_in_the_efield_normalization = 0  ! Counter; updated each time a new point is added
  integer, save :: MAX_number_of_points_stored_in_the_efield_normalization = 100000

  type(rz_boundingbox), save :: plasma_boundingbox

  integer, save :: NRedges_2DgridRZ = 60
  integer, save :: NZedges_2DgridRZ = 100

  logical, save :: ignore_non_serious_error_messages = .TRUE.


contains

  !--------------------------------------------------------------------------------
  subroutine initialise_RFOF_parameters(Rmin_bound,Rmax_bound,Zmin_bound,Zmax_bound)

    ! Input
    real(8), optional, intent(in) :: Rmin_bound,Rmax_bound,Zmin_bound,Zmax_bound

    real(8) :: Rmin, Rmax, Zmin, Zmax

    NAMELIST/simplify_rfof/ &
         simplify__static_resonance_position_during_RF_kick , &
         simplify__drift_velocity_does_not_affect_resonance_condition , &
         simplify__parallel_velocity_does_not_affect_resonance_condition , &
         simplify__assume_zero_larmor_radius_in_KPERPxRHO , &
         simplify__kpar_is_nphi_over_R, &
         width_of_rf_resonance_layer
    NAMELIST/rz_boundingbox/Rmin, Rmax, Zmin, Zmax

    open( io_channel_3872, FILE='input.rfof')
    read( io_channel_3872, simplify_rfof)
    read( io_channel_3872, rz_boundingbox)
    close(io_channel_3872)

    if (present(Rmin_bound) .and. present(Rmax_bound)) then
       plasma_boundingbox%Rmin = Rmin_bound
       plasma_boundingbox%Rmax = Rmax_bound
    endif
    if (present(Zmin_bound) .and. present(Zmax_bound)) then
       plasma_boundingbox%Zmin = Zmin_bound
       plasma_boundingbox%Zmax = Zmax_bound
    endif

  end subroutine initialise_RFOF_parameters

  !--------------------------------------------------------------------------------
  function is_inside_plasma(R,Z) result(is_inside)
    
    ! Input
    real(8), intent(in) :: R, Z

    ! Output
    logical :: is_inside

    is_inside = ( &
         (R.gt.plasma_boundingbox%Rmin) .and. &
         (R.lt.plasma_boundingbox%Rmax) .and. &
         (Z.gt.plasma_boundingbox%Zmin) .and. &
         (Z.lt.plasma_boundingbox%Zmax) )

  end function is_inside_plasma

end module RFOF_parameters
