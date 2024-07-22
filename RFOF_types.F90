#include "config.h"

module RFOF_types

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  use RFOF_constants
  use RFOF_parameters
  use RFOF_waves
  use RFOF_magnetic_field
  use RFOF_markers
  use RFOF_resonance_memory
  use RFOF_resonance_condition
  use RFOF_random_numbers
  use RFOF_diagnostics

  implicit none

  type :: RFOF_state
         type(particle), pointer             :: marker   !< Marker
         type(resonance_memory), pointer     :: mem      !< Resonance memory of the marker
         type(magnetic_field_local), pointer :: Blocal   !< Local magnetic field
         type(rf_wave_local), pointer        :: RFlocal  !< Local RF wave field
         type(rf_wave_global), pointer       :: RFglobal !< Global RF wave field
         real(8)                             :: Iperp    !< The magnetic moment normalised to \f$ Ze/m\omega \f$ in units [MeV]. Useful parameter for cyclotron interaction.
  end type RFOF_state

end module RFOF_types
