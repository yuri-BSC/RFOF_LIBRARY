#include "config.h"

module RFOF_magnetic_field

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif
!   use RFOF_types
    use itm_types    ! canvi Yuri

  implicit none

  type, public :: magnetic_field_local
     real(8) :: R = itm_r8_invalid

     real(8) :: phi = itm_r8_invalid

     real(8) :: z = itm_r8_invalid

     real(8) :: psi = itm_r8_invalid

     real(8) :: theta = itm_r8_invalid

     real(8) :: Bmod = itm_r8_invalid

     real(8) :: F = itm_r8_invalid

     real(8) :: psi_Estatic = itm_r8_invalid

     real(8) :: dBmod_dpsi = itm_r8_invalid

     real(8) :: dF_dpsi = itm_r8_invalid

     real(8) :: dBmod_dtheta = itm_r8_invalid

     real(8) :: dF_dtheta = itm_r8_invalid
  end type magnetic_field_local

end module RFOF_magnetic_field
