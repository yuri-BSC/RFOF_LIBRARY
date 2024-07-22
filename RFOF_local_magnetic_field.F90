#include "config.h"

module RFOF_local_magnetic_field

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  implicit none

contains

  !--------------------------------------------------------------------------------
  ! subroutine get_local_magnetic_field(x1,x2,x3,Bfield,Blocal)
  !
  ! Get the local magetic fields
  !
  !--------------------------------------------------------------------------------
  function get_local_magnetic_field(R,phi,z,marker) result(Blocal)

    use RFOF_magnetic_field, only: magnetic_field_local
    use RFOF_markers, only: particle

    ! Input
    real(8), intent(in) :: R
    real(8), intent(in) :: phi
    real(8), intent(in) :: z
    type(particle), optional, intent(in) :: marker

    ! Output
    type(magnetic_field_local) :: Blocal

#ifdef USE_MAGNETIC_FIELD_INTERFACE
    interface
       subroutine local_magnetic_field_interface_to_RFOF(R,phi,z, &
            psi,theta,Bmod,F,psi_Estatic,dBmod_dpsi,dF_dpsi,dBmod_dtheta,dF_dtheta)
         real(8), intent(in) :: R
         real(8), intent(in) :: phi
         real(8), intent(in) :: z
         real(8), intent(out) :: psi
         real(8), intent(out) :: theta
         real(8), intent(out) :: Bmod
         real(8), intent(out) :: F
         real(8), intent(out) :: psi_Estatic
         real(8), intent(out) :: dBmod_dpsi 
         real(8), intent(out) :: dF_dpsi
         real(8), intent(out) :: dBmod_dtheta
         real(8), intent(out) :: dF_dtheta

       end subroutine local_magnetic_field_interface_to_RFOF
    end interface
#endif  ! USE_MAGNETIC_FIELD_INTERFACE
    Blocal%R = R
    Blocal%phi = phi
    Blocal%z = z
    print *, marker%mass
#ifdef USE_MAGNETIC_FIELD_INTERFACE
    call local_magnetic_field_interface_to_RFOF(R,phi,z, &
         Blocal%psi,Blocal%theta,Blocal%Bmod,Blocal%F,Blocal%psi_Estatic, &
         Blocal%dBmod_dpsi,Blocal%dF_dpsi,Blocal%dBmod_dtheta,Blocal%dF_dtheta)
#else  ! USE_MAGNETIC_FIELD_INTERFACE
    Blocal%Bmod = marker%mass ! (0.5d0 * marker%mass) /  1 !marker%magneticMoment  ! (0.5d0 * marker%mass * marker%vperp**2) /  marker%magneticMoment  
    print *, "test-1"
    Blocal%psi = marker%psi
    Blocal%F = (marker%Pphi - marker%charge * marker%psi) / &
         (marker%mass * marker%vpar / Blocal%Bmod)
  !! REMOVE theta!!!  Blocal%theta = 0d0
    Blocal%psi_Estatic = marker%energy - marker%energy_kinetic

#endif  ! USE_MAGNETIC_FIELD_INTERFACE

  end function get_local_magnetic_field

end module RFOF_local_magnetic_field
