#include "config.h"

module RFOF_local_magnetic_field
use RFOF_markers 

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
  function get_local_magnetic_field(R,phi,z, mass, marker) result(Blocal)

    use RFOF_magnetic_field , only: magnetic_field_local
!    use RFOF_markers , only: particle

    ! Input
    real(8), intent(inout) :: R
    real(8), intent(inout) :: phi
    real(8), intent(inout) :: z
    real(8) :: counter=0
    real(8), target, intent(in) :: mass
    type(particle), optional, intent (inout) :: marker    

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

#ifdef USE_MAGNETIC_FIELD_INTERFACE
    call local_magnetic_field_interface_to_RFOF(R,phi,z, &
         Blocal%psi,Blocal%theta,Blocal%Bmod,Blocal%F,Blocal%psi_Estatic, &
         Blocal%dBmod_dpsi,Blocal%dF_dpsi,Blocal%dBmod_dtheta,Blocal%dF_dtheta)
#else ! USE_MAGNETIC_FIELD_INTERFACE
    
    print *, "The value of marker mass before Blocal%Bmod"
    if (counter > 0) then
      print *, "Now, step_dumorb is using local_magnetic, and the markers issues will arise, it cannot even print marker%mass"
    else 
      print *, "This time init_dumorb is using local_magnetic, no issues occur"   
    end if
    
    print *, marker%mass       
    Blocal%Bmod = (0.5d0 * marker%mass * marker%vperp**2) /  marker%magneticMoment 
    print *, "The value of  marker%psi before Blocal%psi"
    print *, marker%psi
    Blocal%psi = marker%psi  
  
    Blocal%F = (marker%Pphi - marker%charge * marker%psi) / &
        (marker%mass * marker%vpar / Blocal%Bmod)  
  !! REMOVE theta!!!  Blocal%theta = 0d0
    Blocal%theta = 0d0    
    Blocal%psi_Estatic = marker%energy - marker%energy_kinetic  
    counter = counter + 1
#endif  ! USE_MAGNETIC_FIELD_INTERFACE

  end function get_local_magnetic_field

end module RFOF_local_magnetic_field

