module dum_magnetic_field

  implicit none

  type, public :: magnetic_field_global
     real(8) :: B0
     real(8) :: R0
     real(8) :: q
     real(8) :: aminor
  end type magnetic_field_global

  type(magnetic_field_global) :: Bglobal

contains

  subroutine magnetic_field_constructor_simple_field(R0,B0,q,aminor,Bfield)

    real(8), intent(in) :: R0, B0, q, aminor
    type(magnetic_field_global), intent(out) :: Bfield

    Bfield%R0=R0
    Bfield%B0=B0
    Bfield%q =q
    Bfield%aminor=aminor

  end subroutine magnetic_field_constructor_simple_field

  subroutine local_field_from_dum_magnetic_field(R,phi,z, &
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

    real(8) :: rmin, q

    rmin = sqrt( (R-Bglobal%R0)**2 + z**2 )
    q = Bglobal%q
    psi = rmin**2 * Bglobal%B0/(2*q)
    !psia = Bglobal%aminor**2 * Bglobal%B0/(2*q)
    !Btheta = Bglobal%B0 * r / ( q * x1 )

    psi          = psi
    theta        = sign(1d0,z)*acos((R-Bglobal%R0)/rmin)

    Bmod         = Bglobal%B0 * Bglobal%R0 / R
    F            = Bglobal%B0 * Bglobal%R0
    !Btheta       = Btheta
    !BR           =-Btheta * z/rmin
    !Bz           = Btheta * (R-Bglobal%R0)/rmin
    psi_Estatic  = 0.0

    dBmod_dpsi   = 0.0
    dBmod_dtheta = 0.0
    dF_dpsi      = 0.0
    dF_dtheta    = 0.0

  end subroutine local_field_from_dum_magnetic_field

end module dum_magnetic_field
