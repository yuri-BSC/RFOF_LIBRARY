subroutine local_magnetic_field_interface_to_RFOF(R,phi,z, &
     psi,theta,Bmod,F,psi_Estatic,dBmod_dpsi,dF_dpsi,dBmod_dtheta,dF_dtheta)

  use dum_magnetic_field

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

  call local_field_from_dum_magnetic_field(R,phi,z, &
       psi,theta,Bmod,F,psi_Estatic,dBmod_dpsi,dF_dpsi,dBmod_dtheta,dF_dtheta)

end subroutine local_magnetic_field_interface_to_RFOF
