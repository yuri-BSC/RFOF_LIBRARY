Module RFOF_BesselProcess0 
  use RFOF_random_numbers_SG2

contains

  !------------------------------------------------------------------------
  Subroutine Bessel0(x_in,delta_t,C_dif,x_fin)
    use RFOF_random_numbers_SG2
    implicit none
    real(8), intent(in) ::x_in    !< Initial position
    real(8), intent(in) ::delta_t !< The time interval
    real(8), intent(in) ::C_dif   !< From local aproximation: \f[ Diff_coef=2*C_dif*x+O(x**2) \f]
    real(8), intent(out)::x_fin   !< The final position
    real(8):: grw(2)              !  Random normal vector, uncorrelated components
    real(8)::x1,x2,csqrtdt        !  local variables
    csqrtdt=sqrt(delta_t*C_dif)
    call gaussian_random_vectorRM48(grw)
    x1=grw(1)*csqrtdt+sqrt(x_in)
    x2=grw(2)*csqrtdt
    x_fin=x1*x1+x2*x2
  end Subroutine Bessel0

  !------------------------------------------------------------------------
  subroutine exact_moments(x_in, delta_t,C_dif, moments)
    implicit none
    real(8), intent(in)::x_in !starting point of  0 order bessel process
    real(8), intent(in)::delta_t ! duration of Bessel procees
    real(8), intent(in)::C_dif ! Linearity coefficient, see above
    integer, parameter::max_moment_order=2!Returns all moments from 1 to this value
    real(8), intent(out)::moments(max_moment_order)
    ! The first 2 exact moments of the 0 order Bessel Procees, starting at , used for test 
    real(8) :: cdt ! local
    cdt=C_dif*delta_t
    moments(1)=x_in+2.d0*cdt
    moments(2)=8.0d0*cdt*cdt+8.0d0*cdt*x_in+x_in*x_in
  end subroutine exact_moments

End Module RFOF_BesselProcess0


