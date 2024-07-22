Module RFOF_BesselProcess0 

  use RFOF_random_numbers_SG2

contains

  Subroutine Bessel0(x_in,delta_t,C_dif,x_fin)
    use RFOF_random_numbers_SG2
    implicit none
    real(8), intent(in)::x_in    !< initial position
    real(8), intent(in)::delta_t !< the time interval
    real(8),intent(in)::C_dif    !< From local aproximation: \f[ Diff_coef=2*C_dif*x+O(x**2) \f]
    real(8), intent(out)::x_fin  !< the final position
    real(8):: grw(2) !Random normal vector, uncorrelated components
    real(8)::x1,x2,csqrtdt ! local variables
    csqrtdt=sqrt(delta_t*C_dif)
    call gaussian_random_vectorRM48(grw)
    x1=grw(1)*csqrtdt+sqrt(x_in)
    x2=grw(2)*csqrtdt
    x_fin=x1*x1+x2*x2
  end Subroutine Bessel0

  !------------------------------------------------------------------------
  Subroutine Bessel0_General(x_in,delta_t,C_dif,D_0,x_fin)
    use RFOF_random_numbers_SG2
    implicit none
    real(8), intent(in)::x_in !initial position
    real(8), intent(in)::delta_t !the time interval
    real(8),intent(in)::C_dif, D_0 !From local aproximation: Diff_coef=2*C_dif*x +D_0+error
    real(8), intent(out)::x_fin ! the final position
    real(8):: grw(2) !Random normal vector, uncorrelated components
    real(8)::normal_random ! random normal numbers
    real(8)::x1,x2,csqrtdt, sqrtxin,d2c ! local variables
    if(C_dif.NE.0) Then 
       d2c=D_0/(2.0d0*C_dif)
       csqrtdt=sqrt(delta_t*dabs(C_dif) ) 
       if(C_dif>0) Then
          sqrtxin=sqrt(d2c+x_in)
          call gaussian_random_vectorRM48(grw)
          x1=grw(1)*csqrtdt+sqrtxin
          x2=grw(2)*csqrtdt
          x_fin=x1*x1+x2*x2-d2c
       else    ! C_dif<0        
          sqrtxin=sqrt(-d2c-x_in)
          call gaussian_random_vectorRM48(grw)
          x1=grw(1)*csqrtdt+sqrtxin
          x2=grw(2)*csqrtdt
          x_fin=-x1*x1-x2*x2-d2c
       endif
    else
       ! C_dif=0
       call Marsaglia_Bray_RM48(normal_random)
       x_fin=x_in+sqrt(2.0d0*delta_t*D_0)*normal_random
    endif
  end Subroutine Bessel0_General


  !------------------------------------------------------------------------
  subroutine exact_moments(x_in, delta_t,C_dif, moments)
    implicit none
    real(8), intent(in)::x_in !starting point of  0 order bessel process
    real(8), intent(in)::delta_t ! duration of Bessel procees

    real(8), intent(in)::C_dif!From local aproximation: Diff_coef=2*C_dif*x 
    integer, parameter::max_moment_order=2!Returns all moments from 1 to this value
    real(8), intent(out)::moments(max_moment_order)
    ! The first 2 exact moments of the 0 order Bessel Procees, starting at , used for test 
    real(8) :: cdt ! local
    cdt=C_dif*delta_t
    moments(1)=x_in+2.d0*cdt
    moments(2)=8.0d0*cdt*cdt+8.0d0*cdt*x_in+x_in*x_in
  end subroutine exact_moments

  !------------------------------------------------------------------------
  subroutine exact_momentsG(x_in, delta_t,C_dif,D_0 , moments)
    implicit none
    real(8), intent(in)::x_in !starting point of  0 order Bessel process
    real(8), intent(in)::delta_t ! duration of Bessel procees
    real(8), intent(in)::C_dif,D_0 !From local aproximation: Diff_coef=2*C_dif*x +D_0+O(x**2)
    integer, parameter::max_moment_order=2!Returns all moments from 1 to this value
    real(8), intent(out)::moments(max_moment_order)
    ! The first 2 exact moments of the 0 order Bessel Procees, starting at , used for test 
    real(8) :: cdt ! local
    cdt=C_dif*delta_t
    moments(1)=x_in+2.d0*cdt
    moments(2)=8.0d0*cdt*cdt+8.0d0*cdt*x_in+x_in*x_in+2.d0*D_0*delta_t
  end subroutine exact_momentsG

End Module RFOF_BesselProcess0
