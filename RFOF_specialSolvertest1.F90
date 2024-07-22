program specialSolvertest
  use RFOF_SpecialSolver 
  implicit none
  intrinsic CPU_TIME 
  real(8):: ntraject !> Number of MC trajectory sampling
  ntraject=1.0d2
  call multi_step_boundarytest1(ntraject)
end program specialSolvertest

!------------------------------------------------
function exact_moments(x) result(y)
  real(8), intent(in)::x
  real(8)::y
  y=(2.d0+x)**(-2)
  return
end function exact_moments
!---------------------------------------------------------------------------

