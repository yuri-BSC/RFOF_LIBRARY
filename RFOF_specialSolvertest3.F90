program specialSolvertest3 
  use RFOF_SpecialSolver
  implicit none
  intrinsic CPU_TIME 
  real(8):: ntraject !< Number of MC trajectory sampling
  ntraject=1.0d2zd
  call multi_step_generaltest3(ntraject)
end program specialSolvertest3
