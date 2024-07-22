
#include "config.h"

subroutine MC_kick_steinbrecher_integrator( &
     dt, &
     x, &
     calc_diffusion_coeff, &
     state, &
     relative_tolerance_MonteCarlo, &
     errorFlag, &
     MPI_node_Id)

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  use RFOF_types
  use RFOF_random_numbers
  implicit none

  real(8), intent(in) :: dt
  real(8), intent(inout) :: x
  type(RFOF_state), intent(inout) :: state
  real(8) :: relative_tolerance_MonteCarlo
  integer, intent(inout) :: errorFlag
  integer, intent(in) :: MPI_node_Id

  ! Local variables
  real(8) :: dx0, dx1, dW, x0, diffusion

  interface
     subroutine calc_diffusion_coeff(x,state,diffusion,errorFlag)
       use RFOF_types
       real(8), intent(in) :: x                  !< Integration variable
       type(RFOF_state), intent(inout) :: state  !< State of RFOF
       integer, intent(out) :: errorFlag         !< Error flag
       real(8), intent(out) :: diffusion         !< diffusion coefficient
     end subroutine calc_diffusion_coeff
  end interface

  !--------------------------------------------------------------------------------

  ! Calculate the diffusion coefficient
  call calc_diffusion_coeff(x,state,diffusion,errorFlag)

  ! Calculate the length of the kick
  ! (here I've just set the random number to one value just; this is only for illustration)
  dW = rand_uniform_var0mean1()
  dx0 = sqrt(2d0*diffusion*dt)*dW
  dx0 = max( 0.1d0*x , min( 2.0d0*x , dx0 ))
  x0 = x+dx0

  ! Update the diffusion coefficient at the end of the predictor step
  call calc_diffusion_coeff(x0,state,diffusion,errorFlag)

  ! Update position
  dx1 = sqrt(2d0*diffusion*dt)*dW
  dx1 = max( -0.9d0*x , min( 0.9d0*x , dx1 ))
  x = x+dx1

end subroutine MC_kick_steinbrecher_integrator

