subroutine MC_kick_steinbrecher_integrator( &
     dt, &
     x, &
     calc_diffusion_coeff, &
     state, &
     relative_tolerance_MonteCarlo, &
     error_flag, &
     MPI_node_Id)

  use RFOF_types
  implicit none

  real(8), intent(in) :: dt                      !< Time step [s]
  real(8), intent(inout) :: x                    !< Marker position
  type(RFOF_state), intent(inout) :: state       !< Marker state
  real(8) :: relative_tolerance_MonteCarlo       !< Relative tolerance in Monte Carlo stepping
  integer, intent(inout) :: error_flag            !< Error flag
  integer, intent(in) :: MPI_node_Id             !< Number to identify the MPI node

  ! Local variables
  real(8) :: dx, dx0, dx1, dW, x0, x1, dW1, dW2,diffusion, diffusion1,x_a,x_b, x_b1, x_c, x_c1
  type(RFOF_state) :: state0, state1, state_a, state_b,state_b1, state_c, state_c1 
  real(8):: factor
  logical::stop_loop

  !------------------------------------------------------------------------------------
  !
  ! BEGIN: INTERFACES
  !
  !------------------------------------------------------------------------------------

  interface

     subroutine calc_diffusion_coeff(x,state,diffusion,error_flag)
       use RFOF_types
       real(8), intent(in) :: x                  !< Integration variable
       type(RFOF_state), intent(inout) :: state  !< State of RFOF
       integer, intent(out) :: error_flag        !< Flag to be raised if the routine failed. error_flag=0 marker successfully moved; 
                                                 !< error_flag>0 failed to move marker; error_flag=1 new state does not exists; 
                                                 !< error_flag=2 step too long, could not perform move
       real(8), intent(out) :: diffusion         !< diffusion coefficient
     end subroutine calc_diffusion_coeff
  end interface 
  !
  !------------------------------------------------------------------------------------
  !
  !< NEW INTERFACE, TO WRITE THE SUBROUTINE
  !< Subroutine that calculates the diffusion coefficient in the new 'state1', obtained from the 'state'
  !< by the modification of the paramater x (i.e. magnetic moment or equivalent I_perp) only ;
  !< Energy,P_fi remains the same
  interface
     subroutine modifiy_marker(x1,state, state1, error_flag) 
  
       use RFOF_types
       real(8), intent(in) :: x1                   !< The new marker position: magnetic momentum or I_perp  
       type(RFOF_state), intent(in) :: state       !< Marker initial state
       type(RFOF_state), intent(out) :: state1     !< Marker intermediate state; 
                                                   !< only parameter x is modified to x1; rest of parameters are unchanged
       integer, intent(out) :: error_flag          !< Flag to be raised if the routine failed. error_flag=0 marker successfully moved; 
                                                   !< error_flag>0 failed to move marker; error_flag=1 new state does not exists; 
                                                   !<error_flag=2 step too long, could not perform move
       
     end subroutine modifiy_marker
  end interface
  !
  !------------------------------------------------------------------------------------
  ! NEW INTERFACE
  interface
     function timestep_factor(x_a, x_c, state_a) result(timestepFactor)
       !/ compares  the one step relative error with prescribed values and  return a factor<1 , 
       !/if the eror is too large, >1 if the error is close to machine precision, and exactly 1 if OK
       use RFOF_types
       real(8)  ::timestepFactor 
       real(8), intent(in):: x_a, x_c
       type(RFOF_state), intent(in) :: state_a !, state_c ! Gyorgy why is there a state_c? There no state_c in the function declaration.
     end function timestep_factor
  end interface
  !
  !------------------------------------------------------------------------------------
  ! NEW INTERFACE
  interface 
     function close_to_boundary(x_c, state_c) result(closeToBoundary)
       !< return .TRUE. of the distance to the boundary is less then prescribed value
       !< Neumann boundary conditions 
       use RFOF_types
       logical ::closeToBoundary
       real(8), intent(in)::   x_c
       type(RFOF_state), intent(in) ::   state_c
     end function close_to_boundary
  end interface
  !
  !------------------------------------------------------------------------------------
  !/NEW INTERFACE
  interface 
     subroutine Gauss_generator(dw1, dw2)!/ pair of independent N(0,1) variable TO DO!!!!!!!!!
       real(8), intent(out):: dw1, dw2
     end subroutine Gauss_generator
  end interface

  !------------------------------------------------------------------------------------
  !
  ! END: INTERFACES
  !
  !------------------------------------------------------------------------------------



  !------------------------------------------------------------------------------------
  ! Generate a safety copy of the state...e.g. to us to step back to the initial state
  ! somethnings wrong!!        
state0 = copy_RFOF_state(state)
 x0=x 
stop_loop =.false.
  DO while (stop_loop== .false.)
call Gauss_generator(dW1, dW2) 
dW=(dW1+dW2)/sqrt(2.0d0)
 ! Calculate the diffusion coefficient
 diffusion = diffusion_coeff(state)
x1=x0+sqrt(2d0*diffusion*dt)*dW

call modifiy_marker(x1,state, state1,error_flag)
if(error_flag>0) goto 100
diffusion1 = diffusion_coeff(state1)
 ! Calculate the length of the kick
dx=sqrt(2d0*diffusion1*dt)*dW
x_a=x0+dx
 call move_marker(x_a,state_a,error_flag) ! with full step dt
if(error_flag>0) goto 100
! Now we perform 2 steps with dt/2
! First step
  x_b1=x0+sqrt( diffusion*dt)*dW1 
call modifiy_marker(x_b1,state, state1,error_flag)
if(error_flag>0) goto 100
  diffusion1 = diffusion_coeff(state1) 
dx=sqrt( diffusion1*dt)*dW1
x_b=x0+dx
 call move_marker(x_b,state_b,error_flag) 
if(error_flag>0) goto 100
! second step with DW2
 x_c1=x_b+sqrt( diffusion*dt)*dW2 
call modifiy_marker(x_c1,state_b, state1,error_flag)
  diffusion1 = diffusion_coeff(state1) 
dx=sqrt( diffusion1*dt)*dW2
x_c=x_b+dx
 call move_marker(x_c,state_c,error_flag)
if(error_flag>0) goto 100
  ! encoded in timestep_factor error=relative_error(x_a, state_a,x_c, state_c)
  factor=timestep_factor(x_a, state_a,x_c, state_c) 
100  factor=0.5d0
dt=factor*dt 
  if(factor.ge.1) then 
  stop_loop=.true.
endif
   enddo 
  if (close_to_boundary(x_c, state_c)==.true. ) then
  ! somethnings wrong, or Neumann boundary condition, remain at initial x!!     
   state = copy_RFOF_state(state0)
    x=x0;  
! eventually continue subdivision
  else 
     x=x_c
state = copy_RFOF_state(state_c) ! Better: extrapolation in dt
  endif
      
end subroutine MC_kick_steinbrecher_integrator
! TO DO :  
!              Gauss_generator(dW1, dW2) 
! TO VERIFY        subroutine move_marker(dx,state) 
!             modifiy_marker(x_b1,state, state1)
! effecst of boundary of singular domains to treat analytically
! Possible side effects in calling programs, due to change of dt to intent(inout)
