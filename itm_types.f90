module itm_types  
   implicit none

   ! Define the kind parameter
   integer, parameter :: r8 = SELECTED_REAL_KIND(15, 300)

   ! Declare and initialize the variable with the specified kind
   real(r8), parameter :: itm_r8_invalid = -9.0_r8 ** 40


end module itm_types

