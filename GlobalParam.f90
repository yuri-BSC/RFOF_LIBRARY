module global_param
   implicit none
   public :: glbl

   type :: glbl
      real, allocatable :: ntor(:)
      real, allocatable :: p_frac_ntor(:)
      real, allocatable :: pow_i(:)
      real, allocatable :: pow_ntor_i(:,:)
      real, allocatable :: pow_ntor_e(:)
      real, allocatable :: cur_tor_ntor(:)
      real :: power_tot
      real :: pow_e 
      real :: cur_tor 
      real :: frequency    
      character :: name
      character :: type  
!      character(:), allocatable :: name(:)
!      character(:), allocatable :: type(:)
   end type glbl

end module global_param

 