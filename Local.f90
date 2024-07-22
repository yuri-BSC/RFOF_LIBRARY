module local
   implicit none
   public :: lcl
 
   type :: lcl
      real, allocatable :: e_plus(:,:,:)
      real, allocatable :: e_minus(:,:,:)
      real, allocatable :: e_para(:,:,:)
      real, allocatable :: e_plus_ph(:,:,:)
      real, allocatable :: e_minus_ph(:,:,:)
      real, allocatable :: e_para_ph(:,:,:)
   end type lcl
 
end module local
 

