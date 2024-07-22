
module euitm_schemas

   use coherentwave_mod

   type :: type_waves
      type(chw), allocatable :: coherentwave(:)
   end type type_waves 


end module euitm_schemas
