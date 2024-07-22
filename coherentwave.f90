module coherentwave_mod


   use fullwave
   use global_param

   type :: chw
      type(fw), allocatable :: fullwave 
      type(glbl), allocatable :: global_param
   end type chw

end module coherentwave_mod


