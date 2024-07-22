module diffusion_coef5
implicit none

contains

subroutine diffusion_coefficient5(x,d,Error_Flag)
implicit none
real(8), intent(in)::x
real(8), intent(out)::d 
integer, intent(out)::Error_flag
!Error_flag=0 in normal case, Error_flag=11 when x<0

if(x.LT.0) then
    Error_flag=11 
    print *, 'wrong x in subroutine diffusion_coefficient(x,d,Error_Flag). x=',x
     return
      endif
      d= x/(1.d0+x) ! Test the general aspects of subroutines, no singularity test
      return
      end subroutine diffusion_coefficient5

      subroutine D_diffusion_coefficient5(x,Dd)
      implicit none
      real(8), intent(in)::x
      real(8), intent(out)::Dd
      if(x.LE.0) then
          Dd=0; return
      else
          Dd= 1.d0/(1.d0+x)**2  
      endif
   return
 end subroutine D_diffusion_coefficient5 

 !---------------------------------------------------------
subroutine diffusion_coefficient5b(x,d,Error_Flag)
implicit none
real(8), intent(in)::x
real(8), intent(out)::d 
integer, intent(out)::Error_flag
!Error_flag=0 in normal case, Error_flag=11 when x<0
if(x.LE.0) then
    Error_flag=11 
    print *, 'wrong x in subroutine diffusion_coefficient(x,d,Error_Flag). x=',x
     return
      endif
      d= x*x ! Test the general aspects of subroutines, no singularity test
      return
      end subroutine diffusion_coefficient5b

 
      subroutine D_diffusion_coefficient5b(x,Dd)
      implicit none
      real(8), intent(in)::x
      real(8), intent(out)::Dd
      if(x.LE.0) then
          Dd=0; return
      else
          Dd= 2.d0*x  
      endif
   return
 end subroutine D_diffusion_coefficient5b 

 
end module diffusion_coef5
