module diffusion_coef
implicit none
contains

subroutine diffusion_coefficient(x,d,Error_Flag)
real(8), intent(in)::x
real(8), intent(out)::d 
integer, intent(out)::Error_flag



if(x.LE.0) then
Error_flag=11 
print *, 'wrong x in subroutine diffusion_coefficient(x,d,Error_Flag). x=',x
 return
 endif
!x2=x*x
d= x/(1+x**2)
return
end subroutine diffusion_coefficient

subroutine D_diffusion_coefficient(x,Dd)
real(8), intent(in)::x
real(8), intent(out)::Dd

if(x.LE.0) then
Dd=0; return
else
Dd= 2.d0*x/(1+x**2)**2  
endif
return
end subroutine D_diffusion_coefficient

end module diffusion_coef
