                program SDE_Steinb4
                use RFOF_random_numbers_SG2
                !use RFOF_random_generator_init
                implicit none
                intrinsic CPU_TIME 
                integer , parameter:: dim=2, nr_obs=3

                        integer::   ndiv, nout 

                real:: Ntraject, tfinal, xstart(0:dim), observables(nr_obs) 

                integer:: idim ! internal
                integer::IJKLIN, NTOTIN, NTOT2N !< RM48 initialization
                REAL:: time_begin, time_end     !<For speed test
                character(len=13)::outp_fname   !<Output filename 
                common /b1/outp_fname
                outp_fname='sdetest03.txt'
 CALL CPU_TIME ( time_begin )
 PRINT*, 'Time  start= ', time_begin  , ' seconds' 
 call init1_RM48(0); 
  ndiv=20;
 nout=5 ; ! Number of output time values
Ntraject=10000.0;  tfinal=1
                 do idim=0, dim
                 xstart(idim)=idim*0.001 
                 enddo
         
 call init_observables(dim, nr_obs) 
 !call gen_observables(dim, nr_obs, xstart, observables)         
 call IntegratorTest(dim, nr_obs, Ntraject, xstart, tfinal, nout, ndiv)   
CALL CPU_TIME ( time_end )
  PRINT *, 'CPU time= ', time_end - time_begin, ' seconds'
stop ' END program  '
end program SDE_Steinb4

         

                                  
                                subroutine SDEspec_f(dim, x,  f  )

                                implicit none
 integer, intent(in):: dim
 real, dimension(0:dim), intent(in):: x
 real, dimension(dim), intent(out):: f
 real :: s
 integer::idim
 s=x(0)-4;
 do idim =1, dim
 s=s+x(idim)
 end do
 do idim=1, dim
 f(idim)=-x(idim)*sin(-6+s);
 end do
end subroutine SDEspec_f


subroutine SDEspec_a(dim, x, a   )
implicit none
! dx0=a(x)dt+b(x)dw; dx_i=f_i(x)dx0 ;1.le.i .le. dim
  integer, intent(in):: dim
  real, dimension(0:dim), intent(in):: x
  integer idim                     
  real, intent(out):: a 
  real::s
  s=3;
  do idim =0, dim
  s=s+x(idim)
  end do 
  a= cos(s)*sin(s) 
!! $a(x)= partial d b(x) / partial dx_0$, as in FP eq
end subroutine SDEspec_a

subroutine SDEspec_b(dim, x, b )
        implicit none
        ! dx_0=a(x)dt+b(x)dw; dx_i=f_i(x)dx_0 ;1.le.i .le. dim
         integer, intent(in):: dim
         real, dimension(0:dim), intent(in):: x
                                   
          real, intent(out):: b 
          integer idim
          real::s
          s=3;
        do idim =0, dim
        s=s+x(idim)
        end do 
        b= sin(s) 
        end subroutine SDEspec_b 

            Subroutine EulerMaruyama_1step(dim, dt, xold, xnew)
                                use RFOF_random_numbers_SG2
                                implicit none
                                
                                integer, intent(in):: dim
                                real, intent(in):: dt
                        real, dimension(0:dim), intent(in):: xold
                                real, dimension(0:dim), intent(out):: xnew
                                real :: dw, dx0, aa, bb 
                                real, dimension(dim):: ff
                                real(8)::grw(2)
                                integer:: idim1 
                                 
                                call  SDEspec_a(dim, xold, aa   )
                                dx0=aa*dt
                                call    gaussian_random_vector(grw)
                                !here gy unused, will be used at extrapolation
                                 dw=grw(1)*sqrt(dt)
                                call  SDEspec_b(dim, xold, bb   )
                                dx0=dx0+bb*dw 
                                xnew(0)=xold(0)+dx0
                                call SDEspec_f(dim, xold,  ff  )
                                 do idim1=1, dim
                                xnew(idim1)=xold(idim1)+ ff(idim1)*dx0
                                enddo 
                                end Subroutine EulerMaruyama_1step 

                                Subroutine Steinbrecher_1step(dim, dt, xold, xnew)
                                use RFOF_random_numbers_SG2
                                implicit none
                                integer, intent(in):: dim
                                real, intent(in):: dt
                        real, dimension(0:dim), intent(in):: xold
                                real, dimension(0:dim), intent(out):: xnew 
                                 real ,dimension(0:dim):: xinterm
                                real :: dw, dx0, aa, bb1, bb2 
                                real, dimension(dim):: ff1
                                real :: dx0interm, x0interm
                                real(8)::grw(2)
                                integer:: idim1                          
                                call    gaussian_random_vector(grw)
                                !here gy unused, will be used at extrapolation
                                 dw=grw(1)*sqrt(dt)
                                call  SDEspec_b(dim, xold, bb1   )
                                xinterm(0)=xold(0)+bb1*dw 
                                call SDEspec_f(dim, xold,  ff1  )
                                 do idim1=1, dim
                             xinterm(idim1)=xold(idim1)
                                 enddo
                                call  SDEspec_b(dim, xinterm, bb2   )
                                 dx0=bb2*dw 
                                xnew(0)=xold(0)+dx0 
                                do idim1=1, dim
                             xnew(idim1)=xold(idim1)+ ff1(idim1)*dx0
                                 enddo
                                end Subroutine Steinbrecher_1step 




                                subroutine init_observables(dim, nr_obs)
                                implicit none
                                        integer, intent(in):: dim, nr_obs
                                        real, dimension(10, 0:10):: k_obs
                                    real, dimension(10):: obs_phases
                                        common /obs/ k_obs, obs_phases
                                        integer:: idim, iobs            
                                        do iobs=1, nr_obs
                                        obs_phases(iobs)= 0 !5*iobs
                                           do idim=0, dim
                                            k_obs(iobs, idim)=  iobs*2.0+ idim
                                           end do
                                        end do
                                        print *, ' end init observables'
                                        end subroutine init_observables

                                subroutine  gen_observables(dim, nr_obs, x, observables)
                                implicit none
                                integer, intent(in):: dim, nr_obs
                               real, dimension(0:dim), intent(in):: x 
                                real, dimension(nr_obs), intent(out):: observables
                                real, dimension(10, 0:10):: k_obs
                                real, dimension(10):: obs_phases
                                real::obs,ab
                                integer:: idim, iobs
                                common /obs/ k_obs, obs_phases
                                 do  iobs=1, nr_obs
                                   obs=obs_phases(iobs)
                                   do idim=0, dim
                                   obs =obs +x(idim)*k_obs(iobs,idim)
                                   end do
                                    obs=sin(obs)
                                   obs=obs*obs
                                 observables(iobs)=obs 
                        !        print * ,' iobs=' ,iobs, ' observ=',  observables(iobs)
                                 end do 
                                end subroutine  gen_observables


                                subroutine IntegratorTest(dim, nr_obs, Ntraject, xstart, tfinal, nout, ndiv)
                                implicit none
                                integer, intent(in):: dim, nout, ndiv, nr_obs
                                real, dimension(0:dim), intent(in) :: xstart
                                real,   intent(in)::Ntraject,   tfinal 
                                real, dimension(0:dim):: x, xnew
                                real, dimension(nr_obs, nout, 4 ) ::mean_obs
                                 real, dimension(nr_obs) :: observables
                                 real::count_traject,dt, ox(4) 
                                integer:: i4, ix, idim, iobs, iout, idiv
                                  character(len=13)::outp_f ! output file
                                common /b1/outp_f
  
                                 mean_obs=0
                                 dt=(tfinal/ndiv)/nout
                                 do 1000 i4=1, 4
                                 count_traject=0
                                 do 30 while(count_traject .LE. Ntraject) 
                                 count_traject=count_traject+1   
                                 do idim=0, dim
                                 x(idim)  =xstart(idim)
                                 enddo
                                 ! START INTEGRATION
                                
                                    do 10 iout=1, nout
                                   do  20 idiv=1, ndiv 
                                        if (i4.LT.3) then
                                          call EulerMaruyama_1step(dim, dt, x, xnew)
                                        else 
                                          call Steinbrecher_1step(dim, dt, x, xnew)
                                        endif
                                        x=xnew 
        20                              enddo
                                 call   gen_observables(dim, nr_obs, x, observables) 
                                 do iobs=1, nr_obs
                                 mean_obs(iobs, iout, i4 )=mean_obs(iobs, iout, i4 )+observables(iobs)
                                 enddo
        10                              enddo  ! end loop nout
        30                        end do ! end loop ntraject
                                
        1000                   end do ! different integrators
                                mean_obs=mean_obs/Ntraject 
                              ! start output integrator test
                               print *, 'Start output integrator'  
                               print *, ' iout,   x1,  x2,  x3,  x4, (x1+x2)/2, (x3+x4)/2 '

                    open(UNIT=3,file=outp_f,form='formatted',status='new', action='readwrite')
                   print *, ' iout, x1, x2,x3, x4, (x1+x2)/2, (x3+x4)/2 '
                   !write(UNIT=3 ) ' iout, x1, x2,x3, x4, (x1+x2)/2, (x3+x4)/2 '
                                
                                  do iobs=1, nr_obs
                                  print *,''
                                  !write(UNIT=3) ' '
                                   print *,' observable Nr=', iobs
                                   write(UNIT=3,fmt=101) 'observable=', iobs
                                 do iout=1, nout
                                 do ix=1,4
                                 ox(ix)=mean_obs(iobs, iout,ix)
                                 enddo 
                         write (UNIT=3, fmt=100) iout,  ox ,  (ox(1)+ox(2))/2 ,   (ox(3)+ox(4))/2
                          write (*      , fmt=100) iout, ox ,  (ox(1)+ox(2))/2 , (ox(3)+ox(4))/2
                                 enddo; enddo ;     
             
                         close(UNIT=3,  status='keep' )
         100             format( I2, 6(1X, E12.4) )
         101             format( A20,I2)       
                         print *, 'end of   tests'
                        end subroutine IntegratorTest

        
