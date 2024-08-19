module euitm_waves_interface

  use RFOF_parameters
  use euitm_schemas
!  use waves
  implicit none
 
contains

  subroutine construct_dummy_euitm_waves(waves)

    ! Output
    type(type_waves), intent(out) :: waves
    ! Local
    integer :: nfreq, nphi, nnphi, nion
    integer :: npsi, ntheta, max_npsi, max_ntheta
    integer :: jfreq, jnphi, jion, jpsi, jtheta

    real(8) :: RFpower, EfieldNormalisation, freq, kperp

    NAMELIST/input_wavefields/nfreq, nnphi, RFpower, &
    EfieldNormalisation, freq, nphi, kperp

    open( io_channel_3872, FILE='input.rfof')
    read( io_channel_3872, input_wavefields)
    close(io_channel_3872)
!    write(*,*) nfreq, nnphi, nphi , RFpower, EfieldNormalisation, freq, kperp
    nion = 1
    npsi = 1
    ntheta = 1
!    nnphi = 2
    max_npsi = npsi
    max_ntheta = ntheta

    allocate(waves%coherentwave(nfreq))
    
    !---- waves%global_param ----!
   ! allocate(waves%coherentwave(1)%global_param%ntor(1)) !Prova Yuri
     do jfreq = 1, nfreq
       allocate(waves%coherentwave(jfreq)%global_param)
      
       allocate(waves%coherentwave(jfreq)%global_param%ntor(nnphi))
       
       allocate(waves%coherentwave(jfreq)%global_param%p_frac_ntor(nnphi))
      
       allocate(waves%coherentwave(jfreq)%global_param%pow_i(nion))
       allocate(waves%coherentwave(jfreq)%global_param%pow_ntor_i(nnphi,nion))
       allocate(waves%coherentwave(jfreq)%global_param%pow_ntor_e(nnphi))
       allocate(waves%coherentwave(jfreq)%global_param%cur_tor_ntor(nnphi))
       
!       allocate(waves%coherentwave(jfreq)%global_param%name(1))     !Canvi yuri
!       allocate(waves%coherentwave(jfreq)%global_param%type(1))     !Canvi yuri
       
       waves%coherentwave(jfreq)%global_param%frequency = freq * (1d0+0.01d0*dble(jfreq-1))
       waves%coherentwave(jfreq)%global_param%name = 'Antenna made up in RFOF'
       waves%coherentwave(jfreq)%global_param%type = 'IC'
       waves%coherentwave(jfreq)%global_param%ntor = nphi
!       waves%coherentwave(jfreq)%global_param%f_assumption    ! No use here!
       waves%coherentwave(jfreq)%global_param%power_tot = RFpower
       waves%coherentwave(jfreq)%global_param%pow_e = 0d0
       waves%coherentwave(jfreq)%global_param%cur_tor = 0d0
!       waves%coherentwave(jfreq)%global_param%code_type       ! No use here!
!       waves%coherentwave(jfreq)%global_param%toroid_field... ! No use here!
       do jion = 1, nion
          waves%coherentwave(jfreq)%global_param%pow_i(jion) =  RFpower / real(nion)
       enddo

       do jnphi = 1, nnphi
          waves%coherentwave(jfreq)%global_param%ntor(         jnphi) = nphi * jnphi
          waves%coherentwave(jfreq)%global_param%pow_ntor_e(   jnphi) = 0d0
          waves%coherentwave(jfreq)%global_param%p_frac_ntor(  jnphi) = 1d0 / real(nnphi)
          waves%coherentwave(jfreq)%global_param%cur_tor_ntor( jnphi) = 0d0

          do jion = 1, nion
             waves%coherentwave(jfreq)%global_param%pow_ntor_i(jnphi,jion) = RFpower / real(nnphi * nion)
          enddo
       enddo
    enddo
    !print *, "itm_waves: global_param: nnphi", waves%global_param%nntor
    !print *, "itm_waves: global_param: nphi", waves%global_param%ntor
    !print *, "itm_waves: global_param: ptot", waves%global_param%power_tot
    !print *, "itm_waves: global_param: j", waves%global_param%cur_tor_ntor
    !print *, "itm_waves: global_param: freq", waves%global_param%frequency

    !---- waves%fullwave%local ----!
    do jfreq = 1, nfreq
       allocate(waves%coherentwave(jfreq)%fullwave)
       allocate(waves%coherentwave(jfreq)%fullwave%local)
       allocate( waves%coherentwave(jfreq)%fullwave%local%e_plus(     nnphi, max_npsi,max_ntheta) )
       allocate( waves%coherentwave(jfreq)%fullwave%local%e_minus(    nnphi, max_npsi,max_ntheta) )
       allocate( waves%coherentwave(jfreq)%fullwave%local%e_para(     nnphi, max_npsi,max_ntheta) )

       allocate( waves%coherentwave(jfreq)%fullwave%local%e_plus_ph(  nnphi, max_npsi,max_ntheta) )
       allocate( waves%coherentwave(jfreq)%fullwave%local%e_minus_ph( nnphi, max_npsi,max_ntheta) )
       allocate( waves%coherentwave(jfreq)%fullwave%local%e_para_ph(  nnphi, max_npsi,max_ntheta) )
       do jnphi = 1, nnphi
          do jpsi = 1, max_npsi
             do jtheta = 1, max_ntheta
                 waves%coherentwave(jfreq)%fullwave%local%e_plus(    jnphi, jpsi,jtheta) = 0.5d0
                 waves%coherentwave(jfreq)%fullwave%local%e_minus(   jnphi, jpsi,jtheta) = 0.5d0
                 waves%coherentwave(jfreq)%fullwave%local%e_para(    jnphi, jpsi,jtheta) = 0d0
                 waves%coherentwave(jfreq)%fullwave%local%e_plus_ph( jnphi, jpsi,jtheta) = 0d0
                 waves%coherentwave(jfreq)%fullwave%local%e_minus_ph(jnphi, jpsi,jtheta) = 0d0
                 waves%coherentwave(jfreq)%fullwave%local%e_para_ph( jnphi, jpsi,jtheta) = 0d0
              enddo
          enddo
       enddo
    enddo
    
    !print *, "e+ = ",waves%fullwave%local%e_plus
    !print *, "e- = ",waves%fullwave%local%e_minus
    !print *, "ph e+ = ",waves%fullwave%local%e_plus_ph
    !print *, "ph e- = ",waves%fullwave%local%e_minus_ph

  end subroutine construct_dummy_euitm_waves

!!$  !> Fill in a dummy waves CPO.
!!$  !> Using version phase 4.08a,  generated 21/04/2010
!!$  subroutine construct_dummy_euitm_waves(waves)
!!$
!!$    ! Output
!!$    type(type_waves), intent(out) :: waves
!!$    ! Local
!!$    integer :: nfreq, nphi, nnphi, max_nnphi, nion
!!$    integer :: npsi, ntheta, max_npsi, max_ntheta
!!$    integer :: jfreq, jnphi, jion, jpsi, jtheta
!!$
!!$    real(8) :: RFpower, EfieldNormalisation, freq, kperp
!!$
!!$    NAMELIST/input_wavefields/nfreq, nnphi, RFpower, &
!!$         EfieldNormalisation, freq, nphi, kperp
!!$
!!$    open( io_channel_3872, FILE='input.rfof')
!!$    read( io_channel_3872, input_wavefields)
!!$    close(io_channel_3872)
!!$
!!$    max_nnphi = nnphi
!!$
!!$    nion = 1
!!$    npsi = 1
!!$    ntheta = 1
!!$
!!$    max_npsi = npsi
!!$    max_ntheta = ntheta
!!$
!!$    !---- waves%global_param ----!
!!$
!!$    !print *, "ok tryin itm waves stuff...allocate"
!!$
!!$    allocate(waves%global_param%frequency(nfreq))
!!$    allocate(waves%global_param%name(nfreq))
!!$    allocate(waves%global_param%type(nfreq))
!!$    allocate(waves%global_param%nntor(nfreq))
!!$    allocate(waves%global_param%ntor(nfreq,max_nnphi))
!!$    ! ERROR in type description!!!            allocate(waves%global_param%f_assumption(nion+1))
!!$
!!$    allocate(waves%global_param%power_tot(nfreq))
!!$    allocate(waves%global_param%p_frac_ntor(nfreq,max_nnphi))
!!$    allocate(waves%global_param%pow_i(nfreq,nion))
!!$    allocate(waves%global_param%pow_e(nfreq))
!!$
!!$    allocate(waves%global_param%pow_ntor_i(nfreq,nnphi,nion))
!!$    allocate(waves%global_param%pow_ntor_e(nfreq,nnphi))
!!$    allocate(waves%global_param%cur_tor(nfreq))
!!$    allocate(waves%global_param%cur_tor_ntor(nfreq,nnphi))
!!$    allocate(waves%global_param%code_type(nfreq))
!!$    !   Something beam tracing - who care :) !!  allocate(waves%global_param%freq_point(nfreq))
!!$
!!$    do jfreq = 1, nfreq
!!$       waves%global_param%frequency(jfreq) = freq * (1d0+0.01d0*dble(jfreq-1))
!!$       waves%global_param%name(jfreq) = 'No antenna name'
!!$       waves%global_param%type(jfreq) = 'IC'
!!$       waves%global_param%nntor(jfreq) = nnphi
!!$       waves%global_param%power_tot(jfreq) = RFpower
!!$       waves%global_param%pow_e(jfreq) = 0d0
!!$       waves%global_param%cur_tor(jfreq) = 0d0
!!$       waves%global_param%code_type(jfreq) =  0
!!$
!!$       do jion = 1, nion
!!$          waves%global_param%pow_i(        jfreq,jion) =  RFpower / real(nion)
!!$       enddo
!!$
!!$       do jnphi = 1, waves%global_param%nntor(jfreq)
!!$          waves%global_param%ntor(         jfreq,jnphi) = nphi * jnphi
!!$          waves%global_param%pow_ntor_e(   jfreq,jnphi) = 0d0
!!$          waves%global_param%p_frac_ntor(  jfreq,jnphi) = 1d0 / real(nnphi)
!!$          waves%global_param%cur_tor_ntor( jfreq,jnphi) = 0d0
!!$
!!$          do jion = 1, nion
!!$             waves%global_param%pow_ntor_i(jfreq,jnphi,jion) = RFpower / real(nnphi * nion)
!!$          enddo
!!$       enddo
!!$    enddo
!!$    ! ERROR in type description!!!            waves%global_param%f_assumption(nion+1)
!!$
!!$    !print *, "itm_waves: global_param: nnphi", waves%global_param%nntor
!!$    !print *, "itm_waves: global_param: nphi", waves%global_param%ntor
!!$    !print *, "itm_waves: global_param: ptot", waves%global_param%power_tot
!!$    !print *, "itm_waves: global_param: j", waves%global_param%cur_tor_ntor
!!$    !print *, "itm_waves: global_param: freq", waves%global_param%frequency
!!$
!!$
!!$
!!$    !---- waves%fullwave%local ----!
!!$    allocate( waves%fullwave%local%e_plus(nfreq, max_nnphi, max_npsi,max_ntheta) )
!!$    allocate( waves%fullwave%local%e_minus(nfreq, max_nnphi, max_npsi,max_ntheta) )
!!$    allocate( waves%fullwave%local%e_para(nfreq, max_nnphi, max_npsi,max_ntheta) )
!!$
!!$    allocate( waves%fullwave%local%e_plus_ph(nfreq, max_nnphi, max_npsi,max_ntheta) )
!!$    allocate( waves%fullwave%local%e_minus_ph(nfreq, max_nnphi, max_npsi,max_ntheta) )
!!$    allocate( waves%fullwave%local%e_para_ph(nfreq, max_nnphi, max_npsi,max_ntheta) )
!!$
!!$    !---- waves%fullwave%local ----!
!!$    do jfreq = 1, nfreq
!!$       do jnphi = 1, waves%global_param%nntor(jfreq)
!!$          do jpsi = 1, max_npsi
!!$             do jtheta = 1, max_ntheta
!!$                 waves%fullwave%local%e_plus(    jfreq, jnphi, jpsi,jtheta) = 0.5d0
!!$                 waves%fullwave%local%e_minus(   jfreq, jnphi, jpsi,jtheta) = 0.5d0
!!$                 waves%fullwave%local%e_para(    jfreq, jnphi, jpsi,jtheta) = 0d0
!!$                 waves%fullwave%local%e_plus_ph( jfreq, jnphi, jpsi,jtheta) = 0d0
!!$                 waves%fullwave%local%e_minus_ph(jfreq, jnphi, jpsi,jtheta) = 0d0
!!$                 waves%fullwave%local%e_para_ph( jfreq, jnphi, jpsi,jtheta) = 0d0
!!$              enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    !print *, "e+ = ",waves%fullwave%local%e_plus
!!$    !print *, "e- = ",waves%fullwave%local%e_minus
!!$    !print *, "ph e+ = ",waves%fullwave%local%e_plus_ph
!!$    !print *, "ph e- = ",waves%fullwave%local%e_minus_ph
!!$
!!$
!!$  end subroutine construct_dummy_euitm_waves

  subroutine destruct_dummy_euitm_waves(waves)

    type(type_waves), intent(inout) :: waves

    ! Local
    integer :: jfreq
    integer :: nfreq

    nfreq = size(waves%coherentwave)

    !---- waves%global_param ----!
    do jfreq = 1, nfreq
       deallocate(waves%coherentwave(jfreq)%global_param%ntor)
       deallocate(waves%coherentwave(jfreq)%global_param%p_frac_ntor)
       deallocate(waves%coherentwave(jfreq)%global_param%pow_i)
       deallocate(waves%coherentwave(jfreq)%global_param%pow_ntor_i)
       deallocate(waves%coherentwave(jfreq)%global_param%pow_ntor_e)
       deallocate(waves%coherentwave(jfreq)%global_param%cur_tor_ntor)
    enddo


    !---- waves%fullwave%local ----!
    do jfreq = 1, nfreq
       deallocate( waves%coherentwave(jfreq)%fullwave%local%e_plus)
       deallocate( waves%coherentwave(jfreq)%fullwave%local%e_minus)
       deallocate( waves%coherentwave(jfreq)%fullwave%local%e_para)

       deallocate( waves%coherentwave(jfreq)%fullwave%local%e_plus_ph)
       deallocate( waves%coherentwave(jfreq)%fullwave%local%e_minus_ph)
       deallocate( waves%coherentwave(jfreq)%fullwave%local%e_para_ph)
    enddo

    deallocate(waves%coherentwave)


  end subroutine destruct_dummy_euitm_waves

end module euitm_waves_interface
