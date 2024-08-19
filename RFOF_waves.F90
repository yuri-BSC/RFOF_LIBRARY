#include "config.h"

module RFOF_waves

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  use RFOF_constants
  use RFOF_parameters
  use euitm_waves_interface


  !use EZspline_obj  ! at the top
  !use EZspline

  implicit none

  type, public :: rf_wave_mode_global


  end type rf_wave_mode_global

  type, public :: rf_wave_global
     type(type_waves), pointer :: waves

     integer :: nfreq

     integer, allocatable :: nnphi(:)

     integer :: max_nnphi

     integer :: j_freq_present_wave = 1

     integer :: j_nphi_present_wave = 1

     real(8), allocatable :: x1(:,:)

     real(8), allocatable :: x2(:,:)

     real(8), allocatable :: krho(:,:,:,:)

     real(8), allocatable :: kdia(:,:,:,:)

     real(8), allocatable :: kpar(:,:,:,:)


     type( rf_wave_mode_global ), allocatable :: wave_modes(:,:)


     real(8), allocatable :: EnormalisationFactor(:,:)

     real(8), allocatable :: RFpower(:,:)
  end type rf_wave_global

  type, public :: rf_wave_local
     integer :: jfreq
     integer :: jnphi
     complex :: Erho
     complex :: Edia
     complex :: Epar
     complex :: Eplus
     complex :: Eminus
     real(8) :: krho
     real(8) :: kdia
     real(8) :: kpar
     real(8) :: kperp
     real(8) :: omega
     real(8) :: nphi
  end type rf_wave_local


contains

  !--------------------------------------------------------------------------------
  !  subroutine rf_wave_constructor(x1,x2,x3,Erho_norm,Edia_norm,Epar_norm, &
  !             krho,kdia,kpar,EfieldNormalisation,RFpower,omega,nphi,newWave)
  !
  subroutine rf_wave_constructor(RFglobal,x1,x2, &
       EfieldNormalisation,RFpower, &
       freq,nphi, &
       krho,kdia,kpar, &
       Erho_norm,Edia_norm,Epar_norm, &
       itm_waves_in)

    ! Input
    real(8), intent(in) :: x1(:,:)
    real(8), intent(in) :: x2(:,:)
    real(8), intent(in) :: krho(:,:,:,:)
    real(8), intent(in) :: kdia(:,:,:,:)
    real(8), intent(in) :: kpar(:,:,:,:)
    real(8), intent(in) :: EfieldNormalisation(:,:)
    real(8), intent(in) :: RFpower(:,:)
    real(8), intent(in) :: freq(:)
    integer, intent(in) :: nphi(:,:)
    complex, intent(in), optional :: Erho_norm(:,:,:,:)
    complex, intent(in), optional :: Edia_norm(:,:,:,:)
    complex, intent(in), optional :: Epar_norm(:,:,:,:)
    type(type_waves), target, optional :: itm_waves_in

    ! Output
    type(rf_wave_global), intent(out) :: RFglobal

    ! Local
    type(type_waves), save, target :: itm_waves_local
    integer :: n1, n2, nfreq, nnphi
    integer :: j1, j2, jf   , jn
    integer :: AllocateStatus
    real(8), allocatable :: temp_matrix(:,:)


    integer :: boundaryCondition1(2), boundaryCondition2(2)
    integer :: ier

    print *, "Vale"
    n1=size(x1,1)
    n2=size(x1,2)
    nfreq=size(freq)
    nnphi=size(nphi,2)

    allocate( RFglobal%nnphi(nfreq), STAT = AllocateStatus)

    allocate( RFglobal%x1(n1,n2), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    allocate( RFglobal%x2(n1,n2), STAT = AllocateStatus)

    allocate( RFglobal%krho(nfreq,nnphi,n1,n2), STAT = AllocateStatus)
    allocate( RFglobal%kdia(nfreq,nnphi,n1,n2), STAT = AllocateStatus)
    allocate( RFglobal%kpar(nfreq,nnphi,n1,n2), STAT = AllocateStatus)

    allocate( RFglobal%wave_modes(nfreq,nnphi), STAT = AllocateStatus)

    allocate(RFglobal%EnormalisationFactor(nfreq,nnphi), STAT = AllocateStatus)
    allocate(RFglobal%RFpower(nfreq,nnphi), STAT = AllocateStatus)

    RFglobal%nfreq=nfreq
    RFglobal%max_nnphi=nnphi
    do jn = 1, nfreq
       RFglobal%nnphi=nnphi
    enddo

    RFglobal%x1=x1
    RFglobal%x2=x2

    RFglobal%krho = krho
    RFglobal%kdia = kdia
    RFglobal%kpar = kpar

    RFglobal%EnormalisationFactor=EfieldNormalisation
    RFglobal%RFpower=RFpower

    if ( present(Erho_norm) .and.  present(Edia_norm) .and.  present(Epar_norm) ) then
       print *, "WAVES CONSTRUCTOR: make new itm_waves"
       call construct_dummy_euitm_waves(itm_waves_local)
       RFglobal%waves => itm_waves_local

       print *, "WAVE CONSTRUCTOR: ntor from waves%coh...(1):         ", &
            itm_waves_local%coherentwave(1)%global_param%ntor
       print *, "WAVE CONSTRUCTOR: ntor from RFglobal%waves%coh...(1):", &
            RFglobal%waves%coherentwave(1)%global_param%ntor
    endif
    if (present(itm_waves_in)) then
       print *, "WAVES CONSTRUCTOR: use input itm_waves"
       RFglobal%waves => itm_waves_in
    Endif

    allocate( temp_matrix(nfreq,nnphi) )
    temp_matrix(:,:) = Erho_norm(:,:,1,1)
    boundaryCondition1(1) = 1
    boundaryCondition1(2) = 1
    boundaryCondition2(1) = 1
    boundaryCondition2(2) = 1

!    call EZspline_init( RFglobal%wave_modes%Eplus_Real , n1,n2, boundaryCondition1, boundaryCondition2, ier)
!    call EZspline_error(ier)


    deallocate( temp_matrix )

    
  end subroutine rf_wave_constructor


  !--------------------------------------------------------------------------------
  ! subroutine rf_wave_destructor(wave)
  subroutine rf_wave_destructor(wave)

    type(rf_wave_global), intent(inout) :: wave
    integer :: ier

    if (allocated(wave%x1) ) then
       deallocate( wave%x1 )
    endif
    if (allocated(wave%x2) ) then
       deallocate( wave%x2 )
    endif
    if (allocated(wave%krho) ) then
       deallocate( wave%krho )
    endif
    if (allocated(wave%kdia) ) then
       deallocate( wave%kdia )
    endif
    if (allocated(wave%kpar) ) then
       deallocate( wave%kpar )
    endif
    if (allocated(wave%wave_modes)) then
!       call EZspline_free(wave%wave_modes%Eplus_Real, ier)
!       call EZspline_error(ier)
       deallocate( wave%wave_modes )
    endif

    call destruct_dummy_euitm_waves(wave%waves)

  end subroutine rf_wave_destructor

  !--------------------------------------------------------------------------------
  ! function get_rf_wave_local(x1,x2,x3,wi) result(wo)
  function get_rf_wave_local(x1,x2,x3,wi) result(wo)

    ! Input
    real(8), intent(in) :: x1
    real(8), intent(in) :: x2
    real(8), intent(in) :: x3
    type(rf_wave_global), intent(in) :: wi

    ! Output
    type(rf_wave_local) :: wo
    integer :: jf,jn
    real(8) :: ep,em,epar,ep_ph,em_ph,epar_ph

    jf = wi%j_freq_present_wave
    jn = wi%j_nphi_present_wave

    wo%jfreq = jf
    wo%jnphi = jn

    wo%omega = 2d0 * rfof_pi * wi%waves%coherentwave(jf)%global_param%frequency
    wo%nphi  = wi%waves%coherentwave(jf)%global_param%ntor(jn)

    ep      = wi%waves%coherentwave(jf)%fullwave%local%e_plus(    jn,1,1)
    em      = wi%waves%coherentwave(jf)%fullwave%local%e_minus(   jn,1,1)
    epar    = wi%waves%coherentwave(jf)%fullwave%local%e_para(    jn,1,1)
    ep_ph   = wi%waves%coherentwave(jf)%fullwave%local%e_plus_ph( jn,1,1)
    em_ph   = wi%waves%coherentwave(jf)%fullwave%local%e_minus_ph(jn,1,1)
    epar_ph = wi%waves%coherentwave(jf)%fullwave%local%e_para_ph( jn,1,1)

    wo%Eplus  = wi%EnormalisationFactor(jf,jn) * ep   * cmplx(cos(ep_ph  ),sin(ep_ph  ))
    wo%Eminus = wi%EnormalisationFactor(jf,jn) * em   * cmplx(cos(em_ph  ),sin(em_ph  ))
    wo%Epar   = wi%EnormalisationFactor(jf,jn) * epar * cmplx(cos(epar_ph),sin(epar_ph))

    wo%Erho   = (wo%Eplus - wo%Eminus)
    wo%Edia   = (wo%Eplus + wo%Eminus) * cmplx(1d0,0d0)

    wo%krho  = wi%krho(jf,jn,1,1)
    wo%kdia  = wi%kdia(jf,jn,1,1)
    if (simplify__kpar_is_nphi_over_R) then
       wo%kpar  = wo%nphi / x1
    else
       wo%kpar  = wi%kpar(jf,jn,1,1)
    endif
    wo%kperp = sqrt( wo%krho**2 + wo%kdia**2 )

!!$    print *, "get_rf_wave_local: norm=",wi%EnormalisationFactor(jf,jn)
!!$    print *, "get_rf_wave_local: j_nphi=", wi%j_nphi_present_wave
!!$    print *, "get_rf_wave_local: j_freq=", wi%j_freq_present_wave
!!$    print *, "get_rf_wave_local: omega",wi%omega , wo%omega, &
!!$         2.*rfof_pi*wi%waves%coherentwave(wi%j_freq_present_wave)%global_param%frequency
!!$
!!$    print *, "get_rf_wave_local: nphi",wi%nphi , wo%nphi, &
!!$         wi%waves%coherentwave(wi%j_freq_present_wave)%global_param%ntor(wi%j_nphi_present_wave)
!!$
!!$    print *, "get_rf_wave_local: wi%Eplus_norm",abs(wi%Eplus_norm),abs(wo%Eplus) , &
!!$         wi%waves%coherentwave(wi%j_freq_present_wave)%fullwave%local%e_plus(wi%j_nphi_present_wave,1,1)
!!$
!!$    print *, "get_rf_wave_local: wi%Eminus_norm",abs(wi%Eminus_norm), abs(wo%Eminus) , &
!!$         wi%waves%coherentwave(wi%j_freq_present_wave)%fullwave%local%e_minus(wi%j_nphi_present_wave,1,1)
!!$
!!$    print *, "get_rf_wave_local: E_rho", wi%Erho_norm*wi%EnormalisationFactor(jf,jn), &
!!$         wo%Erho
!!$    print *, "get_rf_wave_local: E_dia", wi%Edia_norm*wi%EnormalisationFactor(jf,jn), &
!!$         wo%Edia

  end function get_rf_wave_local


  !--------------------------------------------------------------------------------
  subroutine dummy_rf_wave_field(RFglobal)
!    print *, 'Check1'

    use RFOF_constants

    ! Input/Output
    type(rf_wave_global), intent(inout) :: RFglobal

    ! Local
    type(type_waves), save, target :: itm_waves
    integer,parameter :: N1=1
    integer,parameter :: N2=1

    integer :: j1, j2, jf, jn
    real(8) :: x1(N1,N2), x2(N1,N2)
    complex, allocatable :: Erho(:,:,:,:)
    complex, allocatable :: Edia(:,:,:,:)
    complex, allocatable :: Epar(:,:,:,:)
    real(8), allocatable :: krho(:,:,:,:)
    real(8), allocatable :: kdia(:,:,:,:)
    real(8), allocatable :: kpar(:,:,:,:)
    real(8), allocatable :: freq_vec(:)
    integer, allocatable :: nphi_vec(:,:)

    integer :: nnphi, nfreq, nphi
    real(8) :: RFpower, EfieldNormalisation, freq, kperp
    real(8), allocatable :: EfieldNormalisationMatrix(:,:), RFpowerMatrix(:,:)



    NAMELIST/input_wavefields/nfreq , nnphi , RFpower, &
        EfieldNormalisation, freq, nphi, kperp

    print *, 'Check1'
    open( io_channel_3872, FILE='input.rfof')
    read( io_channel_3872, input_wavefields)
    close(io_channel_3872)
    print *, 'wave field parameters=', nnphi,nfreq,RFpower,EfieldNormalisation,freq,nphi,kperp
    

    allocate( Erho(nfreq,nnphi,N1,N2))
    allocate( Edia(nfreq,nnphi,N1,N2))
    allocate( Epar(nfreq,nnphi,N1,N2))
    allocate( krho(nfreq,nnphi,N1,N2))
    allocate( kdia(nfreq,nnphi,N1,N2))
    allocate( kpar(nfreq,nnphi,N1,N2))

    allocate( freq_vec(nfreq))
    allocate( nphi_vec(nfreq,nnphi))
    allocate( EfieldNormalisationMatrix(nfreq,nnphi))
    allocate( RFpowerMatrix(nfreq,nnphi))
  
    do j1=1,N1
       do j2=1,N2
          x1(j1,j2)=real(j1-1)/real(max(N1-1,1))
          x2(j1,j2)=real(j2-1)/real(max(N2-1,1))
       end do
    end do

    do jf=1,nfreq
       do jn=1,nnphi
          freq_vec = 2d0 * rfof_pi * freq * (1d0+0.1d0*dble(jf-1))
          nphi_vec = nphi * jn
          EfieldNormalisationMatrix(jf,jn) = EfieldNormalisation
          RFpowerMatrix(jf,jn) = RFpower
       enddo
    enddo

    do jf=1,nfreq
       do jn=1,nnphi
          do j1=1,N1
             do j2=1,N2
                Erho(jf,jn,j1,j2)=cmplx(0d0, 0d0)
                Edia(jf,jn,j1,j2)=cmplx(1d0, 0d0)
                Epar(jf,jn,j1,j2)=cmplx(0d0 , 0d0)
                krho(jf,jn,j1,j2)=kperp
                kdia(jf,jn,j1,j2)=0d0
                kpar(jf,jn,j1,j2)=dble(nphi)/3d0
                !Edia(j1,j1,1)=cmplx(1d0+0.01*real(j1),0.01*real(j2))
             enddo
          enddo
       enddo
    enddo

    print *, "Que esta passant"

    call construct_dummy_euitm_waves(itm_waves)
    print *, "Check2"
    IF ( 0 .EQ. 1 ) THEN
       call rf_wave_constructor(RFglobal, x1,x2, &
            EfieldNormalisationMatrix,RFpowerMatrix, &
            freq_vec, nphi_vec, &
            krho,kdia,kpar, &
            itm_waves_in = itm_waves)
       print *, "DUMMY WAVE: ntor from itmwaves", itm_waves%coherentwave(1)%global_param%ntor
     ELSE
       call rf_wave_constructor(RFglobal, x1,x2, &
            EfieldNormalisationMatrix,RFpowerMatrix, &
            freq_vec, nphi_vec, &
            krho,kdia,kpar, &
            Erho_norm = Erho,Edia_norm = Edia,Epar_norm = Epar)
    ENDIF
    
    print *, "DUMMY WAVE: ntor from RFglobal", RFglobal%waves%coherentwave(1)%global_param%ntor

    deallocate( Erho )
    deallocate( Edia )
    deallocate( Epar )
    deallocate( krho )
    deallocate( kdia )
    deallocate( kpar )

    deallocate( freq_vec)
    deallocate( nphi_vec)

  end subroutine dummy_rf_wave_field

end module RFOF_waves
