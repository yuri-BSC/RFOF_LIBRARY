#include "config.h"

module RFOF_kick

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  use RFOF_constants
  use RFOF_parameters
  use RFOF_waves
  use RFOF_magnetic_field
  use RFOF_markers
  use RFOF_resonance_memory
  use RFOF_resonance_condition
  use RFOF_random_numbers
  use RFOF_diagnostics

  implicit none

  type :: coeff_for_guiding_centre_kick
     real(8) :: d_mu_dI            !<  \f$\partial \mu    / \partial I_\perp\f$
     real(8) :: d_vpar_dI          !<  \f$\partial v_\|   / \partial I_\perp\f$
     real(8) :: d_psi_dI           !<  \f$\partial \psi   / \partial I_\perp\f$
     real(8) :: d_theta_dI         !<  \f$\partial \theta / \partial I_\perp\f$
  end type coeff_for_guiding_centre_kick

contains

  !--------------------------------------------------------------------------------
  !
  ! subroutine quasilinear_RF_kick_steinbrecher_integrator()
  !
  subroutine quasilinear_RF_kick_steinbrecher_integrator(marker,mem, &
       Blocal,RFlocal,RFglobal,diagno,MPI_node_Id)

    use RFOF_types

    ! Input
    type(rf_wave_global), target, intent(in) :: RFglobal
    type(resonance_memory), target, intent(in) :: mem
    integer, intent(in) :: MPI_node_Id

    ! Input/Output
    type(particle), target, intent(inout) :: marker
    type(magnetic_field_local), target, intent(inout) :: Blocal
    type(rf_wave_local), target, intent(inout) :: RFlocal
    type(RFOF_cumlative_diagnostics), intent(inout) :: diagno

    ! Local
    type(RFOF_state) :: state
    real(8) :: diffusion, drift
    real(8) :: Iperp, dt
    real(8) :: energy_old, Pphi_old
    character(128) :: RFopType
    real(8) :: relative_tolerance_MonteCarlo = 1d-3
    integer :: validityFlag, errorFlag

    !--------------------------------------------------------------------------------
    interface
       subroutine MC_kick_steinbrecher_integrator( &
            dt, &
            x, &
            calc_diffusion_coeff, &
            state, &
            relative_tolerance_MonteCarlo, &
            errorFlag, &
            MPI_node_Id)
         use RFOF_types
         real(8), intent(in) :: dt                 !> Time step
         real(8), intent(inout) :: x               !> Marker position
         type(RFOF_state), intent(inout) :: state  !> State of RFOF
         real(8) :: relative_tolerance_MonteCarlo  !> Relative tolerance in Monte Carlo stepping
         integer, intent(inout) :: errorFlag       !> Error Flag
         integer, intent(in) :: MPI_node_Id        !> Number to identify the MPI node

         interface
            subroutine calc_diffusion_coeff(x,state,diffusion,errorFlag)
              use RFOF_types
              real(8), intent(in) :: x                  !< Integration variable
              type(RFOF_state), intent(inout) :: state  !< State of RFOF
              integer, intent(out) :: errorFlag         !< Error flag
              real(8), intent(out) :: diffusion         !< diffusion coefficient
            end subroutine calc_diffusion_coeff
         end interface

       end subroutine MC_kick_steinbrecher_integrator
    end interface
    !--------------------------------------------------------------------------------

    if (mem%nharm(1).eq.0) then
       return
    endif

    ! Store initial energy and momentum for diagnostics
    energy_old = marker%energy
    Pphi_old = marker%Pphi

    ! The acceleration is parameterised in terms of Iperp
    Iperp = marker%magneticMoment * &
         ( marker%mass * RFlocal%omega) / ( marker%charge * real(mem%nharm(1)) ) / rfof_Mev
    ! The time step, dt, is unity since the diffusion coefficient give the variance of the total kick
    dt = 1d0

    ! Generate the RFOF-state
    state%marker   => marker
    state%mem      => mem
    state%Blocal   => Blocal
    state%RFlocal  => RFlocal
    state%RFglobal => RFglobal
    state%Iperp    =  Iperp

    ! Integrate stochastic equation for marker. Marker position parameterized by Iperp
    call MC_kick_steinbrecher_integrator( &
         dt, &
         Iperp, &
         wrap_RF_diffusion_and_move_marker, &  ! Function: gives the diffusion coefficient
         state, &
         relative_tolerance_MonteCarlo, &
         errorFlag, &
         MPI_node_Id)

    Iperp = max( abs(Iperp) , 1d-30 )

    call wrap_move_marker_on_characteristic(Iperp,state,errorFlag)

    call validate_new_marker(marker, validityFlag)
    if (validityFlag .gt. 0) then
       write(0,*)'------------------------------------------------------------'
       write(0,*)'WARNING in quasilinear_RF_kick_with_time_acceleration'
       write(0,*)'    Marker invalid; validityFlag=',validityFlag
       write(0,*)'------------------------------------------------------------'
    endif

    ! Update diagnostics
    diffusion = 0d0 ! WARNING
    drift = 0d0 ! WARNING
    RFopType = 'MC_kick_steinbrecher_integrator'
    call add_kick_info_to_diagnostics(diagno, marker, mem, &
         RFlocal%jfreq, RFlocal%jnphi, &
         marker%energy - energy_old, marker%Pphi - Pphi_old, &
         drift, diffusion, RFopType, MPI_node_Id)

  end subroutine quasilinear_RF_kick_steinbrecher_integrator


  !--------------------------------------------------------------------------------
  !
  ! subroutine wrap_RF_diffusion_and_move_marker
  !
  subroutine wrap_RF_diffusion_and_move_marker(Iperp,state,diffusion,errorFlag)

    use RFOF_types

    ! Input
    real(8), intent(in) :: Iperp              !< Magnetic moment normalised for RF
    type(RFOF_state), intent(inout) :: state  !< State of RFOF (Note: may not be latest state - in case state%Iperp is different from Iperp, then state needed updating)
    integer, intent(out) :: errorFlag         !< Error flag
    real(8), intent(out) :: diffusion          !< diffusion coefficient

    call wrap_move_marker_on_characteristic(Iperp,state,errorFlag)

    diffusion = quasilinear_RF_diffusion_coeff(state%marker,state%RFlocal,state%mem)

  end subroutine wrap_RF_diffusion_and_move_marker

  !--------------------------------------------------------------------------------
  ! function quasilinear_RF_diffusion_coeff(marker,RFlocal,mem) result(diffusion)
  !
  function quasilinear_RF_diffusion_coeff(marker,RFlocal,mem) result(diffusion)

    ! Input
    type(particle), intent(in) :: marker
    type(rf_wave_local), intent(in) :: RFlocal
    type(resonance_memory), intent(in) :: mem

    ! Output
    real(8) :: diffusion

    !    real(8) :: Iperp

    diffusion = ( (marker%mass * marker%vperp  / rfof_Mev) * &
         max_RFkick_in_single_pass(marker,RFlocal,mem) )**2

  end function quasilinear_RF_diffusion_coeff

  !--------------------------------------------------------------------------------
  ! subroutine max_RFkick_in_single_pass()
  !
  function max_RFkick_in_single_pass(marker,wave,mem) result(max_kick)

    ! Input
    type(particle), intent(in) :: marker
    type(rf_wave_local), intent(in) :: wave
    type(resonance_memory), intent(in) :: mem

    ! Output
    real(8) :: max_kick

    ! Local
    real(8) :: kperp_rho, Jnm1, Jnp1, e_eff
    real(8) :: tau

    ! Bessel functions
    if (simplify__assume_zero_larmor_radius_in_KPERPxRHO) then
       if (mem%nharm(1) .eq. 1) then
          Jnm1=1d0
       else
          Jnm1=0d0
       endif
       Jnp1=0d0
    else
       kperp_rho = wave%kperp * marker%vperp / marker%omega_gyro
       Jnm1 = BESSELJ( mem%nharm(1)-1 , kperp_rho)
       Jnp1 = BESSELJ( mem%nharm(1)+1 , kperp_rho)
       ! print *, 'bessel function J0,J2,kperp*rho=',Jnm1,Jnp1,kperp_rho
    endif

    ! Electric wave field
    e_eff = abs(wave%Eplus * Jnm1 + wave%Eminus * Jnp1)

    ! Phase integral
    tau = tau_RF_phase_integral(mem%time,mem%omega_res,mem%Number_points_in_memory)

    ! Diffusion coefficient
    max_kick    = (marker%charge / marker%mass) * e_eff * tau

    !print *, "kick param: Eeff[V/m]=", e_eff, ", tau[ns]=", tau*1e9, ", kick=", max_kick

  end function max_RFkick_in_single_pass

  !--------------------------------------------------------------------------------
  ! function tau_RF_phase_integral()
  !
  function tau_RF_phase_integral(time,omega_res,Nelements) result(tau)

    ! Input
    !    type(resonance_memory), intent(in) :: mem
    integer, intent(in) :: Nelements
    real(8), intent(in) :: time(Nelements),omega_res(Nelements)

    ! Output
    real(8) :: tau

    ! Local
    real(8) d_phv, d_phv2, dd_phv, arg, power_third_dd_phv, power_third_2

    if (Nelements .lt. 2) then
       tau = 0.
       return
    endif

    d_phv = abs(omega_res(1) - omega_res(2)) &
         / (abs(time(1)      - time(2)) + 1d-30)

    if (Nelements .lt. 3) then
       d_phv2 = d_phv
       dd_phv = 0.
       tau = sqrt( 2.0 * rfof_pi / abs(d_phv) )
       return
    else
       d_phv2  = (omega_res(2) - omega_res(3)) &
            / (time(2)      - time(3))
       dd_phv = abs((d_phv - d_phv2) / (0.5*(time(1)-time(3))))
    endif

    if ( abs(d_phv)**3 .gt. dd_phv**2 ) then
       tau = sqrt( 2.0 * rfof_pi / d_phv )
    else
       power_third_dd_phv = abs(dd_phv)**(1./3.)
       power_third_2 = 2.0**(1./3.)
       arg = -d_phv**2 / ( power_third_2**2 * power_third_dd_phv**4 )
       tau = 2. * rfof_pi * ( power_third_2 / power_third_dd_phv ) * AIRY_SPECIAL_ARG(arg)
    endif

  end function tau_RF_phase_integral

  !--------------------------------------------------------------------------------
  ! subroutine coeff_RF_characteristic()
  !
  function coeff_RF_characteristic(marker,Blocal,RFlocal,nharm) result(coeff)

    ! Input
    type(particle), intent(in) :: marker
    type(rf_wave_local), intent(in) :: RFlocal
    type(magnetic_field_local), intent(in) :: Blocal
    integer, intent(in) :: nharm

    ! Output
    type(coeff_for_guiding_centre_kick) :: coeff

    ! Local
    real(8) :: vpar

    vpar = marker%vpar
    if (abs(vpar) < marker%velocity * 1e-6) then
       if (vpar>0) then
          vpar = marker%velocity * 1e-6
       else
          vpar = -marker%velocity * 1e-6
       endif
    endif

    coeff%d_mu_dI = rfof_Mev * marker%charge * nharm / ( marker%mass * RFlocal%omega)

    coeff%d_vpar_dI = ( rfof_Mev / ( marker%mass * RFlocal%omega ) ) * &
         ( RFlocal%kpar + ( RFlocal%krho * marker%vDriftRho   + &
                            RFlocal%kdia * marker%vDriftDia ) / vpar &
         )

    if (simplify__static_resonance_position_during_RF_kick) then
       coeff%d_psi_dI = 0.0
       coeff%d_theta_dI = 0.0
    else
       write(0,*) "Not yet implemented!"
       coeff%d_psi_dI = 0.0
       coeff%d_theta_dI = 0.0
    endif

  end function coeff_RF_characteristic

  !--------------------------------------------------------------------------------
  ! subroutine wrap_move_marker_on_characteristic
  !
  subroutine wrap_move_marker_on_characteristic(Iperp,state,errorFlag)

    use RFOF_types
    use RFOF_local_magnetic_field, only: get_local_magnetic_field

    ! Input
    real(8), intent(in) :: Iperp
    type(RFOF_state), intent(inout) :: state
    integer, intent(out) :: errorFlag

    ! Local
    type(coeff_for_guiding_centre_kick) :: coeff
    real(8) :: dIperp
    real(8) :: R,z

    ! The marker is only physical for positive Iperp
    if (Iperp .lt. 0d0) then
       errorFlag = 1
    endif
    errorFlag = 0 ! Default; might be changed later in the routine

    R=state%marker%R
    z=state%marker%z

!!$    write(0,*)"|----------------------------"
!!$    write(0,*)"| MOVE MARKER"
!!$    write(0,*)"|"
!!$    write(0,*)"| dIperp",dIperp
!!$    write(0,*)"| mu", state%marker%magneticMoment

    dIperp = Iperp - state%Iperp
    coeff = coeff_RF_characteristic(state%marker,state%Blocal,state%RFlocal,state%mem%nharm(1))
    call move_marker_on_characteristic(state%marker,dIperp,state%Blocal,coeff)
    state%Iperp = Iperp

    ! If the marker has moved; recalculate the magnetic and the RF wave fields
    if ((abs(R-state%marker%R).gt.1e-8) .or. (abs(z-state%marker%z).gt.1e-8)) then
       state%Blocal  = get_local_magnetic_field(state%marker%R,state%marker%phi,state%marker%z)
       state%RFlocal = get_rf_wave_local(state%marker%R,state%marker%phi,state%marker%z,state%RFglobal)
    endif

!!$    write(0,*)"| dIperp",dIperp
!!$    write(0,*)"| mu", state%marker%magneticMoment
!!$    write(0,*)"|----------------------------"

  end subroutine wrap_move_marker_on_characteristic


  !--------------------------------------------------------------------------------
  ! subroutine move_marker_on_characteristic
  !
  subroutine move_marker_on_characteristic(marker,dIperp,Blocal,coeff)

    ! Input
    real(8) :: dIperp
    type(magnetic_field_local), intent(in) :: Blocal
    type(coeff_for_guiding_centre_kick), intent(in) :: coeff

    ! Input/Output
    type(particle), intent(inout) :: marker

    ! Local
    real(8) :: dvpar, dmu, vpar, mu, vperp2

    ! Steps in \f$ \mu \f$ and \f$ v_\perp \f$
    dvpar = coeff%d_vpar_dI * dIperp
    dmu = coeff%d_mu_dI * dIperp

    if (simplify__static_resonance_position_during_RF_kick .eqv. .false.) then
       write(0,*)"ERROR! Option simplify__static_resonance_position_during_RF_kick==false", &
            " has not been implemented"
       write(0,*)"Abort"
       stop
    endif

    vpar = marker%vpar + dvpar
    mu = marker%magneticMoment + dmu
    if (mu .lt. 0d0) then
       write(0,*) "WARNING in move_marker_on_characteristic! ", &
            "Negative magnetic moment, mu=",mu
       mu=0d0
    endif
    vperp2 = mu * 2d0 * Blocal%Bmod / marker%mass

    ! Update marker
    marker%vpar = vpar
    marker%vperp = sqrt(vperp2)
    marker%velocity = sqrt( vperp2 + vpar**2 )
    marker%energy_kinetic = 0.5d0 * marker%mass * ( vperp2 + vpar**2)
    marker%energy = marker%energy_kinetic + marker%charge*Blocal%psi_Estatic
    marker%magneticMoment = mu
    marker%Pphi = marker%charge * Blocal%psi + &
         (Blocal%F / Blocal%Bmod) * marker%mass * marker%vpar

    !print *, "mv iperp=",dIperp,  marker%magneticMoment * &
    !     ( marker%mass *  337702360.70498121d0) / ( marker%charge * 1d0 ) / rfof_Mev
  end subroutine move_marker_on_characteristic

  !--------------------------------------------------------------------------------
  ! get_dvpar_on_characteristic
  !--------------------------------------------------------------------------------
  subroutine get_dvpar_on_characteristic(marker, Blocal, RFlocal, dvperp2, dvpar)

    ! Input
    type(particle), intent(in) :: marker
    type(magnetic_field_local), intent(in) :: Blocal
    type(rf_wave_local), intent(in) :: RFlocal
    real(8), intent(in) :: dvperp2

    ! Input/Output
    real(8), intent(out) :: dvpar

    ! Local
    real(8) :: coeff_dPphi_over_dvpar
    real(8) :: coeff_denergy_dmvpar
    real(8) :: p_quadratic_eq
    real(8) :: q_quadratic_eq
    real(8) :: Discriminant
    real(8) :: vperp, vperp2

    vperp = marker%vperp
    vperp2=vperp**2

    if (int(RFlocal%nphi+1d-10) .eq. 0) then
       dvpar = 0d0
    else
       coeff_dPphi_over_dvpar = -marker%mass * Blocal%F / Blocal%Bmod
       coeff_denergy_dmvpar = Blocal%F * RFlocal%omega / (Blocal%Bmod * dble(RFlocal%nphi))

       ! Put equation on normal form: x^2 - 2px + q =0, with x=dvpar/vperp (where) 
       p_quadratic_eq = - (coeff_denergy_dmvpar - marker%vpar) / vperp
       q_quadratic_eq = dvperp2 / vperp2
       Discriminant  = p_quadratic_eq**2 - q_quadratic_eq
       if (Discriminant .lt. 0d0) then
          write(0,*)"WARNING in RFkick_single_pass: Discriminant=",Discriminant," in the vpar-equation is negative"
          Discriminant = 0d0
       endif
       dvpar = ( p_quadratic_eq - sign(1d0,p_quadratic_eq) * sqrt(Discriminant) ) ! Quadrativ solution
    endif

  end subroutine get_dvpar_on_characteristic


  !--------------------------------------------------------------------------------
  function timestep_factor(x_a, x_c, state_a) result(timestepFactor)
    use RFOF_types
    real(8), intent(in):: x_a, x_c
    type(RFOF_state), intent(in) :: state_a
    real(8)  :: timestepFactor

    timestepFactor = 1d0

  end function timestep_factor


  !--------------------------------------------------------------------------------
  function close_to_boundary(x, state) result(closeToBoundary)
    use RFOF_types
    real(8), intent(in)::   x
    type(RFOF_state), intent(in) ::   state
    logical ::closeToBoundary

    closeToBoundary = .FALSE.

  end function close_to_boundary


end module RFOF_kick
