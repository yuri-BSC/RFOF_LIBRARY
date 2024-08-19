#include "config.h"

module RFOF_markers

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  use RFOF_constants
  use RFOF_magnetic_field
  use RFOF_parameters
  
  implicit none


  real(8), target :: RFOF_default_weight = 1



  type, public :: particle
     integer, pointer :: Id                  !< Number to identify particle
     real(8), pointer :: weight              !< Number of real particles represented by the marker
     real(8), pointer :: charge              !< Charge (SI units; everything is in SI)
     real(8), pointer :: mass                !< Mass (SI units; everything is in SI)

     real(8), pointer :: R                   !< Major radius
     real(8), pointer :: phi                 !< Toroidal angle
     real(8), pointer :: z                   !< Vertical position
     real(8), pointer :: psi                 !< Poloidal flux function

     real(8), pointer :: energy              !< Total energy, i.e. sum of kinetic and potential energy
     real(8), pointer :: energy_kinetic      !< Kinetic energy
     real(8), pointer :: velocity            !< magnitude of the velocity

     real(8), pointer :: magneticMoment      !< Magnetic moment
     real(8), pointer :: Pphi                !< The canonical toroidal angular momnetum \[f P_{\phi}=Ze\psi+mv_\|F/B \f]

     real(8), pointer :: vpar                !< Parallel velocity
     real(8), pointer :: vperp               !< Perpendicular velocity
     real(8), pointer :: omega_gyro          !< The angular gyro (Larmor) frequency
     real(8), pointer :: tauBounce           !< Estimate for the time of a poloidal orbit (bounce time for trapped / transit time for passing particles)
     real(8), pointer :: vDrift              !< The magnitude of the drift velocity
     real(8), pointer :: vDriftRho           !< The component of the drift velocity w.r.t. the radial      direction
     real(8), pointer :: vDriftDia           !< The component of the drift velocity w.r.t. the diamagnetic direction (diamagnetic=poloidal)
     real(8), pointer :: d_vpar_d_rho        !< Partial derivative of the parallel      velocity w.r.t. the radial      direction
     real(8), pointer :: d_vpar_d_dia        !< Partial derivative of the parallel      velocity w.r.t. the diamagnetic direction (diamagnetic=poloidal)
     real(8), pointer :: d_vperp_d_rho       !< Partial derivative of the perpendicular velocity w.r.t. the radial      direction
     real(8), pointer :: d_vperp_d_dia       !< Partial derivative of the perpendicular velocity w.r.t. the diamagnetic direction (diamagnetic=poloidal)
     real(8), pointer :: d_vDriftRho_d_rho   !< Partial derivative of the radial      component of the drift velocity w.r.t. the radial      direction
     real(8), pointer :: d_vDriftRho_d_dia   !< Partial derivative of the radial      component of the drift velocity w.r.t. the diamagnetic direction (diamagnetic=poloidal)
     real(8), pointer :: d_vDriftDia_d_rho   !< Partial derivative of the diamagnetic component of the drift velocity w.r.t. the radial      direction (diamagnetic=poloidal)
     real(8), pointer :: d_vDriftDia_d_dia   !< Factor by which the orbit integration has been accelerated (number of of real bounce times per calculated bounce in the orbit following.

     real(8), pointer :: time_acceleration
  end type particle

  type, public :: particle_static
     integer :: Id
     real(8) :: weight
     real(8) :: charge
     real(8) :: mass 

     real(8) :: R
     real(8) :: phi
     real(8) :: z
     real(8) :: psi 

     real(8) :: energy
     real(8) :: energy_kinetic
     real(8) :: velocity
     real(8) :: magneticMoment
     real(8) :: Pphi
     real(8) :: vpar
     real(8) :: vperp

     real(8) :: omega_gyro
     real(8) :: tauBounce

     real(8) :: vDrift
     real(8) :: vDriftRho
     real(8) :: vDriftDia
     real(8) :: d_vpar_d_rho
     real(8) :: d_vpar_d_dia
     real(8) :: d_vperp_d_rho
     real(8) :: d_vperp_d_dia
     real(8) :: d_vDriftRho_d_rho
     real(8) :: d_vDriftRho_d_dia
     real(8) :: d_vDriftDia_d_rho
     real(8) :: d_vDriftDia_d_dia
     real(8) :: time_acceleration
  end type particle_static


contains

  function copy_particle_pointer2pointer(pi) result(po)
 
    ! Input
    type(particle), intent(in) :: pi

    ! Output
    type(particle) :: po
    print *, "pointer2pointer"
    po%Id = pi%Id
    po%weight = pi%weight

    po%charge = pi%charge
    po%mass = pi%mass

    po%R = pi%R
    po%phi = pi%phi
    po%z = pi%z
    po%psi = pi%psi

    po%energy = pi%energy
    po%energy_kinetic = pi%energy_kinetic
    po%velocity = pi%velocity
    po%magneticMoment = pi%magneticMoment
    po%Pphi = pi%Pphi
    po%vpar = pi%vpar
    po%vperp = pi%vperp

    po%omega_gyro = pi%omega_gyro
    po%tauBounce = pi%tauBounce

    po%vDrift = pi%vDrift
    po%vDriftRho = pi%vDriftRho
    po%vDriftDia = pi%vDriftDia
    po%d_vpar_d_rho = pi%d_vpar_d_rho
    po%d_vpar_d_dia = pi%d_vpar_d_dia
    po%d_vperp_d_rho = pi%d_vperp_d_rho
    po%d_vperp_d_dia = pi%d_vperp_d_dia
    po%d_vDriftRho_d_rho = pi%d_vDriftRho_d_rho
    po%d_vDriftRho_d_dia = pi%d_vDriftRho_d_dia
    po%d_vDriftDia_d_rho = pi%d_vDriftDia_d_rho
    po%d_vDriftDia_d_dia = pi%d_vDriftDia_d_dia

    po%time_acceleration = pi%time_acceleration

  end function copy_particle_pointer2pointer

  function copy_particle_pointer2static(pi) result(po)

    ! Input
    type(particle), intent(in) :: pi

    ! Output
    type(particle_static) :: po
    print *, "pointer2static"
    po%Id = pi%Id
    po%weight = pi%weight

    po%charge = pi%charge
    po%mass = pi%mass

    po%R = pi%R
    po%phi = pi%phi
    po%z = pi%z
    po%psi = pi%psi

    po%energy = pi%energy
    po%energy_kinetic = pi%energy_kinetic
    po%velocity = pi%velocity
    po%magneticMoment = pi%magneticMoment
    po%Pphi = pi%Pphi
    po%vpar = pi%vpar
    po%vperp = pi%vperp

    po%omega_gyro = pi%omega_gyro
    po%tauBounce = pi%tauBounce

    po%vDrift = pi%vDrift
    po%vDriftRho = pi%vDriftRho
    po%vDriftDia = pi%vDriftDia
    po%d_vpar_d_rho = pi%d_vpar_d_rho
    po%d_vpar_d_dia = pi%d_vpar_d_dia
    po%d_vperp_d_rho = pi%d_vperp_d_rho
    po%d_vperp_d_dia = pi%d_vperp_d_dia
    po%d_vDriftRho_d_rho = pi%d_vDriftRho_d_rho
    po%d_vDriftRho_d_dia = pi%d_vDriftRho_d_dia
    po%d_vDriftDia_d_rho = pi%d_vDriftDia_d_rho
    po%d_vDriftDia_d_dia = pi%d_vDriftDia_d_dia

    po%time_acceleration = pi%time_acceleration

  end function copy_particle_pointer2static


  function copy_particle_static2static(pi) result(po)

    ! Input
    type(particle_static), intent(in) :: pi

    ! Output
    type(particle_static) :: po
    print *, "static2static"
    po%Id = pi%Id
    po%weight = pi%weight

    po%charge = pi%charge
    po%mass = pi%mass

    po%R = pi%R
    po%phi = pi%phi
    po%z = pi%z
    po%psi = pi%psi

    po%energy = pi%energy
    po%energy_kinetic = pi%energy_kinetic
    po%velocity = pi%velocity
    po%magneticMoment = pi%magneticMoment
    po%Pphi = pi%Pphi
    po%vpar = pi%vpar
    po%vperp = pi%vperp

    po%omega_gyro = pi%omega_gyro
    po%tauBounce = pi%tauBounce

    po%vDrift = pi%vDrift
    po%vDriftRho = pi%vDriftRho
    po%vDriftDia = pi%vDriftDia
    po%d_vpar_d_rho = pi%d_vpar_d_rho
    po%d_vpar_d_dia = pi%d_vpar_d_dia
    po%d_vperp_d_rho = pi%d_vperp_d_rho
    po%d_vperp_d_dia = pi%d_vperp_d_dia
    po%d_vDriftRho_d_rho = pi%d_vDriftRho_d_rho
    po%d_vDriftRho_d_dia = pi%d_vDriftRho_d_dia
    po%d_vDriftDia_d_rho = pi%d_vDriftDia_d_rho
    po%d_vDriftDia_d_dia = pi%d_vDriftDia_d_dia

    po%time_acceleration = pi%time_acceleration

  end function copy_particle_static2static

  !--------------------------------------------------------------------------------
  ! subroutine update_marker
  !--------------------------------------------------------------------------------
  subroutine update_marker(marker,Blocal)
 
    ! Input
    type(magnetic_field_local), intent(in) :: Blocal

    ! Input/Output
    type(particle), intent(inout) :: marker
    ! Local
    real(8) :: xi
    print *, "updatemarker"
    xi = sqrt( max(0d0 , (1.0 - Blocal%Bmod*marker%magneticMoment/(marker%energy - marker%charge*Blocal%psi_Estatic)) ) )

    ! Get marker
    marker%R   = Blocal%R
    marker%z   = Blocal%z
    marker%phi = Blocal%phi
    marker%psi = Blocal%psi
    marker%vpar  = xi * marker%velocity
    marker%vperp  = sqrt( max(marker%velocity**2 - marker%vpar**2 , 0d0) )
    print *, "Printeo el marker%mass antes del omega_gyro" , marker%mass
    print *, "Printejo el marker%charge antes del omega_gyro", marker%charge
    marker%omega_gyro = marker%charge * Blocal%Bmod / marker%mass

  end subroutine update_marker


  !--------------------------------------------------------------------------------
  ! subroutine set_marker_pointers_from_marker
  !--------------------------------------------------------------------------------
  subroutine set_marker_pointers_from_marker(m1,m2)

    ! Input
    type(particle_static), target, intent(in) :: m1

    ! Output
    type(particle), intent(out) :: m2
    print *, "set_marker_pointers_from_marker"
    m2%Id = m1%Id
    m2%weight => m1%weight
    m2%R => m1%R
    m2%phi => m1%phi
    m2%z => m1%z
    m2%psi =>m1%psi
    m2%charge => m1%charge
    m2%mass => m1%mass
    m2%energy => m1%energy
    m2%energy_kinetic => m1%energy_kinetic
    m2%velocity => m1%velocity
    m2%magneticMoment => m1%magneticMoment
    m2%Pphi => m1%Pphi
    m2%vpar => m1%vpar
    m2%vperp => m1%vperp
    m2%omega_gyro => m1%omega_gyro
    m2%tauBounce => m1%tauBounce
    m2%vDrift => m1%vDrift
    m2%vDriftRho => m1%vDriftRho
    m2%vDriftDia => m1%vDriftDia
    m2%d_vpar_d_rho => m1%d_vpar_d_rho
    m2%d_vpar_d_dia => m1%d_vpar_d_dia
    m2%d_vperp_d_rho => m1%d_vperp_d_rho
    m2%d_vperp_d_dia => m1%d_vperp_d_dia
    m2%d_vDriftRho_d_rho => m1%d_vDriftRho_d_rho
    m2%d_vDriftRho_d_dia => m1%d_vDriftRho_d_dia
    m2%d_vDriftDia_d_rho => m1%d_vDriftDia_d_rho
    m2%d_vDriftDia_d_dia => m1%d_vDriftDia_d_dia

    m2%time_acceleration => m1%time_acceleration

  end subroutine set_marker_pointers_from_marker


  !--------------------------------------------------------------------------------
  ! subroutine make_marker
  !--------------------------------------------------------------------------------
  subroutine set_marker_pointers(marker , &
       Id, weight, R, phi, z, psi , &
       charge, mass, energy, energy_kinetic, &
       velocity, magneticMoment, Pphi, vpar, &
       vperp, omega_gyro, tauBounce, vDrift, &
       vDriftRho,vDriftDia, d_vpar_d_rho, d_vpar_d_dia, d_vperp_d_rho, d_vperp_d_dia, &
       d_vDriftRho_d_rho, d_vDriftRho_d_dia, d_vDriftDia_d_rho, d_vDriftDia_d_dia,&
       time_acceleration)

    ! Input
    integer, target, intent(in) :: Id
    real(8), target, intent(in) :: weight,R,phi,z, psi, charge
    real(8), target, intent(in) :: mass
    real(8), target, intent(in) :: energy,energy_kinetic,velocity 
    real(8), target, intent(in) :: magneticMoment,Pphi,vpar,vperp,omega_gyro,tauBounce 
    real(8), target, intent(in) :: vDrift,vDriftRho,vDriftDia, d_vpar_d_rho, d_vpar_d_dia 
    real(8), target, intent(in) :: d_vperp_d_rho, d_vperp_d_dia, d_vDriftRho_d_rho, d_vDriftRho_d_dia
    real(8), target, intent(in) :: d_vDriftDia_d_rho, d_vDriftDia_d_dia
    real(8), target, intent(in) :: time_acceleration

    ! Output

    type(particle), intent(out) :: marker
    print *, "set_marker_pointers"
    marker%Id => Id
    marker%weight => weight
    marker%R => R
    marker%phi => phi
    marker%z => z
    marker%psi => psi
    marker%mass => mass
    marker%charge => charge
    marker%energy => energy
    marker%energy_kinetic => energy_kinetic
    marker%velocity => velocity
    marker%magneticMoment => magneticMoment
    marker%Pphi => Pphi
    marker%vpar => vpar
    marker%vperp => vperp
    marker%omega_gyro => omega_gyro
    marker%tauBounce => tauBounce
    marker%vDrift => vDrift
    marker%vDriftRho => vDriftRho
    marker%vDriftDia => vDriftDia
    marker%d_vpar_d_rho => d_vpar_d_rho
    marker%d_vpar_d_dia => d_vpar_d_dia
    marker%d_vperp_d_rho => d_vperp_d_rho
    marker%d_vperp_d_dia => d_vperp_d_dia
    marker%d_vDriftRho_d_rho => d_vDriftRho_d_rho
    marker%d_vDriftRho_d_dia => d_vDriftRho_d_dia
    marker%d_vDriftDia_d_rho => d_vDriftDia_d_rho
    marker%d_vDriftDia_d_dia => d_vDriftDia_d_dia

    marker%time_acceleration => time_acceleration 
    !print *, "In RFOF set_marker_pointers", marker%tauBounce, tauBounce
    !print *, "In RFOF set_marker_pointers", marker%R , R

  end subroutine set_marker_pointers

  !--------------------------------------------------------------------------------
  subroutine make_marker(marker, weight,charge,mass, psi,E,xi,tauBounce,Blocal)

    ! Input
    real(8), intent(inout) :: weight,charge,mass, psi, E,xi,tauBounce
    type(magnetic_field_local), intent(in) :: Blocal
   
    ! Output
    type(particle_static), intent(out) :: marker
    ! Get marker
    print *, "make_marker"
    marker%Id  = 1
    marker%weight = weight
    marker%R   = Blocal%R
    marker%z   = Blocal%z
    marker%phi = Blocal%phi
    marker%psi = psi   ! He definit això perquè no estava definit enlloc. 
    marker%charge  = charge * rfof_ev
    marker%mass  = mass * rfof_amu
    !marker%mass = 1.6d-27 
    marker%energy  = E
    marker%energy_kinetic  = E - marker%charge*Blocal%psi_Estatic
    marker%velocity  = sqrt(2*marker%energy_kinetic/marker%mass)
    marker%vpar  = xi * marker%velocity
    marker%vperp  = sqrt( marker%velocity**2 - marker%vpar**2 )
    marker%magneticMoment  = (1-xi**2) *  marker%energy / Blocal%Bmod
    marker%Pphi  = marker%charge * Blocal%psi + &
         (Blocal%F / Blocal%Bmod) * marker%mass * marker%vpar

    marker%omega_gyro = marker%charge * Blocal%Bmod / marker%mass
    marker%tauBounce = tauBounce

    marker%vDrift = 0.0
    marker%vDriftRho = 0.0
    marker%vDriftDia = 0.0
    marker%d_vpar_d_rho = 0.0
    marker%d_vpar_d_dia = 0.0
    marker%d_vperp_d_rho = 0.0
    marker%d_vperp_d_dia = 0.0
    marker%d_vDriftRho_d_rho = 0.0
    marker%d_vDriftRho_d_dia = 0.0
    marker%d_vDriftDia_d_rho = 0.0
    marker%d_vDriftDia_d_dia = 0.0

    marker%time_acceleration = 1.0
   
  end subroutine make_marker


  !--------------------------------------------------------------------------------
  subroutine validate_new_marker(marker, validityFlag) 
    type(particle), intent(in) :: marker
    integer, intent(out) :: validityFlag
    print *, "validate_new_marker"
    validityFlag = 0

    if ( (marker%energy_kinetic .lt. 0d0) .or. &
         (marker%magneticMoment .lt. 0d0) .or. &
         (marker%vperp          .lt. 0d0) ) then
       validityFlag = 1
       return
    endif

    if ( (marker%R .lt. plasma_boundingbox%Rmin) .or. &
         (marker%R .gt. plasma_boundingbox%Rmax) .or. &
         (marker%Z .lt. plasma_boundingbox%Zmin) .or. & 
         (marker%Z .gt. plasma_boundingbox%Zmax) ) then
       validityFlag = 2
       return
    endif

  end subroutine validate_new_marker

end module RFOF_markers
