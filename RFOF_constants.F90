

#include "config.h"

module RFOF_constants

#ifdef USE_ISO_C_BINDING
  use, intrinsic :: ISO_C_BINDING
#endif

  real (kind = 8), parameter :: rfof_pi = 3.141592653589793238462643383280d0
  real (kind = 8), parameter :: rfof_c = 2.99792458d8         ! speed of light, m/s
  real (kind = 8), parameter :: rfof_me = 9.10938215d-31      ! electron mass, kg
  real (kind = 8), parameter :: rfof_mp = 1.672621637d-27     ! proton mass, kg
  real (kind = 8), parameter :: rfof_md = 3.34358320d-27      ! deuteron mass, kg
  real (kind = 8), parameter :: rfof_mt = 5.00735588d-27      ! triton mass, kg
  real (kind = 8), parameter :: rfof_ma = 6.64465620d-27      ! alpha mass, kg
  real (kind = 8), parameter :: rfof_amu = 1.660538782d-27    ! amu, kg
  real (kind = 8), parameter :: rfof_ev = 1.602176487d-19
  real (kind = 8), parameter :: rfof_Mev = 1.602176487d-13
  real (kind = 8), parameter :: rfof_qe = rfof_ev
  real (kind = 8), parameter :: rfof_mu0 = 4.0d-7 * rfof_pi
  real (kind = 8), parameter :: rfof_eps0 = 1.d0 / (rfof_mu0 * rfof_c * rfof_c)
  real (kind = 8), parameter :: rfof_avogr = 6.02214179d23
  real (kind = 8), parameter :: rfof_KBolt = 1.3806504d-23
  character (len=64), parameter :: rfof_constants_version = '$Id$'

  INTEGER,  PARAMETER :: itm_int_invalid = -999999999
  REAL(8), PARAMETER :: itm_r8_invalid = -9.0D40

end module RFOF_constants

