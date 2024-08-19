#include "../src/config.h"

program rftest

  use dummy_orbit
  implicit none

  print *, " ---dumorb--- (trace dummy orbit and give RF kicks)"
  call run_dumorb
  print *, " ---DONE dumorb--- (trace dummy orbit and give RF kicks)"

  print *, " "
  print *, " --- TEST RFOF ROUTINES ---"
  call my_tests
  print *, " --- DONE: TEST RFOF ROUTINES ---"

  print *, "The end of the program is reached here"

contains

  subroutine my_tests

  use RFOF_Efield_update
!  use mytests
!  use tst_rf_memory
!  use TEST_Monte_Carlo_Stepping_for_RF


  !call enorm_stats_test

  ! call test_poly_stuff
  ! call TEST_Monte_Carlo_Stepping_for_RF_1
  !print *, " ---ok!---"

  end subroutine my_tests
  
end program rftest
