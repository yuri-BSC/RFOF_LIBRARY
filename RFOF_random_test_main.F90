program GaussGeneratorTest 
use RFOF_random_numbers_SG2
use RFOF_random_numbers
implicit none
intrinsic CPU_TIME 
real(8):: Ntest
REAL time_begin, time_end

Print *, 'Start' 
!-------------------------------------------------------------------------------
 print*, 'Accuracy and speed test, by using RM48() uniform random number generator'
CALL CPU_TIME ( time_begin )
call init1_RM48(0)
PRINT*, 'Time  start= ', time_begin  , ' seconds'
Ntest=1.0d4 
  print *,'Start accuracy test of random number generators. Nr iterations=',Ntest
  print *, '  '

 call  PolarGaussAccuracyTest1(Ntest)
CALL CPU_TIME ( time_end )
print *, 'END PolarGaussAccuracyTest1  '
 PRINT *, ' POLAR method1 accuracy test CPU time= for ', time_end - time_begin, ' seconds'
 print *, '  '
 print *, '  '
CALL CPU_TIME ( time_begin )
call  MarsagliaBrayAccuracyTest1(Ntest)
CALL CPU_TIME ( time_end )
   PRINT *, ' MarsagliaBray 1 method CPU time, accuracy test= ', time_end - time_begin, ' seconds'
   print *,' END accuracy test1 of random number generators' 
  print *, '  '
  print *, '  ' 
print *,' Start speed test 1 of random number generators Nr iterations=',Ntest 
CALL CPU_TIME ( time_begin )
call PolarGaussSpeedTest1(Ntest)
CALL CPU_TIME ( time_end )
 print *, 'END Polar Gauss Speed Tes1t  '
PRINT *, ' POLAR method1 CPU time=  ', time_end - time_begin, ' seconds'
print *, '  '
print *, '  '
CALL CPU_TIME ( time_begin )
call  MarsagliaBraySpeedTest1(Ntest)
CALL CPU_TIME ( time_end )
PRINT *, ' MarsagliaBray method1 CPU time= ', time_end - time_begin, ' seconds'
 print *,' END Speed test1 of random number generators'  
 print *, '  '
 print *, '  '
!________________________________________________________________

!-------------------------------------------------------------------------------
 print*, 'Accuracy and speed test 2, by using RANDOM_NUMBER uniform random number generator'
call init_RANDOM_NUMBER(0)
PRINT*, 'Time  start= ', time_begin  , ' seconds'
!Ntest=1.0d6  
 print *,'Start accuracy test 2 of random number generators. Nr iterations=',Ntest
CALL CPU_TIME ( time_begin )
 call  PolarGaussAccuracyTest2(Ntest)
CALL CPU_TIME ( time_end )
 print *, 'END PolarGaussAccuracyTest2 tests'
 PRINT *, ' POLAR method2 CPU time=  ', time_end - time_begin, ' seconds'
  print *,''
 CALL CPU_TIME ( time_begin )
  call  MarsagliaBrayAccuracyTest2(Ntest)
CALL CPU_TIME ( time_end )
   PRINT *, ' MarsagliaBray2 method CPU time, accuracy test= ', time_end - time_begin, ' seconds'
   print *,' END accuracy test2 of random number generators'
   print *,''
   print *,''
   print *,''
!Ntest=1.0d7 
CALL CPU_TIME ( time_begin )
print *,' *********************************'
PRINT *, ' START correlated_GaussTest  '
 call correlated_GaussTest(Ntest)
 CALL CPU_TIME ( time_end )

   PRINT *, ' correlated_GaussTest  CPU time, accuracy test= ', time_end - time_begin, ' seconds'
   print *,' END accuracy test2 of random number generators'  
   print *,' *********************************'
  
  Ntest=1.d5
 print *,' ' 
print *,' Start speed test 2 of random number generators Nr iterations=',Ntest 
print *,' ' 
CALL CPU_TIME ( time_begin )
call PolarGaussSpeedTest2(Ntest) 
CALL CPU_TIME ( time_end )
print *, 'END Polar Gauss Speed Test2  '
PRINT *, ' POLAR method2 CPU time=  ', time_end - time_begin, ' seconds'
CALL CPU_TIME ( time_begin )
call  MarsagliaBraySpeedTest2(Ntest)
CALL CPU_TIME ( time_end )
PRINT *, ' MarsagliaBray method1 CPU time= ', time_end - time_begin, ' seconds'

print *,' END Speed test2 of random number generators'
!________________________________________________________________

stop ' END Gauss Generator tests  '
end program  GaussGeneratorTest 
