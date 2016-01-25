
MODULE ConstantsMod
  IMPLICIT NONE

  INTEGER,  PARAMETER :: dp = KIND(1.0D0)  
  REAL(DP), PARAMETER :: pi = 3.141592653589793_dp 
  !     hbarc = Planck's constant * velocity of light,
  !             given in MeV*fm. 2012 value from NIST.
  REAL(DP), PARAMETER, PUBLIC :: hbarc = 197.3269718 
 

END MODULE ConstantsMod
