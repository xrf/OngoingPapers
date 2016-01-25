
!
!     A general Coupled-Cluster code 
!
PROGRAM coupledCluster
  USE Pt2CorrEnergyMod
  USE CcdCorrEnergyMod
  USE ConstantsMod
  USE SolverMod
  USE TestsMod
  IMPLICIT NONE
  TYPE(Solver) :: mySolver
  
  !     Read the input file and set up
  !     a solver object
  mySolver = readAndSetup()
  
  !     Run the main program
  CALL mySolver%calcEnergy()
  !CALL testSpPotential(mySolver)
  !CALL pt2Sheperd
  
  
END PROGRAM coupledCluster
