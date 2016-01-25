
!
!     A coupled-cluster code for periodic systems.
!
!     Author: Gustav Baardsen, University of Oslo, Norway 
!
PROGRAM coupledCluster
  USE Pt2CorrEnergyMod
  USE CcdCorrEnergyMod
  USE ConstantsMod
  USE SolverMod
  USE TestsMod
  IMPLICIT NONE
  TYPE(Solver), POINTER :: mySolver
  REAL(DP) :: time1, time2, time3
  INCLUDE 'mpif.h'
  
  time1 = MPI_WTIME()
  
  !     Read the input file and set up
  !     a solver object
  mySolver => readAndSetup()
  
  time2 = MPI_WTIME()
  
  !     Run the main program
  CALL mySolver%calcEnergy()
  
  !CALL testSpPotential(mySolver)
  !CALL pt2Sheperd
  ! write(*,*) '## 6'
  ! CALL wait(10)
  
  time3 = MPI_WTIME()
  
  IF (mySolver%energy%hamilt%basisSp%mpiJobs%iAm == 0) THEN
     WRITE(*,*) 'Time for setup: ', time2-time1, ' seconds'
     WRITE(*,*) 'Time for calculations: ', time3-time2, ' seconds'
     WRITE(*,*) 'Total time: ', time3-time1, ' seconds'
  ENDIF
  
  !     Deallocate Solver object.
  CALL Solver_d(mySolver)
  DEALLOCATE(mySolver)
  
  
  
END PROGRAM coupledCluster
