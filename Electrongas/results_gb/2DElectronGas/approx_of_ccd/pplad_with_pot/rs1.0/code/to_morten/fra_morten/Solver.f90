
!
!     Class for Coupled-Cluster solver
!
MODULE SolverMod
  USE ConstantsMod
  USE EnergyMod
  USE MpiMod
  IMPLICIT NONE
  
  TYPE Solver
     TYPE(Energy), POINTER :: energy
   CONTAINS
     PROCEDURE :: calcEnergy
  END TYPE Solver
  
  
CONTAINS
  
  !
  !     Constructor for Solver object
  !
  FUNCTION Solver_(ene) RESULT(this)
    TYPE(Energy), TARGET, INTENT(IN) :: ene
    TYPE(Solver), POINTER :: this
    
    ALLOCATE(this)
    
    this%energy => ene
    
    
  END FUNCTION Solver_
  
  SUBROUTINE calcEnergy(this)
    CLASS(Solver) :: this
    
    !     Calculate the energy in the chosen
    !     approximation
    CALL this%energy%calculateEnergy()
    
    
  END SUBROUTINE calcEnergy
  
  FUNCTION readAndSetup() RESULT(mySolver)
    USE MemoryCounterMod
    USE Pt2CorrEnergyMod
    USE CcdCorrEnergyMod
    USE Yukawa2dHamMod
    USE Yukawa3dHamMod
    USE SpBasisPW2dMod
    USE SpBasisPW3dMod
    TYPE(Solver), POINTER :: mySolver
    CHARACTER(LEN=100) :: interaction
    CHARACTER(LEN=100) :: approximation
    CHARACTER(LEN=100) :: spBasisOption
    CHARACTER(LEN=100) :: spinState
    CLASS(Hamiltonian), POINTER :: myHamiltonian => NULL()
    CLASS(Yukawa2dHam), POINTER :: myYukawa2dHam => NULL()
    CLASS(Yukawa3dHam), POINTER :: myYukawa3dHam => NULL()
    TYPE(RefEnergy), POINTER :: myRefEnergy
    CLASS(CorrEnergy), POINTER :: myCorrEnergy => NULL()
    TYPE(Pt2CorrEnergy), POINTER :: myPt2CorrEnergy => NULL()
    TYPE(CcdCorrEnergy), POINTER :: myCcdCorrEnergy => NULL()
    TYPE(FockMatrix), POINTER :: myFockMatrix
    CLASS(SpBasis), POINTER :: mySpBasis => NULL()
    CLASS(SpBasisPW2d), POINTER :: mySpBasisPW2d => NULL()
    CLASS(SpBasisPW3d), POINTER :: mySpBasisPW3d => NULL()
    TYPE(BlockList), POINTER :: myBlockList
    TYPE(BlockListRe), POINTER :: myBlockListRe
    TYPE(Energy), POINTER :: myEnergy => NULL() 
    TYPE(Mpi), POINTER :: myMpiJobs => NULL()
    REAL(DP) :: tolerance, mu, rs, omega, area, volume, alpha
    INTEGER :: nOccupied, nShells, maxIter, nQNumbers, iError
    CHARACTER(LEN=30) :: inputFile
    INCLUDE 'mpif.h'

    
    CALL GETARG(1, inputFile)
    !     Read the input file
    OPEN(5,FILE=inputFile)
    READ(5,*); READ(5,*)
    READ(5,*) rs
    READ(5,*); READ(5,*); READ(5,*)
    READ(5,*) nOccupied, nShells
    READ(5,*); READ(5,*); READ(5,*); READ(5,*); READ(5,*)
    READ(5,*); READ(5,*)
    READ(5,*) spBasisOption
    READ(5,*); READ(5,*); READ(5,*); READ(5,*); READ(5,*)
    READ(5,*) spinState
    READ(5,*); READ(5,*); READ(5,*); READ(5,*); READ(5,*)
    READ(5,*); READ(5,*)
    READ(5,*) interaction
    READ(5,*); READ(5,*)
    READ(5,*) mu
    READ(5,*); READ(5,*); READ(5,*); READ(5,*)
    READ(5,*) approximation 
    READ(5,*); READ(5,*); READ(5,*)
    READ(5,*) tolerance
    READ(5,*); READ(5,*)
    READ(5,*) maxIter
    READ(5,*); READ(5,*)
    READ(5,*) alpha
    CLOSE(5)
    
    memoryUsed = 0
    memoryMax = 0
    
    !     Create an MPI object
    myMpiJobs => Mpi_()
    
    CALL myMpiJobs%writeMpi(' ')

    !     Create single-particle basis and
    !     single-particle pair objects   
    IF (spBasisOption == 'pw2d') THEN
       !     Plane waves in two dimensions
       ALLOCATE(mySpBasis, SOURCE=mySpBasisPW2d)
       
       area = nOccupied*pi*rs**2
       mySpBasis => SpBasisPW2d_(nOccupied, nShells, &
            spinState, area, myMpiJobs)
       
       nQNumbers = 3
       
    ELSEIF (spBasisOption == 'pw3d') THEN
       !     Plane waves in three dimensions
       ALLOCATE(mySpBasis, SOURCE=mySpBasisPW3d)
       
       volume = nOccupied*4.d0*pi*rs**3/3.d0
       mySpBasis => SpBasisPW3d_(nOccupied, nShells, &
            spinState, volume, myMpiJobs)
       
       nQNumbers = 4
       
    ENDIF
    
    CALL myMpiJobs%writeMpi(' Setting up single-particle basis...')
    !     Set up single-particle basis
    CALL mySpBasis%setupSpBasis()
    !     Set up list of (K_CM, M_S) blocks
    myBlockList => BlockList_(nQNumbers, mySpBasis)
    CALL myBlockList%setupBlockList()
    
    !     Wait until all processes are ready
    !CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
    !write(*,*) 'i am ',myMpiJobs%iAm
    
    myBlockListRe => BlockListRe_(nQNumbers, mySpBasis, myBlockList)
    CALL myBlockListRe%setupPhBlockList()
    
    !     Wait until all processes are ready
    !CALL MPI_BARRIER(MPI_COMM_WORLD, iError)    

    !     Create a Hamiltonian object
    IF (interaction == 'yukawa2d') THEN
       !     Yukawa interaction for two-dimensional
       !     periodic system
       ALLOCATE(myHamiltonian, SOURCE=myYukawa2dHam)
       
       myHamiltonian => Yukawa2dHam_(mySpBasis, &
            myBlockList, myBlockListRe, mu, rs, nOccupied)
       
    ELSEIF (interaction == 'yukawa3d') THEN
       !     Yukawa interaction for three-dimensional
       !     periodic system
       ALLOCATE(myHamiltonian, SOURCE=myYukawa3dHam)
       
       myHamiltonian => Yukawa3dHam_(mySpBasis, &
            myBlockList, myBlockListRe, mu, rs, nOccupied)
    ENDIF
    
    CALL myMpiJobs%writeMpi(' Setting up interaction matrix...')
    !     Set up interaction matrix in the given
    !     single-particle basis
    CALL setupInteractionMatr(myHamiltonian)
    CALL setupInteractionMatrRe(myHamiltonian)

    CALL myMpiJobs%writeMpi(' Setting up reference energy object')
    !     Create a reference energy object
    myRefEnergy => RefEnergy_(myHamiltonian, mySpBasis, &
         myBlockList)

    CALL myMpiJobs%writeMpi(' Setting up Fock matrix...')
    !     Create a single-particle potential object
    myFockMatrix => FockMatrix_(myHamiltonian, mySpBasis)
    !     Set up single-particle potential vector
    CALL setupFockVector(myFockMatrix)
    
    IF (approximation == 'pt2') THEN
       
       CALL myMpiJobs%writeMpi(' Setting up correlation energy object')
       ALLOCATE(myCorrEnergy, SOURCE=myPt2CorrEnergy)
       
       myCorrEnergy => Pt2CorrEnergy_(myHamiltonian,&
            myFockMatrix, mySpBasis)
       
    ELSEIF (approximation == 'ccd') THEN
       
       CALL myMpiJobs%writeMpi(' Setting up correlation energy object')
       ALLOCATE(myCorrEnergy, SOURCE=myCcdCorrEnergy)
       
       myCorrEnergy => CcdCorrEnergy_(myHamiltonian,&
            myFockMatrix, mySpBasis, maxIter, &
            tolerance, alpha)
    ENDIF
    
    
    !     Create an energy object
    myEnergy => Energy_(myHamiltonian, approximation,&
         myRefEnergy, myCorrEnergy)
    !     Create a solver object
    mySolver => Solver_(myEnergy)
    
    !
    !     Deallocation
    !
    
    !IF (ASSOCIATED(myCorrEnergy)) DEALLOCATE(myCorrEnergy)
    
!    IF (spBasisOption == 'pw2d') THEN
       
!       CALL SpBasisPW2d_d(mySpBasisPW2d)
!    ELSEIF (spBasisOption == 'pw3d') THEN
       
!       CALL SpBasisPW3d_d(mySpBasisPW3d)
!    ENDIF
!    DEALLOCATE(mySpBasis)

!    CALL BlockList_d(myBlockList)
!    CALL BlockListRe_d(myBlockListRe)
    
    CALL printMemoryUse(myMpiJobs)
    !myBlockListRe%index2Abij(indexBegin:indexEnd,:)
    
    !     Wait until all processes are ready
    CALL MPI_BARRIER(MPI_COMM_WORLD, iError)

    IF (myMpiJobs%iAm == 0) THEN
       WRITE(*,*) ' '
       WRITE(*,*) '------  Setup is done  ------'
       WRITE(*,*) ' '
    ENDIF
    
    
  END FUNCTION readAndSetup
  
  !
  !     Destructor for Solver objects
  !
  SUBROUTINE Solver_d(this)
    TYPE(Solver), TARGET, INTENT(INOUT) :: this
    CLASS(SpBasis), POINTER :: basis
    
    basis => this%energy%eneRef%spBasis
    CALL printMaxMemoryUsed(basis%mpiJobs)
    
    !     Destruct the Energy object
    CALL Energy_d(this%energy)
    DEALLOCATE(this%energy)
    
    !CALL printMemoryUse(mpiJobs)
    
    !     Complete the MPI jobs
    CALL Mpi_d()
    
    
  END SUBROUTINE Solver_d

  
END MODULE SolverMod
