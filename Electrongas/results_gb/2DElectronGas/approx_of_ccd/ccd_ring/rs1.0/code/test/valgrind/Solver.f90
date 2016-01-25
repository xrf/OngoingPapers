
!
!     Class for Coupled-Cluster solver
!
MODULE SolverMod
  USE ConstantsMod
  USE EnergyMod
  IMPLICIT NONE
  
  TYPE Solver
     TYPE(Energy) :: energy      
   CONTAINS
     PROCEDURE :: calcEnergy
  END TYPE Solver
  
  
CONTAINS
  
  
  SUBROUTINE calcEnergy(this)
    CLASS(Solver) :: this
    
    !     Calculate the energy in the chosen
    !     approximation
    CALL this%energy%calculateEnergy()
    
    
  END SUBROUTINE calcEnergy
  
  FUNCTION readAndSetup() RESULT(mySolver)
    USE Pt2CorrEnergyMod
    USE CcdCorrEnergyMod
    USE Yukawa2dHamMod
    USE Yukawa3dHamMod
    USE SpBasisPW2dMod
    USE SpBasisPW3dMod
    TYPE(Solver) :: mySolver
    CHARACTER(LEN=100) :: interaction
    CHARACTER(LEN=100) :: approximation
    CHARACTER(LEN=100) :: spBasisOption
    CHARACTER(LEN=100) :: spinState
    CLASS(Hamiltonian), POINTER :: myHamiltonian => NULL()
    CLASS(Yukawa2dHam), POINTER :: myYukawa2dHam => NULL()
    CLASS(Yukawa3dHam), POINTER :: myYukawa3dHam => NULL()
    TYPE(RefEnergy) :: myRefEnergy
    CLASS(CorrEnergy), POINTER :: myCorrEnergy => NULL()
    TYPE(Pt2CorrEnergy), POINTER :: myPt2CorrEnergy => NULL()
    TYPE(CcdCorrEnergy), POINTER :: myCcdCorrEnergy => NULL()
    TYPE(FockMatrix) :: myFockMatrix
    CLASS(SpBasis), POINTER :: mySpBasis => NULL()
    CLASS(SpBasisPW2d), POINTER :: mySpBasisPW2d => NULL()
    CLASS(SpBasisPW3d), POINTER :: mySpBasisPW3d => NULL()
    TYPE(BlockList) :: myBlockList
    TYPE(BlockListRe) :: myBlockListRe
    TYPE(Energy) :: myEnergy 
    REAL(DP) :: tolerance, mu, rs, omega, area, volume, alpha
    INTEGER :: nOccupied, nShells, maxIter, nQNumbers
    
    !     Read the input file
    OPEN(5,FILE='input.dat')
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
    
!    WRITE(*,*) ' '
    !     Create single-particle basis and
    !     single-particle pair objects   
    IF (spBasisOption == 'pw2d') THEN
       !     Plane waves in two dimensions
       ALLOCATE(mySpBasis, SOURCE=mySpBasisPW2d)
       
       area = nOccupied*pi*rs**2
       mySpBasis => SpBasisPW2d_(nOccupied, nShells, spinState, area)
       
       nQNumbers = 3
       
    ELSEIF (spBasisOption == 'pw3d') THEN
       !     Plane waves in three dimensions
       ALLOCATE(mySpBasis, SOURCE=mySpBasisPW3d)
       
       volume = nOccupied*4.d0*pi*rs**3/3.d0
    !   write(*,*) 'volume=',volume
       mySpBasis => SpBasisPW3d_(nOccupied, nShells, spinState, volume)
       
       nQNumbers = 4
       
    ENDIF
    
 !   WRITE(*,*) ' Setting up single-particle basis...'
    !     Set up single-particle basis
    CALL mySpBasis%setupSpBasis()
    !     Set up list of (K_CM, M_S) blocks
    myBlockList = BlockList_(nQNumbers, mySpBasis)
    CALL myBlockList%setupBlockList()
    myBlockListRe = BlockListRe_(nQNumbers, mySpBasis, myBlockList)
    CALL myBlockListRe%setupPhBlockList()
    
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

!    WRITE(*,*) ' Setting up interaction matrix...'
    !     Set up interaction matrix in the given
    !     single-particle basis
    CALL setupInteractionMatr(myHamiltonian)
    CALL setupInteractionMatrRe(myHamiltonian)
    
!    WRITE(*,*) ' Setting up Fock matrix...'
    !     Create a reference energy object
    myRefEnergy = RefEnergy_(myHamiltonian, mySpBasis, &
         myBlockList)
    !     Create a single-particle potential object
    myFockMatrix = FockMatrix_(myHamiltonian, mySpBasis)
    !     Set up single-particle potential vector
    CALL setupFockVector(myFockMatrix)
    
    IF (approximation == 'pt2') THEN
       ALLOCATE(myCorrEnergy, SOURCE=myPt2CorrEnergy)
       
       myCorrEnergy => Pt2CorrEnergy_(myHamiltonian,&
            myFockMatrix, mySpBasis)
       
    ELSEIF (approximation == 'ccd') THEN
       ALLOCATE(myCorrEnergy, SOURCE=myCcdCorrEnergy)
       
       myCorrEnergy => CcdCorrEnergy_(myHamiltonian,&
            myFockMatrix, mySpBasis, maxIter, &
            tolerance, alpha)
    ENDIF
    
    !     Create an energy object
    myEnergy = Energy_(myHamiltonian, approximation,&
         myRefEnergy, myCorrEnergy)
    !     Create a solver object
    mySolver = Solver(myEnergy)
    
    !IF (ASSOCIATED(myCorrEnergy)) DEALLOCATE(myCorrEnergy)
    
 !   WRITE(*,*) ' '
 !   WRITE(*,*) '------  Setup is done  ------'
 !   WRITE(*,*) ' '
    
    
  END FUNCTION readAndSetup


  
END MODULE SolverMod
