
!
!     Class for the Fock matrix
!
MODULE FockMatrixMod
  USE HamiltonianMod
  USE ConstantsMod
  USE SpBasisMod
  IMPLICIT NONE
  
  TYPE FockMatrix
     !
     !     hamilt:       Hamiltonian operator
     !     basis:        Single-particle basis
     !     matrix(:,:):  Fock matrix 
     !     vector(:):    Diagonal Fock matrix      
     !
     CLASS(Hamiltonian), POINTER :: hamilt  
     CLASS(SpBasis), POINTER :: basis       
     REAL(DP), ALLOCATABLE :: matrix(:,:)     
     REAL(DP), ALLOCATABLE :: vector(:)
   CONTAINS
     PROCEDURE :: calcFockPotential
     PROCEDURE :: calcFockMatEleD
     PROCEDURE :: setupFockVector
  END TYPE FockMatrix
  
  
CONTAINS
  
  !
  !     Constructor for FockMatrix object.     
  !
  FUNCTION FockMatrix_(ham, basis) RESULT(this)
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    TYPE(FockMatrix), POINTER :: this
    INTEGER :: nOrbits
    
    ALLOCATE(this)
    
    nOrbits = basis%nOccupied + basis%nUnoccupied
    !    ALLOCATE(fockMat%matrix(nOrbits, nOrbits))
    ALLOCATE(this%vector(nOrbits))
    !    this%matrix = 0.d0
    this%vector = 0.d0
    this%hamilt => ham
    this%basis => basis
    
    
  END FUNCTION FockMatrix_

  !
  !     Destructor for FockMatrix object.
  !
  SUBROUTINE FockMatrix_d(this)
    TYPE(FockMatrix), TARGET, INTENT(INOUT) :: this
    
    IF (ALLOCATED(this%vector)) THEN
       DEALLOCATE(this%vector)
    ENDIF

    IF (ASSOCIATED(this%hamilt)) THEN
       CALL Hamiltonian_d(this%hamilt)
       this%hamilt => NULL()
    ENDIF

    IF (ASSOCIATED(this%basis)) THEN
       CALL SpBasis_d(this%basis)
       this%basis => NULL()
    ENDIF
    
    
  END SUBROUTINE FockMatrix_d
  
  !
  !      Routine to calculate the single-particle
  !      potential u(p, q) = \sum_i <p,i|v|q,i>_AS. 
  !
  FUNCTION calcFockPotential(this, orbit1, orbit2) &
       RESULT(potential)
    INTEGER, INTENT(IN) :: orbit1, orbit2
    CLASS(FockMatrix), TARGET, INTENT(IN) :: this
    REAL(DP) :: potential
    REAL(DP) :: temp, interaction
    INTEGER :: hole
    
    temp = 0.d0
    DO hole=1, this%basis%nOccupied
       
       !     Get antisymmetrized two-particle
       !     interaction
       interaction = this%hamilt%getInteraction( &
            orbit1, hole, orbit2, hole)
       
       temp = temp + interaction
    ENDDO
    potential = temp
    
    
  END FUNCTION calcFockPotential
    
  !
  !     Calculate a diagonal Fock matrix element
  !     <p|h0|p> + \sum_i <p,i|v|p,i>_AS. 
  !
  FUNCTION calcFockMatEleD(this, orbit) RESULT(matrixElement)
    INTEGER, INTENT(IN) :: orbit
    CLASS(FockMatrix), TARGET, INTENT(IN) :: this
    REAL(DP) :: nonInteract, potential
    REAL(DP) :: matrixElement
    
    nonInteract = this%basis%spEnergies(orbit)
    potential = this%calcFockPotential(orbit, orbit)
    
    matrixElement = nonInteract !+ potential
    
    
  END FUNCTION calcFockMatEleD
  

  !
  !     Set up vector of single-particle
  !     energies (diagonal elements of Fock
  !     matrix).
  !
  SUBROUTINE setupFockVector(this)
    CLASS(FockMatrix), TARGET, INTENT(INOUT) :: this
    INTEGER :: nOrbits, orbit
    REAL(DP) :: matrixEle
    
    nOrbits = this%basis%nOccupied &
         + this%basis%nUnoccupied
    
    DO orbit=1, nOrbits
       matrixEle = this%calcFockMatEleD(orbit)
       this%vector(orbit) = matrixEle
    ENDDO
    
    
  END SUBROUTINE setupFockVector
  
  
END MODULE FockMatrixMod
