
!
!     Class for the single-particle basis.
!
MODULE SpBasisMod
  USE ConstantsMod
  USE MpiMod
  IMPLICIT NONE
  
  TYPE, ABSTRACT :: SpBasis
     !
     !     nOccupied:         Number of occupied orbitals
     !     nUnoccupied:       Number of unoccupied orbitals
     !     nShells:           Total number of shells
     !     spQuantNr(:,:):    Single-particle quantum numbers
     !     spEnergies(:):     Sinlge-particle energies
     !                        (diagonal matrix elements of h0)
     !
     INTEGER :: nOccupied               
     INTEGER :: nUnoccupied            
     INTEGER :: nShells                 
     INTEGER :: nMax
     INTEGER, ALLOCATABLE :: spQuantNr(:,:) 
     REAL(DP), ALLOCATABLE :: spEnergies(:)
     TYPE(Mpi), POINTER :: mpiJobs
   CONTAINS
     PROCEDURE(setupSp), &
          DEFERRED :: setupSpBasis
  END TYPE SpBasis
  
  ABSTRACT INTERFACE
     SUBROUTINE setupSp(this)
       IMPORT SpBasis
       CLASS(SpBasis), TARGET, INTENT(INOUT) :: this
     END SUBROUTINE setupSp
  END INTERFACE
  
CONTAINS
  
  !
  !     Destructor for polymorphic SpBasis object
  !
  SUBROUTINE SpBasis_d(this)
    CLASS(SpBasis), TARGET, INTENT(INOUT) :: this
    
    IF (ALLOCATED(this%spQuantNr)) DEALLOCATE(this%spQuantNr)
    IF (ALLOCATED(this%spEnergies)) DEALLOCATE(this%spEnergies)

    !IF (ASSOCIATED(this%mpiJobs)) THEN
    !   CALL 
    !ENDIF
    
  END SUBROUTINE SpBasis_d
  
  
END MODULE SpBasisMod
