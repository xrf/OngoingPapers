
!
!     Class for (K_CM, M_S) blocks.
!
MODULE Block1Mod
  USE IntMatrixMod
  IMPLICIT NONE
  
  TYPE :: Block1
     TYPE(IntMatrix) :: hhPairs 
     TYPE(IntMatrix) :: ppPairs 
     TYPE(IntMatrix) :: phPairs
     TYPE(IntMatrix) :: hpPairs
     LOGICAL :: occupiedHh, occupiedPp
     LOGICAL :: occupiedPh, occupiedHp
     INTEGER :: nQNumbers, rank
     INTEGER, ALLOCATABLE :: blockQNrs(:)
  END TYPE Block1
  
CONTAINS
  
  !
  !     Constructor for Block1 objects.
  !
  FUNCTION Block1_(nQNumbers) RESULT(this)
    INTEGER, INTENT(IN) :: nQNumbers
    TYPE(Block1) :: this
    
    this%nQNumbers = nQNumbers
    ALLOCATE(this%blockQNrs(nQNumbers))
    this%blockQNrs = 0
    this%occupiedHh = .FALSE.
    this%occupiedPp = .FALSE.
    this%occupiedPh = .FALSE.
    this%occupiedHp = .FALSE.
     
    CALL add2MemoryInt(nQnumbers)

    
  END FUNCTION Block1_
  
  !
  !     Destrucor for Block1 objects.
  !
  SUBROUTINE Block1_d(this)
    TYPE(Block1), INTENT(INOUT) :: this
    INTEGER :: elements
    
    CALL IntMatrix_d(this%hhPairs)
    CALL IntMatrix_d(this%ppPairs)
    CALL IntMatrix_d(this%hpPairs)
    CALL IntMatrix_d(this%phPairs)
    
    elements = 0 
    IF (ALLOCATED(this%blockQNrs)) THEN
       DEALLOCATE(this%blockQNrs)
       elements = elements - this%nQNumbers
    ENDIF
    CALL add2MemoryInt(elements)
    
    
  END SUBROUTINE Block1_d
  
  
END MODULE Block1Mod
