
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
     TYPE(IntMatrix) :: hh2Index
     TYPE(IntMatrix) :: pp2Index
     TYPE(IntMatrix) :: ph2Index
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

    
  END FUNCTION Block1_
  
  
END MODULE Block1Mod
