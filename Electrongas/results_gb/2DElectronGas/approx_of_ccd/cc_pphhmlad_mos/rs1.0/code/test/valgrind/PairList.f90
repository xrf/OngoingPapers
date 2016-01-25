
!
!     Class for list of pairs of single-particle states.
!
MODULE PairListMod
  IMPLICIT NONE
  
  TYPE :: PairList
     INTEGER :: nPairs, nQNumbers
     INTEGER, ALLOCATABLE :: blockQNumbers(:,:)
     INTEGER, ALLOCATABLE :: pairs(:,:)
     INTEGER, ALLOCATABLE :: rank(:)
     INTEGER, ALLOCATABLE :: blockInd(:)
  END TYPE PairList
  
CONTAINS
  
  !
  !     Constructor for PairList objects.
  !
  FUNCTION PairList_(nPairs, nQNumbers) RESULT(this)
    INTEGER, INTENT(IN) :: nPairs, nQNumbers
    TYPE(PairList) :: this
    
    this%nPairs = nPairs
    this%nQNumbers = nQNumbers
    ALLOCATE(this%blockQNumbers(nPairs, nQNumbers))
    ALLOCATE(this%pairs(nPairs, 2))
    ALLOCATE(this%rank(nPairs))
    ALLOCATE(this%blockInd(nPairs))
    this%blockQNumbers = 0
    this%pairs = 0
    this%rank = 0
    this%blockInd = 0
    
    
  END FUNCTION PairList_
  
  
END MODULE PairListMod
