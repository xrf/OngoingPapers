
!
!     Class for list of pairs of single-particle states.
!
MODULE PairListMod
  USE MemoryCounterMod
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
    INTEGER :: elements
    
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
    
    elements = nPairs*nQNumbers + 4*nPairs 
    CALL add2MemoryInt(elements)
    
    
  END FUNCTION PairList_
  
  !
  !     Destructor for PairList objects
  !
  SUBROUTINE PairList_d(this)
    TYPE(PairList), INTENT(INOUT) :: this
    INTEGER :: elements
    
    elements = 0
    IF (ALLOCATED(this%blockQNumbers)) THEN
       DEALLOCATE(this%blockQNumbers)
       elements = elements - this%nPairs*this%nQNumbers   
    ENDIF
    IF (ALLOCATED(this%pairs)) THEN
       DEALLOCATE(this%pairs)
       elements = elements - 2*this%nPairs
    ENDIF
    IF (ALLOCATED(this%rank)) THEN
       DEALLOCATE(this%rank)
       elements = elements - this%nPairs
    ENDIF
    IF (ALLOCATED(this%blockInd)) THEN
       DEALLOCATE(this%blockInd)
       elements = elements - this%nPairs
    ENDIF
    CALL add2MemoryInt(elements)
    
    
  END SUBROUTINE PairList_d
  
  
END MODULE PairListMod
