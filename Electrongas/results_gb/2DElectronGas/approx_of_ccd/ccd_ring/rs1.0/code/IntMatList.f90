
!
!     Class for an integer matrix list.
!
MODULE IntMatListMod
  USE ConstantsMod
  USE IntMatrixMod
  IMPLICIT NONE
  
  TYPE :: IntMatList
     INTEGER :: nItems
     TYPE(IntMatrix), ALLOCATABLE :: list(:)
  END TYPE IntMatList
  
CONTAINS
  
  !
  !     Constructor for IntMatList object.
  !
  FUNCTION IntMatList_(nItems) RESULT(this)
    INTEGER, INTENT(IN) :: nItems
    TYPE(IntMatList) :: this
    
    this%nItems = nItems
    ALLOCATE(this%list(nItems))
    
    
  END FUNCTION IntMatList_
  
  !
  !     Destructor for an IntMatList object
  !
  SUBROUTINE IntMatList_d(this)
    TYPE(IntMatList), INTENT(INOUT) :: this
    INTEGER :: item
    
    IF (ALLOCATED(this%list)) THEN
       DO item=1, this%nItems
          CALL IntMatrix_d(this%list(item))
       ENDDO
       DEALLOCATE(this%list)
    ENDIF
    
    
  END SUBROUTINE IntMatList_d

END MODULE IntMatListMod
