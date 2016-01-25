
!
!     Class for a general integer matrix.
!
MODULE IntMatrixMod
  USE MemoryCounterMod
  USE ConstantsMod
  USE MatrixMod
  IMPLICIT NONE
  
  TYPE, EXTENDS(Matrix) :: IntMatrix
     INTEGER, ALLOCATABLE :: matr(:,:)
  END TYPE IntMatrix
  
CONTAINS
  
  !
  !     Constructor for IntMatrix object
  !
  FUNCTION IntMatrix_(dim1, dim2) RESULT(this)
    INTEGER, INTENT(IN) :: dim1, dim2
    TYPE(IntMatrix) :: this 
    
    this%dim1 = dim1
    this%dim2 = dim2
    ALLOCATE(this%matr(dim1, dim2))
    this%matr = 0
    
    CALL add2MemoryInt(dim1*dim2)
    
    
  END FUNCTION IntMatrix_
  
  !
  !     Destructor for an IntMatrix object
  !
  SUBROUTINE IntMatrix_d(this)
    TYPE(IntMatrix) :: this
    INTEGER :: dim1

    IF (ALLOCATED(this%matr)) THEN
       DEALLOCATE(this%matr)
       dim1 = MAX(1, this%dim1)
       CALL add2MemoryInt(-dim1*this%dim2)
    ENDIF
    
    
  END SUBROUTINE IntMatrix_d
  
  
END MODULE IntMatrixMod
