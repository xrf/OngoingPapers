
!
!     Class for a general integer matrix.
!
MODULE IntMatrixMod
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
    
    
  END FUNCTION IntMatrix_

  !
  !     Constructor for IntMatrix object
  !
  FUNCTION IntMatrix2_(dim1, dim2) RESULT(this)
    INTEGER, INTENT(IN) :: dim1, dim2
    TYPE(IntMatrix) :: this 
    
    this%dim1 = dim1
    this%dim2 = dim2
    ALLOCATE(this%matr(dim1, dim2))
    !this%matr = 0
    
    
  END FUNCTION IntMatrix2_
  
  !
  !     Destructor for an IntMatrix object
  !
  SUBROUTINE IntMatrix_d(this)
    TYPE(IntMatrix) :: this
    
    DEALLOCATE(this%matr)
    
    
  END SUBROUTINE IntMatrix_d
  
  
END MODULE IntMatrixMod
