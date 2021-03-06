
!
!     Class for a general real matrix.
!
MODULE RealMatrixMod
  USE MemoryCounterMod
  USE ConstantsMod
  USE MatrixMod
  IMPLICIT NONE
  
  TYPE, EXTENDS(Matrix) :: RealMatrix
     REAL(DP), ALLOCATABLE :: matr(:,:)
  END TYPE RealMatrix
  
CONTAINS
  
  !
  !     Constructor for RealMatrix object
  !
  FUNCTION RealMatrix_(dim1, dim2) RESULT(this)
    INTEGER, INTENT(IN) :: dim1, dim2
    TYPE(RealMatrix) :: this
    
    this%dim1 = dim1
    this%dim2 = dim2
    ALLOCATE(this%matr(dim1, dim2))
    this%matr = 0.d0
    
    CALL add2MemoryR8(dim1*dim2)
    
    
  END FUNCTION RealMatrix_
  
  !
  !     Destructor for a RealMatrix object
  !
  SUBROUTINE RealMatrix_d(this)
    TYPE(RealMatrix) :: this
    INTEGER :: dim1
    
    IF (ALLOCATED(this%matr)) THEN
       DEALLOCATE(this%matr)
       dim1 = MAX(this%dim1, 1)
       CALL add2MemoryR8(-dim1*this%dim2)
    ENDIF
    
    
  END SUBROUTINE RealMatrix_d
  
  
END MODULE RealMatrixMod
