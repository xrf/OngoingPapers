
!
!     Abstract class for a general matrix.
!
MODULE MatrixMod
  USE ConstantsMod
  IMPLICIT NONE
  
  TYPE, ABSTRACT :: Matrix
     INTEGER :: dim1, dim2
  END TYPE Matrix
  
  
END MODULE MatrixMod
