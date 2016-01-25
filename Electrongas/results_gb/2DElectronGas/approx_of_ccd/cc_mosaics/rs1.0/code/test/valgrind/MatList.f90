
!
!     Abstract class for a general matrix list
!
MODULE MatListMod
  USE ConstantsMod
  USE MatrixMod
  IMPLICIT NONE

  TYPE, ABSTRACT :: MatList
     INTEGER :: nItems
  END TYPE MatList
  

  
END MODULE MatListMod
