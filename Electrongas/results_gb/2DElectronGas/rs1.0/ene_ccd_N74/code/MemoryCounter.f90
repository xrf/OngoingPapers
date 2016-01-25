
!
!     Class for memory counter
!
MODULE MemoryCounterMod
  USE ConstantsMod
  USE MpiMod
  IMPLICIT NONE
  INTEGER(KIND=8), PUBLIC :: memoryUsed, memoryMax
  
  !TYPE :: MemoryCounter
  !   INTEGER(KIND=8), PUBLIC :: memoryUsed
  !CONTAINS
  !  PROCEDURE :: add2MemoryR8
  !  PROCEDURE :: add2MemoryInt
  !  PROCEDURE :: printMemoryUse
  !END TYPE MemoryCounter
  
CONTAINS
  
  !
  !     Constructor for MemoryCounter objects
  !
  ! FUNCTION MemoryCounter_() RESULT(this)
  !   TYPE(MemoryCounter) :: this
  
  !   memoryUsed = 0
  
  
  ! END FUNCTION MemoryCounter_
  
  !
  !     Add items of type REAL(8) to the memory 
  !     counter.
  ! 
  SUBROUTINE add2MemoryR8(items)
    INTEGER, INTENT(IN) :: items
    
    memoryUsed = memoryUsed + items*sizeReal8
!    write(*,*) 'add:',items*sizeReal8,',memory:',memoryUsed
    
    IF (memoryUsed > memoryMax) memoryMax = memoryUsed
    
    
  END SUBROUTINE add2MemoryR8
  
  !
  !     Add items of type INTEGER to the memory
  !     counter.
  !
  SUBROUTINE add2MemoryInt(items)
    INTEGER, INTENT(IN) :: items
    
    memoryUsed = memoryUsed + items*sizeInt
!    write(*,*) 'add:',items*sizeInt,',memory:',memoryUsed

    IF (memoryUsed > memoryMax) memoryMax = memoryUsed
    
    
  END SUBROUTINE add2MemoryInt
  
  !
  !     Print how much memory has been used
  !
  SUBROUTINE printMemoryUse(mpiJobs)
    TYPE(Mpi), OPTIONAL, TARGET, INTENT(IN) :: mpiJobs
    
    IF (PRESENT(mpiJobs)) THEN
       IF (mpiJobs%iAm == 0) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'Memory in use: ', memoryUsed, ' Bytes'
          WRITE(*,*) ' '
       ENDIF
    ELSE
       WRITE(*,*) ' '
       WRITE(*,*) 'Memory in use: ', memoryUsed, ' Bytes'
       WRITE(*,*) ' '
    ENDIF
    
    
  END SUBROUTINE printMemoryUse
  
  !
  !     Print maximum memory used
  !
  SUBROUTINE printMaxMemoryUsed(mpiJobs)
    TYPE(Mpi), TARGET, INTENT(IN) :: mpiJobs
    
    IF (mpiJobs%iAm == 0) THEN
       WRITE(*,*) ' '
       WRITE(*,*) 'Maximum memory used: ', memoryMax, ' Bytes'
       WRITE(*,*) ' '
    ENDIF
    
    
  END SUBROUTINE printMaxMemoryUsed
  
  
END MODULE MemoryCounterMod
