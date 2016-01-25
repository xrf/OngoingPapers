
PROGRAM testMem
  IMPLICIT NONE

  TYPE type1
     INTEGER :: dim1, dim2
     REAL, ALLOCATABLE :: matrix1(:,:)
  END TYPE type1

  TYPE(type1) :: t

  t%dim1 = 10000
  t%dim2 = 20000
  ALLOCATE(t%matrix1(t%dim1, t%dim2))
  t%matrix1 = 0.0

  CALL wait(5)
  WRITE(*,*) '--'

  CALL routine(t, 5)

  WRITE(*,*) '..'
  CALL wait(5)
  
  DEALLOCATE(t%matrix1)
  
CONTAINS

  
  SUBROUTINE routine(data, time)
    TYPE(type1), INTENT(INOUT) :: data
    INTEGER, INTENT(IN) :: time
    
    CALL wait(time)

  END SUBROUTINE routine


END PROGRAM testMem



!
!     Wait (sleep) in 'interval' seconds.
!
SUBROUTINE wait(interval)
  !USE ConstantsMod
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: interval
  INTEGER :: time(8), second1, second2
  
  CALL DATE_AND_TIME(values=time)
  second1 = 3600*time(5) + 60*time(6) + time(7)
  
  DO
     CALL DATE_AND_TIME(values=time)
     second2 = 3600*time(5) + 60*time(6) + time(7)
     
     IF (second2-second1 >= interval) EXIT
  ENDDO

  
END SUBROUTINE wait
