

PROGRAM testMKL
  IMPLICIT NONE
  INTEGER :: n
  REAL(8) :: time1, time2, time3
  INCLUDE 'mpif.h'
  
  time1 = MPI_WTIME()
  
  n = 4000
  CALL multiply(n)
  
  time2 = MPI_WTIME()
  
  OPEN(6, FILE='out')
  WRITE(6,*) 'Time: ', time2-time1
  WRITE(*,*) 'Time: ', time2-time1
  CLOSE(6)
  
  
END PROGRAM testMKL

SUBROUTINE multiply(n)
  INTEGER, INTENT(IN) :: n
  REAL(8), ALLOCATABLE :: a(:,:), b(:,:), c(:,:)
  
  ALLOCATE(a(n,n), b(n,n), c(n,n))
  CALL RANDOM_NUMBER(a)
  CALL RANDOM_NUMBER(b)
  
  CALL DGEMM('N', 'N', n, n, n, 1.d0, a, n, b, n, 0.d0, c, n)
  
  DEALLOCATE(a, b, c)
  
  
END SUBROUTINE multiply
