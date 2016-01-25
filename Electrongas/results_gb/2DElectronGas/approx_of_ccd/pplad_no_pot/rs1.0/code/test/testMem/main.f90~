
PROGRAM tests
  USE IFPORT
  IMPLICIT NONE
  INTEGER, ALLOCATABLE :: array(:), array2(:)
  INTEGER(2), EXTERNAL :: compare
  INTEGER :: i
  INTEGER(SIZEOF_SIZE_T) :: array_size, int_size
  
  array_size = 5
  int_size = 4
  ALLOCATE(array(array_size), array2(array_size))
  
  array(1) = 2223244
  array(2) = 33435
  array(3) = -12
  array(4) = 2
  array(5) = 120000
  array2 = (/(i, i=1, array_size)/)
  
  DO i=1, array_size
     WRITE(*,*) array(i), array2(i) 
  ENDDO
  
  !CALL qsort(array, array_size, int_size, compare)
  !CALL quicksort(array, array_size)
  CALL quicksort2(array, array2, array_size)
  
  WRITE(*,*) '-------------------------'
  DO i=1, array_size
     WRITE(*,*) array(i), array2(i) 
  ENDDO
  
  !WRITE(*,*) '==========================='
  
  DEALLOCATE(array, array2)
  
  
END PROGRAM tests

INTEGER(2) FUNCTION compare(i1, i2) 
  INTEGER, INTENT(IN) :: i1, i2
  
  compare = INT(i1 - i2, 2)
  write(*,*) 'i1=',i1,',i2=',i2,',c=',compare
  
  
END FUNCTION compare

!
!     Quicksort routine taken from W. H. Press et al.,
!     Numerical Recipes in Fortran 90, 1996.
!
SUBROUTINE quicksort(array, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n 
  INTEGER, INTENT(INOUT) :: array(n)
  INTEGER, PARAMETER :: NN=15, NSTACK=50
  !
  !     Sorts an array 'array' into ascending numerical
  !     order using the Quicksort algorithm. 'array' is
  !     replaced on output by its sorted rearrangement. 
  !
  !     Parameters: NN is the size of subarrays sorted
  !     by straight insertion and NSTACK is the required
  !     auxiliary storage.
  !
  INTEGER :: a
  INTEGER :: k, i, j, jstack, l, r
  INTEGER, DIMENSION(NSTACK) :: istack

  jstack = 0
  l = 1
  r = n
  DO
     
     !     Insertion sort when subarray small enough.
     IF (r-1 < NN) THEN
        DO j=l+1, r
           a = array(j)
           DO i=j-1, l, -1
              IF (array(i) <= a) EXIT
              array(i+1) = array(i)
           ENDDO
           array(i+1) = a
        ENDDO
        IF (jstack == 0) RETURN
        !     Pop stack and begin a new round of 
        !     partitioning.
        r = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack - 2
     ELSE
        !     Choose median of left, center, and 
        !     right elements as partitioning element
        !     a. Also rearrange so that a(l) <= a(l+1)
        !     <= a(r)
        k = (l+r)/2
        CALL swap(array(k), array(l+1))
        CALL swap(array(l), array(r), array(l)>array(r))
        CALL swap(array(l+1), array(r), array(l+1)>array(r))
        CALL swap(array(l), array(l+1), array(l)>array(l+1))
        !     Initialize pointers for partitioning.
        i = l + 1
        j = r
        !     Partitioning element.
        a = array(l+1)
        DO
           !     Here is the meat.
           !     Scan up to find element >= a.
           DO
              i = i+1
              IF (array(i) >= a) EXIT
           ENDDO
           !     Scan down to find element <= a.
           DO
              j = j - 1
              IF (array(j) <= a) EXIT
           ENDDO
           !     Pointers crossed. Exit with 
           !     partitioning complete. 
           IF (j < i) EXIT
           !     Exchange elements.
           CALL swap(array(i), array(j))
        ENDDO
        !     Insert partitioning element.
        array(l+1) = array(j)
        array(j) = a
        jstack = jstack + 2
        !     Push pointers to larger subarray on
        !     stack; process smaller subarray
        !     immediately.
        IF (jstack > NSTACK) THEN
           WRITE(*,*) 'nrerror: sort: NSTACK too small'
           STOP 'program terminated by nrerror'
        ENDIF
        IF (r-i+1 >= j-1) THEN
           istack(jstack) = r
           istack(jstack-1) = i
           r = j - 1
        ELSE
           istack(jstack) = j - 1
           istack(jstack-1) = l
           l = i
        ENDIF
     ENDIF
  ENDDO
  
  
END SUBROUTINE quicksort


!
!     Quicksort routine from W. H. Press et al.,
!     Numerical Recipes in Fortran 90, 1996.
!
SUBROUTINE quicksort2(array, array2, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n 
  INTEGER, INTENT(INOUT) :: array(n), array2(n)
  INTEGER, PARAMETER :: NN=15, NSTACK=50
  !
  !     Sorts an array 'array' into ascending numerical
  !     order using the Quicksort algorithm. 'array' is
  !     replaced on output by its sorted rearrangement. 
  !
  !     Parameters: NN is the size of subarrays sorted
  !     by straight insertion and NSTACK is the required
  !     auxiliary storage.
  !
  INTEGER :: a, b
  INTEGER :: k, i, j, jstack, l, r
  INTEGER, DIMENSION(NSTACK) :: istack
  
  jstack = 0
  l = 1
  r = n
  DO
     
     !     Insertion sort when subarray small enough.
     IF (r-1 < NN) THEN
        DO j=l+1, r
           a = array(j)
           b = array2(j)
           DO i=j-1, l, -1
              IF (array(i) <= a) EXIT
              array(i+1) = array(i)
              array2(i+1) = array2(i)
           ENDDO
           array(i+1) = a
           array2(i+1) = b
        ENDDO
        IF (jstack == 0) RETURN
        !     Pop stack and begin a new round of 
        !     partitioning.
        r = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack - 2
     ELSE
        !     Choose median of left, center, and 
        !     right elements as partitioning element
        !     a. Also rearrange so that a(l) <= a(l+1)
        !     <= a(r)
        k = (l+r)/2
        CALL swap(array(k), array(l+1))
        CALL swap(array2(k), array2(l+1))
        CALL swap(array(l), array(r), array(l)>array(r))
        CALL swap(array2(l), array2(r), array2(l)>array2(r))
        CALL swap(array(l+1), array(r), array(l+1)>array(r))
        CALL swap(array2(l+1), array2(r), array2(l+1)>array2(r))
        CALL swap(array(l), array(l+1), array(l)>array(l+1))
        CALL swap(array2(l), array2(l+1), array2(l)>array2(l+1))
        !     Initialize pointers for partitioning.
        i = l + 1
        j = r
        !     Partitioning element.
        a = array(l+1)
        b = array2(l+1)
        DO
           !     Here is the meat.
           !     Scan up to find element >= a.
           DO
              i = i+1
              IF (array(i) >= a) EXIT
           ENDDO
           !     Scan down to find element <= a.
           DO
              j = j - 1
              IF (array(j) <= a) EXIT
           ENDDO
           !     Pointers crossed. Exit with 
           !     partitioning complete. 
           IF (j < i) EXIT
           !     Exchange elements.
           CALL swap(array(i), array(j))
           CALL swap(array2(i), array2(j))
        ENDDO
        !     Insert partitioning element.
        array(l+1) = array(j)
        array(j) = a
        array2(l+1) = array2(j)
        array2(j) = b
        jstack = jstack + 2
        !     Push pointers to larger subarray on
        !     stack; process smaller subarray
        !     immediately.
        IF (jstack > NSTACK) THEN
           WRITE(*,*) 'nrerror: sort: NSTACK too small'
           STOP 'program terminated by nrerror'
        ENDIF
        IF (r-i+1 >= j-1) THEN
           istack(jstack) = r
           istack(jstack-1) = i
           r = j - 1
        ELSE
           istack(jstack) = j - 1
           istack(jstack-1) = l
           l = i
        ENDIF
     ENDIF
  ENDDO
  
  
END SUBROUTINE quicksort2


SUBROUTINE swap(a, b)
  !     Swap the contents of a and b.
  INTEGER, INTENT(INOUT) :: a, b
  INTEGER :: temp
  
  temp = a
  a = b
  b = temp
  
  
END SUBROUTINE swap


! PROGRAM SORTQ
!   USE IFPORT
!   integer(2), external :: cmp_function
!   integer(2) insort(26), i
!   integer (SIZEOF_SIZE_T) array_len, array_size

!   array_len = 26
!   array_size = 2
  
!   do i=90,65,-1
!      insort(i-64)=91 - i
!   end do
  
!   print *, "Before: "
  
!   print *,insort
  
!   CALL qsort(insort,array_len,array_size,cmp_function)
  
!   print *, 'After: '
!   print *, insort

  
! END PROGRAM SORTQ

! !
! integer(2) function cmp_function(a1, a2) 
!   integer(2) a1, a2
  
!   cmp_function=a1-a2
  
! end function cmp_function
