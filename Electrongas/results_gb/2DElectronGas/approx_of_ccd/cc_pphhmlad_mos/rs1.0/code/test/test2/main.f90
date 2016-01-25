! program showing how to call 'QSORT' on
! a user-defined type.
!
! Define the type to be shared.
!
module share_type
  type element_type
     integer       :: data
     character(10) :: key
  end type element_type
end module share_type

! Main program calls QSORT.
!
program main
  use IFPORT       ! To get QSORT
  use share_type   ! To get shared type
    
  ! Define an overload of the default QSORT signature
  ! with a signature using the shared type.
  !
  interface QSORT      
     subroutine QSORT_element_types(array, len, isize, comp)
       use share_type
       type(element_type) array(len)
       integer len, isize
       integer(2), external :: comp

       !
       ! Hook the overload to the real thing but be careful
       ! to connect to the correct qsort: the Fortran one, not
       ! the C one!
       !
       ! We need to call the _Fortran_ qsort, not the _C_ one, or 
       ! there will be errors from the 1-origin vs. 0-origin indexing
       ! and the row-major vs. column-major ordering.
       !
       ! The symptom is that "OrderCharCI" is called with pointer values
       ! which are outside the bounds of the array to be sorted.
       !
       !DEC$ATTRIBUTES ALIAS:'qsort_' :: QSORT_element_types
       
     end subroutine QSORT_element_types
  end interface QSORT
  
  type(element_type) :: c(7)
  integer(2), external :: OrderCharCI
  integer :: size_of_element, size_of_array
  
  ! Fill in the array to be sorted.  The data value is chosen so
  ! that the sorted array will have the values in numeric order.
  ! Thus we can check the result of the sort.
  !
  
  c(1)%key  = 'aisjdop'
  c(1)%data = 3
  c(2)%key  = '35djf2'
  c(2)%data = 1
  c(3)%key  = 'ss:ss'
  c(3)%data = 6
  c(4)%key  = 'MMhQQ'
  c(4)%data = 4
  c(5)%key  = 'mmHqq'
  c(5)%data = 5
  c(6)%key  = 'aaaa'
  c(6)%data = 2
  c(7)%key  = '["\/'
  c(7)%data = 7
  
  size_of_array   = size(c)         !  7
  size_of_element = sizeof(c(1))    ! 16
  
  write(*,*) '"C" is:'
  do i = 1, 7
     write(*,*) ' "', c(i)%key, '" value ', c(i)%data
  end do
  write(*,*) ' '
  write(*,*) 'size of C is            ', size_of_array, ' elements'
  write(*,*) 'size of element C(1) is ', size_of_element, ' bytes'
  write(*,*) 'len of key in C(1) is   ',   len(c(1)%key)
  write(*,*) ' '

  write(*,*) '###########################################'
  ! Call the overloaded QSORT routine.
  !
  !Call QSort_element_types(C, size_of_array, size_of_element, OrderCharCI)
  Call QSort(C, size_of_array, size_of_element, OrderCharCI)
  
  write(*,*) 'Sorted "C" is '
  do i = 1, 7
     write(*,*) ' "', c(i)%key, '" value ', c(i)%data
  end do

end program

! Computes order of character strings using a case insensitive ordering.
!
! Return -1 if C1 before C2, 0 if C1 = C2, and 1 if C1 after C2.
!
! Called first with the pair (2,3), then (1,2), then (1,3)...when passing
! character strings of length 10.
!
! Passing "element_type" objects, it's called first with the pair (1, <invalid>),
! and the second item has a address well before the beginning of "C".
!
function OrderCharCI(c1, c2)
  use share_type
  implicit none
  type(element_type), intent(in) :: c1 ! Character strings to be ordered.
  type(element_type), intent(in) :: c2 !

  ! Function result:
  !
  integer(2) :: OrderCharCI

  ! Locals:
  !
  character(10) :: c1L !} Local copies of c1 and c2.
  character(10) :: c2L !}
  integer :: i ! Loop index.

  write(*,*)'OrderCharCI, parameter C1 is "', c1%key, '" ', c1%data, ', len is ', len(c1%key)

  write(*,*)' len_trim is ', len_trim(c1%key)
  write(*,*) ' '

  ! SEGV on access to C2
  !
  write(*,*)'OrderCharCI, parameter C2 is "', c2%key, '" ', c2%data, ', len is ', len(c2%key)
  write(*,*)' len_trim is ', len_trim(c2%key)
  write(*,*) ' '

  c1L = c1%key
  c2L = c2%key

  write(*,*) 'about to start do loop'
  do i = 1, len_trim(C1L)
     if ('a' <= C1L(i:i) .and. c1L(i:i) <= 'z') c1L(i:i) = char(ichar(c1L(i:i))- ichar('a') + ichar('A'))
  end do

  do i = 1, len_trim(C2L)
     if ('a' <= c2L(i:i) .and. c2L(i:i) <= 'z') c2L(i:i) = char(ichar(c2L(i:i)) - ichar('a') + ichar('A'))
  end do

  if (c1L == c2L) Then
     OrderCharCI = 0
     write(*,*) ' - equal'
  else if (c1L < c2L) Then
     OrderCharCI = -1
     write(*,*) ' - c1 is less'
  else
     OrderCharCI = 1
     write(*,*) ' - c1 is more'
  end if


end function OrderCharCI
