PROGRAM testSort
  IMPLICIT NONE
  
  type mytype
     integer :: a,b,c
     real :: x,y,z
     real, allocatable :: mydata(:)
  end type mytype

  type listtype
     type(mytype), allocatable :: cell(:)
     type(mytype) :: acell
  end type listtype
  
  integer :: i
  type(listtype) :: this
 
  allocate(this%cell(1000))
  do i = 1,1000
     if (i <= 500) then
        allocate(this%cell(i)%mydata(10))
     else
        allocate(this%cell(i)%mydata(100))
     endif
  enddo
  this%cell(10)%a = 8
  write(*,*) 'a10=',8
  
  
END PROGRAM testSort

!
!     This routine is a modified version of a quicksort
!     routine given in W. H. Press et al., Numerical Recipes
!     in Fortran 90, Second Edition (1996).
!
! SUBROUTINE sort(array)
!   USE nrtype
!   USE nrutil, ONLY :: swap, nrerror
  
  

! END SUBROUTINE sort
