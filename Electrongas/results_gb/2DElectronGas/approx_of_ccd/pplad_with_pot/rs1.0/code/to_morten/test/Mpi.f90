
!
!     Class for an MPI job.
!
MODULE MpiMod
  USE ConstantsMod
  IMPLICIT NONE
  
  TYPE :: Mpi
     INTEGER :: nProcesses
     INTEGER :: iAm, iError
   CONTAINS
     PROCEDURE :: writeMpi
  END TYPE Mpi
  
CONTAINS
  
  !
  !     Constructor for Mpi object
  !
  FUNCTION Mpi_() RESULT(this)
    TYPE(Mpi), POINTER :: this
    INCLUDE 'mpif.h'

    ALLOCATE(this)
    
    !     Initialize MPI
    CALL MPI_INIT(this%iError)
    !     Get the number of processes, nProcesses
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, &
         this%nProcesses, this%iError) 
    !     Get the identification number, iAm, 
    !     of this process
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, this%iAm, &
         this%iError) 
    
    !     Write the number of processes to the 
    !     standard output
    !CALL writeMpi('Number of processes: '+int2str(numprocs))
    
    
  END FUNCTION Mpi_
  
  !
  !     Destructor for Mpi object
  !
  SUBROUTINE Mpi_d()
    INTEGER :: iError
    
    CALL MPI_FINALIZE(iError) 
    
    
  END SUBROUTINE Mpi_d
  
  !
  !     Write only for process number 0
  !
  SUBROUTINE writeMpi(this, text)
    CLASS(Mpi), INTENT(IN) :: this
    CHARACTER (LEN=*) :: text
    
    IF (this%iAm == 0) THEN
       WRITE(*,*) text
    ENDIF
    
    
  END SUBROUTINE writeMpi
  
  
END MODULE MpiMod
