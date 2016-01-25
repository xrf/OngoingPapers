
!
!     Class for the single-particle basis.
!
!     This version is for a plane wave basis
!     of a two-dimensional homogeneous system
!     with periodic boundary conditions
!     (for example the two-dimensional electron
!     gas)
!
MODULE SpBasisPW2dMod
  USE ConstantsMod
  USE SpBasisMod
  IMPLICIT NONE
  
  TYPE, EXTENDS(SpBasis) :: SpBasisPW2d
     !     'normal':          Total spin zero
     !     'polarized':       Spin polarized
     CHARACTER(LEN=100) :: spinState   
     REAL(DP) :: area
   CONTAINS
     PROCEDURE :: setupSpBasis => setupSpBasisPW2d
  END TYPE SpBasisPW2d
  
  
CONTAINS
  
  !
  !     Constructor for SpBasisPW2d object
  !
  FUNCTION SpBasisPW2d_(nOccupied, nShells, spinState, &
       area, mpiJobs) RESULT(this)
    INTEGER, INTENT(IN) :: nOccupied, nShells
    CHARACTER(LEN=*), INTENT(IN) :: spinState
    REAL(DP), INTENT(IN) :: area
    TYPE(Mpi), TARGET, INTENT(IN) :: mpiJobs
    TYPE(SpBasisPW2d), POINTER :: this
    
    ALLOCATE(this)
    
    this%nOccupied = nOccupied
    this%nShells = nShells
    this%spinState = spinState
    this%area = area
    this%nMax = 0
    this%mpiJobs => mpiJobs
    
    
  END FUNCTION SpBasisPW2d_
  
  !
  !     Set up single-particle basis    
  !
  SUBROUTINE setupSpBasisPW2d(this)
    CLASS(SpBasisPW2d), TARGET, INTENT(INOUT) :: this
    INTEGER :: nStates, state, shell, eneInteger
    LOGICAL :: isShell
    INTEGER :: nmax, nx, ny, ms
    REAL(DP) :: factor
    
    !     Factor in the single-particle energy 
    factor = 2.d0*pi**2/this%area
    
    !     Count single-particle states
    nStates = 0
    eneInteger = 0
    shell = 1
    DO WHILE (shell <= this%nShells)
       
       isShell = .FALSE.
       nMax = FLOOR(SQRT(REAL(eneInteger)))
       DO nx=-nMax, nMax
          DO ny=-nMax, nMax
             IF ((nx**2+ny**2) == eneInteger) THEN
                this%nMax = nMax
                isShell = .TRUE.
                nStates = nStates + 1
                
                !     If not spin-polarized
                IF (this%spinState == 'normal') THEN
                   nStates = nStates + 1            
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       eneInteger = eneInteger + 1
       
       IF (isShell) shell = shell + 1
       
    ENDDO
    this%nUnoccupied = nStates - this%nOccupied
    WRITE(*,*) ' Number of orbitals:', nStates
    
    ALLOCATE(this%spQuantNr(nStates, 3))
    ALLOCATE(this%spEnergies(nStates))
    this%spQuantNr = 0.d0
    this%spEnergies = 0.d0
    
    WRITE(*,*) ' '
    WRITE(*,*) ' Sinlge-particle states:'
    WRITE(*,*) ' '
    !     Set up single-particle basis
    eneInteger = 0
    shell = 1
    state = 0
    WRITE(*,*) '    state   nx    ny    ms'
    DO WHILE (shell <= this%nShells)
       
       isShell = .FALSE.
       nMax = FLOOR(SQRT(REAL(eneInteger)))
       DO nx=-nMax, nMax
          DO ny=-nMax, nMax
             
             IF ((nx**2+ny**2) == eneInteger) THEN
                ! IF (.NOT.isShell) THEN
                !    WRITE(*,*) 'eneInt=',eneInteger
                ! ENDIF
                isShell = .TRUE.
                state = state + 1
                
                !     Store (nx, ny, ms)
                this%spQuantNr(state, 1) = nx
                this%spQuantNr(state, 2) = ny
                this%spQuantNr(state, 3) = -1
                !     Store single-particle energy
                this%spEnergies(state) = eneInteger*factor
                
                WRITE(*,'(4(x, I6))') state, nx, ny, -1
                
                !     If not spin-polarized
                IF (this%spinState == 'normal') THEN
                   state = state + 1
                   !     Store (nx, ny, ms)
                   this%spQuantNr(state, 1) = nx
                   this%spQuantNr(state, 2) = ny
                   this%spQuantNr(state, 3) = 1
                   !     Store single-particle energy
                   this%spEnergies(state) = eneInteger*factor
                   
                   WRITE(*,'(4(x, I6))') state, nx, ny, 1
                ENDIF
                
             ENDIF
          ENDDO
       ENDDO
       WRITE(*,*) 'shell: ',shell, ', orbitals: ',state
       
       eneInteger = eneInteger + 1
       
       IF (isShell) shell = shell + 1
       
    ENDDO
    WRITE(*,*) ' '
    
    
  END SUBROUTINE setupSpBasisPW2d
  
  
END MODULE SpBasisPW2dMod
