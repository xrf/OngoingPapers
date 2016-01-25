
!
!     Class for the single-particle basis.
!
!     This version is for a plane wave basis
!     of a three-dimensional homogeneous system
!     with periodic boundary conditions
!     (for example the three-dimensional electron
!     gas)
!
MODULE SpBasisPW3dMod
  USE ConstantsMod
  USE SpBasisMod
  IMPLICIT NONE
  
  TYPE, EXTENDS(SpBasis) :: SpBasisPW3d
     !     'normal': Total spin zero
     !     'polarized': Spin polarized
     CHARACTER(LEN=100) :: spinState
     REAL(DP) :: volume
   CONTAINS
     PROCEDURE :: setupSpBasis => setupSpBasisPW3d
  END TYPE SpBasisPW3d
  
  
CONTAINS
  
  !
  !     Constructor for SpBasisPw3d object
  !
  FUNCTION SpBasisPW3d_(nOccupied, nShells, &
       spinState, volume, mpiJobs) RESULT(this) 
    INTEGER, INTENT(IN) :: nOccupied, nShells
    CHARACTER(LEN=*), INTENT(IN) :: spinState
    REAL(DP), INTENT(IN) :: volume
    TYPE(Mpi), TARGET, INTENT(IN) :: mpiJobs
    TYPE(SpBasisPW3d), POINTER :: this
    
    ALLOCATE(this)
    
    this%nOccupied = nOccupied
    this%nShells = nShells
    this%spinState = spinState
    this%volume = volume
    this%mpiJobs => mpiJobs
    
    
  END FUNCTION SpBasisPW3d_
  
  !
  !     Set up single-particle basis    
  !
  SUBROUTINE setupSpBasisPW3d(this)
    CLASS(SpBasisPW3d), TARGET, INTENT(INOUT) :: this
    INTEGER :: nStates, state, shell, eneInteger
    LOGICAL :: isShell
    INTEGER :: nMax, nx, ny, nz, ms, orbitsShell
    REAL(DP) :: factor, area
    
    !     Factor in the single-particle energy
    area = (this%volume)**(2.d0/3.d0)
    factor = 2.d0*pi**2/area
    
    !     Count single-particle states
    nStates = 0
    eneInteger = 0
    shell = 1
    DO WHILE (shell <= this%nShells)
       
       isShell = .FALSE.
       nMax = FLOOR(SQRT(REAL(eneInteger)))
       DO nx=-nMax, nMax
          DO ny=-nMax, nMax
             DO nz=-nMax, nMax
                IF ((nx**2+ny**2+nz**2) == eneInteger) THEN
                   IF (nMax > this%nMax) this%nMax = nMax
                   isShell = .TRUE.
                   nStates = nStates + 1
                   
                   !     If not spin-polarized
                   IF (this%spinState == 'normal') THEN
                      nStates = nStates + 1            
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       eneInteger = eneInteger + 1
       
       IF (isShell) shell = shell + 1
       
    ENDDO
    this%nUnoccupied = nStates - this%nOccupied
    IF (this%mpiJobs%iAm == 0) THEN
       WRITE(*,*) 'Total number of orbitals: ', nStates
    ENDIF
    
    ALLOCATE(this%spQuantNr(nStates, 4))
    ALLOCATE(this%spEnergies(nStates))
    this%spQuantNr = 0.d0
    this%spEnergies = 0.d0
    
    CALL this%mpiJobs%writeMpi(' ')
    !     Set up single-particle basis
    eneInteger = 0
    shell = 1
    state = 0
    DO WHILE (shell <= this%nShells)
       
       isShell = .FALSE.
       nMax = FLOOR(SQRT(REAL(eneInteger)))
       DO nx=-nMax, nMax
          DO ny=-nMax, nMax
             DO nz=-nMax, nMax
                
                IF ((nx**2+ny**2+nz**2) == eneInteger) THEN
                   !IF (.NOT.isShell) THEN
                   !   WRITE(*,*) 'eneInt=',eneInteger
                   !ENDIF
                   isShell = .TRUE.
                   state = state + 1
                   
                   !     Store (nx, ny, nz, ms)
                   this%spQuantNr(state, 1) = nx
                   this%spQuantNr(state, 2) = ny
                   this%spQuantNr(state, 3) = nz
                   this%spQuantNr(state, 4) = -1
                   !     Store single-particle energy
                   this%spEnergies(state) = eneInteger*factor
                   
                   IF (this%mpiJobs%iAm == 0) THEN
                      WRITE(*,*) 'nx=',nx,',ny=',ny,',nz=',nz,&
                           ',ms=',-1
                   ENDIF
                   
                   !     If not spin-polarized
                   IF (this%spinState == 'normal') THEN
                      state = state + 1
                      !     Store (nx, ny, nz, ms)
                      this%spQuantNr(state, 1) = nx
                      this%spQuantNr(state, 2) = ny
                      this%spQuantNr(state, 3) = nz
                      this%spQuantNr(state, 4) = 1
                      !     Store single-particle energy
                      this%spEnergies(state) = eneInteger*factor
                      
                      IF (this%mpiJobs%iAm == 0) THEN
                         WRITE(*,*) 'nx=',nx,',ny=',ny,',nz=',nz,&
                              ',ms=',1
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       IF (this%mpiJobs%iAm == 0) THEN
          WRITE(*,*) 'shell: ',shell, ', orbitals: ',state
       ENDIF
       
       eneInteger = eneInteger + 1
       
       IF (isShell) shell = shell + 1
       
    ENDDO
    CALL this%mpiJobs%writeMpi(' ')
    
    
  END SUBROUTINE setupSpBasisPW3d
  
  
END MODULE SpBasisPW3dMod
