
!
!     Class for the CCD correlation energy. This class inherits
!     the CorrEnergy class.
!
MODULE CcdCorrEnergyMod
  USE CcdAmplitudeMod
  USE CorrEnergyMod
  IMPLICIT NONE
  
  TYPE, EXTENDS(CorrEnergy) :: CcdCorrEnergy
     TYPE(CcdAmplitude) :: ccdAmpl                ! CCD t2 amplitude
     INTEGER :: maxIter
     REAL(DP) :: tolerance
   CONTAINS
     PROCEDURE :: calcCorrEnergy => calcCcdEnergy
     PROCEDURE :: evaluateEnergy
     PROCEDURE :: evaluateEneTerm1
     PROCEDURE :: evaluateEneTerm2
     PROCEDURE :: evaluateEneTerm3
  END TYPE CcdCorrEnergy
  
  
CONTAINS
  
  !
  !     Constructor for CcdCorrEnergy
  !
  FUNCTION CcdCorrEnergy_(ham, fockMatr, basis, &
       blocks, maxIter, tolerance) RESULT(this)
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    TYPE(FockMatrix), INTENT(IN) :: fockMatr
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    TYPE(BlockList), INTENT(IN) :: blocks
    INTEGER, INTENT(IN) :: maxIter
    REAL(DP), INTENT(IN) :: tolerance
    TYPE(CcdCorrEnergy), POINTER :: this
    
    ALLOCATE(this)
    
    this%hamilt => ham
    this%fockMatr = fockMatr
    this%basisSp => basis
    this%blocks = blocks
    this%eneCorr = 0.d0    
    this%ccdAmpl = CcdAmplitude_(basis, blocks)
    this%maxIter = maxIter
    this%tolerance = tolerance
    
    
  END FUNCTION CcdCorrEnergy_
  
  !
  !     Function to evaluate the CCSD correlation 
  !     energy at a single iteration step.
  !
  FUNCTION evaluateEnergy(this) RESULT(energy)
    CLASS(CcdCorrEnergy), TARGET, INTENT(IN) :: this
    REAL(DP) :: energy
    REAL(DP) :: term1, term2, term3
    
    term1 = this%evaluateEneTerm1()
    term2 = this%evaluateEneTerm2()
    term3 = this%evaluateEneTerm3()
    energy = term1 + term2 + term3
    
    
  END FUNCTION evaluateEnergy
  
  !
  !     Calculate the energy term
  !     \sum_ia <i|f|a><a|t|i>
  !
  FUNCTION evaluateEneTerm1(this) RESULT(energy)
    CLASS(CcdCorrEnergy), TARGET, INTENT(IN) :: this
    REAL(DP) :: energy
    INTEGER :: nHoles, nParticles, hole, particle
    REAL(DP) :: temp, fock, t1
    
    nHoles = this%basisSp%nOccupied 
    nParticles = this%basisSp%nUnoccupied
    
    temp = 0.d0
    DO hole=1, nHoles
       DO particle=1, nParticles
          fock = this%FockMatr%matrix(hole, nHoles+particle)
          t1 = this%CcdAmpl%t1Matrix(particle, hole)
          temp = temp + fock*t1
       ENDDO
    ENDDO
    energy = temp
    
    
  END FUNCTION evaluateEneTerm1
  
  !
  !     Calculate the energy term 
  !     (1/4)*\sum_ijab <ij||ab><ab|t|ij>
  !
  FUNCTION evaluateEneTerm2(this) RESULT(energy)
    CLASS(CcdCorrEnergy), TARGET, INTENT(IN) :: this
    REAL(DP) :: energy
    INTEGER :: nHoles, nParticles
    INTEGER :: hole1, hole2, particle1, particle2
    REAL(DP) :: temp, interaction, tAmplitude
    
    nHoles = this%basisSp%nOccupied 
    nParticles = this%basisSp%nUnoccupied
    
    temp = 0.d0
    DO hole1=1, nHoles
       DO hole2=hole1+1, nHoles
          DO particle1=nHoles+1, nHoles+nParticles             
             DO particle2=particle1+1, nHoles+nParticles
                
                !     Two-particle interaction
!                interaction = this%hamilt%vMatrix(&
!                     hole1, hole2, particle1, particle2)
                
                !     CCD t2-amplitude matrix element
                tAmplitude = this%ccdAmpl%t2Matrix(particle1&
                     -nHoles, particle2-nHoles, hole1, hole2)
                
                temp = temp + interaction*tAmplitude
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    energy = temp
    
    
  END FUNCTION evaluateEneTerm2

  !
  !     Calculate the energy term
  !     (1/2)*\sum_ijab <ij||ab><a|t|i><b|t|j> 
  !
  FUNCTION evaluateEneTerm3(this) RESULT(energy)
    CLASS(CcdCorrEnergy), TARGET, INTENT(IN) :: this
    REAL(DP) :: energy
    INTEGER :: nHoles, nParticles
    INTEGER :: hole1, hole2, particle1, particle2
    REAL(DP) :: temp, interaction, t11, t12
    
    nHoles = this%basisSp%nOccupied 
    nParticles = this%basisSp%nUnoccupied
    
    temp = 0.d0
    DO hole1=1, nHoles
       DO hole2=1, nHoles
          DO particle1=1, nParticles
             DO particle2=1, nParticles
!                interaction = this%hamilt%vMatrix(&
!                     hole1, hole2, &
!                     nHoles+particle1, nHoles+particle2)
                t11 = this%CcdAmpl%t1Matrix(particle1, hole1)
                t12 = this%CcdAmpl%t1Matrix(particle2, hole2)
                
                temp = temp + interaction*t11*t12
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    energy = 0.5d0*temp
    

  END FUNCTION evaluateEneTerm3
  
  SUBROUTINE calcCcdEnergy(this)
    CLASS(CcdCorrEnergy), TARGET, INTENT(INOUT) :: this
    CLASS(CorrEnergy), POINTER :: eneC => NULL()
    REAL(DP) :: energy, energyOld, difference
    INTEGER :: iteration
    
    !     Set up the energy denominator matrix list
    CALL setupEnergyDenomMatrix(this%ccdAmpl, &
         this%fockMatr)
    
    WRITE(*,*) ' '
    WRITE(*,*) 'Self-consistency loop:'
    iteration = 1
    difference = 100.d0
    energyOld = 100000.d0
    DO WHILE ((difference > this%tolerance)&
         .AND.(iteration <= this%maxIter))
       WRITE(*,*) ' '
       WRITE(*,*) 'Iteration number ', iteration
       
       WRITE(*,*) 'Updating intermediate diagrams...'
       CALL updateIntermed1(this%ccdAmpl, this%hamilt)
       CALL updateIntermed2(this%ccdAmpl, this%hamilt)
       CALL updateIntermed3(this%ccdAmpl, this%hamilt)
       CALL updateIntermed4(this%ccdAmpl, this%hamilt)
       
       ! WRITE(*,*) 'Updating the t1-amplitude matrix...'
       ! !     Update the t1-amplitude matrix
       ! CALL updateT1Matrix(this%ccdAmpl, this%hamilt, &
       !      this%fockMatr)
       WRITE(*,*) 'Updating the t2-amplitude matrix...'
       !     Update the t2-amplitude matrix 
       CALL updateT2Matrix(this%ccdAmpl, this%hamilt, &
            this%fockMatr)
       
       this%ccdAmpl%t1MatrixOld = this%ccdAmpl%t1Matrix
       this%ccdAmpl%t2MatrixOld = this%ccdAmpl%t2Matrix
       
       WRITE(*,*) 'Calculating the correlation energy...'
       !     Evaluate the energy
       energy = this%evaluateEnergy()
       this%eneCorr = energy
       
       eneC => this
       !     Print the energy values
       CALL printCorrEnergyOnly(eneC)
       
       !     Get the convergency measure
       difference = ABS(energy - energyOld)
       energyOld = energy
       WRITE(*,*) 'Difference: ', difference
       
       iteration = iteration + 1
    ENDDO
    WRITE(*,*) ' '
    WRITE(*,*) '------  Selfconsistency loop done ------'
    WRITE(*,*) ' '
    
    
  END SUBROUTINE calcCcdEnergy
  
  
  
END MODULE CcdCorrEnergyMod
