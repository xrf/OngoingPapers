
!
!     Class for the CCD correlation energy. This class inherits
!     the CorrEnergy class.
!
MODULE CcdCorrEnergyMod
  USE CcdAmplitudeMod
  USE CorrEnergyMod
  IMPLICIT NONE
  
  TYPE, EXTENDS(CorrEnergy) :: CcdCorrEnergy
     !
     !     ccdAmpl:      CCD t2-amplitude
     !
     TYPE(CcdAmplitude), POINTER :: ccdAmpl               
     INTEGER :: maxIter
     REAL(DP) :: tolerance
     REAL(DP) :: alpha
   CONTAINS
     PROCEDURE :: calcCorrEnergy => calcCcdEnergy
     PROCEDURE :: evaluateEnergy
     PROCEDURE :: evaluateEnergyRe
  END TYPE CcdCorrEnergy
  
CONTAINS
  
  !
  !     Constructor for CcdCorrEnergy
  !
  FUNCTION CcdCorrEnergy_(ham, fockMatr, basis, &
       maxIter, tolerance, alpha) RESULT(this)
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    TYPE(FockMatrix), INTENT(IN) :: fockMatr
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    INTEGER, INTENT(IN) :: maxIter
    REAL(DP), INTENT(IN) :: tolerance, alpha
    TYPE(CcdCorrEnergy), POINTER :: this
    
    ALLOCATE(this)
    
    this%hamilt => ham
    this%fockMatr = fockMatr
    this%basisSp => basis
    this%blocks = ham%blocks
    this%eneCorr = 0.d0    
    this%ccdAmpl => CcdAmplitude_(&
         basis, ham%blocks, ham%blocksRe)
    
    this%maxIter = maxIter
    this%tolerance = tolerance
    this%alpha = alpha
    
    
  END FUNCTION CcdCorrEnergy_
  
  !
  !     Function to evaluate the CCD correlation 
  !     energy at a single iteration step.
  !
  FUNCTION evaluateEnergy(this) RESULT(energy)
    CLASS(CcdCorrEnergy), TARGET, INTENT(IN) :: this
    REAL(DP), ALLOCATABLE :: vHhpp(:,:)
    REAL(DP), ALLOCATABLE :: tPPhh(:,:)
    REAL(DP), ALLOCATABLE :: vt(:,:)
    INTEGER :: blockIndex, nBra, nKet, hh, pphhBlock
    REAL(DP) :: energy, temp
    
    temp = 0.d0
    !     Loop over (K_CM, M_S) blocks.
    DO pphhBlock=1, this%blocks%nPphhBlocks
       blockIndex = this%blocks%pphhBlocks(pphhBlock)
       
       nBra = this%hamilt%vHhpp(blockIndex)%dim1
       nKet = this%hamilt%vHhpp(blockIndex)%dim2
       ALLOCATE(vHhpp(nBra, nKet))
       ALLOCATE(tPPhh(nKet, nBra))
       ALLOCATE(vt(nBra, nBra))
       vHhpp = this%hamilt%vHhpp(blockIndex)%matr
       tPPhh = this%ccdAmpl%tPphh(pphhBlock)%matr

       !     Matrix multiplication vHhpp*tPphh
!       CALL DGEMM('N', 'N', nBra, nBra, nKet, 1.d0, &
!            vHhpp, nBra, tPphh, nKet, 0.d0, vt, nBra)
       
       !     Trace
       DO hh=1, nBra
          temp = temp + vt(hh, hh)
       ENDDO
       
       DEALLOCATE(vHhpp, tPphh, vt)
    ENDDO
    energy = 0.25d0*temp
    
    
  END FUNCTION evaluateEnergy
  
  !
  !     Function to evaluate the CCD correlation 
  !     energy at a single iteration step. Here
  !     evaluated using cross-coupled matrices.
  ! 
  FUNCTION evaluateEnergyRe(this) RESULT(energy)
    CLASS(CcdCorrEnergy), TARGET, INTENT(IN) :: this
    REAL(DP), ALLOCATABLE :: vPphhRe(:,:)
    REAL(DP), ALLOCATABLE :: tPphhRe(:,:)
    REAL(DP), ALLOCATABLE :: vt(:,:)
    INTEGER :: blockIndex, hp, nPhPairs, nHpPairs
    REAL(DP) :: energy, temp
    
    temp = 0.d0
    DO blockIndex=1, this%hamilt%blocksRe%nBlocks
       nPhPairs = this%hamilt%vPphhRe(blockIndex)%dim1
       nHpPairs = this%hamilt%vPphhRe(blockIndex)%dim2
       
       ALLOCATE(vPphhRe(nPhPairs, nPhPairs))
       ALLOCATE(tPphhRe(nPhPairs, nPhPairs))
       ALLOCATE(vt(nPhPairs, nPhPairs))   
       
       vPphhRe = this%hamilt%vPphhRe(blockIndex)%matr
       tPphhRe = this%ccdAmpl%tPphhRe(blockIndex)%matr
       
       !     Matrix multiplication vPhphRe*tPhphRe
 !      CALL DGEMM('N', 'N', nHpPairs, nHpPairs, nPhPairs, &
 !           1.d0, vPphhRe, nHpPairs, tPphhRe, nPhPairs, &
 !           0.d0, vt, nHpPairs)
       
       !     Trace
       DO hp=1, nHpPairs
          temp = temp + vt(hp, hp)
       ENDDO
       
       DEALLOCATE(vPphhRe, tPphhRe, vt)
    ENDDO
    energy = 0.25d0*temp
    
    
  END FUNCTION evaluateEnergyRe
  
  
  SUBROUTINE calcCcdEnergy(this)
    CLASS(CcdCorrEnergy), TARGET, INTENT(INOUT) :: this
    CLASS(CorrEnergy), POINTER :: eneC => NULL()
    REAL(DP) :: energy, energyOld, difference
    INTEGER :: iteration
    
    !     Set up the energy denominator matrix list
    CALL setupEnergyDenomMatrices(this%ccdAmpl, &
         this%fockMatr)
    
!    WRITE(*,*) ' '
!    WRITE(*,*) 'Self-consistency loop:'
    iteration = 1
    difference = 100.d0
    energyOld = 100000.d0
    DO WHILE ((difference > this%tolerance)&
         .AND.(iteration <= this%maxIter))
!       WRITE(*,*) ' '
!       WRITE(*,*) 'Iteration number ', iteration
       
!       WRITE(*,*) 'Updating intermediate diagrams...'
       CALL this%ccdAmpl%updateIntermed1(this%hamilt)
       CALL this%ccdAmpl%updateIntermed2(this%hamilt)
       CALL this%ccdAmpl%updateIntermed3(this%hamilt)
       CALL this%ccdAmpl%updateIntermed4(this%hamilt)
       
!       WRITE(*,*) 'Updating the t2-amplitude matrix...'
       CALL this%ccdAmpl%updateT2Matrix(this%hamilt, &
            this%alpha, iteration)
!       WRITE(*,*) 'Updating the cross-coupled t2-amplitude matrix...'
       CALL this%ccdAmpl%setupT2MatrixRe()
       
!       WRITE(*,*) 'Calculating the correlation energy...'
       !     Evaluate the energy
       energy = this%evaluateEnergy()
       !energy = this%evaluateEnergyRe()
!       write(*,*) 'correne calculated'
       this%eneCorr = energy
       
       eneC => this
       !     Print the energy values
       CALL printCorrEnergyOnly(eneC)
       
       !     Get the convergency measure
       difference = ABS(energy - energyOld)
       energyOld = energy
 !      WRITE(*,*) 'Difference: ', difference
       
       iteration = iteration + 1
       
    ENDDO
 !   WRITE(*,*) ' '
 !   WRITE(*,*) '------  Selfconsistency loop done ------'
 !   WRITE(*,*) ' '
    
    
  END SUBROUTINE calcCcdEnergy
  
  
END MODULE CcdCorrEnergyMod
