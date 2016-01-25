
!
!     Class for the CCD correlation energy. This class inherits
!     the CorrEnergy class.
!
MODULE CcdCorrEnergyMod
  USE CcdAmplitudeMod
  USE CorrEnergyMod
  USE MpiMod
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
     PROCEDURE :: CorrEnergy_d => CcdCorrEnergy_d
     PROCEDURE :: evaluateEnergy
     PROCEDURE :: evaluateEnergyRe
  END TYPE CcdCorrEnergy
  
CONTAINS
  
  !
  !     Constructor for CcdCorrEnergy object.
  !
  FUNCTION CcdCorrEnergy_(ham, fockMatr, basis, &
       maxIter, tolerance, alpha) RESULT(this)
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    TYPE(FockMatrix), TARGET, INTENT(IN) :: fockMatr
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    INTEGER, INTENT(IN) :: maxIter
    REAL(DP), INTENT(IN) :: tolerance, alpha
    TYPE(CcdCorrEnergy), POINTER :: this
    
    ALLOCATE(this)
    
    this%hamilt => ham
    this%fockMatr => fockMatr
    this%basisSp => basis
    this%blocks => ham%blocks
    this%eneCorr = 0.d0    
    this%ccdAmpl => CcdAmplitude_(&
         basis, ham%blocks, ham%blocksRe)
    
    this%maxIter = maxIter
    this%tolerance = tolerance
    this%alpha = alpha
    
    
  END FUNCTION CcdCorrEnergy_
  
  !
  !     Destructor for CcdCorrEnergy object.
  !
  SUBROUTINE CcdCorrEnergy_d(this)
    CLASS(CcdCorrEnergy), TARGET, INTENT(INOUT) :: this
    
    IF (ASSOCIATED(this%hamilt)) THEN
       CALL Hamiltonian_d(this%hamilt)
       this%hamilt => NULL()
    ENDIF
    
    IF (ASSOCIATED(this%fockMatr)) THEN
       CALL FockMatrix_d(this%fockMatr)
       this%fockMatr => NULL()
    ENDIF
    
    IF (ASSOCIATED(this%basisSp)) THEN
       CALL SpBasis_d(this%basisSp)
       this%basisSp => NULL()
    ENDIF
    
    IF (ASSOCIATED(this%blocks)) THEN
       CALL BlockList_d(this%blocks)
       this%blocks => NULL()
    ENDIF
    
    IF (ASSOCIATED(this%ccdAmpl)) THEN
       CALL CcdAmplitude_d(this%ccdAmpl)
       DEALLOCATE(this%ccdAmpl)
       this%ccdAmpl => NULL()
    ENDIF
  
    
  END SUBROUTINE CcdCorrEnergy_d
  
  !
  !     Function to evaluate the CCD correlation 
  !     energy at a single iteration step.
  !
  FUNCTION evaluateEnergy(this) RESULT(energy)
    CLASS(CcdCorrEnergy), TARGET, INTENT(IN) :: this
    REAL(DP), ALLOCATABLE :: vHhpp(:,:)
    REAL(DP), ALLOCATABLE :: tPPhh(:,:)
    REAL(DP), ALLOCATABLE :: vt(:,:)
    INTEGER :: nBra, nKet, hh, pphhBlock
    REAL(DP) :: energy, temp, tempTotal
    INTEGER :: iError
    INCLUDE 'mpif.h'
    
    temp = 0.d0
    !     Loop over (K_CM, M_S) blocks.
    DO pphhBlock=1, this%blocks%nPphhBlocksThis
       
       nBra = this%hamilt%vHhpp(pphhBlock)%dim1
       nKet = this%hamilt%vHhpp(pphhBlock)%dim2
       ALLOCATE(vHhpp(nBra, nKet))
       ALLOCATE(tPPhh(nKet, nBra))
       ALLOCATE(vt(nBra, nBra))
       vHhpp = this%hamilt%vHhpp(pphhBlock)%matr
       tPPhh = this%ccdAmpl%tPphh(pphhBlock)%matr
       
       !     Matrix multiplication vHhpp*tPphh
       CALL DGEMM('N', 'N', nBra, nBra, nKet, 1.d0, &
            vHhpp, nBra, tPphh, nKet, 0.d0, vt, nBra)
       
       !     Trace
       DO hh=1, nBra
          temp = temp + vt(hh, hh)
       ENDDO
       
       DEALLOCATE(vHhpp, tPphh, vt)
    ENDDO
    !     Sum the contributions temp from the different MPI
    !     processes
    CALL MPI_ALLREDUCE(temp, tempTotal, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, iError)
    
    energy = 0.25d0*tempTotal
    
    
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
    REAL(DP) :: energy, temp, tempTotal
    INTEGER :: blockRe, iError
    INCLUDE 'mpif.h'
    
    temp = 0.d0
    DO blockRe=1, this%hamilt%blocksRe%nBlocksReThis
       blockIndex = &
            this%hamilt%blocksRe%blocksThis2BlocksAll(blockRe)
       
       nPhPairs = this%hamilt%vPphhRe(blockRe)%dim1
       nHpPairs = this%hamilt%vPphhRe(blockRe)%dim2
       
       ALLOCATE(vPphhRe(nPhPairs, nPhPairs))
       ALLOCATE(tPphhRe(nPhPairs, nPhPairs))
       ALLOCATE(vt(nPhPairs, nPhPairs))   
       
       vPphhRe = this%hamilt%vPphhRe(blockRe)%matr
       tPphhRe = this%ccdAmpl%tPphhRe(blockIndex)%matr
       
       !     Matrix multiplication vPhphRe*tPhphRe
       CALL DGEMM('N', 'N', nHpPairs, nHpPairs, nPhPairs, &
            1.d0, vPphhRe, nHpPairs, tPphhRe, nPhPairs, &
            0.d0, vt, nHpPairs)
       
       !     Trace
       DO hp=1, nHpPairs
          temp = temp + vt(hp, hp)
       ENDDO
       
       DEALLOCATE(vPphhRe, tPphhRe, vt)
    ENDDO
    
    !     Sum the contributions temp from the different MPI
    !     processes
    CALL MPI_ALLREDUCE(temp, tempTotal, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, iError)

    energy = 0.25d0*tempTotal
    
    
  END FUNCTION evaluateEnergyRe
  
  
  SUBROUTINE calcCcdEnergy(this)
    CLASS(CcdCorrEnergy), TARGET, INTENT(INOUT) :: this
    CLASS(CorrEnergy), POINTER :: eneC => NULL()
    REAL(DP) :: energy, energyOld, difference
    CLASS(SpBasis), POINTER :: basis
    INTEGER :: iteration
    
    !     Set up the energy denominator matrix list
    CALL setupEnergyDenomMatrices(this%ccdAmpl, &
         this%fockMatr)
    
    basis => this%basisSp
    IF (basis%mpiJobs%iAm == 0) THEN
       WRITE(*,*) ' '
       WRITE(*,*) 'Self-consistency loop:'
    ENDIF
    
    iteration = 1
    difference = 100.d0
    energyOld = 100000.d0
    DO WHILE ((difference > this%tolerance)&
         .AND.(iteration <= this%maxIter))
       
       IF (basis%mpiJobs%iAm == 0) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'Iteration number ', iteration
          
          WRITE(*,*) 'Updating intermediate diagrams...'
       ENDIF
       
       CALL this%ccdAmpl%updateIntermed1(this%hamilt)
       CALL this%ccdAmpl%updateIntermed2(this%hamilt)
       CALL this%ccdAmpl%updateIntermed3(this%hamilt)
       CALL this%ccdAmpl%updateIntermed4(this%hamilt)
       
       CALL basis%mpiJobs%writeMpi(&
            'Updating the t2-amplitude matrix...')
       CALL this%ccdAmpl%updateT2Matrix(this%hamilt, &
            this%alpha, iteration)
       CALL basis%mpiJobs%writeMpi(&
            'Updating the cross-coupled t2-amplitude matrix...')
       CALL this%ccdAmpl%setupT2MatrixRe()
       
       CALL basis%mpiJobs%writeMpi(&
            'Calculating the correlation energy...')
       !     Evaluate the energy
       energy = this%evaluateEnergy()
       !energy = this%evaluateEnergyRe()
       this%eneCorr = energy
       
       eneC => this
       !     Print the energy values
       CALL printCorrEnergyOnly(eneC)
       
       !     Get the convergency measure
       difference = ABS(energy - energyOld)
       energyOld = energy
       
       IF (basis%mpiJobs%iAm == 0) THEN
          WRITE(*,*) 'Difference: ', difference
       ENDIF
       
       iteration = iteration + 1
       
    ENDDO

    IF (basis%mpiJobs%iAm == 0) THEN
       WRITE(*,*) ' '
       WRITE(*,*) '------  Selfconsistency loop done ------'
       WRITE(*,*) ' '
    ENDIF
    
    
  END SUBROUTINE calcCcdEnergy
  
  
END MODULE CcdCorrEnergyMod
