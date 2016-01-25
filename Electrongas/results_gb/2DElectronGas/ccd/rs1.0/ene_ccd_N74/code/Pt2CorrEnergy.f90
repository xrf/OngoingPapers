
!
!     Class for the MBPT(2) correlation energy. This class 
!     inherits the general CorrEnergy class.
!
MODULE Pt2CorrEnergyMod
  USE CorrEnergyMod
  IMPLICIT NONE
  
  TYPE, EXTENDS(CorrEnergy) :: Pt2CorrEnergy
   CONTAINS
     PROCEDURE :: calcCorrEnergy => calcPt2Energy 
     PROCEDURE :: CorrEnergy_d => Pt2CorrEnergy_d
  END TYPE Pt2CorrEnergy
  
  
CONTAINS
  
  !
  !     Constructor for Pt2CorrEnergy object.
  !
  FUNCTION Pt2CorrEnergy_(ham, fockMatr, basis) &
       RESULT(this)
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    TYPE(FockMatrix), TARGET, INTENT(IN) :: fockMatr
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    TYPE(Pt2CorrEnergy), POINTER :: this
    
    ALLOCATE(this)
    
    this%hamilt => ham
    this%fockMatr => fockMatr
    this%basisSp => basis
    this%blocks => ham%blocks
    
    
  END FUNCTION Pt2CorrEnergy_
  
  !
  !     Destructor for Pt2CorrEnergy object.
  !
  SUBROUTINE Pt2CorrEnergy_d(this)
    CLASS(Pt2CorrEnergy), TARGET, INTENT(INOUT) :: this
    
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
    
    
  END SUBROUTINE Pt2CorrEnergy_d

  !
  !     Routine to calculate the correlation energy in
  !     second order perturbation theory (MBPT(2)).
  !
  SUBROUTINE calcPt2Energy(this) 
    CLASS(Pt2CorrEnergy), TARGET, INTENT(INOUT) :: this
    REAL(DP), ALLOCATABLE :: vHhpp(:,:)
    REAL(DP), ALLOCATABLE :: vHhpp2(:,:)
    INTEGER :: nBra, nKet, hh, pphhBlock, iError
    REAL(DP) :: temp, tempTotal
    INCLUDE 'mpif.h'

    temp = 0.d0
    !     Loop over (K_CM, M_S) blocks.
    DO pphhBlock=1, this%blocks%nPphhBlocksThis
       
       nBra = this%hamilt%vHhpp(pphhBlock)%dim1
       nKet = this%hamilt%vHhpp(pphhBlock)%dim2
       ALLOCATE(vHhpp(nBra, nKet))
       ALLOCATE(vHhpp2(nBra, nBra))
       vHhpp = this%hamilt%vHhpp(pphhBlock)%matr
       
       !     Matrix multiplication vHhpp*vPphh
       CALL DGEMM('N', 'T', nBra, nBra, nKet, 1.d0, &
            vHhpp, nBra, vHhpp, nBra, 0.d0, vHhpp2, nBra)
       
       !     Trace
       DO hh=1, nBra
          temp = temp + vHhpp2(hh, hh)
       ENDDO
       
       DEALLOCATE(vHhpp, vHhpp2)
    ENDDO
    
    !     Sum the contributions temp from the different MPI
    !     processes
    CALL MPI_ALLREDUCE(temp, tempTotal, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, iError)
    
    this%eneCorr = 0.25d0*tempTotal
    
    
  END SUBROUTINE calcPt2Energy
  
  
!   !
!   !     Routine to calculate the correlation energy in
!   !     second order perturbation theory (MBPT(2)).
!   !
!   SUBROUTINE calcPt2Energy(this) 
!     CLASS(Pt2CorrEnergy), TARGET, INTENT(INOUT) :: this
!     INTEGER :: nHoles, nParticles
!     REAL(DP) :: eneKh1, eneKh2, eneKp1, eneKp2
!     INTEGER :: hole1, hole2, particle1, particle2
!     REAL(DP) :: temp, interaction
    
!     nHoles = this%basisSp%nOccupied 
!     nParticles = this%basisSp%nUnoccupied
    
!     temp = 0.d0
!     DO hole1=1, nHoles
!        DO hole2=hole1+1, nHoles
!           DO particle1=nHoles+1, nHoles+nParticles
!              DO particle2=particle1+1, nHoles+nParticles
                
!                 !     Single-particle potential
!                 eneKh1 = this%fockMatr%matrix(&
!                      hole1, hole1)
!                 eneKh2 = this%fockMatr%matrix(&
!                      hole2, hole2)
!                 eneKp1 = this%fockMatr%matrix(&
!                      particle1, particle1)
!                 eneKp2 = this%fockMatr%matrix(&
!                      particle2, particle2)
                
!                 !     Two-particle interaction
! !                interaction = this%hamilt%vMatrix(&
! !                     hole1, hole2, particle1, particle2)
                
!                 temp = temp + interaction**2&
!                      /(eneKh1 + eneKh2 - eneKp1 - eneKp2)             
!              ENDDO
!           ENDDO
!        ENDDO
!     ENDDO
!     this%eneCorr = temp
    
    
!   END SUBROUTINE calcPt2Energy
  
  
  
END MODULE Pt2CorrEnergyMod
