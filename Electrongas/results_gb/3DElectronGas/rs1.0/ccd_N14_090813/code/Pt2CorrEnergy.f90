
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
  END TYPE Pt2CorrEnergy
  
  
CONTAINS
  
  !
  !     Constructor for Pt2CorrEnergy.
  !
  FUNCTION Pt2CorrEnergy_(ham, fockMatr, basis) &
       RESULT(this)
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    TYPE(FockMatrix), INTENT(IN) :: fockMatr
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    TYPE(Pt2CorrEnergy), POINTER :: this
    
    ALLOCATE(this)
    
    this%hamilt => ham
    this%fockMatr = fockMatr
    this%basisSp => basis
    this%blocks = ham%blocks
    
    
  END FUNCTION Pt2CorrEnergy_
  
  !
  !     Routine to calculate the correlation energy in
  !     second order perturbation theory (MBPT(2)).
  !
  SUBROUTINE calcPt2Energy(this) 
    CLASS(Pt2CorrEnergy), TARGET, INTENT(INOUT) :: this
    REAL(DP), ALLOCATABLE :: vHhpp(:,:)
    REAL(DP), ALLOCATABLE :: vHhpp2(:,:)
    INTEGER :: blockIndex, nBra, nKet, hh, pphhBlock
    REAL(DP) :: temp
    
    temp = 0.d0
    !     Loop over (K_CM, M_S) blocks.
    DO pphhBlock=1, this%blocks%nPphhBlocks
       blockIndex = this%blocks%pphhBlocks(pphhBlock)
       
       nBra = this%hamilt%vHhpp(blockIndex)%dim1
       nKet = this%hamilt%vHhpp(blockIndex)%dim2
       ALLOCATE(vHhpp(nBra, nKet))
       ALLOCATE(vHhpp2(nBra, nBra))
       vHhpp = this%hamilt%vHhpp(blockIndex)%matr
       
       !     Matrix multiplication vHhpp*vPphh
       CALL DGEMM('N', 'T', nBra, nBra, nKet, 1.d0, &
            vHhpp, nBra, vHhpp, nBra, 0.d0, vHhpp2, nBra)
       
       !     Trace
       DO hh=1, nBra
          temp = temp + vHhpp2(hh, hh)
       ENDDO
       
       DEALLOCATE(vHhpp, vHhpp2)
    ENDDO
    this%eneCorr = 0.25d0*temp
    
    
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
