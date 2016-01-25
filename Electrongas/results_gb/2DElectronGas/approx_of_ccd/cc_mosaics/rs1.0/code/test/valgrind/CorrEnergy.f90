
!
!     Generic class for the correlation energy.
!
MODULE CorrEnergyMod
  USE HamiltonianMod
  USE FockMatrixMod
  USE RefEnergyMod
  USE ConstantsMod
  USE SpBasisMod
  IMPLICIT NONE
  
  TYPE, ABSTRACT :: CorrEnergy
     !
     !     hamilt:       Hamiltonian operator
     !     fockMatr:     Fock matrix
     !     spBasis:      Single-particle basis
     !     blocks:       List of (K_CM, M_S) blocks
     !     eneCorr:      Correlation energy
     !
     CLASS(Hamiltonian), POINTER :: hamilt
     TYPE(FockMatrix) :: fockMatr                 
     CLASS(SpBasis), POINTER :: basisSp           
     TYPE(BlockList) :: blocks 
     REAL(DP) :: eneCorr
   CONTAINS
     PROCEDURE(calcCorrEne), &
          DEFERRED :: calcCorrEnergy              ! Calculate correlation
                                                  !   energy
  END TYPE CorrEnergy
  
  ABSTRACT INTERFACE
     SUBROUTINE calcCorrEne(this)
       IMPORT CorrEnergy
       CLASS(CorrEnergy), TARGET, INTENT(INOUT) :: this
     END SUBROUTINE calcCorrEne
  END INTERFACE
  
CONTAINS
  
  SUBROUTINE printCorrEnergy(this, eneRef)
    CLASS(CorrEnergy), POINTER, INTENT(IN) :: this
    CLASS(RefEnergy), INTENT(IN) :: eneRef
    REAL(DP) :: eneTot
    
    !     Energies given in Rydbergs (1 Ry = e^2/(2*a_{0}),
    !     a_{0} = hbar^2/(m*e^2))
    eneTot = eneRef%eneRef + this%eneCorr
    ! WRITE(*,*) ' '
    ! WRITE(*,*) 'In Rydbergs:'
    ! WRITE(*,*) 'eneCorrel/N = ', &
    !      2.d0*this%eneCorr/this%basisSp%nOccupied, &
    !      ', eneTot/N = ', 2.d0*eneTot/this%basisSp%nOccupied
    ! WRITE(*,*) 'In Hartree:'
    ! WRITE(*,*) 'eneCorrel/N = ', this%eneCorr, &
    !      ', eneTot/N = ', eneTot
    ! WRITE(*,*) ' '
    
    ! OPEN(15, FILE='ene_ref_corr_tot')
    ! WRITE(15,'(e14.5,e14.5,e14.5)') &
    !      2.d0*eneRef%eneRef/this%basisSp%nOccupied, &
    !      2.d0*this%eneCorr/this%basisSp%nOccupied, &
    !      2.d0*eneTot/this%basisSp%nOccupied
    ! CLOSE(15)
    
    
  END SUBROUTINE printCorrEnergy
  
  
  SUBROUTINE printCorrEnergyOnly(this)
    CLASS(CorrEnergy), POINTER, INTENT(IN) :: this
    
    ! WRITE(*,*) ' '
    ! WRITE(*,*) 'eneCorrel/N = ', this%eneCorr&
    !      /this%basisSp%nOccupied
    ! WRITE(*,*) ' '
    
    
  END SUBROUTINE printCorrEnergyOnly
  
  
END MODULE CorrEnergyMod
