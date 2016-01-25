
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
     TYPE(FockMatrix), POINTER :: fockMatr                 
     CLASS(SpBasis), POINTER :: basisSp           
     TYPE(BlockList), POINTER :: blocks 
     REAL(DP) :: eneCorr
   CONTAINS
     !     Calculate correlation energy
     PROCEDURE(calcCorrEne), &
          DEFERRED :: calcCorrEnergy  
     !     Destructor
     PROCEDURE(CorrEne_d), &
          DEFERRED :: CorrEnergy_d
  END TYPE CorrEnergy
  
  ABSTRACT INTERFACE
     SUBROUTINE calcCorrEne(this)
       IMPORT CorrEnergy
       CLASS(CorrEnergy), TARGET, INTENT(INOUT) :: this
     END SUBROUTINE calcCorrEne

     SUBROUTINE CorrEne_d(this)
       IMPORT CorrEnergy
       CLASS(CorrEnergy), TARGET, INTENT(INOUT) :: this
     END SUBROUTINE CorrEne_d
  END INTERFACE
  
CONTAINS
  
  SUBROUTINE printCorrEnergy(this, eneRef)
    CLASS(CorrEnergy), POINTER, INTENT(IN) :: this
    CLASS(RefEnergy), INTENT(IN) :: eneRef
    TYPE(Mpi), POINTER :: mpiJobs
    REAL(DP) :: eneTot
    INTEGER :: iError
    INCLUDE 'mpif.h'

    !     Wait until all processes are ready
    CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
    
    IF (this%basisSp%mpiJobs%iAm == 0) THEN
       
       !     Energies given in Rydbergs (1 Ry = e^2/(2*a_{0}),
       !     a_{0} = hbar^2/(m*e^2))
       eneTot = eneRef%eneRef + this%eneCorr
       WRITE(*,*) ' '
       WRITE(*,*) 'In Rydbergs:'
       WRITE(*,*) 'eneCorrel/N = ', &
            2.d0*this%eneCorr/this%basisSp%nOccupied, &
            ', eneTot/N = ', 2.d0*eneTot/this%basisSp%nOccupied
       WRITE(*,*) 'In Hartree:'
       WRITE(*,*) 'eneCorrel/N = ', this%eneCorr, &
            ', eneTot/N = ', eneTot
       WRITE(*,*) ' '
       
       OPEN(15, FILE='ene_ref_corr_tot')
       WRITE(15,'(e14.5,e14.5,e14.5)') &
            2.d0*eneRef%eneRef/this%basisSp%nOccupied, &
            2.d0*this%eneCorr/this%basisSp%nOccupied, &
            2.d0*eneTot/this%basisSp%nOccupied
       CLOSE(15)
       
    ENDIF
    
    
  END SUBROUTINE printCorrEnergy
  
  
  SUBROUTINE printCorrEnergyOnly(this)
    CLASS(CorrEnergy), POINTER, INTENT(IN) :: this
    
    IF (this%basisSp%mpiJobs%iAm == 0) THEN
       
       WRITE(*,*) ' '
       WRITE(*,*) 'eneCorrel/N = ', this%eneCorr&
            /this%basisSp%nOccupied
       WRITE(*,*) ' '
       
    ENDIF
    
    
  END SUBROUTINE printCorrEnergyOnly
  
  
END MODULE CorrEnergyMod
