
!
!     Class for the energy equation.
!
MODULE EnergyMod
  USE HamiltonianMod
  USE CorrEnergyMod
  USE RefEnergyMod
  IMPLICIT NONE
  
  TYPE Energy
     CLASS(Hamiltonian), POINTER :: hamilt        ! Hamiltonian operator
     CHARACTER(30) :: approximation               ! Chosen approximation
     TYPE(RefEnergy), POINTER :: eneRef           ! Reference energy
     CLASS(CorrEnergy), POINTER :: eneCorrel => NULL() ! Correlation energy
   CONTAINS
     PROCEDURE :: calculateEnergy
     PROCEDURE :: printEnergy
  END TYPE Energy
  
CONTAINS
  
  !
  !     Constructor for Energy object
  !
  FUNCTION Energy_(ham, approx, eneRef, eneCorr) RESULT(this)
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    CHARACTER(30), INTENT(IN) :: approx
    TYPE(RefEnergy), TARGET, INTENT(IN) :: eneRef
    CLASS(CorrEnergy), OPTIONAL, TARGET, INTENT(IN) :: eneCorr
    TYPE(Energy), POINTER :: this
    
    ALLOCATE(this)
    
    this%hamilt => ham
    this%approximation = approx
    this%eneRef => eneRef
    IF (PRESENT(eneCorr)) THEN
       this%eneCorrel => eneCorr
    ENDIF
    
    
  END FUNCTION Energy_
  
  !
  !     Destructor for Energy object
  !
  SUBROUTINE Energy_d(this)
    TYPE(Energy), TARGET, INTENT(INOUT) :: this
    
    IF (ASSOCIATED(this%hamilt)) THEN
       CALL Hamiltonian_d(this%hamilt)
       !DEALLOCATE(this%hamilt)
       this%hamilt => NULL()
    ENDIF
    
    IF (ASSOCIATED(this%eneRef)) THEN
       CALL RefEnergy_d(this%eneRef)
       !DEALLOCATE(this%eneRef)
       this%eneRef => NULL()
    ENDIF
    
    IF (ASSOCIATED(this%eneCorrel)) THEN
       CALL this%eneCorrel%CorrEnergy_d()
       DEALLOCATE(this%eneRef)
       this%eneCorrel => NULL()
    ENDIF
    
    
  END SUBROUTINE Energy_d
  
  SUBROUTINE calculateEnergy(this)
    CLASS(Energy), TARGET, INTENT(INOUT) :: this
    
    !     Calculate the reference energy
    !     (here the Hartree-Fock energy)
    CALL this%eneRef%calcRefEnergy()
    !     Print the energy values
    CALL this%eneRef%printRefEnergy()
    
    IF (this%approximation /= 'ref') THEN
       !     Calculate the correlation energy
       CALL this%eneCorrel%calcCorrEnergy()
       !     Print the energy values
       CALL printCorrEnergy(this%eneCorrel, this%eneRef)
    
       CALL this%printEnergy()
    ENDIF
    
    
  END SUBROUTINE calculateEnergy
  
  !
  !     Routine to print energy values
  !
  SUBROUTINE printEnergy(this)
    CLASS(Energy), INTENT(IN) :: this
    TYPE(RefEnergy), POINTER :: refEnergy
    CLASS(CorrEnergy), POINTER :: corrEnergy
    
    refEnergy => this%eneRef
    corrEnergy => this%eneCorrel
    
    IF (refEnergy%spBasis%mpiJobs%iAm == 0) THEN
       
       IF (this%approximation /= 'ref') THEN
          OPEN(17, FILE='nshells_ene')
          !     Energies given in Rydbergs (1 Ry = e^2/(2*a_{0}),
          !     a_{0} = hbar^2/(m*e^2))
          WRITE(17,'(I5,e14.5,e14.5,e14.5)') &
               refEnergy%spBasis%nShells, &
               2.d0*refEnergy%eneRef/refEnergy%spBasis%nOccupied, &
               2.d0*corrEnergy%eneCorr/refEnergy%spBasis%nOccupied
          CLOSE(17)
       ENDIF
    ENDIF
    
    
  END SUBROUTINE printEnergy
  
  
END MODULE EnergyMod
