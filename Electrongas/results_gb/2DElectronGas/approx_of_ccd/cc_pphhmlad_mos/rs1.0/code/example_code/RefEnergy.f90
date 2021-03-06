
!
!     Class for the reference energy (here Hartree-Fock energy).
!
MODULE RefEnergyMod
  USE HamiltonianMod
  USE BlockListMod
  USE ConstantsMod
  USE SpBasisMod
  IMPLICIT NONE
  
  TYPE RefEnergy
     !
     !     hamilt:       Hamiltonian operator
     !     spBasis:      Single-particle basis
     !     blocks:       List of (K_CM, M_S) blocks
     !     eneKin:       Kinetic energy
     !     potRef:       Reference potential energy
     !     eneRef:       Reference energy
     !
     CLASS(Hamiltonian), POINTER :: hamilt  
     CLASS(SpBasis), POINTER :: spBasis      
     TYPE(BlockList), POINTER :: blocks              
     REAL(DP) :: eneKin, potRef, eneRef      
   CONTAINS
     PROCEDURE :: calcRefEnergy
     PROCEDURE :: printRefEnergy
  END TYPE RefEnergy
  
CONTAINS
  
  !
  !     Constructor for RefEnergy.
  !
  FUNCTION RefEnergy_(ham, basis, blocks) RESULT(this)
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    TYPE(BlockList), TARGET, INTENT(IN) :: blocks
    TYPE(RefEnergy), POINTER :: this
    
    ALLOCATE(this)
    
    this%hamilt => ham
    this%spBasis => basis
    this%blocks => blocks
    this%eneKin = 0.d0
    this%potRef = 0.d0
    this%eneRef = 0.d0
    
    
  END FUNCTION RefEnergy_
  
  !
  !     Destructor for RefEnergy object
  !
  SUBROUTINE RefEnergy_d(this)
    TYPE(RefEnergy), TARGET, INTENT(INOUT) :: this
    
    IF (ASSOCIATED(this%hamilt)) THEN
       CALL Hamiltonian_d(this%hamilt)
       this%hamilt => NULL()
    ENDIF
    
    IF (ASSOCIATED(this%spBasis)) THEN
       CALL SpBasis_d(this%spBasis)
       this%spBasis => NULL()
    ENDIF
    
    IF (ASSOCIATED(this%blocks)) THEN
       CALL BlockList_d(this%blocks)
       this%blocks => NULL()
    ENDIF
    
    
  END SUBROUTINE RefEnergy_d
  
  !
  !     Routine to calculate the reference energy,
  !     which is here the same as the Hartree-Fock
  !     energy.
  !  
  SUBROUTINE calcRefEnergy(this)
    CLASS(RefEnergy), TARGET, INTENT(INOUT) :: this
    TYPE(RealMatrix) :: vHhhhMatrix
    INTEGER :: orbit1, hh, nHoles, nParticles
    INTEGER :: hhBlock, blockIndex
    REAL(DP) :: temp, interaction
    integer :: o1,o2
    
    nHoles = this%spBasis%nOccupied 
    nParticles = this%spBasis%nUnoccupied
    
    !     Kinetic energy 
    temp = 0.d0
    DO orbit1=1, nHoles
       temp = temp + this%spBasis%spEnergies(orbit1)
    ENDDO
    this%eneKin = temp

    !     Hartree-Fock potential
    temp = 0.d0
    DO hhBlock=1, this%blocks%nHhBlocks
       blockIndex = this%blocks%hhBlocks(hhBlock)

       !write(*,*) 'hhBlock=',hhBlock
       vHhhhMatrix = this%hamilt%vHhhhRef(blockIndex)
       DO hh=1, vHhhhMatrix%dim1
          interaction = vHhhhMatrix%matr(hh, hh)
          
          temp = temp + interaction
          ! if (abs(interaction) > 0.00000001d0) then
          !    o1 = this%blocks%list(blockIndex)%hhPairs%matr(hh,1)
          !    o2 = this%blocks%list(blockIndex)%hhPairs%matr(hh,2)
             
          !    write(*,*) o1, o2, interaction
          ! endif
       ENDDO
    ENDDO
    !write(*,*) 'temp=',0.5d0*temp
    this%potRef = 0.5d0*temp
    !     Hartree-Fock energy
    this%eneRef = this%eneKin + this%potRef
    
    
  END SUBROUTINE calcRefEnergy
  
  !
  !     Routine to print energy values
  !
  SUBROUTINE printRefEnergy(this)
    CLASS(RefEnergy), INTENT(IN) :: this
    
    IF (this%spBasis%mpiJobs%iAm == 0) THEN
       
       OPEN(14, FILE='N_ene_ref')
       !     Energies given in Rydbergs (1 Ry = e^2/(2*a_{0}),
       !     a_{0} = hbar^2/(m*e^2))
       WRITE(*,*) ' '
       WRITE(*,*) 'eneKin/N = ', 2.d0*this%eneKin/this%spBasis%nOccupied, &
            ', potRef/N = ', 2.d0*this%potRef/this%spBasis%nOccupied
       WRITE(*,*) 'eneRef/N = ', 2.d0*this%eneRef/this%spBasis%nOccupied
       WRITE(*,*) ' '
       WRITE(14,'(I5,e14.5,e14.5,e14.5)') this%spBasis%nOccupied, &
            2.d0*this%eneKin/this%spBasis%nOccupied, &
            2.d0*this%potRef/this%spBasis%nOccupied, &
            2.d0*this%eneRef/this%spBasis%nOccupied
       CLOSE(14)
       
    ENDIF
    
    
  END SUBROUTINE printRefEnergy
  
  
END MODULE RefEnergyMod
