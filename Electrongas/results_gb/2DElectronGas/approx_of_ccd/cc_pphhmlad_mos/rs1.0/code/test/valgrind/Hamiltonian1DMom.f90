
!
!     Class for the Hamiltonian operator.
!
MODULE HamiltonianMod
  USE QuadratureMod
  USE ConstantsMod
  USE SpBasisMod
  IMPLICIT NONE
  
  TYPE Hamiltonian
     !SEQUENCE
     TYPE(SpBasis) :: basisSp           ! Single-particle basis
     CHARACTER(32) :: interaction       ! Two-body interaction 
     TYPE(Quadrature) :: rQuad          ! Coordinate r quadrature 
                                        !   in coordinate to
                                        !   momentum transformation
     REAL(DP) :: boxWidth               ! Box width L
     REAL(DP), ALLOCATABLE :: vMatrix(:,:,:,:) ! Interaction matrix (n1,n2,n3,n4)
   CONTAINS
     PROCEDURE :: setupInteractionMatr
     PROCEDURE :: getInteraction
     PROCEDURE :: setRQuad
     PROCEDURE :: gaussianRelMomDiscr
  END TYPE Hamiltonian
  
CONTAINS
  
  !
  !     Constructor for Hamiltonian object
  !
  FUNCTION Hamiltonian_(basis, interaction, width) RESULT(ham)
    TYPE(SpBasis), INTENT(IN) :: basis
    CHARACTER(*), INTENT(IN) :: interaction
    REAL(DP), INTENT(IN) :: width
    TYPE(Hamiltonian) :: ham
    INTEGER :: nOrbits
    
    ham%basisSp = basis  
    ham%interaction = interaction
    ham%boxWidth = width
    
    nOrbits = basis%nOccupied + basis%nUnoccupied
    ALLOCATE(ham%vMatrix(nOrbits, nOrbits, nOrbits, nOrbits))
    ham%vMatrix = 0.d0
    
    
  END FUNCTION Hamiltonian_
  
  !
  !     Routine to set up an interaction matrix
  !     in the given single-particle basis.
  !
  SUBROUTINE setupInteractionMatr(this)
    CLASS(Hamiltonian), INTENT(INOUT) :: this
    INTEGER :: nOrbits, orbit1, orbit2, orbit3, orbit4
    INTEGER :: n1, n2, n3, n4
    REAL(DP) :: interaction
    
    nOrbits = this%basisSp%nOccupied &
         + this%basisSp%nUnoccupied
    
    DO orbit1=1, nOrbits
       !     Quantum number
       n1 = this%basisSp%spQuantNr(orbit1)
       
       DO orbit2=orbit1, nOrbits  
          !     Quantum number
          n2 = this%basisSp%spQuantNr(orbit2)
          
          DO orbit3=1, nOrbits
             !     Quantum number
             n3 = this%basisSp%spQuantNr(orbit3)
             
             DO orbit4=orbit3, nOrbits
                !     Quantum number
                n4 = this%basisSp%spQuantNr(orbit4)
                
                IF ((n1 + n2) == (n3 + n4)) THEN
                   !     Calculate an antisymmetrized 
                   !     two-particle interaction matrix 
                   !     element
                   interaction = this%getInteraction(&
                        orbit1, orbit2, orbit3, orbit4)
                ELSE
                   interaction = 0.d0
                ENDIF
                !     Store the interaction matrix 
                !     element
                this%vMatrix(orbit1, orbit2, &
                     orbit3, orbit4) = interaction
                this%vMatrix(orbit2, orbit1, &
                     orbit3, orbit4) = -interaction
                this%vMatrix(orbit1, orbit2, &
                     orbit4, orbit3) = -interaction
                this%vMatrix(orbit2, orbit1, &
                     orbit4, orbit3) = interaction
                
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    
  END SUBROUTINE setupInteractionMatr
  
  !
  !     Calculate the interaction matrix element for a  
  !     given set of single-particle states
  !
  FUNCTION getInteraction(this, orbit1, orbit2, &
       orbit3, orbit4) RESULT(interaction)
    INTEGER, INTENT(IN) :: orbit1, orbit2, orbit3, orbit4
    CLASS(Hamiltonian), INTENT(IN) :: this
    REAL(DP) :: k1, k2, k3, k4, krelBra, krelKet
    REAL(DP) :: interaction
    
    !     Laboratory frame momenta
    k1 = this%basisSp%spStates(orbit1)
    k2 = this%basisSp%spStates(orbit2)
    k3 = this%basisSp%spStates(orbit3)
    k4 = this%basisSp%spStates(orbit4)
    
    !     Relative momentum of bra states
    krelBra = 0.5d0*(k1-k2)   
    !     Relative momentum of ket states
    krelKet = 0.5d0*(k3-k4)
    
    interaction = REAL(this%gaussianRelMomDiscr(krelBra, &
         krelKet),DP)
    
    
  END FUNCTION getInteraction
  
  !
  !     Set up coordinate quadrature 
  !
  SUBROUTINE setRQuad(this) 
    CLASS(Hamiltonian), INTENT(INOUT) :: this
    INTEGER :: nr
    TYPE(Quadrature) :: rQuad
    
    nr = INT(150*this%boxWidth/7.d0)
    rQuad = gaussLegendre(-this%boxWidth,this%boxWidth,nr)
    this%rQuad = rQuad
    
    
  END SUBROUTINE setRQuad

  !
  !     A simple Gaussian potential, as used in C. 
  !     Alexandrou et al., Phys. Rev. C, Vol. 39, No. 3,
  !     1076 (1989).
  !
  FUNCTION gaussian(rij) RESULT(potential)
    REAL(DP), INTENT(IN) :: rij
    REAL(DP), PARAMETER :: V1=12.d0,V2=-12.d0
    REAL(DP), PARAMETER :: sigma1=0.2d0,sigma2=0.8d0
    REAL(DP) :: potential
    
    potential = V1*EXP(-rij**2/(sigma1**2))/(sigma1*SQRT(pi)) &
         + V2*EXP(-rij**2/(sigma2**2))/(sigma2*SQRT(pi))
    
    
  END FUNCTION gaussian
  
  !
  !     The Gaussian potential in relative momentum
  !     coordinates. The function gives antisymmetrized
  !     matrix elements.
  !
  FUNCTION gaussianRelMom(k,kp) RESULT(potential)
    USE QuadratureMod
    REAL(DP), INTENT(IN) :: k,kp
    COMPLEX(DP) :: potential
    TYPE(Quadrature) :: rQuad
    COMPLEX(DP) :: tempDir,tempEx
    COMPLEX(DP), PARAMETER :: imag=(0.d0,1.d0)
    REAL(DP) :: rMax,r,wr,potentialR
    INTEGER :: nr,ir
    
    rMax = 5.d0
    nr = 150
    rQuad = gaussLegendre(-rMax,rMax,nr)
    
    tempDir = (0.d0,0.d0)
    tempEx = (0.d0,0.d0)
    DO ir=1, nr
       r = rQuad%xpoints(ir)
       wr = rQuad%weightx(ir)
       potentialR = gaussian(r)
       tempDir = tempDir + wr*potentialR*EXP(imag*(kp-k)*r)
       tempEx = tempEx + wr*potentialR*EXP(imag*(-kp-k)*r)
    ENDDO
    potential = tempDir - tempEx
    
    
  END FUNCTION gaussianRelMom
  
  !
  !     The Gaussian potential in relative momentum
  !     coordinates. The function gives antisymmetrized
  !     matrix elements.
  !
  FUNCTION gaussianRelMomDiscr(this, k, kp) RESULT(potential)
    USE QuadratureMod
    REAL(DP), INTENT(IN) :: k, kp
    CLASS(Hamiltonian), INTENT(IN) :: this
    COMPLEX(DP) :: potential
    TYPE(Quadrature) :: rQuad
    COMPLEX(DP) :: tempDir, tempEx
    COMPLEX(DP), PARAMETER :: imag=(0.d0,1.d0)
    REAL(DP) :: r, wr, potentialR
    INTEGER :: ir
    
    tempDir = (0.d0,0.d0)
    tempEx = (0.d0,0.d0)
    DO ir=1, this%rQuad%n
       r = this%rQuad%xpoints(ir)
       wr = this%rQuad%weightx(ir)
       potentialR = gaussian(r)
       tempDir = tempDir + wr*potentialR*EXP(imag*(kp-k)*r)
       tempEx = tempEx + wr*potentialR*EXP(imag*(-kp-k)*r)
    ENDDO
    potential = (tempDir - tempEx)/this%boxWidth
    
    
  END FUNCTION gaussianRelMomDiscr
  
  
END MODULE HamiltonianMod
