
!
!     Class for Hamiltonian operator with Yukawa interaction
!     for a two-dimensional periodic system.
!
MODULE Yukawa2dHamMod
  USE HamiltonianMod
  USE ConstantsMod
  USE BlockListMod
  USE SpBasisMod
  IMPLICIT NONE
  
  TYPE, EXTENDS(Hamiltonian) :: Yukawa2dHam
     PRIVATE
     REAL(DP) :: mu, rs, area
   CONTAINS 
     PROCEDURE :: getInteraction => getYukawa2d
  END TYPE Yukawa2dHam
  
CONTAINS
  
  !
  !     Constructor for Yukawa2dHam object
  !
  FUNCTION Yukawa2dHam_(basis, blocks, &
       blocksRe, mu, rs, nParticles) RESULT(this)
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    TYPE(BlockList), INTENT(IN) :: blocks
    TYPE(BlockListRe), INTENT(IN) :: blocksRe
    REAL(DP), INTENT(IN) :: mu, rs
    INTEGER, INTENT(IN) :: nParticles
    TYPE(Yukawa2dHam), POINTER :: this
    INTEGER :: nOrbits, nBlocks
    
    ALLOCATE(this)
    this%basisSp => basis  
    this%blocks = blocks
    this%blocksRe = blocksRe
    nBlocks = blocks%nBlocks
    
    ALLOCATE(this%vHhhh(nBlocks))
    ALLOCATE(this%vHhhp(nBlocks))
    ALLOCATE(this%vHhpp(nBlocks))
    ALLOCATE(this%vHppp(nBlocks))
    ALLOCATE(this%vPppp(nBlocks))
    ALLOCATE(this%vPhph(nBlocks))
    ALLOCATE(this%vPphh(nBlocks))
    ALLOCATE(this%vPphhRe(nBlocks))
    ALLOCATE(this%vPhphRe(nBlocks))

    nOrbits = basis%nOccupied + basis%nUnoccupied
    this%mu = mu
    this%rs = rs
    this%area = nParticles*pi*rs**2
    
    
  END FUNCTION Yukawa2dHam_
  
  !
  !     Calculate the antisymmetrized interaction 
  !     matrix element for a given set of 
  !     single-particle states.
  !
  FUNCTION getYukawa2d(this, orbit1, orbit2, &
       orbit3, orbit4) RESULT(interaction)
    CLASS(Yukawa2dHam), TARGET, INTENT(IN) :: this
    INTEGER, INTENT(IN) :: orbit1, orbit2, orbit3, orbit4
    INTEGER :: nx1, ny1, ms1, nx2, ny2, ms2
    INTEGER :: nx3, ny3, ms3, nx4, ny4, ms4
    REAL(DP) :: interaction
    REAL(DP) :: denom2, denom
    
    !     Quantum numbers of orbit 1
    nx1 = this%basisSp%spQuantNr(orbit1, 1)
    ny1 = this%basisSp%spQuantNr(orbit1, 2)
    ms1 = this%basisSp%spQuantNr(orbit1, 3)
    
    !     Quantum numbers of orbit 2
    nx2 = this%basisSp%spQuantNr(orbit2, 1)
    ny2 = this%basisSp%spQuantNr(orbit2, 2)
    ms2 = this%basisSp%spQuantNr(orbit2, 3)
    
    !     Quantum numbers of orbit 3
    nx3 = this%basisSp%spQuantNr(orbit3, 1)
    ny3 = this%basisSp%spQuantNr(orbit3, 2)
    ms3 = this%basisSp%spQuantNr(orbit3, 3)

    !     Quantum numbers of orbit 4
    nx4 = this%basisSp%spQuantNr(orbit4, 1)
    ny4 = this%basisSp%spQuantNr(orbit4, 2)
    ms4 = this%basisSp%spQuantNr(orbit4, 3)
    
    interaction = 0.d0

    IF ((ms1 == ms3).AND.(ms2 == ms4)&
         .AND.((nx1 /= nx3).OR.(ny1 /= ny3))) THEN
       denom2 = this%mu**2 &
            + 4.d0*pi**2*((nx3-nx1)**2 + (ny3-ny1)**2)&
            /this%area
       
       denom = SQRT(denom2)
       interaction = interaction + 1.d0/denom
    ENDIF
    
    IF ((ms1 == ms4).AND.(ms2 == ms3)&
         .AND.((nx1 /= nx4).OR.(ny1 /= ny4))) THEN
       denom2 = this%mu**2 &
            + 4.d0*pi**2*((nx4-nx1)**2 + (ny4-ny1)**2)&
            /this%area
       
       denom = SQRT(denom2)
       interaction = interaction - 1.d0/denom
    ENDIF
    interaction = interaction*2.d0*pi&
         /this%area
    
    
  END FUNCTION getYukawa2d
  
  
END MODULE Yukawa2dHamMod
