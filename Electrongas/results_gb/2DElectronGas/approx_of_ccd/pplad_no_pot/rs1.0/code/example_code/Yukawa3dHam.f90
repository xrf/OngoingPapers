
!
!     Class for Hamiltonian operator with Yukawa interaction
!     for a three-dimensional periodic system.
!
MODULE Yukawa3dHamMod
  USE HamiltonianMod
  USE ConstantsMod
  USE SpBasisMod
  IMPLICIT NONE
  
  TYPE, EXTENDS(Hamiltonian) :: Yukawa3dHam
     PRIVATE
     REAL(DP) :: mu, rs, area, volume
   CONTAINS
     PROCEDURE :: getInteraction => getYukawa3d
  END TYPE Yukawa3dHam
  
CONTAINS
  
  !
  !     Constructor for Yukawa3dHam object
  !
  FUNCTION Yukawa3dHam_(basis, blocks, blocksRe, &
       mu, rs, nParticles) RESULT(this)
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    TYPE(BlockList), TARGET, INTENT(IN) :: blocks
    TYPE(BlockListRe), TARGET, INTENT(IN) :: blocksRe
    REAL(DP), INTENT(IN) :: mu, rs
    INTEGER, INTENT(IN) :: nParticles
    TYPE(Yukawa3dHam), POINTER :: this
    INTEGER :: nOrbits, nPphhBlocks, nBlocksRe, nBlocks
    
    ALLOCATE(this)
    
    this%basisSp => basis
    this%blocks => blocks
    this%blocksRe => blocksRe
    nPphhBlocks = blocks%nPphhBlocksThis
    nBlocksRe = blocksRe%nBlocksReThis
    nBlocks = blocks%nBlocks
    
    ALLOCATE(this%vHhhh(nPphhBlocks))
    ALLOCATE(this%vHhhhRef(nBlocks))
    ALLOCATE(this%vHhpp(nPphhBlocks))
    ALLOCATE(this%vPppp(nPphhBlocks))
    ALLOCATE(this%vPphh(nPphhBlocks))
    ALLOCATE(this%vPphhRe(nBlocksRe))
    ALLOCATE(this%vPhphRe(nBlocksRe))
    
    nOrbits = basis%nOccupied + basis%nUnoccupied
    this%mu = mu
    this%rs = rs
    this%volume = nParticles*4.d0*pi*rs**3/3.d0
    this%area = (this%volume)**(2.d0/3.d0)
    
    
  END FUNCTION Yukawa3dHam_
    
  !
  !     Calculate the antisymmetrized interaction 
  !     matrix element for a given set of 
  !     single-particle states.
  !
  FUNCTION getYukawa3d(this, orbit1, orbit2, &
       orbit3, orbit4) RESULT(interaction)
    CLASS(Yukawa3dHam), TARGET, INTENT(IN) :: this
    INTEGER, INTENT(IN) :: orbit1, orbit2, orbit3, orbit4
    INTEGER :: nx1, ny1, nz1, ms1, nx2, ny2, nz2, ms2
    INTEGER :: nx3, ny3, nz3, ms3, nx4, ny4, nz4, ms4
    REAL(DP) :: interaction
    REAL(DP) :: denom2, denom
    
    !     Quantum numbers of orbit 1
    nx1 = this%basisSp%spQuantNr(orbit1, 1)
    ny1 = this%basisSp%spQuantNr(orbit1, 2)
    nz1 = this%basisSp%spQuantNr(orbit1, 3)
    ms1 = this%basisSp%spQuantNr(orbit1, 4)
    
    !     Quantum numbers of orbit 2
    nx2 = this%basisSp%spQuantNr(orbit2, 1)
    ny2 = this%basisSp%spQuantNr(orbit2, 2)
    nz2 = this%basisSp%spQuantNr(orbit2, 3)
    ms2 = this%basisSp%spQuantNr(orbit2, 4)
    
    !     Quantum numbers of orbit 3
    nx3 = this%basisSp%spQuantNr(orbit3, 1)
    ny3 = this%basisSp%spQuantNr(orbit3, 2)
    nz3 = this%basisSp%spQuantNr(orbit3, 3)
    ms3 = this%basisSp%spQuantNr(orbit3, 4)
    
    !     Quantum numbers of orbit 4
    nx4 = this%basisSp%spQuantNr(orbit4, 1)
    ny4 = this%basisSp%spQuantNr(orbit4, 2)
    nz4 = this%basisSp%spQuantNr(orbit4, 3)
    ms4 = this%basisSp%spQuantNr(orbit4, 4)
    
    interaction = 0.d0
    !    IF ((nx1 + nx2 == nx3 + nx4).AND.&
    !         (ny1 + ny2 == ny3 + ny4).AND.&
    !         (nz1 + nz2 == nz3 + nz4)) THEN
    IF ((ms1 == ms3).AND.(ms2 == ms4)&
         .AND.((nx1 /= nx3).OR.(ny1 /= ny3)&
         .OR.(nz1 /= nz3))) THEN
       denom = this%mu**2 &
            + 4.d0*pi**2*((nx3-nx1)**2 + (ny3-ny1)**2&
            + (nz3-nz1)**2)/this%area

       interaction = interaction + 1.d0/denom
    ENDIF

    IF ((ms1 == ms4).AND.(ms2 == ms3)&
         .AND.((nx1 /= nx4).OR.(ny1 /= ny4)&
         .OR.(nz1 /= nz4))) THEN
       denom = this%mu**2 &
            + 4.d0*pi**2*((nx4-nx1)**2 + (ny4-ny1)**2 &
            + (nz4-nz1)**2)&
            /this%area
       
       interaction = interaction - 1.d0/denom
    ENDIF
    
    interaction = interaction*4.d0*pi&
         /this%volume
    !   ENDIF
    
  END FUNCTION getYukawa3d
  
  
END MODULE Yukawa3dHamMod
