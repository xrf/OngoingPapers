
!
!     Class for pairs of single-particle states.
!
!     This version is for a plane wave basis
!     of a two-dimensional homogeneous system
!     with periodic boundary conditions
!     (for example the two-dimensional electron
!     gas)
!
MODULE SpPairsPW2dMod
  USE ConstantsMod
  USE SpPairsMod
  IMPLICIT NONE
  
  TYPE, EXTENDS(SpPairs) :: SpPairsPW2d
     !
     !     quantToBlock:      Gives block number for a set
     !                        of quantum numbers (K_CM, M_S)
     !     quantToBlockRe:    Gives block number for a set
     !                        of relative (recoupled) quantum 
     !                        numbers (k_rel, m_s_rel)
     !
     INTEGER, ALLOCATABLE  :: quantToBlockHh(:,:,:)
     INTEGER, ALLOCATABLE :: quantToBlockHhRe(:,:,:)
   CONTAINS
     PROCEDURE :: setupSpPairs => setupSpPairsPW2d
  END TYPE SpPairsPW2d
  
  
CONTAINS
  
  FUNCTION SpPairsPw2d_(basis) RESULT(this)
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    TYPE(SpPairsPW2d), POINTER :: this
    INTEGER :: nHoles, nParticles
   
    ALLOCATE(this)
    
    this%nBlocks = 0
    this%basisSp => basis
    
    nHoles = this%basisSp%nOccupied
    nParticles = this%basisSp%nUnoccupied
    this%hhPairs = IntMatList_(nHoles**2)
    this%ppPairs = IntMatList_(nParticles**2)
    this%phPairs = IntMatList_(nHoles*nParticles)
    
    ALLOCATE(this%blockToQuantHh(nHoles**2, 3))
    ALLOCATE(this%blockToQuantPp(nParticles**2, 3))
    ALLOCATE(this%blockToQuantPh(nHoles*nParticles, 3))
    ALLOCATE(this%blockToQuantHhRe(nHoles**2, 3))
    ALLOCATE(this%blockToQuantPpRe(nParticles**2, 3))
    ALLOCATE(this%blockToQuantPhRe(nHoles*nParticles, 3))
    
    
  END FUNCTION SpPairsPw2d_
  
  !
  !     Set up pairs of single-particle states and
  !     block indices.
  !
  SUBROUTINE SetupHhPairsPW2d(this, basis)
    CLASS(SpPairsPW2d), TARGET, INTENT(INOUT) :: this
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    INTEGER :: nHoles, hole1, hole2, pair, nPairs
    INTEGER, ALLOCATABLE :: tempBlockQN(:,:)
    INTEGER, ALLOCATABLE :: tempPairs(:,:)
    INTEGER, ALLOCATABLE :: tempBlockQN2(:,:)
    INTEGER, ALLOCATABLE :: tempPairs2(:,:)
    INTEGER, ALLOCATABLE :: rank(:), order(:)
    INTEGER :: i, newPair, thisBlock, thisPair
    LOGICAL :: equal, isEqual
    
    nHoles = basis%nOccupied
    ALLOCATE(tempBlockQN(nHoles**2, 3))
    ALLOCATE(tempBlockQN2(nHoles**2, 3))
    ALLOCATE(tempPairs(nHoles**2, 2))
    ALLOCATE(tempPairs2(nHoles**2, 2))
    ALLOCATE(rank(nHoles**2), order(nHoles**2))    
    
    !
    !      Get hole-hole pairs
    !
    pair = 1
    DO hole1=1, nHoles
       DO hole2=1, nHoles
          
          !     Store quantum numbers for the
          !     blocks (CM momentum, total spin)
          tempBlockQN(pair,:) = basis%spQuantNr(hole1,:) &
               + basis%spQuantNr(hole2,:)
          !     Store single-particle orbits
          !     corresponding to a block
          tempPairs(pair,:) = (/hole1, hole2/)
          
          !     The variable rank is used as a norm 
          !     to order pair states
          rank(pair) = (2*basis%nMax + 1)**2 &
               *(2+tempBlockQN(pair, 3)) &
               + (2*basis%nMax + 1) &
               *(basis%nMax + 1 + tempBlockQN(pair, 2)) &
               + basis%nMax + 1 + tempBlockQN(pair, 1)
          
          pair = pair + 1
       ENDDO
    ENDDO
    nPairs = pair - 1
    
    order = (/(i, i=1, nHoles**2)/)
    !     Sort the block states using the rank
    !     array. The array order gives the new 
    !     order of the array rank.
    CALL quicksort2(rank, order, nPairs)
    
    !     Reorder the pairs
    DO pair=1, nPairs
       newPair = order(pair)
       tempBlockQN2(pair,:) = tempBlockQN(newPair,:)
       tempPairs2(pair,:) = tempPairs(newPair,:)
    ENDDO
    
    !     Count and allocate blocks
    thisBlock = 0
    thisPair = 0
    DO pair=2, nPairs
       
       thisPair = thisPair + 1
       equal = isEqual(tempBlockQN2(pair, :), &
            tempBlockQN2(pair-1, :), 3)
       IF (.NOT.equal) THEN
          !     Allocate the previous block
          this%hhPairs%list(thisBlock) &
               = IntMatrix_(thisBlock, 2)
          
          thisBlock = thisBlock + 1
          thisPair = 1
          this%blockToQuantHh(thisBlock, :) = &
               tempBlockQN2(pair, :)
          !this%quantToBlock()
       ENDIF
       
       
       this%hhPairs%list(thisBlock)%matrix(thisPair, :) = &
            tempPairs2(pair, :)
    ENDDO
    this%nBlocks = this%nBlocks + thisBlock
    
    !     Assign pairs to blocks
    DO thisBlock=1, this%nBlocks
       
       
    DEALLOCATE(tempBlockQN, tempBlockQN2)
    DEALLOCATE(tempPairs, tempPairs2)
    DEALLOCATE(rank, order)
    
    
  END SUBROUTINE SetupHhPairsPW2d
  
  
END MODULE SpPairsPW2dMod
