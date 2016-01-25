
!
!     Class for list of all recoupled blocks (k_rel, m_s_rel).
!
MODULE BlockListReMod
  USE BlockListMod
  USE PairListMod
  USE SpBasisMod
  USE Block1Mod
  IMPLICIT NONE
  
  TYPE :: BlockListRe
     !
     !     nPhBlocks:    Number of ph blocks
     !     nQNumbers:    Numbers of quantum numbers
     !     phPairs:      Ordered list of all particle-hole
     !                   pairs
     !     phList(:):    List of blocks with recoupled ph
     !                   states
     !     pphhList(:):  List of blocks with recoupled pp
     !                   and hh states
     !
     TYPE(BlockList) :: blocks
     CLASS(SpBasis), POINTER :: basis
     INTEGER :: nPhBlocks, nPphhBlocks, nQNumbers
     TYPE(Block1), ALLOCATABLE :: phList(:)
     TYPE(Block1), ALLOCATABLE :: pphhList(:)
     TYPE(PairList) :: phPairs
   CONTAINS
     PROCEDURE :: setupPhBlockList
     PROCEDURE :: setupPhBlocks
     PROCEDURE :: countPhBlocks
     PROCEDURE :: setupPhPairsRe
  END TYPE BlockListRe
  
CONTAINS
  
  !
  !     Constructor for BlockListRe object.
  !
  FUNCTION BlockListRe_(nQnumbers, basis, blocks) RESULT(this)
    INTEGER, INTENT(IN) :: nQNumbers
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    TYPE(BlockList), INTENT(IN) :: blocks
    TYPE(BlockListRe) :: this
    
    this%basis => basis
    this%blocks = blocks
    this%nQNumbers = nQnumbers
    
    
  END FUNCTION BlockListRe_
  
  SUBROUTINE setupPhBlockList(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    
    IF (this%basis%nUnoccupied > 0) THEN
       CALL setupPhPairsRe(this)
       
       WRITE(*,*) ' '
       WRITE(*,*) ' Counting number of recoupled blocks...'
       WRITE(*,*) ' '
       CALL countPhBlocks(this)
       WRITE(*,*) ' '
       WRITE(*,*) ' Setting up (k_rel, m_s_rel) ph blocks...'
       WRITE(*,*) ' '
       CALL setupPhBlocks(this)
       
    ENDIF
    
    
  END SUBROUTINE setupPhBlockList
  
  !
  !     Set up list of (k_rel, m_s_rel) ph blocks
  !     obtained from a pp and a hh pair
  SUBROUTINE setupPhBlocks(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    TYPE(PairList) :: phPairs
    TYPE(IntMatrix) :: pairsMatr
    INTEGER :: nPairsPh, ph, thisBlock, thisRank, nOcc
    INTEGER :: phThis
    LOGICAL :: occupied
    
    phPairs = this%phPairs
    nPairsPh = phPairs%nPairs

    thisRank = phPairs%rank(1)
    
    ph = 0
    thisBlock = 0
    DO WHILE (thisBlock < this%nPhBlocks)
       occupied = .FALSE.
       
       nOcc = 0
       DO WHILE ((ph + 1 <= nPairsPh) &
            .AND.(phPairs%rank(ph+1) <= thisRank))
          ph = ph + 1
          nOcc = nOcc + 1
          occupied = .TRUE.
       ENDDO
       !     If this block is occupied by particle-hole
       !     states, store them
       IF (occupied) THEN
          thisBlock = thisBLock + 1
          
          this%phList(thisBlock) = Block1_(this%nQNumbers)
          this%phList(thisBlock)%blockQNrs(:) &
               = phPairs%blockQNumbers(ph,:)
          this%phList(thisBlock)%occupiedPh = .TRUE.
          
          WRITE(*,*) ' '
          WRITE(*,*) ' Block nr',thisBlock
          pairsMatr = IntMatrix_(nOcc, 2)
          DO phThis=1, nOcc
             pairsMatr%matr(phThis,:) &
                  = phPairs%pairs(ph-nOcc+phThis,:)
             !WRITE(*,*) phPairs%pairs(ph-nOcc+phThis,:), &
             !     phPairs%rank(ph-nOcc+phThis)
             WRITE(*,*) phPairs%blockQNumbers(ph-nOcc+phThis,:), &
                  phPairs%rank(ph-nOcc+phThis)
          ENDDO
          this%phList(thisBlock)%phPairs = pairsMatr
       ENDIF
       
       thisRank = thisRank + 1
       
    ENDDO
    WRITE(*,*) ' '
    
    
  END SUBROUTINE setupPhBlocks
  
  !
  !     Count the total number of ph blocks obtained
  !     fram a pp and a hh pair
  !
  SUBROUTINE countPhBlocks(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    TYPE(PairList) :: phPairs
    INTEGER :: nPairsPh, ph, blocks, thisRank
    LOGICAL :: occupied
    
    phPairs = this%phPairs
    nPairsPh = phPairs%nPairs
    write(*,*) 'nPairsPh=',nPairsPh

    !     Count number of blocks
    thisRank = phPairs%rank(1)
    
    ph = 0
    blocks = 0
    DO WHILE (ph < nPairsPh)
       occupied = .FALSE.
       !write(*,*) 'ph=',ph
       DO WHILE ((ph + 1 <= nPairsPh) &
            .AND. (phPairs%rank(ph+1) <= thisRank))
          ph = ph + 1
          occupied = .TRUE.
       ENDDO
       
       IF (occupied) blocks = blocks + 1
       
       thisRank = thisRank + 1
    ENDDO
    
    this%nPhBlocks = blocks
    ALLOCATE(this%phList(blocks))
    
    
  END SUBROUTINE countPhBlocks
  
  !
  !     Set up list of particle-hole pairs using
  !     the first states from a pp and a hh pair
  !  
  SUBROUTINE setupPhPairsRe(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    CLASS(SpBasis), POINTER :: basis
    TYPE(Block1) :: thisBlock
    INTEGER, ALLOCATABLE :: blockQNumbers(:,:)
    INTEGER, ALLOCATABLE :: blockQNumbersNew(:,:)
    INTEGER, ALLOCATABLE :: pairs(:,:)
    INTEGER, ALLOCATABLE :: pairsNew(:,:)
    INTEGER, ALLOCATABLE :: rank(:), rankNew(:)
    INTEGER, ALLOCATABLE :: blockInd(:), blockIndNew(:)
    INTEGER, ALLOCATABLE :: order(:)
    INTEGER :: pair, pphhBlock, blockIndex, nPphhBlocks
    INTEGER :: particle, hole, nPphh, nPp, nHh, nPairs
    INTEGER :: i, pp, hh, newPair
    
    nPphhBlocks = this%blocks%nPphhBlocks
    basis => this%basis
    
    !     Count number of pairs
    pair = 0
    DO pphhBlock=1, nPphhBlocks
       blockIndex = this%blocks%pphhBlocks(pphhBlock)
       
       thisBlock = this%blocks%list(blockIndex)
       pair = pair + thisBlock%ppPairs%dim1 &
            *thisBlock%hhPairs%dim1
    ENDDO
    nPairs = pair
    
    this%phPairs = PairList_(nPairs, this%nQnumbers)
    
    ALLOCATE(blockQNumbers(nPairs, this%nQNumbers))
    ALLOCATE(blockQNumbersNew(nPairs, this%nQNumbers))
    ALLOCATE(pairs(nPairs, 2))
    ALLOCATE(pairsNew(nPairs, 2))
    ALLOCATE(rank(nPairs), rankNew(nPairs))
    ALLOCATE(blockInd(nPairs), blockIndNew(nPairs))
    ALLOCATE(order(nPairs))
    
    !     Find pairs
    pair = 0
    DO pphhBlock=1, nPphhBlocks
       blockIndex = this%blocks%pphhBlocks(pphhBlock)
       
       thisBlock = this%blocks%list(blockIndex)
       DO pp=1, thisBlock%ppPairs%dim1 
          DO hh=1, thisBlock%hhPairs%dim1
             pair = pair + 1
             
             particle = thisBlock%ppPairs%matr(pp, 1)
             hole = thisBlock%hhPairs%matr(hh, 1)
             
             !     Store quantum numbers for the
             !     blocks (CM momentum, total spin)
             blockQNumbers(pair,:) = &
                  basis%spQuantNr(particle,:) &
                  - basis%spQuantNr(hole,:)
             
             !     Store single-particle orbits
             !     corresponding to a block
             pairs(pair,:) = (/particle, hole/)
            
             blockInd(pair) = blockIndex 
             
             ! !     The variable rank is used as a norm 
             ! !     to order pair states
             ! IF (this%nQNumbers == 3) THEN
             !    rank(pair) = (4*basis%nMax + 1)**2 &
             !         *(2 + blockQNumbers(pair, 3)) &
             !         + (4*basis%nMax + 1) &
             !         *(2*basis%nMax + blockQNumbers(pair, 2)) &
             !         + 2*basis%nMax + 1 + blockQNumbers(pair, 1)
             ! ELSEIF (this%nQNumbers == 4) THEN
             !    rank(pair) = (4*basis%nMax + 1)**3 &
             !         *(2 + blockQNumbers(pair, 4)) &
             !         + (4*basis%nMax + 1)**2 &
             !         *(2*basis%nMax + blockQNumbers(pair, 3)) &
             !         + (4*basis%nMax + 1) &
             !         *(2*basis%nMax + blockQNumbers(pair, 2)) &
             !         + 2*basis%nMax + 1 + blockQNumbers(pair, 1)
             ! ELSE
             !    WRITE(*,*) 'rank setup only for 2d and 3d'
             !    STOP
             ! ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    rank = blockInd
    order = (/(i, i=1, nPairs)/)
    
    ! !    Sort the pair states using the rank
    ! !    array. The array order gives the new 
    ! !    order of the array rank.
    ! CALL qsorti(order, nPairs, rank)
    
    WRITE(*,*) ' Particle-hole pairs:'
    WRITE(*,*) ' '
    WRITE(*,*) ' particle   hole   ranking number'
    !     Reorder the pairs
    DO pair=1, nPairs
       newPair = order(pair)
       blockQNumbersNew(pair,:) = blockQNumbers(newPair,:)
       pairsNew(pair,:) = pairs(newPair,:)
       rankNew(pair) = rank(newPair)
       blockIndNew(pair) = blockInd(newPair)
       !write(*,'(4(x, I6))')  blockQNumbers(newPair,:), rank(newPair)
       write(*,'(4(x, I6))') pairs(newPair,:), rank(newPair), &
            blockInd(newPair)
    ENDDO
    this%phPairs%blockQNumbers = blockQNumbersNew
    this%phPairs%pairs = pairsNew
    this%phPairs%rank = rankNew
    this%phPairs%blockInd = blockIndNew
    write(*,*) ' '
    
    DEALLOCATE(blockQNumbers, blockQNumbersNew)
    DEALLOCATE(pairs, pairsNew)
    DEALLOCATE(rank, rankNew, order)
    DEALLOCATE(blockInd, blockIndNew)
    
    
  END SUBROUTINE setupPhPairsRe
  
  
  
END MODULE BlockListReMod
