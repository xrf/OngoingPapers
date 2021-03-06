
!
!     Class for list of blocks where the single-particle
!     states are recoupled.
!
MODULE BlockListReMod
  USE IntMatListMod
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
     INTEGER :: nBlocks, nQNumbers
     TYPE(Block1), ALLOCATABLE :: list(:)
     TYPE(PairList) :: phPairs
     TYPE(PairList) :: hpPairs
     TYPE(IntMatList) :: pphh2ReHpph
     TYPE(IntMatList) :: hpphRe2Pphh
     TYPE(IntMatList) :: pphh2Index
     TYPE(IntMatList) :: hpph2Index
     REAL(DP), ALLOCATABLE :: index2ReHpph(:,:)
     REAL(DP), ALLOCATABLE :: index2Pphh(:,:)
     INTEGER, ALLOCATABLE :: rank2Block(:)
     INTEGER :: maxRank
   CONTAINS
     PROCEDURE :: setupPhBlockList
     PROCEDURE :: setupBlocksRe
     PROCEDURE :: setupPhPairsRe
     PROCEDURE :: countBlocksRe
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
    INTEGER :: nHoles, nParticles, nPhPairs
    
    this%basis => basis
    this%blocks = blocks
    this%nQNumbers = nQNumbers
    
    nHoles = basis%nOccupied
    nParticles = basis%nUnoccupied
    nPhPairs = nHoles*nParticles
    this%phPairs = PairList_(nPhPairs, nQNumbers)
    this%hpPairs = PairList_(nPhPairs, nQNumbers)
    
    
  END FUNCTION BlockListRe_
  
  
  SUBROUTINE setupPhBlockList(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    
    IF (this%basis%nUnoccupied > 0) THEN
       
       WRITE(*,*) ' '
       WRITE(*,*) 'Setting up recoupled blocks...'
       WRITE(*,*) ' '
       !     Set up pairs of ph and hp single-particle states
       CALL setupPhPairsRe(this)
       CALL setupHpPairsRe(this)
       !     Count number of (k_rel, m_s_rel) blocks
       CALL countBlocksRe(this)
       !     Set up list of blocks
       CALL setupBlocksRe(this)

       CALL setupHpph2Index(this)
       CALL setupPphh2Index(this)
       CALL setupHpphRe2Pphh(this)
       CALL setupPphh2ReHpph(this)
    ENDIF
    
    
  END SUBROUTINE setupPhBlockList
  
  !
  !     Set up list of (k_rel, m_s_rel) blocks
  !
  SUBROUTINE setupBlocksRe(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    TYPE(PairList) :: phPairs, hpPairs
    TYPE(IntMatrix) :: pairsMatrPh, pairsMatrHp
    INTEGER :: nPairsPh, nPairsHp, ph, hp, thisRank, thisBlock
    INTEGER :: nOccPh, nOccHp, phThis, hpThis
    LOGICAL :: occupied, occupiedPh, occupiedHp
    
    phPairs = this%phPairs
    hpPairs = this%hpPairs
    nPairsPh = phPairs%nPairs
    nPairsHp = hpPairs%nPairs
    
    thisRank = MIN(phPairs%rank(1), &
         hpPairs%rank(1))
    
    ph = 0
    hp = 0
    thisBlock = 0
    DO WHILE (thisBlock < this%nBlocks)
       occupiedPh = .FALSE.
       occupiedHp = .FALSE.
       occupied = .FALSE.
       
       nOccPh = 0
       DO WHILE ((ph + 1 <= nPairsPh) &
            .AND. (phPairs%rank(ph+1) <= thisRank))
          ph = ph + 1
          nOccPh = nOccPh + 1
          occupiedPh = .TRUE.
       ENDDO
       
       nOccHp = 0
       DO WHILE ((hp + 1 <= nPairsHp) &
            .AND. (hpPairs%rank(hp+1) <= thisRank))
          hp = hp + 1
          nOccHp = nOccHp + 1
          occupiedHp = .TRUE.
       ENDDO
       
       occupied = (occupiedPh .AND. occupiedHp)
       
       !     If this block is occupied by particle-hole
       !     states, store them
       IF (occupied) THEN
          thisBlock = thisBlock + 1
          this%rank2Block(thisRank) = thisBlock
          
          this%list(thisBlock) = Block1_(this%nQNumbers)
          this%list(thisBlock)%blockQNrs(:) &
               = phPairs%blockQNumbers(ph, :)
          this%list(thisBlock)%rank = thisRank
          
          WRITE(*,*) ' '
          WRITE(*,*) ' Block nr',thisBlock
          pairsMatrPh = IntMatrix_(nOccPh, 2)
          
          DO phThis=1, nOccPh
             pairsMatrPh%matr(phThis,:) &
                  = phPairs%pairs(ph-nOccPh+phThis,:)
             
             WRITE(*,*) phPairs%pairs(ph-nOccPh+phThis,:), &
                  phPairs%rank(ph-nOccPh+phThis)
             !WRITE(*,*) phPairs%blockQNumbers(ph-nOcc+phThis,:), &
             !    phPairs%rank(ph-nOcc+phThis)
          ENDDO
          this%list(thisBlock)%phPairs = pairsMatrPh
          CALL IntMatrix_d(pairsMatrPh)
          
          pairsMatrHp = IntMatrix_(nOccHp, 2)
          
          WRITE(*,*) ' '
          DO hpThis=1, nOccHp
             pairsMatrHp%matr(hpThis,:) &
                  = hpPairs%pairs(hp-nOccHp+hpThis,:)
             
             WRITE(*,*) hpPairs%pairs(hp-nOccHp+hpThis,:), &
                  hpPairs%rank(hp-nOccHp+hpThis)
          ENDDO
          this%list(thisBlock)%hpPairs = pairsMatrHp
          CALL IntMatrix_d(pairsMatrHp)
          
       ENDIF
       
       thisRank = thisRank + 1
       
    ENDDO
    WRITE(*,*) ' '
    
    
  END SUBROUTINE setupBlocksRe
  
  !
  !     Count the total number of blocks
  !
  SUBROUTINE countBlocksRe(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    TYPE(PairList) :: phPairs, hpPairs
    INTEGER :: nPairsPh, nPairsHp
    INTEGER :: ph, hp, blocks, pairsTot, nPairsTot
    INTEGER :: thisRank, maxRank
    LOGICAL :: occupiedPh, occupiedHp, occupied
    
    phPairs = this%phPairs
    hpPairs = this%hpPairs
    nPairsPh = phPairs%nPairs
    nPairsHp = hpPairs%nPairs
    nPairsTot = nPairsPh + nPairsHp
    WRITE(*,*) ' '
    WRITE(*,*) ' Number of ph pairs:', nPairsPh
    WRITE(*,*) ' Number of hp pairs:', nPairsHp
    
    maxRank = MAX(phPairs%rank(nPairsPh), &
         hpPairs%rank(nPairsHp))
    ALLOCATE(this%rank2Block(maxRank))
    this%rank2Block = 0
    this%maxRank = maxRank
    
    thisRank = MIN(phPairs%rank(1), &
         hpPairs%rank(1))
    
    ph = 0
    hp = 0
    blocks = 0
    pairsTot = 0
    DO WHILE (pairsTot < nPairsTot)
       occupied = .FALSE.
       occupiedPh = .FALSE.
       occupiedHp = .FALSE.
       
       
       DO WHILE ((ph + 1 <= nPairsPh) &
            .AND. (phPairs%rank(ph+1) <= thisRank))
          ph = ph + 1
          occupiedPh = .TRUE.
       ENDDO
       
       DO WHILE ((hp + 1 <= nPairsHp) &
            .AND. (hpPairs%rank(hp+1) <= thisRank))
          hp = hp + 1
          occupiedHp = .TRUE.
       ENDDO
          
       occupied = (occupiedPh .AND. occupiedHp)
       IF (occupied) THEN
          blocks = blocks + 1
          
          this%rank2block(thisRank) = blocks
       ENDIF
       
       pairsTot = ph + hp
       thisRank = thisRank + 1
    ENDDO
    
    !     Allocate block list
    this%nBlocks = blocks
    WRITE(*,*) ' Number of cross-coupled blocks: ', blocks
    ALLOCATE(this%list(this%nBlocks))
    this%hpphRe2Pphh = IntMatList_(this%nBlocks)
    this%pphh2ReHpph = IntMatList_(this%blocks%nPphhBlocks)
    
    
  END SUBROUTINE countBlocksRe
  
  
  SUBROUTINE setupPhPairsRe(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    CLASS(SpBasis), POINTER :: basis
    INTEGER, ALLOCATABLE :: blockQNumbers(:,:)
    INTEGER, ALLOCATABLE :: blockQNumbersNew(:,:)
    INTEGER, ALLOCATABLE :: pairs(:,:)
    INTEGER, ALLOCATABLE :: pairsNew(:,:)
    INTEGER, ALLOCATABLE :: rank(:), rankNew(:)
    INTEGER, ALLOCATABLE :: order(:)
    INTEGER :: nHoles, nParticles, nPairs
    INTEGER :: hole, particle, pair, i
    INTEGER :: newPair
    
    basis => this%basis
    nHoles = basis%nOccupied
    nParticles = basis%nUnoccupied
    nPairs = nHoles*nParticles
    ALLOCATE(blockQNumbers(nPairs, this%nQNumbers))
    ALLOCATE(blockQNumbersNew(nPairs, this%nQNumbers))
    ALLOCATE(pairs(nPairs, 2))
    ALLOCATE(pairsNew(nPairs, 2))
    ALLOCATE(rank(nPairs), rankNew(nPairs))
    ALLOCATE(order(nPairs))
    
    pair = 0
    DO hole=1, nHoles
       DO particle=nHoles+1, nHoles+nParticles
          pair = pair + 1
          
          !     Store quantum numbers for the
          !     blocks (k_rel, m_s_rel)
          blockQNumbers(pair,:) &
               = basis%spQuantNr(particle,:) &
               - basis%spQuantNr(hole,:)
          
          !     Store single-particle orbits
          !     corresponding to a block
          pairs(pair,:) = (/particle, hole/)
          
          !     The variable rank is used as a norm 
          !     to order pair states
          rank(pair) = getRank(basis, blockQNumbers(pair,:), &
               this%nQNumbers)
       ENDDO
    ENDDO
    
    order = (/(i, i=1, nPairs)/)
    !     Sort the pair states using the rank
    !     array. The array order gives the new 
    !     order of the array rank.
    CALL qsorti(order, nPairs, rank)

    WRITE(*,*) ' '
    WRITE(*,*) ' Particle-hole pairs:'
    WRITE(*,*) ' '
    WRITE(*,*) ' particle   hole   ranking number'
    !     Reorder the pairs
    DO pair=1, nPairs
       newPair = order(pair)
       blockQNumbersNew(pair,:) = blockQNumbers(newPair,:)
       pairsNew(pair,:) = pairs(newPair,:)
       rankNew(pair) = rank(newPair)
       !write(*,'(4(x, I6))')  blockQNumbers(newPair,:), rank(newPair)
       write(*,'(3(x, I6))')  pairs(newPair,:), rank(newPair)
    ENDDO
    
    this%phPairs%blockQNumbers = blockQNumbersNew
    this%phPairs%pairs = pairsNew
    this%phPairs%rank = rankNew
    write(*,*) ' '
    
    DEALLOCATE(blockQNumbers)
    DEALLOCATE(blockQNumbersNew)
    DEALLOCATE(pairs, pairsNew)
    DEALLOCATE(rank, rankNew)
    DEALLOCATE(order)
    
    
  END SUBROUTINE setupPhPairsRe
  
  
  SUBROUTINE setupHpPairsRe(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    CLASS(SpBasis), POINTER :: basis
    INTEGER, ALLOCATABLE :: blockQNumbers(:,:)
    INTEGER, ALLOCATABLE :: blockQNumbersNew(:,:)
    INTEGER, ALLOCATABLE :: pairs(:,:)
    INTEGER, ALLOCATABLE :: pairsNew(:,:)
    INTEGER, ALLOCATABLE :: rank(:), rankNew(:)
    INTEGER, ALLOCATABLE :: order(:)
    INTEGER :: nHoles, nParticles, nPairs
    INTEGER :: hole, particle, pair, i
    INTEGER :: newPair
    
    basis => this%basis
    nHoles = basis%nOccupied
    nParticles = basis%nUnoccupied
    nPairs = nHoles*nParticles
    ALLOCATE(blockQNumbers(nPairs, this%nQNumbers))
    ALLOCATE(blockQNumbersNew(nPairs, this%nQNumbers))
    ALLOCATE(pairs(nPairs, 2))
    ALLOCATE(pairsNew(nPairs, 2))
    ALLOCATE(rank(nPairs), rankNew(nPairs))
    ALLOCATE(order(nPairs))
    
    pair = 0
    DO hole=1, nHoles
       DO particle=nHoles+1, nHoles+nParticles
          pair = pair + 1
          
          !     Store quantum numbers for the
          !     blocks (k_rel, m_s_rel)
          blockQNumbers(pair,:) &
               = basis%spQuantNr(hole,:) &
               - basis%spQuantNr(particle,:)
          
          !     Store single-particle orbits
          !     corresponding to a block
          pairs(pair,:) = (/hole, particle/)
          
          !     The variable rank is used as a norm 
          !     to order pair states
          rank(pair) = getRank(basis, blockQNumbers(pair,:), &
               this%nQNumbers)
       ENDDO
    ENDDO

    order = (/(i, i=1, nPairs)/)
    !     Sort the pair states using the rank
    !     array. The array order gives the new 
    !     order of the array rank.
    CALL qsorti(order, nPairs, rank)
    
    WRITE(*,*) ' Hole-particle pairs:'
    WRITE(*,*) ' '
    WRITE(*,*) ' hole   particle   ranking number'
    !     Reorder the pairs
    DO pair=1, nPairs
       newPair = order(pair)
       blockQNumbersNew(pair,:) = blockQNumbers(newPair,:)
       pairsNew(pair,:) = pairs(newPair,:)
       rankNew(pair) = rank(newPair)
       !write(*,'(4(x, I6))')  blockQNumbers(newPair,:), rank(newPair)
       write(*,'(3(x, I6))')  pairs(newPair,:), rank(newPair)
    ENDDO
    
    this%hpPairs%blockQNumbers = blockQNumbersNew
    this%hpPairs%pairs = pairsNew
    this%hpPairs%rank = rankNew
    write(*,*) ' '
    
    DEALLOCATE(blockQNumbers)
    DEALLOCATE(blockQNumbersNew)
    DEALLOCATE(pairs, pairsNew)
    DEALLOCATE(rank, rankNew)
    DEALLOCATE(order)
    
    
  END SUBROUTINE setupHpPairsRe
  
  
  !
  !     Set up the matrix list hpph2Index
  !
  SUBROUTINE setupHpph2Index(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    TYPE(Block1) :: thisBlock
    TYPE(IntMatrix) :: hpphMatrix, phPairs, hpPairs
    INTEGER :: indexHpph, nPhPairs, nHpPairs, hp, ph
    INTEGER :: blockIndex
    
    this%hpph2Index = IntMatList_(this%nBlocks)
    
    indexHpph = 0
    DO blockIndex=1, this%nBlocks
       thisBlock = this%list(blockIndex)
       
       phPairs = thisBlock%phPairs
       nPhPairs = phPairs%dim1
       
       hpPairs = thisBlock%hpPairs
       nHpPairs = hpPairs%dim1
       
       hpphMatrix = IntMatrix_(nHpPairs, nPhPairs)
       DO hp=1, nHpPairs
          DO ph=1, nPhPairs
             indexHpph = indexHpph + 1
             
             hpphMatrix%matr(hp, ph) = indexHpph
          ENDDO
       ENDDO
       this%hpph2Index%list(blockIndex) = hpphMatrix
       
       CALL IntMatrix_d(hpphMatrix)
    ENDDO
    
    ALLOCATE(this%index2ReHpph(indexHpph, 3))
    this%index2ReHpph = 0

    
  END SUBROUTINE setupHpph2Index

  !
  !     Set up the matrix list pphh2Index 
  !
  SUBROUTINE setupPphh2Index(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    TYPE(Block1) :: pphhBlock
    TYPE(IntMatrix) :: ppPairs, hhPairs, pphhMatrix
    INTEGER :: nPphh, pphhBlockIndex, pp, hh
    INTEGER :: nPpPairs, nHhPairs, blockIndex
    INTEGER :: indexPphh
    
    this%pphh2Index = IntMatList_(this%blocks%nPphhBlocks)
    
    indexPphh = 0
    DO pphhBlockIndex=1, this%blocks%nPphhBlocks
       blockIndex = this%blocks%pphhBlocks(pphhBlockIndex) 
       pphhBlock = this%blocks%list(blockIndex)
       
       ppPairs = pphhBlock%ppPairs
       hhPairs = pphhBlock%hhPairs
       nPpPairs = ppPairs%dim1
       nHhPairs = hhPairs%dim1
       
       pphhMatrix = IntMatrix_(nPpPairs, nHhPairs)
       DO pp=1, nPpPairs
          DO hh=1, nHhPairs
             indexPphh = indexPphh + 1
             
             pphhMatrix%matr(pp, hh) = indexPphh
          ENDDO
       ENDDO
       this%pphh2Index%list(pphhBlockIndex) = pphhMatrix
       
       CALL IntMatrix_d(pphhMatrix)
    ENDDO
    
    ALLOCATE(this%index2Pphh(indexPphh, 3))
    this%index2Pphh = 0

    
  END SUBROUTINE setupPphh2Index
  
  !
  !     Set up the list pphh2ReHpph, which gives
  !     a number representing a cross-coupled 
  !     state (blockIndexRe, hp, ph) when a 
  !     normal state (blockIndex, pp, hh) is
  !     the input.
  !     
  SUBROUTINE setupPphh2ReHpph(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    CLASS(SpBasis), POINTER :: basis
    TYPE(BlockList) :: blocks
    TYPE(Block1) :: pphhBlock, hpphBlock
    TYPE(IntMatrix) :: ppPairs, hhPairs
    TYPE(IntMatrix) :: hpphIndices
    INTEGER, ALLOCATABLE :: ij(:), ab(:)
    INTEGER, ALLOCATABLE :: ia(:), bj(:)
    INTEGER, ALLOCATABLE :: qNrsA(:), qNrsB(:)
    INTEGER, ALLOCATABLE :: qNrsI(:), qNrsJ(:)
    INTEGER, ALLOCATABLE :: blockQNrsIa(:)
    INTEGER, ALLOCATABLE :: blockQNrsBj(:)
    INTEGER, ALLOCATABLE :: phThis(:), hpThis(:)
    INTEGER :: indexPphh, blockIndex, pp, hh
    INTEGER :: nPpPairs, nHhPairs, rankIa, rankAi, rankBj
    INTEGER :: hpphBlockIndex, hp, ph
    INTEGER :: pphhBlockIndex, maxPh, maxHp
    LOGICAL :: foundHp, foundPh, occupied
    
    ALLOCATE(ab(2), ij(2))
    ALLOCATE(ia(2), bj(2))
    ALLOCATE(phThis(2), hpThis(2))
    ALLOCATE(qNrsA(this%nQNumbers))
    ALLOCATE(qNrsB(this%nQNumbers))
    ALLOCATE(qNrsI(this%nQNumbers))
    ALLOCATE(qNrsJ(this%nQNumbers))
    ALLOCATE(blockQNrsIa(this%nQNumbers))
    ALLOCATE(blockQNrsBj(this%nQNumbers))
    
    basis => this%basis
    blocks = this%blocks
    
    indexPphh = 0
    DO pphhBlockIndex=1, blocks%nPphhBlocks
       blockIndex = blocks%pphhBlocks(pphhBlockIndex)
       
       pphhBlock = blocks%list(blockIndex)
       ppPairs = pphhBlock%ppPairs
       hhPairs = pphhBlock%hhPairs
       nPpPairs = ppPairs%dim1
       nHhPairs = hhPairs%dim1
       
       hpphIndices = IntMatrix_(nPpPairs, nHhPairs)
       occupied = .FALSE.
       
       !     Loop over pphh pairs
       DO pp=1, nPpPairs
          ab(:) = ppPairs%matr(pp,:)
          DO hh=1, nHhPairs
             ij(:) = hhPairs%matr(hh,:)
             
             indexPphh = indexPphh + 1
             !     Store the pphh pair
             this%index2Pphh(indexPphh,:) &
                  = (/pphhBlockIndex, pp, hh/)
             
             !     Hole-particle and particle-hole
             !     states ia and bj
             ia = (/ij(1), ab(1)/)
             bj = (/ab(2), ij(2)/)
             
             !     Single-particle quantum numbers
             qNrsA(:) = basis%spQuantNr(ab(1),:)
             qNrsB(:) = basis%spQuantNr(ab(2),:)
             qNrsI(:) = basis%spQuantNr(ij(1),:)
             qNrsJ(:) = basis%spQuantNr(ij(2),:)
             
             !     ab ij => ia bj
             blockQNrsIa(:) = qNrsI(:) - qNrsA(:)
             blockQNrsBj(:) = qNrsB(:) - qNrsJ(:)
             
             !     Ranking index for cross-coupled
             !     blocks (k_rel, m_s_rel)
             rankIa = getRank(basis, blockQNrsIa, &
                  this%nQNumbers)
             rankBj = getRank(basis, blockQNrsBj, &
                  this%nQNumbers)
             rankAi = getRank(basis, -blockQNrsIa, &
                  this%nQNumbers)
             
             IF (rankIa == rankBj) THEN
                ! write(*,*) 'a:',qNrsA
                ! write(*,*) 'b:',qNrsB
                ! write(*,*) 'i:',qNrsI
                ! write(*,*) 'j:',qNrsJ
                ! write(*,*) ' '
                
                !     Corresponding (k_rel, m_s_rel)
                !     block index
                hpphBlockIndex = this%rank2Block(rankIa)
                hpphBlock = this%list(hpphBlockIndex)
                
                !     Search for the corresponding hp
                !     and ph pairs
                hp = 0
                foundHp = .FALSE.
                maxHp = hpphBlock%hpPairs%dim1
                
                DO WHILE ((hp < maxHp) .AND. .NOT.foundHp)
                   hp = hp + 1
                   hpThis = hpphBlock%hpPairs%matr(hp,:)
                   IF ((hpThis(1) == ia(1)) .AND. &
                        (hpThis(2) == ia(2))) foundHp = .TRUE.
                ENDDO
                ph = 0
                foundPh = .FALSE.
                maxPh = hpphBlock%phPairs%dim1
                DO WHILE ((ph < maxPh) .AND. .NOT.foundPh)
                   ph = ph + 1
                   phThis = hpphBlock%phPairs%matr(ph,:)
                   IF ((phThis(1) == bj(1)) .AND. &
                        (phThis(2) == bj(2))) foundPh = .TRUE.
                ENDDO
                
                IF (foundHp .AND. foundPh) THEN
                   hpphIndices%matr(pp, hh) &
                        = this%hpph2Index%list(hpphBlockIndex)%matr(hp, ph)
                   occupied = .TRUE.
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       IF (.NOT. occupied) hpphIndices%dim1 = 0
       this%pphh2ReHpph%list(pphhBlockIndex) = hpphIndices
       
       CALL IntMatrix_d(hpphIndices)
    ENDDO
    
    DEALLOCATE(ab, ij)
    DEALLOCATE(ia, bj)
    DEALLOCATE(phThis, hpThis)
    DEALLOCATE(qNrsA, qNrsB)
    DEALLOCATE(qNrsI, qNrsJ)
    DEALLOCATE(blockQNrsIa)
    DEALLOCATE(blockQNrsBj)
    
    
  END SUBROUTINE setupPphh2ReHpph
  
  !
  !     Set up the list hpphRe2Pphh, which gives
  !     a number representing a state 
  !     (blockIndex, pp, hh) when a cross-coupled
  !     state (blockIndexRe, hp, ph) is the input.
  !
  SUBROUTINE setupHpphRe2Pphh(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    CLASS(SpBasis), POINTER :: basis
    TYPE(Block1) :: thisHpphBlock, pphhBlock
    TYPE(IntMatrix) :: phPairs, hpPairs
    TYPE(IntMatrix) :: pphhIndices
    INTEGER, ALLOCATABLE :: bra(:), ket(:)
    INTEGER, ALLOCATABLE :: ij(:), ab(:)
    INTEGER, ALLOCATABLE :: ijThis(:), abThis(:)
    INTEGER, ALLOCATABLE :: qNrsA(:), qNrsB(:)
    INTEGER, ALLOCATABLE :: qNrsI(:), qNrsJ(:)
    INTEGER, ALLOCATABLE :: blockQNrsPp(:)
    INTEGER, ALLOCATABLE :: blockQNrsHh(:)
    INTEGER :: hp, ph, nPhPairs, nHpPairs, a, b, i, j
    INTEGER :: pp, hh, blockIndex, indexHpph
    INTEGER :: rankPp, rankHh, pphhBlockIndex
    LOGICAL :: foundHh, foundPp, occupied
    INTEGER :: maxHh, maxPp, blockIndexCM
    
    ALLOCATE(bra(2), ket(2))
    ALLOCATE(ab(2), ij(2))
    ALLOCATE(abThis(2), ijThis(2))
    ALLOCATE(qNrsA(this%nQNumbers))
    ALLOCATE(qNrsB(this%nQNumbers))
    ALLOCATE(qNrsI(this%nQNumbers))
    ALLOCATE(qNrsJ(this%nQNumbers))
    ALLOCATE(blockQNrsPp(this%nQNumbers))
    ALLOCATE(blockQNrsHh(this%nQNumbers))
    basis => this%basis
    
    indexHpph = 0
    DO blockIndex=1, this%nBlocks
       
       thisHpphBlock = this%list(blockIndex)
       phPairs = thisHpphBlock%phPairs
       hpPairs = thisHpphBlock%hpPairs
       nPhPairs = phPairs%dim1
       nHpPairs = hpPairs%dim1
       pphhIndices = IntMatrix_(nHpPairs, nPhPairs)
       occupied = .FALSE.
       
       !     Loop over hpph pairs
       DO hp=1, nHpPairs
          bra(:) = hpPairs%matr(hp,:)
          DO ph=1, nPhPairs
             ket(:) = phPairs%matr(ph,:)
             
             indexHpph = indexHpph + 1
             !     Store the hpph pair 
             this%index2ReHpph(indexHpph,:) &
                  = (/blockIndex, hp, ph/) 
             
             !     Particle-particle and hole-hole 
             !     states ab and ij
             ab = (/bra(2), ket(1)/)
             ij = (/bra(1), ket(2)/)
             
             !     Single-partcle quantum numbers
             qNrsA(:) = basis%spQuantNr(ab(1),:)
             qNrsB(:) = basis%spQuantNr(ab(2),:)
             qNrsI(:) = basis%spQuantNr(ij(1),:)
             qNrsJ(:) = basis%spQuantNr(ij(2),:)
             !     ia bj => ab ij
             blockQNrsPp(:) = qNrsA(:) + qNrsB(:)
             blockQNrsHh(:) = qNrsI(:) + qNrsJ(:)
             
             !     Ranking index in the normal case
             !     with blocks (K_CM, M_S)
             rankPp = getRank(basis, blockQNrsPp, &
                  this%nQNumbers)
             rankHh = getRank(basis, blockQNrsHh, &
                  this%nQNumbers)
             
             IF (rankPp == rankHh) THEN
                
                !     Corresponding (K_CM, M_S) block index
                pphhBlockIndex &
                     = this%blocks%rank2BlockPphh(rankPp)
                blockIndexCM = this%blocks%pphhBlocks(pphhBlockIndex)
                
                pphhBlock = this%blocks%list(blockIndexCM)
                
                ! if (blockIndex == 1) then
                !    write(*,*) 'hp=',bra
                !    write(*,*) 'ph=',ket
                !    write(*,*) 'ab=',ab
                !    write(*,*) 'ij=',ij
                !    write(*,*) 'occHh=',pphhBlock%occupiedHh
                !    write(*,*) 'occPp=',pphhBlock%occupiedPp
                !    write(*,*) 'pphhBlockIndex=',pphhBlockIndex
                ! endif

                !     Search for the corresponding pp and
                !     hh pairs
                IF (pphhBlock%occupiedHh .AND. &
                     pphhBlock%occupiedPp) THEN
                   
                   hh = 0
                   foundHh = .FALSE.
                   maxHh = pphhBlock%hhPairs%dim1
                   DO WHILE ((hh < maxHh) .AND. .NOT.foundHh) 
                      hh = hh + 1
                      
                      ijThis = pphhBlock%hhPairs%matr(hh,:)
                      IF ((ijThis(1) == ij(1)) .AND. &
                           (ijThis(2) == ij(2))) &
                           foundHh = .TRUE.
                   ENDDO
                   
                   pp = 0
                   foundPp = .FALSE.
                   maxPp = pphhBlock%ppPairs%dim1
                   DO WHILE ((pp < maxPp) .AND. .NOT.foundPp)
                      pp = pp + 1
                      
                      abThis = pphhBlock%ppPairs%matr(pp,:)
                      IF ((abThis(1) == ab(1)) .AND. &
                           (abThis(2) == ab(2))) &
                           foundPp = .TRUE.
                   ENDDO
                   
                   IF (foundHh .AND. foundPp) THEN
                      pphhIndices%matr(hp, ph) &
                           = this%pphh2Index%list(pphhBlockIndex)%matr(pp, hh)        
                      occupied = .TRUE.
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ENDDO
      
       ! if (blockIndex == 1) then
       !    write(*,*) 'pphhIndices%matr=',pphhIndices%matr
       !    write(*,*) 'foundHh=',foundHh,',foundPp=',foundPp 
       ! endif
       IF (.NOT. occupied) pphhIndices%dim1 = 0
       this%hpphRe2Pphh%list(blockIndex) = pphhIndices
       
       CALL IntMatrix_d(pphhIndices)
    ENDDO
    
    DEALLOCATE(bra, ket)
    DEALLOCATE(ab, ij)
    DEALLOCATE(abThis, ijThis)
    DEALLOCATE(qNrsA, qNrsB)
    DEALLOCATE(qNrsI, qNrsJ)
    DEALLOCATE(blockQNrsPp)
    DEALLOCATE(blockQNrsHh)
    
    
  END SUBROUTINE setupHpphRe2Pphh
  
  
END MODULE BlockListReMod
