
!
!     Class for list of blocks where the single-particle
!     states are recoupled.
!
MODULE BlockListReMod
  USE MemoryCounterMod
  USE IntMatListMod
  USE BlockListMod
  USE PairListMod
  USE SpBasisMod
  USE Block1Mod
  IMPLICIT NONE
  
  TYPE :: BlockListRe
     !
     !     nBlocks:    Number of cross-coupled blocks
     !     nQNumbers:    Numbers of quantum numbers
     !     phPairs:      Ordered list of all particle-hole
     !                   pairs
     !     hpPairs:      Ordered list of all hole-particle
     !                   pairs
     !
     TYPE(BlockList), POINTER :: blocks
     CLASS(SpBasis), POINTER :: basis
     INTEGER :: nBlocks, nQNumbers
     TYPE(Block1), ALLOCATABLE :: list(:)
     TYPE(PairList) :: phPairs
     TYPE(PairList) :: hpPairs
     TYPE(IntMatList) :: abij2ReIabj
     TYPE(IntMatList) :: iabjRe2Abij
     TYPE(IntMatList) :: ibajRe2Abij
     TYPE(IntMatList) :: jabiRe2Abij
     TYPE(IntMatList) :: jbaiRe2Abij
     TYPE(IntMatList) :: abij2Index
     TYPE(IntMatList) :: iabj2Index
     INTEGER, ALLOCATABLE :: index2ReIabj(:,:)
     !INTEGER, ALLOCATABLE :: index2ReIbaj(:,:)
     !INTEGER, ALLOCATABLE :: index2ReJabi(:,:)
     !INTEGER, ALLOCATABLE :: index2ReJbai(:,:)
     INTEGER :: nHpph, nPphhAllElements
     INTEGER :: nElementsBeginThis
     INTEGER, ALLOCATABLE :: index2Abij(:,:)
     INTEGER, ALLOCATABLE :: rank2Block(:)
     INTEGER, ALLOCATABLE :: elementsProcs(:)
     INTEGER, ALLOCATABLE :: elementsBeginPs(:)
     INTEGER :: maxRank
     INTEGER :: nBlocksReThis
     INTEGER, ALLOCATABLE :: blocksThis2BlocksAll(:)
   CONTAINS
     PROCEDURE :: setupHpphRe2Pphh
     PROCEDURE :: setupAbij2ReIabj
     PROCEDURE :: setupAbij2Index
     PROCEDURE :: setupIabj2Index
     PROCEDURE :: setupBlocksReThis
     PROCEDURE :: setupPhBlockList
     PROCEDURE :: setupBlocksRe
     PROCEDURE :: setupPhPairsRe
     PROCEDURE :: setupHpPairsRe
     PROCEDURE :: countBlocksRe
  END TYPE BlockListRe
  
  
CONTAINS
  
  !
  !     Constructor for BlockListRe object.
  !
  FUNCTION BlockListRe_(nQnumbers, basis, blocks) RESULT(this)
    INTEGER, INTENT(IN) :: nQNumbers
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    TYPE(BlockList), TARGET, INTENT(IN) :: blocks
    TYPE(BlockListRe), POINTER :: this
    INTEGER :: nHoles, nParticles, nPhPairs
    
    ALLOCATE(this)
    
    this%basis => basis
    this%blocks => blocks
    this%nQNumbers = nQNumbers
    
    nPhPairs = basis%nOccupied*basis%nUnoccupied       
    
    this%phPairs = PairList_(nPhPairs, nQNumbers)
    this%hpPairs = PairList_(nPhPairs, nQNumbers)
    
    
  END FUNCTION BlockListRe_
  
  !
  !     Destructor for BlockListRe object
  !
  SUBROUTINE BlockListRe_d(this)
    TYPE(BlockListRe), INTENT(INOUT) :: this
    INTEGER :: blockRe, pphhBlock
    
    CALL PairList_d(this%phPairs)
    CALL PairList_d(this%hpPairs)
    
    IF (ALLOCATED(this%list)) THEN
       DO blockRe=1, this%nBlocks
          CALL Block1_d(this%list(blockRe))
       ENDDO
       DEALLOCATE(this%list)
    ENDIF
    
    CALL IntMatList_d(this%abij2ReIabj)
    CALL IntMatList_d(this%iabjRe2Abij)
    CALL IntMatList_d(this%ibajRe2Abij)
    CALL IntMatList_d(this%jabiRe2Abij)
    CALL IntMatList_d(this%jbaiRe2Abij)
    CALL IntMatList_d(this%iabj2Index)
    CALL IntMatList_d(this%abij2Index)

    IF (ALLOCATED(this%index2ReIabj)) THEN
       DEALLOCATE(this%index2ReIabj)
       CALL add2MemoryInt(-3*this%nHpph)
    ENDIF
    !IF (ALLOCATED(this%index2ReIbaj)) &
    !     DEALLOCATE(this%index2ReIbaj)
    !IF (ALLOCATED(this%index2ReJabi)) &
    !     DEALLOCATE(this%index2ReJabi)
    !IF (ALLOCATED(this%index2ReJbai)) &
    !     DEALLOCATE(this%index2ReJbai)
    IF (ALLOCATED(this%index2Abij)) THEN
       DEALLOCATE(this%index2Abij)
       CALL add2MemoryInt(-3*this%nPphhAllElements)
    ENDIF
    IF (ALLOCATED(this%rank2Block)) THEN
       DEALLOCATE(this%rank2Block)
       CALL add2MemoryInt(-this%maxRank)
    ENDIF

    IF (ASSOCIATED(this%blocks)) THEN
       CALL BlockList_d(this%blocks)
       this%blocks => NULL()
    ENDIF
    
    IF (ASSOCIATED(this%basis)) THEN
       CALL SpBasis_d(this%basis)
       this%basis => NULL()
    ENDIF
    
    
  END SUBROUTINE BlockListRe_d
  
  SUBROUTINE setupPhBlockList(this)
    CLASS(BlockListRe), TARGET, INTENT(INOUT) :: this
    
    IF (this%basis%nUnoccupied > 0) THEN
       
       IF (this%basis%mpiJobs%iAm == 0) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'Setting up recoupled blocks...'
          WRITE(*,*) ' '
       ENDIF
       
       !     Set up pairs of ph and hp single-particle states
       CALL setupPhPairsRe(this)      
       CALL setupHpPairsRe(this)
       
       !     Count number of (k_rel, m_s_rel) blocks
       CALL countBlocksRe(this)
       !     Set up list of blocks
       CALL setupBlocksRe(this)
       
       !     As part of the MPI parallelization,
       !     assign cross-coupled blocks to the different
       !     processes.
       CALL setupBlocksReThis(this)
       
       CALL setupIabj2Index(this)
       CALL setupAbij2Index(this)
       
       CALL setupHpphRe2Pphh(this, 'iabj')
       CALL setupHpphRe2Pphh(this, 'ibaj')
       CALL setupHpphRe2Pphh(this, 'jabi')
       CALL setupHpphRe2Pphh(this, 'jbai')
       
       !write(*,*) ''
       !if (this%basis%mpiJobs%iAm == 1) then
       !write(*,*) '||>> iam=',this%basis%mpiJobs%iAm 
       !stop
       !endif
       
       CALL setupAbij2ReIabj(this)
       
       
    ENDIF
    
    
  END SUBROUTINE setupPhBlockList
  
  !
  !     Set up list of (k_rel, m_s_rel) blocks
  !
  SUBROUTINE setupBlocksRe(this)
    CLASS(BlockListRe), TARGET, INTENT(INOUT) :: this
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
          
          ! WRITE(*,*) ' '
          ! WRITE(*,*) ' Block nr',thisBlock
          pairsMatrPh = IntMatrix_(nOccPh, 2)
          
          DO phThis=1, nOccPh
             pairsMatrPh%matr(phThis,:) &
                  = phPairs%pairs(ph-nOccPh+phThis,:)
             
             ! WRITE(*,*) phPairs%pairs(ph-nOccPh+phThis,:), &
             !      phPairs%rank(ph-nOccPh+phThis)
             !WRITE(*,*) phPairs%blockQNumbers(ph-nOcc+phThis,:), &
             !    phPairs%rank(ph-nOcc+phThis)
          ENDDO
          this%list(thisBlock)%phPairs = pairsMatrPh
          ! if (thisBlock == 58) then
          !    write(*,*) ':nOccPh=',nOccPh,',rank=',thisRank
          ! endif

          CALL IntMatrix_d(pairsMatrPh)
          
          pairsMatrHp = IntMatrix_(nOccHp, 2)
          
!          WRITE(*,*) ' '
          DO hpThis=1, nOccHp
             pairsMatrHp%matr(hpThis,:) &
                  = hpPairs%pairs(hp-nOccHp+hpThis,:)
             
             ! WRITE(*,*) hpPairs%pairs(hp-nOccHp+hpThis,:), &
             !      hpPairs%rank(hp-nOccHp+hpThis)
          ENDDO
          this%list(thisBlock)%hpPairs = pairsMatrHp
          
          CALL add2MemoryInt(2*(nOccPh + nOccHp))
          CALL IntMatrix_d(pairsMatrHp)
       ENDIF
       
       thisRank = thisRank + 1
       
    ENDDO
    !   WRITE(*,*) ' '
    
    
  END SUBROUTINE setupBlocksRe
  
  !
  !     Count the total number of blocks
  !
  SUBROUTINE countBlocksRe(this)
    CLASS(BlockListRe), TARGET, INTENT(INOUT) :: this
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
    
    IF (this%basis%mpiJobs%iAm == 0) THEN
       WRITE(*,*) ' '
       WRITE(*,*) ' Number of ph pairs:', nPairsPh
       WRITE(*,*) ' Number of hp pairs:', nPairsHp
    ENDIF
    
    maxRank = MAX(phPairs%rank(nPairsPh), &
         hpPairs%rank(nPairsHp))
    ALLOCATE(this%rank2Block(maxRank))
    this%rank2Block = 0
    this%maxRank = maxRank
    CALL add2MemoryInt(maxRank)
    
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
    IF (this%basis%mpiJobs%iAm == 1) THEN
       WRITE(*,*) ' Number of cross-coupled blocks: ', blocks
    ENDIF
    ALLOCATE(this%list(this%nBlocks))
    this%abij2ReIabj = IntMatList_(this%blocks%nPphhBlocksThis)
    
    
  END SUBROUTINE countBlocksRe
  
  !
  !     Divide the cross-coupled blocks to different 
  !     MPI processes.
  !
  SUBROUTINE setupBlocksReThis(this)
    CLASS(BlockListRe), TARGET, INTENT(INOUT) :: this
    CLASS(SpBasis), POINTER :: basis
    INTEGER, ALLOCATABLE :: elementsInBlocks(:)
    TYPE(Block1) :: thisBlock
    TYPE(IntMatrix) :: phPairs
    INTEGER :: element, nPhPairs, elementsTot
    INTEGER :: nProcesses, elementsPerProc
    INTEGER :: process, elementsThis, nBlocksThis
    INTEGER :: reBlockBegin, reBlockEnd, reBlock
    INTEGER :: elementsAll, blockIndex, iError
    INTEGER :: elementsThisSave, nElementsBegin
    INCLUDE 'mpif.h'
    
    basis => this%basis
    ALLOCATE(elementsInBlocks(this%nBlocks))
    
    ! IF (basis%mpiJobs%iAm == 0) THEN
    !    WRITE(*,*) ' '
    !    WRITE(*,*) 'blockRe     elements'
    ! ENDIF
    !     Count elements in the different blocks
    elementsTot = 0
    DO blockIndex=1, this%nBlocks
       
       thisBlock = this%list(blockIndex)
       
       phPairs = thisBlock%phPairs
       nPhPairs = phPairs%dim1
       
       !     Number of elements in this block
       elementsInBlocks(blockIndex) = nPhPairs**2
       ! IF (basis%mpiJobs%iAm == 0) THEN
       !    WRITE(*,*) blockIndex, nPhPairs**2
       ! ENDIF
       
       elementsTot = elementsTot + nPhPairs**2
    ENDDO
    CALL basis%mpiJobs%writeMpi(' ')
    
    !     Number of processes
    nProcesses = basis%mpiJobs%nProcesses
    !     Approximate number of elements per process
    elementsPerProc = &
         CEILING(1.d0*elementsTot/nProcesses)
    
    ! IF (basis%mpiJobs%iAm == 0) THEN
    !    WRITE(*,*) 'Total number of elements: ', elementsTot
    !    WRITE(*,*) 'Elements per process: ', elementsPerProc
    ! ENDIF
    
    process = 0
    elementsThis = 0
    nBlocksThis = 0
    reBlockBegin = 1
    reBlock = 0
    elementsAll = 0
    nElementsBegin = 0
    DO WHILE ((process <= basis%mpiJobs%iAm).AND. &
         (reBlock <= this%nBlocks)) 
       reBlock = reBlock + 1
       
       !     Number of elements assigned to
       !     this process
       elementsThis = elementsThis &
            + elementsInBlocks(reBlock)
       elementsAll = elementsAll & 
            + elementsInBlocks(reBlock)
       nElementsBegin = nElementsBegin &
            + elementsInBlocks(reBlock)
       
       IF (process == basis%mpiJobs%iAm) THEN
          nBlocksThis = nBlocksThis + 1
       ENDIF
       
       IF ((elementsThis >= elementsPerProc) &
            .OR.(elementsAll == elementsTot)) THEN
          process = process + 1
          elementsThisSave = elementsThis
          elementsThis = 0
          this%nElementsBeginThis = nElementsBegin &
               - elementsThisSave
          IF (process - 1 == basis%mpiJobs%iAm) &
               reBlockEnd = reBlock
          IF (process == basis%mpiJobs%iAm) &
               reBlockBegin = reBlock + 1
       ENDIF
    ENDDO
    
    this%nBlocksReThis = nBlocksThis
    ALLOCATE(this%blocksThis2BlocksAll(nBlocksThis))
    
    ! WRITE(*,*) 'iAm:',basis%mpiJobs%iAm, &
    !      ',begin:',reBlockBegin,',end:',reBlockEnd
    ! WRITE(*,*) 'iAm:',basis%mpiJobs%iAm, &
    !      ',elements:',elementsThisSave
    ! WRITE(*,*) 'nBlocksThis=',nBlocksThis, &
    !      ',=',reBlockEnd-reBlockBegin+1
    ! write(*,*) 'iAm:',basis%mpiJobs%iAm, &
    !      ', nEleBegin=',this%nElementsBeginThis
    
    element = 0
    DO reBlock=reBlockBegin, reBlockEnd
       element = element + 1
       this%blocksThis2BlocksAll(element) = reBlock
    ENDDO
    
    ALLOCATE(this%elementsProcs(basis%mpiJobs%nProcesses))
    ALLOCATE(this%elementsBeginPs(basis%mpiJobs%nProcesses))
    
    CALL MPI_ALLGATHER(elementsThisSave, 1, MPI_INTEGER, &
         this%elementsProcs, 1, MPI_INTEGER, &
         MPI_COMM_WORLD, iError)
    
    this%elementsBeginPs(1) = 0
    DO process=2, basis%mpiJobs%nProcesses
       this%elementsBeginPs(process) = &
            this%elementsBeginPs(process-1) &
            + this%elementsProcs(process-1)
    ENDDO
    
    ! write(*,*) 'this%elementsProcs=',&
    !      this%elementsProcs
    ! write(*,*) 'this%elementsBeginPs=',&
    !      this%elementsBeginPs
    
    DEALLOCATE(elementsInBlocks)
    
    
  END SUBROUTINE setupBlocksReThis
  
  
  SUBROUTINE setupPhPairsRe(this)
    CLASS(BlockListRe), TARGET, INTENT(INOUT) :: this
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
    
    ! WRITE(*,*) ' '
    ! WRITE(*,*) ' Particle-hole pairs:'
    ! WRITE(*,*) ' '
    ! WRITE(*,*) ' particle   hole   ranking number'
    !     Reorder the pairs
    DO pair=1, nPairs
       newPair = order(pair)
       blockQNumbersNew(pair,:) = blockQNumbers(newPair,:)
       pairsNew(pair,:) = pairs(newPair,:)
       rankNew(pair) = rank(newPair)
       !write(*,'(4(x, I6))')  blockQNumbers(newPair,:), rank(newPair)
!       write(*,'(3(x, I6))')  pairs(newPair,:), rank(newPair)
    ENDDO
    
    this%phPairs%blockQNumbers = blockQNumbersNew
    this%phPairs%pairs = pairsNew
    this%phPairs%rank = rankNew
!    write(*,*) ' '
    
    DEALLOCATE(blockQNumbers)
    DEALLOCATE(blockQNumbersNew)
    DEALLOCATE(pairs)
    DEALLOCATE(pairsNew)
    DEALLOCATE(rank)
    DEALLOCATE(rankNew)
    DEALLOCATE(order)
    !stop
    
  END SUBROUTINE setupPhPairsRe
  
  
  SUBROUTINE setupHpPairsRe(this)
    CLASS(BlockListRe), TARGET, INTENT(INOUT) :: this
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
    
    ! WRITE(*,*) ' Hole-particle pairs:'
    ! WRITE(*,*) ' '
    ! WRITE(*,*) ' hole   particle   ranking number'
    !     Reorder the pairs
    DO pair=1, nPairs
       newPair = order(pair)
       blockQNumbersNew(pair,:) = blockQNumbers(newPair,:)
       pairsNew(pair,:) = pairs(newPair,:)
       rankNew(pair) = rank(newPair)
       !write(*,'(4(x, I6))')  blockQNumbers(newPair,:), rank(newPair)
 !      write(*,'(3(x, I6))')  pairs(newPair,:), rank(newPair)
    ENDDO
    
    this%hpPairs%blockQNumbers = blockQNumbersNew
    this%hpPairs%pairs = pairsNew
    this%hpPairs%rank = rankNew
!    write(*,*) ' '
    
    DEALLOCATE(blockQNumbers)
    DEALLOCATE(blockQNumbersNew)
    DEALLOCATE(pairs, pairsNew)
    DEALLOCATE(rank, rankNew)
    DEALLOCATE(order)
    
    
  END SUBROUTINE setupHpPairsRe
  
  
  !
  !     Set up the matrix list iabj2Index
  !
  SUBROUTINE setupIabj2Index(this)
    CLASS(BlockListRe), TARGET, INTENT(INOUT) :: this
    TYPE(Block1) :: thisBlock
    TYPE(IntMatrix) :: hpphMatrix, phPairs, hpPairs
    INTEGER :: indexHpph, nPhPairs, nHpPairs, hp, ph
    INTEGER :: blockIndex
    
    this%iabj2Index = IntMatList_(this%nBlocks)
    
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
       this%iabj2Index%list(blockIndex) = hpphMatrix
       
       CALL add2MemoryInt(nHpPairs*nPhPairs)  
       CALL IntMatrix_d(hpphMatrix)
       !CALL IntMatrix_d(phPairs)
       !CALL IntMatrix_d(hpPairs)
    ENDDO
    
    this%nHpph = indexHpph
    ALLOCATE(this%index2ReIabj(3, indexHpph))
    !ALLOCATE(this%index2ReIbaj(3, indexHpph))
    !ALLOCATE(this%index2ReJabi(3, indexHpph))
    !ALLOCATE(this%index2ReJbai(3, indexHpph))
    this%index2ReIabj = 0
    !this%index2ReIbaj = 0
    !this%index2ReJabi = 0
    !this%index2ReJbai = 0
    CALL add2MemoryInt(3*indexHpph)
    
    
  END SUBROUTINE setupIabj2Index
  
  !
  !     Set up the matrix list abij2Index 
  !
  SUBROUTINE setupAbij2Index(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    TYPE(Block1) :: pphhBlock
    TYPE(IntMatrix) :: ppPairs, hhPairs, pphhMatrix
    INTEGER :: nPphh, pphhBlockIndex, pp, hh
    INTEGER :: nPpPairs, nHhPairs, blockIndex
    INTEGER :: indexPphh
    
    this%abij2Index = IntMatList_(this%blocks%nPphhBlocksAll)
    
    indexPphh = 0
    DO pphhBlockIndex=1, this%blocks%nPphhBlocksAll
       blockIndex = this%blocks%pphh2BlocksAll(pphhBlockIndex) 
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
       this%abij2Index%list(pphhBlockIndex) = pphhMatrix
       
       ! if (pphhBlockIndex == 14) then
       !    write(*,*) 'npp:',nPpPairs,',nhh:',nHhPairs 
       !    write(*,*) 'm:',this%abij2Index%list(pphhBlockIndex)%matr(2,1)
       ! endif
       
       CALL add2MemoryInt(nPpPairs*nHhPairs)      
       CALL IntMatrix_d(pphhMatrix)
       !CALL IntMatrix_d(hhPairs)
       !CALL IntMatrix_d(ppPairs)
    ENDDO
    
    this%nPphhAllElements = indexPphh
    ! write(*,*) 'this%nPphhAll=',this%nPphhAllElements
    ! write(*,*) 'this%blocks%nPphhBlocksAll=',&
    !this%blocks%nPphhBlocksAll
    ALLOCATE(this%index2Abij(3, indexPphh))
    ! if (this%basis%mpiJobs%iAm == 1) then
    !    write(*,*) 'size index2Abij=',indexPphh
    ! endif
    this%index2Abij = 0
    CALL add2MemoryInt(3*indexPphh)
    
    
  END SUBROUTINE setupAbij2Index
  
  !
  !     Set up the list abij2ReIabj, which gives
  !     a number representing a cross-coupled 
  !     state (blockIndexRe, hp, ph) when a 
  !     normal state (blockIndex, pp, hh) is
  !     the input.
  !     
  SUBROUTINE setupAbij2ReIabj(this)
    CLASS(BlockListRe), TARGET, INTENT(INOUT) :: this
    CLASS(SpBasis), POINTER :: basis
    TYPE(Block1) :: pphhBlock, hpphBlock
    TYPE(IntMatrix) :: ppPairs, hhPairs
    TYPE(IntMatrix) :: hpphIndices, hpphIndices0
    INTEGER, ALLOCATABLE :: ij(:), ab(:)
    INTEGER, ALLOCATABLE :: ia(:), bj(:)
    INTEGER, ALLOCATABLE :: qNrsA(:), qNrsB(:)
    INTEGER, ALLOCATABLE :: qNrsI(:), qNrsJ(:)
    INTEGER, ALLOCATABLE :: blockQNrsIa(:)
    INTEGER, ALLOCATABLE :: blockQNrsBj(:)
    INTEGER, ALLOCATABLE :: phThis(:), hpThis(:)
    INTEGER, ALLOCATABLE :: displacements(:)
    INTEGER, ALLOCATABLE :: elementsV(:)
    INTEGER, ALLOCATABLE :: index2Abij(:,:)
    INTEGER :: indexPphh, blockIndex, pp, hh
    INTEGER :: nPpPairs, nHhPairs, rankIa, rankAi, rankBj
    INTEGER :: hpphBlockIndex, hp, ph, pphhBlockAll
    INTEGER :: pphhBlockIndex, maxPh, maxHp
    LOGICAL :: foundHp, foundPh, occupied
    INTEGER :: iError, indexBegin, indexEnd, pphhIndices
    INCLUDE 'mpif.h'
    
    ALLOCATE(ab(2), ij(2))
    ALLOCATE(ia(2), bj(2))
    ALLOCATE(phThis(2), hpThis(2))
    ALLOCATE(qNrsA(this%nQNumbers))
    ALLOCATE(qNrsB(this%nQNumbers))
    ALLOCATE(qNrsI(this%nQNumbers))
    ALLOCATE(qNrsJ(this%nQNumbers))
    ALLOCATE(blockQNrsIa(this%nQNumbers))
    ALLOCATE(blockQNrsBj(this%nQNumbers))
    ALLOCATE(index2Abij(3, &
         this%blocks%pphhElementsProcs(this%basis%mpiJobs%iAm+1)))
    
    basis => this%basis
    
    hpphIndices0 = IntMatrix_(1, 1)
    hpphIndices0%dim1 = 0
    
    indexPphh = 0!this%blocks%nPphhElementsBeginThis 
    !write(*,*) 'start iPphh=',indexPphh
    
    !     Wait until all processes are ready
    !CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
    
    ! if (this%basis%mpiJobs%iAm == 0) then
    !    write(*,*) 'aaa'
    !    CALL wait(10)
    !    write(*,*) 'bbb'
    ! endif
    
    DO pphhBlockIndex=1, this%blocks%nPphhBlocksThis
       blockIndex = this%blocks%pphhThis2BlocksAll(&
            pphhBlockIndex)
       pphhBlockAll = this%blocks%pphhThis2PphhAll(&
            pphhBlockIndex)
       
       pphhBlock = this%blocks%list(blockIndex)
       ppPairs = pphhBlock%ppPairs
       hhPairs = pphhBlock%hhPairs
       nPpPairs = ppPairs%dim1
       nHhPairs = hhPairs%dim1
       
       hpphIndices = IntMatrix_(nPpPairs, nHhPairs)
       occupied = .FALSE. 
       
       !     Wait until all processes are ready
       !CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
       
       !     Loop over pphh pairs
       DO pp=1, nPpPairs
          ab(:) = ppPairs%matr(pp,:)
          DO hh=1, nHhPairs
             ij(:) = hhPairs%matr(hh,:)
             
             indexPphh = indexPphh + 1
             
             !     Store the pphh pair
             index2Abij(:, indexPphh) &
                  = (/pphhBlockAll, pp, hh/)
             
             ! if (this%basis%mpiJobs%iAm == 1) then
             !    write(*,*) 'pphhBlockIndex=',pphhBlockIndex
             !    write(*,*) 'pp=',pp,',hh=',hh
             !    write(*,*) 'nq1=',this%list(58)%hpPairs%matr
             !    write(*,*) 'indexPphh=',indexPphh
             !    hpphBlock = this%list(58)
             ! endif
             
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
             !if (this%basis%mpiJobs%iAm == 1) write(*,*) '20'
             IF (rankIa == rankBj) THEN
                !if (this%basis%mpiJobs%iAm == 1) then
                   !write(*,*) 'nq=',this%list(58)%hpPairs%matr
                   !hpphBlock = this%list(58)
                !endif
                
                !if (this%basis%mpiJobs%iAm == 1) write(*,*) '21'
                !     Corresponding (k_rel, m_s_rel)
                !     block index
                hpphBlockIndex = this%rank2Block(rankIa)

                !     Wait until all processes are ready
                !CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
                ! if (this%basis%mpiJobs%iAm == 1) then
                   
                !    write(*,*) '22'
                !    !write(*,*) 'hpph=',hpphBlock
                !    hpphBlock = this%list(59)
                !    write(*,*) 'pphhBlockIndex=',pphhBlockIndex
                !    write(*,*) 'pp=',pp,',hh=',hh
                !    write(*,*) 'ia=',ia
                !    write(*,*) 'bj=',bj
                !    write(*,*) 'hpphBlockIndex=',&
                !         hpphBlockIndex,',nBlocks=',this%nBlocks
                !    write(*,*) 'rankIa=',rankIa
                ! endif
                hpphBlock = this%list(hpphBlockIndex)
                
                !if (this%basis%mpiJobs%iAm == 1) write(*,*) '23'
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
                        = this%iabj2Index%list(hpphBlockIndex)%matr(hp, ph)
                   occupied = .TRUE.
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       
       IF (.NOT. occupied) THEN
          this%abij2ReIabj%list(pphhBlockIndex) = hpphIndices0
          CALL add2MemoryInt(1)
       ELSE
          this%abij2ReIabj%list(pphhBlockIndex) = hpphIndices
          CALL add2MemoryInt(nPpPairs*nHhPairs)
       ENDIF
       
       CALL IntMatrix_d(hpphIndices)
       !CALL IntMatrix_d(ppPairs)
       !CALL IntMatrix_d(hhPairs)
    ENDDO
    !if (this%basis%mpiJobs%iAm == 1) write(*,*) 'indexPphh=',indexPphh
    
    DEALLOCATE(ab, ij)
    DEALLOCATE(ia, bj)
    DEALLOCATE(phThis, hpThis)
    DEALLOCATE(qNrsA, qNrsB)
    DEALLOCATE(qNrsI, qNrsJ)
    DEALLOCATE(blockQNrsIa)
    DEALLOCATE(blockQNrsBj)
    CALL IntMatrix_d(hpphIndices0)
    
    !     Send this slice of the array this%index2Abij
    !     to all other processes
    !   
    ALLOCATE(elementsV(this%basis%mpiJobs%nProcesses))
    ALLOCATE(displacements(this%basis%mpiJobs%nProcesses))
    elementsV = 3*this%blocks%pphhElementsProcs
    displacements = 3*this%blocks%pphhElementsBeginPs

    ! indexBegin = this%blocks%pphhElementsBeginPs(&
    !      this%basis%mpiJobs%iAm+1) + 1
    ! indexEnd = this%blocks%pphhElementsBeginPs(&
    !      this%basis%mpiJobs%iAm+1) + &
    !      this%blocks%pphhElementsProcs(&
    !      this%basis%mpiJobs%iAm+1)
    
    CALL MPI_ALLGATHERV(index2Abij, &
         elementsV(this%basis%mpiJobs%iAm+1), &
         MPI_INTEGER, this%index2Abij, elementsV, &
         displacements, MPI_INTEGER, MPI_COMM_WORLD, iError)
    
    DEALLOCATE(elementsV)
    DEALLOCATE(displacements)
    DEALLOCATE(index2Abij)
    
    ! !     Wait until all processes are ready
    ! CALL MPI_BARRIER(MPI_COMM_WORLD, iError) 
    
    ! write(*,*) 'iam=',this%basis%mpiJobs%iAm,&
    !      'this%index2Abij(:,781)=',this%index2Abij(:,781)
    ! write(*,*) 'iam=',this%basis%mpiJobs%iAm,&
    !      'this%index2Abij(:,20)=',this%index2Abij(:,20)
    ! write(*,*) '-- i am ',this%basis%mpiJobs%iAm  
    
    
  END SUBROUTINE setupAbij2ReIabj
  
    
  !     
  !     If hpph = 'iabj':
  !
  !       Set up the list iabjRe2Abij, which gives
  !       a number representing a state 
  !       (blockIndex, ab, ij) when a cross-coupled
  !       state (blockIndexRe, ia, bj) is given as 
  !       input.
  !
  !     If hpph = 'ibaj':
  !
  !       Set up the list ibajRe2Abij, which gives
  !       a number representing a state
  !       (blockIndex, ab, ij) when a cross-coupled
  !       state (blockIndexRe, ib, aj) is given as
  !       input.
  !
  !     If hpph = 'jabi':
  !
  !       Set up the list jabiRe2Abij, which gives
  !       a number representing a state
  !       (blockIndex, ab, ij) when a cross-coupled
  !       state (blockIndexRe, ja, bi) is given as
  !       input.
  !
  !     If hpph = 'jbai':
  ! 
  !       Set up the list jbaiRe2Abij, which gives
  !       a number representing a state
  !       (blockIndex, ab, ij) when a cross-coupled
  !       state (blockINdexRe, jb, ai) is given as
  !       input.
  !
  SUBROUTINE setupHpphRe2Pphh(this, hpph)
    CLASS(BlockListRe), TARGET, INTENT(INOUT) :: this
    CHARACTER(LEN=4), INTENT(IN) :: hpph 
    CLASS(SpBasis), POINTER :: basis
    TYPE(Block1) :: pphhBlock
    TYPE(IntMatrix) :: phPairs, hpPairs
    TYPE(IntMatrix) :: pphhIndices, pphhIndices0
    INTEGER, ALLOCATABLE :: bra(:), ket(:)
    INTEGER, ALLOCATABLE :: ij(:), ab(:)
    INTEGER, ALLOCATABLE :: ijThis(:), abThis(:)
    INTEGER, ALLOCATABLE :: qNrsA(:), qNrsB(:)
    INTEGER, ALLOCATABLE :: qNrsI(:), qNrsJ(:)
    INTEGER, ALLOCATABLE :: blockQNrsPp(:)
    INTEGER, ALLOCATABLE :: blockQNrsHh(:)
    INTEGER, ALLOCATABLE :: elementsV(:)
    INTEGER, ALLOCATABLE :: displacements(:)
    INTEGER, ALLOCATABLE :: index2ReIabj(:,:)
    INTEGER :: hp, ph, nPhPairs, nHpPairs, a, b, i, j
    INTEGER :: pp, hh, blockIndex, indexHpph
    INTEGER :: rankPp, rankHh, pphhBlockIndex
    LOGICAL :: foundHh, foundPp, occupied
    INTEGER :: maxHh, maxPp, blockIndexCM, iError
    INTEGER :: blockRe, indexBegin, indexEnd
    INCLUDE 'mpif.h'
    
    IF (hpph == 'iabj') THEN
       this%iabjRe2Abij = IntMatList_(this%nBlocksReThis)
    ELSEIF (hpph == 'ibaj') THEN       
       this%ibajRe2Abij = IntMatList_(this%nBlocksReThis)
    ELSEIF (hpph == 'jabi') THEN    
       this%jabiRe2Abij = IntMatList_(this%nBlocksReThis)
    ELSEIF (hpph == 'jbai') THEN
       this%jbaiRe2Abij = IntMatList_(this%nBlocksReThis)
    ENDIF
    
    ALLOCATE(bra(2), ket(2))
    ALLOCATE(ab(2), ij(2))
    ALLOCATE(abThis(2), ijThis(2))
    ALLOCATE(qNrsA(this%nQNumbers))
    ALLOCATE(qNrsB(this%nQNumbers))
    ALLOCATE(qNrsI(this%nQNumbers))
    ALLOCATE(qNrsJ(this%nQNumbers))
    ALLOCATE(blockQNrsPp(this%nQNumbers))
    ALLOCATE(blockQNrsHh(this%nQNumbers))
    ALLOCATE(index2ReIabj(3, &
         this%elementsProcs(this%basis%mpiJobs%iAm+1)))
    basis => this%basis
    
    pphhIndices0 = IntMatrix_(1, 1)
    pphhIndices0%dim1 = 0
    
    indexHpph = 0 !this%nElementsBeginThis
    !write(*,*) '>>> iHpph=',indexHpph
    
    DO blockRe=1, this%nBlocksReThis
       blockIndex = this%blocksThis2BlocksAll(blockRe)
       
       phPairs = this%list(blockIndex)%phPairs
       hpPairs = this%list(blockIndex)%hpPairs 
       nPhPairs = phPairs%dim1
       nHpPairs = hpPairs%dim1
       
       pphhIndices = IntMatrix_(nHpPairs, nPhPairs)
       occupied = .FALSE.
       
       !     Loop over hpph pairs
       DO hp=1, nHpPairs
          bra(1:2) = hpPairs%matr(hp,1:2)
          DO ph=1, nPhPairs
             ket(1:2) = phPairs%matr(ph,1:2)
             
             indexHpph = indexHpph + 1
             
             ! if ((blockIndex == 76).and.(hp==1)&
             !      .and.(ph==3)) then
             !    write(*,*) 'blockRe=',blockRe
             !    write(*,*) 'Yes! iam=',basis%mpiJobs%iAm
             ! endif
             
             IF (hpph == 'iabj') THEN
                ! if (indexHpph == 781) then
                !    write(*,*) '&& ',blockIndex, hp, ph
                ! endif
                ! if (indexHpph == 20) then
                !    write(*,*) '## ',blockIndex, hp, ph
                ! endif
                !     Store the ia-bj pair 
                index2ReIabj(:, indexHpph) &
                     = (/blockIndex, hp, ph/) 
                
                !     Particle-particle and hole-hole 
                !     states ab and ij
                ab = (/bra(2), ket(1)/)
                ij = (/bra(1), ket(2)/)
                
             ELSEIF (hpph == 'ibaj') THEN
                
                !     Store the ib-aj pair
                !this%index2ReIbaj(indexHpph,:) &
                !     = (/blockIndex, hp, ph/)
                
                !     Particle-particle and hole-hole 
                !     states ab and ij
                ab = (/ket(1), bra(2)/)
                ij = (/bra(1), ket(2)/)
                
             ELSEIF (hpph == 'jabi') THEN
                
                !     Store the ja-bi pair
                !this%index2ReJabi(indexHpph,:) &
                !     = (/blockIndex, hp, ph/)

                !     Particle-particle and hole-hole 
                !     states ab and ij
                ab = (/bra(2), ket(1)/)
                ij = (/ket(2), bra(1)/)
                
             ELSEIF (hpph == 'jbai') THEN
                
                !     Store the jb-ai pair
                !this%index2ReJbai(indexHpph,:) &
                !     = (/blockIndex, hp, ph/)

                !     Particle-particle and hole-hole 
                !     states ab and ij
                ab = (/ket(1), bra(2)/)
                ij = (/ket(2), bra(1)/)

             ENDIF
             
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
                blockIndexCM = this%blocks%pphh2BlocksAll(pphhBlockIndex)
                
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
                           = this%abij2Index%list(pphhBlockIndex)%matr(pp, hh)
                      ! if ((blockIndex == 76).and.(hp==1)&
                      !      .and.(ph==1).and.&
                      !      (hpph == 'iabj')) then
                      !    write(*,*) '§> pphhBlockIndex=',&
                      !         pphhBlockIndex,',pp=',pp,',hh=',hh
                      ! endif
                      occupied = .TRUE.
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       
       IF (.NOT. occupied) THEN
          
          IF (hpph == 'iabj') THEN
             this%iabjRe2Abij%list(blockRe) = pphhIndices0
          ELSEIF (hpph == 'ibaj') THEN
             this%ibajRe2Abij%list(blockRe) = pphhIndices0
          ELSEIF (hpph == 'jabi') THEN
             this%jabiRe2Abij%list(blockRe) = pphhIndices0
          ELSEIF (hpph == 'jbai') THEN
             this%jbaiRe2Abij%list(blockRe) = pphhIndices0
          ENDIF
          CALL add2MemoryInt(1)
       ELSE
          
          IF (hpph == 'iabj') THEN
             ! if (blockIndex == 76) then
             !    write(*,*) '.. pphhIndices=',pphhIndices%matr(1,1)
             ! endif
             this%iabjRe2Abij%list(blockRe) = pphhIndices
          ELSEIF (hpph == 'ibaj') THEN
             this%ibajRe2Abij%list(blockRe) = pphhIndices
          ELSEIF (hpph == 'jabi') THEN
             this%jabiRe2Abij%list(blockRe) = pphhIndices
          ELSEIF (hpph == 'jbai') THEN
             this%jbaiRe2Abij%list(blockRe) = pphhIndices
          ENDIF
          CALL add2MemoryInt(nHpPairs*nPhPairs)
       ENDIF
       
       CALL IntMatrix_d(pphhIndices)
       !CALL IntMatrix_d(phPairs)
       !CALL IntMatrix_d(hpPairs)
       
    ENDDO
    
    DEALLOCATE(bra, ket)
    DEALLOCATE(ab, ij)
    DEALLOCATE(abThis, ijThis)
    DEALLOCATE(qNrsA, qNrsB)
    DEALLOCATE(qNrsI, qNrsJ)
    DEALLOCATE(blockQNrsPp)
    DEALLOCATE(blockQNrsHh)
    CALL IntMatrix_d(pphhIndices0)
    
    IF (hpph == 'iabj') THEN
       
       !     Wait until all processes are ready
       CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
       
       !     Send this slice of the array this%index2ReIabj
       !     to all other processes
       ALLOCATE(elementsV(this%basis%mpiJobs%nProcesses))
       ALLOCATE(displacements(this%basis%mpiJobs%nProcesses))
       elementsV = 3*this%elementsProcs
       displacements = 3*this%elementsBeginPs
       
       ! indexBegin = this%elementsBeginPs(&
       !      this%basis%mpiJobs%iAm+1) + 1
       ! indexEnd = this%elementsBeginPs(&
       !      this%basis%mpiJobs%iAm+1) + &
       !      this%elementsProcs(&
       !      this%basis%mpiJobs%iAm+1)
       
       CALL MPI_ALLGATHERV(index2ReIabj, &
            elementsV(this%basis%mpiJobs%iAm+1), &
            MPI_INTEGER, this%index2ReIabj, elementsV, &
            displacements, MPI_INTEGER, MPI_COMM_WORLD, iError)
       
       !!     Wait until all processes are ready
       ! CALL MPI_BARRIER(MPI_COMM_WORLD, iError) 
       ! write(*,*) '== i am ',this%basis%mpiJobs%iAm 
       ! CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
       
       ! write(*,*) '> iam=',this%basis%mpiJobs%iAm,&
       !      'this%index2ReIabj(:,781)=',this%index2ReIabj(:,781)
       ! write(*,*) '> iam=',this%basis%mpiJobs%iAm,&
       !      'this%index2ReIabj(:,20)=',this%index2ReIabj(:,20)
       ! write(*,*) '-- i am ',this%basis%mpiJobs%iAm  
       
       
       DEALLOCATE(elementsV)
       DEALLOCATE(displacements)
    ENDIF
    DEALLOCATE(index2ReIabj)
    
    
  END SUBROUTINE setupHpphRe2Pphh
  
  
END MODULE BlockListReMod
