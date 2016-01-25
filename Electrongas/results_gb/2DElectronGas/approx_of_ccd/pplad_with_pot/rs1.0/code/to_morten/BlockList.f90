
!
!     Class for list of all (K_CM, M_S) blocks.
!
MODULE BlockListMod
  USE MemoryCounterMod
  USE PairListMod
  USE SpBasisMod
  USE Block1Mod
  IMPLICIT NONE
  
  TYPE :: BlockList
     !
     !     nBlocks:      Number of blocks
     !     nQNumbers:    Number of quantum numbers
     !     hhPairs:      Ordered list of all hole-hole pairs
     !     ppPairs:      Ordered list of all particle-particle
     !                   pairs
     !     phPairs:      Ordered list of all particle-hole 
     !                   pairs
     !     nPphhBlocks:  Number of blocks for both pp and hh 
     !                   pairs
     !     nHhBlocks:    Number of blocks for hh pairs
     !     pphhBlocks:   Vector for inidces in the block list 
     !                   corresponding to states for both pp
     !                   and hh pairs
     !     hhBlocks:     Vector for inidices in the block list
     !                   corresponding to states for hh pairs
     !
     INTEGER :: nBlocks, nQNumbers
     TYPE(Block1), ALLOCATABLE :: list(:)
     TYPE(PairList) :: hhPairs
     TYPE(PairList) :: ppPairs
     TYPE(PairList) :: phPairs
     CLASS(SpBasis), POINTER :: basis
     INTEGER, ALLOCATABLE :: rank2BlockPphh(:)
     INTEGER :: nPphhBlocksAll, nHhBlocks
     INTEGER :: nPphhBlocksThis
     INTEGER :: nPphhElementsBeginThis
     INTEGER, ALLOCATABLE :: pphh2BlocksAll(:)
     INTEGER, ALLOCATABLE :: pphhThis2BlocksAll(:)
     INTEGER, ALLOCATABLE :: pphhThis2PphhAll(:)
     INTEGER, ALLOCATABLE :: hhBlocks(:)
     INTEGER, ALLOCATABLE :: elementsPphhBlocks(:)
     INTEGER, ALLOCATABLE :: pphhElementsProcs(:)
     INTEGER, ALLOCATABLE :: pphhElementsBeginPs(:)
     INTEGER :: maxRank
   CONTAINS
     PROCEDURE :: setupPphhBlocksThis
     PROCEDURE :: setupBlockList
     PROCEDURE :: setupBlocks
     PROCEDURE :: countBlocks
     PROCEDURE :: setupHhPairs
     PROCEDURE :: setupPpPairs
     PROCEDURE :: setupPhPairs
  END TYPE BlockList
  
  
CONTAINS
  
  !
  !     Constructor for BlockList object.
  !
  FUNCTION BlockList_(nQNumbers, basis) RESULT(this)
    INTEGER, INTENT(IN) :: nQNumbers
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    TYPE(BlockList), POINTER :: this
    INTEGER :: nHoles, nParticles
    
    ALLOCATE(this)

    this%basis => basis 
    nHoles = basis%nOccupied
    nParticles = basis%nUnoccupied
    this%nQNumbers = nQNumbers

    this%hhPairs = PairList_(nHoles**2, this%nQNumbers)
    this%ppPairs = PairList_(nParticles**2, this%nQNumbers)
    this%phPairs = PairList_(nParticles*nHoles, this%nQnumbers)
    
    
  END FUNCTION BlockList_
  
  !
  !     Destructor for BlockList object.
  !
  SUBROUTINE BlockList_d(this)
    TYPE(BlockList), TARGET, INTENT(INOUT) :: this
    INTEGER :: thisBlock
    
    IF (ASSOCIATED(this%basis)) THEN
       CALL SpBasis_d(this%basis)
       this%basis => NULL()
    ENDIF
   
    CALL PairList_d(this%hhPairs)
    CALL PairList_d(this%ppPairs)
    CALL PairList_d(this%phPairs)
    
    IF (ALLOCATED(this%list)) THEN
       DO thisBlock=1, this%nBlocks
          CALL Block1_d(this%list(thisBlock))
       ENDDO
       DEALLOCATE(this%list)
    ENDIF
    
    IF (ALLOCATED(this%rank2BlockPphh)) &
         DEALLOCATE(this%rank2BlockPphh)
    IF (ALLOCATED(this%pphh2BlocksAll)) &
         DEALLOCATE(this%pphh2BlocksAll)
    IF (ALLOCATED(this%pphhThis2BlocksAll)) &
         DEALLOCATE(this%pphhThis2BlocksAll)
    IF (ALLOCATED(this%pphhThis2PphhAll)) &
         DEALLOCATE(this%pphhThis2PphhAll)
    IF (ALLOCATED(this%hhBlocks)) &
         DEALLOCATE(this%hhBlocks)
    
    
  END SUBROUTINE BlockList_d
  
  SUBROUTINE setupBlockList(this)
    CLASS(BlockList), TARGET, INTENT(INOUT) :: this
    
    !     Set up pairs of single-particle states
    CALL setupHhPairs(this)
    IF (this%basis%nUnoccupied > 0) THEN
       CALL setupPpPairs(this)
       CALL setupPhPairs(this)
    ENDIF

    IF (this%basis%mpiJobs%iAm == 0) THEN
       WRITE(*,*) ' '
       WRITE(*,*) ' Counting number of blocks...'
       WRITE(*,*) ' '
    ENDIF
    !     Count number of (K_CM, M_S) blocks
    CALL countBlocks(this)
    
    IF (this%basis%mpiJobs%iAm == 0) THEN
       WRITE(*,*) ' '
       WRITE(*,*) ' Setting up (K_CM, M_S) blocks...'
       WRITE(*,*) ' '
    ENDIF
    !     Set up list of blocks
    CALL setupBlocks(this)
    
    !     As part of the MPI parallelization,
    !     assign pphh blocks to the different
    !     processes.
    CALL setupPphhBlocksThis(this)
    
    
  END SUBROUTINE setupBlockList
  
  !
  !     Set up list of (K_CM, M_S) blocks
  !
  SUBROUTINE setupBlocks(this)
    CLASS(BlockList), TARGET, INTENT(INOUT) :: this
    TYPE(PairList) :: hhPairs, ppPairs, phPairs
    TYPE(IntMatrix) :: pairsMatr
    LOGICAL :: occupied, occupiedHh, occupiedPp, occupiedPh
    INTEGER :: nPairsHh, nPairsPp, nPairsPh, nPairsTot
    INTEGER :: thisRank, thisBlock, hh, pp, ph
    INTEGER :: hhThis, ppThis, phThis, nOcc
    INTEGER :: hhBlock, pphhBlock
    INTEGER :: hole1, hole2, particle1, particle2
    
    hhPairs = this%hhPairs
    ppPairs = this%ppPairs
    phPairs = this%phPairs
    nPairsHh = hhPairs%nPairs
    nPairsPp = ppPairs%nPairs
    nPairsPh = phPairs%nPairs
    nPairsTot = nPairsHh + nPairsPp + nPairsPh
    
    thisRank = MIN(hhPairs%rank(1), &
         ppPairs%rank(1), &
         phPairs%rank(1))
    
    hh = 0
    pp = 0
    ph = 0
    thisBlock = 0
    hhBlock = 0
    pphhBlock = 0
    DO WHILE (thisBlock < this%nBlocks)
       occupiedHh = .FALSE.
       occupiedPp = .FALSE.
       occupiedPh = .FALSE.
       occupied = .FALSE.
       
       nOcc = 0
       DO WHILE ((hh + 1 <= nPairsHh) &
            .AND.(hhPairs%rank(hh+1) <= thisRank)) 
          hh = hh + 1
          nOcc = nOcc + 1
          occupiedHh = .TRUE.
       ENDDO
       
       !     If this block is occupied by hole-hole 
       !     states, store them
       IF (occupiedHh) THEN
          IF (.NOT.occupied) THEN
             occupied = .TRUE.
             thisBlock = thisBlock + 1
             
             this%list(thisBlock) = Block1_(this%nQNumbers)
             this%list(thisBlock)%blockQNrs(:) &
                  = hhPairs%blockQNumbers(hh,:)
             this%list(thisBlock)%rank = thisRank
             ! WRITE(*,*) ' '
             ! WRITE(*,*) ' Block nr',thisBlock
          ENDIF
          this%list(thisBlock)%occupiedHh = .TRUE.
          
!          WRITE(*,*) ' Hole-hole pairs:'
          pairsMatr = IntMatrix_(nOcc, 2)
          
          DO hhThis=1, nOcc
             pairsMatr%matr(hhThis, :) &
                  = hhPairs%pairs(hh-nOcc+hhThis,:)
             
             hole1 = hhPairs%pairs(hh-nOcc+hhThis, 1)
             hole2 = hhPairs%pairs(hh-nOcc+hhThis, 2)
             
             ! WRITE(*,*) hhPairs%pairs(hh-nOcc+hhThis,:), &
             !      hhPairs%rank(hh-nOcc+hhThis)
             !WRITE(*,*) hhPairs%blockQNumbers(hh-nOcc+hhThis,:), &
             !    hhPairs%rank(hh-nOcc+hhThis)
          ENDDO
          this%list(thisBlock)%hhPairs = pairsMatr
          
          CALL add2MemoryInt(2*nOcc)
          CALL IntMatrix_d(pairsMatr)
        ENDIF
     
       nOcc = 0
       DO WHILE ((pp + 1 <= nPairsPp) &
            .AND.(ppPairs%rank(pp+1) <= thisRank))
          pp = pp + 1
          nOcc = nOcc + 1
          occupiedPp = .TRUE.
       ENDDO
       
       !     If this block is occupied by particle-particle
       !     states, store them
       IF (occupiedPp) THEN
          IF (.NOT.occupied) THEN
             occupied = .TRUE.
             thisBlock = thisBlock + 1
             
             this%list(thisBlock) = Block1_(this%nQNumbers)
             this%list(thisBlock)%blockQNrs(:) &
                  = ppPairs%blockQNumbers(pp,:)
             this%list(thisBlock)%rank = thisRank
 !           WRITE(*,*) ' '
 !            WRITE(*,*) ' Block nr',thisBlock
          ENDIF
          this%list(thisBlock)%occupiedPp = .TRUE.
          
 !          WRITE(*,*) ' Particle-particle pairs:'
          pairsMatr = IntMatrix_(nOcc, 2)
      
          DO ppThis=1, nOcc
             pairsMatr%matr(ppThis, :) &
                  = ppPairs%pairs(pp-nOcc+ppThis,:)
             particle1 = ppPairs%pairs(pp-nOcc+ppThis, 1)
             particle2 = ppPairs%pairs(pp-nOcc+ppThis, 2)
             
             ! WRITE(*,*) ppPairs%pairs(pp-nOcc+ppThis,:), &
             !      ppPairs%rank(pp-nOcc+ppThis)
             !WRITE(*,*) ppPairs%blockQNumbers(pp-nOcc+ppThis,:), &
             !    ppPairs%rank(pp-nOcc+ppThis)
          ENDDO
          
          this%list(thisBlock)%ppPairs = pairsMatr
          
          CALL add2MemoryInt(2*nOcc)
          CALL IntMatrix_d(pairsMatr)
       ENDIF
       
       nOcc = 0
       DO WHILE ((ph + 1 <= nPairsPh) &
            .AND.(phPairs%rank(ph+1) <= thisRank))
          ph = ph + 1
          nOcc = nOcc + 1
          occupiedPh = .TRUE.
       ENDDO
       !     If this block is occupied by particle-hole
       !     states, store them
       IF (occupiedPh) THEN
          IF (.NOT.occupied) THEN
             occupied = .TRUE.
             thisBlock = thisBlock + 1
             
             this%list(thisBlock) = Block1_(this%nQNumbers)
             this%list(thisBlock)%blockQNrs(:) &
                  = phPairs%blockQNumbers(ph,:)
             this%list(thisBlock)%rank = thisRank
!             WRITE(*,*) ' '
!             WRITE(*,*) ' Block nr',thisBlock
          ENDIF
          this%list(thisBlock)%occupiedPh = .TRUE.
          
!          WRITE(*,*) ' Particle-hole pairs:'
          pairsMatr = IntMatrix_(nOcc, 2)
          
          DO phThis=1, nOcc
             pairsMatr%matr(phThis,:) &
                  = phPairs%pairs(ph-nOcc+phThis,:)
             
             particle1 = phPairs%pairs(ph-nOcc+phThis, 1)
             hole1 = phPairs%pairs(ph-nOcc+phThis, 2)
             
             !WRITE(*,*) phPairs%pairs(ph-nOcc+phThis,:), &
             !     phPairs%rank(ph-nOcc+phThis)
             !WRITE(*,*) phPairs%blockQNumbers(ph-nOcc+phThis,:), &
             !    phPairs%rank(ph-nOcc+phThis)
          ENDDO
          this%list(thisBlock)%phPairs = pairsMatr
         
          CALL add2MemoryInt(2*nOcc)
          CALL IntMatrix_d(pairsMatr)
       ENDIF
       
       !     List of blocks that are occupied by
       !     both hh and pp states
       IF (occupiedHh .AND. occupiedPp) THEN
          pphhBlock = pphhBlock + 1
          this%pphh2BlocksAll(pphhBlock) = thisBlock
       ENDIF
       !     List of blocks that are occupied by
       !     hh states
       IF (occupiedHh) THEN
          hhBlock = hhBlock + 1
          this%hhBlocks(hhBlock) = thisBlock
       ENDIF
       
       thisRank = thisRank + 1
       
     ENDDO
!    WRITE(*,*) ' '
    
    
  END SUBROUTINE setupBlocks
  
  !
  !     Count the total number of blocks
  !
  SUBROUTINE countBlocks(this)
    CLASS(BlockList), TARGET, INTENT(INOUT) :: this
    TYPE(PairList) :: hhPairs, ppPairs, phPairs
    INTEGER :: pairsTot, thisRank, hh, pp, ph, blocks
    INTEGER :: nPairsHh, nPairsPp, nPairsPh, nPairsTot
    INTEGER :: thisBlock, nOccPphh, nOccHh, maxRank
    LOGICAL :: occupied, occupiedHh, occupiedPp, occupiedPh
    
    hhPairs = this%hhPairs
    ppPairs = this%ppPairs
    phPairs = this%phPairs
    nPairsHh = hhPairs%nPairs
    nPairsPp = ppPairs%nPairs
    nPairsPh = phPairs%nPairs
    nPairsTot = nPairsHh + nPairsPp + nPairsPh
    
    IF (this%basis%mpiJobs%iAm == 0) THEN
       WRITE(*,*) ' Number of hh pairs:',nPairsHh
       WRITE(*,*) ' Number of pp pairs:',nPairsPp
       WRITE(*,*) ' Number of ph pairs:',nPairsPh
    ENDIF
    
    maxRank = MAX(hhPairs%rank(nPairsHh), &
         ppPairs%rank(nPairsPp), &
         phPairs%rank(nPairsPh))
    ALLOCATE(this%rank2BlockPphh(maxRank))
    this%maxRank = maxRank
    this%rank2BlockPphh = 0
    
    !     Count number of blocks
    thisRank = MIN(hhPairs%rank(1), &
         ppPairs%rank(1), &
         phPairs%rank(1))
    
    hh = 0
    pp = 0
    ph = 0
    blocks = 0
    nOccHh = 0
    nOccPphh = 0
    pairsTot = 0
    DO WHILE (pairsTot < nPairsTot)
       occupied = .FALSE.
       occupiedHh = .FALSE.
       occupiedPp = .FALSE.
       occupiedPh = .FALSE.
       DO WHILE ((hh + 1 <= nPairsHh) &
            .AND. (hhPairs%rank(hh+1) <= thisRank)) 
          hh = hh + 1
          occupiedHh = .TRUE.
       ENDDO
       
       DO WHILE ((pp + 1 <= nPairsPp) &
            .AND. (ppPairs%rank(pp+1) <= thisRank))
          pp = pp + 1
          occupiedPp = .TRUE.
       ENDDO
       
       DO WHILE ((ph + 1 <= nPairsPh) &
            .AND. (phPairs%rank(ph+1) <= thisRank))
          ph = ph + 1
          occupiedPh = .TRUE.
       ENDDO
       
       occupied = (occupiedHh .OR. occupiedPp &
            .OR. occupiedPh)
       IF (occupied) THEN
          blocks = blocks + 1
       ENDIF
       
       IF (occupiedHh .AND. occupiedPp) THEN
          nOccPphh = nOccPphh + 1
          this%rank2BlockPphh(thisRank) = nOccPphh
       ENDIF
       IF (occupiedHh) nOccHh = nOccHh + 1
       
       pairsTot = hh + pp + ph 
       thisRank = thisRank + 1
    ENDDO
    
    !     Allocate block list
    this%nBlocks = blocks
    ALLOCATE(this%list(this%nBlocks))
    
    this%nPphhBlocksAll = nOccPphh
    this%nHhBlocks = nOccHh
    ALLOCATE(this%pphh2BlocksAll(this%nPphhBlocksAll))
    ALLOCATE(this%hhBlocks(this%nHhBlocks))
    
    IF (this%basis%mpiJobs%iAm == 0) THEN
       WRITE(*,*) ' Numbers of blocks: ', blocks
       WRITE(*,*) ' Number of pphh blocks: ', nOccPphh
    ENDIF
    
    
  END SUBROUTINE countBlocks
  
  !
  !     Divide the blocks with pp and hh states to the
  !     different MPI processes.
  !
  SUBROUTINE setupPphhBlocksThis(this)
    CLASS(BlockList), TARGET, INTENT(INOUT) :: this
    CLASS(SpBasis), POINTER :: basis
    INTEGER(KIND=8), ALLOCATABLE :: elementsPpppBlocks(:)
    INTEGER, ALLOCATABLE :: elementsPphhBlocks(:)
    TYPE(Block1) :: thisBlock
    TYPE(IntMatrix) :: ppPairs, hhPairs
    INTEGER(KIND=8) :: elementsPerProc, elementsPpppThis
    INTEGER(KIND=8) :: elementsPpppAll, elementsPpppTot
    INTEGER :: blockIndex, nPpPairs, nHhPairs
    INTEGER :: pphhBlock, nProcesses, iAm
    INTEGER :: process, nBlocksThis, pphhBlockBegin
    INTEGER :: elementsPphhThis, pphhBlockEnd 
    INTEGER :: elementsPphhThisSave, element, iError
    INTEGER :: nPphhElementsBegin, elementsPphhTot
    INCLUDE 'mpif.h'
    
    basis => this%basis
    ALLOCATE(elementsPpppBlocks(this%nPphhBlocksAll))
    ALLOCATE(elementsPphhBlocks(this%nPphhBLocksAll))
    
    ! IF (basis%mpiJobs%iAm == 0) THEN
    !    WRITE(*,*) ' '
    !    WRITE(*,*) 'pphhBlock     elements'
    ! ENDIF
    !     Count pppp elements in the different blocks
    elementsPpppTot = 0
    elementsPphhTOt = 0
    DO pphhBlock=1, this%nPphhBlocksAll
       blockIndex = this%pphh2BlocksAll(pphhBlock)
       
       thisBlock = this%list(blockIndex)
       
       ppPairs = thisBlock%ppPairs
       hhPairs = thisBlock%hhPairs
       nPpPairs = ppPairs%dim1
       nHhPairs = hhPairs%dim1
       
       !     Number of elements in this block
       elementsPpppBlocks(pphhBlock) = nPpPairs**2
       elementsPphhBlocks(pphhBlock) = nPpPairs*nHhPairs
       ! IF (basis%mpiJobs%iAm == 0) THEN
       !    WRITE(*,*) pphhBlock, nPpPairs**2
       ! ENDIF
       
       elementsPpppTot = elementsPpppTot + nPpPairs**2
       elementsPphhTot = elementsPphhTot + nPpPairs*nHhPairs
    ENDDO
 !   CALL basis%mpiJobs%writeMpi(' ')
    
    ALLOCATE(this%elementsPphhBlocks(this%nPphhBlocksAll))
    this%elementsPphhBlocks = elementsPphhBlocks
    
    !     Number of processes
    nProcesses = basis%mpiJobs%nProcesses
    !     Approximate number of elements per process
    elementsPerProc = &
         CEILING(1.d0*elementsPpppTot/nProcesses)
    
    ! IF (basis%mpiJobs%iAm == 0) THEN
    !    WRITE(*,*) 'Total number of elements: ', elementsPpppTot
    !    WRITE(*,*) 'Elements per process: ', elementsPerProc
    ! ENDIF
    
    process = 0
    elementsPpppThis = 0
    elementsPpppAll = 0
    elementsPphhThis = 0
    nPphhElementsBegin = 0
    nBlocksThis = 0
    pphhBlockBegin = 1
    pphhBlock = 0
    DO WHILE ((process <= basis%mpiJobs%iAm).AND. &
         (pphhBlock <= this%nPphhBlocksAll)) 
       pphhBlock = pphhBlock + 1
       
       !     Number of pppp elements assigned to
       !     this process
       elementsPpppThis = elementsPpppThis &
            + elementsPpppBlocks(pphhBlock)
       elementsPpppAll = elementsPpppAll &
            + elementsPpppBlocks(pphhBlock)
       elementsPphhThis = elementsPphhThis &
            + elementsPphhBlocks(pphhBlock)
       nPphhElementsBegin = nPphhElementsBegin &
            + elementsPphhBlocks(pphhBlock)
       
       IF (process == basis%mpiJobs%iAm) THEN
          nBlocksThis = nBlocksThis + 1
       ENDIF
       
       IF ((elementsPpppThis >= elementsPerProc) &
            .OR.(elementsPpppAll == elementsPpppTot)) THEN
          process = process + 1
          elementsPphhThisSave = elementsPphhThis
          elementsPpppThis = 0
          elementsPphhThis = 0
          this%nPphhElementsBeginThis = nPphhElementsBegin &
               - elementsPphhThisSave
          IF (process - 1 == basis%mpiJobs%iAm) &
               pphhBlockEnd = pphhBlock
          IF (process == basis%mpiJobs%iAm) &
               pphhBlockBegin = pphhBlock + 1
       ENDIF
       
    ENDDO
    
    this%nPphhBlocksThis = nBlocksThis
    ALLOCATE(this%pphhThis2BlocksAll(nBlocksThis))
    ALLOCATE(this%pphhThis2PphhAll(nBlocksThis))
    
    ! WRITE(*,*) 'iAm:',basis%mpiJobs%iAm, &
    !      ',begin:',pphhBlockBegin, ',end:',pphhBlockEnd
    ! WRITE(*,*) 'iAm:',basis%mpiJobs%iAm, &
    !      ',elements:',elementsPphhThisSave
    ! WRITE(*,*) 'nBlocksThis=',nBlocksThis,&
    !      ',=',pphhBlockEnd-pphhBlockBegin+1
    ! write(*,*) 'iAm:',basis%mpiJobs%iAm, &
    !      ', nEleBegin=',this%nPphhElementsBeginThis
    ! write(*,*) 'elementsPphhTot=',elementsPphhTot 
    
    element = 0
    DO pphhBlock=pphhBlockBegin, pphhBlockEnd
       element = element + 1
       
       blockIndex = this%pphh2BlocksAll(pphhBlock)
       this%pphhThis2BlocksAll(element) = blockIndex
       
       this%pphhThis2PphhAll(element) = pphhBlock
    ENDDO
    
    ALLOCATE(this%pphhElementsProcs(basis%mpiJobs%nProcesses))
    ALLOCATE(this%pphhElementsBeginPs(basis%mpiJobs%nProcesses))
    
    CALL MPI_ALLGATHER(elementsPphhThisSave, 1, &
         MPI_INTEGER, this%pphhElementsProcs, 1, &
         MPI_INTEGER, MPI_COMM_WORLD, iError)
    
    this%pphhElementsBeginPs(1) = 0
    DO process=2, basis%mpiJobs%nProcesses
       this%pphhElementsBeginPs(process) = &
            this%pphhElementsBeginPs(process-1) &
            + this%pphhElementsProcs(process-1)
    ENDDO

    ! write(*,*) 'this%pphhElementsProcs=',&
    !      this%pphhElementsProcs
    ! write(*,*) 'this%pphhElementsBeginPs=',&
    !      this%pphhElementsBeginPs
    
    DEALLOCATE(elementsPpppBlocks)
    DEALLOCATE(elementsPphhBlocks)
    
    
  END SUBROUTINE setupPphhBlocksThis
  
  
  !
  !     Set up list of hole-hole pairs
  !
  SUBROUTINE setupHhPairs(this)
    CLASS(BlockList), TARGET, INTENT(INOUT) :: this
    CLASS(SpBasis), POINTER :: basis
    INTEGER, ALLOCATABLE :: blockQNumbers(:,:)
    INTEGER, ALLOCATABLE :: blockQNumbersNew(:,:)
    INTEGER, ALLOCATABLE :: pairs(:,:)
    INTEGER, ALLOCATABLE :: pairsNew(:,:)
    INTEGER, ALLOCATABLE :: rank(:), rankNew(:)
    INTEGER, ALLOCATABLE :: order(:)
    INTEGER :: nHoles, nPairs, hole1, hole2
    INTEGER :: pair, i, newPair
    
    basis => this%basis
    nHoles = basis%nOccupied
    nPairs = nHoles**2
    ALLOCATE(blockQNumbers(nPairs, this%nQNumbers))
    ALLOCATE(blockQNumbersNew(nPairs, this%nQNumbers))
    ALLOCATE(pairs(nPairs, 2))
    ALLOCATE(pairsNew(nPairs, 2))
    ALLOCATE(rank(nPairs), rankNew(nPairs))
    ALLOCATE(order(nPairs))
    
!    WRITE(*,*) ' '
    pair = 0
    DO hole1=1, nHoles
       DO hole2=1, nHoles
          pair = pair + 1
          
          !     Store quantum numbers for the
          !     blocks (CM momentum, total spin)
          blockQNumbers(pair,:) &
               = basis%spQuantNr(hole1,:) &
               + basis%spQuantNr(hole2,:)
          
          !     Store single-particle orbits
          !     corresponding to a block
          pairs(pair,:) = (/hole1, hole2/)
          
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
    
!    WRITE(*,*) ' Hole-hole pairs:'
!    WRITE(*,*) ' '
!    WRITE(*,*) ' hole1   hole2   ranking number'
    !     Reorder the pairs
    DO pair=1, nPairs
       newPair = order(pair)
       blockQNumbersNew(pair,:) = blockQNumbers(newPair,:)
!       write(*,'(4(2x, I6))')  blockQNumbers(newPair,:), rank(newPair)
!       write(*,'(3(x, I6))')  pairs(newPair,:), rank(newPair)
       pairsNew(pair,:) = pairs(newPair,:)
       rankNew(pair) = rank(newPair)
    ENDDO
    this%hhPairs%blockQNumbers = blockQNumbersNew
    this%hhPairs%pairs = pairsNew
    this%hhPairs%rank = rankNew
!    WRITE(*,*) ' '
    
    DEALLOCATE(blockQNumbers, blockQNumbersNew)
    DEALLOCATE(pairs, pairsNew)
    DEALLOCATE(rank, rankNew, order)
    basis => NULL()
    
    
  END SUBROUTINE setupHhPairs
      
  !
  !     Set up list of particle-particle pairs
  !
  SUBROUTINE setupPpPairs(this)
    CLASS(BlockList), TARGET, INTENT(INOUT) :: this
    CLASS(SpBasis), POINTER :: basis
    INTEGER, ALLOCATABLE :: blockQNumbers(:,:)
    INTEGER, ALLOCATABLE :: blockQNumbersNew(:,:)
    INTEGER, ALLOCATABLE :: pairs(:,:)
    INTEGER, ALLOCATABLE :: pairsNew(:,:)
    INTEGER, ALLOCATABLE :: rank(:), rankNew(:)
    INTEGER, ALLOCATABLE :: order(:)
    INTEGER :: nHoles, nParticles, nPairs
    INTEGER :: particle1, particle2, pair, i
    INTEGER :: newPair
    
    basis => this%basis
    nHoles = basis%nOccupied
    nParticles = basis%nUnoccupied
    nPairs = nParticles**2
    ALLOCATE(blockQNumbers(nPairs, this%nQNumbers))
    ALLOCATE(blockQNumbersNew(nPairs, this%nQNumbers))
    ALLOCATE(pairs(nPairs, 2))
    ALLOCATE(pairsNew(nPairs, 2))
    ALLOCATE(rank(nPairs), rankNew(nPairs))
    ALLOCATE(order(nPairs))
    
!    WRITE(*,*) ' '
    pair = 0
    DO particle1=nHoles+1, nHoles+nParticles
       DO particle2=nHoles+1, nHoles+nParticles
          pair = pair + 1
          
          !     Store quantum numbers for the
          !     blocks (CM momentum, total spin)
          blockQNumbers(pair,:) &
               = basis%spQuantNr(particle1,:) &
               + basis%spQuantNr(particle2,:)
          
          !     Store single-particle orbits
          !     corresponding to a block
          pairs(pair,:) = (/particle1, particle2/)
          
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

!    WRITE(*,*) ' Particle-particle pairs:'
!    WRITE(*,*) ' '
!    WRITE(*,*) ' particle1   particle2   ranking number'
    !     Reorder the pairs
    DO pair=1, nPairs
       newPair = order(pair)
       blockQNumbersNew(pair,:) = blockQNumbers(newPair,:)
       pairsNew(pair,:) = pairs(newPair,:)
       rankNew(pair) = rank(newPair)
       !write(*,'(4(x, I6))')  blockQNumbers(newPair,:), rank(newPair)
!       write(*,'(3(x, I6))')  pairs(newPair,:), rank(newPair)
    ENDDO
    this%ppPairs%blockQNumbers = blockQNumbersNew
    this%ppPairs%pairs = pairsNew
    this%ppPairs%rank = rankNew
!    WRITE(*,*) ' '

    DEALLOCATE(blockQNumbers, blockQNumbersNew)
    DEALLOCATE(pairs, pairsNew)
    DEALLOCATE(rank, rankNew, order)
    
    
  END SUBROUTINE setupPpPairs
  
  !
  !     Set up list of particle-hole pairs
  !
  SUBROUTINE setupPhPairs(this)
    CLASS(BlockList), TARGET, INTENT(INOUT) :: this
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
    
    write(*,*) ' '
    pair = 0
    DO hole=1, nHoles
       DO particle=nHoles+1, nHoles+nParticles
          pair = pair + 1
          
          !     Store quantum numbers for the
          !     blocks (CM momentum, total spin)
          blockQNumbers(pair,:) &
               = basis%spQuantNr(particle,:) &
               + basis%spQuantNr(hole,:)
          
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
    
!    WRITE(*,*) ' Particle-hole pairs:'
!    WRITE(*,*) ' '
!    WRITE(*,*) ' particle   hole   ranking number'
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
    
    DEALLOCATE(blockQNumbers, blockQNumbersNew)
    DEALLOCATE(pairs, pairsNew)
    DEALLOCATE(rank, rankNew, order)
    
    
  END SUBROUTINE setupPhPairs
  
  !
  !     Get a index corresponding to a given 
  !     set of quantum numbers
  !
  FUNCTION getRank(basis, blockQNumbers, nQNumbers) &
       RESULT(rank)
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    INTEGER, INTENT(IN) :: nQNumbers
    INTEGER, INTENT(IN) :: blockQNumbers(nQNumbers)
    INTEGER :: rank
    
    IF (nQNumbers == 3) THEN
       rank = (4*basis%nMax + 1)**2 &
            *(2 + blockQNumbers(3)) &
            + (4*basis%nMax + 1) &
            *(2*basis%nMax + blockQNumbers(2)) &
            + 2*basis%nMax + 1 + blockQNumbers(1)
    ELSEIF (nQNumbers == 4) THEN
       rank = (4*basis%nMax + 1)**3 &
            *(2 + blockQNumbers(4)) &
            + (4*basis%nMax + 1)**2 &
            *(2*basis%nMax + blockQNumbers(3)) &
            + (4*basis%nMax + 1) &
            *(2*basis%nMax + blockQNumbers(2)) &
            + 2*basis%nMax + 1 + blockQNumbers(1)
    ELSE
       WRITE(*,*) 'rank setup only for 2d and 3d'
       STOP
    ENDIF
    
    
  END FUNCTION getRank
  
  
END MODULE BlockListMod
