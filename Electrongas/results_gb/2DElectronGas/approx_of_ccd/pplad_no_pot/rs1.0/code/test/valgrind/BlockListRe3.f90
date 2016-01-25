
!
!     Class for list of blocks where the single-particle
!     states are recoupled.
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
   CONTAINS
     PROCEDURE :: setupPhBlockList
     PROCEDURE :: setupPhBlocks
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
    this%nPhBlocks = this%blocks%nPphhBlocks
    this%nPphhBlocks = this%blocks%nPphhBlocks
    ALLOCATE(this%phList(this%nPhBlocks))
    ALLOCATE(this%pphhList(this%nPphhBlocks))
    
    
  END FUNCTION BlockListRe_
  
  SUBROUTINE setupPhBlockList(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    
    IF (this%basis%nUnoccupied > 0) THEN
       CALL setupPhBlocks(this)
    ENDIF
    
    
  END SUBROUTINE setupPhBlockList
  
  ! !
  ! !     Set up list of (k_rel, m_s_rel) ph blocks
  ! !     obtained from a pp and a hh pair
  ! SUBROUTINE setupPhBlocks(this)
  !   CLASS(BlockListRe), INTENT(INOUT) :: this
  !   TYPE(PairList) :: phPairs
  !   TYPE(IntMatrix) :: pairsMatr
  !   INTEGER :: nPairsPh, ph, thisBlock, thisRank, nOcc
  !   INTEGER :: phThis
  !   LOGICAL :: occupied
  
  !   phPairs = this%phPairs
  !   nPairsPh = phPairs%nPairs

  !   thisRank = phPairs%rank(1)
    
  !   ph = 0
  !   thisBlock = 0
  !   DO WHILE (thisBlock < this%nPhBlocks)
  !      occupied = .FALSE.
       
  !      nOcc = 0
  !      DO WHILE ((ph + 1 <= nPairsPh) &
  !           .AND.(phPairs%rank(ph+1) <= thisRank))
  !         ph = ph + 1
  !         nOcc = nOcc + 1
  !         occupied = .TRUE.
  !      ENDDO
  !      !     If this block is occupied by particle-hole
  !      !     states, store them
  !      IF (occupied) THEN
  !         thisBlock = thisBLock + 1
          
  !         this%phList(thisBlock) = Block1_(this%nQNumbers)
  !         this%phList(thisBlock)%blockQNrs(:) &
  !              = phPairs%blockQNumbers(ph,:)
  !         this%phList(thisBlock)%occupiedPh = .TRUE.
  
  !         WRITE(*,*) ' '
  !         WRITE(*,*) ' Block nr',thisBlock
  !         pairsMatr = IntMatrix_(nOcc, 2)
  !         DO phThis=1, nOcc
  !            pairsMatr%matr(phThis,:) &
  !                 = phPairs%pairs(ph-nOcc+phThis,:)
  !            !WRITE(*,*) phPairs%pairs(ph-nOcc+phThis,:), &
  !            !     phPairs%rank(ph-nOcc+phThis)
  !            WRITE(*,*) phPairs%blockQNumbers(ph-nOcc+phThis,:), &
  !                 phPairs%rank(ph-nOcc+phThis)
  !         ENDDO
  !         this%phList(thisBlock)%phPairs = pairsMatr
  !      ENDIF
       
  !      thisRank = thisRank + 1
       
  !   ENDDO
  !   WRITE(*,*) ' '
  
  
  ! END SUBROUTINE setupPhBlocks
  
  ! !
  ! !     Count the total number of ph blocks obtained
  ! !     fram a pp and a hh pair
  ! !
  ! SUBROUTINE countPhBlocks(this)
  !   CLASS(BlockListRe), INTENT(INOUT) :: this
  !   TYPE(PairList) :: phPairs
  !   INTEGER :: nPairsPh, ph, blocks, thisRank
  !   LOGICAL :: occupied
    
  !   phPairs = this%phPairs
  !   nPairsPh = phPairs%nPairs
  !   write(*,*) 'nPairsPh=',nPairsPh
  
  !   !     Count number of blocks
  !   thisRank = phPairs%rank(1)
    
  !   ph = 0
  !   blocks = 0
  !   DO WHILE (ph < nPairsPh)
  !      occupied = .FALSE.
  !      !write(*,*) 'ph=',ph
  !      DO WHILE ((ph + 1 <= nPairsPh) &
  !           .AND. (phPairs%rank(ph+1) <= thisRank))
  !         ph = ph + 1
  !         occupied = .TRUE.
  !      ENDDO
       
  !      IF (occupied) blocks = blocks + 1
       
  !      thisRank = thisRank + 1
  !   ENDDO
    
  !   this%nPhBlocks = blocks
  !   ALLOCATE(this%phList(blocks))
    
    
  ! END SUBROUTINE countPhBlocks
  
  !
  !     Set up list of particle-hole pairs using
  !     the first states from a pp and a hh pair
  !  
  SUBROUTINE setupPhBlocks(this)
    CLASS(BlockListRe), INTENT(INOUT) :: this
    CLASS(SpBasis), POINTER :: basis
    TYPE(Block1) :: thisBlock
    TYPE(IntMatrix) :: pairsMatr
    INTEGER, ALLOCATABLE :: blockQNumbers(:,:)
    INTEGER, ALLOCATABLE :: pairs(:,:)
    INTEGER :: pair, pphhBlock, blockIndex, nPphhBlocks
    INTEGER :: particle, hole, pp, hh, nPairsThis
    
    nPphhBlocks = this%blocks%nPphhBlocks
    basis => this%basis
    
!    WRITE(*,*) ' Particle-hole pairs:'
!    WRITE(*,*) ' '
!    WRITE(*,*) ' particle   hole   block'
    
    DO pphhBlock=1, nPphhBlocks
       blockIndex = this%blocks%pphhBlocks(pphhBlock)
       
       thisBlock = this%blocks%list(blockIndex)
       
       nPairsThis = thisBlock%ppPairs%dim1 &
            *thisBlock%hhPairs%dim1
       ALLOCATE(pairs(nPairsThis, 2))
       !ALLOCATE(blockQNumbers(nPairsThis, this%nQNumbers))
       
       pair = 0
       DO pp=1, thisBlock%ppPairs%dim1 
          DO hh=1, thisBlock%hhPairs%dim1
             pair = pair + 1
             
             particle = thisBlock%ppPairs%matr(pp, 1)
             hole = thisBlock%hhPairs%matr(hh, 1)
             
             !     Store quantum numbers for the
             !     blocks (CM momentum, total spin)
             ! blockQNumbers(pair,:) = &
             !      basis%spQuantNr(particle,:) &
             !      - basis%spQuantNr(hole,:)
             
             !     Store single-particle orbits
             !     corresponding to a block
             pairs(pair,:) = (/particle, hole/)
             
             !write(*,'(4(x, I6))')  blockQNumbers(pair,:), rank(pair)
!             write(*,'(4(x, I6))') pairs(pair,:), blockIndex
          ENDDO
       ENDDO
       pairsMatr = IntMatrix_(nPairsThis, 2)
       pairsMatr%matr = pairs
       this%phList(pphhBlock)%phPairs = pairsMatr
       
       !DEALLOCATE(blockQNumbers)
       DEALLOCATE(pairs)
    ENDDO
!    write(*,*) ' '    
    
    
  END SUBROUTINE setupPhBlocks
  
  
  
END MODULE BlockListReMod
