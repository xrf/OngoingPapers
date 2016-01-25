
!
!     Class for the Hamiltonian operator.
!
MODULE HamiltonianMod
  USE RealMatrixMod
  USE IntMatrixMod
  USE BlockListReMod
  USE BlockListMod
  USE ConstantsMod
  USE SpBasisMod
  IMPLICIT NONE
  
  TYPE, ABSTRACT :: Hamiltonian
     !
     !     basisSp:           Single-particle basis
     !     vHhhh etc.:       Interaction matrix lists
     !     vPhphRe etc.:    Matrix lists for recoupled
     !                        interaction matrices
     !
     CLASS(SpBasis), POINTER :: basisSp   
     TYPE(BlockList), POINTER :: blocks
     TYPE(BlockListRe), POINTER :: blocksRe
     TYPE(RealMatrix), ALLOCATABLE :: vHhhh(:)
     TYPE(RealMatrix), ALLOCATABLE :: vHhhhRef(:)
     TYPE(RealMatrix), ALLOCATABLE :: vHhpp(:)
     TYPE(RealMatrix), ALLOCATABLE :: vPphh(:)
     TYPE(RealMatrix), ALLOCATABLE :: vPppp(:)    
     TYPE(RealMatrix), ALLOCATABLE :: vPphhRe(:)
     TYPE(RealMatrix), ALLOCATABLE :: vPhphRe(:)
     
   CONTAINS
     PROCEDURE(getInteract), &
          DEFERRED :: getInteraction
     PROCEDURE :: setupInteractionMatr
     PROCEDURE :: setupInteractionMatrRe
  END TYPE Hamiltonian
  
  ABSTRACT INTERFACE      
     FUNCTION getInteract(this, orbit1, orbit2, orbit3, orbit4) &
          RESULT(interaction)
       USE ConstantsMod
       IMPORT Hamiltonian
       CLASS(Hamiltonian), TARGET, INTENT(IN) :: this
       INTEGER, INTENT(IN) :: orbit1, orbit2, orbit3, orbit4
       REAL(DP) :: interaction
     END FUNCTION getInteract
  END INTERFACE
  
CONTAINS
  
  !
  !     Routine to set up lists of cross-coupled 
  !     interaction matrices in the given single-particle
  !     basis.
  !
  SUBROUTINE setupInteractionMatrRe(this)
    CLASS(Hamiltonian), TARGET, INTENT(INOUT) :: this
    TYPE(RealMatrix) :: vPphhReMatrix, vPhphReMatrix
    TYPE(Block1) :: thisBlock
    TYPE(IntMatrix) :: phPairs, hpPairs
    INTEGER, ALLOCATABLE :: bra(:), ket(:)
    INTEGER :: pair1, pair2, blockIndex, nPhPairs, nHpPairs
    REAL(DP) :: interactionPphhRe, interactionPhphRe
    INTEGER :: blockRe
    
    ALLOCATE(bra(2), ket(2))
    
    DO blockRe=1, this%blocksRe%nBlocksReThis
       blockIndex = this%blocksRe%blocksThis2BlocksAll(blockRe)
       
       thisBlock = this%blocksRe%list(blockIndex)
       phPairs = thisBlock%phPairs
       hpPairs = thisBlock%hpPairs
       nPhPairs = phPairs%dim1
       nHpPairs = hpPairs%dim1
       
       vPphhReMatrix = RealMatrix_(nHpPairs, nPhPairs)
       vPhphReMatrix = RealMatrix_(nPhPairs, nPhPairs)
       
       DO pair1=1, nHpPairs
          bra(:) = hpPairs%matr(pair1,:)
          DO pair2=1, nPhPairs
             ket(:) = phPairs%matr(pair2,:)
             
             !     Cross-coupled abij matrix
             !       => iabj matrix
             !
             !     Get antisymmetrized interaction
             !     matrix element
             interactionPphhRe = this%getInteraction(&
                  bra(2), ket(1), bra(1), ket(2))
             vPphhReMatrix%matr(pair1, pair2) = interactionPphhRe
             
          ENDDO
       ENDDO
       
       DO pair1=1, nPhPairs
          bra(:) = phPairs%matr(pair1,:)
          DO pair2=1, nPhPairs
             ket(:) = phPairs%matr(pair2,:)
             
             !     Cross-coupled aibj matrix
             !       => ajbi matrix
             interactionPhphRe = this%getInteraction(&
                  ket(1), bra(2), bra(1), ket(2))
             vPhphReMatrix%matr(pair1, pair2) = interactionPhphRe
             
          ENDDO
       ENDDO
       
       this%vPphhRe(blockRe) = vPphhReMatrix
       this%vPhphRe(blockRe) = vPhphReMatrix
       
       CALL add2MemoryR8((nHpPairs+nPhPairs)*nPhPairs)
       CALL RealMatrix_d(vPphhReMatrix)
       CALL RealMatrix_d(vPhphReMatrix)
    ENDDO
    
    DEALLOCATE(bra, ket)
    
    
  END SUBROUTINE setupInteractionMatrRe
  
  !
  !     Routine to set up interaction matrix lists 
  !     in the given single-particle basis.
  !
  SUBROUTINE setupInteractionMatr(this)
    CLASS(Hamiltonian), TARGET, INTENT(INOUT) :: this
    INTEGER, ALLOCATABLE :: bra(:), ket(:) 
    TYPE(Block1) :: thisBlock
    TYPE(RealMatrix) :: vHhhhMatrix!, vHhhpMatrix
    TYPE(RealMatrix) :: vHhppMatrix!, vHpppMatrix
    TYPE(RealMatrix) :: vPpppMatrix, vPhPhMatrix
    TYPE(RealMatrix) :: vPphhMatrix, vMatrix0
    TYPE(IntMatrix) :: hhPairs, ppPairs, phPairs
    INTEGER :: nHoles, nParticles, blockIndex
    INTEGER :: nHhPairs, nPpPairs, nPhPairs
    INTEGER :: pair1, pair2, pphhBlock
    REAL(DP) :: interaction
    INTEGER :: memoryVPppp
    
    ALLOCATE(bra(2), ket(2))
    memoryVPppp = 0
    
    nHoles = this%basisSp%nOccupied
    nParticles = this%basisSp%nUnoccupied
    
!    WRITE(*,*) ' '
    !     Loop over (K_CM, M_S) pphh blocks
    DO pphhBlock=1, this%blocks%nPphhBlocksThis
       
       blockIndex = this%blocks%pphhThis2BlocksAll(pphhBlock)
       
!       WRITE(*,*) ' Block number', blockIndex
       thisBlock = this%blocks%list(blockIndex)
       !     hh, pp, and, ph pairs
       hhPairs = thisBlock%hhPairs
       ppPairs = thisBlock%ppPairs
       phPairs = thisBlock%phPairs
       nHhPairs = hhPairs%dim1
       nPpPairs = ppPairs%dim1
       nPhPairs = phPairs%dim1
       
       !          WRITE(*,*) 'hh pairs:', nHhPairs
       !     Hole-hole hole-hole matrix
       vHhhhMatrix = RealMatrix_(nHhPairs, nHhPairs)
       DO pair1=1, nHhPairs
          bra(:) = hhPairs%matr(pair1,:)
          DO pair2=1, nHhPairs
             ket(:) = hhPairs%matr(pair2,:)
             
             !     Get antisymmetrized interaction
             !     matrix element
             interaction = this%getInteraction(&
                  bra(1), bra(2), ket(1), ket(2))
             
             vHhhhMatrix%matr(pair1, pair2) = interaction
          ENDDO
       ENDDO
       this%vHhhh(pphhBlock) = vHhhhMatrix
       
       CALL add2MemoryR8(nHhPairs**2) 
       CALL RealMatrix_d(vHhhhMatrix)
       
       
       !         WRITE(*,*) 'hh pairs:', nHhPairs,', pp pairs:', nPpPairs
       !     Hole-hole particle-particle matrix
       vHhppMatrix = RealMatrix_(nHhPairs, nPpPairs)
       vPphhMatrix = RealMatrix_(nPpPairs, nHhPairs)
       
       DO pair1=1, nHhPairs
          bra(:) = HhPairs%matr(pair1,:)
          DO pair2=1, nPpPairs
             ket(:) = PpPairs%matr(pair2,:)
             
             !     Get antisymmetrized interaction
             !     matrix element
             interaction = this%getInteraction(&
                  bra(1), bra(2), ket(1), ket(2))
             
             vHhppMatrix%matr(pair1, pair2) = interaction
             vPphhMatrix%matr(pair2, pair1) = interaction
             
          ENDDO
       ENDDO
       this%vHhpp(pphhBlock) = vHhppMatrix
       this%vPphh(pphhBlock) = vPphhMatrix
       
       CALL add2MemoryR8(2*nHhPairs*nPpPairs)
       CALL RealMatrix_d(vHhppMatrix)
       CALL RealMatrix_d(vPphhMatrix)
       
       !          WRITE(*,*) 'pp pairs:', nPpPairs
       !     Particle-particle particle-particle matrix
       vPpppMatrix = RealMatrix_(nPpPairs, nPpPairs)
       
       DO pair1=1, nPpPairs
          bra(:) = PpPairs%matr(pair1,:)
          DO pair2=1, nPpPairs
             ket(:) = PpPairs%matr(pair2,:)
             
             !     Get antisymmetrized interaction
             !     matrix element
             interaction = this%getInteraction(&
                  bra(1), bra(2), ket(1), ket(2))
             
             vPpppMatrix%matr(pair1, pair2) = interaction
             
          ENDDO
       ENDDO
       this%vPppp(pphhBlock) = vPpppMatrix
       
       CALL add2MemoryR8(nPpPairs**2)
       memoryVPppp = memoryVPppp + nPpPairs**2*sizeReal8
       CALL RealMatrix_d(vPpppMatrix)
       
    ENDDO
    
    vMatrix0 = RealMatrix_(1, 1)
    vMatrix0%dim1 = 0
    !     Loop over (K_CM, M_S) blocks
    DO blockIndex=1, this%blocks%nBlocks
       thisBlock = this%blocks%list(blockIndex)
       
       hhPairs = thisBlock%hhPairs
       nHhPairs = hhPairs%dim1
       
       IF (thisBlock%occupiedHh) THEN
          
          !     Hole-hole hole-hole matrix
          vHhhhMatrix = RealMatrix_(nHhPairs, nHhPairs)
          DO pair1=1, nHhPairs
             bra(:) = hhPairs%matr(pair1,:)
             DO pair2=1, nHhPairs
                ket(:) = hhPairs%matr(pair2,:)
                
                !     Get antisymmetrized interaction
                !     matrix element
                interaction = this%getInteraction(&
                     bra(1), bra(2), ket(1), ket(2))
                
                vHhhhMatrix%matr(pair1, pair2) = interaction
             ENDDO
          ENDDO
          this%vHhhhRef(blockIndex) = vHhhhMatrix
          
          CALL add2MemoryR8(nHhPairs**2) 
          CALL RealMatrix_d(vHhhhMatrix)
       ELSE
          this%vHhhhRef(blockIndex) = vMatrix0
          CALL add2MemoryR8(1)
          
       ENDIF
    ENDDO
    CALL RealMatrix_d(vMatrix0)
    
!    WRITE(*,*) ' '
    DEALLOCATE(bra, ket)
    write(*,*) 'Memory of vPppp: ', memoryVPppp, ' Bytes'
    
    
  END SUBROUTINE setupInteractionMatr
  
  !
  !     Destructor for polymorphic Hamiltonian object
  !
  SUBROUTINE Hamiltonian_d(this)
    CLASS(Hamiltonian), TARGET, INTENT(INOUT) :: this
    INTEGER :: thisBlock, blockRe

    IF (ALLOCATED(this%vHhhh)) THEN
       DO thisBlock=1, this%blocks%nPphhBlocksThis
          CALL RealMatrix_d(this%vHhhh(thisBlock))
       ENDDO
       DEALLOCATE(this%vHhhh)
    ENDIF
    
    IF (ALLOCATED(this%vHhhhRef)) THEN
       DO thisBlock=1, this%blocks%nBlocks
          CALL RealMatrix_d(this%vHhhhRef(thisBlock))
       ENDDO
       DEALLOCATE(this%vHhhhRef)
    ENDIF
    
    IF (ALLOCATED(this%vHhpp)) THEN
       DO thisBlock=1, this%blocks%nPphhBlocksThis
          CALL RealMatrix_d(this%vHhpp(thisBlock))
       ENDDO
       DEALLOCATE(this%vHhpp)
    ENDIF
    
    IF (ALLOCATED(this%vPppp)) THEN 
       DO thisBlock=1, this%blocks%nPphhBlocksThis
          CALL RealMatrix_d(this%vPppp(thisBlock))
       ENDDO
       DEALLOCATE(this%vPppp)
    ENDIF
   
    IF (ALLOCATED(this%vPphh)) THEN
       DO thisBlock=1, this%blocks%nPphhBlocksThis
          CALL RealMatrix_d(this%vPphh(thisBlock))
       ENDDO
       DEALLOCATE(this%vPphh)
    ENDIF
    
    IF (ALLOCATED(this%vPphhRe)) THEN
       DO blockRe=1, this%blocksRe%nBlocksReThis
          CALL RealMatrix_d(this%vPphhRe(blockRe))
       ENDDO
       DEALLOCATE(this%vPphhRe)
    ENDIF
    
    IF (ALLOCATED(this%vPhphRe)) THEN
       DO blockRe=1, this%blocksRe%nBlocksReThis
          CALL RealMatrix_d(this%vPhphRe(blockRe))
       ENDDO
       DEALLOCATE(this%vPhphRe)
    ENDIF
    
    IF (ASSOCIATED(this%basisSp)) THEN
       CALL SpBasis_d(this%basisSp)
       this%basisSp => NULL()
    ENDIF
    
    IF (ASSOCIATED(this%blocks)) THEN
       CALL BlockList_d(this%blocks)
       this%blocks => NULL()
    ENDIF
    
    IF (ASSOCIATED(this%blocksRe)) THEN
       CALL BlockListRe_d(this%blocksRe)
       this%blocksRe => NULL()
    ENDIF
    
    
  END SUBROUTINE Hamiltonian_d
  
  
END MODULE HamiltonianMod
