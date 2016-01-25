
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
     TYPE(BlockList) :: blocks
     TYPE(BlockListRe) :: blocksRe
     TYPE(RealMatrix), ALLOCATABLE :: vHhhh(:)     
     TYPE(RealMatrix), ALLOCATABLE :: vHhhp(:)
     TYPE(RealMatrix), ALLOCATABLE :: vHhpp(:)
     TYPE(RealMatrix), ALLOCATABLE :: vPphh(:)
     TYPE(RealMatrix), ALLOCATABLE :: vHppp(:)
     TYPE(RealMatrix), ALLOCATABLE :: vPppp(:)
     TYPE(RealMatrix), ALLOCATABLE :: vPhph(:)    
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
    
    ALLOCATE(bra(2), ket(2))
    
    DO blockIndex=1, this%blocksRe%nBlocks
       
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
       
       this%vPphhRe(blockIndex) = vPphhReMatrix
       this%vPhphRe(blockIndex) = vPhphReMatrix
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
    TYPE(RealMatrix) :: vHhhhMatrix, vHhhpMatrix
    TYPE(RealMatrix) :: vHhppMatrix, vHpppMatrix
    TYPE(RealMatrix) :: vPpppMatrix, vPhPhMatrix
    TYPE(RealMatrix) :: vPphhMatrix
    TYPE(IntMatrix) :: hhPairs, ppPairs, phPairs
    INTEGER :: nHoles, nParticles, blockIndex
    INTEGER :: nHhPairs, nPpPairs, nPhPairs
    INTEGER :: pair1, pair2
    REAL(DP) :: interaction
    
    ALLOCATE(bra(2), ket(2))
    
    nHoles = this%basisSp%nOccupied
    nParticles = this%basisSp%nUnoccupied
    
!    WRITE(*,*) ' '
    !     Loop over (K_CM, M_S) blocks
    DO blockIndex=1, this%blocks%nBlocks
       
!       WRITE(*,*) ' Block number', blockIndex
       thisBlock = this%blocks%list(blockIndex)
       !     hh, pp, and, ph pairs
       hhPairs = thisBlock%hhPairs
       ppPairs = thisBlock%ppPairs
       phPairs = thisBlock%phPairs
       nHhPairs = hhPairs%dim1
       nPpPairs = ppPairs%dim1
       nPhPairs = phPairs%dim1
       
       IF (thisBlock%occupiedHh) THEN
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
          this%vHhhh(blockIndex) = vHhhhMatrix
        
       ENDIF
       
       IF (thisBlock%occupiedHh .AND. thisBlock%occupiedPh) THEN   
!          WRITE(*,*) 'hh pairs:', nHhPairs,', ph pairs:', nPhPairs
          !     Hole-hole hole-particle matrix
          vHhhpMatrix = RealMatrix_(nHhPairs, nPhPairs)
          DO pair1=1, nHhPairs
             bra(:) = hhPairs%matr(pair1,:)
             DO pair2=1, nPhPairs
                ket(:) = PhPairs%matr(pair2,:)
                
                !     Get antisymmetrized interaction
                !     matrix element
                interaction = this%getInteraction(&
                     bra(1), bra(2), ket(1), ket(2))
                
                vHhhpMatrix%matr(pair1, pair2) = -interaction
             ENDDO
          ENDDO
          this%vHhhp(blockIndex) = vHhhpMatrix
          
       ENDIF

       IF (thisBlock%occupiedHh .AND. thisBlock%occupiedPp) THEN
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
          this%vHhpp(blockIndex) = vHhppMatrix
          this%vPphh(blockIndex) = vPphhMatrix
         
       ENDIF
     
       IF (thisBlock%occupiedPh .AND. thisBlock%occupiedPp) THEN
!          WRITE(*,*) 'ph pairs:', nPhPairs,', pp pairs:',nPpPairs
          !     Hole-particle particle-particle matrix
          vHpppMatrix = RealMatrix_(nPhPairs, nPpPairs)

          DO pair1=1, nPhPairs
             bra(:) = PhPairs%matr(pair1,:)
             DO pair2=1, nPpPairs
                ket(:) = PpPairs%matr(pair2,:)

                !     Get antisymmetrized interaction
                !     matrix element
                interaction = this%getInteraction(&
                     bra(1), bra(2), ket(1), ket(2))

                vHpppMatrix%matr(pair1, pair2) = -interaction
             ENDDO
          ENDDO
          this%vHppp(blockIndex) = vHpppMatrix

       ENDIF
       
       IF (thisBlock%occupiedPp) THEN
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
          this%vPppp(blockIndex) = vPpppMatrix
       ENDIF
       
       IF (thisBlock%occupiedPh) THEN
!          WRITE(*,*) 'ph pairs:', nPhPairs 
          !     Particle-hole paritcle-hole matrix
          vPhphMatrix = RealMatrix_(nPhPairs, nPhPairs)
          
          DO pair1=1, nPhPairs
             bra(:) = PhPairs%matr(pair1,:)
             DO pair2=1, nPhPairs
                ket(:) = PhPairs%matr(pair2,:)
                
                !     Get antisymmetrized interaction
                !     matrix element
                interaction = this%getInteraction(&
                     bra(1), bra(2), ket(1), ket(2))
                
                vPhphMatrix%matr(pair1, pair2) = interaction
             ENDDO
          ENDDO
          this%vPhph(blockIndex) = vPhphMatrix

       ENDIF
!       WRITE(*,*) ' '
    ENDDO
!    WRITE(*,*) ' '
    DEALLOCATE(bra, ket)
       
    
  END SUBROUTINE setupInteractionMatr
  
  
  ! !
  ! !     Routine to set up an interaction matrix
  ! !     in the given single-particle basis.
  ! !
  ! SUBROUTINE setupInteractionMatr(this)
  !   CLASS(Hamiltonian), TARGET, INTENT(INOUT) :: this
  !   INTEGER :: nOrbits, orbit1, orbit2, orbit3, orbit4  
  !   REAL(DP) :: interaction
    
  !   nOrbits = this%basisSp%nOccupied &
  !        + this%basisSp%nUnoccupied
    
  !   DO orbit1=1, nOrbits
  !      DO orbit2=orbit1, nOrbits
  !         DO orbit3=1, nOrbits
  !            DO orbit4=orbit3, nOrbits
                
  !               !     Get antisymmetrized interaction
  !               !     matrix element
  !               interaction = this%getInteraction(&
  !                    orbit1, orbit2, orbit3, orbit4)
                
  !               !     Store the antisymmetrized 
  !               !     interaction matrix element
  !               this%vMatrix(orbit1, orbit2, &
  !                    orbit3, orbit4) = interaction
  !               this%vMatrix(orbit2, orbit1, &
  !                    orbit3, orbit4) = -interaction
  !               this%vMatrix(orbit1, orbit2, &
  !                    orbit4, orbit3) = -interaction
  !               this%vMatrix(orbit2, orbit1, &
  !                    orbit4, orbit3) = interaction
                
  !            ENDDO
  !         ENDDO
  !      ENDDO
  !   ENDDO
    
    
  ! END SUBROUTINE setupInteractionMatr
  
  
END MODULE HamiltonianMod
