
!
!     Class for pairs of single-particle states.
!
MODULE SpPairsMod
  USE IntMatListMod
  USE ConstantsMod
  USE SpBasisMod
  IMPLICIT NONE
  
  TYPE, ABSTRACT :: SpPairs
     !
     !     nBlocks:           Number of blocks (K_CM, M_S)
     !     blockToQuant:      Gives quantum numbers (K_CM, M_S) 
     !                        for a chosen block index
     !     blockToQuantRe:    Gives relative (recoupled) quantum 
     !                        numbers (k_rel, m_s_rel) for a 
     !                        chosen block index
     !     hhPairs:           Pairs of two hole states
     !     ppPairs:           Pairs of two particle states
     !     phPairs:           Pairs of a particle and a 
     !                        hole state
     !
     INTEGER :: nBlocks
     INTEGER, ALLOCATABLE :: blockToQuantHh(:,:)  
     INTEGER, ALLOCATABLE :: blockToQuantPp(:,:)
     INTEGER, ALLOCATABLE :: blockToQuantPh(:,:)
     INTEGER, ALLOCATABLE :: blockToQuantHhRe(:,:)
     INTEGER, ALLOCATABLE :: blockToQuantPpRe(:,:)
     INTEGER, ALLOCATABLE :: blockToQuantPhRe(:,:)
     TYPE(IntMatList) :: hhPairs 
     TYPE(IntMatList) :: ppPairs 
     TYPE(IntMatList) :: phPairs 
     CLASS(SpBasis), POINTER :: basisSp
   CONTAINS
     PROCEDURE(setupPairs), &
          DEFERRED :: setupSpPairs
  END TYPE SpPairs
  
  ABSTRACT INTERFACE
     SUBROUTINE setupPairs(this, basis) 
       IMPORT SpBasis, SpPairs
       CLASS(SpPairs), TARGET, INTENT(INOUT) :: this
       CLASS(SpBasis), TARGET, INTENT(IN) :: basis
     END SUBROUTINE setupPairs
  END INTERFACE
  
  
END MODULE SpPairsMod
