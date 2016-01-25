
!
!     Class for the Coupled-Cluster doubles (CCD) t-amplitudes.
!
MODULE CcdAmplitudeMod
  USE MemoryCounterMod
  USE HamiltonianMod
  USE FockMatrixMod
  USE ConstantsMod
  IMPLICIT NONE
  
  TYPE CcdAmplitude
     !
     !     basisSp:      Single-particle basis
     !
     CLASS(SpBasis), POINTER :: basisSp
     TYPE(BlockList), POINTER :: blocks 
     TYPE(BlockListRe), POINTER :: blocksRe
     TYPE(RealMatrix), ALLOCATABLE :: tPphh(:)
     TYPE(RealMatrix), ALLOCATABLE :: tPphhOld(:)
     TYPE(RealMatrix), ALLOCATABLE :: tPphhRe(:)
     TYPE(RealMatrix), ALLOCATABLE :: edMatrices(:)
     TYPE(RealMatrix), ALLOCATABLE :: intermediates1(:)
     TYPE(RealMatrix), ALLOCATABLE :: intermediates2(:)
     REAL(DP), ALLOCATABLE :: intermediate3(:)
     REAL(DP), ALLOCATABLE :: intermediate4(:)
   CONTAINS
     PROCEDURE :: setupEnergyDenomMatrices
     PROCEDURE :: updateIntermed1
     PROCEDURE :: updateIntermed2
     PROCEDURE :: updateIntermed3
     PROCEDURE :: updateIntermed4
     PROCEDURE :: updateT2Matrix
     PROCEDURE :: calcT2Matrix
     PROCEDURE :: setupT2MatrixRe
  END TYPE CcdAmplitude
  
CONTAINS
  
  !
  !     Constructor for CcdAmplitude object
  !
  FUNCTION CcdAmplitude_(basis, blocks, blocksRe) &
       RESULT(this)
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    TYPE(BlockList), TARGET, INTENT(IN) :: blocks
    TYPE(BlockListRe), TARGET, INTENT(IN) :: blocksRe
    TYPE(CcdAmplitude), POINTER :: this
    TYPE(Block1) :: thisBlock
    INTEGER :: nPphhBlocksThis, nPphhBlocksAll 
    INTEGER :: nHh, nPp, nPh, nHp, blockIndex
    INTEGER :: elementsR8, pphhBlock, blockRe
    
    ALLOCATE(this)
    this%basisSp => basis
    this%blocks => blocks
    this%blocksRe => blocksRe
    nPphhBlocksThis = blocks%nPphhBlocksThis
    nPphhBlocksAll = blocks%nPphhBlocksAll
    
    ALLOCATE(this%tPphh(nPphhBlocksThis))
    ALLOCATE(this%tPphhOld(nPphhBlocksThis))
    ALLOCATE(this%tPphhRe(blocksRe%nBlocks))
    ALLOCATE(this%edMatrices(nPphhBlocksThis))
    ALLOCATE(this%intermediates1(nPphhBlocksThis))
    ALLOCATE(this%intermediates2(nPphhBlocksAll))
    ALLOCATE(this%intermediate3(basis%nUnoccupied))
    ALLOCATE(this%intermediate4(basis%nOccupied))
    this%intermediate3 = 0.d0
    this%intermediate4 = 0.d0
    
    elementsR8 = basis%nUnoccupied + basis%nOccupied 
    CALL add2MemoryR8(elementsR8)
    
    DO pphhBlock=1, nPphhBlocksThis
       blockIndex = blocks%pphhThis2BlocksAll(pphhBlock)
       nPp = blocks%list(blockIndex)%ppPairs%dim1
       nHh = blocks%list(blockIndex)%hhPairs%dim1
       
       this%tPphh(pphhBlock) = RealMatrix_(nPp, nHh)
       this%tPphhOld(pphhBlock) = RealMatrix_(nPp, nHh)
       
       this%edMatrices(pphhBlock) = RealMatrix_(nPp, nHh)
       this%intermediates1(pphhBlock) = RealMatrix_(nHh, nHh)
    ENDDO
    
    DO pphhBlock=1, nPphhBlocksAll
       blockIndex = blocks%pphh2BlocksAll(pphhBlock)
       nPp = blocks%list(blockIndex)%ppPairs%dim1
       nHh = blocks%list(blockIndex)%hhPairs%dim1
       
       this%intermediates2(pphhBlock) = RealMatrix_(nPp, nHh)
    ENDDO
    
    DO blockRe=1, blocksRe%nBlocks
       nPh = blocksRe%list(blockRe)%phPairs%dim1
       nHp = blocksRe%list(blockRe)%hpPairs%dim1
       
       this%tPphhRe(blockRe) = RealMatrix_(nPh, nHp)       
    ENDDO
    
    
  END FUNCTION CcdAmplitude_
  
  !
  !     Destrucor for CcdAmplitude object
  !
  SUBROUTINE CcdAmplitude_d(this)
    TYPE(CcdAmplitude), INTENT(INOUT) :: this
    INTEGER :: elementsR8, pphhBlock, blockRe
    
    elementsR8 = 0
    IF (ALLOCATED(this%intermediate3)) THEN
       DEALLOCATE(this%intermediate3)
       elementsR8 = elementsR8 - this%basisSp%nUnoccupied
    ENDIF
       
    IF (ALLOCATED(this%intermediate4)) THEN
       DEALLOCATE(this%intermediate4)
       elementsR8 = elementsR8 - this%basisSp%nOccupied
    ENDIF
    CALL add2MemoryR8(elementsR8)
    
    DO pphhBlock=1, this%blocks%nPphhBlocksThis
       
       CALL RealMatrix_d(this%tPphh(pphhBlock))
       CALL RealMatrix_d(this%tPphhOld(pphhBlock))
       CALL RealMatrix_d(this%edMatrices(pphhBlock))
       CALL RealMatrix_d(this%intermediates1(pphhBlock))
    ENDDO
    
    DO pphhBlock=1, this%blocks%nPphhBlocksAll
       
       CALL RealMatrix_d(this%intermediates2(pphhBlock))
    ENDDO
    
    DEALLOCATE(this%tPphh)
    DEALLOCATE(this%tPphhOld)
    DEALLOCATE(this%edMatrices)
    DEALLOCATE(this%intermediates1)
    DEALLOCATE(this%intermediates2)
 
    DO blockRe=1, this%blocksRe%nBlocks
       
       CALL RealMatrix_d(this%tPphhRe(blockRe))
    ENDDO
    
    DEALLOCATE(this%tPphhRe)
    
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
    
    
  END SUBROUTINE CcdAmplitude_d
  
  !
  !     Set up list of energy denominator matrices. 
  !
  SUBROUTINE setupEnergyDenomMatrices(this, fockMatr)
    CLASS(CcdAmplitude), TARGET, INTENT(INOUT) :: this
    TYPE(FockMatrix), TARGET, INTENT(IN) :: fockMatr
    TYPE(IntMatrix) :: ppPairs, hhPairs
    TYPE(RealMatrix) :: edMatrix
    INTEGER :: particle1, particle2, hole1, hole2
    REAL(DP) :: eneP1, eneP2, eneH1, eneH2
    INTEGER :: pphhBlock, blockIndex, pp, hh
    REAL(DP) :: denom
    
    !     Loop over (K_CM, M_S) blocks.
    DO pphhBlock=1, this%blocks%nPphhBlocksThis
       blockIndex = this%blocks%pphhThis2BlocksAll(pphhBlock)
       
       ppPairs = this%blocks%list(blockIndex)%ppPairs
       hhPairs = this%blocks%list(blockIndex)%hhPairs
       
       edMatrix = RealMatrix_(ppPairs%dim1, hhPairs%dim1)
       
       DO pp=1, ppPairs%dim1
          particle1 = ppPairs%matr(pp, 1)
          particle2 = ppPairs%matr(pp, 2)
          
          !     Single-particle energies of particles 
          !     1 and 2
          eneP1 = fockMatr%vector(particle1)
          eneP2 = fockMatr%vector(particle2)
          
          DO hh=1, hhPairs%dim1
             hole1 = hhPairs%matr(hh, 1)
             hole2 = hhPairs%matr(hh, 2)
             
             !     Single-particle energies of holes
             !     1 and 2
             eneH1 = fockMatr%vector(hole1)
             eneH2 = fockMatr%vector(hole2)
             
             !     Energy denominator
             denom = eneH1 + eneH2 - eneP1 - eneP2
             
             edMatrix%matr(pp, hh) = 1.d0/denom
          ENDDO
       ENDDO
       this%edMatrices(pphhBlock) = edMatrix 
       
       CALL RealMatrix_d(edMatrix)
    ENDDO
    
    
  END SUBROUTINE setupEnergyDenomMatrices
  
  !
  !     Intermediate matrix to diagram 3 of
  !     the T2 amplitude equation.
  !
  !     I_ij^kl = <kl||ij> + (1/2)*\sum_cd <kl||cd><cd|t|ij>
  !
  SUBROUTINE updateIntermed1(this, ham)
    CLASS(CcdAmplitude), TARGET, INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    TYPE(RealMatrix) :: vHhhh
    REAL(DP), ALLOCATABLE :: vHhpp(:,:)
    REAL(DP), ALLOCATABLE :: tPphh(:,:)
    INTEGER :: pphhBlock, blockIndex, nPp, nHh
    
    !     Loop over (K_CM, M_S) blocks.
    DO pphhBlock=1, this%blocks%nPphhBlocksThis
       
       nPp = this%tPphh(pphhBlock)%dim1
       nHh = this%tPphh(pphhBlock)%dim2 
       
       ALLOCATE(vHhpp(nHh, nPp))
       ALLOCATE(tPphh(nPp, nHh))
       
       vHhhh = ham%vHhhh(pphhBlock)
       vHhpp = ham%vHhpp(pphhBlock)%matr
       tPphh = this%tPphhOld(pphhBlock)%matr
       CALL DGEMM('N', 'N', nHh, nHh, nPp, 0.5d0, &
            vHhpp, nHh, tPphh, nPp, 1.d0, &
            vHhhh%matr, nHh)
       
       this%intermediates1(pphhBlock) = vHhhh
       
       DEALLOCATE(vHhpp, tPphh)
       !CALL RealMatrix_d(vHhhh)
    ENDDO
    
    
  END SUBROUTINE updateIntermed1
  
  
  !
  !     Intermediate matrix to diagram 4 of
  !     the T2 amplitude equation.
  !
  !     I_cj^kb = <kb||cj> + (1/2)*\sum_ld <kl||cd><db|t|lj>
  !
  SUBROUTINE updateIntermed2(this, ham)
    CLASS(CcdAmplitude), TARGET, INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    TYPE(Block1) :: thisPhBlock, thisPphhBlock
    TYPE(IntMatrix) :: phPairs, thisMatrix
    TYPE(BlockListRe) :: blocksRe
    REAL(DP), ALLOCATABLE :: vPhphRe(:,:)
    REAL(DP), ALLOCATABLE :: vPphhRe(:,:)
    REAL(DP), ALLOCATABLE :: tPphhRe(:,:)
    REAL(DP), ALLOCATABLE :: diagram4ReIabj(:,:)
    TYPE(RealMatrix), ALLOCATABLE :: intermediates2(:)
    INTEGER :: blockIndex, nPhPairs, nHpPairs, hp, ph
    INTEGER :: pphhIndex, pphhBlockIndex, iError
    INTEGER :: pp, hh, blockRe, dim1, dim2
    INCLUDE 'mpif.h'
    
    !     Wait until all processes are ready
    !CALL MPI_BARRIER(MPI_COMM_WORLD, iError)
    
    blocksRe = this%blocksRe
    ALLOCATE(intermediates2(this%blocks%nPphhBlocksAll))
    DO pphhBlockIndex=1, this%blocks%nPphhBlocksAll
       this%intermediates2(pphhBlockIndex)%matr = 0.d0
    ENDDO
    intermediates2 = this%intermediates2
    
    DO blockRe=1, this%blocksRe%nBlocksReThis
       blockIndex = this%blocksRe%blocksThis2BlocksAll(blockRe)
       
       nPhPairs = ham%vPhphRe(blockRe)%dim1
       nHpPairs = ham%vPphhRe(blockRe)%dim1
       
       ALLOCATE(vPhphRe(nPhPairs, nPhPairs))
       ALLOCATE(vPphhRe(nHpPairs, nPhPairs))
       ALLOCATE(tPphhRe(nPhPairs, nHpPairs))
       ALLOCATE(diagram4ReIabj(nPhPairs, nHpPairs))
       
       vPhphRe = ham%vPhphRe(blockRe)%matr
       vPphhRe = ham%vPphhRe(blockRe)%matr
       tPphhRe = this%tPphhRe(blockIndex)%matr
       
       !     Matrix elements cross-coupled from
       !     pphh to hpph are denoted by a star (*),
       !     and from phph to phph by a tilde (~).
       !
       !     <bj|i*|ck> = \sum_ld <bj|t*|ld><ld|v*|ck> 
       !                  + <bj|v~|ck>
       CALL DGEMM('N', 'N', nPhPairs, nPhPairs, nHpPairs, &
            0.5d0, tPphhRe, nPhPairs, vPphhRe, nHpPairs, &
            -1.d0, vPhphRe, nPhPairs)
       ! CALL DGEMM('N', 'N', nPhPairs, nPhPairs, nHpPairs, &
       !      0.5d0, tPphhRe, nPhPairs, vPphhRe, nHpPairs, &
       !      0.d0, vPhphRe, nPhPairs)
       
       !
       !     <bj|t4|ia> = \sum_kc <bj|i*|ck><ck|t*|ia>
       !
       CALL DGEMM('N', 'N', nPhPairs, nPhPairs, nPhPairs, &
            1.d0, vPhphRe, nPhPairs, tPphhRe, nPhPairs, &
            0.d0, diagram4ReIabj, nPhPairs) 
       
       ! CALL DGEMM('N', 'N', nPhPairs, nPhPairs, nPhPairs, &
       !      -1.d0, vPhphRe, nPhPairs, tPphhRe, nPhPairs, &
       !      0.d0, diagram4ReIabj, nPhPairs) 
       
       !     Write back to pphh matrix elements
       thisPhBlock = this%blocksRe%list(blockIndex)
       phPairs = thisPhBlock%phPairs
       DO hp=1, nHpPairs
          DO ph=1, nPhPairs
             
             !
             !     ijab
             !
             thisMatrix = blocksRe%iabjRe2Abij%list(blockRe)
             pphhIndex = thisMatrix%matr(hp, ph)
             
             pphhBlockIndex &
                  = blocksRe%index2Abij(1, pphhIndex)
             pp = blocksRe%index2Abij(2, pphhIndex)
             hh = blocksRe%index2Abij(3, pphhIndex)
             
             intermediates2(pphhBlockIndex)%matr(pp, hh) &
                  = intermediates2(pphhBlockIndex)%matr(pp, hh) &
                  + diagram4ReIabj(ph, hp)             
             
             !
             !     ijba
             !
             thisMatrix = blocksRe%ibajRe2Abij%list(blockRe)
             pphhIndex = thisMatrix%matr(hp, ph)
             
             pphhBlockIndex &
                  = blocksRe%index2Abij(1, pphhIndex)
             pp = blocksRe%index2Abij(2, pphhIndex)
             hh = blocksRe%index2Abij(3, pphhIndex)
             
             intermediates2(pphhBlockIndex)%matr(pp, hh) &
                  = intermediates2(pphhBlockIndex)%matr(pp, hh) &
                  - diagram4ReIabj(ph, hp)
             
             !
             !     jiab
             !
             thisMatrix = blocksRe%jabiRe2Abij%list(blockRe)
             pphhIndex = thisMatrix%matr(hp, ph)
             
             pphhBlockIndex &
                  = blocksRe%index2Abij(1, pphhIndex)
              pp = blocksRe%index2Abij(2, pphhIndex)
             hh = blocksRe%index2Abij(3, pphhIndex)

             intermediates2(pphhBlockIndex)%matr(pp, hh) &
                  = intermediates2(pphhBlockIndex)%matr(pp, hh) &
                  - diagram4ReIabj(ph, hp)
             
             !
             !     jiba
             !
             thisMatrix = blocksRe%jbaiRe2Abij%list(blockRe)
             pphhIndex = thisMatrix%matr(hp, ph)
             
             pphhBlockIndex &
                  = blocksRe%index2Abij(1, pphhIndex)
             pp = blocksRe%index2Abij(2, pphhIndex)
             hh = blocksRe%index2Abij(3, pphhIndex)
             
             intermediates2(pphhBlockIndex)%matr(pp, hh) &
                  = intermediates2(pphhBlockIndex)%matr(pp, hh) &
                  + diagram4ReIabj(ph, hp)
          ENDDO
       ENDDO
       
       DEALLOCATE(vPhphRe)
       DEALLOCATE(vPphhRe)
       DEALLOCATE(tPphhRe)
       DEALLOCATE(diagram4ReIabj)
    ENDDO
    
    DO pphhBlockIndex=1, this%blocks%nPphhBlocksAll
       dim1 = this%intermediates2(pphhBlockIndex)%dim1
       dim2 = this%intermediates2(pphhBlockIndex)%dim2
       
       !     Sum the contribution from all processes
       !     and broadcast the result to all processes
       CALL MPI_ALLREDUCE(intermediates2(pphhBlockIndex)%matr, &
            this%intermediates2(pphhBlockIndex)%matr, &
            dim1*dim2, MPI_DOUBLE_PRECISION, MPI_SUM, &
            MPI_COMM_WORLD, iError)
    ENDDO
    
    DEALLOCATE(intermediates2)
    
    
  END SUBROUTINE updateIntermed2
  
  !
  !     Intermediate matrix to diagram 5 of
  !     the T2 amplitude equation.
  !
  !     I_c^b = -(1/2)*\sum_kld <kl||cd><bd|t|kl> 
  !
  SUBROUTINE updateIntermed3(this, ham)
    CLASS(CcdAmplitude), TARGET, INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    TYPE(Block1) :: thisBlock
    REAL(DP), ALLOCATABLE :: vHhpp(:,:)
    REAL(DP), ALLOCATABLE :: tPphh(:,:)
    REAL(DP), ALLOCATABLE :: tv(:,:)
    REAL(DP), ALLOCATABLE :: intermediate3(:)
    INTEGER :: pphhBlock, blockIndex, nPp, nHh, pp
    INTEGER :: particle1, nHoles, iError
    INCLUDE 'mpif.h'
    
    ALLOCATE(intermediate3(this%basisSp%nUnoccupied))
    intermediate3 = 0.d0
    nHoles = this%basisSp%nOccupied
    
    DO pphhBlock=1, this%blocks%nPphhBlocksThis
       blockIndex = this%blocks%pphhThis2BlocksAll(pphhBlock)
       thisBlock = this%blocks%list(blockIndex)
       
       nPp = this%tPphh(pphhBlock)%dim1
       nHh = this%tPphh(pphhBlock)%dim2
       
       ALLOCATE(vHhpp(nHh, nPp))
       ALLOCATE(tPphh(nPp, nHh))
       ALLOCATE(tv(nPp, nPp))
       
       vHhpp = ham%vHhpp(pphhBlock)%matr
       tPphh = this%tPphhOld(pphhBlock)%matr
       CALL DGEMM('N', 'N', nPp, nPp, nHh, -0.5d0, &
            tPphh, nPp, vHhpp, nHh, 0.d0, &
            tv, nPp)
       
       DO pp=1, nPP
          particle1 = thisBlock%ppPairs%matr(pp, 1) &
               - nHoles
          
          intermediate3(particle1) &
               = intermediate3(particle1) + tv(pp, pp)
       ENDDO
       
       DEALLOCATE(vHhpp, tPphh, tv)
    ENDDO
    
    !     Sum the contribution from all processes
    !     and broadcast the result to all processes
    CALL MPI_ALLREDUCE(intermediate3, this%intermediate3, &
         this%basisSp%nUnoccupied, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_COMM_WORLD, iError)
    
    DEALLOCATE(intermediate3)
    
    
  END SUBROUTINE updateIntermed3
  
  !
  !     Intermediate matrix to diagram 6 of
  !     the T2 amplitude equation.
  !
  !     I_j^k = -(1/2)*\sum_lcd <kl||cd><cd|t|jl>
  !
  SUBROUTINE updateIntermed4(this, ham)
    CLASS(CcdAmplitude), TARGET, INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    TYPE(Block1) :: thisBlock
    REAL(DP), ALLOCATABLE :: vHhpp(:,:)
    REAL(DP), ALLOCATABLE :: tPphh(:,:)
    REAL(DP), ALLOCATABLE :: vt(:,:)
    REAL(DP), ALLOCATABLE :: intermediate4(:)
    INTEGER :: pphhBlock, blockIndex, nPp, nHh, hh
    INTEGER :: hole1, iError
    INCLUDE 'mpif.h'
    
    ALLOCATE(intermediate4(this%basisSp%nOccupied))
    intermediate4 = 0.d0
    
    DO pphhBlock=1, this%blocks%nPphhBlocksThis
       blockIndex = this%blocks%pphhThis2BlocksAll(pphhBlock)
       thisBlock = this%blocks%list(blockIndex)
       
       nPp = this%tPphh(pphhBlock)%dim1
       nHh = this%tPphh(pphhBlock)%dim2
       
       ALLOCATE(vHhpp(nHh, nPp))
       ALLOCATE(tPphh(nPp, nHh))
       ALLOCATE(vt(nHh, nHh))
       
       vHhpp = ham%vHhpp(pphhBlock)%matr
       tPphh = this%tPphhOld(pphhBlock)%matr
       CALL DGEMM('N', 'N', nHh, nHh, nPp, -0.5d0, &
            vHhpp, nHh, tPphh, nPp, 0.d0, &
            vt, nHh)
       
       DO hh=1, nHh
          hole1 = thisBlock%hhPairs%matr(hh, 1) 
          
          intermediate4(hole1) = intermediate4(hole1) &
               + vt(hh, hh)
       ENDDO
       
       DEALLOCATE(vHhpp, tPphh, vt)
    ENDDO

    !     Sum the contribution from all processes
    !     and broadcast the result to all processes
    CALL MPI_ALLREDUCE(intermediate4, this%intermediate4, &
         this%basisSp%nOccupied, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_COMM_WORLD, iError)
    
    DEALLOCATE(intermediate4)
    
    
  END SUBROUTINE updateIntermed4
  
  !
  !     Set up the cross-coupled t-amplitude matrix
  !
  SUBROUTINE setupT2MatrixRe(this)
    CLASS(CcdAmplitude), TARGET, INTENT(INOUT) :: this
    TYPE(BlockListRe), POINTER :: blocksRe
    TYPE(IntMatrix) :: thisMatrix
    TYPE(RealMatrix), ALLOCATABLE :: tPphhRe(:)
    INTEGER :: hpphIndex, hpphBlockIndex
    INTEGER :: pphhBlock, nPp, nHh, pp, hh, hp, ph
    INTEGER :: blockRe, dim1, dim2, iError
    INCLUDE 'mpif.h'
    
    blocksRe => this%blocksRe
    ALLOCATE(tPphhRe(blocksRe%nBlocks))
    DO blockRe=1, blocksRe%nBlocks
       this%tPphhRe(blockRe)%matr = 0.d0
    ENDDO
    tPphhRe = this%tPphhRe
    
    DO pphhBlock=1, this%blocks%nPphhBlocksThis
       
       thisMatrix = blocksRe%abij2ReIabj%list(pphhBlock)
       
       nPp = this%tPphh(pphhBlock)%dim1
       nHh = this%tPphh(pphhBlock)%dim2 
       
       IF (thisMatrix%dim1 == 0) CYCLE
       
       DO pp=1, nPp
          DO hh=1, nHh
             
             hpphIndex = thisMatrix%matr(pp, hh)
             hpphBlockIndex &
                  = blocksRe%index2ReIabj(1, hpphIndex)
             hp = blocksRe%index2ReIabj(2, hpphIndex)
             ph = blocksRe%index2ReIabj(3, hpphIndex)
             
             IF (hpphIndex == 0) CYCLE
             IF (hpphBlockIndex == 0) CYCLE
             IF (hpphBlockIndex > blocksRe%nBlocks) CYCLE
             
             tPphhRe(hpphBlockIndex)%matr(ph, hp) &
                  = this%tPphh(pphhBlock)%matr(pp, hh)
          ENDDO
       ENDDO
    ENDDO
    
    DO hpphBlockIndex=1, blocksRe%nBlocks
       dim1 = this%tPphhRe(hpphBlockIndex)%dim1
       dim2 = this%tPphhRe(hpphBlockIndex)%dim2
       
       !     Sum the contribution from all processes
       !     and broadcast the result to all processes
       CALL MPI_ALLREDUCE(tPphhRe(hpphBlockIndex)%matr, &
            this%tPphhRe(hpphBlockIndex)%matr, &
            dim1*dim2, MPI_DOUBLE_PRECISION, MPI_SUM, &
            MPI_COMM_WORLD, iError)
  
    ENDDO
        
    DEALLOCATE(tPphhRe)
    
    
  END SUBROUTINE setupT2MatrixRe
  
  
  !
  !     Update the CCSD t2-amplitude matrix.
  !
  SUBROUTINE updateT2Matrix(this, ham, alpha, iteration)
    CLASS(CcdAmplitude), TARGET, INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    REAL(DP), INTENT(IN) :: alpha 
    INTEGER, INTENT(IN) :: iteration
    INTEGER :: pphhBlock, blockIndex
    
    !     Loop over (K_CM, M_S) blocks.
    DO pphhBlock=1, this%blocks%nPphhBlocksThis 
       blockIndex = this%blocks%pphhThis2BlocksAll(pphhBlock)
       
       !     Get a new CCD t2-amplitude
       this%tPphh(pphhBlock) = this%calcT2Matrix(&
            ham, pphhBlock, blockIndex)
       
       IF (iteration > 1) THEN
          !     t = alpha*t + (1-alpha)*t_old
          this%tPphh(pphhBlock)%matr = &
               alpha*this%tPphh(pphhBlock)%matr &
               + (1.d0-alpha)*this%tPphhOld(pphhBlock)%matr
       ENDIF
    ENDDO
    this%tPphhOld = this%tPphh
    
    
  END SUBROUTINE updateT2Matrix
  
  !
  !     Calculate a CCD t2-amplitude matrix for a
  !     given (K_CM, M_S) channel (block).
  !
  FUNCTION calcT2Matrix(this, ham, pphhBlock, blockIndex) &
       RESULT(t2Matrix)
    CLASS(CcdAmplitude), TARGET, INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    INTEGER, INTENT(IN) :: pphhBlock, blockIndex
    TYPE(Block1) :: thisBlock
    REAL(DP), ALLOCATABLE :: tempMatrix(:,:)
    REAL(DP), ALLOCATABLE :: tempPphh(:,:), tempPppp(:,:)
    REAL(DP), ALLOCATABLE :: tempHhhh(:,:)
    REAL(DP), ALLOCATABLE :: tPphhMatrix(:,:)
    INTEGER :: nPp, nHh, pp, hh, nHoles, pphhBlockAll
    INTEGER :: particle1, particle2, hole1, hole2
    TYPE(RealMatrix) :: t2Matrix
    REAL(DP) :: intermediate
    
    nPp = this%tPphh(pphhBlock)%dim1
    nHh = this%tPphh(pphhBlock)%dim2  
    ALLOCATE(tempMatrix(nPp, nHh))
    ALLOCATE(tempPphh(nPp, nHh))
    ALLOCATE(tempPppp(nPp, nPP))
    ALLOCATE(tempHhhh(nHh, nHh))
    ALLOCATE(tPphhMatrix(nPp, nHh))
    
    !     Diagram 1: <ab||ij>
    tempMatrix = ham%vPphh(pphhBlock)%matr
    
    !     Diagram 2: (1/2)*\sum_cd <ab||cd><cd|t|ij>
    tempPppp = ham%vPppp(pphhBlock)%matr
    tempPphh = this%tPphhOld(pphhBlock)%matr
    CALL DGEMM('N', 'N', nPp, nHh, nPp, 0.5d0, &
         tempPppp, nPp, tempPphh, nPp, 1.d0, &
         tempMatrix, nPp)
    
    !     Diagram 3: (1/2)*\sum_kl [<kl||ij>
    !                + (1/2)*\sum_cd <kl||cd><cd|t|ij>]
    !                *<ab|t|kl>
    tempHhhh = this%intermediates1(pphhBlock)%matr
    CALL DGEMM('N', 'N', nPp, nHh, nHh, 0.5d0, &
         tempPphh, nPp, tempHhhh, nHh, 1.d0, &
         tempMatrix, nPp)
    
    ! !     Diagram 4: P(ij)P(ab)\sum_kc I2_cj^kb * <ac|t|ik>
    ! pphhBlockAll = this%blocks%pphhThis2PphhAll(pphhBlock)
    ! tempMatrix = tempMatrix &
    !      + this%intermediates2(pphhBlockAll)%matr
    
    ! ! if (this%basisSp%mpiJobs%iAm == 1) then
    ! !    write(*,*) '==================='
    ! !    write(*,*) pphhBlockAll
    ! !    write(*,*) this%intermediates2(pphhBlockAll)%matr
    ! ! endif
    
    !     Diagram 5: P(ab)\sum_c I_c^b * <ac|t|ij>
    thisBlock = this%blocks%list(blockIndex)
    nHoles = this%basisSp%nOccupied
    DO pp=1, nPp
       particle1 = thisBlock%ppPairs%matr(pp, 1)
       particle2 = thisBlock%ppPairs%matr(pp, 2)
       
       intermediate = this%intermediate3(particle2-nHoles) &
            + this%intermediate3(particle1-nHoles)
       tempMatrix(pp,:) = tempMatrix(pp,:) &
            + tempPphh(pp,:)*intermediate
    ENDDO
    
    !     Diagram 6: P(ij)\sum_k I_j^k * <ab|t|ik>
    DO hh=1, nHh
       hole1 = thisBlock%hhPairs%matr(hh, 1)
       hole2 = thisBlock%hhPairs%matr(hh, 2)
       
       intermediate = this%intermediate4(hole2) &
            + this%intermediate4(hole1)
       tempMatrix(:,hh) = tempMatrix(:,hh) &
            + tempPphh(:,hh)*intermediate
    ENDDO
    
    !     Divide all elements by the energy denominator
    CALL elementMatmul(tempMatrix, &
         this%edMatrices(pphhBlock)%matr, &
         nPp, nHh, tPphhMatrix)
    
    t2Matrix = RealMatrix_(nPp, nHh)
    t2Matrix%matr = tPphhMatrix
    
    CALL add2MemoryR8(-Npp*nHh)
    DEALLOCATE(tempMatrix)
    DEALLOCATE(tPphhMatrix)
    DEALLOCATE(tempPphh, tempPppp)
    DEALLOCATE(tempHhhh)
    
    
  END FUNCTION calcT2Matrix
  
  
END MODULE CcdAmplitudeMod
