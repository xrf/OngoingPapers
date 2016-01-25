
!
!     Class for the Coupled-Cluster singles-doubles (CCSD) t-amplitudes.
!
MODULE CcdAmplitudeMod
  USE HamiltonianMod
  USE FockMatrixMod
  USE ConstantsMod
  IMPLICIT NONE
  
  TYPE CcdAmplitude
     !SEQUENCE
     CLASS(SpBasis), POINTER :: basisSp             ! Single-particle basis
     REAL(DP), ALLOCATABLE :: t1Matrix(:,:)         ! t1-amplitude
     REAL(DP), ALLOCATABLE :: t1MatrixOld(:,:)      ! Old t1-amplitude
     REAL(DP), ALLOCATABLE :: t2Matrix(:,:,:,:)     ! t2-amplitude
     REAL(DP), ALLOCATABLE :: t2MatrixOld(:,:,:,:)  ! Old t2-amplitude
     REAL(DP), ALLOCATABLE :: edMatrix(:,:,:,:)     ! energy denominator
     REAL(DP), ALLOCATABLE :: intermediate1(:,:,:,:)! intermediate diagram 
     REAL(DP), ALLOCATABLE :: intermediate2(:,:,:,:)
     REAL(DP), ALLOCATABLE :: intermediate3(:,:)
     REAL(DP), ALLOCATABLE :: intermediate4(:,:)
   CONTAINS
     PROCEDURE :: setupEnergyDenomMatrix
     PROCEDURE :: updateIntermed1
     PROCEDURE :: updateIntermed2
     PROCEDURE :: updateIntermed3
     PROCEDURE :: updateIntermed4
     PROCEDURE :: updateT2Matrix
     PROCEDURE :: updateT1Matrix
     PROCEDURE :: calcT2MatrixElement
     PROCEDURE :: calcT1MatrixElement
  END TYPE CcdAmplitude
  
  
CONTAINS
  
  !
  !     Constructor for CcdAmplitude object
  !
  FUNCTION CcdAmplitude_(basis) RESULT(tAmpl)
    CLASS(SpBasis), TARGET, INTENT(IN) :: basis
    TYPE(CcdAmplitude) :: tAmpl
    INTEGER :: nOcc, nUnocc
    
    nOcc = basis%nOccupied
    nUnocc = basis%nUnoccupied
    tAmpl%basisSp => basis
    ALLOCATE(tAmpl%t1Matrix(nUnocc, nOcc))
    ALLOCATE(tAmpl%t1MatrixOld(nUnocc, nOcc))
    ALLOCATE(tAmpl%t2Matrix(nUnocc, nUnocc, nOcc, nOcc))
    ALLOCATE(tAmpl%t2MatrixOld(nUnocc, nUnocc, nOcc, nOcc))
    ALLOCATE(tAmpl%edMatrix(nUnocc, nUnocc, nOcc, nOcc))
    ALLOCATE(tAmpl%intermediate1(nOcc, nOcc, nOcc, nOcc))
    ALLOCATE(tAmpl%intermediate2(nOcc, nUnocc, nUnocc, nOcc))
    ALLOCATE(tAmpl%intermediate3(nUnocc, nUnocc))
    ALLOCATE(tAmpl%intermediate4(nOcc, nOcc))
    tAmpl%t1Matrix = 0.d0
    tAmpl%t1MatrixOld = 0.d0
    tAmpl%t2Matrix = 0.d0
    tAmpl%t2MatrixOld = 0.d0
    tAmpl%edMatrix = 0.d0
    tAmpl%intermediate1 = 0.d0
    tAmpl%intermediate2 = 0.d0
    tAmpl%intermediate3 = 0.d0
    tAmpl%intermediate4 = 0.d0
    
    
  END FUNCTION CcdAmplitude_
  
  !
  !     Set up the energy denominator matrix. 
  !
  SUBROUTINE setupEnergyDenomMatrix(this, fockMatr)
    CLASS(CcdAmplitude), INTENT(INOUT) :: this
    TYPE(FockMatrix), INTENT(IN) :: fockMatr
    INTEGER :: nHoles, nParticles
    INTEGER :: particle1, particle2, hole1, hole2
    REAL(DP) :: eneP1, eneP2, eneH1, eneH2
    REAL(DP) :: denom
    
    !     Number of hole states
    nHoles = this%basisSp%nOccupied
    !     Number of particle states
    nParticles = this%basisSp%nUnoccupied
    
    DO particle1=nHoles+1, nHoles+nParticles
       !     Single-particle energy of 
       !     particle 1
       eneP1 = fockMatr%matrix(particle1, particle1)
       
       DO particle2=nHoles+1, nHoles+nParticles
          !     Single-particle energy of
          !     particle 2
          eneP2 = fockMatr%matrix(particle2, particle2)
          
          DO hole1=1, nHoles
             !     Single-particle energy of
             !     hole 1
             eneH1 = fockMatr%matrix(hole1, hole1)
             
             DO hole2=1, nHoles
                !     Single-particle energy of
                !     hole 2
                eneH2 = fockMatr%matrix(hole2, hole2)
                
                !     Energy denominator
                denom = eneH1 + eneH2 - eneP1 - eneP2
                
                this%edMatrix(particle1-nHoles,&
                     particle2-nHoles, hole1, hole2) &
                     = 1.d0/denom
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    
  END SUBROUTINE setupEnergyDenomMatrix
  
  
  !
  !     Intermediate matrix to diagram 3 of
  !     the T2 amplitude equation.
  !
  !     I_ij^kl = <kl||ij> + (1/2)*\sum_cd <kl||cd><cd|t|ij>
  !
  SUBROUTINE updateIntermed1(this, ham)
    CLASS(CcdAmplitude), INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    INTEGER :: nHoles, nParticles
    INTEGER :: hole1, hole2, hole3, hole4
    INTEGER :: particle1, particle2
    REAL(DP) :: interaction, tAmplitude, temp
    
    !     Number of hole states
    nHoles = this%basisSp%nOccupied
    !     Number of particle states
    nParticles = this%basisSp%nUnoccupied
    
    DO hole3=1, nHoles
       DO hole4=1, nHoles
          DO hole1=1, nHoles
             DO hole2=1, nHoles
                
                !     Interaction matrix element
!                temp = ham%vMatrix(hole1, hole2, hole3, hole4) 
                
                DO particle1=1, nParticles
                   DO particle2=particle1+1, nParticles
                      
                      !     Interaction matrix element
 !                     interaction = ham%vMatrix(hole1, hole2,&
 !                          nHoles+particle1, nHoles+particle2)
                      !     CCD t2-amplitude matrix element
                      tAmplitude = this%t2MatrixOld(&
                           particle1, particle2, hole3, hole4)
                      temp = temp + interaction*tAmplitude
                   ENDDO
                ENDDO
                
                this%intermediate1(hole1, hole2, &
                     hole3, hole4) = temp
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    
  END SUBROUTINE updateIntermed1
  
  !
  !     Intermediate matrix to diagram 4 of
  !     the T2 amplitude equation.
  !
  !     I_cj^kb = <kb||cj> + (1/2)*\sum_ld <kl||cd><db|t|lj>     
  !
  SUBROUTINE updateIntermed2(this, ham)
    CLASS(CcdAmplitude), INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    INTEGER :: nHoles, nParticles
    INTEGER :: hole1, hole2, hole3
    INTEGER :: particle1, particle2, particle3
    REAL(DP) :: temp, interaction, t2
    
    !     Number of hole states
    nHoles = this%basisSp%nOccupied
    !     Number of particle states
    nParticles = this%basisSp%nUnoccupied
    
    DO hole2=1, nHoles
       DO particle2=1, nParticles
          DO particle1=1, nParticles
             DO hole1=1, nHoles
                
                !     Interaction matrix element
!                temp = ham%vMatrix(hole1, nHoles+particle1, &
!                     nHoles+particle2, hole2)  
                
                DO hole3=1, nHoles
                   DO particle3=1, nParticles
                      !     Interaction matrix element
!                      interaction = ham%vMatrix(hole1, hole3,&
!                           nHoles+particle2, nHoles+particle3)
                      !     t2-amplitude matrix element
                      t2 = this%t2MatrixOld(&
                           particle3, particle1, hole3, hole2)
                      temp = temp + 0.5d0*interaction*t2
                   ENDDO
                ENDDO
                this%intermediate2(hole1, particle1, &
                     particle2, hole2) = temp
                
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    
  END SUBROUTINE updateIntermed2
  
  !
  !     Intermediate matrix to diagram 5 of
  !     the T2 amplitude equation.
  !
  !     I_c^b = -(1/2)*\sum_kld <kl||cd><bd|t|kl> 
  !
  SUBROUTINE updateIntermed3(this, ham)
    CLASS(CcdAmplitude), INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    INTEGER :: nHoles, nParticles
    INTEGER :: particle1, particle2, particle3
    INTEGER :: hole3, hole4
    REAL(DP) :: temp, interaction, t2
    
    !     Number of hole states
    nHoles = this%basisSp%nOccupied
    !     Number of particle states
    nParticles = this%basisSp%nUnoccupied
    
    DO particle1=1, nParticles
       DO particle2=1, nParticles
          
          temp = 0.d0
          DO hole3=1, nHoles
             DO hole4=hole3+1, nHoles
                DO particle3=1, nParticles
                   !     Interaction matrix element
!                   interaction = ham%vMatrix(hole3, hole4, &
!                        nHoles+particle2, nHoles+particle3)
                   !     t2-amplitude matrix element
                   t2 = this%t2MatrixOld(&
                        particle1, particle3, hole3, hole4)
                   temp = temp - interaction*t2
                   
                ENDDO
             ENDDO
          ENDDO
          this%intermediate3(particle1, particle2) &
               = temp

       ENDDO
    ENDDO
    
    
  END SUBROUTINE updateIntermed3
  
  !
  !     Intermediate matrix to diagram 6 of
  !     the T2 amplitude equation.
  !
  !     I_j^k = -(1/2)*\sum_lcd <kl||cd><cd|t|jl>
  !
  SUBROUTINE updateIntermed4(this, ham)
    CLASS(CcdAmplitude), INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    INTEGER :: nHoles, nParticles
    INTEGER :: hole1, hole2, hole3
    INTEGER :: particle3, particle4
    REAL(DP) :: temp, interaction, t2
    
    !     Number of hole states
    nHoles = this%basisSp%nOccupied
    !     Number of particle states
    nParticles = this%basisSp%nUnoccupied
    
    DO hole1=1, nHoles
       DO hole2=1, nHoles
          
          temp = 0.d0
          DO hole3=1, nHoles
             DO particle3=1, nParticles
                DO particle4=particle3+1, nParticles
                   !     Interaction matrix element
!                   interaction = ham%vMatrix(hole1, hole3,&
!                        nHoles+particle3, nHoles+particle4)
                   !     t2-amplitude matrix element
                   t2 = this%t2MatrixOld(&
                        particle3, particle4, hole2, hole3)
                   temp = temp - interaction*t2
                ENDDO
             ENDDO
          ENDDO
          this%intermediate4(hole1, hole2) = temp
       ENDDO
    ENDDO
    
    
  END SUBROUTINE updateIntermed4
  
  !
  !     Update the CCSD t1-ampitude matrix.
  !
  SUBROUTINE updateT1Matrix(this, ham, fockMatr)
    CLASS(CcdAmplitude), INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    TYPE(FockMatrix), INTENT(IN) :: fockMatr
    INTEGER :: nHoles, nParticles
    INTEGER :: particle1, hole1
    REAL(DP) :: t1MatEle
    
    !     Number of hole states
    nHoles = this%basisSp%nOccupied
    !     Number of particle states
    nParticles = this%basisSp%nUnoccupied 
    
    DO particle1=1, nParticles
       DO hole1=1, nHoles
          !     Calculate a t1 matrix element
          t1MatEle = this%calcT1MatrixElement(ham,&
               fockMatr, particle1, hole1)
          
          !     Store new t-amplitude matrix
          this%t1Matrix(particle1, hole1) &
               = t1MatEle
       ENDDO
    ENDDO
    
    
  END SUBROUTINE updateT1Matrix
  
  !
  !     Update the CCSD t2-amplitude matrix.
  !
  SUBROUTINE updateT2Matrix(this, ham, fockMatr)
    CLASS(CcdAmplitude), INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    TYPE(FockMatrix), INTENT(IN) :: fockMatr
    INTEGER :: nHoles, nParticles
    INTEGER :: particle1, particle2, hole1, hole2
    REAL(DP) :: t2MatEle
    integer :: ml1,ml2,ml3,ml4
    
    !     Number of hole states
    nHoles = this%basisSp%nOccupied
    !     Number of particle states
    nParticles = this%basisSp%nUnoccupied
    
    DO particle1=1, nParticles
       DO particle2=particle1, nParticles
          DO hole1=1, nHoles
             DO hole2=hole1, nHoles
                
                !     Calculate a t2 matrix element
                t2MatEle = this%calcT2MatrixElement(&
                     ham, fockMatr, &
                     particle1, particle2, hole1, hole2)             
                
                !     Store new t-amplitude matrix
                this%t2Matrix(particle1, particle2, &
                     hole1, hole2) = t2MatEle
                this%t2Matrix(particle2, particle1, &
                     hole1, hole2) = -t2MatEle
                this%t2Matrix(particle1, particle2, &
                     hole2, hole1) = -t2MatEle
                this%t2Matrix(particle2, particle1, &
                     hole2, hole1) = t2MatEle
                
                ! ml1 = ham%basisSp%spQuantNr(nHoles+particle1, 2)
                ! ml2 = ham%basisSp%spQuantNr(nHoles+particle2, 2)
                ! ml3 = ham%basisSp%spQuantNr(hole1, 2)
                ! ml4 = ham%basisSp%spQuantNr(hole2, 2)
                ! if (abs(t2MatEle) > 0.000001) then
                !    write(*,*) ' '
                !    write(*,*) 'ml1=',ml1,',ml2=',ml2
                !    write(*,*) 'ml3=',ml3,',ml4=',ml4
                !    write(*,*) 't2MatEle=',t2MatEle
                ! endif
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    
  END SUBROUTINE updateT2Matrix

  !
  !     Calculate a CCSD t1-amplitude matrix element.
  !
  FUNCTION calcT1MatrixElement(this, ham, &
       fockMatr, particle1, hole1) RESULT(t1MatEle)
    CLASS(CcdAmplitude), INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    TYPE(FockMatrix), INTENT(IN) :: fockMatr
    INTEGER, INTENT(IN) :: particle1, hole1
    REAL(DP) :: t1MatEle
    INTEGER :: nHoles, nParticles
    REAL(DP) :: diagram1, diagram2, diagram3
    REAL(DP) :: diagram4, diagram5, diagram6
    REAL(DP) :: diagram7, diagram8, diagram9
    REAL(DP) :: diagram10, diagram11, diagram12
    REAL(DP) :: diagram13, diagram14
    REAL(DP) :: fac, tci, fki, tka
    REAL(DP) :: vkaci, tck, fkc, tacik
    REAL(DP) :: vkacd, tcdki, vklci, tcakl
    REAL(DP) :: tak, tal, tdi, vklcd, tdali 
    INTEGER :: hole2, hole3, particle2, particle3
    REAL(DP) :: eneDenom
    
    !     Number of hole states
    nHoles = this%basisSp%nOccupied
    !     Number of particle states
    nParticles = this%basisSp%nUnoccupied
    
    !     Diagram 1: <a|f|i>
    diagram1 = fockMatr%matrix(nHoles+particle1, hole1)
    
    !     Diagram 2: \sum_c (1-\delta_ac)<a|f|c><c|t|i>
    diagram2 = 0.d0
    DO particle2=1, nParticles
       IF (particle2 /= particle1) THEN
          fac = fockMatr%matrix(nHoles+particle1,&
               nHoles+particle2)
          tci = this%t1MatrixOld(particle2, hole1)
          diagram2 = diagram2 + fac*tci
       ENDIF
    ENDDO

    !     Diagram 3: -\sum_k (1-\delta_ik)<k|f|i><a|t|k>
    diagram3 = 0.d0
    DO hole2=1, nHoles
       IF (hole2 /= hole1) THEN
          fki = fockMatr%matrix(hole2, hole1)
          tka = this%t1MatrixOld(particle1, hole2)
          diagram3 = diagram3 - fki*tka
       ENDIF
    ENDDO
    
    !     Diagram 4: \sum_kc <ka||ci><c|t|k>
    diagram4 = 0.d0
    DO hole2=1, nHoles
       DO particle2=1, nParticles
!          vkaci = ham%vMatrix(&
!               hole2, nHoles+particle1, &
!               nHoles+particle2, hole1)
          tck = this%t1MatrixOld(particle2, hole2)
          diagram4 = diagram4 + vkaci*tck
       ENDDO
    ENDDO
    
    !     Diagram 5: \sum_kc <k|f|c><ac|t|ik>
    diagram5 = 0.d0
    DO hole2=1, nHoles
       DO particle2=1, nParticles
          fkc = fockMatr%matrix(hole2, nHoles+particle2)
          tacik = this%t2MatrixOld(particle1, particle2,&
               hole1, hole2)
          diagram5 = diagram5 + fkc*tacik
       ENDDO
    ENDDO
    
    !     Diagram 6: (1/2)*\sum_kcd <ka||cd><cd|t|ki>
    diagram6 = 0.d0
    DO hole2=1, nHoles
       DO particle2=1, nParticles
          DO particle3=particle2+1, nParticles
!             vkacd = ham%vMatrix(hole2, nHoles+particle1,&
!                  nHoles+particle2, nHoles+particle3)
             tcdki = this%t2MatrixOld(particle2, particle3,&
                  hole2, hole1)
             diagram6 = diagram6 + vkacd*tcdki
          ENDDO
       ENDDO
    ENDDO
    
    !     Diagram 7: -(1/2)*\sum_klc <kl||ci><ca|t|kl>
    diagram7 = 0.d0
    DO hole2=1, nHoles
       DO hole3=hole2+1, nHoles
          DO particle2=1, nParticles
!             vklci = ham%vMatrix(hole2, hole3,&
!                  nHoles+particle2, hole1)
             tcakl = this%t2MatrixOld(particle2, particle1,&
                  hole2, hole3)
             diagram7 = diagram7 - vklci*tcakl
          ENDDO
       ENDDO
    ENDDO
    
    !     Diagram 8: -\sum_kc <k|f|c><c|t|i><a|t|k>
    diagram8 = 0.d0
    DO hole2=1, nHoles
       DO particle2=1, nParticles
          fkc = fockMatr%matrix(hole2, nHoles+particle2)
          tci = this%t1MatrixOld(particle2, hole1)
          tak = this%t1MatrixOld(particle1, hole2)
          diagram8 = diagram8 - fkc*tci*tak
       ENDDO
    ENDDO
    
    !     Diagram 9: -\sum_klc <kl||ci><c|t|k><a|t|l>
    diagram9 = 0.d0
    DO hole2=1, nHoles
       DO hole3=1, nHoles
          DO particle2=1, nParticles
!             vklci = ham%vMatrix(hole2, hole3,&
!                  nHoles+particle2, hole1)
             tck = this%t1MatrixOld(particle2, hole2)
             tal = this%t1MatrixOld(particle1, hole3)
             diagram9 = diagram9 - vklci*tck*tal
          ENDDO
       ENDDO
    ENDDO
    
    !     Diagram 10: \sum_kcd <ka||cd><c|t|k><d|t|i>  
    diagram10 = 0.d0
    DO hole2=1, nHoles
       DO particle2=1, nParticles
          DO particle3=1, nParticles
!             vkacd = ham%vMatrix(hole2, nHoles+particle1,&
!                  nHoles+particle2, nHoles+particle3)
             tck = this%t1MatrixOld(particle2, hole2)
             tdi = this%t1MatrixOld(particle3, hole1)
             diagram10 = diagram10 + vkacd*tck*tdi
          ENDDO
       ENDDO
    ENDDO
    
    !     Diagram 11: -\sum_klcd <kl||cd><c|t|k><d|t|i><a|t|l>
    diagram11 = 0.d0
    DO hole2=1, nHoles
       DO hole3=1, nHoles
          DO particle2=1, nParticles
             DO particle3=1, nParticles
!                vklcd = ham%vMatrix(hole2, hole3,&
!                     nHoles+particle2, nHoles+particle3)
                tck = this%t1MatrixOld(particle2, hole2)
                tdi = this%t1MatrixOld(particle3, hole1)
                tal = this%t1MatrixOld(particle1, hole3)
                diagram11 = diagram11 - vklcd*tck*tdi*tal
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !     Diagram 12: \sum_klcd <kl||cd><c|t|k><da|t|li>
    diagram12 = 0.d0
    DO hole2=1, nHoles
       DO hole3=1, nHoles
          DO particle2=1, nParticles
             DO particle3=1, nParticles
!                vklcd = ham%vMatrix(hole2, hole3,&
!                     nHoles+particle2, nHoles+particle3)
                tck = this%t1MatrixOld(particle2, hole2)
                tdali = this%t2MatrixOld(particle3, particle1,&
                     hole3, hole1)
                diagram12 = diagram12 + vklcd*tck*tdali
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    !     Diagram 13: -(1/2)*\sum_klcd <kl||cd><cd|t|ki><a|t|l>
    diagram13 = 0.d0
    DO hole2=1, nHoles
       DO hole3=1, nHoles
          DO particle2=1, nParticles
             DO particle3=particle2+1, nParticles
!                vklcd = ham%vMatrix(hole2, hole3,&
!                     nHoles+particle2, nHoles+particle3)
                tcdki = this%t2MatrixOld(particle2, particle3,&
                     hole2, hole1)
                tal = this%t1MatrixOld(particle1, hole3)
                diagram13 = diagram13 - vklcd*tcdki*tal
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    !     Diagram 14: -(1/2)*\sum_klcd <kl||cd><ca|t|kl><d|t|i>
    diagram14 = 0.d0
    DO hole2=1, nHoles
       DO hole3=hole2+1, nHoles
          DO particle2=1, nParticles
             DO particle3=1, nParticles
!                vklcd = ham%vMatrix(hole2, hole3,&
!                     nHoles+particle2, nHoles+particle3)
                tcakl = this%t2MatrixOld(particle2, particle1,&
                     hole2, hole3)
                tdi = this%t1MatrixOld(particle3, hole1)
                diagram14 = diagram14 - vklcd*tcakl*tdi
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    eneDenom = fockMatr%matrix(hole1, hole1) &
         - fockMatr%matrix(nHoles+particle1, nHoles+particle1)
    
    t1MatEle = diagram1 + diagram2 + diagram3 &
         + diagram4 + diagram5 + diagram6 &
         + diagram7 + diagram8 + diagram9 &
         + diagram10 + diagram11 + diagram12 &
         + diagram13 + diagram14
    
    t1MatEle = t1MatEle/eneDenom
    
    
  END FUNCTION calcT1MatrixElement
  
  !
  !     Calculate a CCSD t2-amplitude matrix element.
  !
  FUNCTION calcT2MatrixElement(this, ham, fockMatr, &
       particle1, particle2, hole1, hole2) RESULT(t2MatEle)
    CLASS(CcdAmplitude), INTENT(INOUT) :: this
    CLASS(Hamiltonian), TARGET, INTENT(IN) :: ham
    TYPE(FockMatrix), INTENT(IN) :: fockMatr
    INTEGER, INTENT(IN) :: particle1, particle2
    INTEGER, INTENT(IN) :: hole1, hole2
    REAL(DP) :: t2MatEle
    INTEGER :: nHoles, nParticles
    INTEGER :: particle3, particle4, hole3, hole4
    REAL(DP) :: tMatrixElement, interaction, tAmplitude
    REAL(DP) :: diagram1, diagram2, diagram3, diagram4
    REAL(DP) :: diagram5, diagram6, diagram7, diagram8
    REAL(DP) :: diagram9, diagram10, diagram11, diagram12
    REAL(DP) :: diagram13, diagram14, diagram15, diagram16
    REAL(DP) :: diagram17, diagram18, diagram19, diagram20
    REAL(DP) :: diagram21, diagram22, diagram23, diagram24
    REAL(DP) :: diagram25, diagram26, diagram27, diagram28
    REAL(DP) :: diagram29, diagram30
    REAL(DP) :: intermediate, energyDenom, t1, t2, fock
    REAL(DP) :: t11, t12, t13, t14
    
    !     Number of hole states
    nHoles = this%basisSp%nOccupied
    !     Number of particle states
    nParticles = this%basisSp%nUnoccupied
    
    !     Diagram 1: <ab||ij>
!    diagram1 = ham%vMatrix(nHoles+particle1,&
!         nHoles+particle2, hole1, hole2)
    
    !     Diagram 2: (1/2)*\sum_cd <ab||cd><cd|t|ij>
    diagram2 = 0.d0
    DO particle3=1, nParticles
       DO particle4=particle3+1, nParticles
          
!          interaction = ham%vMatrix(&
!               nHoles+particle1, nHoles+particle2,&
!               nHoles+particle3, nHoles+Particle4) 
          tAmplitude = this%t2MatrixOld(&
               particle3, particle4, hole1, hole2)
          
          diagram2 = diagram2 + interaction*tAmplitude
       ENDDO
    ENDDO
    
    !     Diagram 3: (1/2)*\sum_kl [<kl||ij>
    !                + (1/2)*\sum_cd <kl||cd><cd|t|ij>]
    !                *<ab|t|kl>
    diagram3 = 0.d0
    DO hole3=1, nHoles
       DO hole4=hole3+1, nHoles       
          
          tAmplitude = this%t2MatrixOld(&
               particle1, particle2, &
               hole3, hole4)
          intermediate = this%intermediate1(&
               hole3, hole4, hole1, hole2)
          
          diagram3 = diagram3 + tAmplitude*intermediate
       ENDDO
    ENDDO
    
    !     Diagram 4: P(ij)P(ab)\sum_kc I2_cj^kb * <ac|t|ik>
    diagram4 = 0.d0
    DO hole3=1, nHoles
       DO particle3=1, nParticles
          !     f(a,b)f(i,j)
          t2 = this%t2MatrixOld(particle1, particle3,&
               hole1, hole3)
          intermediate = this%intermediate2(&
               hole3, particle2, particle3, hole2)
          diagram4 = diagram4 + intermediate*t2
          
          !     -f(b,a)f(i,j)    
          t2 = this%t2MatrixOld(particle2, particle3,&
               hole1, hole3)
          intermediate = this%intermediate2(&
               hole3, particle1, particle3, hole2)
          diagram4 = diagram4 - intermediate*t2
          
          !     -f(a,b)f(j,i)
          t2 = this%t2MatrixOld(particle1, particle3,&
               hole2, hole3)
          intermediate = this%intermediate2(&
               hole3, particle2, particle3, hole1)
          
          diagram4 = diagram4 - intermediate*t2
          
          !     f(b,a)f(j,i)
          t2 = this%t2MatrixOld(particle2, particle3,&
               hole2, hole3)
          intermediate = this%intermediate2(&
               hole3, particle1, particle3, hole1)
          diagram4 = diagram4 + intermediate*t2
       ENDDO
    ENDDO
    
    !     Diagram 5: P(ab)\sum_c I_c^b * <ac|t|ij>
    diagram5 = 0.d0
    DO particle3=1, nParticles
       !     f(a,b)
       intermediate = this%intermediate3(&
            particle2, particle3)
       t2 = this%t2MatrixOld(particle1, particle3,&
            hole1, hole2)
       diagram5 = diagram5 + intermediate*t2
    
       !     -f(b,a)
       intermediate = this%intermediate3(&
            particle1, particle3)
       t2 = this%t2MatrixOld(particle2, particle3,&
            hole1, hole2)
       diagram5 = diagram5 - intermediate*t2
       
    ENDDO
  
    
    !     Diagram 6: P(ij)\sum_k I_j^k * <ab|t|ik>
    diagram6 = 0.d0
    DO hole3=1, nHoles
       !     f(i,j)
       intermediate = this%intermediate4(&
            hole3, hole2)
       t2 = this%t2MatrixOld(particle1, particle2,&
            hole1, hole3)
       diagram6 = diagram6 + intermediate*t2

       !     -f(j,i)
       intermediate = this%intermediate4(&
            hole3, hole1)
       t2 = this%t2MatrixOld(particle1, particle2,&
            hole2, hole3)
       diagram6 = diagram6 - intermediate*t2
       
    ENDDO
    
    ! !     Diagram 7: P(ab)\sum_c (1-\delta_bc)
    ! !                *<b|f|c><ac|t|ij>
    ! diagram7 = 0.d0
    ! DO particle3=1, nParticles
    !    IF (particle3 /= particle2) THEN
    !       !     f(ab)
    !       fock = fockMatr%matrix(nHoles+particle2,&
    !            nHoles+particle3)
    !       t2 = this%t2MatrixOld(particle1, particle3,&
    !            hole1, hole2)
    !       diagram7 = diagram7 + fock*t2
    !    ENDIF
    !    IF (particle3 /= particle1) THEN
    !       !     -f(ba)
    !       fock = fockMatr%matrix(nHoles+particle1,&
    !            nHoles+particle3)
    !       t2 = this%t2MatrixOld(particle2, particle3,&
    !            hole1, hole2)
    !       diagram7 = diagram7 - fock*t2
    !    ENDIF
    ! ENDDO
    
    ! !     Diagram 8: -P(ij)\sum_k (1-\delta_kj)
    ! !                * <k|f|j><ab|t|ik>
    ! diagram8 = 0.d0
    ! DO hole3=1, nHoles
    !    IF (hole3 /= hole2) THEN
    !       !     f(ij)
    !       fock = fockMatr%matrix(hole3, hole2)
    !       t2 = this%t2MatrixOld(particle1, particle2,&
    !            hole1, hole3)
    !       diagram8 = diagram8 - fock*t2
    !    ENDIF
    !    IF (hole3 /= hole1) THEN
    !       !     -f(ji)
    !       fock = fockMatr%matrix(hole3, hole1)
    !       t2 = this%t2MatrixOld(particle1, particle2,&
    !            hole2, hole3)
    !       diagram8 = diagram8 + fock*t2
    !    ENDIF
    ! ENDDO
    
    ! !     Diagram 9: P(ij)\sum_c <ab||cj><c|t|i>
    ! diagram9 = 0.d0
    ! DO particle3=1, nParticles
    !    !     f(ij)
    !    interaction = ham%vMatrix(&
    !         nHoles+particle1, nHoles+particle2,&
    !         nHoles+particle3, hole2)
    !    t1 = this%t1MatrixOld(particle3, hole1)
    !    diagram9 = diagram9 + interaction*t1

    !    !     -f(ji)
    !    interaction = ham%vMatrix(&
    !         nHoles+particle1, nHoles+particle2,&
    !         nHoles+particle3, hole1)
    !    t1 = this%t1MatrixOld(particle3, hole2)
    !    diagram9 = diagram9 - interaction*t1
    ! ENDDO
    
    ! !     Diagram 10: -P(ab)\sum_k <kb||ij><a|t|k>
    ! diagram10 = 0.d0
    ! DO hole3=1, nHoles
    !    !     f(ab)
    !    interaction = ham%vMatrix(&
    !         hole3, nHoles+particle2,&
    !         hole1, hole2)
    !    t1 = this%t1MatrixOld(particle1, hole3)
    !    diagram10 = diagram10 - interaction*t1
       
    !    !     -f(ba)
    !    interaction = ham%vMatrix(&
    !         hole3, nHoles+particle1,&
    !         hole1, hole2)
    !    t1 = this%t1MatrixOld(particle2, hole3)
    !    diagram10 = diagram10 + interaction*t1
    ! ENDDO
    
    ! !     Diagram 11: (1/2)*P(ab)\sum_kl <kl||ij>
    ! !                 *<a|t|k><b|t|l>
    ! diagram11 = 0.d0
    ! DO hole3=1, nHoles
    !    DO hole4=1, nHoles
    !       !     f(ab)
    !       interaction = ham%vMatrix(&
    !            hole3, hole4, hole1, hole2)
    !       t11 = this%t1MatrixOld(particle1, hole3)
    !       t12 = this%t1MatrixOld(particle2, hole4)
    !       diagram11 = diagram11 &
    !            + 0.5d0*interaction*t11*t12
          
    !       !     -f(ba)
    !       interaction = ham%vMatrix(&
    !            hole3, hole4, hole1, hole2)
    !       t11 = this%t1MatrixOld(particle2, hole3)
    !       t12 = this%t1MatrixOld(particle1, hole4)
    !       diagram11 = diagram11 &
    !            - 0.5d0*interaction*t11*t12
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 12: (1/2)*P(ij)\sum_cd <ab||cd>
    ! !                 *<c|t|i><d|t|j>
    ! diagram12 = 0.d0
    ! DO particle3=1, nParticles
    !    DO particle4=1, nParticles
    !       !     f(ij)
    !       interaction = ham%vMatrix(&
    !            nHoles+particle1, nHoles+particle2,&
    !            nHoles+particle3, nHoles+particle4)
    !       t11 = this%t1MatrixOld(particle3, hole1)
    !       t12 = this%t1MatrixOld(particle4, hole2)
    !       diagram12 = diagram12 &
    !            + 0.5d0*interaction*t11*t12 
          
    !       !     -f(ji)
    !       interaction = ham%vMatrix(&
    !            nHoles+particle1, nHoles+particle2,&
    !            nHoles+particle3, nHoles+particle4)
    !       t11 = this%t1MatrixOld(particle3, hole2)
    !       t12 = this%t1MatrixOld(particle4, hole1)
    !       diagram12 = diagram12 &
    !            - 0.5d0*interaction*t11*t12
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 13: -P(ij)P(ab)\sum_kc <kb||ic>
    ! !                 *<a|t|k><c|t|j>
    ! diagram13 = 0.d0
    ! DO hole3=1, nHoles
    !    DO particle3=1, nParticles
    !       !     f(ab)f(ij)
    !       interaction = ham%vMatrix(&
    !            hole3, nHoles+particle2,&
    !            hole1, nHoles+particle3)
    !       t11 = this%t1MatrixOld(particle1, hole3)
    !       t12 = this%t1MatrixOld(particle3, hole2)
    !       diagram13 = diagram13 - interaction*t11*t12
          
    !       !     -f(ba)f(ij)
    !       interaction = ham%vMatrix(&
    !            hole3, nHoles+particle1,&
    !            hole1, nHoles+particle3)
    !       t11 = this%t1MatrixOld(particle2, hole3)
    !       t12 = this%t1MatrixOld(particle3, hole2)
    !       diagram13 = diagram13 + interaction*t11*t12
          
    !       !     -f(ab)f(ji)
    !       interaction = ham%vMatrix(&
    !            hole3, nHoles+particle2,&
    !            hole2, nHoles+particle3)
    !       t11 = this%t1MatrixOld(particle1, hole3)
    !       t12 = this%t1MatrixOld(particle3, hole1)
    !       diagram13 = diagram13 + interaction*t11*t12
          
    !       !     f(ba)f(ji)
    !       interaction = ham%vMatrix(&
    !            hole3, nHoles+particle1,&
    !            hole2, nHoles+particle3)
    !       t11 = this%t1MatrixOld(particle2, hole3)
    !       t12 = this%t1MatrixOld(particle3, hole1)
    !       diagram13 = diagram13 - interaction*t11*t12
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 14: P(ab)\sum_kc <k|f|c><a|t|k>
    ! !                 *<bc|t|ij>
    ! diagram14 = 0.d0
    ! DO hole3=1, nHoles
    !    DO particle3=1, nParticles
    !       !     f(ab)
    !       fock = fockMatr%matrix(hole3, nHoles+particle3)
    !       t1 = this%t1MatrixOld(particle1, hole3)
    !       t2 = this%t2MatrixOld(particle2, particle3,&
    !            hole1, hole2)
    !       diagram14 = diagram14 + fock*t1*t2

    !       !     -f(ba)
    !       fock = fockMatr%matrix(hole3, nHoles+particle3)
    !       t1 = this%t1MatrixOld(particle2, hole3)
    !       t2 = this%t2MatrixOld(particle1, particle3,&
    !            hole1, hole2)
    !       diagram14 = diagram14 - fock*t1*t2
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 15: P(ij)\sum_kc <k|f|c><c|t|i>
    ! !                 *<ab|t|jk>
    ! diagram15 = 0.d0
    ! DO hole3=1, nHoles
    !    DO particle3=1, nParticles
    !       !     f(ij)
    !       fock = fockMatr%matrix(hole3, nHoles+particle3)
    !       t1 = this%t1MatrixOld(particle3, hole1)
    !       t2 = this%t2MatrixOld(particle1, particle2,&
    !            hole2, hole3)
    !       diagram15 = diagram15 + fock*t1*t2
          
    !       !     -f(ji)
    !       fock = fockMatr%matrix(hole3, nHoles+particle3)
    !       t1 = this%t1MatrixOld(particle3, hole2)
    !       t2 = this%t2MatrixOld(particle1, particle2,&
    !            hole1, hole3)
    !       diagram15 = diagram15 - fock*t1*t2
    !    ENDDO
    ! ENDDO

    ! !     Diagram 16: -P(ij)\sum_klc <kl||ci><c|t|k>
    ! !                 *<ab|t|lj>
    ! diagram16 = 0.d0
    ! DO hole3=1, nHoles
    !    DO hole4=1, nHoles
    !       DO particle3=1, nParticles
    !          !     f(ij)
    !          interaction = ham%vMatrix(hole3, hole4,&
    !               nHoles+particle3, hole1)
    !          t1 = this%t1MatrixOld(particle3, hole3)
    !          t2 = this%t2MatrixOld(particle1, particle2,&
    !               hole4, hole2)
    !          diagram16 = diagram16 - interaction*t1*t2
             
    !          !     -f(ji)
    !          interaction = ham%vMatrix(hole3, hole4,&
    !               nHoles+particle3, hole2)
    !          t1 = this%t1MatrixOld(particle3, hole3)
    !          t2 = this%t2MatrixOld(particle1, particle2,&
    !               hole4, hole1)
    !          diagram16 = diagram16 + interaction*t1*t2
    !       ENDDO
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 17: P(ab)\sum_kcd <ka||cd><c|t|k>
    ! !                 * <db|t|ij>
    ! diagram17 = 0.d0
    ! DO hole3=1, nHoles
    !    DO particle3=1, nParticles
    !       DO particle4=1, nParticles
    !          !     f(ab)
    !          interaction = ham%vMatrix(hole3, nHoles+particle1,&
    !               nHoles+particle3, nHoles+particle4)
    !          t1 = this%t1MatrixOld(particle3, hole3)
    !          t2 = this%t2MatrixOld(particle4, particle2,&
    !               hole1, hole2)
    !          diagram17 = diagram17 + interaction*t1*t2
             
    !          !     -f(ab)
    !          interaction = ham%vMatrix(hole3, nHoles+particle2,&
    !               nHOles+particle3, nHoles+particle4)
    !          t1 = this%t1MatrixOld(particle3, hole3)
    !          t2 = this%t2MatrixOld(particle4, particle1,&
    !               hole1, hole2)
    !          diagram17 = diagram17 - interaction*t1*t2
    !       ENDDO
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 18: P(ij)P(ab)\sum_kcd <ak||dc><d|t|i>
    ! !                 *<bc|t|jk>
    ! diagram18 = 0.d0
    ! DO hole3=1, nHoles
    !    DO particle3=1, nParticles
    !       DO particle4=1, nParticles
    !          !     f(ab)f(ij)
    !          interaction = ham%vMatrix(&
    !               nHoles+particle1, hole3,&
    !               nHoles+particle4, nHoles+particle3)
    !          t1 = this%t1MatrixOld(particle4, hole1)
    !          t2 = this%t2MatrixOld(particle2, particle3,&
    !               hole2, hole3)
    !          diagram18 = diagram18 + interaction*t1*t2
             
    !          !     -f(ba)f(ij)
    !          interaction = ham%vMatrix(&
    !               nHoles+particle2, hole3,&
    !               nHoles+particle4, nHoles+particle3)
    !          t1 = this%t1MatrixOld(particle4, hole1)
    !          t2 = this%t2MatrixOld(particle1, particle3,&
    !               hole2, hole3)
    !          diagram18 = diagram18 - interaction*t1*t2
             
    !          !     -f(ab)f(ji)
    !          interaction = ham%vMatrix(&
    !               nHoles+particle1, hole3,&
    !               nHoles+particle4, nHoles+particle3)
    !          t1 = this%t1MatrixOld(particle4, hole2)
    !          t2 = this%t2MatrixOld(particle2, particle3,&
    !               hole1, hole3)
    !          diagram18 = diagram18 - interaction*t1*t2

    !          !     f(ba)f(ji)
    !          interaction = ham%vmatrix(&
    !               nHoles+particle2, hole3,&
    !               nHoles+particle4, nHoles+particle3)
    !          t1 = this%t1MatrixOld(particle4, hole2)
    !          t2 = this%t2MatrixOld(particle1, particle3,&
    !               hole1, hole3)
    !          diagram18 = diagram18 + interaction*t1*t2
    !       ENDDO
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 19: P(ij)P(ab)\sum_klc <kl||ic><a|t|l>
    ! !                 *<bc|t|jk>
    ! diagram19 = 0.d0
    ! DO hole3=1, nHoles
    !    DO hole4=1, nHoles
    !       DO particle3=1, nParticles
    !          !     f(ab)f(ij)
    !          interaction = ham%vMatrix(hole3, hole4,&
    !               hole1, nHoles+particle3)
    !          t1 = this%t1MatrixOld(particle1, hole4)
    !          t2 = this%t2MatrixOld(particle2, particle3,&
    !               hole2, hole3)
    !          diagram19 = diagram19 + interaction*t1*t2
             
    !          !     -f(ba)f(ij)
    !          interaction = ham%vMatrix(hole3, hole4,&
    !               hole1, nHoles+particle3)
    !          t1 = this%t1MatrixOld(particle2, hole4)
    !          t2 = this%t2MatrixOld(particle1, particle3,&
    !               hole2, hole3)
    !          diagram19 = diagram19 - interaction*t1*t2
             
    !          !     -f(ab)f(ji)
    !          interaction = ham%vMatrix(hole3, hole4,&
    !               hole2, nHoles+particle3)
    !          t1 = this%t1MatrixOld(particle1, hole4)
    !          t2 = this%t2MatrixOld(particle2, particle3,&
    !               hole1, hole3)
    !          diagram19 = diagram19 - interaction*t1*t2
             
    !          !     f(ba)f(ji)
    !          interaction = ham%vMatrix(hole3, hole4,&
    !               hole2, nHoles+particle3)
    !          t1 = this%t1MatrixOld(particle2, hole4)
    !          t2 = this%t2MatrixOld(particle1, particle3,&
    !               hole1, hole3)
    !          diagram19 = diagram19 + interaction*t1*t2
    !       ENDDO
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 20: (1/2)*P(ij)\sum_klc <kl||cj><c|t|i>
    ! !                 *<ab|t|kl>
    ! diagram20 = 0.d0
    ! DO hole3=1, nHoles
    !    DO hole4=1, nHoles
    !       DO particle3=1, nParticles
    !          !     f(ij)
    !          interaction = ham%vMatrix(hole3, hole4,&
    !               nHoles+particle3, hole2)
    !          t1 = this%t1MatrixOld(particle3, hole1)
    !          t2 = this%t2MatrixOld(particle1, particle2,&
    !               hole3, hole4)
    !          diagram20 = diagram20 + 0.5d0*interaction*t1*t2
             
    !          !     -f(ji)
    !          interaction = ham%vMatrix(hole3, hole4,&
    !               nHoles+particle3, hole1)
    !          t1 = this%t1MatrixOld(particle3, hole2)
    !          t2 = this%t2MatrixOld(particle1, particle2,&
    !               hole3, hole4)
    !          diagram20 = diagram20 - 0.5d0*interaction*t1*t2
    !       ENDDO
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 21: -(1/2)*P(ab)\sum_kcd <kb||cd><a|t|k>
    ! !                 *<cd|t|ij>
    ! diagram21 = 0.d0
    ! DO hole3=1, nHoles
    !    DO particle3=1, nParticles
    !       DO particle4=particle3+1, nParticles
    !          !     f(ab)
    !          interaction = ham%vMatrix(hole3, nHoles+particle2,&
    !               nHoles+particle3, nHoles+particle4)
    !          t1 = this%t1MatrixOld(particle1, hole3)
    !          t2 = this%t2MatrixOld(particle3, particle4,&
    !               hole1, hole2)
    !          diagram21 = diagram21 - interaction*t1*t2
             
    !          !     -f(ba)
    !          interaction= ham%vMatrix(hole3, nHoles+particle1,&
    !               nHoles+particle3, nHoles+particle4)
    !          t1 = this%t1MatrixOld(particle2, hole3)
    !          t2 = this%t2MatrixOld(particle3, particle4,&
    !               hole1, hole2)
    !          diagram21 = diagram21 + interaction*t1*t2
    !       ENDDO
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 22: -(1/2)*P(ij)P(ab)\sum_kcd <kb||cd><c|t|i>
    ! !                *<a|t|k><d|t|j>
    ! diagram22 = 0.d0
    ! DO hole3=1, nHoles
    !    DO particle3=1, nParticles
    !       DO particle4=1, nParticles
    !          !     f(ab)f(ij)
    !          interaction = ham%vMatrix(hole3, nHoles+particle2,&
    !               nHoles+particle3, nHoles+particle4)
    !          t11 = this%t1MatrixOld(particle3, hole1)
    !          t12 = this%t1MatrixOld(particle1, hole3)
    !          t13 = this%t1MatrixOld(particle4, hole2)
    !          diagram22 = diagram22 &
    !               - 0.5d0*interaction*t11*t12*t13
             
    !          !     -f(ba)f(ij)
    !          interaction = ham%vMatrix(hole3, nHoles+particle1,&
    !               nHoles+particle3, nHoles+particle4)
    !          t11 = this%t1MatrixOld(particle3, hole1)
    !          t12 = this%t1MatrixOld(particle2, hole3)
    !          t13 = this%t1MatrixOld(particle4, hole2)
    !          diagram22 = diagram22 &
    !               + 0.5d0*interaction*t11*t12*t13
             
    !          !     -f(ab)f(ji)
    !          interaction = ham%vMatrix(hole3, nHoles+particle2,&
    !               nHoles+particle3, nHoles+particle4)
    !          t11 = this%t1MatrixOld(particle3, hole2)
    !          t12 = this%t1MatrixOld(particle1, hole3)
    !          t13 = this%t1MatrixOld(particle4, hole1)
    !          diagram22 = diagram22 &
    !               + 0.5d0*interaction*t11*t12*t13

    !          !     f(ba)f(ji)
    !          interaction = ham%vMatrix(hole3, nHoles+particle1,&
    !               nHoles+particle3, nHoles+particle4)
    !          t11 = this%t1MatrixOld(particle3, hole2)
    !          t12 = this%t1MatrixOld(particle2, hole3)
    !          t13 = this%t1MatrixOld(particle4, hole1)
    !          diagram22 = diagram22 &
    !               - 0.5d0*interaction*t11*t12*t13
    !       ENDDO
    !    ENDDO
    ! ENDDO

    ! !     Diagram 23: (1/2)*P(ij)P(ab)\sum_klc <kl||cj>
    ! !                 *<c|t|i><a|t|k><b|t|l>
    ! diagram23 = 0.d0
    ! DO hole3=1, nHoles
    !    DO hole4=1, nHoles
    !       DO particle3=1, nParticles
    !          !     f(ab)f(ij)
    !          interaction = ham%vMatrix(hole3, hole4,&
    !               nHoles+particle3, hole2)
    !          t11 = this%t1MatrixOld(particle3, hole1)
    !          t12 = this%t1MatrixOld(particle1, hole3)
    !          t13 = this%t1MatrixOld(particle2, hole4)
    !          diagram23 = diagram23 &
    !               + 0.5d0*interaction*t11*t12*t13
             
    !          !     -f(ba)f(ij)
    !          interaction = ham%vMatrix(hole3, hole4,&
    !               nHoles+particle3, hole2)
    !          t11 = this%t1MatrixOld(particle3, hole1)
    !          t12 = this%t1MatrixOld(particle2, hole3)
    !          t13 = this%t1MatrixOld(particle1, hole4)
    !          diagram23 = diagram23 &
    !               - 0.5d0*interaction*t11*t12*t13
             
    !          !     -f(ab)f(ji)
    !          interaction = ham%vMatrix(hole3, hole4,&
    !               nHoles+particle3, hole1)
    !          t11 = this%t1MatrixOld(particle3, hole2)
    !          t12 = this%t1MatrixOld(particle1, hole3)
    !          t13 = this%t1MatrixOld(particle2, hole4)
    !          diagram23 = diagram23 &
    !               - 0.5d0*interaction*t11*t12*t13
             
    !          !     f(ba)f(ji)
    !          interaction = ham%vMatrix(hole3, hole4,&
    !               nHoles+particle3, hole1)
    !          t11 = this%t1MatrixOld(particle3, hole2)
    !          t12 = this%t1MatrixOld(particle2, hole3)
    !          t13 = this%t1MatrixOld(particle1, hole4)
    !          diagram23 = diagram23 &
    !               + 0.5d0*interaction*t11*t12*t13
    !       ENDDO
    !    ENDDO
    ! ENDDO

    ! !     Diagram 24: -P(ij)\sum_klcd <kl||cd>
    ! !                 *<c|t|k><d|t|i><ab|t|lj>
    ! diagram24 = 0.d0
    ! DO hole3=1, nHoles
    !    DO hole4=1, nHoles
    !       DO particle3=1, nParticles
    !          DO particle4=1, nParticles
    !             !     f(ij)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole3)
    !             t12 = this%t1MatrixOld(particle4, hole1)
    !             t2 = this%t2MatrixOld(particle1, particle2,&
    !                  hole4, hole2)
    !             diagram24 = diagram24 &
    !                  - interaction*t11*t12*t2
                
    !             !     -f(ji)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole3)
    !             t12 = this%t1MatrixOld(particle4, hole2)
    !             t2 = this%t2MatrixOld(particle1, particle2,&
    !                  hole4, hole1)
    !             diagram24 = diagram24 &
    !                  + interaction*t11*t12*t2
    !          ENDDO
    !       ENDDO
    !    ENDDO
    ! ENDDO

    ! !     Diagram 25: -P(ab)\sum_klcd <kl||cd><c|t|k>
    ! !                 *<a|t|l><db|t|ij>
    ! diagram25 = 0.d0
    ! DO hole3=1, nHoles
    !    DO hole4=1, nHoles
    !       DO particle3=1, nParticles
    !          DO particle4=1, nParticles
    !             !     f(ab)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole3)
    !             t12 = this%t1MatrixOld(particle1, hole4)
    !             t2 = this%t2MatrixOld(particle4, particle2,&
    !                  hole1, hole2)
    !             diagram25 = diagram25 &
    !                  - interaction*t11*t12*t2
                
    !             !     -f(ba)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole3)
    !             t12 = this%t1MatrixOld(particle2, hole4)
    !             t2 = this%t2MatrixOld(particle4, particle1,&
    !                  hole1, hole2)
    !             diagram25 = diagram25 &
    !                  + interaction*t11*t12*t2
    !          ENDDO
    !       ENDDO
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 26: (1/4)*P(ij)\sum_klcd <kl||cd>
    ! !                 *<c|t|i><d|t|j><ab|t|kl>
    ! diagram26 = 0.d0
    ! DO hole3=1, nHoles
    !    DO hole4=1, nHoles
    !       DO particle3=1, nParticles
    !          DO particle4=1, nParticles
    !             !     f(ij)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole1)
    !             t12 = this%t1MatrixOld(particle4, hole2)
    !             t2 = this%t2MatrixOld(particle1, particle2,&
    !                  hole3, hole4)
    !             diagram26 = diagram26 &
    !                  + 0.25d0*interaction*t11*t12*t2
                
    !             !     -f(ji)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole2)
    !             t12 = this%t1MatrixOld(particle4, hole1)
    !             t2 = this%t2MatrixOld(particle1, particle2,&
    !                  hole3, hole4)
    !             diagram26 = diagram26 &
    !                  - 0.25d0*interaction*t11*t12*t2
    !          ENDDO
    !       ENDDO
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 27: (1/4)*P(ab)\sum_klcd <kl||cd>
    ! !                 *<a|t|k><b|t|l><cd|t|ij>
    ! diagram27 = 0.d0
    ! DO hole3=1, nHoles
    !    DO hole4=1, nHoles
    !       DO particle3=1, nParticles
    !          DO particle4=particle3+1, nParticles
    !             !     f(ab)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle1, hole3)
    !             t12 = this%t1MatrixOld(particle2, hole4)
    !             t2 = this%t2MatrixOld(particle3, particle4,&
    !                  hole1, hole2)
    !             diagram27 = diagram27 &
    !                  + 0.5d0*interaction*t11*t12*t2
                
    !             !     -f(ba)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle2, hole3)
    !             t12 = this%t1MatrixOld(particle1, hole4)
    !             t2 = this%t2MatrixOld(particle3, particle4,&
    !                  hole1, hole2)
    !             diagram27 = diagram27 &
    !                  - 0.5d0*interaction*t11*t12*t2
    !          ENDDO
    !       ENDDO
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 28: P(ij)P(ab)\sum_klcd <kl||cd>
    ! !                 *<c|t|i><b|t|l><ad|t|kj>
    ! diagram28 = 0.d0
    ! DO hole3=1, nHoles
    !    DO hole4=1, nHoles
    !       DO particle3=1, nParticles
    !          DO particle4=1, nParticles
    !             !     f(ab)f(ij)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole1)
    !             t12 = this%t1MatrixOld(particle2, hole4)
    !             t2 = this%t2MatrixOld(particle1, particle4,&
    !                  hole3, hole2)
    !             diagram28 = diagram28 &
    !                  + interaction*t11*t12*t2
                
    !             !     -f(ba)f(ij)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole1)
    !             t12 = this%t1MatrixOld(particle1, hole4)
    !             t2 = this%t2MatrixOld(particle2, particle4,&
    !                  hole3, hole2)
    !             diagram28 = diagram28 &
    !                  - interaction*t11*t12*t2

    !             !     -f(ab)f(ji)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole2)
    !             t12 = this%t1MatrixOld(particle2, hole4)
    !             t2 = this%t2MatrixOld(particle1, particle4,&
    !                  hole3, hole1)
    !             diagram28 = diagram28 &
    !                  - interaction*t11*t12*t2

    !             !     f(ba)f(ji)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole2)
    !             t12 = this%t1MatrixOld(particle1, hole4)
    !             t2 = this%t2MatrixOld(particle2, particle4,&
    !                  hole3, hole1)
    !             diagram28 = diagram28 &
    !                  + interaction*t11*t12*t2
    !          ENDDO
    !       ENDDO
    !    ENDDO
    ! ENDDO
    
    ! !     Diagram 29: (1/4)*P(ij)P(ab)\sum_klcd <kl||cd>
    ! !                 *<c|t|i><a|t|k><d|t|j><b|t|l>
    ! diagram29 = 0.d0
    ! DO hole3=1, nHoles
    !    DO hole4=1, nHoles
    !       DO particle3=1, nParticles
    !          DO particle4=1, nParticles
    !             !     f(ab)f(ij)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole1)
    !             t12 = this%t1MatrixOld(particle1, hole3)
    !             t13 = this%t1MatrixOld(particle4, hole2)
    !             t14 = this%t1MatrixOld(particle2, hole4)
    !             diagram29 = diagram29 &
    !                  + 0.25d0*interaction*t11*t12*t13*t14
                
    !             !     -f(ba)f(ij)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole1)
    !             t12 = this%t1MatrixOld(particle2, hole3)
    !             t13 = this%t1MatrixOld(particle4, hole2)
    !             t14 = this%t1MatrixOld(particle1, hole4)
    !             diagram29 = diagram29 &
    !                  - 0.25d0*interaction*t11*t12*t13*t14
                
    !             !     -f(ab)f(ji)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole2)
    !             t12 = this%t1MatrixOld(particle1, hole3)
    !             t13 = this%t1MatrixOld(particle4, hole1)
    !             t14 = this%t1MatrixOld(particle2, hole4)
    !             diagram29 = diagram29 &
    !                  - 0.25d0*interaction*t11*t12*t13*t14
                
    !             !     f(ba)f(ji)
    !             interaction = ham%vMatrix(hole3, hole4,&
    !                  nHoles+particle3, nHoles+particle4)
    !             t11 = this%t1MatrixOld(particle3, hole2)
    !             t12 = this%t1MatrixOld(particle2, hole3)
    !             t13 = this%t1MatrixOld(particle4, hole1)
    !             t14 = this%t1MatrixOld(particle1, hole4)
    !             diagram29 = diagram29 &
    !                  + 0.25d0*interaction*t11*t12*t13*t14
    !          ENDDO
    !       ENDDO
    !    ENDDO
    ! ENDDO
    
    !     Collect different diagrams to the
    !     t-amplitude matrix element
    ! tMatrixElement = diagram1 + diagram2 + diagram3 &
    !      + diagram4 + diagram5 + diagram6 + diagram7 &
    !      + diagram8 + diagram9 + diagram10 + diagram11 &
    !      + diagram12 + diagram13 + diagram14 + diagram15 &
    !      + diagram16 + diagram17 + diagram18 + diagram19 &
    !      + diagram20 + diagram21 + diagram22 + diagram23 &
    !      + diagram24 + diagram25 + diagram26 + diagram27 &
    !      + diagram28 + diagram29
    
    tMatrixElement = diagram1 + diagram2 + diagram3 &
         + diagram4 + diagram5 + diagram6 !+ diagram7 !&
    !+ diagram8 + diagram9 + diagram10 
    !tMatrixElement = diagram1 + diagram2 + diagram3 &
    !     + diagram4 + diagram5 + diagram6
    
    !    write(*,*) 'diagram1=',diagram1
    !    write(*,*) 'diagram2=',diagram2
    !    write(*,*) 'diagram3=',diagram3
    !    write(*,*) 'diagram4=',diagram4
    !write(*,*) 'diagram5=',diagram5
    !write(*,*) 'diagram6=',diagram6
    !write(*,*) 'diagram7=',diagram7
    !write(*,*) 'diagram8=',diagram8
    !write(*,*) 'diagram9=',diagram9
    !write(*,*) 'diagram10=',diagram10
    
    !     1/(energy denominator)
    energyDenom = this%edMatrix(particle1, particle2, &
         hole1, hole2)
    
    t2MatEle = tMatrixElement*energyDenom
    
    
  END FUNCTION calcT2MatrixElement
  
  
  
  
  
END MODULE CcdAmplitudeMod


  
