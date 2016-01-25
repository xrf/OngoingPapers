
!
!     Module for tests.
!
MODULE TestsMod
  USE ConstantsMod
  IMPLICIT NONE
  
  
  
CONTAINS


  SUBROUTINE pt2Sheperd
    USE SpBasisPW3dMod
    USE SpBasisMod
    CLASS(SpBasis), POINTER :: basis
    CLASS(SpBasisPW3d), POINTER :: mySpBasisPW3d => NULL()
    CHARACTER(LEN=100) :: spinState
    REAL(DP), ALLOCATABLE :: fockVector(:)
    INTEGER :: nOccupied, nShells, nOrbits
    INTEGER :: nix, niy, niz, njx, njy, njz
    INTEGER :: nax, nay, naz, nbx, nby, nbz
    INTEGER :: i, j, a, b, nxkf, nykf, nzkf
    INTEGER :: p, npx, npy, npz
    REAL(DP) :: volume, area, rs, via, vja, temp
    REAL(DP) :: factor, kf2, kf
    REAL(DP) :: ki2, ki, kj2, kj, ka2, ka, kb2, kb
    REAL(DP) :: denom, length, v
    
    write(*,*) 'in pt2'
    nOccupied = 14
    nShells = 6
    rs = 1.d0
    volume = nOccupied*4.d0*pi*rs**3/3.d0
    area = volume**(2.d0/3.d0)
    length = volume**(1.d0/3.d0)
    factor = 2.d0*pi/length
    spinState = 'normal'
    ALLOCATE(basis, SOURCE=mySpBasisPW3d)
    basis => SpBasisPW3d_(nOccupied, nShells, &
         spinState, volume)
    
    CALL basis%setupSpBasis()
    nOrbits = nOccupied + basis%nUnoccupied
    write(*,*) 'Total number of orbits: ', nOrbits

    nxkf = basis%spQuantNr(nOccupied, 1)
    nykf = basis%spQuantNr(nOccupied, 2)
    nzkf = basis%spQuantNr(nOccupied, 3)
    kf2 = factor**2*(nxkf**2 + nykf**2 + nzkf**2)
    kf = SQRT(kf2)

    ALLOCATE(fockVector(nOrbits))
    fockVector = 0.d0
    DO p=1, nOrbits, 2
       npx = basis%spQuantNr(p, 1)
       npy = basis%spQuantNr(p, 2)
       npz = basis%spQuantNr(p, 3)
       
       temp = 0.d0
       DO j=1, nOccupied, 2
          njx = basis%spQuantNr(j, 1)
          njy = basis%spQuantNr(j, 2)
          njz = basis%spQuantNr(j, 3)
          
          v = interaction3d(npx-njx, npy-njy, npz-njz, &
               factor, volume) 
          
          temp = temp - v
       ENDDO
       fockVector(p) = temp
       !write(*,*) 'p=',p,',fockV=',temp
    ENDDO
    
    temp = 0.d0
    DO i=1, nOccupied, 2
       nix = basis%spQuantNr(i, 1)
       niy = basis%spQuantNr(i, 2)
       niz = basis%spQuantNr(i, 3)
       
       ki2 = factor**2*(nix**2 + niy**2 + niz**2)
       ki = SQRT(ki2)
       
       DO j=1, nOccupied, 2
          njx = basis%spQuantNr(j, 1)
          njy = basis%spQuantNr(j, 2)
          njz = basis%spQuantNr(j, 3)
          
          kj2 = factor**2*(njx**2 + njy**2 + njz**2)
          kj = SQRT(kj2)
          
          DO a=nOccupied+1, nOrbits, 2
             nax = basis%spQuantNr(a, 1)
             nay = basis%spQuantNr(a, 2)
             naz = basis%spQuantNr(a, 3)
             
             ka2 = factor**2*(nax**2 + nay**2 + naz**2)
             ka = SQRT(ka2)
             
             DO b=nOccupied+1, nOrbits, 2
                nbx = basis%spQuantNr(b, 1)
                nby = basis%spQuantNr(b, 2)
                nbz = basis%spQuantNr(b, 3)

                kb2 = factor**2*(nbx**2 + nby**2 + nbz**2)
                kb = SQRT(kb2)
                
                !write(*,*) i,j,a,b

                IF ((nix+njx == nax+nbx).AND.&
                     (niy+njy == nay+nby).AND.&
                     (niz+njz == naz+nbz)) THEN
                   via = interaction3d(&
                        nix-nax, niy-nay, niz-naz, factor, volume)
                   vja = interaction3d(&
                        njx-nax, njy-nay, njz-naz, factor, volume)
                   ! denom = 0.5d0*(ki2 + kj2 - ka2 - kb2) &
                   !      + (fockExact(ki/kf) + fockExact(kj/kf) &
                   !      - fockExact(ka/kf) - fockExact(kb/kf))
                   denom = 0.5d0*(ki2 + kj2 - ka2 - kb2) &
                        + (fockVector(i) + fockVector(j) &
                        - fockVector(a) - fockVector(b))
                   
                   temp = temp + (2.d0*via**2 &
                        - via*vja)/denom
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    write(*,*) 'ene_corr=',temp, ' Ha'
    
    DEALLOCATE(fockVector)
    write(*,*) '--'
    
  END SUBROUTINE pt2Sheperd
  
  FUNCTION fockExact(x) RESULT(fock)
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: fock

    fock = 1.d0 + 0.5d0*(1-x**2)*LOG((1+x)/(1-x))/x

  END FUNCTION fockExact
  
  FUNCTION interaction3d(nx, ny, nz, factor, volume) RESULT(v)
    INTEGER, INTENT(IN) :: nx, ny, nz
    REAL(DP), INTENT(IN) :: factor, volume
    REAL(DP) :: q2, v
    
    IF ((nx /= 0).OR.(ny /= 0).OR.(nz /= 0)) THEN
       q2 = factor**2*(nx**2 + ny**2 + nz**2)
       v = 4.d0*pi/(q2*volume)
       !v = 4.d0*pi/q2
    ELSE
       v = 0.d0
    ENDIF
    
    
  END FUNCTION interaction3d
    
  SUBROUTINE testQuicksort2
    REAL(DP), ALLOCATABLE :: rand1(:), rand2(:)
    INTEGER, ALLOCATABLE :: array1(:), array2(:)
    INTEGER :: n, i
    
    n = 100
    ALLOCATE(rand1(n), rand2(n))
    ALLOCATE(array1(n), array2(n))
    
    CALL RANDOM_NUMBER(rand2)
    CALL RANDOM_NUMBER(rand1)
    array1 = INT(10.d0*rand1)
    array2 = INT(10.d0*rand2)
    !array1 = (/8, 2, -1, 4, 3, 6, 5/)
    !array2 = array1
    !array2 = (/(i, i=1, 7)/)
    rand1 = 1.d0*array1
    rand2 = rand1

    DO i=1, n
       write(*,*) array1(i), array2(i)
       !write(*,*) rand1(i), array2(i)
    ENDDO
    WRITE(*,*) '--------------------------'
    
    !CALL quicksort2(array1, array2, n)
    !CALL quicksort(array1, n)
    !CALL sort2(n,rand1,rand2)
!    CALL sortrx(n,rand1,array2)
    CALL qsorti(array2,n,array1)

    DO i=1, n
       write(*,*) array1(array2(i)), array2(i)
       !write(*,*) rand1(i), array2(i)
    ENDDO
    
    DEALLOCATE(array1, array2)
    
    
  END SUBROUTINE testQuicksort2
  
  ! SUBROUTINE testSpPotential(solv)
  !   USE SpPotentialMod
  !   USE SolverMod
  !   TYPE(SpPotential) :: spPot
  !   TYPE(Solver) :: solv
  !   REAL(DP) :: potential,k

  !   k = 1.d0
  !   spPot = solv%energy%eneCorrel%spPot
  !   potential = calcSpPotential(k,spPot)
  !   write(*,*) 'spPot = ', potential
    
  ! END SUBROUTINE testSpPotential

  ! SUBROUTINE testGaussian
  !   USE QuadratureMod
  !   USE HamiltonianMod
  !   TYPE(Quadrature) :: quad
  !   REAL(DP) :: xmin,xmax,rij,potential,k,kp
  !   COMPLEX(DP) :: potentialC
  !   INTEGER :: nx,i

  !   nx = 200
  !   xmin = 0.d0
  !   xmax = 2.d0
  !   !     Set up Gausss-Legendre quadrature
  !   quad = gaussLegendre(xmin,xmax,nx)
  
  !   OPEN(12,FILE='r_potential')

  !   DO i=1, nx
  !      rij = quad%xpoints(i)
  !      potential = gaussian(rij)
  !      write(12,*) rij, potential
  !   ENDDO
  !   CLOSE(12)

  !   k = 2.d0
  !   kp = 0.3d0
  !   potentialC = gaussianRelMom(k,kp)
  !   write(*,*) 'potential = ', potentialC

  !   !     Deallocate quadrature object
  !   CALL deallocateQuad(quad)
    
    
  ! END SUBROUTINE testGaussian

  
  ! SUBROUTINE testHF
  !   USE HamiltonianMod
  !   INTEGER, ALLOCATABLE :: spStates(:)
  !   INTEGER :: nParticles,n,spState,spState1,spState2
  !   REAL(DP) :: boxWidth,density,eneKin,potHF,eneHF
  !   REAL(DP) :: temp,k1,k2,krel,potential,interaction
  !   REAL(DP) :: spPotential
    
  !   density = 1.27323950930405d0
  !   nParticles = 181
  !   IF (MOD(nParticles,2)==0) THEN
  !      WRITE(*,*) 'The number of particles must be odd!'
  !      STOP
  !   ENDIF
    
  !   ALLOCATE(spStates(nParticles))
  !   spStates = 0
  !   n = 1
  !   DO spState=2, nParticles, 2
  !      !write(*,*) 'spState=',spState
  !      spStates(spState) = -n
  !      spStates(spState+1) = n
  !      n = n + 1
  !   ENDDO
    
  !   !write(*,*) 'spStates = ', spStates
    
  !   boxWidth = 1.d0*nParticles/density

  !   temp = 0.d0
  !   DO spState=1, nParticles
  !      temp = temp + 1.d0*spStates(spState)**2
  !   ENDDO
  !   !     Kinetic energy per particle
  !   eneKin = 2.d0*pi**2*temp/(nParticles*mass*boxWidth**2)
    
  !   temp = 0.d0
  !   DO spState1=1, nParticles
  !      !     Lab frame momentum of state 1
  !      k1 = 2.d0*pi*spStates(spState1)/boxWidth
  !      DO spState2=spState1+1, nParticles
  !         !     Lab frame momentum of state 2
  !         k2 = 2.d0*pi*spStates(spState2)/boxWidth
  !         !     Relative momentum
  !         krel = 0.5d0*(k1-k2)
          
  !         potential = REAL(gaussianRelMomDiscr(krel,krel,&
  !              boxWidth),DP)
  !         temp = temp + potential
  !      ENDDO
  !   ENDDO
  !   !      HF potential energy per particle
  !   potHF = temp/nParticles
    
  !   !      HF energy per particle
  !   eneHF = eneKin + potHF
    
  !   write(*,*) 'boxWidth = ', boxWidth
  !   write(*,*) 'eneKin = ', eneKin, ', potHF = ', potHF
  !   write(*,*) 'eneHF = ', eneHF, ', density = ', density
    
    
  !   !     Test the single-particle potential
  !   k1 = 1.d0
    
  !   temp = 0.d0
  !   DO spState=1, nParticles
  !      !     Lab frame momentum of state 1
  !      k2 = 2.d0*pi*spStates(spState)/boxWidth
  !      !     Relative momentum
  !      krel = 0.5d0*(k1-k2)
       
  !      !     Get antisymmetrized two-particle
  !      !     interaction
  !      interaction = REAL(gaussianRelMomDiscr(krel,krel,&
  !           boxWidth),DP)
       
  !      temp = temp + interaction
  !   ENDDO
  !   spPotential = temp
  !   write(*,*) 'spPotential = ', spPotential
    
  !   DEALLOCATE(spStates)
    
    
  ! END SUBROUTINE testHF

  
  
END MODULE TestsMod
