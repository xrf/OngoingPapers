
!
!     Class for one-dimensional mesh
!
MODULE OneDimMeshMod
  USE ConstantsMod
  IMPLICIT NONE

  TYPE OneDimMesh
     SEQUENCE
     INTEGER :: nMeshp                            ! Number of mesh points
     REAL(DP), ALLOCATABLE :: points(:)           ! Mesh points
     REAL(DP), ALLOCATABLE :: values(:)           ! Vaules at the given
                                                  !   mesh points
  END TYPE OneDimMesh
  
CONTAINS
  
  !
  !     Constructor for OneDimMesh object
  !
  FUNCTION OneDimMesh_(nMeshp) RESULT(mesh)
    INTEGER, INTENT(IN) :: nMeshp
    TYPE(OneDimMesh) :: mesh
    
    mesh%nMeshp = nMeshp
    ALLOCATE(mesh%points(nMeshp))
    ALLOCATE(mesh%values(nMeshp))
    mesh%points = 0.d0
    mesh%values = 0.d0
    
    
  END FUNCTION OneDimMesh_
  
  !
  !     Initialize OneDimMesh object
  !
  FUNCTION initOneDimMesh(points,values,nPoints) RESULT(mesh)
    INTEGER, INTENT(IN) :: nPoints
    REAL(DP), INTENT(IN) :: points(nPoints)
    REAL(DP), INTENT(IN) :: values(nPoints)
    TYPE(OneDimMesh) :: mesh
    
    mesh%nMeshp = nPoints
    mesh%points = points
    mesh%values = values
  
    
  END FUNCTION initOneDimMesh
  
  !
  !     Interface for one-dimensional Lagrange interpolation
  !
  SUBROUTINE lagrInterpol1dim(mesh,evaluationPoint,interpValue)
    TYPE(OneDimMesh), INTENT(IN) :: mesh
    REAL(DP), INTENT(IN) :: evaluationPoint
    REAL(DP), INTENT(OUT) :: interpValue
    REAL(DP) :: error
    
    CALL polint_1dim(mesh%points,mesh%values,mesh%nMeshp,&
         evaluationPoint,interpValue,error)
    
    
  END SUBROUTINE lagrInterpol1dim
  
  !
  !     Lagrange interpolation using the Neville's algorithm. This 
  !     algorithm is implemented as suggested in Numerical Recipes.
  !     Given arrays xa and ya, each of length n, and given a value
  !     x, the routine returns a value y, and an error estimate dy.
  !
  !     This version does a four point interpolation, i.e. inter-
  !     polation with a polynomial of third degree. When the vector
  !     dimension is large, the low order polynomial interpolation 
  !     of polint_1dim() is much more reliable than a higher order 
  !     interpolation, such as that obtained from lagrange_1dim().
  !
  !     The routine is a modification of an example in Numerical
  !     Recipes. Modifications and implementation in Fortran 95 
  !     by Gustav Baardsen in March 2011. 
  !
  !     July 26, 2011: Modified to handle interpolation close to 
  !     the end points correctly, using linear interpolation there. 
  !
  !     Tested and seems to work.
  !
  SUBROUTINE polint_1dim(xa,ya,n,x,y,dy)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: xa(n),ya(n),x
    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(OUT) :: y,dy
    REAL(DP) :: c(4),d(4),x0(4)
    REAL(DP) :: ho,hp,w,den
    INTEGER :: ns,i,m,j_lower,j_upper,j_mid
    
    IF (n < 4) THEN
       write(*,*) ' '
       write(*,*) 'WARNING! Dimension too small for interpolation'
       write(*,*) ' '
       STOP
    ENDIF
    IF ((x < xa(1)).OR.(x > xa(n))) THEN
       write(*,*) ' '
       write(*,*) 'WARNING! The point is outside the allowed interpolation domain'
       write(*,*) ' '
       STOP
    ENDIF

    !     Do a bisection search to find the index j such that
    !     xa(j) <= x <= xa(j+1).
    j_lower = 0
    j_upper = n+1
    DO WHILE (j_upper-j_lower > 1)
       j_mid = (j_upper+j_lower)/2
       IF ((xa(n)>xa(1)).EQV.(x>xa(j_mid))) THEN
          j_lower = j_mid
       ELSE
          j_upper = j_mid
       ENDIF
    ENDDO
    ns = j_lower
    !     Four point interpolation, with two points smaller 
    !     than x and two points larger than x.
    IF ((ns >= 2).AND.(ns+2 <= n)) THEN
       
       DO i=1, 4
          c(i) = ya(ns-2+i)
          d(i) = ya(ns-2+i)
          x0(i) = xa(ns-2+i)
       ENDDO
       
       y = ya(ns)        ! This is the initial approximation to Y.
       ns = 1
       DO m=1, 3         ! For each column of the tableau,
          DO i=1, 4-m    ! we loop over the current c's and d's 
             ! and update them. 
             ho = x0(i) - x
             hp = x0(i+m) - x
             w = c(i+1) - d(i)
             den = ho - hp         
             !     This error can occur only if two input xa's are 
             !     to within roundoff identical.
             IF (den.EQ.0.) THEN
                write(*,*) '     OBS! Two input x values are identical in'
                write(*,*) '     interpolation routine.'
                PAUSE
             ENDIF
             den = w/den
             d(i) = hp*den         ! Here the c's and d's are updated.
             c(i) = ho*den
          ENDDO
          IF (2*ns.LT.4-m) THEN
             dy = c(ns+1)
          ELSE
             dy = d(ns)
             ns = ns-1
          ENDIF
          y = y+dy
       ENDDO
       
       !     Only first order interpolation at the borders
    ELSEIF (ns < 2) THEN
       y = ya(ns) + (x-xa(ns))*(ya(ns+1)-ya(ns))/(xa(ns+1)-xa(ns))
       dy = ya(ns+1) - ya(ns)
       
    ELSEIF (ns+2 > n) THEN
       if (ns+1 > n) write(*,*) 'ns + 1 = ', ns, ', n = ', n
       y = ya(ns) + (x-xa(ns))*(ya(ns+1)-ya(ns))/(xa(ns+1)-xa(ns))
       dy = ya(ns+1) - ya(ns)
       
    ENDIF
    
    
  END SUBROUTINE polint_1dim

  
  
  
END MODULE OneDimMeshMod
