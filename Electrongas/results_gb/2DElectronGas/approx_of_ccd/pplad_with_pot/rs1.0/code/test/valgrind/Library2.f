C From HDK@psuvm.psu.edu Thu Dec  8 15:27:16 MST 1994
C 
C The following was converted from Algol recursive to Fortran iterative
C by a colleague at Penn State (a long time ago - Fortran 66, please
C excuse the GoTo's). The following code also corrects a bug in the
C Quicksort algorithm published in the ACM (see Algorithm 402, CACM,
C Sept. 1970, pp 563-567; also you younger folks who weren't born at
C that time might find interesting the history of the Quicksort
C algorithm beginning with the original published in CACM, July 1961,
C pp 321-322, Algorithm 64). Note that the following algorithm sorts
C integer data; actual data is not moved but sort is affected by sorting
C a companion index array (see leading comments). The data type being
C sorted can be changed by changing one line; see comments after
C declarations and subsequent one regarding comparisons(Fortran
C 77 takes care of character comparisons of course, so that comment
C is merely historical from the days when we had to write character
C compare subprograms, usually in assembler language for a specific
C mainframe platform at that time). But the following algorithm is
C good, still one of the best available.


      SUBROUTINE QSORTI (ORD,N,A)
C
C==============SORTS THE ARRAY A(I),I=1,2,...,N BY PUTTING THE
C   ASCENDING ORDER VECTOR IN ORD.  THAT IS ASCENDING ORDERED A
C   IS A(ORD(I)),I=1,2,...,N; DESCENDING ORDER A IS A(ORD(N-I+1)),
C   I=1,2,...,N .  THIS SORT RUNS IN TIME PROPORTIONAL TO N LOG N .
C
C
C     ACM QUICKSORT - ALGORITHM #402 - IMPLEMENTED IN FORTRAN 66 BY
C                                 WILLIAM H. VERITY, WHV@PSUVM.PSU.EDU
C                                 CENTER FOR ACADEMIC COMPUTING
C                                 THE PENNSYLVANIA STATE UNIVERSITY
C                                 UNIVERSITY PARK, PA.  16802
C
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION ORD(N),POPLST(2,20)
      INTEGER X,XX,Z,ZZ,Y
C
C     TO SORT DIFFERENT INPUT TYPES, CHANGE THE FOLLOWING
C     SPECIFICATION STATEMENTS; FOR EXAMPLE, FOR FORTRAN CHARACTER
C     USE THE FOLLOWING:  CHARACTER *(*) A(N)
C
      INTEGER A(N)
C
      NDEEP=0
      U1=N
      L1=1
      DO 1  I=1,N
    1 ORD(I)=I
    2 IF (U1.LE.L1) RETURN
C
    3 L=L1
      U=U1
C
C PART
C
    4 P=L
      Q=U
C     FOR CHARACTER SORTS, THE FOLLOWING 3 STATEMENTS WOULD BECOME
C     X = ORD(P)
C     Z = ORD(Q)
C     IF (A(X) .LE. A(Z)) GO TO 2
C
C     WHERE "CLE" IS A LOGICAL FUNCTION WHICH RETURNS "TRUE" IF THE
C     FIRST ARGUMENT IS LESS THAN OR EQUAL TO THE SECOND, BASED ON "LEN"
C     CHARACTERS.
C
      X=A(ORD(P))
      Z=A(ORD(Q))
      IF (X.LE.Z) GO TO 5
      Y=X
      X=Z
      Z=Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
    5 IF (U-L.LE.1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q
C
C LEFT
C
    6 P=P+1
      IF (P.GE.Q) GO TO 7
      X=A(ORD(P))
      IF (X.GE.XX) GO TO 8
      GO TO 6
    7 P=Q-1
      GO TO 13
C
C RIGHT
C
    8 Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=A(ORD(Q))
      IF (Z.LE.ZZ) GO TO 10
      GO TO 8
    9 Q=P
      P=P-1
      Z=X
      X=A(ORD(P))
C
C DIST
C
   10 IF (X.LE.Z) GO TO 11
      Y=X
      X=Z
      Z=Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
   11 IF (X.LE.XX) GO TO 12
      XX=X
      IX=P
   12 IF (Z.GE.ZZ) GO TO 6
      ZZ=Z
      IZ=Q
      GO TO 6
C
C OUT
C
   13 CONTINUE
      IF (.NOT.(P.NE.IX.AND.X.NE.XX)) GO TO 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
   14 CONTINUE
      IF (.NOT.(Q.NE.IZ.AND.Z.NE.ZZ)) GO TO 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
   15 CONTINUE
      IF (U-Q.LE.P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
   16 U1=U
      L1=Q+1
      U=P-1
   17 CONTINUE
      IF (U1.LE.L1) GO TO 18
C
C START RECURSIVE CALL
C
      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
   18 IF (U.GT.L) GO TO 4
C
C POP BACK UP IN THE RECURSION LIST
C
      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18
C
C END SORT
C END QSORT
C
      END



c$$$
c$$$C From Leonard J. Moss of SLAC:
c$$$
c$$$C Here's a hybrid QuickSort I wrote a number of years ago.  It's
c$$$C based on suggestions in Knuth, Volume 3, and performs much better
c$$$C than a pure QuickSort on short or partially ordered input arrays.  
c$$$
c$$$      SUBROUTINE SORTRX(N,DATA,INDEX)
c$$$C===================================================================
c$$$C
c$$$C     SORTRX -- SORT, Real input, indeX output
c$$$C
c$$$C
c$$$C     Input:  N     INTEGER
c$$$C             DATA  REAL
c$$$C
c$$$C     Output: INDEX INTEGER (DIMENSION N)
c$$$C
c$$$C This routine performs an in-memory sort of the first N elements of
c$$$C array DATA, returning into array INDEX the indices of elements of
c$$$C DATA arranged in ascending order.  Thus,
c$$$C
c$$$C    DATA(INDEX(1)) will be the smallest number in array DATA;
c$$$C    DATA(INDEX(N)) will be the largest number in DATA.
c$$$C
c$$$C The original data is not physically rearranged.  The original order
c$$$C of equal input values is not necessarily preserved.
c$$$C
c$$$C===================================================================
c$$$C
c$$$C SORTRX uses a hybrid QuickSort algorithm, based on several
c$$$C suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
c$$$C "pivot key" [my term] for dividing each subsequence is chosen to be
c$$$C the median of the first, last, and middle values of the subsequence;
c$$$C and the QuickSort is cut off when a subsequence has 9 or fewer
c$$$C elements, and a straight insertion sort of the entire array is done
c$$$C at the end.  The result is comparable to a pure insertion sort for
c$$$C very short arrays, and very fast for very large arrays (of order 12
c$$$C micro-sec/element on the 3081K for arrays of 10K elements).  It is
c$$$C also not subject to the poor performance of the pure QuickSort on
c$$$C partially ordered data.
c$$$C
c$$$C Created:  15 Jul 1986  Len Moss
c$$$C
c$$$C===================================================================
c$$$ 
c$$$      INTEGER   N,INDEX(N)
c$$$      REAL      DATA(N)
c$$$ 
c$$$      INTEGER   LSTK(31),RSTK(31),ISTK
c$$$      INTEGER   L,R,I,J,P,INDEXP,INDEXT
c$$$      REAL      DATAP
c$$$ 
c$$$C     QuickSort Cutoff
c$$$C
c$$$C     Quit QuickSort-ing when a subsequence contains M or fewer
c$$$C     elements and finish off at end with straight insertion sort.
c$$$C     According to Knuth, V.3, the optimum value of M is around 9.
c$$$ 
c$$$      INTEGER   M
c$$$      PARAMETER (M=9)
c$$$ 
c$$$C===================================================================
c$$$C
c$$$C     Make initial guess for INDEX
c$$$      
c$$$      DO 50 I=1,N
c$$$         INDEX(I)=I
c$$$ 50   CONTINUE
c$$$ 
c$$$C     If array is short, skip QuickSort and go directly to
c$$$C     the straight insertion sort.
c$$$ 
c$$$      IF (N.LE.M) GOTO 900
c$$$ 
c$$$C===================================================================
c$$$C
c$$$C     QuickSort
c$$$C
c$$$C     The "Qn:"s correspond roughly to steps in Algorithm Q,
c$$$C     Knuth, V.3, PP.116-117, modified to select the median
c$$$C     of the first, last, and middle elements as the "pivot
c$$$C     key" (in Knuth's notation, "K").  Also modified to leave
c$$$C     data in place and produce an INDEX array.  To simplify
c$$$C     comments, let DATA[I]=DATA(INDEX(I)).
c$$$ 
c$$$C Q1: Initialize
c$$$      ISTK=0
c$$$      L=1
c$$$      R=N
c$$$ 
c$$$ 200  CONTINUE
c$$$      
c$$$C Q2: Sort the subsequence DATA[L]..DATA[R].
c$$$C
c$$$C     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
c$$$C     r > R, and L <= m <= R.  (First time through, there is no
c$$$C     DATA for l < L or r > R.)
c$$$ 
c$$$      I=L
c$$$      J=R
c$$$ 
c$$$C Q2.5: Select pivot key
c$$$C
c$$$C     Let the pivot, P, be the midpoint of this subsequence,
c$$$C     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
c$$$C     so the corresponding DATA values are in increasing order.
c$$$C     The pivot key, DATAP, is then DATA[P].
c$$$ 
c$$$      P=(L+R)/2
c$$$      INDEXP=INDEX(P)
c$$$      DATAP=DATA(INDEXP)
c$$$      
c$$$      IF (DATA(INDEX(L)) .GT. DATAP) THEN
c$$$         INDEX(P)=INDEX(L)
c$$$         INDEX(L)=INDEXP
c$$$         INDEXP=INDEX(P)
c$$$         DATAP=DATA(INDEXP)
c$$$      ENDIF
c$$$      
c$$$      IF (DATAP .GT. DATA(INDEX(R))) THEN
c$$$         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
c$$$            INDEX(P)=INDEX(L)
c$$$            INDEX(L)=INDEX(R)
c$$$         ELSE
c$$$            INDEX(P)=INDEX(R)
c$$$         ENDIF
c$$$         INDEX(R)=INDEXP
c$$$         INDEXP=INDEX(P)
c$$$         DATAP=DATA(INDEXP)
c$$$      ENDIF
c$$$      
c$$$C     Now we swap values between the right and left sides and/or
c$$$C     move DATAP until all smaller values are on the left and all
c$$$C     larger values are on the right.  Neither the left or right
c$$$C     side will be internally ordered yet; however, DATAP will be
c$$$C     in its final position.
c$$$ 
c$$$ 300  CONTINUE
c$$$ 
c$$$C Q3: Search for datum on left >= DATAP
c$$$C
c$$$C     At this point, DATA[L] <= DATAP.  We can therefore start scanning
c$$$C     up from L, looking for a value >= DATAP (this scan is guaranteed
c$$$C     to terminate since we initially placed DATAP near the middle of
c$$$C     the subsequence).
c$$$ 
c$$$      I=I+1
c$$$      IF (DATA(INDEX(I)).LT.DATAP) GOTO 300
c$$$         
c$$$ 400  CONTINUE
c$$$      
c$$$C Q4: Search for datum on right <= DATAP
c$$$C
c$$$C     At this point, DATA[R] >= DATAP.  We can therefore start scanning
c$$$C     down from R, looking for a value <= DATAP (this scan is guaranteed
c$$$C     to terminate since we initially placed DATAP near the middle of
c$$$C     the subsequence).
c$$$ 
c$$$      J=J-1
c$$$      IF (DATA(INDEX(J)).GT.DATAP) GOTO 400
c$$$ 
c$$$C     Q5: Have the two scans collided?
c$$$ 
c$$$      IF (I.LT.J) THEN
c$$$ 
c$$$C Q6: No, interchange DATA[I] <--> DATA[J] and continue
c$$$ 
c$$$         INDEXT=INDEX(I)
c$$$         INDEX(I)=INDEX(J)
c$$$         INDEX(J)=INDEXT
c$$$         GOTO 300
c$$$      ELSE
c$$$         
c$$$C Q7: Yes, select next subsequence to sort
c$$$C
c$$$C     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
c$$$C     for all L <= l < I and J < r <= R.  If both subsequences are
c$$$C     more than M elements long, push the longer one on the stack and
c$$$C     go back to QuickSort the shorter; if only one is more than M
c$$$C     elements long, go back and QuickSort it; otherwise, pop a
c$$$C     subsequence off the stack and QuickSort it.
c$$$ 
c$$$         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
c$$$            ISTK=ISTK+1
c$$$            LSTK(ISTK)=J+1
c$$$            RSTK(ISTK)=R
c$$$            R=I-1
c$$$         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
c$$$            ISTK=ISTK+1
c$$$            LSTK(ISTK)=L
c$$$            RSTK(ISTK)=I-1
c$$$            L=J+1
c$$$         ELSE IF (R-J .GT. M) THEN
c$$$            L=J+1
c$$$         ELSE IF (I-L .GT. M) THEN
c$$$            R=I-1
c$$$         ELSE
c$$$C Q8: Pop the stack, or terminate QuickSort if empty
c$$$            IF (ISTK.LT.1) GOTO 900
c$$$            L=LSTK(ISTK)
c$$$            R=RSTK(ISTK)
c$$$            ISTK=ISTK-1
c$$$         ENDIF
c$$$         GOTO 200
c$$$      ENDIF
c$$$      
c$$$ 900  CONTINUE
c$$$      
c$$$C===================================================================
c$$$C
c$$$C Q9: Straight Insertion sort
c$$$      
c$$$      DO 950 I=2,N
c$$$         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
c$$$            INDEXP=INDEX(I)
c$$$            DATAP=DATA(INDEXP)
c$$$            P=I-1
c$$$ 920        CONTINUE
c$$$            INDEX(P+1) = INDEX(P)
c$$$            P=P-1
c$$$            IF (P.GT.0) THEN
c$$$               IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
c$$$            ENDIF
c$$$            INDEX(P+1) = INDEXP
c$$$         ENDIF
c$$$ 950  CONTINUE
c$$$ 
c$$$C===================================================================
c$$$C
c$$$C     All done
c$$$      
c$$$      END
c$$$      
c$$$      
c$$$
c$$$!
c$$$!     Quicksort routine from W. H. Press et al.,
c$$$!     Numerical Recipes in Fortran 77, second ed.,
c$$$!     1992.
c$$$!
c$$$      SUBROUTINE sort2(n,arr,brr)
c$$$      INTEGER n,M,NSTACK
c$$$      REAL arr(n),brr(n)
c$$$      PARAMETER (M=1,NSTACK=1000)
c$$$!     Sorts an array arr(1:n) into ascending order using 
c$$$!     Quicksort, while making the corresponding rearrangement
c$$$!     of the array brr(1:n).
c$$$      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
c$$$      REAL a,b,temp
c$$$
c$$$!      write(*, 950) arr
c$$$ 950  format (7(2X, I3))
c$$$      jstack=0
c$$$      l=1
c$$$      ir=n
c$$$!     Insertion sort when subarray small enough.
c$$$ 1    if(ir-l.lt.M)then
c$$$         do 12 j=l+1,ir
c$$$            a=arr(j)
c$$$            b=brr(j)
c$$$            do 11 i=j-1,l,-1
c$$$               if(arr(i).le.a)goto 2
c$$$               arr(i+1)=arr(i)
c$$$               brr(i+1)=brr(i)
c$$$ 11         enddo 
c$$$            i=l-1
c$$$ 2          arr(i+1)=a
c$$$            brr(i+1)=b
c$$$ 12      enddo 
c$$$         if(jstack.eq.0)return
c$$$!     Pop stack and begin a new round of partitioning.
c$$$         ir=istack(jstack)
c$$$         l=istack(jstack-1)
c$$$         jstack=jstack-2
c$$$      else
c$$$!     Coose median of left, center and right elements as
c$$$!     partitioning element a. Also rearrange so that 
c$$$!     a(1) <= a(l+1) <= a(ir).
c$$$         k=(l+ir)/2
c$$$         temp=arr(k)
c$$$         arr(k)=arr(l+1)
c$$$         arr(l+1)=temp
c$$$         temp=brr(k)
c$$$         brr(k)=brr(l+1)
c$$$         brr(l+1)=temp
c$$$         if(arr(l).gt.arr(ir))then
c$$$            temp=arr(l)
c$$$            arr(l)=arr(ir)
c$$$            arr(ir)=temp
c$$$            temp=brr(l)
c$$$            brr(l)=brr(ir)
c$$$            brr(ir)=temp
c$$$         endif
c$$$         if(arr(l+1).gt.arr(ir))then
c$$$            temp=arr(l+1)
c$$$            arr(l+1)=arr(ir)
c$$$            arr(ir)=temp
c$$$            temp=brr(l+1)
c$$$            brr(l+1)=brr(ir)
c$$$            brr(ir)=temp
c$$$         endif
c$$$         if(arr(l).gt.arr(l+1))then
c$$$            temp=arr(l)
c$$$            arr(l)=arr(l+1)
c$$$            arr(l+1)=temp
c$$$            temp=brr(l)
c$$$            brr(l)=brr(l+1)
c$$$            brr(l+1)=temp
c$$$         endif
c$$$!     Initialize pointers for partitioning.
c$$$         i=l+1
c$$$         j=ir
c$$$!     Partitioning element.
c$$$         a=arr(l+1)
c$$$         b=brr(l+1)
c$$$!     Beginning of innermost loop.
c$$$ 3       continue
c$$$!     Scan up to find element > a.
c$$$         i=i+1
c$$$         if(arr(i).lt.a)goto 3
c$$$ 4       continue
c$$$!     Scan down to find element < a.
c$$$         j=j-1
c$$$         if(arr(j).gt.a)goto 4
c$$$!     Pointers crossed. Exit with partitioning 
c$$$!     complete. Exchange elements of both arrrays.
c$$$         if(j.lt.i)goto 5
c$$$         temp=arr(i)
c$$$         arr(i)=arr(j)
c$$$         arr(j)=temp
c$$$         temp=brr(i)
c$$$         brr(i)=brr(j)
c$$$         brr(j)=temp
c$$$!     End of innermost loop.
c$$$         goto 3
c$$$!     Insert partitioning element in both arrays.
c$$$ 5       arr(l+1)=arr(j)
c$$$         arr(j)=a
c$$$         brr(l+1)=brr(j)
c$$$         brr(j)=b
c$$$         jstack=jstack+2
c$$$!     Push pointers to larger subarray on stack, 
c$$$!     process smaller subarray immediately.
c$$$         if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
c$$$         if(ir-i+1.ge.j-l)then
c$$$            istack(jstack)=ir
c$$$            istack(jstack-1)=i
c$$$            ir=j-1
c$$$         else
c$$$            istack(jstack)=j-1
c$$$            istack(jstack-1)=l
c$$$            l=i
c$$$         endif
c$$$         write(*, 950) arr
c$$$      endif
c$$$      goto 1  
c$$$      END 
