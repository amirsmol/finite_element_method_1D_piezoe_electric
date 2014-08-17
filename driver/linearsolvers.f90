module linearsolvers
   implicit none
      integer,parameter :: ndf=2 !number of degrees of freedom
      integer,parameter :: dimen=1 !dimention of problem

      integer,parameter::inp=1 !input file
      integer,parameter::out=2 !output file
      integer,parameter::msh=3 !mesh file
      integer,parameter::bnd=4 !boundary file
      integer,parameter::tec=5 !the techplot file
      integer,parameter::dat=6 !mesh file itech
      integer,parameter::cur=8 !curve file
      integer,parameter::csv=11 !comma seperated file
      integer,parameter::vtu=13 !output file
      integer,parameter::OUTPUT_UNIT=7

   ! the default value for the smallest pivot that will be accepted
   ! using the linearsolvers subroutines.  pivots smaller than this
   ! threshold will cause premature termination of the linear equation
   ! solver and return false as the return value of the function.

   real*8, parameter :: default_smallest_pivot = 1.0e-6

contains




    FUNCTION SAWTOOTH(TIME)
    IMPLICIT NONE
    REAL*8::TIME,SAWTOOTH
    SAWTOOTH=4*(ABS(TIME+0.25-FLOOR(TIME+0.75))-0.25)
    ENDFUNCTION SAWTOOTH



      subroutine lapack_gesv(a,b)
      implicit none
      integer, parameter :: nmax=1000

! ============================================================================
! name        : lapack solver compact
! author      : amir
! version     :
! copyright   : your copyright notice
! description : this solves a * x=b and put the solution
! into the b vector
! ============================================================================


real*8::a(:,:),b(:)
integer::n, nrhs, lda, ldb, info
integer,allocatable::ipiv(:)

n=size(b)
nrhs=1
lda=n
ldb=n;
allocate(ipiv(n))


call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)

endsubroutine lapack_gesv



    SUBROUTINE CONJGRAD_SOLVER(A,B,X)
    IMPLICIT NONE
    REAL*8,DIMENSION(:),INTENT(IN)::B
    REAL*8,DIMENSION(:)::X
    REAL*8,DIMENSION(:,:),INTENT(IN)::A
    REAL*8,DIMENSION(:),ALLOCATABLE :: P,R,AP
    REAL*8::RSOLD,RSNEW,ALPHA
    INTEGER::DIMEN
! ========================CONVERGENCE VARIABLES
    INTEGER::CG_ITERS,CG_LIMIT
    LOGICAL::CG_CONVERGED
    REAL*8::CG_TOL
! ================================================================

    DIMEN=SIZE(B)
    CG_TOL=1E-5
    CG_LIMIT=500
    CG_ITERS=0

    ALLOCATE(R(DIMEN),P(DIMEN),AP(DIMEN))
    R=B-MATMUL(A,X)
    P=R
    RSOLD=DOT_PRODUCT(R,R)

 PCG: DO; CG_ITERS=CG_ITERS+1;
    AP=MATMUL(A,P)
    ALPHA=RSOLD/DOT_PRODUCT(P,AP)
    X=X+ALPHA*P
    R=R-ALPHA*AP
    RSNEW=DOT_PRODUCT(R,R)
    CG_CONVERGED=(SQRT(RSNEW).LT.CG_TOL)
    IF(CG_CONVERGED.OR.CG_ITERS==CG_LIMIT)EXIT
       P=R+(RSNEW/RSOLD)*P;
       RSOLD=RSNEW;
     END DO PCG
     WRITE(1,*) CG_ITERS
    ENDSUBROUTINE CONJGRAD_SOLVER



    SUBROUTINE CONJGRAD_SOLVER_K(K,F,U)
    IMPLICIT NONE
    REAL*8,DIMENSION(:),INTENT(IN)::F
    REAL*8,DIMENSION(:)::U
    REAL*8,DIMENSION(:,:),INTENT(IN)::K
    REAL*8,DIMENSION(:),ALLOCATABLE :: P,R,Q
    REAL*8::RSOLD,RSNEW,ALPHA,BETA,RTR
    INTEGER::DIMEN
! ========================CONVERGENCE VARIABLES
    INTEGER::CG_ITERS,CG_LIMIT
    LOGICAL::CG_CONVERGED
    REAL*8::CG_TOL
! ================================================================

    DIMEN=SIZE(F)
    CG_TOL=1E-5
    CG_LIMIT=DIMEN*2
    CG_ITERS=0


    ALLOCATE(P(DIMEN),Q(DIMEN),R(DIMEN))
    R=F-MATMUL(K,U)
    P=R

 PCG: DO;
    CG_ITERS=CG_ITERS+1;
!     WRITE(1,*)CG_ITERS
    Q=MATMUL(K,P)
    RTR=DOT_PRODUCT(R,R)
    ALPHA=RTR/DOT_PRODUCT(P,Q)

    U=U+ALPHA*P
    R=R-ALPHA*Q
    BETA=DOT_PRODUCT(R,R)/RTR
    P=R+BETA*P


    CG_CONVERGED=(SQRT(DOT_PRODUCT(R,R)).LT.CG_TOL)
    IF(CG_CONVERGED.OR.CG_ITERS==CG_LIMIT)EXIT
     END DO PCG


     WRITE(1,*)CG_ITERS
    ENDSUBROUTINE CONJGRAD_SOLVER_k


   ! Use Gaussian elimination to calculate the solution to the linear
   ! system, A x = b.  No partial pivoting is done.  If the threshold
   ! argument is present, it is used as the smallest allowable pivot
   ! encountered in the computation; otherwise, DEFAULT_SMALLEST_PIVOT,
   ! defined in this module, is used as the default threshold.  The status
   ! of the computation is a logical returned by the function indicating
   ! the existence of a unique solution (.true.), or the nonexistence of
   ! a unique solution or threshold passed (.false.).

   ! Note that this is an inappropriate method for some linear systems.
   ! In particular, the linear system, M x = b, where M = 10e-12 I, will
   ! cause this routine to fail due to the presence of small pivots.
   ! However, this system is perfectly conditioned, with solution x = b.

   function gaussianElimination( A, b, x, threshold )
      implicit none
      logical gaussianElimination
      real*8, dimension( :, : ), intent( in ) :: A   ! Assume the shape of A.
      real*8, dimension( : ), intent( in ) ::  b     ! Assume the shape of b.
      real*8, dimension( : ), intent( out ) :: x     ! Assume the shape of x.

      ! The optional attribute specifies that the indicated argument
      ! is not required to be present in a call to the function.  The
      ! presence of optional arguments, such as threshold, may be checked
      ! using the intrinsic logical function, present (see below).

      REAL*8, optional, intent( in ) :: threshold

      integer i, j   ! Local index variables.
      integer N      ! Order of the linear system.
      real m         ! Multiplier.
      real :: smallestPivot = DEFAULT_SMALLEST_PIVOT

      ! Pointers to the appropriate rows of the matrix during the elmination.
      REAL*8, dimension( : ), pointer :: pivotRow
      REAL*8, dimension( : ), pointer :: currentRow

      ! Copies of the input arguments.  These copies are modified during
      ! the computation.
      ! The target attribute is used to indicate that the specified
      ! variable may be the target of a pointer.  Rows of ACopy are targets
      ! of pivotRow and currentRow, defined above.

      REAL*8, dimension( size( A, 1 ), size( A, 2) ), target :: ACopy
      REAL*8, dimension( size( b ) ) :: bCopy

      ! Status of the computation.  The return value of the function.
      logical successful
      ! Change the smallestPivot if the threshold argument was included.
      if ( present( threshold ) ) smallestPivot = abs( threshold )

      ! Setup the order of the system by using the intrinsic function size.
      ! size returns the number of elements in the specified dimension of
      ! an array or the total number of elements if the dimension is not
      ! specified.  Also assume that a unique solution exists initially.

      N = size( b )
      ACopy = A
      bCopy = b
      successful = .true.

      ! Begin the Gaussian elimination algorithm.
      ! Note the use of array sections in the following loops.  These
      ! eliminate the need for many do loops that are common in Fortran
      ! 77 code.
      ! Pointers are also used below and enhance the readability of the
      ! elimination process.

      ! Begin with the first row.

      i = 1

      ! Reduce the system to upper triangular.
      do while ( ( successful ) .and. ( i <= N-1 ) )

         ! The following statement is called pointer assignment and uses
         ! the pointer assignment operator `=>'.  This causes pivotRow
         ! to be an alias for the ith row of ACopy.  Note that this does
         ! not cause any movement of data.

         ! Assign the pivot row.
         pivotRow => ACopy( i, : )

         ! Verify that the current pivot is not smaller than smallestPivot.
         successful = abs( pivotRow( i ) ) >= smallestPivot

         if ( successful ) then

        ! Eliminate the entries in the pivot column below the pivot row.

        do j = i+1, N
           ! Assign the current row.
           currentRow => ACopy( j, : )

               ! Calculate the multiplier.
               m = currentRow( i ) / pivotRow( i )

               ! Perform the elimination step on currentRow and right
               ! hand side, bCopy.
               currentRow = m * pivotRow - currentRow
               bCopy( j ) = m * bCopy( i ) - bCopy( j )
        end do

         end if

         ! Move to the next row.
         i = i + 1

      end do

      ! Check the last pivot.
      pivotRow => ACopy( N, : )
      if ( successful ) successful = abs( pivotRow( N ) ) >= smallestPivot

      if ( successful ) then
         do i = N, 2, -1   ! Backward substitution.

            ! Determine the ith unknown, x( i ).
            x( i ) = bCopy( i ) / ACopy( i, i )

            ! Substitute the now known value of x( i ), reducing the order of
            ! the system by 1.
            bCopy = bCopy - x( i ) * ACopy( :, i )

         end do
      end if

      ! Determine the value of x( 1 ) as a special case.
      if ( successful ) x( 1 ) = bCopy( 1 ) / ACopy( 1, 1 )

      ! Prepare the return value of the function.
      gaussianElimination = successful

   end function gaussianElimination


   ! The LU decomposition of a matrix may be represented in a compact form
   ! existing in a single matrix, M,  if the assignments M=L and M=U are
   ! done (in that order).  The diagonal entries in L are assumed to be
   ! unity so that no storage space is necessary.  Instead, the diagonal
   ! of M is used to hold the diagonal entries of U.  This is a common
   ! method of storing the LU decomposition of a matrix.

   ! The algorithm belows makes an additional assumption concerning the
   ! pivots or diagonal elements of U.  Computation terminates if one of
   ! these pivots is smaller than the given or default threshold.  In this
   ! case, the LU decomposition is not formed.  Note that this algorithm
   ! successfully terminates if such an LU can be computed.  In this case
   ! the coefficient matrix, A, is nonsingular.  (No attempt for recovery,
   ! such as permutation of rows, is done.)

   ! Compute the LU decomposition of A, storing the result in LU so that
   ! A is not overwritten.  If the threshold argument is present, it is used
   ! as the smallest allowable pivot encountered in the computation;
   ! otherwise, DEFAULT_SMALLEST_PIVOT, defined in this module, is used as
   ! the default threshold during the computation.  The status of the
   ! computation is a logical returned by the function indicating the
   ! success (.true.) or failure (.false.) of the factorization
   ! After the computation, LU will contain the multipliers below the main
   ! diagonal (L) and the result after elimination on and above the main
   ! diagonal (U), so that A = L * U.

   function LUFactor ( A, LU, threshold )
      implicit none
      logical LUFactor
      REAL*8, dimension( :, : ), intent( in ) :: A
      REAL*8, dimension( :, : ), intent( out ) :: LU
      REAL*8, optional, intent( in ) :: threshold

      integer k, i
      integer N
      logical successful   ! Status of the computation.
      real :: smallestPivot = DEFAULT_SMALLEST_PIVOT

      ! Reassign the smallestPivot, set the order of the system, and
      ! copy A into LU as it will be written to during the factorization.

      if ( present( threshold ) ) smallestPivot = abs( threshold )
      N = size( A, 1 )
      LU = A

      ! Begin the LU factorization algorithm.
      ! The status of the computation is initially successful.
      successful = .true.

      k = 1   ! Begin with the first column.
      do while ( ( successful ) .and. ( k <= N-1 ) )

         ! Verify that the kth pivot is not smaller than smallestPivot.
         successful = abs( LU( k, k ) ) >= smallestPivot

         if ( successful ) then
            ! Calculate the multipliers (L) for the current column.
            LU( k+1:N, k ) = LU( k+1:N, k ) / LU( k, k )

            ! Perform elimination on the upper portion of the matrix (U).
            do i = k+1, N
               LU( i, k+1:N ) = LU( i, k+1:N ) - LU( i, k ) * LU( k, k+1:N )
            enddo

            k = k + 1   ! Move to the next column.
         end if

      enddo

      ! Prepare the return value of the function.
      LUFactor = successful

   end function LUFactor


   ! Let A = L*U where LU represents the LU decomposition of A stored in the
   ! format produced by LUFactor, A, L, U in R**(NxN).
   ! Solve the linear system, A x = b, using the LU decomposition of A stored
   ! in LU.  Since LU is the LU decomposition of A, A is nonsingular.
   ! Consequently, the columns of A constitute a basis for R**N.   So, there
   ! must exist a unique solution to the linear system A x = b.
   ! LUSolve returns the solution to this linear system.

   function LUSolve( LU, b ) result( x )
      implicit none
      REAL*8, dimension( :, : ), intent( in ) :: LU
      REAL*8, dimension( : ), intent( in ) :: b
      REAL*8, dimension( size( b ) ) :: x

      integer k
      integer N
      REAL*8, dimension( size( b ) ) :: bCopy

      ! Determine the order of the system and store a copy of b in bCopy
      ! as it is written during the computation.
      N = size( b )
      bCopy = b

      ! Assume LU is in the form of LU and solve the system in two steps.
      ! First, using forward elmination to solve L y = b, then using
      ! backward elmination to solve U x = y.  In both cases, the right
      ! hand side is overwritten with the solution as it is computed.

      ! Forward elimination.  Store the solution into the right hand side.
      do k = 1, N-1
         bCopy( k+1:N ) = bCopy( k+1:N ) - bCopy( k ) * LU( k+1:N, k )
      end do

      ! Backward elimination.  Store the solution into the right hand side.
      do k = N, 2, -1
         bCopy( k ) = bcopy( k ) / LU( k, k )
         bCopy( 1:k-1 ) = bCopy( 1:k-1 ) - bCopy( k ) * LU( 1:k-1, k )
      end do

      ! Solve for the 1st unknown as a special case.
      bCopy( 1 ) = bCopy( 1 ) / LU( 1, 1 )

      ! Assign a return value for the function via its result variable, x.
      x = bCopy

   end function LUSolve


   ! Output A in Matlab format, using name in the Matlab assignment statement.
   subroutine printMatrix( A, name )
      implicit none
      real*8, dimension( :, : ) :: A   ! Assume the shape of A.
      character name  ! Name for use in assignment, ie, name = ......

      integer n, m, i, j


      OPEN (OUTPUT_UNIT,file='matrixes_out.out')
      n = size( A, 1 )
      m = size( A, 2 )

      write(OUTPUT_UNIT, fmt="(a1,a5)", advance = "no" ) name, ' = [ '

      ! Output the matrix, except for the last row, which needs no `;'.
      do i = 1, n-1

         ! Output current row.
         do j = 1, m-1
            write( OUTPUT_UNIT, fmt="(f10.6,a2)", advance = "no" ) A( i, j ), ', '
         end do

         ! Output last element in row and end current row.
         write( OUTPUT_UNIT, fmt="(f10.6,a1)" ) A( i, m ), ';'

      end do

      ! Output the last row.
      do j = 1, m-1
         write( OUTPUT_UNIT, fmt="(f10.6,a2)", advance = "no" ) A( i, j ), ', '
      end do

      ! Output last element in row and end.
      write( OUTPUT_UNIT, fmt="(f10.6,a1)" ) A( i, m ), ']'

   end subroutine printMatrix


   subroutine printMatrix_int( A, name )
      implicit none
      integer, dimension( :, : ) :: A   ! Assume the shape of A.
      character name  ! Name for use in assignment, ie, name = ......

      integer n, m, i, j

      n = size( A, 1 )
      m = size( A, 2 )
    OPEN (OUTPUT_UNIT,file='matrixes_out.out')
      write( OUTPUT_UNIT, fmt="(a1,a5)", advance = "no" ) name, ' = [ '

      ! Output the matrix, except for the last row, which needs no `;'.
      do i = 1, n-1

         ! Output current row.
         do j = 1, m-1
            write( OUTPUT_UNIT, fmt="(I3,a2)", advance = "no" ) A( i, j ), ', '
         end do

         ! Output last element in row and end current row.
         write( OUTPUT_UNIT, fmt="(I3,a2)" ) A( i, m ), ';'

      end do

      ! Output the last row.
      do j = 1, m-1
         write( OUTPUT_UNIT, fmt="(I3,a2)", advance = "no" ) A( i, j ), ', '
      end do

      ! Output last element in row and end.
      write( 7, fmt="(I3,a2)" ) A( i, m ), ']'

   end subroutine printMatrix_int


   ! Output b in Matlab format, using name in the Matlab assignment statement.
   subroutine printVector_int( b, name )
      implicit none
      integer, dimension( : ) :: b   ! Assume the shape of b.
      character name   ! Name for use in assignment, ie, name = ......

      integer n, i
    OPEN (OUTPUT_UNIT,file='matrixes_out.out')
      n = size( b )

      write( OUTPUT_UNIT, fmt="(a1,a5)", advance = "no" ) name, ' = [ '

      do i = 1, n-1
         write( OUTPUT_UNIT, fmt = "(i3,a2)", advance = "no" ) b( i ), ', '
      end do

      write( OUTPUT_UNIT, fmt = "(i3,a2)" ) b( n ), ']'''

   end subroutine printVector_int

     subroutine printVector( b, name )
      implicit none
      REAL*8, dimension( : ) :: b   ! Assume the shape of b.
      character name   ! Name for use in assignment, ie, name = ......

      integer n, i

      n = size( b )
    OPEN (OUTPUT_UNIT,file='matrixes_out.out')

      write( OUTPUT_UNIT, fmt="(a1,a5)", advance = "no" ) name, ' = [ '

      do i = 1, n-1
         write( OUTPUT_UNIT, fmt = "(E14.5,a2)", advance = "no" ) b( i ), ', '
      end do

      write( OUTPUT_UNIT, fmt = "(E14.5,a2)" ) b( n ), ']'''

   end subroutine printVector



      SUBROUTINE SLVUNSYM(A,NRMAX,NCMAX,N,ITERM)
!C     _________________________________________________________________
!C
!C       Solver for BANDED UNSYMMETRIC system of algebraic equations
!C     ______________________________________________________________
!C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER::NRMAX,NCMAX,N,ITERM
      REAL*8::A(NRMAX,NCMAX)
      CERO=1.0D-15
      PARE=CERO**2
      NBND=2*ITERM
      NBM=NBND-1
!C
!C     Begin elimination of the lower left
!C
      DO 80 I=1,N
      IF (DABS(A(I,ITERM)).LT.CERO) GO TO 10
      GO TO 20
   10 IF (DABS(A(I,ITERM)).LT.PARE) GO TO 110
   20 JLAST=MIN0(I+ITERM-1,N)
      L=ITERM+1
      DO 40 J=I,JLAST
      L=L-1
      IF (DABS(A(J,L)).LT.PARE) GO TO 40
      B=A(J,L)
      DO 30 K=L,NBND
   30 A(J,K)=A(J,K)/B
      IF (I.EQ.N) GO TO 90
   40 CONTINUE
      L=0
      JFIRST=I+1
      IF (JLAST.LE.I) GO TO 80
      DO 70 J=JFIRST,JLAST
      L=L+1
      IF (DABS(A(J,ITERM-L)).LT.PARE) GO TO 70
      DO 50 K=ITERM,NBM
   50 A(J,K-L)=A(J-L,K)-A(J,K-L)
      A(J,NBND)=A(J-L,NBND)-A(J,NBND)
      IF (I.GE.N-ITERM+1) GO TO 70
      DO 60 K=1,L
   60 A(J,NBND-K)=-A(J,NBND-K)
   70 CONTINUE
   80 CONTINUE
   90 L=ITERM-1
      DO 100 I=2,N
      DO 100 J=1,L
      IF (N+1-I+J.GT.N) GO TO 100
      A(N+1-I,NBND)=A(N+1-I,NBND)-A(N+1-I+J,NBND)*A(N+1-I,ITERM+J)
  100 CONTINUE
      RETURN
  110 WRITE (*,140) I,A(I,ITERM)
  140 FORMAT (/,2X,'Computation stopped in SLVUNSYM because zero appears on the main diagonal *** Eqn no. and value:',I5,E12.4)
      ENDSUBROUTINE SLVUNSYM

! -------------------------- MODULE cg.f90 ----------------------------

!************************************************************************
!*                                                                      *
!* Conjugate Gradient Method (CG Method)                                *
!* -------------------------------------                                *
!*                                                                      *
!* Programming language: ANSI C                                         *
!* Compiler:             Turbo C 2.0                                    *
!* Computer:             IBM PS/2 70 with 80387                         *
!* Sources:              [BUNS85], [SCHW], [MAES84]                     *
!* Author:               Juergen Dietel, Computer Center, RWTH Aachen   *
!* Date:                 7.31.1992                                      *
!*                                                                      *
!*             F90 version by J-P Moreau (without dynamic allocations). *
!*                                (www.jpmoreau.fr)                     *
!************************************************************************
Subroutine cg_method (    &     ! Conjugate Gradient Method
                       n, &     ! Size of the linear system
                       a, &     ! System matrix
                       y, &     ! right hand side
                       x, &     ! solution vector
                       fehler & ! error code
                     )
! original name: cg_verfahren()
integer::n
real*8, parameter :: ZERO=0.d0, MACH_EPS=2.d-16
real*8  a(0:n-1,0:n-1),x(0:n-1),y(0:n-1)
integer fehler
!************************************************************************
!* CG solves the linear system                                          *
!*                         A * X = Y                                    *
!* for a symmetric, positive definite matrix A via the conjugate        *
!* gradient method.                                                     *
!*                                                                      *
!* Input parameters:                                                    *
!* =================                                                    *
!* n  Size of the linear system                                         *
!* a  [0..n-1,0..n-1] system matrix A. Only the upper triangle of A is  *
!*    used.                                                             *
!* y  [0..n-1] vector of the right hand side                            *
!*                                                                      *
!* Output parameters:                                                   *
!* ==================                                                   *
!* x  [0..n-1] vector giving the solution                               *
!*                                                                      *
!* Return value:                                                        *
!* =============                                                        *
!* = 0: all is ok                                                       *
!* = 1: n < 2 or other disallowed input parameters                      *
!* = 2: memory exceeded                                                 *
!*                                                                      *
!************************************************************************
  real*8 d(0:n-1), &   ! (0..n-1) auxiliary vectors d and g
         g(0:n-1), &
         AmalD(0:n-1)  ! (0..n-1) auxiliary vector A * d
  real*8 alpha,    &   ! coefficient
         beta,     &   ! coefficient
         dividend, &   ! numerator and denominator of a fraction
         divisor,  &   ! respectively, used to compute alpha, beta
         hilf,     &   ! auxiliary variables
         hilf2,    &
         abstand,  &   ! distance of two successive approximations
                       ! for the solution vector x (taken in the
                       ! euclidean norm)
         xnorm         ! euklidean norm of x
  integer k, i, j      ! loop variables

  if (n < 2) then      ! invalid parameter?
    fehler=1
    return
  end if

  !------------------------------------------------------------------
  ! start with x at the origin
  !------------------------------------------------------------------
  do i = n - 1, 0, -1
    x(i) = ZERO
  end do

  !------------------------------------------------------------------
  ! initialize  d and g :
  ! d = -g = -(a*x - y) = y (since x = 0)
  !------------------------------------------------------------------
  do i = n - 1, 0, -1
    hilf = y(i)
    d(i) = hilf
    g(i) = -hilf
  end do


  !------------------------------------------------------------------
  ! perform at most n steps of the CG Method
  !------------------------------------------------------------------
  do k = n, 0, -1

    !----------------------------------------------------------------
    ! compute new alpha:
    ! alpha = -(d(transp) * g) / (d(transp) * (a * d))
    !----------------------------------------------------------------

    dividend = ZERO
    divisor  = ZERO

    do i = n - 1, 0, -1
      dividend = dividend + d(i) * g(i)
      hilf = ZERO
      do j = 0, i-1
        hilf = hilf + a(j,i) * d(j)
      end do
      do j = i, n-1
        hilf = hilf + a(i,j) * d(j)
      end do
      AmalD(i) = hilf
      divisor = divisor + d(i) * hilf
    end do

    if (divisor.eq.ZERO) then
      fehler=0
      return
    end if

    alpha = -dividend / divisor

    !----------------------------------------------------------------
    ! compute the norm of x und  alpha * d  and find a new x:
    ! x  =  x + alpha * d, then check whether x is close enough,
    ! in order to stop the process before n complete steps
    !----------------------------------------------------------------
    xnorm   = ZERO
    abstand = ZERO

    do i = n - 1, 0, -1
      hilf =  x(i)
      xnorm   = xnorm + hilf*hilf
      hilf2   =  alpha * d(i)
      abstand = abstand + hilf2*hilf2
      x(i)    =  hilf + hilf2
    end do

    if (abstand < MACH_EPS * xnorm) then
      fehler=0
      return
    end if


    !----------------------------------------------------------------
    ! compute new g:   g  =  g + alpha * (a * d)
    !----------------------------------------------------------------
    do i = n - 1, 0, -1
      g(i) = g(i) + alpha * AmalD(i)
    end do

    !----------------------------------------------------------------
    ! compute new beta :
    ! beta = (g(transp) * (a * d)) / (d(transp) * (a * d))
    !----------------------------------------------------------------
    dividend = ZERO
    do i = n - 1, 0, -1
      dividend = dividend + g(i) * AmalD(i)
    end do

    beta = dividend / divisor

    !----------------------------------------------------------------
    ! compute new d :   d  =  - g + beta * d
    !----------------------------------------------------------------
    do i = n - 1, 0, -1
      d(i) = -g(i) + beta * d(i)
    end do

  end do  !k loop

  fehler=0
  return
endSubroutine cg_method

! ---------------------------- END cg.cpp --------------------------

subroutine Gaussian_Elimination_Solver(a,y)
    integer :: k,i
    real*8::a(:,:),y(:)
    integer::n
    n=size(a,dim=1)

    do k = 1, n-1
        a(k+1: n, k) = a(k+1: n, k) / a(k, k)

        a(k+1: n, k+1: n) = a(k+1: n, k+1: n) - &
                matmul(a(k+1: n, k: k), a(k: k, k+1: n))
    end do


    ! L x = f  =>  x = L \ f
    do i = 1, n
        y(i) = y(i) - dot_product(a(i, 1: i-1), y(1: i-1))
    end do

    ! U y = x  =>  y = U \ x
    do i = n, 1, -1
        y(i) = y(i) - dot_product(a(i, i+1: n), y(i+1: n))
        y(i) = y(i) / a(i, i)
    end do

end subroutine Gaussian_Elimination_Solver


      FUNCTION NORM_VECT(A)
      REAL*8::A(:),NORM_VECT
      NORM_VECT=SQRT(DOT_PRODUCT(A,A))
      END FUNCTION NORM_VECT


      SUBROUTINE TENSOR_TRANSFORM_RANK_1(TRANSFORMATION,TENSOR_LOCAL_COORD,TENSOR_GLOBAL_COORD)
!-------------------------------------------------------------------
!This will transform the tensor from local into the global coordinate
!it used the transformation tensor to present the tensor that is
!represented in local coordinate in the local coordinate
!-------------------------------------------------------------------
      IMPLICIT NONE
! ----------------------------------------------------------
      REAL*8,INTENT(IN) :: TRANSFORMATION(3,3) ! TRANSFORMATION TENSOR
      REAL*8,INTENT(IN) :: TENSOR_LOCAL_COORD(3) ! LOCAL COORDINATE TENSOR
      REAL*8,INTENT(OUT) :: TENSOR_GLOBAL_COORD(3) ! GLOBAL COORDINATE TENSOR
      INTEGER :: I,J ! The counter integers
! ----------------------------------------------------------
      TENSOR_GLOBAL_COORD=0.0D0 !INITIALIZING TENSOR
! ----------------------------------------------------------
!   IMPLEMENTING THE SYMMETRY
    DO  I=1,3
    DO  J=1,3
    TENSOR_GLOBAL_COORD(I)=TENSOR_GLOBAL_COORD(I)+TRANSFORMATION(I,J)*TENSOR_LOCAL_COORD(J)
    ENDDO;ENDDO
      ENDSUBROUTINE TENSOR_TRANSFORM_RANK_1
!===================================================================
!
      SUBROUTINE TENSOR_TRANSFORM_RANK_2(TRANSFORMATION,TENSOR_LOCAL_COORD,TENSOR_GLOBAL_COORD)
!-------------------------------------------------------------------
!This will transform the tensor from local into the global coordinate
!it used the transformation tensor to present the tensor that is
!represented in local coordinate in the local coordinate
!-------------------------------------------------------------------
      IMPLICIT NONE
! ----------------------------------------------------------
      REAL*8,INTENT(IN) :: TRANSFORMATION(3,3) ! TRANSFORMATION TENSOR
      REAL*8,INTENT(IN) :: TENSOR_LOCAL_COORD(3,3) ! LOCAL COORDINATE TENSOR
      REAL*8,INTENT(OUT) :: TENSOR_GLOBAL_COORD(3,3) ! GLOBAL COORDINATE TENSOR
      INTEGER :: I,J,K,L! The counter integers
! ----------------------------------------------------------
      TENSOR_GLOBAL_COORD=0.0D0 !INITIALIZING TENSOR
! ----------------------------------------------------------
    DO  I=1,3
    DO  J=1,3
    DO  K=1,3
    DO  L=1,3
    TENSOR_GLOBAL_COORD(I,J)=TENSOR_GLOBAL_COORD(I,J)+TRANSFORMATION(I,K)*TRANSFORMATION(J,L)*TENSOR_LOCAL_COORD(K,L)
    ENDDO;ENDDO;ENDDO;ENDDO

      ENDSUBROUTINE TENSOR_TRANSFORM_RANK_2


!===================================================================
      SUBROUTINE ROTATION(XLOCAL,QTRANSFORM)
!====================================================================
      ! THE PhD THESIS OF AMIR SOHRABI MOLLAYOUSEF
      ! Programmer(s):Amir Sohrabi
      ! Objective: to find the rotation tensor for a given point with respect
      ! to the origin of the problem. This transformation matrix is to be used
      ! in non linear electromechancial analysis of telescopic actoator
!=====================================================================
      IMPLICIT NONE
!=====================================================================
      REAL*8 :: XLOCAL(3) ! THIS IS THE COORDINATE OF THE POINT OF INTEREST
      REAL*8 :: QTRANSFORM(3,3) ! THIS IS THE TRANSFORMTION MATRIX
      REAL*8 :: RADIOUS,EUCLI_NORM ! THIS IS THE TRANSFORMTION MATRIX
      INTEGER::DIMENSIONV,INPUT_UNIT
!=====================================================================
      INPUT_UNIT=2
!      OPEN(INPUT_UNIT,file='INPUT_UNIT.TXT')
!      READ(INPUT_UNIT,*)XLOCAL(:)
      QTRANSFORM=0.0D0
      DIMENSIONV=3
      RADIOUS=SQRT(DOT_PRODUCT(XLOCAL,XLOCAL))

      IF (RADIOUS>0) THEN
      QTRANSFORM(1,1)=XLOCAL(1)/RADIOUS
      QTRANSFORM(1,2)=XLOCAL(2)/RADIOUS
      QTRANSFORM(2,1)=-XLOCAL(2)/RADIOUS
      QTRANSFORM(2,2)=XLOCAL(1)/RADIOUS
      QTRANSFORM(3,3)=1
      ENDIF


!     ================================================================
!                           FORMATS
!     ================================================================
  100 FORMAT(3E14.5)
      ENDSUBROUTINE ROTATION

!     ===============================================================
!                       TRACE OF A 3X3 TENSOR
!     ===============================================================

      FUNCTION TRACE(A)
      REAL*8::A(:,:),TRACE
      INTEGER::I
      TRACE=0.0D0;
      DO I=1,SIZE(A,DIM=1);TRACE=TRACE+A(I,I);ENDDO
      END FUNCTION TRACE

      FUNCTION MATMUL_EPZ(EPZ,ELECT)
      REAL*8::ELECT(3),EPZ(3,3,3),MATMUL_EPZ(3,3)
      INTEGER::I,J,K
      MATMUL_EPZ=0.0D0;
      DO I=1,3;
      DO J=1,3;
      DO K=1,3;
      MATMUL_EPZ(I,J)=MATMUL_EPZ(I,J)+EPZ(K,I,J)*ELECT(K)
      ENDDO;ENDDO;ENDDO
      END FUNCTION MATMUL_EPZ

      FUNCTION MATMUL_D_EPZ(EPZ,STRAN)
      REAL*8::MATMUL_D_EPZ(3),EPZ(3,3,3),STRAN(3,3)
      INTEGER::I,J,K
      MATMUL_D_EPZ=0.0D0;
      DO I=1,3;
      DO J=1,3;
      DO K=1,3;
      MATMUL_D_EPZ(K)=MATMUL_D_EPZ(K)+EPZ(K,I,J)*STRAN(I,J)
      ENDDO;ENDDO;ENDDO
      END FUNCTION MATMUL_D_EPZ

!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************

      SUBROUTINE M33INV (A, AINV,DET, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)


      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET
      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE M33INV

!     ================================================================
!                 NORM OF A VECTOR
!     ================================================================
      FUNCTION NORM_MATRIX(A)
      INTEGER::I,J
      REAL*8::A(:,:),NORM_MATRIX


      NORM_MATRIX=0.0D0
      DO I=1,SIZE(A,DIM=1);DO J=1,SIZE(A,DIM=2);
      NORM_MATRIX=NORM_MATRIX+A(I,J)*A(I,J)
      ENDDO;ENDDO


      END FUNCTION NORM_MATRIX

!     ================================================================
!                 MACAULAY BRACKETS
!     ================================================================
      FUNCTION MACAULAY(A)
      INTEGER::MACAULAY
      REAL*8::A
      MACAULAY=0
      IF (A>=0)THEN
       MACAULAY=1
      ENDIF
      END FUNCTION MACAULAY


!     ================================================================
!                 UNIT VECTOR IN THE DIRECTION OF A VECTOR
!     ================================================================
      FUNCTION UNIT_VECT(A)
      REAL*8::A(:)
      REAL*8,ALLOCATABLE :: UNIT_VECT(:)

      ALLOCATE(UNIT_VECT(SIZE(A)))
      UNIT_VECT=0.0D0

      UNIT_VECT=A/NORM_VECT(A)

      IF(NORM_VECT(A).EQ.0)UNIT_VECT=0.0D0


      END FUNCTION UNIT_VECT



      subroutine elapsedtime (END_STAMP_V,START_STAMP_V)
      integer ( kind = 4 )START_STAMP_V(8),END_STAMP_V(8),values(8),SECONDS,MINUTES,HOURS,MILISECONDS,START_STAMP,END_STAMP,ELAPSED_TIME

      MILISECONDS=mod(END_STAMP_V(8)-START_STAMP_V(8),100)
      SECONDS=mod(END_STAMP_V(7)-START_STAMP_V(7),60)
      MINUTES=mod(END_STAMP_V(6)-START_STAMP_V(6),60)
      HOURS=END_STAMP_V(5)-START_STAMP_V(5)
      write ( *,'(1X,a,2x,i2,a1,i4.4,a1,i4.4,a1,i4.4)' )'Elapsed Time=', HOURS,':',MINUTES,':',SECONDS,'.',ABS(MILISECONDS)
      endsubroutine elapsedtime


subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

    end module linearsolvers
