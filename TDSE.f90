!Modern Computational Techniques Midterm Coursework - Part 2
!URN: 6457454

!This program aims to solve the Time Independant Schrondinger Equation in 1D

PROGRAM TDSE
    IMPLICIT NONE
 
    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: psi
    COMPLEX, DIMENSION(:), ALLOCATABLE ::  A, B, C, D 
    REAL :: dx, dt, xrange, tperiod, x, t, expv
    INTEGER, PARAMETER :: xmax = 50 !Number of space steps divided by 2.
    INTEGER, PARAMETER :: tmax = 25 !Number of time steps.
    REAL, PARAMETER :: pi = 4 * ATAN(1.0)
    COMPLEX, PARAMETER :: i = CMPLX(0.0,1.0)
    INTEGER :: j, k, info

    xrange = 5
    tperiod = 10
  
    dx = xrange/xmax !Step length.
    dt = tperiod/tmax !Time length per step.
  
    !For tridiagonal matrix solution, we lose elements of A and C diagonals.
    ALLOCATE(A(-xmax+1:xmax))
    ALLOCATE(B(-xmax:xmax))
    ALLOCATE(C(-xmax+1:xmax))
    ALLOCATE(D(-xmax:xmax))
  
    !Generate an x by t array.
    ALLOCATE(psi(-xmax:xmax,0:tmax))

    !Open files for output.
    OPEN(10, file='wavefunction.dat')
    OPEN(11, file='expectation.dat')
  
    !Initial conditions: DO loop to generate list of values psi at time t = 0 for all x in dx intervals.
    DO j = -xmax, xmax
        x = dx*j
        psi(j,0) = ((2/pi)**1/4)*EXP(-(x-1.0)**2.0) !Initial condition provided multiplied by analytical determined N.
    END DO 

    !DO for all timesteps.
    DO j = 1, tmax 
        !Define constant coefficients for A and C array for each iteration as LAPACK destroys input arguments 
        A = dt 
        C = dt 
        !DO for all values of x (within edge bounds of array).
        DO k = -xmax+1, xmax-1
            x = dx*k 
            !Generate the coefficient for B for each iteration as LAPACK destroys input arguments.
            B(k) = 4.0*i*(dx**2.0)-2.0*dt-(dt*(dx**2)*((x**2))) !check if you can take out set of brackets here.
            !D is a function of B and psi values; defined for each iteration as LAPACK detroys input arguments.
            !Crank-Nicholson Method used.
            D(k) = CMPLX(-REAL(B(k)), AIMAG(B(K)))*psi(k,j-1) + psi(k+1,j-1)*(-dt) + psi(k-1,j-1)*(-dt) 
        END DO 

        !LAPACK library to solve tridiagonal matrix. (xmax*2)+1 is used as argument to account for negative and positive sides of symmetrical space x, adding 1 to account for zero.
        CALL CGTSV((xmax*2)+1, 1, A, B, C, D, (xmax*2)+1, info)

        !Variable INFO outputs non-zero if subroutine failed to execute.
        IF (info /= 0) THEN
            WRITE(6,*) 'LAPACK routine CGSTV failed to execute.'
            EXIT 
        END IF 

        !Write output to D to be reused in the next iteration.
        psi(:,j) = D 

    END DO 

    !Write results to file 'wavefunction.dat'.
    DO j = 0, tmax
        t = dt*j 
        DO k = -xmax,xmax 
            x = dx*k
            WRITE(10,*) x, t, REAL(psi(k,j))
        END DO
        WRITE(10,*)
    END DO
    CLOSE(10)

    !Generate expectation value across all x - using Riemann approximation.
    DO j = 0, tmax 
        expv = 0
        t = dt*j 
        DO k = -xmax,xmax 
            x = dx*k  
            expv = expv + x*dx*ABS(psi(k,j)*CONJG(psi(k,j)))
        END DO
        WRITE(11,*) t, expv
    END DO
    WRITE(11,*) 
    CLOSE(11)

END PROGRAM TDSE

