!This program will read any single column file named 'signal.dat' of length 2**n and compute the Fast Fourier
!Transform of that data provided. Parallelisation via OpenMP is implemented within the iteration process in
!order to speed up the matrix multiplication of the recursion relation. The notes provided contain more
!information about how I parallelised. Once the FFT is finished, the program then uses the FFTW library
!to perform the inverse FFT on the parallelised FFT output in order to check if the program Fourier
!transformed the data correctly. If so, the output of the inverse FFT will be equal to the original input data.

PROGRAM fft

    !These lines tell the compiler to: use the module omp_lib which define the OpenMP functions, use the
    !intrinsic module for C-Fortran in order to use the FFTW library.
    USE omp_lib
    USE, INTRINSIC :: iso_c_binding

    IMPLICIT NONE

    !This line includes 'fftw3.f03', a file which contains declarations, function interfaces and subroutine
    !interfaces for FFTW.
    INCLUDE 'fftw3.f03'

    TYPE(C_PTR) :: plan
    REAL(C_DOUBLE), DIMENSION(:), ALLOCATABLE :: fftw_out
    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), ALLOCATABLE :: fftw_in
    INTEGER :: me, nprocs, NperCPU, lbeg, lend
    INTEGER :: logN, N, i, m, kplus, kmax, j, k, l
    REAL :: pi, x, fsum, test
    REAL, DIMENSION(:,:), ALLOCATABLE :: a, b, an, bn
    REAL, DIMENSION(:), ALLOCATABLE :: f


!---READ IN DATA-----------------------------------------------------------------------------------------------------------------

    !This will count the number of lines in the data file signal.dat and set N equal to the number of data points.
    !REWIND is used to reset the data file so that it can be read into the f array without error.
    OPEN(42,FILE="signal.dat")
    N = 0
    DO
        READ(42,*,END=7)
        N = N + 1
    END DO
7   WRITE(6,*) '       File: signal.dat successfully read.'
    WRITE(6,*) N, "lines of data in this file."
    REWIND(42)

    !Ends program if data file does not meet 2^n lines od data rule.
    test = log(real(N))/log(2.0)
    IF(mod(test*2,2.0)/=0) THEN
        WRITE(6,*) "The file provided does not have N = 2^n data points, program will now finish."
        STOP
    END IF

    !Sets up f array with data point provided from signal.dat
    ALLOCATE(f(0:N-1))
    DO i=0,N-1
        READ(42,*) f(i)
    END DO
    CLOSE(42)


!---SET UP CONDITIONS-----------------------------------------------------------------------------------------------------------------

    !Arrays to hold the values for the a and b coefficients. Two arrays are needed such that a holds the running
    !value and an holds the new value through each iteration, likewise for b and bn.
    ALLOCATE(a(0:N-1,0:N-1))
    ALLOCATE(b(0:N-1,0:N-1))
    ALLOCATE(an(0:N-1,0:N-1))
    ALLOCATE(bn(0:N-1,0:N-1))

    !Initialise arrays
    an=0.0; bn=0.0; b=0.0
    a(:,0) = 2*f
    a(:,1) = a(:,0)

    !Set contants and initial variables.
    pi = acos(-1.0)
    logN  = log(real(N))/log(2.0)
    m = 1
    kmax = N/2 -1
    kplus = N/2


!---BEGIN FFT ITERATION-----------------------------------------------------------------------------------------------------------------

    DO i=1,logN
       WRITE(6,*) 'Iteration no.',i

       !$omp PARALLEL
       nprocs = omp_get_num_threads()
       NperCPU = ((kmax+1)*(m+1)) / nprocs
       !$omp END PARALLEL

       !$omp PARALLEL PRIVATE(me, lbeg, lend, k, l, j), SHARED(NperCPU, an, bn)

       me = omp_get_thread_num()
       lbeg = NperCPU*me
       lend = lbeg + NperCPU -1

          DO l=lbeg,lend

            !The following lines act as rules to allocate the correct values to j and k.
            j = floor(real(l)/(kmax+1))
            k = l - j*(kmax+1)

             an(k,j) = 0.5*(a(k,j) + cos(j*pi/m)*a(k+kplus,j) &
                                   - sin(j*pi/m)*b(k+kplus,j))
             bn(k,j) = 0.5*(b(k,j) + sin(j*pi/m)*a(k+kplus,j) &
                                   + cos(j*pi/m)*b(k+kplus,j))

          END DO

       !$omp END PARALLEL

       !Symmetry relations.
       IF(i /= logN) THEN
          DO j=m+1,2*m
             an(:,j) = an(:,2*m-j)
             bn(:,j) = -bn(:,2*m-j)
          END DO
       END IF

       kplus = kplus/2
       m = 2*m
       kmax = (kmax+1)/2-1
       a=an; b=bn
    END DO


!---FFTW-------------------------------------------------------------------------------------------------------------------------

    ALLOCATE(fftw_in(0:(N-1)/2))
    ALLOCATE(fftw_out(0:(N-1)))

    a(0,0) = a(0,0)*2

    fftw_in=cmplx(a(0,:)*N/2,-b(0,:)*N/2)

    plan = fftw_plan_dft_c2r_1d(size(fftw_out),fftw_in, fftw_out,FFTW_ESTIMATE);

    CALL fftw_execute_dft_c2r(plan,fftw_in,fftw_out)

    !Write FFTW reverse and pre-transform data to file to check if they are the same.
    DO i=0,N-1
        WRITE(42,'(2(f8.4,1x))') REAL(fftw_out(i)/N), f(i)
    END DO

    CALL fftw_destroy_plan(plan)


!---Output-----------------------------------------------------------------------------------------------------------------------

    !Writes out data for function before transform.
    DO i=0,N-1
        WRITE(11,'(2(f8.4,1x))') pi/(N/2)*i,f(i)
    END DO

    !Output Fourier series.
    a(0,0) = a(0,0)/2
    DO i=0,1000
       x= i*2*pi/1001
       fsum = 0.0
       DO j=0,N/2-1
          fsum = fsum + a(0,j)*cos(j*x) + b(0,j)*sin(j*x)
       END DO
       WRITE(77,*) x,fsum
    END DO

    !Output coefficients
    WRITE(97,*) a(0,0:N/2-1)
    WRITE(97,*) b(0,0:N/2-1)

    WRITE(6,*) 'SUCCESS!'
    WRITE(6,*) 'Data files contain output as follows:'
    WRITE(6,*) 'fort.42 - FFTW reverse data and pre-transform data'
    WRITE(6,*) 'fort.11 - Data for function prior transformation'
    WRITE(6,*) 'fort.97 - A and B coefficients'

END PROGRAM fft
