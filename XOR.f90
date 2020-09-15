!Sigmoid Function
FUNCTION f(k,net)

    IMPLICIT NONE

    REAL, INTENT(in) :: net, k
    REAL :: f
    f = 1.0/(1.0+exp(-k*net))
END FUNCTION f


!This program trains a neural network to reproduce the XOR gate. It uses a feed-forward multilayer perceptron
!with 2 nodes in the hidden layer. It uses the Metropolis algorithm as a training rule, minimising the
!'temperature' through each iteration (with exception).
PROGRAM XOR
    IMPLICIT NONE

    REAL, DIMENSION(1:8,1:4) :: truth
    REAL, DIMENSION(1:16) :: w, wn
    REAL, DIMENSION(1:8) :: error
    REAL, DIMENSION(4:7) :: o
    REAL :: globerror, globerrorn, beta, k, random, m, p, targeterror
    INTEGER :: RandomR, i, count, reset, stopper, totalcount

    INTERFACE
        REAL FUNCTION f(k,net)
            REAL, INTENT(in) :: net, k
        END FUNCTION f
    END INTERFACE

!---SET UP CONDITIONS-------------------------------------------------------------------------------------------------------------------------

    k = 42
    beta = 0.1
    count = 0
    totalcount = 0
    stopper = 0
    reset = 0
    targeterror = 1E-7


    !Initialise with with a random number between +-0.5
    CALL random_seed
    CALL random_number(w)
    w(:) = w(:)-0.5
    wn(:) = w(:)

    !Set up truth table.
    truth(1,:) = (/0.0, 0.0, 0.0, 0.0/)
    truth(2,:) = (/0.0, 0.0, 1.0, 1.0/)
    truth(3,:) = (/0.0, 1.0, 0.0, 1.0/)
    truth(4,:) = (/0.0, 1.0, 1.0, 0.0/)
    truth(5,:) = (/1.0, 0.0, 0.0, 1.0/)
    truth(6,:) = (/1.0, 0.0, 1.0, 0.0/)
    truth(7,:) = (/1.0, 1.0, 0.0, 0.0/)
    truth(8,:) = (/1.0, 1.0, 1.0, 1.0/)

    !Calculate initial value for global error
    DO i=1,8

        o(4) = f(k,w(1) + w(5)*truth(i,1) + w(8)*truth(i,2) + w(11)*truth(i,3))
        o(5) = f(k,w(2) + w(6)*truth(i,1) + w(9)*truth(i,2) + w(12)*truth(i,3))
        o(6) = f(k,w(3) + w(7)*truth(i,1) + w(10)*truth(i,2) + w(13)*truth(i,3))
        o(7) = f(k,w(4) + w(14)*o(4) + w(15)*o(5) + w(16)*o(6))

        error(i) = (o(7)-truth(i,4))**2

    END DO

    globerrorn = 0.5*sum(error)

!---BEGIN NEURAL NETWORK-------------------------------------------------------------------------------------------------------------------

    DO WHILE(stopper < 1)

        !Choose weight in array.
        CALL random_number(random)
        Random = (random*20.0)+1.0
        RandomR = random

        !Allocate random value between +-0.5 to weight.
        CALL random_number(m)
        m = m-0.5
        w(RandomR) = m

        !Computes new outputs and error with new random weight
        DO i=1,8

            o(4) = f(k,w(1) + w(5)*truth(i,1) + w(8)*truth(i,2) + w(11)*truth(i,3))
            o(5) = f(k,w(2) + w(6)*truth(i,1) + w(9)*truth(i,2) + w(12)*truth(i,3))
            o(6) = f(k,w(3) + w(7)*truth(i,1) + w(10)*truth(i,2) + w(13)*truth(i,3))
            o(7) = f(k,w(4) + w(14)*o(4) + w(15)*o(5) + w(16)*o(6))

            error(i) = (o(7)-truth(i,4))**2

        END DO

        globerror = 0.5*sum(error)


!-------METROPOLIS ALGORITHM---------------------------------------------------------------------------------------------------------

        CALL random_number(p)

        !Keep new weight if it reduces the total error.
        IF (globerror - globerrorn < 0) THEN
            wn(RandomR) = w(RandomR)
            globerrorn = globerror
        !Random probability of accepting weight even if it increases error to prevent it getting stuck.
        ELSE IF (p < 1.0*exp(beta*(globerrorn-globerror))) THEN
            wn(RandomR) = w(RandomR)
            globerrorn = globerror
        END IF

        !This resets the algorithm with new initial weights to prevent the program getting stuck.
        IF (count > 1000) THEN
            IF (globerrorn > 0.01) THEN
                CALL random_number(w)
                count = 0
                reset = reset + 1
                beta = 0.1
            END IF
        END IF

        !Finish Neural Network when error is less than predefined level.
        IF (globerrorn < targeterror) THEN
            stopper = 1
        END IF

        beta = beta + 0.1
        totalcount = totalcount + 1
        count = count + 1
        w = wn

    END DO


!---OUTPUT TO SCREEN-------------------------------------------------------------------------------------------------------------------------------

    WRITE(6,*)
    WRITE(6,*) '                            XOR Table'
    WRITE(6,*)
    WRITE(6,*) '       Input 1    Input 2     Input 3   Ouput-Table Output-Network'
    WRITE(42,*) 'Actual Output-Network'
    WRITE(42,*)


    DO i = 1,8

        o(4) = f(k,w(1) + w(5)*truth(i,1) + w(8)*truth(i,2) + w(11)*truth(i,3))
        o(5) = f(k,w(2) + w(6)*truth(i,1) + w(9)*truth(i,2) + w(12)*truth(i,3))
        o(6) = f(k,w(3) + w(7)*truth(i,1) + w(10)*truth(i,2) + w(13)*truth(i,3))
        o(7) = f(k,w(4) + w(14)*o(4) + w(15)*o(5) + w(16)*o(6))

        error(i) = (o(7)-truth(i,4))**2

        WRITE(6,*) nint(truth(i,1)), nint(truth(i,2)), nint(truth(i,3)), nint(truth(i,4)), nint(o(7))
        WRITE(42,*) o(7)


    END DO

    globerror = 0.5*sum(error)

    WRITE(6,*)
    WRITE(6,*) 'Check fort.42 for actual network output for XOR gate and end values.'
    WRITE(6,*)
    WRITE(42,*)


!---OUTPUT TO FILE-----------------------------------------------------------------------------------------------------------------------------------

    WRITE(42,*) 'END VALUES'
    WRITE(42,*)
    WRITE(42,*) 'Bias (4) =', w(1)
    WRITE(42,*) 'Bias (5) =', w(2)
    WRITE(42,*) 'Bias (6) =', w(3)
    WRITE(42,*) 'Bias (7) =', w(4)
    WRITE(42,*)

    WRITE(42,*) 'Weight(41)=',w(5), 'Weight(42)=',w(8), 'Weight(43)=',w(11)
    WRITE(42,*) 'Weight(51)=',w(6), 'Weight(52)=',w(9), 'Weight(53)=',w(12)
    WRITE(42,*) 'Weight(61)=',w(7),'Weight(62)=',w(10), 'Weight(63)=',w(13)
    WRITE(42,*) 'Weight(74)=',w(14), 'Weight(75)=',w(15), 'Weight(76)=',w(16)

    WRITE(42,*)
    WRITE(42,*) 'Global error =', globerror
    WRITE(42,*)
    WRITE(42,*) 'Number of iterations =', totalcount
    WRITE(42,*)
    WRITE(42,*) 'Number of resets =', reset

END PROGRAM xor
