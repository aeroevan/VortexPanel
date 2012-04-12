PROGRAM PPANEL
    use mpi
    use workingprecision
    IMPLICIT NONE
    INTEGER :: dest, i, j, k, ierr, my_rank, p, n
    integer, parameter :: Nseg = 500, Nact=2*Nseg-2
    real(kind=wp), parameter :: xx = 0.12_wp
    REAL, allocatable, dimension(:,:) :: A, AT
    REAL, allocatable, dimension(:) :: buf, me, temp2
    REAL, allocatable, dimension(:) :: ans, ans2, b, x1, p1, q, r
    REAL, allocatable, dimension(:,:) :: tempW, tempAT
    INTEGER count_rows, rows, row_index, sender, status(MPI_STATUS_SIZE)
    ! set up tolerances
    REAL :: alpha, rho, rho0, tol = .001
    REAL et, dtime, t(2)

    ! get info on processors/cores
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,p,ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
    
    ! ****get A and b matrix for given problem (call Evan's code). need a consistent 'n'
    ! CALL SUBROUTINE... should return A, b, n (****change next three lines)
    n = Nact + 1
    allocate(b(n))
    allocate(A(n,n))
    n = 100 ! number of panels
    
    ! set up problem
    dest = 0
    rows = n/(p-1)
    allocate(AT(n,n))

    ! transpose A matrix and store in A_T since for CG, we must solve AT*A*x=AT*b
    AT = transpose(A)
    CALL TRANS(n,A,AT)

    ! have master node allocate the problem to different nodes
    IF (my_rank == dest) THEN
        allocate(p1(n))
        allocate(q(n))
        allocate(r(n))
        allocate(temp2(n))
        allocate(x1(n))
        allocate(ans(rows))
        allocate(ans2(rows))

        ! CG method (see http://en.wikipedia.org/wiki/Conjugate_gradient_method)
        x1 = 0.0 ! An initial guess
        ! since we're solving AT*A*x=AT*b
        b = MATMUL(AT, b(1:n))
        ! r = residual
        r = b(1:n) - MATMUL(AT, MATMUL(A,x1))
        p1 = r(1:n)
        
        rho = DOT_PRODUCT( r(1:n), r(1:n) )
        q = MATMUL(AT, MATMUL(A,p1(1:n)))
        alpha = rho / DOT_PRODUCT( p1(1:n), q(1:n) )
        
        ! update solution and residual
        x1 = x1(1:n) + alpha * p1(1:n) !saxpy
        r = r(1:n) - alpha * q(1:n) !saxpy

        ! divide and conquer matrix-vector multiplications
        DO k = 2, n
                rho0 = rho
                rho = DOT_PRODUCT( r(1:n), r(1:n) )
                p1 = r(1:n) + ( rho/rho0 ) * p1(1:n)

                !start
                !multiply A*p1 = temp2
                DO i=1,p-1
                        CALL MPI_SEND(p1,n,MPI_DOUBLE_PRECISION,i,0,MPI_COMM_WORLD,ierr)
                END DO
                DO i=1,p-1
                        CALL MPI_RECV(ans,rows,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                        row_index = status(MPI_TAG)
                        temp2(row_index:row_index+rows-1) = ans
                END DO

                !multiply AT*(temp2 = A*p1) = q
                DO i=1,p-1
                        CALL MPI_SEND(temp2,n,MPI_DOUBLE_PRECISION,i,1,MPI_COMM_WORLD,ierr)

                END DO
                DO i=1,p-1
                        CALL MPI_RECV(ans2,rows,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                        row_index = status(MPI_TAG)
                        q(row_index:row_index+rows-1) = ans2
                END DO
                
                ! use updated q to set up next iteration
                alpha = rho / DOT_PRODUCT( p1(1:n), q(1:n) )
                ! update solution and residual
                x1 = x1(1:n) + alpha * p1(1:n) !saxpy
                r = r(1:n) - alpha * q(1:n) !saxpy
                
                ! check to see if residual is within tolerance
                IF ( rho < tol ) EXIT
        END DO

        !print solution x1
        print *, 'X =', x1
        et = dtime(t)
        print *, 'Elapsed time =', et

    ELSE !other nodes receive work from master node
        allocate(buf(n))
        allocate(tempA(rows,n))
        allocate(tempAT(rows,n))

        row_index = (my_rank-1)*rows+1

        tempA(1:rows,1:n) = A(row_index:row_index+rows-1,1:n)
        tempAT(1:rows,1:n) = AT(row_index:row_index+rows-1,1:n)

        DO k = 1, 2*(n-1)
                CALL MPI_RECV(buf,n,MPI_DOUBLE_PRECISION,dest,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                IF (status(MPI_TAG) == 0) THEN
                        ans = MATMUL(tempA(1:rows,1:n),buf(1:n))
                        CALL MPI_SEND(ans,rows,MPI_DOUBLE_PRECISION,dest,row_index,MPI_COMM_WORLD,ierr)
                ELSE
                        ans = MATMUL(tempAT(1:rows,1:n),buf(1:n))
                        CALL MPI_SEND(ans,rows,MPI_DOUBLE_PRECISION,dest,row_index,MPI_COMM_WORLD,ierr)
                END IF
        END DO
    END IF

    CALL MPI_FINALIZE(ierr)
    
    contains
    !just a normal transpose function
    SUBROUTINE TRANS(n,A,B)
        Integer n, i, j
        REAL, dimension(:,:) :: A,B
        DO i = 1,n
                DO j=1,n
                        B(j,i) = A(i,j)
                END DO
        END DO
    END SUBROUTINE TRANS

END PROGRAM PPANEL                                                                                                                         
