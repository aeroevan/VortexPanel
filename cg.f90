module cg
    use workingprecision
    use matrixtools
    use mpi
    implicit none
contains
    function congrad(A, b) result(x)
        real(kind=wp) :: A(:,:), b(:), x(size(b)), r(size(b))
        real(kind=wp) :: AtA(size(A,1), size(A,2)), bt(size(b)), p(size(b))
        real(kind=wp) :: Ap(size(A,1))
        real(kind=wp) :: rsold, rsnew, alph
        real(kind=wp), parameter :: tol = 1e-10
        integer :: i, n

        n = size(A,1)

        x = 0_wp

        AtA = matmul(transpose(A), A)
        bt = matmul(transpose(A), b)
        r = bt - matmul(AtA,x)
        p = r
        rsold = dot_product(r,r)
        do i = 1, n
            Ap = matmul(AtA,p)
            alph = rsold/(sum(p*Ap))
            x = x + alph*p
            r = r - alph*Ap
            rsnew = dot_product(r,r)
            if (rsnew < epsilon(1.0)**2) then
                print *, "Converged!"
                exit
            end if
            p = r + rsnew/rsold*p
            rsold = rsnew
        end do
    end function congrad
    subroutine solve_cg(A, b, x)
        ! Solve a linear system using QR decomposition
        real(kind=wp) :: A(:,:), b(:), x(size(b))
        real(kind=wp) :: AT(size(A,1),size(A,2))
        integer :: dest, i, j, k, ierr, my_rank, p, n
        real(kind=wp), allocatable, dimension(:) :: buf, me, temp2
        real(kind=wp), allocatable, dimension(:) :: ans, ans2, x1, p1, q, r
        real(kind=wp), allocatable, dimension(:,:) :: tempA, tempAT
        real(kind=wp) :: alpha, rho, rho0
        real(kind=wp), parameter :: tol = 0.001_wp
        integer :: count_rows, rows, row_index, sender, status(MPI_STATUS_SIZE)

        call MPI_COMM_SIZE(MPI_COMM_WORLD, p, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

        n = size(A,1)

        dest = 0
        rows = n/(p-1)

        AT = transpose(A)

        if (my_rank == dest) then
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
                    !temp2 = matmul(A,p1)
                    DO i=1,p-1
                            CALL MPI_SEND(p1,n,MPI_REAL,i,0,MPI_COMM_WORLD,ierr)
                    END DO
                    DO i=1,p-1
                            CALL MPI_RECV(ans,rows,MPI_REAL,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                            row_index = status(MPI_TAG)
                            temp2(row_index:row_index+rows-1) = ans
                    END DO

                    !multiply AT*(temp2 = A*p1) = q
                    !q = matmul(AT, temp2)
                    DO i=1,p-1
                            CALL MPI_SEND(temp2,n,MPI_REAL,i,1,MPI_COMM_WORLD,ierr)

                    END DO
                    DO i=1,p-1
                            CALL MPI_RECV(ans2,rows,MPI_REAL,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                            row_index = status(MPI_TAG)
                            q(row_index:row_index+rows-1) = ans2
                    END DO
                    
                    ! use updated q to set up next iteration
                    alpha = rho / DOT_PRODUCT( p1(1:n), q(1:n) )
                    ! update solution and residual
                    x1 = x1(1:n) + alpha * p1(1:n) !saxpy
                    r = r(1:n) - alpha * q(1:n) !saxpy
                    
                    ! check to see if residual is within tolerance
                    IF ( rho < epsilon(1.0)**2 ) EXIT
            END DO

            !print solution x1
            !print *, 'X =', x1
            !et = dtime(t)
            !print *, 'Elapsed time =', et
            x = x1

        ELSE !other nodes receive work from master node
            allocate(buf(n))
            allocate(tempA(rows,n))
            allocate(tempAT(rows,n))

            row_index = (my_rank-1)*rows+1

            tempA(1:rows,1:n) = A(row_index:row_index+rows-1,1:n)
            tempAT(1:rows,1:n) = AT(row_index:row_index+rows-1,1:n)

            DO k = 1, 2*(n-1)
                    CALL MPI_RECV(buf,n,MPI_REAL,dest,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                    IF (status(MPI_TAG) == 0) THEN
                            ans = MATMUL(tempA(1:rows,1:n),buf(1:n))
                            CALL MPI_SEND(ans,rows,MPI_REAL,dest,row_index,MPI_COMM_WORLD,ierr)
                    ELSE
                            ans = MATMUL(tempAT(1:rows,1:n),buf(1:n))
                            CALL MPI_SEND(ans,rows,MPI_REAL,dest,row_index,MPI_COMM_WORLD,ierr)
                    END IF
            END DO
            CALL MPI_FINALIZE(ierr)
            stop
        END IF

    end subroutine solve_cg
end module cg
