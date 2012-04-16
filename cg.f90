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
    subroutine solve_cg_old(A, b, x)
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
        
        ! rows indicates the number of rows of A each child node will focus
        ! on when multiplying by a vector. the idea is that the results from
        ! each child node can be stacked together to form a complete result.
        ! for example, in Ax=b, if we have 2 child nodes, we would form two
        ! systems: A1*x1=b1 and A2*x2=b2, where A = [A1;A2], x=[x1;x2], b=[b1;b2]
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

            ! divide and conquer matrix-vector multiplications:
            ! this sends the current p1 vector to all nodes. the
            ! nodes each have access the A and AT and based on
            ! the number of these nodes (processors), they each
            ! multiply their own set of rows of A and AT by p1.
            ! the results from each node are compiled into a single
            ! result by the master node. this procedure is iterated
            ! at most n times, stopping when a particular tolerance
            ! is reached.
            DO k = 2, n
                    rho0 = rho
                    rho = DOT_PRODUCT( r(1:n), r(1:n) )
                    p1 = r(1:n) + ( rho/rho0 ) * p1(1:n)

                    !start
                    !multiply A*p1 = temp2
                    !temp2 = matmul(A,p1)
                    DO i=1,p-1
                            ! send current p1 vector to all nodes to be multiplied:
                            CALL MPI_SEND(p1,n,MPI_REAL,i,0,MPI_COMM_WORLD,ierr)
                    END DO
                    DO i=1,p-1
                            CALL MPI_RECV(ans,rows,MPI_REAL,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                            row_index = status(MPI_TAG)
                            ! store vector received from child node into the relevant rows
                            ! of temp2 (since each node only multiples a certain set of rows
                            ! of A...)
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
                    ! ****if it is, close out all child nodes by sending a message...
                    ! (read below --see ****-- for more details)
                    IF ( rho < epsilon(1.0)**2 ) THEN
                        EXIT
                    END IF
            END DO
            
            ! return result
            x = x1

        ELSE ! other nodes receive work from master node
            allocate(buf(n))
            allocate(tempA(rows,n))
            allocate(tempAT(rows,n))

            ! the starting index for this node to work on...
            ! so this node will multiply (row_index:row_index+rows-1) 
            ! rows in A
            row_index = (my_rank-1)*rows+1

            tempA(1:rows,1:n) = A(row_index:row_index+rows-1,1:n)
            tempAT(1:rows,1:n) = AT(row_index:row_index+rows-1,1:n)

            ! ****currently, the child nodes iterate the default maximum number of times
            ! (remember, CG converges in n steps; I use 2*(n-1) below since the child nodes
            ! are called TWICE--once to multiply A and once for AT. HOWEVER, since the
            ! procedure can hault before n steps due to tolerance, the root node probably
            ! needs to send a message in this case telling the child nodes to stop
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
        END IF

        ! this gets called when either the root node or a child node falls out of the do-loops above
        !CALL MPI_FINALIZE(ierr)

    end subroutine solve_cg_old
    subroutine solve_cg(A, b, x)
        real(kind=wp), intent(in) :: A(:,:), b(:)
        real(kind=wp), intent(out) :: x(size(b,1))
        real(kind=wp), allocatable :: At(:,:), AtA(:,:), Atb(:)
        real(kind=wp), allocatable :: Ap(:), bt(:), r(:), p(:)
        real(kind=wp) :: rsnew, rsold, alph
        integer, parameter :: from_master = 1, from_worker = 2
        integer :: numtasks, id, numworkers, source, dest
        integer :: m, n, rows, avgrow, extra, offset, i, k, ierr
        integer :: status(MPI_STATUS_SIZE), mdata(2)

        call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
        numworkers = numtasks - 1 ! Need one for master
        if (numtasks < 2) then
            print *, "Number of processors must be at least 2"
            call MPI_FINALIZE(ierr)
            stop
        end if

        m = size(A,1)
        n = size(A,2)
        if (size(b,1) /= m) then
            print *, "Dimensions of A and b must agree"
            call MPI_FINALIZE(ierr)
            stop
        end if

        if (id == 0) then

            allocate(At(n,m))
            At = transpose(A)

            ! compute: AtA = matmul(At,A) and bt = matmul(At,b)
            avgrow = n/numworkers ! Since we are working with transpose(A), rows = n
            extra = mod(n,numworkers)
            ! Send data to workers:
            offset = 1
            do dest = 1, numworkers
                if (dest <= extra) then
                    rows = avgrow + 1
                else
                    rows = avgrow
                end if
                mdata(1) = offset
                mdata(2) = rows
                call MPI_SEND(mdata, 2, MPI_INTEGER, dest, from_master, &
                    MPI_COMM_WORLD, ierr)
                if (wp == sp) then
                    call MPI_SEND(At(offset:offset+rows-1,:), rows*m, MPI_REAL, &
                        dest, from_master, MPI_COMM_WORLD, ierr)
                else
                    call MPI_SEND(At(offset:offset+rows-1,:), rows*m, &
                        MPI_DOUBLE_PRECISION, dest, from_master, &
                        MPI_COMM_WORLD, ierr)
                end if
                offset = offset + rows
            end do
            ! Get results
            allocate(AtA(n,n),Atb(n))
            do source = 1, numworkers
                call MPI_RECV(mdata, 2, MPI_INTEGER, source, from_worker, &
                    MPI_COMM_WORLD, status, ierr)
                offset = mdata(1)
                rows = mdata(2)
                if (wp == sp) then
                    call MPI_RECV(AtA(offset:offset+rows-1,:), rows*n, MPI_REAL, &
                        source, from_worker, MPI_COMM_WORLD, status, ierr)
                    call MPI_RECV(Atb(offset:offset+rows-1), rows, MPI_REAL, &
                        source, from_worker, MPI_COMM_WORLD, status, ierr)
                else
                    call MPI_RECV(AtA(offset:offset+rows-1,:), rows*n, &
                        MPI_DOUBLE_PRECISION, source, from_worker, MPI_COMM_WORLD, &
                        status, ierr)
                    call MPI_RECV(Atb(offset:offset+rows-1), rows, &
                        MPI_DOUBLE_PRECISION, source, from_worker, MPI_COMM_WORLD, &
                        status, ierr)
                end if
            end do

            allocate(Ap(n), r(n), p(n))
            x = 0_wp
            r = Atb
            p = r
            rsold = dot_product(r,r)
            do i = 1, n
                do dest = 1, numworkers
                    call MPI_SEND(1, 1, MPI_INTEGER, dest, from_master, &
                        MPI_COMM_WORLD, ierr)
                    if (wp == sp) then
                        call MPI_SEND(p, n, MPI_REAL, dest, from_master, &
                            MPI_COMM_WORLD, ierr)
                    else
                        call MPI_SEND(p, n, MPI_DOUBLE_PRECISION, dest, &
                            from_master, MPI_COMM_WORLD, ierr)
                    end if
                end do
                do source = 1, numworkers
                    call MPI_RECV(mdata, 2, MPI_INTEGER, source, from_worker, &
                        MPI_COMM_WORLD, status, ierr)
                    offset = mdata(1)
                    rows = mdata(2)
                    if (wp == sp) then
                        call MPI_RECV(Ap(offset:offset+rows-1), rows, MPI_REAL, &
                            source, from_worker, MPI_COMM_WORLD, status, ierr)
                    else
                        call MPI_RECV(Ap(offset:offset+rows-1), rows, &
                            MPI_DOUBLE_PRECISION, source, from_worker, &
                            MPI_COMM_WORLD, status, ierr)
                    end if
                end do
                alph = rsold/dot_product(p,Ap)
                x = x + alph*p
                r = r - alph*Ap
                rsnew = dot_product(r,r)
                if (rsnew < epsilon(1.0)**2) then
                    do dest = 1, numworkers
                        call MPI_SEND(0, 1, MPI_INTEGER, dest, from_master, &
                            MPI_COMM_WORLD, ierr)
                    end do
                    print *, "Converged!"
                    exit
                end if
                p = r + rsnew/rsold*p
                rsold = rsnew
            end do

        else ! workers

            ! Compute AtA = matmul(At,A) and bt = matmul(At,b)
            call MPI_RECV(mdata, 2, MPI_INTEGER, 0, from_master, &
                MPI_COMM_WORLD, status, ierr)
            offset = mdata(1)
            rows = mdata(2)

            allocate(At(rows,m))
            if (wp == sp) then
                call MPI_RECV(At, rows*m, MPI_REAL, 0, from_master, &
                    MPI_COMM_WORLD, status, ierr)
            else
                call MPI_RECV(At, rows*m, MPI_DOUBLE_PRECISION, 0, from_master, &
                    MPI_COMM_WORLD, status, ierr)
            end if
            allocate(AtA(rows,n), bt(rows))
            AtA = matmul(At,A)
            bt = matmul(At,b)

            ! Send results back
            call MPI_SEND(mdata, 2, MPI_INTEGER, 0, from_worker, MPI_COMM_WORLD, &
                ierr)
            if (wp == sp) then
                call MPI_SEND(AtA, rows*n, MPI_REAL, 0, from_worker, &
                    MPI_COMM_WORLD, ierr)
                call MPI_SEND(bt, rows, MPI_REAL, 0, from_worker, &
                    MPI_COMM_WORLD, ierr)
            else
                call MPI_SEND(AtA, rows*n, MPI_DOUBLE_PRECISION, 0, from_worker, &
                    MPI_COMM_WORLD, ierr)
                call MPI_SEND(bt, rows, MPI_DOUBLE_PRECISION, 0, from_worker, &
                    MPI_COMM_WORLD, ierr)
            end if

            allocate(p(n))

            do i = 1, n
                ! If we send 0, shut down this process.
                call MPI_RECV(k, 1, MPI_INTEGER, 0, from_master, &
                    MPI_COMM_WORLD, status, ierr)
                if (k == 0) then
                    call MPI_FINALIZE(ierr)
                    stop
                end if
                if (wp == sp) then
                    call MPI_RECV(p, n, MPI_REAL, 0, from_master, &
                        MPI_COMM_WORLD, status, ierr)
                else
                    call MPI_RECV(p, n, MPI_DOUBLE_PRECISION, 0, from_master, &
                        MPI_COMM_WORLD, status, ierr)
                end if
                bt = matmul(AtA,p)
                call MPI_SEND(mdata, 2, MPI_INTEGER, 0, from_worker, MPI_COMM_WORLD, &
                    ierr)
                if (wp == sp) then
                    call MPI_SEND(bt, rows, MPI_REAL, 0, from_worker, &
                        MPI_COMM_WORLD, ierr)
                else
                    call MPI_SEND(bt, rows, MPI_DOUBLE_PRECISION, 0, from_worker, &
                        MPI_COMM_WORLD, ierr)
                end if
            end do
            call MPI_FINALIZE(ierr)
            stop
        end if
    end subroutine solve_cg
end module cg
