module cg
    use workingprecision
    use matrixtools
    use mpi
    implicit none
    ! This module contains the conjugate gradient related functions and
    ! subroutines.
contains
    function congrad(A, b) result(x)
        ! Simple sequential algorithm for the CG method.
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
        ! Use MPI to distribute the matrix multiplication
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
            avgrow = n/numworkers ! Since we are working with transpose(A),
                                  ! rows = n
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
                    call MPI_RECV(AtA(offset:offset+rows-1,:), rows*n, &
                        MPI_REAL, source, from_worker, MPI_COMM_WORLD, &
                        status, ierr)
                    call MPI_RECV(Atb(offset:offset+rows-1), rows, MPI_REAL, &
                        source, from_worker, MPI_COMM_WORLD, status, ierr)
                else
                    call MPI_RECV(AtA(offset:offset+rows-1,:), rows*n, &
                        MPI_DOUBLE_PRECISION, source, from_worker, &
                        MPI_COMM_WORLD, status, ierr)
                    call MPI_RECV(Atb(offset:offset+rows-1), rows, &
                        MPI_DOUBLE_PRECISION, source, from_worker, &
                        MPI_COMM_WORLD, status, ierr)
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
                        call MPI_RECV(Ap(offset:offset+rows-1), rows, &
                            MPI_REAL, source, from_worker, MPI_COMM_WORLD, &
                            status, ierr)
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
                call MPI_RECV(At, rows*m, MPI_DOUBLE_PRECISION, 0, &
                    from_master, MPI_COMM_WORLD, status, ierr)
            end if
            allocate(AtA(rows,n), bt(rows))
            AtA = matmul(At,A)
            bt = matmul(At,b)

            ! Send results back
            call MPI_SEND(mdata, 2, MPI_INTEGER, 0, from_worker, &
                MPI_COMM_WORLD, ierr)
            if (wp == sp) then
                call MPI_SEND(AtA, rows*n, MPI_REAL, 0, from_worker, &
                    MPI_COMM_WORLD, ierr)
                call MPI_SEND(bt, rows, MPI_REAL, 0, from_worker, &
                    MPI_COMM_WORLD, ierr)
            else
                call MPI_SEND(AtA, rows*n, MPI_DOUBLE_PRECISION, 0, &
                    from_worker, MPI_COMM_WORLD, ierr)
                call MPI_SEND(bt, rows, MPI_DOUBLE_PRECISION, 0, &
                    from_worker, MPI_COMM_WORLD, ierr)
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
                call MPI_SEND(mdata, 2, MPI_INTEGER, 0, from_worker, &
                    MPI_COMM_WORLD, ierr)
                if (wp == sp) then
                    call MPI_SEND(bt, rows, MPI_REAL, 0, from_worker, &
                        MPI_COMM_WORLD, ierr)
                else
                    call MPI_SEND(bt, rows, MPI_DOUBLE_PRECISION, 0, &
                        from_worker, MPI_COMM_WORLD, ierr)
                end if
            end do
            call MPI_FINALIZE(ierr)
            stop
        end if
    end subroutine solve_cg
end module cg
