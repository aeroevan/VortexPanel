program test_qr
    use workingprecision
    use matrixtools
    use qr
    implicit none
    integer, parameter :: n = 500
    real(kind=wp) :: A(n,n), Q(n,n), R(n,n)
    integer :: i, j

    do i = 1, n
        do j = 1, n
            A(i,j) = i+j
        end do
    end do

    print *, "A:"
    !call print_matrix(A)
    print *, ""


    call house_qr(A, Q, R)

    print *, "Q:"
    !call print_matrix(Q)
    print *, ""

    print *, "Q'*Q:"
    !call print_matrix(matmul(transpose(Q),Q))
    print *, ""

    print *, "R:"
    !call print_matrix(R)
    print *, ""

    print *, "Q*R:"
    !call print_matrix(matmul(Q,R))
end program test_qr
