program powermethod
    use workingprecision
    use matrixtools
    implicit none
    integer, parameter :: n = 500
    integer :: i
    real(kind=wp), allocatable :: A(:,:), z(:), q(:)

    allocate(A(n,n), z(n), q(n))

    call init_random_seed()

    call random_number(A)
    call random_number(q)

    q = q / sqrt(dot_product(q,q))

    do i = 1, 1000
        z = matmul(A,q)
        q = z / sqrt(dot_product(z,z))
    end do
    z = matmul(A,q)
    print *, z(1)/q(1)

    print *, sum(z/q)/n

end program powermethod
