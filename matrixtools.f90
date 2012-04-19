module matrixtools
    use workingprecision
    implicit none
    ! This module contains some common routines to both CG and QR methods.
contains

    subroutine print_matrix(A)
        real(kind=wp), intent(in) :: A(:,:)
        integer :: i, n

        n = size(A,1)
        do i = 1, n
            print *, A(i,:)
        end do
    end subroutine print_matrix

    function backsolve(R, b) result(x)
        ! Backsolve right triangular systems as is common with QR
        ! decomposition
        real(kind=wp) :: R(:, :), b(:), x(size(b))
        integer :: i, j, n
        n = size(R, 1)

        ! Solve Q x = b
        do j = n, 1, -1
            x(j) = b(j) / R(j, j)
            !$OMP parallel do
            do i = j-1, 1, -1
                b(i) = b(i) - R(i, j)*x(j)
            end do
            !$OMP end parallel do
        end do
    end function backsolve


    subroutine tril(A, B)
        real(kind=wp), intent(in) :: A(:,:)
        real(kind=wp), intent(out) :: B(size(A,1),size(A,2))
        integer :: m, n, i, j

        m = size(A,1)
        n = size(A,2)

        B = 0_wp
        !$OMP parallel do private(j)
        do i = 1, m
            do j = 1, i
                B(i,j) = A(i,j)
            end do
        end do
        !$OMP end parallel do
    end subroutine

    subroutine triu(A, B)
        real(kind=wp), intent(in) :: A(:,:)
        real(kind=wp), intent(out) :: B(size(A,1),size(A,2))
        integer :: m, n, i, j

        m = size(A,1)
        n = size(A,2)

        B = 0_wp
        !$OMP parallel do private(j)
        do i = 1, m
            do j = i, n
                B(i,j) = A(i,j)
            end do
        end do
        !$OMP end parallel do
    end subroutine

    subroutine eye(n, A)
        integer, intent(in) :: n
        real(kind=wp), intent(out) :: A(n,n)
        integer :: i
        A = 0_wp
        !$OMP parallel do
        do i = 1, n
            A(i,i) = 1_wp
        end do
        !$OMP end parallel do
    end subroutine eye

    subroutine outer_product(x, y, A)
        real(kind=wp), intent(in) :: x(:), y(:)
        real(kind=wp), intent(out) :: A(size(x,1),size(y,1))
        integer :: m, n, i, j

        m = size(x,1)
        n = size(y,1)

        !$OMP parallel do private(j)
        do i = 1, m
            do j = 1, n
                A(i,j) = x(i)*y(j)
            end do
        end do
        !$OMP end parallel do
    end subroutine outer_product

    subroutine init_random_seed()
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size=n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed = clock + 37 * [(i -1, i=1, n)]

        call random_seed(put = seed)
        deallocate(seed)
    end subroutine init_random_seed

end module matrixtools
