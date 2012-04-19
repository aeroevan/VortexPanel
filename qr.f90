module qr
    use workingprecision
    use matrixtools
    implicit none
    ! This module contains the QR related functions and subroutines for both
    ! Gives and Householder methods.
contains
    ! {{{ Householder QR subroutines
    ! Given in input vector, compute the Householder vector.
    pure subroutine house(x, v, b)
        real(kind=wp), intent(in) :: x(:)
        real(kind=wp), intent(out) :: v(:), b
        real(kind=wp) :: s, u
        integer n

        n = size(x,1)
        s = dot_product(x(2:n), x(2:n))
        v(1) = 1_wp
        v(2:n) = x(2:n)

        if (s == 0) then
            b = 0
        else
            u = sqrt(x(1)*x(1)+s)
            if (x(1) <= 0) then
                v(1) = x(1) - u
            else
                v(1) = -s/(x(1)+u)
            end if
            b = 2*v(1)*v(1)/(s+v(1)*v(1))
            v = v / v(1)
        end if
    end subroutine house

    ! Use Householder vectors to do QR decomposition.
    subroutine house_qr(A, Q, R)
        real(kind=wp), intent(in) :: A(:,:)
        real(kind=wp), intent(out) :: R(size(A,1),size(A,2)), &
            Q(size(A,1),size(A,1))
        real(kind=wp) :: W(size(A,1),size(A,2)), b(size(A,2)), &
            v(size(A,1)), v_outer(size(A,1),size(A,1))
        integer :: m, n, i

        ! I coule just overwrite A inplace...
        W = A

        m = size(A,1)
        n = size(A,2)

        b = 0_wp

        do i = 1, n
            call house(W(i:m,i), v(1:m-i+1), b(i))
            call outer_product(v(1:m-i+1), v(1:m-i+1), v_outer(1:m-i+1,1:m-i+1))
            W(i:m,i:n) = W(i:m,i:n) - b(i)*matmul(v_outer(1:m-i+1,1:m-i+1), W(i:m,i:n))
            if (i < m) then
                W(i+1:m,i) = v(2:m-i+1)
            end if
        end do

        ! We only want the upper triangular parts of W.
        call triu(W, R)

        call eye(m, Q)
        do i = n, 1, -1
            v(1) = 1_wp
            v(2:m-i+1) = W(i+1:m,i)
            call outer_product(v(1:m-i+1), v(1:m-i+1), v_outer(1:m-i+1,1:m-i+1))
            Q(i:m,i:m) = Q(i:m,i:m) - b(i)*matmul(v_outer(1:m-i+1,1:m-i+1), Q(i:m,i:m))
        end do
    end subroutine house_qr
    ! }}}
    ! {{{ Givens QR subroutines
    pure subroutine givens(a, b, c, s)
        real(kind=wp), intent(in) :: a, b
        real(kind=wp), intent(out) :: c, s
        real(kind=wp) :: tau

        if (b == 0.0_wp) then
            c = 1.0_wp
            s = 0.0_wp
        else if (abs(b) > abs(a)) then
            tau = -a/b
            s = 1.0_wp/(sqrt(1.0_wp+tau*tau))
            c = s*tau
        else
            tau = -b/a
            c = 1.0_wp/(sqrt(1.0_wp+tau*tau))
            s = c*tau
        end if
    end subroutine givens
    subroutine givens_row(A, c, s)
        real(kind=wp), intent(inout) :: A(:,:)
        real(kind=wp), intent(in) :: c, s
        real(kind=wp) :: tau, sigma
        integer :: k, q

        q = size(A,2)

        do k = 1, q
            tau = A(1,k)
            sigma = A(2,k)
            A(1,k) = c*tau - s*sigma
            A(2,k) = s*tau + c*sigma
        end do
    end subroutine givens_row
    subroutine givens_qr(A, Q, R)
        real(kind=wp), intent(in) :: A(:,:)
        real(kind=wp), intent(out) :: Q(size(A,1), size(A,2)), R(size(A,1), size(A,2))
        integer :: n, m, i, j, T
        real(kind=wp) :: c, s
        logical :: updated

        n = size(A,1)
        m = size(A,2)
        
        Q = 0.0_wp
        do i = 1, n
            Q(i,i) = 1.0_wp
        end do

        updated = .true.
        T = 1

        R = A

        do while(updated)
            updated = .false.
            !$OMP parallel do private(j,c,s)
            do i = n, 2, -1
                do j = 1, i-1
                    if (i-2*j==n-1-T) then
                        updated = .true.
                        call givens(R(i-1,j), R(i,j), c, s)
                        ! Only need to update j:m here (rest should be zero)
                        call givens_row(R(i-1:i,j:m), c, s)
                        call givens_row(Q(i-1:i,1:m), c, s)
                    end if
                end do
            end do
            !$OMP end parallel do
            T = T + 1
        end do
        Q = transpose(Q)
    end subroutine givens_qr
    ! }}}
    ! {{{ Linear system drivers
    function solve_givens_qr(A, b) result(x)
        ! Solve a linear system using QR decomposition
        real(kind=wp) :: A(:,:), b(:), x(size(b))
        real(kind=wp) :: Q(size(A,1), size(A,2)), R(size(A,1), size(A,2))
        integer :: l, m, n
        n = size(A, 1)
        m = size(A, 2)
        if (n /= m) then
            print *, "n /= m, A is not square"
            stop
        end if
        l = size(b, 1)
        if (n /= l) then
            print *, "n /= l, A and b are of different sizes"
            stop
        end if

        call givens_qr(A, Q, R)

        ! Ax = QRx = b => Rx = Q'b
        Q = transpose(Q)
        b = matmul(Q, b)
        x = backsolve(R, b)
    end function solve_givens_qr
    function solve_house_qr(A, b) result(x)
        ! Solve a linear system using QR decomposition
        real(kind=wp) :: A(:,:), b(:), x(size(b))
        real(kind=wp) :: Q(size(A,1), size(A,2)), R(size(A,1), size(A,2))
        integer :: l, m, n
        n = size(A, 1)
        m = size(A, 2)
        if (n /= m) then
            print *, "n /= m, A is not square"
            stop
        end if
        l = size(b, 1)
        if (n /= l) then
            print *, "n /= l, A and b are of different sizes"
            stop
        end if

        call house_qr(A, Q, R)

        ! Ax = QRx = b => Rx = Q'b
        Q = transpose(Q)
        b = matmul(Q, b)
        x = backsolve(R, b)
    end function solve_house_qr
    ! }}}
end module qr
! vim: set foldmethod=marker :
