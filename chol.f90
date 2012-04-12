module cholesky
    use workingprecision
    use matrixtools
    implicit none
contains
    subroutine gaxpychol(A, G)
        real(kind=wp), intent(in) :: A(:,:)
        real(kind=wp), intent(out) :: G(size(A,1),size(A,2))
        real(kind=wp) :: W(size(A,1),size(A,2))
        integer :: m, n, i

        m = size(A,1)
        n = size(A,2)
        if (n /= m) then
            print *, "A must be square (and spd)!"
            stop
        end if

        W = A

        do i = 1, n
            if (i > 1) then
                W(i:n,i) = W(i:n,i) - matmul(W(i:n,1:i-1),W(i,1:i-1))
            end if
            W(i:n,i) = W(i:n,i) / sqrt(W(i,i))
        end do
        call tril(A, G)
    end subroutine gaxpychol
end module cholesky
