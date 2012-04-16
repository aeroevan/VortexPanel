program project
    use workingprecision
    use matrixtools
    use cg
    implicit none
    real(kind=wp), parameter :: pi = 3.14159265358979323846264338327950288_wp
    integer, parameter :: Nseg = 500, N=2*Nseg-2
    ! NACA 0012 => 12% thickness (0.12)
    real(kind=wp), parameter :: xx = 0.12_wp
    real(kind=wp) :: alpha, dx, dy, t1, t2, cy, cx, cm, cl, cd, xarm
    real(kind=wp), dimension(2*Nseg-1) :: x, y
    real(kind=wp), dimension(N+1,N+1) :: A
    real(kind=wp), dimension(N) :: ds, xmid, ymid, cp
    real(kind=wp), dimension(N+1) :: rhs, gam
    character(len=100) :: buffer
    integer :: i

    call getarg(1,buffer)
    if (buffer(1:2) == "-h") then
        print *, "./project [alpha]"
        stop
    elseif (trim(buffer) == '') then
        print *, "Need alpha..."
        print *, "./project [alpha]"
        stop
    end if
    read(buffer,*) alpha
    ! deg => rad
    alpha = pi/180_wp*alpha

    call build_panel(alpha, Nseg, N, xx, x, y, ds, xmid, ymid, A, rhs)

    gam = 0_wp

    gam = congrad(A, rhs)

    !$OMP parallel do
    do i = 2, N-1
        cp(i) = 1_wp - gam(i)*gam(i)
    end do
    !$OMP end parallel do
    !cp(1) = -cp(1)
    !cp(N) = -cp(N)
    do i = 2, N-1
        print *, xmid(i), cp(i), ymid(i), gam(i)
    end do

    cy = 0_wp
    cx = 0_wp
    cm = 0_wp
    !$OMP parallel do private(xarm)
    do i = 1, N
        dx = x(i+1) - x(i)
        dy = y(i+1) - y(i)
        ! moment arm = midpoint of current point to the quarter chord.
        xarm = xmid(i)-x(Nseg)-0.25_wp
        cy = cy - cp(i)*dx
        cx = cx + cp(i)*dy
        cm = cm - cp(i)*dx*xarm
    end do
    !$OMP end parallel do

    cl = cy*cos(alpha) - cx*sin(alpha)
    cd = cy*sin(alpha) + cx*cos(alpha)

    print *, "#", cl, cd, cm
contains
    subroutine build_panel(alpha, Nseg, N, xx, x, y, ds, xmid, ymid, A, rhs)
        real(kind=wp), intent(in) :: alpha, xx
        integer, intent(in) :: Nseg, N
        real(kind=wp), intent(out), dimension(2*Nseg-1) :: x, y
        real(kind=wp), intent(out), dimension(N+1,N+1) :: A
        real(kind=wp), intent(out), dimension(N) :: ds, xmid, ymid
        real(kind=wp), intent(out), dimension(N+1) :: rhs
        integer :: i, j
        ! Upper surface
        !$OMP parallel do
        do i = Nseg, 2*Nseg-1
            x(i) = real(i-Nseg)/Nseg
            y(i) = naca00xx(xx, x(i))
        end do
        !$OMP end parallel do
        ! Lower surface is symmetric... index so bottom then top
        !$OMP parallel do
        do i = 1, Nseg
            x(Nseg+1-i) = x(Nseg-1+i)
            y(Nseg+1-i) = -y(Nseg-1+i)
        end do
        !$OMP end parallel do

        ! Compute panel sizes
        !$OMP parallel do
        do i = 1, N
            t1 = x(i+1) - x(i)
            t2 = y(i+1) - y(i)
            ds(i) = sqrt(t1*t1 + t2*t2)
        end do
        !$OMP end parallel do

        ! Compute RHS
        rhs = 0_wp
        xmid = 0_wp
        ymid = 0_wp

        !$OMP parallel do
        do i = 1, N
            xmid(i) = 0.5_wp * (x(i) + x(i+1))
            ymid(i) = 0.5_wp * (y(i) + y(i+1))
            rhs(i) = ymid(i) * cos(alpha) - xmid(i) * sin(alpha)
        end do
        !$OMP end parallel do

        ! Parallelize this...
        A = 0_wp
        !$OMP parallel do private(i,j)
        do i = 1, N
            A(i,N+1) = 1_wp
            do j = 1, N
                A(i, j) = make_A(x, y, ds, xmid, ymid, i, j)
            end do
        end do
        !$OMP end parallel do
        ! Kutta condition
        A(N+1,1) = 1_wp
        A(N+1,N) = 1_wp
    end subroutine build_panel

    pure function make_A(x, y, ds, xmid, ymid, i, j) result(aij)
        real(kind=wp), dimension(2*Nseg-1), intent(in) :: x, y
        real(kind=wp), dimension(N), intent(in) :: ds, xmid, ymid
        integer, intent(in) :: i, j
        real(kind=wp) :: aij
        real(kind=wp) :: dx, dy, t1, t2, t3, t4, t5, t6, t7
        if (i == j) then
            aij = ds(i)/(2_wp*pi) * (log(0.5_wp*ds(i)) - 1_wp)
        else
            dx  = (x(j+1)-x(j))/ds(j);
            dy  = (y(j+1)-y(j))/ds(j);
            t1  = x(j) - xmid(i);
            t2  = y(j) - ymid(i);
            t3  = x(j+1) - xmid(i);
            t4  = y(j+1) - ymid(i);
            t5  = t1 * dx + t2 * dy;
            t6  = t3 * dx + t4 * dy;
            t7  = t2 * dx - t1 * dy;
            t1  = t6 * log(t6*t6+t7*t7) - t5 * log(t5*t5+t7*t7);
            t2  = atan2(t7,t5)-atan2(t7,t6);
            aij = (0.5_wp * t1-t6+t5+t7*t2)/(2_wp*pi);
        end if
    end function make_A
    pure function naca00xx(xx, x, c) result(y)
        real(kind=wp), intent(in) :: xx, x
        real(kind=wp), intent(in), optional :: c
        real(kind=wp) :: y
        if (.not.present(c)) then
            ! Assume c = 1
            y = xx/0.2_wp*(0.2969_wp*sqrt(x) &
                - 0.1260_wp*(x) - 0.3516_wp*(x)**2 &
                + 0.2843_wp*(x)**3 - 0.1015_wp*(x)**4)
        else
            y = xx/0.2_wp*c*(0.2969_wp*sqrt(x/c) &
                - 0.1260_wp*(x/c) - 0.3516_wp*(x/c)**2 &
                + 0.2843_wp*(x/c)**3 - 0.1015_wp*(x/c)**4)
        end if
    end function naca00xx
end program project

