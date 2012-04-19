module workingprecision
    ! This module defines the working precision of the routines.
    implicit none
    integer, parameter :: sp = kind(1.0), dp = kind(1.0d0)
    ! Working precision is currently double precision.
    integer, parameter :: wp = dp
end module
