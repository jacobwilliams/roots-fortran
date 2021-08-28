!*****************************************************************************************
!>
!  Test of the object-oriented interface.

program root_tests_oo

    use stdlib_root_module
    use iso_fortran_env, only: wp => real64, output_unit

    implicit none

    real(wp) :: root,ax,bx,xzero,fzero,error
    integer :: iflag
    character(len=1000) :: line

    real(wp),parameter :: ftol = 1.0e-15_wp
    real(wp),parameter :: rtol = 1.0e-13_wp
    real(wp),parameter :: atol = 1.0e-16_wp
    integer,parameter  :: maxiter = 1000

    character(len=*),parameter :: fmt  = '(A25,   1X,A25,   1X,A25,  1X,A5,1X,A5)' !! format for header
    character(len=*),parameter :: dfmt = '(E25.10,1X,E25.10,1X,E25.6,1X,I5,1X,I5)' !! format for results

    type,extends(muller_solver_real64) :: my_solver
        integer :: ifunc = 0
        integer :: n = 1
    end type my_solver
    type(my_solver) :: solver

    ! case 19
    ax = 0.0_wp
    bx = 1.0_wp
    root = 0.22370545765466_wp

    call solver%initialize(f,ftol,rtol,atol,maxiter); solver%n = 3
    call solver%solve(ax,bx,xzero,fzero,iflag)

    error = xzero-root

    write(output_unit,fmt) repeat('-',25),repeat('-',25),repeat('-',25),repeat('-',5),repeat('-',5)
    write(output_unit,fmt) 'error','x','f','evals','iflag'
    write(output_unit,fmt) repeat('-',25),repeat('-',25),repeat('-',25),repeat('-',5),repeat('-',5)

    write(line, dfmt) error,xzero,fzero,solver%ifunc,iflag
    write(output_unit, '(A)') trim(line)

    contains

    function f(me,x)

        !! function to find the root of.

        class(root_solver_real64),intent(inout) :: me
        real(wp),intent(in) :: x
        real(wp) :: f

        real(wp) :: dn

        select type(me)
        class is (my_solver)

            me%ifunc = me%ifunc + 1

            dn = real(me%n, wp)
            f = 2.0_wp*x*exp(-dn) - 2.0_wp*exp(-dn*x) + 1.0_wp

        end select

    end function f

    end program root_tests_oo
!*****************************************************************************************
