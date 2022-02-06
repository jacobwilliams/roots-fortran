!*****************************************************************************************
!>
!  Test of the object-oriented interface.

program root_tests_oo

    use root_module, wp => root_module_rk
    use iso_fortran_env, only: output_unit

    implicit none

    real(wp) :: root,ax,bx,xzero,fzero,error
    integer :: iflag
    character(len=1000) :: line

    real(wp),parameter :: ftol = epsilon(1.0_wp) * 10    !1.0e-15_wp
    real(wp),parameter :: rtol = epsilon(1.0_wp) * 100   !1.0e-13_wp
    real(wp),parameter :: atol = epsilon(1.0_wp) * 1     !1.0e-16_wp
    integer,parameter  :: maxiter = 1000

    character(len=*),parameter :: fmt  = '(A25,   1X,A25,   1X,A25,  1X,A5,1X,A5)' !! format for header
    character(len=*),parameter :: dfmt = '(E25.10,1X,E25.10,1X,E25.6,1X,I5,1X,I5)' !! format for results

    type,extends(muller_solver) :: my_solver
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

        class(root_solver),intent(inout) :: me
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
