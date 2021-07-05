!*****************************************************************************************
!>
!  Tests for scalar root functions without derivatives.

    program root_tests

    use iso_fortran_env,    only: wp => real64, ip => int32
    use iso_fortran_env,    only: error_unit, output_unit
    use stdlib_root_module, only: root_scalar
    use face, only: colorize

    implicit none

    integer  :: imeth !! method counter
    integer  :: icase !! case counter
    integer  :: ifunc !! function counter
    real(wp) :: x1    !! left endpoint
    real(wp) :: x2    !! right endpoint
    real(wp) :: root  !! actual root
    real(wp) :: xzero !! root found
    real(wp) :: fzero !! value of function at found root
    real(wp) :: error !! xzero-root
    integer  :: iflag !! status flag
    character(len=1000) :: line

    integer,parameter :: number_of_cases = 17 !! number of functions to test
    integer,parameter :: number_of_methods = 10 !! number of methods to test
    character(len=100),dimension(number_of_methods),parameter :: methods = [ &
        'brent          ', &
        'brentq         ', &
        'brenth         ', &
        'bisection      ', &
        'anderson_bjorck', &
        'ridders        ', &
        'pegasus        ', &
        'bdqrf          ', &
        'muller         ', &
        'chandrupatla   ' ] !! method names

    write(output_unit,*) ''

    do imeth = 1, number_of_methods

        write(output_unit, '(A20,1X,A3,1X,1X,A25,1X,A25,1X,A25,1X,A5,1X,A5)') &
        repeat('-',20),repeat('-',20),repeat('-',25),repeat('-',25),repeat('-',25),repeat('-',5),repeat('-',5)

        write(output_unit, '(A20,1X,A3,1X,1X,A25,1X,A25,1X,A25,1X,A5,1X,A5)') &
        'method','function','error','x','f','evals','iflag'

        write(output_unit, '(A20,1X,A3,1X,1X,A25,1X,A25,1X,A25,1X,A5,1X,A5)') &
        repeat('-',20),repeat('-',20),repeat('-',25),repeat('-',25),repeat('-',25),repeat('-',5),repeat('-',5)

        do icase = 1,number_of_cases

            ifunc = 0

            select case (icase)

            case(1)

                x1  = 0.0_wp
                x2  = 1.0_wp
                root = asin(2.0_wp/3.0_wp)

            case(2)

                x1  = -1.0_wp
                x2  = 1.0_wp
                root = 0.5671432904_wp   ! approximate

            case(3)

                x1  = 0.1_wp
                x2  = 0.9_wp
                root = 0.8041330975_wp   ! approximate

            case(4)

                x1  = 2.8_wp
                x2  = 3.1_wp
                root = 3.0_wp

            case(5)

                x1  = -1.3_wp
                x2  = -0.5_wp
                root = -0.6294464841_wp  ! approximate

            case(6)

                x1  = 2.0_wp
                x2  = 3.0_wp
                root = 2.0945514815_wp   ! approximate

            case(7)

                x1  = 0.5_wp
                x2  = 2.0_wp
                root = 1.0_wp

            case(8)

                x1  = -10.0_wp
                x2  = 10.0_wp
                root = 1.0_wp

            case(9)

                x1  = 0.0_wp
                x2  = 10.0_wp
                root = 1.0_wp

            case(10)

                x1  = 0.0_wp
                x2  = 10.0_wp
                root = 1.0_wp

            case(11)

                x1  = 0.0_wp
                x2  = 10.0_wp
                root = 1.0_wp

            case(12)

                !From: http://www.mth.pdx.edu/~daescu/mth451_551/Muller_bisection.pdf
                x1  = -2.0_wp
                x2  = 3.0_wp
                root = 0.0_wp

            case(13)

                !From: http://www.mth.pdx.edu/~daescu/mth451_551/Muller_bisection.pdf
                x1  = 0.0_wp     ! x1 is the location of the root for this case
                x2  = 3.0_wp
                root = 0.0_wp

            case(14)

                x1  = 0.0_wp
                x2  = 0.7_wp
                root = 0.677712926329_wp       ! approximate

            case(15)

                x1  = -0.7_wp
                x2  = 0.4_wp
                root = -0.640083696085_wp      ! approximate

            case(16)

                !case 11 from "Algorithm 748: Enclosing Zeros of Continuous Functions" (n=20)
                x1  = 0.01_wp
                x2  = 1.0_wp
                root = 0.05_wp

            case(17)

                x1  = -10.0_wp
                x2  = 5.0_wp
                root = 0.0_wp

            end select

            ifunc = 0 ! reset function eval counter

            call root_scalar(methods(imeth),test_func,x1,x2,xzero,fzero,iflag,&
                                atol = 1.0e-16_wp, rtol = 1.0e-12_wp)

            error = xzero-root
            write(line, '(A20,1X,I3,1X,E25.10,1X,E25.10,1X,E25.6,1X,I5,1X,I5)') &
                    trim(methods(imeth)),icase,xzero-root,xzero,fzero,ifunc,iflag

            if (abs(fzero) <= 1.0e-9_wp .or. abs(error) <= 1.0e-9_wp) then
                write(output_unit, '(A)') trim(line)
            else
                write(output_unit, '(A)') colorize(trim(line), color_fg='red')
            end if

        end do

    end do

    write(output_unit,*) ''

    contains
!*****************************************************************************************

!**********************************************************
    function test_func(x) result(f)

    implicit none

    real(wp),intent(in) :: x
    real(wp) :: f

    integer :: n !! for case 16

    ifunc = ifunc + 1   ! update function counter

    ! Note: cases 1-7 are the functions from Gottlieb's paper [table 1 & 2]

    select case (icase)
    case(0)
        f = sind(x) * 123456.789_wp
    case(1)
        f = 3.0_wp * sin(x) - 2.0_wp
    case(2)
        f = x*exp(x) - 1.0_wp
    case(3)
        f = 11.0_wp * x**11 - 1.0_wp
    case(4)
        f = exp(x**2 + 7.0_wp*x - 30.0_wp) - 1.0_wp
    case(5)
        f = 1.0_wp/x - sin(x) + 1.0_wp
    case(6)
        f = x**3 - 2.0_wp*x - 5.0_wp
    case(7)
        f = 1.0_wp/x - 1.0_wp
    case(8)
        f = x**3 - 5.0_wp*x**2 + 12.0_wp*x - 8.0_wp
    case(9)
        f = x**5 - 61.0_wp*x**4 + 1368.0_wp*x**3 - 13548.0_wp*x**2 + 54256.0_wp*x - 42016.0_wp ! x in [0, 10]
    case(10)
        f = x**5-61.0_wp*x**4+6801.0_wp/5_wp*x**3-66531.0_wp/5_wp*x**2+5205601.0_wp/100.0_wp*x-4005001.0_wp/100.0_wp
    case(11)
        f = (x-1.0_wp)**5
    case(12)
        f = atan(x)
    case(13)
        f = exp(x) - 2.0_wp*x - 1.0_wp
    case(14)
        f = 3.0_wp * x*sin(x*20.0_wp)*cos(x) + x - 2.0_wp
    case(15)
        f = 3.0_wp * x**2*sin(x*20.0_wp) + cos(x*2.0_wp)
    case(16)
        ! case 11 from "Algorithm 748: Enclosing Zeros of Continuous Functions" (n=20)
        n = 20
        f = (real(n,wp)*x-1.0_wp)/((real(n,wp)-1.0_wp)*x)
    case(17)
        f = sin(x*10.0_wp) + sin(x*2.0_wp) + atan(x/3.0_wp)*x**2 + exp(x) - 1.0_wp
    case default
        error stop 'Invalid selection'
    end select

    end function test_func
!**********************************************************

!*****************************************************************************************
    end program root_tests
!*****************************************************************************************