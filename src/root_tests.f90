!*****************************************************************************************
!>
!  Tests for scalar root functions without derivatives.
!
!### Reference
!  * G.E. Alefeld, F.A. Potra and Yixun Shi,
!    "Algorithm 748: Enclosing Zeros of Continuous Functions",
!    ACM Transactions on Mathematica1 Software,
!    Vol. 21. No. 3. September 1995. Pages 327-344.

program root_tests

    use stdlib_root_module
    use iso_fortran_env, only: wp => real64, output_unit
    use face, only: colorize

    implicit none

    integer :: nprob
    integer :: n
    integer :: imeth !! method counter
    integer :: ifunc !! number of function evaluations
    real(wp) :: root !! known root value
    integer :: iunit,istat,iunit_failed
    real(wp) :: ax,bx
    integer :: i !! counter

    character(len=*),parameter :: fmt  = '(A20,1X,A3,1X,A4,1X,A25,   1X,A25,   1X,A25,  1X,A5,1X,A5)' !! format for header
    character(len=*),parameter :: dfmt = '(A20,1X,I3,1X,I4,1X,E25.10,1X,E25.10,1X,E25.6,1X,I5,1X,I5)' !! format for results

    integer,parameter :: number_of_methods = 14 !! number of methods to test
    character(len=100),dimension(number_of_methods),parameter :: methods = [ &
        'brent               ', &
        'brentq              ', &
        'brenth              ', &
        'bisection           ', &
        'regula_falsi        ', &
        'anderson_bjorck     ', &
        'anderson_bjorck_king', &
        'ridders             ', &
        'pegasus             ', &
        'bdqrf               ', &
        'muller              ', &
        'chandrupatla        ', &
        'toms748             ', &
        'zhang               '] !! method names

    integer,dimension(number_of_methods) :: number_of_wins, ivec, number_of_failures, ivec2

    number_of_wins = 0
    number_of_failures = 0

    open(newunit=iunit,        file='root_report_best.txt',     status='REPLACE', iostat=istat)
    open(newunit=iunit_failed, file='root_report_failures.txt', status='REPLACE', iostat=istat)

    write(iunit, '(A5,1X,A5,1X,A5,1X,A)')  'nprob', 'n', 'evals', 'Best Methods'
    write(iunit_failed, '(A5,1X,A5,1X,A)') 'nprob', 'n', 'Failed Methods'

    ! run the tests:
    nprob = 101; n = 1;   call test()
    nprob = 102; n = 1;   call test()
    nprob = 103; n = 1;   call test()
    nprob = 104; n = 1;   call test()
    nprob = 105; n = 1;   call test()
    nprob = 106; n = 1;   call test()
    nprob = 107; n = 1;   call test()
    nprob = 108; n = 1;   call test()
    nprob = 109; n = 1;   call test()
    nprob = 110; n = 1;   call test()
    nprob = 111; n = 1;   call test()
    nprob = 112; n = 1;   call test()
    nprob = 113; n = 1;   call test()
    nprob = 114; n = 1;   call test()
    nprob = 115; n = 1;   call test()
    nprob = 116; n = 20;  call test()
    nprob = 117; n = 1;   call test()
    nprob = 1;  n = 1;    call test()
    nprob = 2;  n = 1;    call test()
    nprob = 3;  n = 1;    call test()
    nprob = 4;  n = 1;    call test()
    nprob = 5;  n = 1;    call test()
    nprob = 6;  n = 1;    call test()
    nprob = 7;  n = 1;    call test()
    nprob = 8;  n = 1;    call test()
    nprob = 9;  n = 1;    call test()
    nprob = 10; n = 1;    call test()
    nprob = 11; n = 1;    call test()
    nprob = 12; n = 1;    call test()
    nprob = 13; n = 1;    call test()
    nprob = 14; n = 1;    call test()
    nprob = 15; n = 4;    call test()
    nprob = 15; n = 6;    call test()
    nprob = 15; n = 8;    call test()
    nprob = 15; n = 10;   call test()
    nprob = 15; n = 12;   call test()
    nprob = 16; n = 4;    call test()
    nprob = 16; n = 6;    call test()
    nprob = 16; n = 8;    call test()
    nprob = 16; n = 10;   call test()
    nprob = 16; n = 12;   call test()
    nprob = 17; n = 8;    call test()
    nprob = 17; n = 10;   call test()
    nprob = 17; n = 12;   call test()
    nprob = 17; n = 14;   call test()
    nprob = 18; n = 1;    call test()
    nprob = 19; n = 1;    call test()
    nprob = 19; n = 2;    call test()
    nprob = 19; n = 3;    call test()
    nprob = 19; n = 4;    call test()
    nprob = 19; n = 5;    call test()
    nprob = 19; n = 20;   call test()
    nprob = 19; n = 40;   call test()
    nprob = 19; n = 60;   call test()
    nprob = 19; n = 80;   call test()
    nprob = 19; n = 100;  call test()
    nprob = 20; n = 5;    call test()
    nprob = 20; n = 10;   call test()
    nprob = 20; n = 20;   call test()
    nprob = 21; n = 2;    call test()
    nprob = 21; n = 5;    call test()
    nprob = 21; n = 10;   call test()
    nprob = 21; n = 15;   call test()
    nprob = 21; n = 20;   call test()
    nprob = 22; n = 1;    call test()
    nprob = 22; n = 2;    call test()
    nprob = 22; n = 4;    call test()
    nprob = 22; n = 5;    call test()
    nprob = 22; n = 8;    call test()
    nprob = 22; n = 15;   call test()
    nprob = 22; n = 20;   call test()
    nprob = 23; n = 1;    call test()
    nprob = 23; n = 5;    call test()
    nprob = 23; n = 10;   call test()
    nprob = 23; n = 15;   call test()
    nprob = 23; n = 20;   call test()
    nprob = 24; n = 2;    call test()
    nprob = 24; n = 5;    call test()
    nprob = 24; n = 15;   call test()
    nprob = 24; n = 20;   call test()
    nprob = 25; n = 2;    call test()
    nprob = 25; n = 3;    call test()
    nprob = 25; n = 4;    call test()
    nprob = 25; n = 5;    call test()
    nprob = 25; n = 6;    call test()
    nprob = 25; n = 7;    call test()
    nprob = 25; n = 9;    call test()
    nprob = 25; n = 11;   call test()
    nprob = 25; n = 13;   call test()
    nprob = 25; n = 15;   call test()
    nprob = 25; n = 17;   call test()
    nprob = 25; n = 19;   call test()
    nprob = 25; n = 21;   call test()
    nprob = 25; n = 23;   call test()
    nprob = 25; n = 25;   call test()
    nprob = 25; n = 27;   call test()
    nprob = 25; n = 29;   call test()
    nprob = 25; n = 31;   call test()
    nprob = 25; n = 33;   call test()
    !nprob = 26; n = 1;    call test()  ! this is just almost impossible so commented out for now
    nprob = 27; n = 1;    call test()
    nprob = 27; n = 2;    call test()
    nprob = 27; n = 3;    call test()
    nprob = 27; n = 4;    call test()
    nprob = 27; n = 5;    call test()
    nprob = 27; n = 6;    call test()
    nprob = 27; n = 7;    call test()
    nprob = 27; n = 8;    call test()
    nprob = 27; n = 9;    call test()
    nprob = 27; n = 10;   call test()
    nprob = 27; n = 11;   call test()
    nprob = 27; n = 12;   call test()
    nprob = 27; n = 13;   call test()
    nprob = 27; n = 14;   call test()
    nprob = 27; n = 15;   call test()
    nprob = 27; n = 16;   call test()
    nprob = 27; n = 17;   call test()
    nprob = 27; n = 18;   call test()
    nprob = 27; n = 19;   call test()
    nprob = 27; n = 20;   call test()
    nprob = 27; n = 21;   call test()
    nprob = 27; n = 22;   call test()
    nprob = 27; n = 23;   call test()
    nprob = 27; n = 24;   call test()
    nprob = 27; n = 25;   call test()
    nprob = 27; n = 26;   call test()
    nprob = 27; n = 27;   call test()
    nprob = 27; n = 28;   call test()
    nprob = 27; n = 29;   call test()
    nprob = 27; n = 30;   call test()
    nprob = 27; n = 31;   call test()
    nprob = 27; n = 32;   call test()
    nprob = 27; n = 33;   call test()
    nprob = 27; n = 34;   call test()
    nprob = 27; n = 35;   call test()
    nprob = 27; n = 36;   call test()
    nprob = 27; n = 37;   call test()
    nprob = 27; n = 38;   call test()
    nprob = 27; n = 39;   call test()
    nprob = 27; n = 40;   call test()
    nprob = 28; n = 20;   call test()
    nprob = 28; n = 21;   call test()
    nprob = 28; n = 22;   call test()
    nprob = 28; n = 23;   call test()
    nprob = 28; n = 24;   call test()
    nprob = 28; n = 25;   call test()
    nprob = 28; n = 26;   call test()
    nprob = 28; n = 27;   call test()
    nprob = 28; n = 28;   call test()
    nprob = 28; n = 29;   call test()
    nprob = 28; n = 30;   call test()
    nprob = 28; n = 31;   call test()
    nprob = 28; n = 32;   call test()
    nprob = 28; n = 33;   call test()
    nprob = 28; n = 34;   call test()
    nprob = 28; n = 35;   call test()
    nprob = 28; n = 36;   call test()
    nprob = 28; n = 37;   call test()
    nprob = 28; n = 38;   call test()
    nprob = 28; n = 39;   call test()
    nprob = 28; n = 40;   call test()
    nprob = 28; n = 100;  call test()
    nprob = 28; n = 200;  call test()
    nprob = 28; n = 300;  call test()
    nprob = 28; n = 400;  call test()
    nprob = 28; n = 500;  call test()
    nprob = 28; n = 600;  call test()
    nprob = 28; n = 700;  call test()
    nprob = 28; n = 800;  call test()
    nprob = 28; n = 900;  call test()
    nprob = 28; n = 1000; call test()
    nprob = 29; n = 1;    call test()

    !nprob = 202; n = 1;  call test()
    nprob = 203; n = 1;   call test()
    nprob = 204; n = 1;   call test()
    nprob = 205; n = 1;   call test()
    nprob = 206; n = 1;   call test()
    nprob = 207; n = 1;   call test()
    nprob = 208; n = 1;   call test()
    nprob = 209; n = 1;   call test()
    nprob = 210; n = 1;   call test()
    nprob = 211; n = 1;   call test()
    nprob = 212; n = 1;   call test()
    nprob = 213; n = 1;   call test()
    nprob = 214; n = 1;   call test()
    nprob = 215; n = 1;   call test()
    nprob = 216; n = 1;   call test()
    nprob = 217; n = 1;   call test()
    nprob = 218; n = 1;   call test()

    nprob = 300; n = 1;   call test()
    nprob = 301; n = 1;   call test()
    nprob = 303; n = 1;   call test()
    nprob = 304; n = 1;   call test()
    nprob = 305; n = 1;   call test()
    nprob = 306; n = 1;   call test()

    nprob = 996; n = 1;   call test()
    nprob = 997; n = 1;   call test()
    nprob = 998; n = 1;   call test()
    nprob = 999; n = 1;   call test()

    close(iunit)
    close(iunit_failed)

    ! another summary:
    ivec  = [(i, i = 1, number_of_methods)]
    ivec2 = ivec

    call insertion_sort(number_of_wins, ivec, number_of_failures)
   ! call insertion_sort(number_of_failures, ivec2)

    write(*,*) ''
    write(*,'(A25,1X,A5,1X,A5)') repeat('-',25), repeat('-',5), repeat('-',5)
    write(*,'(A25,1X,A5,1X,A5)') 'Method', 'Win', 'Fail'
    write(*,'(A25,1X,A5,1X,A5)') repeat('-',25), repeat('-',5), repeat('-',5)
    do imeth = 1, number_of_methods
        write(*,'(A25,1X,I5,1X,I5)') trim(methods(ivec(imeth))), number_of_wins(imeth), number_of_failures(imeth)
    end do
    write(*,*) ''

    contains
!*****************************************************************************************

    subroutine test()

        implicit none

        real(wp) :: xzero,fzero,error
        integer :: iflag
        character(len=1000) :: line
        integer,dimension(number_of_methods) :: fevals !! func evals for each case
        integer :: best_feval
        character(len=:),allocatable :: best,failures
        logical :: root_found

        real(wp),parameter :: tol_for_check = 1.0e-7_wp

        write(output_unit,fmt) &
            repeat('-',20),repeat('-',3),repeat('-',4),repeat('-',25),&
            repeat('-',25),repeat('-',25),repeat('-',5),repeat('-',5)

        write(output_unit,fmt) &
            'method','function','n','error','x','f','evals','iflag'

        write(output_unit,fmt) &
            repeat('-',20),repeat('-',3),repeat('-',4),repeat('-',25),&
            repeat('-',25),repeat('-',25),repeat('-',5),repeat('-',5)

        do imeth = 1, number_of_methods

            ifunc = 0 ! reset func evals counter
            call problems(ax=ax, bx=bx, xroot=root)
            call root_scalar(methods(imeth),func,ax,bx,xzero,fzero,iflag,&
                             atol = 1.0e-15_wp, rtol = 1.0e-13_wp, ftol = 1.0e-15_wp, maxiter = 1000)

            error = xzero-root
            write(line, dfmt) trim(methods(imeth)),nprob,n,error,xzero,fzero,ifunc,iflag

            root_found = abs(fzero) <= tol_for_check !.and. abs(error) <= tol_for_check

            if (root_found) then
                write(output_unit, '(A)') trim(line)
            else
                write(output_unit, '(A)') colorize(trim(line), color_fg='red') ! failed case
            end if

            ! save results for this case:
            if (iflag==0 .and. root_found) then
                fevals(imeth) = ifunc
            else
                fevals(imeth) = huge(1)
            end if

        end do

        ! get the best and the failures for the output files:
        best_feval = minval(fevals)
        if (best_feval==huge(1)) best_feval = -1 ! if non of them converged
        best = ''
        failures = ''
        do imeth = 1, number_of_methods
            if (fevals(imeth) == best_feval) then
                if (best=='') then
                    best = trim(methods(imeth))
                else
                    best = best//', '//trim(methods(imeth))
                end if
                number_of_wins(imeth) = number_of_wins(imeth) + 1
            elseif (fevals(imeth) == huge(1)) then
                if (failures=='') then
                    failures = trim(methods(imeth))
                else
                    failures = failures//', '//trim(methods(imeth))
                end if
                number_of_failures(imeth) = number_of_failures(imeth) + 1
            end if
        end do

        ! output the report:
        if (best /= '')     write(iunit, '(I5,1X,I5,1X,I5,1X,A)')  nprob, n, best_feval, 'Best: '//trim(best)
        if (failures /= '') write(iunit_failed, '(I5,1X,I5,1X,A)') nprob, n, 'Failures: '//trim(failures)

    end subroutine test

    subroutine problems(x, ax, bx, fx, xroot)

    implicit none

    real(wp),intent(in),optional  :: x       !! indep variable
    real(wp),intent(out),optional :: ax, bx  !! bounds
    real(wp),intent(out),optional :: fx      !! function value
    real(wp),intent(out),optional :: xroot   !! value of indep variable at the root

    real(wp),parameter :: pi = acos(-1.0_wp)

    real(wp) :: a, b, root, f, dn, di, xi, t1, emx, ex
    integer :: i !! counter

    dn = real(n, wp)

    select case (nprob)
    case (1)
        a = pi/2.0_wp
        b = pi
        root = 1.8954942670340_wp
        if (present(x)) f = sin(x) - x/2.0_wp
    case (2)
        a = 1.0_wp + 1.0e-9_wp
        b = (2.0_wp)*(2.0_wp) - 1.0e-9_wp
        root = 3.0229153472731_wp
        if (present(x)) then
            f = 0.0_wp
            do i = 1, 20
                di = real(i,wp)
                f = f + ((2.0_wp*di - 5.0_wp)**2)/(x - di*di)**3
            end do
            f = -2.0_wp*f
        end if
    case (3)
        a = (2.0_wp)*(2.0_wp) + 1.0e-9_wp
        b = (3.0_wp)*(3.0_wp) - 1.0e-9_wp
        root = 6.6837535608081_wp
        if (present(x)) then
            f = 0.0_wp
            do i = 1, 20
                di = real(i,wp)
                f = f + ((2.0_wp*di - 5.0_wp)**2)/(x - di*di)**3
            end do
            f = -2.0_wp*f
        end if
    case (4)
        a = (3.0_wp)*(3.0_wp) + 1.0e-9_wp
        b = (4.0_wp)*(4.0_wp) - 1.0e-9_wp
        root = 11.238701655002_wp
        if (present(x)) then
            f = 0.0_wp
            do i = 1, 20
                di = real(i,wp)
                f = f + ((2.0_wp*di - 5.0_wp)**2)/(x - di*di)**3
            end do
            f = -2.0_wp*f
        end if
    case (5)
        a = (4.0_wp)*(4.0_wp) + 1.0e-9_wp
        b = (5.0_wp)*(5.0_wp) - 1.0e-9_wp
        root = 19.676000080623_wp
        if (present(x)) then
            f = 0.0_wp
            do i = 1, 20
                di = real(i,wp)
                f = f + ((2.0_wp*di - 5.0_wp)**2)/(x - di*di)**3
            end do
            f = -2.0_wp*f
        end if
    case (6)
        a = (5.0_wp)*(5.0_wp) + 1.0e-9_wp
        b = (6.0_wp)*(6.0_wp) - 1.0e-9_wp
        root = 29.828227326505_wp
        if (present(x)) then
            f = 0.0_wp
            do i = 1, 20
                di = real(i,wp)
                f = f + ((2.0_wp*di - 5.0_wp)**2)/(x - di*di)**3
            end do
            f = -2.0_wp*f
        end if
    case (7)
        a = (6.0_wp)*(6.0_wp) + 1.0e-9_wp
        b = (7.0_wp)*(7.0_wp) - 1.0e-9_wp
        root = 41.906116195289_wp
        if (present(x)) then
            f = 0.0_wp
            do i = 1, 20
                di = real(i,wp)
                f = f + ((2.0_wp*di - 5.0_wp)**2)/(x - di*di)**3
            end do
            f = -2.0_wp*f
        end if
    case (8)
        a = (7.0_wp)*(7.0_wp) + 1.0e-9_wp
        b = (8.0_wp)*(8.0_wp) - 1.0e-9_wp
        root = 55.953595800143_wp
        if (present(x)) then
            f = 0.0_wp
            do i = 1, 20
                di = real(i,wp)
                f = f + ((2.0_wp*di - 5.0_wp)**2)/(x - di*di)**3
            end do
            f = -2.0_wp*f
        end if
    case (9)
        a = (8.0_wp)*(8.0_wp) + 1.0e-9_wp
        b = (9.0_wp)*(9.0_wp) - 1.0e-9_wp
        root = 71.985665586588_wp
        if (present(x)) then
            f = 0.0_wp
            do i = 1, 20
                di = real(i,wp)
                f = f + ((2.0_wp*di - 5.0_wp)**2)/(x - di*di)**3
            end do
            f = -2.0_wp*f
        end if
    case (10)
        a = (9.0_wp)*(9.0_wp) + 1.0e-9_wp
        b = (10.0_wp)*(10.0_wp) - 1.0e-9_wp
        root = 90.008868539167_wp
        if (present(x)) then
            f = 0.0_wp
            do i = 1, 20
                di = real(i,wp)
                f = f + ((2.0_wp*di - 5.0_wp)**2)/(x - di*di)**3
            end do
            f = -2.0_wp*f
        end if
    case (11)
        a = (10.0_wp)*(10.0_wp) + 1.0e-9_wp
        b = (11.0_wp)*(11.0_wp) - 1.0e-9_wp
        root = 110.02653274833_wp
        if (present(x)) then
            f = 0.0_wp
            do i = 1, 20
                di = real(i,wp)
                f = f + ((2.0_wp*di - 5.0_wp)**2)/(x - di*di)**3
            end do
            f = -2.0_wp*f
        end if
    case (12)
        a = -9.0_wp
        b = 31.0_wp
        root = 0.0_wp
        if (present(x)) f = -40.0_wp*x*exp(-1.0_wp*x)
    case (13)
        a = -9.0_wp
        b = 31.0_wp
        root = 0.0_wp
        if (present(x)) f = -100.0_wp*x*exp(-2.0_wp*x)
    case (14)
        a = -9.0_wp
        b = 31.0_wp
        root = 0.0_wp
        if (present(x)) f = -200.0_wp*x*exp(-3.0_wp*x)
    case (15)
        a = 0.0_wp
        b = 5.0_wp
        select case (n)
        case(4);  root = 0.66874030497642_wp
        case(6);  root = 0.76472449133173_wp
        case(8);  root = 0.81776543395794_wp
        case(10); root = 0.85133992252078_wp
        case(12); root = 0.87448527222117_wp
        case default
            error stop 'invalid n'
        end select
        if (present(x)) f = x**n - 0.2_wp
    case (16)
        a = 0.0_wp
        b = 5.0_wp
        select case (n)
        case(4);    root = 1.0_wp
        case(6);    root = 1.0_wp
        case(8);    root = 1.0_wp
        case(10);   root = 1.0_wp
        case(12);   root = 1.0_wp
        case default; error stop 'invalid n'
        end select
        if (present(x)) f = x**n - 1.0_wp
    case (17)
        a = -0.95_wp
        b = 4.05_wp
        root = 1.0_wp
        if (present(x)) f = x**n - 1.0_wp
    case (18)
        a = 0.0_wp
        b = 1.5_wp
        root = 0.52359877559830_wp
        if (present(x)) f = sin(x) - 0.5_wp
    case (19)
        a = 0.0_wp
        b = 1.0_wp
        select case (n)
        case(1);    root = 0.42247770964124_wp
        case(2);    root = 0.30669941048320_wp
        case(3);    root = 0.22370545765466_wp
        case(4);    root = 0.17171914751951_wp
        case(5);    root = 0.13825715505682_wp
        case(20);   root = 3.4657359020854e-02_wp
        case(40);   root = 1.7328679513999e-02_wp
        case(60);   root = 1.1552453009332e-02_wp
        case(80);   root = 8.6643397569993e-03_wp
        case(100);  root = 6.9314718055995e-03_wp
        case default; error stop 'invalid n'
        end select
        if (present(x)) f = 2.0_wp*x*exp(-dn) - 2.0_wp*exp(-dn*x) + 1.0_wp
    case (20)
        a = 0.0_wp
        b = 1.0_wp
        select case (n)
        case(5);    root = 3.8402551840622e-02_wp
        case(10);   root = 9.9000099980005e-03_wp
        case(20);   root = 2.4937500390620e-03_wp
        case default; error stop 'invalid n'
        end select
        if (present(x)) f = (1.0_wp + (1.0_wp - dn)**2)*x - (1.0_wp - dn*x)**2
    case (21)
        a = 0.0_wp
        b = 1.0_wp
        select case (n)
        case(2);    root = 0.5_wp
        case(5);    root = 0.34595481584824_wp
        case(10);   root = 0.24512233375331_wp
        case(15);   root = 0.19554762353657_wp
        case(20);   root = 0.16492095727644_wp
        case default; error stop 'invalid n'
        end select
        if (present(x)) f = x**2 - (1.0_wp - x)**n
    case (22)
        a = 0.0_wp
        b = 1.0_wp
        select case (n)
        case(1);    root = 0.27550804099948_wp
        case(2);    root = 0.13775402049974_wp
        case(4);    root = 1.0305283778156e-02_wp
        case(5);    root = 3.6171081789041e-03_wp
        case(8);    root = 4.1087291849640e-04_wp
        case(15);   root = 2.5989575892908e-05_wp
        case(20);   root = 7.6685951221853e-06_wp
        case default; error stop 'invalid n'
        end select
        if (present(x)) f = (1.0_wp + (1.0_wp - dn)**4)*x - (1.0_wp - dn*x)**4
    case (23)
        a = 0.0_wp
        b = 1.0_wp
        select case (n)
        case(1);    root = 0.40105813754155_wp
        case(5);    root = 0.51615351875793_wp
        case(10);   root = 0.53952222690842_wp
        case(15);   root = 0.54818229434066_wp
        case(20);   root = 0.55270466667849_wp
        case default; error stop 'invalid n'
        end select
        if (present(x)) f = (x - 1.0_wp)*exp(-dn*x) + x**n
    case (24)
        a = 1.0e-2_wp
        b = 1.0_wp
        select case (n)
        case(2);    root = 0.5_wp
        case(5);    root = 0.2_wp
        case(15);   root = 6.6666666666667e-02_wp
        case(20);   root = 5.0e-02_wp
        case default; error stop 'invalid n'
        end select
        if (present(x)) f = (dn*x - 1.0_wp)/((dn - 1.0_wp)*x)
    case (25)
        a = 1.0_wp
        b = 100.0_wp
        root = real(n,wp)
        if (present(x)) f = x**(1.0_wp/dn) - dn**(1.0_wp/dn)
    case (26)
        a = -1.0_wp
        b = 4.0_wp
        root = 0.0_wp    ! this is just almost impossible so commented out for now
        if (present(x)) then
            if (x == 0.0_wp) then
                f = 0.0_wp
            else
                f = x/exp(1.0_wp/(x*x))
            end if
        end if
    case (27)
        a = -10000.0_wp
        b = pi/2.0_wp
        root = 0.62380651896161_wp
        if (present(x)) then
            if (x >= 0.0_wp) then
                f = (x/1.5_wp + sin(x) - 1.0_wp)*dn/20.0_wp
            else
                f = (-1.0_wp*dn)/20.0_wp
            end if
        end if
    case (28)
        a = -10000.0_wp
        b = 1.0e-4_wp
        select case (n)
        case(20);   root = 5.9051305594220e-05_wp
        case(21);   root = 5.6367155339937e-05_wp
        case(22);   root = 5.3916409455592e-05_wp
        case(23);   root = 5.1669892394942e-05_wp
        case(24);   root = 4.9603096699145e-05_wp
        case(25);   root = 4.7695285287639e-05_wp
        case(26);   root = 4.5928793239949e-05_wp
        case(27);   root = 4.4288479195665e-05_wp
        case(28);   root = 4.2761290257883e-05_wp
        case(29);   root = 4.1335913915954e-05_wp
        case(30);   root = 4.0002497338020e-05_wp
        case(31);   root = 3.8752419296207e-05_wp
        case(32);   root = 3.7578103559958e-05_wp
        case(33);   root = 3.6472865219959e-05_wp
        case(34);   root = 3.5430783356532e-05_wp
        case(35);   root = 3.4446594929961e-05_wp
        case(36);   root = 3.3515605877800e-05_wp
        case(37);   root = 3.2633616249437e-05_wp
        case(38);   root = 3.1796856858426e-05_wp
        case(39);   root = 3.1001935436965e-05_wp
        case(40);   root = 3.0245790670210e-05_wp
        case(100);  root = 1.2277994232462e-05_wp
        case(200);  root = 6.1695393904409e-06_wp
        case(300);  root = 4.1198585298293e-06_wp
        case(400);  root = 3.0924623877272e-06_wp
        case(500);  root = 2.4752044261050e-06_wp
        case(600);  root = 2.0633567678513e-06_wp
        case(700);  root = 1.7690120078154e-06_wp
        case(800);  root = 1.5481615698859e-06_wp
        case(900);  root = 1.3763345366022e-06_wp
        case(1000); root = 1.2388385788997e-06_wp
        case default; error stop 'invalid n'
        end select
        if (present(x)) then
            if (x >= (1.0e-3_wp)*2.0_wp/(dn + 1.0_wp)) then
                f = exp(1.0_wp) - 1.859_wp
            elseif (x >= 0.0_wp) then
                f = exp((dn + 1.0_wp)*0.5_wp*x*1000.0_wp) - 1.859_wp
            else
                f = -0.859_wp
            end if
        end if
    case (29)
        ! Zhang test case
        a = 0.0_wp
        b = 4.0_wp
        root = 0.8654740331e+00_wp
        if (present(x)) f = cos(x) - x**3

    ! 101-107 : Gottlieb's paper [table 1 & 2]
    case (101)
        a = 0.0_wp
        b = 1.0_wp
        root = asin(2.0_wp/3.0_wp)
        if (present(x)) f = 3.0_wp * sin(x) - 2.0_wp
    case (102)
        a = -1.0_wp
        b = 1.0_wp
        root = 0.5671432904_wp
        if (present(x)) f = x*exp(x) - 1.0_wp
    case (103)
        a = 0.1_wp
        b = 0.9_wp
        root = 0.8041330975_wp
        if (present(x)) f = 11.0_wp * x**11 - 1.0_wp
    case (104)
        a = 2.8_wp
        b = 3.1_wp
        root = 3.0_wp
        if (present(x)) f = exp(x**2 + 7.0_wp*x - 30.0_wp) - 1.0_wp
    case (105)
        a = -1.3_wp
        b = -0.5_wp
        root = -0.6294464841_wp
        if (present(x)) f = 1.0_wp/x - sin(x) + 1.0_wp
    case (106)
        a = 2.0_wp
        b = 3.0_wp
        root = 2.0945514815_wp
        if (present(x)) f = x**3 - 2.0_wp*x - 5.0_wp
    case (107)
        a = 0.5_wp
        b = 2.0_wp
        root = 1.0_wp
        if (present(x)) f = 1.0_wp/x - 1.0_wp

    case (108)
        a = -10.0_wp
        b = 10.0_wp
        root = 1.0_wp
        if (present(x)) f = x**3 - 5.0_wp*x**2 + 12.0_wp*x - 8.0_wp
    case (109)
        a = 0.0_wp
        b = 10.0_wp
        root = 1.0_wp
        if (present(x)) f = x**5 - 61.0_wp*x**4 + 1368.0_wp*x**3 - &
                            13548.0_wp*x**2 + 54256.0_wp*x - 42016.0_wp ! x in [0, 10]
    case (110)
        a = 0.0_wp
        b = 10.0_wp
        root = 1.0_wp
        if (present(x)) f = x**5-61.0_wp*x**4+6801.0_wp/5.0_wp*x**3-66531.0_wp/&
                            5.0_wp*x**2+5205601.0_wp/100.0_wp*x-4005001.0_wp/100.0_wp
    case (111)
        a = 0.0_wp
        b = 10.0_wp
        root = 1.0_wp
        if (present(x)) f = (x-1.0_wp)**5
    case (112)
        !From: http://www.mth.pdx.edu/~daescu/mth451_551/Muller_bisection.pdf
        a = -2.0_wp
        b = 3.0_wp
        root = 0.0_wp
        if (present(x)) f = atan(x)
    case (113)
        !From: http://www.mth.pdx.edu/~daescu/mth451_551/Muller_bisection.pdf
        a = 0.0_wp     ! x1 is the location of the root for this case
        b = 3.0_wp
        root = 0.0_wp
        if (present(x)) f = exp(x) - 2.0_wp*x - 1.0_wp
    case (114)
        a = 0.0_wp
        b = 0.7_wp
        root = 0.677712926329_wp
        if (present(x)) f = 3.0_wp * x*sin(x*20.0_wp)*cos(x) + x - 2.0_wp
    case (115)
        a = -0.7_wp
        b = 0.4_wp
        root = -0.640083696085_wp
        if (present(x)) f = 3.0_wp * x**2*sin(x*20.0_wp) + cos(x*2.0_wp)
    case (116)
        !case 11 from "Algorithm 748: Enclosing Zeros of Continuous Functions" (n=20)
        a = 0.01_wp
        b = 1.0_wp
        root = 0.05_wp
        if (present(x)) f = (real(n,wp)*x-1.0_wp)/((real(n,wp)-1.0_wp)*x)
    case (117)
        a = -10.0_wp
        b = 5.0_wp
        root = 0.0_wp
        if (present(x)) f = sin(x*10.0_wp) + sin(x*2.0_wp) + &
                            atan(x/3.0_wp)*x**2 + exp(x) - 1.0_wp

    ! 2xx functions are from:
    !   T. R. Chandrupatla, "A new hybrid quadratic/bisection algorithm for finding the zero of
    !   a nonlinear function without using derivatives", Advances in Engineering Software
    !   Volume 28, Issue 3, April 1997, Pages 145-149

    !case (202)
    !   a  = 2.0_wp
    !   b  = 3.0_wp
    !   root = 2.0945515_wp
    !   if (present(x)) f = x**3 - 2.0_wp*x - 5.0_wp  ! same as #106 above
    case (203)
        a = 0.5_wp
        b = 1.51_wp
        root = 1.0_wp
        if (present(x)) f = 1.0_wp - 1.0_wp/(x**2)
    case (204)
        a =  1.0e-12_wp
        b =  1.0e12_wp
        root = 1.0_wp
        if (present(x)) f = 1.0_wp - 1.0_wp/(x**2)
    case (205)
        a = 0.0_wp
        b = 5.0_wp
        root = 3.0_wp
        if (present(x)) f = (x-3.0_wp)**3
    case (206)
        a = -1.0e10_wp
        b =  1.0e10_wp
        root = 3.0_wp
        if (present(x)) f = (x-3.0_wp)**3
    case (207)
        a = 0.0_wp
        b = 5.0_wp
        root = 2.0_wp
        if (present(x)) f = 6.0_wp*(x-2.0_wp)**5
    case (208)
        a = -1.0e10_wp
        b =  1.0e10_wp
        root = 2.0_wp
        if (present(x)) f = 6.0_wp*(x-2.0_wp)**5
    case (209)
        a = -1.0_wp
        b =  4.0_wp
        root = 0.0_wp
        if (present(x)) f = x**9
    case (210)
        a = -10.0_wp
        b =  100.0_wp
        root = 0.0_wp
        if (present(x)) f = x**9
    case (211)
        a = -1.0_wp
        b =  4.0_wp
        root = 0.0_wp
        if (present(x)) f = x**19
    case (212)
        a = -10.0_wp
        b =  100.0_wp
        root = 0.0_wp
        if (present(x)) f = x**19
    case (213)
        a = -1.0_wp
        b =  4.0_wp
        root = 0.0_wp
        if (present(x)) then
            if (abs(x) < 3.8e-4_wp) then  ! same as #26 above except for this
            f = 0.0_wp
            else
            f = x*exp(-x**(-2))
            endif
        end if
    case (214)
        a = -10.0_wp
        b =  100.0_wp
        root = 0.0_wp
        if (present(x)) then
            if (abs(x) < 3.8e-4_wp) then  ! same as #26 above except for this
            f = 0.0_wp
            else
            f = x*exp(-x**(-2))
            endif
        end if
    case (215)
        a = 2.0e-4_wp
        b = 2.0_wp
        root = 1.037536033_wp
        if (present(x)) then
            xi = 0.61489_wp
            t1 = 1.0_wp-xi
            emx = exp(-x)
            f = -(3062.0_wp*t1*emx)/(xi + t1*emx) - 1013.0_wp + 1628.0_wp/x
        end if
    case (216)
        a = 2.0e-4_wp
        b = 81.0_wp
        root = 1.037536033_wp
        if (present(x)) then
            xi = 0.61489_wp
            t1 = 1.0_wp-xi
            emx = exp(-x)
            f = -(3062.0_wp*t1*emx)/(xi + t1*emx) - 1013.0_wp + 1628.0_wp/x
        end if
    case (217)
        a = 2.0e-4_wp
        b = 1.0_wp
        root = 0.7032048404_wp
        if (present(x)) then
            ex = exp(x)
            f = ex - 2.0_wp - 0.01_wp/(x*x) + 2.0e-6_wp/(x*x*x)
        end if
    case (218)
        a = 2.0e-4_wp
        b = 81.0_wp
        root = 0.7032048404_wp
        if (present(x)) then
            ex = exp(x)
            f = ex - 2.0_wp - 0.01_wp/(x*x) + 2.0e-6_wp/(x*x*x)
        end if

    ! from A modified three-point Secant method with improved
    !      rate and characteristics of convergence
    !      July 2019

    case(300)
        a = 0.0_wp
        b = 1.5_wp
        root = 1.365230013414100_wp
        if (present(x)) f = x**3 + 4.0_wp*x**2 - 10.0_wp
    case(301)
        a = -2.0_wp
        b = -1.0_wp
        root = -1.404491648215340_wp
        if (present(x)) f = sin(x)**2 - x**2 + 1.0_wp
    case(303)
        a = 1.0_wp
        b = 2.1_wp
        root = 2.0_wp
        if (present(x)) f = (x - 1.0_wp)**6 - 1.0_wp
    case(304)
        a = -0.9_wp
        b = -0.1_wp
        root = -0.603231971557215_wp
        if (present(x)) f = sin(x)*exp(x) + log(x**2+1.0_wp)
    case(305)
        a = 1.0_wp
        b = 3.1_wp
        root = 3.0_wp
        if (present(x)) f = exp(x**2+7.0_wp*x-30.0_wp) - 1.0_wp
    case(306)
        a = 1.0_wp
        b = 2.1_wp
        root = 1.857183860207840_wp
        if (present(x)) f = x - 3.0_wp * log(x)

    ! some simple functions
    case(996)
        a = 0.1001_wp
        b = 10.0001_wp
        root = 1.0_wp
        if (present(x)) f = x**2 - 1.0_wp
    case(997)
        a = 0.0001_wp
        b = 10.0001_wp
        root = 0.1001861840E+01_wp
        if (present(x)) f = x**2 - 1.0038273894_wp + sin(x/10000.0_wp)
    case(998)
        a = -2.0001_wp
        b = 1.0001_wp
        root = 0.0_wp
        if (present(x)) f = x
    case(999)
        a = -0.0001_wp  ! root very close to a
        b = 1.0001_wp
        root = 0.0_wp
        if (present(x)) f = x

    case default
        error stop 'invalid case'
    end select

    ! outputs:
    if (present(ax)) ax = a
    if (present(bx)) bx = b
    if (present(xroot)) xroot = root
    if (present(fx)) fx = f

    end subroutine problems

    function func(x) result(f)

        implicit none

        real(wp),intent(in) :: x
        real(wp) :: f

        call problems(x=x, fx=f)
        ifunc = ifunc + 1

    end function func

    pure subroutine insertion_sort(vec, vec2, vec3)

    implicit none

    integer,dimension(:),intent(inout) :: vec  !! sort ascending
    integer,dimension(:),intent(inout) :: vec2 !! carry this one along
    integer,dimension(:),intent(inout) :: vec3 !! carry this one along

    integer :: i,j,n

    n = size(vec)

    do i = 1+1, n
        do j = i, 1+1, -1
            if ( vec(j) < vec(j-1) ) then
                call swap(vec(j), vec(j-1))
                call swap(vec2(j),vec2(j-1))
                call swap(vec3(j),vec3(j-1))
            else
                exit
            end if
        end do
    end do

    end subroutine insertion_sort

    pure elemental subroutine swap(a,b)

    implicit none

    integer,intent(inout) :: a
    integer,intent(inout) :: b

    integer :: tmp

    tmp = a
    a   = b
    b   = tmp

    end subroutine swap

!*****************************************************************************************
    end program root_tests
!*****************************************************************************************