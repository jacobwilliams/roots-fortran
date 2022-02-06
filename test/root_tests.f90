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

    use iso_fortran_env
    use root_module, wp => root_module_rk
    use face,            only: colorize
    use pyplot_module,   only: pyplot

    implicit none

    integer :: nprob
    integer :: n
    integer :: imeth !! method counter
    integer :: ifunc !! number of function evaluations
    real(wp) :: root !! known root value
    integer :: iunit,istat,iunit_failed
    real(wp) :: ax,bx
    integer :: i !! counter
    integer :: ic !! case counter
    integer :: num_of_problems
    integer,dimension(:),allocatable :: cases_to_run
    type(pyplot) :: stats_plot   !! pyplot handler

    character(len=*),parameter :: fmt  = '(   A20,1X,A3,1X,A4,1X,A16,  1X,A25,   1X,A16,  1X,A5,1X,A5,1X,A8  )' !! format for header
    character(len=*),parameter :: dfmt = '(1P,A20,1X,I3,1X,I4,1X,E16.6,1X,E25.16,1X,E16.6,1X,I5,1X,I5,1X,E8.1)' !! format for results

    integer,parameter :: number_of_methods = 16 !! number of methods to test
    character(len=100),dimension(number_of_methods),parameter :: methods = [ &
        'bisection           ', &
        'brent               ', &
        'brentq              ', &
        'brenth              ', &
        'regula_falsi        ', &
        'illinois            ', &
        'anderson_bjorck     ', &
        'anderson_bjorck_king', &
        'ridders             ', &
        'pegasus             ', &
        'bdqrf               ', &
        'muller              ', &
        'chandrupatla        ', &
        'toms748             ', &
        'zhang               ', &
        'blendtf             '] !! method names

    integer,dimension(number_of_methods) :: number_of_wins, ivec, number_of_failures, ivec2

    number_of_wins = 0
    number_of_failures = 0

    ! open output files:
    open(newunit=iunit,        file='root_report_best.txt',     status='REPLACE', iostat=istat)
    open(newunit=iunit_failed, file='root_report_failures.txt', status='REPLACE', iostat=istat)

    ! output file headers:
    write(iunit, '(A5,1X,A5,1X,A5,1X,A)')  'nprob', 'n', 'evals', 'Best Methods'
    write(iunit_failed, '(A5,1X,A5,1X,A)') 'nprob', 'n', 'Failed Methods'

    ! plot:
    call stats_plot%initialize(grid=.true.,xlabel='Methods',ylabel='Function Evaluations',figsize=[10,5],&
                                title='Performance of each method',&
                                tight_layout=.true.)

    ! number of problems to test:
    nprob = 1 ! initialize
    n = 1
    call problems(num_of_problems=num_of_problems)

    ! run the tests:
    do nprob = 1, num_of_problems
        ! write(*,*) 'nprob: ', nprob
        ! write(*,*) 'cases_to_run: ', cases_to_run
        call problems(cases = cases_to_run)
        do ic = 1, size(cases_to_run)
            n = cases_to_run(ic)
            call test()
        end do
    end do

    ! close output files:
    close(iunit)
    close(iunit_failed)

    call stats_plot%savefig('stats_plot.pdf',istat=istat,pyfile='stats_plot.py')

    ! another summary:
    ivec  = [(i, i = 1, number_of_methods)]
    ivec2 = ivec

    ! sort them based on which had the most wins:
    call insertion_sort(number_of_wins, ivec, number_of_failures)

    write(*,*) ''
    write(*,'(A25,1X,A5,1X,A5)') repeat('-',25), repeat('-',5), repeat('-',5)
    write(*,'(A25,1X,A5,1X,A5)') 'Method', 'Win', 'Fail'
    write(*,'(A25,1X,A5,1X,A5)') repeat('-',25), repeat('-',5), repeat('-',5)
    do imeth = 1, number_of_methods
        write(*,'(A25,1X,I5,1X,I5)') trim(methods(ivec(imeth))), &
                                     number_of_wins(imeth), &
                                     number_of_failures(imeth)
    end do
    write(*,*) ''

   ! call generate_plots()  ! comment out for now...

    contains
!*****************************************************************************************

    subroutine test()

        implicit none

        real(wp) :: xzero,fzero,error
        integer :: iflag
        character(len=1000) :: line
        integer,dimension(number_of_methods) :: fevals !! func evals for each case
        real(wp),dimension(number_of_methods) :: cases_per_sec !! speed of each case
        integer :: best_feval
        real(wp) :: best_cases_per_second
        character(len=:),allocatable :: best,failures
        logical :: root_found
        real(wp) :: tstart, tfinish  !! for `cpu_time`
        integer :: irepeat !! test repeat counter
        integer :: i
        real(wp) :: atol, rtol, ftol
        real(wp) :: tol_for_check !! for pass/fail check

        integer,parameter :: n_repeat = 1  !! number of times to repeat each test for timing purposes
        integer,parameter :: maxiter = 1000 !! maximum number of iterations

        select case (wp)
        case(real32)
            atol = 1.0e-5_wp
            rtol = 1.0e-5_wp
            ftol = 1.0e-5_wp
            tol_for_check = 1.0e-4_wp
        case(real64)
            atol = 1.0e-15_wp
            rtol = 1.0e-13_wp
            ftol = 1.0e-15_wp
            tol_for_check = 1.0e-7_wp
        case(real128)
            atol = 1.0e-25_wp
            rtol = 1.0e-23_wp
            ftol = 1.0e-25_wp
            tol_for_check = 1.0e-16_wp
        case default
            error stop 'unknown real kind'
        end select

        write(output_unit,fmt) &
            repeat('-',20),repeat('-',3),repeat('-',4),repeat('-',16),&
            repeat('-',25),repeat('-',16),repeat('-',5),repeat('-',5),repeat('-',8)

        write(output_unit,fmt) &
            'method','function','n','error','x','f','evals','iflag','runs/sec'

        write(output_unit,fmt) &
            repeat('-',20),repeat('-',3),repeat('-',4),repeat('-',16),&
            repeat('-',25),repeat('-',16),repeat('-',5),repeat('-',5),repeat('-',8)

        do imeth = 1, number_of_methods

            call problems(ax=ax, bx=bx, xroot=root)

            call cpu_time(tstart)
            do irepeat = 1, n_repeat
                ifunc = 0 ! reset func evals counter
                call root_scalar(methods(imeth),func,ax,bx,xzero,fzero,iflag,bisect_on_failure=.true., &
                                 atol = atol, rtol = rtol, ftol = ftol, maxiter = maxiter)
            end do
            call cpu_time(tfinish)

            error = xzero-root
            write(line, dfmt) trim(methods(imeth)),nprob,n,error,xzero,fzero,ifunc,iflag,n_repeat / (tfinish-tstart)

            root_found = abs(fzero) <= tol_for_check !.and. abs(error) <= tol_for_check

            if (root_found) then
                write(output_unit, '(A)') trim(line)
            else
                write(output_unit, '(A)') colorize(trim(line), color_fg='red') ! failed case
            end if

            ! save results for this case:
            if (iflag==0 .and. root_found) then
                fevals(imeth) = ifunc
                cases_per_sec(imeth) = n_repeat / (tfinish-tstart)
            else
                fevals(imeth) = huge(1)
                cases_per_sec(imeth) = huge(1.0_wp)
            end if

        end do

        ! update plot:
        call stats_plot%add_plot([(real(i,wp), i=1,number_of_methods)],real(fevals,wp),&
                                    label='',linestyle='-',linewidth=2,istat=istat)

        ! get the best and the failures for the output files:

        ! best in terms of lowest function evals:
        best_feval = minval(fevals)
        if (best_feval==huge(1)) best_feval = -1 ! if none of them converged
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

    !     ! best in terms of fastest cpu time:
    !     best_cases_per_second = maxval(cases_per_sec)
    !     if (best_cases_per_second==huge(1.0_wp)) best_cases_per_second = -1 ! if none of them converged
    !     best = ''
    !     failures = ''
    !     do imeth = 1, number_of_methods
    !         if (cases_per_sec(imeth) == best_cases_per_second) then
    !             if (best=='') then
    !                 best = trim(methods(imeth))
    !             else
    !                 best = best//', '//trim(methods(imeth))
    !             end if
    !             number_of_wins(imeth) = number_of_wins(imeth) + 1
    !         elseif (cases_per_sec(imeth) == huge(1.0_wp)) then
    !             if (failures=='') then
    !                 failures = trim(methods(imeth))
    !             else
    !                 failures = failures//', '//trim(methods(imeth))
    !             end if
    !             number_of_failures(imeth) = number_of_failures(imeth) + 1
    !         end if
    !     end do

    ! output the report:
    if (best /= '')     write(iunit, '(I5,1X,I5,1X,I5,1X,A)')  nprob, n, best_feval, 'Best: '//trim(best)
    if (failures /= '') write(iunit_failed, '(I5,1X,I5,1X,A)') nprob, n, 'Failures: '//trim(failures)

    end subroutine test

    subroutine problems(x, ax, bx, fx, xroot, cases, num_of_problems, latex)

    !! returns all the information about the test case.

    implicit none

    real(wp),intent(in),optional  :: x       !! indep variable
    real(wp),intent(out),optional :: ax, bx  !! bounds
    real(wp),intent(out),optional :: fx      !! function value `f(x)`
    real(wp),intent(out),optional :: xroot   !! value of indep variable at the root
    integer,dimension(:),allocatable,intent(out),optional :: cases !! list of `n` cases to test
    integer,intent(out),optional :: num_of_problems  !! total number of problems
    character(len=:),allocatable,intent(out),optional :: latex !! LaTeX string of `f(x)`

    real(wp),parameter :: pi = acos(-1.0_wp)

    real(wp) :: a, b, root, f, dn, di, xi, t1, emx, ex
    integer :: i !! counter
    integer,dimension(:),allocatable :: ns

    dn = real(n, wp)
    ns = [1] ! default case
    if (present(latex)) latex = ''

    select case (nprob)
    case (1)
        a = pi/2.0_wp
        b = pi
        root = 1.8954942670339809E+00_wp
        if (present(x)) f = sin(x) - x/2.0_wp
        if (present(latex)) latex = '\sin x - x/2'
    case (2:11)

        select case (nprob)
        case(2)
            a = 1.0_wp + 1.0e-9_wp
            b = (2.0_wp)*(2.0_wp) - 1.0e-9_wp
            root = 3.0229153472730570E+00_wp
        case (3)
            a = (2.0_wp)*(2.0_wp) + 1.0e-9_wp
            b = (3.0_wp)*(3.0_wp) - 1.0e-9_wp
            root = 6.6837535608080781E+00_wp
        case (4)
            a = (3.0_wp)*(3.0_wp) + 1.0e-9_wp
            b = (4.0_wp)*(4.0_wp) - 1.0e-9_wp
            root = 1.1238701655002212E+01_wp
        case (5)
            a = (4.0_wp)*(4.0_wp) + 1.0e-9_wp
            b = (5.0_wp)*(5.0_wp) - 1.0e-9_wp
            root = 1.9676000080623409E+01_wp
        case (6)
            a = (5.0_wp)*(5.0_wp) + 1.0e-9_wp
            b = (6.0_wp)*(6.0_wp) - 1.0e-9_wp
            root = 2.9828227326504754E+01_wp
        case (7)
            a = (6.0_wp)*(6.0_wp) + 1.0e-9_wp
            b = (7.0_wp)*(7.0_wp) - 1.0e-9_wp
            root = 4.1906116195289413E+01_wp
        case (8)
            a = (7.0_wp)*(7.0_wp) + 1.0e-9_wp
            b = (8.0_wp)*(8.0_wp) - 1.0e-9_wp
            root = 5.5953595800143094E+01_wp
        case (9)
            a = (8.0_wp)*(8.0_wp) + 1.0e-9_wp
            b = (9.0_wp)*(9.0_wp) - 1.0e-9_wp
            root = 7.1985665586587795E+01_wp
        case (10)
            a = (9.0_wp)*(9.0_wp) + 1.0e-9_wp
            b = (10.0_wp)*(10.0_wp) - 1.0e-9_wp
            root = 9.0008868539166666E+01_wp
        case (11)
            a = (10.0_wp)*(10.0_wp) + 1.0e-9_wp
            b = (11.0_wp)*(11.0_wp) - 1.0e-9_wp
            root = 1.1002653274833019E+02_wp
        end select
        if (present(x)) then
            f = 0.0_wp
            do i = 1, 20
                di = real(i,wp)
                f = f + ((2.0_wp*di - 5.0_wp)**2)/(x - di*di)**3
            end do
            f = -2.0_wp*f
        end if
        if (present(latex)) latex = '-2 \sum_{i=1}^{20} ((2i - 5)^2)/(x - i^2)^3'
    case (12)
        a = -9.0_wp
        b = 31.0_wp
        root = 0.0_wp
        if (present(x)) f = -40.0_wp*x*exp(-1.0_wp*x)
        if (present(latex)) latex = '-40 x \mathrm{e}^{-x}'
    case (13)
        a = -9.0_wp
        b = 31.0_wp
        root = 0.0_wp
        if (present(x)) f = -100.0_wp*x*exp(-2.0_wp*x)
        if (present(latex)) latex = '-100 x \mathrm{e}^{-2 x}'
    case (14)
        a = -9.0_wp
        b = 31.0_wp
        root = 0.0_wp
        if (present(x)) f = -200.0_wp*x*exp(-3.0_wp*x)
        if (present(latex)) latex = '-200 x \mathrm{e}^{-3 x}'
    case (15)
        a = 0.0_wp
        b = 5.0_wp
        ns = [4, 6, 8, 10, 12]
        if (present(xroot)) then
            select case (n)
            case(4);  root = 6.6874030497642202E-01_wp
            case(6);  root = 7.6472449133173001E-01_wp
            case(8);  root = 8.1776543395794250E-01_wp
            case(10); root = 8.5133992252078460E-01_wp
            case(12); root = 8.7448527222116784E-01_wp
            case default
                write(*,*) 'invalid n: ', n ; error stop
            end select
        end if
        if (present(x)) f = x**n - 0.2_wp
        if (present(latex)) latex = 'x^n - 0.2'
    case (16:17)
        select case (nprob)
        case(16)
            a = 0.0_wp
            b = 5.0_wp
            ns = [4, 6, 8, 10, 12]
        case (17)
            a = -0.95_wp
            b = 4.05_wp
            root = 1.0_wp
            ns = [8, 10, 12, 14]
        end select
        root = 1.0_wp
        if (present(x)) f = x**n - 1.0_wp
        if (present(latex)) latex = 'x^n - 1'
    case (18)  ! note: 18-24 are also in Dowell & Jarratt
        a = 0.0_wp
        b = 1.5_wp
        root = pi / 6.0_wp
        if (present(x)) f = sin(x) - 0.5_wp
        if (present(latex)) latex = '\sin x - 0.5'
    case (19)
        a = 0.0_wp
        b = 1.0_wp
        ns = [1,2,3,4,5,20,40,60,80,100]
        if (present(xroot)) then
            select case (n)
            case(1);    root = 4.2247770964123666E-01_wp
            case(2);    root = 3.0669941048320373E-01_wp
            case(3);    root = 2.2370545765466297E-01_wp
            case(4);    root = 1.7171914751950839E-01_wp
            case(5);    root = 1.3825715505682408E-01_wp
            case(20);   root = 3.4657359020853851E-02_wp
            case(40);   root = 1.7328679513998633E-02_wp
            case(60);   root = 1.1552453009332422E-02_wp
            case(80);   root = 8.6643397569993164E-03_wp
            case(100);  root = 6.9314718055994531E-03_wp
            case default; write(*,*) 'invalid n: ', n ; error stop
            end select
        end if
        if (present(x)) f = 2.0_wp*x*exp(-dn) - 2.0_wp*exp(-dn*x) + 1.0_wp
        if (present(latex)) latex = '2 x \mathrm{e}^{-n} - 2 \mathrm{e}^{-n x} + 1'
    case (20)
        a = 0.0_wp
        b = 1.0_wp
        ns = [5,10,20]
        if (present(xroot)) then
            select case (n)
            case(5);    root = 3.8402551840621900E-02_wp
            case(10);   root = 9.9000099980004999E-03_wp
            case(20);   root = 2.4937500390620117E-03_wp
            case default; write(*,*) 'invalid n: ', n ; error stop
            end select
        end if
        if (present(x)) f = (1.0_wp + (1.0_wp - dn)**2)*x - (1.0_wp - dn*x)**2
        if (present(latex)) latex = '(1 + (1 - n)^2)x - (1 - n x)^2'
    case (21)
        a = 0.0_wp
        b = 1.0_wp
        ns = [2,5,10,15,20]
        if (present(xroot)) then
            select case (n)
            case(2);    root = 0.5_wp
            case(5);    root = 3.4595481584824202E-01_wp
            case(10);   root = 2.4512233375330724E-01_wp
            case(15);   root = 1.9554762353656561E-01_wp
            case(20);   root = 1.6492095727644095E-01_wp
            case default; write(*,*) 'invalid n: ', n ; error stop
            end select
        end if
        if (present(x)) f = x**2 - (1.0_wp - x)**n
        if (present(latex)) latex = 'x^2 - (1-x)^n'
    case (22)
        a = 0.0_wp
        b = 1.0_wp
        ns = [1,2,4,5,8,15,20]
        if (present(xroot)) then
            select case (n)
            case(1);    root = 2.7550804099948439E-01_wp
            case(2);    root = 1.3775402049974219E-01_wp
            case(4);    root = 1.0305283778156444E-02_wp
            case(5);    root = 3.6171081789040635E-03_wp
            case(8);    root = 4.1087291849639540E-04_wp
            case(15);   root = 2.5989575892907627E-05_wp
            case(20);   root = 7.6685951221853367E-06_wp
            case default; write(*,*) 'invalid n: ', n ; error stop
            end select
        end if
        if (present(x)) f = (1.0_wp + (1.0_wp - dn)**4)*x - (1.0_wp - dn*x)**4
        if (present(latex)) latex = '(1 + (1 - n)^4) x - (1 - n x)^4'
    case (23)
        a = 0.0_wp
        b = 1.0_wp
        ns = [1,5,10,15,20]
        if (present(xroot)) then
            select case (n)
            case(1);    root = 4.0105813754154704E-01_wp
            case(5);    root = 5.1615351875793357E-01_wp
            case(10);   root = 5.3952222690841584E-01_wp
            case(15);   root = 5.4818229434065527E-01_wp
            case(20);   root = 5.5270466667848779E-01_wp
            case default; write(*,*) 'invalid n: ', n ; error stop
            end select
        end if
        if (present(x)) f = (x - 1.0_wp)*exp(-dn*x) + x**n
        if (present(latex)) latex = '(x - 1) \mathrm{e}^{-n x} + x^n'
    case (24)
        a = 0.01_wp
        b = 1.0_wp
        ns = [2,5,15,20]
        if (present(xroot)) then
            select case (n)
            case(2);    root = 0.5_wp
            case(5);    root = 0.2_wp
            case(15);   root = 1.0_wp / 15.0_wp !6.6666666666667e-02_wp
            case(20);   root = 5.0e-02_wp
            case default; write(*,*) 'invalid n: ', n ; error stop
            end select
        end if
        if (present(x)) f = (dn*x - 1.0_wp)/((dn - 1.0_wp)*x)
        if (present(latex)) latex = '\frac{n x - 1}{(n - 1)x}'
    case (25)
        a = 1.0_wp
        b = 100.0_wp
        ns = [2,3,4,5,6,7,9,11,13,15,17,19,21,23,25,27,29,31,33]
        root = real(n,wp)
        if (present(x)) f = x**(1.0_wp/dn) - dn**(1.0_wp/dn)
        if (present(latex)) latex = 'x^{1/n} - n^{1/n}'
    case (26)
        a = -1.0_wp
        b = 4.0_wp
        root = 0.0_wp    ! this is just almost impossible to get
        if (present(x)) then
            if (x == 0.0_wp) then
                f = 0.0_wp
            else
                f = x/exp(1.0_wp/(x*x))
            end if
        end if
        if (present(latex)) latex = '\left\{ \begin{array}{cl} 0 & \mathrm{if~}x=0 \\'//&
                                    ' x \mathrm{e}^{-1/x^2} & \mathrm{otherwise} \end{array} \right.'
    case (27)
        a = -10000.0_wp
        b = pi/2.0_wp
        root = 6.2380651896161232E-01_wp
        ns = [(i, i = 1, 40)]
        if (present(x)) then
            if (x >= 0.0_wp) then
                f = (x/1.5_wp + sin(x) - 1.0_wp)*dn/20.0_wp
            else
                f = (-1.0_wp*dn)/20.0_wp
            end if
        end if
        if (present(latex)) latex = '\left\{ \begin{array}{cl} (x/1.5 + \sin(x) - 1) n/20 & \mathrm{if~}x \ge 0 \\'//&
                                    ' -n/20 & \mathrm{otherwise} \end{array} \right.'
    case (28)
        a = -10000.0_wp
        b = 1.0e-4_wp
        ns = [ [(i, i = 20,40)], [(i*100, i = 1, 10)] ]
        if (present(xroot)) then
            select case (n)
            case(20);   root = 5.9051305594219711E-05_wp
            case(21);   root = 5.6367155339936997E-05_wp
            case(22);   root = 5.3916409455591910E-05_wp
            case(23);   root = 5.1669892394942247E-05_wp
            case(24);   root = 4.9603096699144557E-05_wp
            case(25);   root = 4.7695285287638997E-05_wp
            case(26);   root = 4.5928793239948664E-05_wp
            case(27);   root = 4.4288479195664783E-05_wp
            case(28);   root = 4.2761290257883239E-05_wp
            case(29);   root = 4.1335913915953798E-05_wp
            case(30);   root = 4.0002497338019804E-05_wp
            case(31);   root = 3.8752419296206685E-05_wp
            case(32);   root = 3.7578103559957998E-05_wp
            case(33);   root = 3.6472865219959233E-05_wp
            case(34);   root = 3.5430783356531827E-05_wp
            case(35);   root = 3.4446594929961498E-05_wp
            case(36);   root = 3.3515605877800376E-05_wp
            case(37);   root = 3.2633616249437209E-05_wp
            case(38);   root = 3.1796856858425998E-05_wp
            case(39);   root = 3.1001935436965348E-05_wp
            case(40);   root = 3.0245790670210096E-05_wp
            case(100);  root = 1.2277994232461524E-05_wp
            case(200);  root = 6.1695393904408653E-06_wp
            case(300);  root = 4.1198585298292822E-06_wp
            case(400);  root = 3.0924623877272168E-06_wp
            case(500);  root = 2.4752044261050178E-06_wp
            case(600);  root = 2.0633567678512711E-06_wp
            case(700);  root = 1.7690120078154264E-06_wp
            case(800);  root = 1.5481615698859100E-06_wp
            case(900);  root = 1.3763345366022352E-06_wp
            case(1000); root = 1.2388385788997142E-06_wp
            case default; write(*,*) 'invalid n: ', n ; error stop
            end select
        end if
        if (present(x)) then
            if (x >= (1.0e-3_wp)*2.0_wp/(dn + 1.0_wp)) then
                f = exp(1.0_wp) - 1.859_wp
            elseif (x >= 0.0_wp) then
                f = exp((dn + 1.0_wp)*0.5_wp*x*1000.0_wp) - 1.859_wp
            else
                f = -0.859_wp
            end if
        end if
        if (present(latex)) latex = '\left\{ \begin{array}{cl} '//&
                                    ' \mathrm{e} - 1.859 & \mathrm{if~} x \ge 0.002/(n + 1) \\'//&
                                    ' \mathrm{e}^{500 (n + 1) x} - 1.859 & x \ge 0\\'//&
                                    ' -0.859 & \mathrm{otherwise}'//&
                                    ' \end{array} \right.'
    case (29)
        ! Zhang test case
        a = 0.0_wp
        b = 4.0_wp
        root = 8.6547403310161445E-01_wp
        if (present(x)) f = cos(x) - x**3
        if (present(latex)) latex = '\cos x - x^3'

    ! 30-36 : Gottlieb's paper [table 1 & 2]
    case (30)
        a = 0.0_wp
        b = 1.0_wp
        root = asin(2.0_wp/3.0_wp)
        if (present(x)) f = 3.0_wp * sin(x) - 2.0_wp
        if (present(latex)) latex = '3 \sin x - 2'
    case (31)
        a = -1.0_wp
        b = 1.0_wp
        root = 5.6714329040978387E-01_wp
        if (present(x)) f = x*exp(x) - 1.0_wp
        if (present(latex)) latex = 'x \mathrm{e}^x - 1'
    case (32)
        a = 0.1_wp
        b = 0.9_wp
        root = 8.0413309750366432E-01_wp
        if (present(x)) f = 11.0_wp * x**11 - 1.0_wp
        if (present(latex)) latex = '11 x^{11} - 1'
    case (33)
        a = 2.8_wp
        b = 3.1_wp
        root = 3.0_wp
        if (present(x)) f = exp(x**2 + 7.0_wp*x - 30.0_wp) - 1.0_wp
        if (present(latex)) latex = '\mathrm{e}^{x^2 + 7x - 30} - 1'
    case (34)
        a = -1.3_wp
        b = -0.5_wp
        root = -6.2944648407333333E-01_wp
        if (present(x)) f = 1.0_wp/x - sin(x) + 1.0_wp
        if (present(latex)) latex = '1/x - \sin x + 1'
    case (35)
        a = 2.0_wp
        b = 3.0_wp
        root = 2.0945514815423266E+00_wp
        if (present(x)) f = x**3 - 2.0_wp*x - 5.0_wp
        if (present(latex)) latex = 'x^3 - 2 x - 5'
    case (36)
        a = 0.5_wp
        b = 2.0_wp
        root = 1.0_wp
        if (present(x)) f = 1.0_wp/x - 1.0_wp
        if (present(latex)) latex = '1/x - 1'
    case (37)
        a = -10.0_wp
        b = 10.0_wp
        root = 1.0_wp
        if (present(x)) f = x**3 - 5.0_wp*x**2 + 12.0_wp*x - 8.0_wp
        if (present(latex)) latex = 'x^3 - 5 x^2 + 12 x - 8'
    case (38)
        a = 0.0_wp
        b = 10.0_wp
        root = 1.0_wp
        if (present(x)) f = x**5 - 61.0_wp*x**4 + 1368.0_wp*x**3 - &
                            13548.0_wp*x**2 + 54256.0_wp*x - 42016.0_wp ! x in [0, 10]
        if (present(latex)) latex = 'x^5 - 61 x^4 + 1368 x^3 - 13548 x^2 + 54256 x - 42016'
    case (39)
        a = 0.0_wp
        b = 10.0_wp
        root = 1.0_wp
        if (present(x)) f = x**5-61.0_wp*x**4+6801.0_wp/5.0_wp*x**3-66531.0_wp/&
                            5.0_wp*x**2+5205601.0_wp/100.0_wp*x-4005001.0_wp/100.0_wp
        if (present(latex)) latex = 'x^5 - 61 x^4 + 6801 x^3 / 5 - 66531 x^2 / 5 + 5205601 x / 100 - 4005001/100'
    case (40)
        a = 0.0_wp
        b = 10.0_wp
        root = 1.0_wp
        if (present(x)) f = (x-1.0_wp)**5
        if (present(latex)) latex = '(x-1)^5'
    case (41)
        !From: http://www.mth.pdx.edu/~daescu/mth451_551/Muller_bisection.pdf
        a = -2.0_wp
        b = 3.0_wp
        root = 0.0_wp
        if (present(x)) f = atan(x)
        if (present(latex)) latex = '\tan^{-1} x'
    case (42)
        !From: http://www.mth.pdx.edu/~daescu/mth451_551/Muller_bisection.pdf
        a = 0.0_wp     ! x1 is the location of the root for this case
        b = 3.0_wp
        root = 0.0_wp
        if (present(x)) f = exp(x) - 2.0_wp*x - 1.0_wp
        if (present(latex)) latex = '\mathrm{e}^{x} - 2 x - 1'
    case (43)
        a = 0.0_wp
        b = 0.7_wp
        root = 6.7771292632868404E-01_wp
        if (present(x)) f = 3.0_wp * x * sin(x*20.0_wp)*cos(x) + x - 2.0_wp
        if (present(latex)) latex = '3 x \sin(20 x) \cos x + x - 2'
    case (44)
        a = -0.7_wp
        b = 0.4_wp
        root = -6.4008369608468403E-01_wp
        if (present(x)) f = 3.0_wp * x**2 * sin(x*20.0_wp) + cos(x*2.0_wp)
        if (present(latex)) latex = '3 x^2 \sin (20x) + \cos(2x)'
    case (45)
        !case 11 from "Algorithm 748: Enclosing Zeros of Continuous Functions" (n=20)
        a = 0.01_wp
        b = 1.0_wp
        root = 0.05_wp
        ns = [20]
        if (present(x)) f = (real(n,wp)*x-1.0_wp)/((real(n,wp)-1.0_wp)*x)
        if (present(latex)) latex = '\frac{n x - 1}{(n - 1) x}'
    case (46)
        a = -10.0_wp
        b = 5.0_wp
        root = 0.0_wp
        if (present(x)) f = sin(x*10.0_wp) + sin(x*2.0_wp) + &
                            atan(x/3.0_wp)*x**2 + exp(x) - 1.0_wp
        if (present(latex)) latex = '\sin(10 x) + \sin(2 x) + \atan(x/3) x^2 + \mathrm{e}^x - 1'

    ! functions are from:
    !   T. R. Chandrupatla, "A new hybrid quadratic/bisection algorithm for finding the zero of
    !   a nonlinear function without using derivatives", Advances in Engineering Software
    !   Volume 28, Issue 3, April 1997, Pages 145-149

    case (47:48)
        select case (nprob)
        case (47)
            a = 0.5_wp
            b = 1.51_wp
        case (48)
            a =  1.0e-12_wp
            b =  1.0e12_wp
        end select
        root = 1.0_wp
        if (present(x)) f = 1.0_wp - 1.0_wp/(x**2)
        if (present(latex)) latex = '1 - 1/x^2'
    case (49:50)
        select case (nprob)
        case (49)
            a = 0.0_wp
            b = 5.0_wp
        case (50)
            a = -1.0e10_wp
            b =  1.0e10_wp
        end select
        root = 3.0_wp
        if (present(x)) f = (x-3.0_wp)**3
        if (present(latex)) latex = '(x - 3)^3'
    case (51:52)
        select case (nprob)
        case (51)
            a = 0.0_wp
            b = 5.0_wp
        case (52)
            a = -1.0e10_wp
            b =  1.0e10_wp
        end select
        root = 2.0_wp
        if (present(x)) f = 6.0_wp*(x-2.0_wp)**5
        if (present(latex)) latex = '6 (x-2)^5'
    case (53:54)
        select case (nprob)
        case (53)
            a = -1.0_wp
            b =  4.0_wp
        case (54)
            a = -10.0_wp
            b =  100.0_wp
        end select
        root = 0.0_wp
        if (present(x)) f = x**9
        if (present(latex)) latex = 'x^9'
    case (55:56)
        select case (nprob)
        case (55)
            a = -1.0_wp
            b =  4.0_wp
        case (56)
            a = -10.0_wp
            b =  100.0_wp
        end select
        root = 0.0_wp
        if (present(x)) f = x**19
        if (present(latex)) latex = 'x^{19}'
    case (57:58)
        select case (nprob)
        case (57)
            a = -1.0_wp
            b =  4.0_wp
        case (58)
            a = -10.0_wp
            b =  100.0_wp
        end select
        root = 0.0_wp
        if (present(x)) then
            if (abs(x) < 3.8e-4_wp) then  ! same as one above except for this
                f = 0.0_wp
            else
                f = x*exp(-x**(-2))
            endif
        end if
        if (present(latex)) latex = '\left\{ \begin{array}{cl} 0 & \mathrm{if~} | x | < 0.00038 \\'//&
                                    ' x \mathrm{e}^{-1/x^2} & \mathrm{otherwise} \end{array} \right.'
    case (59:60)
        select case (nprob)
        case (59)
            a = 2.0e-4_wp
            b = 2.0_wp
        case (60)
            a = 2.0e-4_wp
            b = 81.0_wp
        end select
        root = 1.0375360332870403E+00_wp
        if (present(x)) then
            xi = 0.61489_wp
            t1 = 1.0_wp-xi
            emx = exp(-x)
            f = -(3062.0_wp*t1*emx)/(xi + t1*emx) - 1013.0_wp + 1628.0_wp/x
        end if
        if (present(latex)) latex = '-\frac{3062 * 0.38511 * \mathrm{e}^{-x}}{0.61489 + 0.38511 * \mathrm{e}^{-x}} - 1013 + 1628/x'
    case (61:62)
        select case (nprob)
        case (61)
            a = 2.0e-4_wp
            b = 1.0_wp
        case (62)
            a = 2.0e-4_wp
            b = 81.0_wp
        end select
        root = 7.0320484036313581E-01_wp
        if (present(x)) then
            ex = exp(x)
            f = ex - 2.0_wp - 0.01_wp/(x*x) + 2.0e-6_wp/(x*x*x)
        end if
        if (present(latex)) latex = '\mathrm{e}^x - 2 - 0.01/x^2 + 0.000002/x^3'

    ! from A modified three-point Secant method with improved
    !      rate and characteristics of convergence
    !      July 2019

    case(63)
        a = 0.0_wp
        b = 1.5_wp
        root = 1.3652300134140968E+00_wp
        if (present(x)) f = x**3 + 4.0_wp*x**2 - 10.0_wp
        if (present(latex)) latex = 'x^3 + 4 * x^2 - 10'
    case(64)
        a = -2.0_wp
        b = -1.0_wp
        root = -1.4044916482153412E+00_wp
        if (present(x)) f = sin(x)**2 - x**2 + 1.0_wp
        if (present(latex)) latex = '\sin^2 x - x^2 + 1'
    case(65)
        a = 1.0_wp
        b = 2.1_wp
        root = 2.0_wp
        if (present(x)) f = (x - 1.0_wp)**6 - 1.0_wp
        if (present(latex)) latex = '(x - 1)^6 - 1'
    case(66)
        a = -0.9_wp
        b = -0.1_wp
        root = -6.0323197155721517E-01_wp
        if (present(x)) f = sin(x)*exp(x) + log(x**2+1.0_wp)
        if (present(latex)) latex = '\sin(x) \mathrm{e}^x + \log(x^2 + 1)'
    case(67)
        a = 1.0_wp
        b = 3.1_wp
        root = 3.0_wp
        if (present(x)) f = exp(x**2+7.0_wp*x-30.0_wp) - 1.0_wp
        if (present(latex)) latex = '\mathrm{e}^{x^2 + 7 x - 30} - 1'
    case(68)
        a = 1.0_wp
        b = 2.1_wp
        root = 1.8571838602078353E+00_wp
        if (present(x)) f = x - 3.0_wp * log(x)
        if (present(latex)) latex = 'x - 3 \log x'

    ! some simple functions
    case(69)
        a = 0.1001_wp
        b = 10.0001_wp
        root = 1.0_wp
        if (present(x)) f = x**2 - 1.0_wp
        if (present(latex)) latex = 'x^2 - 1'
    case(70)
        a = 0.0001_wp
        b = 10.0001_wp
        root = 1.0018618683298285E+00_wp
        if (present(x)) f = x**2 - 1.0038273894_wp + sin(x/10000.0_wp)
        if (present(latex)) latex = 'x^2 - 1.0038273894 + \sin (x/10000)'
    case(71:72)
        select case (nprob)
        case (71)
            a = -2.0001_wp
            b = 1.0001_wp
        case(72)
            a = -0.0001_wp  ! root very close to a
            b = 1.0001_wp
        end select
        root = 0.0_wp
        if (present(x)) f = x
        if (present(latex)) latex = 'x'

    ! D. Popovski, "A note on King's method F for finding a bracketed root", 1982.
    ! https://www.semanticscholar.org/paper/A-note-on-King's-method-F-for-finding-a-bracketed-Popovski/513ffdb3befbf48c38bb70dbec60934472cc9903
    case(73)
        a = 0.0_wp
        b = 5.0_wp
        root =  1.4655712318767680E+00_wp
        if (present(x)) f = (x - 1.0_wp)*x*x - 1.0_wp
        if (present(latex)) latex = '(x - 1) x^2 - 1'

    case(74)
        a = 0.0_wp
        b = 8.0_wp
        root =  5.0275246628429326E+00_wp
        if (present(x)) f = ((x - 3.0_wp)*x - 9.0_wp)*x - 6.0_wp
        if (present(latex)) latex = '((x - 3) x - 9) x - 6'
    case(75:78)
        select case (nprob)
        case (75)
            a = -5.0_wp
            b = -3.0_wp
            root = -3.5437609487522587E+00_wp
        case (76)
            a = -3.0_wp
            b = -2.0_wp
            root = -2.2701871608581255E+00_wp
        case (77)
            a = -2.0_wp
            b = 0.0_wp
            root = -2.0646554491023637E-01_wp
        case (78)
            a = 0.0_wp
            b = 9.0_wp
            root = 6.0204136545206206E+00_wp
        end select
        if (present(x)) f = ((x*x - 27.0_wp)*x - 54.0_wp)*x - 10.0_wp
        if (present(latex)) latex = '((x^2 - 27) x - 54 ) x - 10'
    case(79:82)
        select case (nprob)
        case (79)
            a = 0.0_wp
            b = 2.0_wp
            root =  1.1893836076828531E+00_wp
        case (80)
            a = 2.0_wp
            b = 4.0_wp
            root = 3.4254515851826428E+00_wp
        case (81)
            a = 4.0_wp
            b = 8.0_wp
            root = 6.5745484148173572E+00_wp
        case (82)
            a = 8.0_wp
            b = 9.0_wp
            root = 8.8106163923171469E+00_wp
        end select
        if (present(x)) f = (((x-20.0_wp)*x + 133.0_wp)*x - 330.0_wp)*x + 236.0_wp
        if (present(latex)) latex = '(((x-20) x + 133) x - 330) x + 236'
    case(83)
        a = 0.0_wp
        b = 3.0_wp
        root =  1.2554228710768465E+00_wp
        if (present(x)) f = (x - 1.0_wp)*x**6 - 1.0_wp
        if (present(latex)) latex = '(x - 1) x^6 - 1'
    case(84:85)
        select case (nprob)
        case(84)
            a = -5.0_wp
            b = 0.0_wp
            root = -1.3926261305968417E+00_wp
        case(85)
            a = 0.0_wp
            b = 5.0_wp
            root = 3.0001928237608510E+00_wp
        end select
        if (present(x)) f = (((x*x - 6.0_wp)*x - 8.0_wp)*x - 3.0_wp)*x**4 - 1.0_wp
        if (present(latex)) latex = '(((x^2 - 6) x - 8) x - 3) x^4 - 1'

    case(86:88)
        select case (nprob)
        case(86)
            a = 0.0_wp
            b = 2.0_wp
            root =  1.0883595391412838E+00_wp
        case(87)
            a = 2.0_wp
            b = 6.0_wp
            root = 4.6333772273659109E+00_wp
        case(88)
            a = 6.0_wp
            b = 8.0_wp
            root = 7.2609923862369942E+00_wp
        end select
        if (present(x)) f = ((x - 13.0_wp)*x + 47.0_wp)*x - 36.0_wp - sqrt(x)
        if (present(latex)) latex = '((x - 13) x + 47 ) x - 36 - \sqrt{x}'
    case(89:90)
        select case (nprob)
        case (89)
            a = 0.0_wp
            b = 2.0_wp
            root =  1.1118325591589630E+00_wp
        case (90)
            a = 2.0_wp
            b = 5.0_wp
            root = 4.5771520639572972E+00_wp
        end select
        if (present(x)) f = exp(1.0_wp-x)*(x-1.0_wp)*10.0_wp - 1.0_wp
        if (present(latex)) latex = '10 \mathrm{e}^{1-x} (x-1) - 1'

    case(91)
        a = 2.0_wp
        b = 3.0_wp
        root =  2.8424389537844471E+00_wp
        if (present(x)) f = exp(x) + x - 20.0_wp
        if (present(latex)) latex = '\mathrm{e}^x + x - 20'

    case(92)
        a = 0.0_wp
        b = 2.0_wp
        root =  1.3065586410393502E+00_wp   ! the root in the paper is wrong for this one !
        if (present(x)) f = exp(x) + x - 5.0_wp
        if (present(latex)) latex = '\mathrm{e}^x + x - 5'

    case(93:94)
        select case (nprob)
        case (93)
            a = 0.0_wp
            b = 6.0_wp
            root =  3.1094324994726575E+00_wp
        case (94)
            a = 6.0_wp
            b = 8.0_wp
            root = 7.0453708632631685E+00_wp
        end select
        if (present(x)) f = (x - 10.0_wp)*x + 23.0_wp - x**0.4_wp
        if (present(latex)) latex = '(x - 10) x + 23 - x^{0.4}'

    case(95:98)
        select case (nprob)
        case (95)
            a = 0.0_wp
            b = 1.0_wp
            root =  3.6682419352992085E-01_wp
        case (96)
            a = 1.0_wp
            b = 5.0_wp
            root = 3.5735001943758262E+00_wp
        case (97)
            a = 5.0_wp
            b = 7.0_wp
            root = 5.7880397474352980E+00_wp
        case (98)
            a = 7.0_wp
            b = 10.0_wp
            root = 9.2131701574933741E+00_wp
        end select
        if (present(x)) f = (0.04_wp*x - 0.4_wp)*x + 0.5_wp - sin(x)
        if (present(latex)) latex = '(0.04 x - 0.4) x + 0.5 - \sin x'

    ! cases from Stage paper
    ! see also #26
    case(99)
        a = -10.0_wp
        b = 10.0_wp
        root = -7.3908513321516067E-01_wp
        if (present(x)) f = cos(x) + x
        if (present(latex)) latex = '\cos x + x'
    case(100)
        a = -10.0_wp
        b = 10.0_wp
        root = 4.0_wp / 3.0_wp
        if (present(x)) f = 1.0_wp - 0.75_wp * x
        if (present(latex)) latex = '1 - 0.75 x'

    case(101)
        a = -10.0_wp
        b = 10.0_wp
        root = 2.0_wp / 3.0_wp - 0.01_wp
        if (present(x)) then
            if (x <= 2.0_wp / 3.0_wp) then
                f =  (abs(x - 2.0_wp / 3.0_wp))**0.5_wp - 0.1_wp
            else
                f = -(abs(x - 2.0_wp / 3.0_wp))**0.5_wp - 0.1_wp
            end if
        end if
        if (present(latex)) latex = '\left\{ \begin{array}{cl} |x - 2/3|^{0.5} - 0.1 & \mathrm{if~} x \le 2/3 \\'//&
                                    ' -|x - 2/3|^{0.5} - 0.1& \mathrm{otherwise} \end{array} \right.'
    ! case(102)
    !     a = -10.0_wp
    !     b = 10.0_wp
    !     root = 2.0_wp / 3.0_wp
    !     if (present(x)) then
    !         if (x <= 2.0_wp / 3.0_wp) then
    !             f = (abs(x - 2.0_wp / 3.0_wp))**0.5_wp
    !         else
    !             f = -(abs(x - 2.0_wp / 3.0_wp))**0.5_wp
    !         end if
    !     end if
    !     if (present(latex)) latex = 'TODO'
    ! case(103)
    !     a = -10.0_wp
    !     b = 10.0_wp
    !     root = 2.0_wp / 3.0_wp
    !     if (present(x)) then
    !         if (x <= 2.0_wp / 3.0_wp) then
    !             f = (abs(x - 2.0_wp / 3.0_wp))**0.2_wp
    !         else
    !             f = -(abs(x - 2.0_wp / 3.0_wp))**0.2_wp
    !         end if
    !     end if
    !     if (present(latex)) latex = 'TODO'
    case(102)
        a = -10.0_wp
        b = 10.0_wp
        root = 7.0_wp / 9.0_wp
        if (present(x)) f = (x - 7.0_wp/9.0_wp)**3 + 0.001_wp * (x - 7.0_wp/9.0_wp)
        if (present(latex)) latex = '(x - 7/9)^3 + 0.001 (x - 7/9)'

    ! P1 in BlendTF paper
    case(103)
        a = 1.0_wp
        b = 2.0_wp
        root = 1.7320508075688773E+00_wp
        if (present(x)) f = x**2 - 3.0_wp
        if (present(latex)) latex = 'x^2 - 3'
    case(104)
        a = 2.0_wp
        b = 7.0_wp
        root = 2.2360679774997897E+00_wp
        if (present(x)) f = x**2 - 5.0_wp
        if (present(latex)) latex = 'x^2 - 5'
    case(105)
        a = 3.0_wp
        b = 7.0_wp
        root = 3.1622776601683793E+00_wp
        if (present(x)) f = x**2 - 10.0_wp
        if (present(latex)) latex = 'x^2 - 10'
    case(106)
        a = 1.0_wp
        b = 4.0_wp
        root = 2.0_wp
        if (present(x)) f = x**2 - x - 2.0_wp
        if (present(latex)) latex = 'x^2 - x - 2'
    case(107)
        a = 1.0_wp
        b = 3.0_wp
        root = 1.8284271247461901E+00_wp
        if (present(x)) f = x**2 + 2.0_wp * x - 7.0_wp
        if (present(latex)) latex = 'x^2 + 2 x - 7'
    case(108)
        a = 0.0_wp
        b = 2.0_wp
        root = 1.2599210498948732E+00_wp
        if (present(x)) f = x**3 - 2.0_wp
        if (present(latex)) latex = 'x^3 - 2'
    case(109)
        a = 0.0_wp
        b = 2.0_wp
        root = 1.5243452049841444E+00_wp
        if (present(x)) f = x * exp(x) - 7.0_wp
        if (present(latex)) latex = 'x \mathrm{e}^x - 7'
    case(110)
        a = 0.0_wp
        b = 1.0_wp
        root = 7.3908513321516064E-01_wp
        if (present(x)) f = x  - cos(x)
        if (present(latex)) latex = 'x  - \cos(x)'
    case(111)
        a = 0.0_wp
        b = 2.0_wp
        root = 1.1141571408719301E+00_wp
        if (present(x)) f = x*sin(x) - 1.0_wp
        if (present(latex)) latex = 'x \sin(x) - 1'
    case(112)
        a = -2.0_wp
        b = 4.0_wp
        root = 2.0739328090912149E+00_wp
        if (present(x)) f = x*cos(x) + 1.0_wp
        if (present(latex)) latex = 'x \cos(x) + 1'
    case(113)
        a = 0.0_wp
        b = 1.3_wp
        root = 1.0_wp
        if (present(x)) f = x**10 - 1.0_wp
        if (present(latex)) latex = 'x^{10} - 1'
    case(114)
        a = 1.0_wp
        b = 2.0_wp
        root = 1.6490132683031901E+00_wp
        if (present(x)) f = x**2 + exp(x/2.0_wp) - 5.0_wp
        if (present(latex)) latex = 'x^2 + \mathrm{e}^{x/2} - 5'
    case(115)
        a = 3.0_wp
        b = 4.0_wp
        root = 3.2215883990939420E+00_wp
        if (present(x)) f = sin(x) * sinh(x) + 1.0_wp
        if (present(latex)) latex = '\sin(x) \sinh(x) + 1'
    case(116)
        a = 2.0_wp
        b = 3.0_wp
        root = 2.1253911988111299E+00_wp
        if (present(x)) f = exp(x) - 3.0_wp * x - 2.0_wp
        if (present(latex)) latex = '\mathrm{e}^x - 3 x - 2'
    case(117)
        a = 0.5_wp
        b = 1.0_wp
        root = 8.7672621539506245E-01_wp
        if (present(x)) f = sin(x) - x**2
        if (present(latex)) latex = '\sin(x) - x^2'

    case default
        write(*,*) 'invalid case: ', nprob
        error stop 'invalid case'
    end select

    if (present(num_of_problems)) num_of_problems = 117

    ! outputs:
    if (present(ax))    ax = a
    if (present(bx))    bx = b
    if (present(xroot)) xroot = root
    if (present(fx))    fx = f
    if (present(cases)) cases = ns

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

    subroutine generate_plots()

    !! generate a plot of each function.

    implicit none

    integer,parameter :: n_points = 1000
    character(len=*),parameter :: ext = '.pdf' !! plot file type

    type(pyplot) :: plt   !! pytplot handler
    integer :: num_of_problems
    integer :: istat,i,ic
    character(len=10) :: s1,s2
    real(wp) :: x,ax,bx,xroot,froot,delta
    integer,dimension(:),allocatable :: cases_to_run
    real(wp),dimension(n_points) :: xvec
    real(wp),dimension(n_points) :: yvec
    character(len=:),allocatable :: filename, title, latex

    write(*,*) 'generate_plots...'

    n = 1  ! initialize
    nprob = 1
    call problems(num_of_problems=num_of_problems)

    do nprob = 1, num_of_problems

        write(*,*) 'case', nprob, '/', num_of_problems

        call problems(cases=cases_to_run,latex=latex)
        write(s1,'(I10)') nprob; s1 = trim(adjustl(s1)) ! int to str

        do ic = 1, size(cases_to_run)

            n = cases_to_run(ic)

            write(s2,'(I10)') n; s2 = trim(adjustl(s2)) ! int to str

            call problems(ax=ax, bx=bx, xroot=xroot)
            call problems(x=xroot, fx=froot)
            delta = (bx-ax)/(n_points-1)

            do i = 1, n_points
                xvec(i) = ax + delta * (i-1)
                call problems(x=xvec(i), fx=yvec(i))
            end do

            if (size(cases_to_run)==1) then
                title = 'Case '//trim(s1)//' : $f(x) = '//latex//'$'
                filename = 'plottest_'//'CASE='//trim(s1)//ext
            else
                title = 'Case '//trim(s1)//' (n='//trim(s2)//')'//' : $f(x) = '//latex//'$'
                filename = 'plottest_'//'CASE='//trim(s1)//'_N='//trim(s2)//ext
            end if

            call plt%initialize(grid=.true.,xlabel='x',ylabel='f(x)',figsize=[10,5],&
                                title=title,legend=.true.,&
                                tight_layout=.true.,&
                                real_fmt = '(E30.16E3)',&
                                usetex=.true.,raw_strings=.true.)
            call plt%add_plot(xvec,yvec,label='f(x)',linestyle='-',markersize=5,linewidth=2,istat=istat)
            call plt%add_plot([xroot],[froot],label='root',linestyle='.',markersize=10,linewidth=2,istat=istat)
            call plt%savefig(filename,istat=istat)

        end do

    end do

    end subroutine generate_plots

!*****************************************************************************************
    end program root_tests
!*****************************************************************************************