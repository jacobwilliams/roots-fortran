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

    character(len=*),parameter :: fmt  = '(A20,1X,A3,1X,A4,1X,A25,   1X,A25,   1X,A25,  1X,A5,1X,A5)' !! format for header
    character(len=*),parameter :: dfmt = '(A20,1X,I3,1X,I4,1X,E25.10,1X,E25.10,1X,E25.6,1X,I5,1X,I5)' !! format for results

    integer,parameter :: number_of_methods = 11 !! number of methods to test
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
        'chandrupatla   ', &
        'toms748        ' ] !! method names

    do imeth = 1, number_of_methods

        write(output_unit,fmt) &
            repeat('-',20),repeat('-',3),repeat('-',4),repeat('-',25),&
            repeat('-',25),repeat('-',25),repeat('-',5),repeat('-',5)

        write(output_unit,fmt) &
            'method','function','n','error','x','f','evals','iflag'

        write(output_unit,fmt) &
            repeat('-',20),repeat('-',3),repeat('-',4),repeat('-',25),&
            repeat('-',25),repeat('-',25),repeat('-',5),repeat('-',5)

        ! run the tests:
        nprob = 101; n = 1;   root = asin(2.0_wp/3.0_wp);    call test()
        nprob = 102; n = 1;   root = 0.5671432904_wp;        call test()
        nprob = 103; n = 1;   root = 0.8041330975_wp;        call test()
        nprob = 104; n = 1;   root = 3.0_wp;                 call test()
        nprob = 105; n = 1;   root = -0.6294464841_wp;       call test()
        nprob = 106; n = 1;   root = 2.0945514815_wp;        call test()
        nprob = 107; n = 1;   root = 1.0_wp;                 call test()
        nprob = 108; n = 1;   root = 1.0_wp;                 call test()
        nprob = 109; n = 1;   root = 1.0_wp;                 call test()
        nprob = 110; n = 1;   root = 1.0_wp;                 call test()
        nprob = 111; n = 1;   root = 1.0_wp;                 call test()
        nprob = 112; n = 1;   root = 0.0_wp;                 call test()
        nprob = 113; n = 1;   root = 0.0_wp;                 call test()
        nprob = 114; n = 1;   root = 0.677712926329_wp;      call test()
        nprob = 115; n = 1;   root = -0.640083696085_wp;     call test()
        nprob = 116; n = 20;  root = 0.05_wp;                call test()
        nprob = 117; n = 1;   root = 0.0_wp;                 call test()
        nprob = 1;  n = 1;    root = 1.8954942670340_wp;     call test()
        nprob = 2;  n = 1;    root = 3.0229153472731_wp;     call test()
        nprob = 3;  n = 1;    root = 6.6837535608081_wp;     call test()
        nprob = 4;  n = 1;    root = 11.238701655002_wp;     call test()
        nprob = 5;  n = 1;    root = 19.676000080623_wp;     call test()
        nprob = 6;  n = 1;    root = 29.828227326505_wp;     call test()
        nprob = 7;  n = 1;    root = 41.906116195289_wp;     call test()
        nprob = 8;  n = 1;    root = 55.953595800143_wp;     call test()
        nprob = 9;  n = 1;    root = 71.985665586588_wp;     call test()
        nprob = 10; n = 1;    root = 90.008868539167_wp;     call test()
        nprob = 11; n = 1;    root = 110.02653274833_wp;     call test()
        nprob = 12; n = 1;    root = 0.0_wp;                 call test()
        nprob = 13; n = 1;    root = 0.0_wp;                 call test()
        nprob = 14; n = 1;    root = 0.0_wp;                 call test()
        nprob = 15; n = 4;    root = 0.66874030497642_wp;    call test()
        nprob = 15; n = 6;    root = 0.76472449133173_wp;    call test()
        nprob = 15; n = 8;    root = 0.81776543395794_wp;    call test()
        nprob = 15; n = 10;   root = 0.85133992252078_wp;    call test()
        nprob = 15; n = 12;   root = 0.87448527222117_wp;    call test()
        nprob = 16; n = 4;    root = 1.0000000000000_wp;     call test()
        nprob = 16; n = 6;    root = 1.0000000000000_wp;     call test()
        nprob = 16; n = 8;    root = 1.0000000000000_wp;     call test()
        nprob = 16; n = 10;   root = 1.0000000000000_wp;     call test()
        nprob = 16; n = 12;   root = 1.0000000000000_wp;     call test()
        nprob = 17; n = 8;    root = 1.0000000000000_wp;     call test()
        nprob = 17; n = 10;   root = 1.0000000000000_wp;     call test()
        nprob = 17; n = 12;   root = 1.0000000000000_wp;     call test()
        nprob = 17; n = 14;   root = 1.0000000000000_wp;     call test()
        nprob = 18; n = 1;    root = 0.52359877559830_wp;    call test()
        nprob = 19; n = 1;    root = 0.42247770964124_wp;    call test()
        nprob = 19; n = 2;    root = 0.30669941048320_wp;    call test()
        nprob = 19; n = 3;    root = 0.22370545765466_wp;    call test()
        nprob = 19; n = 4;    root = 0.17171914751951_wp;    call test()
        nprob = 19; n = 5;    root = 0.13825715505682_wp;    call test()
        nprob = 19; n = 20;   root = 3.4657359020854e-02_wp; call test()
        nprob = 19; n = 40;   root = 1.7328679513999e-02_wp; call test()
        nprob = 19; n = 60;   root = 1.1552453009332e-02_wp; call test()
        nprob = 19; n = 80;   root = 8.6643397569993e-03_wp; call test()
        nprob = 19; n = 100;  root = 6.9314718055995e-03_wp; call test()
        nprob = 20; n = 5;    root = 3.8402551840622e-02_wp; call test()
        nprob = 20; n = 10;   root = 9.9000099980005e-03_wp; call test()
        nprob = 20; n = 20;   root = 2.4937500390620e-03_wp; call test()
        nprob = 21; n = 2;    root = 0.50000000000000_wp;    call test()
        nprob = 21; n = 5;    root = 0.34595481584824_wp;    call test()
        nprob = 21; n = 10;   root = 0.24512233375331_wp;    call test()
        nprob = 21; n = 15;   root = 0.19554762353657_wp;    call test()
        nprob = 21; n = 20;   root = 0.16492095727644_wp;    call test()
        nprob = 22; n = 1;    root = 0.27550804099948_wp;    call test()
        nprob = 22; n = 2;    root = 0.13775402049974_wp;    call test()
        nprob = 22; n = 4;    root = 1.0305283778156e-02_wp; call test()
        nprob = 22; n = 5;    root = 3.6171081789041e-03_wp; call test()
        nprob = 22; n = 8;    root = 4.1087291849640e-04_wp; call test()
        nprob = 22; n = 15;   root = 2.5989575892908e-05_wp; call test()
        nprob = 22; n = 20;   root = 7.6685951221853e-06_wp; call test()
        nprob = 23; n = 1;    root = 0.40105813754155_wp;    call test()
        nprob = 23; n = 5;    root = 0.51615351875793_wp;    call test()
        nprob = 23; n = 10;   root = 0.53952222690842_wp;    call test()
        nprob = 23; n = 15;   root = 0.54818229434066_wp;    call test()
        nprob = 23; n = 20;   root = 0.55270466667849_wp;    call test()
        nprob = 24; n = 2;    root = 0.50000000000000_wp;    call test()
        nprob = 24; n = 5;    root = 0.20000000000000_wp;    call test()
        nprob = 24; n = 15;   root = 6.6666666666667e-02_wp; call test()
        nprob = 24; n = 20;   root = 5.0000000000000e-02_wp; call test()
        nprob = 25; n = 2;    root = 2.0000000000000_wp;     call test()
        nprob = 25; n = 3;    root = 3.0000000000000_wp;     call test()
        nprob = 25; n = 4;    root = 4.0000000000000_wp;     call test()
        nprob = 25; n = 5;    root = 5.0000000000000_wp;     call test()
        nprob = 25; n = 6;    root = 6.0000000000000_wp;     call test()
        nprob = 25; n = 7;    root = 7.0000000000000_wp;     call test()
        nprob = 25; n = 9;    root = 9.0000000000000_wp;     call test()
        nprob = 25; n = 11;   root = 11.000000000000_wp;     call test()
        nprob = 25; n = 13;   root = 13.000000000000_wp;     call test()
        nprob = 25; n = 15;   root = 15.000000000000_wp;     call test()
        nprob = 25; n = 17;   root = 17.000000000000_wp;     call test()
        nprob = 25; n = 19;   root = 19.000000000000_wp;     call test()
        nprob = 25; n = 21;   root = 21.000000000000_wp;     call test()
        nprob = 25; n = 23;   root = 23.000000000000_wp;     call test()
        nprob = 25; n = 25;   root = 25.000000000000_wp;     call test()
        nprob = 25; n = 27;   root = 27.000000000000_wp;     call test()
        nprob = 25; n = 29;   root = 29.000000000000_wp;     call test()
        nprob = 25; n = 31;   root = 31.000000000000_wp;     call test()
        nprob = 25; n = 33;   root = 33.000000000000_wp;     call test()
        nprob = 26; n = 1;    root = 2.2317679157465e-02_wp; call test()
        nprob = 27; n = 1;    root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 2;    root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 3;    root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 4;    root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 5;    root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 6;    root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 7;    root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 8;    root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 9;    root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 10;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 11;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 12;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 13;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 14;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 15;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 16;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 17;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 18;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 19;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 20;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 21;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 22;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 23;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 24;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 25;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 26;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 27;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 28;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 29;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 30;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 31;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 32;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 33;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 34;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 35;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 36;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 37;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 38;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 39;   root = 0.62380651896161_wp;    call test()
        nprob = 27; n = 40;   root = 0.62380651896161_wp;    call test()
        nprob = 28; n = 20;   root = 5.9051305594220e-05_wp; call test()
        nprob = 28; n = 21;   root = 5.6367155339937e-05_wp; call test()
        nprob = 28; n = 22;   root = 5.3916409455592e-05_wp; call test()
        nprob = 28; n = 23;   root = 5.1669892394942e-05_wp; call test()
        nprob = 28; n = 24;   root = 4.9603096699145e-05_wp; call test()
        nprob = 28; n = 25;   root = 4.7695285287639e-05_wp; call test()
        nprob = 28; n = 26;   root = 4.5928793239949e-05_wp; call test()
        nprob = 28; n = 27;   root = 4.4288479195665e-05_wp; call test()
        nprob = 28; n = 28;   root = 4.2761290257883e-05_wp; call test()
        nprob = 28; n = 29;   root = 4.1335913915954e-05_wp; call test()
        nprob = 28; n = 30;   root = 4.0002497338020e-05_wp; call test()
        nprob = 28; n = 31;   root = 3.8752419296207e-05_wp; call test()
        nprob = 28; n = 32;   root = 3.7578103559958e-05_wp; call test()
        nprob = 28; n = 33;   root = 3.6472865219959e-05_wp; call test()
        nprob = 28; n = 34;   root = 3.5430783356532e-05_wp; call test()
        nprob = 28; n = 35;   root = 3.4446594929961e-05_wp; call test()
        nprob = 28; n = 36;   root = 3.3515605877800e-05_wp; call test()
        nprob = 28; n = 37;   root = 3.2633616249437e-05_wp; call test()
        nprob = 28; n = 38;   root = 3.1796856858426e-05_wp; call test()
        nprob = 28; n = 39;   root = 3.1001935436965e-05_wp; call test()
        nprob = 28; n = 40;   root = 3.0245790670210e-05_wp; call test()
        nprob = 28; n = 100;  root = 1.2277994232462e-05_wp; call test()
        nprob = 28; n = 200;  root = 6.1695393904409e-06_wp; call test()
        nprob = 28; n = 300;  root = 4.1198585298293e-06_wp; call test()
        nprob = 28; n = 400;  root = 3.0924623877272e-06_wp; call test()
        nprob = 28; n = 500;  root = 2.4752044261050e-06_wp; call test()
        nprob = 28; n = 600;  root = 2.0633567678513e-06_wp; call test()
        nprob = 28; n = 700;  root = 1.7690120078154e-06_wp; call test()
        nprob = 28; n = 800;  root = 1.5481615698859e-06_wp; call test()
        nprob = 28; n = 900;  root = 1.3763345366022e-06_wp; call test()
        nprob = 28; n = 1000; root = 1.2388385788997e-06_wp; call test()

    end do

    contains
!*****************************************************************************************

    subroutine test()

        implicit none

        real(wp) :: ax,bx,xzero,fzero,error
        integer :: iflag
        character(len=1000) :: line

        ifunc = 0 ! reset func evals counter
        call get_bounds(ax, bx)
        call root_scalar(methods(imeth),func,ax,bx,xzero,fzero,iflag,&
                            atol = 1.0e-16_wp, rtol = 1.0e-12_wp, maxiter = 1000)

        error = xzero-root
        write(line, dfmt) trim(methods(imeth)),nprob,n,xzero-root,xzero,fzero,ifunc,iflag

        if (abs(fzero) <= 1.0e-9_wp .or. abs(error) <= 1.0e-9_wp) then
            write(output_unit, '(A)') trim(line)
        else
            write(output_unit, '(A)') colorize(trim(line), color_fg='red') ! failed case
        end if

    end subroutine test

    subroutine get_bounds(a, b)

        implicit none

        real(wp),intent(out) :: a, b    !! bounds

        real(wp),parameter :: pi = acos(-1.0_wp)

        select case (nprob)
        case (1)
            a = pi/2.0_wp
            b = pi
        case (2)
            a = 1.0_wp + 1.0e-9_wp
            b = (2.0_wp)*(2.0_wp) - 1.0e-9_wp
        case (3)
            a = (2.0_wp)*(2.0_wp) + 1.0e-9_wp
            b = (3.0_wp)*(3.0_wp) - 1.0e-9_wp
        case (4)
            a = (3.0_wp)*(3.0_wp) + 1.0e-9_wp
            b = (4.0_wp)*(4.0_wp) - 1.0e-9_wp
        case (5)
            a = (4.0_wp)*(4.0_wp) + 1.0e-9_wp
            b = (5.0_wp)*(5.0_wp) - 1.0e-9_wp
        case (6)
            a = (5.0_wp)*(5.0_wp) + 1.0e-9_wp
            b = (6.0_wp)*(6.0_wp) - 1.0e-9_wp
        case (7)
            a = (6.0_wp)*(6.0_wp) + 1.0e-9_wp
            b = (7.0_wp)*(7.0_wp) - 1.0e-9_wp
        case (8)
            a = (7.0_wp)*(7.0_wp) + 1.0e-9_wp
            b = (8.0_wp)*(8.0_wp) - 1.0e-9_wp
        case (9)
            a = (8.0_wp)*(8.0_wp) + 1.0e-9_wp
            b = (9.0_wp)*(9.0_wp) - 1.0e-9_wp
        case (10)
            a = (9.0_wp)*(9.0_wp) + 1.0e-9_wp
            b = (10.0_wp)*(10.0_wp) - 1.0e-9_wp
        case (11)
            a = (10.0_wp)*(10.0_wp) + 1.0e-9_wp
            b = (11.0_wp)*(11.0_wp) - 1.0e-9_wp
        case (12, 13, 14)
            a = -9.0_wp
            b = 31.0_wp
        case (15, 16)
            a = 0.0_wp
            b = 5.0_wp
        case (17)
            a = -0.95_wp
            b = 4.05_wp
        case (18)
            a = 0.0_wp
            b = 1.5_wp
        case (19, 20, 21, 22, 23)
            a = 0.0_wp
            b = 1.0_wp
        case (24)
            a = 1.0e-2_wp
            b = 1.0_wp
        case (25)
            a = 1.0_wp
            b = 100.0_wp
        case (26)
            a = -1.0_wp
            b = 4.0_wp
        case (27)
            a = -10000.0_wp
            b = pi/2.0_wp
        case (28)
            a = -10000.0_wp
            b = 1.0e-4_wp
        case (101)
            a = 0.0_wp
            b = 1.0_wp
        case (102)
            a = -1.0_wp
            b = 1.0_wp
        case (103)
            a = 0.1_wp
            b = 0.9_wp
        case (104)
            a = 2.8_wp
            b = 3.1_wp
        case (105)
            a = -1.3_wp
            b = -0.5_wp
        case (106)
            a = 2.0_wp
            b = 3.0_wp
        case (107)
            a = 0.5_wp
            b = 2.0_wp
        case (108)
            a = -10.0_wp
            b = 10.0_wp
        case (109)
            a = 0.0_wp
            b = 10.0_wp
        case (110)
            a = 0.0_wp
            b = 10.0_wp
        case (111)
            a = 0.0_wp
            b = 10.0_wp
        case (112)
            !From: http://www.mth.pdx.edu/~daescu/mth451_551/Muller_bisection.pdf
            a = -2.0_wp
            b = 3.0_wp
        case (113)
            !From: http://www.mth.pdx.edu/~daescu/mth451_551/Muller_bisection.pdf
            a = 0.0_wp     ! x1 is the location of the root for this case
            b = 3.0_wp
        case (114)
            a = 0.0_wp
            b = 0.7_wp
        case (115)
            a = -0.7_wp
            b = 0.4_wp
        case (116)
            !case 11 from "Algorithm 748: Enclosing Zeros of Continuous Functions" (n=20)
            a = 0.01_wp
            b = 1.0_wp
        case (117)
            a = -10.0_wp
            b = 5.0_wp

        case default
            error stop 'invalid case'
        end select

    end subroutine get_bounds

    function func(x) result(f)

        implicit none

        real(wp),intent(in) :: x
        real(wp) :: f

        integer :: i !! counter
        real(wp) :: dn, di

        dn = real(n, wp)

        ifunc = ifunc + 1

        select case (nprob)
        case (1)
            f = sin(x) - x/2.0_wp
        case (2:11)
            f = 0.0_wp
            do i = 1, 20
                di = real(i,wp)
                f = f + ((2.0_wp*di - 5.0_wp)**2)/(x - di*di)**3
            end do
            f = -2.0_wp*f
        case (12)
            f = -40.0_wp*x*exp(-1.0_wp*x)
        case (13)
            f = -100.0_wp*x*exp(-2.0_wp*x)
        case (14)
            f = -200.0_wp*x*exp(-3.0_wp*x)
        case (15)
            f = x**n - 0.2_wp
        case (16, 17)
            f = x**n - 1.0_wp
        case (18)
            f = sin(x) - 0.5_wp
        case (19)
            f = 2.0_wp*x*exp(-dn) - 2.0_wp*exp(-dn*x) + 1.0_wp
        case (20)
            f = (1.0_wp + (1.0_wp - dn)**2)*x - (1.0_wp - dn*x)**2
        case (21)
            f = x**2 - (1.0_wp - x)**n
        case (22)
            f = (1.0_wp + (1.0_wp - dn)**4)*x - (1.0_wp - dn*x)**4
        case (23)
            f = (x - 1.0_wp)*exp(-dn*x) + x**n
        case (24)
            f = (dn*x - 1.0_wp)/((dn - 1.0_wp)*x)
        case (25)
            f = x**(1.0_wp/dn) - dn**(1.0_wp/dn)
        case (26)
            if (x == 0.0_wp) then
                f = 0.0_wp
            else
                f = x/exp(1.0_wp/(x*x))
            end if
        case (27)
            if (x >= 0.0_wp) then
                f = (x/1.5_wp + sin(x) - 1.0_wp)*dn/20.0_wp
            else
                f = (-1.0_wp*dn)/20.0_wp
            end if
        case (28)
            if (x >= (1.0e-3_wp)*2.0_wp/(dn + 1.0_wp)) then
                f = exp(1.0_wp) - 1.859_wp
            elseif (x >= 0.0_wp) then
                f = exp((dn + 1.0_wp)*0.5_wp*x*1000.0_wp) - 1.859_wp
            else
                f = -0.859_wp
            end if
        case (101)
           ! Gottlieb's paper [table 1 & 2]
           f = 3.0_wp * sin(x) - 2.0_wp
        case (102)
           ! Gottlieb's paper [table 1 & 2]
           f = x*exp(x) - 1.0_wp
        case (103)
           ! Gottlieb's paper [table 1 & 2]
           f = 11.0_wp * x**11 - 1.0_wp
        case (104)
           ! Gottlieb's paper [table 1 & 2]
           f = exp(x**2 + 7.0_wp*x - 30.0_wp) - 1.0_wp
        case (105)
           ! Gottlieb's paper [table 1 & 2]
           f = 1.0_wp/x - sin(x) + 1.0_wp
        case (106)
           ! Gottlieb's paper [table 1 & 2]
           f = x**3 - 2.0_wp*x - 5.0_wp
        case (107)
           ! Gottlieb's paper [table 1 & 2]
           f = 1.0_wp/x - 1.0_wp
        case (108)
           f = x**3 - 5.0_wp*x**2 + 12.0_wp*x - 8.0_wp
        case (109)
           f = x**5 - 61.0_wp*x**4 + 1368.0_wp*x**3 - 13548.0_wp*x**2 + 54256.0_wp*x - 42016.0_wp ! x in [0, 10]
        case (110)
           f = x**5-61.0_wp*x**4+6801.0_wp/5_wp*x**3-66531.0_wp/5_wp*x**2+5205601.0_wp/100.0_wp*x-4005001.0_wp/100.0_wp
        case (111)
           f = (x-1.0_wp)**5
        case (112)
           f = atan(x)
        case (113)
           f = exp(x) - 2.0_wp*x - 1.0_wp
        case (114)
           f = 3.0_wp * x*sin(x*20.0_wp)*cos(x) + x - 2.0_wp
        case (115)
           f = 3.0_wp * x**2*sin(x*20.0_wp) + cos(x*2.0_wp)
        case (116)
           ! case 11 from "Algorithm 748: Enclosing Zeros of Continuous Functions" (n=20)
           f = (real(n,wp)*x-1.0_wp)/((real(n,wp)-1.0_wp)*x)
        case (117)
           f = sin(x*10.0_wp) + sin(x*2.0_wp) + atan(x/3.0_wp)*x**2 + exp(x) - 1.0_wp

        case default
            error stop 'invalid case'
        end select

    end function func

!*****************************************************************************************
    end program root_tests
!*****************************************************************************************