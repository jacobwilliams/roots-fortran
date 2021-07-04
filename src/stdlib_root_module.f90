!*****************************************************************************************
!>
!  Root solver methods for:
!  * Bracked interval
!  * Without derivatives

    module stdlib_root_module

    use iso_fortran_env, only: wp => real64, ip => int32
    use iso_fortran_env, only: error_unit

    implicit none

    type,abstract,public :: root_solver
    !! abstract class for the root solver methods
    private
    procedure(func),pointer :: f => null()  !! user function to find the root of
    real(wp) :: ftol = 0.0_wp   !! absolute tolerance for `f=0`
    contains
    private
    procedure,public :: solve !! main routine for finding the root
    procedure(root_f),deferred :: find_root !! root solver function
    procedure :: get_fa_fb
    end type root_solver

    type,extends(root_solver),public :: zeroin_solver
    !! Classic zeroin root solver
    private
    real(wp) :: tol = 1.0e-6_wp     !! desired length of the interval of
                                    !! uncertainty of the final result (>=0)
    contains
    private
    procedure,public :: find_root => zeroin
    end type zeroin_solver

    type,extends(root_solver),public :: bisection_solver
    !! Classic bisection root solver
    private
    real(wp) :: rtol = 1.0e-6_wp    !! relative tol for x
    real(wp) :: atol = 1.0e-12_wp   !! absolute tol for x
    integer  :: maxiter = 2000      !! maximum number of iterations
    contains
    private
    procedure,public :: find_root => bisection
    end type bisection_solver

    type,extends(root_solver),public :: anderson_bjorck_solver
    !! anderson bjorck root solver
    private
    real(wp) :: rtol = 1.0e-6_wp    !! relative tol for x
    real(wp) :: atol = 1.0e-12_wp   !! absolute tol for x
    integer  :: maxiter = 2000      !! maximum number of iterations
    contains
    private
    procedure,public :: find_root => anderson_bjorck
    end type anderson_bjorck_solver

    type,extends(root_solver),public :: ridders_solver
    !! anderson bjorck root solver
    private
    real(wp) :: tol = 1.0e-6_wp     !! relative tol for x
    integer  :: maxiter = 2000      !! maximum number of iterations
    contains
    private
    procedure,public :: find_root => ridders
    end type ridders_solver

    type,extends(root_solver),public :: pegasus_solver
    !! anderson bjorck root solver
    private
    real(wp) :: rtol = 1.0e-6_wp    !! relative tol for x
    real(wp) :: atol = 1.0e-12_wp   !! absolute tol for x
    integer  :: maxiter = 2000      !! maximum number of iterations
    contains
    private
    procedure,public :: find_root => pegasus
    end type pegasus_solver

    type,extends(root_solver),public :: bdqrf_solver
    !! anderson bjorck root solver
    private
    integer  :: maxiter = 2000      !! maximum number of iterations
    contains
    private
    procedure,public :: find_root => bdqrf
    end type bdqrf_solver

    type,extends(root_solver),public :: muller_solver
    !! anderson bjorck root solver
    private
    real(wp) :: rtol = 1.0e-6_wp    !! relative tol for x
    real(wp) :: atol = 1.0e-12_wp   !! absolute tol for x
    integer  :: maxiter = 2000      !! maximum number of iterations
    contains
    private
    procedure,public :: find_root => muller
    end type muller_solver

    abstract interface
        function func(me,x) result(f)
            !! Interface to the function to be minimized
            !! (Object-oriented version).
            !! It should evaluate f(x) for any x in the interval (ax,bx)
            import :: root_solver, wp
            implicit none
            class(root_solver),intent(inout) :: me
            real(wp),intent(in) :: x
            real(wp) :: f
        end function func
        function func2(x) result(f)
            !! Interface to the function to be minimized
            !! (Functional version).
            !! It should evaluate f(x) for any x in the interval (ax,bx)
            import :: wp
            implicit none
            real(wp),intent(in) :: x
            real(wp) :: f
        end function func2
        subroutine root_f(me,ax,bx,fax,fbx,xzero,fzero,iflag)
            !! Root solver function interface
            import :: root_solver, wp, ip
            implicit none
            class(root_solver),intent(inout) :: me
            real(wp),intent(in)       :: ax
            real(wp),intent(in)       :: bx
            real(wp),intent(in)       :: fax
            real(wp),intent(in)       :: fbx
            real(wp),intent(out)      :: xzero
            real(wp),intent(out)      :: fzero
            integer(ip),intent(out)   :: iflag
        end subroutine root_f
    end interface

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Non-object-oriented wrapper.

    subroutine root_scalar(method,fun,ax,bx,xzero,fzero,iflag,ftol,tol,rtol,atol,maxiter,fax,fbx)

    implicit none

    character(len=*),intent(in)   :: method   !! the method to use
    procedure(func2)              :: fun      !! user function to find the root of
    real(wp),intent(in)           :: ax       !! left endpoint of initial interval
    real(wp),intent(in)           :: bx       !! right endpoint of initial interval
    real(wp),intent(out)          :: xzero    !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)          :: fzero    !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)           :: iflag    !! status flag (`-1`=error, `0`=root found, `-999`=invalid method)
    real(wp),intent(in),optional  :: ftol     !! absolute tolerance for `f=0`
    real(wp),intent(in),optional  :: tol      !! desired length of the interval of uncertainty of the final result (>=0)
    real(wp),intent(in),optional  :: rtol     !! relative tol for x
    real(wp),intent(in),optional  :: atol     !! absolute tol for x
    integer,intent(in),optional   :: maxiter  !! maximum number of iterations
    real(wp),intent(in),optional  :: fax      !! if `f(ax)` is already known, it can be input here
    real(wp),intent(in),optional  :: fbx      !! if `f(ax)` is already known, it can be input here

    select case (lowercase(method))

    case('zeroin')

        block
            type(zeroin_solver) :: s
            if (present(ftol)) s%ftol = ftol
            if (present(tol))  s%tol  = tol
            s%f => func_wrapper
            call s%solve(ax,bx,xzero,fzero,iflag,fax,fbx)
        end block

    case('bisection')

        block
            type(bisection_solver) :: s
            if (present(ftol))    s%ftol    = ftol
            if (present(rtol))    s%rtol    = rtol
            if (present(atol))    s%atol    = atol
            if (present(maxiter)) s%maxiter = maxiter
            s%f => func_wrapper
            call s%solve(ax,bx,xzero,fzero,iflag,fax,fbx)
        end block

    case('anderson_bjorck', 'anderson-bjorck')

        block
            type(anderson_bjorck_solver) :: s
            if (present(ftol))    s%ftol    = ftol
            if (present(rtol))    s%rtol    = rtol
            if (present(atol))    s%atol    = atol
            if (present(maxiter)) s%maxiter = maxiter
            s%f => func_wrapper
            call s%solve(ax,bx,xzero,fzero,iflag,fax,fbx)
        end block

    case('ridders')

        block
            type(ridders_solver) :: s
            if (present(ftol))    s%ftol    = ftol
            if (present(tol))     s%tol     = tol
            if (present(maxiter)) s%maxiter = maxiter
            s%f => func_wrapper
            call s%solve(ax,bx,xzero,fzero,iflag,fax,fbx)
        end block

    case('pegasus')

        block
            type(pegasus_solver) :: s
            if (present(ftol))    s%ftol    = ftol
            if (present(rtol))    s%rtol    = rtol
            if (present(atol))    s%atol    = atol
            if (present(maxiter)) s%maxiter = maxiter
            s%f => func_wrapper
            call s%solve(ax,bx,xzero,fzero,iflag,fax,fbx)
        end block

    case('bdqrf')

        block
            type(bdqrf_solver) :: s
            if (present(ftol))    s%ftol    = ftol
            if (present(maxiter)) s%maxiter = maxiter
            s%f => func_wrapper
            call s%solve(ax,bx,xzero,fzero,iflag,fax,fbx)
        end block

    case('muller')

        block
            type(muller_solver) :: s
            if (present(ftol))    s%ftol    = ftol
            if (present(rtol))    s%rtol    = rtol
            if (present(atol))    s%atol    = atol
            if (present(maxiter)) s%maxiter = maxiter
            s%f => func_wrapper
            call s%solve(ax,bx,xzero,fzero,iflag,fax,fbx)
        end block

    case default
        ! invalid method
        iflag = -999
    end select

    contains

        function func_wrapper(me,x) result(f)
            implicit none
            class(root_solver),intent(inout) :: me
            real(wp),intent(in) :: x
            real(wp) :: f
            f = fun(x)
        end function func_wrapper

        pure function lowercase(str) result(s_lower)

        implicit none

        character(len=*),intent(in) :: str      !! input string
        character(len=(len(str)))   :: s_lower  !! lowercase version of the string

        integer :: i  !! counter
        integer :: j  !! index of uppercase character

        character(len=*),parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' !! uppercase characters
        character(len=*),parameter :: lower = 'abcdefghijklmnopqrstuvwxyz' !! lowercase characters

        s_lower = str

        do i = 1, len_trim(str)
            j = index(upper,s_lower(i:i))
            if (j>0) s_lower(i:i) = lower(j:j)
        end do

        end function lowercase

    end subroutine root_scalar
!*****************************************************************************************

!*****************************************************************************************
!>
!  Main wrapper routine for all the methods.

    subroutine solve(me,ax,bx,xzero,fzero,iflag,fax,fbx)

    implicit none

    class(root_solver),intent(inout) :: me
    real(wp),intent(in)              :: ax      !! left endpoint of initial interval
    real(wp),intent(in)              :: bx      !! right endpoint of initial interval
    real(wp),intent(out)             :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)             :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)              :: iflag   !! status flag (`-1`=error, `0`=root found, `-4`=ax must be /= bx)
    real(wp),intent(in),optional     :: fax     !! if `f(ax)` is already known, it can be input here
    real(wp),intent(in),optional     :: fbx     !! if `f(ax)` is already known, it can be input here

    real(wp) :: fa, fb

    if (ax==bx) then
        ! ax must be /= bx
        iflag = -4
        xzero = ax  ! just to return something
        fzero = fa  !
    else

        call me%get_fa_fb(ax,bx,fax,fbx,fa,fb)

        ! check trivial cases first:
        if (abs(fa)<=me%ftol) then

            iflag = 0
            xzero = ax
            fzero = fa

        elseif (abs(fb)<=me%ftol) then

            iflag = 0
            xzero = bx
            fzero = fb

        elseif (fa*fb>0.0_wp) then

            ! f(ax) and f(bx) do not have different signs
            iflag = -1
            xzero = ax  ! just to return something
            fzero = fa  !

        else

            ! call the root solver.
            ! make sure order is correct.
            if (ax<bx) then
                call me%find_root(ax,bx,fa,fb,xzero,fzero,iflag)
            else
                call me%find_root(bx,ax,fb,fa,xzero,fzero,iflag)
            end if

        end if

    end if

    end subroutine solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the function values at `ax` and `bx` to start the root finding algorithm.

    subroutine get_fa_fb(me,ax,bx,fax,fbx,fa,fb)

    implicit none

    class(root_solver),intent(inout) :: me
    real(wp),intent(in)              :: ax      !! left endpoint of initial interval
    real(wp),intent(in)              :: bx      !! right endpoint of initial interval
    real(wp),intent(in),optional     :: fax     !! if `f(ax)` is already known, it can be input here
    real(wp),intent(in),optional     :: fbx     !! if `f(ax)` is already known, it can be input here
    real(wp),intent(out)             :: fa      !! `f(ax)` to use
    real(wp),intent(out)             :: fb      !! `f(ax)` to use

    if (present(fax)) then
        fa = fax
    else
        fa = me%f(ax)
    end if

    if (present(fbx)) then
        fb = fbx
    else
        fb = me%f(bx)
    end if

    end subroutine get_fa_fb
!*****************************************************************************************

!*****************************************************************************************
!>
!  Find a zero of the function \( f(x) \) in the given interval
!  \( [a_x,b_x] \) to within a tolerance \( 4 \epsilon |x| + tol \),
!  where \( \epsilon \) is the relative machine precision defined as
!  the smallest representable number such that \( 1.0 + \epsilon > 1.0 \).
!
!  It is assumed that \( f(a_x) \) and \( f(b_x) \) have opposite signs.
!
!### References
!  * R. P. Brent, "[An algorithm with guaranteed convergence for
!    finding a zero of a function](http://maths-people.anu.edu.au/~brent/pd/rpb005.pdf)",
!    The Computer Journal, Vol 14, No. 4., 1971.
!  * R. P. Brent, "[Algorithms for minimization without derivatives](http://maths-people.anu.edu.au/~brent/pub/pub011.html)",
!    Prentice-Hall, Inc., 1973.
!
!### See also
!  1. [zeroin.f](http://www.netlib.org/go/zeroin.f) from Netlib

    subroutine zeroin(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(zeroin_solver),intent(inout) :: me
    real(wp),intent(in)    :: ax    !! left endpoint of initial interval
    real(wp),intent(in)    :: bx    !! right endpoint of initial interval
    real(wp),intent(in)    :: fax   !! `f(ax)`
    real(wp),intent(in)    :: fbx   !! `f(ax)`
    real(wp),intent(out)   :: xzero !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)   :: fzero !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)    :: iflag !! status flag (`0`=root found)

    real(wp),parameter :: eps = epsilon(1.0_wp)  !! d1mach(4) in original code
    real(wp) :: a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s

    ! initialize:
    iflag = 0
    tol1  = eps+1.0_wp
    a     = ax
    b     = bx
    fa    = fax
    fb    = fbx
    c     = a
    fc    = fa
    d     = b-a
    e     = d

    do

        if (abs(fc)<abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
        end if

        tol1 = 2.0_wp*eps*abs(b)+0.5_wp*me%tol
        xm = 0.5_wp*(c-b)
        if ((abs(xm)<=tol1).or.(fb==0.0_wp)) exit

        ! see if a bisection is forced
        if ((abs(e)>=tol1).and.(abs(fa)>abs(fb))) then
            s=fb/fa
            if (a/=c) then
                ! inverse quadratic interpolation
                q=fa/fc
                r=fb/fc
                p=s*(2.0_wp*xm*q*(q-r)-(b-a)*(r-1.0_wp))
                q=(q-1.0_wp)*(r-1.0_wp)*(s-1.0_wp)
            else
                ! linear interpolation
                p=2.0_wp*xm*s
                q=1.0_wp-s
            end if
            if (p<=0.0_wp) then
                p=-p
            else
                q=-q
            end if
            s=e
            e=d
            if (((2.0_wp*p)>=(3.0_wp*xm*q-abs(tol1*q))) .or. (p>=abs(0.5_wp*s*q))) then
                d=xm
                e=d
            else
                d=p/q
            end if
        else
            d=xm
            e=d
        end if

        a=b
        fa=fb
        if (abs(d)<=tol1) then
            if (xm<=0.0_wp) then
                b=b-tol1
            else
                b=b+tol1
            end if
        else
            b=b+d
        end if
        fb=me%f(b)
        if (abs(fb)<=me%ftol) exit  ! absolute convergence in f
        if ((fb*(fc/abs(fc)))>0.0_wp) then
            c=a
            fc=fa
            d=b-a
            e=d
        end if

    end do

    xzero = b
    fzero = fb

    end subroutine zeroin
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the zero of the function f(x) in the interval ax,bx using the bisection method.
!
!  It is assumed that \( f(a_x) \) and \( f(b_x) \) have opposite signs.
!
!### See also
!  * G.E. Mullges & F. Uhlig, "Numerical Algorithms with Fortran",
!    Springer, 1996. Section 2.8.1, p 32-34.

    subroutine bisection(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(bisection_solver),intent(inout) :: me
    real(wp),intent(in)    :: ax      !! left endpoint of initial interval
    real(wp),intent(in)    :: bx      !! right endpoint of initial interval
    real(wp),intent(in)    :: fax     !! `f(ax)`
    real(wp),intent(in)    :: fbx     !! `f(ax)`
    real(wp),intent(out)   :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)   :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)    :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    real(wp) :: x1,x2,x3,f1,f2,f3
    integer :: i !! iteration counter
    logical :: root_found !! convergence in x

    ! initialize:
    iflag = 0
    x1    = ax
    x2    = bx
    f1    = fax
    f2    = fbx

    ! main loop
    do i=1,me%maxiter

        ! bisection of the inclusion interval:
        !  x1------x3------x2
        x3 = x2 + (x1 - x2) / 2.0_wp

        ! calculate the new function value:
        f3 = me%f(x3)
        ! check for root:
        if (abs(f3)<=me%ftol) then
            xzero = x3
            fzero = f3
            return
        end if

        ! determine new inclusion interval:
        if (f2*f3<0.0_wp) then
            ! root lies between x2 and x3
            x1 = x3
            x2 = x2
            f1 = f3
            f2 = f2
        else
            ! root lies between x1 and x3
            x2 = x3
            f2 = f3
        end if

        ! check for convergence:
        root_found = abs(x2-x1) <= abs(x2) * me%rtol + me%atol
        if (root_found .or. i==me%maxiter) then
            xzero = x2
            fzero = f2
            if (.not. root_found) iflag = -2  ! max iterations reached
            exit
        end if

    end do

    end subroutine bisection
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the zero of the function f(x) in the interval ax,bx using the Anderson-Bjorck method.
!
!  It is assumed that \( f(a_x) \) and \( f(b_x) \) have opposite signs.
!
!### See also
!  * G.E. Mullges & F. Uhlig, "Numerical Algorithms with Fortran",
!    Springer, 1996. Section 2.8.2, p 36.

    subroutine anderson_bjorck(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(anderson_bjorck_solver),intent(inout) :: me
    real(wp),intent(in)    :: ax      !! left endpoint of initial interval
    real(wp),intent(in)    :: bx      !! right endpoint of initial interval
    real(wp),intent(in)    :: fax     !! `f(ax)`
    real(wp),intent(in)    :: fbx     !! `f(ax)`
    real(wp),intent(out)   :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)   :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)    :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    integer :: i !! counter
    logical :: root_found !! convergence in x
    real(wp) :: x1,x2,x3,f1,f2,f3,s12,g

    ! initialize:
    iflag = 0
    x1    = ax
    x2    = bx
    f1    = fax
    f2    = fbx

    ! main loop:
    do i = 1,me%maxiter

        ! secant step:
        s12 = (f2 - f1) / (x2 - x1)

        ! intersection of this secant with the x-axis:
        x3 = x2 - f2 / s12

        ! calculate f3:
        f3 = me%f(x3)

        if (abs(f3)<=me%ftol)  then  ! f3 is a root
            xzero = x3
            fzero = f3
            exit
        end if

        ! determine a new inclusion interval:
        if (f2*f3<0.0_wp) then
            ! zero lies between x2 and x3
            x1 = x2
            x2 = x3
            f1 = f2
            f2 = f3
        else
            ! zero lies between x1 and x3
            g = 1.0_wp-f3/f2
            if (g<=0.0_wp) g = 0.5_wp
            x2 = x3
            f1 = g*f1
            f2 = f3
        end if

        ! check for convergence:
        root_found = abs(x2-x1) <= abs(x2)*me%rtol + me%atol
        if (root_found .or. i == me%maxiter) then
            xzero = x2
            fzero = f2
            if (.not. root_found) iflag = -2  ! max iterations reached
            exit
        end if

    end do

    end subroutine anderson_bjorck
!*****************************************************************************************

!*****************************************************************************************
!>
!  Ridders method to find a root of f(x).
!
!### See also
!  * Ridders, C., "A new algorithm for computing a single root of a real continuous function",
!    IEEE Trans. on Circuits and Systems, Vol 26, Issue 11.

    subroutine ridders(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(ridders_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached, `-3`=singularity in the algorithm)

    integer  :: i !! counter
    real(wp) :: fh,fl,fm,fnew,denom,xh,xl,xm,xnew

    ! initialize:
    iflag = 0
    fl    = fax
    fh    = fbx
    xl    = ax
    xh    = bx

    do i = 1, me%maxiter

        xm = (xl+xh)/2.0_wp
        fm = me%f(xm)
        if (abs(fm) <= me%ftol) then
            ! abs convergence in f
            xzero = xm
            fzero = fm
            exit
        end if

        denom = sqrt(fm**2-fl*fh)
        if (denom == 0.0_wp) then
            xzero = xm
            fzero = fm
            iflag = -3        ! can't proceed: denominator is zero
            exit
        end if

        xnew = xm+(xm-xl)*(sign(1.0_wp,fl-fh)*fm/denom)
        if (abs(xnew-xzero) <= me%tol) then
            ! relative convergence in x
            exit
        end if

        xzero = xnew
        fnew = me%f(xzero)
        fzero = fnew
        if (abs(fnew) <= me%ftol) then
            ! abs convergence in f
            exit
        end if

        ! to keep the root bracketed:
        if (sign(fm,fnew) /= fm) then
            xl = xm
            fl = fm
            xh = xzero
            fh = fnew
        else if (sign(fl,fnew) /= fl) then
            xh = xzero
            fh = fnew
        else if (sign(fh,fnew) /= fh) then
            xl = xzero
            fl = fnew
        end if

        if (abs(xh-xl) <= me%tol) then
            ! relative convergence in x
            exit
        else if (i == me%maxiter) then
            iflag = -2    ! max iterations exceeded
        end if

    end do

    end subroutine ridders
!*****************************************************************************************

!*****************************************************************************************
!>
!  Pegasus method to find a root of f(x).
!
!### See also
!  * G.E. Mullges & F. Uhlig, "Numerical Algorithms with Fortran",
!    Springer, 1996. Section 2.8.2, p 35.

    subroutine pegasus(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(pegasus_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    integer :: i !! counter
    real(wp) :: x1,x2,x3,f1,f2,f3,s12

    ! initialize:
    iflag = 0
    x1    = ax
    x2    = bx
    f1    = fax
    f2    = fbx

    ! main loop:
    do i = 1, me%maxiter

        s12 = (f2 - f1) / (x2 - x1) ! secant step
        x3  = x2 - f2 / s12         ! intersection of this secant with the x-axis
        f3  = me%f(x3)              ! calculate f3

        if (f3==0.0_wp)  then ! f3 is a root
            fzero = f3
            xzero = x3
            iflag = 0
            return
        end if

        ! determine a new inclusion interval:
        if (f2*f3<=0.0_wp) then
            x1 = x2
            f1 = f2
        else
            f1 = f1 * f2 / (f2 + f3)
        end if

        x2 = x3
        f2 = f3

        ! Check for break-off condition:
        if (abs(f2)<me%ftol) exit
        if (abs(x2-x1)<=abs(x2)*me%rtol + me%atol) exit
        if (i == me%maxiter) iflag = -2   ! max iterations exceeded

    end do

    fzero = f2
    xzero = x2

    end subroutine pegasus
!*****************************************************************************************

!*****************************************************************************************
!>
!  Bisected Direct Quadratic Regula Falsi (BDQRF) root solver method
!  to find the root of a 1D function.
!
!### See also
!  * R. G. Gottlieb, B. F. Thompson, "Bisected Direct Quadratic Regula Falsi",
!    Applied Mathematical Sciences, Vol. 4, 2010, no. 15, 709-718.

    subroutine bdqrf(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(bdqrf_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached, `-3`=limit of precision reached)

    real(wp) :: xdn,ydn,xup,yup,xlast,d,xm,ym,a,b,y2
    integer :: i !! counter

    ! initialize:
    iflag = 0
    xzero = ax
    fzero = fax
    y2    = fbx
    xlast = huge(1.0_wp)

    if (fzero<0.0_wp) then
        xdn = ax
        ydn = fzero
        xup = bx
        yup = y2
    else
        xup = ax
        yup = fzero
        xdn = bx
        ydn = y2
    end if

    ! main loop:
    do i = 1, me%maxiter

        d = (xup - xdn) / 2.0_wp
        xm = (xup + xdn) / 2.0_wp
        ym = me%f(xm)
        a = (yup + ydn - 2.0_wp*ym)/(2.0_wp*d**2)
        b = (yup - ydn)/(2.0_wp*d)
        xzero = xm - 2.0_wp*ym / (b * (1.0_wp + sqrt(1.0_wp - 4.0_wp*a*ym/b**2)))

        if (xzero==xlast) exit ! limit of computing precision has been reached.

        xlast = xzero
        fzero = me%f(xzero)

        if (abs(fzero)<=me%ftol) exit ! Convergence

        if (fzero>0.0_wp) then
            yup = fzero
            xup = xzero
            if (ym<0.0_wp) then
                ydn = ym
                xdn = xm
            end if
        else
            ydn = fzero
            xdn = xzero
            if (ym>0.0_wp) then
                yup = ym
                xup = xm
            end if
        end if

        if (i==me%maxiter) iflag = -2 ! maximum number of iterations

    end do

    end subroutine bdqrf
!*****************************************************************************************

!*****************************************************************************************
!>
!   Muller's method to find a real root of f(x).

    ! this one seems to fail for some cases ...

    subroutine muller (me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(muller_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached, `-3`=singularity in the algorithm)

    real(wp)  :: a,b,c,a2,d
    real(wp)  :: fminus,fplus,fxmid,fxnew,fxold
    real(wp)  :: x_ave,x_inc,xlast,xmid,xminus,xold,xplus
    integer   :: i !! iteration counter

    iflag = 0
    xzero = ax
    xold  = bx
    fxnew = fax
    fxold = fbx
    fzero = fxnew

    xmid  = (ax + bx) / 2.0_wp   ! pick a third point in the middle
    fxmid = me%f(xmid)
    if (abs(fxmid)<me%ftol) then
        xzero = xmid
        fzero = fxmid
        return
    end if

    ! main loop:
    do i = 1, me%maxiter

        if ( abs(fxnew) >= abs(fxmid) ) then
            call swap ( xzero, xmid )
            call swap ( fxnew, fxmid )
        end if

        xlast = xzero

        a = (xmid-xzero)*(fxold-fxnew)-(xold-xzero)*(fxmid-fxnew)
        b = (xold-xzero)**2*(fxmid-fxnew)-(xmid-xzero)**2*(fxold-fxnew)
        c = (xold-xzero)*(xmid-xzero)*(xold-xmid)*fxnew

        if ( a == 0.0_wp ) then
            iflag = -3
            exit
        end if

        xold = xmid
        xmid = xzero

        !  Apply the quadratic formula to get roots xplus and xminus.
        d = b**2 - 4.0_wp * a * c
        if ( d < 0.0_wp ) then
            d = 0.0_wp  ! to avoid complex roots
        else
            d = sqrt(d)
        end if
        a2 = 2.0_wp*a

        xplus  = xzero + ( - b + d ) / a2
        xminus = xzero + ( - b - d ) / a2

        fplus  = me%f(xplus)
        if ( abs(fplus) <= me%ftol ) then
            ! Absolute convergence in f
            xzero = xplus
            fzero = fplus
            exit
        end if

        fminus = me%f(xminus)
        if ( abs(fminus) <= me%ftol ) then
            ! Absolute convergence in f
            xzero = xminus
            fzero = fminus
            exit
        end if

        !  Take whichever of the two quadratic roots is closest to a root of the function.
        if ( abs(fminus) < abs(fplus) ) then
            xzero = xminus
            fzero = fminus
        else
            xzero = xplus
            fzero = fplus
        end if

        fxold = fxmid
        fxmid = fxnew
        fxnew = fzero

        x_inc = xzero - xmid
        if ( abs ( x_inc ) <= me%atol ) exit ! Absolute convergence in X

        x_ave = ( abs ( xzero ) + abs ( xmid ) + abs ( xold ) ) / 3.0_wp
        if ( abs ( x_inc ) <= me%rtol * x_ave ) exit ! Relative convergence in X

        if ( i == me%maxiter ) iflag = -2 ! max iterations exceeded

    end do

    contains

        pure elemental subroutine swap(a,b)

        !! Swap two real(wp) values.
        implicit none

        real(wp),intent(inout) :: a
        real(wp),intent(inout) :: b

        real(wp) :: tmp

        tmp = a
        a   = b
        b   = tmp

        end subroutine swap

    end subroutine muller
!*****************************************************************************************

!*****************************************************************************************
    end module stdlib_root_module
!*****************************************************************************************
