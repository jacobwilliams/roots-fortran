program modab_edge_tests
!! Edge-case tests for the shared modAB safeguards, driven through [[root_scalar]].
!!
!! The 100-problem benchmark suite in root_tests never produces a non-finite or
!! overflowing residual, so none of it reaches the overflow/NaN branches of
!! same_sign, safe_midpoint and safe_secant, or the overflow cases of the
!! switching test. The brackets below do.
!!
!! Those helpers are private to root_module, so this program reaches them only
!! through the solver. The first case is decisive on its own: before
!! safe_midpoint it returned +Inf. The rest pin down the behaviour of brackets
!! and residuals at the extremes of binary64, where the old inline arithmetic
!! overflowed even though the mathematical answer is an ordinary finite number.

    use root_module, only: wp => root_module_rk, root_scalar
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan, &
                               ieee_positive_inf, ieee_negative_inf

    implicit none

    integer :: passed, total

    passed = 0
    total = 0

    ! x1 + x2 overflows although both endpoints are finite. The naive midpoint
    ! (x1+x2)*0.5 becomes +Inf, the tolerance test then passes trivially and the
    ! solver returns +Inf. safe_midpoint computes 0.5*x1 + 0.5*x2 instead.
    call ck_root('x1 + x2 overflows',           f_mid,     1.0e308_wp,  1.7e308_wp,  1.2e308_wp)
    call ck_root('x1 + x2 overflows, negative', f_mid_neg, -1.7e308_wp, -1.0e308_wp, -1.2e308_wp)

    ! |f1| + |f2| overflows, so the old symmetry factor divided by an infinite
    ! f2 - f1 and always saw a perfectly symmetric bracket. The switching test
    ! is now skipped while f2 - f1 overflows, so the solver keeps bisecting.
    call ck_root('|f1| + |f2| overflows',       f_sym,     0.0_wp,      1.0_wp,      1.7_wp/2.2_wp)
    call ck_root('|f1| + |f2| overflows, wide', f_sym2,    -0.6_wp,     1.0_wp,      0.6_wp/1.62_wp)

    ! Residuals large enough that the textbook secant numerator x1*y2 - y1*x2
    ! overflows while the intersection itself is an ordinary interior point.
    call ck_root('step-like huge residuals',    f_step,    0.0_wp,      1.0_wp,      0.25_wp)

    ! A NaN residual has no usable sign, so the bracket cannot be updated and
    ! the solver must not report success.
    call ck_rejects('NaN residual mid-solve',   f_nan,     0.0_wp,      1.0_wp)

    ! Overflow and infinite residuals in the switching test.
    ! |ym| + |y3| overflows at the first midpoint while |ym - y3| is finite.
    call ck_root('hump',                        f_hump,      -1.0_wp, 1.0_wp, -5.0_wp/6.0_wp)
    ! f2 - f1 overflows: switching is disabled until the residuals shrink.
    call ck_root('f2 - f1 overflows',           f_tanh,      -1.0_wp, 1.0_wp, 0.3_wp)
    ! Infinite residuals at one or both ends of the bracket.
    call ck_root('infinite left residual',      f_inf_left,  -1.0_wp, 1.0_wp, 0.1_wp)
    call ck_root('infinite residuals',          f_inf_both,  -1.0_wp, 1.0_wp, 0.1_wp)
    ! Subnormal residuals: A&B corrections underflow towards zero.
    call ck_root('subnormal residuals',         f_subnormal, -1.0_wp, 1.0_wp, 0.2_wp**(1.0_wp/3.0_wp))
    ! f101, L. Tomov counterexample: A&B keeps halving the residual without shrinking
    ! the bracket; the residual test must not keep it on A&B, but fall back to bisection.
    call ck_root('residual test',               f_tomov,     -acos(-1.0_wp), 1.0_wp, 0.0_wp)  ! [-pi, 1]

    write(*,'(A,I0,A,I0,A)') 'Fortran edge-cases: ', passed, '/', total, &
        merge(' PASS', ' FAIL', passed == total)
    if (passed /= total) error stop 1

contains

    !! f(x) = x*1e-308 - 1.2; the root is 1.2e308 and x1 + x2 overflows.
    function f_mid(x) result(f)
    real(wp),intent(in) :: x
    real(wp) :: f
    f = x*1.0e-308_wp - 1.2_wp
    end function f_mid

    !! Mirror of f_mid on the negative axis; the root is -1.2e308.
    function f_mid_neg(x) result(f)
    real(wp),intent(in) :: x
    real(wp) :: f
    f = x*1.0e-308_wp + 1.2_wp
    end function f_mid_neg

    !! Endpoint residuals -1.7e308 and 0.5e308: finite apart, overflowing summed.
    !! Written so that no intermediate product overflows. The root is 1.7/2.2.
    function f_sym(x) result(f)
    real(wp),intent(in) :: x
    real(wp) :: f
    f = 1.7e308_wp*(x - 1.0_wp) + 0.5e308_wp*x
    end function f_sym

    !! Residuals -0.6e308 and 1.02e308 over an asymmetric bracket.
    function f_sym2(x) result(f)
    real(wp),intent(in) :: x
    real(wp) :: f
    f = 1.62e308_wp*x - 0.6e308_wp
    end function f_sym2

    !! A near-discontinuity with residuals at the top of the exponent range.
    function f_step(x) result(f)
    real(wp),intent(in) :: x
    real(wp) :: f
    if (x == 0.25_wp) then
        f = 0.0_wp
    else
        f = sign(1.0e300_wp, x - 0.25_wp)
    end if
    end function f_step

    !! Finite endpoints of opposite sign, but NaN across the middle of the
    !! bracket, so the first bisection step produces a NaN residual.
    function f_nan(x) result(f)
    real(wp),intent(in) :: x
    real(wp) :: f
    if (x > 0.4_wp .and. x < 0.6_wp) then
        f = ieee_value(f, ieee_quiet_nan)
    else
        f = x - 0.5_wp
    end if
    end function f_nan

    !! Linear pieces up to huge(): |ym| + |y3| overflows at the first midpoint.
    function f_hump(x) result(f)
    real(wp),intent(in) :: x
    real(wp) :: f
    if (x <= 0.0_wp) then
        f = huge(1.0_wp)*(-0.2_wp + 1.2_wp*(x + 1.0_wp))
    else
        f = huge(1.0_wp)*(1.0_wp - 0.2_wp*x)
    end if
    end function f_hump

    !! Residuals near +-1.7e308, so f2 - f1 overflows; the root is 0.3.
    function f_tanh(x) result(f)
    real(wp),intent(in) :: x
    real(wp) :: f
    f = 1.7e308_wp*tanh(10.0_wp*(x - 0.3_wp))
    end function f_tanh

    !! -Inf left of -0.5; the root is 0.1.
    function f_inf_left(x) result(f)
    real(wp),intent(in) :: x
    real(wp) :: f
    if (x < -0.5_wp) then
        f = ieee_value(f, ieee_negative_inf)
    else
        f = x - 0.1_wp
    end if
    end function f_inf_left

    !! -Inf left of -0.5 and +Inf right of 0.9; the root is 0.1.
    function f_inf_both(x) result(f)
    real(wp),intent(in) :: x
    real(wp) :: f
    if (x < -0.5_wp) then
        f = ieee_value(f, ieee_negative_inf)
    else if (x > 0.9_wp) then
        f = ieee_value(f, ieee_positive_inf)
    else
        f = x - 0.1_wp
    end if
    end function f_inf_both

    !! Residuals of order 1e-300; the root is 0.2**(1/3).
    function f_subnormal(x) result(f)
    real(wp),intent(in) :: x
    real(wp) :: f
    f = 1.0e-300_wp*(x*x*x - 0.2_wp)
    end function f_subnormal

    !! f101: -x**2 for x < 0, 7x**2 otherwise; the root is 0.
    function f_tomov(x) result(f)
    real(wp),intent(in) :: x
    real(wp) :: f
    if (x < 0.0_wp) then
        f = -x*x
    else
        f = 7.0_wp*x*x
    end if
    end function f_tomov

    !! Asserts the solver reports success and lands on the expected root.
    subroutine ck_root(name, fun, ax, bx, want)
    character(len=*),intent(in) :: name
    procedure(fun_if) :: fun
    real(wp),intent(in) :: ax, bx, want

    real(wp) :: root, froot, scale
    integer :: iflag

    total = total + 1
    call root_scalar('modab', fun, ax, bx, root, froot, iflag, &
                     atol=1.0e-14_wp, rtol=1.0e-14_wp)
    scale = max(abs(want), 1.0_wp)
    if (iflag == 0 .and. .not. ieee_is_nan(root) .and. &
        abs(root - want) <= 1.0e-12_wp*scale) then
        passed = passed + 1
    else
        write(*,'(A,A,A,ES24.17,A,ES24.17,A,I0)') 'FAIL ', name, &
            ': got ', root, ' want ', want, ' iflag ', iflag
    end if
    end subroutine ck_root

    !! Asserts the solver refuses the bracket rather than returning a bogus root.
    subroutine ck_rejects(name, fun, ax, bx)
    character(len=*),intent(in) :: name
    procedure(fun_if) :: fun
    real(wp),intent(in) :: ax, bx

    real(wp) :: root, froot
    integer :: iflag

    total = total + 1
    call root_scalar('modab', fun, ax, bx, root, froot, iflag, &
                     atol=1.0e-14_wp, rtol=1.0e-14_wp)
    if (iflag /= 0 .or. ieee_is_nan(root)) then
        passed = passed + 1
    else
        write(*,'(A,A,A,ES24.17,A,I0)') 'FAIL ', name, &
            ': reported success at ', root, ' iflag ', iflag
    end if
    end subroutine ck_rejects

    !! Shape of the user functions above, for the dummy-procedure declarations.
    function fun_if(x) result(f)
    real(wp),intent(in) :: x
    real(wp) :: f
    f = x
    end function fun_if

end program modab_edge_tests
