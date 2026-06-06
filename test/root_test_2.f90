program root_test_2

    use root_module, wp => root_module_rk

    implicit none
    integer,parameter :: n = 100000
    real(wp) :: xzero
    integer :: k, iflag
    real(wp),dimension(n) :: p, out
    type(brenth_solver) :: solver

    write(*,*) ''
    write(*,*) '-------------------------------------------------'
    write(*,*) 'root_test_2'
    write(*,*) '-------------------------------------------------'
    write(*,*) ''

    call random_number(p); p = 1.5_wp * p
    call solver%initialize(f, rtol=epsilon(1.0_wp), atol=epsilon(1.0_wp))
    do k = 1, n
       call solver%solve(0.0_wp,2.0_wp,xzero,out(k),iflag)
    end do
    print*, 'norm   =', norm2(out)
    print*, 'maxval =', maxval(abs(out))
    print*, 'minval =', minval(abs(out))

contains

    function f(me,x)
        class(root_solver),intent(inout) :: me
        real(wp),intent(in) :: x
        real(wp) :: f
        f = x * sin(x) - p(k)
    end function f

end program root_test_2