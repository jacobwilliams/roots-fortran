![roots-fortran](media/logo.png)
============

**roots-fortran**: root solvers for modern Fortran

[![GitHub release](https://img.shields.io/github/release/jacobwilliams/roots-fortran.svg?style=plastic)](https://github.com/jacobwilliams/roots-fortran/releases/latest)
![Build Status](https://github.com/jacobwilliams/roots-fortran/actions/workflows/CI.yml/badge.svg)
[![codecov](https://codecov.io/gh/jacobwilliams/roots-fortran/branch/master/graph/badge.svg?token=43HK33CSMY)](https://codecov.io/gh/jacobwilliams/roots-fortran)

## Description

A modern Fortran library for finding the roots of continuous scalar functions of a single real variable.

## Compiling

A `fmp.toml` file is provided for compiling roots-fortran with the [Fortran Package Manager](https://github.com/fortran-lang/fpm). For example, to build:

```
fpm build --profile release
```

By default, the library is built with double precision (`real64`) real values. Explicitly specifying the real kind can be done using the following processor flags:

Preprocessor flag | Kind  | Number of bytes
----------------- | ----- | ---------------
`REAL32`  | `real(kind=real32)`  | 4
`REAL64`  | `real(kind=real64)`  | 8
`REAL128` | `real(kind=real128)` | 16

For example, to build a single precision version of the library, use:

```
fpm build --profile release --flag "-DREAL32"
```

To run the unit tests:

```
fpm test
```

To use `roots-fortran` within your fpm project, add the following to your `fpm.toml` file:
```toml
[dependencies]
roots-fortran = { git="https://github.com/jacobwilliams/roots-fortran.git" }
```

or, to use a specific version:
```toml
[dependencies]
roots-fortran = { git="https://github.com/jacobwilliams/roots-fortran.git", tag = "1.0.0"  }
```

To generate the documentation using [ford](https://github.com/Fortran-FOSS-Programmers/ford), run: ```ford roots-fortran.md```

## Usage

### Methods

The module contains the following methods (in alphabetical order):

 * `anderson_bjorck`
 * `anderson_bjorck_king` (a variant of `anderson_bjorck`)
 * `barycentric`
 * `bdqrf`
 * `bisection`
 * `blendtf`
 * `brent`
 * `brenth` (a variant of `brent`)
 * `brentq` (a variant of `brent`)
 * `chandrupatla`
 * `illinois`
 * `muller`
 * `pegasus`
 * `regula_falsi`
 * `ridders`
 * `toms748`
 * `zhang`

In general, all the methods are guaranteed to converge. Some will be more efficient (in terms of number of function evaluations) than others for various problems. The methods can be broadly classified into three groups:

 * Simple classical methods (`bisection`, `regula_falsi`, `illinois`, `ridders`).
 * Newfangled methods (`zhang`, `barycentric`, `blendtf`, `bdqrf`, `anderson_bjorck_king`). These rarely or ever seem to be better than the best methods.
 * Best methods (`anderson_bjorck`, `muller`, `pegasus`, `toms748`, `brent`, `brentq`, `brenth`, `chandrupatla`). Generally, one of these will be the most efficient method.

### Functional Interface Example

```fortran
program main

  use root_module, wp => root_module_rk

  implicit none

  real(wp) :: x, f
  integer :: iflag

  call root_scalar('bisection',func,-9.0_wp,31.0_wp,x,f,iflag)

  write(*,*) 'f(',x,') = ', f
  write(*,*) 'iflag    = ', iflag

contains

  function func(x) result(f)

  implicit none

  real(wp),intent(in) :: x
  real(wp) :: f

  f = -200.0_wp * x * exp(-3.0_wp*x)

  end function func

end program main
```

### Object Oriented Interface Example

```fortran
program main

  use root_module, wp => root_module_rk

  implicit none

  type,extends(bisection_solver) :: my_solver
  end type my_solver

  real(wp) :: x, f
  integer :: iflag
  type(my_solver) :: solver

  call solver%initialize(func)
  call solver%solve(-9.0_wp,31.0_wp,x,f,iflag)

  write(*,*) 'f(',x,') = ', f
  write(*,*) 'iflag    = ', iflag

contains

  function func(me,x)

    class(root_solver),intent(inout) :: me
    real(wp),intent(in) :: x
    real(wp) :: f

    f = -200.0_wp * x * exp(-3.0_wp*x)

  end function func

end program main
```

### Result

```
 f( -2.273736754432321E-013 ) =   4.547473508867743E-011
 iflag    =            0
```

## Notes

* Originally this was an idea for the Fortran [stdlib](https://github.com/fortran-lang/stdlib). See: [#87: Optimization, Root finding, and Equation Solvers](https://github.com/fortran-lang/stdlib/issues/87). Eventually, it may be merged into the other one in one form or another.

## Documentation

The latest API documentation for the `master` branch can be found [here](https://jacobwilliams.github.io/roots-fortran/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

## License

The roots-fortran source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/roots-fortran/blob/master/LICENSE.md) (BSD-style).

## Similar libraries in other programming languages

* Julia: [Roots.jl](https://github.com/JuliaMath/Roots.jl)
* Python: [scipy.optimize.root_scalar](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root_scalar.html)
* C: [GSL](https://www.gnu.org/software/gsl/doc/html/roots.html)

## References
  * D. E. Muller, "[A Method for Solving Algebraic Equations Using an Automatic Computer](https://www.ams.org/journals/mcom/1956-10-056/S0025-5718-1956-0083822-0/S0025-5718-1956-0083822-0.pdf)", Mathematical Tables and Other Aids to Computation, 10 (1956), 208-215.
  * M. Dowell, P. Jarratt, "[A modified regula falsi method for computing the root of an equation](https://personal.math.ubc.ca/~loew/mech2/Dowell+Jarratt.pdf)", BIT 11 (1971), 168-174.
  * R. P. Brent, "[An algorithm with guaranteed convergence for finding a zero of a function](http://maths-people.anu.edu.au/~brent/pd/rpb005.pdf)", The Computer Journal, Vol 14, No. 4., 1971.
  * R. P. Brent, "[Algorithms for minimization without derivatives](http://maths-people.anu.edu.au/~brent/pub/pub011.html)", Prentice-Hall, Inc., 1973.
  * Ridders, C., "[A new algorithm for computing a single root of a real continuous function](https://cs.fit.edu/~dmitra/SciComp/Resources/RidderMethod.pdf)", IEEE Trans. on Circuits and Systems, Vol 26, Issue 11, Nov 1979.
  * G. E. Alefeld, F. A. Potra and Yixun Shi, "[Algorithm 748: Enclosing Zeros of Continuous Functions](https://dl.acm.org/doi/abs/10.1145/210089.210111)", ACM Transactions on Mathematical Software, Vol. 21. No. 3. September 1995. Pages 327-344.
  * G.E. Mullges & F. Uhlig, "Numerical Algorithms with Fortran", Springer, 1996. Section 2.8.1, p 32-34.
  * T.R. Chandrupatla, "[A new hybrid quadratic/bisection algorithm for finding the zero of a nonlinear function without derivatives](https://dl.acm.org/doi/10.1016/S0965-9978%2896%2900051-8)", Advances in Engineering Software, Vol 28, 1997, pp. 145-149.
  * R. G. Gottlieb, B. F. Thompson, "[Bisected Direct Quadratic Regula Falsi](https://www.researchgate.net/publication/228712261_Bisected_Direct_Quadratic_Regula_Falsi)", Applied Mathematical Sciences, Vol. 4, 2010, no. 15, 709-718.
  * A. Zhang, "[An Improvement to the Brent's Method](https://www.cscjournals.org/download/issuearchive/IJEA/Volume2/IJEA_V2_I1.pdf)", International Journal of Experimental Algorithms (IJEA), Volume (2) : Issue (1) : 2011.
  * E Badr, S Almotairi, A El Ghamry, "[A Comparative Study among New Hybrid Root Finding Algorithms and Traditional Methods](https://www.mdpi.com/2227-7390/9/11/1306)", Mathematics 2021, 9, 1306.

## See also

 * [Code coverage statistics](https://app.codecov.io/gh/jacobwilliams/roots-fortran) [codecov.io]
