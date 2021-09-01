![roots-fortran](media/logo.png)
============

**roots-fortran**: root solvers for modern Fortran

[![GitHub release](https://img.shields.io/github/release/jacobwilliams/roots-fortran.svg?style=plastic)](https://github.com/jacobwilliams/roots-fortran/releases/latest)
![Build Status](https://github.com/jacobwilliams/roots-fortran/actions/workflows/CI.yml/badge.svg)

## Description

A library for finding the roots of continuous scalar functions of a single real variable.

## Compiling

### FPM

A `fmp.toml` file is provided for compiling roots-fortran with the [Fortran Package Manager](https://github.com/fortran-lang/fpm). For example, to build:

```
  fpm build --profile release
```

And to run the unit tests:

```
  fpm test
```

### FoBiS

A [FoBiS](https://github.com/szaghi/FoBiS) configuration file (`roots-fortran.fobis`) is also provided that can also build the library and tests. Use the `mode` flag to indicate what to build. For example:

  * To build all the examples using gfortran: `FoBiS.py build -f roots-fortran.fobis -mode tests-gnu`
  * To build all the examples using ifort: `FoBiS.py build -f roots-fortran.fobis -mode tests-intel`
  * To build a static library using gfortran: `FoBiS.py build -f roots-fortran.fobis -mode static-gnu`
  * To build a static library using ifort: `FoBiS.py build -f roots-fortran.fobis -mode static-intel`

  The full set of modes are: `static-gnu`, `static-gnu-debug`, `static-intel`, `static-intel-debug`, `shared-gnu`, `shared-gnu-debug`, `shared-intel`, `shared-intel-debug`, `tests-gnu`, `tests-gnu-debug`, `tests-intel`, `tests-intel-debug`

  To generate the documentation using [ford](https://github.com/Fortran-FOSS-Programmers/ford), run: ```FoBis.py rule --execute makedoc -f roots-fortran.fobis```


## Usage

### Methods

The module contains the following methods:

 * anderson_bjorck
 * anderson_bjorck_king
 * bdqrf
 * bisection
 * blendtf
 * brent
 * brenth
 * brentq
 * chandrupatla
 * illinois
 * muller
 * pegasus
 * regula_falsi
 * ridders
 * toms748
 * zhang

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

The latest API documentation can be found [here](http://jacobwilliams.github.io/roots-fortran/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

## License

The roots-fortran source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/roots-fortran/blob/master/LICENSE.md) (BSD-style).

## Similar libraries in other programming languages

* Julia: [Roots.jl](https://github.com/JuliaMath/Roots.jl)
* Python: [scipy.optimize.root_scalar](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root_scalar.html)
* C: [GSL](https://www.gnu.org/software/gsl/doc/html/roots.html)
