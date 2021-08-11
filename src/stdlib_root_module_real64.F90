!*****************************************************************************************
!>
!  64 bit root solver [[stdlib_root_module]].

    module stdlib_root_module_real64

    use iso_fortran_env, only: wp => real64

#define INCLUDE_STDLIB_ROOT_CORE
#include "stdlib_root_core.F90"

    end module stdlib_root_module_real64
!*****************************************************************************************
