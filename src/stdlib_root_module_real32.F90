!*****************************************************************************************
!>
!  32 bit root solver for [[stdlib_root_module]].

    module stdlib_root_module_real32

    use iso_fortran_env, only: wp => real32

#define INCLUDE_STDLIB_ROOT_CORE
#include "stdlib_root_core.F90"

    end module stdlib_root_module_real32
!*****************************************************************************************
