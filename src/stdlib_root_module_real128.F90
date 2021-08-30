!*****************************************************************************************
!>
!  128 bit root solver for [[stdlib_root_module]].

    module stdlib_root_module_real128

#ifdef HAS_REAL128
#define REAL_KIND real128
#ifdef __GFORTRAN__
#define RKIND(__PROCEDURE__) __PROCEDURE__/**/_real128
#else
#define RKIND(__PROCEDURE__) __PROCEDURE__##_real128
#endif
#include "stdlib_root_core.inc"
#undef REAL_KIND
#endif

    end module stdlib_root_module_real128
!*****************************************************************************************
