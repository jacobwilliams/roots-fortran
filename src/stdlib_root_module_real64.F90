!*****************************************************************************************
!>
!  64 bit root solver for [[stdlib_root_module]].

    module stdlib_root_module_real64

#define REAL_KIND real64
#ifdef __GFORTRAN__
#define RKIND(__PROCEDURE__) __PROCEDURE__/**/_real64
#else
#define RKIND(__PROCEDURE__) __PROCEDURE__##_real64
#endif
#include "stdlib_root_core.inc"
#undef REAL_KIND

    end module stdlib_root_module_real64
!*****************************************************************************************
