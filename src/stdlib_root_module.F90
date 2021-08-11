!*****************************************************************************************
!>
!  Root solver methods for:
!  * Bracked interval
!  * Without derivatives

    module stdlib_root_module

    use stdlib_root_module_real32,  only: root_solver_real32 => root_solver, &
                                          root_scalar_real32 => root_scalar
    use stdlib_root_module_real64,  only: root_solver, &
                                          root_scalar_real64 => root_scalar, &
                                          root_solver_real64 => root_solver
#ifdef HAS_REAL128
    use stdlib_root_module_real128, only: root_solver_real128 => root_solver, &
                                          root_scalar_real128 => root_scalar
#endif

    implicit none

    private

    interface root_scalar
        module procedure :: root_scalar_real32
        module procedure :: root_scalar_real64
#ifdef HAS_REAL128
        module procedure :: root_scalar_real128
#endif
    end interface root_scalar
    public :: root_scalar

    public :: root_solver_real32
    public :: root_solver_real64, root_solver
#ifdef HAS_REAL128
    public :: root_solver_real128
#endif

    end module stdlib_root_module
!*****************************************************************************************