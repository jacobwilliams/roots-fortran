!*****************************************************************************************
!>
!  Root solver methods for:
!  * Bracked interval
!  * Without derivatives

    module stdlib_root_module

    use stdlib_root_module_real32,  only: root_solver_real32                 => root_solver, &
                                          root_scalar_real32                 => root_scalar, &
                                          brent_solver_real32                => brent_solver, &
                                          bisection_solver_real32            => bisection_solver, &
                                          regula_falsi_solver_real32         => regula_falsi_solver, &
                                          illinois_solver_real32             => illinois_solver, &
                                          anderson_bjorck_solver_real32      => anderson_bjorck_solver, &
                                          ridders_solver_real32              => ridders_solver, &
                                          pegasus_solver_real32              => pegasus_solver, &
                                          bdqrf_solver_real32                => bdqrf_solver, &
                                          muller_solver_real32               => muller_solver, &
                                          brenth_solver_real32               => brenth_solver, &
                                          brentq_solver_real32               => brentq_solver, &
                                          chandrupatla_solver_real32         => chandrupatla_solver, &
                                          toms748_solver_real32              => toms748_solver, &
                                          zhang_solver_real32                => zhang_solver, &
                                          anderson_bjorck_king_solver_real32 => anderson_bjorck_king_solver, &
                                          blendtf_solver_real32              => blendtf_solver

    use stdlib_root_module_real64,  only: root_solver_real64                 => root_solver, &
                                          root_scalar_real64                 => root_scalar, &
                                          brent_solver_real64                => brent_solver, &
                                          bisection_solver_real64            => bisection_solver, &
                                          regula_falsi_solver_real64         => regula_falsi_solver, &
                                          illinois_solver_real64             => illinois_solver, &
                                          anderson_bjorck_solver_real64      => anderson_bjorck_solver, &
                                          ridders_solver_real64              => ridders_solver, &
                                          pegasus_solver_real64              => pegasus_solver, &
                                          bdqrf_solver_real64                => bdqrf_solver, &
                                          muller_solver_real64               => muller_solver, &
                                          brenth_solver_real64               => brenth_solver, &
                                          brentq_solver_real64               => brentq_solver, &
                                          chandrupatla_solver_real64         => chandrupatla_solver, &
                                          toms748_solver_real64              => toms748_solver, &
                                          zhang_solver_real64                => zhang_solver, &
                                          anderson_bjorck_king_solver_real64 => anderson_bjorck_king_solver, &
                                          blendtf_solver_real64              => blendtf_solver

#ifdef HAS_REAL128
    use stdlib_root_module_real128, only: root_solver_real128                 => root_solver, &
                                          root_scalar_real128                 => root_scalar, &
                                          brent_solver_real128                => brent_solver, &
                                          bisection_solver_real128            => bisection_solver, &
                                          regula_falsi_solver_real128         => regula_falsi_solver, &
                                          illinois_solver_real128             => illinois_solver, &
                                          anderson_bjorck_solver_real128      => anderson_bjorck_solver, &
                                          ridders_solver_real128              => ridders_solver, &
                                          pegasus_solver_real128              => pegasus_solver, &
                                          bdqrf_solver_real128                => bdqrf_solver, &
                                          muller_solver_real128               => muller_solver, &
                                          brenth_solver_real128               => brenth_solver, &
                                          brentq_solver_real128               => brentq_solver, &
                                          chandrupatla_solver_real128         => chandrupatla_solver, &
                                          toms748_solver_real128              => toms748_solver, &
                                          zhang_solver_real128                => zhang_solver, &
                                          anderson_bjorck_king_solver_real128 => anderson_bjorck_king_solver, &
                                          blendtf_solver_real128              => blendtf_solver
#endif

    implicit none

    public  ! everything is public

    interface root_scalar
        !! Function interface to the root solvers.
        module procedure :: root_scalar_real32
        module procedure :: root_scalar_real64
#ifdef HAS_REAL128
        module procedure :: root_scalar_real128
#endif
    end interface root_scalar
    public :: root_scalar

    public :: root_solver_real32
    public :: root_solver_real64
#ifdef HAS_REAL128
    public :: root_solver_real128
#endif

    end module stdlib_root_module
!*****************************************************************************************