module constants_mod

    implicit none

   !MPI variables
#ifdef USE_MPI
    integer, parameter :: MPI_GOOD_TAG = 1
    logical, parameter :: reorder=.false.   !<--- probably shouldn't reoder processors
#endif

    integer         , parameter :: dp = kind(1.0d0)

    real(dp)    , parameter :: m_pml = 3.0d0
    real(dp)    , parameter :: ma_pml = 1.0d0
    real(dp)    , parameter :: kappa_max_pml = 5.0d0
    real(dp)    , parameter :: c0      = 1.0d0 !137.035999679
    real(dp)    , parameter :: eps0_SI = 8.8541878188d-12
    real(dp)    , parameter :: c0_SI   = 299792458.0d0
    real(dp)    , parameter :: pi      = 3.141592653589793238462643383279502884197d0
    real(dp)    , parameter :: eps0    = 1.0d0 !0.079577471545947672804111050481878919527
    real(dp)    , parameter :: mu0     = 1.0_dp
    real(dp)    , parameter :: R_0     = 0.0_dp
    complex(dp) , parameter :: Z_I     = (0.0_dp, 1.0_dp)
    complex(dp) , parameter :: Z_ONE   = (1.0_dp, 0.0_dp)
    complex(dp) , parameter :: Z_0     = (0.0_dp, 0.0_dp)

    integer, public, parameter :: &
    CLOSE_BOUNDARIES     = 1,                &
    PERIODIC_BOUNDARIES  = 2,                &
    PML_BOUNDARIES       = 3

    integer, public, parameter :: &
    Jx_SOURCE       = 1,                &
    Jy_SOURCE       = 2,                &
    Jz_SOURCE       = 3,                &
    Re_Ex_TARGET    = 4,                &
    Im_Ex_TARGET    = 5,                &
    Re_Ey_TARGET    = 6,                &
    Im_Ey_TARGET    = 7,                &
    Re_Ez_TARGET    = 8,                &
    Im_Ez_TARGET    = 9,                &
    Abs_Ex_TARGET   = 10,               &
    Abs_Ey_TARGET   = 11,               &
    Abs_Ez_TARGET   = 12

end module constants_mod