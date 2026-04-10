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
    real(dp)    , parameter :: mu0     = 1.0d0
    real(dp)    , parameter :: M_0     = 0.0d0
    complex(dp) , parameter :: Z_I     = (0.0d0, 1.0d0)
    complex(dp) , parameter :: Z_ONE   = (1.0d0, 0.0d0)
    complex(dp) , parameter :: Z_0     = (0.0d0, 0.0d0)

    integer, public, parameter :: &
    CLOSE_BOUNDARIES     = 1,                &
    PERIODIC_BOUNDARIES  = 2,                &
    PML_BOUNDARIES       = 3

end module constants_mod