module constants_mod
    
    implicit none
   
    !MPI variables
#ifdef USE_MPI
    integer, parameter :: MPI_GOOD_TAG = 1
    logical, parameter :: reorder=.false.   !<--- probably shouldn't reoder processors
#endif

    !costants
    
    integer         , parameter :: dp = kind(1.0d0)
    double precision, parameter :: pi0   = 3.141592653589793d0
    double precision, parameter :: c0    = 137.0                 !299792458.0d0
    double precision, parameter :: sqrt2 = 1.414213562373095d0
    double precision, parameter :: sqrt3 = 1.732050807568877d0
    double precision, parameter :: sqrt6 = 2.449489742783178d0
    double precision, parameter :: eps0  = 1.0/(4.0*pi0)         !1.0d0/(c*c*mu0)
    double precision, parameter :: mu0   = 1.0/(eps0*c0**2)      !4.0d-7*pi
    double precision, parameter :: hbar  = 1.0_dp                  !1.054571628d-34
    double precision, parameter :: M_ZERO= 0.0_dp
    double precision, parameter :: M_ONE = 1.0_dp
    double complex  , parameter :: Z_I   = (0.0_dp, 1.0_dp)
    double complex  , parameter :: Z_0   = (0.0_dp, 0.0_dp)
    double complex  , parameter :: Z_ONE = (1.0_dp, 0.0_dp)
    double precision, parameter :: alphaCPML=0.05,kappaCPML=5.0
    double precision, parameter, dimension (4) :: aBH=(/0.353222222d0,-0.488d0,0.145d0,-0.010222222d0/) !Blackman-Harris window 
    integer         , parameter :: m=3,ma=1

    !unit convertions

    double precision, parameter :: sec_to_au       = 1.0/2.4188843265864D-17
    double precision, parameter :: C_to_au         = 1.0/1.602176634D-19
    double precision, parameter :: kg_to_au        = 1.0/9.1093837139D-31
    double precision, parameter :: m_to_au         = 1.0/5.29177210544D-11
    double precision, parameter :: nm_to_au        = 18.897261259077823
    double precision, parameter :: au_to_nm        = 1.0/nm_to_au
    double precision, parameter :: ev_to_radsec    = 2.0*pi0*2.418d14
    double precision, parameter :: ev_to_au        = 1.0/27.2114
    double precision, parameter :: Hz_to_ev        = 4.1356691d-15
    double precision, parameter :: Debye_to_Cm     = 3.33564d-30
    double precision, parameter :: Debye_to_au     = 3.0/7.63
    double precision, parameter :: fs_to_au        = 1.0d0/2.4188843265864D-2  !0.0241888432650516d0
    double precision, parameter :: au_to_fs        = 2.4188843265864D-2

    double precision, parameter :: dipole_au_to_SI = 8.47835281d-30

    !options

    integer, public, parameter :: &
    CLOSE_BOUNDARIES     = 1,                &
    PERIODIC_BOUNDARIES  = 2,                &
    CPML_BOUNDARIES      = 3

    integer, public, parameter :: &
    TMZ_2D_MODE   = 1,                &
    TEZ_2D_MODE   = 2,                &
    FULL_2D_MODE  = 3

    integer, public, parameter :: &
    NONE_MEDIUM        = 0,                &
    DIELECTRIC_MEDIUM  = 1,                &
    DRUDE_MEDIUM       = 2,                &
    DL_MEDIUM          = 3

    integer, public, parameter :: &
    DFTB_ELEC_DYN    = 1,                &
    DFTB_EHREN_DYN   = 2,                &
    DFTB_BO_DYN      = 3

    integer, public, parameter :: &
    Q_MATERIAL    = 1,            &  
    Q_SINGLE      = 2
    
    integer, public, parameter :: &
    Q_SYS_DFTB    = 1

    integer, public, parameter :: &
    POINT_DETECTOR  = 1,                &
    LINE_X_DETECTOR   = 2,                &
    LINE_Y_DETECTOR  = 3,                &
    LINE_Z_DETECTOR  = 4,                &
    PLANE_XY_DETECTOR  = 5,                &
    PLANE_YZ_DETECTOR  = 6,                &
    PLANE_ZX_DETECTOR  = 7,                &
    VOLUME_DETECTOR = 8

    integer, public, parameter :: &
    Ex_FIELD = 1,                &
    Ey_FIELD = 2,                &
    Ez_FIELD = 3,                &
    Hx_FIELD = 4,                &
    Hy_FIELD = 5,                &
    Hz_FIELD = 6

end module constants_mod