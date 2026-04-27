module operator_mod

    use constants_mod
    use rs_vec_base_mod
    use rs_vec_dimensions_mod
    use medium_mod
    use parallel_subs_mod

    implicit none

    type :: TOperator
 
        logical  :: cpml_pos(6) = .false.  ! [left, right, front, back, top, bottom]
        integer  :: dimensions
        integer  :: boundaries(3)
        integer  :: n_pml
        integer  :: n_der
        integer  :: nx
        integer  :: ny
        integer  :: nz
        real(dp) :: dr
        real(dp) :: w0 !frequency

        complex(dp), allocatable :: s_inv_x(:)
        complex(dp), allocatable :: s_inv_y(:)
        complex(dp), allocatable :: s_inv_z(:)
        real(dp)   , allocatable :: coef_der(:)

    contains
        procedure :: init_operator
        procedure :: kill_operator
        procedure :: apply_operator
        
    end type TOperator

contains

!###################################################################################################

subroutine init_operator(this, dr, freq, dimensions, grid_Ndims, n_pml, n_der, &
                         boundaries,mpi_dims, mpi_coords)
    
    class(TOperator), intent(inout) :: this
    real(dp)    , intent(in)        :: dr
    real(dp)    , intent(in)        :: freq
    integer     , intent(in)        :: dimensions
    integer     , intent(in)        :: n_pml
    integer     , intent(in)        :: grid_Ndims(3)
    integer     , intent(in)        :: boundaries(3)
    integer     , intent(in)        :: mpi_dims(3)
    integer     , intent(in)        :: mpi_coords(3)
    integer     , intent(in)        :: n_der
    complex(dp) :: s_i
    real(dp)    :: sig_max
    real(dp)    :: kappa
    real(dp)    :: sig
    integer     :: i, ii

    this%dimensions     = dimensions
    this%boundaries     = boundaries
    this%n_pml          = n_pml
    this%dr             = dr
    this%w0             = freq
    this%n_der          = n_der

    if (dimensions == 1) then
        this%nx         = grid_Ndims(1)
        this%ny         = 0
        this%nz         = 0
    else if (dimensions == 2) then
        this%nx         = grid_Ndims(1)
        this%ny         = grid_Ndims(2)
        this%nz         = 0
    else if (dimensions == 3) then
        this%nx         = grid_Ndims(1)
        this%ny         = grid_Ndims(2)
        this%nz         = grid_Ndims(3)
    end if

    if (.not. allocated(this%s_inv_x)) allocate(this%s_inv_x(this%nx))
    if (.not. allocated(this%s_inv_y)) allocate(this%s_inv_y(this%ny))
    if (.not. allocated(this%s_inv_z)) allocate(this%s_inv_z(this%nz))
    
    this%s_inv_x = Z_ONE
    this%s_inv_y = Z_ONE
    this%s_inv_z = Z_ONE
    
    if (.not. allocated(this%coef_der)) allocate(this%coef_der(-this%n_der:this%n_der))
    
    select case (this%n_der)
    case (1)
        this%coef_der(-1) = -0.5d0
        this%coef_der(0)  =  0.0d0
        this%coef_der(1)  =  0.5d0
    case(2)
        this%coef_der(-2) =  1.0d0 / 12.0d0
        this%coef_der(-1) = -2.0d0 / 3.0d0
        this%coef_der(0)  =  0.0d0
        this%coef_der(1)  =  2.0d0 / 3.0d0
        this%coef_der(2)  = -1.0d0 / 12.0d0
    case (3)
        this%coef_der(-3) = -1.0d0 / 60.0d0
        this%coef_der(-2) =  3.0d0 / 20.0d0
        this%coef_der(-1) = -3.0d0 / 4.0d0
        this%coef_der(0)  =  0.0d0
        this%coef_der(1)  =  3.0d0 / 4.0d0
        this%coef_der(2)  = -3.0d0 / 20.0d0
        this%coef_der(3)  =  1.0d0 / 60.0d0
    case(4)
        this%coef_der(-4) =  1.0d0 / 280.0d0
        this%coef_der(-3) = -4.0d0 / 105.0d0
        this%coef_der(-2) =  1.0d0 / 5.0d0
        this%coef_der(-1) = -4.0d0 / 5.0d0
        this%coef_der(0)  =  0.0d0
        this%coef_der(1)  =  4.0d0 / 5.0d0
        this%coef_der(2)  = -1.0d0 / 5.0d0
        this%coef_der(3)  =  4.0d0 / 105.0d0
        this%coef_der(4)  = -1.0d0 / 280.0d0

    end select

    this%coef_der = this%coef_der / this%dr

    select case (this%dimensions)
    case (1)

        if (boundaries(1) == PML_BOUNDARIES) then
            this%cpml_pos(1) = .true.  ! left
            this%cpml_pos(2) = .true.  ! right
        end if

    case (2)

        if (boundaries(1) == PML_BOUNDARIES) then
            if (mpi_coords(1) == 0) then
                this%cpml_pos(1) = .true.  ! left
            end if
            if (mpi_coords(1) == mpi_dims(1)-1) then
                this%cpml_pos(2) = .true.  ! right
            end if
        end if

        if (boundaries(2) == PML_BOUNDARIES) then
            if (mpi_coords(2) == 0) then
                this%cpml_pos(3) = .true.  ! front
            end if
            if (mpi_coords(2) == mpi_dims(2)-1) then
                this%cpml_pos(4) = .true.  ! back
            end if
        end if

    case (3)

        if (boundaries(1) == PML_BOUNDARIES) then
            if (mpi_coords(1) == 0) then
                this%cpml_pos(1) = .true.  ! left
            end if
            if (mpi_coords(1) == mpi_dims(1)-1) then
                this%cpml_pos(2) = .true.  ! right
            end if
        end if

        if (boundaries(2) == PML_BOUNDARIES) then
            if (mpi_coords(2) == 0) then
                this%cpml_pos(3) = .true.  ! front
            end if
            if (mpi_coords(2) == mpi_dims(2)-1) then
                this%cpml_pos(4) = .true.  ! back
            end if
        end if

        if (boundaries(3) == PML_BOUNDARIES) then
            if (mpi_coords(3) == 0) then
                this%cpml_pos(5) = .true.  ! top
            end if
            if (mpi_coords(3) == mpi_dims(3)-1) then
                this%cpml_pos(6) = .true.  ! bottom
            end if
        end if

    end select

    sig_max = 0.8d0 * (m_pml + 1.0d0) / (this%dr * DSQRT(mu0/eps0))

    if (this%cpml_pos(1)) then
        do i = 1, n_pml
            sig   = sig_max * (DBLE(n_pml-i+1)/DBLE(n_pml))**m_pml
            kappa = 1.0d0 + (kappa_max_pml - 1.0d0) * ((DBLE(n_pml-i+1)/DBLE(n_pml)))**m_pml
            s_i   = Z_ONE*kappa + sig / (Z_I * this%w0 * eps0)
            this%s_inv_x(i) = (1.0d0 / s_i)
        end do
    end if

    if (this%cpml_pos(2)) then
        ii = this%nx - n_pml + 1
        do i = n_pml, 1, -1
            sig   = sig_max * ((DBLE(n_pml-i+1)/DBLE(n_pml)))**m_pml
            kappa = 1.0d0 + (kappa_max_pml - 1.0d0) * ((DBLE(n_pml-i+1)/DBLE(n_pml)))**m_pml
            s_i   = Z_ONE*kappa + sig / (Z_I * this%w0 * eps0)
            this%s_inv_x(ii) = (1.0d0 / s_i)
            ii = ii + 1
        end do
    end if

    if (this%cpml_pos(3)) then
        do i = 1, n_pml
            sig   = sig_max * (DBLE(n_pml-i+1)/DBLE(n_pml))**m_pml
            kappa = 1.0d0 + (kappa_max_pml - 1.0d0) * (DBLE(n_pml-i+1)/DBLE(n_pml))**m_pml
            s_i   = Z_ONE*kappa + sig / (Z_I * this%w0 * eps0)
            this%s_inv_y(i) = (1.0d0 / s_i)
        end do
    end if

    if (this%cpml_pos(4)) then
        ii = this%ny - n_pml + 1
        do i = n_pml, 1, -1
            sig   = sig_max * (DBLE(n_pml-i+1)/DBLE(n_pml))**m_pml
            kappa = 1.0d0 + (kappa_max_pml - 1.0d0) * (DBLE(n_pml-i+1)/DBLE(n_pml))**m_pml
            s_i   = Z_ONE*kappa + sig / (Z_I * this%w0 * eps0)
            this%s_inv_y(ii) = (1.0d0 / s_i)
            ii = ii + 1
        end do
    end if

    if (this%cpml_pos(5)) then
        do i = 1, n_pml
            sig   = sig_max * (DBLE(n_pml-i+1)/DBLE(n_pml))**m_pml
            kappa = 1.0d0 + (kappa_max_pml - 1.0d0) * (DBLE(n_pml-i+1)/DBLE(n_pml))**m_pml
            s_i   = Z_ONE*kappa + sig / (Z_I * this%w0 * eps0)
            this%s_inv_z(i) = (1.0d0 / s_i)
        end do
    end if

    if (this%cpml_pos(6)) then
        ii = this%nz - n_pml + 1
        do i = n_pml, 1, -1
            sig   = sig_max * (DBLE(n_pml-i+1)/DBLE(n_pml))**m_pml
            kappa = 1.0d0 + (kappa_max_pml - 1.0d0) * (DBLE(n_pml-i+1)/DBLE(n_pml))**m_pml
            s_i   = Z_ONE*kappa + sig / (Z_I * this%w0 * eps0)
            this%s_inv_z(ii) = (1.0d0 / s_i)
            ii = ii + 1
        end do
    end if

end subroutine init_operator

!###################################################################################################

subroutine kill_operator(this)
    class(TOperator), intent(inout) :: this

    if (allocated(this%s_inv_x)) deallocate(this%s_inv_x)
    if (allocated(this%s_inv_y)) deallocate(this%s_inv_y)
    if (allocated(this%s_inv_z)) deallocate(this%s_inv_z)
    if (allocated(this%coef_der)) deallocate(this%coef_der)

end subroutine kill_operator

!###################################################################################################

subroutine apply_operator(this, f_vec, Af_vec, eps_r, transpose)
    class(TOperator)    , intent(inout) :: this
    class(TRSvec)       , intent(inout) :: f_vec
    class(TRSvec)       , intent(inout) :: Af_vec
    class(TMedium_eps_r), intent(in)    :: eps_r
    logical             , intent(in)    :: transpose

    call extend_fvec_to_ranks(f_vec)

    select type(f_vec)
    type is (TRSvec_1D)
        select type(Af_vec)
        type is (TRSvec_1D)
            call apply_operator_1D(this, f_vec, Af_vec, eps_r, transpose)
        end select
    type is (TRSvec_2D)
        select type(Af_vec)
        type is (TRSvec_2D)
            call apply_operator_2D(this, f_vec, Af_vec, eps_r, transpose)
        end select

    type is (TRSvec_3D)
        select type(Af_vec)
        type is (TRSvec_3D)
            call apply_operator_3D(this, f_vec, Af_vec, eps_r, transpose)
        end select
    end select

end subroutine apply_operator

!###################################################################################################

subroutine apply_operator_1D(Aop, f_vec, Af_vec, eps_r, transpose)  
    type(TOperator)     , intent(inout) :: Aop
    class(TRSvec_1D)    , intent(inout) :: Af_vec
    class(TRSvec_1D)    , intent(inout) :: f_vec
    class(TMedium_eps_r), intent(in)    :: eps_r   
    logical             , intent(in)    :: transpose

    complex(dp) :: dfpydz, dfpxdz, dfmydz, dfmxdz
    complex(dp) :: C1, C2, C3, C4
    complex(dp) :: Ex
    integer     :: sign = 1
    integer     :: i, ii

    if (transpose) sign = -1

    C1 =  sign * c0 * Z_I
    C2 = -Aop%w0 * Z_I
    C3 = -Aop%w0 * Z_I / (2.0d0 * eps0)
    C4 =  Z_ONE * eps0

#ifdef USE_MPI
    ! If MPI is being used, the periodic condition were already applied.
#else
    if (Aop%boundaries(1) == PERIODIC_BOUNDARIES) then

        f_vec%pl_x(-3:0)                  = f_vec%pl_x(f_vec%nx-3:f_vec%nx)
        f_vec%pl_x(f_vec%nx+1:f_vec%nx+4) = f_vec%pl_x(1:4)
        f_vec%pl_y(-3:0)                  = f_vec%pl_y(f_vec%nx-3:f_vec%nx)
        f_vec%pl_y(f_vec%nx+1:f_vec%nx+4) = f_vec%pl_y(1:4)
        f_vec%mi_x(-3:0)                  = f_vec%mi_x(f_vec%nx-3:f_vec%nx)
        f_vec%mi_x(f_vec%nx+1:f_vec%nx+4) = f_vec%mi_x(1:4)
        f_vec%mi_y(-3:0)                  = f_vec%mi_y(f_vec%nx-3:f_vec%nx)
        f_vec%mi_y(f_vec%nx+1:f_vec%nx+4) = f_vec%mi_y(1:4)

    end if
#endif

    do i = 1, f_vec%nx
        
        dfpxdz = Z_0
        dfpydz = Z_0
        dfmxdz = Z_0
        dfmydz = Z_0

        do ii = -Aop%n_der, Aop%n_der
            dfpxdz = dfpxdz + Aop%coef_der(ii) * f_vec%pl_x(i+ii)
            dfpydz = dfpydz + Aop%coef_der(ii) * f_vec%pl_y(i+ii)
            dfmxdz = dfmxdz + Aop%coef_der(ii) * f_vec%mi_x(i+ii)
            dfmydz = dfmydz + Aop%coef_der(ii) * f_vec%mi_y(i+ii)
        end do

        Ex = f_vec%pl_x(i) + f_vec%mi_x(i)

        Af_vec%pl_x(i) = -C1 * dfpydz * Aop%s_inv_x(i) + C2 * f_vec%pl_x(i) + &
                            C3 * (eps_r%mat1D(i) - C4) * Ex
        Af_vec%pl_y(i) =  C1 * dfpxdz * Aop%s_inv_x(i) +  C2 * f_vec%pl_y(i)
        Af_vec%mi_x(i) =  C1 * dfmydz * Aop%s_inv_x(i) + C2 * f_vec%mi_x(i) + C3 * (eps_r%mat1D(i) - C4) * Ex
        Af_vec%mi_y(i) = -C1 * dfmxdz * Aop%s_inv_x(i) + C2 * f_vec%mi_y(i)

    end do
    
end subroutine apply_operator_1D

!###################################################################################################
subroutine apply_operator_2D(Aop, f_vec, Af_vec, eps_r, transpose)  
    type(TOperator)     , intent(inout) :: Aop  
    class(TRSvec_2D)    , intent(inout) :: Af_vec
    class(TRSvec_2D)    , intent(inout) :: f_vec
    class(TMedium_eps_r), intent(in)    :: eps_r   
    logical             , intent(in)    :: transpose

    complex(dp) :: dfpzdy, dfpzdx, dfpxdy, dfpydx
    complex(dp) :: dfmzdy, dfmzdx, dfmxdy, dfmydx
    complex(dp) :: C1, C2, C3, C4
    complex(dp) :: Ex, Ey, Ez
    integer     :: sign = 1
    integer     :: i, ii
    integer     :: j, jj

    if (transpose) sign = -1

    C1 =  sign * c0 * Z_I
    C2 = -Aop%w0 * Z_I
    C3 = -Aop%w0 * Z_I / (2.0d0 * eps0)
    C4 =  Z_ONE * eps0

#ifdef USE_MPI
    ! If MPI is being used, the periodic condition were already applied.
#else
    if (Aop%boundaries(1) == PERIODIC_BOUNDARIES) then

        f_vec%pl_x(-3:0, :)                  = f_vec%pl_x(f_vec%nx-3:f_vec%nx, :)
        f_vec%pl_x(f_vec%nx+1:f_vec%nx+4, :) = f_vec%pl_x(1:4, :)
        f_vec%pl_y(-3:0, :)                  = f_vec%pl_y(f_vec%nx-3:f_vec%nx, :)
        f_vec%pl_y(f_vec%nx+1:f_vec%nx+4, :) = f_vec%pl_y(1:4, :)
        f_vec%pl_z(-3:0, :)                  = f_vec%pl_z(f_vec%nx-3:f_vec%nx, :)
        f_vec%pl_z(f_vec%nx+1:f_vec%nx+4, :) = f_vec%pl_z(1:4, :)
        f_vec%mi_x(-3:0, :)                  = f_vec%mi_x(f_vec%nx-3:f_vec%nx, :)
        f_vec%mi_x(f_vec%nx+1:f_vec%nx+4, :) = f_vec%mi_x(1:4, :)
        f_vec%mi_y(-3:0, :)                  = f_vec%mi_y(f_vec%nx-3:f_vec%nx, :)
        f_vec%mi_y(f_vec%nx+1:f_vec%nx+4, :) = f_vec%mi_y(1:4, :)
        f_vec%mi_z(-3:0, :)                  = f_vec%mi_z(f_vec%nx-3:f_vec%nx, :)
        f_vec%mi_z(f_vec%nx+1:f_vec%nx+4, :) = f_vec%mi_z(1:4, :)

    end if

    if (Aop%boundaries(2) == PERIODIC_BOUNDARIES) then

        f_vec%pl_x(:, -3:0)                  = f_vec%pl_x(:, f_vec%ny-3:f_vec%ny)
        f_vec%pl_x(:, f_vec%ny+1:f_vec%ny+4) = f_vec%pl_x(:, 1:4)
        f_vec%pl_y(:, -3:0)                  = f_vec%pl_y(:, f_vec%ny-3:f_vec%ny)
        f_vec%pl_y(:, f_vec%ny+1:f_vec%ny+4) = f_vec%pl_y(:, 1:4)
        f_vec%pl_z(:, -3:0)                  = f_vec%pl_z(:, f_vec%ny-3:f_vec%ny)
        f_vec%pl_z(:, f_vec%ny+1:f_vec%ny+4) = f_vec%pl_z(:, 1:4)
        f_vec%mi_x(:, -3:0)                  = f_vec%mi_x(:, f_vec%ny-3:f_vec%ny)
        f_vec%mi_x(:, f_vec%ny+1:f_vec%ny+4) = f_vec%mi_x(:, 1:4)
        f_vec%mi_y(:, -3:0)                  = f_vec%mi_y(:, f_vec%ny-3:f_vec%ny)
        f_vec%mi_y(:, f_vec%ny+1:f_vec%ny+4) = f_vec%mi_y(:, 1:4)
        f_vec%mi_z(:, -3:0)                  = f_vec%mi_z(:, f_vec%ny-3:f_vec%ny)
        f_vec%mi_z(:, f_vec%ny+1:f_vec%ny+4) = f_vec%mi_z(:, 1:4)

    end if
#endif


    do j = 1, f_vec%ny
    do i = 1, f_vec%nx
    
        dfpxdy = Z_0
        dfpydx = Z_0
        dfpzdx = Z_0
        dfpzdy = Z_0

        dfmxdy = Z_0
        dfmydx = Z_0
        dfmzdx = Z_0
        dfmzdy = Z_0

        do ii = -Aop%n_der, Aop%n_der
            dfpydx = dfpydx + Aop%coef_der(ii) * f_vec%pl_y(i+ii,j)
            dfpzdx = dfpzdx + Aop%coef_der(ii) * f_vec%pl_z(i+ii,j)
            dfmydx = dfmydx + Aop%coef_der(ii) * f_vec%mi_y(i+ii,j)
            dfmzdx = dfmzdx + Aop%coef_der(ii) * f_vec%mi_z(i+ii,j)
        end do

        do jj = -Aop%n_der, Aop%n_der
            dfpxdy = dfpxdy + Aop%coef_der(jj) * f_vec%pl_x(i,j+jj)
            dfpzdy = dfpzdy + Aop%coef_der(jj) * f_vec%pl_z(i,j+jj)
            dfmxdy = dfmxdy + Aop%coef_der(jj) * f_vec%mi_x(i,j+jj)
            dfmzdy = dfmzdy + Aop%coef_der(jj) * f_vec%mi_z(i,j+jj)
        end do

        Ex = f_vec%pl_x(i,j) + f_vec%mi_x(i,j)
        Ey = f_vec%pl_y(i,j) + f_vec%mi_y(i,j)
        Ez = f_vec%pl_z(i,j) + f_vec%mi_z(i,j)

        dfpxdy = dfpxdy * Aop%s_inv_y(j)
        dfpydx = dfpydx * Aop%s_inv_x(i)
        dfpzdx = dfpzdx * Aop%s_inv_x(i)
        dfpzdy = dfpzdy * Aop%s_inv_y(j)

        dfmxdy = dfmxdy * Aop%s_inv_y(j)
        dfmydx = dfmydx * Aop%s_inv_x(i)
        dfmzdx = dfmzdx * Aop%s_inv_x(i)
        dfmzdy = dfmzdy * Aop%s_inv_y(j)

        Af_vec%pl_x(i,j) =  -C1 * dfpzdy + C2 * f_vec%pl_x(i,j) + &
                                C3 * (eps_r%mat2D(i,j) - C4) * Ex
        Af_vec%pl_y(i,j) =   C1 * dfpzdx + C2 * f_vec%pl_y(i,j) + &
                                C3 * (eps_r%mat2D(i,j) - C4) * Ey
        Af_vec%pl_z(i,j) =  -C1 * (dfpydx - dfpxdy) + C2 * f_vec%pl_z(i,j) + &
                                C3 * (eps_r%mat2D(i,j) - C4) * Ez
        Af_vec%mi_x(i,j) =   C1 * dfmzdy + C2 * f_vec%mi_x(i,j) + &
                                C3 * (eps_r%mat2D(i,j) - C4) * Ex
        Af_vec%mi_y(i,j) =  -C1 * dfmzdx + C2 * f_vec%mi_y(i,j) + &
                                C3 * (eps_r%mat2D(i,j) - C4) * Ey
        Af_vec%mi_z(i,j) =   C1 * (dfmydx - dfmxdy) + C2 * f_vec%mi_z(i,j) + &
                                C3 * (eps_r%mat2D(i,j) - C4) * Ez

    end do
    end do

end subroutine apply_operator_2D

!###################################################################################################

subroutine apply_operator_3D(Aop, f_vec, Af_vec, eps_r, transpose)  
    type(TOperator)     , intent(inout) :: Aop  
    class(TRSvec_3D)    , intent(inout) :: Af_vec
    class(TRSvec_3D)    , intent(inout) :: f_vec
    class(TMedium_eps_r), intent(in)    :: eps_r   
    logical             , intent(in)    :: transpose

    complex(dp) :: dfpydz, dfpzdx, dfpxdy, dfpydx, dfpzdy, dfpxdz
    complex(dp) :: dfmydz, dfmzdx, dfmxdy, dfmydx, dfmzdy, dfmxdz
    complex(dp) :: C1, C2, C3, C4
    complex(dp) :: Ex, Ey, Ez
    integer     :: sign = 1
    integer     :: i, ii
    integer     :: j, jj
    integer     :: k, kk

    if (transpose) sign = -1

    C1 =  sign * c0 * Z_I
    C2 = -Aop%w0 * Z_I
    C3 = -Aop%w0 * Z_I / (2.0d0 * eps0)
    C4 =  Z_ONE * eps0

#ifdef USE_MPI
    ! If MPI is being used, the periodic condition were already applied.
#else
    if (Aop%boundaries(1) == PERIODIC_BOUNDARIES) then
        f_vec%pl_x(-3:0, :, :)                  = f_vec%pl_x(f_vec%nx-3:f_vec%nx, :, :)
        f_vec%pl_x(f_vec%nx+1:f_vec%nx+4, :, :) = f_vec%pl_x(1:4, :, :)
        f_vec%pl_y(-3:0, :, :)                  = f_vec%pl_y(f_vec%nx-3:f_vec%nx, :, :)
        f_vec%pl_y(f_vec%nx+1:f_vec%nx+4, :, :) = f_vec%pl_y(1:4, :, :)
        f_vec%pl_z(-3:0, :, :)                  = f_vec%pl_z(f_vec%nx-3:f_vec%nx, :, :)
        f_vec%pl_z(f_vec%nx+1:f_vec%nx+4, :, :) = f_vec%pl_z(1:4, :, :)
        f_vec%mi_x(-3:0, :, :)                  = f_vec%mi_x(f_vec%nx-3:f_vec%nx, :, :)
        f_vec%mi_x(f_vec%nx+1:f_vec%nx+4, :, :) = f_vec%mi_x(1:4, :, :)
        f_vec%mi_y(-3:0, :, :)                  = f_vec%mi_y(f_vec%nx-3:f_vec%nx, :, :)
        f_vec%mi_y(f_vec%nx+1:f_vec%nx+4, :, :) = f_vec%mi_y(1:4, :, :)
        f_vec%mi_z(-3:0, :, :)                  = f_vec%mi_z(f_vec%nx-3:f_vec%nx, :, :)
        f_vec%mi_z(f_vec%nx+1:f_vec%nx+4, :, :) = f_vec%mi_z(1:4, :, :)
    end if

    if (Aop%boundaries(2) == PERIODIC_BOUNDARIES) then
        f_vec%pl_x(:, -3:0, :)                  = f_vec%pl_x(:, f_vec%ny-3:f_vec%ny, :)
        f_vec%pl_x(:, f_vec%ny+1:f_vec%ny+4, :) = f_vec%pl_x(:, 1:4, :)
        f_vec%pl_y(:, -3:0, :)                  = f_vec%pl_y(:, f_vec%ny-3:f_vec%ny, :)
        f_vec%pl_y(:, f_vec%ny+1:f_vec%ny+4, :) = f_vec%pl_y(:, 1:4, :)
        f_vec%pl_z(:, -3:0, :)                  = f_vec%pl_z(:, f_vec%ny-3:f_vec%ny, :)
        f_vec%pl_z(:, f_vec%ny+1:f_vec%ny+4, :) = f_vec%pl_z(:, 1:4, :)
        f_vec%mi_x(:, -3:0, :)                  = f_vec%mi_x(:, f_vec%ny-3:f_vec%ny, :)
        f_vec%mi_x(:, f_vec%ny+1:f_vec%ny+4, :) = f_vec%mi_x(:, 1:4, :)
        f_vec%mi_y(:, -3:0, :)                  = f_vec%mi_y(:, f_vec%ny-3:f_vec%ny, :)
        f_vec%mi_y(:, f_vec%ny+1:f_vec%ny+4, :) = f_vec%mi_y(:, 1:4, :)
        f_vec%mi_z(:, -3:0, :)                  = f_vec%mi_z(:, f_vec%ny-3:f_vec%ny, :)
        f_vec%mi_z(:, f_vec%ny+1:f_vec%ny+4, :) = f_vec%mi_z(:, 1:4, :)
    end if

    if (Aop%boundaries(3) == PERIODIC_BOUNDARIES) then
        f_vec%pl_x(:, :, -3:0)                  = f_vec%pl_x(:, :, f_vec%nz-3:f_vec%nz)
        f_vec%pl_x(:, :, f_vec%nz+1:f_vec%nz+4) = f_vec%pl_x(:, :, 1:4)
        f_vec%pl_y(:, :, -3:0)                  = f_vec%pl_y(:, :, f_vec%nz-3:f_vec%nz)
        f_vec%pl_y(:, :, f_vec%nz+1:f_vec%nz+4) = f_vec%pl_y(:, :, 1:4)
        f_vec%pl_z(:, :, -3:0)                  = f_vec%pl_z(:, :, f_vec%nz-3:f_vec%nz)
        f_vec%pl_z(:, :, f_vec%nz+1:f_vec%nz+4) = f_vec%pl_z(:, :, 1:4)
        f_vec%mi_x(:, :, -3:0)                  = f_vec%mi_x(:, :, f_vec%nz-3:f_vec%nz)
        f_vec%mi_x(:, :, f_vec%nz+1:f_vec%nz+4) = f_vec%mi_x(:, :, 1:4)
        f_vec%mi_y(:, :, -3:0)                  = f_vec%mi_y(:, :, f_vec%nz-3:f_vec%nz)
        f_vec%mi_y(:, :, f_vec%nz+1:f_vec%nz+4) = f_vec%mi_y(:, :, 1:4)
        f_vec%mi_z(:, :, -3:0)                  = f_vec%mi_z(:, :, f_vec%nz-3:f_vec%nz)
        f_vec%mi_z(:, :, f_vec%nz+1:f_vec%nz+4) = f_vec%mi_z(:, :, 1:4)
    end if

#endif


    do k = 1, f_vec%nz
    do j = 1, f_vec%ny
    do i = 1, f_vec%nx
    
        dfpxdy = Z_0
        dfpxdz = Z_0
        dfpydx = Z_0
        dfpydz = Z_0
        dfpzdx = Z_0
        dfpzdy = Z_0

        dfmxdz = Z_0
        dfmxdy = Z_0
        dfmydz = Z_0
        dfmydx = Z_0
        dfmzdx = Z_0
        dfmzdy = Z_0

        do ii = -Aop%n_der, Aop%n_der
            dfpydx = dfpydx + Aop%coef_der(ii) * f_vec%pl_y(i+ii,j,k)
            dfpzdx = dfpzdx + Aop%coef_der(ii) * f_vec%pl_z(i+ii,j,k)
            dfmydx = dfmydx + Aop%coef_der(ii) * f_vec%mi_y(i+ii,j,k)
            dfmzdx = dfmzdx + Aop%coef_der(ii) * f_vec%mi_z(i+ii,j,k)
        end do

        do jj = -Aop%n_der, Aop%n_der
            dfpxdy = dfpxdy + Aop%coef_der(jj) * f_vec%pl_x(i,j+jj,k)
            dfpzdy = dfpzdy + Aop%coef_der(jj) * f_vec%pl_z(i,j+jj,k)
            dfmxdy = dfmxdy + Aop%coef_der(jj) * f_vec%mi_x(i,j+jj,k)
            dfmzdy = dfmzdy + Aop%coef_der(jj) * f_vec%mi_z(i,j+jj,k)
        end do

        do kk = -Aop%n_der, Aop%n_der
            dfpxdz = dfpxdz + Aop%coef_der(kk) * f_vec%pl_x(i,j,k+kk)
            dfpydz = dfpydz + Aop%coef_der(kk) * f_vec%pl_y(i,j,k+kk)
            dfmxdz = dfmxdz + Aop%coef_der(kk) * f_vec%mi_x(i,j,k+kk)
            dfmydz = dfmydz + Aop%coef_der(kk) * f_vec%mi_y(i,j,k+kk)
        end do
    
        Ex = f_vec%pl_x(i,j,k) + f_vec%mi_x(i,j,k)
        Ey = f_vec%pl_y(i,j,k) + f_vec%mi_y(i,j,k)
        Ez = f_vec%pl_z(i,j,k) + f_vec%mi_z(i,j,k)

        dfpxdy = dfpxdy * Aop%s_inv_y(j)
        dfpxdz = dfpxdz * Aop%s_inv_z(k)
        dfpydx = dfpydx * Aop%s_inv_x(i)
        dfpydz = dfpydz * Aop%s_inv_z(k)
        dfpzdx = dfpzdx * Aop%s_inv_x(i)
        dfpzdy = dfpzdy * Aop%s_inv_y(j)
        dfmxdy = dfmxdy * Aop%s_inv_y(j)
        dfmxdz = dfmxdz * Aop%s_inv_z(k)
        dfmydz = dfmydz * Aop%s_inv_z(k)
        dfmydx = dfmydx * Aop%s_inv_x(i)
        dfmzdx = dfmzdx * Aop%s_inv_x(i)
        dfmzdy = dfmzdy * Aop%s_inv_y(j)

        Af_vec%pl_x(i,j,k) =  C1 * (dfpzdy - dfpydz) + C2 * f_vec%pl_x(i,j,k) + &
                                C3 * (eps_r%mat3D(i,j,k) - C4) * Ex
        Af_vec%pl_y(i,j,k) =  C1 * (dfpxdz - dfpzdx) + C2 * f_vec%pl_y(i,j,k) + &
                                C3 * (eps_r%mat3D(i,j,k) - C4) * Ey
        Af_vec%pl_z(i,j,k) =  C1 * (dfpydx - dfpxdy) + C2 * f_vec%pl_z(i,j,k) + &
                                C3 * (eps_r%mat3D(i,j,k) - C4) * Ez
        Af_vec%mi_x(i,j,k) = -C1 * (dfmzdy - dfmydz) + C2 * f_vec%mi_x(i,j,k) + &
                                C3 * (eps_r%mat3D(i,j,k) - C4) * Ex
        Af_vec%mi_y(i,j,k) = -C1 * (dfmxdz - dfmzdx) + C2 * f_vec%mi_y(i,j,k) + &
                                C3 * (eps_r%mat3D(i,j,k) - C4) * Ey
        Af_vec%mi_z(i,j,k) = -C1 * (dfmydx - dfmxdy) + C2 * f_vec%mi_z(i,j,k) + &
                                C3 * (eps_r%mat3D(i,j,k) - C4) * Ez

    end do
    end do
    end do


end subroutine apply_operator_3D

!###################################################################################################

end module operator_mod