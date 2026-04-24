module optimization_problem_mod

    use constants_mod
    use operator_mod
    use target_source_mod
    use medium_mod
    use rs_vec_base_mod
    use rs_vec_dimensions_mod
    use allocator_multidim_mod
    use bicgstab_mod

    implicit none

    type TOptPrblm  

        logical             :: is_complex
        type(TOperator)     :: A_op
        type(TMedium_eps_r) :: eps_r
        integer             :: id
        integer             :: n_src
        integer             :: n_trg
        integer             :: n_src_trg
        integer             :: dimensions
        integer             :: nx
        integer             :: ny
        integer             :: nz
        real(dp)            :: eps_Re
        real(dp)            :: eps_Im
        real(dp)            :: kappa_max
        real(dp)            :: eta_max
        real(dp)            :: eta_min
        real(dp)            :: sigma_ker
        real(dp)            :: eta_rho
        real(dp)            :: beta_rho
        real(dp)            :: w_total

        class(TRSvec), allocatable :: f_vec
        class(TRSvec), allocatable :: f_vec_new
        class(TRSvec), allocatable :: f_adj_vec(:)
        class(TRSvec), allocatable :: f_adj_vec_new(:)
        class(TRSvec), allocatable :: Af_vec
        class(TRSvec), allocatable :: j_src
        class(TRSvec), allocatable :: j_trg
        type(TSrcTrg), allocatable :: src_trg(:)
        logical      , allocatable :: opt_region(:,:,:)
        complex(dp)  , allocatable :: dA(:,:,:)
        real(dp)     , allocatable :: grad(:,:,:)
        real(dp)     , allocatable :: w_dL(:) ! Weigh of the objective 
                                              ! function gradient

    contains

        procedure :: init_optprblm
        procedure :: kill_optprblm
        procedure :: solve_forward
        procedure :: solve_adjoint
        procedure :: compute_gradient
        procedure :: update_fields
        procedure :: update_permittivity
        procedure :: update_targets

    end type TOptPrblm

contains

!###################################################################################################

subroutine init_optprblm(this, id,dimensions, dr, freq, eps_Re, eps_Im, boundaries, restart, &
                         n_pml, grid_Ndims, mpi_cart_comm, mpi_coords, mpi_dims)

    class(TOptPrblm) , intent(inout) :: this
    logical          , intent(in)    :: restart
    integer          , intent(in)    :: dimensions
    integer          , intent(in)    :: id
    integer          , intent(in)    :: boundaries(3)
    integer          , intent(in)    :: grid_Ndims(3)
    integer          , intent(in)    :: mpi_cart_comm
    integer          , intent(in)    :: mpi_coords(3)
    integer          , intent(in)    :: mpi_dims(3)
    integer          , intent(in)    :: n_pml
    real(dp)         , intent(in)    :: dr
    real(dp)         , intent(in)    :: freq
    real(dp)         , intent(in)    :: eps_Re
    real(dp)         , intent(in)    :: eps_Im

    integer  :: i
    real(dp) :: eps_amp

    this%dimensions = dimensions
    this%nx         = grid_Ndims(1)
    this%ny         = grid_Ndims(2)
    this%nz         = grid_Ndims(3)

    if (dimensions == 1) this%nz = this%nx

    this%id         = id
    this%eps_Re     = eps_Re
    this%eps_Im     = eps_Im

    if (eps_Im /= 0.0d0) then
        this%is_complex  = .true.
        eps_amp          = SQRT(eps_Re**2 + eps_Im**2)
        this%eta_max     = DSQRT((eps_amp + eps_Re)/(2.0d0))
        this%eta_min     = DSQRT(eps0)
        this%kappa_max   = DSQRT((eps_amp - eps_Re)/(2.0d0))
    else
        this%is_complex = .false.
    end if
    
    this%f_vec     = rs_vec_factory(dimensions)
    this%f_vec_new = rs_vec_factory(dimensions)
    this%Af_vec    = rs_vec_factory(dimensions)
    this%j_src     = rs_vec_factory(dimensions)
    this%j_trg     = rs_vec_factory(dimensions)

    call this%f_vec%init(grid_Ndims, dr, dimensions, freq)
    call this%f_vec_new%init(grid_Ndims, dr, dimensions, freq)
    call this%Af_vec%init(grid_Ndims, dr, dimensions, freq)
    call this%j_src%init(grid_Ndims, dr, dimensions, freq)
    call this%j_trg%init(grid_Ndims, dr, dimensions, freq)

    call this%A_op%init_operator(dr, freq, dimensions, grid_Ndims, n_pml, boundaries, &
                                 mpi_dims, mpi_coords)

    call allocate_multidim(array = this%opt_region, dim = dimensions, i_max = this%nx, &
                           j_max = this%ny, k_max = this%nz)

    call allocate_multidim(array = this%grad, dim = dimensions, i_max = this%nx, &
                           j_max = this%ny, k_max = this%nz)

    call allocate_multidim(array = this%dA, dim = dimensions, i_max = this%nx, &
                           j_max = this%ny, k_max = this%nz)

    this%opt_region = .false.
    this%grad = 0.0_dp
    this%dA   = Z_0
    
    call this%eps_r%init_medium(dimensions, grid_Ndims)
    call this%eps_r%read_medium(id, dr, grid_Ndims, mpi_coords, mpi_dims)

    call read_opt_region(this, dr, grid_Ndims, mpi_coords, mpi_dims)
    call clean_opt_region(this)

    call read_init_src_trg(this%src_trg, this%n_src_trg, this%n_src, this%n_trg, dimensions, id, dr, &
                           grid_Ndims, mpi_cart_comm, mpi_coords, mpi_dims)

    do i = 1, this%n_src_trg
        call set_source_J(this%src_trg(i), this%j_trg)
    end do

    this%f_adj_vec     = rs_vec_factory(dimensions, this%n_trg)
    this%f_adj_vec_new = rs_vec_factory(dimensions, this%n_trg)

    do i =1, this%n_trg
        call this%f_adj_vec(i)%init(grid_Ndims, dr, dimensions, freq)
        call this%f_adj_vec_new(i)%init(grid_Ndims, dr, dimensions, freq)
    end do

    if (.not. allocated(this%w_dL)) allocate(this%w_dL(this%n_trg))

end subroutine init_optprblm

!###################################################################################################

subroutine kill_optprblm(this)

    class(TOptPrblm) , intent(inout) :: this

    integer :: i

    call this%A_op%kill_operator()

    call this%f_vec%kill()
    call this%f_vec_new%kill()
    call this%Af_vec%kill()
    call this%j_src%kill()
    call this%j_trg%kill()

    if (allocated(this%opt_region)) deallocate(this%opt_region)
    if (allocated(this%w_dL))       deallocate(this%w_dL)
    if (allocated(this%grad))      deallocate(this%grad)
    if (allocated(this%dA))         deallocate(this%dA)
    call this%eps_r%kill_medium()

    do i = 1, this%n_src_trg
        call this%src_trg(i)%kill_src_trg()
    end do

    if (allocated(this%src_trg)) deallocate(this%src_trg)
    
    do i = 1, this%n_trg
        call this%f_adj_vec(i)%kill()
        call this%f_adj_vec_new(i)%kill()
    end do

    if (allocated(this%f_adj_vec))     deallocate(this%f_adj_vec)
    if (allocated(this%f_adj_vec_new)) deallocate(this%f_adj_vec_new)

end subroutine kill_optprblm

!###################################################################################################

subroutine solve_forward(this, bicgstab_tol, bicgstab_max_iter, bicgstab_L_term, &
                         delta_p, opt_step, converged)

    class(TOptPrblm) , intent(inout) :: this
    integer          , intent(in)    :: bicgstab_max_iter
    integer          , intent(in)    :: bicgstab_L_term
    logical          , intent(out)   :: converged
    real(dp)         , intent(in)    :: bicgstab_tol
    integer          , intent(in)    :: opt_step
    real(dp)         , intent(in)    :: delta_p
    
    logical    :: transpose
    real(dp)   :: bicgstab_error

    transpose = .false.

    call BICGStab_L(this%A_op, this%f_vec, this%j_src, this%f_vec_new, this%Af_vec, &
                    this%eps_r, converged, transpose, bicgstab_tol, bicgstab_max_iter, &
                    bicgstab_L_term, bicgstab_error)

    if (.not. converged) then
        write(*, *) "Warning, bicgstab-L has not converged F vector"
        write(*, *) "Problem:",this%id,"drho:",delta_p,"Step:", opt_step
        write(*, *) "Total Error =", bicgstab_error
    end if

end subroutine solve_forward

!###################################################################################################

subroutine solve_adjoint(this, bicgstab_tol, bicgstab_max_iter, bicgstab_L_term, &
                         delta_p, opt_step, converged)

    class(TOptPrblm) , intent(inout) :: this
    integer          , intent(in)    :: bicgstab_max_iter
    integer          , intent(in)    :: bicgstab_L_term
    logical          , intent(out)   :: converged
    real(dp)         , intent(in)    :: bicgstab_tol
    integer          , intent(in)    :: opt_step
    real(dp)         , intent(in)    :: delta_p
    
    integer    :: n
    logical    :: transpose
    real(dp)   :: bicgstab_error

    transpose = .true.

    do n = 1, this%n_src_trg
        if (this%src_trg(n)%is_a_target) then

            call set_jtrg(this%src_trg(n), this%j_trg)

            call BICGStab_L(this%A_op, this%f_adj_vec(n), this%j_trg, this%f_adj_vec_new(n), &
                            this%Af_vec, this%eps_r, converged, transpose, bicgstab_tol,     &
                            bicgstab_max_iter, bicgstab_L_term, bicgstab_error)

            if (.not. converged) then
                write(*, *) "Warning, bicgstab-L has not converged adjoint F vector"
                write(*, *) "Problem:",this%id,"drho:",delta_p,"Step:", opt_step
                write(*, *) "Target:", n
                write(*, *) "Total Error =", bicgstab_error
            end if

        end if
    end do
        
end subroutine solve_adjoint

!###################################################################################################

subroutine compute_gradient(this, rho_conv)

    class(TOptPrblm) , intent(inout) :: this
    real(dp)         , intent(in)    :: rho_conv(:,:,:)

    complex(dp)   :: deps_drho
    real(dp)      :: eps_r
    real(dp)      :: func_rho
    real(dp)      :: df_drho
    real(dp)      :: eta
    real(dp)      :: kappa
    real(dp)      :: beta_p
    real(dp)      :: eta_p
    real(dp)      :: C1
    real(dp)      :: C2
    real(dp)      :: w0
    integer       :: i, j, k, n
    integer       :: nx, ny, nz

    beta_p = this%beta_rho
    eta_p  = this%eta_rho
    eps_r  = this%eps_Re

    w0 = this%f_vec%freq

    C1 = 1.0d0 / (DTANH(beta_p*eta_p) + DTANH(beta_p*(1.0-eta_p)))
    C2 = DTANH(beta_p*eta_p)*C1

    select type (f_vec => this%f_vec)
    type is (TRSvec_1D)
    select type (f_adj_vec => this%f_adj_vec)
    type is (TRSvec_1D)
        nz = this%nz

        if(this%is_complex) then
            do k =1, nz
                df_drho   = C1*beta_p*(1.0/DCOSH(beta_p*(rho_conv(k,1,1)-eta_p))**2)
                func_rho  = C1*DTANH(beta_p*(rho_conv(k,1,1)-eta_p)) + C2
                eta       = (this%eta_max - this%eta_min)*func_rho + this%eta_min
                kappa     = this%kappa_max*func_rho
                deps_drho = 2.0d0*eps0*df_drho*((this%eta_max-this%eta_min)*(Z_ONE*eta-Z_I*kappa)+ &
                                                 this%kappa_max*(-Z_ONE*kappa-Z_I*eta))
                this%dA(k,1,1) = -Z_I * w0 * deps_drho 
            end do
        else
            do k = 1, nz
                df_drho   =   C1*beta_p*(1.0/DCOSH(beta_p*(rho_conv(k,1,1)-eta_p))**2)
                deps_drho =   Z_ONE * df_drho * (eps_r*eps0 - eps0)
                this%dA(k,1,1) = -Z_I * w0 * deps_drho
            end do
        end if

        do k = 1, nz
            this%grad(k,1,1) = 0.0d0
            do n = 1, this%n_trg

                this%grad(k,1,1) = this%grad(k,1,1) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%pl_x(k) * this%dA(k,1,1) * &
                (f_vec%pl_x(k) + f_vec%mi_x(k)), kind=dp)

                this%grad(k,1,1) = this%grad(k,1,1) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%pl_y(k) * this%dA(k,1,1) * &
                (f_vec%pl_y(k) + f_vec%mi_y(k)), kind=dp)

                this%grad(k,1,1) = this%grad(k,1,1) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%mi_x(k) * this%dA(k,1,1) * &
                (f_vec%pl_x(k) + f_vec%mi_x(k)), kind=dp)

                this%grad(k,1,1) = this%grad(k,1,1) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%mi_y(k) * this%dA(k,1,1) * &
                (f_vec%pl_y(k) + f_vec%mi_y(k)), kind=dp)

            end do  
        end do
    end select

    type is (TRSvec_2D)
    select type (f_adj_vec => this%f_adj_vec)
    type is (TRSvec_2D)

        nx = this%nx
        ny = this%ny

        if (this%is_complex) then
            do j = 1, ny
            do i = 1, nx
                df_drho   = C1*beta_p*(1.0/DCOSH(beta_p*(rho_conv(i,j,1)-eta_p))**2)
                func_rho  = C1*DTANH(beta_p*(rho_conv(i,j,1)-eta_p)) + C2
                eta       = (this%eta_max - this%eta_min)*func_rho + this%eta_min
                kappa     = this%kappa_max*func_rho
                deps_drho = 2.0d0*eps0*df_drho*((this%eta_max-this%eta_min)*(Z_ONE*eta-Z_I*kappa)+ &
                                                    this%kappa_max*(-Z_ONE*kappa-Z_I*eta))
                this%dA(i,j,1) = -Z_I * w0 * deps_drho 
            end do
            end do
        else
            do j = 1, ny
            do i = 1, nx
                df_drho   =   C1*beta_p*(1.0/DCOSH(beta_p*(rho_conv(i,j,1)-eta_p))**2)
                deps_drho =   Z_ONE * df_drho * (eps_r*eps0 - eps0)
                this%dA(i,j,1) = -Z_I * w0 * deps_drho 
            end do
            end do

        end if

        do j = 1, ny
        do i = 1, nx
            this%grad(i,j,1) = 0.0d0
            do n = 1, this%n_trg

                this%grad(i,j,1) = this%grad(i,j,1) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%pl_x(i,j) * this%dA(i,j,1) * &
                (f_vec%pl_x(i,j) + f_vec%mi_x(i,j)), kind=dp)

                this%grad(i,j,1) = this%grad(i,j,1) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%pl_y(i,j) * this%dA(i,j,1) * &
                (f_vec%pl_y(i,j) + f_vec%mi_y(i,j)), kind=dp)

                this%grad(i,j,1) = this%grad(i,j,1) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%pl_z(i,j) * this%dA(i,j,1) * &
                (f_vec%pl_z(i,j) + f_vec%mi_z(i,j)), kind=dp)

                this%grad(i,j,1) = this%grad(i,j,1) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%mi_x(i,j) * this%dA(i,j,1) * &
                (f_vec%pl_x(i,j) + f_vec%mi_x(i,j)), kind=dp)

                this%grad(i,j,1) = this%grad(i,j,1) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%mi_y(i,j) * this%dA(i,j,1) * &
                (f_vec%pl_y(i,j) + f_vec%mi_y(i,j)), kind=dp)

                this%grad(i,j,1) = this%grad(i,j,1) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%mi_z(i,j) * this%dA(i,j,1) * &
                (f_vec%pl_z(i,j) + f_vec%mi_z(i,j)), kind=dp)

            end do
        end do
        end do
    end select

    type is (TRSvec_3D)
    select type (f_adj_vec => this%f_adj_vec)
    type is (TRSvec_3D)

        nx = this%nx
        ny = this%ny
        nz = this%nz

        if (this%is_complex) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                df_drho   = C1*beta_p*(1.0/DCOSH(beta_p*(rho_conv(i,j,k)-eta_p))**2)
                func_rho  = C1*DTANH(beta_p*(rho_conv(i,j,k)-eta_p)) + C2
                eta       = (this%eta_max - this%eta_min)*func_rho + this%eta_min
                kappa     = this%kappa_max*func_rho
                deps_drho = 2.0d0*eps0*df_drho*((this%eta_max-this%eta_min)*(Z_ONE*eta-Z_I*kappa)+ &
                                                    this%kappa_max*(-Z_ONE*kappa-Z_I*eta))
                this%dA(i,j,k) = -Z_I * w0 * deps_drho 
            end do
            end do
            end do
        else
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                df_drho   =   C1*beta_p*(1.0/DCOSH(beta_p*(rho_conv(i,j,k)-eta_p))**2)
                deps_drho =   Z_ONE * df_drho * (eps_r*eps0 - eps0)
                this%dA(i,j,k) = -Z_I * w0 * deps_drho 
            end do
            end do
            end do
        end if
    
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            this%grad(i,j,k) = 0.0d0
            do n = 1, this%n_trg

                this%grad(i,j,k) = this%grad(i,j,k) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%pl_x(i,j,k) * this%dA(i,j,k) * &
                (f_vec%pl_x(i,j,k) + f_vec%mi_x(i,j,k)), kind=dp)

                this%grad(i,j,k) = this%grad(i,j,k) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%pl_y(i,j,k) * this%dA(i,j,k) * &
                (f_vec%pl_y(i,j,k) + f_vec%mi_y(i,j,k)), kind=dp)

                this%grad(i,j,k) = this%grad(i,j,k) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%pl_z(i,j,k) * this%dA(i,j,k) * &
                (f_vec%pl_z(i,j,k) + f_vec%mi_z(i,j,k)), kind=dp)

                this%grad(i,j,k) = this%grad(i,j,k) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%mi_x(i,j,k) * this%dA(i,j,k) * &
                (f_vec%pl_x(i,j,k) + f_vec%mi_x(i,j,k)), kind=dp)

                this%grad(i,j,k) = this%grad(i,j,k) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%mi_y(i,j,k) * this%dA(i,j,k) * &
                (f_vec%pl_y(i,j,k) + f_vec%mi_y(i,j,k)), kind=dp)

                this%grad(i,j,k) = this%grad(i,j,k) + real( &
                1.0d0/this%w_dL(n) * f_adj_vec(n)%mi_z(i,j,k) * this%dA(i,j,k) * &
                (f_vec%pl_z(i,j,k) + f_vec%mi_z(i,j,k)), kind=dp)

            end do
        end do
        end do
        end do

    end select

    end select

end subroutine compute_gradient

!###################################################################################################

subroutine update_permittivity(this, rho_conv)

    class(TOptPrblm) , intent(inout) :: this
    real(dp)         , intent(in)    :: rho_conv(:,:,:)

    integer :: i, j, k
    real(dp) :: beta_p
    real(dp) :: eta_p
    real(dp) :: eta, eta_max, eta_min
    real(dp) :: kappa, kappa_max
    real(dp) :: eps_Re, eps_Im
    real(dp) :: C1, C2
    real(dp) :: func_rho

    beta_p    = this%beta_rho
    eta_p     = this%eta_rho
    eta_max   = this%eta_max
    eta_min   = this%eta_min
    kappa_max = this%kappa_max
    eps_Re    = this%eps_Re
    eps_Im    = this%eps_Im

    C1 = 1.0d0 / (DTANH(beta_p*eta_p) + DTANH(beta_p*(1.0-eta_p)))
    C2 = DTANH(beta_p*eta_p)*C1

    select case (this%dimensions)
    case (1)
        if (this%is_complex) then
            do k = 1, this%nz
                func_rho  = C1*DTANH(beta_p*(rho_conv(k,1,1)-eta_p)) + C2
                eta       = (eta_max - eta_min)*func_rho + eta_min
                kappa     = kappa_max*func_rho
                this%eps_r%mat1D(k) = Z_ONE*eps0*(eta**2 - kappa**2) - Z_I*2.0d0*eps0*eta*kappa
            end do
        else
            do k = 1, this%nz
                func_rho  = C1*DTANH(beta_p*(rho_conv(k,1,1)-eta_p)) + C2
                this%eps_r%mat1D(k) = Z_ONE*(func_rho*(eps_Re*eps0 - eps0) + eps0)
            end do
        end if
    case (2)
        if (this%is_complex) then
            do j = 1, this%ny
            do i = 1, this%nx
                func_rho  = C1*DTANH(beta_p*(rho_conv(i,j,1)-eta_p)) + C2
                eta       = (eta_max - eta_min)*func_rho + eta_min
                kappa     = kappa_max*func_rho
                this%eps_r%mat2D(i,j) = Z_ONE*eps0*(eta**2 - kappa**2) - Z_I*2.0d0*eps0*eta*kappa
            end do
            end do
        else
            do j = 1, this%ny
            do i = 1, this%nx
                func_rho  = C1*DTANH(beta_p*(rho_conv(i,j,1)-eta_p)) + C2
                this%eps_r%mat2D(i,j) = Z_ONE*(func_rho*(eps_Re*eps0 - eps0) + eps0)
            end do
            end do
        end if
    case (3)

        if (this%is_complex) then
            do k = 1, this%nz
            do j = 1, this%ny
            do i = 1, this%nx
                func_rho  = C1*DTANH(beta_p*(rho_conv(i,j,k)-eta_p)) + C2
                eta       = (eta_max - eta_min)*func_rho + eta_min
                kappa     = kappa_max*func_rho
                this%eps_r%mat3D(i,j,k) = Z_ONE*eps0*(eta**2 - kappa**2) - Z_I*2.0d0*eps0*eta*kappa
            end do
            end do
            end do
        else
            do k = 1, this%nz
            do j = 1, this%ny
            do i = 1, this%nx
                func_rho  = C1*DTANH(beta_p*(rho_conv(i,j,k)-eta_p)) + C2
                this%eps_r%mat3D(i,j,k) = Z_ONE*(func_rho*(eps_Re*eps0 - eps0) + eps0)
            end do
            end do
            end do
        end if

    end select

end subroutine update_permittivity

!###################################################################################################

subroutine update_fields(this)

    class(TOptPrblm) , intent(inout) :: this

    integer :: n

    select type (f_vec => this%f_vec)
    
    type is (TRSvec_1D)
    select type (f_vec_new => this%f_vec_new)
    type is (TRSvec_1D)
    select type ( f_adj_vec => this%f_adj_vec)
    type is (TRSvec_1D)
    select type (f_adj_vec_new => this%f_adj_vec_new)
    type is (TRSvec_1D)
        f_vec%pl_x = f_vec_new%pl_x
        f_vec%mi_x = f_vec_new%mi_x
        f_vec%pl_y = f_vec_new%pl_y
        f_vec%mi_y = f_vec_new%mi_y

        do n = 1, this%n_trg
            f_adj_vec(n)%pl_x = f_adj_vec_new(n)%pl_x
            f_adj_vec(n)%mi_x = f_adj_vec_new(n)%mi_x
            f_adj_vec(n)%pl_y = f_adj_vec_new(n)%pl_y
            f_adj_vec(n)%mi_y = f_adj_vec_new(n)%mi_y
        end do        

    end select
    end select
    end select

    type is (TRSvec_2D)
    select type (f_vec_new => this%f_vec_new)
    type is (TRSvec_2D)
    select type (f_adj_vec => this%f_adj_vec)
    type is (TRSvec_2D)
    select type (f_adj_vec_new => this%f_adj_vec_new)
    type is (TRSvec_2D)

        f_vec%pl_x = f_vec_new%pl_x
        f_vec%mi_x = f_vec_new%mi_x
        f_vec%pl_y = f_vec_new%pl_y
        f_vec%mi_y = f_vec_new%mi_y
        f_vec%pl_z = f_vec_new%pl_z
        f_vec%mi_z = f_vec_new%mi_z
    
        do n = 1, this%n_trg
            f_adj_vec(n)%pl_x = f_adj_vec_new(n)%pl_x
            f_adj_vec(n)%mi_x = f_adj_vec_new(n)%mi_x
            f_adj_vec(n)%pl_y = f_adj_vec_new(n)%pl_y
            f_adj_vec(n)%mi_y = f_adj_vec_new(n)%mi_y
            f_adj_vec(n)%pl_z = f_adj_vec_new(n)%pl_z
            f_adj_vec(n)%mi_z = f_adj_vec_new(n)%mi_z
        end do

    end select
    end select
    end select

    type is (TRSvec_3D)
    select type (f_vec_new => this%f_vec_new)
    type is (TRSvec_3D)
    select type (f_adj_vec => this%f_adj_vec)
    type is (TRSvec_3D)
    select type (f_adj_vec_new => this%f_adj_vec_new)
    type is (TRSvec_3D)

        f_vec%pl_x = f_vec_new%pl_x
        f_vec%mi_x = f_vec_new%mi_x
        f_vec%pl_y = f_vec_new%pl_y
        f_vec%mi_y = f_vec_new%mi_y
        f_vec%pl_z = f_vec_new%pl_z
        f_vec%mi_z = f_vec_new%mi_z

        do n = 1, this%n_trg
            f_adj_vec(n)%pl_x = f_adj_vec_new(n)%pl_x
            f_adj_vec(n)%mi_x = f_adj_vec_new(n)%mi_x
            f_adj_vec(n)%pl_y = f_adj_vec_new(n)%pl_y
            f_adj_vec(n)%mi_y = f_adj_vec_new(n)%mi_y
            f_adj_vec(n)%pl_z = f_adj_vec_new(n)%pl_z
            f_adj_vec(n)%mi_z = f_adj_vec_new(n)%mi_z
        end do

    end select
    end select
    end select

    end select

end subroutine update_fields

!###################################################################################################

subroutine update_targets(this)

    class(TOptPrblm) , intent(inout) :: this

    integer :: n

    this%w_total = 0.0d0

    do n = 1, this%n_src_trg
        call update_target(this%src_trg(n), this%f_vec_new, this%w_dL(n), this%w_total)
    end do

end subroutine update_targets

!###################################################################################################

subroutine read_opt_region(prblm, dr, grid_Ndims, mpi_coords, mpi_dims)

    type(TOptPrblm)  , intent(inout) :: prblm
    real(dp)         , intent(in)    :: dr
    integer          , intent(in)    :: grid_Ndims(3)
    integer          , intent(in)    :: mpi_coords(3)
    integer          , intent(in)    :: mpi_dims(3)

    character(len=20) :: file_name = "opt_region_"
    character(len=20) :: file_exten = ".in"
    character(len=20) :: file_number
    character(len=20) :: input_name
    integer           :: ierr, funit
    integer           :: i_ndx, j_ndx, k_ndx
    integer           :: i_min, i_max, j_min, j_max, k_min, k_max
    real(dp)          :: x, y, z


    i_min = mpi_coords(1)*prblm%nx + 1
    i_max = (mpi_coords(1)+1)*prblm%nx
    j_min = mpi_coords(2)*prblm%ny + 1
    j_max = (mpi_coords(2)+1)*prblm%ny
    k_min = mpi_coords(3)*prblm%nz + 1
    k_max = (mpi_coords(3)+1)*prblm%nz


    write(file_number, '(I3.3)') prblm%id

    input_name = trim(file_name) // trim(file_number) // trim(file_exten)

    inquire(file=trim(input_name), iostat=ierr)
    if (ierr /= 0) then
        write(*,*) 'Error: Optimization problem', prblm%id, 'does not count'
        write(*,*) 'with an optimization region file:'
        write(*,*) trim(input_name)
        error stop
    end if

    open(newunit=funit, file=trim(input_name), status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        write(*,*) 'Error opening optimization region file:'
        write(*,*) trim(input_name)
        error stop
    end if

    select case (prblm%dimensions)
    case (1)

        do
            read(funit, *, iostat=ierr) z
            if (ierr /= 0) exit

            i_ndx = int(z/dr) + int(grid_Ndims(1)/2)
            
            prblm%opt_region(i_ndx,1,1) = .true.

        end do

    case (2)
    
        do
            read(funit, *, iostat=ierr) x, y
            if (ierr /= 0) exit

            i_ndx = int(x/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
            j_ndx = int(y/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)

            if ((i_ndx >= i_min .and. i_ndx <= i_max) .and. &
                (j_ndx >= j_min .and. j_ndx <= j_max)) then
                prblm%opt_region(i_ndx,j_ndx,1) = .true.
            end if

        end do 

    case (3)
   
        do
            read(funit, *, iostat=ierr) x, y, z
            if (ierr /= 0) exit

            i_ndx = int(x/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
            j_ndx = int(y/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
            k_ndx = int(z/dr) + int(grid_Ndims(3)*mpi_dims(3)/2)

            if ((i_ndx >= i_min .and. i_ndx <= i_max) .and. &
                (j_ndx >= j_min .and. j_ndx <= j_max) .and. &
                (k_ndx >= k_min .and. k_ndx <= k_max)) then
                prblm%opt_region(i_ndx,j_ndx,k_ndx) = .true.
            end if
        end do

    end select

end subroutine read_opt_region

!###################################################################################################

subroutine clean_opt_region(prblm)

    type(TOptPrblm)    , intent(inout) :: prblm

    integer :: i, j, k

    select case (prblm%dimensions)
    case (1)
        do k = 1, prblm%nz
            if (prblm%eps_r%mat1D(k) /= Z_ONE) then
                prblm%opt_region(k,1,1) = .false.
            end if
        end do
    case (2)
        do j = 1, prblm%ny
        do i = 1, prblm%nx
            if (prblm%eps_r%mat2D(i,j) /= Z_ONE) then
                prblm%opt_region(i,j,1) = .false.
            end if
        end do
        end do
    case (3)
        do k = 1, prblm%nz
        do j = 1, prblm%ny
        do i = 1, prblm%nx
            if (prblm%eps_r%mat3D(i,j,k) /= Z_ONE) then
                prblm%opt_region(i,j,k) = .false.
            end if
        end do
        end do
        end do
    end select

end subroutine clean_opt_region

!###################################################################################################

end module optimization_problem_mod