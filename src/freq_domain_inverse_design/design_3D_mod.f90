module design_3D_mod

#ifdef USE_MPI
    use mpi
#endif

    use design_base_mod

    implicit none

    type, extends(TDesign) :: TDesign_3D

        real(dp) , allocatable :: ker_mat(:,:,:)
        real(dp) , allocatable :: rho(:,:,:)
        real(dp) , allocatable :: rho_conv(:,:,:)
        real(dp) , allocatable :: rho_old(:,:,:)
        real(dp) , allocatable :: grad(:,:,:)
        real(dp) , allocatable :: grad_old(:,:,:)
        real(dp) , allocatable :: grad_conv(:,:,:)
        real(dp) , allocatable :: grad_new(:,:,:)
        logical  , allocatable :: opt_region(:,:,:)

        contains
            procedure :: init_design             => init_3D_design
            procedure :: kill_design             => kill_3D_design
            procedure :: collect_opt_regions     => collect_3D_opt_regions
            procedure :: collect_FOM             => collect_3D_FOM
            procedure :: collect_gradients       => collect_3D_gradients
            procedure :: apply_kernel_on_rho     => apply_kernel_on_rho_3D
            procedure :: apply_kernel_on_grad    => apply_kernel_on_grad_3D
            procedure :: displace_rho            => displace_rho_3D
            procedure :: reset_rho_one_step_back => reset_rho_one_step_back_3D
            procedure :: reset_grad              => reset_grad_3D
            procedure :: update_rho              => update_rho_3D
            procedure :: update_grad             => update_grad_3D

    end type TDesign_3D

contains
!###################################################################################################

subroutine init_3D_design(this, dimensions, dr, sigma, grid_Ndims)

    class(TDesign_3D), intent(inout) :: this
    integer          , intent(in)     :: dimensions
    real(dp)         , intent(in)     :: dr
    real(dp)         , intent(in)     :: sigma
    integer          , intent(in)     :: grid_Ndims(3)

    integer       :: i, j, k
    integer       :: n_ker, nx, ny, nz
    real(dp)      :: ker_sum

    this%dimensions = dimensions
    this%nx         = grid_Ndims(1)
    this%ny         = grid_Ndims(2)
    this%nz         = grid_Ndims(3)
    this%n_ker      = int(3.0_dp*sigma/dr)

    nx    = this%nx
    ny    = this%ny
    nz    = this%nz
    n_ker = this%n_ker

    if (.not. allocated(this%rho))       allocate(this%rho(-n_ker+1:nx+n_ker, -n_ker+1:ny+n_ker, &
                                                  -n_ker+1:nz+n_ker))
    if (.not. allocated(this%rho_conv))  allocate(this%rho_conv(nx, ny, nz))
    if (.not. allocated(this%rho_old))   allocate(this%rho_old(nx, ny, nz))
    if (.not. allocated(this%grad))      allocate(this%grad(-n_ker+1:nx+n_ker, -n_ker+1:ny+n_ker, &
                                                  -n_ker+1:nz+n_ker))
    if (.not. allocated(this%grad_conv)) allocate(this%grad_conv(nx, ny, nz))
    if (.not. allocated(this%grad_old))  allocate(this%grad_old(nx, ny, nz))
    if (.not. allocated(this%grad_new))  allocate(this%grad_new(nx, ny, nz))
    if (.not. allocated(this%opt_region)) &
        allocate(this%opt_region(-n_ker+1:nx+n_ker, -n_ker+1:ny+n_ker, -n_ker+1:nz+n_ker))
    if (.not. allocated(this%ker_mat))   allocate(this%ker_mat(-n_ker:n_ker, -n_ker:n_ker, &
                                                  -n_ker:n_ker))

    ker_sum = 0.0_dp
    do k = -n_ker, n_ker
    do j = -n_ker, n_ker
    do i = -n_ker, n_ker
        this%ker_mat(i, j, k) = EXP(-REAL(i**2+j**2+k**2, kind=dp)/(2.0_dp*sigma**2))
        ker_sum = ker_sum + this%ker_mat(i, j, k)
    end do
    end do
    end do

    this%ker_mat = this%ker_mat / ker_sum

    this%rho        = 0.0_dp
    this%rho_conv   = 0.0_dp
    this%rho_old    = 0.0_dp
    this%opt_region = .false.
    this%grad_max   = 0.0_dp

end subroutine init_3D_design

!###################################################################################################

subroutine kill_3D_design(this)

    class(TDesign_3D), intent(inout) :: this

    if (allocated(this%rho))       deallocate(this%rho)
    if (allocated(this%rho_conv))  deallocate(this%rho_conv)
    if (allocated(this%rho_old))   deallocate(this%rho_old)
    if (allocated(this%grad))      deallocate(this%grad)
    if (allocated(this%grad_conv)) deallocate(this%grad_conv)
    if (allocated(this%grad_old))  deallocate(this%grad_old)
    if (allocated(this%grad_new))  deallocate(this%grad_new)
    if (allocated(this%opt_region)) deallocate(this%opt_region)
    if (allocated(this%ker_mat))   deallocate(this%ker_mat)

end subroutine kill_3D_design

!###################################################################################################

subroutine collect_3D_opt_regions(this, opt_region_i, rho_init)

    class(TDesign_3D), intent(inout) :: this
    logical          , intent(in)    :: opt_region_i(:,:,:)
    real(dp)         , intent(in)    :: rho_init

    integer :: i, j, k

    do k = 1, this%nz
    do j = 1, this%ny
    do i = 1, this%nx
        if (opt_region_i(i, j, k)) then
            this%opt_region(i, j, k) = .true.
            this%rho(i, j, k) = rho_init
        else
            this%opt_region(i, j, k) = .false.
            this%rho(i, j, k) = 0.0_dp
        end if
    end do
    end do
    end do

end subroutine collect_3D_opt_regions

!###################################################################################################

subroutine collect_3D_FOM(this, w_p, p, n_opt_problems)

    class(TDesign_3D), intent(inout) :: this
    real(dp)         , intent(in)    :: w_p
    integer          , intent(in)    :: p
    integer          , intent(in)    :: n_opt_problems

    integer  :: ierr
    real(dp) :: fom_loc = 0.0_dp
    real(dp) :: fom_sum = 0.0_dp

    if (p == 1) this%fom = 0.0_dp

    this%fom = this%fom + w_p

    !MPI is already used to share w_p across ranks.

end subroutine collect_3D_FOM

!###################################################################################################

subroutine collect_3D_gradients(this, grad_in, p, n_opt_problems)

    class(TDesign_3D), intent(inout) :: this
    real(dp)         , intent(in)    :: grad_in(:,:,:)
    integer          , intent(in)    :: p
    integer          , intent(in)    :: n_opt_problems

    if  (p == 1) this%grad = 0.0_dp

    this%grad(1:this%nx,1:this%ny,1:this%nz) = this%grad(1:this%nx,1:this%ny,1:this%nz) + &
                                           grad_in(1:this%nx,1:this%ny,1:this%nz)

    if (p == n_opt_problems) then
        this%grad = 2.0_dp * this%fom * this%grad
    end if

end subroutine collect_3D_gradients

!###################################################################################################

subroutine apply_kernel_on_rho_3D(this)

    class(TDesign_3D), intent(inout) :: this

    integer :: i, j, k, ii, jj, kk

    this%rho_conv = 0.0_dp

    do k = 1, this%nz
    do j = 1, this%ny
    do i = 1, this%nx
        if (this%opt_region(i,j,k)) then
            do kk = -this%n_ker, this%n_ker
            do jj = -this%n_ker, this%n_ker
            do ii = -this%n_ker, this%n_ker
                if (this%opt_region(i+ii,j+jj,k+kk)) then
                    this%rho_conv(i,j,k) = this%rho_conv(i,j,k) + &
                        this%ker_mat(ii,jj,kk)*this%rho(i+ii,j+jj,k+kk)
                end if
            end do
            end do
            end do
        end if
    end do
    end do
    end do

end subroutine apply_kernel_on_rho_3D

!###################################################################################################

subroutine apply_kernel_on_grad_3D(this)

    class(TDesign_3D), intent(inout) :: this

    integer :: i, j, k, ii, jj, kk

    this%grad_conv = 0.0_dp

    do k = 1, this%nz
    do j = 1, this%ny
    do i = 1, this%nx
        if (this%opt_region(i,j,k)) then
            do kk = -this%n_ker, this%n_ker
            do jj = -this%n_ker, this%n_ker
            do ii = -this%n_ker, this%n_ker
                if (this%opt_region(i+ii,j+jj,k+kk)) then
                    this%grad_conv(i,j,k) = this%grad_conv(i,j,k) + &
                        this%ker_mat(ii,jj,kk)*this%grad(i+ii,j+jj,k+kk)
                end if
            end do
            end do
            end do
        end if
    end do
    end do
    end do

end subroutine apply_kernel_on_grad_3D

!###################################################################################################

subroutine displace_rho_3D(this)

    class(TDesign_3D), intent(inout) :: this

    integer  :: ierr
    real(dp) :: gamma
    real(dp) :: deno_loc
    real(dp) :: nume_loc
    real(dp) :: deno
    real(dp) :: nume
    real(dp) :: norm_loc
    real(dp) :: norm_global

    deno_loc  = SUM(this%grad_old*this%grad_old)
    nume_loc = SUM((this%grad_conv - this%grad_old)*this%grad_conv)

#ifdef USE_MPI
    call MPI_Allreduce(deno_loc, deno, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(nume_loc, nume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
    deno = deno_loc
    nume = nume_loc
#endif

    if (deno > 0.0_dp) then
        gamma = nume/deno
    else
        gamma = 0.0_dp
    end if

   this%grad_new = this%grad_conv + gamma*this%grad_old

    norm_loc = MAXVAL(ABS(this%grad_new))

#ifdef USE_MPI
    call MPI_Allreduce(norm_loc, norm_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                       MPI_COMM_WORLD, ierr)
#else
    norm_global = norm_loc
#endif

    this%grad_max = norm_global
    
    this%rho(1:this%nx,1:this%ny,1:this%nz) = this%rho_old(1:this%nx,1:this%ny,1:this%nz) + &
            this%grad_new(1:this%nx,1:this%ny,1:this%nz) * this%drho / norm_global

end subroutine displace_rho_3D

!###################################################################################################

subroutine reset_rho_one_step_back_3D(this)

    class(TDesign_3D), intent(inout) :: this

    this%rho(1:this%nx,1:this%ny,1:this%nz) = this%rho_old(1:this%nx,1:this%ny,1:this%nz)

end subroutine reset_rho_one_step_back_3D

!###################################################################################################

subroutine reset_grad_3D(this)

    class(TDesign_3D), intent(inout) :: this

    this%grad(1:this%nx,1:this%ny,1:this%nz) = 0.0_dp

end subroutine reset_grad_3D

!###################################################################################################

subroutine update_rho_3D(this)

    class(TDesign_3D), intent(inout) :: this

    this%rho_old(1:this%nx,1:this%ny,1:this%nz) = this%rho(1:this%nx,1:this%ny,1:this%nz)

end subroutine update_rho_3D

!###################################################################################################

subroutine update_grad_3D(this)

    class(TDesign_3D), intent(inout) :: this

    this%grad_old(1:this%nx,1:this%ny,1:this%nz) = this%grad_conv(1:this%nx,1:this%ny,1:this%nz)

end subroutine update_grad_3D

!###################################################################################################

end module design_3D_mod