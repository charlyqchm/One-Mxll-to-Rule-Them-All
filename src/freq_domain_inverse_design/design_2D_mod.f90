module design_2D_mod

#ifdef USE_MPI
    use mpi
#endif

    use design_base_mod

    implicit none

    type, extends(TDesign) :: TDesign_2D

        real(dp) , allocatable :: ker_mat(:, :)
        real(dp) , allocatable :: rho(:, :)
        real(dp) , allocatable :: rho_conv(:, :)
        real(dp) , allocatable :: rho_old(:, :)
        real(dp) , allocatable :: grad(:, :)
        real(dp) , allocatable :: grad_old(:, :)
        real(dp) , allocatable :: grad_conv(:, :)
        real(dp) , allocatable :: grad_new(:, :)
        logical  , allocatable :: opt_region(:, :)

        contains
            procedure :: init_design             => init_2D_design
            procedure :: kill_design             => kill_2D_design
            procedure :: collect_opt_regions     => collect_2D_opt_regions
            procedure :: collect_FOM             => collect_2D_FOM
            procedure :: collect_gradients       => collect_2D_gradients
            procedure :: apply_kernel_on_rho     => apply_kernel_on_rho_2D
            procedure :: apply_kernel_on_grad    => apply_kernel_on_grad_2D
            procedure :: displace_rho            => displace_rho_2D
            procedure :: reset_rho_one_step_back => reset_rho_one_step_back_2D
            procedure :: reset_grad              => reset_grad_2D
            procedure :: update_rho              => update_rho_2D
            procedure :: update_grad             => update_grad_2D

    end type TDesign_2D

contains
!###################################################################################################

subroutine init_2D_design(this, dimensions, dr, sigma, grid_Ndims)

    class(TDesign_2D), intent(inout) :: this
    integer          , intent(in)     :: dimensions
    real(dp)         , intent(in)     :: dr
    real(dp)         , intent(in)     :: sigma
    integer          , intent(in)     :: grid_Ndims(3)

    integer       :: i, j
    integer       :: n_ker, nx, ny
    real(dp)      :: ker_sum

    this%dimensions = dimensions
    this%nx         = grid_Ndims(1)
    this%ny         = grid_Ndims(2)
    this%n_ker      = int(3.0_dp*sigma/dr)

    nx    = this%nx
    ny    = this%ny
    n_ker = this%n_ker

    if (.not. allocated(this%rho))       allocate(this%rho(-n_ker+1:nx+n_ker, -n_ker+1:ny+n_ker))
    if (.not. allocated(this%rho_conv))  allocate(this%rho_conv(nx, ny))
    if (.not. allocated(this%rho_old))   allocate(this%rho_old(nx, ny))
    if (.not. allocated(this%grad))      allocate(this%grad(-n_ker+1:nx+n_ker, -n_ker+1:ny+n_ker))
    if (.not. allocated(this%grad_conv)) allocate(this%grad_conv(nx, ny))
    if (.not. allocated(this%grad_old))  allocate(this%grad_old(nx, ny))
    if (.not. allocated(this%grad_new))  allocate(this%grad_new(nx, ny))
    if (.not. allocated(this%opt_region)) allocate(this%opt_region(-n_ker+1:nx+n_ker, &
                                                                   -n_ker+1:ny+n_ker))
    if (.not. allocated(this%ker_mat))   allocate(this%ker_mat(-n_ker:n_ker, -n_ker:n_ker))

    ker_sum = 0.0_dp
    do j = -n_ker, n_ker
        do i = -n_ker, n_ker
            this%ker_mat(i, j) = EXP(-REAL(i**2+j**2, kind=dp)/(2.0_dp*sigma**2))
            ker_sum = ker_sum + this%ker_mat(i, j)
        end do
    end do

    this%ker_mat = this%ker_mat / ker_sum

    this%rho        = 0.0_dp
    this%rho_conv   = 0.0_dp
    this%rho_old    = 0.0_dp
    this%opt_region = .false.
    this%grad_max   = 0.0_dp

end subroutine init_2D_design

!###################################################################################################

subroutine kill_2D_design(this)

    class(TDesign_2D), intent(inout) :: this

    if (allocated(this%rho))       deallocate(this%rho)
    if (allocated(this%rho_conv))  deallocate(this%rho_conv)
    if (allocated(this%rho_old))   deallocate(this%rho_old)
    if (allocated(this%grad))      deallocate(this%grad)
    if (allocated(this%grad_conv)) deallocate(this%grad_conv)
    if (allocated(this%grad_old))  deallocate(this%grad_old)
    if (allocated(this%grad_new))  deallocate(this%grad_new)
    if (allocated(this%opt_region)) deallocate(this%opt_region)
    if (allocated(this%ker_mat))   deallocate(this%ker_mat)

end subroutine kill_2D_design

!###################################################################################################

subroutine collect_2D_opt_regions(this, opt_region_i, rho_init)

    class(TDesign_2D), intent(inout) :: this
    logical          , intent(in)    :: opt_region_i(:,:,:)
    real(dp)         , intent(in)    :: rho_init

    integer :: i, j

    do j = 1, this%ny
    do i = 1, this%nx
        if (opt_region_i(i, j, 1)) then
            this%opt_region(i, j) = .true.
            this%rho(i, j) = rho_init
        else
            this%opt_region(i, j) = .false.
            this%rho(i, j) = 0.0_dp
        end if
    end do
    end do

end subroutine collect_2D_opt_regions

!###################################################################################################

subroutine collect_2D_FOM(this, w_p, p, n_opt_problems)

    class(TDesign_2D), intent(inout) :: this
    real(dp)         , intent(in)    :: w_p
    integer          , intent(in)    :: p
    integer          , intent(in)    :: n_opt_problems

    integer  :: ierr
    real(dp) :: fom_loc = 0.0_dp
    real(dp) :: fom_sum = 0.0_dp

    if (p == 1) this%fom = 0.0_dp

    this%fom = this%fom + w_p

    !MPI is already used to share w_p across ranks.

end subroutine collect_2D_FOM

!###################################################################################################

subroutine collect_2D_gradients(this, grad_in, p, n_opt_problems)

    class(TDesign_2D), intent(inout) :: this
    real(dp)         , intent(in)    :: grad_in(:,:,:)
    integer          , intent(in)    :: p
    integer          , intent(in)    :: n_opt_problems

    if  (p == 1) this%grad = 0.0_dp

    this%grad(1:this%nx,1:this%ny) = this%grad(1:this%nx,1:this%ny) + &
                                           grad_in(1:this%nx,1:this%ny,1)

    if (p == n_opt_problems) then
        this%grad = 2.0_dp * this%fom * this%grad
    end if

end subroutine collect_2D_gradients

!###################################################################################################

subroutine apply_kernel_on_rho_2D(this)

    class(TDesign_2D), intent(inout) :: this

    integer :: i, j, ii, jj

    this%rho_conv = 0.0_dp

    do j = 1, this%ny
    do i = 1, this%nx
        if (this%opt_region(i,j)) then
            do ii = -this%n_ker, this%n_ker
            do jj = -this%n_ker, this%n_ker
                if (this%opt_region(i+ii,j+jj)) then
                    this%rho_conv(i,j) = this%rho_conv(i,j) + &
                        this%ker_mat(ii,jj)*this%rho(i+ii,j+jj)
                end if
            end do
            end do
        end if
    end do
    end do

end subroutine apply_kernel_on_rho_2D

!###################################################################################################

subroutine apply_kernel_on_grad_2D(this)

    class(TDesign_2D), intent(inout) :: this

    integer :: i, j, ii, jj

    this%grad_conv = 0.0_dp

    do i = 1, this%nx
    do j = 1, this%ny
        if (this%opt_region(i,j)) then
            do ii = -this%n_ker, this%n_ker
            do jj = -this%n_ker, this%n_ker
                if (this%opt_region(i+ii,j+jj)) then
                    this%grad_conv(i,j) = this%grad_conv(i,j) + &
                        this%ker_mat(ii,jj)*this%grad(i+ii,j+jj)
                end if
            end do
            end do
        end if
    end do
    end do

end subroutine apply_kernel_on_grad_2D

!###################################################################################################

subroutine displace_rho_2D(this)

    class(TDesign_2D), intent(inout) :: this

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
    
    this%rho(1:this%nx,1:this%ny) = this%rho_old(1:this%nx,1:this%ny) + &
            this%grad_new(1:this%nx,1:this%ny) * this%drho / norm_global

end subroutine displace_rho_2D

!###################################################################################################

subroutine reset_rho_one_step_back_2D(this)

    class(TDesign_2D), intent(inout) :: this

    this%rho(1:this%nx,1:this%ny) = this%rho_old(1:this%nx,1:this%ny)

end subroutine reset_rho_one_step_back_2D

!###################################################################################################

subroutine reset_grad_2D(this)

    class(TDesign_2D), intent(inout) :: this

    this%grad(1:this%nx,1:this%ny) = 0.0_dp

end subroutine reset_grad_2D

!###################################################################################################

subroutine update_rho_2D(this)

    class(TDesign_2D), intent(inout) :: this

    this%rho_old(1:this%nx,1:this%ny) = this%rho(1:this%nx,1:this%ny)

end subroutine update_rho_2D

!###################################################################################################

subroutine update_grad_2D(this)

    class(TDesign_2D), intent(inout) :: this

    this%grad_old(1:this%nx,1:this%ny) = this%grad_conv(1:this%nx,1:this%ny)

end subroutine update_grad_2D

!###################################################################################################

end module design_2D_mod