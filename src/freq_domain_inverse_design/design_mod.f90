module design_mod

    use constants_mod
    use medium_mod
    use allocator_multidim_mod
    use parallel_subs_mod

    implicit none

    type TDesign

        integer  :: dimensions
        integer  :: nx
        integer  :: ny
        integer  :: nz
        integer  :: n_ker
        real(dp) :: sigma
        real(dp) :: fom
        real(dp) :: drho

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
            procedure :: init_design
            procedure :: kill_design
            procedure :: collect_opt_regions
            procedure :: collect_FOM
            procedure :: collect_gradients
            procedure :: apply_kernel_on_rho
            procedure :: apply_kernel_on_grad
            procedure :: displace_rho
            procedure :: reset_rho_one_step_back
            procedure :: reset_grad
            procedure :: update_rho
            procedure :: update_grad

    end type TDesign

contains
!###################################################################################################

subroutine init_design(this, dimensions, dr, sigma, grid_Ndims)

    class(TDesign) :: this
    integer        :: dimensions
    real(dp)       :: dr
    real(dp)       :: sigma
    integer        :: grid_Ndims(3)

    this%dimensions = dimensions
    this%nx = grid_Ndims(1)
    this%ny = grid_Ndims(2)
    this%nz = grid_Ndims(3)
    this%n_ker = int(3.0_dp*sigma/dr)

    call allocate_multidim(array=this%rho, dim=dimensions, i_max=this%nx+this%n_ker, &
                           i_min=-this%n_ker+1, j_max=this%ny+this%n_ker, j_min=-this%n_ker+1, &
                           k_max=this%nz+this%n_ker, k_min=-this%n_ker+1)

    call allocate_multidim(array=this%rho_conv, dim=dimensions, i_max=this%nx, &
                           j_max=this%ny, k_max=this%nz)

    call allocate_multidim(array=this%rho_old, dim=dimensions, i_max=this%nx, &
                            j_max=this%ny, k_max=this%nz)

    call allocate_multidim(array=this%grad, dim=dimensions, i_max=this%nx+this%n_ker, &
                           i_min=-this%n_ker+1, j_max=this%ny+this%n_ker, j_min=-this%n_ker+1, &
                           k_max=this%nz+this%n_ker, k_min=-this%n_ker+1)

    call allocate_multidim(array=this%grad_conv, dim=dimensions, i_max=this%nx, &
                           j_max=this%ny, k_max=this%nz)

    call allocate_multidim(array=this%grad_new, dim=dimensions, i_max=this%nx, &
                           j_max=this%ny, k_max=this%nz)

    call allocate_multidim(array=this%grad_old, dim=dimensions, i_max=this%nx, &
                            j_max=this%ny, k_max=this%nz)

    call allocate_multidim(array=this%opt_region, dim=dimensions, i_max=this%nx+this%n_ker, &
                           i_min=-this%n_ker+1, j_max=this%ny+this%n_ker, j_min=-this%n_ker+1, &
                           k_max=this%nz+this%n_ker, k_min=-this%n_ker+1)

    call allocate_multidim(array=this%ker_mat, dim=dimensions, i_max=this%n_ker, i_min=-this%n_ker, &
                           j_max=this%n_ker, j_min=-this%n_ker, k_max=this%n_ker, k_min=-this%n_ker)

    this%rho        = R_0
    this%rho_conv   = R_0
    this%rho_old    = R_0
    this%opt_region = .false.

end subroutine init_design

!###################################################################################################

subroutine kill_design(this)
    class(TDesign) :: this

    if (allocated(this%ker_mat))     deallocate(this%ker_mat)
    if (allocated(this%rho))         deallocate(this%rho)
    if (allocated(this%rho_conv))    deallocate(this%rho_conv)
    if (allocated(this%rho_old))     deallocate(this%rho_old)
    if (allocated(this%grad))        deallocate(this%grad)
    if (allocated(this%grad_conv))   deallocate(this%grad_conv)
    if (allocated(this%grad_new))    deallocate(this%grad_new)
    if (allocated(this%grad_old))    deallocate(this%grad_old)
    if (allocated(this%opt_region))  deallocate(this%opt_region)

end subroutine kill_design

!###################################################################################################

subroutine collect_opt_regions(this, opt_region_i, rho_init, p_id, p_tot)

    class(TDesign) :: this
    logical        :: opt_region_i(:,:,:)
    real(dp)       :: rho_init
    integer        :: p_id
    integer        :: p_tot

    integer :: i, j, k

    ! this%opt_region is .false. from init_design.

    select case (this%dimensions)
    case (1)
        this%opt_region(1:this%nx,1,1) = this%opt_region(1:this%nx,1,1) .or. opt_region_i
        do i = 1, this%nx
            if (this%opt_region(i,1,1)) this%rho(i,1,1) = rho_init
        end do
    case (2)
        this%opt_region(1:this%nx, 1:this%ny, 1) = this%opt_region(1:this%nx, 1:this%ny, 1) .or. &
                                                   opt_region_i
        do j = 1, this%ny
        do i = 1, this%nx
            if (this%opt_region(i,j,1)) this%rho(i,j,1) = rho_init
        end do
        end do
    case (3)
        this%opt_region(1:this%nx, 1:this%ny, 1:this%nz) = &
            this%opt_region(1:this%nx, 1:this%ny, 1:this%nz) .or. opt_region_i
    
        do k = 1, this%nz
        do j = 1, this%ny
        do i = 1, this%nx
            if (this%opt_region(i,j,k)) this%rho(i,j,k) = rho_init
        end do
        end do
        end do
    end select 

    if (p_id == p_tot) then
        call extend_array_to_ranks(this%rho, this%dimensions, this%nx, this%ny, this%nz, this%n_ker)
        call extend_array_to_ranks(this%opt_region, this%dimensions, this%nx, this%ny, this%nz, &
                                   this%n_ker)
    end if



end subroutine collect_opt_regions

!###################################################################################################

subroutine collect_FOM(this, w_p, p, n_opt_problems)

    class(TDesign)  :: this
    real(dp)        :: w_p
    integer         :: p
    integer         :: n_opt_problems

    real(dp) :: fom_loc = R_0
    real(dp) :: fom_sum = R_0
#ifdef USE_MPI
    integer :: ierr
#endif

    if (p == 1) this%fom = R_0

    this%fom = this%fom + w_p

#ifdef USE_MPI

    if (p == n_opt_problems) then
        fom_loc = this%fom
        call MPI_Allreduce(fom_loc, fom_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        this%fom = fom_sum
    end if

#endif

end subroutine collect_FOM

!###################################################################################################

subroutine collect_gradients(this, grad_in, p, n_opt_problems)

    class(TDesign)  :: this
    real(dp)        :: grad_in(:,:,:)
    integer         :: p
    integer         :: n_opt_problems

    if (p == 1) this%grad = R_0

    this%grad = this%grad + grad_in

    if (p == n_opt_problems) then
        this%grad = 2.0_dp* this%fom *this%grad

        call extend_array_to_ranks(this%grad, this%dimensions, this%nx, this%ny, this%nz, &
                                   this%n_ker)

    end if

end subroutine collect_gradients

!###################################################################################################

subroutine apply_kernel_on_rho(this)

    class(TDesign) :: this

    integer :: i, j, k, ii, jj, kk

    this%rho_conv = R_0

    select case (this%dimensions)
    case (1)
        do i = 1, this%nx
        do ii = -this%n_ker, this%n_ker
            if (this%opt_region(i+ii,1,1)) then
                this%rho_conv(i,1,1) = this%rho_conv(i,1,1) + &
                                       this%ker_mat(ii, 1, 1)*this%rho(i+ii,1,1)
            end if
        end do
        end do
    case (2)
        do i = 1, this%nx
        do j = 1, this%ny
        do ii = -this%n_ker, this%n_ker
        do jj = -this%n_ker, this%n_ker
            if (this%opt_region(i+ii,j+jj,1)) then
                this%rho_conv(i,j,1) = this%rho_conv(i,j,1) + &
                    this%ker_mat(ii,jj,1)*this%rho(i+ii,j+jj,1)
            end if
        end do
        end do
        end do
        end do
    case (3)
        do i = 1, this%nx
        do j = 1, this%ny
        do k = 1, this%nz
        do ii = -this%n_ker, this%n_ker
        do jj = -this%n_ker, this%n_ker
        do kk = -this%n_ker, this%n_ker
            if (this%opt_region(i+ii,j+jj,k+kk)) then
                this%rho_conv(i,j,k) = this%rho_conv(i,j,k) + &
                    this%ker_mat(ii,jj,kk)*this%rho(i+ii,j+jj,k+kk)
            end if
        end do
        end do
        end do
        end do
        end do
        end do
    end select

end subroutine apply_kernel_on_rho

!###################################################################################################

subroutine apply_kernel_on_grad(this)

    class(TDesign) :: this

    integer :: i, j, k, ii, jj, kk

    this%grad_conv = R_0

    select case (this%dimensions)
    case (1)
        do i = 1, this%nx
        do ii = -this%n_ker, this%n_ker
            if (this%opt_region(i+ii,1,1)) then
                this%grad_conv(i+ii,1,1) = this%grad_conv(i+ii,1,1) + &
                                       this%ker_mat(-ii, 1, 1)*this%grad(i,1,1)
            end if
        end do
        end do
    case (2)
        do i = 1, this%nx
        do j = 1, this%ny
        do ii = -this%n_ker, this%n_ker
        do jj = -this%n_ker, this%n_ker
            if (this%opt_region(i+ii,j+jj,1)) then
                this%grad_conv(i+ii,j+jj,1) = this%grad_conv(i+ii,j+jj,1) + &
                    this%ker_mat(-ii,-jj,1)*this%grad(i,j,1)
            end if
        end do
        end do
        end do
        end do
    case (3)
        do i = 1, this%nx
        do j = 1, this%ny
        do k = 1, this%nz
        do ii = -this%n_ker, this%n_ker
        do jj = -this%n_ker, this%n_ker
        do kk = -this%n_ker, this%n_ker
            if (this%opt_region(i+ii,j+jj,k+kk)) then
                this%grad_conv(i+ii,j+jj,k+kk) = this%grad_conv(i+ii,j+jj,k+kk) + &
                    this%ker_mat(-ii,-jj,-kk)*this%grad(i,j,k)
            end if
        end do
        end do
        end do
        end do
        end do
        end do
    end select

end subroutine apply_kernel_on_grad

!###################################################################################################

subroutine displace_rho(this)

    class(TDesign) :: this

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

    select case (this%dimensions)
    case (1)
        this%rho(1:this%nx,1,1) = this%rho_old(1:this%nx,1,1) + &
            this%grad_new(1:this%nx,1,1) * this%drho / norm_global
    case (2)
        this%rho(1:this%nx,1:this%ny,1) = this%rho_old(1:this%nx,1:this%ny,1) + &
            this%grad_new(1:this%nx,1:this%ny,1) * this%drho / norm_global
    case (3)
        this%rho(1:this%nx,1:this%ny,1:this%nz) = this%rho_old(1:this%nx,1:this%ny,1:this%nz) + &
            this%grad_new(1:this%nx,1:this%ny,1:this%nz) * this%drho / norm_global
    end select

    call extend_array_to_ranks(this%rho, this%dimensions, this%nx, this%ny, this%nz, this%n_ker)

end subroutine displace_rho

!###################################################################################################

subroutine reset_rho_one_step_back(this)

    class(TDesign) :: this

    select case (this%dimensions)
    case (1)
        this%rho(1:this%nx,1,1) = this%rho_old(1:this%nx,1,1)
    case (2)
        this%rho(1:this%nx,1:this%ny,1) = this%rho_old(1:this%nx,1:this%ny,1)
    case (3)
        this%rho(1:this%nx,1:this%ny,1:this%nz) = this%rho_old(1:this%nx,1:this%ny,1:this%nz)
    end select

    call extend_array_to_ranks(this%rho, this%dimensions, this%nx, this%ny, this%nz, this%n_ker)

end subroutine reset_rho_one_step_back

!###################################################################################################

subroutine reset_grad(this)

    class(TDesign) :: this

    select case (this%dimensions)
    case (1)
        this%grad(1:this%nx,1,1) = R_0
    case (2)
        this%grad(1:this%nx,1:this%ny,1) = R_0
    case (3)
        this%grad(1:this%nx,1:this%ny,1:this%nz) = R_0
    end select

    call extend_array_to_ranks(this%grad, this%dimensions, this%nx, this%ny, this%nz, this%n_ker)

end subroutine reset_grad

!###################################################################################################

subroutine update_rho(this)

    class(TDesign) :: this

    select case (this%dimensions)
    case (1)
        this%rho_old(1:this%nx,1,1) = this%rho(1:this%nx,1,1)
    case (2)
        this%rho_old(1:this%nx,1:this%ny,1) = this%rho(1:this%nx,1:this%ny,1)
    case (3)
        this%rho_old(1:this%nx,1:this%ny,1:this%nz) = this%rho(1:this%nx,1:this%ny,1:this%nz)
    end select

end subroutine update_rho

!###################################################################################################

subroutine update_grad(this)

    class(TDesign) :: this

    select case (this%dimensions)
    case (1)
        this%grad_old(1:this%nx,1,1) = this%grad_conv(1:this%nx,1,1)
    case (2)
        this%grad_old(1:this%nx,1:this%ny,1) = this%grad_conv(1:this%nx,1:this%ny,1)
    case (3)
        this%grad_old(1:this%nx,1:this%ny,1:this%nz) = this%grad_conv(1:this%nx,1:this%ny,1:this%nz)
    end select

end subroutine update_grad

!###################################################################################################

end module design_mod