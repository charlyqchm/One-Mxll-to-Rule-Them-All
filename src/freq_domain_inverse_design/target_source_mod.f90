module target_source_mod

#ifdef USE_MPI
    use mpi
#endif

    use constants_mod
    use rs_vec_base_mod
    use rs_vec_dimensions_mod

    implicit none

    type TSrcTrg

        complex(dp) :: J_amp
        integer     :: type
        integer     :: dimensions
        integer     :: n_sig
        logical     :: is_a_source
        logical     :: is_a_target
        real(dp)    :: sigma
        real(dp)    :: r0(3)
        
        complex(dp), allocatable :: trg_amp(:,:,:)
        integer    , allocatable :: i_ndx(:,:,:)
        integer    , allocatable :: j_ndx(:,:,:)
        integer    , allocatable :: k_ndx(:,:,:)
        integer    , allocatable :: rank_id(:,:,:)
        logical    , allocatable :: in_this_rank(:,:,:)
        real(dp)   , allocatable :: ker_mat(:,:,:)

        contains

            procedure :: init_src_trg
            procedure :: kill_src_trg

    end type TSrcTrg

contains
!###################################################################################################

subroutine init_src_trg(this, dim, r0, dr, type_ch, sigma, J_amp, grid_Ndims, &
                        mpi_cart_comm, mpi_coords, mpi_dims)

    class(TSrcTrg)   , intent(inout) :: this
    character(len=10), intent(in)    :: type_ch
    complex(dp)      , intent(in)    :: J_amp
    real(dp)         , intent(in)    :: r0(3)
    real(dp)         , intent(in)    :: dr
    real(dp)         , intent(in)    :: sigma
    integer          , intent(in)    :: dim
    integer          , intent(in)    :: grid_Ndims(3)
    integer          , intent(in)    :: mpi_cart_comm
    integer          , intent(in)    :: mpi_coords(3)
    integer          , intent(in)    :: mpi_dims(3)

    integer  :: i, j, k
    integer  :: i0, j0, k0
    integer  :: rank_id
    integer  :: rank_x, rank_y, rank_z
    integer  :: rank_vec(3)
    integer  :: ierr
    real(dp) :: norm

    this%sigma      = sigma
    this%r0         = r0
    this%J_amp      = J_amp
    this%dimensions = dim

    select case (type_ch)
        case ("Jx")
            this%is_a_source = .true.
            this%is_a_target = .false.
            this%type = Jx_SOURCE
        case ("Jy")
            this%is_a_source = .true.
            this%is_a_target = .false.
            this%type = Jy_SOURCE
            if (dim == 1) then
                write (*, '("Error: Jy source can only be defined in 2D/3D simulations")')
                error stop
            end if
        case ("Jz")
            this%is_a_source = .true.
            this%is_a_target = .false.
            this%type = Jz_SOURCE
            if (dim == 1) then
                write (*, '("Error: Jz source can only be defined in 2D/3D simulations")')
                error stop
            end if
        case ("Re_Ex")
            this%is_a_source = .false.
            this%is_a_target = .true.
            this%type = Re_Ex_TARGET
        case ("Im_Ex")
            this%is_a_source = .false.
            this%is_a_target = .true.
            this%type = Im_Ex_TARGET
        case ("Re_Ey")
            this%is_a_source = .false.
            this%is_a_target = .true.
            this%type = Re_Ey_TARGET
            if (dim == 1) then
                write (*, '("Error: Re_Ey target can only be defined in 2D/3D simulations")')
                error stop
            end if
        case ("Im_Ey")
            this%is_a_source = .false.
            this%is_a_target = .true.
            this%type = Im_Ey_TARGET
            if (dim == 1) then
                write (*, '("Error: Im_Ey target can only be defined in 2D/3D simulations")')
                error stop
            end if
        case ("Re_Ez")
            this%is_a_source = .false.
            this%is_a_target = .true.
            this%type = Re_Ez_TARGET
            if (dim == 1) then
                write (*, '("Error: Re_Ez target can only be defined in 2D/3D simulations")')
                error stop
            end if
        case ("Im_Ez")
            this%is_a_source = .false.
            this%is_a_target = .true.
            this%type = Im_Ez_TARGET
            if (dim == 1) then
                write (*, '("Error: Im_Ez target can only be defined in 2D/3D simulations")')
                error stop
            end if
        case ("Abs_Ex")
            this%is_a_source = .false.
            this%is_a_target = .true.
            this%type = Abs_Ex_TARGET
        case ("Abs_Ey")
            this%is_a_source = .false.
            this%is_a_target = .true.
            this%type = Abs_Ey_TARGET
            if (dim == 1) then
                write (*, '("Error: Abs_Ey target can only be defined in 2D/3D simulations")')
                error stop
            end if
        case ("Abs_Ez")
            this%is_a_source = .false.
            this%is_a_target = .true.
            this%type = Abs_Ez_TARGET
            if (dim == 1) then
                write (*, '("Error: Abs_Ez target can only be defined in 2D/3D simulations")')
                error stop
            end if
        case default
            write (*, '("Error: invalid source/target type")')
            error stop
    end select

    this%n_sig = int(3*this%sigma/dr)

    select case (this%dimensions)

    case (1)
        if (.not.allocated(this%ker_mat))      allocate(this%ker_mat(-this%n_sig:this%n_sig, 1, 1))
        if (.not.allocated(this%in_this_rank)) &
        allocate(this%in_this_rank(-this%n_sig:this%n_sig, 1, 1))
        if (.not.allocated(this%rank_id))      allocate(this%rank_id(-this%n_sig:this%n_sig, 1, 1))
        if (.not.allocated(this%trg_amp))      allocate(this%trg_amp(-this%n_sig:this%n_sig, 1, 1))
        if (.not.allocated(this%i_ndx))        allocate(this%i_ndx(-this%n_sig:this%n_sig, 1, 1))

        this%in_this_rank = .false.

        i0 = int(this%r0(1)/dr) + int(grid_Ndims(1)/2)        
        
        this%i_ndx = 0
        norm = 0.0_dp
        do i = -this%n_sig, this%n_sig
            this%ker_mat(i, 1, 1) = exp(-(i*dr)**2/(2*this%sigma**2))
            norm = norm + this%ker_mat(i, 1, 1)
            this%in_this_rank(i, 1, 1) = .true.
            this%rank_id(i, 1, 1) = 0
            this%i_ndx(i, 1, 1) = i + i0
        end do

        this%ker_mat = this%ker_mat/norm

    case (2)
        if (.not.allocated(this%ker_mat)) allocate(this%ker_mat(-this%n_sig:this%n_sig, &
                                                   -this%n_sig:this%n_sig, 1))
        if (.not.allocated(this%in_this_rank)) &
        allocate(this%in_this_rank(-this%n_sig:this%n_sig, -this%n_sig:this%n_sig, 1))
        if (.not.allocated(this%rank_id)) &
        allocate(this%rank_id(-this%n_sig:this%n_sig, -this%n_sig:this%n_sig, 1))
        if (.not.allocated(this%trg_amp)) &
        allocate(this%trg_amp(-this%n_sig:this%n_sig, -this%n_sig:this%n_sig, 1))
        if (.not.allocated(this%i_ndx)) &
        allocate(this%i_ndx(-this%n_sig:this%n_sig, -this%n_sig:this%n_sig, 1))
        if (.not.allocated(this%j_ndx)) &
        allocate(this%j_ndx(-this%n_sig:this%n_sig, -this%n_sig:this%n_sig, 1))

        this%in_this_rank = .false.

        this%i_ndx = 0
        this%j_ndx = 0

        i0 = int(this%r0(1)/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
        j0 = int(this%r0(2)/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)

        do j = -this%n_sig, this%n_sig
        do i = -this%n_sig, this%n_sig

            rank_x = int((i0+i+1)/grid_Ndims(1))
            rank_y = int((j0+j+1)/grid_Ndims(2))

            rank_vec = (/rank_x, rank_y, 0/)

#ifdef USE_MPI
            call MPI_Cart_rank(mpi_cart_comm, rank_vec(1:2), rank_id, ierr)
            this%rank_id(i, j, 1) = rank_id
#else
            this%rank_id(i, j, 1) = 0
#endif

            if (rank_x == mpi_coords(1) .and. rank_y == mpi_coords(2)) then
                this%i_ndx(i, j, 1) = (i+i0) - rank_x*grid_Ndims(1)
                this%j_ndx(i, j, 1) = (j+j0) - rank_y*grid_Ndims(2)
                this%in_this_rank(i,j,1) = .true.
            end if

        end do
        end do

        norm = 0.0_dp
        do j = -this%n_sig, this%n_sig
        do i = -this%n_sig, this%n_sig
            this%ker_mat(i, j, 1) = exp(-(i*dr)**2/(2*this%sigma**2)) * &
                                       exp(-(j*dr)**2/(2*this%sigma**2))
            norm = norm + this%ker_mat(i, j, 1)
        end do
        end do

        this%ker_mat = this%ker_mat/norm

    case (3)
        if (.not.allocated(this%ker_mat)) &
        allocate(this%ker_mat(-this%n_sig:this%n_sig, -this%n_sig:this%n_sig, &
                              -this%n_sig:this%n_sig))
        if (.not.allocated(this%in_this_rank)) &
        allocate(this%in_this_rank(-this%n_sig:this%n_sig, -this%n_sig:this%n_sig, &
                                   -this%n_sig:this%n_sig))
        if (.not.allocated(this%rank_id)) &
        allocate(this%rank_id(-this%n_sig:this%n_sig, -this%n_sig:this%n_sig, &
                              -this%n_sig:this%n_sig))
        if (.not.allocated(this%trg_amp)) &
        allocate(this%trg_amp(-this%n_sig:this%n_sig, -this%n_sig:this%n_sig, &
                              -this%n_sig:this%n_sig))
        if (.not.allocated(this%i_ndx)) &
        allocate(this%i_ndx(-this%n_sig:this%n_sig, -this%n_sig:this%n_sig, &
                            -this%n_sig:this%n_sig))
        if (.not.allocated(this%j_ndx)) &
        allocate(this%j_ndx(-this%n_sig:this%n_sig, -this%n_sig:this%n_sig, &
                            -this%n_sig:this%n_sig))
        if (.not.allocated(this%k_ndx)) &
        allocate(this%k_ndx(-this%n_sig:this%n_sig, -this%n_sig:this%n_sig, &
                 -this%n_sig:this%n_sig))

        this%in_this_rank = .false.

        this%i_ndx = 0
        this%j_ndx = 0
        this%k_ndx = 0

        i0 = int(this%r0(1)/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
        j0 = int(this%r0(2)/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
        k0 = int(this%r0(3)/dr) + int(grid_Ndims(3)*mpi_dims(3)/2)

        do k = -this%n_sig, this%n_sig
        do j = -this%n_sig, this%n_sig
        do i = -this%n_sig, this%n_sig

            rank_x = int((i0+i+1)/grid_Ndims(1))
            rank_y = int((j0+j+1)/grid_Ndims(2))
            rank_z = int((k0+k+1)/grid_Ndims(3))

            rank_vec = (/rank_x, rank_y, rank_z/)

#ifdef USE_MPI
            call MPI_Cart_rank(mpi_cart_comm, rank_vec, rank_id, ierr)
            this%rank_id(i, j, k) = rank_id
#else
            this%rank_id(i, j, k) = 0
#endif
            if (rank_x == mpi_coords(1) .and. rank_y == mpi_coords(2) .and. &
                rank_z == mpi_coords(3)) then
                this%i_ndx(i, j, k) = (i+i0) - rank_x*grid_Ndims(1)
                this%j_ndx(i, j, k) = (j+j0) - rank_y*grid_Ndims(2)
                this%k_ndx(i, j, k) = (k+k0) - rank_z*grid_Ndims(3)
                this%in_this_rank(i,j,k) = .true.
            end if

        end do 
        end do
        end do

        norm = 0.0_dp
        do k = -this%n_sig, this%n_sig
        do j = -this%n_sig, this%n_sig
        do i = -this%n_sig, this%n_sig
            this%ker_mat(i, j, k) = exp(-(i*dr)**2/(2*this%sigma**2)) * &
                                       exp(-(j*dr)**2/(2*this%sigma**2)) * &
                                       exp(-(k*dr)**2/(2*this%sigma**2))
            norm = norm + this%ker_mat(i, j, k)
        end do
        end do
        end do

        this%ker_mat = this%ker_mat/norm

    end select

end subroutine init_src_trg

!###################################################################################################   
subroutine kill_src_trg(this)

    class(TSrcTrg), intent(inout) :: this

    if (allocated(this%ker_mat))      deallocate(this%ker_mat)
    if (allocated(this%in_this_rank)) deallocate(this%in_this_rank)
    if (allocated(this%rank_id))      deallocate(this%rank_id)
    if (allocated(this%trg_amp))      deallocate(this%trg_amp)

end subroutine kill_src_trg
!###################################################################################################

subroutine read_init_src_trg(src_trg_list, n_src_trg, n_src, n_trg, dimensions, id, dr, &
                             grid_Ndims, mpi_cart_comm, mpi_coords, mpi_dims)

    type(TSrcTrg), allocatable, intent(inout) :: src_trg_list(:)
    integer                   , intent(out)   :: n_src_trg
    integer                   , intent(out)   :: n_src
    integer                   , intent(out)   :: n_trg
    integer                   , intent(in)    :: dimensions
    integer                   , intent(in)    :: id
    real(dp)                  , intent(in)    :: dr
    integer                   , intent(in)    :: grid_Ndims(3)
    integer                   , intent(in)    :: mpi_coords(3)
    integer                   , intent(in)    :: mpi_dims(3)
    integer                   , intent(in)    :: mpi_cart_comm

    character(len=20) :: file_name = "src_trg_list_"
    character(len=20) :: file_exten = ".in"
    character(len=20) :: file_number
    character(len=30) :: input_name
    integer           :: ierr, funit
    integer           :: i, j, k
    real(dp)          :: r0(3)
    real(dp)          :: sigma
    complex(dp)       :: J_amp
    character(len=10) :: type_ch
    real(dp)          :: amp_Re, amp_Im

    write (file_number, '(I3.3)') id

    input_name = trim(file_name) // trim(file_number) // trim(file_exten)

    inquire (file=input_name, iostat=ierr)
    
    if (ierr /= 0) then
        write (*, '("Error: target/sources list file ", A, " does not exist")') input_name
        error stop
    end if

    open (action='read', file=input_name, iostat=ierr, newunit=funit)
    if (ierr /= 0) then
        write (*, '("Error: cannot open target/sources list file ", A)') input_name
        error stop
    end if

    do
        read (funit, *, iostat=ierr)
        if (ierr /= 0) exit
        n_src_trg = n_src_trg + 1 
    end do
    
    rewind(funit)
    
    if (.not. allocated(src_trg_list)) allocate(src_trg_list(n_src_trg))
    
    n_src = 0
    n_trg = 0
    do i = 1, n_src_trg
        read (funit, *, iostat=ierr) type_ch, amp_Re, amp_Im, r0(1), r0(2), r0(3), sigma
        J_amp = amp_Re*Z_ONE + amp_Im*Z_I
        call src_trg_list(i)%init_src_trg(dimensions, r0, dr, type_ch, sigma, J_amp, &
                                           grid_Ndims, mpi_cart_comm, mpi_coords, mpi_dims)
        
        if (src_trg_list(i)%is_a_source) n_src = n_src + 1
        if (src_trg_list(i)%is_a_target) n_trg = n_trg + 1

    end do

end subroutine read_init_src_trg
!###################################################################################################

subroutine set_source_J(src, j_vec)
    
    type(TSrcTrg) , intent(in)    :: src
    class(TRSvec)  , intent(inout) :: j_vec

    integer :: i, j, k
    integer :: i_ndx, j_ndx, k_ndx

    if (.not. src%is_a_source) return

    select type (j_vec)
    type is (TRSvec_1D)
        i_ndx = src%i_ndx(0,1,1)
        select case (src%type)
        case (Jx_SOURCE)
            do i = -src%n_sig, src%n_sig
                j_vec%pl_x(i_ndx + i) = src%J_amp * src%ker_mat(i, 1, 1)
            end do
            j_vec%mi_x = j_vec%pl_x
        end select
    type is (TRSvec_2D)
        i_ndx = src%i_ndx(0,0,1)
        j_ndx = src%j_ndx(0,0,1)

        select case (src%type)
        case (Jx_SOURCE)

            do j = -src%n_sig, src%n_sig
            do i = -src%n_sig, src%n_sig
                if (src%in_this_rank(i, j, 1)) then
                    j_vec%pl_x(i_ndx + i, j_ndx + j) = src%J_amp * src%ker_mat(i, j, 1)
                end if
            end do
            end do

            j_vec%mi_x = j_vec%pl_x

        case (Jy_SOURCE)

            do j = -src%n_sig, src%n_sig
            do i = -src%n_sig, src%n_sig
                if (src%in_this_rank(i, j, 1)) then
                    j_vec%pl_y(i_ndx + i, j_ndx + j) = src%J_amp * src%ker_mat(i, j, 1)
                end if
            end do
            end do

            j_vec%mi_y = j_vec%pl_y

        case (Jz_SOURCE)

            do j = -src%n_sig, src%n_sig
            do i = -src%n_sig, src%n_sig
                if (src%in_this_rank(i, j, 1)) then
                    j_vec%pl_z(i_ndx + i, j_ndx + j) = src%J_amp * src%ker_mat(i, j, 1)
                end if
            end do
            end do

            j_vec%mi_z = j_vec%pl_z
            
        end select
    type is (TRSvec_3D)

        i_ndx = src%i_ndx(0,0,0)
        j_ndx = src%j_ndx(0,0,0)
        k_ndx = src%k_ndx(0,0,0)

        select case (src%type)

        case (Jx_SOURCE)
            do k = -src%n_sig, src%n_sig
            do j = -src%n_sig, src%n_sig
            do i = -src%n_sig, src%n_sig
                if (src%in_this_rank(i, j, k)) then
                    j_vec%pl_x(i_ndx + i, j_ndx + j, k_ndx + k) = &
                        src%J_amp * src%ker_mat(i, j, k)
                end if
            end do
            end do
            end do
            j_vec%mi_x = j_vec%pl_x

        case (Jy_SOURCE)
            do k = -src%n_sig, src%n_sig
            do j = -src%n_sig, src%n_sig
            do i = -src%n_sig, src%n_sig
                if (src%in_this_rank(i, j, k)) then
                    j_vec%pl_y(i_ndx + i, j_ndx + j, k_ndx + k) = &
                        src%J_amp * src%ker_mat(i, j, k)
                end if
            end do
            end do
            end do
            j_vec%mi_y = j_vec%pl_y

        case (Jz_SOURCE)
            do k = -src%n_sig, src%n_sig
            do j = -src%n_sig, src%n_sig
            do i = -src%n_sig, src%n_sig
                if (src%in_this_rank(i, j, k)) then
                    j_vec%pl_z(i_ndx + i, j_ndx + j, k_ndx + k) = &
                        src%J_amp * src%ker_mat(i, j, k)
                end if
            end do
            end do
            end do
            j_vec%mi_z = j_vec%pl_z

        end select
    end select

end subroutine set_source_J

!###################################################################################################

subroutine update_target(trg, f_vec, w_dL, w_total)

    type(TSrcTrg)  , intent(inout) :: trg
    class(TRSvec)  , intent(inout) :: f_vec
    real(dp)       , intent(inout) :: w_dL
    real(dp)       , intent(inout) :: w_total

    complex(dp) :: E_field
    integer     :: i0, j0, k0

    if (.not. trg%is_a_target) return

    select type (f_vec)
    type is (TRSvec_1D)

        if (.not. trg%in_this_rank(0, 1, 1)) return

        i0 = trg%i_ndx(0, 1, 1)

        E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_x(i0)+f_vec%mi_x(i0))

        select case (trg%type)
        case(Re_Ex_TARGET)
            w_dL = DABS(DREAL(E_field))
            trg%J_amp = DSQRT(1.0_dp/(2.0_dp*eps0)) * DREAL(E_field)
        case(Im_Ex_TARGET)
            w_dL  = DABS(DIMAG(E_field))
            trg%J_amp = Z_I*DSQRT(1.0_dp/(2.0_dp*eps0)) * DIMAG(E_field)
        case(Abs_Ex_TARGET)
            w_dL  = ABS(E_field)
            trg%J_amp = -DSQRT(1.0_dp/(2.0_dp*eps0)) * DCONJG(E_field)
        end select

    type is (TRSvec_2D)

        if (.not. trg%in_this_rank(0, 0, 1)) return

        i0 = trg%i_ndx(0, 0, 1)
        j0 = trg%j_ndx(0, 0, 1)

        select case (trg%type)
        case(Re_Ex_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_x(i0, j0)+f_vec%mi_x(i0, j0))
            w_dL  = DABS(DREAL(E_field))
            trg%J_amp = DSQRT(1.0_dp/(2.0_dp*eps0)) * DREAL(E_field)
        case(Im_Ex_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_x(i0, j0)+f_vec%mi_x(i0, j0))
            w_dL  = DABS(DIMAG(E_field))
            trg%J_amp = Z_I*DSQRT(1.0_dp/(2.0_dp*eps0)) * DIMAG(E_field)
        case(Abs_Ex_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_x(i0, j0)+f_vec%mi_x(i0, j0))
            w_dL  = ABS(E_field)
            trg%J_amp = -DSQRT(1.0_dp/(2.0_dp*eps0)) * DCONJG(E_field)
        case(Re_Ey_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_y(i0, j0)+f_vec%mi_y(i0, j0))
            w_dL  = DABS(DREAL(E_field))
            trg%J_amp = DSQRT(1.0_dp/(2.0_dp*eps0)) * DREAL(E_field)
        case(Im_Ey_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_y(i0, j0)+f_vec%mi_y(i0, j0))
            w_dL  = DABS(DIMAG(E_field))
            trg%J_amp = Z_I*DSQRT(1.0_dp/(2.0_dp*eps0)) * DIMAG(E_field)
        case(Abs_Ey_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_y(i0, j0)+f_vec%mi_y(i0, j0))
            w_dL  = ABS(E_field)
            trg%J_amp = -DSQRT(1.0_dp/(2.0_dp*eps0)) * DCONJG(E_field)
        case(Re_Ez_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_z(i0, j0)+f_vec%mi_z(i0, j0))
            w_dL  = DABS(DREAL(E_field))
            trg%J_amp = DSQRT(1.0_dp/(2.0_dp*eps0)) * DREAL(E_field)
        case(Im_Ez_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_z(i0, j0)+f_vec%mi_z(i0, j0))
            w_dL  = DABS(DIMAG(E_field))
            trg%J_amp = Z_I*DSQRT(1.0_dp/(2.0_dp*eps0)) * DIMAG(E_field)
        case(Abs_Ez_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_z(i0, j0)+f_vec%mi_z(i0, j0))
            w_dL  = ABS(E_field)
            trg%J_amp = -DSQRT(1.0_dp/(2.0_dp*eps0)) * DCONJG(E_field)
        end select
    type is (TRSvec_3D)

        if (.not. trg%in_this_rank(0, 0, 0)) return

        i0 = trg%i_ndx(0, 0, 0)
        j0 = trg%j_ndx(0, 0, 0)
        k0 = trg%k_ndx(0, 0, 0)

        select case (trg%type)
        case(Re_Ex_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_x(i0, j0, k0)+f_vec%mi_x(i0, j0, k0))
            w_dL  = DABS(DREAL(E_field))
            trg%J_amp = DSQRT(1.0_dp/(2.0_dp*eps0)) * DREAL(E_field)
        case(Im_Ex_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_x(i0, j0, k0)+f_vec%mi_x(i0, j0, k0))
            w_dL  = DABS(DIMAG(E_field))
            trg%J_amp = Z_I*DSQRT(1.0_dp/(2.0_dp*eps0)) * DIMAG(E_field)
        case(Abs_Ex_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_x(i0, j0, k0)+f_vec%mi_x(i0, j0, k0))
            w_dL  = ABS(E_field)
            trg%J_amp = -DSQRT(1.0_dp/(2.0_dp*eps0)) * DCONJG(E_field)
        case(Re_Ey_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_y(i0, j0, k0)+f_vec%mi_y(i0, j0, k0))
            w_dL  = DABS(DREAL(E_field))
            trg%J_amp = DSQRT(1.0_dp/(2.0_dp*eps0)) * DREAL(E_field)
        case(Im_Ey_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_y(i0, j0, k0)+f_vec%mi_y(i0, j0, k0))
            w_dL  = DABS(DIMAG(E_field))
            trg%J_amp = Z_I*DSQRT(1.0_dp/(2.0_dp*eps0)) * DIMAG(E_field)
        case(Abs_Ey_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_y(i0, j0, k0)+f_vec%mi_y(i0, j0, k0))
            w_dL  = ABS(E_field)
            trg%J_amp = -DSQRT(1.0_dp/(2.0_dp*eps0)) * DCONJG(E_field)
        case(Re_Ez_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_z(i0, j0, k0)+f_vec%mi_z(i0, j0, k0))
            w_dL  = DABS(DREAL(E_field))
            trg%J_amp = DSQRT(1.0_dp/(2.0_dp*eps0)) * DREAL(E_field)
        case(Im_Ez_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_z(i0, j0, k0)+f_vec%mi_z(i0, j0, k0))
            w_dL  = DABS(DIMAG(E_field))
            trg%J_amp = Z_I*DSQRT(1.0_dp/(2.0_dp*eps0)) * DIMAG(E_field)
        case(Abs_Ez_TARGET)
            E_field = DSQRT(1.0_dp/(2.0_dp*eps0)) * (f_vec%pl_z(i0, j0, k0)+f_vec%mi_z(i0, j0, k0))
            w_dL  = ABS(E_field)
            trg%J_amp = -DSQRT(1.0_dp/(2.0_dp*eps0)) * DCONJG(E_field)
        end select
    end select

end subroutine update_target

!###################################################################################################

subroutine set_jtrg(trg, j_vec)

    type(TSrcTrg)  , intent(in)    :: trg
    class(TRSvec)  , intent(inout) :: j_vec

    logical :: source_in_this_rank
    integer :: i, j, k
    integer :: i_ndx, j_ndx, k_ndx
    complex(dp) :: J_recv, J_sent

#ifdef USE_MPI
    integer :: istatus(MPI_STATUS_SIZE)
    integer :: ierr
#endif

    if (.not. trg%is_a_target) return

    select type (j_vec)
    type is (TRSvec_1D)

        j_vec%pl_x = 0.0_dp
        j_vec%mi_x = 0.0_dp
        j_vec%pl_y = 0.0_dp
        j_vec%mi_y = 0.0_dp

        select case (trg%type)
        case (Re_Ex_TARGET, Im_Ex_TARGET, Abs_Ex_TARGET)

            do i = -trg%n_sig, trg%n_sig
                
                i_ndx = trg%i_ndx(i,1,1)
                
                j_vec%pl_x(i_ndx) = trg%J_amp * trg%ker_mat(i, 1, 1)

            end do
        end select

    type is (TRSvec_2D)

        j_vec%pl_x = 0.0_dp
        j_vec%mi_x = 0.0_dp
        j_vec%pl_y = 0.0_dp
        j_vec%mi_y = 0.0_dp
        j_vec%pl_z = 0.0_dp
        j_vec%mi_z = 0.0_dp

        source_in_this_rank = trg%in_this_rank(0, 0, 1)

        do j = -trg%n_sig, trg%n_sig
        do i = -trg%n_sig, trg%n_sig

            i_ndx = trg%i_ndx(i,j,1)
            j_ndx = trg%j_ndx(i,j,1)

            select case (trg%type)
            case (Re_Ex_TARGET, Im_Ex_TARGET, Abs_Ex_TARGET)
                if (trg%in_this_rank(i, j, 1) .and. source_in_this_rank) then
                    j_vec%pl_x(i_ndx, j_ndx) = trg%J_amp * trg%ker_mat(i, j, 1)
#ifdef USE_MPI
                else if (trg%in_this_rank(i, j, 1) .and. .not. source_in_this_rank) then
                    call mpi_recv(J_recv,1,mpi_double_complex, &
                                  trg%rank_id(0,0,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)
                    trg%J_amp = J_recv
                    j_vec%pl_x(i_ndx, j_ndx) = trg%J_amp * trg%ker_mat(i, j, 1)

                else if (.not. trg%in_this_rank(i, j, 1) .and. source_in_this_rank) then
                    J_sent = trg%J_amp
                    call mpi_send(J_sent,1,mpi_double_complex, &
                                  trg%rank_id(i, j, 1), MPI_GOOD_TAG, MPI_COMM_WORLD, ierr)
#endif
                end if

            case (Re_Ey_TARGET, Im_Ey_TARGET, Abs_Ey_TARGET)
                if (trg%in_this_rank(i, j, 1) .and. source_in_this_rank) then
                    j_vec%pl_y(i_ndx, j_ndx) = trg%J_amp * trg%ker_mat(i, j, 1)
#ifdef USE_MPI
                else if (trg%in_this_rank(i, j, 1) .and. .not. source_in_this_rank) then
                    call mpi_recv(J_recv,1,mpi_double_complex, &
                                  trg%rank_id(0,0,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)
                    trg%J_amp = J_recv
                    j_vec%pl_y(i_ndx, j_ndx) = trg%J_amp * trg%ker_mat(i, j, 1)

                else if (.not. trg%in_this_rank(i, j, 1) .and. source_in_this_rank) then
                    J_sent = trg%J_amp
                    call mpi_send(J_sent,1,mpi_double_complex, &
                                  trg%rank_id(i, j, 1), MPI_GOOD_TAG, MPI_COMM_WORLD, ierr)
#endif
                end if

            case (Re_Ez_TARGET, Im_Ez_TARGET, Abs_Ez_TARGET)
                if (trg%in_this_rank(i, j, 1) .and. source_in_this_rank) then
                    j_vec%pl_z(i_ndx, j_ndx) = trg%J_amp * trg%ker_mat(i, j, 1)
#ifdef USE_MPI
                else if (trg%in_this_rank(i, j, 1) .and. .not. source_in_this_rank) then
                    call mpi_recv(J_recv,1,mpi_double_complex, &
                                  trg%rank_id(0,0,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)
                    trg%J_amp = J_recv
                    j_vec%pl_z(i_ndx, j_ndx) = trg%J_amp * trg%ker_mat(i, j, 1)

                else if (.not. trg%in_this_rank(i, j, 1) .and. source_in_this_rank) then
                    J_sent = trg%J_amp
                    call mpi_send(J_sent,1,mpi_double_complex, &
                                  trg%rank_id(i, j, 1), MPI_GOOD_TAG, MPI_COMM_WORLD, ierr)
#endif
                end if

            end select
        end do
        end do

    type is (TRSvec_3D)

        j_vec%pl_x = 0.0_dp
        j_vec%mi_x = 0.0_dp
        j_vec%pl_y = 0.0_dp
        j_vec%mi_y = 0.0_dp
        j_vec%pl_z = 0.0_dp
        j_vec%mi_z = 0.0_dp

         source_in_this_rank = trg%in_this_rank(0, 0, 0)

        do k = -trg%n_sig, trg%n_sig
        do j = -trg%n_sig, trg%n_sig
        do i = -trg%n_sig, trg%n_sig

            i_ndx = trg%i_ndx(i,j,k)
            j_ndx = trg%j_ndx(i,j,k)
            k_ndx = trg%k_ndx(i,j,k)

            select case (trg%type)
            case (Re_Ex_TARGET, Im_Ex_TARGET, Abs_Ex_TARGET)
                if (trg%in_this_rank(i, j, k) .and. source_in_this_rank) then
                    j_vec%pl_x(i_ndx, j_ndx, k_ndx) = trg%J_amp * trg%ker_mat(i, j, k)
#ifdef USE_MPI
                else if (trg%in_this_rank(i, j, k) .and. .not. source_in_this_rank) then
                    call mpi_recv(J_recv,1,mpi_double_complex, &
                                  trg%rank_id(0,0,0), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)
                    trg%J_amp = J_recv
                    j_vec%pl_x(i_ndx, j_ndx, k_ndx) = trg%J_amp * trg%ker_mat(i, j, k)
                else if (.not. trg%in_this_rank(i, j, k) .and. source_in_this_rank) then
                    J_sent = trg%J_amp
                    call mpi_send(J_sent,1,mpi_double_complex, &
                                  trg%rank_id(i, j, k), MPI_GOOD_TAG, MPI_COMM_WORLD, ierr)
#endif
                end if

            case (Re_Ey_TARGET, Im_Ey_TARGET, Abs_Ey_TARGET)
                if (trg%in_this_rank(i, j, k) .and. source_in_this_rank) then
                    j_vec%pl_y(i_ndx, j_ndx, k_ndx) = trg%J_amp * trg%ker_mat(i, j, k)
#ifdef USE_MPI
                else if (trg%in_this_rank(i, j, k) .and. .not. source_in_this_rank) then
                    call mpi_recv(J_recv,1,mpi_double_complex, &
                                  trg%rank_id(0,0,0), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)
                    trg%J_amp = J_recv
                    j_vec%pl_y(i_ndx, j_ndx, k_ndx) = trg%J_amp * trg%ker_mat(i, j, k)
                else if (.not. trg%in_this_rank(i, j, k) .and. source_in_this_rank) then
                    J_sent = trg%J_amp
                    call mpi_send(J_sent,1,mpi_double_complex, &
                                  trg%rank_id(i, j, k), MPI_GOOD_TAG, MPI_COMM_WORLD, ierr)
#endif
                end if

            case (Re_Ez_TARGET, Im_Ez_TARGET, Abs_Ez_TARGET)
                if (trg%in_this_rank(i, j, k) .and. source_in_this_rank) then
                    j_vec%pl_z(i_ndx, j_ndx, k_ndx) = trg%J_amp * trg%ker_mat(i, j, k)
#ifdef USE_MPI
                else if (trg%in_this_rank(i, j, k) .and. .not. source_in_this_rank) then
                    call mpi_recv(J_recv,1,mpi_double_complex, &
                                  trg%rank_id(0,0,0), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)
                    trg%J_amp = J_recv
                    j_vec%pl_z(i_ndx, j_ndx, k_ndx) = trg%J_amp * trg%ker_mat(i, j, k)
                else if (.not. trg%in_this_rank(i, j, k) .and. source_in_this_rank) then
                    J_sent = trg%J_amp
                    call mpi_send(J_sent,1,mpi_double_complex, &
                                  trg%rank_id(i, j, k), MPI_GOOD_TAG, MPI_COMM_WORLD, ierr)
#endif
                end if

            end select

        end do
        end do
        end do
    end select

end subroutine set_jtrg 
!###################################################################################################

end module target_source_mod