module target_source_mod

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
        
        integer,  allocatable :: i_ndx
        integer,  allocatable :: j_ndx
        integer,  allocatable :: k_ndx
        logical,  allocatable :: in_this_rank(:,:,:)
        real(dp), allocatable :: ker_mat(:,:,:)

        contains

            procedure :: init_src_trg
            procedure :: kill_src_trg

    end type TSrcTrg

contains
!###################################################################################################

subroutine init_src_trg(this, dim, r0, dr, type_ch, sigma, J_amp, grid_Ndims, mpi_coords, mpi_dims)

    class(TSrcTrg)   , intent(inout) :: this
    character(len=10), intent(in)    :: type_ch
    complex(dp)      , intent(in)    :: J_amp
    real(dp)         , intent(in)    :: r0(3)
    real(dp)         , intent(in)    :: dr
    real(dp)         , intent(in)    :: sigma
    integer          , intent(in)    :: dim
    integer          , intent(in)    :: grid_Ndims(3)
    integer          , intent(in)    :: mpi_coords(3)
    integer          , intent(in)    :: mpi_dims(3)

    integer  :: i, j, k
    integer  :: i0, j0, k0
    integer  :: rank_x, rank_y, rank_z
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

        i0 = int(this%r0(1)/dr) + int(grid_Ndims(1)/2)        
        
        this%i_ndx = i0
        norm = R_0
        do i = -this%n_sig, this%n_sig
            this%ker_mat(i, 1, 1) = exp(-(i*dr)**2/(2*this%sigma**2))
            norm = norm + this%ker_mat(i, 1, 1)
            this%in_this_rank(i, 1, 1) = .true.
        end do

        this%ker_mat = this%ker_mat/norm

    case (2)
        if (.not.allocated(this%ker_mat)) allocate(this%ker_mat(-this%n_sig:this%n_sig, &
                                                   -this%n_sig:this%n_sig, 1))
        if (.not.allocated(this%in_this_rank)) &
        allocate(this%in_this_rank(-this%n_sig:this%n_sig, -this%n_sig:this%n_sig, 1))

        i0 = int(this%r0(1)/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
        j0 = int(this%r0(2)/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)

        rank_x = int((i0+1)/grid_Ndims(1))
        rank_y = int((j0+1)/grid_Ndims(2))

        if (rank_x == mpi_coords(1) .and. rank_y == mpi_coords(2)) then
            this%i_ndx = i0 - rank_x*grid_Ndims(1)
            this%j_ndx = j0 - rank_y*grid_Ndims(2)
            this%in_this_rank(this%i_ndx, this%j_ndx, 1) = .true.
        end if

        norm = R_0
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

        norm = R_0
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

end subroutine kill_src_trg
!###################################################################################################

subroutine read_init_src_trg(src_trg_list, n_src_trg, n_src, n_trg, dimensions, id, dr, &
                             grid_Ndims, mpi_coords, mpi_dims)

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
                                           grid_Ndims, mpi_coords, mpi_dims)
        
        if (src_trg_list(i)%is_a_source) n_src = n_src + 1
        if (src_trg_list(i)%is_a_target) n_trg = n_trg + 1

    end do

end subroutine read_init_src_trg
!###################################################################################################

subroutine set_source_J(src_trg, j_vec)
    
    type(TSrcTrg) , intent(in)    :: src_trg
    class(TRSvec)  , intent(inout) :: j_vec

    integer :: i, j, k
    integer :: i_ndx, j_ndx, k_ndx

    if (.not. src_trg%is_a_source) return

    select type (j_vec)
    type is (TRSvec_1D)
        i_ndx = src_trg%i_ndx
        select case (src_trg%type)
        case (Jx_SOURCE)
            do i = -src_trg%n_sig, src_trg%n_sig
                j_vec%pl_x(i_ndx + i) = src_trg%J_amp * src_trg%ker_mat(i, 1, 1)
            end do
            j_vec%mi_x = j_vec%pl_x
        end select
    type is (TRSvec_2D)
        i_ndx = src_trg%i_ndx
        j_ndx = src_trg%j_ndx

        select case (src_trg%type)
        case (Jx_SOURCE)

            do j = -src_trg%n_sig, src_trg%n_sig
            do i = -src_trg%n_sig, src_trg%n_sig
                if (src_trg%in_this_rank(i, j, 1)) then
                    j_vec%pl_x(i_ndx + i, j_ndx + j) = src_trg%J_amp * src_trg%ker_mat(i, j, 1)
                end if
            end do
            end do

            j_vec%mi_x = j_vec%pl_x

        case (Jy_SOURCE)

            do j = -src_trg%n_sig, src_trg%n_sig
            do i = -src_trg%n_sig, src_trg%n_sig
                if (src_trg%in_this_rank(i, j, 1)) then
                    j_vec%pl_y(i_ndx + i, j_ndx + j) = src_trg%J_amp * src_trg%ker_mat(i, j, 1)
                end if
            end do
            end do

            j_vec%mi_y = j_vec%pl_y

        case (Jz_SOURCE)

            do j = -src_trg%n_sig, src_trg%n_sig
            do i = -src_trg%n_sig, src_trg%n_sig
                if (src_trg%in_this_rank(i, j, 1)) then
                    j_vec%pl_z(i_ndx + i, j_ndx + j) = src_trg%J_amp * src_trg%ker_mat(i, j, 1)
                end if
            end do
            end do

            j_vec%mi_z = j_vec%pl_z
            
        end select
    type is (TRSvec_3D)

        i_ndx = src_trg%i_ndx
        j_ndx = src_trg%j_ndx
        k_ndx = src_trg%k_ndx

        select case (src_trg%type)

        case (Jx_SOURCE)
            do k = -src_trg%n_sig, src_trg%n_sig
            do j = -src_trg%n_sig, src_trg%n_sig
            do i = -src_trg%n_sig, src_trg%n_sig
                if (src_trg%in_this_rank(i, j, k)) then
                    j_vec%pl_x(i_ndx + i, j_ndx + j, k_ndx + k) = &
                        src_trg%J_amp * src_trg%ker_mat(i, j, k)
                end if
            end do
            end do
            end do
            j_vec%mi_x = j_vec%pl_x

        case (Jy_SOURCE)
            do k = -src_trg%n_sig, src_trg%n_sig
            do j = -src_trg%n_sig, src_trg%n_sig
            do i = -src_trg%n_sig, src_trg%n_sig
                if (src_trg%in_this_rank(i, j, k)) then
                    j_vec%pl_y(i_ndx + i, j_ndx + j, k_ndx + k) = &
                        src_trg%J_amp * src_trg%ker_mat(i, j, k)
                end if
            end do
            end do
            end do
            j_vec%mi_y = j_vec%pl_y

        case (Jz_SOURCE)
            do k = -src_trg%n_sig, src_trg%n_sig
            do j = -src_trg%n_sig, src_trg%n_sig
            do i = -src_trg%n_sig, src_trg%n_sig
                if (src_trg%in_this_rank(i, j, k)) then
                    j_vec%pl_z(i_ndx + i, j_ndx + j, k_ndx + k) = &
                        src_trg%J_amp * src_trg%ker_mat(i, j, k)
                end if
            end do
            end do
            end do
            j_vec%mi_z = j_vec%pl_z

        end select
    end select

end subroutine set_source_J

!###################################################################################################

end module target_source_mod