module source_mod

    use constants_mod

    implicit none

    type TSource

        character(len=2) :: dir
        integer          :: dim
        integer          :: n_ker
        real(dp)         :: j_amp
        real(dp)         :: w0
        real(dp)         :: t0
        real(dp)         :: tau
        real(dp)         :: t_init
        real(dp)         :: t_final
        real(dp)         :: phase
        real(dp)         :: r0(3)
        real(dp)         :: rad

        integer , allocatable  :: ind_i(:, :, :)
        integer , allocatable  :: ind_j(:, :, :)
        integer , allocatable  :: ind_k(:, :, :)
        real(dp), allocatable  :: ker_mat(:,:,:)
        real(dp), allocatable  :: J_mat(:,:,:)
        logical , allocatable  :: in_this_rank(:,:,:)

        contains
            procedure :: init_source, kill_source, update_source

    end type TSource   

contains

!###################################################################################################
    subroutine init_source(this, direction, dimensions, j_amp, freq, t0, tau, &
                           r0, radius, t_init, t_final, phase, dr, grid_Ndims, mpi_coords, mpi_dims)

        class(TSource)  ,intent(inout) :: this
        character(len=2),intent(in)    :: direction
        integer         ,intent(in)    :: dimensions
        real(dp)        ,intent(in)    :: j_amp
        real(dp)        ,intent(in)    :: freq
        real(dp)        ,intent(in)    :: t0
        real(dp)        ,intent(in)    :: tau
        real(dp)        ,intent(in)    :: t_init
        real(dp)        ,intent(in)    :: t_final
        real(dp)        ,intent(in)    :: phase
        real(dp)        ,intent(in)    :: r0(3)
        real(dp)        ,intent(in)    :: radius
        integer         ,intent(in)    :: grid_Ndims(3)
        real(dp)        ,intent(in)    :: dr
        integer         ,intent(in)    :: mpi_coords(3)
        integer         ,intent(in)    :: mpi_dims(3)
        
        integer :: n_ker

        this%dim     = dimensions
        this%dir     = direction
        this%j_amp   = j_amp
        this%w0      = freq
        this%t0      = t0
        this%tau     = tau
        this%r0      = r0
        this%rad     = radius
        this%t_init  = t_init
        this%t_final = t_final
        this%phase   = phase

        n_ker = int(this%rad/dr)

        this%n_ker = n_ker

        select case (this%dim)
        case (1)
            if (.not. allocated(this%ind_i))        allocate(this%ind_i(-n_ker:n_ker, 1, 1))
            if (.not. allocated(this%ker_mat))      allocate(this%ker_mat(-n_ker:n_ker, 1, 1))
            if (.not. allocated(this%J_mat))        allocate(this%J_mat(-n_ker:n_ker, 1, 1))
            if (.not. allocated(this%in_this_rank)) allocate(this%in_this_rank(-n_ker:n_ker, 1, 1))
        case (2)
            if (.not. allocated(this%ind_i))        allocate(this%ind_i(-n_ker:n_ker, -n_ker:n_ker, 1))
            if (.not. allocated(this%ind_j))        allocate(this%ind_j(-n_ker:n_ker, -n_ker:n_ker, 1))
            if (.not. allocated(this%ker_mat))      allocate(this%ker_mat(-n_ker:n_ker, -n_ker:n_ker, 1))
            if (.not. allocated(this%J_mat))        allocate(this%J_mat(-n_ker:n_ker, -n_ker:n_ker, 1))
            if (.not. allocated(this%in_this_rank)) allocate(this%in_this_rank(-n_ker:n_ker, -n_ker:n_ker, 1))
        case (3)
            if (.not. allocated(this%ind_i))        allocate(this%ind_i(-n_ker:n_ker, -n_ker:n_ker, -n_ker:n_ker))
            if (.not. allocated(this%ind_j))        allocate(this%ind_j(-n_ker:n_ker, -n_ker:n_ker, -n_ker:n_ker))
            if (.not. allocated(this%ind_k))        allocate(this%ind_k(-n_ker:n_ker, -n_ker:n_ker, -n_ker:n_ker))   
            if (.not. allocated(this%ker_mat))      allocate(this%ker_mat(-n_ker:n_ker, -n_ker:n_ker, -n_ker:n_ker))
            if (.not. allocated(this%J_mat))        allocate(this%J_mat(-n_ker:n_ker, -n_ker:n_ker, -n_ker:n_ker))
            if (.not. allocated(this%in_this_rank)) allocate(this%in_this_rank(-n_ker:n_ker, -n_ker:n_ker, -n_ker:n_ker))
        end select

        call compute_kernel(this, dr)


        this%in_this_rank=.true.
        call determine_indx_and_ranks(this, dr, grid_Ndims, mpi_coords, mpi_dims)

end subroutine init_source

!###################################################################################################

subroutine update_source(this, time)

    class(TSource), intent(inout) :: this
    real(dp)     , intent(in)    :: time

    real(dp) :: envelope
    real(dp) :: cos_t

    if (time >= this%t_init .and. time <= this%t_final) then

        envelope = DEXP( -((time - this%t0)/this%tau)**2 )
        cos_t    = DCOS( this%w0*(time - this%t0) + this%phase)

        this%J_mat = this%j_amp * envelope * cos_t * this%ker_mat

        ! this%J_mat = this%j_amp * DCOS( this%w0*time)* &
        !              (aBH(1)+ &
        !              aBH(2)*DCOS(2.0*pi0*(time)/this%tau)+ &
        !              aBH(3)*DCOS(2.0*pi0*2.0*(time)/this%tau)+ &
        !              aBH(4)*DCOS(2.0*pi0*3.0*(time)/this%tau))  !Blackman-Harris window in time

    else
        this%J_mat = 0.0d0

    end if

end subroutine update_source

!###################################################################################################

subroutine kill_source(this)

    class(TSource), intent(inout) :: this

    if (allocated(this%ind_i))        deallocate(this%ind_i)
    if (allocated(this%ind_j))        deallocate(this%ind_j)
    if (allocated(this%ind_k))        deallocate(this%ind_k)
    if (allocated(this%ker_mat))      deallocate(this%ker_mat)
    if (allocated(this%J_mat))        deallocate(this%J_mat)
    if (allocated(this%in_this_rank)) deallocate(this%in_this_rank)

end subroutine kill_source

!###################################################################################################

subroutine read_init_sources(source_list, n_src, dimensions, dr, grid_Ndims, mpi_coords, mpi_dims)
    type(TSource), allocatable, intent(out) :: source_list(:)
    integer      , intent(in)  :: n_src
    integer      , intent(in)  :: dimensions
    integer      , intent(in)  :: grid_Ndims(3)
    integer      , intent(in)  :: mpi_coords(3)
    integer      , intent(in)  :: mpi_dims(3)

    character(len=2) :: direction
    real(dp)         :: j_amp
    real(dp)         :: freq
    real(dp)         :: t0
    real(dp)         :: tau
    real(dp)         :: r0(3)
    real(dp)         :: radius
    real(dp)         :: t_init
    real(dp)         :: t_final
    real(dp)         :: phase
    real(dp)         :: dr

    character(len=20) :: src_file = "sources.in"
    integer           :: ierr, funit
    integer           :: i

    if (n_src == 0) return

    if (.not. allocated(source_list)) allocate(source_list(n_src))

    ! Check whether file exists.
    inquire (file=src_file, iostat=ierr)
    if (ierr /= 0) then
        write (*, '("Error: source input file ", A, " does not exist")') trim(src_file)
        error stop
    end if

    ! Open and read source input file.
    open (action='read', file=src_file, iostat=ierr, newunit=funit)
    if (ierr /= 0) then
        write (*, '("Error: could not open source input file ", A)') trim(src_file)
        error stop
    end if

    do i=1, n_src
        read (unit=funit, fmt=*) direction, j_amp, freq, t0, tau, &
                                 r0(1), r0(2), r0(3), radius, t_init, t_final, phase

        r0      = r0 * nm_to_au
        radius  = radius * nm_to_au
        freq    = freq * ev_to_au
        t_init  = t_init * fs_to_au
        t_final = t_final * fs_to_au
        t0      = t0 * fs_to_au
        tau     = tau * fs_to_au

        call source_list(i)%init_source(direction, dimensions, j_amp, freq, t0, tau, &
                                        r0, radius, t_init, t_final, phase, dr, &
                                        grid_Ndims, mpi_coords, mpi_dims)
    end do

    close(funit)

end subroutine read_init_sources

!###################################################################################################

subroutine kill_sources(source_list, n_src)
    type(TSource), allocatable, intent(inout) :: source_list(:)
    integer      , intent(in)  :: n_src

    integer :: i

    if (n_src == 0) return

    do i=1, n_src
        call source_list(i)%kill_source()
    end do

    if (allocated(source_list)) deallocate(source_list)

end subroutine kill_sources 

!###################################################################################################

subroutine update_sources(source_list, n_src, time)
    type(TSource), allocatable, intent(inout) :: source_list(:)
    integer      , intent(in)  :: n_src
    real(dp)     , intent(in)  :: time

    integer :: i

    if (n_src == 0) return

    do i=1, n_src
        call source_list(i)%update_source(time)
    end do

end subroutine update_sources

!###################################################################################################

subroutine compute_kernel(this, dr)

    class(TSource)  , intent(inout) :: this
    real(dp)       , intent(in)    :: dr

    integer  :: i, j, k
    integer  :: n_ker
    real(dp) :: x, y, z
    real(dp) :: r_max, r, r0
    real(dp) :: norm

    n_ker = this%n_ker

    r_max = (this%n_ker*2 + 2)*dr
    r0    = r_max / 2.0d0

    norm = 0.0d0
    select case (this%dim)
    case (1)
        do i = -n_ker, n_ker

            x = i*dr

            this%ker_mat(i,1,1) = aBH(1)+ &
                          aBH(2)*DCOS(2.0*pi0*(x-r0)/r_max)+ &
                          aBH(3)*DCOS(2.0*pi0*2.0*(x-r0)/r_max)+ &
                          aBH(4)*DCOS(2.0*pi0*3.0*(x-r0)/r_max)

            norm = norm + this%ker_mat(i,1,1)
        end do
        
    case (2)
        do j = -n_ker, n_ker
        do i = -n_ker, n_ker
                
            x = i*dr
            y = j*dr
            
            r = SQRT(x**2 + y**2)
            
            this%ker_mat(i,j,1) = aBH(1)+ &
            aBH(2)*DCOS(2.0*pi0*(r-r0)/r_max)+ &
            aBH(3)*DCOS(2.0*pi0*2.0*(r-r0)/r_max)+ &
            aBH(4)*DCOS(2.0*pi0*3.0*(r-r0)/r_max)
            
            if (this%ker_mat(i,j,1) < 0.0d0) this%ker_mat(i,j,1) = 0.0d0
            
            norm = norm + this%ker_mat(i,j,1)
        end do
        end do
    case (3)
        do k = -n_ker, n_ker
        do j = -n_ker, n_ker
        do i = -n_ker, n_ker

            x = i*dr
            y = j*dr
            z = k*dr

            r = SQRT(x**2 + y**2 + z**2)

            this%ker_mat(i,j,k) = aBH(1)+ &
                          aBH(2)*DCOS(2.0*pi0*(r-r0)/r_max)+ &
                          aBH(3)*DCOS(2.0*pi0*2.0*(r-r0)/r_max)+ &
                          aBH(4)*DCOS(2.0*pi0*3.0*(r-r0)/r_max)

            if (this%ker_mat(i,j,k) < 0.0d0) this%ker_mat(i,j,k) = 0.0d0

            norm = norm + this%ker_mat(i,j,k)
        end do
        end do
        end do
    end select

    this%ker_mat = this%ker_mat / norm

end subroutine compute_kernel

!###################################################################################################

subroutine determine_indx_and_ranks(this, dr, grid_Ndims, mpi_coords, mpi_dims)

    class(TSource) , intent(inout) :: this
    real(dp)       , intent(in)    :: dr
    integer        , intent(in)    :: grid_Ndims(3)
    integer        , intent(in)    :: mpi_coords(3)
    integer        , intent(in)    :: mpi_dims(3)

    integer :: i, j, k
    integer :: rank_x, rank_y, rank_z

    !By default, we consider that the origin is where global_i, global_j and/or global_k
    !reach the half of the global grid size.

    select case (this%dim)
    case (1)
        do i = -this%n_ker, this%n_ker
            this%ind_i(i,1,1) = i+int(this%r0(1)/dr) + int(grid_Ndims(1)/2)
        end do
    case (2)
        do j = -this%n_ker, this%n_ker
        do i = -this%n_ker, this%n_ker
            this%ind_i(i,j,1) = i+int(this%r0(1)/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
            this%ind_j(i,j,1) = j+int(this%r0(2)/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
        end do
        end do
    case (3)
        do k = -this%n_ker, this%n_ker
        do j = -this%n_ker, this%n_ker
        do i = -this%n_ker, this%n_ker
            this%ind_i(i,j,k) = i+int(this%r0(1)/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
            this%ind_j(i,j,k) = j+int(this%r0(2)/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
            this%ind_k(i,j,k) = k+int(this%r0(3)/dr) + int(grid_Ndims(3)*mpi_dims(3)/2)
        end do
        end do
        end do
    end select

#ifdef USE_MPI

    select case (this%dim)
    case (1)
    case (2)
        do j = -this%n_ker, this%n_ker
        do i = -this%n_ker, this%n_ker
            rank_x = int((this%ind_i(i,j,1)-1)/grid_Ndims(1))
            rank_y = int((this%ind_j(i,j,1)-1)/grid_Ndims(2))

            if (rank_x == mpi_coords(1) .and. rank_y == mpi_coords(2)) then
                this%in_this_rank(i,j,1) = .true.
                this%ind_i(i,j,1) = this%ind_i(i,j,1) - rank_x*grid_Ndims(1)
                this%ind_j(i,j,1) = this%ind_j(i,j,1) - rank_y*grid_Ndims(2)
            else
                this%in_this_rank(i,j,1) = .false.
            end if

        end do
        end do
    case (3)
        do k = -this%n_ker, this%n_ker
        do j = -this%n_ker, this%n_ker
        do i = -this%n_ker, this%n_ker
            rank_x = int((this%ind_i(i,j,k)-1)/grid_Ndims(1))
            rank_y = int((this%ind_j(i,j,k)-1)/grid_Ndims(2))
            rank_z = int((this%ind_k(i,j,k)-1)/grid_Ndims(3))

            if (rank_x == mpi_coords(1) .and. rank_y == mpi_coords(2) .and. rank_z == mpi_coords(3)) then
                this%in_this_rank(i,j,k) = .true.
                this%ind_i(i,j,k) = this%ind_i(i,j,k) - rank_x*grid_Ndims(1)
                this%ind_j(i,j,k) = this%ind_j(i,j,k) - rank_y*grid_Ndims(2)
                this%ind_k(i,j,k) = this%ind_k(i,j,k) - rank_z*grid_Ndims(3)
            else
                this%in_this_rank(i,j,k) = .false.
            end if

        end do
        end do
        end do
    end select

#endif

end subroutine determine_indx_and_ranks

!###################################################################################################

end module source_mod