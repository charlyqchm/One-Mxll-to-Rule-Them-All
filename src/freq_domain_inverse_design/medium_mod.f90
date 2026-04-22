module medium_mod

    use constants_mod
    
    implicit none

    type :: TMedium_eps_r
    
        integer :: dimensions

        complex(dp), allocatable :: mat1D(:)
        complex(dp), allocatable :: mat2D(:,:)
        complex(dp), allocatable :: mat3D(:,:,:)

    contains
        procedure :: init_medium
        procedure :: kill_medium
        procedure :: read_medium

    end type TMedium_eps_r

contains

!###################################################################################################

    subroutine init_medium(this, dimensions, grid_Ndims)
        
        class(TMedium_eps_r), intent(inout) :: this
        integer             , intent(in)    :: dimensions
        integer             , intent(in)    :: grid_Ndims(3)

        integer :: nx, ny, nz

        this%dimensions = dimensions

        select case (dimensions)
        case (1)
            nz = grid_Ndims(1)
            if (.not. allocated(this%mat1D)) allocate(this%mat1D(nz))
            this%mat1D = Z_ONE
        case (2)
            nx = grid_Ndims(1)
            ny = grid_Ndims(2)
            if (.not. allocated(this%mat2D)) allocate(this%mat2D(nx, ny))
            this%mat2D = Z_ONE
        case (3)
            nx = grid_Ndims(1)
            ny = grid_Ndims(2)
            nz = grid_Ndims(3)
            if (.not. allocated(this%mat3D)) allocate(this%mat3D(nx, ny, nz))
            this%mat3D = Z_ONE
        end select
        
    end subroutine init_medium

!###################################################################################################

    subroutine kill_medium(this)
        class(TMedium_eps_r), intent(inout) :: this

        if (allocated(this%mat1D))     deallocate(this%mat1D)
        if (allocated(this%mat2D))     deallocate(this%mat2D)
        if (allocated(this%mat3D))     deallocate(this%mat3D)

    end subroutine kill_medium

!###################################################################################################

subroutine read_medium(this, id, dr, grid_Ndims, mpi_coords, mpi_dims)
    class(TMedium_eps_r), intent(inout) :: this
    integer             , intent(in)    :: id
    real(dp)            , intent(in)    :: dr
    integer             , intent(in)    :: grid_Ndims(3)
    integer             , intent(in)    :: mpi_coords(3)
    integer             , intent(in)    :: mpi_dims(3)

    
    character(len=20) :: file_name = "medium_"
    character(len=20) :: file_exten = ".in"
    character(len=20) :: file_number
    character(len=20) :: input_name
    integer           :: ierr, funit
    integer           :: i, j, k
    integer           :: rank_x, rank_y, rank_z
    real(dp)          :: x, y, z
    real(dp)          :: eps_Re, eps_Im


    write(file_number, '(I3.3)') id

    input_name = trim(file_name) // trim(file_number) // trim(file_exten)

    inquire (file=input_name, iostat=ierr)
    
    if (ierr /= 0) then
        write (*, '("Warning: medium file ", A, " does not exist")') input_name
        write(*, '("No medium is concidered in problem ", I0)') id

    else

        open (action='read', file=input_name, iostat=ierr, newunit=funit)
        if (ierr /= 0) then
            write (*, '("Error: cannot open medium file ", A)') input_name
            error stop
        end if

        select case (this%dimensions)
        case (1)
            do
                read (funit, *, iostat=ierr) z, eps_Re, eps_Im
                if (ierr /= 0) exit
                i = int(z/dr) + int(grid_Ndims(1)/2)
                this%mat1D(i) = eps_Re*Z_ONE +  eps_Im*Z_I

            end do
        case (2)

            do
                read (funit, *, iostat=ierr) x, y, eps_Re, eps_Im
                if (ierr /= 0) exit
                i = int(x/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
                j = int(y/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
                rank_x = int((i-1)/grid_Ndims(1))
                rank_y = int((j-1)/grid_Ndims(2))
                if (rank_x == mpi_coords(1) .and. rank_y == mpi_coords(2)) then
                    i = i - rank_x*grid_Ndims(1)
                    j = j - rank_y*grid_Ndims(2)
                    this%mat2D(i, j) = eps_Re*Z_ONE +  eps_Im*Z_I
                end if
            end do

        case (3)

            do
                read (funit, *, iostat=ierr) x, y, z, eps_Re, eps_Im
                if (ierr /= 0) exit
                i = int(x/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
                j = int(y/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
                k = int(z/dr) + int(grid_Ndims(3)*mpi_dims(3)/2)
                rank_x = int((i-1)/grid_Ndims(1))
                rank_y = int((j-1)/grid_Ndims(2))
                rank_z = int((k-1)/grid_Ndims(3))
                if (rank_x == mpi_coords(1) .and. rank_y == mpi_coords(2) .and. &
                    rank_z == mpi_coords(3)) then
                    i = i - rank_x*grid_Ndims(1)
                    j = j - rank_y*grid_Ndims(2)
                    k = k - rank_z*grid_Ndims(3)
                    this%mat3D(i, j, k) = eps_Re*Z_ONE +  eps_Im*Z_I
                end if 
            end do

        end select

        close (funit)

    end if

end subroutine read_medium

!###################################################################################################

end module medium_mod