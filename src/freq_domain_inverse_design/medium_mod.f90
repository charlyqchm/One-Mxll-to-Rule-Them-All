module medium_mod

    use constants_mod
    
    implicit none

    type :: TMedium_eps_r
    
        integer :: dimensions

        complex(dp), allocatable :: mat1D(:)
        complex(dp), allocatable :: mat2D(:,:)
        complex(dp), allocatable :: mat3D(:,:,:)

        complex(dp), allocatable :: old_mat1D(:)
        complex(dp), allocatable :: old_mat2D(:,:)
        complex(dp), allocatable :: old_mat3D(:,:,:)

    contains
        procedure :: init_medium
        procedure :: kill_medium

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
            if (.not. allocated(this%old_mat1D)) allocate(this%old_mat1D(nz))
            this%mat1D = Z_ONE
            this%old_mat1D = Z_ONE
        case (2)
            nx = grid_Ndims(1)
            ny = grid_Ndims(2)
            if (.not. allocated(this%mat2D)) allocate(this%mat2D(nx, ny))
            if (.not. allocated(this%old_mat2D)) allocate(this%old_mat2D(nx, ny))
            this%mat2D = Z_ONE
            this%old_mat2D = Z_ONE
        case (3)
            nx = grid_Ndims(1)
            ny = grid_Ndims(2)
            nz = grid_Ndims(3)
            if (.not. allocated(this%mat3D)) allocate(this%mat3D(nx, ny, nz))
            if (.not. allocated(this%old_mat3D)) allocate(this%old_mat3D(nx, ny, nz))
            this%mat3D = Z_ONE
            this%old_mat3D = Z_ONE  
        end select
        
    end subroutine init_medium

!###################################################################################################

    subroutine kill_medium(this)
        class(TMedium_eps_r), intent(inout) :: this

        if (allocated(this%mat1D))     deallocate(this%mat1D)
        if (allocated(this%mat2D))     deallocate(this%mat2D)
        if (allocated(this%mat3D))     deallocate(this%mat3D)
        if (allocated(this%old_mat1D)) deallocate(this%old_mat1D)
        if (allocated(this%old_mat2D)) deallocate(this%old_mat2D)
        if (allocated(this%old_mat3D)) deallocate(this%old_mat3D)

    end subroutine kill_medium

!###################################################################################################
end module medium_mod