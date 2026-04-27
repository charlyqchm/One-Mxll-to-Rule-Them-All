module rs_vec_dimensions_mod

    use constants_mod
    use rs_vec_base_mod

    implicit none

    type, extends(TRSvec) :: TRSvec_1D

        integer :: nx

        complex(dp), allocatable :: pl_x(:)
        complex(dp), allocatable :: pl_y(:)
        complex(dp), allocatable :: mi_x(:)
        complex(dp), allocatable :: mi_y(:)

    contains

        procedure :: init       => init_1Dgrid
        procedure :: kill       => kill_1Dgrid

    end type TRSvec_1D

    type, extends(TRSvec) :: TRSvec_2D

        integer :: nx
        integer :: ny

        complex(dp), allocatable :: pl_x(:, :)
        complex(dp), allocatable :: pl_y(:, :)
        complex(dp), allocatable :: pl_z(:, :)
        complex(dp), allocatable :: mi_x(:, :)
        complex(dp), allocatable :: mi_y(:, :)
        complex(dp), allocatable :: mi_z(:, :)

    contains

        procedure :: init       => init_2Dgrid
        procedure :: kill       => kill_2Dgrid

    end type TRSvec_2D

    type, extends(TRSvec) :: TRSvec_3D

        integer :: nx
        integer :: ny
        integer :: nz

        complex(dp), allocatable :: pl_x(:, :, :)
        complex(dp), allocatable :: pl_y(:, :, :)
        complex(dp), allocatable :: pl_z(:, :, :)
        complex(dp), allocatable :: mi_x(:, :, :)
        complex(dp), allocatable :: mi_y(:, :, :)
        complex(dp), allocatable :: mi_z(:, :, :)

    contains

        procedure :: init       => init_3Dgrid
        procedure :: kill       => kill_3Dgrid

    end type TRSvec_3D

contains

    subroutine init_1Dgrid(this, grid_Ndims, dr, dimensions, freq, n_der)

        class(TRSvec_1D), intent(inout) :: this
        integer     , intent(in)     :: grid_Ndims(3)
        integer     , intent(in)     :: dimensions
        real(dp)    , intent(in)     :: dr
        real(dp)    , intent(in)     :: freq
        integer     , intent(in)     :: n_der

        if (grid_Ndims(1) < 1) then
            write(*,*) 'Error: grid_Ndims(1) must be >= 1, got ', grid_Ndims(1)
            error stop
        end if

        this%nx    = grid_Ndims(1)
        this%n_der = n_der

        if (.not. allocated (this%pl_x))     allocate(this%pl_x((-n_der+1):(this%nx+n_der)))
        if (.not. allocated (this%pl_y))     allocate(this%pl_y((-n_der+1):(this%nx+n_der)))
        if (.not. allocated (this%mi_x))     allocate(this%mi_x((-n_der+1):(this%nx+n_der)))
        if (.not. allocated (this%mi_y))     allocate(this%mi_y((-n_der+1):(this%nx+n_der)))

        this%pl_x = Z_0
        this%pl_y = Z_0
        this%mi_x = Z_0
        this%mi_y = Z_0

        this%dimensions = dimensions
        this%freq       = freq
        this%dr         = dr

    end subroutine init_1Dgrid

    subroutine kill_1Dgrid(this)
        class(TRSvec_1D), intent(inout) :: this

        if (allocated(this%pl_x))  deallocate(this%pl_x)
        if (allocated(this%pl_y))  deallocate(this%pl_y)
        if (allocated(this%mi_x))  deallocate(this%mi_x)
        if (allocated(this%mi_y))  deallocate(this%mi_y)

    end subroutine kill_1Dgrid

    subroutine init_2Dgrid(this, grid_Ndims, dr, dimensions, freq, n_der)
        class(TRSvec_2D), intent(inout) :: this
        integer     , intent(in)     :: grid_Ndims(3)
        integer     , intent(in)     :: dimensions
        real(dp)    , intent(in)     :: dr
        real(dp)    , intent(in)     :: freq
        integer     , intent(in)     :: n_der

        if (grid_Ndims(1) < 1 .or. grid_Ndims(2) < 1) then
            write(*,*) 'Error: grid_Ndims(1:2) must be >= 1, got ', grid_Ndims(1), grid_Ndims(2)
            error stop
        end if

        this%nx = grid_Ndims(1)
        this%ny = grid_Ndims(2)
        this%n_der = n_der

        if (.not. allocated (this%pl_x))     allocate(this%pl_x((-n_der+1):(this%nx+n_der), &
                                                                (-n_der+1):(this%ny+n_der)))
        if (.not. allocated (this%pl_y))     allocate(this%pl_y((-n_der+1):(this%nx+n_der), &
                                                                (-n_der+1):(this%ny+n_der)))
        if (.not. allocated (this%pl_z))     allocate(this%pl_z((-n_der+1):(this%nx+n_der), &
                                                                (-n_der+1):(this%ny+n_der)))
        if (.not. allocated (this%mi_x))     allocate(this%mi_x((-n_der+1):(this%nx+n_der), &
                                                                (-n_der+1):(this%ny+n_der)))
        if (.not. allocated (this%mi_y))     allocate(this%mi_y((-n_der+1):(this%nx+n_der), &
                                                                (-n_der+1):(this%ny+n_der)))
        if (.not. allocated (this%mi_z))     allocate(this%mi_z((-n_der+1):(this%nx+n_der), &
                                                                (-n_der+1):(this%ny+n_der)))

        this%pl_x = Z_0
        this%pl_y = Z_0
        this%pl_z = Z_0
        this%mi_x = Z_0
        this%mi_y = Z_0
        this%mi_z = Z_0

        this%dimensions = dimensions
        this%freq       = freq
        this%dr         = dr

    end subroutine init_2Dgrid

    subroutine kill_2Dgrid(this)
        class(TRSvec_2D), intent(inout) :: this

        if (allocated(this%pl_x))  deallocate(this%pl_x)
        if (allocated(this%pl_y))  deallocate(this%pl_y)
        if (allocated(this%pl_z))  deallocate(this%pl_z)
        if (allocated(this%mi_x))  deallocate(this%mi_x)
        if (allocated(this%mi_y))  deallocate(this%mi_y)
        if (allocated(this%mi_z))  deallocate(this%mi_z)

    end subroutine kill_2Dgrid

    subroutine init_3Dgrid(this, grid_Ndims, dr, dimensions, freq, n_der)
        class(TRSvec_3D), intent(inout) :: this
        integer     , intent(in)     :: grid_Ndims(3)
        integer     , intent(in)     :: dimensions
        real(dp)    , intent(in)     :: dr
        real(dp)    , intent(in)     :: freq
        integer     , intent(in)     :: n_der

        if (any(grid_Ndims(1:3) < 1)) then
            write(*,*) 'Error: grid_Ndims(1:3) must be >= 1, got ', grid_Ndims(1), grid_Ndims(2), grid_Ndims(3)
            error stop
        end if

        this%nx = grid_Ndims(1)
        this%ny = grid_Ndims(2)
        this%nz = grid_Ndims(3)
        this%n_der = n_der

        if (.not. allocated (this%pl_x))     allocate(this%pl_x((-n_der+1):(this%nx+n_der), &
                                                                (-n_der+1):(this%ny+n_der), &
                                                                (-n_der+1):(this%nz+n_der)))
        if (.not. allocated (this%pl_y))     allocate(this%pl_y((-n_der+1):(this%nx+n_der), &
                                                                (-n_der+1):(this%ny+n_der), &
                                                                (-n_der+1):(this%nz+n_der)))
        if (.not. allocated (this%pl_z))     allocate(this%pl_z((-n_der+1):(this%nx+n_der), &
                                                                (-n_der+1):(this%ny+n_der), &
                                                                (-n_der+1):(this%nz+n_der)))
        if (.not. allocated (this%mi_x))     allocate(this%mi_x((-n_der+1):(this%nx+n_der), &
                                                                (-n_der+1):(this%ny+n_der), &
                                                                (-n_der+1):(this%nz+n_der)))
        if (.not. allocated (this%mi_y))     allocate(this%mi_y((-n_der+1):(this%nx+n_der), &
                                                                (-n_der+1):(this%ny+n_der), &
                                                                (-n_der+1):(this%nz+n_der)))
        if (.not. allocated (this%mi_z))     allocate(this%mi_z((-n_der+1):(this%nx+n_der), &
                                                                (-n_der+1):(this%ny+n_der), &
                                                                (-n_der+1):(this%nz+n_der)))

        this%pl_x = Z_0
        this%pl_y = Z_0
        this%pl_z = Z_0
        this%mi_x = Z_0
        this%mi_y = Z_0
        this%mi_z = Z_0

        this%dimensions = dimensions
        this%freq       = freq
        this%dr         = dr

    end subroutine init_3Dgrid

    subroutine kill_3Dgrid(this)
        class(TRSvec_3D), intent(inout) :: this

        if (allocated(this%pl_x))  deallocate(this%pl_x)
        if (allocated(this%pl_y))  deallocate(this%pl_y)
        if (allocated(this%pl_z))  deallocate(this%pl_z)
        if (allocated(this%mi_x))  deallocate(this%mi_x)
        if (allocated(this%mi_y))  deallocate(this%mi_y)
        if (allocated(this%mi_z))  deallocate(this%mi_z)

    end subroutine kill_3Dgrid

end module rs_vec_dimensions_mod