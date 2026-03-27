module detector_mod

    use constants_mod

    implicit none
    
    type :: TDetector
        
        integer :: field
        integer :: detector_type
        integer :: dimensions
        integer :: i_min
        integer :: i_max
        integer :: j_min
        integer :: j_max
        integer :: k_min
        integer :: k_max
        integer :: nd
        logical :: detect_rank

        integer, allocatable :: indx_list(:,:)

        character(len=2) :: f_ch

        contains
            procedure :: init_detector, kill_detector

    end type TDetector

contains

!###################################################################################################

subroutine init_detector(this, field_ch, detector_ch, dimensions, &
                         grid_Ndims, dr, mpi_dims, mpi_coords, &
                         x_min, x_max, y_min, y_max, z_min, z_max)

    class(TDetector), intent(inout) :: this

    character(len=4) , intent(in)   :: field_ch
    character(len=10), intent(in)   :: detector_ch
    integer          , intent(in)   :: dimensions
    integer          , intent(in)   :: grid_Ndims(3)
    integer          , intent(in)   :: mpi_dims(3)
    integer          , intent(in)   :: mpi_coords(3)
    real(dp)         , intent(in)   :: x_min
    real(dp)         , intent(in)   :: x_max
    real(dp)         , intent(in)   :: y_min
    real(dp)         , intent(in)   :: y_max
    real(dp)         , intent(in)   :: z_min
    real(dp)         , intent(in)   :: z_max
    real(dp)         , intent(in)   :: dr

    integer :: i, ii
    integer :: j, jj
    integer :: k, kk
    integer :: n
    integer :: nd
    integer :: nx, ny, nz
    integer :: nx_tot, ny_tot, nz_tot
    integer :: i_loc_min, j_loc_min, k_loc_min
    integer :: i_loc_max, j_loc_max, k_loc_max

    nx_tot = grid_Ndims(1)*mpi_dims(1)
    ny_tot = grid_Ndims(2)*mpi_dims(2)
    nz_tot = grid_Ndims(3)*mpi_dims(3)

    nx = grid_Ndims(1)
    ny = grid_Ndims(2)
    nz = grid_Ndims(3)

    i_loc_min = mpi_coords(1)*nx + 1
    j_loc_min = mpi_coords(2)*ny + 1
    k_loc_min = mpi_coords(3)*nz + 1

    i_loc_max = mpi_coords(1)*nx + nx 
    j_loc_max = mpi_coords(2)*ny + ny
    k_loc_max = mpi_coords(3)*nz + nz

    nd = 0

    this%dimensions = dimensions

    this%detect_rank = .false.

    select case (field_ch)
        case ("Ex", "ex", "EX")
            this%field = Ex_FIELD
            this%f_ch = "Ex"
        case ("Ey", "ey", "EY")
            this%field = Ey_FIELD
            this%f_ch = "Ey"
            if (this%dimensions == 1) then
                write (*, '("Error: Ey field cannot be detected in 1D")')
                write (*, '("Use Ex or Hy fields instead")')
                error stop
            end if
        case ("Ez", "ez", "EZ")
            this%field = Ez_FIELD
            this%f_ch = "Ez"
            if (this%dimensions == 1) then
                write (*, '("Error: Ez field cannot be detected in 1D")')
                write (*, '("Use Ex or Hy fields instead")')
                error stop
            end if
        case ("Hx", "hx", "HX")
            this%field = Hx_FIELD
            this%f_ch = "Hx"
            if (this%dimensions == 1) then
                write (*, '("Error: Hx field cannot be detected in 1D")')
                write (*, '("Use Ex or Hy fields instead")')
                error stop
            end if
        case ("Hy", "hy", "HY")
            this%field = Hy_FIELD
            this%f_ch = "Hy"
        case ("Hz", "hz", "HZ")
            this%field = Hz_FIELD
            this%f_ch = "Hz"
            if (this%dimensions == 1) then
                write (*, '("Error: Hz field cannot be detected in 1D")')
                write (*, '("Use Ex or Hy fields instead")')
                error stop
            end if
        case default
            write (*, '("Error: invalid field type for detector")')
            error stop
    end select

    select case (detector_ch)
        case ("point", "POINT", "Point")
            this%detector_type = POINT_DETECTOR
        case ("line_x", "LINE_X", "Line_x")
            this%detector_type = LINE_X_DETECTOR
            if (this%dimensions == 1) then
                write (*, '("Error: line_x detector cannot be used in 1D")')
                write (*, '("Use point or line_z detectors instead")')
                error stop
            end if
        case ("line_y", "LINE_Y", "Line_y")
            this%detector_type = LINE_Y_DETECTOR
            if (this%dimensions == 1) then
                write (*, '("Error: line_y detector cannot be used in 1D")')
                write (*, '("Use point or line_z detectors instead")')
                error stop
            end if
        case ("line_z", "LINE_Z", "Line_z")
            this%detector_type = LINE_Z_DETECTOR
        case ("plane_xy", "PLANE_XY", "Plane_xy")
            this%detector_type = PLANE_XY_DETECTOR
            if (this%dimensions == 1) then
                write (*, '("Error: plane_xy detector cannot be used in 1D")')
                write (*, '("Use point or line_z detectors instead")')
                error stop
            end if
        case ("plane_yz", "PLANE_YZ", "Plane_yz")
            this%detector_type = PLANE_YZ_DETECTOR
            if (this%dimensions == 1) then
                write (*, '("Error: plane_yz detector cannot be used in 1D")')
                write (*, '("Use point or line_z detectors instead")')
                error stop
            end if

            if (this%dimensions == 2) then
                write (*, '("Error: plane_yz detector cannot be used in 2D")')
                write (*, '("Use plane_xy detectors instead")')
                error stop
            end if

        case ("plane_zx", "PLANE_ZX", "Plane_zx")
            this%detector_type = PLANE_ZX_DETECTOR

            if (this%dimensions == 1) then
                write (*, '("Error: plane_zx detector cannot be used in 1D")')
                write (*, '("Use point or line_z detectors instead")')
                error stop
            end if

            if (this%dimensions == 2) then
                write (*, '("Error: plane_zx detector cannot be used in 2D")')
                write (*, '("Use plane_xy detectors instead")')
                error stop
            end if

        case ("volume", "VOLUME", "Volume")
            this%detector_type = VOLUME_DETECTOR

            if (this%dimensions == 1) then
                write (*, '("Error: volume detector cannot be used in 1D")')
                write (*, '("Use point or line_z detectors instead")')
                error stop
            end if 

             if (this%dimensions == 2) then
                write (*, '("Error: volume detector cannot be used in 2D")')
                write (*, '("Use point, line, or plane_xy detectors instead")')
                error stop
            end if

        case default
            write (*, '("Error: invalid detector type")')
            error stop
    end select


    select case (this%dimensions)
    case (1)
        select case (this%detector_type)
        case (POINT_DETECTOR)
            !Some change of notations since in 1D the z-axis is the only axis,
            !and just nx is declared.
            this%k_min = int(z_min/dr) + int(nx_tot/2)
            nd = 1
            if (.not. allocated(this%indx_list)) allocate(this%indx_list(nd, dimensions))
            this%indx_list(1,1) = this%k_min
        case (LINE_Z_DETECTOR)
            this%k_min = int(z_min/dr) + int(nx_tot/2)
            this%k_max = int(z_max/dr) + int(nx_tot/2)
            do ii = 1, nx_tot
                if (ii >= this%k_min .and. ii <= this%k_max) then
                    nd = nd + 1
                end if
            end do

            if (.not. allocated(this%indx_list)) allocate(this%indx_list(nd, dimensions))

            n = 1
            do ii = 1, nx_tot
                if (ii >= this%k_min .and. ii <= this%k_max) then
                    this%indx_list(n, 1) = ii
                    n = n + 1
                end if
            end do

        case default
            write (*, '("Error: invalid combination of detector type and dimensions")')
            error stop
        end select
    case (2)
        select case (this%detector_type)
        case (POINT_DETECTOR)
            this%i_min = int(x_min/dr) + int(nx_tot/2)
            this%j_min = int(y_min/dr) + int(ny_tot/2)

            ii = this%i_min
            jj = this%j_min
            if (ii >= i_loc_min .and. ii <= i_loc_max .and. &
                jj >= j_loc_min .and. jj <= j_loc_max) then
                nd = 1
                this%detect_rank = .true.
            end if

            if (.not. allocated(this%indx_list) .and. this%detect_rank) &
            allocate(this%indx_list(nd, dimensions))

            
            if (this%detect_rank) then
                ii = this%i_min - i_loc_min + 1
                jj = this%j_min - j_loc_min + 1
                this%indx_list(1,1) = ii
                this%indx_list(1,2) = jj
            end if

        case (LINE_X_DETECTOR)
            this%i_min = int(x_min/dr) + int(nx_tot/2)
            this%i_max = int(x_max/dr) + int(nx_tot/2)
            this%j_min = int(y_min/dr) + int(ny_tot/2)

            jj = this%j_min
            if (jj >= j_loc_min .and. jj <= j_loc_max) then
                do i = 1, nx
                    ii = i_loc_min - 1 + i
                    if (ii >= this%i_min .and. ii <= this%i_max) then
                        nd = nd + 1
                        this%detect_rank = .true.
                    end if
                end do
            end if

            if (.not. allocated(this%indx_list) .and. this%detect_rank) &
            allocate(this%indx_list(nd, dimensions))

            if (this%detect_rank) then
                jj = this%j_min - j_loc_min + 1
                n = 1
                do i = 1, nx
                    ii = i_loc_min - 1 + i
                    if (ii >= this%i_min .and. ii <= this%i_max) then
                        ii  = ii - i_loc_min + 1
                        this%indx_list(n,1) = ii
                        this%indx_list(n,2) = jj
                        n = n + 1
                    end if
                end do
            end if

        case (LINE_Y_DETECTOR)
            this%j_min = int(y_min/dr) + int(ny_tot/2)
            this%j_max = int(y_max/dr) + int(ny_tot/2)
            this%i_min = int(x_min/dr) + int(nx_tot/2)

            ii = this%i_min
            if (ii >= i_loc_min .and. ii <= i_loc_max) then
                do j = 1, ny
                    jj = j_loc_min - 1 + j
                    if (jj >= this%j_min .and. jj <= this%j_max) then
                        nd = nd + 1
                        this%detect_rank = .true.
                    end if
                end do
            end if

            if (.not. allocated(this%indx_list) .and. this%detect_rank) &
            allocate(this%indx_list(nd, dimensions))

            if (this%detect_rank) then
                ii = this%i_min - i_loc_min + 1
                n = 1
                do j = 1, ny
                    jj = j_loc_min - 1 + j
                    if (jj >= this%j_min .and. jj <= this%j_max) then
                        jj  = jj - j_loc_min + 1
                        this%indx_list(n,1) = ii
                        this%indx_list(n,2) = jj
                        n = n + 1
                    end if
                end do
            end if

        case (PLANE_XY_DETECTOR)
            this%i_min = int(x_min/dr) + int(nx_tot/2)
            this%i_max = int(x_max/dr) + int(nx_tot/2)
            this%j_min = int(y_min/dr) + int(ny_tot/2)
            this%j_max = int(y_max/dr) + int(ny_tot/2)

            do j = 1, ny
            do i = 1, nx
                ii = i_loc_min - 1 + i
                jj = j_loc_min - 1 + j
                if (ii >= this%i_min .and. ii <= this%i_max .and. &
                    jj >= this%j_min .and. jj <= this%j_max) then
                    nd = nd + 1
                    this%detect_rank = .true.
                end if
            end do
            end do

            if (.not. allocated(this%indx_list) .and. this%detect_rank) &
            allocate(this%indx_list(nd, dimensions))

            if (this%detect_rank) then
                n = 1
                do j = 1, ny
                do i = 1, nx
                    ii = i_loc_min - 1 + i
                    jj = j_loc_min - 1 + j
                    if (ii >= this%i_min .and. ii <= this%i_max .and. &
                        jj >= this%j_min .and. jj <= this%j_max) then
                        ii  = ii - i_loc_min + 1
                        jj  = jj - j_loc_min + 1
                        this%indx_list(n,1) = ii
                        this%indx_list(n,2) = jj
                        n = n + 1
                    end if
                end do
                end do
            end if

        case default
            write (*, '("Error: invalid combination of detector type and dimensions")')
            error stop
        end select
    case (3)
        select case (this%detector_type)
        case (POINT_DETECTOR)
            this%i_min = int(x_min/dr) + int(nx_tot/2)
            this%j_min = int(y_min/dr) + int(ny_tot/2)
            this%k_min = int(z_min/dr) + int(nz_tot/2)

            ii = this%i_min
            jj = this%j_min
            kk = this%k_min
            if (ii >= i_loc_min .and. ii <= i_loc_max .and. &
                jj >= j_loc_min .and. jj <= j_loc_max .and. &
                kk >= k_loc_min .and. kk <= k_loc_max) then
                nd = 1
                this%detect_rank = .true.
            end if

            if (.not. allocated(this%indx_list) .and. this%detect_rank) &
            allocate(this%indx_list(nd, dimensions))

            if (this%detect_rank) then
                ii = this%i_min - i_loc_min + 1
                jj = this%j_min - j_loc_min + 1
                kk = this%k_min - k_loc_min + 1
                this%indx_list(1,1) = ii
                this%indx_list(1,2) = jj
                this%indx_list(1,3) = kk
            end if

        case (LINE_X_DETECTOR)
            this%i_min = int(x_min/dr) + int(nx_tot/2)
            this%i_max = int(x_max/dr) + int(nx_tot/2)
            this%j_min = int(y_min/dr) + int(ny_tot/2)
            this%k_min = int(z_min/dr) + int(nz_tot/2)

            jj = this%j_min
            kk = this%k_min

            if (jj >= j_loc_min .and. jj <= j_loc_max .and. &
                kk >= k_loc_min .and. kk <= k_loc_max) then
                do i = 1, nx
                    ii = i_loc_min - 1 + i
                    if (ii >= this%i_min .and. ii <= this%i_max) then
                        nd = nd + 1
                        this%detect_rank = .true.
                    end if
                end do
            end if

            if (.not. allocated(this%indx_list) .and. this%detect_rank) &
            allocate(this%indx_list(nd, dimensions))

            if (this%detect_rank) then
                jj = this%j_min - j_loc_min + 1
                kk = this%k_min - k_loc_min + 1
                n = 1
                do i = 1, nx
                    ii = i_loc_min - 1 + i
                    if (ii >= this%i_min .and. ii <= this%i_max) then
                        ii  = ii - i_loc_min + 1
                        this%indx_list(n,1) = ii
                        this%indx_list(n,2) = jj
                        this%indx_list(n,3) = kk
                        n = n + 1
                    end if
                end do
            end if

        case (LINE_Y_DETECTOR)
            this%j_min = int(y_min/dr) + int(ny_tot/2)
            this%j_max = int(y_max/dr) + int(ny_tot/2)
            this%k_min = int(z_min/dr) + int(nz_tot/2)
            this%i_min = int(x_min/dr) + int(nx_tot/2)

            ii = this%i_min
            kk = this%k_min
            if (ii >= i_loc_min .and. ii <= i_loc_max .and. &
                kk >= k_loc_min .and. kk <= k_loc_max) then
                do j = 1, ny
                    jj = j_loc_min - 1 + j
                    if (jj >= this%j_min .and. jj <= this%j_max) then
                        nd = nd + 1
                        this%detect_rank = .true.
                    end if
                end do
            end if

            if (.not. allocated(this%indx_list) .and. this%detect_rank) &
            allocate(this%indx_list(nd, dimensions))

            if (this%detect_rank) then
                ii = this%i_min - i_loc_min + 1
                kk = this%k_min - k_loc_min + 1
                n = 1
                do j = 1, ny
                    jj = j_loc_min - 1 + j
                    if (jj >= this%j_min .and. jj <= this%j_max) then
                        jj  = jj - j_loc_min + 1
                        this%indx_list(n,1) = ii
                        this%indx_list(n,2) = jj
                        this%indx_list(n,3) = kk
                        n = n + 1
                    end if
                end do
            end if

        case (LINE_Z_DETECTOR)
            this%k_min = int(z_min/dr) + int(nz_tot/2)
            this%k_max = int(z_max/dr) + int(nz_tot/2)
            this%i_min = int(x_min/dr) + int(nx_tot/2)
            this%j_min = int(y_min/dr) + int(ny_tot/2)

            ii = this%i_min
            jj = this%j_min
            if (ii >= i_loc_min .and. ii <= i_loc_max .and. &
                jj >= j_loc_min .and. jj <= j_loc_max) then
                do k = 1, nz
                    kk = k_loc_min - 1 + k
                    if (kk >= this%k_min .and. kk <= this%k_max) then
                        nd = nd + 1
                        this%detect_rank = .true.
                    end if
                end do
            end if

            if (.not. allocated(this%indx_list) .and. this%detect_rank) &
            allocate(this%indx_list(nd, dimensions))

            if (this%detect_rank) then
                ii = this%i_min - i_loc_min + 1
                jj = this%j_min - j_loc_min + 1
                n = 1
                do k = 1, nz
                    kk = k_loc_min - 1 + k
                    if (kk >= this%k_min .and. kk <= this%k_max) then
                        kk  = kk - k_loc_min + 1
                        this%indx_list(n,1) = ii
                        this%indx_list(n,2) = jj
                        this%indx_list(n,3) = kk
                        n = n + 1
                    end if
                end do
            end if

        case (PLANE_XY_DETECTOR)
            this%i_min = int(x_min/dr) + int(nx_tot/2)
            this%i_max = int(x_max/dr) + int(nx_tot/2)
            this%j_min = int(y_min/dr) + int(ny_tot/2)
            this%j_max = int(y_max/dr) + int(ny_tot/2)
            this%k_min = int(z_min/dr) + int(nz_tot/2)

            kk = this%k_min
            if (kk >= k_loc_min .and. kk <= k_loc_max) then
                do j = 1, ny
                do i = 1, nx
                    ii = i_loc_min - 1 + i
                    jj = j_loc_min - 1 + j
                    if (ii >= this%i_min .and. ii <= this%i_max .and. &
                        jj >= this%j_min .and. jj <= this%j_max) then
                        nd = nd + 1
                        this%detect_rank = .true.
                    end if
                end do
                end do
            end if

            if (.not. allocated(this%indx_list) .and. this%detect_rank) &
            allocate(this%indx_list(nd, dimensions))

            if (this%detect_rank) then
                kk = this%k_min - k_loc_min + 1
                n = 1
                do j = 1, ny
                do i = 1, nx
                    ii = i_loc_min - 1 + i
                    jj = j_loc_min - 1 + j
                    if (ii >= this%i_min .and. ii <= this%i_max .and. &
                        jj >= this%j_min .and. jj <= this%j_max) then
                        ii  = ii - i_loc_min + 1
                        jj  = jj - j_loc_min + 1
                        this%indx_list(n,1) = ii
                        this%indx_list(n,2) = jj
                        this%indx_list(n,3) = kk
                        n = n + 1
                    end if
                end do
                end do
            end if

        case (PLANE_YZ_DETECTOR)
            this%j_min = int(y_min/dr) + int(ny_tot/2)
            this%j_max = int(y_max/dr) + int(ny_tot/2)
            this%k_min = int(z_min/dr) + int(nz_tot/2)
            this%k_max = int(z_max/dr) + int(nz_tot/2)
            this%i_min = int(x_min/dr) + int(nx_tot/2)

            ii = this%i_min
            if (ii >= i_loc_min .and. ii <= i_loc_max) then
                do k = 1, nz
                do j = 1, ny
                    jj = j_loc_min - 1 + j
                    kk = k_loc_min - 1 + k
                    if (jj >= this%j_min .and. jj <= this%j_max .and. &
                        kk >= this%k_min .and. kk <= this%k_max) then
                        nd = nd + 1
                        this%detect_rank = .true.
                    end if
                end do
                end do
            end if

            if (.not. allocated(this%indx_list) .and. this%detect_rank) &
            allocate(this%indx_list(nd, dimensions))

            if (this%detect_rank) then
                ii = this%i_min - i_loc_min + 1
                n = 1
                do k = 1, nz
                do j = 1, ny
                    jj = j_loc_min - 1 + j
                    kk = k_loc_min - 1 + k
                    if (jj >= this%j_min .and. jj <= this%j_max .and. &
                        kk >= this%k_min .and. kk <= this%k_max) then
                        jj  = jj - j_loc_min + 1
                        kk  = kk - k_loc_min + 1
                        this%indx_list(n,1) = ii
                        this%indx_list(n,2) = jj
                        this%indx_list(n,3) = kk
                        n = n + 1
                    end if
                end do
                end do
            end if

        case (PLANE_ZX_DETECTOR)
            this%k_min = int(z_min/dr) + int(nz_tot/2)
            this%k_max = int(z_max/dr) + int(nz_tot/2)
            this%i_min = int(x_min/dr) + int(nx_tot/2)
            this%i_max = int(x_max/dr) + int(nx_tot/2)
            this%j_min = int(y_min/dr) + int(ny_tot/2)

            jj = this%j_min
            if (jj >= j_loc_min .and. jj <= j_loc_max) then
                do k = 1, nz
                do i = 1, nx
                    ii = i_loc_min - 1 + i
                    kk = k_loc_min - 1 + k
                    if (ii >= this%i_min .and. ii <= this%i_max .and. &
                        kk >= this%k_min .and. kk <= this%k_max) then
                        nd = nd + 1
                        this%detect_rank = .true.
                    end if
                end do
                end do
            end if

            if (.not. allocated(this%indx_list) .and. this%detect_rank) &
            allocate(this%indx_list(nd, dimensions))

            if (this%detect_rank) then
                jj = this%j_min - j_loc_min + 1
                n = 1
                do k = 1, nz
                do i = 1, nx
                    ii = i_loc_min - 1 + i
                    kk = k_loc_min - 1 + k
                    if (ii >= this%i_min .and. ii <= this%i_max .and. &
                        kk >= this%k_min .and. kk <= this%k_max) then
                        ii  = ii - i_loc_min + 1
                        kk  = kk - k_loc_min + 1
                        this%indx_list(n,1) = ii
                        this%indx_list(n,2) = jj
                        this%indx_list(n,3) = kk
                        n = n + 1
                    end if
                end do
                end do
            end if

        case (VOLUME_DETECTOR)
            this%i_min = int(x_min/dr) + int(nx_tot/2)
            this%i_max = int(x_max/dr) + int(nx_tot/2)
            this%j_min = int(y_min/dr) + int(ny_tot/2)
            this%j_max = int(y_max/dr) + int(ny_tot/2)
            this%k_min = int(z_min/dr) + int(nz_tot/2)
            this%k_max = int(z_max/dr) + int(nz_tot/2)

            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                ii = i_loc_min - 1 + i
                jj = j_loc_min - 1 + j
                kk = k_loc_min - 1 + k
                if (ii >= this%i_min .and. ii <= this%i_max .and. &
                    jj >= this%j_min .and. jj <= this%j_max .and. &
                    kk >= this%k_min .and. kk <= this%k_max) then
                    nd = nd + 1
                    this%detect_rank = .true.
                end if
            end do
            end do
            end do

            if (.not. allocated(this%indx_list) .and. this%detect_rank) &
            allocate(this%indx_list(nd, dimensions))

            if (this%detect_rank) then
                n = 1
                do k = 1, nz
                do j = 1, ny
                do i = 1, nx
                    ii = i_loc_min - 1 + i
                    jj = j_loc_min - 1 + j
                    kk = k_loc_min - 1 + k
                    if (ii >= this%i_min .and. ii <= this%i_max .and. &
                        jj >= this%j_min .and. jj <= this%j_max .and. &
                        kk >= this%k_min .and. kk <= this%k_max) then
                        ii  = ii - i_loc_min + 1
                        jj  = jj - j_loc_min + 1
                        kk  = kk - k_loc_min + 1
                        this%indx_list(n,1) = ii
                        this%indx_list(n,2) = jj
                        this%indx_list(n,3) = kk
                        n = n + 1
                    end if
                end do
                end do
                end do
            end if

        case default
            write (*, '("Error: invalid combination of detector type and dimensions")')
            error stop
        end select
    end select

    this%nd = nd

end subroutine init_detector

!###################################################################################################

subroutine kill_detector(this)

    class(TDetector), intent(inout) :: this

    if (allocated(this%indx_list)) deallocate(this%indx_list)

end subroutine kill_detector

end module detector_mod