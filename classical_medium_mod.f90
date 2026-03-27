module classical_medium_mod

    use constants_mod
    implicit none

    type TClassicalMedium

        character(len=10)  :: material_type
        real(dp)           :: A1   !Temporally this variables will be equal    
        real(dp)           :: A2   ! for all the cases. They should be indpendent
        real(dp)           :: C1   ! to cover parametrical changes
        real(dp)           :: C2 
        real(dp)           :: C3
        real(dp)           :: C4
        real(dp)           :: eps_r
        integer            :: medium_type
        integer            :: nx,ny

        real(dp), allocatable :: AP_i(:)
        real(dp), allocatable :: omegaP_i(:)
        real(dp), allocatable :: GammaP_i(:)
        real(dp), allocatable :: alpha_k(:)
        real(dp), allocatable :: beta_k(:)
        real(dp), allocatable :: gamma_k(:)
        real(dp), allocatable :: B1_k(:)
        real(dp), allocatable :: B2_k(:)

            
        integer          :: n_poles
        double precision :: omegaD
        double precision :: GammaD
        double precision :: omegaP

        contains
            procedure :: init_classical_medium, kill_classical_medium

    end type TClassicalMedium

    interface read_init_media
        module procedure read_init_media_1D, read_init_media_2D, read_init_media_3D
    end interface read_init_media

contains

!###################################################################################################
    subroutine init_classical_medium(this, medium_type_ch, material_type, omegaD, GammaD, &
                                     eps_r, dt)

        class(TClassicalMedium), intent(inout) :: this
        character(len=20)      , intent(in)    :: medium_type_ch
        
        character(len=10)      , intent(in), optional :: material_type
        real(dp)               , intent(in), optional :: omegaD
        real(dp)               , intent(in), optional :: GammaD
        real(dp)               , intent(in), optional :: eps_r
        real(dp)               , intent(in), optional :: dt

        integer  :: i

        select case(trim(medium_type_ch))

        case ("dielectric")
            this%medium_type = DIELECTRIC_MEDIUM
        case ("drude")
            this%medium_type = DRUDE_MEDIUM
            this%omegaD      = omegaD
            this%GammaD      = GammaD
            this%eps_r       = eps_r

            this%A1 = (2.0d0-this%GammaD*dt)/(2.0+this%GammaD*dt)
            this%A2 = eps0*this%omegaD**2*dt/(2.0+this%GammaD*dt)
            this%C1 = (this%eps_r*eps0/dt-0.5*this%A2)/(this%eps_r*eps0/dt+0.5*this%A2)
            this%C3 = 1.0d0/(this%eps_r*eps0/dt + 0.5d0*this%A2)
            this%C4 = 0.5*(this%A1+1.0)/(this%eps_r*eps0/dt+0.5*this%A2)

        case ("drude-lorentz")    
            this%medium_type   = DL_MEDIUM
            this%material_type = material_type
            this%eps_r         = eps_r

            call get_DL_library_data(this)

            do i =1, this%n_poles
                this%alpha_k(i) = (2.0-this%omegaP_i(i)**2*dt**2)/(1.0+0.5*this%GammaP_i(i)*dt)
                this%beta_k(i)  = (this%GammaP_i(i)*dt-2.0)/(this%GammaP_i(i)*dt+2.0)
                this%gamma_k(i) = eps0*this%AP_i(i)*dt/(this%GammaP_i(i)*dt+2.0)
            end do

            this%A1 = (2.0-this%GammaD*dt)/(2.0+this%GammaD*dt)
            this%A2 = eps0*this%omegaD**2*dt/(2.0+this%GammaD*dt)

            this%C1=(this%eps_r*eps0/dt-0.5*this%A2)/(this%eps_r*eps0/dt+0.5*this%A2+0.5*SUM(this%gamma_k))
            this%C2=0.5*SUM(this%gamma_k)/(this%eps_r*eps0/dt+0.5*this%A2+0.5*SUM(this%gamma_k))
            this%C3=1.0/(this%eps_r*eps0/dt+0.5*this%A2+0.5*SUM(this%gamma_k))
            this%C4=0.5*(this%A1+1.0)/(this%eps_r*eps0/dt+0.5*this%A2+0.5*SUM(this%gamma_k))
            
            do i = 1, this%n_poles
                this%B1_k(i)=0.5*(1.0+this%alpha_k(i))/(this%eps_r*eps0/dt+0.5*this%A2+0.5*SUM(this%gamma_k))
                this%B2_k(i)=0.5*this%beta_k(i)/(this%eps_r*eps0/dt+0.5*this%A2+0.5*SUM(this%gamma_k))
            enddo

        end select

    end subroutine init_classical_medium

!###################################################################################################
 
    subroutine kill_classical_medium(this)

        class(TClassicalMedium), intent(inout) :: this

        if (allocated(this%AP_i))       deallocate(this%AP_i)
        if (allocated(this%omegaP_i))   deallocate(this%omegaP_i)
        if (allocated(this%GammaP_i))   deallocate(this%GammaP_i)
        if (allocated(this%alpha_k))    deallocate(this%alpha_k)
        if (allocated(this%beta_k))     deallocate(this%beta_k)
        if (allocated(this%gamma_k))    deallocate(this%gamma_k)
        if (allocated(this%B1_k))       deallocate(this%B1_k)
        if (allocated(this%B2_k))       deallocate(this%B2_k)

    end subroutine kill_classical_medium

!###################################################################################################

    subroutine get_DL_library_data(medium)
    
        type(TClassicalMedium), intent(inout) :: medium

        select case(medium%material_type)
        case("Ag")

            medium%n_poles = 5

            if (.not. allocated(medium%omegaP_i)) allocate(medium%omegaP_i(medium%n_poles))
            if (.not. allocated(medium%GammaP_i)) allocate(medium%GammaP_i(medium%n_poles))
            if (.not. allocated(medium%AP_i))     allocate(medium%AP_i(medium%n_poles))
            
            medium%GammaD = 0.048 * ev_to_au ; medium%omegaP = 9.01 * ev_to_au 
            medium%omegaD = (9.01*0.919238815542512) * ev_to_au  !omega_p*sqrt(f0 = 0.845)
            medium%omegaP_i(1) = 0.816 * ev_to_au ; medium%omegaP_i(2) = 4.481 * ev_to_au 
            medium%omegaP_i(3) = 8.185 * ev_to_au ; medium%omegaP_i(4) = 9.083 * ev_to_au
            medium%omegaP_i(5) = 20.29 * ev_to_au
            medium%GammaP_i(1) = 3.886 * ev_to_au ; medium%GammaP_i(2) = 0.452 * ev_to_au
            medium%GammaP_i(3) = 0.065 * ev_to_au ; medium%GammaP_i(4) = 0.916 * ev_to_au
            medium%GammaP_i(5) = 2.419* ev_to_au
            medium%AP_i(1)     = 0.065 * medium%omegaP**2 ; medium%AP_i(2)     = 0.124 * medium%omegaP**2
            medium%AP_i(3)     = 0.011 * medium%omegaP**2 ; medium%AP_i(4)     = 0.840 * medium%omegaP**2
            medium%AP_i(5)     = 5.646 * medium%omegaP**2

        case("Al")
            medium%n_poles = 4

            if (.not. allocated(medium%omegaP_i)) allocate(medium%omegaP_i(medium%n_poles))
            if (.not. allocated(medium%GammaP_i)) allocate(medium%GammaP_i(medium%n_poles))
            if (.not. allocated(medium%AP_i))     allocate(medium%AP_i(medium%n_poles))
            
            medium%GammaD = 0.047 * ev_to_au ; medium%omegaP = 14.98 * ev_to_au 
            medium%omegaD = (14.98 * 0.723187389270582) * ev_to_au  !omega_p*sqrt(f0 = 0.523)
            medium%omegaP_i(1) = 0.162 * ev_to_au ; medium%omegaP_i(2) = 1.544 * ev_to_au 
            medium%omegaP_i(3) = 1.808 * ev_to_au ; medium%omegaP_i(4) = 3.473 * ev_to_au
            medium%GammaP_i(1) = 0.333 * ev_to_au ; medium%GammaP_i(2) = 0.312 * ev_to_au
            medium%GammaP_i(3) = 1.351 * ev_to_au ; medium%GammaP_i(4) = 3.382 * ev_to_au
            medium%AP_i(1)     = 0.227 * medium%omegaP**2 ; medium%AP_i(2)     = 0.050 * medium%omegaP**2
            medium%AP_i(3)     = 0.166 * medium%omegaP**2 ; medium%AP_i(4)     = 0.030 * medium%omegaP**2
            
        case("Au")

            medium%n_poles = 5

            if (.not. allocated(medium%omegaP_i)) allocate(medium%omegaP_i(medium%n_poles))
            if (.not. allocated(medium%GammaP_i)) allocate(medium%GammaP_i(medium%n_poles))
            if (.not. allocated(medium%AP_i))     allocate(medium%AP_i(medium%n_poles))
            
            medium%GammaD = 0.053 * ev_to_au ; medium%omegaP = 9.03 * ev_to_au 
            medium%omegaD = (9.03*0.871779788708135) * ev_to_au  !omega_p*sqrt(f0 = 0.760)
            medium%omegaP_i(1) = 0.415 * ev_to_au ; medium%omegaP_i(2) = 0.830 * ev_to_au 
            medium%omegaP_i(3) = 2.969 * ev_to_au ; medium%omegaP_i(4) = 4.304 * ev_to_au
            medium%omegaP_i(5) = 13.32 * ev_to_au
            medium%GammaP_i(1) = 0.241 * ev_to_au ; medium%GammaP_i(2) = 0.345 * ev_to_au
            medium%GammaP_i(3) = 0.870 * ev_to_au ; medium%GammaP_i(4) = 2.494 * ev_to_au
            medium%GammaP_i(5) = 2.214 * ev_to_au
            medium%AP_i(1)     = 0.024 * medium%omegaP**2 ; medium%AP_i(2) = 0.010 * medium%omegaP**2
            medium%AP_i(3)     = 0.071 * medium%omegaP**2 ; medium%AP_i(4) = 0.601 * medium%omegaP**2
            medium%AP_i(5)     = 4.384 * medium%omegaP**2

        case default

            write(*, *) 'Error: The material', medium%medium_type, "is not an option for drude-lorentz medium."

        end select

        if (.not. allocated(medium%alpha_k))  allocate(medium%alpha_k(medium%n_poles))
        if (.not. allocated(medium%beta_k))   allocate(medium%beta_k(medium%n_poles))
        if (.not. allocated(medium%gamma_k))  allocate(medium%gamma_k(medium%n_poles))
        if (.not. allocated(medium%B1_k))     allocate(medium%B1_k(medium%n_poles))
        if (.not. allocated(medium%B2_k))     allocate(medium%B2_k(medium%n_poles))


    end subroutine get_DL_library_data

!###################################################################################################

    subroutine read_init_media_1D(media, n_media, eps_x, grid_Ndims, dr, media_map, nz, dt)

        type(TClassicalMedium), allocatable, intent(inout) :: media(:)
        integer               , allocatable, intent(inout) :: media_map(:)
        integer , intent(in)  :: n_media
        integer , intent(in)  :: nz
        integer , intent(in)  :: grid_Ndims(3)
        real(dp), intent(in)  :: dr
        real(dp), intent(out) :: eps_x(nz)
        real(dp), intent(in)  :: dt
        character(len=20) :: file_name = "medium_"
        character(len=20) :: file_exten = ".dat"
        character(len=20) :: file_number
        character(len=20) :: input_name
        character(len=20) :: medium_type_ch
        character(len=10) :: material_type
        integer           :: ierr, funit
        integer           :: n
        integer           :: i_ndx
        real(dp)          :: x
        real(dp)          :: eps_read
        real(dp)          :: omegaD
        real(dp)          :: GammaD

        if (.not. allocated(media_map)) allocate(media_map(nz))
        media_map = 0
        
        if (n_media == 0) return
        
        if (.not. allocated(media))     allocate(media(n_media))

        do n = 1, n_media

            write(file_number,'(I4.4)') n
            input_name = trim(file_name)//trim(file_number)//trim(file_exten)

            open (action='read', file=input_name, iostat=ierr, newunit=funit)
            if (ierr /= 0) then
                write (*, '("Error: medium file ",A," could not be opened")') &
                       input_name
                error stop
            end if

            read (unit=funit, fmt=*, iostat=ierr) medium_type_ch

             if (ierr /= 0) then
                write (*, '("Error: medium file ",A," has invalid format")') &
                       input_name
                error stop
            end if

            select case (trim(medium_type_ch))
                case ("dielectric")

                    call media(n)%init_classical_medium(medium_type_ch=medium_type_ch)
                    read (unit=funit, fmt=*) eps_read

                    do
                        read (unit=funit, fmt=*, iostat=ierr) x
                        if (ierr /= 0) exit

                        x = x*nm_to_au

                        i_ndx = int(x/dr) + int(grid_Ndims(1)/2)

                        eps_x(i_ndx) = 1.0d0/eps_read ! Store the inverse of eps_r for faster calculations

                        if (eps_read /= M_ONE) then
                            if (media_map(i_ndx)/= 0) then
                                write (*, '("Error: medium ",I4.4," overlaps medium ",I4.4)') &
                                media_map(i_ndx), n
                            else
                                media_map(i_ndx) = n
                            end if
                        end if

                    end do                    
                case ("drude")
                    read (unit=funit, fmt=*) omegaD, GammaD, eps_read

                    omegaD = omegaD*ev_to_au
                    GammaD = GammaD*ev_to_au

                    call media(n)%init_classical_medium(medium_type_ch=medium_type_ch, &
                                                        omegaD=omegaD, GammaD=GammaD,  &
                                                        eps_r=eps_read, dt=dt)
                    do
                        read (unit=funit, fmt=*, iostat=ierr) x
                        if (ierr /= 0) exit

                        x = x*nm_to_au
                        
                        i_ndx = int(x/dr) + int(grid_Ndims(1)/2)
                        
                        if (media_map(i_ndx)/= 0) then
                            write (*, '("Error: medium ",I4.4," overlaps medium ",I4.4)') &
                            media_map(i_ndx), n
                        else
                            media_map(i_ndx) = n
                        end if
                        
                    end do
                case ("drude-lorentz")
                    read (unit=funit, fmt=*) material_type, eps_read

                    call media(n)%init_classical_medium(medium_type_ch=medium_type_ch, &
                                                        material_type=material_type,   &
                                                        eps_r=eps_read, dt=dt)

                    do 
                        read (unit=funit, fmt=*, iostat=ierr) x
                        if (ierr /= 0) exit

                        x = x*nm_to_au

                        i_ndx = int(x/dr) + int(grid_Ndims(1)/2)

                        if (media_map(i_ndx)/= 0) then
                            write (*, '("Error: medium ",I4.4," overlaps medium ",I4.4)') &
                            media_map(i_ndx), n
                        else
                            media_map(i_ndx) = n
                        end if
                    end do
                case default
                    write (*, '("Error: medium type ",A," not recognized in",A)') &
                           trim(medium_type_ch), input_name
                    error stop
            end select

            close (funit)

        end do

    end subroutine read_init_media_1D

!###################################################################################################   

subroutine read_init_media_2D(media, n_media, eps_x, eps_y, eps_z, grid_Ndims, dr, media_map, &
                              nx, ny, dt, mpi_coords, mpi_dims)

        type(TClassicalMedium), allocatable, intent(inout) :: media(:)
        integer               , allocatable, intent(inout) :: media_map(:, :, :)
        integer , intent(in)  :: n_media
        integer , intent(in)  :: nx
        integer , intent(in)  :: ny
        integer , intent(in)  :: grid_Ndims(3)
        integer , intent(in)  :: mpi_coords(3) !must be equal to 0, when mpi is not used.
        integer , intent(in)  :: mpi_dims(3)   !must be equal to 1, when mpi is not used.
        real(dp), intent(in)  :: dt
        real(dp), intent(in)  :: dr
        real(dp), allocatable, intent(inout) :: eps_x(:,:) !if allocated the dimensions are (nx, ny)
        real(dp), allocatable, intent(inout) :: eps_y(:,:) !if allocated the dimensions are (nx, ny) 
        real(dp), allocatable, intent(inout) :: eps_z(:,:) !if allocated the dimensions are (nx, ny) 

        character(len=20) :: file_name = "medium_"
        character(len=20) :: file_exten = ".in"
        character(len=20) :: file_number
        character(len=20) :: input_name
        character(len=20) :: medium_type_ch
        character(len=10) :: material_type
        integer           :: ierr, funit
        integer           :: i,j,n
        integer           :: i_ndx, j_ndx
        integer           :: i_min, i_max, j_min, j_max
        integer           :: ii, jj
        integer           :: nx_tot, ny_tot
        integer           :: one_zero_i, one_zero_j, one_zero_k
        real(dp)          :: eps_x_read, eps_y_read, eps_z_read, eps_read
        real(dp)          :: x, y
        real(dp)          :: omegaD
        real(dp)          :: GammaD

        integer, allocatable :: media_map_aux(:, :, :)

        if (.not. allocated(media_map)) allocate(media_map(nx, ny, 3))
        media_map = 0

        allocate(media_map_aux(nx*mpi_dims(1), ny*mpi_dims(2), 3))

        media_map_aux = 0

        if (n_media == 0) return
        if (.not. allocated(media))     allocate(media(n_media))


        i_min = mpi_coords(1)*nx + 1
        i_max = (mpi_coords(1)+1)*nx
        j_min = mpi_coords(2)*ny + 1
        j_max = (mpi_coords(2)+1)*ny

        nx_tot = mpi_dims(1)*nx
        ny_tot = mpi_dims(2)*ny

        do n = 1, n_media

            write(file_number,'(I4.4)') n
            input_name = trim(file_name)//trim(file_number)//trim(file_exten)

            open (action='read', file=input_name, iostat=ierr, newunit=funit)
            if (ierr /= 0) then
                write (*, '("Error: medium file ",A," could not be opened")') &
                       input_name
                error stop
            end if

            read (unit=funit, fmt=*, iostat=ierr) medium_type_ch

             if (ierr /= 0) then
                write (*, '("Error: medium file ",A," has invalid format")') &
                       input_name
                error stop
            end if

            select case (trim(medium_type_ch))
                case ("dielectric")

                    read (unit=funit, fmt=*) eps_read
                    call media(n)%init_classical_medium(medium_type_ch=medium_type_ch,eps_r=eps_read)

                case ("drude")
                    read (unit=funit, fmt=*) omegaD, GammaD, eps_read

                    omegaD = omegaD*ev_to_au
                    GammaD = GammaD*ev_to_au

                    call media(n)%init_classical_medium(medium_type_ch=medium_type_ch, &
                                                        omegaD=omegaD, GammaD=GammaD, &
                                                        eps_r=eps_read, dt=dt)

                case ("drude-lorentz")
                    read (unit=funit, fmt=*) material_type, eps_read

                    call media(n)%init_classical_medium(medium_type_ch=medium_type_ch, &
                                                        material_type=material_type, &
                                                        eps_r=eps_read, dt=dt)

                case default
                    write (*, '("Error: medium type ",A," not recognized in",A)') &
                           trim(medium_type_ch), input_name
                    error stop
            end select

            do

                read (unit=funit, fmt=*, iostat=ierr) x, y
                if (ierr /= 0) exit
                
                x = x*nm_to_au
                y = y*nm_to_au

                ii = int(x/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
                jj = int(y/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)

                if (media_map_aux(ii,jj,1)/= 0) then
                    write (*, '("Error: medium ",I4.4," overlaps medium ",I4.4)') &
                    media_map_aux(ii,jj,1), n
                else
                    media_map_aux(ii,jj,1) = n
                end if

                if (media_map_aux(ii,jj,2)/= 0) then
                    write (*, '("Error: medium ",I4.4," overlaps medium ",I4.4)') &
                    media_map_aux(ii,jj,2), n
                else
                    media_map_aux(ii,jj,2) = n
                end if

                if (media_map_aux(ii,jj,3)/= 0) then
                    write (*, '("Error: medium ",I4.4," overlaps medium ",I4.4)') &
                    media_map_aux(ii,jj,3), n
                else
                    media_map_aux(ii,jj,3) = n
                end if

            end do

            close (funit)

        end do

        do n=1, n_media
        do j=1, ny_tot
        do i=1, nx_tot
            if (i >= i_min .and. i <= i_max .and. &
                j >= j_min .and. j <= j_max) then

                i_ndx = i - i_min + 1
                j_ndx = j - j_min + 1

                if (i == nx_tot .and. media_map_aux(i,j,1) == n) then
                
                    media_map(i_ndx,j_ndx,1) = n
                
                else if (j == ny_tot .and. media_map_aux(i,j,1) == n) then
                
                    media_map(i_ndx,j_ndx,1) = n
                
                else

                    if (media_map_aux(i,j,1) == n .and. media_map_aux(i+1,j,1) == n) then
                        media_map(i_ndx,j_ndx,1) = n
                    else if (media_map_aux(i,j,1) /= 0 .and. media_map_aux(i+1,j,1) == n) then
                        media_map(i_ndx,j_ndx,1) = n
                    end if

                    if (media_map_aux(i,j,2) == n .and. media_map_aux(i,j+1,2) == n) then
                        media_map(i_ndx,j_ndx,2) = n
                    else if (media_map_aux(i,j,2) /= 0 .and. media_map_aux(i,j+1,2) == n) then
                        media_map(i_ndx,j_ndx,2) = n
                    end if

                    if (media_map_aux(i,j,3) == n) then
                        media_map(i_ndx,j_ndx,3) = n
                    end if

                end if

            end if
        end do
        end do
        end do
        
        do n=1, n_media
            if (media(n)%medium_type == DIELECTRIC_MEDIUM) then
                do j=1, ny
                do i=1, nx
                    if (media_map(i,j,1) == n ) then
                        eps_x(i,j) = 1.0d0/media(n)%eps_r
                        eps_y(i,j) = 1.0d0/media(n)%eps_r
                        eps_z(i,j) = 1.0d0/media(n)%eps_r
                    end if
                end do
                end do
            end if
        end do

        deallocate(media_map_aux)

    end subroutine read_init_media_2D

!###################################################################################################

    subroutine read_init_media_3D(media, n_media, eps_x, eps_y, eps_z, grid_Ndims, dr , media_map, &
                                  nx, ny, nz, dt, mpi_coords, mpi_dims)

        type(TClassicalMedium), allocatable, intent(inout) :: media(:)
        integer               , allocatable, intent(inout) :: media_map(:, :, :, :)
        integer , intent(in)  :: n_media
        integer , intent(in)  :: nx
        integer , intent(in)  :: ny
        integer , intent(in)  :: nz
        integer , intent(in)  :: mpi_coords(3) !must be equal to 0, when mpi is not used.
        integer , intent(in)  :: mpi_dims(3)   !must be equal to 1, when mpi is not used.
        integer , intent(in)  :: grid_Ndims(3)
        real(dp), intent(in)  :: dr
        real(dp), intent(in)  :: dt
        real(dp), intent(out) :: eps_x(nx, ny, nz)
        real(dp), intent(out) :: eps_y(nx, ny, nz)
        real(dp), intent(out) :: eps_z(nx, ny, nz)
        
        character(len=20) :: file_name = "medium_"
        character(len=20) :: file_exten = ".dat"
        character(len=20) :: file_number
        character(len=20) :: input_name
        character(len=20) :: medium_type_ch
        character(len=10) :: material_type
        integer           :: ierr, funit
        integer           :: i,j,k,n
        integer           :: i_ndx, j_ndx, k_ndx
        integer           :: i_min, i_max, j_min, j_max, k_min, k_max
        integer           :: ii, jj, kk
        integer           :: nx_tot, ny_tot, nz_tot
        integer           :: one_zero_i, one_zero_j, one_zero_k
        real(dp)          :: eps_x_read, eps_y_read, eps_z_read, eps_read
        real(dp)          :: omegaD
        real(dp)          :: GammaD
        real(dp)          :: x, y, z

        integer, allocatable :: media_map_aux(:, :, :, :)

        if (.not. allocated(media_map)) allocate(media_map(nx, ny, nz, 3))
        media_map = 0
        
        if (n_media == 0) return
        if (.not. allocated(media))     allocate(media(n_media))
        allocate(media_map_aux(nx*mpi_dims(1), ny*mpi_dims(2), nz*mpi_dims(3), 3))
        media_map_aux = 0

        i_min = mpi_coords(1)*nx + 1
        i_max = (mpi_coords(1)+1)*nx
        j_min = mpi_coords(2)*ny + 1
        j_max = (mpi_coords(2)+1)*ny
        k_min = mpi_coords(3)*nz + 1
        k_max = (mpi_coords(3)+1)*nz

        nx_tot = mpi_dims(1)*nx
        ny_tot = mpi_dims(2)*ny
        nz_tot = mpi_dims(3)*nz

        do n = 1, n_media

            write(file_number,'(I4.4)') n
            input_name = trim(file_name)//trim(file_number)//trim(file_exten)

            open (action='read', file=input_name, iostat=ierr, newunit=funit)
            if (ierr /= 0) then
                write (*, '("Error: medium file ",A," could not be opened")') &
                       input_name
                error stop
            end if

            read (unit=funit, fmt=*, iostat=ierr) medium_type_ch

             if (ierr /= 0) then
                write (*, '("Error: medium file ",A," has invalid format")') &
                       input_name
                error stop
            end if

            select case (trim(medium_type_ch))
                case ("dielectric")

                    read (unit=funit, fmt=*) eps_read

                    call media(n)%init_classical_medium(medium_type_ch=medium_type_ch)

                case ("drude")
                    read (unit=funit, fmt=*) omegaD, GammaD, eps_read

                    omegaD = omegaD*ev_to_au
                    GammaD = GammaD*ev_to_au

                    call media(n)%init_classical_medium(medium_type_ch=medium_type_ch, &
                                                        omegaD=omegaD, GammaD=GammaD, &
                                                        eps_r=eps_read, dt=dt)

                case ("drude-lorentz")
                    read (unit=funit, fmt=*) material_type, eps_read

                    call media(n)%init_classical_medium(medium_type_ch=medium_type_ch, &
                                                        material_type=material_type, &
                                                        eps_r=eps_read, dt=dt)

                case default
                    write (*, '("Error: medium type ",A," not recognized in",A)') &
                           trim(medium_type_ch), input_name
                    error stop
            end select
        

            do

                read (unit=funit, fmt=*, iostat=ierr) x, y, z
                if (ierr /= 0) exit
                
                x = x*nm_to_au
                y = y*nm_to_au
                z = z*nm_to_au

                ii = int(x/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
                jj = int(y/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
                kk = int(z/dr) + int(grid_Ndims(3)*mpi_dims(3)/2)

                if (media_map_aux(ii,jj,kk,1)/= 0) then
                    write (*, '("Error: medium ",I4.4," overlaps medium ",I4.4)') &
                    media_map_aux(ii,jj,kk,1), n
                else
                    media_map_aux(ii,jj,kk,1) = n
                end if

                if (media_map_aux(ii,jj,kk,2)/= 0) then
                    write (*, '("Error: medium ",I4.4," overlaps medium ",I4.4)') &
                    media_map_aux(ii,jj,kk,2), n
                else
                    media_map_aux(ii,jj,kk,2) = n
                end if

                if (media_map_aux(ii,jj,kk,3)/= 0) then
                    write (*, '("Error: medium ",I4.4," overlaps medium ",I4.4)') &
                    media_map_aux(ii,jj,kk,3), n
                else
                    media_map_aux(ii,jj,kk,3) = n
                end if

            end do

        end do

        close (funit)

        do n=1, n_media
        do k=1, nz_tot
        do j=1, ny_tot
        do i=1, nx_tot
            if (i >= i_min .and. i <= i_max .and. &
                j >= j_min .and. j <= j_max .and. &
                k >= k_min .and. k <= k_max) then

                i_ndx = i - i_min + 1
                j_ndx = j - j_min + 1
                k_ndx = k - k_min + 1

                if (i == nx_tot .and. media_map_aux(i,j,k,1) == n) then
                
                    media_map(i_ndx,j_ndx,k_ndx,1) = n
                
                else if (j == ny_tot .and. media_map_aux(i,j,k,1) == n) then
                
                    media_map(i_ndx,j_ndx,k_ndx,1) = n
                
                else if (k == nz_tot .and. media_map_aux(i,j,k,1) == n) then
                
                    media_map(i_ndx,j_ndx,k_ndx,1) = n
                
                else


                    if (media_map_aux(i,j,k,1) == n .and. media_map_aux(i+1,j,k,1) == n) then
                        media_map(i_ndx,j_ndx, k_ndx,1) = n
                    else if (media_map_aux(i,j,k,1) /= 0 .and. media_map_aux(i+1,j,k,1) == n) then
                        media_map(i_ndx,j_ndx, k_ndx,1) = n
                    end if

                    if (media_map_aux(i,j,k,2) == n .and. media_map_aux(i,j+1,k,2) == n) then
                        media_map(i_ndx,j_ndx, k_ndx,2) = n
                    else if (media_map_aux(i,j,k,2) /= 0 .and. media_map_aux(i,j+1,k,2) == n) then
                        media_map(i_ndx,j_ndx, k_ndx,2) = n
                    end if

                    if (media_map_aux(i,j,k,3) == n .and. media_map_aux(i,j,k+1,3) == n) then
                        media_map(i_ndx,j_ndx, k_ndx,3) = n
                    else if (media_map_aux(i,j,k,3) /= 0 .and. media_map_aux(i,j,k+1,3) == n) then
                        media_map(i_ndx,j_ndx, k_ndx,3) = n
                    end if

                end if

            end if
        end do
        end do
        end do
        end do

        do n=1, n_media
            if (media(n)%medium_type == DIELECTRIC_MEDIUM) then
                do k=1, nz
                do j=1, ny
                do i=1, nx
                    if (media_map(i,j,k,1) == n ) then
                        eps_x(i,j,k) = 1.0d0/media(n)%eps_r
                    end if
                    if (media_map(i,j,k,2) == n ) then
                        eps_y(i,j,k) = 1.0d0/media(n)%eps_r
                    end if
                    if (media_map(i,j,k,3) == n ) then
                        eps_z(i,j,k) = 1.0d0/media(n)%eps_r
                    end if
                end do
                end do
                end do
            end if
        end do

        deallocate(media_map_aux)

    end subroutine read_init_media_3D

!###################################################################################################

    subroutine get_medium_polarization(media, idx, PDi, PLi, PLi_old, &
                                        Ei, Ei_old, rotH, no_medium)

        !If there is no medium, media is not allocated and idx=0 always.
        !TO-DO: improve this condition.
        type(TClassicalMedium), allocatable, intent(in)    :: media(:)
        integer               , intent(in)    :: idx
        logical               , intent(inout) :: no_medium
        real(dp)              , intent(inout) :: PDi
        real(dp)              , intent(inout) :: PLi(:)
        real(dp)              , intent(inout) :: PLi_old(:)
        real(dp)              , intent(inout) :: Ei
        real(dp)              , intent(inout) :: Ei_old
        real(dp)              , intent(in)    :: rotH
        
        integer  :: i
        real(dp) :: tmpE
        real(dp) :: tmpPL(size(PLi))

        if (idx==0) then
            no_medium = .true.
            return
        end if

        select case (media(idx)%medium_type)

        case(DRUDE_MEDIUM)
            tmpE = media(idx)%C1 * Ei + media(idx)%C3 * rotH - media(idx)%C4 * PDi
            PDi  = media(idx)%A1 * PDi + media(idx)%A2 * (tmpE + Ei)
            Ei   = tmpE

        case(DL_MEDIUM)
            tmpE = media(idx)%C1 * Ei + media(idx)%C2 * Ei_old + media(idx)%C3 * rotH &
                   - media(idx)%C4 * PDi

            do i =1, media(idx)%n_poles
                tmpE = tmpE - media(idx)%B1_k(i) * PLi(i) - media(idx)%B2_k(i) * PLi_old(i)
            end do

            PDi  = media(idx)%A1 * PDi + media(idx)%A2 * (tmpE + Ei)

            do i= 1, media(idx)%n_poles
                tmpPL(i) = media(idx)%alpha_k(i) * PLi(i) + &
                           media(idx)%beta_k(i)  * PLi_old(i) + &
                           media(idx)%gamma_k(i) * (tmpE - Ei_old)
            end do

            Ei_old = Ei
            
            do i=1, media(idx)%n_poles
                PLi_old(i) = PLi(i)
                PLi(i)     = tmpPL(i)
            end do
            
            Ei     = tmpE

        case default
            no_medium = .true. !Dielectric case
        end select

    end subroutine get_medium_polarization

!###################################################################################################

end module classical_medium_mod