module outputs_mod

    use constants_mod
    use mxll_base_mod
    use mxll_1D_mod
    use mxll_2D_mod
    use mxll_3D_mod
    use detector_mod
#ifdef USE_MPI
    use mpi
#endif
    implicit none

contains

!###################################################################################################

subroutine init_outputs(n_detectors, t_det_print, dt_det_print, detectors, dimensions, &
                        grid_Ndims, mxll_mode, dr, dt, mpi_dims, mpi_coords, myrank)

    type(TDetector), allocatable, intent(inout) :: detectors(:)
    integer , intent(out) :: t_det_print
    integer , intent(in)  :: n_detectors
    integer , intent(in)  :: dimensions
    integer , intent(in)  :: grid_Ndims(3)
    integer , intent(in)  :: mxll_mode
    integer , intent(in)  :: mpi_dims(3)
    integer , intent(in)  :: mpi_coords(3)
    integer , intent(in)  :: myrank
    real(dp), intent(in)  :: dt_det_print
    real(dp), intent(in)  :: dr
    real(dp), intent(in)  :: dt

    integer           :: ierr, funit
    integer           :: i
    real(dp)          :: x_min, x_max, y_min, y_max, z_min, z_max
    real(dp)          :: x, y, z
    character(len=4)  :: field_name
    character(len=10) :: detector_type
    character(len=20) :: input_name = "detectors.in"
    character(len=20) :: dirname = "./detector"
    character(len=99) :: full_dirname
    character(len=20) :: number

    if (n_detectors == 0) return

    if (.not. allocated(detectors)) allocate(detectors(n_detectors))

    t_det_print = int(dt_det_print/dt)

    if (t_det_print <= 0) then
        write(*, *) "Error: dt_det_print must be greater than Maxwell dt."
        error stop
    end if

    open (action='read', file=input_name, iostat=ierr, newunit=funit)

    if (ierr /= 0) then
        write(*, *) "Error: Could not open file ", trim(input_name)
        stop
    end if

    do i = 1, n_detectors
        read(funit, *) detector_type, field_name, x_min, x_max, y_min, y_max, z_min, z_max

        x_min = x_min*nm_to_au
        x_max = x_max*nm_to_au
        y_min = y_min*nm_to_au
        y_max = y_max*nm_to_au
        z_min = z_min*nm_to_au
        z_max = z_max*nm_to_au


        call detectors(i)%init_detector(field_name, detector_type, dimensions, &
                           grid_Ndims, dr, mpi_dims, mpi_coords, &
                           x_min, x_max, y_min, y_max, z_min, z_max)
    end do

    close(funit)

!Checking the 2D modes. These are not checked in the initialization of the detectors.

    do i = 1, n_detectors
        if (mxll_mode == TMZ_2D_MODE .and. detectors(i)%field == Ey_FIELD) then
            write(*, *) "Error: Ey field cannot be detected in TMZ 2D mode."
            stop
        else if (mxll_mode == TMZ_2D_MODE .and. detectors(i)%field == Ex_FIELD) then
            write(*, *) "Error: Ex field cannot be detected in TMZ 2D mode."
            stop
        else if (mxll_mode == TMZ_2D_MODE .and. detectors(i)%field == Hz_FIELD) then
            write(*, *) "Error: Hz field cannot be detected in TMZ 2D mode."
            stop
        else if (mxll_mode == TEZ_2D_MODE .and. detectors(i)%field == Ez_FIELD) then
            write(*, *) "Error: Ez field cannot be detected in TEZ 2D mode."
            stop
        else if (mxll_mode == TEZ_2D_MODE .and. detectors(i)%field == Hx_FIELD) then
            write(*, *) "Error: Hx field cannot be detected in TEZ 2D mode."
            stop
        else if (mxll_mode == TEZ_2D_MODE .and. detectors(i)%field == Hy_FIELD) then
            write(*, *) "Error: Hy field cannot be detected in TEZ 2D mode."
            stop
        end if
    end do

! Just the rank 0 will create directories and files for the detectors, 
! to avoid conflicts in parallel runs.
! This header is not included in another subroutine because 
! it is needed to be created only once.
    if (myrank == 0) then
        do i=1, n_detectors
            write(number, '(I7.7)') i
            full_dirname = trim(dirname) // "_" // trim(number)

            call execute_command_line("mkdir -p " // trim(full_dirname), exitstat=ierr)
        
            if (ierr /= 0) then
                write(*, *) "Error: Could not create directory ", trim(full_dirname)
                stop
            end if

            if (detectors(i)%detector_type == POINT_DETECTOR) then         
                open(newunit=funit, file=trim(full_dirname) // "/point.dat", &
                     status='replace', action='write', iostat=ierr)

                if (ierr /= 0) then
                    write(*, *) "Error: Could not open file for writing: ", &
                            trim(full_dirname) // "/point.dat"
                    stop
                end if
                
                if (dimensions == 1 ) then
                    x = 0.0
                    y = 0.0
                    z = detectors(i)%k_min*dr - int(grid_Ndims(1)/2)*dr
                else
                    x = detectors(i)%i_min*dr - int(mpi_dims(1)*grid_Ndims(1)/2)*dr
                    y = detectors(i)%j_min*dr - int(mpi_dims(2)*grid_Ndims(2)/2)*dr
                    z = detectors(i)%k_min*dr - int(mpi_dims(3)*grid_Ndims(3)/2)*dr
                end if

                x = x*au_to_nm
                y = y*au_to_nm
                z = z*au_to_nm

                select case (detectors(i)%field)
                case (Ex_FIELD)
                    write(funit, '("# Time (a.u.)              Ex (a.u.),      &
                    x = ", F10.5, " nm,   y = ", F10.5, " nm,   z = ", F10.5, " nm")') x, y, z
                case (Ey_FIELD)
                    write(funit, '("# Time (a.u.)              Ey (a.u.),      &
                    x = ", F10.5, " nm,   y = ", F10.5, " nm,   z = ", F10.5, " nm")') x, y, z
                case (Ez_FIELD)
                    write(funit, '("# Time (a.u.)              Ez (a.u.),      &
                    x = ", F10.5, " nm,   y = ", F10.5, " nm,   z = ", F10.5, " nm")') x, y, z
                case (Hx_FIELD)
                    write(funit, '("# Time (a.u.)              Hx (a.u.),      &
                    x = ", F10.5, " nm,   y = ", F10.5, " nm,   z = ", F10.5, " nm")') x, y, z
                case (Hy_FIELD)
                    write(funit, '("# Time (a.u.)              Hy (a.u.),      &
                    x = ", F10.5, " nm,   y = ", F10.5, " nm,   z = ", F10.5, " nm")') x, y, z
                case (Hz_FIELD)
                    write(funit, '("# Time (a.u.)              Hz (a.u.),      &
                    x = ", F10.5, " nm,   y = ", F10.5, " nm,   z = ", F10.5, " nm")') x, y, z
                end select
                close(funit)
            end if
       end do 
    end if

end subroutine init_outputs

!###################################################################################################

subroutine kill_outputs(detectors, n_detectors)

    type(TDetector), allocatable, intent(inout) :: detectors(:)
    integer , intent(in)  :: n_detectors

    integer :: i

    if (.not. allocated(detectors)) return

    do i = 1, n_detectors
        call kill_detector(detectors(i))
    end do

    deallocate(detectors)

end subroutine kill_outputs

!###################################################################################################

subroutine write_detectors_outputs(detectors, mxll, n_detectors, t_step, print_det_step, &
                                   t_det_print, time, dr, grid_Ndims, mpi_dims, mpi_coords, myrank)

    type(TDetector), intent(in)    :: detectors(n_detectors)
    class(TMxll)   , intent(in)    :: mxll
    integer        , intent(in)    :: n_detectors
    integer        , intent(inout) :: print_det_step
    integer        , intent(in)    :: t_step
    integer        , intent(in)    :: myrank
    integer        , intent(in)    :: grid_Ndims(3)
    integer        , intent(in)    :: mpi_dims(3)
    integer        , intent(in)    :: mpi_coords(3)
    integer        , intent(in)    :: t_det_print
    real(dp)       , intent(in)    :: dr
    real(dp)       , intent(in)    :: time

    if (MOD(t_step, t_det_print) /= 0) return
    if (n_detectors == 0) return

    select type (mxll)
    class is (TMxll_1D)
        call write_1D_header(detectors, n_detectors, time, print_det_step, dr, grid_Ndims, mpi_dims, myrank) 
        call write_1D_field(detectors, mxll, n_detectors, print_det_step, time, dr, &
                            grid_Ndims, mpi_dims, mpi_coords, myrank)
    class is (TMxll_2D)

        call write_2D_headers(detectors, n_detectors, print_det_step, time, dr, grid_Ndims, mpi_dims, myrank)
        call write_2D_field(detectors, mxll, n_detectors, print_det_step, time, dr, &
                            grid_Ndims, mpi_dims, mpi_coords, myrank)

    class is (TMxll_3D)

        call write_3D_headers(detectors, n_detectors, print_det_step, time, dr, grid_Ndims, mpi_dims, myrank)
        call write_3D_field(detectors, mxll, n_detectors, print_det_step, time, dr, &
                            grid_Ndims, mpi_dims, mpi_coords, myrank)

    end select

    print_det_step = print_det_step + 1

end subroutine write_detectors_outputs
!###################################################################################################
subroutine write_1D_header(detectors, n_detectors, time, print_det_step, &
                            dr, grid_Ndims, mpi_dims, myrank)

    type(TDetector) , intent(in) :: detectors(n_detectors)
    integer         , intent(in) :: n_detectors
    integer         , intent(in) :: print_det_step
    integer         , intent(in) :: myrank
    integer         , intent(in) :: grid_Ndims(3)
    integer         , intent(in) :: mpi_dims(3)
    real(dp)        , intent(in) :: time
    real(dp)        , intent(in) :: dr

    integer  :: i
    real(dp) :: z
    integer  :: nz

    character(len=20) :: directory = "./detector"
    character(len=20) :: number
    character(len=20) :: print_number
    character(len=99) :: full_dirname
    integer           :: ierr, funit

    if (myrank /= 0) return

    nz = grid_Ndims(1)


    do i = 1, n_detectors

        select case (detectors(i)%detector_type)

        case (LINE_Z_DETECTOR)

            write(number, '(I7.7)') i
            write(print_number, '(I9.9)') print_det_step

            select case (detectors(i)%field)

            case (Ex_FIELD)

                full_dirname = trim(directory) // "_" // trim(number) // "/Ex_line_" // &
                                        trim(print_number) // ".dat"

                open(newunit=funit, file=trim(full_dirname), status='replace', &
                        action='write', iostat=ierr)

                write(funit, '("# Time = ", ES18.8, " (a.u.)")') time    
                write(funit, '("# z (nm)                   Ex (a.u.)")')

                close (funit)

        case (Hy_FIELD)

            full_dirname = trim(directory) // "_" // trim(number) // "/Hy_line_z_" // &
                                    trim(print_number) // ".dat"

            open(newunit=funit, file=trim(full_dirname), status='replace', &
                    action='write', iostat=ierr)

            write(funit, '("# Time = ", ES18.8, " (a.u.)")') time    
            write(funit, '("# z (nm)                   Hy (a.u.)")')

            close (funit)

        end select

    end select
    end do

end subroutine write_1D_header

!###################################################################################################
subroutine write_1D_field(detectors, mxll, n_detectors, print_det_step, time, &
                          dr, grid_Ndims, mpi_dims, mpi_coords, myrank)

    type(TDetector) , intent(in) :: detectors(n_detectors)
    class(TMxll_1D) , intent(in) :: mxll
    integer         , intent(in) :: n_detectors
    integer         , intent(in) :: print_det_step
    integer         , intent(in) :: myrank
    integer         , intent(in) :: grid_Ndims(3)
    integer         , intent(in) :: mpi_dims(3)
    integer         , intent(in) :: mpi_coords(3)
    real(dp)        , intent(in) :: dr
    real(dp)        , intent(in) :: time

    integer  :: n
    integer  :: i
    integer  :: k_ndx
    real(dp) :: z

    character(len=20) :: directory = "./detector"
    character(len=20) :: number
    character(len=20) :: print_number
    character(len=4)  :: field_name
    character(len=99) :: full_dirname
    integer           :: ierr, funit

    if (myrank /= 0) return

    do i = 1, n_detectors

        write(number, '(I7.7)') i
        write(print_number, '(I9.9)') print_det_step

        select case(detectors(i)%detector_type)

        case (POINT_DETECTOR)

            select case (detectors(i)%field)
                
            case (Ex_FIELD)

                full_dirname = trim(directory) // "_" // trim(number) // "/point" // ".dat"

                open(newunit=funit, file=trim(full_dirname), status='old', &
                    action='write', position='append', iostat=ierr)

                k_ndx = detectors(i)%indx_list(1,1)
                write(funit, *) time, mxll%Ex(k_ndx)
                close(funit)

            case (Hy_FIELD)

                full_dirname = trim(directory) // "_" // trim(number) // "/point.dat"
                open(newunit=funit, file=trim(full_dirname), status='old', &
                    action='write', position='append', iostat=ierr)

                k_ndx = detectors(i)%indx_list(1,1)
                write(funit, *) time, 0.5*(mxll%Hy(k_ndx-1)+mxll%Hy(k_ndx))
                close(funit)

            end select

        case (LINE_Z_DETECTOR)
            select case (detectors(i)%field)
            case (Ex_FIELD)
                full_dirname = trim(directory) // "_" // trim(number) // "/Ex_line_" // &
                               trim(print_number) // ".dat"

                open(newunit=funit, file=trim(full_dirname), status='replace', &
                     action='write', iostat=ierr)

                write(funit, '("# Time = ", ES18.8, " (a.u.)")') time    
                write(funit, '("# z (nm)                   Ex (a.u.)")')

                do n = 1, detectors(i)%nd
                    k_ndx = detectors(i)%indx_list(n,1)
                    z    = k_ndx*dr - int(grid_Ndims(1)/2)*dr
                    write(funit, *) z , mxll%Ex(k_ndx)
                end do

                close(funit)

            case (Hy_FIELD)
                full_dirname = trim(directory) // "_" // trim(number) // "/Hy_line_z_" // &
                               trim(print_number) // ".dat"
                open(newunit=funit, file=trim(full_dirname), status='replace', &
                     action='write', iostat=ierr)

                write(funit, '("# Time = ", ES18.8, " (a.u.)")') time    
                write(funit, '("# z (nm)                   Hy (a.u.)")')

                do n = 1, detectors(i)%nd
                    k_ndx = detectors(i)%indx_list(n,1)
                    z    = k_ndx*dr - int(grid_Ndims(1)/2)*dr
                    write(funit, *) z , 0.5*(mxll%Hy(k_ndx-1)+mxll%Hy(k_ndx))
                end do

                close(funit)
            end select
        end select

    end do
                        
end subroutine write_1D_field

!###################################################################################################

subroutine write_2D_headers(detectors, n_detectors, print_det_step, time, &
                            dr, grid_Ndims, mpi_dims, myrank)

    type(TDetector) , intent(in) :: detectors(n_detectors)
    integer         , intent(in) :: n_detectors
    integer         , intent(in) :: myrank
    integer         , intent(in) :: grid_Ndims(3)
    integer         , intent(in) :: mpi_dims(3)
    integer         , intent(in) :: print_det_step
    real(dp)        , intent(in) :: dr
    real(dp)        , intent(in) :: time

    integer  :: n
    integer  :: i
    integer  :: i_ndx, j_ndx
    integer  :: nx, ny
    integer  :: nx_tot, ny_tot
    real(dp) :: x, y
    real(dp) :: coor_val

    character(len=20) :: directory = "./detector"
    character(len=20) :: number
    character(len=20) :: print_number
    character(len=2)  :: field_name
    character(len=1)  :: coor_head
    character(len=1) :: coor_write
    character(len=99) :: full_dirname
    integer           :: ierr, funit


    if (myrank /= 0) return

    do i = 1, n_detectors

        field_name = detectors(i)%f_ch
        write(number, '(I7.7)') i
        write(print_number, '(I9.9)') print_det_step

        select case(detectors(i)%detector_type)
        case (LINE_X_DETECTOR, LINE_Y_DETECTOR)
            
            if (detectors(i)%detector_type == LINE_X_DETECTOR) then
                coor_head  = "y"
                coor_write = "x"
                y = detectors(i)%j_min*dr - int(mpi_dims(2)*grid_Ndims(2)/2)*dr
                coor_val  = y*au_to_nm
            else if (detectors(i)%detector_type == LINE_Y_DETECTOR) then
                coor_head  = "x"
                coor_write = "y"
                x = detectors(i)%i_min*dr - int(mpi_dims(1)*grid_Ndims(1)/2)*dr
                coor_val  = x*au_to_nm
            end if

            full_dirname = trim(directory) // "_" // trim(number) // &
                                "/" // trim(field_name) // "_line_" // trim(print_number) // ".dat"
    
            open(newunit=funit, file=trim(full_dirname), status='replace', &
                action='write', iostat=ierr)
            write(funit, '("# Time = ", ES18.8, " (a.u.), ' // trim(coor_head) // &
                            ' = ", F12.6, " (nm)")') time, coor_val
            write(funit, '("# ' // trim(coor_write) // ' (nm)                   ' &
                            // trim(field_name) // ' (a.u.)")')
            close(funit)

        case (PLANE_XY_DETECTOR)
            
            full_dirname = trim(directory) // "_" // trim(number) // &
                                "/" // trim(field_name) // "_plane_xy_" // trim(print_number) // ".dat"
           
            open(newunit=funit, file=trim(full_dirname), status='replace', &
                action='write', iostat=ierr)
            write(funit, '("# Time = ", ES18.8, " (a.u.)")') time
            write(funit, '("# x (nm)                   y (nm)                   '&
                           // trim(field_name) // ' (a.u.)")')
        
            close(funit)
        end select
    end do

end subroutine write_2D_headers

!###################################################################################################
subroutine write_2D_field(detectors, mxll, n_detectors, print_det_step, time, &
                          dr, grid_Ndims, mpi_dims, mpi_coords, myrank)

    type(TDetector) , intent(in) :: detectors(n_detectors)
    class(TMxll_2D) , intent(in) :: mxll
    integer         , intent(in) :: n_detectors
    integer         , intent(in) :: print_det_step
    integer         , intent(in) :: myrank
    integer         , intent(in) :: grid_Ndims(3)
    integer         , intent(in) :: mpi_dims(3)
    integer         , intent(in) :: mpi_coords(3)
    real(dp)        , intent(in) :: dr
    real(dp)        , intent(in) :: time

    integer  :: n
    integer  :: i
    integer  :: i_ndx, j_ndx
    integer  :: nx, ny
    integer  :: nx_tot, ny_tot
    real(dp) :: x, y

    character(len=20) :: directory = "./detector"
    character(len=20) :: number
    character(len=20) :: print_number
    character(len=4)  :: field_name
    character(len=99) :: full_dirname
    integer           :: ierr, funit

#ifdef USE_MPI
    integer :: rank_iter, nprocs_mpi
#endif

    nx = grid_Ndims(1)
    ny = grid_Ndims(2)

    nx_tot = grid_Ndims(1)*mpi_dims(1)
    ny_tot = grid_Ndims(2)*mpi_dims(2)

#ifdef USE_MPI

    nprocs_mpi = mpi_dims(1)*mpi_dims(2)*mpi_dims(3)
    do rank_iter = 0, nprocs_mpi-1
        if (myrank == rank_iter) then
#endif

    do i = 1, n_detectors
    
        write(number, '(I7.7)') i
        write(print_number, '(I9.9)') print_det_step

        select case(detectors(i)%detector_type)

        case (POINT_DETECTOR)

            full_dirname = trim(directory) // "_" // trim(number) // "/point.dat"
            
            open(newunit=funit, file=trim(full_dirname), status='old', &
            action='write', position='append', iostat=ierr)
            
            i_ndx = detectors(i)%indx_list(1,1)
            j_ndx = detectors(i)%indx_list(1,2)

            select case (detectors(i)%field)
            
            case (Ex_FIELD)
                
                write(funit, *) time, 0.5*(mxll%Ex(i_ndx-1,j_ndx)+mxll%Ex(i_ndx,j_ndx))

            case (Ey_FIELD)

                write(funit, *) time, 0.5*(mxll%Ey(i_ndx,j_ndx-1)+mxll%Ey(i_ndx,j_ndx))

            case (Ez_FIELD)
                
                write(funit, *) time, mxll%Ez(i_ndx,j_ndx)

            case (Hx_FIELD)
                
                write(funit, *) time, 0.5*(mxll%Hx(i_ndx,j_ndx) + mxll%Hx(i_ndx,j_ndx-1))
             
            case (Hy_FIELD)

                write(funit, *) time, 0.5*(mxll%Hy(i_ndx,j_ndx) + mxll%Hy(i_ndx-1,j_ndx))
             
            case (Hz_FIELD)

                write(funit, *) time, 0.25*(mxll%Hz(i_ndx,j_ndx) + mxll%Hz(i_ndx-1,j_ndx-1) +&
                                            mxll%Hz(i_ndx-1,j_ndx) + mxll%Hz(i_ndx,j_ndx-1))
             
            end select

            close(funit)

        case (LINE_X_DETECTOR)

            field_name = detectors(i)%f_ch
            full_dirname = trim(directory) // "_" // trim(number) // &
                                "/" // trim(field_name) // "_line_" // trim(print_number) // ".dat"
            
            open(newunit=funit, file=trim(full_dirname), status='old', &
                     action='write', position='append', iostat=ierr)

            select case (detectors(i)%field)

            case (Ex_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    x    = x*au_to_nm
                    write(funit, *) x , 0.5*(mxll%Ex(i_ndx-1,j_ndx)+mxll%Ex(i_ndx,j_ndx))
                end do

            case (Ey_FIELD)
                
                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    x    = x*au_to_nm
                    write(funit, *) x , 0.5*(mxll%Ey(i_ndx,j_ndx-1)+mxll%Ey(i_ndx,j_ndx))
                end do

            case (Ez_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    x    = x*au_to_nm
                    write(funit, *) x , mxll%Ez(i_ndx,j_ndx)
                end do
    
            case (Hx_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    x    = x*au_to_nm
                    write(funit, *) x , 0.5*(mxll%Hx(i_ndx,j_ndx) + mxll%Hx(i_ndx,j_ndx-1))
                end do

            case (Hy_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    x    = x*au_to_nm
                    write(funit, *) x , 0.5*(mxll%Hy(i_ndx,j_ndx) + mxll%Hy(i_ndx-1,j_ndx))
                end do

            case (Hz_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    x    = x*au_to_nm
                    write(funit, *) x , 0.25*(mxll%Hz(i_ndx,j_ndx) + mxll%Hz(i_ndx-1,j_ndx-1) +&
                                                  mxll%Hz(i_ndx-1,j_ndx) + mxll%Hz(i_ndx,j_ndx-1))
                end do
    
            end select
            close(funit)

        case (LINE_Y_DETECTOR)

                field_name = detectors(i)%f_ch
                full_dirname = trim(directory) // "_" // trim(number) // &
                                "/" // trim(field_name) // "_line_" // trim(print_number) // ".dat"

                open(newunit=funit, file=trim(full_dirname), status='old', &
                         action='write', position='append', iostat=ierr)

                select case (detectors(i)%field)

                case (Ex_FIELD)

                    do n = 1, detectors(i)%nd
                        i_ndx = detectors(i)%indx_list(n,1)
                        j_ndx = detectors(i)%indx_list(n,2)
                        y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                        y    = y*au_to_nm
                        write(funit, *) y , 0.5*(mxll%Ex(i_ndx-1,j_ndx)+mxll%Ex(i_ndx,j_ndx))
                    end do

                case (Ey_FIELD)

                    do n = 1, detectors(i)%nd
                        i_ndx = detectors(i)%indx_list(n,1)
                        j_ndx = detectors(i)%indx_list(n,2)
                        y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                        y    = y*au_to_nm
                        write(funit, *) y , 0.5*(mxll%Ey(i_ndx,j_ndx-1)+mxll%Ey(i_ndx,j_ndx))
                    end do

                case (Ez_FIELD)

                    do n = 1, detectors(i)%nd
                        i_ndx = detectors(i)%indx_list(n,1)
                        j_ndx = detectors(i)%indx_list(n,2)
                        y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                        y    = y*au_to_nm
                        write(funit, *) y , mxll%Ez(i_ndx,j_ndx)
                    end do

                case (Hx_FIELD)

                    do n = 1, detectors(i)%nd
                        i_ndx = detectors(i)%indx_list(n,1)
                        j_ndx = detectors(i)%indx_list(n,2)
                        y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                        y    = y*au_to_nm
                        write(funit, *) y , 0.5*(mxll%Hx(i_ndx,j_ndx) + mxll%Hx(i_ndx,j_ndx-1))
                    end do

                case (Hy_FIELD)

                    do n = 1, detectors(i)%nd
                        i_ndx = detectors(i)%indx_list(n,1)
                        j_ndx = detectors(i)%indx_list(n,2)
                        y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                        y    = y*au_to_nm
                        write(funit, *) y , 0.5*(mxll%Hy(i_ndx,j_ndx) + mxll%Hy(i_ndx-1,j_ndx))
                    end do

                case (Hz_FIELD)

                    do n = 1, detectors(i)%nd
                        i_ndx = detectors(i)%indx_list(n,1)
                        j_ndx = detectors(i)%indx_list(n,2)
                        y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                        y    = y*au_to_nm
                        write(funit, *) y , 0.25*(mxll%Hz(i_ndx,j_ndx) + mxll%Hz(i_ndx-1,j_ndx-1) +&
                                                  mxll%Hz(i_ndx-1,j_ndx) + mxll%Hz(i_ndx,j_ndx-1))
                    end do

                end select
                
                close(funit)

        case (PLANE_XY_DETECTOR)

            field_name = detectors(i)%f_ch

            full_dirname = trim(directory) // "_" // trim(number) // &
                               "/" // trim(field_name) // "_plane_xy_" // trim(print_number) // ".dat"

            open(newunit=funit, file=trim(full_dirname), status='old', &
                        action='write', position='append', iostat=ierr)

            select case (detectors(i)%field)
            case (Ex_FIELD)
                
                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    write(funit, *) x , y, 0.5*(mxll%Ex(i_ndx-1,j_ndx)+mxll%Ex(i_ndx,j_ndx))
                end do

            case (Ey_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    write(funit, *) x , y, 0.5*(mxll%Ey(i_ndx,j_ndx-1)+mxll%Ey(i_ndx,j_ndx))
                end do

            case (Ez_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    write(funit, *) x , y, mxll%Ez(i_ndx,j_ndx)
                end do

            case (Hx_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    write(funit, *) x , y, 0.5*(mxll%Hx(i_ndx,j_ndx) + mxll%Hx(i_ndx,j_ndx-1))
                end do

                close(funit)

            case (Hy_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    write(funit, *) x , y, 0.5*(mxll%Hy(i_ndx,j_ndx) + mxll%Hy(i_ndx-1,j_ndx))
                end do

            case (Hz_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    write(funit, *) x , y, 0.25*(mxll%Hz(i_ndx,j_ndx) + mxll%Hz(i_ndx-1,j_ndx-1) +&
                                                  mxll%Hz(i_ndx-1,j_ndx) + mxll%Hz(i_ndx,j_ndx-1))
                end do
                
            end select

            close(funit)

        end select

    end do

#ifdef USE_MPI
        end if
        call mpi_barrier(MPI_COMM_WORLD, ierr)
    end do
#endif

end subroutine write_2D_field
!###################################################################################################

subroutine write_3D_headers(detectors, n_detectors, print_det_step, time, &
                            dr, grid_Ndims, mpi_dims, myrank)

    type(TDetector) , intent(in) :: detectors(n_detectors)
    integer         , intent(in) :: n_detectors
    integer         , intent(in) :: myrank
    integer         , intent(in) :: grid_Ndims(3)
    integer         , intent(in) :: mpi_dims(3)
    integer         , intent(in) :: print_det_step
    real(dp)        , intent(in) :: dr
    real(dp)        , intent(in) :: time

    integer  :: n
    integer  :: i
    integer  :: i_ndx, j_ndx, k_ndx
    integer  :: nx, ny, nz
    integer  :: nx_tot, ny_tot, nz_tot
    real(dp) :: x, y, z
    real(dp) :: coor_val1, coor_val2

    character(len=20) :: directory = "./detector"
    character(len=20) :: number
    character(len=20) :: print_number
    character(len=2)  :: field_name
    character(len=1)  :: coor_head1
    character(len=1)  :: coor_head2
    character(len=1)  :: coor_write1
    character(len=1)  :: coor_write2
    character(len=99) :: full_dirname
    integer           :: ierr, funit

    if (myrank /= 0) return

    do i = 1, n_detectors

        field_name = detectors(i)%f_ch
        write(number, '(I7.7)') i
        write(print_number, '(I9.9)') print_det_step

        select case(detectors(i)%detector_type)
        case (LINE_X_DETECTOR, LINE_Y_DETECTOR, LINE_Z_DETECTOR)
            if (detectors(i)%detector_type == LINE_X_DETECTOR) then
                coor_head1 = "y" ; coor_head2 = "z"
                coor_write1 = "x"
                y = detectors(i)%j_min*dr - int(mpi_dims(2)*grid_Ndims(2)/2)*dr
                z = detectors(i)%k_min*dr - int(mpi_dims(3)*grid_Ndims(3)/2)*dr
                coor_val1  = y*au_to_nm
                coor_val2  = z*au_to_nm
            else if (detectors(i)%detector_type == LINE_Y_DETECTOR) then
                coor_head1 = "z" ; coor_head2 = "x"
                coor_write1 = "y"
                x = detectors(i)%i_min*dr - int(mpi_dims(1)*grid_Ndims(1)/2)*dr
                coor_val1  = z*au_to_nm
                coor_val2  = x*au_to_nm
            else if (detectors(i)%detector_type == LINE_Z_DETECTOR) then
                coor_head1 = "x" ; coor_head2 = "y"
                coor_write1 = "z"
                z = detectors(i)%k_min*dr - int(mpi_dims(3)*grid_Ndims(3)/2)*dr
                coor_val1  = x*au_to_nm
                coor_val2  = y*au_to_nm
            end if

            full_dirname = trim(directory) // "_" // trim(number) // &
                                "/" // trim(field_name) // "_line_" // trim(print_number) // ".dat"
    
            open(newunit=funit, file=trim(full_dirname), status='replace', &
                action='write', iostat=ierr)
            write(funit, '("# Time = ", ES18.8, " (a.u.), ' // trim(coor_head1) // &
                            ' = ", F12.6, " (nm), ' // trim(coor_head2) // &
                            ' = ", F12.6, " (nm)")') time, coor_val1, coor_val2
            write(funit, '("# ' // trim(coor_write1) // ' (nm)                   ' &
                            // trim(coor_write2) // ' (nm)                   ' &
                            // trim(field_name) // ' (a.u.)")')
            close(funit)

        case (PLANE_XY_DETECTOR, PLANE_YZ_DETECTOR, PLANE_ZX_DETECTOR)
            
            if (detectors(i)%detector_type == PLANE_XY_DETECTOR) then
                coor_head1 = "z"
                coor_write1 = "x" ; coor_write2 = "y"
            else if (detectors(i)%detector_type == PLANE_YZ_DETECTOR) then
                coor_head1 = "x"
                coor_write1 = "y" ; coor_write2 = "z"
            else if (detectors(i)%detector_type == PLANE_ZX_DETECTOR) then
                coor_head1 = "y"
                coor_write1 = "z" ; coor_write2 = "x"
            end if

            full_dirname = trim(directory) // "_" // trim(number) // &
                                "/" // trim(field_name) // "_plane_" // trim(print_number) // ".dat"
            
            open(newunit=funit, file=trim(full_dirname), status='replace', &
                action='write', iostat=ierr)
            write(funit, '("# Time = ", ES18.8, " (a.u.), ' // trim(coor_head1) // ' = ", F12.6, " (nm)")') time, coor_val1
            write(funit, '("# ' // trim(coor_write1) // ' (nm)                   ' &
                           // trim(coor_write2) // ' (nm)                   ' &
                           // trim(field_name) // ' (a.u.)")')
        
            close(funit)

        case (VOLUME_DETECTOR)

            full_dirname = trim(directory) // "_" // trim(number) // &
                                "/" // trim(field_name) // "_volume_" // trim(print_number) // ".dat"
            
            open(newunit=funit, file=trim(full_dirname), status='replace', &
                action='write', iostat=ierr)
            write(funit, '("# Time = ", ES18.8, " (a.u.)")') time
            write(funit, '("# x (nm)                   y (nm)                   z (nm)' &
                          //'                    '// trim(field_name) // ' (a.u.)")')
        
            close(funit)

        end select
    end do

end subroutine write_3D_headers

!###################################################################################################

subroutine write_3D_field(detectors, mxll, n_detectors, print_det_step, time, &
                          dr, grid_Ndims, mpi_dims, mpi_coords, myrank)

    type(TDetector) , intent(in) :: detectors(n_detectors)
    class(TMxll_3D) , intent(in) :: mxll
    integer         , intent(in) :: n_detectors
    integer         , intent(in) :: print_det_step
    integer         , intent(in) :: myrank
    integer         , intent(in) :: grid_Ndims(3)
    integer         , intent(in) :: mpi_dims(3)
    integer         , intent(in) :: mpi_coords(3)
    real(dp)        , intent(in) :: dr
    real(dp)        , intent(in) :: time

    integer  :: n
    integer  :: i
    integer  :: i_ndx, j_ndx, k_ndx
    integer  :: nx, ny, nz
    integer  :: nx_tot, ny_tot, nz_tot
    real(dp) :: x, y, z

    character(len=20) :: directory = "./detector"
    character(len=20) :: number
    character(len=20) :: print_number
    character(len=4)  :: field_name
    character(len=99) :: full_dirname
    character(len=32) :: position1, position2
    integer           :: ierr, funit

#ifdef USE_MPI
    integer :: rank_iter, nprocs_mpi
#endif

    nx = grid_Ndims(1)
    ny = grid_Ndims(2)
    nz = grid_Ndims(3)
    
    nx_tot = grid_Ndims(1)*mpi_dims(1)
    ny_tot = grid_Ndims(2)*mpi_dims(2)
    nz_tot = grid_Ndims(3)*mpi_dims(3)

#ifdef USE_MPI

    nprocs_mpi = mpi_dims(1)*mpi_dims(2)*mpi_dims(3)
    do rank_iter = 0, nprocs_mpi-1
        if (myrank == rank_iter) then
#endif


    do i = 1, n_detectors ! Main loop over detectors
   
        if (detectors(i)%detect_rank) then !---- Here it checks if the rank has a detector
        
        write(number, '(I7.7)') i
        write(print_number, '(I9.9)') print_det_step

        select case(detectors(i)%detector_type)

        case (POINT_DETECTOR)

            full_dirname = trim(directory) // "_" // trim(number) // "/point.dat"
                
            open(newunit=funit, file=trim(full_dirname), status='old', &
                action='write', position='append', iostat=ierr)

            i_ndx = detectors(i)%indx_list(1,1)
            j_ndx = detectors(i)%indx_list(1,2)
            k_ndx = detectors(i)%indx_list(1,3)

            select case(detectors(i)%field)
            case (Ex_FIELD)

                write(funit, *) time, 0.5*(mxll%Ex(i_ndx-1,j_ndx,k_ndx)+mxll%Ex(i_ndx,j_ndx,k_ndx))
             
            case (Ey_FIELD)

                write(funit, *) time, 0.5*(mxll%Ey(i_ndx,j_ndx-1,k_ndx)+mxll%Ey(i_ndx,j_ndx,k_ndx))
             
            case (Ez_FIELD)

                write(funit, *) time, mxll%Ez(i_ndx,j_ndx,k_ndx)
             
            case (Hx_FIELD)

                write(funit, *) time, 0.25*(mxll%Hx(i_ndx,j_ndx-1,k_ndx) + mxll%Hx(i_ndx,j_ndx,k_ndx) +&
                                            mxll%Hx(i_ndx,j_ndx,k_ndx-1) + mxll%Hx(i_ndx,j_ndx-1,k_ndx-1))
             
            case (Hy_FIELD)

                write(funit, *) time, 0.25*(mxll%Hy(i_ndx-1,j_ndx,k_ndx) + mxll%Hy(i_ndx,j_ndx,k_ndx) +&
                                            mxll%Hy(i_ndx,j_ndx,k_ndx-1) + mxll%Hy(i_ndx-1,j_ndx,k_ndx-1))
             
            case (Hz_FIELD)

                write(funit, *) time, 0.25*(mxll%Hz(i_ndx,j_ndx,k_ndx) + mxll%Hz(i_ndx-1,j_ndx,k_ndx) + &
                                            mxll%Hz(i_ndx,j_ndx-1,k_ndx) + mxll%Hz(i_ndx-1,j_ndx-1,k_ndx))

            end select

            close(funit)

        case (LINE_X_DETECTOR)

            field_name = detectors(i)%f_ch

            full_dirname = trim(directory) // "_" // trim(number) // &
                               "/" // trim(field_name) // "_line_" // trim(print_number) // ".dat"

            open(newunit=funit, file=trim(full_dirname), status='old', &
                action='write', position='append', iostat=ierr)

            select case(detectors(i)%field)
            case (Ex_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    x    = x*au_to_nm
                    write(funit, *) x , 0.5*(mxll%Ex(i_ndx-1,j_ndx,k_ndx)+mxll%Ex(i_ndx,j_ndx,k_ndx))
                end do

            case (Ey_FIELD)
                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    x    = x*au_to_nm
                    write(funit, *) x , 0.5*(mxll%Ey(i_ndx,j_ndx-1,k_ndx)+mxll%Ey(i_ndx,j_ndx,k_ndx))
                end do

            case (Ez_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    x    = x*au_to_nm
                    write(funit, *) x , 0.5*(mxll%Ez(i_ndx,j_ndx,k_ndx-1)+mxll%Ez(i_ndx,j_ndx,k_ndx))
                end do

            case (Hx_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    x    = x*au_to_nm
                    write(funit, *) x , 0.25*(mxll%Hx(i_ndx,j_ndx-1,k_ndx) + mxll%Hx(i_ndx,j_ndx,k_ndx) + &
                                                mxll%Hx(i_ndx,j_ndx,k_ndx-1) + mxll%Hx(i_ndx,j_ndx-1,k_ndx-1))
                end do                

            case (Hy_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    x    = x*au_to_nm
                    write(funit, *) x , 0.25*(mxll%Hy(i_ndx-1,j_ndx,k_ndx) + mxll%Hy(i_ndx,j_ndx,k_ndx) +&
                                                mxll%Hy(i_ndx,j_ndx,k_ndx-1) + mxll%Hy(i_ndx-1,j_ndx,k_ndx-1))
                end do

            case (Hz_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    x    = x*au_to_nm
                    write(funit, *) x , 0.25*(mxll%Hz(i_ndx,j_ndx,k_ndx) + mxll%Hz(i_ndx-1,j_ndx,k_ndx) + &
                                                mxll%Hz(i_ndx,j_ndx-1,k_ndx) + mxll%Hz(i_ndx-1,j_ndx-1,k_ndx))
                end do

            end select

            close(funit)

        case (LINE_Y_DETECTOR)

            field_name = detectors(i)%f_ch
            full_dirname = trim(directory) // "_" // trim(number) // &
                           "/" // trim(field_name) // "_line_" // trim(print_number) // ".dat"

            open(newunit=funit, file=trim(full_dirname), status='old', &
                     action='write', position='append', iostat=ierr)

            select case(detectors(i)%field)

            case (Ex_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    y    = y*au_to_nm
                    write(funit, *) y , 0.5*(mxll%Ex(i_ndx-1,j_ndx,k_ndx)+mxll%Ex(i_ndx,j_ndx,k_ndx))
                end do

            case (Ey_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    y    = y*au_to_nm
                    write(funit, *) y , 0.5*(mxll%Ey(i_ndx,j_ndx-1,k_ndx)+mxll%Ey(i_ndx,j_ndx,k_ndx))
                end do

            case (Ez_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    y    = y*au_to_nm
                    write(funit, *) y , 0.5*(mxll%Ez(i_ndx,j_ndx,k_ndx-1)+mxll%Ez(i_ndx,j_ndx,k_ndx))
                end do

            case (Hx_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    y    = y*au_to_nm
                    write(funit, *) y , 0.25*(mxll%Hx(i_ndx,j_ndx-1,k_ndx) + mxll%Hx(i_ndx,j_ndx,k_ndx) +&
                                                mxll%Hx(i_ndx,j_ndx,k_ndx-1) + mxll%Hx(i_ndx,j_ndx-1,k_ndx-1))
                end do

            case (Hy_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    y    = y*au_to_nm
                    write(funit, *) y , 0.25*(mxll%Hy(i_ndx-1,j_ndx,k_ndx) + mxll%Hy(i_ndx,j_ndx,k_ndx) +&
                                                mxll%Hy(i_ndx,j_ndx,k_ndx-1) + mxll%Hy(i_ndx-1,j_ndx,k_ndx-1))
                end do

            case (Hz_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    y    = y*au_to_nm
                    write(funit, *) y , 0.25*(mxll%Hz(i_ndx,j_ndx,k_ndx) + mxll%Hz(i_ndx-1,j_ndx,k_ndx) + &
                                                mxll%Hz(i_ndx,j_ndx-1,k_ndx) + mxll%Hz(i_ndx-1,j_ndx-1,k_ndx))
                end do

            end select

            close(funit)

        case (LINE_Z_DETECTOR)

            field_name = detectors(i)%f_ch
            full_dirname = trim(directory) // "_" // trim(number) // &
                           "/" // trim(field_name) // "_line_" // trim(print_number) // ".dat"

            open(newunit=funit, file=trim(full_dirname), status='old', &
                     action='write', position='append', iostat=ierr)

            select case(detectors(i)%field)

            case (Ex_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    z    = z*au_to_nm
                    write(funit, *) z , 0.5*(mxll%Ex(i_ndx-1,j_ndx,k_ndx)+mxll%Ex(i_ndx,j_ndx,k_ndx))
                end do

            case (Ey_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    z    = z*au_to_nm
                    write(funit, *) z , 0.5*(mxll%Ey(i_ndx,j_ndx-1,k_ndx)+mxll%Ey(i_ndx,j_ndx,k_ndx))
                end do

            case (Ez_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    z    = z*au_to_nm
                    write(funit, *) z , 0.5*(mxll%Ez(i_ndx,j_ndx,k_ndx-1)+mxll%Ez(i_ndx,j_ndx,k_ndx))
                end do

            case (Hx_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    z    = z*au_to_nm
                    write(funit, *) z , 0.25*(mxll%Hx(i_ndx,j_ndx-1,k_ndx) + mxll%Hx(i_ndx,j_ndx,k_ndx) + &
                                                mxll%Hx(i_ndx,j_ndx,k_ndx-1) + mxll%Hx(i_ndx,j_ndx-1,k_ndx-1))
                end do

            case (Hy_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    z    = z*au_to_nm
                    write(funit, *) z , 0.25*(mxll%Hy(i_ndx-1,j_ndx,k_ndx) + mxll%Hy(i_ndx,j_ndx,k_ndx) +&
                                                mxll%Hy(i_ndx,j_ndx,k_ndx-1) + mxll%Hy(i_ndx-1,j_ndx,k_ndx-1))
                end do

            case (Hz_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    z    = z*au_to_nm
                    write(funit, *) z , 0.25*(mxll%Hz(i_ndx,j_ndx,k_ndx) + mxll%Hz(i_ndx-1,j_ndx,k_ndx) + &
                                                mxll%Hz(i_ndx,j_ndx-1,k_ndx) + mxll%Hz(i_ndx-1,j_ndx-1,k_ndx))
                end do

            end select

            close(funit)

        case (PLANE_XY_DETECTOR)

            field_name = detectors(i)%f_ch
            full_dirname = trim(directory) // "_" // trim(number) // &
                               "/" // trim(field_name) // "_plane_" // trim(print_number) // ".dat"

            open(newunit=funit, file=trim(full_dirname), status='old', &
                action='write', position='append', iostat=ierr)

            select case(detectors(i)%field)

            case (Ex_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    write(funit, *) x , y, 0.5*(mxll%Ex(i_ndx-1,j_ndx,k_ndx)+mxll%Ex(i_ndx,j_ndx,k_ndx))
                end do

            case (Ey_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    write(funit, *) x , y, 0.5*(mxll%Ey(i_ndx,j_ndx-1,k_ndx)+mxll%Ey(i_ndx,j_ndx,k_ndx))
                end do

            case (Ez_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    write(funit, *) x , y, 0.5*(mxll%Ez(i_ndx,j_ndx,k_ndx-1)+mxll%Ez(i_ndx,j_ndx,k_ndx))
                end do

            case (Hx_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    write(funit, *) x , y, 0.25*(mxll%Hx(i_ndx,j_ndx-1,k_ndx) + mxll%Hx(i_ndx,j_ndx,k_ndx) + &
                                                mxll%Hx(i_ndx,j_ndx,k_ndx-1) + mxll%Hx(i_ndx,j_ndx-1,k_ndx-1))
                end do

            case (Hy_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    write(funit, *) x , y, 0.25*(mxll%Hy(i_ndx-1,j_ndx,k_ndx) + mxll%Hy(i_ndx,j_ndx,k_ndx) +&
                                                mxll%Hy(i_ndx,j_ndx,k_ndx-1) + mxll%Hy(i_ndx-1,j_ndx,k_ndx-1))
                end do

            case (Hz_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    write(funit, *) x , y, 0.25*(mxll%Hz(i_ndx,j_ndx,k_ndx) + mxll%Hz(i_ndx-1,j_ndx,k_ndx) + &
                                                mxll%Hz(i_ndx,j_ndx-1,k_ndx) + mxll%Hz(i_ndx-1,j_ndx-1,k_ndx))
                end do

            end select

            close(funit)

        case (PLANE_YZ_DETECTOR)

            field_name = detectors(i)%f_ch
            full_dirname = trim(directory) // "_" // trim(number) // &
                               "/" // trim(field_name) // "_plane_" // trim(print_number) // ".dat"

            open(newunit=funit, file=trim(full_dirname), status='old', &
                action='write', position='append', iostat=ierr)

            select case(detectors(i)%field)

            case (Ex_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    y    = y*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) y , z, 0.5*(mxll%Ex(i_ndx-1,j_ndx,k_ndx)+mxll%Ex(i_ndx,j_ndx,k_ndx))
                end do

            case (Ey_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    y    = y*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) y , z, 0.5*(mxll%Ey(i_ndx,j_ndx-1,k_ndx)+mxll%Ey(i_ndx,j_ndx,k_ndx))
                end do

            case (Ez_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    y    = y*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) y , z, 0.5*(mxll%Ez(i_ndx,j_ndx,k_ndx-1)+mxll%Ez(i_ndx,j_ndx,k_ndx))
                end do

            case (Hx_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    y    = y*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) y , z, 0.25*(mxll%Hx(i_ndx,j_ndx-1,k_ndx) + mxll%Hx(i_ndx,j_ndx,k_ndx) + &
                                                mxll%Hx(i_ndx,j_ndx,k_ndx-1) + mxll%Hx(i_ndx,j_ndx-1,k_ndx-1))
                end do

            case (Hy_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    y    = y*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) y , z, 0.25*(mxll%Hy(i_ndx-1,j_ndx,k_ndx) + mxll%Hy(i_ndx,j_ndx,k_ndx) +&
                                                mxll%Hy(i_ndx,j_ndx,k_ndx-1) + mxll%Hy(i_ndx-1,j_ndx,k_ndx-1))
                end do

            case (Hz_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    y    = y*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) y , z, 0.25*(mxll%Hz(i_ndx,j_ndx,k_ndx) + mxll%Hz(i_ndx-1,j_ndx,k_ndx) + &
                                                mxll%Hz(i_ndx,j_ndx-1,k_ndx) + mxll%Hz(i_ndx-1,j_ndx-1,k_ndx))
                end do

            end select

            close(funit)

        case (PLANE_ZX_DETECTOR)

            field_name = detectors(i)%f_ch
            full_dirname = trim(directory) // "_" // trim(number) // &
                               "/" // trim(field_name) // "_plane_" // trim(print_number) // ".dat"

            open(newunit=funit, file=trim(full_dirname), status='old', &
                action='write', position='append', iostat=ierr)

            select case(detectors(i)%field)

            case (Ex_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    x    = x*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) z , x, 0.5*(mxll%Ex(i_ndx-1,j_ndx,k_ndx)+mxll%Ex(i_ndx,j_ndx,k_ndx))
                end do

            case (Ey_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    x    = x*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) z , x, 0.5*(mxll%Ey(i_ndx,j_ndx-1,k_ndx)+mxll%Ey(i_ndx,j_ndx,k_ndx))
                end do

            case (Ez_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    x    = x*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) z , x, 0.5*(mxll%Ez(i_ndx,j_ndx,k_ndx-1)+mxll%Ez(i_ndx,j_ndx,k_ndx))
                end do

            case (Hx_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    x    = x*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) z , x, 0.25*(mxll%Hx(i_ndx,j_ndx-1,k_ndx) + mxll%Hx(i_ndx,j_ndx,k_ndx) + &
                                                mxll%Hx(i_ndx,j_ndx,k_ndx-1) + mxll%Hx(i_ndx,j_ndx-1,k_ndx-1))
                end do

            case (Hy_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    x    = x*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) z , x, 0.25*(mxll%Hy(i_ndx-1,j_ndx,k_ndx) + mxll%Hy(i_ndx,j_ndx,k_ndx) +&
                                                mxll%Hy(i_ndx,j_ndx,k_ndx-1) + mxll%Hy(i_ndx-1,j_ndx,k_ndx-1))
                end do

            case (Hz_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    x    = x*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) z , x, 0.25*(mxll%Hz(i_ndx-1,j_ndx-1,k_ndx) + mxll%Hz(i_ndx,j_ndx-1,k_ndx) +&
                                                mxll%Hz(i_ndx-1,j_ndx,k_ndx) + mxll%Hz(i_ndx,j_ndx,k_ndx))
                end do

            end select

            close(funit)

        case (VOLUME_DETECTOR)

            field_name = detectors(i)%f_ch
            full_dirname = trim(directory) // "_" // trim(number) // &
                            "/" // trim(field_name) // "_volume_" // trim(print_number) // ".dat"

            open(newunit=funit, file=trim(full_dirname), status='old', &
                action='write', position='append', iostat=ierr)

            select case(detectors(i)%field)

            case (Ex_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) x , y, z, 0.5*(mxll%Ex(i_ndx-1,j_ndx,k_ndx)+mxll%Ex(i_ndx,j_ndx,k_ndx))
                end do

            case (Ey_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) x , y, z, 0.5*(mxll%Ey(i_ndx,j_ndx-1,k_ndx)+mxll%Ey(i_ndx,j_ndx,k_ndx))
                end do

            case (Ez_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) x , y, z, 0.5*(mxll%Ez(i_ndx,j_ndx,k_ndx-1)+mxll%Ez(i_ndx,j_ndx,k_ndx))
                end do

            case (Hx_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) x , y, z, 0.25*(mxll%Hx(i_ndx,j_ndx-1,k_ndx) + mxll%Hx(i_ndx,j_ndx,k_ndx) + &
                                                mxll%Hx(i_ndx,j_ndx,k_ndx-1) + mxll%Hx(i_ndx,j_ndx-1,k_ndx-1))
                end do

            case (Hy_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) x , y, z, 0.25*(mxll%Hy(i_ndx-1,j_ndx,k_ndx) + mxll%Hy(i_ndx,j_ndx,k_ndx) +&
                                                mxll%Hy(i_ndx,j_ndx,k_ndx-1) + mxll%Hy(i_ndx-1,j_ndx,k_ndx-1))
                end do

            case (Hz_FIELD)

                do n = 1, detectors(i)%nd
                    i_ndx = detectors(i)%indx_list(n,1)
                    j_ndx = detectors(i)%indx_list(n,2)
                    k_ndx = detectors(i)%indx_list(n,3)
                    x    = (i_ndx+nx*mpi_coords(1))*dr - int(nx_tot/2)*dr
                    y    = (j_ndx+ny*mpi_coords(2))*dr - int(ny_tot/2)*dr
                    z    = (k_ndx+nz*mpi_coords(3))*dr - int(nz_tot/2)*dr
                    x    = x*au_to_nm
                    y    = y*au_to_nm
                    z    = z*au_to_nm
                    write(funit, *) x , y, z, 0.25*(mxll%Hz(i_ndx,j_ndx,k_ndx) + mxll%Hz(i_ndx-1,j_ndx,k_ndx) + &
                                                mxll%Hz(i_ndx,j_ndx-1,k_ndx) + mxll%Hz(i_ndx-1,j_ndx-1,k_ndx))
                end do

            end select

            close(funit)
        
        end select
    
        end if !----Here finishes the condition for detect_rank.

    end do !----Here finishes the loop over detectors.
#ifdef USE_MPI
        end if
        call mpi_barrier(MPI_COMM_WORLD, ierr)
    end do
#endif
end subroutine write_3D_field
!###################################################################################################

end module outputs_mod 