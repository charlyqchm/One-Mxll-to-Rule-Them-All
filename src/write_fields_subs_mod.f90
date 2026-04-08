module write_fields_subs_mod

    use constants_mod
    use mxll_1D_mod
    use mxll_2D_mod
    use mxll_3D_mod
    use detector_mod
#ifdef USE_MPI
    use mpi
#endif

    implicit none

contains

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

    character(len=20) :: directory = "./output_detector"
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

    character(len=20) :: directory = "./output_detector"
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
                    z    = z*au_to_nm
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
                    z    = z*au_to_nm
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

    character(len=20) :: directory = "./output_detector"
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

    character(len=20) :: directory = "./output_detector"
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
        
        if (detectors(i)%detect_rank) then  !---- Here it checks if the rank has a detector

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

        end if !-------- Here it finishes the check if the rank has a detector

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

    character(len=20) :: directory = "./output_detector"
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

    character(len=20) :: directory = "./output_detector"
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

end module write_fields_subs_mod