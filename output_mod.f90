module outputs_mod

    use constants_mod
    use mxll_base_mod
    use mxll_1D_mod
    use mxll_2D_mod
    use mxll_3D_mod
    use write_fields_subs_mod
    use detector_mod
    use q_group_mod
    
    implicit none

contains

!###################################################################################################

subroutine init_detectors_outputs(n_detectors, t_det_print, dt_det_print, detectors, dimensions, &
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
    character(len=20) :: dirname = "./output_detector"
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
                    &x = ", F10.5, " nm,   y = ", F10.5, " nm,   z = ", F10.5, " nm")') x, y, z
                case (Ey_FIELD)
                    write(funit, '("# Time (a.u.)              Ey (a.u.),      &
                    &x = ", F10.5, " nm,   y = ", F10.5, " nm,   z = ", F10.5, " nm")') x, y, z
                case (Ez_FIELD)
                    write(funit, '("# Time (a.u.)              Ez (a.u.),      &
                    &x = ", F10.5, " nm,   y = ", F10.5, " nm,   z = ", F10.5, " nm")') x, y, z
                case (Hx_FIELD)
                    write(funit, '("# Time (a.u.)              Hx (a.u.),      &
                    &x = ", F10.5, " nm,   y = ", F10.5, " nm,   z = ", F10.5, " nm")') x, y, z
                case (Hy_FIELD)
                    write(funit, '("# Time (a.u.)              Hy (a.u.),      &
                    &x = ", F10.5, " nm,   y = ", F10.5, " nm,   z = ", F10.5, " nm")') x, y, z
                case (Hz_FIELD)
                    write(funit, '("# Time (a.u.)              Hz (a.u.),      &
                    &x = ", F10.5, " nm,   y = ", F10.5, " nm,   z = ", F10.5, " nm")') x, y, z
                end select
                close(funit)
            end if
       end do 
    end if

end subroutine init_detectors_outputs

!###################################################################################################

subroutine kill_detectors_outputs(detectors, n_detectors)

    type(TDetector), allocatable, intent(inout) :: detectors(:)
    integer , intent(in)  :: n_detectors

    integer :: i

    if (.not. allocated(detectors)) return

    do i = 1, n_detectors
        call kill_detector(detectors(i))
    end do

    deallocate(detectors)

end subroutine kill_detectors_outputs

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

subroutine init_q_groups_outputs(q_groups, n_q_groups, dt_q_print, t_q_print, dt, myrank)

    type(TQ_group), intent(in)  :: q_groups(:)
    integer       , intent(in)  :: n_q_groups
    integer       , intent(in)  :: myrank
    integer       , intent(out) :: t_q_print
    real(dp)      , intent(in)  :: dt_q_print
    real(dp)      , intent(in)  :: dt

    character(len=20) :: dirname = "./output_q_group_"
    character(len=4)  :: number
    integer :: i
    integer :: ierr

    if (n_q_groups == 0) return

    t_q_print = int(dt_q_print/dt)

    if (myrank == 0) then
        
        do i = 1, n_q_groups
            if (q_groups(i)%print_option /= PRINT_NONE) then

                write(number, '(I4.4)') i

                call execute_command_line("mkdir -p " // trim(dirname) // trim(number), exitstat=ierr)
                
                if (ierr /= 0) then
                    write(*, *) "Error: Could not create directory ", trim(dirname) // trim(number)
                    stop
                end if
            else
                write(*, *) "Output for q_group ", i, " is turned off."
            end if
        end do

    end if

end subroutine init_q_groups_outputs

!###################################################################################################

subroutine write_q_groups_outputs(q_groups, n_q_groups, t_q_print, t_step, time, print_q_step,&
                                  myrank, mpi_dims)

    type(TQ_group), intent(in)    :: q_groups(:)
    integer       , intent(in)    :: n_q_groups
    integer       , intent(in)    :: myrank
    integer       , intent(in)    :: t_step
    integer       , intent(in)    :: t_q_print
    integer       , intent(in)    :: mpi_dims(3)
    real(dp)      , intent(in)    :: time
    integer       , intent(inout) :: print_q_step
    
    character(len=20) :: dirname = "./output_q_group_"
    character(len=4)  :: number
    character(len=30) :: file_name = "energy_and_dipole.dat"
    character(len=99) :: full_dirname
    integer :: i
    integer :: n
    integer :: ierr
    integer :: funit
    
#ifdef USE_MPI
    integer :: rank_iter, nprocs_mpi

    nprocs_mpi = mpi_dims(1)*mpi_dims(2)*mpi_dims(3)
#endif

    if (n_q_groups == 0) return
    
    if (MOD(t_step, t_q_print) /= 0) return
   
    if (myrank == 0) then
        do n = 1, n_q_groups
            if (q_groups(n)%print_option /= PRINT_NONE) then
                
                write(number, '(I4.4)') n
                
                full_dirname = trim(dirname) // trim(number) // "/" // trim(file_name)
                
                if (print_q_step == 0) then
                    open(newunit=funit, file=full_dirname, status='replace', action='write', iostat=ierr)
                else
                    open(newunit=funit, file=full_dirname, status='old', action='write', position='append', iostat=ierr)
                end if

                if (ierr /= 0) then
                    write(*, *) "Error: Could not open file ", trim(full_dirname)
                    stop
                end if

                write(funit, '("# Time = ", ES18.8, " (a.u.)")') time
                write(funit, '("# N_mol        Energy (eV)                   Dipole_x (a.u.)             &
                            &Dipole_y (a.u.)             Dipole_z (a.u.)")')

                close(funit)

            end if
        end do
    end if

#ifdef USE_MPI
    do rank_iter = 0, nprocs_mpi-1
        if (myrank == rank_iter) then
#endif

    do n = 1, n_q_groups
        if (q_groups(n)%print_option /= PRINT_NONE) then
            
            write(number, '(I4.4)') n
            
            full_dirname = trim(dirname) // trim(number) // "/" // trim(file_name)
            
            open(newunit=funit, file=trim(full_dirname), status='old', &
                action='write', position='append', iostat=ierr)

            if (ierr /= 0) then
                write(*, *) "Error: Could not open file ", trim(full_dirname)
                stop
            end if

            do i = 1, q_groups(n)%n_sys_loc
                if (q_groups(n)%q_sys(i)%print_on) then
                    write(funit, *) q_groups(n)%q_sys(i)%id, &
                                    q_groups(n)%q_sys(i)%Et*au_to_ev, &
                                    q_groups(n)%q_sys(i)%dipole(1), &
                                    q_groups(n)%q_sys(i)%dipole(2), &
                                    q_groups(n)%q_sys(i)%dipole(3) 
                end if
            end do

            close(funit)

        end if
    end do      

#ifdef USE_MPI
        end if
        call mpi_barrier(MPI_COMM_WORLD, ierr)
    end do
#endif

    do n = 1, n_q_groups
    do i = 1, q_groups(n)%n_sys_loc
        call q_groups(n)%q_sys(i)%write_output(time, n, print_q_step)
    end do
    end do

    print_q_step = print_q_step + 1

end subroutine write_q_groups_outputs

!###################################################################################################

end module outputs_mod 