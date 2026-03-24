module input_mod

    use constants_mod

    implicit none

    ! Mxll boundary conditions: "close", "periodic", "cpml", "none"
    character(len = 10):: mxll_boundaries(3)    = "none"
    ! Mxll 2D mode: "TMz", "TEz", "full", "none"
    character(len = 10):: mxll_2D_mode          = "none"
    ! Mxll dimensions: 1, 2, 3
    integer            :: mxll_dimensions       = 1
    ! Mxll PML thickness in number of grid points.
    integer            :: mxll_npml             = 19
    ! Number of sources to check in the "sources.in" file.
    integer            :: mxll_n_src            = 0
    ! Number of media to check in the "medium_xxxx.in" file.
    integer            :: mxll_n_media          = 0
    ! Number of quantum groups to check in the "mol_group_xxxxxxx.in" file.
    integer            :: mxll_n_q_groups       = 0
    ! Number of MPI processes per axis, when running in parallel.
    integer            :: mpi_procs_per_axis(3) = 1
    ! Size of each dimension in nm for the Mxll simulation box. The origin is always moved at the center of the box.
    real(dp)           :: mxll_box_size(3)      = 0.0d0
    ! Total simulation time in fs.
    real(dp)           :: mxll_total_time       = 1.0
    ! Spatial step in nm.
    real(dp)           :: mxll_dr               = 1.0
    ! Time step in fs. It will be automatically reduced if it is larger than the stability limit for the Maxwell solver (dr/(2*c)).
    real(dp)           :: mxll_dt               = 10.0E10
    ! Quantum time step in fs.
    real(dp)           :: mxll_dt_q             = 0.01
    ! Density factor in units of 1/nm^3.
    real(dp)           :: mxll_density_factor   = 1.0e-3
    ! Relative permittivity of the medium (for now, it is not used).
    real(dp)           :: mxll_eps_r            = 1.0d0

contains

subroutine read_input_file(boundaries, mode_2D, dimensions, npml, grid_Ndims, &
                           Nt, dr, dt, dt_q, density_factor, mpi_dims, eps_r, n_src, n_media, &
                           n_q_groups)
    
    integer   , intent(out) :: boundaries(3)
    integer   , intent(out) :: mode_2D
    integer   , intent(out) :: dimensions
    integer   , intent(out) :: npml
    integer   , intent(out) :: grid_Ndims(3)
    integer   , intent(out) :: Nt
    integer   , intent(out) :: n_src
    integer   , intent(out) :: n_media
    integer   , intent(out) :: n_q_groups
    real(dp)  , intent(out) :: dr
    real(dp)  , intent(out) :: dt
    real(dp)  , intent(out) :: dt_q
    real(dp)  , intent(out) :: density_factor
    real(dp)  , intent(out) :: eps_r
    integer   , intent(out) :: mpi_dims(3)
    integer :: ierr, funit
    integer :: i 

    namelist /OMxRTA/ mxll_boundaries, mxll_2D_mode, mxll_dimensions, mxll_npml, &
                         mxll_n_src, mxll_box_size, mxll_total_time, mxll_dr, mxll_dt, &
                         mxll_density_factor, mxll_eps_r, mxll_n_media, mxll_n_q_groups, &
                         mxll_dt_q, mpi_procs_per_axis

    
    ! Check whether file exists.
    inquire (file="inp", iostat=ierr)
    
    if (ierr /= 0) then
        write (*, '("Error: input file inp does not exist")')
        error stop
    end if
    
    ! Open and read Namelist file.
    open (action='read', file="inp", iostat=ierr, newunit=funit)
    read (nml=OMxRTA, iostat=ierr, unit=funit)
    if (ierr /= 0) write (*, '("Error: invalid Namelist format")')
    
    close (funit)

#ifdef USE_MPI

    mpi_dims(1) = mpi_procs_per_axis(1)
    mpi_dims(2) = mpi_procs_per_axis(2)
    mpi_dims(3) = mpi_procs_per_axis(3)

    if (mxll_dimensions > 1) then
        grid_Ndims(1)  = int(mxll_box_size(1)/mxll_dr/mpi_dims(1))
        grid_Ndims(2)  = int(mxll_box_size(2)/mxll_dr/mpi_dims(2))
        grid_Ndims(3)  = int(mxll_box_size(3)/mxll_dr/mpi_dims(3))
    else if (mxll_dimensions == 1) then
        grid_Ndims(1)  = int(mxll_box_size(1)/mxll_dr)
        grid_Ndims(2)  = int(mxll_box_size(2)/mxll_dr)
        grid_Ndims(3)  = int(mxll_box_size(3)/mxll_dr)
    end if

#else
    grid_Ndims(1)  = int(mxll_box_size(1)/mxll_dr)
    grid_Ndims(2)  = int(mxll_box_size(2)/mxll_dr)
    grid_Ndims(3)  = int(mxll_box_size(3)/mxll_dr)
#endif

    dimensions     = mxll_dimensions
    npml           = mxll_npml
    n_src          = mxll_n_src
    n_media        = mxll_n_media
    n_q_groups     = mxll_n_q_groups
    dr             = mxll_dr * nm_to_au
    dt             = MIN(dr/(2.0*c0), mxll_dt)
    dt_q           = mxll_dt_q
    Nt             = INT((DBLE(mxll_total_time)*fs_to_au)/(dt))
    density_factor = mxll_density_factor/(nm_to_au)**3
    eps_r          = mxll_eps_r

    select case (mxll_2D_mode)
    case ("TMz")
        mode_2D = TMZ_2D_MODE
    case ("TEz")
        mode_2D = TEZ_2D_MODE
    case ("full")
        mode_2D = FULL_2D_MODE
    case ("none")
        write (*, '("none option selected for mxll_2D_mode")')  
    case default
        write (*, '("mxll_2D_mode: ", A)') mxll_2D_mode
        write (*, '("Error: invalid 2D mode")')
        error stop
    end select

    do i=1,3
        select case (mxll_boundaries(i))
        case ("close")
            boundaries(i) = CLOSE_BOUNDARIES
        case ("periodic")
            boundaries(i) = PERIODIC_BOUNDARIES
        case ("cpml")
            boundaries(i) = CPML_BOUNDARIES
        case ("none")
            write (*, '("none option selected for boundary condition in axis ", I0)') i
        case default
            write (*, '("Error: invalid boundary condition in axis ", I0)') i
            error stop
        end select
    end do

end subroutine read_input_file

end module input_mod