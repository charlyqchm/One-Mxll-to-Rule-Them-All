module input_mod

    use constants_mod
    implicit none

    !Boundaries of the simulation box.
    !This can be: "closed", "periodic", "pml"
    character(len=10) :: mxll_boundaries(3)
    
    !Number of dimensions of the simulation box. This can be 1, 2 or 3.
    integer           :: mxll_dimensions = 1
    
    !Number of PML layers in every direction of the simulation box.
    integer           :: mxll_n_pml=20

    !Grid steps in every direction.     
    real(dp)          :: mxll_dr
    
    !Dimensions of the simulation box in units of dr.
    real(dp)          :: mxll_box_size(3)

    !List of frequencies of each optimization problem.
    real(dp)          :: mxll_freq_list(100)

    !List of permittivity values of each optimization problem.
    real(dp)          :: mxll_eps_Re(100) = 1.0_dp
    real(dp)          :: mxll_eps_Im(100) = 0.0_dp

    !Number of optimization problems to solve within the same
    !function of merit.
    integer           :: design_n_opt_problems
    
    !Maximum number of iteration steps for the optimization.
    integer           :: design_max_iter_steps = 100

    !If true, the optimization will scan every term from
    !delta_rho. Once this does not improve the objective function,
    !the optimization will continue with the next beta value or
    !stop if there are no more beta values to scan.
    logical           :: design_converge_optimization = .true.

    !If true, the optimization will restart from a restart file.
    logical           :: design_restart = .false.

    !List of delta_rho values to scan in the optimization. Being rho
    !the design variable and delta_rho the step size of the optimization.
    real(dp)          :: design_delta_rho(50)

    !List of beta values to scan in the optimization. 
    real(dp)          :: design_beta(50)
    
    !Initial value of the design variable rho for the optimization.
    real(dp)          :: design_rho_init = 0.0_dp

    !This variable control the position of the binary function in the optimization.
    !It is not recomended to change its value. 
    real(dp)          :: design_eta = 0.5_dp

    !This variable control the gaussian function used over rho, to give
    !more dimension to the optimization variable.
    real(dp)          :: design_sigma_rho = 0.0_dp

    !Maximum number of iteration steps for the BiCGSTAB_L solver 
    !used in the optimization.
    integer           :: design_bicgstab_max_iter = 1000

    !Number of L parameter of the BiCGSTAB_L solver.
    !L=1 corresponds to the standard BiCGSTAB solver.
    integer           :: design_bicgstab_l_term = 5

    !Tolerance for the BiCGSTAB_L solver used in the optimization.
    real(dp)          :: design_bicgstab_tol = 1.0d-6

    !Number of MPI processes per axis. This is only used 
    !if the code is compiled with MPI support.
    integer           :: mpi_procs_per_axis(3) = 1

contains

!###################################################################################################

subroutine read_input_file(boundaries, restart, converge_optimization, n_opt_problems,        &
                           max_iter_steps, dimensions, n_pml, grid_Ndims, dr, freq_list,      &
                           eps_Re, eps_Im, delta_rho, beta, eta, rho_init, sigma_rho, mpi_dims, &
                           n_delta_rho_steps, n_beta_steps, bicgstab_max_iter, bicgstab_l_term, &
                           bicgstab_tol)
    
    integer   , intent(out) :: boundaries(3)
    logical   , intent(out) :: restart
    logical   , intent(out) :: converge_optimization
    integer   , intent(out) :: n_opt_problems
    integer   , intent(out) :: max_iter_steps
    integer   , intent(out) :: dimensions
    integer   , intent(out) :: n_pml
    integer   , intent(out) :: grid_Ndims(3)
    real(dp)  , intent(out) :: dr
    real(dp)  , intent(out) :: freq_list(100)
    real(dp)  , intent(out) :: eps_Re(100)
    real(dp)  , intent(out) :: eps_Im(100)
    real(dp)  , intent(out) :: delta_rho(50)
    real(dp)  , intent(out) :: beta(50)
    real(dp)  , intent(out) :: eta
    real(dp)  , intent(out) :: rho_init
    real(dp)  , intent(out) :: sigma_rho
    real(dp)  , intent(out) :: bicgstab_tol
    integer   , intent(out) :: mpi_dims(3)
    integer   , intent(out) :: n_delta_rho_steps
    integer   , intent(out) :: n_beta_steps
    integer   , intent(out) :: bicgstab_max_iter
    integer   , intent(out) :: bicgstab_l_term


    integer :: ierr, funit
    integer :: i

    namelist /OMxRTA/ mxll_boundaries, mxll_dimensions, mxll_n_pml, mxll_box_size, mxll_dr, mxll_freq_list, &
                         mxll_eps_Re, mxll_eps_Im, design_n_opt_problems, design_max_iter_steps, &
                         design_converge_optimization, design_restart, design_delta_rho, design_beta, &
                         design_eta, design_rho_init, design_sigma_rho, design_bicgstab_max_iter,        &
                         design_bicgstab_l_term, design_bicgstab_tol, mpi_procs_per_axis

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
    end if else if (mxll_dimensions == 1) then
        grid_Ndims(1)  = int(mxll_box_size(1)/mxll_dr)
        grid_Ndims(2)  = int(mxll_box_size(2)/mxll_dr)
        grid_Ndims(3)  = int(mxll_box_size(3)/mxll_dr)
    end if
#else
    
    grid_Ndims(1)  = int(mxll_box_size(1)/mxll_dr)
    grid_Ndims(2)  = int(mxll_box_size(2)/mxll_dr)
    grid_Ndims(3)  = int(mxll_box_size(3)/mxll_dr)

#endif

    dimensions            = mxll_dimensions
    n_pml                 = mxll_n_pml
    dr                    = mxll_dr
    freq_list             = mxll_freq_list
    eps_Re                = mxll_eps_Re
    eps_Im                = mxll_eps_Im
    n_opt_problems        = design_n_opt_problems
    max_iter_steps        = design_max_iter_steps
    converge_optimization = design_converge_optimization
    restart               = design_restart
    delta_rho             = design_delta_rho
    beta                  = design_beta
    eta                   = design_eta
    rho_init              = design_rho_init
    sigma_rho             = design_sigma_rho
    bicgstab_max_iter     = design_bicgstab_max_iter
    bicgstab_L_term       = design_bicgstab_l_term
    bicgstab_tol          = design_bicgstab_tol

    do i = 1, 3
        select case (mxll_boundaries(i))
        case ("closed")
            boundaries(i) = CLOSE_BOUNDARIES
        case ("periodic")
            boundaries(i) = PERIODIC_BOUNDARIES
        case ("pml")
            boundaries(i) = PML_BOUNDARIES
        case default
            write (*, '("Error: invalid boundary condition in axis ", I0)') i
            error stop
        end select
    end do

    i=1
    n_delta_rho_steps = 0
    do while (delta_rho(i) /= 0.0_dp .and. i < size(delta_rho))
        n_delta_rho_steps = n_delta_rho_steps + 1
        i = i + 1
    end do    

    i=1
    n_beta_steps = 0
    do while (beta(i) /= 0.0_dp .and. i < size(beta))
        n_beta_steps = n_beta_steps + 1
        i = i + 1
    end do

end subroutine read_input_file

!###################################################################################################

end module input_mod