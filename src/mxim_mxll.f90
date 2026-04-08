program mxim_mxll

    use constants_mod
    use mxll_base_mod
    use factory_mod
    use input_mod
    use source_mod
    use interactions_mod
    use parallel_subs_mod
    use q_group_mod
    use detector_mod
    use outputs_mod

    implicit none

    class(TMxll),    allocatable :: mxll
    type(TSource),   allocatable :: source_list(:)
    type(TQ_group),  allocatable :: q_groups(:)
    type(TDetector), allocatable :: detectors(:)
    integer                     :: boundaries(3)
    integer                     :: mode_2D
    integer                     :: mpi_coords(3) = 0
    integer                     :: mpi_dims(3)   = 1
    integer                     :: myrank = 0
    integer                     :: dimensions
    integer                     :: npml
    integer                     :: n_src
    integer                     :: n_media
    integer                     :: n_q_groups
    integer                     :: grid_Ndims(3)
    integer                     :: Nt
    integer                     :: Nt_q
    integer                     :: tt_0
    integer                     :: n_detectors
    integer                     :: t_det_print
    integer                     :: t_q_print
    integer                     :: print_det_step = 0 ! Counter of the printing seps of the detectors.
    integer                     :: print_q_step = 0 ! Counter of the printing seps of the quantum system.
    logical                     :: move_q_system
    real(dp)                    :: dr 
    real(dp)                    :: dt
    real(dp)                    :: dt_q
    real(dp)                    :: density_factor
    real(dp)                    :: eps_r
    real(dp)                    :: dt_det_print
    real(dp)                    :: dt_q_print

    integer  :: tt
    integer  :: i, j
    real(dp) :: time

    !TO-DO: With a huge number of MPI_ranks, this subroutine may be a bottleneck.
    call read_input_file(boundaries, mode_2D, dimensions, npml, grid_Ndims, &
                         Nt, dr, dt, dt_q, density_factor, mpi_dims, eps_r, n_src, n_media, &
                         n_q_groups, n_detectors, dt_det_print, dt_q_print)

    n_procs = mpi_dims(1)*mpi_dims(2)*mpi_dims(3)

    call init_parallelization(dimensions, mpi_coords, mpi_dims, n_procs, boundaries, myrank)

    mxll = maxwell_factory(dimensions)

    if (.not. allocated(q_groups)) allocate(q_groups(n_q_groups))

    call read_init_sources(source_list, n_src, dimensions, dr, grid_Ndims, mpi_coords, mpi_dims)

    call mxll%init(grid_Ndims, npml, boundaries, dt, dr, mode_2D, n_media, mpi_coords, mpi_dims) 
    
    call init_detectors_outputs(n_detectors, t_det_print, dt_det_print, detectors, dimensions, &
                                grid_Ndims, mxll%mode, dr, dt, mpi_dims, mpi_coords, myrank)

    !Define the real time step for the quantum system (smaller than the mxll dt)
    mxll%n_skip_steps = 1
    if (dt_q > dt) mxll%n_skip_steps = NINT(dt_q/dt)

    tt_0 = 0
    if (mxll%n_skip_steps == 1) tt_0 = 1
    
    mxll%tq_step_old = 0

    dt_q = dble(mxll%n_skip_steps)*dt
    Nt_q = 1 + int(Nt/mxll%n_skip_steps)

    do i = 1, n_q_groups
        call q_groups(i)%init_q_group(i, dimensions, mpi_dims, mpi_coords, & 
                                      cartesian_comm, myrank, dr, dt_q, Nt_q, &
                                      density_factor, grid_Ndims)
    end do

    call init_q_groups_outputs(q_groups, n_q_groups, dt_q_print, t_q_print, dt, myrank)
    
    do tt = 1, Nt
        time = dt*(DBLE(tt)-0.5d0)
        
        mxll%time = time
        
        mxll%tq_step = 1 + int((tt-tt_0)/mxll%n_skip_steps)
        move_q_system = .false.
        if (mxll%tq_step > mxll%tq_step_old) then
            mxll%t_skip       = time
            mxll%tq_step_old  = mxll%tq_step
            move_q_system     = .true.
            if (n_q_groups == 0) move_q_system = .false.
        end if
       
        call exchange_E_field_between_ranks(mxll)
        call mxll%td_propagate_H_field()   
       
        call update_sources(source_list, n_src, time)
        call source_interactions(mxll, source_list, n_src)
        
        call exchange_H_field_between_ranks(mxll)
        call mxll%td_propagate_E_field(tt)
        
        call expand_E_field_between_ranks(mxll, move_q_system)
        
        do i = 1, n_q_groups
            call send_E_to_J_ranks(mxll, q_groups(i), move_q_system, myrank)
            call q_groups(i)%td_propagate_q_group(mxll%tq_step, move_q_system)
            call send_J_to_E_ranks(mxll, q_groups(i), move_q_system, myrank)
        end do 
        call expand_J_field_between_ranks(mxll, move_q_system)

        call extend_fields_to_detectors(mxll, detectors, n_detectors, tt, t_det_print)

        call write_detectors_outputs(detectors, mxll, n_detectors, tt, print_det_step, &
                                   t_det_print, time, dr, grid_Ndims, mpi_dims, mpi_coords, myrank)

        call write_q_groups_outputs(q_groups, n_q_groups, t_q_print, tt, time, print_q_step,&
                                    myrank, mpi_dims)

    end do

    call kill_sources(source_list, n_src)
    call mxll%kill()
    do i = 1, n_q_groups
        call q_groups(i)%kill_q_group()
    end do

    call kill_detectors_outputs(detectors, n_detectors)

    if (allocated(mxll))     deallocate(mxll)
    if (allocated(q_groups)) deallocate(q_groups)
    call finalize_parallelization()
end program mxim_mxll