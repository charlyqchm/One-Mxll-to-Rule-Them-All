module parallel_subs_mod

#ifdef USE_MPI
    use mpi
#endif
    use constants_mod
    use rs_vec_base_mod
    use rs_vec_dimensions_mod
    use design_base_mod
    use design_1D_mod
    use design_2D_mod
    use design_3D_mod

    implicit none

    integer :: nprocs
    integer :: Xprev
    integer :: Xnext
    integer :: Ynext
    integer :: Yprev
    integer :: Znext
    integer :: Zprev
    integer :: old_comm
    integer :: cartesian_comm
    integer :: ndims
    integer :: n_procs
    
    logical, allocatable :: periods(:)
    integer, allocatable :: proc_coords(:,:)
    integer, allocatable :: dims(:)
    integer, allocatable :: coords(:)

    interface extend_array_to_x_ranks_2D
        module procedure extend_array_to_x_ranks_2D_R
        module procedure extend_array_to_x_ranks_2D_C
        module procedure extend_array_to_x_ranks_2D_L
    end interface extend_array_to_x_ranks_2D

    interface extend_array_to_y_ranks_2D
        module procedure extend_array_to_y_ranks_2D_R
        module procedure extend_array_to_y_ranks_2D_C
        module procedure extend_array_to_y_ranks_2D_L
    end interface extend_array_to_y_ranks_2D

    interface extend_array_to_x_ranks_3D
        module procedure extend_array_to_x_ranks_3D_R
        module procedure extend_array_to_x_ranks_3D_C
        module procedure extend_array_to_x_ranks_3D_L
    end interface extend_array_to_x_ranks_3D

    interface extend_array_to_y_ranks_3D
        module procedure extend_array_to_y_ranks_3D_R
        module procedure extend_array_to_y_ranks_3D_C
        module procedure extend_array_to_y_ranks_3D_L
    end interface extend_array_to_y_ranks_3D

    interface extend_array_to_z_ranks_3D
        module procedure extend_array_to_z_ranks_3D_R
        module procedure extend_array_to_z_ranks_3D_C
        module procedure extend_array_to_z_ranks_3D_L
    end interface extend_array_to_z_ranks_3D

contains
!###################################################################################################
subroutine init_parallelization(dimensions, mpi_coords, mpi_dims, mpi_nprocs, mpi_cart_comm, &
                                boundaries, myrank)
        integer, intent(in)    :: dimensions
        integer, intent(in)    :: mpi_nprocs
        integer, intent(in)    :: boundaries(3)   
        integer, intent(inout) :: myrank
        integer, intent(inout) :: mpi_dims(3)
        integer, intent(inout) :: mpi_coords(3)
        integer, intent(inout) :: mpi_cart_comm

        integer :: i
        integer :: coords_aux(dimensions)
        integer :: ierr
        integer :: nprocs_world
        integer :: requested_nprocs

#ifdef USE_MPI

        requested_nprocs = mpi_nprocs
        ndims = dimensions

        !~~~~=== setting and initializing MPI ===~~~~!
        !MPI is initialized before defining grid sizes. This should be separated in a different module
        !and also consdier non-parallelized cases.
        call mpi_init(ierr)
        call mpi_comm_size(MPI_COMM_WORLD, nprocs_world, ierr)      !<--- nprocs=total number of processors

        if (ierr /= 0) then
            write (*, '("Error: MPI_Comm_size failed")')
            error stop
        end if

        if (nprocs_world /= requested_nprocs) then
            write (*, '("Error: mpirun -np ", I0, " does not match mpi_dims product (", I0, ")")') &
                nprocs_world, requested_nprocs
            error stop
        end if

        nprocs = nprocs_world

        if (.not. allocated(periods)) allocate(periods(ndims))
        if (.not. allocated(dims)) allocate(dims(ndims))
        if (.not. allocated(coords)) allocate(coords(ndims))
        if (.not. allocated(proc_coords)) allocate(proc_coords(nprocs, ndims))

        dims   = 1
        coords = 0

        periods = .false.
        do i=1, ndims
            if (boundaries(i) == PERIODIC_BOUNDARIES) periods(i) = .true.
        end do
        
        call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)      !<--- myrank now is set to unique integer for each processor
       
        if (ndims == 1) then

            mpi_coords(1) = myrank
            cartesian_comm = MPI_COMM_WORLD

        else if (ndims > 1) then

            do i=1, ndims
                dims(i) = mpi_dims(i)
            end do

            old_comm=MPI_COMM_WORLD

            !~~~~=== creating Cartesian topology ===~~~~!
            call mpi_cart_create( &
            old_comm,         &        !<--- original communicator
            ndims,            &        !<--- ndims = integer, number of dimensions
            dims,             &        !<--- dims(ndims) = number of processors in each dimension
            periods,          &        !<--- periods(ndims) = logical array defining boundary conditions
            reorder,          &        !<--- reorder = logical, reorder or not processors [set to .false.]
            cartesian_comm,   &        !<--- new topology defined
            ierr              &        !<--- ierror = integer, error
                        )
            !<---- now I should have a Cartesian grid, first index is for X and the second one is for Y

            !~~~~=== get Cartesian coordinates ===~~~~!
            call mpi_cart_coords( &
            cartesian_comm, &
            myrank,         &
            ndims,          &
            coords,         &     !<--- coordinates(ndims) = gives integers identifying local coordinates of a given block for myrank
            ierr            &     !<--- coords(1) = x, coords(2) = y, coords(3) = z
                        )

            !write(*,*) 'processors and coordinates (X,Y,Z)'
            !write(*,*) myrank,coords(1),coords(2),coords(3)

            !~~~~=== map of neighbors for X ===~~~~!
            call mpi_cart_shift(&
            cartesian_comm, &
            0             , &     !<~~~ direction is X
            1             , &     !<~~~ displacement along X
            Xprev         , &     !<~~~ neighbor to the Xprev
            Xnext         , &     !<~~~ neighbor to the Xnext
            ierr            &
            )

            if (Xprev<0) then
                Xprev = MPI_PROC_NULL    !<~~~ making sure we don't send/receive outside of the grid
            endif
            if (Xnext<0) then
                Xnext = MPI_PROC_NULL    !<~~~ making sure we don't send/receive outside of the grid
            endif



            !write(*,*) 'processors and neighbors along X'
            !write(*,*) myrank,Xprev,Xnext

            !~~~~=== map of neighbors for Y ===~~~~!
            call mpi_cart_shift(&
            cartesian_comm, &
            1             , &     !<~~~ direction is Y
            1             , &     !<~~~ displacement along Y
            Yprev         , &     !<~~~ neighbor Yprev
            Ynext         , &     !<~~~ neighbor Ynext
            ierr            &
            )

            if (Yprev<0) then
                Yprev = MPI_PROC_NULL    !<~~~ making sure we don't send/receive outside of the grid
            endif
            if (Ynext<0) then
                Ynext = MPI_PROC_NULL    !<~~~ making sure we don't send/receive outside of the grid
            endif

            !write(*,*) 'processors and neighbors along Y'
            !write(*,*) myrank,Yprev,Ynext

            if (ndims == 3) then

                !~~~~=== map of neighbors for Z ===~~~~!
                call mpi_cart_shift(&
                cartesian_comm, &
                2             , &     !<~~~ direction is Z
                1             , &     !<~~~ displacement along Z
                Zprev         , &     !<~~~ neighbor Zprev
                Znext         , &     !<~~~ neighbor Znext
                ierr            &
                )

                if (Zprev<0) then
                    Zprev = MPI_PROC_NULL    !<~~~ making sure we don't send/receive outside of the grid
                endif
                if (Znext<0) then
                    Znext = MPI_PROC_NULL    !<~~~ making sure we don't send/receive outside of the grid
                endif

                !write(*,*) 'processors and neighbors along Z'
                !write(*,*) myrank,Zprev,Znext

            end if
            
            do i= 0, nprocs-1
                call mpi_cart_coords(cartesian_comm, i, ndims, coords_aux, ierr)
                proc_coords(i+1, :) = coords_aux(:)
            end do
        
            if (ndims == 2) then
                mpi_coords(1) = coords(1)
                mpi_coords(2) = coords(2)

            else if (ndims == 3) then
                mpi_coords(1) = coords(1)
                mpi_coords(2) = coords(2)
                mpi_coords(3) = coords(3)
            end if

        end if

        mpi_cart_comm = cartesian_comm

#else
        return
#endif

end subroutine init_parallelization
!###################################################################################################

subroutine finalize_parallelization()

    integer :: ierr

#ifdef USE_MPI
    if (allocated(periods)) deallocate(periods)
    if (allocated(dims)) deallocate(dims)
    if (allocated(coords)) deallocate(coords)
    if (allocated(proc_coords)) deallocate(proc_coords)

    call MPI_FINALIZE(ierr)
#else
    return
#endif

end subroutine finalize_parallelization

!###################################################################################################

subroutine extend_fvec_to_ranks(f_vec)

    class(TRSvec), intent(inout) :: f_vec

    select type (f_vec)
    type is (TRSvec_1D)
        return
    type is (TRSvec_2D)
        call extend_array_to_x_ranks_2D(f_vec%pl_y,f_vec%nx, f_vec%ny, f_vec%n_der)
        call extend_array_to_y_ranks_2D(f_vec%pl_x,f_vec%nx, f_vec%ny, f_vec%n_der)
        call extend_array_to_x_ranks_2D(f_vec%mi_y,f_vec%nx, f_vec%ny, f_vec%n_der)
        call extend_array_to_y_ranks_2D(f_vec%mi_x,f_vec%nx, f_vec%ny, f_vec%n_der)
        call extend_array_to_x_ranks_2D(f_vec%pl_z,f_vec%nx, f_vec%ny, f_vec%n_der)
        call extend_array_to_y_ranks_2D(f_vec%pl_z,f_vec%nx, f_vec%ny, f_vec%n_der)
        call extend_array_to_x_ranks_2D(f_vec%mi_z,f_vec%nx, f_vec%ny, f_vec%n_der)
        call extend_array_to_y_ranks_2D(f_vec%mi_z,f_vec%nx, f_vec%ny, f_vec%n_der)
    type is (TRSvec_3D)
        call extend_array_to_x_ranks_3D(f_vec%pl_y, f_vec%nx, f_vec%ny, f_vec%nz, f_vec%n_der)
        call extend_array_to_y_ranks_3D(f_vec%pl_x, f_vec%nx, f_vec%ny, f_vec%nz, f_vec%n_der)
        call extend_array_to_z_ranks_3D(f_vec%pl_x, f_vec%nx, f_vec%ny, f_vec%nz, f_vec%n_der)
        call extend_array_to_x_ranks_3D(f_vec%pl_z, f_vec%nx, f_vec%ny, f_vec%nz, f_vec%n_der)
        call extend_array_to_y_ranks_3D(f_vec%pl_z, f_vec%nx, f_vec%ny, f_vec%nz, f_vec%n_der)
        call extend_array_to_z_ranks_3D(f_vec%pl_y, f_vec%nx, f_vec%ny, f_vec%nz, f_vec%n_der)    
        call extend_array_to_x_ranks_3D(f_vec%mi_y, f_vec%nx, f_vec%ny, f_vec%nz, f_vec%n_der)
        call extend_array_to_y_ranks_3D(f_vec%mi_x, f_vec%nx, f_vec%ny, f_vec%nz, f_vec%n_der)
        call extend_array_to_z_ranks_3D(f_vec%mi_x, f_vec%nx, f_vec%ny, f_vec%nz, f_vec%n_der)
        call extend_array_to_x_ranks_3D(f_vec%mi_z, f_vec%nx, f_vec%ny, f_vec%nz, f_vec%n_der)
        call extend_array_to_y_ranks_3D(f_vec%mi_z, f_vec%nx, f_vec%ny, f_vec%nz, f_vec%n_der)
        call extend_array_to_z_ranks_3D(f_vec%mi_y, f_vec%nx, f_vec%ny, f_vec%nz, f_vec%n_der)    
    end select

end subroutine extend_fvec_to_ranks
!###################################################################################################

subroutine extend_rho_to_ranks(design)

    class(TDesign), intent(inout) :: design

    select type (design)
    type is (TDesign_1D)
        return
    type is (TDesign_2D)
        call extend_array_to_x_ranks_2D(design%rho, design%nx, design%ny, design%n_ker)
        call extend_array_to_y_ranks_2D(design%rho, design%nx, design%ny, design%n_ker)
    type is (TDesign_3D)
        call extend_array_to_x_ranks_3D(design%rho, design%nx, design%ny, design%nz, design%n_ker)
        call extend_array_to_y_ranks_3D(design%rho, design%nx, design%ny, design%nz, design%n_ker)
        call extend_array_to_z_ranks_3D(design%rho, design%nx, design%ny, design%nz, design%n_ker)
    end select

end subroutine extend_rho_to_ranks

!###################################################################################################

subroutine extend_grad_to_ranks(design)

    class(TDesign), intent(inout) :: design

    select type (design)
    type is (TDesign_1D)
        return
    type is (TDesign_2D)
        call extend_array_to_x_ranks_2D(design%grad, design%nx, design%ny, design%n_ker)
        call extend_array_to_y_ranks_2D(design%grad, design%nx, design%ny, design%n_ker)
    type is (TDesign_3D)
        call extend_array_to_x_ranks_3D(design%grad, design%nx, design%ny, design%nz, design%n_ker)
        call extend_array_to_y_ranks_3D(design%grad, design%nx, design%ny, design%nz, design%n_ker)
        call extend_array_to_z_ranks_3D(design%grad, design%nx, design%ny, design%nz, design%n_ker)
    end select

end subroutine extend_grad_to_ranks
!###################################################################################################

subroutine extend_opt_region_to_ranks(design)

    class(TDesign), intent(inout) :: design

    select type (design)
    type is (TDesign_1D)
        return
    type is (TDesign_2D)
        call extend_array_to_x_ranks_2D(design%opt_region, design%nx, design%ny, design%n_ker)
        call extend_array_to_y_ranks_2D(design%opt_region, design%nx, design%ny, design%n_ker)
    type is (TDesign_3D)
        call extend_array_to_x_ranks_3D(design%opt_region, design%nx, design%ny, &
                                        design%nz, design%n_ker)
        call extend_array_to_y_ranks_3D(design%opt_region, design%nx, design%ny, &
                                        design%nz, design%n_ker)
        call extend_array_to_z_ranks_3D(design%opt_region, design%nx, design%ny, &
                                        design%nz, design%n_ker)
    end select

end subroutine extend_opt_region_to_ranks

!###################################################################################################

subroutine extend_array_to_x_ranks_2D_C(array, nx, ny, n_ghost)

    complex(dp), intent(inout) :: array(-n_ghost+1:nx+n_ghost,-n_ghost+1:ny+n_ghost)
    integer    , intent(in)    :: nx
    integer    , intent(in)    :: ny
    integer    , intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
    integer     :: i, j, n
    complex(dp) :: send_buffer(ny*n_ghost)
    complex(dp) :: recv_buffer(ny*n_ghost)

    recv_buffer = Z_0

    n = 1
    do j = 1, ny
    do i = 1, n_ghost
        send_buffer(n) = array(i, j)
        n = n + 1
    end do
    end do

    call mpi_sendrecv(send_buffer, &  !<=== sending
    ny*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Xprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    recv_buffer, &               !<=== receiving
    ny*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Xnext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    n = 1
    do j = 1, ny
    do i = nx+1, nx+n_ghost
        array(i, j) = recv_buffer(n)
        n = n + 1
    end do
    end do

    n = 1
    do j = 1, ny
    do i = nx-n_ghost+1, nx
        send_buffer(n) = array(i, j)
        n = n + 1
    end do
    end do

    recv_buffer = Z_0

    !i+1,j
    call mpi_sendrecv(send_buffer, &  !<=== sending 
    ny*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Xnext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    recv_buffer, &               !<=== receiving
    ny*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Xprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    n = 1
    do j = 1, ny
    do i = -n_ghost+1, 0
        array(i, j) = recv_buffer(n)
        n = n + 1
    end do
    end do

#else

    return

#endif

end subroutine extend_array_to_x_ranks_2D_C

!###################################################################################################

subroutine extend_array_to_x_ranks_2D_R(array, nx, ny, n_ghost)

    real(dp), intent(inout) :: array(-n_ghost+1:nx+n_ghost,-n_ghost+1:ny+n_ghost)
    integer , intent(in)    :: nx
    integer , intent(in)    :: ny
    integer , intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
    integer     :: i, j, n
    real(dp)    :: send_buffer(ny*n_ghost)
    real(dp)    :: recv_buffer(ny*n_ghost)

    recv_buffer = 0.0_dp

    n = 1
    do j = 1, ny
    do i = 1, n_ghost
        send_buffer(n) = array(i, j)
        n = n + 1
    end do
    end do

    call mpi_sendrecv(send_buffer, &  !<=== sending
    ny*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Xprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    recv_buffer, &               !<=== receiving
    ny*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Xnext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    n = 1
    do j = 1, ny
    do i = nx+1, nx+n_ghost
        array(i, j) = recv_buffer(n)
        n = n + 1
    end do
    end do

    n = 1
    do j = 1, ny
    do i = nx-n_ghost+1, nx
        send_buffer(n) = array(i, j)
        n = n + 1
    end do
    end do

    recv_buffer = 0.0_dp


    !i+1,j
    call mpi_sendrecv(send_buffer, &  !<=== sending
    ny*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Xnext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    recv_buffer, &               !<=== receiving
    ny*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Xprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    n = 1
    do j = 1, ny
    do i = -n_ghost+1, 0
        array(i, j) = recv_buffer(n)
        n = n + 1
    end do
    end do

#else

    return

#endif

end subroutine extend_array_to_x_ranks_2D_R

!###################################################################################################

subroutine extend_array_to_x_ranks_2D_L(array, nx, ny, n_ghost)

    logical, intent(inout) :: array(-n_ghost+1:nx+n_ghost,-n_ghost+1:ny+n_ghost)
    integer, intent(in)    :: nx
    integer, intent(in)    :: ny
    integer, intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
    integer     :: i, j, n
    logical     :: send_buffer(ny*n_ghost)
    logical     :: recv_buffer(ny*n_ghost)


    recv_buffer = .false.

    n = 1
    do j = 1, ny
    do i = 1, n_ghost
        send_buffer(n) = array(i, j)
        n = n + 1
    end do
    end do

    call mpi_sendrecv(send_buffer, &  !<=== sending
    ny*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Xprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    recv_buffer, &               !<=== receiving
    ny*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Xnext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    n = 1
    do j = 1, ny
    do i = nx+1, nx+n_ghost
        array(i, j) = recv_buffer(n)
        n = n + 1
    end do
    end do

    n = 1
    do j = 1, ny
    do i = nx-n_ghost+1, nx
        send_buffer(n) = array(i, j)
        n = n + 1
    end do
    end do

    recv_buffer = .false.

    call mpi_sendrecv(send_buffer, &  !<=== sending
    ny*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Xnext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    recv_buffer, &                                   !<=== receiving
    ny*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Xprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    n = 1
    do j = 1, ny
    do i = -n_ghost+1, 0
        array(i, j) = recv_buffer(n)
        n = n + 1
    end do
    end do

#else
    return
#endif

end subroutine extend_array_to_x_ranks_2D_L
!###################################################################################################

subroutine extend_array_to_y_ranks_2D_C(array, nx, ny, n_ghost)

    complex(dp), intent(inout) :: array(-n_ghost+1:nx+n_ghost,-n_ghost+1:ny+n_ghost)
    integer    , intent(in)    :: nx
    integer    , intent(in)    :: ny
    integer    , intent(in)    :: n_ghost

#ifdef USE_MPI
    integer     :: ierr
    integer     :: istatus(MPI_STATUS_SIZE)
    integer     :: i, j, n
    complex(dp) :: send_buffer(nx*n_ghost)
    complex(dp) :: recv_buffer(nx*n_ghost)

    recv_buffer = Z_0

    n = 1
    do j = 1, n_ghost
    do i = 1, nx
        send_buffer(n) = array(i, j)
        n = n + 1
    end do
    end do

    call mpi_sendrecv(send_buffer, &  !<=== sending
    nx*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Yprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    recv_buffer, &                                   !<=== receiving
    nx*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Ynext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    n = 1
    do j = ny+1, ny+n_ghost
    do i = 1, nx
        array(i, j) = recv_buffer(n)
        n = n + 1
    end do
    end do

    n = 1
    do j = ny-n_ghost+1, ny
    do i = 1, nx
        send_buffer(n) = array(i, j)
        n = n + 1
    end do
    end do

    recv_buffer = Z_0

    call mpi_sendrecv(send_buffer, &  !<=== sending
    nx*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Ynext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    recv_buffer, &               !<=== receiving
    nx*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Yprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    n = 1
    do j = -n_ghost+1, 0
    do i = 1, nx
        array(i, j) = recv_buffer(n)
        n = n + 1
    end do
    end do

#else
    return
#endif
end subroutine extend_array_to_y_ranks_2D_C

!###################################################################################################

subroutine extend_array_to_y_ranks_2D_R(array, nx, ny, n_ghost)

    real(dp), intent(inout) :: array(-n_ghost+1:nx+n_ghost,-n_ghost+1:ny+n_ghost)
    integer , intent(in)    :: nx
    integer , intent(in)    :: ny
    integer , intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
    integer     :: i, j, n
    real(dp)    :: send_buffer(nx*n_ghost)
    real(dp)    :: recv_buffer(nx*n_ghost)

    recv_buffer = 0.0_dp

    n = 1
    do j = 1, n_ghost
    do i = 1, nx
        send_buffer(n) = array(i, j)
        n = n + 1
    end do
    end do

    call mpi_sendrecv(send_buffer, &  !<=== sending
    nx*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Yprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    recv_buffer, &                                   !<=== receiving
    nx*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Ynext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    n = 1
    do j = ny+1, ny+n_ghost
    do i = 1, nx
        array(i, j) = recv_buffer(n)
        n = n + 1
    end do
    end do

    n = 1
    do j = ny-n_ghost+1, ny
    do i = 1, nx
        send_buffer(n) = array(i, j)
        n = n + 1
    end do
    end do

    recv_buffer = 0.0_dp

    call mpi_sendrecv(send_buffer, &  !<=== sending
    nx*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Ynext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    recv_buffer, &                                   !<=== receiving
    nx*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Yprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    n = 1
    do j = -n_ghost+1, 0
    do i = 1, nx
        array(i, j) = recv_buffer(n)
        n = n + 1
    end do
    end do

#else
    return
#endif
end subroutine extend_array_to_y_ranks_2D_R

!###################################################################################################

subroutine extend_array_to_y_ranks_2D_L(array, nx, ny, n_ghost)

    logical, intent(inout) :: array(-n_ghost+1:nx+n_ghost,-n_ghost+1:ny+n_ghost)
    integer, intent(in)    :: nx
    integer, intent(in)    :: ny
    integer, intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
    integer     :: i, j, n
    logical     :: send_buffer(nx*n_ghost)
    logical     :: recv_buffer(nx*n_ghost)

    recv_buffer = .false.

    n = 1
    do j = 1, n_ghost
    do i = 1, nx
        send_buffer(n) = array(i, j)
        n = n + 1
    end do
    end do

    call mpi_sendrecv(send_buffer, &  !<=== sending
    nx*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Yprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    recv_buffer, &                                   !<=== receiving
    nx*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Ynext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    n = 1
    do j = ny+1, ny+n_ghost
    do i = 1, nx
        array(i, j) = recv_buffer(n)
        n = n + 1
    end do
    end do

    n = 1
    do j = ny-n_ghost+1, ny
    do i = 1, nx
        send_buffer(n) = array(i, j)
        n = n + 1
    end do
    end do

    recv_buffer = .false.

    call mpi_sendrecv(send_buffer, &  !<=== sending
    nx*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Ynext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    recv_buffer, &                                   !<=== receiving
    nx*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Yprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    n = 1
    do j = -n_ghost+1, 0
    do i = 1, nx
        array(i, j) = recv_buffer(n)
        n = n + 1
    end do
    end do

#else
    return
#endif

end subroutine extend_array_to_y_ranks_2D_L

!###################################################################################################

subroutine extend_array_to_x_ranks_3D_C(array, nx, ny, nz, n_ghost)

    complex(dp), intent(inout) :: array(:,:,:)
    integer    , intent(in)    :: nx
    integer    , intent(in)    :: ny
    integer    , intent(in)    :: nz
    integer    , intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    array(-n_ghost+1:0, 1:ny, 1:nz)     = Z_0
    array(nx+1:nx+n_ghost, 1:ny, 1:nz)  = Z_0

    call mpi_sendrecv(array(1:n_ghost, 1:ny, 1:nz), &  !<=== sending
    ny*nz*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Xprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    array(nx+1:nx+n_ghost, 1:ny, 1:nz), &               !<=== receiving
    ny*nz*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Xnext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    call mpi_sendrecv(array(nx-n_ghost+1:nx, 1:ny, 1:nz), &  !<=== sending 
    ny*nz*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Xnext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    array(-n_ghost+1:0, 1:ny, 1:nz), &               !<=== receiving
    ny*nz*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Xprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

#else
    return
#endif

end subroutine extend_array_to_x_ranks_3D_C

!###################################################################################################

subroutine extend_array_to_x_ranks_3D_R(array, nx, ny, nz, n_ghost)

    real(dp), intent(inout) :: array(:,:,:)
    integer , intent(in)    :: nx
    integer , intent(in)    :: ny
    integer , intent(in)    :: nz
    integer , intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    array(-n_ghost+1:0, 1:ny, 1:nz)     = 0.0_dp
    array(nx+1:nx+n_ghost, 1:ny, 1:nz)  = 0.0_dp

    call mpi_sendrecv(array(1:n_ghost, 1:ny, 1:nz), &  !<=== sending
    ny*nz*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Xprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    array(nx+1:nx+n_ghost, 1:ny, 1:nz), &               !<=== receiving
    ny*nz*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Xnext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    call mpi_sendrecv(array(nx-n_ghost+1:nx, 1:ny, 1:nz), &  !<=== sending 
    ny*nz*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Xnext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    array(-n_ghost+1:0, 1:ny, 1:nz), &               !<=== receiving
    ny*nz*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Xprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

#else
    return
#endif

end subroutine extend_array_to_x_ranks_3D_R

!###################################################################################################

subroutine extend_array_to_x_ranks_3D_L(array, nx, ny, nz, n_ghost)

    logical, intent(inout) :: array(:,:,:)
    integer, intent(in)    :: nx
    integer, intent(in)    :: ny
    integer, intent(in)    :: nz
    integer, intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    array(-n_ghost+1:0, 1:ny, 1:nz)     = .false.
    array(nx+1:nx+n_ghost, 1:ny, 1:nz)  = .false.

    call mpi_sendrecv(array(1:n_ghost, 1:ny, 1:nz), &  !<=== sending
    ny*nz*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Xprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    array(nx+1:nx+n_ghost, 1:ny, 1:nz), &               !<=== receiving
    ny*nz*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Xnext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    call mpi_sendrecv(array(nx-n_ghost+1:nx, 1:ny, 1:nz), &  !<=== sending 
    ny*nz*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Xnext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    array(-n_ghost+1:0, 1:ny, 1:nz), &               !<=== receiving
    ny*nz*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Xprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

#else
    return
#endif

end subroutine extend_array_to_x_ranks_3D_L

!###################################################################################################

subroutine extend_array_to_y_ranks_3D_C(array, nx, ny, nz, n_ghost)

    complex(dp), intent(inout) :: array(:,:,:)
    integer    , intent(in)    :: nx
    integer    , intent(in)    :: ny
    integer    , intent(in)    :: nz
    integer    , intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    array(1:nx, -n_ghost+1:0, 1:nz)     = Z_0
    array(1:nx, ny+1:ny+n_ghost, 1:nz)  = Z_0

    call mpi_sendrecv(array(1:nx, 1:n_ghost, 1:nz), &  !<=== sending
    nx*nz*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Yprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    array(1:nx, ny+1:ny+n_ghost, 1:nz), &               !<=== receiving
    nx*nz*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Ynext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    call mpi_sendrecv(array(1:nx, ny-n_ghost+1:ny, 1:nz), &  !<=== sending
    nx*nz*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Ynext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    array(1:nx, -n_ghost+1:0, 1:nz), &               !<=== receiving
    nx*nz*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Yprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

#else
    return
#endif

end subroutine extend_array_to_y_ranks_3D_C

!###################################################################################################

subroutine extend_array_to_y_ranks_3D_R(array, nx, ny, nz, n_ghost)

    real(dp), intent(inout) :: array(:,:,:)
    integer , intent(in)    :: nx
    integer , intent(in)    :: ny
    integer , intent(in)    :: nz
    integer , intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    array(1:nx, -n_ghost+1:0, 1:nz)     = 0.0_dp
    array(1:nx, ny+1:ny+n_ghost, 1:nz)  = 0.0_dp

    call mpi_sendrecv(array(1:nx, 1:n_ghost, 1:nz), &  !<=== sending
    nx*nz*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Yprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    array(1:nx, ny+1:ny+n_ghost, 1:nz), &               !<=== receiving
    nx*nz*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Ynext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    call mpi_sendrecv(array(1:nx, ny-n_ghost+1:ny, 1:nz), &  !<=== sending
    nx*nz*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Ynext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    array(1:nx, -n_ghost+1:0, 1:nz), &               !<=== receiving
    nx*nz*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Yprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

#else
    return
#endif

end subroutine extend_array_to_y_ranks_3D_R

!###################################################################################################

subroutine extend_array_to_y_ranks_3D_L(array, nx, ny, nz, n_ghost)

    logical, intent(inout) :: array(:,:,:)
    integer, intent(in)    :: nx
    integer, intent(in)    :: ny
    integer, intent(in)    :: nz
    integer, intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    array(1:nx, -n_ghost+1:0, 1:nz)     = .false.
    array(1:nx, ny+1:ny+n_ghost, 1:nz)  = .false.

    call mpi_sendrecv(array(1:nx, 1:n_ghost, 1:nz), &  !<=== sending
    nx*nz*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Yprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    array(1:nx, ny+1:ny+n_ghost, 1:nz), &               !<=== receiving
    nx*nz*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Ynext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    call mpi_sendrecv(array(1:nx, ny-n_ghost+1:ny, 1:nz), &  !<=== sending
    nx*nz*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Ynext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    array(1:nx, -n_ghost+1:0, 1:nz), &               !<=== receiving
    nx*nz*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Yprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

#else
    return
#endif

end subroutine extend_array_to_y_ranks_3D_L

!###################################################################################################

subroutine extend_array_to_z_ranks_3D_C(array, nx, ny, nz, n_ghost)

    complex(dp), intent(inout) :: array(:,:,:)
    integer    , intent(in)    :: nx
    integer    , intent(in)    :: ny
    integer    , intent(in)    :: nz
    integer    , intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    array(1:nx, 1:ny, -n_ghost+1:0)     = Z_0
    array(1:nx, 1:ny, nz+1:nz+n_ghost)  = Z_0

    call mpi_sendrecv(array(1:nx, 1:ny, 1:n_ghost), &  !<=== sending
    nx*ny*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Zprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    array(1:nx, 1:ny, nz+1:nz+n_ghost), &               !<=== receiving
    nx*ny*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Znext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    call mpi_sendrecv(array(1:nx, 1:ny, nz-n_ghost+1:nz), &  !<=== sending
    nx*ny*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Znext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    array(1:nx, 1:ny, -n_ghost+1:0), &               !<=== receiving
    nx*ny*n_ghost, &                                       !<=== the size
    mpi_double_complex, &                          !<=== type
    Zprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

#else
    return
#endif

end subroutine extend_array_to_z_ranks_3D_C
!###################################################################################################

subroutine extend_array_to_z_ranks_3D_R(array, nx, ny, nz, n_ghost)

    real(dp), intent(inout) :: array(:,:,:)
    integer , intent(in)    :: nx
    integer , intent(in)    :: ny
    integer , intent(in)    :: nz
    integer , intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    array(1:nx, 1:ny, -n_ghost+1:0)     = 0.0_dp
    array(1:nx, 1:ny, nz+1:nz+n_ghost)  = 0.0_dp

    call mpi_sendrecv(array(1:nx, 1:ny, 1:n_ghost), &  !<=== sending
    nx*ny*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Zprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    array(1:nx, 1:ny, nz+1:nz+n_ghost), &               !<=== receiving
    nx*ny*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Znext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    call mpi_sendrecv(array(1:nx, 1:ny, nz-n_ghost+1:nz), &  !<=== sending
    nx*ny*n_ghost, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Znext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    array(1:nx, 1:ny, -n_ghost+1:0), &               !<=== receiving
    nx*ny*n_ghost, &                                      !<=== the size
    mpi_double_precision, &                          !<=== type
    Zprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

#else
    return
#endif

end subroutine extend_array_to_z_ranks_3D_R
!###################################################################################################

subroutine extend_array_to_z_ranks_3D_L(array, nx, ny, nz, n_ghost)

    logical, intent(inout) :: array(:,:,:)
    integer , intent(in)    :: nx
    integer , intent(in)    :: ny
    integer , intent(in)    :: nz
    integer , intent(in)    :: n_ghost

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    array(1:nx, 1:ny, -n_ghost+1:0)     = .false.
    array(1:nx, 1:ny, nz+1:nz+n_ghost)  = .false.

    call mpi_sendrecv(array(1:nx, 1:ny, 1:n_ghost), &  !<=== sending
    nx*ny*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Zprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    array(1:nx, 1:ny, nz+1:nz+n_ghost), &               !<=== receiving
    nx*ny*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Znext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

    call mpi_sendrecv(array(1:nx, 1:ny, nz-n_ghost+1:nz), &  !<=== sending
    nx*ny*n_ghost, &                                       !<=== the size
    mpi_logical, &                          !<=== type
    Znext, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    array(1:nx, 1:ny, -n_ghost+1:0), &               !<=== receiving
    nx*ny*n_ghost, &                                      !<=== the size
    mpi_logical, &                          !<=== type
    Zprev, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code

#else
    return
#endif

end subroutine extend_array_to_z_ranks_3D_L

!###################################################################################################
end module parallel_subs_mod