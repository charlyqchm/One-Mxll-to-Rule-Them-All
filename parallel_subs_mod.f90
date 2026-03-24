module parallel_subs_mod

#ifdef USE_MPI
    use mpi
#endif
    use constants_mod
    use mxll_base_mod
    use mxll_1D_mod
    use mxll_2D_mod
    use mxll_3D_mod
    use q_group_mod

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

contains

!###################################################################################################

    subroutine init_parallelization(dimensions, mpi_coords, mpi_dims, mpi_nprocs, boundaries, myrank)
        integer, intent(in)    :: dimensions
        integer, intent(in)    :: mpi_nprocs
        integer, intent(in)    :: boundaries(3)   
        integer, intent(inout) :: myrank
        integer, intent(inout) :: mpi_dims(3)
        integer, intent(inout) :: mpi_coords(3)

        integer :: i
        integer :: coords_aux(dimensions)
        integer :: ierr

#ifdef USE_MPI

        nprocs = mpi_nprocs
        ndims = dimensions
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

        !~~~~=== setting and initializing MPI ===~~~~!
        !MPI is initialized before defining grid sizes. This should be separated in a different module
        !and also consdier non-parallelized cases.
        call mpi_init(ierr)
        call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr)      !<--- nprocs=total number of processors
        
        if (ierr /= 0) then
            write (*, '("Error: number of processors does not match the number specified in the input file")')
            error stop
        end if
        
        call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)      !<--- myrank now is set to unique integer for each processor
       
        if (ndims == 1) then

            mpi_coords(1) = myrank

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

subroutine exchange_E_field_between_ranks(mxll)

    class(TMxll), intent(inout) :: mxll

#ifdef USE_MPI

    select type(mxll)
    class is(TMxll_1D)

        return

    class is(TMxll_2D)

        call exchange_E_field_2D(mxll)

    class is(TMxll_3D)

        call exchange_E_field_3D(mxll)

    end select

#else
    return
#endif

end subroutine exchange_E_field_between_ranks

!###################################################################################################

subroutine exchange_H_field_between_ranks(mxll)

    class(TMxll), intent(inout) :: mxll

#ifdef USE_MPI

    select type(mxll)
    class is(TMxll_1D)

        return

    class is(TMxll_2D)

        call exchange_H_field_2D(mxll)

    class is(TMxll_3D)

        call exchange_H_field_3D(mxll)

    end select
#else
    return
#endif
end subroutine exchange_H_field_between_ranks

!###################################################################################################

subroutine expand_E_field_between_ranks(mxll, move_q_system)

    class(TMxll), intent(inout) :: mxll
    logical     , intent(in)    :: move_q_system

    if (.not. move_q_system) return

#ifdef USE_MPI

    select type(mxll)
    class is(TMxll_1D)

        return

    class is(TMxll_2D)

        call expand_E_field_2D(mxll)

    class is(TMxll_3D)

        call expand_E_field_3D(mxll)
    end select

#else
    return
#endif

end subroutine expand_E_field_between_ranks

!###################################################################################################

subroutine expand_J_field_between_ranks(mxll, move_q_system)

    class(TMxll), intent(inout) :: mxll
    logical     , intent(in)    :: move_q_system

    if (.not. move_q_system) return

#ifdef USE_MPI

    select type(mxll)
    class is(TMxll_1D)

        return

    class is(TMxll_2D)

        call expand_J_field_2D(mxll)

    class is(TMxll_3D)

        call expand_J_field_3D(mxll)
    end select

#else
    return
#endif

end subroutine expand_J_field_between_ranks

!###################################################################################################

subroutine exchange_E_field_3D(mxll)
    class(TMxll_3D), intent(inout) :: mxll

    integer :: nx
    integer :: ny
    integer :: nz
    
#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    nx = mxll%nx
    ny = mxll%ny
    nz = mxll%nz

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Expanding arrays of E  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
    !i,j+1,k
    call mpi_sendrecv(mxll%Ex(1:nx, 1, 1:nz), &  !<=== sending
    nx * nz, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Yprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    mxll%Ex(1:nx, ny + 1, 1:nz), &               !<=== receiving
    nx * nz, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Ynext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code
    !i,j,k+1
    call mpi_sendrecv(mxll%Ex(1:nx, 1:ny, 1), &  !<=== sending
    nx * ny, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Zprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag
    mxll%Ex(1:nx, 1:ny, nz + 1), &               !<=== receiving
    nx * ny, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Znext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code
    !i+1,j,k
    call mpi_sendrecv(mxll%Ey(1, 1:ny, 1:nz), &  !<=== sending 
    ny * nz, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Xprev, &                                         !<=== where sending
    MPI_GOOD_TAG, &                                  !<=== sending tag   
    mxll%Ey(nx + 1, 1:ny, 1:nz), &               !<=== receiving
    ny * nz, &                                       !<=== the size
    mpi_double_precision, &                          !<=== type
    Xnext, &                                         !<=== receiving from where
    MPI_GOOD_TAG, &                                  !<=== receiving tag
    cartesian_comm, &                                !<=== handle of Cartesian coordinates
    istatus, &                                       !<=== istatus
    ierr)                                            !<=== error code
    !i,j,k+1
    call mpi_sendrecv(mxll%Ey(1:nx, 1:ny, 1), &     !<=== sending
    nx * ny, &     !<=== the size
    mpi_double_precision, &     !<=== type
    Zprev, &     !<=== where sending
    MPI_GOOD_TAG, &     !<=== sending tag
    mxll%Ey(1:nx, 1:ny, nz + 1), &     !<=== receiving
    nx * ny, &     !<=== the size
    mpi_double_precision, &     !<=== type
    Znext, &     !<=== receiving from where
    MPI_GOOD_TAG, &     !<=== receiving tag
    cartesian_comm, &     !<=== handle of Cartesian coordinates
    istatus, &     !<=== istatus
    ierr)                       !<=== error code
    !i+1,j,k
    call mpi_sendrecv(mxll%Ez(1, 1:ny, 1:nz), &     !<=== sending
    ny * nz, &     !<=== the size
    mpi_double_precision, &     !<=== type
    Xprev, &     !<=== where sending
    MPI_GOOD_TAG, &     !<=== sending tag
    mxll%Ez(nx + 1, 1:ny, 1:nz), &     !<=== receiving
    ny * nz, &     !<=== the size
    mpi_double_precision, &     !<=== type
    Xnext, &     !<=== receiving from where
    MPI_GOOD_TAG, &     !<=== receiving tag
    cartesian_comm, &     !<=== handle of Cartesian coordinates
    istatus, &     !<=== istatus
    ierr)                       !<=== error code
    !i,j+1,k
    call mpi_sendrecv(mxll%Ez(1:nx, 1, 1:nz), &     !<=== sending
    nx * nz, &     !<=== the size
    mpi_double_precision, &     !<=== type
    Yprev, &     !<=== where sending
    MPI_GOOD_TAG, &     !<=== sending tag
    mxll%Ez(1:nx, ny + 1, 1:nz), &     !<=== receiving
    nx * nz, &     !<=== the size
    mpi_double_precision, &     !<=== type
    Ynext, &     !<=== receiving from where
    MPI_GOOD_TAG, &     !<=== receiving tag
    cartesian_comm, &     !<=== handle of Cartesian coordinates
    istatus, &     !<=== istatus
    ierr)                       !<=== error code

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

#else

    return

#endif

end subroutine exchange_E_field_3D

!###################################################################################################

subroutine exchange_E_field_2D(mxll)
    class(TMxll_2D), intent(inout) :: mxll

    integer :: nx
    integer :: ny
    
#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    nx = mxll%nx
    ny = mxll%ny

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Expanding arrays of E  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then

        !i,j+1
        call mpi_sendrecv(mxll%Ex(1:nx, 1), &  !<=== sending
        nx, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Yprev, &                                         !<=== where sending
        MPI_GOOD_TAG, &                                  !<=== sending tag
        mxll%Ex(1:nx, ny + 1), &               !<=== receiving
        nx, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Ynext, &                                         !<=== receiving from where
        MPI_GOOD_TAG, &                                  !<=== receiving tag
        cartesian_comm, &                                !<=== handle of Cartesian coordinates
        istatus, &                                       !<=== istatus
        ierr)                                            !<=== error code

        !i+1,j
        call mpi_sendrecv(mxll%Ey(1, 1:ny), &  !<=== sending 
        ny, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Xprev, &                                         !<=== where sending
        MPI_GOOD_TAG, &                                  !<=== sending tag   
        mxll%Ey(nx + 1, 1:ny), &               !<=== receiving
        ny, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Xnext, &                                         !<=== receiving from where
        MPI_GOOD_TAG, &                                  !<=== receiving tag
        cartesian_comm, &                                !<=== handle of Cartesian coordinates
        istatus, &                                       !<=== istatus
        ierr)                                            !<=== error code

    end if

    if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
        !i+1,j,k
        call mpi_sendrecv(mxll%Ez(1, 1:ny), &     !<=== sending
        ny, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Xprev, &     !<=== where sending
        MPI_GOOD_TAG, &     !<=== sending tag
        mxll%Ez(nx + 1, 1:ny), &     !<=== receiving
        ny, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Xnext, &     !<=== receiving from where
        MPI_GOOD_TAG, &     !<=== receiving tag
        cartesian_comm, &     !<=== handle of Cartesian coordinates
        istatus, &     !<=== istatus
        ierr)                       !<=== error code
        !i,j+1
        call mpi_sendrecv(mxll%Ez(1:nx, 1), &     !<=== sending
        nx, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Yprev, &     !<=== where sending
        MPI_GOOD_TAG, &     !<=== sending tag
        mxll%Ez(1:nx, ny + 1), &     !<=== receiving
        nx, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Ynext, &     !<=== receiving from where
        MPI_GOOD_TAG, &     !<=== receiving tag
        cartesian_comm, &     !<=== handle of Cartesian coordinates
        istatus, &     !<=== istatus
        ierr)                       !<=== error code
    
    end if

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

#else

    return

#endif

end subroutine exchange_E_field_2D

!###################################################################################################

subroutine exchange_H_field_3D(mxll)
    class(TMxll_3D), intent(inout) :: mxll

    integer :: nx
    integer :: ny
    integer :: nz
    
#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
    
    nx = mxll%nx
    ny = mxll%ny
    nz = mxll%nz

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Expanding arrays of H  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
    !i,j-1,k
        call mpi_sendrecv(mxll%Hx(1:nx, ny, 1:nz), &     !<=== sending
        nx * nz, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Ynext, &     !<=== where sending
        MPI_GOOD_TAG, &     !<=== sending tag
        mxll%Hx(1:nx, 0, 1:nz), &     !<=== receiving
        nx * nz, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Yprev, &     !<=== receiving from where
        MPI_GOOD_TAG, &     !<=== receiving tag
        cartesian_comm, &     !<=== handle of Cartesian coordinates
        istatus, &     !<=== istatus
        ierr)                       !<=== error code
    !i,j,k-1
    call mpi_sendrecv(mxll%Hx(1:nx, 1:ny, nz), &     !<=== sending
            nx * ny, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Znext, &     !<=== where sending
            MPI_GOOD_TAG, &     !<=== sending tag
            mxll%Hx(1:nx, 1:ny, 0), &     !<=== receiving
            nx * ny, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Zprev, &     !<=== receiving from where
            MPI_GOOD_TAG, &     !<=== receiving tag
            cartesian_comm, &     !<=== handle of Cartesian coordinates
            istatus, &     !<=== istatus
            ierr)                       !<=== error code
    !i-1,j,k
    call mpi_sendrecv(mxll%Hy(nx, 1:ny, 1:nz), &     !<=== sending
            ny * nz, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Xnext, &     !<=== where sending
            MPI_GOOD_TAG, &     !<=== sending tag
            mxll%Hy(0, 1:ny, 1:nz), &     !<=== receiving
            ny * nz, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Xprev, &     !<=== receiving from where
            MPI_GOOD_TAG, &     !<=== receiving tag
            cartesian_comm, &     !<=== handle of Cartesian coordinates
            istatus, &     !<=== istatus
            ierr)                       !<=== error code
    !i,j,k-1
    call mpi_sendrecv(mxll%Hy(1:nx, 1:ny, nz), &     !<=== sending
            nx * ny, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Znext, &     !<=== where sending
            MPI_GOOD_TAG, &     !<=== sending tag
            mxll%Hy(1:nx, 1:ny, 0), &     !<=== receiving
            nx * ny, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Zprev, &     !<=== receiving from where
            MPI_GOOD_TAG, &     !<=== receiving tag
            cartesian_comm, &     !<=== handle of Cartesian coordinates
            istatus, &     !<=== istatus
            ierr)                       !<=== error code
    !i-1,j,k
    call mpi_sendrecv(mxll%Hz(nx, 1:ny, 1:nz), &     !<=== sending
            ny * nz, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Xnext, &     !<=== where sending
            MPI_GOOD_TAG, &     !<=== sending tag
            mxll%Hz(0, 1:ny, 1:nz), &     !<=== receiving
            ny * nz, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Xprev, &     !<=== receiving from where
            MPI_GOOD_TAG, &     !<=== receiving tag
            cartesian_comm, &     !<=== handle of Cartesian coordinates
            istatus, &     !<=== istatus
            ierr)                       !<=== error code
    !i,j-1,k
    call mpi_sendrecv(mxll%Hz(1:nx, ny, 1:nz), &     !<=== sending
            nx * nz, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Ynext, &     !<=== where sending
            MPI_GOOD_TAG, &     !<=== sending tag
            mxll%Hz(1:nx, 0, 1:nz), &     !<=== receiving
            nx * nz, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Yprev, &     !<=== receiving from where
            MPI_GOOD_TAG, &     !<=== receiving tag
            cartesian_comm, &     !<=== handle of Cartesian coordinates
            istatus, &     !<=== istatus
            ierr)                       !<=== error code       
    
    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

#else 
    return 
#endif    

end subroutine exchange_H_field_3D

!###################################################################################################

subroutine exchange_H_field_2D(mxll)
    class(TMxll_2D), intent(inout) :: mxll

    integer :: nx
    integer :: ny
    
#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    nx = mxll%nx
    ny = mxll%ny

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Expanding arrays of H  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
    

    if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
        !i,j-1
            call mpi_sendrecv(mxll%Hx(1:nx, ny), &     !<=== sending
            nx, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Ynext, &     !<=== where sending
            MPI_GOOD_TAG, &     !<=== sending tag
            mxll%Hx(1:nx, 0), &     !<=== receiving
            nx, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Yprev, &     !<=== receiving from where
            MPI_GOOD_TAG, &     !<=== receiving tag
            cartesian_comm, &     !<=== handle of Cartesian coordinates
            istatus, &     !<=== istatus
            ierr)                       !<=== error code
        !i-1,j
        call mpi_sendrecv(mxll%Hy(nx, 1:ny), &     !<=== sending
                ny, &     !<=== the size
                mpi_double_precision, &     !<=== type
                Xnext, &     !<=== where sending
                MPI_GOOD_TAG, &     !<=== sending tag
                mxll%Hy(0, 1:ny), &     !<=== receiving
                ny, &     !<=== the size
                mpi_double_precision, &     !<=== type
                Xprev, &     !<=== receiving from where
                MPI_GOOD_TAG, &     !<=== receiving tag
                cartesian_comm, &     !<=== handle of Cartesian coordinates
                istatus, &     !<=== istatus
                ierr)                       !<=== error code
    end if

    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
        !i-1,j
        call mpi_sendrecv(mxll%Hz(nx, 1:ny), &     !<=== sending
                ny, &     !<=== the size
                mpi_double_precision, &     !<=== type
                Xnext, &     !<=== where sending
                MPI_GOOD_TAG, &     !<=== sending tag
                mxll%Hz(0, 1:ny), &     !<=== receiving
                ny, &     !<=== the size
                mpi_double_precision, &     !<=== type
                Xprev, &     !<=== receiving from where
                MPI_GOOD_TAG, &     !<=== receiving tag
                cartesian_comm, &     !<=== handle of Cartesian coordinates
                istatus, &     !<=== istatus
                ierr)                       !<=== error code
        !i,j-1
        call mpi_sendrecv(mxll%Hz(1:nx, ny), &     !<=== sending
                nx, &     !<=== the size
                mpi_double_precision, &     !<=== type
                Ynext, &     !<=== where sending
                MPI_GOOD_TAG, &     !<=== sending tag
                mxll%Hz(1:nx, 0), &     !<=== receiving
                nx, &     !<=== the size
                mpi_double_precision, &     !<=== type
                Yprev, &     !<=== receiving from where
                MPI_GOOD_TAG, &     !<=== receiving tag
                cartesian_comm, &     !<=== handle of Cartesian coordinates
                istatus, &     !<=== istatus
                ierr)                       !<=== error code      
    end if
    
    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

#else 
    return 
#endif    

end subroutine exchange_H_field_2D

!###################################################################################################

subroutine expand_E_field_2D(mxll)
    class(TMxll_2D), intent(inout) :: mxll

    integer :: nx
    integer :: ny
    
#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    nx = mxll%nx
    ny = mxll%ny

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Expanding arrays of E fields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
        
        call mpi_sendrecv(mxll%Ex(nx, 1:ny), &  !<=== sending
        ny, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Xnext, &                                         !<=== where sending
        MPI_GOOD_TAG, &                                  !<=== sending tag
        mxll%Ex(0, 1:ny), &               !<=== receiving
        ny, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Xprev, &                                         !<=== receiving from where
        MPI_GOOD_TAG, &                                  !<=== receiving tag
        cartesian_comm, &                                !<=== handle of Cartesian coordinates
        istatus, &                                       !<=== istatus
        ierr)                                            !<=== error code

        call mpi_sendrecv(mxll%Ey(1:nx, ny), &  !<=== sending
        nx, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Ynext, &                                         !<=== where sending
        MPI_GOOD_TAG, &                                  !<=== sending tag
        mxll%Ex(1:nx, 0), &                       !<=== receiving
        nx, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Yprev, &                                         !<=== receiving from where
        MPI_GOOD_TAG, &                                  !<=== receiving tag
        cartesian_comm, &                                !<=== handle of Cartesian coordinates
        istatus, &                                       !<=== istatus
        ierr)                                            !<=== error code

    end if

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

#else
    return
#endif   

end subroutine expand_E_field_2D

!###################################################################################################

subroutine expand_E_field_3D(mxll)
    class(TMxll_3D), intent(inout) :: mxll

    integer :: nx
    integer :: ny
    integer :: nz
    
#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

    nx = mxll%nx
    ny = mxll%ny
    nz = mxll%nz

    call mpi_sendrecv(mxll%Ex(nx, 1:ny, 1:nz), &  !<=== sending
        ny*nz, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Xnext, &                                         !<=== where sending
        MPI_GOOD_TAG, &                                  !<=== sending tag
        mxll%Ex(0, 1:ny, 1:nz), &               !<=== receiving
        ny*nz, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Xprev, &                                         !<=== receiving from where
        MPI_GOOD_TAG, &                                  !<=== receiving tag
        cartesian_comm, &                                !<=== handle of Cartesian coordinates
        istatus, &                                       !<=== istatus
        ierr)                                            !<=== error code

    call mpi_sendrecv(mxll%Ey(1:nx, ny, 1:nz), &  !<=== sending
        nx*nz, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Ynext, &                                         !<=== where sending
        MPI_GOOD_TAG, &                                  !<=== sending tag
        mxll%Ey(1:nx, 0, 1:nz), &               !<=== receiving
        nx*nz, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Yprev, &                                         !<=== receiving from where
        MPI_GOOD_TAG, &                                  !<=== receiving tag
        cartesian_comm, &                                !<=== handle of Cartesian coordinates
        istatus, &                                       !<=== istatus
        ierr)                                            !<=== error code

    call mpi_sendrecv(mxll%Ez(1:nx, 1:ny, nz), &  !<=== sending
        nx*ny, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Znext, &                                         !<=== where sending
        MPI_GOOD_TAG, &                                  !<=== sending tag
        mxll%Ez(1:nx, 1:ny, 0), &               !<=== receiving
        nx*ny, &                                       !<=== the size
        mpi_double_precision, &                          !<=== type
        Zprev, &                                         !<=== receiving from where
        MPI_GOOD_TAG, &                                  !<=== receiving tag
        cartesian_comm, &                                !<=== handle of Cartesian coordinates
        istatus, &                                       !<=== istatus
        ierr)                                            !<=== error code

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

#else
    return
#endif

end subroutine expand_E_field_3D

!###################################################################################################

subroutine expand_J_field_2D(mxll)
    class(TMxll_2D), intent(inout) :: mxll

    integer :: nx
    integer :: ny

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
#endif

    nx = mxll%nx
    ny = mxll%ny

#ifdef USE_MPI

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
    
        call mpi_sendrecv(mxll%Jx_old(1,1:ny), &     !<=== sending
            ny, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Xprev, &     !<=== where sending
            MPI_GOOD_TAG, &     !<=== sending tag
            mxll%Jx_old(nx+1,1:ny), &     !<=== receiving
            ny, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Xnext, &     !<=== receiving from where
            MPI_GOOD_TAG, &     !<=== receiving tag
            cartesian_comm, &     !<=== handle of Cartesian coordinates
            istatus, &     !<=== istatus
            ierr)                       !<=== error code

        call mpi_sendrecv(mxll%Jy_old(1:nx,1), &     !<=== sending
            nx, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Yprev, &     !<=== where sending
            MPI_GOOD_TAG, &     !<=== sending tag
            mxll%Jy_old(1:nx,ny+1), &     !<=== receiving
            nx, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Ynext, &     !<=== receiving from where
            MPI_GOOD_TAG, &     !<=== receiving tag
            cartesian_comm, &     !<=== handle of Cartesian coordinates
            istatus, &     !<=== istatus
            ierr)
    
        call mpi_sendrecv(mxll%dJx(1,1:ny), &     !<=== sending
            ny, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Xprev, &     !<=== where sending
            MPI_GOOD_TAG, &     !<=== sending tag
            mxll%dJx(nx+1,1:ny), &     !<=== receiving
            ny, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Xnext, &     !<=== receiving from where
            MPI_GOOD_TAG, &     !<=== receiving tag
            cartesian_comm, &     !<=== handle of Cartesian coordinates
            istatus, &     !<=== istatus
            ierr)                       !<=== error code

        call mpi_sendrecv(mxll%dJy(1:nx,1), &     !<=== sending
            nx, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Yprev, &     !<=== where sending
            MPI_GOOD_TAG, &     !<=== sending tag
            mxll%dJy(1:nx,ny+1), &     !<=== receiving
            nx, &     !<=== the size
            mpi_double_precision, &     !<=== type
            Ynext, &     !<=== receiving from where
            MPI_GOOD_TAG, &     !<=== receiving tag
            cartesian_comm, &     !<=== handle of Cartesian coordinates
            istatus, &     !<=== istatus
            ierr)

    end if

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

#else
    return
#endif

end subroutine expand_J_field_2D

!###################################################################################################

subroutine expand_J_field_3D(mxll)
    class(TMxll_3D), intent(inout) :: mxll

    integer :: nx
    integer :: ny
    integer :: nz

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
#endif

    nx = mxll%nx
    ny = mxll%ny
    nz = mxll%nz

#ifdef USE_MPI

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

    call mpi_sendrecv(mxll%Jx_old(1,1:ny,1:nz), &     !<=== sending
        ny*nz, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Xprev, &     !<=== where sending
        MPI_GOOD_TAG, &     !<=== sending tag
        mxll%Jx_old(nx+1,1:ny,1:nz), &     !<=== receiving
        ny*nz, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Xnext, &     !<=== receiving from where
        MPI_GOOD_TAG, &     !<=== receiving tag
        cartesian_comm, &     !<=== handle of Cartesian coordinates
        istatus, &     !<=== istatus
        ierr)                       !<=== error code

    call mpi_sendrecv(mxll%Jy_old(1:nx,1,1:nz), &     !<=== sending
        nx*nz, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Yprev, &     !<=== where sending
        MPI_GOOD_TAG, &     !<=== sending tag
        mxll%Jy_old(1:nx,ny+1,1:nz), &     !<=== receiving
        nx*nz, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Ynext, &     !<=== receiving from where
        MPI_GOOD_TAG, &     !<=== receiving tag
        cartesian_comm, &     !<=== handle of Cartesian coordinates
        istatus, &     !<=== istatus
        ierr)

    call mpi_sendrecv(mxll%Jz_old(1:nx,1:ny,1), &     !<=== sending
        nx*ny, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Zprev, &     !<=== where sending
        MPI_GOOD_TAG, &    !<=== sending tag
        mxll%Jz_old(1:nx,1:ny,nz+1), &     !<=== receiving
        nx*ny, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Znext, &     !<=== receiving from where
        MPI_GOOD_TAG, &     !<=== receiving tag
        cartesian_comm, &     !<=== handle of Cartesian coordinates
        istatus, &     !<=== istatus
        ierr)                       !<=== error code

    call mpi_sendrecv(mxll%dJx(1,1:ny,1:nz), &     !<=== sending
        ny*nz, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Xprev, &     !<=== where sending
        MPI_GOOD_TAG, &     !<=== sending tag
        mxll%dJx(nx+1,1:ny,1:nz), &     !<=== receiving
        ny*nz, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Xnext, &     !<=== receiving from where
        MPI_GOOD_TAG, &     !<=== receiving tag
        cartesian_comm, &     !<=== handle of Cartesian coordinates
        istatus, &     !<=== istatus
        ierr)                       !<=== error code

    call mpi_sendrecv(mxll%dJy(1:nx,1,1:nz), &     !<=== sending
        nx*nz, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Yprev, &     !<=== where sending
        MPI_GOOD_TAG, &     !<=== sending tag
        mxll%dJy(1:nx,ny+1,1:nz), &     !<=== receiving
        nx*nz, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Ynext, &     !<=== receiving from where
        MPI_GOOD_TAG, &     !<=== receiving tag
        cartesian_comm, &     !<=== handle of Cartesian coordinates
        istatus, &     !<=== istatus
        ierr)

    call mpi_sendrecv(mxll%dJz(1:nx,1:ny,1), &     !<=== sending
        nx*ny, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Zprev, &     !<=== where sending
        MPI_GOOD_TAG, &    !<=== sending tag
        mxll%dJz(1:nx,1:ny,nz+1), &     !<=== receiving
        nx*ny, &     !<=== the size
        mpi_double_precision, &     !<=== type
        Znext, &     !<=== receiving from where
        MPI_GOOD_TAG, &     !<=== receiving tag
        cartesian_comm, &     !<=== handle of Cartesian coordinates
        istatus, &     !<=== istatus
        ierr)                       !<=== error code

    call MPI_BARRIER( MPI_COMM_WORLD, ierr)

#else
    return
#endif

end subroutine expand_J_field_3D

!###################################################################################################
end module parallel_subs_mod