module interactions_mod

#ifdef USE_MPI
    use mpi
#endif

    use constants_mod
    use source_mod
    use mxll_base_mod
    use mxll_1D_mod
    use mxll_2D_mod
    use mxll_3D_mod
    use q_group_mod

    implicit none

contains

!###################################################################################################

subroutine source_interactions(mxll, source_list, n_sources)
    class(TMxll) ,intent(inout) :: mxll
    type(TSource),intent(in)    :: source_list(n_sources)
    integer      ,intent(in)    :: n_sources

    integer  :: i, j, k, s
    integer  :: i_id, j_id, k_id
    integer  :: n_ker
    real(dp) :: J_av
    real(dp) :: c_src

    c_src = mxll%dt_eps0/mxll%dr/c0/2.0d0

    select type(mxll)
    class is(TMxll_1D)

#ifdef USE_MPI
        if (.not. allocated(mxll%Ex)) return
#endif

        do s=1, n_sources
            n_ker = source_list(s)%n_ker
            select case (source_list(s)%dir)
            case ('x')
                do i=-n_ker, n_ker
                    i_id = source_list(s)%ind_i(i,1,1)
                    J_av = source_list(s)%J_mat(i,1,1)
                    mxll%Ex(i_id) = mxll%Ex(i_id) + c_src * J_av
                end do
            case default
                print *, "Error: Unknown source direction in 1D mxll."
                stop
            end select
        end do

    class is(TMxll_2D)
    
        do s=1, n_sources
            n_ker = source_list(s)%n_ker
            select case (source_list(s)%dir)
            case ('x')
                do i=-n_ker, n_ker
                do j=-n_ker, n_ker
                    if (source_list(s)%in_this_rank(i,j,1)) then
                        i_id = source_list(s)%ind_i(i,j,1)
                        j_id = source_list(s)%ind_j(i,j,1)
                        J_av = 0.5d0*(source_list(s)%J_mat(i,j,1) + &
                                        source_list(s)%J_mat(i+1,j,1))
                        mxll%Ex(i_id,j_id) = mxll%Ex(i_id,j_id) + c_src * J_av
                    end if
                end do
                end do
            case ('y')
                do i=-n_ker, n_ker
                do j=-n_ker, n_ker
                    if (source_list(s)%in_this_rank(i,j,1)) then
                        i_id = source_list(s)%ind_i(i,j,1)
                        j_id = source_list(s)%ind_j(i,j,1)
                        J_av = 0.5d0*(source_list(s)%J_mat(i,j,1) + &
                                        source_list(s)%J_mat(i,j+1,1))
                        mxll%Ey(i_id,j_id) = mxll%Ey(i_id,j_id) + c_src * J_av
                    end if
                end do
                end do
            case ('z')
                do i=-n_ker, n_ker
                do j=-n_ker, n_ker
                    if (source_list(s)%in_this_rank(i,j,1)) then
                        i_id = source_list(s)%ind_i(i,j,1)
                        j_id = source_list(s)%ind_j(i,j,1)
                        J_av = source_list(s)%J_mat(i,j,1)
                        mxll%Ez(i_id,j_id) = mxll%Ez(i_id,j_id) + c_src * J_av
                    end if
                end do
                end do
            case default
                print *, "Error: Unknown source direction in 2D mxll."
                stop
            end select
        end do

    class is(TMxll_3D)
    
        do s=1, n_sources
            n_ker = source_list(s)%n_ker
            select case (source_list(s)%dir)
            case ('x')
                do i=-n_ker, n_ker
                do j=-n_ker, n_ker
                do k=-n_ker, n_ker
                    if (source_list(s)%in_this_rank(i,j,k)) then
                        i_id = source_list(s)%ind_i(i,j,k)
                        j_id = source_list(s)%ind_j(i,j,k)
                        k_id = source_list(s)%ind_k(i,j,k)
                        J_av = 0.5d0*(source_list(s)%J_mat(i,j,k) + &
                                        source_list(s)%J_mat(i+1,j,k))
                        mxll%Ex(i_id,j_id,k_id) = mxll%Ex(i_id,j_id,k_id) + c_src * J_av
                    end if
                end do
                end do
                end do
            case ('y')
                do i=-n_ker, n_ker
                do j=-n_ker, n_ker
                do k=-n_ker, n_ker
                    if (source_list(s)%in_this_rank(i,j,k)) then
                        i_id = source_list(s)%ind_i(i,j,k)
                        j_id = source_list(s)%ind_j(i,j,k)
                        k_id = source_list(s)%ind_k(i,j,k)
                        J_av = 0.5d0*(source_list(s)%J_mat(i,j,k) + &
                                        source_list(s)%J_mat(i,j+1,k))
                        mxll%Ey(i_id,j_id,k_id) = mxll%Ey(i_id,j_id,k_id) + c_src * J_av
                    end if
                end do
                end do
                end do
            case ('z')
                do i=-n_ker, n_ker
                do j=-n_ker, n_ker
                do k=-n_ker, n_ker
                    if (source_list(s)%in_this_rank(i,j,k)) then
                        i_id = source_list(s)%ind_i(i,j,k)
                        j_id = source_list(s)%ind_j(i,j,k)
                        k_id = source_list(s)%ind_k(i,j,k)
                        J_av = 0.5d0*(source_list(s)%J_mat(i,j,k) + &
                                        source_list(s)%J_mat(i,j,k+1))
                        mxll%Ez(i_id,j_id,k_id) = mxll%Ez(i_id,j_id,k_id) + c_src * J_av
                    end if
                end do
                end do
                end do
            case default
                print *, "Error: Unknown source direction in 3D mxll."
                stop
            end select
        end do

    class default
        print *, "Error: Unknown mxll type in source_interactions."
        stop
    end select

end subroutine source_interactions

!###################################################################################################
    
subroutine send_E_to_J_ranks(mxll, q_group, move_q_system, myrank)

    class(TMxll)  , intent(inout) :: mxll
    type(TQ_Group), intent(inout) :: q_group
    logical       , intent(in)    :: move_q_system
    integer       , intent(in)    :: myrank
    integer :: ierr

    if (.not. move_q_system) return

#ifdef USE_MPI
    call MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif

    select type(mxll)
    class is(TMxll_1D)

        call send_E_1D_to_J_ranks(mxll, q_group, myrank)

    class is(TMxll_2D)

        call send_E_2D_to_J_ranks(mxll, q_group, myrank)
    class is(TMxll_3D)

        call send_E_3D_to_J_ranks(mxll, q_group, myrank)
    end select

#ifdef USE_MPI
    call MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif

end subroutine send_E_to_J_ranks

!###################################################################################################

subroutine send_J_to_E_ranks(mxll, q_group, move_q_system, myrank)

    class(TMxll),  intent(inout)  :: mxll
    type(TQ_Group), intent(inout) :: q_group
    logical       , intent(in)    :: move_q_system
    integer       , intent(in)    :: myrank
    integer :: ierr

    if (.not. move_q_system) return

#ifdef USE_MPI
    call MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif

    select type(mxll)
    class is(TMxll_1D)

        call send_J_to_E_1D_ranks(mxll, q_group, myrank)

    class is(TMxll_2D)

        call send_J_to_E_2D_ranks(mxll, q_group, myrank)
    class is(TMxll_3D)

        call send_J_to_E_3D_ranks(mxll, q_group, myrank)
    end select

#ifdef USE_MPI
    call MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif

end subroutine send_J_to_E_ranks

!###################################################################################################

subroutine send_E_1D_to_J_ranks(mxll, q_group, myrank)
    class(TMxll_1D),  intent(inout)  :: mxll
    type(TQ_Group),   intent(inout)  :: q_group
    integer       , intent(in)       :: myrank

    integer :: n
    integer :: i_idx
    integer :: n_mol
    real(dp) :: E_field_send
    real(dp) :: E_field_get

#ifdef USE_MPI

    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
    
    select case(q_group%group_type)

    case(Q_MATERIAL)

        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%map(n,3)

            if(q_group%map(n,1) == myrank .and. q_group%map(n,2) == myrank) then
                
                q_group%E_field_list(n_mol, 1) = mxll%Ex(i_idx)
                n_mol = n_mol + 1

            else

                if (q_group%map(n,2) == myrank) then

                    E_field_send = mxll%Ex(i_idx)
                    call mpi_send(E_field_send,1,mpi_double_precision, &
                                  q_group%map(n,1), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)


                else if (q_group%map(n,1) == myrank) then

                    call mpi_recv(E_field_get,1,mpi_double_precision, &
                                  q_group%map(n,2), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    q_group%E_field_list(n_mol, 1) = E_field_get
                    n_mol = n_mol + 1

                end if
            end if
        end do

    case(Q_SINGLE)

        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%kernel_map(n,0,0,0,3)

            if(q_group%kernel_map(n,0,0,0,1) == myrank .and. q_group%kernel_map(n,0,0,0,2) == myrank) then
            
                q_group%E_field_list(n_mol, 1) = mxll%Ex(i_idx)
                n_mol = n_mol + 1

            else

                if (q_group%kernel_map(n,0,0,0,2) == myrank) then

                    E_field_send = mxll%Ex(i_idx)

                    call mpi_send(E_field_send,1,mpi_double_precision, &
                                  q_group%kernel_map(n,0,0,0,1), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)


                else if (q_group%kernel_map(n,0,0,0,1) == myrank) then

                    call mpi_recv(E_field_get,1,mpi_double_precision, &
                                  q_group%kernel_map(n,0,0,0,2), MPI_GOOD_TAG, MPI_COMM_WORLD, &
                                  istatus,ierr)

                    q_group%E_field_list(n_mol, 1) = E_field_get
                    n_mol = n_mol + 1

                end if
            end if
        end do

    end select

#else

    select case(q_group%group_type)

    case(Q_MATERIAL)

        do n = 1, q_group%n_systems
            n_mol = n
            i_idx = q_group%map(n,3)
            q_group%E_field_list(n_mol, 1) = mxll%Ex(i_idx)
        end do

    case(Q_SINGLE)

        do n = 1, q_group%n_systems
            n_mol = n
            i_idx = q_group%kernel_map(n,0,0,0,3)
            q_group%E_field_list(n_mol, 1) = mxll%Ex(i_idx)
        end do
    end select

#endif

end subroutine send_E_1D_to_J_ranks

!###################################################################################################

subroutine send_E_2D_to_J_ranks(mxll, q_group, myrank)
    class(TMxll_2D),  intent(inout) :: mxll
    type(TQ_Group) ,  intent(inout) :: q_group
    integer        , intent(in)     :: myrank
 
    integer :: n
    integer :: i_idx
    integer :: j_idx
    integer :: n_mol
    real(dp) :: E_field_send(3)
    real(dp) :: E_field_get(3)

#ifdef USE_MPI

    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)

    E_field_send = M_ZERO
    E_field_get  = M_ZERO

    select case(q_group%group_type)
    case(Q_MATERIAL)

        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            
            if(q_group%map(n,1) == myrank .and. q_group%map(n,2) == myrank) then
            
                if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                
                    q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx)+ &
                                                            mxll%Ex(i_idx, j_idx))
                    q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1)+ &
                                                            mxll%Ey(i_idx, j_idx))
                else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                    q_group%E_field_list(n_mol, 3) = mxll%Ez(i_idx, j_idx)

                end if
                
                n_mol = n_mol + 1

            else
                if (q_group%map(n,2) == myrank) then

                    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        E_field_send(1) =  0.5d0*(mxll%Ex(i_idx-1, j_idx)+ &
                                                    mxll%Ex(i_idx, j_idx))
                        E_field_send(2) =  0.5d0*(mxll%Ey(i_idx, j_idx-1)+ &
                                                    mxll%Ey(i_idx, j_idx))
                    else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        E_field_send(3) = mxll%Ez(i_idx, j_idx)
                    end if

                    call mpi_send(E_field_send,3,mpi_double_precision, &
                                    q_group%map(n,1), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)

                else if (q_group%map(n,1) == myrank) then

                    call mpi_recv(E_field_get,3,mpi_double_precision, &
                                    q_group%map(n,2), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    q_group%E_field_list(n_mol, 1) = E_field_get(1)
                    q_group%E_field_list(n_mol, 2) = E_field_get(2)
                    q_group%E_field_list(n_mol, 3) = E_field_get(3) 
                    n_mol = n_mol + 1

                end if
            end if
        end do

    case(Q_SINGLE)

        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%kernel_map(n,0,0,0,3)
            j_idx = q_group%kernel_map(n,0,0,0,4)
            
                if(q_group%kernel_map(n,0,0,0,1) == myrank .and. &
                    q_group%kernel_map(n,0,0,0,2) == myrank) then
            
                if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                
                    q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx)+ &
                                                            mxll%Ex(i_idx, j_idx))
                    q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1)+ &
                                                            mxll%Ey(i_idx, j_idx))

                else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                    q_group%E_field_list(n_mol, 3) = mxll%Ez(i_idx, j_idx)
                end if
                
                n_mol = n_mol + 1

            else
                if (q_group%kernel_map(n,0,0,0,2) == myrank) then

                    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        E_field_send(1) =  0.5d0*(mxll%Ex(i_idx-1, j_idx)+ &
                                                    mxll%Ex(i_idx, j_idx))
                        E_field_send(2) =  0.5d0*(mxll%Ey(i_idx, j_idx-1)+ &
                                                    mxll%Ey(i_idx, j_idx))
                    else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        E_field_send(3) = mxll%Ez(i_idx, j_idx)
                    end if
                
                    call mpi_send(E_field_send,3,mpi_double_precision, &
                                    q_group%kernel_map(n,0,0,0,1), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)
                
                else if (q_group%kernel_map(n,0,0,0,1) == myrank) then
                    
                    call mpi_recv(E_field_get,3,mpi_double_precision, &
                                    q_group%kernel_map(n,0,0,0,2), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    q_group%E_field_list(n_mol, 1) = E_field_get(1)
                    q_group%E_field_list(n_mol, 2) = E_field_get(2)
                    q_group%E_field_list(n_mol, 3) = E_field_get(3) 
                    n_mol = n_mol + 1

                end if
            end if
        end do

    end select

#else

    select case(q_group%group_type)
    case(Q_MATERIAL)

        do n = 1, q_group%n_systems
            
            n_mol = n
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            
            if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                
                q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx)+ &
                                                        mxll%Ex(i_idx, j_idx))
                q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1)+ &
                                                        mxll%Ey(i_idx, j_idx))
            else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                q_group%E_field_list(n_mol, 3) = mxll%Ez(i_idx, j_idx)

            end if
        end do

    case(Q_SINGLE)

        do n = 1, q_group%n_systems
            
            n_mol = n
            i_idx = q_group%kernel_map(n,0,0,0,3)
            j_idx = q_group%kernel_map(n,0,0,0,4)
            
            if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                
                q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx)+ &
                                                        mxll%Ex(i_idx, j_idx))
                q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1)+ &
                                                        mxll%Ey(i_idx, j_idx))
            else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                q_group%E_field_list(n_mol, 3) = mxll%Ez(i_idx, j_idx)
            end if
        end do
    
    end select

#endif

end subroutine send_E_2D_to_J_ranks

!###################################################################################################

subroutine send_E_3D_to_J_ranks(mxll, q_group, myrank)
    class(TMxll_3D),  intent(inout) :: mxll
    type(TQ_Group) ,  intent(inout) :: q_group
    integer        , intent(in)     :: myrank
 
    integer :: n
    integer :: i_idx
    integer :: j_idx
    integer :: k_idx
    integer :: n_mol
    real(dp) :: E_field_send(3)
    real(dp) :: E_field_get(3)
    
#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
    

    E_field_send = M_ZERO
    E_field_get  = M_ZERO

    select case(q_group%group_type)
    case(Q_MATERIAL)

        n_mol = 1
        do n = 1, q_group%n_systems
            
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            k_idx = q_group%map(n,5)
            
            if(q_group%map(n,1) == myrank .and. q_group%map(n,2) == myrank) then
            
                q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx, k_idx)+ &
                                                        mxll%Ex(i_idx, j_idx, k_idx))
                q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1, k_idx)+ &
                                                        mxll%Ey(i_idx, j_idx, k_idx))
                q_group%E_field_list(n_mol, 3) = 0.5d0*(mxll%Ez(i_idx, j_idx, k_idx-1)+ &
                                                        mxll%Ez(i_idx, j_idx, k_idx))
                n_mol = n_mol + 1

            else
                
                if (q_group%map(n,2) == myrank) then

                    E_field_send(1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx, k_idx)+ &
                                             mxll%Ex(i_idx, j_idx, k_idx))
                    E_field_send(2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1, k_idx)+ &
                                             mxll%Ey(i_idx, j_idx, k_idx))
                    E_field_send(3) = 0.5d0*(mxll%Ez(i_idx, j_idx, k_idx-1)+ &
                                             mxll%Ez(i_idx, j_idx, k_idx))

                    call mpi_send(E_field_send,3,mpi_double_precision, &
                                  q_group%map(n,1), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)

                else if (q_group%map(n,1) == myrank) then

                    call mpi_recv(E_field_get,3,mpi_double_precision, &
                                  q_group%map(n,2), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    q_group%E_field_list(n_mol, 1) = E_field_get(1)
                    q_group%E_field_list(n_mol, 2) = E_field_get(2)
                    q_group%E_field_list(n_mol, 3) = E_field_get(3) 
                    n_mol = n_mol + 1

                end if

            end if
        end do

    case(Q_SINGLE)

        n_mol = 1
        do n = 1, q_group%n_systems
            
            i_idx = q_group%kernel_map(n,0,1,1,3)
            j_idx = q_group%kernel_map(n,0,1,1,4)
            k_idx = q_group%kernel_map(n,0,1,1,5)
            
                if(q_group%kernel_map(n,0,0,0,1) == myrank .and. &
                    q_group%kernel_map(n,0,0,0,2) == myrank) then
            
                q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx, k_idx)+ &
                                                        mxll%Ex(i_idx, j_idx, k_idx))
                q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1, k_idx)+ &
                                                        mxll%Ey(i_idx, j_idx, k_idx))
                q_group%E_field_list(n_mol, 3) = 0.5d0*(mxll%Ez(i_idx, j_idx, k_idx-1)+ &
                                                        mxll%Ez(i_idx, j_idx, k_idx))
                n_mol = n_mol + 1

            else
                
                if (q_group%kernel_map(n,0,0,0,2) == myrank) then

                    E_field_send(1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx, k_idx)+ &
                                             mxll%Ex(i_idx, j_idx, k_idx))
                    E_field_send(2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1, k_idx)+ &
                                             mxll%Ey(i_idx, j_idx, k_idx))
                    E_field_send(3) = 0.5d0*(mxll%Ez(i_idx, j_idx, k_idx-1)+ &
                                             mxll%Ez(i_idx, j_idx, k_idx))

                    call mpi_send(E_field_send,3,mpi_double_precision, &
                                  q_group%kernel_map(n,0,0,0,1), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)

                else if (q_group%kernel_map(n,0,0,0,1) == myrank) then

                    call mpi_recv(E_field_get,3,mpi_double_precision, &
                                  q_group%kernel_map(n,0,0,0,2), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    q_group%E_field_list(n_mol, 1) = E_field_get(1)
                    q_group%E_field_list(n_mol, 2) = E_field_get(2)
                    q_group%E_field_list(n_mol, 3) = E_field_get(3) 
                    n_mol = n_mol + 1

                end if
            end if
        end do

    end select  

#else

    select case(q_group%group_type)
    case(Q_MATERIAL)

        do n = 1, q_group%n_systems
            
            n_mol = n
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            k_idx = q_group%map(n,5)
            
            q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx, k_idx)+ &
                                                    mxll%Ex(i_idx, j_idx, k_idx))
            q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1, k_idx)+ &
                                                    mxll%Ey(i_idx, j_idx, k_idx))
            q_group%E_field_list(n_mol, 3) = 0.5d0*(mxll%Ez(i_idx, j_idx, k_idx-1)+ &
                                                    mxll%Ez(i_idx, j_idx, k_idx))
        end do

    case(Q_SINGLE)

        do n = 1, q_group%n_systems
            
            n_mol = n
            i_idx = q_group%kernel_map(n,0,0,0,3)
            j_idx = q_group%kernel_map(n,0,0,0,4)
            k_idx = q_group%kernel_map(n,0,0,0,5)
            
            q_group%E_field_list(n_mol, 1) = 0.5d0*(mxll%Ex(i_idx-1, j_idx, k_idx)+ &
                                                    mxll%Ex(i_idx, j_idx, k_idx))
            q_group%E_field_list(n_mol, 2) = 0.5d0*(mxll%Ey(i_idx, j_idx-1, k_idx)+ &
                                                    mxll%Ey(i_idx, j_idx, k_idx))
            q_group%E_field_list(n_mol, 3) = 0.5d0*(mxll%Ez(i_idx, j_idx, k_idx-1)+ &
                                                    mxll%Ez(i_idx, j_idx, k_idx))
        end do
    
    end select

#endif

end subroutine send_E_3D_to_J_ranks

!###################################################################################################

subroutine send_J_to_E_1D_ranks(mxll, q_group, myrank)
    class(TMxll_1D),  intent(inout)  :: mxll
    type(TQ_Group),   intent(inout)  :: q_group
    integer       , intent(in)       :: myrank

    integer :: n, i
    integer :: i_idx
    integer :: n_mol
    logical  :: move_mol
    real(dp) :: J_field_send
    real(dp) :: J_field_get

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)


    if (myrank == 0) then
        mxll%Jx_old = mxll%Jx
    end if

#else
    mxll%Jx_old = mxll%Jx
#endif

    
#ifdef USE_MPI    
    

    select case(q_group%group_type)

    case(Q_MATERIAL)

        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%map(n,3)


            if(q_group%map(n,1) == myrank .and. q_group%map(n,2) == myrank) then

                mxll%Jx(i_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
                n_mol = n_mol + 1

            else

                if (q_group%map(n,1) == myrank) then
                    J_field_send = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)

                    call mpi_send(J_field_send,1,mpi_double_precision, &
                                  q_group%map(n,2), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)

                    n_mol = n_mol + 1

                else if (q_group%map(n,2) == myrank) then
                    call mpi_recv(J_field_get,1,mpi_double_precision, &
                                  q_group%map(n,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    mxll%Jx(i_idx) = J_field_get

                end if
            end if
        end do

    case(Q_SINGLE)

        n_mol = 1

        do n = 1, q_group%n_systems
            move_mol = .false.
            do i = -q_group%n_ker, q_group%n_ker
                
                i_idx = q_group%kernel_map(n,i,0,0,3)

                if(q_group%kernel_map(n,i,0,0,1) == myrank .and. q_group%kernel_map(n,i,0,0,2) == myrank) then
                
                    mxll%Jx(i_idx) = q_group%density*q_group%mat_kernel(i,1,1)*  &
                                     q_group%q_sys(n_mol)%dPt_dt(1)
                    move_mol = .true.
                else

                    if (q_group%kernel_map(n,i,0,0,1) == myrank) then

                        J_field_send = q_group%density*q_group%mat_kernel(i,1,1)*  &
                                       q_group%q_sys(n_mol)%dPt_dt(1)

                        call mpi_send(J_field_send,1,mpi_double_precision, &
                                    q_group%kernel_map(n,i,0,0,2), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)
                        move_mol = .true.

                    else if (q_group%kernel_map(n,i,0,0,2) == myrank) then
                        
                        call mpi_recv(J_field_get,1,mpi_double_precision, &
                                    q_group%kernel_map(n,i,0,0,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)
                        mxll%Jx(i_idx) = J_field_get
                        
                    end if
                end if
            end do
            if (move_mol) n_mol = n_mol + 1
        end do
    end select

    if (myrank == 0) then
        mxll%dJx = (mxll%Jx - mxll%Jx_old)/q_group%dt
    end if

#else

    select case(q_group%group_type)

    case(Q_MATERIAL)

        do n = 1, q_group%n_systems
            n_mol = n
            i_idx = q_group%map(n,3)
            mxll%Jx(i_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
        end do

    case(Q_SINGLE)

        do n = 1, q_group%n_systems
        do i = -q_group%n_ker, q_group%n_ker
            n_mol = n
            i_idx = q_group%kernel_map(n,i,0,0,3)
            mxll%Jx(i_idx) = q_group%density*q_group%mat_kernel(i,1,1)*  &
                             q_group%q_sys(n_mol)%dPt_dt(1)
        end do
        end do
    end select
    
    mxll%dJx = (mxll%Jx - mxll%Jx_old)/q_group%dt

#endif


end subroutine send_J_to_E_1D_ranks

!###################################################################################################

subroutine send_J_to_E_2D_ranks(mxll, q_group, myrank)
    class(TMxll_2D),  intent(inout) :: mxll
    type(TQ_Group) ,  intent(inout) :: q_group
    integer        , intent(in)     :: myrank
 
    integer :: n
    integer :: i, j
    integer :: i_idx
    integer :: j_idx
    integer :: n_mol
    logical  :: move_mol
    real(dp) :: J_field_send(3)
    real(dp) :: J_field_get(3)

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
#endif

    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
        mxll%Jx_old = mxll%Jx
        mxll%Jy_old = mxll%Jy
    end if
    if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
        mxll%Jz_old = mxll%Jz
    end if

#ifdef USE_MPI

    select case(q_group%group_type)
    case(Q_MATERIAL)
        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            
            if(q_group%map(n,1) == myrank .and. q_group%map(n,2) == myrank) then
            
                if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                    mxll%Jx(i_idx, j_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
                    mxll%Jy(i_idx, j_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(2)
                end if
                if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                    mxll%Jz(i_idx, j_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(3)
                end if
                
                n_mol = n_mol + 1

            else
            
                if (q_group%map(n,1) == myrank) then

                    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        J_field_send(1) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
                        J_field_send(2) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(2)
                    end if
                    if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        J_field_send(3) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(3)
                    end if

                    call mpi_send(J_field_send,3,mpi_double_precision, &
                                    q_group%map(n,2), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)

                    n_mol = n_mol + 1

                else if (q_group%map(n,2) == myrank) then

                    call mpi_recv(J_field_get,3,mpi_double_precision, &
                                    q_group%map(n,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        mxll%Jx(i_idx, j_idx) = J_field_get(1)
                        mxll%Jy(i_idx, j_idx) = J_field_get(2)
                    end if
                    if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        mxll%Jz(i_idx, j_idx) = J_field_get(3)
                    end if
                end if

            end if
        end do

    case(Q_SINGLE)

        n_mol = 1

        do n = 1, q_group%n_systems
            move_mol = .false.
            do j = -q_group%n_ker, q_group%n_ker
            do i = -q_group%n_ker, q_group%n_ker
                
                i_idx = q_group%kernel_map(n,i,j,0,3)
                j_idx = q_group%kernel_map(n,i,j,0,4)

                     if(q_group%kernel_map(n,i,j,0,1) == myrank .and. &
                         q_group%kernel_map(n,i,j,0,2) == myrank) then
                
                    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        mxll%Jx(i_idx, j_idx) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                                q_group%q_sys(n_mol)%dPt_dt(1)
                        mxll%Jy(i_idx, j_idx) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                                q_group%q_sys(n_mol)%dPt_dt(2)
                    end if
                    if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                        mxll%Jz(i_idx, j_idx) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                                q_group%q_sys(n_mol)%dPt_dt(3)
                    end if
                    move_mol = .true.
                else

                    if (q_group%kernel_map(n,i,j,0,1) == myrank) then

                        if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                            J_field_send(1) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                              q_group%q_sys(n_mol)%dPt_dt(1)
                            J_field_send(2) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                              q_group%q_sys(n_mol)%dPt_dt(2)
                        end if
                        if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                            J_field_send(3) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                              q_group%q_sys(n_mol)%dPt_dt(3)
                        end if
                        call mpi_send(J_field_send,3,mpi_double_precision, &
                                    q_group%kernel_map(n,i,j,0,2), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)
                        move_mol = .true.

                    else if (q_group%kernel_map(n,i,j,0,2) == myrank) then
                        
                         call mpi_recv(J_field_get,3,mpi_double_precision, &
                                    q_group%kernel_map(n,i,j,0,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                        if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                            mxll%Jx(i_idx, j_idx) = J_field_get(1)
                            mxll%Jy(i_idx, j_idx) = J_field_get(2)
                        end if
                        if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                            mxll%Jz(i_idx, j_idx) = J_field_get(3)
                        end if
                        
                    end if
                end if
            end do
            end do
            if (move_mol) n_mol = n_mol + 1
        end do
    end select

#else

    select case(q_group%group_type)
    case(Q_MATERIAL)

        do n = 1, q_group%n_systems
            n_mol = n
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            
            if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                mxll%Jx(i_idx, j_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
                mxll%Jy(i_idx, j_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(2)
            end if
            if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                mxll%Jz(i_idx, j_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(3)
            end if
                
        end do

    case(Q_SINGLE)

        do n = 1, q_group%n_systems
        do j = -q_group%n_ker, q_group%n_ker
        do i = -q_group%n_ker, q_group%n_ker
            n_mol = n
            i_idx = q_group%kernel_map(n,i,j,0,3)
            j_idx = q_group%kernel_map(n,i,j,0,4)

            if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                mxll%Jx(i_idx, j_idx) = q_group%density*q_group%mat_kernel(i,j,1)* &
                                        q_group%q_sys(n_mol)%dPt_dt(1)
                mxll%Jy(i_idx, j_idx) = q_group%density*q_group%mat_kernel(i,j,1)* &
                                        q_group%q_sys(n_mol)%dPt_dt(2)
            end if
            if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
                mxll%Jz(i_idx, j_idx) = q_group%density*q_group%mat_kernel(i,j,1)*  &
                                        q_group%q_sys(n_mol)%dPt_dt(3)
            end if
        end do
        end do
        end do
    end select

#endif

    if (mxll%mode == TEZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
        mxll%dJx = (mxll%Jx - mxll%Jx_old)/q_group%dt
        mxll%dJy = (mxll%Jy - mxll%Jy_old)/q_group%dt
    else if (mxll%mode == TMZ_2D_MODE .or. mxll%mode == FULL_2D_MODE) then
        mxll%dJz = (mxll%Jz - mxll%Jz_old)/q_group%dt
    end if

end subroutine send_J_to_E_2D_ranks

!###################################################################################################

subroutine send_J_to_E_3D_ranks(mxll, q_group, myrank)
    class(TMxll_3D),  intent(inout) :: mxll
    type(TQ_Group) ,  intent(inout) :: q_group
    integer        , intent(in)     :: myrank
 
    integer :: n
    integer :: i, j, k
    integer :: i_idx
    integer :: j_idx
    integer :: k_idx
    integer :: n_mol
    logical  :: move_mol
    real(dp) :: J_field_send(3)
    real(dp) :: J_field_get(3)

#ifdef USE_MPI
    integer :: ierr
    integer :: istatus(MPI_STATUS_SIZE)
#endif


    mxll%Jx_old = mxll%Jx
    mxll%Jy_old = mxll%Jy
    mxll%Jz_old = mxll%Jz

#ifdef USE_MPI

    select case(q_group%group_type)
    case(Q_MATERIAL)

        n_mol = 1

        do n = 1, q_group%n_systems
            
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            k_idx = q_group%map(n,5)
            
            if(q_group%map(n,1) == myrank .and. q_group%map(n,2) == myrank) then
            
                mxll%Jx(i_idx, j_idx, k_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
                mxll%Jy(i_idx, j_idx, k_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(2)
                mxll%Jz(i_idx, j_idx, k_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(3)
                
                n_mol = n_mol + 1

            else
            
                if (q_group%map(n,1) == myrank) then

                    J_field_send(1) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
                    J_field_send(2) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(2)
                    J_field_send(3) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(3)

                    call mpi_send(J_field_send,3,mpi_double_precision, &
                                  q_group%map(n,2), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)

                    n_mol = n_mol + 1

                else if (q_group%map(n,2) == myrank) then

                    call mpi_recv(J_field_get,3,mpi_double_precision, &
                                  q_group%map(n,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                    mxll%Jx(i_idx, j_idx, k_idx) = J_field_get(1)
                    mxll%Jy(i_idx, j_idx, k_idx) = J_field_get(2)
                    mxll%Jz(i_idx, j_idx, k_idx) = J_field_get(3)

                end if

            end if
        end do

    case(Q_SINGLE)

        n_mol = 1

        do n = 1, q_group%n_systems
            move_mol = .false.
            do k = -q_group%n_ker, q_group%n_ker
            do j = -q_group%n_ker, q_group%n_ker
            do i = -q_group%n_ker, q_group%n_ker
                
                i_idx = q_group%kernel_map(n,i,j,k,3)
                j_idx = q_group%kernel_map(n,i,j,k,4)
                k_idx = q_group%kernel_map(n,i,j,k,5)

                     if(q_group%kernel_map(n,i,j,k,1) == myrank .and. &
                         q_group%kernel_map(n,i,j,k,2) == myrank) then
                
                    mxll%Jx(i_idx, j_idx, k_idx) = q_group%density*&
                                                   q_group%mat_kernel(i,j,k)*&
                                                   q_group%q_sys(n_mol)%dPt_dt(1)
                    mxll%Jy(i_idx, j_idx, k_idx) = q_group%density*&
                                                   q_group%mat_kernel(i,j,k)*&
                                                   q_group%q_sys(n_mol)%dPt_dt(2)
                    mxll%Jz(i_idx, j_idx, k_idx) = q_group%density*&
                                                   q_group%mat_kernel(i,j,k)*&
                                                   q_group%q_sys(n_mol)%dPt_dt(3)
                    move_mol = .true.
                else

                    if (q_group%kernel_map(n,i,j,k,1) == myrank) then

                        J_field_send(1) = q_group%density*q_group%mat_kernel(i,j,k)* &
                                          q_group%q_sys(n_mol)%dPt_dt(1)
                        J_field_send(2) = q_group%density*q_group%mat_kernel(i,j,k)* &
                                          q_group%q_sys(n_mol)%dPt_dt(2)
                        J_field_send(3) = q_group%density*q_group%mat_kernel(i,j,k)* &
                                          q_group%q_sys(n_mol)%dPt_dt(3)

                        call mpi_send(J_field_send,3,mpi_double_precision, &
                                    q_group%kernel_map(n,i,j,k,2), MPI_GOOD_TAG, MPI_COMM_WORLD,ierr)
                        move_mol = .true.

                    else if (q_group%kernel_map(n,i,j,k,2) == myrank) then
                        
                         call mpi_recv(J_field_get,3,mpi_double_precision, &
                                    q_group%kernel_map(n,i,j,k,1), MPI_GOOD_TAG, MPI_COMM_WORLD,istatus,ierr)

                        mxll%Jx(i_idx, j_idx, k_idx) = J_field_get(1)
                        mxll%Jy(i_idx, j_idx, k_idx) = J_field_get(2)
                        mxll%Jz(i_idx, j_idx, k_idx) = J_field_get(3)
                    end if
                end if
            end do
            end do
            end do
            if (move_mol) n_mol = n_mol + 1
        end do
    end select
#else
    select case(q_group%group_type)
    case(Q_MATERIAL)

        do n = 1, q_group%n_systems
            n_mol = n
            i_idx = q_group%map(n,3)
            j_idx = q_group%map(n,4)
            k_idx = q_group%map(n,5)
            
            mxll%Jx(i_idx, j_idx, k_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(1)
            mxll%Jy(i_idx, j_idx, k_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(2)
            mxll%Jz(i_idx, j_idx, k_idx) = q_group%density*q_group%q_sys(n_mol)%dPt_dt(3)
        end do

    case(Q_SINGLE)

        do n = 1, q_group%n_systems
        do k = -q_group%n_ker, q_group%n_ker
        do j = -q_group%n_ker, q_group%n_ker
        do i = -q_group%n_ker, q_group%n_ker
            n_mol = n
            i_idx = q_group%kernel_map(n,i,j,k,3)
            j_idx = q_group%kernel_map(n,i,j,k,4)
            k_idx = q_group%kernel_map(n,i,j,k,5)

            mxll%Jx(i_idx, j_idx, k_idx) = q_group%density*q_group%mat_kernel(i,j,k)*  &
                                           q_group%q_sys(n_mol)%dPt_dt(1)
            mxll%Jy(i_idx, j_idx, k_idx) = q_group%density*q_group%mat_kernel(i,j,k)*  &
                                           q_group%q_sys(n_mol)%dPt_dt(2)
            mxll%Jz(i_idx, j_idx, k_idx) = q_group%density*q_group%mat_kernel(i,j,k)*  &
                                           q_group%q_sys(n_mol)%dPt_dt(3)
        end do
        end do
        end do
        end do
    end select
#endif

    mxll%dJx = (mxll%Jx - mxll%Jx_old)/q_group%dt
    mxll%dJy = (mxll%Jy - mxll%Jy_old)/q_group%dt
    mxll%dJz = (mxll%Jz - mxll%Jz_old)/q_group%dt

end subroutine send_J_to_E_3D_ranks

!###################################################################################################

end module interactions_mod