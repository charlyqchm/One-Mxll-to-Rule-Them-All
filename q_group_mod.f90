module q_group_mod

#ifdef USE_MPI
    use mpi
#endif
    use constants_mod
    use q_sys_base_mod
    use factory_mod, only: q_system_factory

    implicit none

    type TQ_group

        class(TQ_sys_base), allocatable :: q_sys(:)
        integer :: group_type
        integer :: q_sys_type
        integer :: n_systems
        integer :: n_sys_loc
        integer :: n_ker=3 !this value is harcoded for avoiding memeory issues.
        integer :: Nt_steps
        integer :: print_option
        real(dp):: dt
        real(dp):: density

        integer,  allocatable :: map(:,:)
        integer,  allocatable :: kernel_map(:,:,:,:,:)
        real(dp), allocatable :: mat_kernel(:,:,:)
        real(dp), allocatable :: E_field_list(:,:)

!Depending on the dimensions, map stores the following:
!           1D     |     2D      |    3D
!map(:,1) : J_rank |   J_rank    |   J_rank
!map(:,2) : E_rank |   E_rank    |   E_rank
!map(:,3) : i_loc  |   i_loc     |   i_loc
!map(:,4) :        |   j_loc     |   j_loc
!map(:,5) :        |             |   k_loc
        
!J_rank: Rank that calculates the molecule's J field.
!E_rank: Rank where the molecule's E field is located.
!i_loc, j_loc, k_loc: Local grid indices where the molecule is located in the E_rank.

!The last index of kernel_map stores the same information for each point in the kernel.

    contains
        procedure :: init_q_group
        procedure :: init_all_q_systems
        procedure :: kill_q_group
        procedure :: do_mapping
        procedure :: do_kernel_mapping 
        procedure :: td_propagate_q_group

    end type TQ_group

contains
    
    subroutine init_q_group(this, group_id, dimensions, mpi_dims, mpi_coords,  &
                             mpi_cart_comm, myrank, dr, dt, Nt_steps, density, &
                             grid_Ndims)
        class(TQ_group), intent(inout) :: this
        integer        , intent(in)    :: group_id
        integer        , intent(in)    :: dimensions
        integer        , intent(in)    :: mpi_dims(3)
        integer        , intent(in)    :: mpi_coords(3)
        integer        , intent(in)    :: mpi_cart_comm
        integer        , intent(in)    :: myrank
        integer        , intent(in)    :: Nt_steps
        integer        , intent(in)    :: grid_Ndims(3)
        real(dp)       , intent(in)    :: dr
        real(dp)       , intent(in)    :: dt
        real(dp)       , intent(in)    :: density
        character(len=20)  :: group_type_str
        character(len=20)  :: q_sys_type_str
        character(len=20)  :: file_name = "mol_group_"
        character(len=20)  :: file_exten = ".in"
        character(len=20)  :: file_number
        character(len=50)  :: print_option_ch
        character(len=50)  :: input_name
        character(len=20)  :: print_on_ch
        real(dp)           :: r0
        real(dp)           :: r_max
        real(dp)           :: x, y, z
        real(dp)           :: norm
        real(dp)           :: r
        integer            :: ierr
        integer            :: funit
        integer            :: n_systems 
        integer            :: i, j, k
        integer            :: n_left
        integer            :: n_total_ranks
        integer            :: n_ker
        integer , allocatable :: mol_id(:)
        integer , allocatable :: id_file(:)
        logical , allocatable :: print_on(:)
        real(dp), allocatable :: grid_coord(:,:)

        this%dt = dt
        this%Nt_steps = Nt_steps
        this%density = density

        n_ker = this%n_ker

        write(file_number,'(I7.7)') group_id
        input_name = trim(file_name)//trim(file_number)//trim(file_exten)
   
        open (action='read', file=input_name, iostat=ierr, newunit=funit)

        if (ierr /= 0) then
            write(*,*) "Error: Could not open file ", trim(input_name)
            stop
        end if      

        read (unit=funit,fmt=*, iostat=ierr) group_type_str, q_sys_type_str, n_systems, &
                                             print_option_ch

        this%n_systems = n_systems
        n_total_ranks  = mpi_dims(1)*mpi_dims(2)*mpi_dims(3)
        this%n_sys_loc = int(n_systems / n_total_ranks)
        n_left         = mod(n_systems, n_total_ranks)


        if (myrank < n_left) then
            this%n_sys_loc = this%n_sys_loc + 1
        end if

        select case(trim(print_option_ch))
        case ("print_none")
            this%print_option = PRINT_NONE
        case ("print_all")
            this%print_option = PRINT_ALL
        case ("print_selected")
            this%print_option = PRINT_SELECTED
        case default
            write(*,*) "Error: Unknown print option for the q_system ", trim(print_option_ch)
            stop
        end select

        if (.not. allocated(print_on)) allocate(print_on(n_systems))
        print_on = .false.

        select case(trim(group_type_str))
        case ("material")
            this%group_type = Q_MATERIAL
            select case (dimensions)
            case (1)
                if (.not. allocated(this%map)) allocate(this%map(n_systems, 3))
            case (2)
                if (.not. allocated(this%map)) allocate(this%map(n_systems, 4))
            case (3)
                if (.not. allocated(this%map)) allocate(this%map(n_systems, 5))
            case default
                write(*,*) "Error: Unknown dimensions ", dimensions
                stop
            end select
        case ("single")
            this%group_type = Q_SINGLE
            select case (dimensions)
            case (1)
                if (.not. allocated(this%kernel_map)) then 
                    allocate(this%kernel_map(n_systems, -n_ker:n_ker, 0:0,0:0, 3))
                end if
                if (.not. allocated(this%mat_kernel)) then
                    allocate(this%mat_kernel(-n_ker:n_ker, 1, 1))
                end if

                r_max = dble(n_ker*2+4)
                r0    = r_max/2.0d0
                norm  = M_ZERO
                do i = -n_ker, n_ker
                    x = dble(i)
                    this%mat_kernel(i,1,1) = aBH(1)+ &
                          aBH(2)*DCOS(2.0*pi0*(x-r0)/r_max)+ &
                          aBH(3)*DCOS(2.0*pi0*2.0*(x-r0)/r_max)+ &
                          aBH(4)*DCOS(2.0*pi0*3.0*(x-r0)/r_max)

                    norm = norm + this%mat_kernel(i,1,1)
                end do

            case (2)
                if (.not. allocated(this%kernel_map)) then 
                    allocate(this%kernel_map(n_systems, -n_ker:n_ker, -n_ker:n_ker, 0:0, 4))
                end if
                if (.not. allocated(this%mat_kernel)) then
                    allocate(this%mat_kernel(-n_ker:n_ker, -n_ker:n_ker, 1))
                end if

                r_max = dble(n_ker*2 + 4)
                r0    = r_max/2.0d0
                norm  = M_ZERO

                do i = -n_ker, n_ker
                do j = -n_ker, n_ker
                    x = dble(i)
                    y = dble(j)

                    r = SQRT(x**2 + y**2)
                    this%mat_kernel(i,j,1) = aBH(1)+ &
                          aBH(2)*DCOS(2.0*pi0*(r-r0)/r_max)+ &
                          aBH(3)*DCOS(2.0*pi0*2.0*(r-r0)/r_max)+ &
                          aBH(4)*DCOS(2.0*pi0*3.0*(r-r0)/r_max)

                    if (this%mat_kernel(i,j,1) < M_ZERO) this%mat_kernel(i,j,1) = M_ZERO

                    norm = norm + this%mat_kernel(i,j,1)

                end do
                end do

            case (3)
                if (.not. allocated(this%kernel_map)) then 
                    allocate(this%kernel_map(n_systems,-n_ker:n_ker,-n_ker:n_ker,-n_ker:n_ker,5))                                                 
                end if
                if (.not. allocated(this%mat_kernel)) then
                    allocate(this%mat_kernel(-n_ker:n_ker,-n_ker:n_ker,-n_ker:n_ker))
                end if

                r_max = dble(n_ker*2+4)
                r0    = r_max/2.0d0
                norm  = M_ZERO

                do i = -n_ker, n_ker
                do j = -n_ker, n_ker
                do k = -n_ker, n_ker
                    x = dble(i)
                    y = dble(j)
                    z = dble(k)

                    r = SQRT(x**2 + y**2 + z**2)
                    this%mat_kernel(i,j,k) = aBH(1)+ &
                          aBH(2)*DCOS(2.0*pi0*(r-r0)/r_max)+ &
                          aBH(3)*DCOS(2.0*pi0*2.0*(r-r0)/r_max)+ &
                          aBH(4)*DCOS(2.0*pi0*3.0*(r-r0)/r_max)

                    if (this%mat_kernel(i,j,k) < M_ZERO) this%mat_kernel(i,j,k) = M_ZERO

                    norm = norm + this%mat_kernel(i,j,k)

                end do
                end do
                end do

            case default
                write(*,*) "Error: Unknown dimensions ", dimensions
                stop
            end select
        case default
            write(*,*) "Error: Unknown Q_group type ", trim(group_type_str)
            stop
        end select

        select case(trim(q_sys_type_str))
        case ("dftb")
            this%q_sys_type = Q_SYS_DFTB
        case default
            write(*,*) "Error: Unknown Q_sys type ", trim(q_sys_type_str)
            stop
        end select


        if (.not. allocated(mol_id))            allocate(mol_id(n_systems))
        if (.not. allocated(id_file))           allocate(id_file(n_systems))
        if (.not. allocated(grid_coord))        allocate(grid_coord(n_systems, dimensions))
        if (.not. allocated(this%E_field_list)) allocate(this%E_field_list(this%n_sys_loc, 3))

        this%E_field_list = M_ZERO

        if (this%print_option == PRINT_NONE) then
            print_on = .false.
            do i = 1, n_systems
                read (unit=funit,fmt=*, iostat=ierr) mol_id(i), id_file(i), grid_coord(i, :)
            end do
        else if (this%print_option == PRINT_ALL) then
            print_on = .true.
            do i = 1, n_systems
                read (unit=funit,fmt=*, iostat=ierr) mol_id(i), id_file(i), grid_coord(i, :)
            end do
        else if (this%print_option == PRINT_SELECTED) then
            do i = 1, n_systems
                read (unit=funit,fmt=*, iostat=ierr) mol_id(i), id_file(i), grid_coord(i, :), &
                                                     print_on_ch
                if (trim(print_on_ch) == "print_on") then
                    print_on(i) = .true.
                else if (trim(print_on_ch) == "print_off") then
                    print_on(i) = .false.
                else
                    write(*,*) "Error: Unknown print option for the q_system ", i
                    stop
                end if
            end do
        end if


        grid_coord = grid_coord * nm_to_au

        select case (this%group_type)
            
        case (Q_MATERIAL)
            call this%do_mapping(grid_coord, dimensions, mpi_dims, mpi_coords, &
            mpi_cart_comm, grid_Ndims, dr)
        case (Q_SINGLE)
            call this%do_kernel_mapping(grid_coord, dimensions, mpi_dims, mpi_coords, &
            mpi_cart_comm, grid_Ndims, dr)
        end select
       
        call this%init_all_q_systems(mol_id, id_file, myrank, print_on)
       
        if (allocated(mol_id))     deallocate(mol_id)
        if (allocated(id_file))    deallocate(id_file)
        if (allocated(grid_coord)) deallocate(grid_coord)
        if (allocated(print_on))   deallocate(print_on)

        close(funit)

    end subroutine init_q_group
!###################################################################################################
    subroutine kill_q_group(this)
        class(TQ_group), intent(inout) :: this

        integer :: i

        do i = 1, this%n_sys_loc
            select case (this%q_sys_type)
            case (Q_SYS_DFTB)
                call this%q_sys(i)%kill()
            case default
                write(*,*) "Error: Unknown Q_sys type in kill_q_group."
                stop
            end select
        end do
        
        if (allocated(this%q_sys))        deallocate(this%q_sys)
        if (allocated(this%map))          deallocate(this%map)
        if (allocated(this%kernel_map))   deallocate(this%kernel_map)
        if (allocated(this%mat_kernel))   deallocate(this%mat_kernel)
        if (allocated(this%E_field_list)) deallocate(this%E_field_list)

    end subroutine kill_q_group

    
!###################################################################################################

    subroutine do_mapping(this, grid_coord, dimensions, mpi_dims, mpi_coords, mpi_cart_comm, &
                          grid_Ndims,dr)
        class(TQ_group), intent(inout) :: this
        real(dp)       , intent(in)    :: dr
        real(dp)       , intent(in)    :: grid_coord(:,:)
        integer        , intent(in)    :: dimensions
        integer        , intent(in)    :: mpi_dims(3)
        integer        , intent(in)    :: mpi_coords(3)
        integer        , intent(in)    :: mpi_cart_comm
        integer        , intent(in)    :: grid_Ndims(3)

        integer :: rank_counter
        integer :: n_total_ranks
        integer :: n_base_loc
        integer :: n_left
        integer :: i, j, k, n
        integer :: global_i, global_j, global_k
        integer :: ierr
        integer :: rank_coor(3)
        integer :: aux_rank

        n_total_ranks = mpi_dims(1)*mpi_dims(2)*mpi_dims(3)
        n_base_loc    = int(this%n_systems / n_total_ranks)
        n_left        = mod(this%n_systems, n_total_ranks)

#ifdef USE_MPI

        !By default, we consider that the origin is where global_i, global_j and/or global_k
        !reach the half of the global grid size.

       select case(dimensions)
        case (1)
            rank_counter = 0
            n = 1 
            do i = 1, this%n_systems
                global_i  = int(grid_coord(i, 1)/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
                
                this%map(i,1) = rank_counter
                this%map(i,2) = 0
                this%map(i,3) = int(grid_coord(i, 1)/dr) + int(grid_Ndims(1)/2)
                
                if (rank_counter < n_left) then
                    if (n >= n_base_loc + 1) then
                        rank_counter = rank_counter + 1
                        n = 0
                    end if
                else
                    if (n >= n_base_loc) then
                        rank_counter = rank_counter + 1
                        n = 0
                    end if
                end if
                
                n = n + 1
            end do

        case (2)

            rank_counter = 0
            n = 1 

            do i = 1, this%n_systems
                global_i = int(grid_coord(i, 1)/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
                global_j = int(grid_coord(i, 2)/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)

                rank_coor(1) = int((global_i-1)/grid_Ndims(1))
                rank_coor(2) = int((global_j-1)/grid_Ndims(2))
                rank_coor(3) = 0

                call MPI_Cart_rank(mpi_cart_comm, rank_coor(:2), aux_rank, ierr)

                this%map(i,1) =  rank_counter
                this%map(i,2) =  aux_rank

                this%map(i,3) =  MOD(global_i, grid_Ndims(1))
                this%map(i,4) =  MOD(global_j, grid_Ndims(2))
                if (MOD(global_i, grid_Ndims(1)) == 0) this%map(i,3) =  grid_Ndims(1)
                if (MOD(global_j, grid_Ndims(2)) == 0) this%map(i,4) =  grid_Ndims(2)

                if (rank_counter < n_left) then
                    if (n >= n_base_loc + 1) then
                        rank_counter = rank_counter + 1
                        n = 0
                    end if
                else
                    if (n >= n_base_loc) then
                        rank_counter = rank_counter + 1
                        n = 0
                    end if
                end if

                n = n + 1
            end do

        case (3)

            rank_counter = 0
            n = 1
            do i = 1, this%n_systems
                global_i = int(grid_coord(i, 1)/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
                global_j = int(grid_coord(i, 2)/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
                global_k = int(grid_coord(i, 3)/dr) + int(grid_Ndims(3)*mpi_dims(3)/2)

                rank_coor(1) = int((global_i-1)/grid_Ndims(1))
                rank_coor(2) = int((global_j-1)/grid_Ndims(2))
                rank_coor(3) = int((global_k-1)/grid_Ndims(3))

                call MPI_Cart_rank(mpi_cart_comm, rank_coor, aux_rank, ierr)

                this%map(i,1) =  rank_counter
                this%map(i,2) =  aux_rank
                this%map(i,3) =  MOD(global_i, grid_Ndims(1))
                this%map(i,4) =  MOD(global_j, grid_Ndims(2))
                this%map(i,5) =  MOD(global_k, grid_Ndims(3))

                if (this%map(i,3) == 0) this%map(i,3) = grid_Ndims(1)
                if (this%map(i,4) == 0) this%map(i,4) = grid_Ndims(2)
                if (this%map(i,5) == 0) this%map(i,5) = grid_Ndims(3)

                if (rank_counter < n_left) then
                    if (n >= n_base_loc + 1) then
                        rank_counter = rank_counter + 1
                        n = 0
                    end if
                else
                    if (n >= n_base_loc) then
                        rank_counter = rank_counter + 1
                        n = 0
                    end if
                end if

                n = n + 1
            end do

        case default
            write(*,*) "Error: Unknown dimensions in do_mapping."
            stop
        end select

#else

        select case(dimensions)
        case (1)
            rank_counter = 0

            do i = 1, this%n_systems
                this%map(i,1) = 0
                this%map(i,2) = 0
                this%map(i,3) = int(grid_coord(i, 1)/dr) + int(grid_Ndims(1)/2)
            end do

        case (2)

            do i = 1, this%n_systems
                this%map(i,1) = 0
                this%map(i,2) = 0
                this%map(i,3) = int(grid_coord(i, 1)/dr) + int(grid_Ndims(1)/2)
                this%map(i,4) = int(grid_coord(i, 2)/dr) + int(grid_Ndims(2)/2)
            end do

        case (3)

            do i = 1, this%n_systems
                this%map(i,1) = 0
                this%map(i,2) = 0
                this%map(i,3) = int(grid_coord(i, 1)/dr) + int(grid_Ndims(1)/2)
                this%map(i,4) = int(grid_coord(i, 2)/dr) + int(grid_Ndims(2)/2)
                this%map(i,5) = int(grid_coord(i, 3)/dr) + int(grid_Ndims(3)/2)
            end do

        case default
            write(*,*) "Error: Unknown dimensions in do_mapping."
            stop
        end select

#endif

    end subroutine do_mapping
    
!###################################################################################################
 
    subroutine do_kernel_mapping(this, grid_coord, dimensions, mpi_dims, mpi_coords, mpi_cart_comm, &
                                 grid_Ndims, dr)
        class(TQ_group), intent(inout) :: this
        real(dp)       , intent(in)    :: dr
        real(dp)       , intent(in)    :: grid_coord(:,:)
        integer        , intent(in)    :: dimensions
        integer        , intent(in)    :: mpi_dims(3)
        integer        , intent(in)    :: mpi_coords(3)
        integer        , intent(in)    :: mpi_cart_comm
        integer        , intent(in)    :: grid_Ndims(3)

        integer :: rank_counter
        integer :: n_total_ranks
        integer :: n_base_loc
        integer :: n_left
        integer :: i, n
        integer :: ii, jj, kk
        integer :: global_i, global_j, global_k
        integer :: ierr
        integer :: rank_coor(3)
        integer :: aux_rank

        n_total_ranks = mpi_dims(1)*mpi_dims(2)*mpi_dims(3)
        n_base_loc    = int(this%n_systems / n_total_ranks)
        n_left        = mod(this%n_systems, n_total_ranks)

#ifdef USE_MPI

       select case(dimensions)
        case (1)
            rank_counter = 0
            n = 1 
            do i = 1, this%n_systems
                
                do ii = -this%n_ker, this%n_ker
                    this%kernel_map(i,ii,0,0,1) = rank_counter
                    this%kernel_map(i,ii,0,0,1) = 0
                    this%kernel_map(i,ii,0,0,3) = ii + global_i
                end do

                if (rank_counter < n_left) then
                    if (n >= n_base_loc + 1) then
                        rank_counter = rank_counter + 1
                        n = 0
                    end if
                else
                    if (n >= n_base_loc) then
                        rank_counter = rank_counter + 1
                        n = 0
                    end if
                end if

                n = n + 1

            end do

        case (2)

            rank_counter = 0
            n = 1 

            do i = 1, this%n_systems
                global_i = int(grid_coord(i, 1)/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
                global_j = int(grid_coord(i, 2)/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)


                do jj = -this%n_ker, this%n_ker
                do ii = -this%n_ker, this%n_ker

                    rank_coor(1) = int((ii+global_i-1)/grid_Ndims(1))
                    rank_coor(2) = int((jj+global_j-1)/grid_Ndims(2))

                    call MPI_Cart_rank(mpi_cart_comm, rank_coor(:2), aux_rank, ierr)

                    this%kernel_map(i,ii,jj,0,1) =  rank_counter
                    this%kernel_map(i,ii,jj,0,2) =  aux_rank
                    this%kernel_map(i,ii,jj,0,3) =  MOD(ii+global_i, grid_Ndims(1))
                    this%kernel_map(i,ii,jj,0,4) =  MOD(jj+global_j, grid_Ndims(2))
                    if (this%kernel_map(i,ii,jj,0,3) == 0) this%kernel_map(i,ii,jj,0,3) = grid_Ndims(1)
                    if (this%kernel_map(i,ii,jj,0,4) == 0) this%kernel_map(i,ii,jj,0,4) = grid_Ndims(2)
                    
                end do
                end do

                if (rank_counter < n_left) then
                    if (n >= n_base_loc + 1) then
                        rank_counter = rank_counter + 1
                        n = 0
                    end if
                else
                    if (n >= n_base_loc) then
                        rank_counter = rank_counter + 1
                        n = 0
                    end if
                end if

                n = n + 1

            end do

        case (3)

            rank_counter = 0
            n = 1
            do i = 1, this%n_systems
                global_i = int(grid_coord(i, 1)/dr) + int(grid_Ndims(1)*mpi_dims(1)/2)
                global_j = int(grid_coord(i, 2)/dr) + int(grid_Ndims(2)*mpi_dims(2)/2)
                global_k = int(grid_coord(i, 3)/dr) + int(grid_Ndims(3)*mpi_dims(3)/2)

                do kk = -this%n_ker, this%n_ker
                do jj = -this%n_ker, this%n_ker
                do ii = -this%n_ker, this%n_ker

                    rank_coor(1) = int((ii+global_i-1)/grid_Ndims(1))
                    rank_coor(2) = int((jj+global_j-1)/grid_Ndims(2))
                    rank_coor(3) = int((kk+global_k-1)/grid_Ndims(3))

                    call MPI_Cart_rank(mpi_cart_comm, rank_coor, aux_rank, ierr)

                    this%kernel_map(i,ii,jj,kk,1) =  rank_counter
                    this%kernel_map(i,ii,jj,kk,2) =  aux_rank
                    this%kernel_map(i,ii,jj,kk,3) =  MOD(ii+global_i, grid_Ndims(1))
                    this%kernel_map(i,ii,jj,kk,4) =  MOD(jj+global_j, grid_Ndims(2))
                    this%kernel_map(i,ii,jj,kk,5) =  MOD(kk+global_k, grid_Ndims(3))
                    if (this%kernel_map(i,ii,jj,kk,3) == 0) this%kernel_map(i,ii,jj,kk,3) = grid_Ndims(1)
                    if (this%kernel_map(i,ii,jj,kk,4) == 0) this%kernel_map(i,ii,jj,kk,4) = grid_Ndims(2)
                    if (this%kernel_map(i,ii,jj,kk,5) == 0) this%kernel_map(i,ii,jj,kk,5) = grid_Ndims(3)

                end do
                end do
                end do

                if (rank_counter < n_left) then
                    if (n >= n_base_loc + 1) then
                        rank_counter = rank_counter + 1
                        n = 0
                    end if
                else
                    if (n >= n_base_loc) then
                        rank_counter = rank_counter + 1
                        n = 0
                    end if
                end if

                n = n + 1

            end do

        case default
            write(*,*) "Error: Unknown dimensions in do_mapping."
            stop
        end select

#else

        select case(dimensions)
        case (1)
            rank_counter = 0

            do i = 1, this%n_systems
                
                global_i = int(grid_coord(i, 1)/dr) + int(grid_Ndims(1)/2)
                
                do ii = -this%n_ker, this%n_ker
                    this%kernel_map(i,ii,0,0,1) = 0
                    this%kernel_map(i,ii,0,0,2) = 0
                    this%kernel_map(i,ii,0,0,3) = ii + global_i
                end do
            end do

        case (2)

            do i = 1, this%n_systems

                global_i = int(grid_coord(i, 1)/dr) + int(grid_Ndims(1)/2)
                global_j = int(grid_coord(i, 2)/dr) + int(grid_Ndims(2)/2)

                do jj = -this%n_ker, this%n_ker
                do ii = -this%n_ker, this%n_ker
                    this%kernel_map(i,ii,jj,0,1) =  0
                    this%kernel_map(i,ii,jj,0,2) =  0
                    this%kernel_map(i,ii,jj,0,3) =  global_i + ii
                    this%kernel_map(i,ii,jj,0,4) =  global_j + jj
                end do
                end do
            end do

        case (3)

            do i = 1, this%n_systems

                global_i = int(grid_coord(i, 1)/dr) + int(grid_Ndims(1)/2)
                global_j = int(grid_coord(i, 2)/dr) + int(grid_Ndims(2)/2)
                global_k = int(grid_coord(i, 3)/dr) + int(grid_Ndims(3)/2)

                do kk = -this%n_ker, this%n_ker
                do jj = -this%n_ker, this%n_ker
                do ii = -this%n_ker, this%n_ker
                    this%kernel_map(i,ii,jj,kk,1) =  0
                    this%kernel_map(i,ii,jj,kk,2) =  0
                    this%kernel_map(i,ii,jj,kk,3) =  global_i + ii
                    this%kernel_map(i,ii,jj,kk,4) =  global_j + jj
                    this%kernel_map(i,ii,jj,kk,5) =  global_k + kk
                end do
                end do          
                end do

            end do

        case default
            write(*,*) "Error: Unknown dimensions in do_mapping."
            stop
        end select

#endif

    end subroutine do_kernel_mapping

!###################################################################################################

    subroutine init_all_q_systems(this, mol_id, id_file, myrank, print_on)
        class(TQ_group), intent(inout) :: this
        integer        , intent(in)    :: mol_id(:)
        integer        , intent(in)    :: id_file(:)
        integer        , intent(in)    :: myrank
        logical        , intent(in)    :: print_on(:)

        integer  :: i
        integer  :: n_mol
        integer  :: Nt_steps
        real(dp) :: dt

        dt = this%dt
        Nt_steps = this%Nt_steps

        if (.not. allocated(this%q_sys)) then
            select case (this%q_sys_type)
            case (Q_SYS_DFTB)
                this%q_sys = q_system_factory(this%q_sys_type, this%n_sys_loc)
            case default
                write(*,*) "Error: Unknown Q_sys type in init_all_q_systems."
                stop
            end select
        end if

        n_mol = 0
        if (this%group_type == Q_MATERIAL) then
            do i = 1, this%n_systems
                if (myrank == this%map(i,1)) then
                    n_mol = n_mol + 1
                    call this%q_sys(n_mol)%init(mol_id(i), id_file(i), dt, Nt_steps, myrank, print_on(i))
                end if
            end do
        else if (this%group_type == Q_SINGLE) then
            do i = 1, this%n_systems
                if (myrank == this%kernel_map(i,0,0,0,1)) then
                    n_mol = n_mol + 1
                    call this%q_sys(n_mol)%init(mol_id(i), id_file(i), dt, Nt_steps, myrank, print_on(i))
                end if
            end do
        else
            write(*,*) "Error: Unknown Q_group type in init_all_q_systems."
            stop
        end if

        if (n_mol /= this%n_sys_loc) then
            write(*,*) "Error: n_sys_loc mismatch in init_all_q_systems.", n_mol, this%n_sys_loc
            stop
        end if

    !Calculate initial ground state
    !TO-DO: this should be somewhere else.
        do i = 1, this%n_sys_loc
            call this%q_sys(i)%gs_calculate()
        end do

    end subroutine init_all_q_systems

!###################################################################################################
 
    subroutine td_propagate_q_group(this, tq_step, move_q_system)

        class(TQ_group), intent(inout) :: this
        integer        , intent(in)    :: tq_step
        logical        , intent(in)    :: move_q_system

        integer  :: i
        real(dp) :: E_field(3)

        if (.not. move_q_system) return

        do i = 1, this%n_sys_loc
            E_field = this%E_field_list(i, :)
            call this%q_sys(i)%td_propagate(tq_step, E_field)
        end do

    end subroutine td_propagate_q_group

!###################################################################################################
    
end module q_group_mod