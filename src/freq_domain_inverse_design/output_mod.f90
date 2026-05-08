module output_mod

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
    use medium_mod

    implicit none

    type :: TOutput

        character(len=20) :: output_file          ='output_gral.out'
        logical           :: print_restart        = .false.
        integer           :: unit_general_output  = 100
        integer           :: unit_final_structure = 101
        integer           :: unit_final_fields    = 102
        integer           :: unit_restart         = 200
        integer           :: n_fields             = 0

        integer, allocatable :: field_list(:)

        contains
    
            procedure :: init_output
            procedure :: write_gral_output
            procedure :: write_restart
            procedure :: write_final_structure
            procedure :: write_final_fields
            procedure :: close_output

    end type TOutput

contains

!###################################################################################################

subroutine init_output(this, print_restart, print_field_list, myrank)

    class(TOutput), intent(inout) :: this
    logical, intent(in) :: print_restart
    integer, intent(in) :: print_field_list(6)
    integer, intent(in) :: myrank

    integer :: n_fields
    integer :: ierr
    integer :: i

    this%print_restart = print_restart

    n_fields = 0
    do i = 1, 6
        if (print_field_list(i) /= PRINT_NONE) then
            n_fields = n_fields + 1
        end if
    end do

    this%n_fields = n_fields

    if (n_fields > 0) then
        if (.not. allocated(this%field_list)) allocate(this%field_list(n_fields))
        n_fields = 0
        do i = 1, 6
            if (print_field_list(i) /= PRINT_NONE) then
                n_fields = n_fields + 1
                this%field_list(n_fields) = print_field_list(i)
            end if
        end do
    end if

    if (myrank == 0) then

        open(unit=this%unit_general_output, file=this%output_file, status='replace', action='write')

        write(this%unit_general_output, '(a)') '# iter_step  FOM             delta_rho &
                                            &        beta            |grad_max|'

        if (this%print_restart) then
            call execute_command_line("mkdir -p " // trim("restart"), exitstat=ierr)
            if (ierr /= 0) then
                write (*, '("Error: could not create restart directory")')
                error stop
            end if
        end if

        call execute_command_line("mkdir -p " // trim("./outputs"), exitstat=ierr)
        if (ierr /= 0) then
            write (*, '("Error: could not create outputs directory")')
            error stop
        end if

        this%unit_restart = 200 + myrank

    end if

end subroutine init_output

!###################################################################################################

subroutine write_gral_output(this, iter_step, fom, delta_rho, beta, grad_max, myrank)

    class(TOutput), intent(inout) :: this
    integer, intent(in)  :: iter_step
    integer, intent(in)  :: myrank
    real(dp), intent(in) :: fom
    real(dp), intent(in) :: delta_rho
    real(dp), intent(in) :: beta
    real(dp), intent(in) :: grad_max

    if (myrank == 0) then

        write(this%unit_general_output, *) iter_step, fom, delta_rho, beta, grad_max

    end if

end subroutine write_gral_output

!###################################################################################################

subroutine write_restart(this, f_vec, f_adj_list, n_trg, design, p_prblm,&
                         mpi_dims, mpi_coords, beta)

    class(TOutput), intent(inout) :: this
    class(TRSvec) , intent(in)    :: f_vec
    class(TRSvec) , intent(in)    :: f_adj_list(n_trg)
    class(TDesign), intent(in)    :: design
    integer       , intent(in)    :: n_trg
    integer       , intent(in)    :: p_prblm
    integer       , intent(in)    :: mpi_dims(3)
    integer       , intent(in)    :: mpi_coords(3)
    real(dp)      , intent(in)    :: beta

    character(len=30)  :: restart_file="./restart/restart"
    character(len=5)   :: extension=".bin"
    character(len=10)  :: beta_number
    character(len=100) :: full_restart_file
    integer            :: ierr
    integer            :: beta_step
    integer            :: t

#ifdef USE_MPI
    character(len=10)  :: mpi_name = "mpi_coords"
    character(len=10)  :: x_ranks, y_ranks, z_ranks
#endif
    
    if (.not. this%print_restart) return

    beta_step = int(beta)
    write(beta_number, '(I3.3)') beta_step

#ifdef USE_MPI

    write (x_ranks, '(A,I3.3)') "_x", mpi_dims(1)
    write (y_ranks, '(A,I3.3)') "_y", mpi_dims(2)
    write (z_ranks, '(A,I3.3)') "_z", mpi_dims(3)

    mpi_name = mpi_name // "_" // trim(x_ranks) // "_" //trim(y_ranks) // "_" // trim(z_ranks)
    full_restart_file = trim(restart_file) // "_" // trim(mpi_name) // "_" // &
                       "_beta_" //trim(beta_number) // trim(extension)
#else
    full_restart_file = trim(restart_file) // "_beta_" // trim(beta_number) // trim(extension)
#endif

    if (p_prblm == 1) then

        open(unit=this%unit_restart, file=trim(full_restart_file), status='replace', &
            access='stream', form='unformatted', &
            action='write', iostat=ierr)
        if (ierr /= 0) then
            write(*,*) 'Error opening restart file for write:'
            write(*,*) trim(full_restart_file)
            error stop
        end if

        select type(design)
        class is (TDesign_1D)
            write(this%unit_restart) design%rho(:)
        class is (TDesign_2D)
            write(this%unit_restart) design%rho(:,:)
        class is (TDesign_3D)
            write(this%unit_restart) design%rho(:,:,:)
        end select

        close(this%unit_restart)
            
    end if

    open(unit=this%unit_restart, file=trim(full_restart_file), status='old', &
        action='write', position='append', access='stream', form='unformatted',     &
        iostat=ierr)
    if (ierr /= 0) then
        write(*,*) 'Error opening restart file for append:'
        write(*,*) trim(full_restart_file)
        error stop
    end if
    
    select type(f_vec)

        class is (TRSvec_1D)
        select type(f_adj_list)
        class is (TRSvec_1D)

            write(this%unit_restart) f_vec%pl_x
            write(this%unit_restart) f_vec%pl_y
            write(this%unit_restart) f_vec%mi_x
            write(this%unit_restart) f_vec%mi_y

            do t = 1, n_trg
                write(this%unit_restart) f_adj_list(t)%pl_x
                write(this%unit_restart) f_adj_list(t)%pl_y
                write(this%unit_restart) f_adj_list(t)%mi_x
                write(this%unit_restart) f_adj_list(t)%mi_y
            end do

        end select
        class is (TRSvec_2D)
        select type(f_adj_list)
        class is (TRSvec_2D)

            write(this%unit_restart) f_vec%pl_x
            write(this%unit_restart) f_vec%pl_y
            write(this%unit_restart) f_vec%pl_z
            write(this%unit_restart) f_vec%mi_x
            write(this%unit_restart) f_vec%mi_y
            write(this%unit_restart) f_vec%mi_z

            do t = 1, n_trg
                write(this%unit_restart) f_adj_list(t)%pl_x
                write(this%unit_restart) f_adj_list(t)%pl_y
                write(this%unit_restart) f_adj_list(t)%pl_z
                write(this%unit_restart) f_adj_list(t)%mi_x
                write(this%unit_restart) f_adj_list(t)%mi_y
                write(this%unit_restart) f_adj_list(t)%mi_z
            end do

        end select
        class is (TRSvec_3D)
        select type(f_adj_list)
        class is (TRSvec_3D)

            write(this%unit_restart) f_vec%pl_x
            write(this%unit_restart) f_vec%pl_y
            write(this%unit_restart) f_vec%pl_z
            write(this%unit_restart) f_vec%mi_x
            write(this%unit_restart) f_vec%mi_y
            write(this%unit_restart) f_vec%mi_z

            do t = 1, n_trg
                write(this%unit_restart) f_adj_list(t)%pl_x
                write(this%unit_restart) f_adj_list(t)%pl_y
                write(this%unit_restart) f_adj_list(t)%pl_z
                write(this%unit_restart) f_adj_list(t)%mi_x
                write(this%unit_restart) f_adj_list(t)%mi_y
                write(this%unit_restart) f_adj_list(t)%mi_z
            end do

        end select

    end select

    close(this%unit_restart)

end subroutine write_restart
!###################################################################################################

subroutine close_output(this)

    class(TOutput), intent(inout) :: this
    logical :: is_open

    inquire(unit=this%unit_general_output, opened=is_open)
    if (is_open) close(this%unit_general_output)

    inquire(unit=this%unit_restart, opened=is_open)
    if (is_open) close(this%unit_restart)

end subroutine close_output

!###################################################################################################

subroutine write_final_structure(this, eps_r, grid_Ndims, dr, p_prblm, mpi_coords, mpi_dims)

    class(TOutput)      , intent(inout) :: this
    class(TMedium_eps_r), intent(in)    :: eps_r
    integer             , intent(in)    :: grid_Ndims(3)
    integer             , intent(in)    :: mpi_coords(3)
    integer             , intent(in)    :: mpi_dims(3)
    integer             , intent(in)    :: p_prblm
    real(dp)            , intent(in)    :: dr


    character(len=60)  :: structure_file="./outputs/final_permitivity_prblm_"
    character(len=10)  :: prblm_number
    character(len=5)   :: extension=".out"
    integer            :: ierr
    integer            :: i, j, k
    integer            :: i_rank, j_rank, k_rank
    integer            :: nx, ny, nz
    integer            :: nx_rank, ny_rank, nz_rank
    logical            :: is_root
    real(dp)           :: x,y,z
    real(dp)           :: Re_eps, Im_eps

    nx_rank = mpi_dims(1)
    ny_rank = mpi_dims(2)
    nz_rank = mpi_dims(3)

    nx = grid_Ndims(1)
    ny = grid_Ndims(2)
    nz = grid_Ndims(3)

    write(prblm_number, '(I3.3)') p_prblm

    structure_file = trim(structure_file) // trim(prblm_number) // trim(extension)

    is_root = all(mpi_coords == 0)

    if (is_root) then
        open(unit=this%unit_final_structure, file=trim(structure_file), status='replace', &
             action='write', iostat=ierr)
        if (ierr /= 0) then
            write (*, '("Error: could not open final structure output file")')
            error stop
        end if

        select case (eps_r%dimensions)
        case (1)
            write(this%unit_final_structure, '(a)') '#    z [Arb.u.]          Re(eps_r) &
                                                    &           Im(eps_r)'
        case (2)
            write(this%unit_final_structure, '(a)') '#    x [Arb.u.]          y [Arb.u.] &
                                                    &         Re(eps_r)           Im(eps_r)'
        case (3)
            write(this%unit_final_structure, '(a)') '#    x [Arb.u.]          y [Arb.u.] &
                                                    &          z [Arb.u.] &
                                                    &         Re(eps_r)           Im(eps_r)'
        end select

        close(this%unit_final_structure)

    end if

#ifdef USE_MPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

    do k_rank = 0, nz_rank-1
    do j_rank = 0, ny_rank-1
    do i_rank = 0, nx_rank-1
        if (mpi_coords(1) == i_rank .and. mpi_coords(2) == j_rank .and. & 
            mpi_coords(3) == k_rank) then
        
            open(unit=this%unit_final_structure, file=trim(structure_file), status='old', &
                 action='write', position='append', iostat=ierr)
            if (ierr /= 0) then
                write (*, '("Error: could not open final structure output file for appending")')
                error stop
            end if

            select case (eps_r%dimensions)
            case (1)
                do i = 1, nx
                    Re_eps = REAL(eps_r%mat1D(i))
                    Im_eps = IMAG(eps_r%mat1D(i))
                    x = (i - int(grid_Ndims(1)*mpi_dims(1)/2)) * dr
                    write(this%unit_final_structure, *) x,  Re_eps, Im_eps
                end do
            case (2)
                do j = 1, ny
                    y = (j + ny*mpi_coords(2) - int(ny*mpi_dims(2)/2)) * dr
                    do i = 1, nx
                        x = (i + nx*mpi_coords(1) - int(nx*mpi_dims(1)/2)) * dr
                        Re_eps = REAL(eps_r%mat2D(i,j))
                        Im_eps = IMAG(eps_r%mat2D(i,j))
                        write(this%unit_final_structure, *) x, y, Re_eps, Im_eps
                    end do
                end do
            case (3)
                do k = 1, nz
                    z = (k + nz*mpi_coords(3) - int(nz*mpi_dims(3)/2)) * dr
                    do j = 1, ny
                        y = (j + ny*mpi_coords(2) - int(ny*mpi_dims(2)/2)) * dr
                        do i = 1, nx
                            x = (i + nx*mpi_coords(1) - int(nx*mpi_dims(1)/2)) * dr
                            Re_eps = REAL(eps_r%mat3D(i,j,k))
                            Im_eps = IMAG(eps_r%mat3D(i,j,k))
                            write(this%unit_final_structure, *) x, y, z, Re_eps, Im_eps
                        end do
                    end do
                end do
            end select

            close(this%unit_final_structure)
        
        end if

#ifdef USE_MPI
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

    end do
    end do
    end do

end subroutine write_final_structure

!###################################################################################################

subroutine write_final_fields(this, fvec, grid_Ndims, p_prblm, dr, mpi_coords, mpi_dims)

    class(TOutput)      , intent(inout) :: this
    class(TRSvec)       , intent(in)    :: fvec
    integer             , intent(in)    :: grid_Ndims(3)
    integer             , intent(in)    :: p_prblm
    integer             , intent(in)    :: mpi_coords(3)
    integer             , intent(in)    :: mpi_dims(3)
    real(dp)            , intent(in)    :: dr

    character(len=100) :: field_file
    character(len=20)  :: final_file="./outputs/final"
    character(len=5)   :: extension=".out"
    character(len=10)  :: prblm_name = "prblm"
    character(len=10)  :: prblm_number
    character(len=10)  :: field_name
    integer            :: ierr

    complex(dp)        :: field
    integer            :: i, j, k, f
    integer            :: i_rank, j_rank, k_rank
    integer            :: nx, ny, nz
    integer            :: nx_rank, ny_rank, nz_rank
    logical            :: is_root
    real(dp)           :: x,y,z
    real(dp)           :: Re_field, Im_field
    real(dp)           :: C1, C2

    write(prblm_number, '(I3.3)') p_prblm
    
    C1 = 1.0_dp / DSQRT(2.0_dp * eps0)
    C2 = 1.0_dp / DSQRT(1 / (2.0_dp * mu0))

    nx_rank = mpi_dims(1)
    ny_rank = mpi_dims(2)
    nz_rank = mpi_dims(3)
    
    nx = grid_Ndims(1)
    ny = grid_Ndims(2)
    nz = grid_Ndims(3)
    
    is_root = all(mpi_coords == 0)
    
    if (is_root) then
        do f =1, this%n_fields
            
            select case (this%field_list(f))
            case (PRINT_Ex)
                field_name = "Ex"
            case (PRINT_Ey)
                field_name = "Ey"
            case (PRINT_Ez)
                field_name = "Ez"
            case (PRINT_Hx)
                field_name = "Hx"
            case (PRINT_Hy)
                field_name = "Hy"
            case (PRINT_Hz)
                field_name = "Hz"
            end select
            
            field_file = trim(final_file) // "_" // trim(field_name)  // "_" // &
                         trim(prblm_name) // trim(prblm_number) // trim(extension)

            open(unit=this%unit_final_fields, file=trim(field_file), status='replace', &
                action='write', iostat=ierr)
            if (ierr /= 0) then
                write (*, '("Error: could not open final E-field output file")')
                error stop
            end if

            select case (fvec%dimensions)
            case (1)
                write(this%unit_final_fields, '(a)') '#    z [Arb.u.] &
                                                            &           Re('//trim(field_name)//') &
                                                            &           Im('//trim(field_name)//')'
            case (2)
                write(this%unit_final_fields, '(a)') '#    x [Arb.u.]          y [Arb.u.] &
                                                            &         Re('//trim(field_name)//') &
                                                            &         Im('//trim(field_name)//')'
            case (3)
                write(this%unit_final_fields, '(a)') '#    x [Arb.u.]          y [Arb.u.] &
                                                            &          z [Arb.u.] &
                                                            &         Re('//trim(field_name)//') &
                                                            &         Im('//trim(field_name)//')'
            end select

            close(this%unit_final_fields)
            
        end do
    end if

#ifdef USE_MPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

    do f =1, this%n_fields

        select case (this%field_list(f))
        case (PRINT_Ex)
            field_name = "Ex"
        case (PRINT_Ey)
            field_name = "Ey"
        case (PRINT_Ez)
            field_name = "Ez"
        case (PRINT_Hx)
            field_name = "Hx"
        case (PRINT_Hy)
            field_name = "Hy"
        case (PRINT_Hz)
            field_name = "Hz"
        end select

            
        field_file = trim(final_file) // "_" // trim(field_name)  // "_" // &
                     trim(prblm_name) // trim(prblm_number) // trim(extension)

        do k_rank = 0, nz_rank-1
        do j_rank = 0, ny_rank-1
        do i_rank = 0, nx_rank-1

            if (mpi_coords(1) == i_rank .and. mpi_coords(2) == j_rank .and. & 
                mpi_coords(3) == k_rank) then
           
                open(unit=this%unit_final_fields, file=trim(field_file), status='old', &
                    action='write', position='append', iostat=ierr)
                if (ierr /= 0) then
                    write (*, '("Error: could not open final field output file for appending")')
                    error stop
                end if

                select type (fvec)
                class is (TRSvec_1D)
                    do i = 1, nx
                        x = (i - int(grid_Ndims(1)*mpi_dims(1)/2)) * dr

                        select case (this%field_list(f))
                        case (PRINT_Ex)
                                field = C1 *(fvec%pl_x(i) + fvec%mi_x(i))
                        case (PRINT_Hy)
                                field = -Z_I*C2 *(fvec%pl_y(i) - fvec%mi_y(i))
                        end select

                        Re_field = REAL(field)
                        Im_field = IMAG(field)
                        write(this%unit_final_fields, *) x, Re_field, Im_field

                    end do

                class is (TRSvec_2D)
                    do j = 1, ny
                        y = (j + ny * mpi_coords(2) - int(ny*mpi_dims(2)/2)) * dr
                        do i = 1, nx
                            x = (i + nx * mpi_coords(1) - int(nx*mpi_dims(1)/2)) * dr

                            select case (this%field_list(f))
                            case (PRINT_Ex)
                                    field = C1 *(fvec%pl_x(i,j) + fvec%mi_x(i,j))
                            case (PRINT_Ey)
                                    field = C1 *(fvec%pl_y(i,j) + fvec%mi_y(i,j))
                            case (PRINT_Ez)
                                    field = C1 *(fvec%pl_z(i,j) + fvec%mi_z(i,j))
                            case (PRINT_Hx)
                                    field = -Z_I*C2 *(fvec%pl_x(i,j) - fvec%mi_x(i,j))
                            case (PRINT_Hy)
                                    field = -Z_I*C2 *(fvec%pl_y(i,j) - fvec%mi_y(i,j))
                            case (PRINT_Hz)
                                    field = -Z_I*C2 *(fvec%pl_z(i,j) - fvec%mi_z(i,j))
                            end select

                            Re_field = REAL(field)
                            Im_field = IMAG(field)
                            write(this%unit_final_fields, *) x, y, Re_field, Im_field

                        end do
                    end do

                class is (TRSvec_3D)
                    do k = 1, nz
                        z = (k + nz * mpi_coords(3) - int(nz*mpi_dims(3)/2)) * dr
                        do j = 1, ny
                            y = (j + ny * mpi_coords(2) - int(ny*mpi_dims(2)/2)) * dr
                            do i = 1, nx
                                x = (i + nx * mpi_coords(1) - int(nx*mpi_dims(1)/2)) * dr

                                select case (this%field_list(f))
                                case (PRINT_Ex)
                                        field = C1 *(fvec%pl_x(i,j,k) + fvec%mi_x(i,j,k))
                                case (PRINT_Ey)
                                        field = C1 *(fvec%pl_y(i,j,k) + fvec%mi_y(i,j,k))
                                case (PRINT_Ez)
                                        field = C1 *(fvec%pl_z(i,j,k) + fvec%mi_z(i,j,k))
                                case (PRINT_Hx)
                                        field = -Z_I*C2 *(fvec%pl_x(i,j,k) - fvec%mi_x(i,j,k))
                                case (PRINT_Hy)
                                        field = -Z_I*C2 *(fvec%pl_y(i,j,k) - fvec%mi_y(i,j,k))
                                case (PRINT_Hz)
                                        field = -Z_I*C2 *(fvec%pl_z(i,j,k) - fvec%mi_z(i,j,k))
                                end select

                                Re_field = REAL(field)
                                Im_field = IMAG(field)

                                write(this%unit_final_fields, *) x, y, z, Re_field, Im_field

                            end do
                        end do
                    end do
                end select

            end if

#ifdef USE_MPI
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

        end do
        end do
        end do
    end do

end subroutine write_final_fields

!###################################################################################################

end module output_mod