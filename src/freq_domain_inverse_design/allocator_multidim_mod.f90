module allocator_multidim_mod
 
     use constants_mod
     use rs_vec_base_mod
     use rs_vec_dimensions_mod

    interface rs_vec_factory
        module procedure rs_vec_factory_single
        module procedure rs_vec_factory_array
    end interface

    interface allocate_multidim
        module procedure allocate_multidim_C
        module procedure allocate_multidim_R
        module procedure allocate_multidim_L
    end interface

contains

!###################################################################################################

function rs_vec_factory_single(dim) result(vec)
    integer, intent(in) :: dim
    class(TRSvec), allocatable :: vec

    select case (dim)
    case (1)
        allocate(TRSvec_1D :: vec)
    case (2)
        allocate(TRSvec_2D :: vec)        
    case (3)
        allocate(TRSvec_3D :: vec)
    case default
        write (*, '("Error: invalid number of dimensions: ", I0)') dim
        error stop
    end select
end function rs_vec_factory_single

!###################################################################################################

function rs_vec_factory_array(dim, n_max,n_min) result(vec)
    integer, intent(in) :: dim
    integer, intent(in) :: n_max
    integer, optional, intent(in) :: n_min
    class(TRSvec), allocatable :: vec(:)

    integer :: n

    if (present(n_min)) then
        n = n_min
    else
        n = 1
    end if

    select case (dim)
    case (1)
        allocate(TRSvec_1D :: vec(n:n_max))
    case (2)
        allocate(TRSvec_2D :: vec(n:n_max))        
    case (3)
        allocate(TRSvec_3D :: vec(n:n_max))
    case default
        write (*, '("Error: invalid number of dimensions: ", I0)') dim
        error stop
    end select
end function rs_vec_factory_array

!###################################################################################################

subroutine allocate_multidim_C(array, dim, i_max, i_min, j_max, j_min, k_max, k_min)
    
    complex(dp), allocatable :: array(:,:,:)
    integer, intent(in) :: dim
    integer, intent(in) :: i_max
    integer, intent(in) :: j_max
    integer, intent(in) :: k_max
    integer, optional, intent(in) :: i_min
    integer, optional, intent(in) :: j_min
    integer, optional, intent(in) :: k_min

    integer :: i_min_loc, j_min_loc, k_min_loc

    if (present(i_min)) then
        i_min_loc = i_min
    else
        i_min_loc = 1
    end if

    if (present(j_min)) then
        j_min_loc = j_min
    else
        j_min_loc = 1
    end if

    if (present(k_min)) then
        k_min_loc = k_min
    else
        k_min_loc = 1
    end if

    select case (dim)
    case (1)
        if (.not. allocated(array)) allocate(array(i_min_loc:i_max,1,1))
    case (2)
        if (.not. allocated(array)) allocate(array(i_min_loc:i_max, j_min_loc:j_max,1))
    case (3)
        if (.not. allocated(array)) allocate(array(i_min_loc:i_max, j_min_loc:j_max, &
                                             k_min_loc:k_max))
    end select

end subroutine allocate_multidim_C

!###################################################################################################

subroutine allocate_multidim_R(array, dim, i_max, i_min, j_max, j_min, k_max, k_min)
    
    real(dp), allocatable :: array(:,:,:)
    integer, intent(in) :: dim
    integer, intent(in) :: i_max
    integer, intent(in) :: j_max
    integer, intent(in) :: k_max
    integer, optional, intent(in) :: i_min
    integer, optional, intent(in) :: j_min
    integer, optional, intent(in) :: k_min

    integer :: i_min_loc, j_min_loc, k_min_loc

    if (present(i_min)) then
        i_min_loc = i_min
    else
        i_min_loc = 1
    end if

    if (present(j_min)) then
        j_min_loc = j_min
    else
        j_min_loc = 1
    end if

    if (present(k_min)) then
        k_min_loc = k_min
    else
        k_min_loc = 1
    end if

    select case (dim)
    case (1)
        if (.not. allocated(array)) allocate(array(i_min_loc:i_max,1,1))
    case (2)
        if (.not. allocated(array)) allocate(array(i_min_loc:i_max, j_min_loc:j_max,1))
    case (3)
        if (.not. allocated(array)) allocate(array(i_min_loc:i_max, j_min_loc:j_max, &
                                             k_min_loc:k_max))
    end select

end subroutine allocate_multidim_R

!###################################################################################################

subroutine allocate_multidim_L(array, dim, i_max, i_min, j_max, j_min, k_max, k_min)
    
    logical, allocatable :: array(:,:,:)
    integer, intent(in) :: dim
    integer, intent(in) :: i_max
    integer, intent(in) :: j_max
    integer, intent(in) :: k_max
    integer, optional, intent(in) :: i_min
    integer, optional, intent(in) :: j_min
    integer, optional, intent(in) :: k_min

    integer :: i_min_loc, j_min_loc, k_min_loc

    if (present(i_min)) then
        i_min_loc = i_min
    else
        i_min_loc = 1
    end if

    if (present(j_min)) then
        j_min_loc = j_min
    else
        j_min_loc = 1
    end if

    if (present(k_min)) then
        k_min_loc = k_min
    else
        k_min_loc = 1
    end if

    select case (dim)
    case (1)
        if (.not. allocated(array)) allocate(array(i_min_loc:i_max,1,1))
    case (2)
        if (.not. allocated(array)) allocate(array(i_min_loc:i_max, j_min_loc:j_max,1))
    case (3)
        if (.not. allocated(array)) allocate(array(i_min_loc:i_max, j_min_loc:j_max, &
                                             k_min_loc:k_max))
    end select

end subroutine allocate_multidim_L

!###################################################################################################

end module allocator_multidim_mod