module factory_mod

    use mxll_base_mod
    use mxll_1D_mod
    use mxll_2D_mod
    use mxll_3D_mod

    use q_sys_base_mod
    use q_sys_dftb_mod

    implicit none

contains

function maxwell_factory(dim) result(maxwell_obj)
    integer, intent(in) :: dim
    class(TMxll), allocatable :: maxwell_obj

    select case (dim)
    case (1)
        allocate(TMxll_1D :: maxwell_obj)
    case (2)
        allocate(TMxll_2D :: maxwell_obj)
    case (3)
        allocate(TMxll_3D :: maxwell_obj)    
    case default
        error stop "The selected dim, must be 1, 2, or 3."
    end select

end function maxwell_factory

function q_system_factory(sys_type, n_sys) result(q_sys_obj)
    integer           , intent(in)  :: sys_type
    integer           , intent(in)  :: n_sys
    class(TQ_sys_base), allocatable :: q_sys_obj(:)

    select case (sys_type)
    case (Q_SYS_DFTB)
        allocate(TQ_sys_dftb :: q_sys_obj(n_sys))
    case default
        error stop "The selected quantum system type is not implemented."
    end select

end function q_system_factory


end module factory_mod