module q_sys_base_mod

    use constants_mod

    implicit none

    type, abstract :: TQ_sys_base

        integer  :: id
        integer  :: id_file
        integer  :: t_steps
        integer  :: rank
        logical  :: print_on
        real(dp) :: dt
        real(dp) :: E0
        real(dp) :: Et

        real(dp), allocatable :: dipole(:)
        real(dp), allocatable :: dip_old(:)
        real(dp), allocatable :: dPt_dt(:)
        real(dp), allocatable :: dPt_dt_old(:)

        contains
            procedure(init_interface)                , deferred :: init
            procedure(kill_interface)                , deferred :: kill
            procedure(gs_calculate_interface)        , deferred :: gs_calculate
            procedure(td_propagate_interface)        , deferred :: td_propagate
            procedure(write_output_interface)        , deferred :: write_output

    end type TQ_sys_base

    abstract interface

        subroutine init_interface(this, id, id_file, dt, t_steps, rank, print_on)
            ! This subroutine alse read the needed variables from the molecule files.
            import :: TQ_sys_base, dp
            class(TQ_sys_base), intent(inout) :: this
            integer           , intent(in)    :: id
            integer           , intent(in)    :: id_file
            real(dp)          , intent(in)    :: dt
            integer           , intent(in)    :: t_steps
            integer           , intent(in)    :: rank
            logical           , intent(in)    :: print_on
        end subroutine init_interface

        subroutine kill_interface(this)
            import :: TQ_sys_base
            class(TQ_sys_base), intent(inout) :: this
        end subroutine kill_interface

        subroutine gs_calculate_interface(this)
            import :: TQ_sys_base
            class(TQ_sys_base), intent(inout) :: this
        end subroutine gs_calculate_interface

        subroutine td_propagate_interface(this, tq_step, E_field)
            import :: TQ_sys_base, dp
            class(TQ_sys_base), intent(inout) :: this
            real(dp)     , intent(in)    :: E_field(3)
            integer      , intent(in)    :: tq_step
        end subroutine td_propagate_interface

        subroutine write_output_interface(this, time, i_group, print_q_step)
             import :: TQ_sys_base, dp
             class(TQ_sys_base), intent(in) :: this
             real(dp)          , intent(in) :: time
             integer           , intent(in) :: i_group
             integer           , intent(in) :: print_q_step
        end subroutine write_output_interface

    end interface

end module q_sys_base_mod