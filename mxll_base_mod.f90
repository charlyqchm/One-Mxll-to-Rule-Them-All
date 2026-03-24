module mxll_base_mod

    use constants_mod

    implicit none

    type, abstract :: TMxll

        integer  :: dimensions
        integer  :: tq_step
        integer  :: tq_step_old
        integer  :: n_skip_steps
        real(dp) :: time
        real(dp) :: t_skip
!Temporal variable for testing. Remove when not needed anymore.
        integer  :: chunk_coor(3)


    contains
        procedure(init_interface)      , deferred :: init
        procedure(kill_interface)      , deferred :: kill
        procedure(td_prop_H_interface) , deferred :: td_propagate_H_field
        procedure(td_prop_E_interface) , deferred :: td_propagate_E_field

    end type TMxll

    abstract interface

        subroutine init_interface(this, grid_Ndims, npml, boundaries, dt, dr, mode, n_media, &
                                  mpi_coords, mpi_dims)
        
            import :: TMxll, dp
            class(TMxll), intent(inout) :: this
            integer     , intent(in)    :: grid_Ndims(3)
            integer     , intent(in)    :: npml
            integer     , intent(in)    :: boundaries(3)
            integer     , intent(in)    :: mode
            integer     , intent(in)    :: n_media
            integer     , intent(in)    :: mpi_coords(3)
            integer     , intent(in)    :: mpi_dims(3)
            real(dp)    , intent(in)    :: dt
            real(dp)    , intent(in)    :: dr

        end subroutine init_interface

        subroutine kill_interface(this)
            import :: TMxll
            class(TMxll), intent(inout) :: this
        end subroutine kill_interface

        subroutine td_prop_H_interface(this)
            import :: TMxll
            class(TMxll), intent(inout) :: this
        end subroutine td_prop_H_interface

        subroutine td_prop_E_interface(this, t)
            import :: TMxll
            class(TMxll), intent(inout) :: this
            integer     , intent(in)    :: t
        end subroutine td_prop_E_interface

    end interface

end module mxll_base_mod