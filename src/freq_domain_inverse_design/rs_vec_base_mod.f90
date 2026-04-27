module rs_vec_base_mod

    use constants_mod

    implicit none

    type, abstract :: TRSvec

        integer  :: dimensions
        integer  :: n_der
        real(dp) :: freq
        real(dp) :: dr
        
    contains

        procedure(init_interface)      , deferred :: init
        procedure(kill_interface)      , deferred :: kill

    end type TRSvec

    abstract interface

        subroutine init_interface(this, grid_Ndims, dr, dimensions, freq, n_der)
        
            import :: TRSvec, dp
            class(TRSvec), intent(inout) :: this
            integer     , intent(in)     :: grid_Ndims(3)
            integer     , intent(in)     :: dimensions
            real(dp)    , intent(in)     :: dr
            real(dp)    , intent(in)     :: freq
            integer     , intent(in)     :: n_der

        end subroutine init_interface

        subroutine kill_interface(this)
            import :: TRSvec
            class(TRSvec), intent(inout) :: this
        end subroutine kill_interface

    end interface

end module rs_vec_base_mod