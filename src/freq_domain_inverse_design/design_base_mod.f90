module design_base_mod

    use constants_mod

    implicit none

    type, abstract :: TDesign

        integer  :: dimensions
        integer  :: nx
        integer  :: ny
        integer  :: nz
        integer  :: n_ker
        real(dp) :: sigma
        real(dp) :: fom
        real(dp) :: drho
        real(dp) :: grad_max

    contains

        procedure(init_design_interface)             ,deferred :: init_design
        procedure(kill_design_interface)             ,deferred :: kill_design
        procedure(collect_opt_regions_interface)     ,deferred :: collect_opt_regions
        procedure(collect_FOM_interface)             ,deferred :: collect_FOM
        procedure(collect_gradients_interface)       ,deferred :: collect_gradients
        procedure(apply_kernel_on_rho_interface)     ,deferred :: apply_kernel_on_rho
        procedure(apply_kernel_on_grad_interface)    ,deferred :: apply_kernel_on_grad
        procedure(displace_rho_interface)            ,deferred :: displace_rho
        procedure(reset_rho_one_step_back_interface) ,deferred :: reset_rho_one_step_back
        procedure(reset_grad_interface)              ,deferred :: reset_grad
        procedure(update_rho_interface)              ,deferred :: update_rho
        procedure(update_grad_interface)             ,deferred :: update_grad
    end type TDesign

    abstract interface

        subroutine init_design_interface(this, dimensions, dr, sigma, grid_Ndims)

            import :: TDesign, dp
            class(TDesign), intent(inout) :: this
            integer     , intent(in)     :: dimensions
            real(dp)    , intent(in)     :: dr
            real(dp)    , intent(in)     :: sigma
            integer     , intent(in)     :: grid_Ndims(3)

        end subroutine init_design_interface

        subroutine kill_design_interface(this)
            import :: TDesign
            class(TDesign), intent(inout) :: this
        end subroutine kill_design_interface

        subroutine collect_opt_regions_interface(this, opt_region_i, rho_init)
            import :: TDesign, dp
            class(TDesign), intent(inout) :: this
            logical       , intent(in)    :: opt_region_i(:,:,:)
            real(dp)      , intent(in)    :: rho_init
        end subroutine collect_opt_regions_interface

        subroutine collect_FOM_interface(this, w_p, p, n_opt_problems)
            import :: TDesign, dp
            class(TDesign), intent(inout) :: this
            real(dp)      , intent(in)    :: w_p
            integer       , intent(in)    :: p
            integer       , intent(in)    :: n_opt_problems
        end subroutine collect_FOM_interface

        subroutine collect_gradients_interface(this, grad_in, p, n_opt_problems)
            import :: TDesign, dp
            class(TDesign), intent(inout) :: this
            real(dp)      , intent(in)    :: grad_in(:,:,:)
            integer       , intent(in)    :: p
            integer       , intent(in)    :: n_opt_problems
        end subroutine collect_gradients_interface

        subroutine apply_kernel_on_rho_interface(this)
            import :: TDesign
            class(TDesign), intent(inout) :: this
        end subroutine apply_kernel_on_rho_interface

        subroutine apply_kernel_on_grad_interface(this)
            import :: TDesign
            class(TDesign), intent(inout) :: this
        end subroutine apply_kernel_on_grad_interface

        subroutine displace_rho_interface(this)
            import :: TDesign
            class(TDesign), intent(inout) :: this
        end subroutine displace_rho_interface

        subroutine reset_rho_one_step_back_interface(this)
            import :: TDesign
            class(TDesign), intent(inout) :: this
        end subroutine reset_rho_one_step_back_interface

        subroutine reset_grad_interface(this)
            import :: TDesign
            class(TDesign), intent(inout) :: this
        end subroutine reset_grad_interface

        subroutine update_rho_interface(this)
            import :: TDesign
            class(TDesign), intent(inout) :: this
        end subroutine update_rho_interface

        subroutine update_grad_interface(this)
            import :: TDesign
            class(TDesign), intent(inout) :: this
        end subroutine update_grad_interface

    end interface

end module design_base_mod