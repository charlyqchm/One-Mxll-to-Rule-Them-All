module bicgstab_mod

    use constants_mod
    use rs_vec_base_mod
    use rs_vec_dimensions_mod
    use allocator_multidim_mod

    implicit none

    type(TRSvec), allocatable :: r0_vec
    type(TRSvec), allocatable :: r_vec(:)
    type(TRSvec), allocatable :: u_vec(:)

    complex(dp) , allocatable :: tau(:,:)
    complex(dp) , allocatable :: gam(:)
    complex(dp) , allocatable :: gam_p(:)
    complex(dp) , allocatable :: gam_pp(:)
    complex(dp) , allocatable :: sig(:)

contains
!###################################################################################################

subroutine init_BICGStab_L_variables(grid_Ndims, dr, dimensionsm, L_term)

    real(dp), intent(in) :: dr(3)
    integer , intent(in) :: grid_Ndims(3)
    integer , intent(in) :: dimensionsm
    integer , intent(in) :: L_term

    r0_vec = rs_vec_factory(dimensionsm)
    r_vec  = rs_vec_factory_array(dim= dimensionsm, n_max= L_term, n_min= 0)
    u_vec  = rs_vec_factory_array(dim= dimensionsm, n_max= L_term, n_min= 0)

    if (.not. allocated(tau))    allocate(tau(L_term-1, L_term))
    if (.not. allocated(gam))    allocate(gam(L_term+1))
    if (.not. allocated(gam_p))  allocate(gam_p(L_term+1))
    if (.not. allocated(gam_pp)) allocate(gam_pp(L_term+1))
    if (.not. allocated(sig))    allocate(sig(L_term))

end subroutine init_BICGStab_L_variables

!###################################################################################################

subroutine kill_BICGStab_L_variables()

    if (allocated(r0_vec)) deallocate(r0_vec)
    if (allocated(r_vec))  deallocate(r_vec)
    if (allocated(u_vec))  deallocate(u_vec)

    if (allocated(tau))    deallocate(tau)
    if (allocated(gam))    deallocate(gam)
    if (allocated(gam_p))  deallocate(gam_p)
    if (allocated(gam_pp)) deallocate(gam_pp)
    if (allocated(sig))    deallocate(sig)

end subroutine kill_BICGStab_L_variables

!###################################################################################################



end module