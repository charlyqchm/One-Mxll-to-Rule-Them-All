module bicgstab_mod

    use constants_mod
    use rs_vec_base_mod
    use rs_vec_dimensions_mod
    use medium_mod
    use operator_mod
    use allocator_multidim_mod
    use rs_operations_subs_mod

    implicit none

    class(TRSvec), allocatable :: r0_vec
    class(TRSvec), allocatable :: r_vec(:)
    class(TRSvec), allocatable :: u_vec(:)

    complex(dp) , allocatable :: tau(:,:)
    complex(dp) , allocatable :: gam(:)
    complex(dp) , allocatable :: gam_p(:)
    complex(dp) , allocatable :: gam_pp(:)
    complex(dp) , allocatable :: sig(:)

contains
!###################################################################################################

subroutine init_BICGStab_L_variables(grid_Ndims, dimensions, n_der, L_term)

    integer , intent(in) :: grid_Ndims(3)
    integer , intent(in) :: dimensions
    integer , intent(in) :: n_der
    integer , intent(in) :: L_term

    integer :: i

    if (allocated(r0_vec)) deallocate(r0_vec)
    if (allocated(r_vec))  deallocate(r_vec)
    if (allocated(u_vec))  deallocate(u_vec)

    allocate(r0_vec, mold=rs_vec_factory(dimensions))
    allocate(r_vec(0:L_term), mold=rs_vec_factory(dimensions))
    allocate(u_vec(0:L_term), mold=rs_vec_factory(dimensions))

    call r0_vec%init(grid_Ndims, 0.0_dp, dimensions, 0.0_dp, n_der)

    do i = 0, L_term
        call r_vec(i)%init(grid_Ndims, 0.0_dp, dimensions, 0.0_dp, n_der)
        call u_vec(i)%init(grid_Ndims, 0.0_dp, dimensions, 0.0_dp, n_der)
    end do

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

subroutine BICGStab_L(A_op, f_vec, j_vec, f_vec_out, Af_vec, eps_r, converged, transpose, &
                      tol, n_max, L_iter, abs_err)

    type(TOperator)    , intent(inout) :: A_op
    type(TMedium_eps_r), intent(inout) :: eps_r
    class(TRSvec)      , intent(inout) :: f_vec 
    class(TRSvec)      , intent(inout) :: j_vec
    class(TRSvec)      , intent(inout) :: f_vec_out
    class(TRSvec)      , intent(inout) :: Af_vec
    logical            , intent(out)   :: converged
    logical            , intent(in)    :: transpose
    real(dp)           , intent(in)    :: tol
    real(dp)           , intent(out)   :: abs_err
    integer            , intent(in)    :: n_max
    integer            , intent(in)    :: L_iter

    complex(dp) :: rho_0
    complex(dp) :: rho_1
    complex(dp) :: alpha
    complex(dp) :: omega
    complex(dp) :: beta
    complex(dp) :: dot_product
    integer     :: n, l, i
    real(dp)    :: norm_t
    real(dp)    :: norm_inv 
    real(dp)    :: err_r

    omega = Z_ONE
    rho_0 = Z_ONE
    rho_1 = Z_ONE
    alpha = Z_ONE
    sig   = Z_0

    call reset_V1_to_zero(r0_vec)

    do i = 0, L_iter
        call reset_V1_to_zero(r_vec(i))
        call reset_V1_to_zero(u_vec(i))
    end do

    call copy_V2_on_V1(f_vec_out, f_vec)

    call A_op%apply_operator(f_vec, Af_vec, eps_r, transpose)

    call linear_op_V1_aV2(j_vec, Af_vec, -Z_ONE, r0_vec)

    call copy_V2_on_V1(r_vec(0), r0_vec)

    call dot_product_V1_V2(r0_vec, r0_vec, dot_product)

    norm_t = DSQRT(DBLE(dot_product))

    norm_inv = 1.0_dp/norm_t
    call self_product_aV1(r0_vec, Z_ONE*norm_inv)

    call dot_product_V1_V2(j_vec, j_vec, dot_product)
    norm_t = DSQRT(DBLE(dot_product))

    err_r = 1.0_dp
    converged = .false.

    n = 0
    !Beginning of the BICGStab iteration---------------------------------------!
    do while((n <= n_max) .and. (.not. converged))

        rho_0 = - omega * rho_0

        !Beginning of the L-term BICGStab iteration----------------------------!
        do l = 0, L_iter-1

            call dot_product_V1_V2(r0_vec, r_vec(l), dot_product)

            rho_1 = dot_product
            beta  = alpha * rho_1 / rho_0
            rho_0 = rho_1

            do i = 0, l
                call self_linear_op_aV1_V2(u_vec(i), r_vec(i), -beta)
            end do

            call A_op%apply_operator(u_vec(l), u_vec(l+1), eps_r, transpose)

            call dot_product_V1_V2(r0_vec, u_vec(l+1), dot_product)

            alpha = rho_0 / dot_product

            do i = 0, l
                call self_linear_op_V1_aV2(r_vec(i), u_vec(i+1), -alpha)
            end do

            call A_op%apply_operator(r_vec(l), r_vec(l+1), eps_r, transpose)

            call self_linear_op_V1_aV2(f_vec_out, u_vec(0), alpha)

        end do
        !End of the L-term BICGStab iteration--------------------------------!

        !Beginning of the L-term BiCGStab iteration to compute gammas and tau!
        do l = 1, L_iter
            do i = 1, l-1

                call dot_product_V1_V2(r_vec(i), r_vec(l), dot_product)

                tau(i,l) = dot_product / sig(i)

                call self_linear_op_V1_aV2(r_vec(l), r_vec(i), -tau(i,l))

            end do

            call dot_product_V1_V2(r_vec(l), r_vec(l), dot_product)

            sig(l) = dot_product

            call dot_product_V1_V2(r_vec(l), r_vec(0), dot_product)

            gam_p(l) = dot_product / sig(l)

        end do
        !End of the L-term BiCGStab iteration to compute gammas and tau------!

        gam(L_iter) = gam_p(L_iter)
        omega       = gam(L_iter)

        do l = L_iter-1, 1, -1
            gam(l) = gam_p(l)
            do i = l+1, L_iter
                gam(l) = gam(l) - tau(l,i)*gam(i)
            end do
        end do

        do l =1, L_iter-1
            gam_pp(l) = gam(l+1)
            do i = l+1, L_iter-1
                gam_pp(l) = gam_pp(l) + tau(l,i)*gam(i+1)
            end do
        end do

        call self_linear_op_V1_aV2(f_vec_out, r_vec(0), gam(1))

        call self_linear_op_V1_aV2(r_vec(0), r_vec(L_iter), -gam_p(L_iter))

        call self_linear_op_V1_aV2(u_vec(0), u_vec(L_iter), -gam(L_iter))

        do l = 1, L_iter-1
            call self_linear_op_V1_aV2(u_vec(0), u_vec(l), -gam(l))
            call self_linear_op_V1_aV2(f_vec_out, r_vec(l), gam_pp(l))
            call self_linear_op_V1_aV2(r_vec(0), r_vec(l), -gam_p(l))
        end do

        call dot_product_V1_V2(r_vec(0), r_vec(0), dot_product)

        err_r = DSQRT(DBLE(dot_product))

        n = n + 1

        if (err_r <= tol*norm_t) converged = .true.

    end do
    !End of the BICGStab iteration-------------------------------------------!

    abs_err = err_r/norm_t

end subroutine BICGStab_L

!###################################################################################################
end module