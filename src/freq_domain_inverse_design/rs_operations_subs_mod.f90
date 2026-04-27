module rs_operations_subs_mod

#ifdef USE_MPI
    use mpi
#endif

    use constants_mod
    use rs_vec_base_mod
    use rs_vec_dimensions_mod
    
    implicit none

contains

!###################################################################################################

subroutine linear_op_V1_aV2(vec1, vec2, a, vec_out)

    class(TRSvec), intent(in) :: vec1
    class(TRSvec), intent(in) :: vec2
    class(TRSvec), intent(inout) :: vec_out
    complex(dp)       , intent(in) :: a

    select type(vec1)
    type is (TRSvec_1D)
    select type(vec2)
    type is (TRSvec_1D)
    select type(vec_out)
    type is (TRSvec_1D)

        vec_out%nx         = vec1%nx
        vec_out%dimensions = vec1%dimensions
        vec_out%freq       = vec1%freq
        vec_out%dr         = vec1%dr
    
        vec_out%pl_x = vec1%pl_x + a * vec2%pl_x
        vec_out%pl_y = vec1%pl_y + a * vec2%pl_y
        vec_out%mi_x = vec1%mi_x + a * vec2%mi_x
        vec_out%mi_y = vec1%mi_y + a * vec2%mi_y

    end select
    end select


    type is (TRSvec_2D)
    select type(vec2)
    type is (TRSvec_2D)
    select type(vec_out)
    type is (TRSvec_2D)

        vec_out%nx         = vec1%nx
        vec_out%ny         = vec1%ny
        vec_out%dimensions = vec1%dimensions
        vec_out%freq       = vec1%freq
        vec_out%dr         = vec1%dr

        vec_out%pl_x = vec1%pl_x + a * vec2%pl_x
        vec_out%pl_y = vec1%pl_y + a * vec2%pl_y
        vec_out%pl_z = vec1%pl_z + a * vec2%pl_z
        vec_out%mi_x = vec1%mi_x + a * vec2%mi_x
        vec_out%mi_y = vec1%mi_y + a * vec2%mi_y
        vec_out%mi_z = vec1%mi_z + a * vec2%mi_z

    end select
    end select

    type is (TRSvec_3D)
    select type(vec2)
    type is (TRSvec_3D)
    select type(vec_out)
    type is (TRSvec_3D)

        vec_out%nx         = vec1%nx
        vec_out%ny         = vec1%ny
        vec_out%nz         = vec1%nz
        vec_out%dimensions = vec1%dimensions
        vec_out%freq       = vec1%freq
        vec_out%dr         = vec1%dr

        vec_out%pl_x = vec1%pl_x + a * vec2%pl_x
        vec_out%pl_y = vec1%pl_y + a * vec2%pl_y
        vec_out%pl_z = vec1%pl_z + a * vec2%pl_z
        vec_out%mi_x = vec1%mi_x + a * vec2%mi_x
        vec_out%mi_y = vec1%mi_y + a * vec2%mi_y
        vec_out%mi_z = vec1%mi_z + a * vec2%mi_z
    
    end select
    end select

    end select

end subroutine linear_op_V1_aV2

!###################################################################################################


subroutine self_linear_op_V1_aV2(vec1, vec2, a)

    class(TRSvec), intent(inout) :: vec1
    class(TRSvec), intent(in) :: vec2
    complex(dp)       , intent(in) :: a

    select type(vec1)
    type is (TRSvec_1D)
    select type(vec2)
    type is (TRSvec_1D)

        vec1%pl_x = vec1%pl_x + a * vec2%pl_x
        vec1%pl_y = vec1%pl_y + a * vec2%pl_y
        vec1%mi_x = vec1%mi_x + a * vec2%mi_x
        vec1%mi_y = vec1%mi_y + a * vec2%mi_y

    end select

    type is (TRSvec_2D)
    select type(vec2)
    type is (TRSvec_2D)

        vec1%pl_x = vec1%pl_x + a * vec2%pl_x
        vec1%pl_y = vec1%pl_y + a * vec2%pl_y
        vec1%pl_z = vec1%pl_z + a * vec2%pl_z
        vec1%mi_x = vec1%mi_x + a * vec2%mi_x
        vec1%mi_y = vec1%mi_y + a * vec2%mi_y
        vec1%mi_z = vec1%mi_z + a * vec2%mi_z

    end select

    type is (TRSvec_3D)
    select type(vec2)
    type is (TRSvec_3D)

        vec1%pl_x = vec1%pl_x + a * vec2%pl_x
        vec1%pl_y = vec1%pl_y + a * vec2%pl_y
        vec1%pl_z = vec1%pl_z + a * vec2%pl_z
        vec1%mi_x = vec1%mi_x + a * vec2%mi_x
        vec1%mi_y = vec1%mi_y + a * vec2%mi_y
        vec1%mi_z = vec1%mi_z + a * vec2%mi_z
    
    end select

    end select


end subroutine self_linear_op_V1_aV2

!###################################################################################################

subroutine self_linear_op_aV1_V2(vec1, vec2, a)

    class(TRSvec), intent(inout) :: vec1
    class(TRSvec), intent(in)    :: vec2
    complex(dp)  , intent(in)    :: a

    select type(vec1)
    type is (TRSvec_1D)
    select type(vec2)
    type is (TRSvec_1D)

        vec1%pl_x = a * vec1%pl_x + vec2%pl_x
        vec1%pl_y = a * vec1%pl_y + vec2%pl_y
        vec1%mi_x = a * vec1%mi_x + vec2%mi_x
        vec1%mi_y = a * vec1%mi_y + vec2%mi_y

    end select

    type is (TRSvec_2D)
    select type(vec2)
    type is (TRSvec_2D)

        vec1%pl_x = a * vec1%pl_x + vec2%pl_x
        vec1%pl_y = a * vec1%pl_y + vec2%pl_y
        vec1%pl_z = a * vec1%pl_z + vec2%pl_z
        vec1%mi_x = a * vec1%mi_x + vec2%mi_x
        vec1%mi_y = a * vec1%mi_y + vec2%mi_y
        vec1%mi_z = a * vec1%mi_z + vec2%mi_z

    end select

    type is (TRSvec_3D)
    select type(vec2)
    type is (TRSvec_3D)

        vec1%pl_x = a * vec1%pl_x + vec2%pl_x
        vec1%pl_y = a * vec1%pl_y + vec2%pl_y
        vec1%pl_z = a * vec1%pl_z + vec2%pl_z
        vec1%mi_x = a * vec1%mi_x + vec2%mi_x
        vec1%mi_y = a * vec1%mi_y + vec2%mi_y
        vec1%mi_z = a * vec1%mi_z + vec2%mi_z
    
    end select

    end select

end subroutine self_linear_op_aV1_V2

!###################################################################################################

subroutine dot_product_V1_V2(vec1, vec2, dot_prod)

    class(TRSvec), intent(in) :: vec1
    class(TRSvec), intent(in) :: vec2
    complex(dp)  , intent(out) :: dot_prod

    integer     :: nx, ny, nz
    complex(dp) :: local_dot_prod
#ifdef USE_MPI
    integer     :: ierr
    complex(dp) :: global_dot_prod
#endif

    dot_prod = 0.0_dp

    select type(vec1)
    type is (TRSvec_1D)
    select type(vec2)
    type is (TRSvec_1D)

        nx = vec1%nx

        local_dot_prod = SUM(DCONJG(vec1%pl_x(1:nx))*vec2%pl_x(1:nx))
        local_dot_prod = local_dot_prod + SUM(DCONJG(vec1%pl_y(1:nx))*vec2%pl_y(1:nx))
        local_dot_prod = local_dot_prod + SUM(DCONJG(vec1%mi_x(1:nx))*vec2%mi_x(1:nx))
        local_dot_prod = local_dot_prod + SUM(DCONJG(vec1%mi_y(1:nx))*vec2%mi_y(1:nx))

    end select

    type is (TRSvec_2D)
    select type(vec2)
    type is (TRSvec_2D)

        nx = vec1%nx
        ny = vec1%ny

        local_dot_prod = SUM(DCONJG(vec1%pl_x(1:nx,1:ny))*vec2%pl_x(1:nx,1:ny))
        local_dot_prod = local_dot_prod + SUM(DCONJG(vec1%pl_y(1:nx,1:ny))*vec2%pl_y(1:nx,1:ny))
        local_dot_prod = local_dot_prod + SUM(DCONJG(vec1%pl_z(1:nx,1:ny))*vec2%pl_z(1:nx,1:ny))
        local_dot_prod = local_dot_prod + SUM(DCONJG(vec1%mi_x(1:nx,1:ny))*vec2%mi_x(1:nx,1:ny))
        local_dot_prod = local_dot_prod + SUM(DCONJG(vec1%mi_y(1:nx,1:ny))*vec2%mi_y(1:nx,1:ny))
        local_dot_prod = local_dot_prod + SUM(DCONJG(vec1%mi_z(1:nx,1:ny))*vec2%mi_z(1:nx,1:ny))

    end select

    type is (TRSvec_3D)
    select type(vec2)
    type is (TRSvec_3D)
        nx = vec1%nx
        ny = vec1%ny
        nz = vec1%nz

        local_dot_prod = SUM(DCONJG(vec1%pl_x(1:nx,1:ny,1:nz))*vec2%pl_x(1:nx,1:ny,1:nz))
        local_dot_prod = local_dot_prod + SUM(DCONJG(vec1%pl_y(1:nx,1:ny,1:nz))*vec2%pl_y(1:nx,1:ny,1:nz))
        local_dot_prod = local_dot_prod + SUM(DCONJG(vec1%pl_z(1:nx,1:ny,1:nz))*vec2%pl_z(1:nx,1:ny,1:nz))
        local_dot_prod = local_dot_prod + SUM(DCONJG(vec1%mi_x(1:nx,1:ny,1:nz))*vec2%mi_x(1:nx,1:ny,1:nz))
        local_dot_prod = local_dot_prod + SUM(DCONJG(vec1%mi_y(1:nx,1:ny,1:nz))*vec2%mi_y(1:nx,1:ny,1:nz))
        local_dot_prod = local_dot_prod + SUM(DCONJG(vec1%mi_z(1:nx,1:ny,1:nz))*vec2%mi_z(1:nx,1:ny,1:nz))

    end select

    end select

#ifdef USE_MPI
    call MPI_Allreduce(local_dot_prod, global_dot_prod, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    dot_prod = global_dot_prod
#else
    dot_prod = local_dot_prod
#endif

end subroutine dot_product_V1_V2

!###################################################################################################

subroutine copy_V2_on_V1(vec1, vec2)

    class(TRSvec), intent(inout) :: vec1
    class(TRSvec), intent(in)    :: vec2

    select type(vec1)
    type is (TRSvec_1D)
    select type(vec2)
    type is (TRSvec_1D)

        vec1%pl_x = vec2%pl_x
        vec1%pl_y = vec2%pl_y
        vec1%mi_x = vec2%mi_x
        vec1%mi_y = vec2%mi_y

    end select

    type is (TRSvec_2D)
    select type(vec2)
    type is (TRSvec_2D)

        vec1%pl_x = vec2%pl_x
        vec1%pl_y = vec2%pl_y
        vec1%pl_z = vec2%pl_z
        vec1%mi_x = vec2%mi_x
        vec1%mi_y = vec2%mi_y
        vec1%mi_z = vec2%mi_z

    end select

    type is (TRSvec_3D)
    select type(vec2)
    type is (TRSvec_3D)

        vec1%pl_x = vec2%pl_x
        vec1%pl_y = vec2%pl_y
        vec1%pl_z = vec2%pl_z
        vec1%mi_x = vec2%mi_x
        vec1%mi_y = vec2%mi_y
        vec1%mi_z = vec2%mi_z
    
    end select

    end select

end subroutine copy_V2_on_V1

!###################################################################################################

subroutine self_product_aV1(vec1, a)

    class(TRSvec), intent(inout) :: vec1
    complex(dp)  , intent(in)    :: a

    select type(vec1)
    type is (TRSvec_1D)

        vec1%pl_x = a * vec1%pl_x
        vec1%pl_y = a * vec1%pl_y
        vec1%mi_x = a * vec1%mi_x
        vec1%mi_y = a * vec1%mi_y

    type is (TRSvec_2D)


        vec1%pl_x = a * vec1%pl_x
        vec1%pl_y = a * vec1%pl_y
        vec1%pl_z = a * vec1%pl_z
        vec1%mi_x = a * vec1%mi_x
        vec1%mi_y = a * vec1%mi_y
        vec1%mi_z = a * vec1%mi_z

    type is (TRSvec_3D)

        vec1%pl_x = a * vec1%pl_x
        vec1%pl_y = a * vec1%pl_y
        vec1%pl_z = a * vec1%pl_z
        vec1%mi_x = a * vec1%mi_x
        vec1%mi_y = a * vec1%mi_y
        vec1%mi_z = a * vec1%mi_z
    
    end select

end subroutine self_product_aV1

!###################################################################################################
subroutine reset_V1_to_zero(vec1)

    class(TRSvec), intent(inout) :: vec1

    select type(vec1)
    type is (TRSvec_1D)

        vec1%pl_x = Z_0
        vec1%pl_y = Z_0
        vec1%mi_x = Z_0
        vec1%mi_y = Z_0
    
    type is (TRSvec_2D)

        vec1%pl_x = Z_0
        vec1%pl_y = Z_0
        vec1%pl_z = Z_0
        vec1%mi_x = Z_0
        vec1%mi_y = Z_0
        vec1%mi_z = Z_0

    type is (TRSvec_3D)

        vec1%pl_x = Z_0
        vec1%pl_y = Z_0
        vec1%pl_z = Z_0
        vec1%mi_x = Z_0
        vec1%mi_y = Z_0
        vec1%mi_z = Z_0
    
    end select

end subroutine reset_V1_to_zero

!###################################################################################################

end module rs_operations_subs_mod