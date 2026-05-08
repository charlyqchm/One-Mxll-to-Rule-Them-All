program ch_inv_design

    use constants_mod
    use medium_mod
    use optimization_problem_mod
    use allocator_multidim_mod
    use design_base_mod
    use design_1D_mod
    use design_2D_mod
    use design_3D_mod
    use bicgstab_mod
    use parallel_subs_mod
    use input_mod
    use output_mod

    implicit none

    type(TOptPrblm), allocatable :: opt_problems(:)
    class(TDesign) , allocatable :: design
    type(TOutput)   :: output
    
    logical         :: restart
    logical         :: print_restart
    logical         :: converge_optimization
    logical         :: iterate_optimization
    logical         :: change_delta_rho
    logical         :: converged_forward
    logical         :: converged_adjoint
    
    integer         :: grid_Ndims(3)
    integer         :: max_iter_steps
    integer         :: n_opt_problems
    integer         :: boundaries(3)
    integer         :: mpi_cart_comm
    integer         :: mpi_coords(3) = 0
    integer         :: mpi_dims(3)   = 1
    integer         :: myrank        = 0
    integer         :: mpi_nprocs
    integer         :: dimensions
    integer         :: n_pml
    integer         :: bicgstab_max_iter
    integer         :: bicgstab_L_term
    integer         :: print_field_list(6)
    integer         :: unit_restart = 300
    integer         :: iter_step
    integer         :: n_beta_steps
    integer         :: n_delta_rho_steps
    integer         :: n_accuracy_der
    integer         :: i, j, k, p
    
    real(dp)        :: fom, fom_old
    real(dp)        :: dr
    real(dp)        :: freq_list(100)
    real(dp)        :: eps_Re(100)
    real(dp)        :: eps_Im(100)
    real(dp)        :: delta_rho(50)
    real(dp)        :: beta(50)
    real(dp)        :: eta
    real(dp)        :: rho_init
    real(dp)        :: sigma_rho
    real(dp)        :: bicgstab_tol


    call read_input_file(boundaries, restart, converge_optimization, n_opt_problems,        &
                         max_iter_steps, dimensions, n_pml, grid_Ndims, dr, freq_list,      &
                         eps_Re, eps_Im, delta_rho, beta, eta, rho_init, sigma_rho, mpi_dims, &
                         n_delta_rho_steps, n_beta_steps, n_accuracy_der, bicgstab_max_iter, bicgstab_L_term,    &
                         bicgstab_tol, print_field_list, print_restart)

    mpi_nprocs = mpi_dims(1) * mpi_dims(2) * mpi_dims(3)

    call init_parallelization(dimensions, mpi_coords, mpi_dims, mpi_nprocs, mpi_cart_comm, &
                              boundaries, myrank)

    unit_restart = unit_restart + myrank

    call output%init_output(print_restart, print_field_list, myrank)

    allocate(opt_problems(n_opt_problems))

    do i = 1, n_opt_problems
        call opt_problems(i)%init_optprblm(i, dimensions, dr, freq_list(i), eps_Re(i),    &
                                           eps_Im(i), eta, boundaries, n_pml, grid_Ndims, &
                                           n_accuracy_der, mpi_cart_comm,  mpi_coords, mpi_dims)
        opt_problems(i)%beta_rho = beta(1)
    end do

    design = design_factory(dimensions)

    call design%init_design(dimensions, dr, sigma_rho, grid_Ndims)

    do i = 1, n_opt_problems
        call design%collect_opt_regions(opt_problems(i)%opt_region, rho_init)
    end do

    call extend_rho_to_ranks(design)
    call extend_opt_region_to_ranks(design)

    call init_BICGStab_L_variables(grid_Ndims, dimensions, n_accuracy_der, bicgstab_L_term)

    if (restart) then
        do p = 1, n_opt_problems
            call read_problem_from_restart_file(opt_problems(p), design, mpi_dims, &
                                                mpi_coords, p, unit_restart, n_opt_problems)
        end do
        
    else 
        
        do p = 1, n_opt_problems

            call opt_problems(p)%solve_forward(bicgstab_tol, bicgstab_max_iter, bicgstab_L_term, &
            0.0_dp, -1, converged_forward)
            
            call opt_problems(p)%update_targets()
            
            call opt_problems(p)%solve_adjoint(bicgstab_tol, bicgstab_max_iter, bicgstab_L_term, &
            0.0_dp, -1, converged_adjoint)
            
            call opt_problems(p)%update_fields()
            
        end do

    end if
        
    call design%apply_kernel_on_rho()
   
    do p = 1, n_opt_problems

        call opt_problems(p)%update_permittivity(design)

        call opt_problems(p)%solve_forward(bicgstab_tol, bicgstab_max_iter, bicgstab_L_term, &
                                           0.0_dp, 0, converged_forward)

        call opt_problems(p)%update_targets()

        call opt_problems(p)%solve_adjoint(bicgstab_tol, bicgstab_max_iter, bicgstab_L_term, &
                                           0.0_dp, 0, converged_adjoint)

        call opt_problems(p)%update_fields()

        call design%collect_FOM(opt_problems(p)%w_total, p, n_opt_problems)

    end do

    fom     = design%fom**2
    fom_old = fom

    call output%write_gral_output(0, fom, 0.0_dp, beta(1), design%grad_max, myrank)

    iterate_optimization = .true.
    iter_step = 1
    k = 1
    !Init optimization loop. --------------------------------------------------!
    do while(iterate_optimization .and. iter_step <= max_iter_steps)
        
        !Calculate the gradient ---------------------------------------!
        do p = 1, n_opt_problems
            opt_problems(p)%beta_rho = beta(k)
            call opt_problems(p)%compute_gradient(design)
            call design%collect_gradients(opt_problems(p)%grad, &
            p, n_opt_problems)
        end do
        
        call extend_grad_to_ranks(design)

        call design%apply_kernel_on_grad()
        !--------------------------------------------------------------!
        
        j = 1
        design%drho = delta_rho(j)
        change_delta_rho = .true.
        do while(change_delta_rho)
          
            call design%displace_rho()

            call extend_rho_to_ranks(design)

            call design%apply_kernel_on_rho()

            do p = 1, n_opt_problems
                
                call opt_problems(p)%update_permittivity(design)

                call opt_problems(p)%solve_forward(bicgstab_tol, bicgstab_max_iter, bicgstab_L_term, &
                                                   delta_rho(j), iter_step, converged_forward)
                
                call opt_problems(p)%update_targets()

                call design%collect_FOM(opt_problems(p)%w_total, p, n_opt_problems)                                   
            end do

            fom = design%fom**2

            if (converge_optimization) then
                change_delta_rho = (fom < fom_old) .or. (.not. converged_forward)
            else
                change_delta_rho = .not. converged_forward
            end if

            if (change_delta_rho .and. converge_optimization)then
                j = j + 1
                design%drho = delta_rho(j)
            end if

            !At this point we have used all the elements of delta_rho.
            !Since the FOM is not maximized, we change beta from the
            !binary function.
            if (j == n_delta_rho_steps .and. converge_optimization) then
                change_delta_rho = .false.
                design%drho = 0.0d0

                do p = 1, n_opt_problems

                    call output%write_restart(opt_problems(p)%f_vec, &
                        opt_problems(p)%f_adj_vec, opt_problems(p)%n_trg, &
                        design, p, mpi_dims, mpi_coords, beta(k))

                end do

                k = k + 1
                do p = 1, n_opt_problems
                    opt_problems(p)%beta_rho = beta(k)
                end do

                if (k > n_beta_steps) then
                    iterate_optimization = .false.
                else

                    call design%reset_rho_one_step_back()

                    call extend_rho_to_ranks(design)

                    call design%reset_grad()

                    call extend_grad_to_ranks(design)

                    call design%apply_kernel_on_rho()
                    do p = 1, n_opt_problems

                        call opt_problems(p)%update_permittivity(design)

                        call opt_problems(p)%solve_forward(bicgstab_tol, bicgstab_max_iter, &
                                                           bicgstab_L_term, delta_rho(j),   &
                                                           iter_step, converged_forward)
                    
                        call opt_problems(p)%update_targets()
                    end do
                end if

            else if (.not. converge_optimization .and. change_delta_rho) then
                iterate_optimization = .false.
                change_delta_rho     = .false.
            end if

        end do

        do p = 1, n_opt_problems

            call opt_problems(p)%solve_adjoint(bicgstab_tol, bicgstab_max_iter, bicgstab_L_term, &
                                               0.0_dp, iter_step, converged_adjoint)

            call design%collect_FOM(opt_problems(p)%w_total, p, n_opt_problems)

        end do

        fom = design%fom**2

        if (k <= n_beta_steps) then

            do p = 1, n_opt_problems
                call opt_problems(p)%update_fields()
            end do

            call design%update_grad()
            call design%update_rho()
            fom_old = fom

            call output%write_gral_output(iter_step, fom, delta_rho(j), beta(k), &
                                          design%grad_max, myrank)

        end if

        iter_step = iter_step + 1

    end do
    ! End of optimization loop. --------------------------------------------------!


    do p = 1, n_opt_problems

        call output%write_final_structure(opt_problems(p)%eps_r, grid_Ndims, dr, p, &
                                          mpi_coords, mpi_dims)

    end do

    do p = 1, n_opt_problems

        call output%write_final_fields(opt_problems(p)%f_vec, grid_Ndims, p, &
                                       dr, mpi_coords, mpi_dims) 

    end do

    call kill_BICGStab_L_variables()

    call design%kill_design()

    do p = 1, n_opt_problems
        call opt_problems(p)%kill_optprblm()
    end do

    if (allocated(opt_problems)) deallocate(opt_problems)

    call output%close_output()

    call finalize_parallelization()

end program ch_inv_design