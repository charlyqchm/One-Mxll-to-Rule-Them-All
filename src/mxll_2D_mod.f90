module mxll_2D_mod

    use constants_mod
    use mxll_base_mod
    use classical_medium_mod

    implicit none

    type, extends(TMxll) :: TMxll_2D

        type(TClassicalMedium), allocatable :: media(:)

        integer               :: nx,ny
        integer               :: boundaries(2)
        integer               :: n_cpml_sections
        logical               :: cpml_pos(4) !Indicates if the rank has some of the
                                               ! 4 possible CPML boundaries.
        integer      , allocatable :: media_map(:,:,:)
        real(dp)     , allocatable :: Ex(:,:), Ey(:,:), Ez(:,:)
        real(dp)     , allocatable :: Ex_old(:,:), Ey_old(:,:), Ez_old(:,:)
        real(dp)     , allocatable :: Hx(:,:), Hy(:,:), Hz(:,:)
        real(dp)     , allocatable :: eps_z(:,:), eps_y(:,:), eps_x(:,:)
        real(dp)     , allocatable :: Jx(:,:), Jy(:,:), Jz(:,:)
        real(dp)     , allocatable :: Jx_old(:,:), Jy_old(:,:), Jz_old(:,:)
        real(dp)     , allocatable :: dJx(:,:), dJy(:,:), dJz(:,:)
        real(dp)     , allocatable :: den_ex(:), den_ey(:)
        real(dp)     , allocatable :: den_hx(:), den_hy(:)
        !PML STUFF
        real(dp)     , allocatable :: psix_Ezx(:,:)
        real(dp)     , allocatable :: psix_Eyx(:,:)
        real(dp)     , allocatable :: psiy_Exy(:,:)
        real(dp)     , allocatable :: psiy_Ezy(:,:)
        real(dp)     , allocatable :: psix_Hyx(:,:)
        real(dp)     , allocatable :: psix_Hzx(:,:)
        real(dp)     , allocatable :: psiy_Hxy(:,:)
        real(dp)     , allocatable :: psiy_Hzy(:,:)
        real(dp)     , allocatable :: be(:)
        real(dp)     , allocatable :: ce(:)
        real(dp)     , allocatable :: bh(:)
        real(dp)     , allocatable :: ch(:)

        !Medium stuff
        real(dp)     , allocatable :: PDx(:,:), PDy(:,:), PDz(:,:)
        real(dp)     , allocatable :: PLx(:,:,:), PLy(:,:,:), PLz(:,:,:)
        real(dp)     , allocatable :: PLx_old(:,:,:), PLy_old(:,:,:), PLz_old(:,:,:)

    contains
        procedure :: init => init_2Dgrid
        procedure :: kill => kill_2Dgrid
        procedure :: td_propagate_H_field => td_propagate_H_2Dfield
        procedure :: td_propagate_E_field => td_propagate_E_2Dfield 

    end type TMxll_2D

contains

!###################################################################################################

    subroutine init_2Dgrid(this, grid_Ndims, npml, boundaries, dt, dr, mode, n_media, &
                           mpi_coords, mpi_dims)
        
        class(TMxll_2D), intent(inout) :: this
        integer        , intent(in)    :: grid_Ndims(3)
        integer        , intent(in)    :: npml
        integer        , intent(in)    :: boundaries(3)
        integer        , intent(in)    :: mode
        integer        , intent(in)    :: n_media
        integer        , intent(in)    :: mpi_coords(3)
        integer        , intent(in)    :: mpi_dims(3)
        real(dp)       , intent(in)    :: dt
        real(dp)       , intent(in)    :: dr

        integer  :: nx, ny
        integer  :: kk, k
        integer  :: n_max_poles
        real(dp) :: sig_to_au
        real(dp) :: sigmaCPML
        integer  :: n_sec
        logical  :: rank_with_cpml

        real(dp), allocatable :: alphae(:), sige(:), kappae(:)
        real(dp), allocatable :: alphah(:), sigh(:), kappah(:)

        sig_to_au =  C_to_au**2*sec_to_au/(kg_to_au*m_to_au**3)

        nx = grid_Ndims(1)
        ny = grid_Ndims(2)
      
        this%nx            = nx
        this%ny            = ny
        this%npml          = npml
        this%mode          = mode
        this%n_media       = n_media
        this%boundaries(1) = boundaries(1)
        this%boundaries(2) = boundaries(2)
        this%dr            = dr
        this%dt_eps0       = dt/eps0
        this%dt_mu0        = dt/mu0
        this%dt            = dt

        this%chunk_coor = mpi_coords

        if (.not. allocated(this%den_ex))   allocate(this%den_ex(nx))
        if (.not. allocated(this%den_hx))   allocate(this%den_hx(nx))
        if (.not. allocated(this%den_ey))   allocate(this%den_ey(ny))
        if (.not. allocated(this%den_hy))   allocate(this%den_hy(ny))

        this%den_ex = 1.0d0 / dr
        this%den_ey = 1.0d0 / dr
        this%den_hx = 1.0d0 / dr
        this%den_hy = 1.0d0 / dr

        if (this%mode == TEZ_2D_MODE .or. this%mode == FULL_2D_MODE) then

            if (.not. allocated(this%Ex))      allocate(this%Ex(0:nx, ny+1))
            if (.not. allocated(this%Ey))      allocate(this%Ey(nx+1, 0:ny))
            if (.not. allocated(this%Ex_old))  allocate(this%Ex_old(0:nx, ny+1))
            if (.not. allocated(this%Ey_old))  allocate(this%Ey_old(nx+1, 0:ny))
            if (.not. allocated(this%Hz))      allocate(this%Hz(0:nx, 0:ny))
            if (.not. allocated(this%eps_x))   allocate(this%eps_x(nx,ny))
            if (.not. allocated(this%eps_y))   allocate(this%eps_y(nx,ny))
            if (.not. allocated(this%Jx))      allocate(this%Jx(nx+1, ny))
            if (.not. allocated(this%Jx_old))  allocate(this%Jx_old(nx+1, ny))
            if (.not. allocated(this%dJx))     allocate(this%dJx(nx+1, ny))
            if (.not. allocated(this%Jy))      allocate(this%Jy(nx, ny+1))
            if (.not. allocated(this%Jy_old))  allocate(this%Jy_old(nx, ny+1))
            if (.not. allocated(this%dJy))     allocate(this%dJy(nx, ny+1))

            this%Ex      = 0.0d0
            this%Ex_old  = 0.0d0
            this%Ey      = 0.0d0
            this%Ey_old  = 0.0d0
            this%Hz      = 0.0d0
            this%Jx      = 0.0d0 
            this%Jx_old  = 0.0d0  
            this%dJx     = 0.0d0  
            this%Jy      = 0.0d0  
            this%Jy_old  = 0.0d0  
            this%dJy     = 0.0d0  
            this%eps_x   = 1.0d0
            this%eps_y   = 1.0d0

        elseif (this%mode == TMZ_2D_MODE .or. this%mode == FULL_2D_MODE) then

            if (.not. allocated(this%Ez))      allocate(this%Ez(nx+1, ny+1))
            if (.not. allocated(this%Ez_old))  allocate(this%Ez_old(nx+1, ny+1))
            if (.not. allocated(this%Hx))      allocate(this%Hx(nx  , 0:ny))
            if (.not. allocated(this%Hy))      allocate(this%Hy(0:nx, ny  ))
            if (.not. allocated(this%eps_z))   allocate(this%eps_z(nx,ny))
            if (.not. allocated(this%Jz))      allocate(this%Jz(nx, ny))
            if (.not. allocated(this%Jz_old))  allocate(this%Jz_old(nx, ny))
            if (.not. allocated(this%dJz))     allocate(this%dJz(nx, ny))

            this%Ez      = 0.0d0
            this%Ez_old  = 0.0d0
            this%Hx      = 0.0d0
            this%Hy      = 0.0d0
            this%Jz      = 0.0d0
            this%Jz_old  = 0.0d0  
            this%dJz     = 0.0d0
            this%eps_z   = 1.0d0
        else
            write (*, '("The selected mode for 2D Maxwell does not exist")')
            stop
        end if

        call read_init_media(this%media, n_media, this%eps_x, this%eps_y, this%eps_z, grid_Ndims, dr, &
                             this%media_map, this%nx, this%ny, dt, mpi_coords, mpi_dims)


!Polarization grids are allocated even in the absence of a medium.
!This avoids unnecessary checks during the time propagation loop.
!TO-DO: Optimize memory usage later if needed.
        n_max_poles = 1

        do k=1,this%n_media
            if (this%media(k)%medium_type == DL_MEDIUM) then
                if (this%media(k)%n_poles > n_max_poles) then
                    n_max_poles = this%media(k)%n_poles
                end if
            end if
        end do                        

        if (.not. allocated(this%PDx))    allocate(this%PDx(nx, ny))
        if (.not. allocated(this%PDy))    allocate(this%PDy(nx, ny))
        if (.not. allocated(this%PDz))    allocate(this%PDz(nx, ny))
        this%PDx = 0.0d0
        this%PDy = 0.0d0
        this%PDz = 0.0d0

        if (.not. allocated(this%PLx))      allocate(this%PLx(nx, ny, n_max_poles))
        if (.not. allocated(this%PLy))      allocate(this%PLy(nx, ny, n_max_poles))
        if (.not. allocated(this%PLz))      allocate(this%PLz(nx, ny, n_max_poles))
        if (.not. allocated(this%PLx_old))  allocate(this%PLx_old(nx, ny, n_max_poles))
        if (.not. allocated(this%PLy_old))  allocate(this%PLy_old(nx, ny, n_max_poles))
        if (.not. allocated(this%PLz_old))  allocate(this%PLz_old(nx, ny, n_max_poles))
        this%PLx = 0.0d0
        this%PLy = 0.0d0
        this%PLz = 0.0d0
        this%PLx_old = 0.0d0
        this%PLy_old = 0.0d0
        this%PLz_old = 0.0d0

        !In case of using MPI, just the ranks on the borders of the box can have CPML
        !boundaries conditions.

#ifdef USE_MPI
        this%cpml_pos(1) = (mpi_coords(1) == 0) .and.&
                            this%boundaries(1) == CPML_BOUNDARIES
        this%cpml_pos(2) = (mpi_coords(1) == mpi_dims(1)-1) .and.&
                            this%boundaries(1) == CPML_BOUNDARIES
        this%cpml_pos(3) = (mpi_coords(2) == 0) .and.&
                            this%boundaries(2) == CPML_BOUNDARIES
        this%cpml_pos(4) = (mpi_coords(2) == mpi_dims(2)-1) .and.&
                            this%boundaries(2) == CPML_BOUNDARIES
        
        n_sec                = 1
        this%n_cpml_sections = n_sec
#else
        this%cpml_pos(1) = this%boundaries(1) == CPML_BOUNDARIES
        this%cpml_pos(2) = this%boundaries(1) == CPML_BOUNDARIES
        this%cpml_pos(3) = this%boundaries(2) == CPML_BOUNDARIES
        this%cpml_pos(4) = this%boundaries(2) == CPML_BOUNDARIES
    
        n_sec                = 2
        this%n_cpml_sections = n_sec
#endif
        
        if (this%cpml_pos(1) .or. this%cpml_pos(2)) then
            
            if (this%mode == TEZ_2D_MODE .or. this%mode == FULL_2D_MODE) then    
                if (.not. allocated(this%psix_Eyx)) allocate(this%psix_Eyx(n_sec*npml, ny))
                if (.not. allocated(this%psix_Hzx)) allocate(this%psix_Hzx(n_sec*(npml-1), ny))
                this%psix_Eyx = 0.0d0
                this%psix_Hzx = 0.0d0
            end if
            
            if (this%mode == TMZ_2D_MODE .or. this%mode == FULL_2D_MODE) then 
                if (.not. allocated(this%psix_Hyx)) allocate(this%psix_Hyx(n_sec*(npml-1), ny))
                if (.not. allocated(this%psix_Ezx)) allocate(this%psix_Ezx(n_sec*npml, ny))
                this%psix_Ezx = 0.0d0
                this%psix_Hyx = 0.0d0
            end if
            
        end if
        
        if (this%cpml_pos(3) .or. this%cpml_pos(4)) then      
            
            if (this%mode == TEZ_2D_MODE .or. this%mode == FULL_2D_MODE) then
                if (.not. allocated(this%psiy_Exy)) allocate(this%psiy_Exy(nx, n_sec*npml))
                if (.not. allocated(this%psiy_Hzy)) allocate(this%psiy_Hzy(nx, n_sec*(npml-1)))
                this%psiy_Exy = 0.0d0
                this%psiy_Hzy = 0.0d0
            end if 
            
            if (this%mode == TMZ_2D_MODE .or. this%mode == FULL_2D_MODE) then
                if (.not. allocated(this%psiy_Ezy)) allocate(this%psiy_Ezy(nx, n_sec*npml))
                if (.not. allocated(this%psiy_Hxy)) allocate(this%psiy_Hxy(nx, n_sec*(npml-1)))
                this%psiy_Ezy = 0.0d0
                this%psiy_Hxy = 0.0d0
            end if 
            
        end if
        
        rank_with_cpml = this%cpml_pos(1) .or. this%cpml_pos(2) .or. &
                         this%cpml_pos(3) .or. this%cpml_pos(4)

        if (rank_with_cpml) then
            if (.not. allocated(this%ce))        allocate(this%ce(npml))
            if (.not. allocated(this%be))        allocate(this%be(npml))
            if (.not. allocated(this%ch))        allocate(this%ch(npml))
            if (.not. allocated(this%bh))        allocate(this%bh(npml))
            if (.not. allocated(alphae))         allocate(alphae(npml))
            if (.not. allocated(alphah))         allocate(alphah(npml-1))
            if (.not. allocated(sige))           allocate(sige(npml))
            if (.not. allocated(sigh))           allocate(sigh(npml-1))
            if (.not. allocated(kappae))         allocate(kappae(npml))
            if (.not. allocated(kappah))         allocate(kappah(npml-1))

            sigmaCPML = 0.8 * (m + 1) / (dr * (mu0 / eps0)**0.5)

            do kk = 1, npml
                sige(kk)   = sigmaCPML * ((npml - kk) / (npml - 1.0))**m
                alphae(kk) = alphaCPML * sig_to_au * ((kk - 1) / (npml - 1.0))**ma
                kappae(kk) = 1.0 + (kappaCPML - 1.0) * ((npml - kk) / (npml - 1.0))**m
                this%be(kk)     = exp(-(sige(kk) / kappae(kk) + alphae(kk)) * dt / eps0)
                if ((sige(kk)==M_ZERO) .and.  &
                        (alphae(kk)==M_ZERO) .and.  &
                        (kk==npml)) then
                    this%ce(kk) = 0.0
                else
                    this%ce(kk) = sige(kk) * (this%be(kk) - 1.0) / &
                            (sige(kk) + kappae(kk) * alphae(kk)) / kappae(kk)
                end if
            end do
            
            do kk = 1, npml - 1
                sigh(kk)   = sigmaCPML * ((npml - kk - 0.5) / (npml - 1.0))**m
                alphah(kk) = alphaCPML * sig_to_au * ((kk - 0.5) / (npml - 1.0))**ma
                kappah(kk) = 1.0 + (kappaCPML - 1.0) * ((npml - kk - 0.5) / (npml - 1.0))**m
                this%bh(kk) = exp(-(sigh(kk) / kappah(kk) + alphah(kk)) * dt / eps0)
                this%ch(kk) = sigh(kk) * (this%bh(kk) - 1.0) / &
                        (sigh(kk) + kappah(kk) * alphah(kk)) / kappah(kk)
            end do

            if(this%cpml_pos(1)) then
                do k = 1, nx
                    if (k<=npml) then
                        this%den_ex(k) = 1.0 / (kappae(k) * dr)
                    endif
                end do

                do k = 1, nx
                    if (k<=(npml - 1)) then
                        this%den_hx(k) = 1.0 / (kappah(k) * dr)
                    endif
                end do
            end if

            if(this%cpml_pos(2)) then    
                kk = npml
                do k = 1, nx-1
                    if (k>=(nx + 1 - npml)) then
                        this%den_ex(k) = 1.0 / (kappae(kk) * dr)
                        kk = kk - 1
                    endif
                end do

                kk = npml - 1
                do k = 1, nx-1
                    if (k>=(nx + 1 - npml)) then
                        this%den_hx(k) = 1.0 / (kappah(kk) * dr)
                        kk = kk - 1
                    endif
                end do
            end if

            if(this%cpml_pos(3)) then

                do k = 1, ny
                    if (k<=npml) then
                        this%den_ey(k) = 1.0 / (kappae(k) * dr)
                    endif
                end do
                
                do k = 1, ny
                    if (k<=(npml - 1)) then
                        this%den_hy(k) = 1.0 / (kappah(k) * dr)
                    endif
                end do

            end if

            if (this%cpml_pos(4)) then
            
                kk = npml
                do k = 1, ny-1
                    if (k>=(ny + 1 - npml)) then
                        this%den_ey(k) = 1.0 / (kappae(kk) * dr)
                        kk = kk - 1
                    endif
                enddo

                kk = npml - 1
                do k = 1, ny-1
                    if (k>=(ny + 1 - npml)) then
                        this%den_hy(k) = 1.0 / (kappah(kk) * dr)
                        kk = kk - 1
                    endif
                enddo
            end if

        end if

        if (allocated(alphae)) deallocate(alphae)
        if (allocated(alphah)) deallocate(alphah)
        if (allocated(sige))   deallocate(sige)
        if (allocated(sigh))   deallocate(sigh)
        if (allocated(kappae)) deallocate(kappae)
        if (allocated(kappah)) deallocate(kappah)

    end subroutine init_2Dgrid

!###################################################################################################

    subroutine kill_2Dgrid(this)

        class(TMxll_2D), intent(inout) :: this
        integer :: i

        if (allocated(this%media)) then
            do i=1,size(this%media)
                call this%media(i)%kill_classical_medium()
            end do
            deallocate(this%media)
        end if

        if (allocated(this%Ex))       deallocate(this%Ex)
        if (allocated(this%Ey))       deallocate(this%Ey)
        if (allocated(this%Ez))       deallocate(this%Ez)
        if (allocated(this%Ex_old))   deallocate(this%Ex_old)
        if (allocated(this%Ey_old))   deallocate(this%Ey_old)
        if (allocated(this%Ez_old))   deallocate(this%Ez_old)
        if (allocated(this%Hx))       deallocate(this%Hx)
        if (allocated(this%Hy))       deallocate(this%Hy)
        if (allocated(this%Hz))       deallocate(this%Hz)
        if (allocated(this%eps_x))    deallocate(this%eps_x)
        if (allocated(this%eps_y))    deallocate(this%eps_y)
        if (allocated(this%eps_z))    deallocate(this%eps_z)
        if (allocated(this%psix_Eyx)) deallocate(this%psix_Eyx)
        if (allocated(this%psix_Ezx)) deallocate(this%psix_Ezx)
        if (allocated(this%psix_Hyx)) deallocate(this%psix_Hyx)
        if (allocated(this%psix_Hzx)) deallocate(this%psix_Hzx)
        if (allocated(this%psiy_Exy)) deallocate(this%psiy_Exy)
        if (allocated(this%psiy_Ezy)) deallocate(this%psiy_Ezy)
        if (allocated(this%psiy_Hxy)) deallocate(this%psiy_Hxy)
        if (allocated(this%psiy_Hzy)) deallocate(this%psiy_Hzy)
        if (allocated(this%Jx))       deallocate(this%Jx)
        if (allocated(this%Jy))       deallocate(this%Jy)
        if (allocated(this%Jz))       deallocate(this%Jz)
        if (allocated(this%Jx_old))   deallocate(this%Jx_old)
        if (allocated(this%Jy_old))   deallocate(this%Jy_old)
        if (allocated(this%Jz_old))   deallocate(this%Jz_old)
        if (allocated(this%dJx))      deallocate(this%dJx)
        if (allocated(this%dJy))      deallocate(this%dJy)
        if (allocated(this%dJz))      deallocate(this%dJz)
        if (allocated(this%den_ex))   deallocate(this%den_ex)
        if (allocated(this%den_ey))   deallocate(this%den_ey)
        if (allocated(this%den_hx))   deallocate(this%den_hx)
        if (allocated(this%den_hy))   deallocate(this%den_hy)     
        if (allocated(this%ce))       deallocate(this%ce)
        if (allocated(this%be))       deallocate(this%be)
        if (allocated(this%ch))       deallocate(this%ch)
        if (allocated(this%bh))       deallocate(this%bh)
        if (allocated(this%media_map))deallocate(this%media_map)
        if (allocated(this%PDx))      deallocate(this%PDx)
        if (allocated(this%PDy))      deallocate(this%PDy)
        if (allocated(this%PDz))      deallocate(this%PDz)
        if (allocated(this%PLx))      deallocate(this%PLx)
        if (allocated(this%PLy))      deallocate(this%PLy)
        if (allocated(this%PLz))      deallocate(this%PLz)
        if (allocated(this%PLx_old))  deallocate(this%PLx_old)
        if (allocated(this%PLy_old))  deallocate(this%PLy_old)
        if (allocated(this%PLz_old))  deallocate(this%PLz_old)


    end subroutine kill_2Dgrid   

!###################################################################################################
    subroutine td_propagate_H_2Dfield(this)

        class(TMxll_2D), intent(inout) :: this

        integer  :: nx, ny
        integer  :: npml
        integer  :: n_base
        integer  :: n_sec
        integer  :: i,j, ii, jj
        real(dp) :: rotE

        nx    = this%nx
        ny    = this%ny
        npml  = this%npml
        n_sec = this%n_cpml_sections

        if (n_sec==1) n_base = 0
        if (n_sec==2) n_base = npml-1

        if (this%mode == TMZ_2D_MODE .or. this%mode == FULL_2D_MODE) then


#ifdef USE_MPI
       !MPI subroutines will handle the periodic boundaries.
#else
            if (this%boundaries(1) == PERIODIC_BOUNDARIES) then
                this%Ez(nx+1, :)     = this%Ez(1, :)
            end if

            if (this%boundaries(2) == PERIODIC_BOUNDARIES) then
                this%Ez(:, ny+1)     = this%Ez(:, 1)
            end if
#endif 

            !$omp parallel default(shared) private(i, j, rotE)
            !$omp do collapse(2) schedule(static)
            do i = 1, nx
            do j = 1, ny
                rotE          = (this%Ez(i, j) - this%Ez(i, j + 1)) * this%den_hy(j)
                this%Hx(i, j) = this%Hx(i, j) + this%dt_mu0 * rotE
            end do
            end do
            !$omp end do nowait

            !$omp do collapse(2) schedule(static)
            do i = 1, nx
            do j = 1, ny
                rotE          = (this%Ez(i + 1, j) - this%Ez(i, j)) * this%den_hx(i)
                this%Hy(i, j) = this%Hy(i, j) + this%dt_mu0 * rotE
            end do
            end do
            !$omp end do
            !$omp end parallel

            if (this%cpml_pos(1)) then
                do j = 1, ny
                do i = 1, npml - 1
                    rotE                = (this%Ez(i+1, j)-this%Ez(i, j))/this%dr
                    this%psix_Hyx(i, j) = this%bh(i)*this%psix_Hyx(i, j) + &
                                        this%ch(i)*rotE
                    this%Hy(i, j)       = this%Hy(i, j) + this%dt_mu0*this%psix_Hyx(i, j)
                end do
                end do
            end if

            if (this%cpml_pos(2)) then
                do j = 1, ny
                    ii = n_sec*(npml - 1)
                    do i = nx + 1 - npml, nx - 1
                        rotE                 = (this%Ez(i+1, j) - this%Ez(i, j)) / this%dr
                        this%psix_Hyx(ii, j) = this%bh(ii-n_base)*this%psix_Hyx(ii, j) + &
                                            this%ch(ii-n_base)*rotE

                        this%Hy(i, j) = this%Hy(i, j) + this%dt_mu0*this%psix_Hyx(ii, j)
                        ii = ii - 1
                    end do
                end do
            end if

            if (this%cpml_pos(3)) then
                do i = 1, nx
                do j = 1, npml - 1
                    rotE                = (this%Ez(i, j) - this%Ez(i, j+1))/this%dr
                    this%psiy_Hxy(i, j) = this%bh(j)*this%psiy_Hxy(i, j) + this%ch(j)*rotE

                    this%Hx(i, j) = this%Hx(i, j) + this%dt_mu0 * this%psiy_Hxy(i, j)
                end do
                end do
            end if

            if (this%cpml_pos(4)) then
                do i = 1, nx
                    jj = n_sec*(npml-1)
                    do j = ny + 1 - npml, ny - 1
                        rotE                 = (this%Ez(i, j) - this%Ez(i, j+1)) / this%dr
                        this%psiy_Hxy(i, jj) = this%bh(jj-n_base)*this%psiy_Hxy(i, jj) + &
                                               this%ch(jj-n_base)*rotE

                        this%Hx(i, j) = this%Hx(i, j) + this%dt_mu0*this%psiy_Hxy(i, jj)
                        jj = jj - 1
                    enddo
                enddo  
            end if

        end if

        if (this%mode == TEZ_2D_MODE .or. this%mode == FULL_2D_MODE) then

#ifdef USE_MPI
       !MPI subroutines will handle the periodic boundaries.
#else
            if (this%boundaries(1) == PERIODIC_BOUNDARIES) then
                this%Ey(nx+1, :)     = this%Ey(1, :)
            end if

            if (this%boundaries(2) == PERIODIC_BOUNDARIES) then
                this%Ex(:, ny+1)     = this%Ex(:, 1)
            end if
#endif  

            !$omp parallel default(shared) private(i, j, rotE)
            !$omp do collapse(2) schedule(static)
            do i = 1, nx
            do j = 1, ny
                rotE          = (this%Ey(i, j)-this%Ey(i + 1, j))*this%den_hx(i) + &
                                (this%Ex(i, j + 1)-this%Ex(i, j))*this%den_hy(j)
                this%Hz(i, j) = this%Hz(i, j) + this%dt_mu0 * rotE
            end do
            end do
            !$omp end do
            !$omp end parallel

            if (this%cpml_pos(1)) then
                do j = 1, ny
                do i = 1, npml - 1
                    rotE                = (this%Ey(i, j)-this%Ey(i+1, j))/this%dr
                    this%psix_Hzx(i, j) = this%bh(i) * this%psix_Hzx(i, j) + &
                                          this%ch(i) * rotE 

                    this%Hz(i, j) = this%Hz(i, j) + this%dt_mu0 * this%psix_Hzx(i, j)
                end do
                end do
            end if

            if (this%cpml_pos(2)) then
                do j = 1, ny
                    ii = n_sec*(npml - 1)
                    do i = nx + 1 - npml, nx - 1
                        rotE                 = (this%Ey(i, j) - this%Ey(i+1, j))/this%dr 
                        this%psix_Hzx(ii, j) = this%bh(ii-n_base) * this%psix_Hzx(ii, j) + &
                                               this%ch(ii-n_base) * rotE

                        this%Hz(i, j) = this%Hz(i, j) + this%dt_mu0 * this%psix_Hzx(ii, j)
                        ii = ii - 1

                    end do
                end do
            end if

            if (this%cpml_pos(3)) then
                do i = 1, nx
                do j = 1, npml - 1
                    rotE                = (this%Ex(i, j+1) - this%Ex(i, j))/this%dr
                    this%psiy_Hzy(i, j) = this%bh(j) * this%psiy_Hzy(i, j) + &
                                          this%ch(j) * rotE

                    this%Hz(i, j) = this%Hz(i, j) + this%dt_mu0 * this%psiy_Hzy(i, j)
                enddo
                enddo
            end if

            if (this%cpml_pos(4)) then
                do i = 1, nx
                    jj = n_sec*(npml-1)
                    do j = ny + 1 - npml, ny - 1
                        rotE                 = (this%Ex(i, j+1) - this%Ex(i, j)) / this%dr
                        this%psiy_Hzy(i, jj) = this%bh(jj-n_base) * this%psiy_Hzy(i, jj) + &
                                               this%ch(jj-n_base) * rotE

                        this%Hz(i, j) = this%Hz(i, j) + this%dt_mu0 * this%psiy_Hzy(i, jj)
                        jj = jj - 1
                    enddo
                enddo
            end if
        end if        


    end subroutine td_propagate_H_2Dfield
!###################################################################################################

    subroutine td_propagate_E_2Dfield(this, t)

        class(TMxll_2D), intent(inout) :: this
        integer, intent(in)            :: t

        integer  :: nx, ny
        integer  :: npml
        integer  :: n_base
        integer  :: n_sec
        integer  :: i,j, ii, jj
        logical  :: no_medium
        real(dp) :: rotH, Jav_x, Jav_y, Jav_z

        nx    = this%nx
        ny    = this%ny
        npml  = this%npml
        n_sec = this%n_cpml_sections

        if (n_sec==1) n_base = 0
        if (n_sec==2) n_base = npml

        if (this%mode == TMZ_2D_MODE .or. this%mode == FULL_2D_MODE) then

#ifdef USE_MPI
       !MPI subroutines will handle the periodic boundaries.
#else
            if (this%boundaries(1) == PERIODIC_BOUNDARIES) then
                this%Hy(0, :)     = this%Hy(nx, :)
            end if

            if (this%boundaries(2) == PERIODIC_BOUNDARIES) then
                this%Hx(:, 0)     = this%Hx(:, ny)
            end if
#endif

            !$omp parallel default(shared) private(i, j, rotH, Jav_z, no_medium)
            !$omp do collapse(2) schedule(static)
            do j = 1, ny
            do i = 1, nx
                rotH = (this%Hy(i, j) - this%Hy(i - 1, j)) * this%den_ex(i) +  &
                       (this%Hx(i, j - 1) - this%Hx(i, j)) * this%den_ey(j)

                no_medium = .false.
                call get_medium_polarization(this%media, this%media_map(i,j,3), &
                                            this%PDz(i,j), this%PLz(i,j,:), this%PLz_old(i,j,:), &
                                            this%Ez(i,j), this%Ez_old(i,j), rotH, no_medium)

                if (no_medium) then

                    Jav_z         = this%Jz_old(i,j) + (this%time - this%t_skip) * this%dJz(i,j)
                    this%Ez(i, j) = this%Ez(i, j) + this%dt_eps0 * this%eps_z(i, j) * rotH - &
                                    this%dt_eps0 * Jav_z
                end if

            end do
            end do
            !$omp end do
            !$omp end parallel

            if (this%cpml_pos(1)) then
                do j = 1, ny
                do i = 2, npml
                    rotH  = (this%Hy(i, j) - this%Hy(i-1, j)) / this%dr 
                    this%psix_Ezx(i, j) = this%be(i) * this%psix_Ezx(i, j) + this%ce(i) * rotH
                    this%Ez(i, j) = this%Ez(i, j) + this%dt_eps0 * this%psix_Ezx(i, j)
                end do
                end do
            end if

            if (this%cpml_pos(2)) then
                do j = 1, ny
                    ii = n_sec*npml
                    do i = nx + 1 - npml, nx - 1
                        rotH = (this%Hy(i, j) - this%Hy(i-1, j)) / this%dr
                        this%psix_Ezx(ii, j) = this%be(ii-n_base) * this%psix_Ezx(ii, j) + &
                                               this%ce(ii-n_base) * rotH 

                        this%Ez(i, j) = this%Ez(i, j) + this%dt_eps0 * this%psix_Ezx(ii, j)
                        ii = ii - 1
                    enddo
                enddo
            end if

            if (this%cpml_pos(3)) then
                do i = 1, nx
                do j = 2, npml
                    rotH = (this%Hx(i, j-1) - this%Hx(i, j)) / this%dr
                    this%psiy_Ezy(i, j) = this%be(j) * this%psiy_Ezy(i, j) + &
                                          this%ce(j) * rotH

                    this%Ez(i, j) = this%Ez(i, j) + this%dt_eps0 * this%psiy_Ezy(i, j)
                end do
                end do
            end if

            if (this%cpml_pos(4)) then
                do i = 1, nx
                    jj = n_sec*npml
                    do j = ny + 1 - npml, ny - 1
                        rotH = (this%Hx(i, j-1) - this%Hx(i, j)) / this%dr
                        this%psiy_Ezy(i, jj) = this%be(jj-n_base) * this%psiy_Ezy(i, jj) + &
                                this%ce(jj-n_base) * rotH

                        this%Ez(i, j) = this%Ez(i, j) + this%dt_eps0 * this%psiy_Ezy(i, jj)
                        jj = jj - 1

                    end do
                end do    
            end if

        end if

        if (this%mode == TEZ_2D_MODE .or. this%mode == FULL_2D_MODE) then

#ifdef USE_MPI
       !MPI subroutines will handle the periodic boundaries.
#else
            if (this%boundaries(1) == PERIODIC_BOUNDARIES) then
                this%Hz(0, :)     = this%Hz(nx, :)
            end if

            if (this%boundaries(2) == PERIODIC_BOUNDARIES) then
                this%Hz(:, 0)     = this%Hz(:, ny)
            end if
#endif

            !$omp parallel default(shared) private(i, j, rotH, Jav_x, Jav_y, no_medium)
            !$omp do collapse(2) schedule(static)
            do j=1,ny
            do i=1,nx
                rotH         = (this%Hz(i, j) - this%Hz(i, j - 1)) * this%den_ey(j)

                no_medium = .false.
                call get_medium_polarization(this%media, this%media_map(i,j,1), &
                                            this%PDx(i,j), this%PLx(i,j,:), this%PLx_old(i,j,:), &
                                            this%Ex(i,j), this%Ex_old(i,j), rotH, no_medium)

                if (no_medium) then
                    Jav_x        = 0.5*(this%Jx_old(i,j) + (this%time - this%t_skip)*this%dJx(i,j) + &
                                       this%Jx_old(i+1,j) + (this%time - this%t_skip)*this%dJx(i+1,j))
                    this%Ex(i,j) = this%Ex(i,j) + this%dt_eps0*this%eps_x(i, j)*rotH - &
                                this%dt_eps0 * Jav_x
                end if
            end do
            end do
            !$omp end do nowait

            !$omp do collapse(2) schedule(static)
            do j = 1, ny
            do i = 1, nx
                rotH          = (this%Hz(i - 1, j) - this%Hz(i, j)) * this%den_ex(i)

                no_medium = .false.
                call get_medium_polarization(this%media, this%media_map(i,j,2), &
                                            this%PDy(i,j), this%PLy(i,j,:), this%PLy_old(i,j,:), &
                                            this%Ey(i,j), this%Ey_old(i,j), rotH, no_medium)

                if (no_medium) then
                    Jav_y         = 0.5*(this%Jy_old(i,j) + (this%time - this%t_skip)*this%dJy(i,j) + &
                                        this%Jy_old(i,j+1) + (this%time - this%t_skip)*this%dJy(i,j+1))
                    this%Ey(i, j) = this%Ey(i, j) + this%dt_eps0 * this%eps_y(i, j) * rotH -&
                                    this%dt_eps0 * Jav_y
                end if
            enddo
            enddo
            !$omp end do
            !$omp end parallel

            if (this%cpml_pos(1)) then
                do j = 1, ny
                do i = 2, npml
                    rotH = (this%Hz(i-1, j) - this%Hz(i, j)) / this%dr
                    this%psix_Eyx(i, j) = this%be(i) * this%psix_Eyx(i, j) + &
                                          this%ce(i) * rotH

                    this%Ey(i, j) = this%Ey(i, j) + this%dt_eps0 * this%psix_Eyx(i, j)
                enddo
                enddo
            end if

            if (this%cpml_pos(2)) then
                do j = 1, ny
                    ii = n_sec*npml
                    do i = nx + 1 - npml, nx - 1
                        rotH = (this%Hz(i-1, j) - this%Hz(i, j)) / this%dr
                        this%psix_Eyx(ii, j) = this%be(ii-n_base) * this%psix_Eyx(ii, j) + &
                                               this%ce(ii-n_base) * rotH

                        this%Ey(i, j) = this%Ey(i, j) + this%dt_eps0 * this%psix_Eyx(ii, j)
                        ii = ii - 1

                    enddo
                enddo
            end if

            if (this%cpml_pos(3)) then
                do i = 1, nx
                do j = 2, npml
                    rotH = (this%Hz(i, j) - this%Hz(i, j-1)) / this%dr
                    this%psiy_Exy(i, j) = this%be(j) * this%psiy_Exy(i, j) + &
                            this%ce(j) * rotH

                    this%Ex(i, j) = this%Ex(i, j) + this%dt_eps0 * this%psiy_Exy(i, j)

                enddo
                enddo  
            end if

            if (this%cpml_pos(4)) then
                do i = 1, nx
                    jj = 2*npml
                    do j = ny + 1 - npml, ny - 1
                        rotH = (this%Hz(i, j) - this%Hz(i, j-1)) / this%dr
                        this%psiy_Exy(i, jj) = this%be(jj-n_base) * this%psiy_Exy(i, jj) + &
                                this%ce(jj-n_base) * rotH

                        this%Ex(i, j) = this%Ex(i, j) + this%dt_eps0 * this%psiy_Exy(i, jj)
                        jj = jj - 1
                    end do
                end do                
            end if

        end if

    end subroutine td_propagate_E_2Dfield

!###################################################################################################    
end module mxll_2D_mod