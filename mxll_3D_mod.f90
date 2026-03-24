module mxll_3D_mod

    use constants_mod
    use mxll_base_mod
    use classical_medium_mod

    implicit none

    type, extends(TMxll) :: TMxll_3D

        type(TClassicalMedium), allocatable :: media(:)

        integer               :: nx,ny, nz
        integer               :: npml
        integer               :: n_media
        integer               :: boundaries(3)
        integer               :: n_cpml_sections
        real(dp)              :: dt_mu0
        real(dp)              :: dt_eps0
        real(dp)              :: dr
        logical               :: cpml_pos(6) !Indicates if the rank has some of the
                                               ! 6 possible CPML boundaries.
        integer     , allocatable :: media_map(:,:,:,:)
        real(dp)    , allocatable :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
        real(dp)    , allocatable :: Ex_old(:,:,:), Ey_old(:,:,:), Ez_old(:,:,:)
        real(dp)    , allocatable :: Hx(:,:,:), Hy(:,:,:), Hz(:,:,:)
        real(dp)    , allocatable :: eps_z(:,:,:), eps_y(:,:,:), eps_x(:,:,:)
        real(dp)    , allocatable :: Jx(:,:,:), Jy(:,:,:), Jz(:,:,:)
        real(dp)    , allocatable :: Jx_old(:,:,:), Jy_old(:,:,:), Jz_old(:,:,:)
        real(dp)    , allocatable :: dJx(:,:,:), dJy(:,:,:), dJz(:,:,:)
        real(dp)    , allocatable :: den_ex(:), den_ey(:), den_ez(:)
        real(dp)    , allocatable :: den_hx(:), den_hy(:), den_hz(:)
        !PML STUFF
        real(dp)    , allocatable :: psix_Ezx(:,:,:)
        real(dp)    , allocatable :: psix_Eyx(:,:,:)
        real(dp)    , allocatable :: psiy_Exy(:,:,:)
        real(dp)    , allocatable :: psiy_Ezy(:,:,:)
        real(dp)    , allocatable :: psiz_Exz(:,:,:)
        real(dp)    , allocatable :: psiz_Eyz(:,:,:)
        real(dp)    , allocatable :: psix_Hyx(:,:,:)
        real(dp)    , allocatable :: psix_Hzx(:,:,:)
        real(dp)    , allocatable :: psiy_Hxy(:,:,:)
        real(dp)    , allocatable :: psiy_Hzy(:,:,:)
        real(dp)    , allocatable :: psiz_Hxz(:,:,:)
        real(dp)    , allocatable :: psiz_Hyz(:,:,:)
        real(dp)    , allocatable :: be(:)
        real(dp)    , allocatable :: ce(:)
        real(dp)    , allocatable :: bh(:)
        real(dp)    , allocatable :: ch(:)

        !Medium stuff
        real(dp)    , allocatable :: PDx(:,:,:), PDy(:,:,:), PDz(:,:,:)
        real(dp)    , allocatable :: PLx(:,:,:,:), PLy(:,:,:,:), PLz(:,:,:,:)
        real(dp)    , allocatable :: PLx_old(:,:,:,:), PLy_old(:,:,:,:), PLz_old(:,:,:,:)

    contains
        procedure :: init => init_3Dgrid
        procedure :: kill => kill_3Dgrid
        procedure :: td_propagate_H_field => td_propagate_H_3Dfield
        procedure :: td_propagate_E_field => td_propagate_E_3Dfield 

    end type TMxll_3D

contains

!###################################################################################################

    subroutine init_3Dgrid(this, grid_Ndims, npml, boundaries, dt, dr, mode, n_media, &
                           mpi_coords, mpi_dims)
        
        class(TMxll_3D), intent(inout) :: this
        integer        , intent(in)    :: grid_Ndims(3)
        integer        , intent(in)    :: npml
        integer        , intent(in)    :: boundaries(3)
        integer        , intent(in)    :: mode
        integer        , intent(in)    :: n_media
        integer        , intent(in)    :: mpi_coords(3)
        integer        , intent(in)    :: mpi_dims(3)
        real(dp)       , intent(in)    :: dt
        real(dp)       , intent(in)    :: dr

        integer  :: nx, ny, nz
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
        nz = grid_Ndims(3)
     
        this%chunk_coor = mpi_coords
        
        this%nx            = nx
        this%ny            = ny
        this%nz            = nz
        this%npml          = npml
        this%n_media       = n_media
        this%boundaries(1) = boundaries(1)
        this%boundaries(2) = boundaries(2)
        this%boundaries(3) = boundaries(3)
        this%dr            = dr
        this%dt_eps0       = dt/eps0
        this%dt_mu0        = dt/mu0

        if (.not. allocated(this%den_ex))   allocate(this%den_ex(nx))
        if (.not. allocated(this%den_ey))   allocate(this%den_ey(ny))
        if (.not. allocated(this%den_ez))   allocate(this%den_ez(nz))
        if (.not. allocated(this%den_hx))   allocate(this%den_hx(nx))
        if (.not. allocated(this%den_hy))   allocate(this%den_hy(ny))
        if (.not. allocated(this%den_hz))   allocate(this%den_hz(nz))

        this%den_ex = 1.0d0 / dr
        this%den_ey = 1.0d0 / dr
        this%den_ez = 1.0d0 / dr
        this%den_hx = 1.0d0 / dr
        this%den_hy = 1.0d0 / dr
        this%den_hz = 1.0d0 / dr

        if (.not. allocated(this%Ex))      allocate(this%Ex(0:nx, ny+1, nz+1))
        if (.not. allocated(this%Ey))      allocate(this%Ey(nx+1, 0:ny, nz+1))
        if (.not. allocated(this%Ez))      allocate(this%Ez(nx+1, ny+1, 0:nz))
        if (.not. allocated(this%Ex_old))  allocate(this%Ex_old(0:nx, ny+1, nz+1))
        if (.not. allocated(this%Ey_old))  allocate(this%Ey_old(nx+1, 0:ny, nz+1))
        if (.not. allocated(this%Ez_old))  allocate(this%Ez_old(nx+1, ny+1, 0:nz))
        if (.not. allocated(this%Hx))      allocate(this%Hx(nx, 0:ny, 0:nz))
        if (.not. allocated(this%Hy))      allocate(this%Hy(0:nx, ny, 0:nz))
        if (.not. allocated(this%Hz))      allocate(this%Hz(0:nx, 0:ny, nz))
        if (.not. allocated(this%eps_x))   allocate(this%eps_x(nx,ny,nz))
        if (.not. allocated(this%eps_y))   allocate(this%eps_y(nx,ny,nz))
        if (.not. allocated(this%eps_z))   allocate(this%eps_z(nx,ny,nz))
        if (.not. allocated(this%Jx))      allocate(this%Jx(nx+1, ny, nz))
        if (.not. allocated(this%Jy))      allocate(this%Jy(nx, ny+1, nz))
        if (.not. allocated(this%Jz))      allocate(this%Jz(nx, ny, nz+1))
        if (.not. allocated(this%Jx_old))  allocate(this%Jx_old(nx+1, ny, nz))
        if (.not. allocated(this%Jy_old))  allocate(this%Jy_old(nx, ny+1, nz))
        if (.not. allocated(this%Jz_old))  allocate(this%Jz_old(nx, ny, nz+1))
        if (.not. allocated(this%dJx))     allocate(this%dJx(nx+1, ny, nz))
        if (.not. allocated(this%dJy))     allocate(this%dJy(nx, ny+1, nz))
        if (.not. allocated(this%dJz))     allocate(this%dJz(nx, ny, nz+1))
        
        this%Ex      = 0.0d0
        this%Ey      = 0.0d0
        this%Ez      = 0.0d0
        this%Ex_old  = 0.0d0
        this%Ey_old  = 0.0d0
        this%Ez_old  = 0.0d0
        this%Hz      = 0.0d0
        this%Hx      = 0.0d0
        this%Hy      = 0.0d0
        this%Jx      = 0.0d0 
        this%Jy      = 0.0d0  
        this%Jz      = 0.0d0
        this%Jx_old  = 0.0d0  
        this%Jy_old  = 0.0d0  
        this%Jz_old  = 0.0d0  
        this%dJx     = 0.0d0  
        this%dJy     = 0.0d0  
        this%dJz     = 0.0d0
        this%eps_x   = 1.0d0
        this%eps_y   = 1.0d0
        this%eps_z   = 1.0d0


        call read_init_media(this%media, n_media, this%eps_x, this%eps_y, this%eps_z, grid_Ndims, dr, &
                             this%media_map, this%nx, this%ny, this%nz, dt, mpi_coords, mpi_dims)


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

        if (.not. allocated(this%PDx))   allocate(this%PDx(nx,ny,nz))
        if (.not. allocated(this%PDy))   allocate(this%PDy(nx,ny,nz))
        if (.not. allocated(this%PDz))   allocate(this%PDz(nx,ny,nz))
        this%PDx = 0.0d0
        this%PDy = 0.0d0
        this%PDz = 0.0d0

        if (.not. allocated(this%PLx))       allocate(this%PLx(nx,ny,nz,n_max_poles))
        if (.not. allocated(this%PLy))       allocate(this%PLy(nx,ny,nz,n_max_poles))
        if (.not. allocated(this%PLz))       allocate(this%PLz(nx,ny,nz,n_max_poles))
        if (.not. allocated(this%PLx_old))   allocate(this%PLx_old(nx,ny,nz,n_max_poles))
        if (.not. allocated(this%PLy_old))   allocate(this%PLy_old(nx,ny,nz,n_max_poles))
        if (.not. allocated(this%PLz_old))   allocate(this%PLz_old(nx,ny,nz,n_max_poles))
        this%PLx    = 0.0d0
        this%PLy    = 0.0d0
        this%PLz    = 0.0d0
        this%PLx_old  = 0.0d0
        this%PLy_old  = 0.0d0
        this%PLz_old  = 0.0d0

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
        this%cpml_pos(5) = (mpi_coords(3) == 0) .and.&
                            this%boundaries(3) == CPML_BOUNDARIES
        this%cpml_pos(6) = (mpi_coords(3) == mpi_dims(3)-1) .and.&
                            this%boundaries(3) == CPML_BOUNDARIES
        
        n_sec                = 1
        this%n_cpml_sections = n_sec
#else
        this%cpml_pos(1) = this%boundaries(1) == CPML_BOUNDARIES
        this%cpml_pos(2) = this%boundaries(1) == CPML_BOUNDARIES
        this%cpml_pos(3) = this%boundaries(2) == CPML_BOUNDARIES
        this%cpml_pos(4) = this%boundaries(2) == CPML_BOUNDARIES
        this%cpml_pos(5) = this%boundaries(3) == CPML_BOUNDARIES
        this%cpml_pos(6) = this%boundaries(3) == CPML_BOUNDARIES
    
        n_sec                = 2
        this%n_cpml_sections = n_sec
#endif
        
        if (this%cpml_pos(1) .or. this%cpml_pos(2)) then
            
            if (.not. allocated(this%psix_Eyx)) allocate(this%psix_Eyx(n_sec*npml, ny, nz))
            if (.not. allocated(this%psix_Ezx)) allocate(this%psix_Ezx(n_sec*npml, ny, nz))
            if (.not. allocated(this%psix_Hyx)) allocate(this%psix_Hyx(n_sec*(npml-1), ny, nz))
            if (.not. allocated(this%psix_Hzx)) allocate(this%psix_Hzx(n_sec*(npml-1), ny, nz))
           
            this%psix_Eyx = 0.0d0
            this%psix_Hzx = 0.0d0
            this%psix_Ezx = 0.0d0
            this%psix_Hyx = 0.0d0

        end if
        
        if (this%cpml_pos(3) .or. this%cpml_pos(4)) then      
            
            if (.not. allocated(this%psiy_Exy)) allocate(this%psiy_Exy(nx, n_sec*npml, nz))
            if (.not. allocated(this%psiy_Ezy)) allocate(this%psiy_Ezy(nx, n_sec*npml, nz))
            if (.not. allocated(this%psiy_Hxy)) allocate(this%psiy_Hxy(nx, n_sec*(npml-1), nz))
            if (.not. allocated(this%psiy_Hzy)) allocate(this%psiy_Hzy(nx, n_sec*(npml-1), nz))
            
            this%psiy_Exy = 0.0d0
            this%psiy_Hzy = 0.0d0
            this%psiy_Ezy = 0.0d0
            this%psiy_Hxy = 0.0d0
            
        end if

        if (this%cpml_pos(5) .or. this%cpml_pos(6)) then      
            
            if (.not. allocated(this%psiz_Exz)) allocate(this%psiz_Exz(nx, ny, n_sec*npml))
            if (.not. allocated(this%psiz_Eyz)) allocate(this%psiz_Eyz(nx, ny, n_sec*npml))
            if (.not. allocated(this%psiz_Hxz)) allocate(this%psiz_Hxz(nx, ny, n_sec*(npml-1)))
            if (.not. allocated(this%psiz_Hyz)) allocate(this%psiz_Hyz(nx, ny, n_sec*(npml-1)))
            
            this%psiz_Exz = 0.0d0
            this%psiz_Eyz = 0.0d0
            this%psiz_Hxz = 0.0d0
            this%psiz_Hyz = 0.0d0
            
        end if
        
        rank_with_cpml = this%cpml_pos(1) .or. this%cpml_pos(2) .or. &
                         this%cpml_pos(3) .or. this%cpml_pos(4) .or. &
                         this%cpml_pos(5) .or. this%cpml_pos(6)

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

            if(this%cpml_pos(5)) then

                do k = 1, nz
                    if (k<=npml) then
                        this%den_ez(k) = 1.0 / (kappae(k) * dr)
                    endif
                end do
                
                do k = 1, nz
                    if (k<=(npml - 1)) then
                        this%den_hz(k) = 1.0 / (kappah(k) * dr)
                    endif
                end do

            end if

            if (this%cpml_pos(6)) then
            
                kk = npml
                do k = 1, nz-1
                    if (k>=(nz + 1 - npml)) then
                        this%den_ez(k) = 1.0 / (kappae(kk) * dr)
                        kk = kk - 1
                    endif
                enddo

                kk = npml - 1
                do k = 1, nz-1
                    if (k>=(nz + 1 - npml)) then
                        this%den_hz(k) = 1.0 / (kappah(kk) * dr)
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

    end subroutine init_3Dgrid

!###################################################################################################

    subroutine kill_3Dgrid(this)

        class(TMxll_3D), intent(inout) :: this
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
        if (allocated(this%psiy_Exy)) deallocate(this%psiy_Exy)
        if (allocated(this%psiy_Ezy)) deallocate(this%psiy_Ezy)
        if (allocated(this%psiz_Exz)) deallocate(this%psiz_Exz)
        if (allocated(this%psiz_Eyz)) deallocate(this%psiz_Eyz)
        if (allocated(this%psix_Hyx)) deallocate(this%psix_Hyx)
        if (allocated(this%psix_Hzx)) deallocate(this%psix_Hzx)
        if (allocated(this%psiy_Hxy)) deallocate(this%psiy_Hxy)
        if (allocated(this%psiy_Hzy)) deallocate(this%psiy_Hzy)
        if (allocated(this%psiz_Hxz)) deallocate(this%psiz_Hxz)
        if (allocated(this%psiz_Hyz)) deallocate(this%psiz_Hyz)
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
        if (allocated(this%den_ez))   deallocate(this%den_ez)
        if (allocated(this%den_hx))   deallocate(this%den_hx)
        if (allocated(this%den_hy))   deallocate(this%den_hy)     
        if (allocated(this%den_hz))   deallocate(this%den_hz)
        if (allocated(this%ce))       deallocate(this%ce)
        if (allocated(this%be))       deallocate(this%be)
        if (allocated(this%ch))       deallocate(this%ch)
        if (allocated(this%bh))       deallocate(this%bh)
        if (allocated(this%media_map))deallocate(this%media_map)
        if (allocated(this%PLx_old))  deallocate(this%PLx_old)
        if (allocated(this%PLy_old))  deallocate(this%PLy_old)
        if (allocated(this%PLz_old))  deallocate(this%PLz_old)
        if (allocated(this%PDx))      deallocate(this%PDx)
        if (allocated(this%PDy))      deallocate(this%PDy)
        if (allocated(this%PDz))      deallocate(this%PDz)
        if (allocated(this%PLx))      deallocate(this%PLx)
        if (allocated(this%PLy))      deallocate(this%PLy)
        if (allocated(this%PLz))      deallocate(this%PLz)


    end subroutine kill_3Dgrid   

!###################################################################################################
    subroutine td_propagate_H_3Dfield(this)

        class(TMxll_3D), intent(inout) :: this

        integer  :: nx, ny, nz
        integer  :: npml
        integer  :: n_base
        integer  :: n_sec
        integer  :: i,j,k, ii, jj, kk
        real(dp) :: rotE

        nx    = this%nx
        ny    = this%ny
        nz    = this%nz
        npml  = this%npml
        n_sec = this%n_cpml_sections

        if (n_sec==1) n_base = 0
        if (n_sec==2) n_base = npml-1

#ifdef USE_MPI
        !MPI subroutines will handle the periodic boundaries.
#else
        if (this%boundaries(1) == PERIODIC_BOUNDARIES) then
            this%Ez(nx+1,:,:) = this%Ez(1,:,:)
            this%Ey(nx+1,:,:) = this%Ey(1,:,:)
        end if

        if (this%boundaries(2) == PERIODIC_BOUNDARIES) then
            this%Ex(:,ny+1,:) = this%Ex(:,1,:)
            this%Ez(:,ny+1,:) = this%Ez(:,1,:)
        end if

        if (this%boundaries(3) == PERIODIC_BOUNDARIES) then
            this%Ex(:,:,nz+1) = this%Ex(:,:,1)
            this%Ey(:,:,nz+1) = this%Ey(:,:,1)
        end if

#endif


        !$omp parallel default(shared) private(i, j, k, rotE)
        !$omp do collapse(3) schedule(static)
        do i = 1, nx
        do j = 1, ny
        do k = 1, nz
            rotE             = (this%Ez(i, j, k) - this%Ez(i, j + 1, k)) * this%den_hy(j) + &
            (this%Ey(i, j, k + 1) - this%Ey(i, j, k)) * this%den_hz(k)
            this%Hx(i, j, k) = this%Hx(i, j, k) + this%dt_mu0 * rotE
        end do
        end do
        end do
        !$omp end do nowait

        !$omp do collapse(3) schedule(static)
        do i = 1, nx
        do j = 1, ny
        do k = 1, nz
            rotE             = (this%Ez(i + 1, j, k) - this%Ez(i, j, k)) * this%den_hx(i) + &
            (this%Ex(i, j, k) - this%Ex(i, j, k + 1)) * this%den_hz(k)
            this%Hy(i, j, k) = this%Hy(i, j, k) + this%dt_mu0 * rotE
        end do
        end do
        end do
        !$omp end do nowait

        !$omp do collapse(3) schedule(static)
        do i = 1, nx
        do j = 1, ny
        do k = 1, nz
            rotE             = (this%Ey(i, j, k)-this%Ey(i + 1, j, k))*this%den_hx(i) + &
                               (this%Ex(i, j + 1, k)-this%Ex(i, j, k))*this%den_hy(j)
            this%Hz(i, j, k) = this%Hz(i, j, k) + this%dt_mu0 * rotE
        end do
        end do
        end do
        !$omp end do
        !$omp end parallel
        

        if (this%cpml_pos(1)) then
            do j = 1, ny
            do k = 1, nz
            do i = 1, npml - 1
                this%psix_Hyx(i, j, k) = this%bh(i)*this%psix_Hyx(i, j, k) + &
                                       this%ch(i)*(this%Ez(i+1, j, k)-this%Ez(i, j, k))/this%dr
                this%psix_Hzx(i, j, k) = this%bh(i)*this%psix_Hzx(i, j, k) + &
                                       this%ch(i)*(this%Ey(i, j, k)-this%Ey(i+1, j, k))/this%dr

                this%Hy(i, j, k) = this%Hy(i, j, k) + this%dt_mu0*this%psix_Hyx(i, j, k)
                this%Hz(i, j, k) = this%Hz(i, j, k) + this%dt_mu0*this%psix_Hzx(i, j, k)
            
            end do 
            end do
            end do
        end if

        if (this%cpml_pos(2)) then
            do j = 1, ny
            do k = 1, nz
                ii = n_sec*(npml - 1)
                do i = nx + 1 - npml, nx - 1
                    this%psix_Hyx(ii, j, k) = this%bh(ii-n_base)*this%psix_Hyx(ii, j, k) + &
                                        this%ch(ii-n_base)*(this%Ez(i+1, j, k) - this%Ez(i, j, k))/this%dr
                    this%psix_Hzx(ii, j, k) = this%bh(ii-n_base)*this%psix_Hzx(ii, j, k) + &
                                        this%ch(ii-n_base)*(this%Ey(i, j, k) - this%Ey(i+1, j, k))/this%dr 

                    this%Hy(i, j, k) = this%Hy(i, j, k) + this%dt_mu0*this%psix_Hyx(ii, j, k)
                    this%Hz(i, j, k) = this%Hz(i, j, k) + this%dt_mu0*this%psix_Hzx(ii, j, k)

                    ii = ii - 1

                end do
            end do
            end do
        end if

        if (this%cpml_pos(3)) then
            do k = 1, nz
            do i = 1, nx
            do j = 1, npml - 1
                this%psiy_Hxy(i, j, k) = this%bh(j)*this%psiy_Hxy(i, j, k) + &
                                this%ch(j)*(this%Ez(i, j, k) - this%Ez(i, j+1, k))/this%dr
                this%psiy_Hzy(i, j, k) = this%bh(j) * this%psiy_Hzy(i, j, k) + &
                                this%ch(j)*(this%Ex(i, j+1, k) - this%Ex(i, j, k))/this%dr

                this%Hx(i, j, k) = this%Hx(i, j, k) + this%dt_mu0 * this%psiy_Hxy(i, j, k)
                this%Hz(i, j, k) = this%Hz(i, j, k) + this%dt_mu0 * this%psiy_Hzy(i, j, k)

            end do
            end do
            end do
        end if

        if (this%cpml_pos(4)) then
            do k = 1, nz
            do i = 1, nx
                jj = n_sec*(npml-1)
                do j = ny + 1 - npml, ny - 1
                    
                    this%psiy_Hxy(i, jj, k) = this%bh(jj-n_base)*this%psiy_Hxy(i, jj, k) + &
                                        this%ch(jj-n_base)*(this%Ez(i, j, k) - this%Ez(i, j+1, k)) / this%dr
                    this%psiy_Hzy(i, jj, k) = this%bh(jj-n_base) * this%psiy_Hzy(i, jj, k) + &
                                        this%ch(jj-n_base)*(this%Ex(i, j+1, k) - this%Ex(i, j, k)) / this%dr

                    this%Hx(i, j, k) = this%Hx(i, j, k) + this%dt_mu0*this%psiy_Hxy(i, jj, k)
                    this%Hz(i, j, k) = this%Hz(i, j, k) + this%dt_mu0*this%psiy_Hzy(i, jj, k)

                    jj = jj - 1
                end do
            end do
            end do  
        end if

        if (this%cpml_pos(5)) then
            do i = 1, nx
            do j = 1, ny
            do k = 1, npml-1
                this%psiz_Hxz(i, j, k) = this%bh(k)*this%psiz_Hxz(i, j, k) + &
                                this%ch(k)*(this%Ey(i, j, k+1) - this%Ey(i, j, k))/this%dr
                this%psiz_Hyz(i, j, k) = this%bh(k) * this%psiz_Hyz(i, j, k) + &
                                this%ch(k)*(this%Ex(i, j, k) - this%Ex(i, j, k+1))/this%dr

                this%Hx(i, j, k) = this%Hx(i, j, k) + this%dt_mu0 * this%psiz_Hxz(i, j, k)
                this%Hy(i, j, k) = this%Hy(i, j, k) + this%dt_mu0 * this%psiz_Hyz(i, j, k)

            end do
            end do
            end do
        end if

        if (this%cpml_pos(6)) then
            do i = 1, nx
            do j = 1, ny
                kk = n_sec*(npml-1)
                do k = nz + 1 - npml, nz - 1
                    
                    this%psiz_Hxz(i, j, kk) = this%bh(kk-n_base)*this%psiz_Hxz(i, j, kk) + &
                                        this%ch(kk-n_base)*(this%Ey(i, j, k+1) - this%Ey(i, j, k)) / this%dr
                    this%psiz_Hyz(i, j, kk) = this%bh(kk-n_base)*this%psiz_Hyz(i, j, kk) + &
                                        this%ch(kk-n_base)*(this%Ex(i, j, k) - this%Ex(i, j, k+1)) / this%dr

                    this%Hx(i, j, k) = this%Hx(i, j, k) + this%dt_mu0*this%psiz_Hxz(i, j, kk)
                    this%Hy(i, j, k) = this%Hy(i, j, k) + this%dt_mu0*this%psiz_Hyz(i, j, kk)

                    kk = kk - 1
                end do
            end do
            end do  
        end if

    end subroutine td_propagate_H_3Dfield

!###################################################################################################

    subroutine td_propagate_E_3Dfield(this, t)

        class(TMxll_3D), intent(inout) :: this
        integer, intent(in)            :: t

        integer  :: nx, ny, nz
        integer  :: npml
        integer  :: n_base
        integer  :: n_sec
        integer  :: i,j,k, ii, jj, kk
        logical  :: no_medium
        real(dp) :: rotH, Jav_x, Jav_y, Jav_z

        nx    = this%nx
        ny    = this%ny
        nz    = this%nz
        npml  = this%npml
        n_sec = this%n_cpml_sections

        if (n_sec==1) n_base = 0
        if (n_sec==2) n_base = npml

#ifdef USE_MPI
        !MPI subroutines will handle the periodic boundaries.
#else
        if (this%boundaries(1) == PERIODIC_BOUNDARIES) then
            this%Hy(0,:,:) = this%Hx(nx,:,:)
            this%Hz(0,:,:) = this%Hz(nx,:,:)
        end if

        if (this%boundaries(2) == PERIODIC_BOUNDARIES) then
            this%Hx(:,0,:) = this%Hx(:,ny,:)
            this%Hz(:,0,:) = this%Hz(:,ny,:)
        end if

        if (this%boundaries(3) == PERIODIC_BOUNDARIES) then
            this%Hx(:,:,0) = this%Hx(:,:,nz)
            this%Hy(:,:,0) = this%Hy(:,:,nz)
        end if

#endif

        !$omp parallel default(shared) private(i, j, k, rotH, Jav_x, Jav_y, Jav_z, no_medium)

        !$omp do collapse(3) schedule(static)
        do k=1,nz
        do j=1,ny
        do i=1,nx
            rotH  = (this%Hz(i, j, k) - this%Hz(i, j - 1, k)) * this%den_ey(j) + &
                    (this%Hy(i, j, k-1) - this%Hy(i, j, k)) * this%den_ez(k)

            no_medium = .false.
            call get_medium_polarization(this%media, this%media_map(i,j,k,1), &
                                            this%PDx(i,j,k), this%PLx(i,j,k,:), this%PLx_old(i,j,k,:), &
                                            this%Ex(i,j,k), this%Ex_old(i,j,k), rotH, no_medium)

            if (no_medium) then
                Jav_x = 0.5*(this%Jx_old(i,j,k) + (this%time - this%t_skip)*this%dJx(i,j,k) + &
                             this%Jx_old(i+1,j,k) + (this%time - this%t_skip)*this%dJx(i+1,j,k))
                this%Ex(i,j, k) = this%Ex(i,j, k) + this%dt_eps0*this%eps_x(i, j, k)*rotH - &
                                this%dt_eps0 * Jav_x
            end if
        end do
        end do
        end do
        !$omp end do nowait

        !$omp do collapse(3) schedule(static)
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
            rotH  = (this%Hz(i - 1, j, k) - this%Hz(i, j, k)) * this%den_ex(i) + &
                    (this%Hx(i, j, k) - this%Hx(i, j, k - 1)) * this%den_ez(k)

            no_medium = .false.
            call get_medium_polarization(this%media, this%media_map(i,j,k,2), &
                                            this%PDy(i,j,k), this%PLy(i,j,k,:), this%PLy_old(i,j,k,:), &
                                            this%Ey(i,j,k), this%Ey_old(i,j,k), rotH, no_medium)

            if (no_medium) then
                Jav_y = 0.5*(this%Jy_old(i,j,k) + (this%time - this%t_skip)*this%dJy(i,j,k) + &
                             this%Jy_old(i,j+1,k) + (this%time - this%t_skip)*this%dJy(i,j+1,k))
                this%Ey(i, j, k) = this%Ey(i, j, k) + this%dt_eps0*this%eps_y(i, j, k)*rotH - &
                                this%dt_eps0 * Jav_y
            end if

        end do
        end do
        end do
        !$omp end do nowait

        !$omp do collapse(3) schedule(static)
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx

            rotH  = (this%Hy(i, j, k) - this%Hy(i - 1, j, k)) * this%den_ex(i) +  &
                    (this%Hx(i, j - 1, k) - this%Hx(i, j, k)) * this%den_ey(j)

            no_medium = .false.
            call get_medium_polarization(this%media, this%media_map(i,j,k,3), &
                                            this%PDz(i,j,k), this%PLz(i,j,k,:), this%PLz_old(i,j,k,:), &
                                            this%Ez(i,j,k), this%Ez_old(i,j,k), rotH, no_medium)

            if (no_medium) then
                Jav_z = 0.5*(this%Jz_old(i,j,k) + (this%time - this%t_skip)*this%dJz(i,j,k) + &
                             this%Jz_old(i,j,k+1) + (this%time - this%t_skip)*this%dJz(i,j,k+1))
                this%Ez(i, j, k) = this%Ez(i, j, k) + this%dt_eps0 * this%eps_z(i, j, k) * rotH - &
                                this%dt_eps0 * Jav_z
            end if
        end do
        end do
        end do
        !$omp end do

        !$omp end parallel

        if (this%cpml_pos(1)) then
            do j = 1, ny
            do k = 1, nz
            do i = 2, npml
                this%psix_Eyx(i, j, k) = this%be(i) * this%psix_Eyx(i, j, k) + &
                                    this%ce(i) * (this%Hz(i-1, j, k) - this%Hz(i, j, k)) / this%dr
                this%psix_Ezx(i, j, k) = this%be(i) * this%psix_Ezx(i, j, k) + &
                                      this%ce(i) * (this%Hy(i, j, k) - this%Hy(i-1, j, k)) / this%dr 

                this%Ey(i, j, k) = this%Ey(i, j, k) + this%dt_eps0 * this%psix_Eyx(i, j, k)                
                this%Ez(i, j, k) = this%Ez(i, j, k) + this%dt_eps0 * this%psix_Ezx(i, j, k)
            end do
            end do
            end do
        end if

        if (this%cpml_pos(2)) then
            do j = 1, ny
            do k = 1, nz
                ii = n_sec*npml
                do i = nx + 1 - npml, nx - 1
                    this%psix_Eyx(ii, j, k) = this%be(ii-n_base) * this%psix_Eyx(ii, j, k) + &
                                            this%ce(ii-n_base) * (this%Hz(i-1, j, k) - this%Hz(i, j, k)) / this%dr
                    this%psix_Ezx(ii, j, k) = this%be(ii-n_base) * this%psix_Ezx(ii, j, k) + &
                                            this%ce(ii-n_base) * (this%Hy(i, j, k) - this%Hy(i-1, j, k)) / this%dr

                    this%Ey(i, j, k) = this%Ey(i, j, k) + this%dt_eps0 * this%psix_Eyx(ii, j, k)
                    this%Ez(i, j, k) = this%Ez(i, j, k) + this%dt_eps0 * this%psix_Ezx(ii, j, k)
                    ii = ii - 1
                end do
            end do
            end do
        end if

        if (this%cpml_pos(3)) then
            do k = 1, nz
            do i = 1, nx
            do j = 2, npml
                this%psiy_Exy(i, j, k) = this%be(j) * this%psiy_Exy(i, j, k) + &
                        this%ce(j) * (this%Hz(i, j, k) - this%Hz(i, j-1, k)) / this%dr
                this%psiy_Ezy(i, j, k) = this%be(j) * this%psiy_Ezy(i, j, k) + &
                                        this%ce(j) * (this%Hx(i, j-1, k) - this%Hx(i, j, k)) / this%dr


                this%Ex(i, j, k) = this%Ex(i, j, k) + this%dt_eps0 * this%psiy_Exy(i, j, k)                       
                this%Ez(i, j, k) = this%Ez(i, j, k) + this%dt_eps0 * this%psiy_Ezy(i, j, k)
            end do
            end do
            end do
        end if

        if (this%cpml_pos(4)) then
            do k = 1, nz
            do i = 1, nx
                jj = n_sec*npml
                do j = ny + 1 - npml, ny - 1
                    this%psiy_Exy(i, jj, k) = this%be(jj-n_base) * this%psiy_Exy(i, jj, k) + &
                            this%ce(jj-n_base) * (this%Hz(i, j, k) - this%Hz(i, j-1, k)) / this%dr
                    this%psiy_Ezy(i, jj, k) = this%be(jj-n_base) * this%psiy_Ezy(i, jj, k) + &
                            this%ce(jj-n_base) * (this%Hx(i, j-1, k) - this%Hx(i, j, k)) / this%dr

                    this%Ex(i, j, k) = this%Ex(i, j, k) + this%dt_eps0 * this%psiy_Exy(i, jj, k)       
                    this%Ez(i, j, k) = this%Ez(i, j, k) + this%dt_eps0 * this%psiy_Ezy(i, jj, k)
                    jj = jj - 1

                end do
            end do
            end do    
        end if

        if (this%cpml_pos(5)) then
            do i = 1, nx
            do j = 1, ny
            do k = 2, npml
                this%psiz_Exz(i, j, k) = this%be(k) * this%psiz_Exz(i, j, k) + &
                                    this%ce(k) * (this%Hy(i, j, k-1) - this%Hy(i, j, k)) / this%dr
                this%psiz_Eyz(i, j, k) = this%be(k) * this%psiz_Eyz(i, j, k) + &
                                    this%ce(k) * (this%Hx(i, j, k) - this%Hx(i, j, k-1)) / this%dr


                this%Ex(i, j, k) = this%Ex(i, j, k) + this%dt_eps0 * this%psiz_Exz(i, j, k)                       
                this%Ey(i, j, k) = this%Ey(i, j, k) + this%dt_eps0 * this%psiz_Eyz(i, j, k)
            end do
            end do
            end do
        end if

        if (this%cpml_pos(6)) then
            do i = 1, nx
            do j = 1, ny
                kk = n_sec*npml
                do k = nz + 1 - npml, nz - 1
                    this%psiz_Exz(i, j, kk) = this%be(kk-n_base) * this%psiz_Exz(i, j, kk) + &
                            this%ce(kk-n_base) * (this%Hy(i, j, k-1) - this%Hy(i, j, k)) / this%dr
                    this%psiz_Eyz(i, j, kk) = this%be(kk-n_base) * this%psiz_Eyz(i, j, kk) + &
                            this%ce(kk-n_base) * (this%Hx(i, j, k) - this%Hx(i, j, k-1)) / this%dr

                    this%Ex(i, j, k) = this%Ex(i, j, k) + this%dt_eps0 * this%psiz_Exz(i, j, kk)       
                    this%Ey(i, j, k) = this%Ey(i, j, k) + this%dt_eps0 * this%psiz_Eyz(i, j, kk)
                    kk = kk - 1

                end do
            end do
            end do    
        end if

        if (this%chunk_coor(1) == 0 .and. this%chunk_coor(2) == 0 .and. this%chunk_coor(3) == 0) then
            write(700,*) this%Ez(30,70,70)
        end if

    end subroutine td_propagate_E_3Dfield

!###################################################################################################    
end module mxll_3D_mod