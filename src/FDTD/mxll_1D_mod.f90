module mxll_1D_mod

    use constants_mod
    use mxll_base_mod
    use classical_medium_mod

    implicit none

    type, extends(TMxll) :: TMxll_1D

        type(TClassicalMedium), allocatable :: media(:)

        integer               :: nz
        integer               :: boundaries

        integer      , allocatable :: media_map(:)
        real(dp)     , allocatable :: Ex(:)
        real(dp)     , allocatable :: Ex_old(:)
        real(dp)     , allocatable :: Hy(:)
        real(dp)     , allocatable :: rotH_x(:)
        real(dp)     , allocatable :: eps_x(:)
        real(dp)     , allocatable :: Jx(:)
        real(dp)     , allocatable :: Jx_old(:)
        real(dp)     , allocatable :: dJx(:)
        real(dp)     , allocatable :: den_ez(:)
        real(dp)     , allocatable :: den_hz(:)
        
        !CPML variables
        real(dp), allocatable :: psiz_Exz(:)
        real(dp), allocatable :: psiz_Hyz(:)
        real(dp), allocatable :: be(:)
        real(dp), allocatable :: ce(:)
        real(dp), allocatable :: bh(:)
        real(dp), allocatable :: ch(:)

        !Medium stuff
        real(dp), allocatable :: PDx(:)
        real(dp), allocatable :: PLx(:,:)
        real(dp), allocatable :: PLx_old(:,:)

        contains
            procedure :: init => init_1Dgrid
            procedure :: kill => kill_1Dgrid
            procedure :: td_propagate_H_field => td_propagate_H_1Dfield
            procedure :: td_propagate_E_field => td_propagate_E_1Dfield

    end type TMxll_1D

contains

!###################################################################################################

    subroutine init_1Dgrid(this, grid_Ndims, npml, boundaries, dt, dr, mode, n_media,&
                         mpi_coords, mpi_dims)

        class(TMxll_1D), intent(inout) :: this
        integer        , intent(in)    :: grid_Ndims(3)
        integer        , intent(in)    :: npml
        integer        , intent(in)    :: boundaries(3)
        integer        , intent(in)    :: mode
        integer        , intent(in)    :: n_media
        integer        , intent(in)    :: mpi_coords(3)
        integer        , intent(in)    :: mpi_dims(3)
        real(dp)       , intent(in)    :: dt
        real(dp)       , intent(in)    :: dr

        real(dp), allocatable :: alphae(:), sige(:), kappae(:)
        real(dp), allocatable :: alphah(:), sigh(:), kappah(:)

        integer  :: nz
        integer  :: kk, k
        integer  :: n_max_poles
        real(dp) :: sig_to_au
        real(dp) :: sigmaCPML

#ifdef USE_MPI
        if (mpi_coords(1)/=0) return
#endif
        sig_to_au =  C_to_au**2*sec_to_au/(kg_to_au*m_to_au**3)

        nz              = grid_Ndims(1)
        this%nz         = nz
        this%n_media    = n_media
        this%boundaries = boundaries(1)
        this%dr         = dr
        this%dt_eps0    = dt/eps0
        this%dt_mu0     = dt/mu0
        this%dt         = dt
        this%npml       = npml

        if (.not. allocated(this%Ex))       allocate(this%Ex(nz))
        if (.not. allocated(this%Ex_old))   allocate(this%Ex_old(nz))       
        if (.not. allocated(this%Hy))       allocate(this%Hy(nz-1))
        if (.not. allocated(this%Jx))       allocate(this%Jx(nz))
        if (.not. allocated(this%eps_x))    allocate(this%eps_x(nz))
        if (.not. allocated(this%Jx_old))   allocate(this%Jx_old(nz))
        if (.not. allocated(this%dJx))      allocate(this%dJx(nz))
        if (.not. allocated(this%den_ez))   allocate(this%den_ez(nz))
        if (.not. allocated(this%den_hz))   allocate(this%den_hz(nz))
        
        this%Ex        = M_ZERO
        this%Ex_old    = M_ZERO      
        this%Hy        = M_ZERO
        this%Jx        = M_ZERO
        this%Jx_old    = M_ZERO
        this%dJx       = M_ZERO
        this%den_ez    = M_ONE/dr
        this%den_hz    = M_ONE/dr
        this%eps_x     = M_ONE

        call read_init_media(this%media, this%n_media, this%eps_x, grid_Ndims, dr, &
                             this%media_map, this%nz, dt)

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

        if (.not. allocated(this%PDx)) allocate(this%PDx(nz))
        this%PDx   = M_ZERO

        if (.not. allocated(this%PLx))     allocate(this%PLx(nz, n_max_poles))
        if (.not. allocated(this%PLx_old)) allocate(this%PLx_old(nz, n_max_poles))
        this%PLx     = M_ZERO
        this%PLx_old   = M_ZERO

        if (this%boundaries == CPML_BOUNDARIES) then
            if (.not. allocated(this%psiz_Exz))  allocate(this%psiz_Exz(2*npml))
            if (.not. allocated(this%psiz_Hyz))  allocate(this%psiz_Hyz(2*(npml-1)))
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
            
            this%psiz_Exz  = 0.0d0
            this%psiz_Hyz  = 0.0d0

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
            
            kk=npml
            do k=1,nz-1
                if(k<=npml)then
                    this%den_ez(k)=1.0/(kappae(k)*dr)
                else if(k>=(nz+1-npml))then
                    this%den_ez(k)=1.0/(kappae(kk)*dr)
                    kk=kk-1
                else
                    this%den_ez(k)=1.0/dr
                end if
            end do

            kk=npml-1
            do k=1, nz-1
                if(k<=(npml-1))then
                    this%den_hz(k)=1.0/(kappah(k)*dr)
                elseif(k>=(nz+1-npml))then
                    this%den_hz(k)=1.0/(kappah(kk)*dr)
                    kk=kk-1
                else
                    this%den_hz(k)=1.0/dr
                endif
            enddo

        end if


        if (allocated(alphae)) deallocate(alphae)
        if (allocated(alphah)) deallocate(alphah)
        if (allocated(sige))   deallocate(sige)
        if (allocated(sigh))   deallocate(sigh)
        if (allocated(kappae)) deallocate(kappae)
        if (allocated(kappah)) deallocate(kappah)

    end subroutine init_1Dgrid

!###################################################################################################

    subroutine kill_1Dgrid(this)

        class(TMxll_1D), intent(inout) :: this
        integer :: i

        if (allocated(this%media)) then
            do i=1,size(this%media)
                call this%media(i)%kill_classical_medium()
            end do
            deallocate(this%media)
        end if

        if (allocated(this%Ex))        deallocate(this%Ex)
        if (allocated(this%Ex_old))    deallocate(this%Ex_old)       
        if (allocated(this%Hy))        deallocate(this%Hy)
        if (allocated(this%Jx))        deallocate(this%Jx)
        if (allocated(this%eps_x))     deallocate(this%eps_x)
        if (allocated(this%Jx_old))    deallocate(this%Jx_old)
        if (allocated(this%dJx))       deallocate(this%dJx)
        if (allocated(this%psiz_Exz))  deallocate(this%psiz_Exz)
        if (allocated(this%psiz_Hyz))  deallocate(this%psiz_Hyz)
        if (allocated(this%den_ez))    deallocate(this%den_ez)
        if (allocated(this%den_hz))    deallocate(this%den_hz)
        if (allocated(this%ce))        deallocate(this%ce)
        if (allocated(this%be))        deallocate(this%be)
        if (allocated(this%ch))        deallocate(this%ch)
        if (allocated(this%bh))        deallocate(this%bh)
        if (allocated(this%media_map)) deallocate(this%media_map)
        if (allocated(this%PDx))       deallocate(this%PDx)
        if (allocated(this%PLx))       deallocate(this%PLx)
        if (allocated(this%PLx_old))   deallocate(this%PLx_old)

    end subroutine kill_1Dgrid

!###################################################################################################

    subroutine td_propagate_H_1Dfield(this)

        class(TMxll_1D), intent(inout) :: this

        integer :: nz
        integer :: npml
        integer :: i, kk

#ifdef USE_MPI
        !This checks that only myrank==0 propagates the field when using MPI.
        !TO-DO: consider a more elgant way to do this.
        if (.not. allocated(this%Ex)) return
#endif

        nz   = this%nz
        npml = this%npml 

        if (this%boundaries == PERIODIC_BOUNDARIES) then
            this%Ex(nz+1) = this%Ex(1)
        end if

        !~~~~~~ Hy ~~~~~~~~!
        do i=1, nz-1
            this%Hy(i) = this%Hy(i) + this%dt_mu0 * (this%Ex(i)-this%Ex(i+1))*this%den_hz(i)
        enddo
        
        if (this%boundaries == CPML_BOUNDARIES) then
            !  PML for the left side Hy
            do i=1,npml-1
                this%psiz_Hyz(i) = this%bh(i)*this%psiz_Hyz(i)+ &
                                    this%ch(i)*(this%Ex(i)-this%Ex(i+1))/this%dr
                this%Hy(i) = this%Hy(i)+this%dt_mu0*this%psiz_Hyz(i)
            enddo
            !  PML for the right side Hy
            kk=2*(npml-1)
            do i=nz+1-npml,nz-1
                this%psiz_Hyz(kk) = this%bh(kk-(npml-1))*this%psiz_Hyz(kk)+&
                                this%ch(kk-(npml-1))*(this%Ex(i)-this%Ex(i+1))/this%dr
                this%Hy(i) = this%Hy(i)+this%dt_mu0*this%psiz_Hyz(kk)
                kk=kk-1
            enddo
        end if

    end subroutine td_propagate_H_1Dfield

!###################################################################################################

    subroutine td_propagate_E_1Dfield(this, t)

        class(TMxll_1D), intent(inout) :: this
        integer     , intent(in)    :: t

        integer  :: nz
        integer  :: npml
        integer  :: i, kk
        real(dp) :: rotH
        real(dp) :: Jx
        logical  :: no_medium

#ifdef USE_MPI
        !This checks that only myrank==0 propagates the field when using MPI.
        !TO-DO: consider a more elgant way to do this.
        if (.not. allocated(this%Ex)) return
#endif

        nz   = this%nz
        npml = this%npml

        do i=2,nz-1
            rotH       = (this%Hy(i-1)-this%Hy(i))*this%den_ez(i)

            no_medium = .false.

            call get_medium_polarization(this%media, this%media_map(i), &
                                            this%PDx(i), this%PLx(i,:), this%PLx_old(i,:), &
                                            this%Ex(i), this%Ex_old(i), rotH, no_medium)
            if (no_medium) then

                Jx = this%Jx_old(i) + (this%time - this%t_skip) * this%dJx(i)

                this%Ex(i) = this%Ex(i) + this%dt_eps0*this%eps_x(i)*rotH - this%dt_eps0*Jx
            end if

        end do             

        !  PML for the left side Ex
        if (this%boundaries == CPML_BOUNDARIES) then
            do i=2,npml
                rotH             = (this%Hy(i-1)-this%Hy(i))/this%dr
                this%psiz_Exz(i) = this%be(i)*this%psiz_Exz(i)+this%ce(i)*rotH
                this%Ex(i)       = this%Ex(i) + this%dt_eps0*this%psiz_Exz(i)
            enddo
            !  PML for right side Ex
            kk=2*npml
            do i=nz+1-npml,nz-1
                rotH              = (this%Hy(i-1)-this%Hy(i))/this%dr
                this%psiz_Exz(kk) = this%be(kk-npml)*this%psiz_Exz(kk)+this%ce(kk-npml)*rotH
                this%Ex(i)       = this%Ex(i)+this%dt_eps0*this%psiz_Exz(kk)
                kk=kk-1
            enddo
        end if


    end subroutine td_propagate_E_1Dfield

!###################################################################################################

end module mxll_1D_mod