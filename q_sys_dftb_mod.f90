module q_sys_dftb_mod

    use constants_mod
    use, intrinsic :: iso_fortran_env, only : output_unit
    use dftbp_common_constants, only : AA__Bohr, V_m__au, eV__Hartree, fs__au, imag
    use dftbplus
    ! Only needed for the internal test system
    !    use testhelpers, only : writeAutotestTag
    use q_sys_base_mod

    implicit none

    type, extends(TQ_sys_base) :: TQ_sys_dftb

        character(len=20)    :: temp_file
        integer              :: dyn_type
        integer              :: n_atoms
        integer              :: n_at_typ
        integer              :: euler_steps
        integer              :: std_unit=100
        real(dp)             :: Kinetic_Energy 
        type(TDftbPlus)      :: dftbp
        type(TDftbPlusInput) :: dftbp_input

        character(len=2)    , allocatable  :: atom_names(:)
        real(dp)            , allocatable  :: coor(:,:)
        real(dp)            , allocatable  :: coor_old(:,:)
        real(dp)            , allocatable  :: coor_new(:,:)
        real(dp)            , allocatable  :: forces(:,:)
        real(dp)            , allocatable  :: vel(:,:)
        real(dp)            , allocatable  :: at_masses(:)
        real(dp)            , allocatable  :: at_charges(:)

        type(fnode), pointer :: pRoot
        type(fnode), pointer :: pGeo 
        type(fnode), pointer :: pHam
        type(fnode), pointer :: pDftb
        type(fnode), pointer :: pMaxAng
        type(fnode), pointer :: pSlakos
        type(fnode), pointer :: pType2Files
        type(fnode), pointer :: pElecDyn
        type(fnode), pointer :: pElecStatic
        type(fnode), pointer :: pExternal
        type(fnode), pointer :: pPerturb
        type(fnode), pointer :: pLaser
        type(fnode), pointer :: pAnalysis

        contains
            procedure :: init => init_dftb
            procedure :: kill => kill_dftb
            procedure :: gs_calculate => gs_calculate_dftb
            procedure :: td_propagate => td_propagate_dftb
    end type TQ_sys_dftb

contains

!###################################################################################################

subroutine init_dftb(this, id, id_file, dt, t_steps, rank)

    class(TQ_sys_dftb), intent(inout) :: this
    integer           , intent(in)    :: id
    integer           , intent(in)    :: id_file
    integer           , intent(in)    :: t_steps
    integer           , intent(in)    :: rank
    real(dp)          , intent(in)    :: dt


    character(len=20)  :: file_name = "molecule_"
    character(len=20)  :: file_exten = ".dat"
    character(len=20)  :: file_number
    character(len=20)  :: input_name
    character(len=100) :: dyn_type_str
    integer            :: n_atoms
    integer            :: n_at_typ
    integer            :: euler_steps
    integer            :: i, ii, j
    integer            :: ierr
    integer            :: funit
    logical            :: atom_type_exists
    logical            :: periodic = .false. !This could be read from the file in the future
    logical            :: ion_dyn = .false.
    logical            :: scc = .true. !This could be read from the file in the future
    real(dp)           :: scc_tol

    character(len = 2), allocatable :: atom_type_list(:)
    character(len = 2), allocatable :: max_ang_orb(:)
    real(dp)          , allocatable :: atomNetCarges(:,:)
    real(dp)          , allocatable :: coor_aux(:,:)
    integer           , allocatable :: atom_type(:)

    this%id      = id
    this%id_file = id_file
    this%dt      = dt
    this%t_steps = t_steps
    this%rank    = rank

    this%std_unit = this%std_unit + this%rank

    write(file_number,'(I7.7)') id_file
    input_name = trim(file_name)//trim(file_number)//trim(file_exten)

    write(file_number,'(I7.7)') this%rank
    this%temp_file = "temp_output_"//trim(file_number)

    open (action='read', file=input_name, iostat=ierr, newunit=funit)

    if (ierr /= 0) then
        write(*,*) "Error: Could not open molecule file ", trim(input_name)
        stop
    end if

    read (unit=funit, fmt=*, iostat=ierr) n_atoms, n_at_typ
    read (unit=funit, fmt=*, iostat=ierr) dyn_type_str, euler_steps

    this%n_atoms = n_atoms
    this%n_at_typ = n_at_typ
    this%euler_steps = euler_steps

    select case (trim(dyn_type_str))
    case ('electrons')
        this%dyn_type = DFTB_ELEC_DYN
    case ('ehrenfest')
        this%dyn_type = DFTB_EHREN_DYN
    case ('born-oppenheimer')
        this%dyn_type = DFTB_BO_DYN
    case default
        write(*,*) "Error: Unknown DFTB dynamics type ", trim(dyn_type_str)
        stop
    end select

    if (.not. allocated(this%atom_names)) allocate(this%atom_names(n_atoms))
    if (.not. allocated(this%coor))       allocate(this%coor(3, n_atoms))
    if (.not. allocated(this%at_charges)) allocate(this%at_charges(n_atoms))
    if (.not. allocated(this%dipole))     allocate(this%dipole(3))
    if (.not. allocated(this%dip_old))    allocate(this%dip_old(3))
    if (.not. allocated(this%dPt_dt))     allocate(this%dPt_dt(3))
    if (.not. allocated(this%dPt_dt_old)) allocate(this%dPt_dt_old(3))

    this%coor       = M_ZERO
    this%at_charges = M_ZERO
    this%dipole     = M_ZERO
    this%dip_old    = M_ZERO
    this%dPt_dt     = M_ZERO
    this%dPt_dt_old = M_ZERO    

    if (.not. allocated(atom_type_list))  allocate(atom_type_list(n_at_typ))
    if (.not. allocated(max_ang_orb))     allocate(max_ang_orb(n_at_typ))
    if (.not. allocated(atom_type))       allocate(atom_type(n_atoms))
    if (.not. allocated(coor_aux))        allocate(coor_aux(3, n_atoms))
    
    if (this%dyn_type == DFTB_BO_DYN) then
        if (.not. allocated(this%coor_old))   allocate(this%coor_old(3, n_atoms))
        if (.not. allocated(this%coor_new))   allocate(this%coor_new(3, n_atoms))
        if (.not. allocated(this%forces))     allocate(this%forces(3, n_atoms))
        if (.not. allocated(this%vel))        allocate(this%vel(3, n_atoms))
        if (.not. allocated(this%at_masses))  allocate(this%at_masses(n_atoms))

        this%at_masses  = M_ZERO
        this%forces     = M_ZERO
        this%coor_old   = M_ZERO
        this%coor_new   = M_ZERO
        this%vel        = M_ZERO

    end if 

    read (unit=funit, fmt=*, iostat=ierr) (atom_type_list(i), i=1,n_at_typ)
    read (unit=funit, fmt=*, iostat=ierr) (max_ang_orb(i), i=1,n_at_typ)
    read (unit=funit, fmt=*, iostat=ierr) scc_tol

    read (unit=funit, fmt=*, iostat=ierr) 

    do i=1, n_atoms
        read (unit=funit, fmt=*, iostat=ierr) this%atom_names(i), &
                                              this%coor(1,i), &
                                              this%coor(2,i), &
                                              this%coor(3,i)

        do ii=1, n_at_typ
            atom_type_exists = .false.
            if (this%atom_names(i) == atom_type_list(ii)) then
                atom_type(i) = ii
                atom_type_exists = .true.
            end if
        end do
        if (.not. atom_type_exists) then
            write(*,*) "ERROR.", this%atom_names(i), "is not an atom type."
            stop
        end if
    end do

    !This file should save data printed by DFTB+. The program will overwrite it
    !in different stages to not accumulate too much data.

    
    open(newunit=this%std_unit, file=trim(this%temp_file), status='replace', action="write")

    call TDftbPlus_init(this%dftbp, outputUnit=this%std_unit)

    call this%dftbp%getEmptyInput(this%dftbp_input)
    call this%dftbp_input%getRootNode(this%pRoot)


    call setChild(this%pRoot, "Geometry", this%pGeo)
    call setChildValue(this%pGeo, "Periodic", periodic)
    call setChildValue(this%pGeo, "TypeNames", atom_type_list)

    coor_aux = M_ZERO
    call setChildValue(this%pGeo, "TypesAndCoordinates", reshape(atom_type, &
                       [1, size(atom_type)]), coor_aux)
    call setChild(this%pRoot, "Hamiltonian", this%pHam)
    call setChild(this%pHam, "Dftb", this%pDftb)
    call setChildValue(this%pDftb, "Scc", scc)
    call setChildValue(this%pDftb, "SccTolerance", scc_tol)

    call setChild(this%pDftb, "MaxAngularMomentum", this%pMaxAng)
    do j=1, n_at_typ
        call setChildValue(this%pMaxAng, atom_type_list(j), max_ang_orb(j))
    end do

    call setChild(this%pDftb, "SlaterKosterFiles", this%pSlakos)
    call setChild(this%pSlakos, "Type2FileNames", this%pType2Files)
    call setChildValue(this%pType2Files, "Prefix", "./")
    call setChildValue(this%pType2Files, "Separator", "-")
    call setChildValue(this%pType2Files, "Suffix", ".skf")

    if (this%dyn_type == DFTB_BO_DYN) then
        call setChild(this%pRoot, "Analysis", this%pAnalysis)
        call setChildValue(this%pAnalysis, "PrintForces", .true.)
        
        call setChild(this%pDftb, "ElectricField", this%pElecStatic)
        call setChild(this%pElecStatic, "External", this%pExternal)
        call setChildValue(this%pExternal, "Strength", 0.0_dp)
        call setChildValue(this%pExternal, "Direction", [1.0_dp, 1.0_dp, 1.0_dp])
    else
        !  set up electron dynamics options
        call setChild(this%pRoot, "ElectronDynamics", this%pElecDyn)
        call setChildValue(this%pElecDyn, "Steps", this%t_steps)
        call setChildValue(this%pElecDyn, "TimeStep", this%dt)
        call setChildValue(this%pElecDyn, "FieldStrength", 1.0_dp)
     
        if (this%dyn_type == DFTB_EHREN_DYN) ion_dyn = .true.
     
        call setChildValue(this%pElecDyn, "IonDynamics", ion_dyn)
        call setChildValue(this%pElecDyn, "InitialTemperature", 0.0_dp)
        call setChildValue(this%pElecDyn, "EulerFrequency", this%euler_steps)
        call setChildValue(this%pElecDyn, "VerboseDynamics", .false.)

        call setChild(this%pElecDyn, "Perturbation", this%pPerturb)
        call setChild(this%pPerturb, "Laser", this%pLaser)
    ! these twovalues will be overriden
        call setChildValue(this%pLaser, "PolarisationDirection", [1.0_dp, 1.0_dp, 1.0_dp])
        call setChildValue(this%pLaser, "LaserEnergy", 1.0_dp)
    end if

    call this%dftbp%setupCalculator(this%dftbp_input)

    close(funit)

    if (allocated(atom_type_list)) deallocate(atom_type_list)
    if (allocated(max_ang_orb))    deallocate(max_ang_orb)
    if (allocated(atom_type))      deallocate(atom_type)
    if (allocated(atomNetCarges))  deallocate(atomNetCarges)

end subroutine init_dftb

!###################################################################################################

subroutine kill_dftb(this)

    class(TQ_sys_dftb), intent(inout) :: this

    call TDftbPlus_destruct(this%dftbp)

    if (allocated(this%atom_names)) deallocate(this%atom_names)
    if (allocated(this%coor))       deallocate(this%coor)
    if (allocated(this%coor_old))   deallocate(this%coor_old)
    if (allocated(this%coor_new))   deallocate(this%coor_new)
    if (allocated(this%forces))     deallocate(this%forces)
    if (allocated(this%vel))        deallocate(this%vel)
    if (allocated(this%at_masses))  deallocate(this%at_masses)
    if (allocated(this%at_charges)) deallocate(this%at_charges)
    if (allocated(this%dipole))     deallocate(this%dipole)
    if (allocated(this%dip_old))    deallocate(this%dip_old)
    if (allocated(this%dPt_dt))     deallocate(this%dPt_dt)
    if (allocated(this%dPt_dt_old)) deallocate(this%dPt_dt_old)

    close(this%std_unit)

end subroutine kill_dftb

!###################################################################################################

subroutine gs_calculate_dftb(this)

    class(TQ_sys_dftb), intent(inout) :: this

    integer  :: i, kk

    this%E0      = M_ZERO
    this%Et      = M_ZERO
    this%dipole  = M_ZERO
    this%dip_old = M_ZERO

    do kk = 1, this%n_atoms
        this%coor(:, kk) = this%coor(:, kk) * AA__Bohr
    end do

    call this%dftbp%setGeometry(this%coor)
    call this%dftbp%getEnergy(this%E0)

    call this%dftbp%getGrossCharges(this%at_charges)

    do i = 1, this%n_atoms
        this%dipole(1) = this%dipole(1) + this%at_charges(i) * this%coor(1,i)
        this%dipole(2) = this%dipole(2) + this%at_charges(i) * this%coor(2,i)
        this%dipole(3) = this%dipole(3) + this%at_charges(i) * this%coor(3,i)
    end do

    if (this%dyn_type == DFTB_BO_DYN) then
        call this%dftbp%getAtomicMasses(this%at_masses)
        this%Kinetic_Energy = M_ZERO
    else
    !This nitialize time propagation for Ehrenfest dynamics as well
    !TO-DO: maybe this should be somewhere else.
        call this%dftbp%initializeTimeProp(this%dt, .true., .false.)
    end if

    if (allocated(this%coor_old)) this%coor_old = this%coor
    this%dip_old = this%dipole

end subroutine gs_calculate_dftb

!###################################################################################################

subroutine td_propagate_dftb(this, tq_step, E_field)

    class(TQ_sys_dftb), intent(inout) :: this
    integer           , intent(in)    :: tq_step
    real(dp)          , intent(in)    :: E_field(3)
    
    integer  :: i
    real(dp) :: E_amp
    real(dp) :: vv

    real(dp), allocatable :: aux_at_charges(:,:)
    real(dp), allocatable :: dip_aux(:,:)

    if(.not. allocated(aux_at_charges)) allocate(aux_at_charges(this%n_atoms, 1))
    if(.not. allocated(dip_aux)) allocate(dip_aux(3, 1))

    E_amp = dsqrt(E_field(1)**2 + E_field(2)**2 + E_field(3)**2)

    this%dPt_dt_old = this%dPt_dt

    if (this%dyn_type == DFTB_BO_DYN) then
        this%Kinetic_Energy = M_ZERO

        call this%dftbp%setGeometry(this%coor)
        call this%dftbp%setExternalEfield(E_amp, E_field)
        call this%dftbp%getGradients(this%forces) !The force has the opposite sign of the gradient.
        call this%dftbp%getEnergy(this%Et)
        call this%dftbp%getGrossCharges(this%at_charges)

        this%dip_old = this%dipole
        this%dipole  = M_ZERO

        do i = 1, this%n_atoms
            this%coor_new(:, i) = 2.0d0*this%coor(:, i) - this%coor_old(:, i) - &
                                 (this%forces(:, i) / this%at_masses(i)) * this%dt**2
            this%dipole(:) = this%dipole(:) + this%at_charges(i) * this%coor(:,i)

            this%vel(:, i) = (this%coor_new(:, i) - this%coor_old(:, i)) / (2.0d0 * this%dt)
            vv = this%vel(1,i)**2 + this%vel(2,i)**2 + this%vel(3,i)**2

            this%Kinetic_Energy = this%Kinetic_Energy + 0.5d0 * this%at_masses(i) * vv

        end do
    
        this%coor_old = this%coor
        this%coor     = this%coor_new

        this%dPt_dt = (this%dipole - this%dip_old) / this%dt

    else ! electron dynamics and Ehrenfest dynamics.

        this%dip_old = this%dipole
        call this%dftbp%setTdElectricField(E_field)

        call this%dftbp%doOneTdStep(tq_step, dipole=dip_aux, energy=this%Et,  &
                                   atomNetCharges=aux_at_charges, coord=this%coor, &
                                   force=this%forces)

        this%dipole     = dip_aux(:,1)
        this%at_charges = aux_at_charges(:,1)

        this%dPt_dt = (this%dipole - this%dip_old) / this%dt

    end if

    if (allocated(aux_at_charges)) deallocate(aux_at_charges)
    if (allocated(dip_aux))        deallocate(dip_aux)

end subroutine td_propagate_dftb

!###################################################################################################
end module q_sys_dftb_mod