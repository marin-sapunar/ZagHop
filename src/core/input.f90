!--------------------------------------------------------------------------------------------------
! MODULE: input_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date December, 2016
!
! DESCRIPTION: 
!> @brief Input subroutines.
!> @details
!! Subroutines for setting default program options, reading command line arguments and reading 
!! program input files. 
!! Arrays of dimension natom are allocated in the ReadGeom subroutine. 
!! Arrays of dimension nstate are allocated in the ReadMainInp subroutine.
!--------------------------------------------------------------------------------------------------
module input_mod
    use global_defs
    use system_var
    use control_var
    use string_mod
    use file_mod, only : reader
    use constants

    implicit none

    private
    public :: readinput

    character(len=:), allocatable :: maininp !< Name of main input file.
    character(len=:), allocatable :: geominp !< Name of geometry input file.
    character(len=:), allocatable :: veloinp !< Name of velocity input file.
    character(len=:), allocatable :: pcinp !< Name of partial charges input file.
    real(dp) :: tol_cns

contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: DefaultOptions
    !
    ! DESCRIPTION:
    !> @brief Set program options to their default values.
    !> @details
    !! See manual for details concerning the input.
    !> @note Default values should be changeable without causing problems elsewhere in the program.
    !----------------------------------------------------------------------------------------------
    subroutine defaultoptions()
        ! File names.
        maininp = 'dynamics.in'
        geominp = 'geom'
        veloinp = 'veloc'
        pcinp = 'pcharge'
        ctrl%output_dir = 'Results'
        ctrl%bufile = 'backup.dat'
        ctrl%buinterval = 1
        ! Method options.
        ctrl%qext = .true.
        ! Initial state.
        t(1)%nstate = 1
        t(1)%cstate = 1
        ! Set start step at 0
        t(1)%step = 0
        t(1)%time = 0.0_dp
        ! Print options.
        ctrl%print = .false.
        ctrl%print(1:6) = .true.
        ctrl%print(7:10) = .true.
        ctrl%print(11) = .true.
        ! Program flow options.
        ctrl%min01gap = -1.0_dp
        ctrl%max_tot_en_change = 1.0_dp / eh_eV
        ctrl%max_tot_en_change_step = 0.2_dp / eh_eV
        ctrl%target_state = 0
        ctrl%target_state_time = 0.0_dp
        ! Nuclear dynamics options.
        ctrl%max_time = 10.0_dp / aut_fs
        ctrl%dt_0 = 0.5_dp / aut_fs
        ctrl%dt = ctrl%dt_0
        ctrl%lz_min_dt = ctrl%dt_0 / 20.0_dp
        ctrl%lz_prob_conv = 0.05_dp
        ctrl%orientlvl = 0
        tol_cns = 0.000000000001_dp ! Default rattle tolerance.
        ! Surface hopping options.
        ctrl%sh = 2
        ctrl%shnstep = 10000
        ctrl%decohlvl = 1
        ctrl%couplvl = 0
        ctrl%coupndiff = 1
        ctrl%coupediff = 1 / eh_eV
        ctrl%tdc_type = 1
        ctrl%tdc_interpolate = 2
        ctrl%ene_interpolate = 2
        ctrl%vrescale = 1
        ctrl%fhop = 1
        ctrl%phaselvl = 1
        ctrl%oscill = .false.
        ! MM options.
        ctrl%mmcut = 10000.0_dp
        ctrl%pbc = .false.
        ! Point charges
        ctrl%pcharge = .false.
    end subroutine defaultoptions


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ReadInput
    !
    ! DESCRIPTION:
    !> @brief Read all program input.
    !> @details 
    !! Read the command line arguments passed to the program, the main program input file and the 
    !! files containing the initial geometry and velocity of the system.
    !! At the moment, only the name of the main input file can be passed through the command line 
    !! arguments, all other options are set in this file.
    !----------------------------------------------------------------------------------------------
    subroutine readinput
        use nuclear_dyn_mod
        integer :: narg
        character(len=1000) :: temp
        integer :: i

        ! Get work directory.
        call getcwd(temp)
        ctrl%maindir = trim(adjustl(temp))//'/'

        ! Set default options.
        call defaultoptions()

        ! Get name of the main input file.
        narg = command_argument_count()
        if (narg > 0) then
            call get_command_argument(1, temp)
            maininp = temp
        end if

         ! Read main input file.
        write(stdout, *)
        write(stdout, '(1x,a,a,a)') 'Reading dynamics input file ', trim(maininp), '.'
        call readmaininp()

        ! Set initial directory.
        ctrl%qmdir = 'qmdir'
        call get_environment_variable('QMDIR', temp)
        if (temp /= '') ctrl%qmdir = trim(adjustl(temp))
        write(stdout, *) 
        write(stdout, '(1x,a,a)') "Work directory for QM calculations: ", ctrl%qmdir
        call system('mkdir -p '//ctrl%qmdir)
       
        ! If restarting don't read initial conditions.
        if (ctrl%restart) return

        ! Read initial conditions.
        write(stdout, *)
        write(stdout, '(1x,a,a,a)') 'Reading geometry file ', geominp, '.'
        call readgeom()
        write(stdout, '(1x,a,a,a)') 'Reading velocity file ', veloinp, '.'
        call readvelo()

        ! Reorient the molecule based on orientlvl
        select case(ctrl%orientlvl)
        case(1)
            call set_geom_center_of_mass(t(1)%mass, t(1)%geom)
        case(2)
            call set_geom_center_of_mass(t(1)%mass, t(1)%geom)
            call project_translation_rotation_from_velocity(t(1)%mass, t(1)%geom, t(1)%velo)
        end select


        ! Read list of partial charges.
        if (ctrl%pcharge) then
            write(stdout, '(1x,a,a,a)') 'Reading partial charges file ', pcinp, '.'
            call readpc()
        end if

        ! Print initial information about the full system.
        write(stdout, *)
        write(stdout, '(1x,a,i0,a)') 'Starting calculation for a system with ', t(1)%natom, ' atoms.'
        write(stdout, '(3x,a,i0,a)') 'QM system consists of ', t(1)%qnatom, '.'
        write(stdout, '(3x,a,i0,a)') 'MM system consists of ', t(1)%mnatom, '.'

        if (ctrl%qm) then
            write(stdout, *)
            write(stdout, '(3x,a)') 'Initial geometry of the QM system: '
            do i = 1, t(1)%qnatom
                write(stdout, '(2x,a5,1000e18.10)') t(1)%sym(t(1)%qind(i)), t(1)%geom(:, t(1)%qind(i))
            end do
            write(stdout, '(3x,a)') 'Initial velocity of the QM system: '
            do i = 1, t(1)%qnatom
                write(stdout, '(2x,a5,1000e18.10)') t(1)%sym(t(1)%qind(i)), t(1)%velo(:, t(1)%qind(i))
            end do
        end if
        if (ctrl%mm) then
            write(stdout, *)
            write(stdout, '(3x,a)') 'Initial geometry of the MM system: '
            do i = 1, t(1)%mnatom
                write(stdout, '(2x,a5,1000e18.10)') t(1)%sym(t(1)%mind(i)), t(1)%geom(:, t(1)%mind(i))
            end do
            write(stdout, '(3x,a)') 'Initial velocity of the MM system: '
            do i = 1, t(1)%mnatom
                write(stdout, '(2x,a5,1000e18.10)') t(1)%sym(t(1)%mind(i)), t(1)%velo(:, t(1)%mind(i))
            end do
        end if

    end subroutine readinput


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ReadMainInp
    !
    ! DESCRIPTION:
    !> @brief Read the main program control variables.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine readmaininp()
        use matrix_mod, only : unit_mat
        logical :: check
        logical :: buffer
        type(reader) :: readf
        integer, allocatable :: skipstates(:) !< States for which no couplings are calculated.
        integer :: i
        buffer = .false.

        inquire(file=maininp, exist=check)
        if (.not. check) then
            write(stderr,*) 'Error in Input module, ReadMainInp subroutine.'
            write(stderr,*) ' Main input file (', trim(maininp),') not found.'
            stop
        end if

        call readf%open(maininp, comment='#')
        ! Mandatory keywords are read first:
        call readmethod(readf)
        call readsystem(readf, t(1)%nstate, t(1)%cstate, geominp, veloinp, t(1)%ndim)

        ! Followed by optional keywords.
        call readdynamics(readf)
        call readoutput(readf)
        call readsurfhop(readf, skipstates)
        call readmm(readf)
        call readconstraints(readf, tol_cns)
        call readrestart(readf)
        call readf%close()

        ! Modify options based on type of calculation.
        if (t(1)%nstate == 1) then
            ctrl%sh = 0
            ctrl%print(11) = .false.
        end if
        if (ctrl%print(11)) ctrl%oscill = .true.
        if (ctrl%sh == 0) then
            ctrl%print(6:9) = .false.
        else if (ctrl%sh == 1) then
            ctrl%print(6:9) = .false.
        else if (ctrl%sh >= 2) then
            allocate(ctrl%couple(t(1)%nstate), source=.true.)
            if (allocated(skipstates)) then
                do i = 1, size(skipstates, 1)
                    if (skipstates(i) > t(1)%nstate) exit
                    ctrl%couple(skipstates(i)) = .false.
                end do
            end if
        end if

        ! Allocate all arrays of size t(1)%nstate.
        if (.not. ctrl%restart) then
            allocate(t(1)%qe(t(1)%nstate))
            if (ctrl%oscill) allocate(t(1)%qo(t(1)%nstate - 1))
            if (ctrl%sh == 1) then
                allocate(t(1)%prob(t(1)%nstate), source=0.0_dp)
            else if (ctrl%sh >= 2) then
                allocate(t(1)%cwf(t(1)%nstate))
                allocate(t(1)%prob(t(1)%nstate), source=0.0_dp)
                select case (ctrl%tdc_type)
                case(1)
                    allocate(t(1)%olap(t(1)%nstate, t(1)%nstate))
                    t(1)%olap = unit_mat(t(1)%nstate)
                case(2)
                    ! Nonadiabatic coupling vecotrs are allocated after reading number of atoms.
                end select
                t(1)%cwf = cmplx((0.0_dp, 0.0_dp), kind = dp)
                t(1)%cwf(t(1)%cstate) = cmplx((1.0_dp, 0.0_dp), kind = dp)
            end if
        end if


    end subroutine readmaininp

 
    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ReadGeom
    !
    ! DESCRIPTION:
    !> @brief Read the list of symbols, masses and positions of all atoms in the full system.
    !> @details 
    !! Read number of atoms in the full system and the QM and MM subsystems. Allocate all arrays
    !! of dimensions natom, qnatom or mnatom. Read the symbol, mass and position arrays for
    !! the full system.
    !! Assign QM and MM pointers to the appropriate sections of the main arrays.
    !----------------------------------------------------------------------------------------------
    subroutine readgeom()
        logical :: check
        character(len=2) :: tsym
        real(dp) :: tmass
        character :: part
        integer :: i
        integer :: d
        type(reader) :: readf

        inquire(file=geominp, exist=check)
        if (.not. check) then
            write(stderr,*) 'Error in Input module, ReadGeom subroutine.'
            write(stderr,*) ' Geometry file (', trim(geominp),') not found.'
            stop
        end if

        ! First run through the file to check number of atoms in the QM and MM part of the system.
        call readf%open(geominp, abort_on_eof=.false., comment='#')
        do
            call readf%next()
            if (is_iostat_end(readf%iostat)) exit
            call readf%parseline(' ,')
            if (readf%narg < t(1)%ndim + 3) then
                write(stderr, '(a)') 'Error in Input module, ReadGeom subroutine.'
                write(stderr, '(a,i0,a)') ' Need ', 3+t(1)%ndim, ' arguments per line:'
                write(stderr, '(a)') '  symbol mass coords q/m'
                write(stderr, '(a,a)') '  Line: ', readf%line
                stop 
            end if
            read(readf%args(3+t(1)%ndim)%s, *) part
            select case (tolower(part))
            case('q')
                t(1)%qnatom = t(1)%qnatom + 1
            case('m')
                t(1)%mnatom = t(1)%mnatom + 1
            case default
                write(stderr, *) 'Error in Input module, ReadGeom subroutine.'
                write(stderr, *) ' Lines should end with "q" for QM atoms or "m" for MM atoms.'
                stop
            end select
        end do
        call readf%rewind()

        t(1)%natom = t(1)%qnatom + t(1)%mnatom
        if (t(1)%natom == 0) then
            write(stderr, *) 'Error in Input module, ReadGeom subroutine.'
            write(stderr, *) ' No atoms found in geometry file!'
            stop
        end if

        ! Allocate arrays of natom dimensions.
        allocate(t(1)%sym(t(1)%natom))
        allocate(t(1)%mass(t(1)%natom))
        allocate(t(1)%chrg(t(1)%natom))
        allocate(t(1)%geom(t(1)%ndim, t(1)%natom))
        allocate(t(1)%velo(t(1)%ndim, t(1)%natom))
        allocate(t(1)%grad(t(1)%ndim, t(1)%natom))
        if (t(1)%qnatom > 0) then
            ctrl%qm = .true.
            allocate(t(1)%qind(t(1)%qnatom))
        else
            ctrl%qm = .false.
        end if
        if (t(1)%mnatom > 0) then
            ctrl%mm = .true.
            ctrl%print(12) = .true.
            allocate(t(1)%mind(t(1)%mnatom))
        else
            ctrl%mm = .false.
        end if
        if (ctrl%tdc_type == 2) then
            allocate(t(1)%nadv(t(1)%ndim*t(1)%natom, t(1)%nstate, t(1)%nstate))
        end if

        ! Second run through the file to read the atom, mass, position and index of each atom.
        t(1)%qnatom = 0
        t(1)%mnatom = 0
        do i = 1, t(1)%natom
            call readf%next()
            call readf%parseline(' ,')

            ! Read symbol.
            read(readf%args(1)%s, *) tsym
            t(1)%sym(i) = tolower(tsym)

            ! Read mass.
            read(readf%args(2)%s, *) tmass
            t(1)%mass(i) = tmass * Da_me

            ! Read initial geometry.
            do d = 1, t(1)%ndim
                read(readf%args(2+d)%s, *) t(1)%geom(d, i)
            end do

            ! Check part (q = qm, m = mm)
            read(readf%args(3+t(1)%ndim)%s, *) part
            select case (tolower(part))
            case('q')
                t(1)%qnatom = t(1)%qnatom + 1
                t(1)%qind(t(1)%qnatom) = i
            case('m')
                t(1)%mnatom = t(1)%mnatom + 1
                t(1)%mind(t(1)%mnatom) = i
            end select
        end do
        call readf%close()
    end subroutine readgeom

 
    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ReadVelo
    !
    ! DESCRIPTION:
    !> @brief Read the velocities of the atoms in the system.
    !----------------------------------------------------------------------------------------------
    subroutine readvelo()
        logical :: check
        integer :: i
        integer :: d
        type(reader) :: readf

        inquire(file=veloinp, exist=check)
        if (.not. check) then
            write(stderr, *) 'Error in Input module, ReadVelo subroutine.'
            write(stderr, *) ' Velocity file (', trim(veloinp),') not found.'
            stop
        end if

        call readf%open(veloinp, comment='#')
        do i = 1, t(1)%natom
            call readf%next()
            call readf%parseline(' ,')
            if (readf%narg < t(1)%ndim) then
                write(stderr, *) 'Error in Input module, ReadVelo subroutine.'
                write(stderr, '(a,i0,a)') ' Premature end of line, expect ', t(1)%ndim, ' arguments.'
                write(stderr, '(a,a)') ' Line: ', readf%line
                stop
            end if
            do d = 1, t(1)%ndim
                read(readf%args(d)%s, *) t(1)%velo(d, i)
            end do
        end do
        call readf%close()
    end subroutine readvelo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ReadPC
    !
    ! DESCRIPTION:
    !> @brief Read the partial charges of the atoms in the system.
    !----------------------------------------------------------------------------------------------
    subroutine readpc()
        logical :: check
        integer :: i
        type(reader) :: readf

        inquire(file=pcinp, exist=check)
        if (.not. check) then
            write(stderr, *) 'Error in Input module, ReadPC subroutine.'
            write(stderr, *) ' Point charges file (', trim(pcinp),') not found.'
            stop
        end if

        call readf%open(pcinp)
        do i = 1, t(1)%natom
            call readf%next()
            read(readf%line, *) t(1)%chrg(i)
        end do
        call readf%close()
    end subroutine readpc
 

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ReadMethod
    !
    ! DESCRIPTION:
    !> @brief Read the $method section of the main input file.
    !> @details 
    !! In this section, programs to run the QM and MM calculations are selected.
    !----------------------------------------------------------------------------------------------
    subroutine readmethod(readf)
        type(reader), intent(inout) :: readf

        if (.not. allocated(ctrl%qprog)) ctrl%qprog = ''
        if (.not. allocated(ctrl%mprog)) ctrl%mprog = ''

        write(stdout, '(3x,a)') 'Reading $method section of the dynamics input file.'
        call readf%rewind()
        call readf%go_to_keyword('$method')
        do 
            call readf%next()
            if (index(readf%line, '$') == 1) exit
            call readf%parseline(' ')
            select case(readf%args(1)%s)
            case('qm')
                ctrl%qprog = readf%args(2)%s
                write(stdout, '(5x,a)') '"'//ctrl%qprog//'" will be used for QM calculations.'
            case('mm')
                ctrl%mprog = readf%args(2)%s
                write(stdout, '(5x,a)') '"'//ctrl%mprog//'" will be used for MM calculations.'
            case('overlap')
                ctrl%oprog = readf%args(2)%s
                write(stdout, '(5x,a)') '"'//ctrl%oprog//'" will be used for overlap calculations.'
            case('model')
                ctrl%qext = .false.
            case('qm_en_error')
                read(readf%args(2)%s, *) ctrl%qm_en_err
            end select
        end do
    end subroutine readmethod


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ReadSystem
    !
    ! DESCRIPTION:
    !> @brief Read information about the system.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine readsystem(readf, nstate, cstate, geominp, veloinp, ndim)
        type(reader), intent(inout) :: readf
        integer, intent(inout) :: nstate
        integer, intent(out) :: cstate
        character(len=:), allocatable :: geominp
        character(len=:), allocatable :: veloinp
        integer, intent(out) :: ndim

        call readf%rewind()
        call readf%go_to_keyword('$system')
        do
            call readf%next()
            if (index(readf%line, '$') == 1) exit
            call readf%parseline(' =')
            select case(readf%args(1)%s)
            case('nstate')
                read(readf%args(2)%s, *) nstate
            case('istate')
                read(readf%args(2)%s, *) cstate
                cstate = cstate
            case('geometry')
                geominp = readf%args(2)%s
            case('velocity')
                veloinp = readf%args(2)%s
            case('ndim')
                read(readf%args(2)%s, *) ndim
            case default
                write(stderr, *) 'Warning in InputSection module, ReadSystem subroutine.' 
                write(stderr, *) '  Skipping line with unrecognized keyword:' 
                write(stderr, *) readf%line
            end select
        end do

        if (nstate == 0) then
            write(stderr,*) 'Error in InputSection module, ReadSystem subroutine.'
            write(stderr,'(a)') ' Number of states not defined.'
            stop
        end if
        if (cstate == 0) then
            write(stderr,*) 'Warning in InputSection module, ReadSystem subroutine.'
            write(stderr,'(a)') ' Initial state not defined, setting initial state = ', nstate, '.'
            cstate = nstate
        end if

    end subroutine readsystem


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ReadOutput
    !
    ! DESCRIPTION:
    !> @brief Read output options of the program.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine readoutput(readf)
        type(reader), intent(inout) :: readf
        integer, allocatable :: printopts(:)
        logical :: check
        integer :: i

        call readf%rewind()
        call readf%go_to_keyword('$output', found=check)
        if (.not. check) return
        do
            call readf%next()
            if (index(readf%line, '$') == 1) exit
            call readf%parseline(' =')
            select case(readf%args(1)%s)
            case('results_dir')
                ctrl%output_dir = readf%args(2)%s
            case('print')
                call readf%parseline('t')
                call read_index_list(readf%args(2)%s, printopts)
                do i = 1, size(printopts)
                    ctrl%print(printopts(i)) = .true.
                end do
            case('noprint')
                call readf%parseline('t')
                call read_index_list(readf%args(2)%s, printopts)
                do i = 1, size(printopts)
                    ctrl%print(printopts(i)) = .false.
                end do
            case('bufile')
                ctrl%bufile = readf%args(2)%s
            case('buinterval')
                read(readf%args(2)%s, *) ctrl%buinterval
            case default
                write(stderr, *) 'Warning in InputSection module, ReadOutput subroutine.' 
                write(stderr, *) '  Skipping line with unrecognized keyword:' 
                write(stderr, *) readf%line
            end select
        end do 
    end subroutine readoutput


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ReadDynamics
    !
    ! DESCRIPTION:
    !> @brief Read information about the nuclear dynamics.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine readdynamics(readf)
        type(reader), intent(inout) :: readf
        logical :: check

        call readf%rewind()
        call readf%go_to_keyword('$dynamics', found=check)
        if (.not. check) return
        do
            call readf%next()
            if (index(readf%line, '$') == 1) exit
            call readf%parseline(' =')
            select case(readf%args(1)%s)
            case('time', 'max_time')
                read(readf%args(2)%s, *) ctrl%max_time
                ctrl%max_time = ctrl%max_time / aut_fs
            case('tstep', 'step', 'dt')
                read(readf%args(2)%s, *) ctrl%dt
                ctrl%dt = ctrl%dt / aut_fs
                ctrl%dt_0 = ctrl%dt
            case('end_state')
                read(readf%args(2)%s, *) ctrl%target_state
            case('end_state_time')
                read(readf%args(2)%s, *) ctrl%target_state_time
                ctrl%target_state_time = ctrl%target_state_time / aut_fs
            case('max_toten_d')
                read(readf%args(2)%s, *) ctrl%max_tot_en_change
                ctrl%max_tot_en_change = ctrl%max_tot_en_change / eh_eV
            case('max_toten_d_step')
                read(readf%args(2)%s, *) ctrl%max_tot_en_change_step
                ctrl%max_tot_en_change_step = ctrl%max_tot_en_change_step / eh_eV
            case('orient')
                read(readf%args(2)%s, *) ctrl%orientlvl
            case default
                write(stderr, *) 'Warning in InputSection module, ReadDynamics subroutine.' 
                write(stderr, *) '  Skipping line with unrecognized keyword:' 
                write(stderr, *) readf%line
            end select
        end do
    end subroutine readdynamics


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ReadConstraints
    !
    ! DESCRIPTION:
    !> @brief Read information about constraints.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine readconstraints(readf, default_tol)
        type(reader), intent(inout) :: readf
        real(dp), intent(in) :: default_tol
        logical :: check
        integer :: i

        ctrl%n_cns = 0
        call readf%rewind()
        call readf%go_to_keyword('$constraints', found=check)
        if (.not. check) return
        ctrl%n_cns = readf%count_lines_to_keyword('$')
        if (ctrl%n_cns == 0) return

        ctrl%constrain_dyn = .true.
        allocate(ctrl%cns(ctrl%n_cns))
        ctrl%cns(:)%tol = default_tol
        do i = 1, ctrl%n_cns
            call readf%next()
            call readf%parseline(' =,')
            select case(readf%args(1)%s)
            case('b')
                read(readf%args(2)%s, *) ctrl%cns(i)%i(1)
                read(readf%args(3)%s, *) ctrl%cns(i)%i(2)
                read(readf%args(4)%s, *) ctrl%cns(i)%cval
                if (readf%narg > 4) read(readf%args(5)%s, *) ctrl%cns(i)%tol
            case default
                write(stderr, *) 'Error in InputSection module, ReadConstraints subroutine.' 
                write(stderr, *) '  Unrecognized constraint type: ', readf%args(1)%s
                stop
            end select
        end do
    end subroutine readconstraints


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ReadSurfHop
    !
    ! DESCRIPTION:
    !> @brief Read the $surfhop section of the main input file.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine readsurfhop(readf, skipstates)
        type(reader), intent(inout) :: readf
        integer, allocatable, intent(out) :: skipstates(:)
        logical :: check

        call readf%rewind()
        call readf%go_to_keyword('$surfhop', found=check)
        if (.not. check) return

        call readf%parseline(' =')
        if (readf%narg > 1) then
            select case(readf%args(2)%s)
            case('off')
                ctrl%sh = 0
                ctrl%tdc_type = 0
                return
            end select
        end if

        do
            call readf%next()
            if (index(readf%line, '$') == 1) exit
            call readf%parseline(' =')
            select case(readf%args(1)%s)
            case('lz')
                ctrl%sh = 1
                ctrl%tdc_type = 0
            case('fssh')
                ctrl%sh = 2
                if (readf%narg > 1) then
                    select case(readf%args(2)%s)
                    case('adiabatic')
                        ctrl%sh = 2
                    case('diabatic', 'ldiab', 'diab')
                        ctrl%sh = 3
                    end select
                end if
            case('tdse_steps')
                read(readf%args(2)%s, *) ctrl%shnstep
            case('decoherence')
                select case(readf%args(2)%s)
                case('off')
                    ctrl%decohlvl = 0
                case('nldm')
                    ctrl%decohlvl = 1
                case default
                    write(stderr, *) 'Error in InputSection module, ReadSurfHop subroutine.'
                    write(stderr, '(2x,a)') '  Unrecognized decoherence method: '
                    write(stderr, '(2x,a)') readf%line
                    stop
                end select
            case('overlap')
                ctrl%tdc_type = 1
                if (readf%narg > 1) then
                    select case(readf%args(2)%s)
                    case('constant')
                        ctrl%tdc_interpolate = 1
                    case('linear')
                        ctrl%tdc_interpolate = 2
                    case('npi')
                        ctrl%tdc_interpolate = 4
                    end select
                end if
            case('nadvec')
                ctrl%tdc_type = 2
                ctrl%print(9) = .false. !< @todo Move this.
            case('energy')
                select case(readf%args(2)%s)
                case('constant')
                    ctrl%ene_interpolate = 0
                case('step')
                    ctrl%ene_interpolate = 1
                case('linear')
                    ctrl%ene_interpolate = 2
                end select
            case('velocity')
                select case(readf%args(2)%s)
                case('off')
                    ctrl%vrescale = 0
                case('vel')
                    ctrl%vrescale = 1
                case('gdif')
                    ctrl%vrescale = 2
                case('nadvec')
                    ctrl%vrescale = 3
                end select
            case('frustrated')
                select case(readf%args(2)%s)
                case('continue')
                    ctrl%fhop = 1
                case('reverse')
                    ctrl%fhop = 2
                end select
            case('phase')
                select case(readf%args(2)%s)
                case('off')
                    ctrl%phaselvl = 0
                case('diagonal')
                    ctrl%phaselvl = 1
                case('assigned')
                    ctrl%phaselvl = 2
                case default
                    write(stderr, *) 'Error in InputSection module, ReadSurfHop subroutine.'
                    write(stderr, '(2x,a)') '  Unrecognized phase matching algorithm: '
                    write(stderr, '(2x,a)') readf%line
                    stop
                end select
            case('lz_prob_conv')
                read(readf%args(2)%s, *) ctrl%lz_prob_conv
            case('s0_ci')
                read(readf%args(2)%s, *) ctrl%min01gap
                ctrl%min01gap = ctrl%min01gap / eh_eV
            case('skip_state')
                call read_index_list(readf%line(11:), skipstates)
            case('seed')
                read(readf%args(2)%s, *) ctrl%seed
            case default
                write(stderr, *) 'Warning in InputSection module, ReadSurfHop subroutine.'
                write(stderr, *) '  Skipping line with unrecognized keyword:' 
                write(stderr, *) readf%line
            end select
        end do
    end subroutine readsurfhop

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ReadMM
    !
    ! DESCRIPTION:
    !> @brief Read the $mm section of the main input file.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine readmm(readf)
        type(reader), intent(inout) :: readf
        logical :: check

        call readf%rewind()
        call readf%go_to_keyword('$mm', found=check)
        if (.not. check) return
        do
            call readf%next()
            if (index(readf%line, '$') == 1) exit
            call readf%parseline(' =')
            select case(readf%args(1)%s)
            case('cutoff')
                read(readf%args(2)%s, *) ctrl%mmcut
            case('pcharge')
                ctrl%pcharge = .true.
                if (readf%narg > 1) then
                    if (readf%args(2)%s == 'off') ctrl%pcharge = .false.
                end if
            case('pcharge_file')
                pcinp = readf%args(2)%s
            case('pbc')
                ctrl%pbc = .true.
                if (readf%narg > 1) then
                    if (readf%args(2)%s == 'off') ctrl%pbc = .false.
                end if
                if (ctrl%pbc) then
                    call readf%next()
                    read (readf%line, *) t(1)%pbcbox
                end if
            case default
                write(stderr, *) 'Warning in InputSection module, ReadMM subroutine.' 
                write(stderr, *) '  Skipping line with unrecognized keyword:' 
                write(stderr, *) readf%line
            end select
        end do
    end subroutine readmm


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ReadRestart
    !
    ! DESCRIPTION:
    !> @brief Check if the $restart keyword is present in the main input file.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine readrestart(readf)
        type(reader), intent(inout) :: readf
        logical :: check

        call readf%rewind()
        call readf%go_to_keyword('$restart', found=check)
        if (check) ctrl%restart = .true.
    end subroutine readrestart

 
end module input_mod
