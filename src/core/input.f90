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
!! Arrays of dimension natom are allocated in the read_geom subroutine.
!! Arrays of dimension nstate are allocated in the read_main subroutine.
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
    public :: set_defaults
    public :: read_input
    public :: command_line_interface

    character(len=:), allocatable :: maininp !< Name of main input file.
    character(len=:), allocatable :: geominp !< Name of geometry input file.
    character(len=:), allocatable :: veloinp !< Name of velocity input file.
    character(len=:), allocatable :: pcinp !< Name of partial charges input file.
    real(dp) :: tol_cns
    real(dp) :: mb_temperature = -1.0_dp

contains



    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: command_line_interface
    !> @brief Read command line options.
    !----------------------------------------------------------------------------------------------
    subroutine command_line_interface()
        integer :: narg, i
        character(len=1000) :: temp
        integer :: print_level

        narg = command_argument_count()
        if (narg == 0) return
        i = 0
        do while (i < narg)
            i = i + 1
            call get_command_argument(i, temp)
            if ((i == narg) .and. (temp(1:1) /= '-')) then
                maininp = trim(adjustl(temp))
                continue
            end if
            select case(temp)
            case('--help', '-h')
                call print_help()
            case('--print-level', '-p')
                i = i + 1
                call get_command_argument(i, temp)
                read(temp, *) print_level
                stdp1 = .false.
                stdp2 = .false.
                stdp3 = .false.
                stdp4 = .false.
                if (print_level >= 1) stdp1 = .true.
                if (print_level >= 2) stdp2 = .true.
                if (print_level >= 3) stdp3 = .true.
                if (print_level >= 4) stdp4 = .true.
            end select
        end do
    end subroutine command_line_interface


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: set_defaults
    !
    ! DESCRIPTION:
    !> @brief Set program options to their default values.
    !> @details
    !! See manual for details concerning the input.
    !> @note Default values should be changeable without causing problems elsewhere in the program.
    !----------------------------------------------------------------------------------------------
    subroutine set_defaults()
        use random_pcg_mod, only : rng_pcg_xsh_rr
        ! File names.
        maininp = 'dynamics.in'
        geominp = 'geom'
        veloinp = 'veloc'
        pcinp = 'pcharge'
        ctrl%output_dir = 'Results'
        ctrl%bufile = 'backup.dat'
        ctrl%printerval = 1
        ctrl%buinterval = 1
        ! RNG.
        allocate(rng_pcg_xsh_rr::ctrl%rng)
        ctrl%rng%seed = -1
        ! Method options.
        ctrl%qlib = 0
        ctrl%noise = 0.0_dp
        ! Set start step at 0
        tr1%step = 0
        tr1%time = 0.0_dp
        ! Print options.
        ctrl%print = .false.
        ctrl%print(1:7) = .true.
        ctrl%print(10) = .true.
        ! Program flow options.
        ctrl%stop_s0s1_ci = -1.0_dp
        ctrl%max_tot_en_change = 1.0_dp / eh_eV
        ctrl%max_tot_en_change_step = 0.2_dp / eh_eV
        ctrl%target_state = 0
        ctrl%target_state_time = 0.0_dp
        ! Nuclear dynamics options.
        ctrl%max_time = 10.0_dp / aut_fs
        ctrl%dt_0 = 0.5_dp / aut_fs
        ctrl%dt = ctrl%dt_0
        ctrl%lz_min_dt = 0.05_dp / aut_fs
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
        ctrl%vrescale = 2
        ctrl%fhop = 1
        ctrl%phaselvl = 1
        ctrl%variable_nstate = 0
        ctrl%oscill = .false.
        ! Thermostat
        ctrl%thermostat = 0
        ctrl%target_t = -1.0_dp
        ctrl%tau_t = -1.0_dp
        ! MM options.
        ctrl%mmcut = 10000.0_dp
        ctrl%pbc = .false.
        ! Point charges
        ctrl%pcharge = .false.
    end subroutine set_defaults


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_input
    !
    ! DESCRIPTION:
    !> @brief Read all program input.
    !> @details
    !! Read the command line arguments passed to the program, the main program input file and the
    !! files containing the initial geometry and velocity of the system.
    !! At the moment, only the name of the main input file can be passed through the command line
    !! arguments, all other options are set in this file.
    !----------------------------------------------------------------------------------------------
    subroutine read_input()
        use nuclear_dyn_mod
        character(len=1000) :: temp
        integer :: i

        ! Get work directory.
        call getcwd(temp)
        ctrl%maindir = trim(adjustl(temp))//'/'

        ! Read main input file.
        if (stdp1) then
            write(stdout, *)
            write(stdout, '(1x,a,a,a)') 'Reading dynamics input file ', trim(maininp), '.'
        end if
        call read_main()

        ! Start the random number generator.
        call ctrl%rng%init(ctrl%rng%seed)
        if (stdp1) then
            write(stdout, *)
            write(stdout, '(1x,a,i0)') 'Random number generator started with seed: ', ctrl%rng%seed
        end if

        ! Set initial directory.
        if (ctrl%qlib == 0) then
            ctrl%qmdir = 'qmdir'
            call get_environment_variable('QMDIR', temp)
            if (temp /= '') ctrl%qmdir = trim(adjustl(temp))
            if (stdp1) then
                write(stdout, *)
                write(stdout, '(1x,a,a)') "Work directory for QM calculations: ", ctrl%qmdir
            end if
            call system('mkdir -p '//ctrl%qmdir)
        end if

        ! If restarting don't read initial conditions.
        if (ctrl%restart) return

        ! Read initial conditions.
        if (stdp2) then
            write(stdout, *)
            write(stdout, '(1x,a,a,a)') 'Reading geometry file ', geominp, '.'
        end if
        call read_geom()
        select case(veloinp)
        case('maxwell-boltzmann', 'mb')
            if (stdp2) then
                write(stdout, '(1x,a,a)') 'Generating random velocities following ', &
                &                         'Maxwell-Boltzmann distribution.'
            end if
            call maxwell_boltzmann_velo(ctrl%rng, tr1%mass, mb_temperature, tr1%velo)
        case default
            if (stdp2) write(stdout, '(1x,a,a,a)') 'Reading velocity file ', veloinp, '.'
            call read_velo()
        end select

        ! Reorient the molecule based on orientlvl
        select case(ctrl%orientlvl)
        case(1)
            call set_geom_center_of_mass(tr1%mass, tr1%geom)
        case(2)
            call set_geom_center_of_mass(tr1%mass, tr1%geom)
            call project_translation_rotation_from_velocity(tr1%mass, tr1%geom, tr1%velo)
        end select


        ! Read list of partial charges.
        if (ctrl%pcharge) then
            if (stdp2) write(stdout, '(1x,a,a,a)') 'Reading partial charges file ', pcinp, '.'
            call read_pc()
        end if

        ! Print initial information about the full system.
        if (stdp1) then
            write(stdout, *)
            write(stdout, '(1x,a,i0,a)') 'Starting calculation with ', tr1%natom, ' atoms.'
            write(stdout, '(3x,a,i0,a)') 'QM system consists of ', tr1%qnatom, '.'
            write(stdout, '(3x,a,i0,a)') 'MM system consists of ', tr1%mnatom, '.'
        end if

        if (ctrl%qm .and. stdp2) then
            write(stdout, *)
            write(stdout, '(3x,a)') 'Initial geometry of the QM system: '
            do i = 1, tr1%qnatom
                write(stdout, '(2x,a5,1000e18.10)') tr1%sym(tr1%qind(i)), tr1%geom(:, tr1%qind(i))
            end do
            write(stdout, '(3x,a)') 'Initial velocity of the QM system: '
            do i = 1, tr1%qnatom
                write(stdout, '(2x,a5,1000e18.10)') tr1%sym(tr1%qind(i)), tr1%velo(:, tr1%qind(i))
            end do
        end if
        if (ctrl%mm .and. stdp2) then
            write(stdout, *)
            write(stdout, '(3x,a)') 'Initial geometry of the MM system: '
            do i = 1, tr1%mnatom
                write(stdout, '(2x,a5,1000e18.10)') tr1%sym(tr1%mind(i)), tr1%geom(:, tr1%mind(i))
            end do
            write(stdout, '(3x,a)') 'Initial velocity of the MM system: '
            do i = 1, tr1%mnatom
                write(stdout, '(2x,a5,1000e18.10)') tr1%sym(tr1%mind(i)), tr1%velo(:, tr1%mind(i))
            end do
        end if

    end subroutine read_input


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_main
    !
    ! DESCRIPTION:
    !> @brief Read the main program control variables.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine read_main()
        use matrix_mod, only : unit_mat
        logical :: check
        logical :: buffer
        type(reader) :: readf
        integer, allocatable :: uncouple_states(:) !< States for which no couplings are calculated.
        integer :: i
        buffer = .false.

        inquire(file=maininp, exist=check)
        if (.not. check) then
            write(stderr,*) 'Error in Input module, read_main subroutine.'
            write(stderr,*) ' Main input file (', trim(maininp),') not found.'
            stop
        end if

        call readf%open(maininp, comment='#')
        ! Mandatory keywords are read first:
        call read_method(readf)
        call read_system(readf)

        ! Followed by optional keywords.
        call read_dynamics(readf)
        call read_output(readf)
        call read_surfhop(readf, uncouple_states)
        call read_mm(readf)
        call read_constraints(readf, tol_cns)
        call read_thermostat(readf)
        call read_restart(readf)
        call readf%close()

        ! Modify options based on type of calculation.
        if (tr1%max_nstate == 1) then
            ctrl%sh = 0
            ctrl%print(10) = .false.
        end if
        if (ctrl%print(10)) ctrl%oscill = .true.
        select case(ctrl%sh)
        case(0)
            ctrl%print(6:9) = .false.
        case(1)
            ctrl%print(6:9) = .false.
        case(2)
            if (ctrl%tdc_type == 2) then
                ctrl%print(7) = .false.
            end if
        end select
        if (ctrl%thermostat /= 0) then
            ctrl%max_tot_en_change = huge(ctrl%max_tot_en_change)
            ctrl%max_tot_en_change_step = huge(ctrl%max_tot_en_change_step)
        end if

        ! Allocate all arrays of size tr1%max_nstate.
        if (.not. ctrl%restart) then
            allocate(tr1%qe(tr1%max_nstate))
            if (ctrl%oscill) allocate(tr1%qo(tr1%max_nstate - 1))
            select case(ctrl%sh)
            case(1)
                allocate(tr1%prob(tr1%max_nstate), source=0.0_dp)
                allocate(tr1%gap_2deriv(3, tr1%max_nstate), source=0.0_dp)
            case(2, 3)
                allocate(tr1%cwf(tr1%max_nstate))
                allocate(tr1%prob(tr1%max_nstate), source=0.0_dp)
                if (ctrl%phaselvl > 0) allocate(tr1%phase(tr1%max_nstate), source=1)
                select case (ctrl%tdc_type)
                case(1)
                    allocate(tr1%olap(tr1%max_nstate, tr1%max_nstate))
                    tr1%olap = unit_mat(tr1%max_nstate)
                case(2)
                    ! Nonadiabatic coupling vecotrs are allocated after reading number of atoms.
                case(3)
                    allocate(tr1%olap(tr1%max_nstate, tr1%max_nstate))
                    allocate(tr1%adt(tr1%max_nstate, tr1%max_nstate))
                    tr1%olap = unit_mat(tr1%max_nstate)
                    tr1%adt = unit_mat(tr1%max_nstate)
                end select
                tr1%cwf = cmplx((0.0_dp, 0.0_dp), kind = dp)
                tr1%cwf(tr1%cstate) = cmplx((1.0_dp, 0.0_dp), kind = dp)
                if (ctrl%print(8) .and. (.not. allocated(tr1%adt))) then
                    allocate(tr1%adt(tr1%max_nstate, tr1%max_nstate))
                end if
            end select
        end if
        select case(ctrl%sh)
        case(2, 3)
            allocate(ctrl%couple(tr1%max_nstate), source=.true.)
            if (allocated(uncouple_states)) then
                do i = 1, size(uncouple_states, 1)
                    if (uncouple_states(i) > tr1%max_nstate) exit
                    ctrl%couple(uncouple_states(i)) = .false.
                end do
            end if
        end select

    end subroutine read_main


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_geom
    !
    ! DESCRIPTION:
    !> @brief Read the list of symbols, masses and positions of all atoms in the full system.
    !> @details
    !! Read number of atoms in the full system and the QM and MM subsystems. Allocate all arrays
    !! of dimensions natom, qnatom or mnatom. Read the symbol, mass and position arrays for
    !! the full system.
    !! Assign QM and MM pointers to the appropriate sections of the main arrays.
    !----------------------------------------------------------------------------------------------
    subroutine read_geom()
        logical :: check
        character(len=2) :: tsym
        real(dp) :: tmass
        character :: part
        integer :: i
        integer :: d
        type(reader) :: readf

        inquire(file=geominp, exist=check)
        if (.not. check) then
            write(stderr,*) 'Error in Input module, read_geom subroutine.'
            write(stderr,*) ' Geometry file (', trim(geominp),') not found.'
            stop
        end if

        ! First run through the file to check number of atoms in the QM and MM part of the system.
        call readf%open(geominp, abort_on_eof=.false., comment='#')
        do
            call readf%next()
            if (is_iostat_end(readf%iostat)) exit
            call readf%parseline(' ,')
            if (readf%narg < tr1%ndim + 3) then
                write(stderr, '(a)') 'Error in Input module, read_geom subroutine.'
                write(stderr, '(a,i0,a)') ' Need ', 3+tr1%ndim, ' arguments per line:'
                write(stderr, '(a)') '  symbol mass coords q/m'
                write(stderr, '(a,a)') '  Line: ', readf%line
                stop
            end if
            read(readf%args(3+tr1%ndim)%s, *) part
            select case (tolower(part))
            case('q')
                tr1%qnatom = tr1%qnatom + 1
            case('m')
                tr1%mnatom = tr1%mnatom + 1
            case default
                write(stderr, *) 'Error in Input module, read_geom subroutine.'
                write(stderr, *) ' Lines should end with "q" for QM atoms or "m" for MM atoms.'
                stop
            end select
        end do
        call readf%rewind()

        tr1%natom = tr1%qnatom + tr1%mnatom
        if (tr1%natom == 0) then
            write(stderr, *) 'Error in Input module, read_geom subroutine.'
            write(stderr, *) ' No atoms found in geometry file!'
            stop
        end if

        ! Allocate arrays of natom dimensions.
        allocate(tr1%sym(tr1%natom))
        allocate(tr1%mass(tr1%natom))
        allocate(tr1%chrg(tr1%natom))
        allocate(tr1%geom(tr1%ndim, tr1%natom))
        allocate(tr1%velo(tr1%ndim, tr1%natom))
        allocate(tr1%grad(tr1%ndim, tr1%natom))
        if (tr1%qnatom > 0) then
            ctrl%qm = .true.
            allocate(tr1%qind(tr1%qnatom))
        else
            ctrl%qm = .false.
        end if
        if (tr1%mnatom > 0) then
            ctrl%mm = .true.
            ctrl%print(12) = .true.
            allocate(tr1%mind(tr1%mnatom))
        else
            ctrl%mm = .false.
        end if
        if ((ctrl%tdc_type == 2) .or. (ctrl%vrescale == 3)) then
            allocate(tr1%nadv(tr1%ndim*tr1%qnatom, tr1%max_nstate, tr1%max_nstate))
        end if

        ! Second run through the file to read the atom, mass, position and index of each atom.
        tr1%qnatom = 0
        tr1%mnatom = 0
        do i = 1, tr1%natom
            call readf%next()
            call readf%parseline(' ,')

            ! Read symbol.
            read(readf%args(1)%s, *) tsym
            tr1%sym(i) = tolower(tsym)

            ! Read mass.
            read(readf%args(2)%s, *) tmass
            tr1%mass(i) = tmass * Da_me

            ! Read initial geometry.
            do d = 1, tr1%ndim
                read(readf%args(2+d)%s, *) tr1%geom(d, i)
            end do

            ! Check part (q = qm, m = mm)
            read(readf%args(3+tr1%ndim)%s, *) part
            select case (tolower(part))
            case('q')
                tr1%qnatom = tr1%qnatom + 1
                tr1%qind(tr1%qnatom) = i
            case('m')
                tr1%mnatom = tr1%mnatom + 1
                tr1%mind(tr1%mnatom) = i
            end select
        end do
        call readf%close()
    end subroutine read_geom


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_velo
    !
    ! DESCRIPTION:
    !> @brief Read the velocities of the atoms in the system.
    !----------------------------------------------------------------------------------------------
    subroutine read_velo()
        logical :: check
        integer :: i
        integer :: d
        type(reader) :: readf

        inquire(file=veloinp, exist=check)
        if (.not. check) then
            write(stderr, *) 'Error in Input module, read_velo subroutine.'
            write(stderr, *) ' Velocity file (', trim(veloinp),') not found.'
            stop
        end if

        call readf%open(veloinp, comment='#')
        do i = 1, tr1%natom
            call readf%next()
            call readf%parseline(' ,')
            if (readf%narg < tr1%ndim) then
                write(stderr, *) 'Error in Input module, read_velo subroutine.'
                write(stderr, '(a,i0,a)') ' Premature end of line, expect ', tr1%ndim, ' arguments.'
                write(stderr, '(a,a)') ' Line: ', readf%line
                stop
            end if
            do d = 1, tr1%ndim
                read(readf%args(d)%s, *) tr1%velo(d, i)
            end do
        end do
        call readf%close()
    end subroutine read_velo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_pc
    !
    ! DESCRIPTION:
    !> @brief Read the partial charges of the atoms in the system.
    !----------------------------------------------------------------------------------------------
    subroutine read_pc()
        logical :: check
        integer :: i
        type(reader) :: readf

        inquire(file=pcinp, exist=check)
        if (.not. check) then
            write(stderr, *) 'Error in Input module, read_pc subroutine.'
            write(stderr, *) ' Point charges file (', trim(pcinp),') not found.'
            stop
        end if

        call readf%open(pcinp)
        do i = 1, tr1%natom
            call readf%next()
            read(readf%line, *) tr1%chrg(i)
        end do
        call readf%close()
    end subroutine read_pc


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_method
    !
    ! DESCRIPTION:
    !> @brief Read the $method section of the main input file.
    !> @details
    !! In this section, programs to run the QM and MM calculations are selected.
    !----------------------------------------------------------------------------------------------
    subroutine read_method(readf)
        use model_mod, only : qmodel
        type(reader), intent(inout) :: readf

        if (.not. allocated(ctrl%qprog)) ctrl%qprog = ''
        if (.not. allocated(ctrl%mprog)) ctrl%mprog = ''

        if (stdp3) write(stdout, '(3x,a)') 'Reading $method section of the dynamics input file.'
        call readf%rewind()
        call readf%go_to_keyword('$method')
        do
            call readf%next()
            if (index(readf%line, '$') == 1) exit
            call readf%parseline(' =')
            select case(readf%args(1)%s)
            case('qm')
                ctrl%qprog = readf%args(2)%s
                if (stdp3) then
                    write(stdout, '(5x,a)') '"'//ctrl%qprog//&
                    &                       '" will be used for QM calculations.'
                end if
            case('mm')
                ctrl%mprog = readf%args(2)%s
                if (stdp3) then
                    write(stdout, '(5x,a)') '"'//ctrl%mprog//&
                    &                       '" will be used for QM calculations.'
                end if
            case('overlap')
                ctrl%oprog = readf%args(2)%s
                if (stdp3) then
                    write(stdout, '(5x,a)') '"'//ctrl%oprog//&
                    &                       '" will be used for overlap calculations.'
                end if
            case('qlib')
                select case(readf%args(2)%s)
                case('quantics')
                    ctrl%qlib = 1
                case('model')
                    ctrl%qlib = 2
                    call qmodel%init(readf%args(3:))
                case default
                    write(stderr, *) 'Error in Input module, read_method subroutine.'
                    write(stderr, *) '  Unrecognized qlib keyword: ', readf%args(2)%s
                    stop
                end select
            case('noise')
                read(readf%args(2)%s, *) ctrl%noise
            case('seed')
                read(readf%args(2)%s, *) ctrl%rng%seed
            case('qm_en_error')
                read(readf%args(2)%s, *) ctrl%qm_en_err
            end select
        end do

        if ((ctrl%qlib > 0) .and. (ctrl%qprog /= '')) then
            write(stderr, *) 'Error in Input module, read_method subroutine.'
            write(stderr, *) '  Both qlib and qm keywords given in input.'
            stop
        end if
    end subroutine read_method


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_system
    !
    ! DESCRIPTION:
    !> @brief Read information about the system.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine read_system(readf)
        type(reader), intent(inout) :: readf

        call readf%rewind()
        call readf%go_to_keyword('$system')
        do
            call readf%next()
            if (index(readf%line, '$') == 1) exit
            call readf%parseline(' =')
            select case(readf%args(1)%s)
            case('nstate')
                read(readf%args(2)%s, *) tr1%nstate
            case('variable_nstate')
                read(readf%args(2)%s, *) ctrl%variable_nstate
            case('max_nstate')
                read(readf%args(2)%s, *) tr1%max_nstate
            case('min_nstate')
                read(readf%args(2)%s, *) tr1%min_nstate
            case('istate')
                read(readf%args(2)%s, *) tr1%cstate
            case('geometry')
                geominp = readf%args(2)%s
            case('velocity')
                veloinp = readf%args(2)%s
                if (veloinp == 'maxwell-boltzmann' .or. veloinp == 'mb') then
                    if (readf%narg < 3) then
                        write(stderr, *) 'Error in Input module, read_system subroutine.'
                        write(stderr, '(2x,a,a)') 'Maxwell-Boltzmann velocities requested, ', &
                                                  'but no temperature given'
                        call abort()
                    end if
                    read(readf%args(3)%s, *) mb_temperature
                end if
            case('ndim')
                read(readf%args(2)%s, *) tr1%ndim
            case default
                write(stderr, *) 'Warning in Input module, read_system subroutine.'
                write(stderr, *) '  Skipping line with unrecognized keyword:'
                write(stderr, *) readf%line
            end select
        end do

        select case (ctrl%variable_nstate)
        case(0)
            if (tr1%nstate == 0) then
                write(stderr, *) 'Error in Input module, read_system subroutine.'
                write(stderr, *) '  Number of states not defined.'
                stop
            end if
            if ((tr1%max_nstate /= 0) .or. (tr1%min_nstate /= 0)) then
                write(stderr, *) 'Warning in Input module, read_system subroutine.'
                write(stderr, *) '  Options max_nstate and min_nstate are ignored when '//&
                                'variable_nstate = 0.'
            end if
            tr1%max_nstate = tr1%nstate
            tr1%min_nstate = tr1%nstate
        case(1)
            if ((tr1%max_nstate == 0) .or. (tr1%min_nstate == 0)) then
                write(stderr, *) 'Error in Input module, read_system subroutine.'
                write(stderr, *) '  Maximum/minimum number of states not defined.'
                stop
            end if
            if (tr1%nstate == 0) then
                write(stderr, *) 'Warning in Input module, read_system subroutine.'
                write(stderr, *) '  Option nstate is ignored when variable_nstate = 1'
            end if
            tr1%nstate = tr1%max_nstate
        end select
        if (tr1%cstate == 0) then
            write(stderr, *) 'Warning in Input module, read_system subroutine.'
            write(stderr, '(3x,a,i0,a)') 'Initial state not defined, setting initial state = ', &
            &                            tr1%nstate, '.'
            tr1%cstate = tr1%nstate
        end if

    end subroutine read_system


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_output
    !
    ! DESCRIPTION:
    !> @brief Read output options of the program.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine read_output(readf)
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
            case('printerval')
                read(readf%args(2)%s, *) ctrl%printerval
            case('buinterval')
                read(readf%args(2)%s, *) ctrl%buinterval
            case default
                write(stderr, *) 'Warning in Input module, read_output subroutine.'
                write(stderr, *) '  Skipping line with unrecognized keyword:'
                write(stderr, *) readf%line
            end select
        end do
    end subroutine read_output


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_dynamics
    !
    ! DESCRIPTION:
    !> @brief Read information about the nuclear dynamics.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine read_dynamics(readf)
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
                write(stderr, *) 'Warning in Input module, read_dynamics subroutine.'
                write(stderr, *) '  Skipping line with unrecognized keyword:'
                write(stderr, *) readf%line
            end select
        end do
    end subroutine read_dynamics


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_constraints
    !
    ! DESCRIPTION:
    !> @brief Read information about constraints.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine read_constraints(readf, default_tol)
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
                write(stderr, *) 'Error in Input module, read_constraints subroutine.'
                write(stderr, *) '  Unrecognized constraint type: ', readf%args(1)%s
                stop
            end select
        end do
    end subroutine read_constraints


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_surfhop
    !
    ! DESCRIPTION:
    !> @brief Read the $surfhop section of the main input file.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine read_surfhop(readf, uncouple_states)
        type(reader), intent(inout) :: readf
        integer, allocatable, intent(out) :: uncouple_states(:)
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
                case('edc')
                    ctrl%decohlvl = 1
                case default
                    write(stderr, *) 'Error in Input module, read_surfhop subroutine.'
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
                if (readf%narg > 1) then
                    select case(readf%args(2)%s)
                    case('constant')
                        ctrl%tdc_interpolate = 1
                    case('linear')
                        ctrl%tdc_interpolate = 2
                    end select
                end if
            case('adt')
                ctrl%tdc_type = 3
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
                    write(stderr, *) 'Error in Input module, read_surfhop subroutine.'
                    write(stderr, '(2x,a)') '  Unrecognized phase matching algorithm: '
                    write(stderr, '(2x,a)') readf%line
                    stop
                end select
            case('lz_prob_conv')
                read(readf%args(2)%s, *) ctrl%lz_prob_conv
            case('lz_min_dt')
                read(readf%args(2)%s, *) ctrl%lz_min_dt
                ctrl%lz_min_dt = ctrl%lz_min_dt / aut_fs
            case('stop_s0s1_ci')
                read(readf%args(2)%s, *) ctrl%stop_s0s1_ci
                ctrl%stop_s0s1_ci = ctrl%stop_s0s1_ci / eh_eV
            case('uncouple_state')
                call compact(readf%line)
                call read_index_list(readf%line(15:), uncouple_states)
            case default
                write(stderr, *) 'Warning in Input module, read_surfhop subroutine.'
                write(stderr, *) '  Skipping line with unrecognized keyword:'
                write(stderr, *) readf%line
            end select
        end do
    end subroutine read_surfhop

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_mm
    !
    ! DESCRIPTION:
    !> @brief Read the $mm section of the main input file.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine read_mm(readf)
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
                    read (readf%line, *) tr1%pbcbox
                end if
            case default
                write(stderr, *) 'Warning in Input module, read_mm subroutine.'
                write(stderr, *) '  Skipping line with unrecognized keyword:'
                write(stderr, *) readf%line
            end select
        end do
    end subroutine read_mm


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_thermostat
    !
    ! DESCRIPTION:
    !> @brief Read the $thermostat section of the main input file.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine read_thermostat(readf)
        type(reader), intent(inout) :: readf
        logical :: check

        call readf%rewind()
        call readf%go_to_keyword('$thermostat', found=check)
        if (.not. check) return

        call readf%parseline(' =')
        if (readf%narg > 1) then
            select case(readf%args(2)%s)
            case('off')
                ctrl%thermostat = 0
                return
            case('berendsen')
                ctrl%thermostat = 1
            case default
                write(stderr, *) 'Warning in Input module, read_thermostat subroutine.'
                write(stderr, *) '  Unrecognized thermostat option:'
                write(stderr, *) readf%line
            end select
        end if


        do
            call readf%next()
            if (index(readf%line, '$') == 1) exit
            call readf%parseline(' =')
            select case(readf%args(1)%s)
            case('target_t')
                read(readf%args(2)%s, *) ctrl%target_t
            case('tau_t')
                if (ctrl%thermostat /= 1) then
                    write(stderr, *) 'Warning in Input module, read_thermostat subroutine.'
                    write(stderr, *) '  Skipping tau_t, only relevant for Berendsen thermostat.'
                end if
                read(readf%args(2)%s, *) ctrl%tau_t
                ctrl%tau_t = ctrl%tau_t / aut_fs
            case default
                write(stderr, *) 'Warning in Input module, read_thermostat subroutine.'
                write(stderr, *) '  Skipping line with unrecognized keyword:'
                write(stderr, *) readf%line
            end select
        end do
    end subroutine read_thermostat


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_restart
    !
    ! DESCRIPTION:
    !> @brief Check if the $restart keyword is present in the main input file.
    !> @details
    !! See manual for details concerning the input.
    !----------------------------------------------------------------------------------------------
    subroutine read_restart(readf)
        type(reader), intent(inout) :: readf
        logical :: check

        call readf%rewind()
        call readf%go_to_keyword('$restart', found=check)
        if (check) ctrl%restart = .true.
    end subroutine read_restart


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: maxwell_boltzmann_velo
    !
    ! DESCRIPTION:
    !> @brief Generate random velocities based on Maxwell Boltzmann distribution
    !----------------------------------------------------------------------------------------------
    subroutine maxwell_boltzmann_velo(rng, mass, temperature, veloc)
        use random_mod, only : rng_type
        class(rng_type), allocatable, intent(in) :: rng
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: temperature
        real(dp), intent(out) :: veloc(:, :)
        real(dp) :: sigma
        integer :: i

        do i = 1, size(veloc, 2)
            sigma = sqrt(num2 * au_boltzmann * temperature / mass(i))
            call rng%gaussian(0.0_dp, sigma, veloc(:, i))
        end do
    end subroutine maxwell_boltzmann_velo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: print_help
    !> @brief Print help message and exit program.
    !----------------------------------------------------------------------------------------------
    subroutine print_help()
        write(stdout, '(a)') 'usage: zaghop [options] input-file'
        write(stdout, '(a)')
        write(stdout, '(a)') 'Run MQC trajectory calculations.'
        write(stdout, '(a)')
        write(stdout, '(a)') 'positional arguments:'
        write(stdout, '(a)') '  input-file                  name of input file                                   '
        write(stdout, '(32x,a,a)') 'default: ', maininp
        write(stdout, '(a)')
        write(stdout, '(a)') 'optional arguments:'
        write(stdout, '(a)') '  -h, --help                  show this message and exit                      '
        write(stdout, '(a)') '  -p, --print-level p         control output level of program (0 = quiet)          '
        write(stdout, '(32x,a,i0)') 'default: ', count([stdp1, stdp2, stdp3, stdp4])
        stop
    end subroutine print_help


end module input_mod
