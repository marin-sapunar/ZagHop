!--------------------------------------------------------------------------------------------------
!> @mainpage Dynamics
!> Trajectory surface-hopping nonadiabatic dynamics program
!--------------------------------------------------------------------------------------------------
! PROGRAM: zaghop
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2016
!
! DESCRIPTION:
!> @brief Trajectory surface-hopping nonadiabatic dynamics program.
!> @details
!! Main program file. Calls input and output subroutines and contains the main dynamics loop.
!! Interfaces to QM and MM programs are called from within the loop to calculate energies,
!! gradients and couplings of the system. Once the calculations are performed, the electronic
!! wavefunction and nuclear coordinates are propagated.
!--------------------------------------------------------------------------------------------------
program zaghop
    ! Import variables
    use global_defs
    use constants
    use control_var
    use system_var

    ! Import subroutines
    use timing_mod, only : timer
    use input_mod, only : read_input, command_line_interface, set_defaults
    use interface_mod
    use file_mod, only : check_is_dir
    use hopping_mod
    use nuclear_dyn_mod
    implicit none

    logical :: check
    logical :: abort_flag = .false.
    type(timer) :: mainclock
    type(timer) :: stepclock
    real(dp), allocatable :: hop_grad(:, :)
    real(dp) :: tinydp = 1.0e-8_dp
    integer :: i

    ! Allocate main trajectory data array.
    call allocate_trajectory_data()
    ! Set default values for most options and read command line arguments.
    call set_defaults()
    call command_line_interface()

    if (stdp1) then
        call mainclock%start()
        write(stdout, *) '-------------------------------------------------------------------------------'
        write(stdout, *) '                               ZagHop, v0.93                                   '
        write(stdout, *) 'Program compiled '//__DATE__//' at '//__TIME__//'.'
        write(stdout, *) '-------------------------------------------------------------------------------'
    end if

    ! Read input file(s)
    call read_input()

    ! Restart old trajectory or open new output file.
    if (ctrl%restart) then
        ! Restart old trajectory.
        call trajectory_read_backup(ctrl%bufile, trajectory_data)
    else
        ! Create output directory.
        if (check_is_dir(ctrl%output_dir)) then
            write(stderr,*) 'Warning. Results directory already exists.'
            write(stderr,*) '  New results will be appended.'
            call tr1%open_files(ctrl%print, ctrl%print_units, ctrl%output_dir)
        else
            call system('mkdir -p '//ctrl%output_dir)
            call tr1%open_files(ctrl%print, ctrl%print_units, ctrl%output_dir)
            call tr1%writeheader(ctrl%print, ctrl%print_units, ctrl%output_dir)
        end if

        ! Run energy/gradient calculation for initial geometry.
        if (ctrl%mm) then
            if (stdp1) write(stdout, '(a)') ' Running initial MM calculation: '
            if (stdp1) call stepclock%start()
            call run_mm(tr1)
            if (stdp1) call stepclock%print(stdout, '  MM run time:')
        end if
        if (stdp1) write(stdout, '(a)') ' Running initial QM calculation: '
        if (stdp1) call stepclock%start()
        call run_qm(tr1, .false.)
        if (stdp1) call stepclock%print(stdout, '  QM run time:')
        call tr1%writestep(ctrl%print, ctrl%print_units, ctrl%output_dir)
        call trajectory_next(ctrl%dt)
    end if
    ctrl%t0_tot_en = tr1%tote()

    ! Stop the program if max_time already reached.
    if (tr1%time > ctrl%max_time) then
        if (stdp1) write(stdout, *) 'No steps to be done, ending calculation.'
        if (stdp1) call mainclock%print(stdout, ' Total run time:')
        stop
    end if

    if (stdp1) then
        write(stdout, *)
        write(stdout, *) '-------------------------------------------------------------------------------'
        write(stdout, *) '                          STARTING MAIN PROGRAM LOOP                           '
        write(stdout, *)
    end if

    main: do
        if (stdp2) call stepclock%start()
        if (stdp1) write(stdout, '(1x,a,i0)') 'Starting step: ', tr1%step

        ! Get new geometry.
        call dyn_updategeom(ctrl%dt, tr1%mass, tr2%geom, tr2%grad, tr2%velo, tr1%geom)

        ! Get new gradients.
        if (ctrl%mm) call run_mm(tr1)
        if (ctrl%variable_nstate == 1) then
            do i = tr1%nstate, max(tr1%cstate + 2, -1, tr1%min_nstate), -1
                ! Remove top states whose contribution to the total wave function is below the
                ! threshold until the first state which should be kept is reached.
                if (abs(tr1%cwf(i)) > 0.01) exit
                if (stdp2) then
                    write(stdout, '(3x,a, i0, a)') 'Excluding state ', i, ' from calculation.'
                end if
                tr1%nstate = i-1
                tr1%cwf(i) = 0.0_dp
            end do
        end if
        if (stdp2) write(stdout, *) '  Running QM calculation.'
        call run_qm(tr1, .false.)

        ! Get new velocity.
        call dyn_updatevelo(ctrl%dt, tr1%mass, tr1%geom, tr2%grad, tr1%grad, tr2%velo,        &
        &                   tr1%velo)

        ! Integrate TDSE.
        call hopping()
        if (ctrl%hop) then
            if (stdp2) then
                write(stdout, '(3x,a)') 'Change of state: '
                write(stdout, '(5x,a,i0)') 'Previous state: ', tr2%cstate
                write(stdout, '(5x,a,i0)') 'Current state: ', tr1%cstate
                write(stdout, '(5x,a)') 'Running QM gradient calculation for new state.'
            end if
            allocate(hop_grad, source=tr1%grad)
            call run_qm(tr1, .true.)
            call sh_rescalevelo(ctrl%vrescale, ctrl%fhop, tr1%qind, tr2%cstate, tr1%cstate,     &
            &                   tr1%mass, tr1%qe, hop_grad, tr1%grad, tr1%nadv, tr1%velo)
            deallocate(hop_grad)
        end if

        ! Stop the program after max_time was reached. Add tinydp to time for precision.
        if (tr1%time + tinydp >= ctrl%max_time) then
            if (ctrl%target_state == -2) then
                if (stdp1) write(stdout, *) '  Target state max_t reached.'
            else
                if (stdp1) write(stdout, *) '  Reached max_time.'
            end if
            abort_flag = .true.
        end if

        ! Stop the program if the total energy changed above tolerance levels.
        if (abs(tr1%tote() - tr2%tote()) > ctrl%max_tot_en_change_step) then
            if (stdp1) write(stdout, *) '  Change in total energy too large.'
            abort_flag = .true.
        else if (abs(tr1%tote() - ctrl%t0_tot_en) > ctrl%max_tot_en_change) then
            if (stdp1) write(stdout, *) '  Drift in total energy too large.'
            abort_flag = .true.
        end if

        ! Stop the program at S0/S1 conical intersection.
        if (tr1%nstate > 1) then
            if (tr1%qe(2) - tr1%qe(1) < ctrl%stop_s0s1_ci) then
                if (stdp1) write(stdout, *) '  Intersection with ground state.'
                abort_flag = .true.
            end if
        end if

        ! Stop the program after reaching the target state.
        if (tr1%cstate == ctrl%target_state) then
            if (ctrl%target_state_time == 0.0_dp) then
                abort_flag = .true.
                if (stdp1) write(stdout, *) '  Target state reached.'
            else
                if (stdp1) write(stdout, *) '  Target state reached, modifying max_time.'
                if (tr1%time + ctrl%target_state_time < ctrl%max_time) then
                    ctrl%max_time = tr1%time + ctrl%target_state_time
                    ctrl%target_state = -2 ! Signal max_time was changed because of target state.
                else
                    ctrl%target_state = -1 ! Ignore target_state, trajectory will end anyways.
                end if
            end if
        end if

        ! The program can be stopped after completing a step by creating a 'dynamics.stop' file
        ! in the working directory.
        inquire(file='dynamics.stop', exist=check)
        if (check) then
            if (stdp1) write(stdout, *) '  Stop file found in working directory.'
            abort_flag = .true.
            call system('rm dynamics.stop')
        end if

        ! Write output and prepare next step.
        if ((mod(tr1%step, ctrl%printerval) == 0)) then
            call tr1%writestep(ctrl%print, ctrl%print_units, ctrl%output_dir)
        end if
        if ((mod(tr1%step, ctrl%buinterval) == 0) .or. (abort_flag)) then
            call trajectory_write_backup(ctrl%bufile, trajectory_data)
        end if

        call trajectory_next(ctrl%dt)
        if (stdp2) call stepclock%print(stdout, '   Step run time:')
        if (abort_flag) then
            if (stdp1) write(stdout, *) '  Ending calculation.'
            if (stdp1) write(stdout, *)
            exit main
        end if

    end do main

    if (stdp1) then
        write(stdout, *)
        write(stdout, *) '                         FINISHING MAIN PROGRAM LOOP                           '
        write(stdout, *) '-------------------------------------------------------------------------------'
        write(stdout, *)
        call mainclock%print(stdout, ' Total run time:')
    end if



end program zaghop
