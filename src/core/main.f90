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
    use system_type_mod
    use tsh_method_mod, only : tsh_method, tsh_request, request_goto, request_insert

    ! Import subroutines
    use timing_mod, only : timer
    use input_mod, only : read_input, command_line_interface, set_defaults
    use interface_mod
    use file_mod, only : check_is_dir
    use nuclear_dyn_mod
    implicit none

    logical :: check
    logical :: insert
    integer :: abort_flag = 0
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
        write(stdout, *) '                               ZagHop, v0.94                                   '
        write(stdout, *) 'Program compiled '//__DATE__//' at '//__TIME__//'.'
        write(stdout, *) '-------------------------------------------------------------------------------'
    end if

    ! Read input file(s)
    call read_input()

    ! Restart old trajectory or open new output file.
    if (ctrl%restart) then
        ! Restart old trajectory.
      !  call trajectory_read_backup(ctrl%bufile, trajectory_data)
        call write_status(abort_flag)
    else
        ! Create output directory.
        if (check_is_dir(ctrl%output_dir)) then
            write(stderr,*) 'Error. Results directory already exists.'
            write(stderr,*) 'Remove or rename the directory and restart the program.'
            stop
        else
            call system('mkdir -p '//ctrl%output_dir)
            call tr1%open_files(ctrl%print, ctrl%print_units, ctrl%output_dir)
            call tr1%writeheader(ctrl%print, ctrl%print_units)
            call write_status(abort_flag)
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
        call tr1%wf%start(tr1%time, tr1%pot)
        if (stdp1) call stepclock%print(stdout, '  QM run time:')
        call trajectory_write_final(tr1%wf%request%t_final, ctrl%print, ctrl%print_units,         &
        &                           ctrl%printerval)
        call trajectory_next(ctrl%dt_0)
    end if
    ctrl%t0_tot_en = tr1%tote()

    ! Stop the program if max_time already reached.
    if (tr1%time > ctrl%max_time) then
        if (stdp1) write(stdout, *) 'No steps to be done, ending calculation.'
        if (stdp1) call mainclock%print(stdout, ' Total run time:')
        call write_status(1)
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
        call dyn_updategeom(tr1%dt, tr1%mass, tr2%geom, tr2%grad, tr2%velo, tr1%geom)

        ! Get new gradients.
        if (ctrl%mm) call run_mm(tr1)
        !> @todo Re-check and re-activate variable_nstate option.
        ! if (ctrl%variable_nstate == 1) then
        !     do i = tr1%nstate, max(tr1%cstate + 2, -1, tr1%min_nstate), -1
        !         ! Remove top states whose contribution to the total wave function is below the
        !         ! threshold until the first state which should be kept is reached.
        !         if (abs(tr1%cwf(i)) > 0.01) exit
        !         if (stdp2) then
        !             write(stdout, '(3x,a, i0, a)') 'Excluding state ', i, ' from calculation.'
        !         end if
        !         tr1%nstate = i-1
        !         tr1%cwf(i) = 0.0_dp
        !     end do
        ! end if
        if (stdp2) write(stdout, *) '  Running QM calculation.'
        call run_qm(tr1, .false.)

        ! Get new velocity.
        call dyn_updatevelo(tr1%dt, tr1%mass, tr1%geom, tr2%grad, tr1%grad, tr2%velo,         &
        &                   tr1%velo)

        ! Propagate the surface hopping method and carry out its request.
        call hopping(insert)

        ! Checks and output wait until the method stops inserting points.
        if (.not. insert) then
            call end_of_step_checks()

            ! Write output.
            if (abort_flag /= 0) then
                call trajectory_write_final(huge(1.0_dp), ctrl%print, ctrl%print_units,           &
                &                           ctrl%printerval)
            else
                call trajectory_write_final(tr1%wf%request%t_final, ctrl%print, ctrl%print_units, &
                &                           ctrl%printerval)
            end if
            if ((mod(tr1%step, ctrl%buinterval) == 0) .or. (abort_flag /= 0)) then
                call trajectory_write_backup(ctrl%bufile, trajectory_data)
            end if
        end if

        ! Prepare next step.
        if (insert) then
            call trajectory_next(ctrl%dt_0, tr1%wf%request%dt)
        else
            call trajectory_next(ctrl%dt_0)
        end if
        if (stdp2) call stepclock%print(stdout, '   Step run time:')
        if (abort_flag /= 0) then
            if (stdp1) write(stdout, *) '  Ending calculation.'
            if (stdp1) write(stdout, *)
            call write_status(abort_flag)
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


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: write_status
    !
    ! DESCRIPTION:
    !> @brief Write the abort flag to the status file in the results directory.
    !----------------------------------------------------------------------------------------------
    subroutine write_status(flag)
        integer, intent(in) :: flag
        integer :: sunit

        open(newunit=sunit, file=ctrl%output_dir//'/status', action='write', status='replace')
        write(sunit, '(i0)') flag
        close(sunit)
    end subroutine write_status


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: hopping
    !
    ! DESCRIPTION:
    !> @brief Propagate the surface hopping method to the current point and carry out its request.
    !> @details
    !! If the request makes a stored point current, the method object is moved to that point, so
    !! that it keeps its own history. The previous current point can stay stored before the new
    !! one, so it keeps a copy of the method object with its active state from before the hop.
    !! After a hop, the gradient of the new state is calculated and the velocity is rescaled. A
    !! frustrated hop is undone by the method.
    !----------------------------------------------------------------------------------------------
    subroutine hopping(insert)
        logical, intent(out) :: insert !< The next step inserts a point before later stored points.
        class(tsh_method), allocatable :: wf
        type(tsh_request) :: request
        logical :: frustrated

        call tr1%wf%decohere(tr1%qkine(), tr1%dt, tr1%pot)
        call tr1%wf%propagate(tr1%time, tr1%pot)
        request = tr1%wf%request

        select case(request%action)
        case(request_goto, request_insert)
            allocate(wf, source=tr1%wf)
            if (request%hop) tr1%wf%active_state = request%previous_state
            call trajectory_goto(request%time, keep_later=(request%action == request_insert))
            call move_alloc(wf, tr1%wf)
        end select
        insert = (request%action == request_insert)

        if (request%hop) then
            if (stdp2) then
                write(stdout, '(3x,a)') 'Change of state: '
                write(stdout, '(5x,a,i0)') 'Previous state: ', request%previous_state
                write(stdout, '(5x,a,i0)') 'Current state: ', tr1%wf%active_state
                write(stdout, '(5x,a)') 'Running QM gradient calculation for new state.'
            end if
            call run_qm(tr1, .true.)
            call tr1%wf%rescale_velocity(ctrl%vrescale, ctrl%fhop, tr1%qind,                      &
            &                            request%previous_state, tr1%pot, tr1%mass, tr1%velo,    &
            &                            frustrated)
            if (frustrated) call tr1%wf%hop_rejected()
            tr1%grad(:, tr1%qind) = tr1%pot%get_gradient(tr1%wf%basis, tr1%wf%active_state)
        end if
    end subroutine hopping


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: end_of_step_checks
    !
    ! DESCRIPTION:
    !> @brief Set abort_flag if the trajectory should end after the current point.
    !----------------------------------------------------------------------------------------------
    subroutine end_of_step_checks()
        ! Stop the program after max_time was reached.
        if (tr1%time + tinydp >= ctrl%max_time) then
            if (ctrl%target_state == -2) then
                if (stdp1) write(stdout, *) '  Target state max_t reached.'
                abort_flag = 3
            else
                if (stdp1) write(stdout, *) '  Reached max_time.'
                abort_flag = 1
            end if
        end if

        ! Stop the program if the total energy changed above tolerance levels.
        if (abs(tr1%tote() - tr2%tote()) > ctrl%max_tot_en_change_step) then
            if (stdp1) write(stdout, *) '  Change in total energy too large.'
            abort_flag = -1
        else if (abs(tr1%tote() - ctrl%t0_tot_en) > ctrl%max_tot_en_change) then
            if (stdp1) write(stdout, *) '  Drift in total energy too large.'
            abort_flag = -2
        end if

        ! Stop the program at S0/S1 conical intersection.
        if (tr1%wf%n_state > 1) then
            if (tr1%pot%get_energy(tr1%wf%basis, 2) - tr1%pot%get_energy(tr1%wf%basis, 1)        &
            &   < ctrl%stop_s0s1_ci) then
                if (stdp1) write(stdout, *) '  Intersection with ground state.'
                abort_flag = 4
            end if
        end if

        ! Stop the program after reaching the target state.
        if (tr1%wf%active_state == ctrl%target_state) then
            if (ctrl%target_state_time == 0.0_dp) then
                abort_flag = 2
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
            abort_flag = -3
            call system('rm dynamics.stop')
        end if
    end subroutine end_of_step_checks


end program zaghop
