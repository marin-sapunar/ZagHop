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
    use control_var
    use system_var
    use constants

    ! Import subroutines
    use timing_mod, only : timer
    use input_mod, only : read_input
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

    call mainclock%start()

    write(stdout, *) '-------------------------------------------------------------------------------'
    write(stdout, *) '                               ZagHop, v0.92                                   '
    write(stdout, *) 'Program compiled '//__DATE__//' at '//__TIME__//'.'
    write(stdout, *) '-------------------------------------------------------------------------------'

    ! Read input file(s)
    allocate(t(memory)) ! Allocate number of previous steps to keep in memory.
    call read_input()

    ! Restart old trajectory or open new output file.
    if (ctrl%restart) then
        ! Restart old trajectory.
        call trajectory_read_backup(ctrl%bufile, t)
      ! if (t(1)%substep /= -1) ctrl%dt = ctrl%dt_0 / ctrl%lz_substep
    else
        ! Create output directory.
        if (check_is_dir(ctrl%output_dir)) then
            write(stderr,*) 'Warning. Results directory already exists.'
            write(stderr,*) '  New results will be appended.'
        else
            call system('mkdir -p '//ctrl%output_dir)
            call t(1)%writeheader(ctrl%print, ctrl%output_dir)
        end if

        ! Run energy/gradient calculation for initial geometry.
        if (ctrl%mm) then
            write(stdout, '(a)') ' Running initial MM calculation: '
            call stepclock%start()
            call run_mm(t(1))
            call stepclock%print(stdout, '  MM run time:')
        end if
        write(stdout, '(a)') ' Running initial QM calculation: '
        call stepclock%start()
        call run_qm(t(1), .false.)
        call stepclock%print(stdout, '  QM run time:')
        call t(1)%writestep(ctrl%print, ctrl%output_dir)
        call trajectory_next(1, t, ctrl%dt, .true.)
    end if
    ctrl%t0_tot_en = t(1)%tote()

    ! Stop the program if max_time already reached.
    if (t(1)%time > ctrl%max_time) then
        write(stdout, *) ' NO STEPS TO BE DONE, ending calculation.'
        call mainclock%print(stdout, ' Total run time:')
        stop
    end if

    write(stdout, *)
    write(stdout, *) '-------------------------------------------------------------------------------'
    write(stdout, *) '                          STARTING MAIN PROGRAM LOOP                           '
    write(stdout, *)

    main: do
        call stepclock%start()
        write(stdout, '(1x,a,i0)') 'Starting step: ', t(1)%step

        ! Get new geometry.
        call dyn_updategeom(ctrl%dt, t(1)%mass, t(2)%geom, t(2)%grad, t(2)%velo, t(1)%geom)

        ! Get new gradients.
        if (ctrl%mm) call run_mm(t(1))
        if (ctrl%variable_nstate == 1) then
            do i = t(1)%nstate, max(t(1)%cstate + 2, -1, t(1)%min_nstate), -1
                ! Remove top states whose contribution to the total wave function is below the
                ! threshold until the first state which should be kept is reached.
                if (abs(t(1)%cwf(i)) > 0.01) exit
                write(stdout, '(3x,a, i0, a)') 'Excluding state ', i, ' from calculation.'
                t(1)%nstate = i-1
                t(1)%cwf(i) = 0.0_dp
            end do
        end if
        call run_qm(t(1), .false.)

        ! Get new velocity.
        call dyn_updatevelo(ctrl%dt, t(1)%mass, t(1)%geom, t(2)%grad, t(1)%grad, t(2)%velo,        &
        &                   t(1)%velo)

        ! Integrate TDSE.
        call hopping()
        if (ctrl%hop) then
            write(stdout, '(3x,a)') 'Change of state: '
            write(stdout, '(5x,a,i0)') 'Previous state: ', t(2)%cstate
            write(stdout, '(5x,a,i0)') 'Current state: ', t(1)%cstate
            write(stdout, '(5x,a)') 'Running QM gradient calculation for new state.'
            allocate(hop_grad, source=t(1)%grad)
            call run_qm(t(1), .true.)
            write(stdout, '(5x,a)') 'Ensuring energy conservation.'
            call sh_rescalevelo(ctrl%vrescale, ctrl%fhop, t(1)%qind, t(2)%cstate, t(1)%cstate,     &
            &                   t(1)%mass, t(1)%qe, hop_grad, t(1)%grad, t(1)%nadv, t(1)%velo)
            deallocate(hop_grad)
        end if

        ! Stop the program after max_time was reached. Add tinydp to time for precision.
        if (t(1)%time + tinydp >= ctrl%max_time) then
            if (ctrl%target_state == -2) then
                write(stdout, *) 'Target state max_t reached, ending calculation.'
            else
                write(stdout, *) 'Reached max_time, ending calculation.'
            end if
            abort_flag = .true.
        end if

        ! Stop the program if the total energy changed above tolerance levels.
        if (abs(t(1)%tote() - t(2)%tote()) > ctrl%max_tot_en_change_step) then
            write(stdout, *) 'Change in total energy too large, ending calculation.'
            write(stdout, *)
            abort_flag = .true.
        else if (abs(t(1)%tote() - ctrl%t0_tot_en) > ctrl%max_tot_en_change) then
            write(stdout, *) 'Drift in total energy too large, ending calculation.'
            write(stdout, *)
            abort_flag = .true.
        end if

        ! Stop the program at S0/S1 conical intersection.
        if (t(1)%nstate > 1) then
            if (t(1)%qe(2) - t(1)%qe(1) < ctrl%stop_s0s1_ci) then
                write(stdout, *) 'Intersection with ground state detected, ending calculation.'
                write(stdout, *)
                abort_flag = .true.
            end if
        end if

        ! Stop the program after reaching the target state.
        if (t(1)%cstate == ctrl%target_state) then
            if (ctrl%target_state_time == 0.0_dp) then
                abort_flag = .true.
                write(stdout, *) 'Target state reached, ending calculation.'
                write(stdout, *)
            else
                write(stdout, '(1x,a,f10.5,a)') 'Target state reached, ending calculation after',   &
                &                 ctrl%target_state_time * aut_fs, 'fs.'
                write(stdout, *)
                if (t(1)%time + ctrl%target_state_time < ctrl%max_time) then
                    ctrl%max_time = t(1)%time + ctrl%target_state_time
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
            write(stdout, *) 'Stop file found in working directory, ending calculation.'
            write(stdout, *)
            abort_flag = .true.
            call system('rm dynamics.stop')
        end if

        ! Write output and prepare next step.
        if ((mod(t(1)%step, ctrl%printerval) == 0)) then
            call t(1)%writestep(ctrl%print, ctrl%output_dir)
        end if
        if ((mod(t(1)%step, ctrl%buinterval) == 0) .or. (abort_flag)) then
            call trajectory_write_backup(ctrl%bufile, t)
        end if

        call trajectory_next(1, t, ctrl%dt, .true.)
        call stepclock%print(stdout, '   Step run time:')
        if (abort_flag) exit main

    end do main


    write(stdout, *)
    write(stdout, *) '                         FINISHING MAIN PROGRAM LOOP                           '
    write(stdout, *) '-------------------------------------------------------------------------------'
    write(stdout, *)

    call mainclock%print(stdout, ' Total run time:')


end program zaghop
