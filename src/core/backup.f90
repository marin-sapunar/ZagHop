!--------------------------------------------------------------------------------------------------
! MODULE: backup_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2026
!
! DESCRIPTION:
!> @brief Backup file for restarting a trajectory.
!> @details
!! The backup is a JSON file with everything needed to continue the trajectory from the point at
!! which it was written, giving the same results as an uninterrupted run:
!! - format_version
!! - trajectory - The stored points of the trajectory (see trajectory_save_state), each with the
!!   states of its potential and surface hopping method, the last id of the sets of states
!!   (last_id), the total energy at the start of the trajectory (t0_tot_en), the time at which
!!   the target state was reached (t_target), and the size of each output file (output_files).
!! - rng - The state of the random number generator.
!!
!! The setup of the system (atoms, potential, method and options) is not saved. It is built from
!! the input, which must be the same as for the saved run, except for the options that do not
!! change the dynamics up to the saved point (such as max_time).
!!
!! The output files can contain lines written after the backup. On restart, each file is cut back
!! to its size at the time of the backup before the trajectory continues.
!--------------------------------------------------------------------------------------------------
module backup_mod
    use global_defs
    use constants, only : aut_fs
    use control_var, only : ctrl
    use system_type_mod, only : tr1, trajectory_save_state, trajectory_load_state,              &
    &                           output_options, output_file_name
    use json_module, only : json_core, json_value
    use json_util_mod, only : json_get_required, json_get_child, json_check_value
    implicit none

    private
    public :: write_backup
    public :: read_backup


    integer, parameter :: format_version = 1 !< Version of the backup file format.


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: write_backup
    !> @brief Write the backup file ctrl%bufile.
    !> @details
    !! The output files are flushed and their sizes are saved. The file is written as
    !! <bufile>.tmp and then renamed, so the previous backup stays complete if the program stops
    !! while writing.
    !----------------------------------------------------------------------------------------------
    subroutine write_backup()
        use state_overlap_mod, only : last_requested_id
        type(json_core) :: json
        type(json_value), pointer :: root, traj, files, entry, rng
        integer, allocatable :: iopt(:)
        integer(int64) :: n_byte
        character(len=24) :: n_byte_str
        character(len=:), allocatable :: tmp_file
        integer :: i, unit, stat

        call json%initialize()
        call json%create_object(root, '')
        call json%add(root, 'format_version', format_version)

        call json%create_object(traj, 'trajectory')
        call trajectory_save_state(json, traj)
        call json%add(traj, 'last_id', last_requested_id())
        call json%add(traj, 't0_tot_en', ctrl%t0_tot_en)
        call json%add(traj, 't_target', ctrl%t_target)
        call json%create_array(files, 'output_files')
        iopt = output_options(tr1, ctrl%print)
        do i = 1, size(iopt)
            unit = ctrl%print_units(iopt(i))
            flush(unit)
            inquire(unit=unit, size=n_byte)
            write(n_byte_str, '(i0)') n_byte
            call json%create_object(entry, '')
            call json%add(entry, 'print', iopt(i))
            call json%add(entry, 'file', output_file_name(iopt(i)))
            call json%add(entry, 'size', trim(n_byte_str))
            call json%add(files, entry)
        end do
        call json%add(traj, files)
        call json%add(root, traj)

        call json%create_object(rng, 'rng')
        call ctrl%rng%save_state(json, rng)
        call json%add(root, rng)

        tmp_file = ctrl%bufile//'.tmp'
        call json%print(root, tmp_file)
        call json%destroy(root)
        if (json%failed()) then
            call json%print_error_message(stderr)
            call errstop('write_backup', 'Failed to write '//tmp_file//'.')
        end if
        call execute_command_line('mv -f "'//tmp_file//'" "'//ctrl%bufile//'"', exitstat=stat)
        if (stat /= 0) call errstop('write_backup', 'Failed to rename '//tmp_file//'.', stat)
    end subroutine write_backup


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_backup
    !> @brief Continue the trajectory from the backup file ctrl%bufile.
    !> @details
    !! The input must be read first, so that the system, potential and method are set up. The
    !! stored points, the ids, the control values and the random number generator are restored,
    !! and the output files in the results directory are cut back to their size at the time of the
    !! backup. They must be the same files as in the saved run.
    !!
    !! If the target state was reached, the end time of the trajectory is set again from the time
    !! at which it was reached (t_target) and the current max_time, as in end_of_step_checks.
    !!
    !! If the saved point already reached this end time, there is nothing to continue: finished
    !! is set and nothing else is restored. The output files and the status file are left as the
    !! saved run wrote them.
    !----------------------------------------------------------------------------------------------
    subroutine read_backup(finished)
        use state_overlap_mod, only : restore_last_id
        use file_mod, only : need_file, need_dir
        logical, intent(out) :: finished !< The saved point is at or after the end time.
        real(dp), parameter :: t_tol = 1.0e-8_dp !< Tolerance of the end time, as tinydp in the
                                                 !! end_of_step_checks of the main program.
        type(json_core) :: json
        type(json_value), pointer :: root, traj
        integer :: last_id

        call need_file(ctrl%bufile, 'The backup file is needed to restart the trajectory.')
        call need_dir(ctrl%output_dir, 'The results directory is needed to restart the trajectory.')
        if (stdp1) then
            write(stdout, *)
            write(stdout, '(1x,3a)') 'Restarting trajectory from backup file ', ctrl%bufile, '.'
        end if

        call json%initialize()
        call json%load(file=ctrl%bufile, p=root)
        if (json%failed()) then
            call json%print_error_message(stderr)
            call errstop('read_backup', 'Failed to read '//ctrl%bufile//'.')
        end if
        call json_check_value(json, root, 'format_version', format_version)

        traj => json_get_child(json, root, 'trajectory')
        call trajectory_load_state(json, traj)
        call json_get_required(json, traj, 't_target', ctrl%t_target)
        if (ctrl%t_target >= 0.0_dp) then
            if (ctrl%t_target + ctrl%target_state_time < ctrl%max_time) then
                ctrl%max_time = ctrl%t_target + ctrl%target_state_time
                ctrl%target_state = -2
            else
                ctrl%target_state = -1
            end if
        end if
        finished = (tr1%time + t_tol >= ctrl%max_time)
        if (finished) then
            call json%destroy(root)
            return
        end if

        call json_get_required(json, traj, 'last_id', last_id)
        call restore_last_id(last_id)
        call json_get_required(json, traj, 't0_tot_en', ctrl%t0_tot_en)
        call restore_output_files(json, json_get_child(json, traj, 'output_files'))
        call ctrl%rng%load_state(json, json_get_child(json, root, 'rng'))
        call json%destroy(root)

        if (stdp1) then
            write(stdout, '(3x,a,f0.4,a)') 'Continuing from t=', tr1%time * aut_fs, ' fs.'
        end if
    end subroutine read_backup


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: restore_output_files
    !> @brief Cut the output files back to their sizes in the backup.
    !> @details
    !! Stops if the output files of the current print options are not those of the backup.
    !----------------------------------------------------------------------------------------------
    subroutine restore_output_files(json, files)
        use file_mod, only : truncate_file
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: files !< Array of output files.
        type(json_value), pointer :: entry
        integer, allocatable :: iopt(:)
        character(len=:), allocatable :: size_str
        integer(int64) :: n_byte
        integer :: i, ios

        iopt = output_options(tr1, ctrl%print)
        if (json%count(files) /= size(iopt)) then
            write(stderr, '(1x,a)') 'Error in backup_mod, restore_output_files subroutine.'
            write(stderr, '(3x,a,i0,a,i0,a)') 'The backup has ', json%count(files),               &
            &   ' output files, the print options give ', size(iopt), '.'
            stop 1
        end if
        do i = 1, size(iopt)
            call json%get_child(files, i, entry)
            call json_check_value(json, entry, 'print', iopt(i))
            call json_check_value(json, entry, 'file', output_file_name(iopt(i)))
            call json_get_required(json, entry, 'size', size_str)
            read(size_str, *, iostat=ios) n_byte
            if (ios /= 0) call errstop('restore_output_files', 'Bad size: '//size_str)
            call truncate_file(ctrl%output_dir//'/'//output_file_name(iopt(i)), n_byte)
        end do
    end subroutine restore_output_files


end module backup_mod
