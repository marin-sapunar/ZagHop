!--------------------------------------------------------------------------------------------------
! MODULE: system_type_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date October, 2017
!
! DESCRIPTION:
!> @brief Defines type for holding all variables about the system in the dynamics calculation.
!--------------------------------------------------------------------------------------------------
module system_type_mod
    ! Import variables
    use global_defs
    use constants
    use string_mod, only : pad_right
    use tsh_method_mod, only : tsh_method
    use vc_evaluator_mod
    use vc_model_mod
    use coord_transform_mod, only : normal_modes
    implicit none

    private
    public :: system_type
    public :: vc_data
    public :: vc_nmode
    public :: trajectory_data, tr1, tr2
    public :: memory
    public :: allocate_trajectory_data
    public :: trajectory_n_point
    public :: trajectory_point
    public :: trajectory_goto
    public :: trajectory_next
    public :: trajectory_write_final
    public :: trajectory_save_state
    public :: trajectory_load_state
    public :: output_options
    public :: output_file_name
    public :: ekin

    !----------------------------------------------------------------------------------------------
    ! TYPE: system_type
    !> @brief Holds all information about the system at a single time step.
    !----------------------------------------------------------------------------------------------
    type system_type
        integer :: step = -1 !< Index of the last point of the output grid at or before this point.
        real(dp) :: time = -1.0_dp !< Current time.
        real(dp) :: dt = 0.0_dp !< Step from the previous point to this point.
        logical :: on_grid = .true. !< The point lies on the output grid.
        logical :: written = .false. !< The point was passed to the output.

        ! Full system:
        integer :: ndim = 3 !< Number of dimensions.
        character(len=:), allocatable :: classical_coordinates !< Coordinates of the classical
                                                               !! propagation (geom/velo).
        integer :: natom = 0 !< Total number of atoms.
        character(len=16), allocatable :: sym(:) !< List of atom symbols (or DOF labels).
        real(dp), allocatable :: chrg(:) !< List of atom partial charges.
        real(dp), allocatable :: mass(:) !< List of atom masses.
        real(dp), allocatable :: geom(:, :) !< Geometry.
        real(dp), allocatable :: velo(:, :) !< Velocity.
        real(dp), allocatable :: grad(:, :) !< Gradient.

        ! QM system:
        integer :: qnatom = 0 !< Number of QM atoms.
        integer, allocatable :: qind(:) !< Indexes of the QM atoms in the full system.
        class(potential_evaluator), allocatable :: pot
        class(tsh_method), allocatable :: wf !< Surface hopping method and its electronic data.

        ! MM system:
        integer :: mnatom = 0 !< Number of MM atoms.
        integer, allocatable :: mind(:) !< Indexes of the MM atoms in the full system.
        real(dp) :: men(2) = 0.0_dp !< MM energies (1 - QM model sys, 2 - MM sys).

        real(dp) :: pbcbox(1:6) = 0.0_dp
    contains
        procedure :: writestep => traj_writestep
        procedure :: writeheader => traj_writeheader
        procedure :: write_state_layout => traj_write_state_layout
        procedure :: open_files => traj_open_files
        procedure :: tote => traj_tote !< Total energy of the system.
        procedure :: kine => traj_kine !< Kinetic energy of the system.
        procedure :: pote => traj_pote !< Potential energy of the system.
        procedure :: qkine => traj_qkine !< QM kinetic energy.
        procedure :: mkine => traj_mkine !< MM kinetic energy.
    end type system_type


    ! The stored points of the trajectory, in time order, are the slots point_slot(1:n_point) of
    ! trajectory_data. Usually the current point is the latest one. After a step back to an
    ! earlier point, the later points are kept only when a new point is inserted between them.
    integer :: memory = 5 !< Maximum number of stored points.
    type(system_type), allocatable, target :: trajectory_data(:)
    type(system_type), pointer :: tr1 !< Pointer to trajectory data for current time step.
    type(system_type), pointer :: tr2 !< Pointer to trajectory data for previous time step.
    integer, allocatable :: point_slot(:) !< Slots of the stored points in time order.
    integer :: n_point = 0 !< Number of stored points.
    integer :: i_current = 0 !< Position of the current point in point_slot.

    type(vc_model), target :: vc_data !< LVC model for the system.
    type(normal_modes), target :: vc_nmode !< Normal modes of the LVC model.


contains

    !----------------------------------------------------------------------------------------------
    ! FUNCTION: ekin
    !> @brief Calculate the kinetic energy for a set of particles.
    !----------------------------------------------------------------------------------------------
    pure function ekin(mass, velo) result (en)
        real(dp), intent(in) :: mass(:)
        real(dp), intent(in) :: velo(:, :)
        real(dp) :: en
        integer :: i

        en = 0.0_dp
        do i = 1, size(mass)
            en = en + mass(i) * dot_product(velo(:, i), velo(:, i)) * 0.5_dp
        end do
    end function ekin

    pure function traj_tote(t) result(tote)
        class(system_type), intent(in) :: t
        real(dp) :: tote
        tote = t%pote() + t%kine()
    end function traj_tote

    pure function traj_pote(t) result(pote)
        class(system_type), intent(in) :: t
        real(dp) :: pote
        pote = t%pot%get_energy(t%wf%basis, t%wf%active_state) + t%men(2) - t%men(1)
    end function traj_pote

    pure function traj_kine(t) result(kine)
        class(system_type), intent(in) :: t
        real(dp) :: kine
        kine = ekin(t%mass, t%velo)
    end function traj_kine

    pure function traj_qkine(t) result(qkine)
        class(system_type), intent(in) :: t
        real(dp) :: qkine
        qkine = ekin(t%mass(t%qind), t%velo(:, t%qind))
    end function traj_qkine

    pure function traj_mkine(t) result(mkine)
        class(system_type), intent(in) :: t
        real(dp) :: mkine
        mkine = ekin(t%mass(t%mind), t%velo(:, t%mind))
    end function traj_mkine
            

    !---------------------------------------------------------------------------------------------
    ! SUBROUTINE: allocate_trajectory_data
    !> @brief Allocate trajectory_data array and associate pointers.
    !> @details
    !! On the first call, the initial point is stored in the first slot. Later calls can increase
    !! the number of slots, keeping the stored points.
    !---------------------------------------------------------------------------------------------
    subroutine allocate_trajectory_data(n)
        integer, intent(in), optional :: n !< Number of slots.
        type(system_type), allocatable, target :: new_data(:)
        integer :: i

        if (present(n)) memory = n
        if (.not. allocated(trajectory_data)) then
            allocate(trajectory_data(memory))
            allocate(point_slot(memory))
            n_point = 1
            i_current = 1
            point_slot(1) = 1
        else
            if (memory <= size(trajectory_data)) then
                memory = size(trajectory_data)
                return
            end if
            allocate(new_data(memory))
            do i = 1, n_point
                new_data(i) = trajectory_data(point_slot(i))
            end do
            call move_alloc(new_data, trajectory_data)
            deallocate(point_slot)
            allocate(point_slot(memory))
            point_slot(1:n_point) = [(i, i = 1, n_point)]
        end if
        call set_pointers()
    end subroutine allocate_trajectory_data


    !---------------------------------------------------------------------------------------------
    ! SUBROUTINE: set_pointers
    !> @brief Point tr1 to the current point and tr2 to the point before it.
    !> @details
    !! Without an earlier point, tr2 points to the current point.
    !---------------------------------------------------------------------------------------------
    subroutine set_pointers()
        tr1 => trajectory_data(point_slot(i_current))
        tr2 => trajectory_data(point_slot(max(i_current - 1, 1)))
    end subroutine set_pointers


    !---------------------------------------------------------------------------------------------
    ! FUNCTION: trajectory_n_point
    !> @brief Number of stored points.
    !---------------------------------------------------------------------------------------------
    integer function trajectory_n_point() result(n)
        n = n_point
    end function trajectory_n_point


    !---------------------------------------------------------------------------------------------
    ! FUNCTION: trajectory_point
    !> @brief Pointer to the i-th stored point, in time order.
    !---------------------------------------------------------------------------------------------
    function trajectory_point(i) result(t)
        integer, intent(in) :: i
        type(system_type), pointer :: t

        if (i < 1 .or. i > n_point) call errstop('trajectory_point', 'Point not stored.', i)
        t => trajectory_data(point_slot(i))
    end function trajectory_point


    !---------------------------------------------------------------------------------------------
    ! SUBROUTINE: trajectory_goto
    !> @brief Make the stored point at time current.
    !> @details
    !! Later points are discarded, unless keep_later is true. They are then kept for the next call
    !! of trajectory_next, which inserts a new point before them.
    !---------------------------------------------------------------------------------------------
    subroutine trajectory_goto(time, keep_later)
        real(dp), intent(in) :: time
        logical, intent(in) :: keep_later
        integer :: i

        do i = 1, n_point
            if (trajectory_data(point_slot(i))%time == time) exit
        end do
        if (i > n_point) then
            write(stderr, '(1x,a)') 'Error in trajectory_goto. No stored point at the requested time.'
            write(stderr, '(3x,a,e24.16)') 'Requested time: ', time
            write(stderr, '(3x,a,1000e24.16)') 'Stored times: ',                                   &
            &   (trajectory_data(point_slot(i))%time, i = 1, n_point)
            call errstop('trajectory_goto')
        end if
        i_current = i
        if (.not. keep_later) n_point = i
        call set_pointers()
    end subroutine trajectory_goto


    !---------------------------------------------------------------------------------------------
    ! SUBROUTINE: trajectory_next
    !> @brief Add a new point after the current one, as a copy of it, and make it current.
    !> @details
    !! Without dt, the new point is the next point of the output grid with spacing dt_0. If the
    !! current point is on the grid, the step is dt_0. Otherwise, the step ends at the time of the
    !! last grid point before the current point plus dt_0. The current point must then be the
    !! latest stored point.
    !!
    !! With dt, the new point lies off the output grid, dt after the current point. Stored points
    !! after the current point are kept after the new point.
    !!
    !! If all slots are used, the oldest point is overwritten. This is an error if that point was
    !! not yet passed to the output, or if it is the current point or the one before it.
    !---------------------------------------------------------------------------------------------
    subroutine trajectory_next(dt_0, dt)
        real(dp), intent(in) :: dt_0 !< Time step of the output grid.
        real(dp), intent(in), optional :: dt !< Step to a new point off the output grid.
        real(dp) :: time, new_dt
        integer :: step, slot, i
        logical :: on_grid

        if (present(dt)) then
            time = tr1%time + dt
            new_dt = dt
            step = tr1%step
            on_grid = .false.
        else
            if (i_current /= n_point) then
                call errstop('trajectory_next', 'Stored points after the current point.')
            end if
            if (tr1%on_grid) then
                time = tr1%time + dt_0
                new_dt = dt_0
            else
                do i = i_current - 1, 1, -1
                    if (trajectory_data(point_slot(i))%on_grid) exit
                end do
                if (i < 1) call errstop('trajectory_next', 'Last point of output grid not stored.')
                time = trajectory_data(point_slot(i))%time + dt_0
                new_dt = time - tr1%time
            end if
            step = tr1%step + 1
            on_grid = .true.
        end if

        ! Find a free slot.
        if (n_point < memory) then
            do slot = 1, memory
                if (all(point_slot(1:n_point) /= slot)) exit
            end do
        else
            associate(oldest => trajectory_data(point_slot(1)))
                if (i_current <= 2 .or. (oldest%on_grid .and. .not. oldest%written)) then
                    write(stderr, '(1x,a)') 'Error in trajectory_next. All stored points are needed.'
                    write(stderr, '(3x,a,i0,a)') 'Increase the number of stored points (', memory, ').'
                    call errstop('trajectory_next')
                end if
            end associate
            slot = point_slot(1)
            point_slot(1:n_point - 1) = point_slot(2:n_point)
            n_point = n_point - 1
            i_current = i_current - 1
        end if

        trajectory_data(slot) = trajectory_data(point_slot(i_current))
        point_slot(i_current + 2:n_point + 1) = point_slot(i_current + 1:n_point)
        point_slot(i_current + 1) = slot
        n_point = n_point + 1
        i_current = i_current + 1
        call set_pointers()
        tr1%time = time
        tr1%dt = new_dt
        tr1%step = step
        tr1%on_grid = on_grid
        tr1%written = .false.
    end subroutine trajectory_next


    !---------------------------------------------------------------------------------------------
    ! SUBROUTINE: trajectory_write_final
    !> @brief Write the stored points of the output grid up to time t_final, which were not yet
    !! written, in time order.
    !> @details
    !! Only every printerval-th point of the grid is written to the output files.
    !---------------------------------------------------------------------------------------------
    subroutine trajectory_write_final(t_final, popt, punit, printerval)
        real(dp), intent(in) :: t_final !< Points up to this time are final.
        logical, intent(in) :: popt(:) !< Print options.
        integer, intent(in) :: punit(:) !< Units for output files.
        integer, intent(in) :: printerval !< Number of grid steps between written points.
        integer :: i

        do i = 1, n_point
            associate(t => trajectory_data(point_slot(i)))
                if (t%time > t_final) exit
                if (t%on_grid .and. .not. t%written) then
                    if (mod(t%step, printerval) == 0) call t%writestep(popt, punit)
                    t%written = .true.
                end if
            end associate
        end do
    end subroutine trajectory_write_final


    !---------------------------------------------------------------------------------------------
    ! SUBROUTINE: trajectory_save_state
    !> @brief Add the stored points of the trajectory to a JSON object, for a backup.
    !> @details
    !! The points are saved in time order, in the array 'points', with the position of the current
    !! point (i_current). Each point holds its time and grid data, geometry, velocity, gradient
    !! and MM data, and the states of its potential ('potential') and surface hopping method
    !! ('method'). The atoms and options of the system are not saved; they come from the input.
    !---------------------------------------------------------------------------------------------
    subroutine trajectory_save_state(json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_add_array
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the points are added.
        type(json_value), pointer :: arr, pt, sub
        integer :: i

        call json%add(obj, 'i_current', i_current)
        call json%create_array(arr, 'points')
        do i = 1, n_point
            associate(t => trajectory_data(point_slot(i)))
                call json%create_object(pt, '')
                call json%add(pt, 'step', t%step)
                call json%add(pt, 'time', t%time)
                call json%add(pt, 'dt', t%dt)
                call json%add(pt, 'on_grid', t%on_grid)
                call json%add(pt, 'written', t%written)
                call json_add_array(json, pt, 'geom', t%geom)
                call json_add_array(json, pt, 'velo', t%velo)
                call json_add_array(json, pt, 'grad', t%grad)
                call json%add(pt, 'men', t%men)
                call json%add(pt, 'pbcbox', t%pbcbox)
                call json%create_object(sub, 'potential')
                call t%pot%save_state(json, sub)
                call json%add(pt, sub)
                call json%create_object(sub, 'method')
                call t%wf%save_state(json, sub)
                call json%add(pt, sub)
                call json%add(arr, pt)
            end associate
        end do
        call json%add(obj, arr)
    end subroutine trajectory_save_state


    !---------------------------------------------------------------------------------------------
    ! SUBROUTINE: trajectory_load_state
    !> @brief Restore the stored points of the trajectory saved by trajectory_save_state.
    !> @details
    !! The current point must be set up from the input (atoms, potential and method), and serves
    !! as the template of all loaded points, whose saved state then replaces that of the template.
    !! Stops if there are more points than slots, or if the shape of the geometry differs.
    !---------------------------------------------------------------------------------------------
    subroutine trajectory_load_state(json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_get_required, json_get_child
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object holding the points.
        type(json_value), pointer :: arr, pt
        type(system_type) :: template
        real(dp), allocatable :: men(:), pbcbox(:)
        integer :: i, n, cur

        template = tr1
        arr => json_get_child(json, obj, 'points')
        n = json%count(arr)
        call json_get_required(json, obj, 'i_current', cur)
        if (n < 1 .or. n > memory .or. cur < 1 .or. cur > n) then
            write(stderr, '(1x,a)') 'Error in system_type_mod, trajectory_load_state subroutine.'
            write(stderr, '(3x,a,i0,a,i0,a,i0,a)') 'Backup has ', n, ' points (current: ', cur,  &
            &   '), the trajectory can store ', memory, '.'
            stop 1
        end if
        do i = 1, n
            call json%get_child(arr, i, pt)
            trajectory_data(i) = template
            associate(t => trajectory_data(i))
                call json_get_required(json, pt, 'step', t%step)
                call json_get_required(json, pt, 'time', t%time)
                call json_get_required(json, pt, 'dt', t%dt)
                call json_get_required(json, pt, 'on_grid', t%on_grid)
                call json_get_required(json, pt, 'written', t%written)
                call json_get_required(json, pt, 'geom', t%geom, shape(template%geom))
                call json_get_required(json, pt, 'velo', t%velo, shape(template%velo))
                call json_get_required(json, pt, 'grad', t%grad, shape(template%grad))
                call json_get_required(json, pt, 'men', men, 2)
                t%men = men
                call json_get_required(json, pt, 'pbcbox', pbcbox, 6)
                t%pbcbox = pbcbox
                call t%pot%load_state(json, json_get_child(json, pt, 'potential'))
                call t%wf%load_state(json, json_get_child(json, pt, 'method'))
            end associate
        end do
        point_slot(1:n) = [(i, i = 1, n)]
        n_point = n
        i_current = cur
        call set_pointers()
    end subroutine trajectory_load_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: traj_open_files
    !
    ! DESCRIPTION:
    !> @brief Prepare output files for writing dynamics outputs.
    !----------------------------------------------------------------------------------------------
    subroutine traj_open_files(t, popt, punit, res_dir)
        class(system_type) :: t
        logical, intent(in) :: popt(:) !< Print options.
        integer, intent(out) :: punit(:) !< Units for output files.
        character(len=*), intent(in) :: res_dir ! Directory for writing results.
        integer, allocatable :: iopt(:)
        integer :: i

        iopt = output_options(t, popt)
        do i = 1, size(iopt)
            open(newunit=punit(iopt(i)), file=res_dir//'/'//output_file_name(iopt(i)),           &
            &    action='write', position='append')
        end do
    end subroutine traj_open_files


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: output_options
    !> @brief Indices of the output files opened by open_files for the print options popt.
    !> @details
    !! The MM energies (mm.dat, index 50) are written with the energies, if there are MM atoms.
    !----------------------------------------------------------------------------------------------
    function output_options(t, popt) result(iopt)
        class(system_type), intent(in) :: t
        logical, intent(in) :: popt(:) !< Print options.
        integer, allocatable :: iopt(:)
        integer :: i

        allocate(iopt(0))
        do i = 1, size(popt)
            if (i == 50) cycle
            if (popt(i) .and. output_file_name(i) /= '') iopt = [iopt, i]
        end do
        if (popt(1) .and. t%mnatom > 0) iopt = [iopt, 50]
    end function output_options


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: output_file_name
    !> @brief Name of the output file with index i (in the results directory), or '' if there is
    !! no file with this index.
    !----------------------------------------------------------------------------------------------
    pure function output_file_name(i) result(fname)
        integer, intent(in) :: i
        character(len=:), allocatable :: fname

        select case(i)
        case(1)
            fname = 'energy.dat'
        case(2)
            fname = 'trajectory.xyz'
        case(3)
            fname = 'geometry'
        case(4)
            fname = 'velocity'
        case(5)
            fname = 'gradient'
        case(6)
            fname = 'cwf.dat'
        case(7)
            fname = 'prob'
        case(10)
            fname = 'oscill.dat'
        case(20)
            fname = 'trans_ldiab'
        case(21)
            fname = 'trans_diab'
        case(22)
            fname = 'trans_adiab'
        case(30)
            fname = 'cstate_ldiab'
        case(31)
            fname = 'cstate_diab'
        case(32)
            fname = 'cstate_adiab'
        case(50)
            fname = 'mm.dat'
        case default
            fname = ''
        end select
    end function output_file_name



    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: traj_writestep
    !
    ! DESCRIPTION:
    !> @brief Write requested output from each step.
    !----------------------------------------------------------------------------------------------
    subroutine traj_writestep(t, popt, punit)
        class(system_type) :: t
        logical, intent(in) :: popt(:) !< Print options.
        integer, intent(in) :: punit(:) !< Units for output files.
        integer :: i
        real(dp) :: time_fs
        real(dp), allocatable :: wrk(:, :)
        real(dp), allocatable :: cstate(:)

        time_fs = t%time * aut_fs
        if (popt(1)) then
            write(punit(1), 1001, advance='no') time_fs, t%wf%active_state
            write(punit(1), 1002, advance='no') t%tote(), t%pote()
            write(punit(1), 1002, advance='no') t%pot%get_energy(t%wf%basis)
            write(punit(1), *)

            if (t%mnatom > 0) then
                write(punit(50), 1001, advance='no') time_fs
                write(punit(50), 1002, advance='no') t%mkine(), t%men
                write(punit(50), 1002, advance='no') t%pbcbox
                write(punit(50), *)
            end if
        end if

        if (popt(2)) then
            write(punit(2), *) t%natom
            write(punit(2), *) '# t= ', time_fs, 'fs, state=', t%wf%active_state
            do i = 1, t%natom
                write(punit(2), 1003) pad_right(t%sym(i), 2), t%geom(:, i) * a0_A
            end do
        end if

        if (popt(3)) then
            do i = 1, t%natom
                write(punit(3), 1002) t%geom(:, i)
            end do
        end if

        if (popt(4)) then
            do i = 1, t%natom
                write(punit(4), 1002) t%velo(:, i)
            end do
        end if

        if (popt(5)) then
            do i = 1, t%natom
                write(punit(5), 1002) t%grad(:, i)
            end do
        end if

        if (popt(6)) then
            call t%wf%write_coefficients(punit(6))
        end if

        if (popt(7)) then
            write(punit(7), 1006) t%wf%prob
        end if

        ! The locally diabatic basis is that of the propagation: the orthogonalised overlaps with
        ! the states of the previous point, from the link of the potential.
        if (popt(20) .and. t%pot%link%valid()) then
            write(punit(20), *) '# t= ', time_fs, 'fs, state=', t%wf%active_state
            wrk = t%pot%link%t
            do i = 1, size(wrk, 2)
                write(punit(20), 1006) wrk(:, i)
            end do
        end if

        if (popt(21)) then
           write(punit(21), *) '# t= ', time_fs, 'fs, state=', t%wf%active_state
           wrk = t%pot%get_transformation('diabatic', 'adiabatic')
           do i = 1, t%wf%n_state
               write(punit(21), 1006) wrk(:, i)
           end do
       end if

        if (popt(22)) then
            write(punit(22), *) '# t= ', time_fs, 'fs, state=', t%wf%active_state
            wrk = t%pot%get_transformation('adiabatic', 'adiabatic')
            do i = 1, t%wf%n_state
                write(punit(22), 1006) wrk(:, i)
            end do
        end if

        if (popt(30) .and. t%pot%link%valid()) then
            if (t%pot%link%basis == t%wf%basis) then
                write(punit(30), 1006) t%pot%link%t(:, t%wf%active_state)
            end if
        end if

        if (popt(31)) then
            wrk = t%pot%get_transformation('diabatic', t%wf%basis)
            write(punit(31), 1006) wrk(:, t%wf%active_state)
        end if

        if (popt(32)) then
            wrk = t%pot%get_transformation('adiabatic', t%wf%basis)
            write(punit(32), 1006) wrk(:, t%wf%active_state)
        end if
        

1001 format (f12.5,1x,i4)
1002 format (1x,1000f18.10)
1003 format (1x,a,2x,1000f18.10)
1006 format (1x,1000e22.12)
    end subroutine traj_writestep


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Traj_writeheader
    !
    ! DESCRIPTION:
    !> @brief Initialize output files and write column headers.
    !----------------------------------------------------------------------------------------------
    subroutine traj_writeheader(t, popt, punit)
        class(system_type) :: t
        logical, intent(in) :: popt(:) !< Print options.
        integer, intent(in) :: punit(:) !< Units for output files.
        integer :: i
        real(dp) :: time_fs

        time_fs = t%time * aut_fs
        if (popt(1)) then
            write(punit(1), '(a)', advance='no') '#'
            write(punit(1), '(a12,1x,a4)', advance='no') 't,', 'cst,'
            write(punit(1), '(1x,2a18)', advance='no') 'etot,', 'epot,'
            write(punit(1), '(1x)', advance='no')
            do i = 1, t%wf%n_state - 1
                write(punit(1), '(12x,a2,i3.3,a)', advance='no') 'qe', i, ','
            end do
            write(punit(1), '(12x,a2,i3.3)', advance='no') 'qe', i
            write(punit(1), *)

            if (t%mnatom > 0) then
                write(punit(50), '(a)', advance='no') '#'
                write(punit(50), '(a12,1x)', advance='no') 't,'
                write(punit(50), '(2x,3a18)', advance='no') 'mekin,', 'me1,', 'me2,'
                write(punit(50), '(1x,5a18,a17)', advance='no') 'pbc1,', 'pbc2,', 'pbc3,', 'pbc4,', 'pbc5,', 'pbc6'
                write(punit(50), *)
            end if
        end if
    end subroutine traj_writeheader


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: traj_write_state_layout
    !> @brief Write the layout of the electronic states to the file states.json.
    !> @details
    !! The file gives the format version, the basis of the dynamics (dynamics_basis), in which the
    !! cst and qe columns of energy.dat are given, and the layout of the states of the potential
    !! (see add_state_layout of the potential).
    !----------------------------------------------------------------------------------------------
    subroutine traj_write_state_layout(t, res_dir)
        use json_module, only : json_core, json_value
        class(system_type), intent(in) :: t
        character(len=*), intent(in) :: res_dir !< Directory for writing results.
        integer, parameter :: format_version = 1
        type(json_core) :: json
        type(json_value), pointer :: root

        call json%initialize()
        call json%create_object(root, '')
        call json%add(root, 'format_version', format_version)
        call json%add(root, 'dynamics_basis', t%wf%basis)
        call t%pot%add_state_layout(json, root)
        call json%print(root, res_dir//'/states.json')
        call json%destroy(root)
        if (json%failed()) then
            write(stderr, *) 'Warning in system_type_mod, traj_write_state_layout subroutine.'
            write(stderr, *) '  Failed to write '//res_dir//'/states.json.'
            call json%clear_exceptions()
        end if
    end subroutine traj_write_state_layout


end module system_type_mod
