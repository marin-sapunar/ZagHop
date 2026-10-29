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
    use tsh_method_mod, only : tsh_method
    use vc_evaluator_mod
    use vc_model_mod
    use coord_transform_mod, only : normal_modes
    implicit none

    private
    public :: system_type
    public :: vc_potential !< TODO move to more appropriate place
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
    public :: trajectory_write_backup
    public :: trajectory_read_backup
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
        character(len=2), allocatable :: sym(:) !< List of atom symbols.
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
    type(vc_evaluator) :: vc_potential !< LVC potential interface for the system.


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

        if (popt(1)) then
            open(newunit=punit(1), file=res_dir//'/energy.dat', action='write', position='append')

            if (t%mnatom > 0) then
                open(newunit=punit(50), file=res_dir//'/mm.dat', action='write', position='append')
            end if
        end if

        if (popt(2)) then
            open(newunit=punit(2), file=res_dir//'/trajectory.xyz', action='write', position='append')
        end if

        if (popt(3)) then
            open(newunit=punit(3), file=res_dir//'/geometry', action='write', position='append')
        end if

        if (popt(4)) then
            open(newunit=punit(4), file=res_dir//'/velocity', action='write', position='append')
        end if

        if (popt(5)) then
            open(newunit=punit(5), file=res_dir//'/gradient', action='write', position='append')
        end if

        if (popt(6)) then
            open(newunit=punit(6), file=res_dir//'/cwf.dat', action='write', position='append')
        end if
        
        if (popt(7)) then
            open(newunit=punit(7), file=res_dir//'/prob', action='write', position='append')
        end if
        
        if (popt(10)) then
            open(newunit=punit(10), file=res_dir//'/oscill.dat', action='write', position='append')
        end if

        if (popt(20)) then
            open(newunit=punit(20), file=res_dir//'/trans_ldiab', action='write', position='append')
        end if

        if (popt(21)) then
            open(newunit=punit(21), file=res_dir//'/trans_diab', action='write', position='append')
        end if

        if (popt(22)) then
            open(newunit=punit(22), file=res_dir//'/trans_adiab', action='write', position='append')
        end if

        if (popt(30)) then
            open(newunit=punit(30), file=res_dir//'/cstate_ldiab', action='write', position='append')
        end if

        if (popt(31)) then
            open(newunit=punit(31), file=res_dir//'/cstate_diab', action='write', position='append')
        end if

        if (popt(32)) then
            open(newunit=punit(32), file=res_dir//'/cstate_adiab', action='write', position='append')
        end if

    end subroutine traj_open_files



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
                write(punit(2), 1003) t%sym(i), t%geom(:, i) * a0_A
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

        if (popt(20) .and. t%step > 0) then
            write(punit(20), *) '# t= ', time_fs, 'fs, state=', t%wf%active_state
            wrk = t%pot%get_transformation('locally_diabatic', 'adiabatic')
            do i = 1, t%wf%n_state
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

        if (popt(30) .and. t%step > 0) then
            wrk = t%pot%get_transformation('locally_diabatic', t%wf%basis)
            write(punit(30), 1006) wrk(:, t%wf%active_state)
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
1003 format (1x,a2,2x,1000f18.10)
1006 format (1x,1000e22.12)
    end subroutine traj_writestep


    subroutine trajectory_write_backup(backup_file, t)
        character(len=*), intent(in) :: backup_file
        class(system_type), intent(in) :: t(:)
        integer :: i, ounit

!         open(newunit=ounit, file=backup_file, action='write')
!         write(ounit, *) memory
!         do i = 1, memory
!             write(ounit, *) t(i)%step
!             write(ounit, *) t(i)%substep
!             write(ounit, *) t(i)%time
!             write(ounit, *) t(i)%ndim
!             write(ounit, *) t(i)%natom
!             write(ounit, *) allocated(t(i)%sym)
!             if (allocated(t(i)%sym)) write(ounit, *) t(i)%sym(:)
!             write(ounit, *) allocated(t(i)%chrg)
!             if (allocated(t(i)%chrg)) write(ounit, *) t(i)%chrg(:)
!             write(ounit, *) allocated(t(i)%mass)
!             if (allocated(t(i)%mass)) write(ounit, *) t(i)%mass(:)
!             write(ounit, *) allocated(t(i)%geom)
!             if (allocated(t(i)%geom)) write(ounit, *) t(i)%geom(:, :)
!             write(ounit, *) allocated(t(i)%velo)
!             if (allocated(t(i)%velo)) write(ounit, *) t(i)%velo(:, :)
!             write(ounit, *) allocated(t(i)%grad)
!             if (allocated(t(i)%grad)) write(ounit, *) t(i)%grad(:, :)
!             write(ounit, *) t(i)%qnatom
!             write(ounit, *) allocated(t(i)%qind)
!             if (allocated(t(i)%qind)) write(ounit, *) t(i)%qind(:)
!             write(ounit, *) t(i)%max_nstate
!             write(ounit, *) t(i)%cstate
!             write(ounit, *) allocated(t(i)%qe)
!             if (allocated(t(i)%qe)) write(ounit, *) t(i)%qe(:)
!             write(ounit, *) allocated(t(i)%qo)
!             if (allocated(t(i)%qo)) write(ounit, *) t(i)%qo(:)
!             write(ounit, *) t(i)%mnatom
!             write(ounit, *) allocated(t(i)%mind)
!             if (allocated(t(i)%mind)) write(ounit, *) t(i)%mind(:)
!             write(ounit, *) t(i)%me(2)
!             write(ounit, *) allocated(t(i)%cwf)
!             if (allocated(t(i)%cwf)) write(ounit, *) t(i)%cwf(:)
!             write(ounit, *) allocated(t(i)%prob)
!             if (allocated(t(i)%prob)) write(ounit, *) t(i)%prob(:)
!             write(ounit, *) allocated(t(i)%olap)
!             if (allocated(t(i)%olap)) write(ounit, *) t(i)%olap(:, :)
!             write(ounit, *) allocated(t(i)%nadv)
!             if (allocated(t(i)%nadv)) write(ounit, *) t(i)%nadv(:, :, :)
!             write(ounit, *) allocated(t(i)%sov)
!             if (allocated(t(i)%sov)) write(ounit, *) t(i)%sov(:, :)
! !            write(ounit, *) allocated(t(i)%spinv)
! !            if (allocated(t(i)%spinv)) write(ounit, *) t(i)%spinv(:)
!             write(ounit, *) t(i)%pbcbox
!             write(ounit, *) allocated(t(i)%phase)
!             if (allocated(t(i)%phase)) write(ounit, *) t(i)%phase
!         end do
!         close(ounit)
    end subroutine trajectory_write_backup


    subroutine trajectory_read_backup(backup_file, t)
        character(len=*), intent(in) :: backup_file
        class(system_type), intent(inout) :: t(:)
        integer :: i, iunit
        logical :: check

!         open(newunit=iunit, file=backup_file, status='old', action='read')
!         read(iunit, *) i
!         if (i /= memory) stop 'Error. Bad memory parameter in system_var.'
!         do i = 1, memory
!             read(iunit, *) t(i)%step
!             read(iunit, *) t(i)%substep
!             read(iunit, *) t(i)%time
!             read(iunit, *) t(i)%ndim
!             read(iunit, *) t(i)%natom
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%sym(t(i)%natom))
!                 read(iunit, *) t(i)%sym(:)
!             end if
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%chrg(t(i)%natom))
!                 read(iunit, *) t(i)%chrg(:)
!             end if
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%mass(t(i)%natom))
!                 read(iunit, *) t(i)%mass(:)
!             end if
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%geom(t(i)%ndim, t(i)%natom))
!                 read(iunit, *) t(i)%geom(:, :)
!             end if
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%velo(t(i)%ndim, t(i)%natom))
!                 read(iunit, *) t(i)%velo(:, :)
!             end if
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%grad(t(i)%ndim, t(i)%natom))
!                 read(iunit, *) t(i)%grad(:, :)
!             end if
!             read(iunit, *) t(i)%qnatom
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%qind(t(i)%qnatom))
!                 read(iunit, *) t(i)%qind(:)
!             end if
!             read(iunit, *) t(i)%max_nstate
!             read(iunit, *) t(i)%cstate
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%qe(t(i)%max_nstate))
!                 read(iunit, *) t(i)%qe(:)
!             end if
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%qo(t(i)%max_nstate-1))
!                 read(iunit, *) t(i)%qo(:)
!             end if
!             read(iunit, *) t(i)%mnatom
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%mind(t(i)%mnatom))
!                 read(iunit, *) t(i)%mind(:)
!             end if
!             read(iunit, *) t(i)%me(2)
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%cwf(t(i)%max_nstate))
!                 read(iunit, *) t(i)%cwf(:)
!             end if
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%prob(t(i)%max_nstate))
!                 read(iunit, *) t(i)%prob(:)
!             end if
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%olap(t(i)%max_nstate, t(i)%max_nstate))
!                 read(iunit, *) t(i)%olap(:,:)
!             end if
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%nadv(t(i)%ndim*t(i)%qnatom, t(i)%max_nstate, t(i)%max_nstate))
!                 read(iunit, *) t(i)%nadv(:, :, :)
!             end if
!             if (check) then
!                 allocate(t(i)%sov(t(i)%max_nstate, t(i)%max_nstate))
!                 read(iunit, *) t(i)%sov(:, :)
!             end if
! !            if (check) then
! !                allocate(t(i)%spinv(t(i)%max_nstate))
! !                read(iunit, *) t(i)%spinv(:)
! !            end if
!             read(iunit, *) t(i)%pbcbox
!             read(iunit, *) check
!             if (check) then
!                 allocate(t(i)%phase(t(i)%max_nstate))
!                 read(iunit, *) t(i)%phase(:)
!             end if
!         end do
!         close(iunit)
    end subroutine trajectory_read_backup


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


end module system_type_mod
