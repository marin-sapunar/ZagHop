!--------------------------------------------------------------------------------------------------
! MODULE: system_var
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date October, 2017
!
! DESCRIPTION: 
!> @brief Defines type for holding all variables about the system in the dynamics calculation.
!--------------------------------------------------------------------------------------------------
module system_var
    ! Import variables
    use global_defs
    use constants
    implicit none

    private
    public :: t, memory
    public :: trajtype
    public :: trajectory_next
    public :: trajectory_rewind
    public :: trajectory_write_backup
    public :: trajectory_read_backup

    !----------------------------------------------------------------------------------------------
    ! TYPE: TrajType
    !> @brief Holds dynamics variables.
    !----------------------------------------------------------------------------------------------
    type trajtype
        integer :: step = -1 !< Current step.
        integer :: substep = -1 !< Signal that substeps are currently being performed.
        real(dp) :: time = -1.0_dp !< Current time.
        real(dp) :: last_hop = -1.0_dp !< Time at which the last hop occured.

        ! Full system:
        integer :: ndim = 3 !< Number of dimensions.
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
        integer :: nstate = 0
        integer :: nspin = 1 !< Number of different spins. 1 or 2 supported.
        integer :: nstate_spin(2) = 0 !< Number of singlet/triplet states.
        integer :: cstate = 0 !< Current state.
      ! integer :: cspin = 1 !< Current spin.
        real(dp), allocatable :: qe(:) !< Potential energies of the QM system.
        real(dp), allocatable :: qo(:) !< Oscillator strengths of the QM system.

        ! MM system:
        integer :: mnatom = 0 !< Number of MM atoms.
        integer, allocatable :: mind(:) !< Indexes of the MM atoms in the full system.
        real(dp) :: me(2) = 0.0_dp !< MM energies (1 - QM model sys, 2 - MM sys).

        ! Surface hopping variables:
        complex(dp), allocatable :: cwf(:) !< Coeffs of the el. states in the total wf.
        real(dp), allocatable :: prob(:) !< Probabilities of hopping during current step.
        real(dp), allocatable :: olap(:, :) !< Overlaps between wfs between this and previous step.
        real(dp), allocatable :: nadv(:, :, :) !< Nonadiabatic coupling vectors.
        real(dp), allocatable :: soc(:, :) !< Spin orbit couplings.

        real(dp) :: pbcbox(6) = 0.0_dp
    contains
        procedure :: writestep => traj_writestep
        procedure :: writeheader => traj_writeheader
        procedure :: tote => traj_tote !< Total energy of the system.
        procedure :: kine => traj_kine !< Kinetic energy of the system.
        procedure :: pote => traj_pote !< Potential energy of the system.
        procedure :: qkine => traj_qkine !< QM kinetic energy.
        procedure :: mkine => traj_mkine !< MM kinetic energy.
    end type trajtype


    integer, parameter :: memory = 5
    type(trajtype), allocatable :: t(:)


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
        class(trajtype), intent(in) :: t
        real(dp) :: tote
        tote = t%pote() + t%kine()
    end function traj_tote

    pure function traj_pote(t) result(pote)
        class(trajtype), intent(in) :: t
        real(dp) :: pote
        pote = t%qe(t%cstate) + t%me(2) - t%me(1)
    end function traj_pote

    pure function traj_kine(t) result(kine)
        class(trajtype), intent(in) :: t
        real(dp) :: kine
        kine = ekin(t%mass, t%velo)
    end function traj_kine

    pure function traj_qkine(t) result(qkine)
        class(trajtype), intent(in) :: t
        real(dp) :: qkine
        qkine = ekin(t%mass(t%qind), t%velo(:, t%qind))
    end function traj_qkine

    pure function traj_mkine(t) result(mkine)
        class(trajtype), intent(in) :: t
        real(dp) :: mkine
        mkine = ekin(t%mass(t%mind), t%velo(:, t%mind))
    end function traj_mkine


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: trajectory_next
    !> @brief Move array of trajectory values forward by one step at position pos.
    !----------------------------------------------------------------------------------------------
    subroutine trajectory_next(pos, t, dt, increment_time)
        integer, intent(in) :: pos
        type(trajtype), allocatable :: t(:)
        real(dp), intent(in) :: dt
        logical, intent(in) :: increment_time
        integer :: i

        do i = size(t), pos+1, -1
            t(i) = t(i-1)
        end do
        t(pos)%step = t(pos+1)%step + 1
        if (increment_time) then
            t(pos)%time = t(pos)%time + dt
        end if
    end subroutine trajectory_next


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: trajectory_rewind
    !> @brief Move array of trajectory values forward by one or more steps.
    !----------------------------------------------------------------------------------------------
    subroutine trajectory_rewind(t, nstep, increment_step)
        type(trajtype), allocatable :: t(:)
        integer, intent(in) :: nstep
        logical :: increment_step !< Increment step instead of copying value.
        integer :: i, step

        step = t(1)%step + 1
        do i = 1, size(t) - nstep
            t(i) = t(i+nstep)
        end do
        do i = size(t) - nstep + 1, size(t)
            if (allocated(t(i)%geom)) deallocate(t(i)%geom)
            if (allocated(t(i)%velo)) deallocate(t(i)%velo)
            if (allocated(t(i)%grad)) deallocate(t(i)%grad)
            if (allocated(t(i)%olap)) deallocate(t(i)%olap)
            if (allocated(t(i)%nadv)) deallocate(t(i)%nadv)
        end do
        if (increment_step) t(1)%step = step
    end subroutine trajectory_rewind


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Traj_WriteStep
    !
    ! DESCRIPTION:
    !> @brief Write requested output from each step.
    !----------------------------------------------------------------------------------------------
    subroutine traj_writestep(t, popt, res_dir)
        class(trajtype) :: t
        logical, intent(in) :: popt(:) !< Print options.
        character(len=*), intent(in) :: res_dir ! Directory for writing results.
        integer :: ounit, i
        real(dp) :: time_fs

        time_fs = t%time * aut_fs
        if (popt(1)) then
            open(newunit=ounit, file=res_dir//'/energy.dat', action='write', position='append')
            write(ounit, 1001, advance='no') time_fs, t%cstate
            write(ounit, 1002, advance='no') t%tote(), t%qe(t%cstate)
            write(ounit, 1002, advance='no') t%qe(:)
            write(ounit, *)
            close(ounit)

            if (t%mnatom > 0) then
                open(newunit=ounit, file=res_dir//'/mm.dat', action='write', position='append')
                write(ounit, 1001, advance='no') time_fs
                write(ounit, 1002, advance='no') t%mkine(), t%me
                write(ounit, 1002, advance='no') t%pbcbox
                write(ounit, *)
                close(ounit)
            end if
        end if

        if (popt(2)) then
            open(newunit=ounit, file=res_dir//'/trajectory.xyz', action='write', position='append')
            write(ounit, *) t%natom
            write(ounit, *) 't= ', time_fs, 'fs, state=', t%cstate
            do i = 1, t%natom
                write(ounit, 1003) t%sym(i), t%geom(:, i) * a0_A
            end do
            close(ounit)
        end if

        if (popt(3)) then
            open(newunit=ounit, file=res_dir//'/geometry', action='write', position='append')
            do i = 1, t%natom
                write(ounit, 1002) t%geom(:, i)
            end do
            close(ounit)
        end if
         
        if (popt(4)) then
            open(newunit=ounit, file=res_dir//'/velocity', action='write', position='append')
            do i = 1, t%natom
                write(ounit, 1002) t%velo(:, i)
            end do
            close(ounit)
        end if

        if (popt(5)) then
            open(newunit=ounit, file=res_dir//'/gradient', action='write', position='append')
            do i = 1, t%natom
                write(ounit, 1002) t%grad(:, i)
            end do
            close(ounit)
        end if

        if (popt(6)) then
            open(newunit=ounit, file=res_dir//'/cwf.dat', action='write', position='append')
            write(ounit, 1006) abs(t%cwf)
            close(ounit)
        end if

      ! if (popt(8)) then
      !     open(newunit=ounit, file=res_dir//'/prob.dat', action='write', position='append')
      !     write(ounit, 1006) t%prob
      !     close(ounit)
      ! end if

        if (popt(7)) then
            open(newunit=ounit, file=res_dir//'/overlap', action='write', position='append')
            write(ounit, *) 't= ', time_fs, 'fs, state=', t%cstate
            do i = 1, t%nstate
                write(ounit, 1006) t%olap(i, 1:t%nstate)
            end do
            close(ounit)
        end if

        if (popt(10)) then
            open(newunit=ounit, file=res_dir//'/oscill.dat', action='write', position='append')
            write(ounit, 1006) t%qo
            close(ounit)
        end if

        if (popt(11)) then
            open(newunit=ounit, file=res_dir//'/qm_traj.xyz', action='write', position='append')
            write(ounit, *) t%qnatom
            write(ounit, *) 't= ', time_fs, 'fs, state=', t%cstate
            do i = 1, t%qnatom
                write(ounit, 1003) t%sym(t%qind(i)), t%geom(:, t%qind(i)) * a0_A
            end do
            close(ounit)
        end if


1001 format (f12.5,2x,i3) 
1002 format (1x,1000f18.10)
1003 format (1x,a2,2x,1000f18.10) 
1006 format (1x,1000e22.12)
    end subroutine traj_writestep


    subroutine trajectory_write_backup(backup_file, t)
        character(len=*), intent(in) :: backup_file
        class(trajtype), intent(in) :: t(:)
        integer :: i, ounit

        open(newunit=ounit, file=backup_file, action='write')
        write(ounit, *) memory
        do i = 1, memory
            write(ounit, *) t(i)%step
            write(ounit, *) t(i)%substep
            write(ounit, *) t(i)%time
            write(ounit, *) t(i)%ndim
            write(ounit, *) t(i)%natom
            write(ounit, *) allocated(t(i)%sym)
            if (allocated(t(i)%sym)) write(ounit, *) t(i)%sym(:)
            write(ounit, *) allocated(t(i)%chrg)
            if (allocated(t(i)%chrg)) write(ounit, *) t(i)%chrg(:)
            write(ounit, *) allocated(t(i)%mass)
            if (allocated(t(i)%mass)) write(ounit, *) t(i)%mass(:)
            write(ounit, *) allocated(t(i)%geom)
            if (allocated(t(i)%geom)) write(ounit, *) t(i)%geom(:, :)
            write(ounit, *) allocated(t(i)%velo)
            if (allocated(t(i)%velo)) write(ounit, *) t(i)%velo(:, :)
            write(ounit, *) allocated(t(i)%grad)
            if (allocated(t(i)%grad)) write(ounit, *) t(i)%grad(:, :)
            write(ounit, *) t(i)%qnatom
            write(ounit, *) allocated(t(i)%qind)
            if (allocated(t(i)%qind)) write(ounit, *) t(i)%qind(:)
            write(ounit, *) t(i)%nspin
            write(ounit, *) t(i)%nstate
            write(ounit, *) t(i)%nstate_spin
            write(ounit, *) t(i)%cstate
            write(ounit, *) allocated(t(i)%qe)
            if (allocated(t(i)%qe)) write(ounit, *) t(i)%qe(:)
            write(ounit, *) allocated(t(i)%qo)
            if (allocated(t(i)%qo)) write(ounit, *) t(i)%qo(:)
            write(ounit, *) t(i)%mnatom
            write(ounit, *) allocated(t(i)%mind)
            if (allocated(t(i)%mind)) write(ounit, *) t(i)%mind(:)
            write(ounit, *) t(i)%me(2)
            write(ounit, *) allocated(t(i)%cwf)
            if (allocated(t(i)%cwf)) write(ounit, *) t(i)%cwf(:)
            write(ounit, *) allocated(t(i)%prob)
            if (allocated(t(i)%prob)) write(ounit, *) t(i)%prob(:)
            write(ounit, *) allocated(t(i)%olap)
            if (allocated(t(i)%olap)) write(ounit, *) t(i)%olap(:, :)
            write(ounit, *) allocated(t(i)%nadv)
            if (allocated(t(i)%nadv)) write(ounit, *) t(i)%nadv(:, :, :)
            write(ounit, *) t(i)%pbcbox(6)
        end do
        close(ounit)
    end subroutine trajectory_write_backup


    subroutine trajectory_read_backup(backup_file, t)
        character(len=*), intent(in) :: backup_file
        class(trajtype), intent(inout) :: t(:)
        integer :: i, iunit
        logical :: check

        open(newunit=iunit, file=backup_file, status='old', action='read')
        read(iunit, *) i
        if (i /= memory) stop 'Error. Bad memory parameter in system_var.'
        do i = 1, memory
            read(iunit, *) t(i)%step
            read(iunit, *) t(i)%substep
            read(iunit, *) t(i)%time
            read(iunit, *) t(i)%ndim
            read(iunit, *) t(i)%natom
            read(iunit, *) check
            if (check) then
                allocate(t(i)%sym(t(i)%natom))
                read(iunit, *) t(i)%sym(:)
            end if
            read(iunit, *) check
            if (check) then
                allocate(t(i)%chrg(t(i)%natom))
                read(iunit, *) t(i)%chrg(:)
            end if
            read(iunit, *) check
            if (check) then
                allocate(t(i)%mass(t(i)%natom))
                read(iunit, *) t(i)%mass(:)
            end if
            read(iunit, *) check
            if (check) then
                allocate(t(i)%geom(t(i)%ndim, t(i)%natom))
                read(iunit, *) t(i)%geom(:, :)
            end if
            read(iunit, *) check
            if (check) then
                allocate(t(i)%velo(t(i)%ndim, t(i)%natom))
                read(iunit, *) t(i)%velo(:, :)
            end if
            read(iunit, *) check
            if (check) then
                allocate(t(i)%grad(t(i)%ndim, t(i)%natom))
                read(iunit, *) t(i)%grad(:, :)
            end if
            read(iunit, *) t(i)%qnatom
            read(iunit, *) check
            if (check) then
                allocate(t(i)%qind(t(i)%qnatom))
                read(iunit, *) t(i)%qind(:)
            end if
            read(iunit, *) t(i)%nspin
            read(iunit, *) t(i)%nstate
            read(iunit, *) t(i)%nstate_spin
            read(iunit, *) t(i)%cstate
            read(iunit, *) check
            if (check) then
                allocate(t(i)%qe(t(i)%nstate))
                read(iunit, *) t(i)%qe(:)
            end if
            read(iunit, *) check
            if (check) then
                allocate(t(i)%qo(t(i)%nstate_spin(1)-1))
                read(iunit, *) t(i)%qo(:)
            end if
            read(iunit, *) t(i)%mnatom
            read(iunit, *) check
            if (check) then
                allocate(t(i)%mind(t(i)%mnatom))
                read(iunit, *) t(i)%mind(:)
            end if
            read(iunit, *) t(i)%me(2)
            read(iunit, *) check
            if (check) then
                allocate(t(i)%cwf(t(i)%nstate))
                read(iunit, *) t(i)%cwf(:)
            end if
            read(iunit, *) check
            if (check) then
                allocate(t(i)%prob(t(i)%nstate))
                read(iunit, *) t(i)%prob(:)
            end if
            read(iunit, *) check
            if (check) then
                allocate(t(i)%olap(t(i)%nstate, t(i)%nstate))
                read(iunit, *) t(i)%olap(:,:)
            end if
            read(iunit, *) check
            if (check) then
                allocate(t(i)%nadv(t(i)%ndim*t(i)%natom, t(i)%nstate, t(i)%nstate))
                read(iunit, *) t(i)%nadv(:, :, :)
            end if
            read(iunit, *) t(i)%pbcbox(6)
        end do
        close(iunit)
    end subroutine trajectory_read_backup


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Traj_writeheader
    !
    ! DESCRIPTION:
    !> @brief Initialize output files and write column headers.
    !----------------------------------------------------------------------------------------------
    subroutine traj_writeheader(t, popt, res_dir)
        class(trajtype) :: t
        logical, intent(in) :: popt(:) !< Print options.
        character(len=*), intent(in) :: res_dir ! Directory for writing results.
        integer :: ounit, i
        real(dp) :: time_fs

        time_fs = t%time * aut_fs
        if (popt(1)) then
            open(newunit=ounit, file=res_dir//'/energy.dat', action='write', position='append')
            write(ounit, '(a)', advance='no') '#'
            write(ounit, '(a12,1x,a4)', advance='no') 't,', 'cst,' 
            write(ounit, '(1x,2a18)', advance='no') 'etot,', 'epot,'
            write(ounit, '(1x)', advance='no') 
            do i = 1, size(t%qe)
                write(ounit, '(12x,a2,i3.3)', advance='no') 'qe', i
                if (i < size(t%qe)) write(ounit, '(a)', advance='no') ','
            end do
            write(ounit, *)
            close(ounit)

            if (t%mnatom > 0) then
                open(newunit=ounit, file=res_dir//'/mm.dat', action='write', position='append')
                write(ounit, '(a)', advance='no') '#'
                write(ounit, '(a12,1x)', advance='no') 't,'
                write(ounit, '(2x,3a18)', advance='no') 'mekin,', 'me1,', 'me2,'
                write(ounit, '(1x,5a18,a17)', advance='no') 'pbc1,', 'pbc2,', 'pbc3,', 'pbc4,', 'pbc5,', 'pbc6'
                write(ounit, *)
                close(ounit)
            end if
        end if
    end subroutine traj_writeheader



end module system_var
