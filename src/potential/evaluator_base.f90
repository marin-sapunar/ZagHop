module evaluator_base_mod
    use global_defs
    implicit none

    private
    public :: potential_evaluator
    public :: set_group_nstate
    public :: state_location

    type, abstract :: potential_evaluator
        integer :: n_state = 0
        integer :: n_group = 0
        integer, allocatable :: group_nstate(:)
        integer, allocatable :: group_i0(:)
        integer, allocatable :: group_mult(:) !< Multiplicity of each group.
        integer, allocatable :: group_block(:) !< Multiplicity block of each group. The groups
                                               !! holding the Ms components of the same states
                                               !! belong to the same block.
        logical :: soc_in_group = .false. !< SOC couples states within a group, so the
                                          !! group_adiabatic states are not spin-free.
        character(len=:), allocatable :: classical_coordinates !< Coordinates expected by
                                                               !! update_geometry.
    contains
        procedure(eval), deferred :: eval
        procedure(update_geometry), deferred :: update_geometry
        procedure(get_Hamiltonian), deferred :: get_Hamiltonian
        procedure(get_transformation), deferred :: get_transformation
        procedure(get_energy_single), deferred :: get_energy_single
        procedure(get_energy_all), deferred :: get_energy_all
        procedure(get_gradient), deferred :: get_gradient
        procedure(get_nadv), deferred :: get_nadv
        procedure :: set_group_nstate
        procedure :: locate_state
        procedure :: get_reference_geometry
        generic :: get_energy => get_energy_single, get_energy_all
    end type potential_evaluator


    !> Position of a state in the groups of states of a potential_evaluator.
    type :: state_location
        integer :: group = 0 !< Group containing the state.
        integer :: block = 0 !< Multiplicity block of the group.
        integer :: mult = 0 !< Multiplicity of the state.
        integer :: component = 0 !< Ms component (position of the group within its block).
        integer :: index = 0 !< Index of the state within its group.
    end type state_location

    
    abstract interface
        subroutine eval(self)
            import potential_evaluator
            class(potential_evaluator), intent(inout) :: self
        end subroutine eval
    end interface


    abstract interface
        pure function get_energy_single(self, basis, state) result(energy)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            character(len=*), intent(in) :: basis
            integer, intent(in) :: state
            real(dp) :: energy
        end function get_energy_single
    end interface


    abstract interface
        pure function get_energy_all(self, basis) result(energy)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            character(len=*), intent(in) :: basis
            real(dp), allocatable :: energy(:)
        end function get_energy_all
    end interface


    abstract interface
        !> Coordinates gives the coordinate system of the geometry. Gradients and nonadiabatic
        !! couplings are returned in the same coordinate system.
        subroutine update_geometry(self, geometry, coordinates)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(inout) :: self
            real(dp), intent(in) :: geometry(:, :)
            character(len=*), intent(in) :: coordinates
        end subroutine update_geometry
    end interface


    abstract interface
        function get_Hamiltonian(self, basis) result(H)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            character(len=*), intent(in) :: basis
            real(dp), allocatable :: H(:, :)
        end function get_Hamiltonian
    end interface


    abstract interface
        function get_transformation(self, from_basis, to_basis) result(trans)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            character(len=*), intent(in) :: from_basis, to_basis
            real(dp), allocatable :: trans(:, :)
        end function get_transformation
    end interface


    abstract interface
        function get_gradient(self, basis, istate) result(grad)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            character(len=*), intent(in) :: basis
            integer, intent(in) :: istate
            real(dp), allocatable :: grad(:, :)
        end function get_gradient
    end interface


    abstract interface
        function get_nadv(self, basis, istate1, istate2) result(nadv)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            character(len=*), intent(in) :: basis
            integer, intent(in) :: istate1, istate2
            real(dp), allocatable :: nadv(:)
        end function get_nadv
    end interface


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: set_group_nstate
    !> @brief Set the groups of states from the number of states of each multiplicity.
    !> @details
    !! Each block i of nstate(i) states with multiplicity multiplicity(i) gives multiplicity(i)
    !! consecutive groups of nstate(i) states, one for each Ms component. Blocks without states
    !! are skipped, but keep their index in group_block.
    !----------------------------------------------------------------------------------------------
    subroutine set_group_nstate(self, nstate, multiplicity)
        class(potential_evaluator), intent(inout) :: self
        integer, intent(in) :: nstate(:)
        integer, intent(in) :: multiplicity(:)
        integer :: i, j, ci, cgroup

        self%n_group = sum(multiplicity, mask=(nstate > 0))
        allocate(self%group_nstate(self%n_group))
        allocate(self%group_i0(self%n_group))
        allocate(self%group_mult(self%n_group))
        allocate(self%group_block(self%n_group))
        cgroup = 1
        ci = 1
        do i = 1, size(nstate)
            if (nstate(i) == 0) cycle
            do j = 1, multiplicity(i)
                self%group_nstate(cgroup) = nstate(i)
                self%group_i0(cgroup) = ci
                self%group_mult(cgroup) = multiplicity(i)
                self%group_block(cgroup) = i
                cgroup = cgroup + 1
                ci = ci + nstate(i)
            end do
        end do
        self%n_state = ci - 1
    end subroutine set_group_nstate


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: locate_state
    !> @brief Group, multiplicity block, multiplicity, Ms component and index within the group of
    !! a state.
    !----------------------------------------------------------------------------------------------
    pure function locate_state(self, state) result(loc)
        class(potential_evaluator), intent(in) :: self
        integer, intent(in) :: state
        type(state_location) :: loc
        integer :: g

        if (state < 1 .or. state > self%n_state) then
            error stop 'Error in evaluator_base_mod, locate_state function. State out of range.'
        end if
        do g = 1, self%n_group
            if (state < self%group_i0(g) + self%group_nstate(g)) exit
        end do
        loc%group = g
        loc%block = self%group_block(g)
        loc%mult = self%group_mult(g)
        loc%component = g - findloc(self%group_block, loc%block, dim=1) + 1
        loc%index = state - self%group_i0(g) + 1
    end function locate_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: get_reference_geometry
    !> @brief Reference geometry which defines the orientation of the potential.
    !> @details
    !! The geometry is not allocated (default) if the potential is invariant to overall rotation
    !! of the system given in the coordinates.
    !----------------------------------------------------------------------------------------------
    subroutine get_reference_geometry(self, coordinates, geom)
        class(potential_evaluator), intent(in) :: self
        character(len=*), intent(in) :: coordinates
        real(dp), allocatable, intent(out) :: geom(:, :)
    end subroutine get_reference_geometry


end module evaluator_base_mod