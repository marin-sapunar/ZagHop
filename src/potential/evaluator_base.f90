!--------------------------------------------------------------------------------------------------
! MODULE: evaluator_base_mod
!
! DESCRIPTION:
!> @brief Base class of all potentials.
!> @details
!! The states of a potential are given in three bases:
!! - diabatic - For model systems, the basis in which the Hamiltonian is provided.
!! - mch - Eigenstates of the spin-free Hamiltonian (molecular Coulomb Hamiltonian) within each
!!   multiplicity. One state for each spin-free state.
!! - adiabatic - Eigenstates of the full Hamiltonian (with spin-orbit coupling).
!--------------------------------------------------------------------------------------------------
module evaluator_base_mod
    use global_defs
    use state_overlap_mod, only : state_overlap, request_id
    use state_layout_mod, only : state_layout, spin_manifold, state_location
    implicit none

    private
    public :: potential_evaluator
    public :: eval_request
    public :: new_request
    public :: state_layout
    public :: spin_manifold
    public :: state_location
    public :: add_state_layout
    public :: soc_group_tol
    public :: evaluator_set_model_states
    public :: evaluator_save_state
    public :: evaluator_load_state


    real(dp), parameter :: soc_group_tol = 1.0e-10_dp !< Largest spin-orbit coupling within a group
                                                      !! of states without a warning.


    !> Derivatives that an evaluation of the potential must provide, besides the energies.
    !> @details
    !! Model potentials compute all derivatives at every point and ignore requests. Potentials
    !! computed by external programs compute only what is requested.
    type :: eval_request
        character(len=:), allocatable :: basis !< Basis of the states.
        integer, allocatable :: gradient(:) !< States whose gradients are needed.
        integer, allocatable :: nadv(:, :) !< Pairs of states (2, n_pair) whose nonadiabatic
                                           !! coupling vectors are needed.
    end type eval_request


    !> Base class of all potentials.
    !> @details
    !! Specific evaluators implement update_geometry and the evaluation of the potential at the
    !! current geometry, and fill the variables of the base classes. The base classes do the rest
    !!
    !! An evaluator holds the quantities of a single point. The only dependence on the path is
    !! through the previous point passed to eval: the phases of the states (and any other sign
    !! convention) follow those of the previous point, and link holds the overlaps between the
    !! states of the two points. Each evaluation gets a new id, which identifies the point in
    !! link.
    type, abstract :: potential_evaluator
        type(state_layout) :: model_states !< States of the model of the potential (diabatic
                                           !! basis).
        type(state_layout) :: states !< States of the dynamics (adiabatic and mch bases): all
                                     !! states of the model, or a subset of them.
        integer, allocatable :: sel(:) !< MCH state of the model, with all spin components in the
                                       !! order of model_states, selected for each state of
                                       !! states.
        integer, allocatable :: sel_mch(:) !< MCH state of the model for each mch state of the
                                           !! dynamics.
        integer, allocatable :: uns(:) !< MCH states of the model, with all spin components,
                                       !! which are not selected.
        character(len=:), allocatable :: classical_coordinates !< Coordinates expected by
                                                               !! update_geometry.
        character(len=12) :: component_order = 'by_component' !< Order of the spin components in
        !! a full Hamiltonian given by the source of the potential:
        !! - by_component - All states of the first component, then of the second, ... (ZagHop).
        !! - by_state - All components of the first state, then of the second, ... (Quantics).
        integer, allocatable :: state_map(:) !< Index in the source of each state of the model.
        integer :: id = 0 !< Id of the evaluated point. 0 before the first evaluation.
        type(state_overlap) :: link !< Overlaps between the states of the previous point and of
                                   !! this point. Not set without a previous point.
        type(eval_request) :: request !< Derivatives needed from the next evaluation.
    contains
        procedure(eval), deferred :: eval
        procedure(update_geometry), deferred :: update_geometry
        procedure :: new_point
        procedure :: set_request
        procedure :: complete
        procedure :: get_derivative_target_coord
        procedure :: set_cartesian_frame
        procedure :: check_atoms
        procedure(get_Hamiltonian), deferred :: get_Hamiltonian
        procedure(get_transformation), deferred :: get_transformation
        procedure(get_energy_single), deferred :: get_energy_single
        procedure(get_energy_all), deferred :: get_energy_all
        procedure(get_gradient), deferred :: get_gradient
        procedure(get_nadv), deferred :: get_nadv
        procedure(get_soc2), deferred :: get_soc2
        procedure :: set_model_states => evaluator_set_model_states
        procedure :: select_states
        procedure :: basis_size
        procedure :: diabatic_from_source
        procedure :: get_reference_geometry
        procedure :: default_reference
        procedure :: print_summary
        procedure :: add_state_layout
        procedure :: save_state => evaluator_save_state
        procedure :: load_state => evaluator_load_state
        generic :: get_energy => get_energy_single, get_energy_all
    end type potential_evaluator


    abstract interface
        !> Evaluate the potential at the current geometry. prev is the evaluator of the previous
        !! point of the trajectory, if there is one. The phases of the states then follow those of
        !! prev, and link holds the overlaps between the states of prev and of this point.
        !! Implementations start with new_point. prev must not be the evaluator itself.
        subroutine eval(self, prev)
            import potential_evaluator
            class(potential_evaluator), intent(inout) :: self
            class(potential_evaluator), intent(in), optional :: prev
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
        !> Set the geometry at which the potential is evaluated next. Coordinates gives the
        !! coordinate system of the geometry. Anything that get_derivative_target_coord needs to
        !! convert derivatives to these coordinates must be stored here. Gradients and
        !! nonadiabatic couplings are returned in the coordinates of this geometry.
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
        !> Gradient of a state, flat, in the coordinates of the last geometry passed to
        !! update_geometry. Callers reshape it to the shape of their geometry.
        function get_gradient(self, basis, istate) result(grad)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            character(len=*), intent(in) :: basis
            integer, intent(in) :: istate
            real(dp), allocatable :: grad(:)
        end function get_gradient
    end interface


    abstract interface
        !> Nonadiabatic coupling vector between two states, flat, in the coordinates of the last
        !! geometry passed to update_geometry.
        function get_nadv(self, basis, istate1, istate2) result(nadv)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            character(len=*), intent(in) :: basis
            integer, intent(in) :: istate1, istate2
            real(dp), allocatable :: nadv(:)
        end function get_nadv
    end interface


    abstract interface
        !> Squared magnitudes of the spin-orbit couplings between all states of a basis. In the
        !! 'mch' basis, the squares are summed over the components of each pair of
        !! spin-free states.
        function get_soc2(self, basis) result(s2)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(in) :: self
            character(len=*), intent(in) :: basis
            real(dp), allocatable :: s2(:, :)
        end function get_soc2
    end interface


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: new_point
    !> @brief Start the evaluation of a new point: give it a new id and remove the link to the
    !! previous point.
    !> @details
    !! The link must be removed because the evaluator of a new point of the trajectory is a copy
    !! of the evaluator of the previous point, with its link.
    !----------------------------------------------------------------------------------------------
    subroutine new_point(self)
        class(potential_evaluator), intent(inout) :: self

        self%id = request_id()
        call self%link%clear()
    end subroutine new_point


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: new_request
    !> @brief Request for the gradients of some states and the nonadiabatic coupling vectors
    !! between some pairs of states of a basis.
    !----------------------------------------------------------------------------------------------
    pure function new_request(basis, gradient, nadv) result(req)
        character(len=*), intent(in) :: basis !< Basis of the states.
        integer, intent(in), optional :: gradient(:) !< States whose gradients are needed.
        integer, intent(in), optional :: nadv(:, :) !< Pairs of states (2, n_pair).
        type(eval_request) :: req

        req%basis = basis
        if (present(gradient)) then
            req%gradient = gradient
        else
            allocate(req%gradient(0))
        end if
        if (present(nadv)) then
            req%nadv = nadv
        else
            allocate(req%nadv(2, 0))
        end if
    end function new_request


    !> Set the derivatives that the next evaluation must provide (see eval_request).
    subroutine set_request(self, req)
        class(potential_evaluator), intent(inout) :: self
        type(eval_request), intent(in) :: req

        self%request = req
    end subroutine set_request


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: complete
    !> @brief Make the requested derivatives available at the evaluated point.
    !> @details
    !! Used when more derivatives are needed at a point after its evaluation, for example the
    !! gradient of a new active state after a hop. Nothing is done in the base class, for
    !! potentials which compute all derivatives at every point. Potentials computed by external
    !! programs override this to compute the derivatives which are missing.
    !----------------------------------------------------------------------------------------------
    subroutine complete(self, req)
        class(potential_evaluator), intent(inout) :: self
        type(eval_request), intent(in) :: req
    end subroutine complete


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: get_derivative_target_coord
    !> @brief Convert a derivative with respect to the coordinates of the potential to the
    !! coordinates of the last geometry passed to update_geometry.
    !> @details
    !! The result is flat. In the base class the two coordinate systems are the same. An
    !! evaluator that accepts geometries in other coordinates overrides this function. Everything
    !! it needs must be stored in update_geometry.
    !----------------------------------------------------------------------------------------------
    function get_derivative_target_coord(self, d) result(dr)
        class(potential_evaluator), intent(in) :: self
        real(dp), intent(in) :: d(:) !< Derivative with respect to the coordinates of the potential.
        real(dp), allocatable :: dr(:)

        dr = d
    end function get_derivative_target_coord


    !> Set the frame in which Cartesian geometries are transformed to the coordinates of the
    !! potential. Not available in the base class.
    subroutine set_cartesian_frame(self, frame)
        class(potential_evaluator), intent(inout) :: self
        character(len=*), intent(in) :: frame

        write(stderr, *) 'Error in evaluator_base_mod, set_cartesian_frame subroutine.'
        write(stderr, *) '  Option cartesian_frame is not available for this potential: ', frame
        stop
    end subroutine set_cartesian_frame


    !> Check that the atoms of the system match the potential. Nothing is checked in the base
    !! class.
    subroutine check_atoms(self, sym, mass)
        class(potential_evaluator), intent(in) :: self
        character(len=*), intent(in) :: sym(:) !< Atom symbols.
        real(dp), intent(in) :: mass(:) !< Atom masses (amu).
    end subroutine check_atoms


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: evaluator_set_model_states
    !> @brief Set the states of the model from the number of states of each multiplicity.
    !> @details
    !! See state_layout for nstate, multiplicity and spin_states. state_map is built from
    !! component_order. The dynamics runs on all states of the model, until select_states is
    !! called. Evaluators whose source fixes the number of states override this to check it.
    !----------------------------------------------------------------------------------------------
    subroutine evaluator_set_model_states(self, nstate, multiplicity, spin_states)
        class(potential_evaluator), intent(inout) :: self
        integer, intent(in) :: nstate(:) !< Number of spin-free states of each multiplicity.
        integer, intent(in) :: multiplicity(:) !< Multiplicity of each block of states.
        character(len=*), intent(in) :: spin_states
        integer :: m, ci, c, k

        call self%model_states%set(nstate, multiplicity, spin_states)

        ! Index in the source of each state. With one component, both orders are the same.
        if (allocated(self%state_map)) deallocate(self%state_map)
        allocate(self%state_map(self%model_states%n_state))
        do m = 1, size(self%model_states%manifold)
            associate(mf => self%model_states%manifold(m))
                do c = 1, mf%n_comp
                    do k = 1, mf%n
                        ci = mf%i0 + (c - 1) * mf%n + k - 1
                        select case(self%component_order)
                        case('by_component')
                            self%state_map(ci) = ci
                        case('by_state')
                            self%state_map(ci) = mf%i0 + (k - 1) * mf%n_comp + c - 1
                        case default
                            write(stderr, *) 'Error in evaluator_base_mod, set_model_states ',    &
                            &                'subroutine.'
                            write(stderr, *) '  Unknown component_order: ', trim(self%component_order)
                            stop
                        end select
                    end do
                end do
            end associate
        end do

        self%states = self%model_states
        call self%model_states%subset_index(self%states, self%sel, self%sel_mch, self%uns)
    end subroutine evaluator_set_model_states


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: select_states
    !> @brief Choose the states of the dynamics: the nstate(i) lowest MCH states of each
    !! multiplicity of the model.
    !> @details
    !! The states of the model must be set. A multiplicity which is not given has no states in the
    !! dynamics. Stops if the states are not a subset of the states of the model (see subset of
    !! state_layout). Must be called before the first evaluation of the potential.
    !----------------------------------------------------------------------------------------------
    subroutine select_states(self, nstate, multiplicity)
        class(potential_evaluator), intent(inout) :: self
        integer, intent(in) :: nstate(:) !< Number of spin-free states of each multiplicity.
        integer, intent(in) :: multiplicity(:) !< Multiplicity of each block of states.

        if (.not. self%model_states%is_set()) then
            write(stderr, *) 'Error in evaluator_base_mod, select_states subroutine.'
            write(stderr, *) '  The states of the model are not set.'
            stop
        end if
        self%states = self%model_states%subset(nstate, multiplicity)
        call self%model_states%subset_index(self%states, self%sel, self%sel_mch, self%uns)
    end subroutine select_states


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: basis_size
    !> @brief Number of states of a basis.
    !> @details
    !! The 'diabatic' basis has all states of the model. The 'mch' basis has one state for each
    !! spin-free state of the dynamics, and the other bases one state for each state of the
    !! dynamics.
    !----------------------------------------------------------------------------------------------
    pure function basis_size(self, basis) result(n)
        class(potential_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer :: n

        select case(basis)
        case('diabatic')
            n = self%model_states%n_state
        case('mch')
            n = self%states%n_sf
        case default
            n = self%states%n_state
        end select
    end function basis_size


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: diabatic_from_source
    !> @brief Diabatic state given by its index in the source of the potential.
    !> @details
    !! The index in the source follows component_order (see state_map). The result is 0 if the
    !! source has no state with this index.
    !----------------------------------------------------------------------------------------------
    pure function diabatic_from_source(self, i_source) result(state)
        class(potential_evaluator), intent(in) :: self
        integer, intent(in) :: i_source !< Index of the state in the source.
        integer :: state

        if (allocated(self%state_map)) then
            state = findloc(self%state_map, i_source, dim=1)
        else if (i_source >= 1 .and. i_source <= self%model_states%n_state) then
            state = i_source
        else
            state = 0
        end if
    end function diabatic_from_source


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


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: default_reference
    !> @brief Default reference geometry of the potential, in its classical coordinates.
    !> @details
    !! Used when a quantity must be computed at the same geometry for all trajectories, and no
    !! reference geometry is given in the input. The geometry is not allocated (default) if the
    !! potential has no natural reference.
    !----------------------------------------------------------------------------------------------
    subroutine default_reference(self, geom)
        class(potential_evaluator), intent(in) :: self
        real(dp), allocatable, intent(out) :: geom(:, :)
    end subroutine default_reference


    !> Print a summary of the potential at the end of the trajectory. Nothing is printed in the
    !! base class.
    subroutine print_summary(self, unit)
        class(potential_evaluator), intent(in) :: self
        integer, intent(in) :: unit
    end subroutine print_summary


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: add_state_layout
    !> @brief Add the layout of the states of the potential to a JSON object.
    !> @details
    !! Adds the number of states of each basis (n_state), and spin_states and the spin manifolds
    !! of the states of the dynamics (see add_to_json of state_layout). The object 'model' holds
    !! spin_states and the spin manifolds of the states of the model, which are the states of the
    !! diabatic basis. Evaluators that are built from another potential override this to describe
    !! it as well.
    !----------------------------------------------------------------------------------------------
    subroutine add_state_layout(self, json, obj)
        use json_module, only : json_core, json_value
        class(potential_evaluator), intent(in) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the layout is added.
        type(json_value), pointer :: nst, model

        call json%create_object(nst, 'n_state')
        call json%add(nst, 'diabatic', self%basis_size('diabatic'))
        call json%add(nst, 'adiabatic', self%basis_size('adiabatic'))
        call json%add(nst, 'mch', self%basis_size('mch'))
        call json%add(obj, nst)
        call self%states%add_to_json(json, obj)
        call json%create_object(model, 'model')
        call self%model_states%add_to_json(json, model)
        call json%add(obj, model)
    end subroutine add_state_layout


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: evaluator_save_state
    !> @brief Add the state of the evaluated point to a JSON object, for a backup.
    !> @details
    !! Only what depends on the point is saved: the quantities evaluated at it, its id and its
    !! link to the previous point. The setup (model, spin manifolds, options) is built from the
    !! input when the backup is loaded. Each type adds its own object to obj (here 'base', with
    !! the numbers of states of the dynamics and of the model for checking), so a backup of
    !! another type of potential is detected
    !! on loading. Extensions override this and call the subroutine of their parent first.
    !----------------------------------------------------------------------------------------------
    subroutine evaluator_save_state(self, json, obj)
        use json_module, only : json_core, json_value
        class(potential_evaluator), intent(in) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the state is added.
        type(json_value), pointer :: base, link

        call json%create_object(base, 'base')
        call json%add(base, 'n_state', self%states%n_state)
        call json%add(base, 'n_sf', self%states%n_sf)
        call json%add(base, 'n_model_state', self%model_states%n_state)
        call json%add(base, 'n_model_sf', self%model_states%n_sf)
        call json%add(base, 'spin_states', self%states%spin_states)
        call json%add(base, 'id', self%id)
        call json%create_object(link, 'link')
        call self%link%save_state(json, link)
        call json%add(base, link)
        call json%add(obj, base)
    end subroutine evaluator_save_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: evaluator_load_state
    !> @brief Restore the state of an evaluated point saved by save_state.
    !> @details
    !! The evaluator must be set up from the same input as the one that was saved. Stops if the
    !! numbers of states differ.
    !----------------------------------------------------------------------------------------------
    subroutine evaluator_load_state(self, json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_get_required, json_get_child, json_check_value
        class(potential_evaluator), intent(inout) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object holding the state.
        type(json_value), pointer :: base

        base => json_get_child(json, obj, 'base')
        call json_check_value(json, base, 'n_state', self%states%n_state)
        call json_check_value(json, base, 'n_sf', self%states%n_sf)
        call json_check_value(json, base, 'n_model_state', self%model_states%n_state)
        call json_check_value(json, base, 'n_model_sf', self%model_states%n_sf)
        call json_check_value(json, base, 'spin_states', self%states%spin_states)
        call json_get_required(json, base, 'id', self%id)
        call self%link%load_state(json, json_get_child(json, base, 'link'))
    end subroutine evaluator_load_state


end module evaluator_base_mod
