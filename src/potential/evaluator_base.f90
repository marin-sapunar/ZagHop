!--------------------------------------------------------------------------------------------------
! MODULE: evaluator_base_mod
!
! DESCRIPTION:
!> @brief Base class of all potentials.
!> @details
!! The states of a potential are given in five bases, in two spaces (see basis_mod): 'diabatic',
!! 'mch' and 'adiabatic' with the spin components of the states, and 'diabatic_sf' and 'mch_sf'
!! with one state for each multiplet. The spin-orbit couplings are the difference between the
!! Hamiltonians of the two spaces (get_soc2).
!!
!! An evaluator relates its sets of states by three quantities, which it fills when a point is
!! evaluated:
!! - diab - The spin-free MCH states of the dynamics in the spin-free diabatic states
!!   ('diabatic_sf' -> 'mch_sf'). The same vectors hold for every spin component. Not set for
!!   potentials without diabatic states (has_diabatic).
!! - adiab_vec - The adiabatic states in the MCH states with their spin components
!!   ('mch' -> 'adiabatic'): U of Mai, Marquetand and Gonzalez, WIREs Comput. Mol. Sci. 8, e1370
!!   (2018), doi: 10.1002/wcms.1370, eqs. 7 and 8.
!! - link - The overlaps between the states of one basis at the previous point and at this point.
!!
!! All transformations (get_transformation, transform) are built from these in the base class and
!! follow one convention. T has a row for each state of from_basis and a column for each state of
!! to_basis: T(i, j) = <from_i|to_j>, so its columns are the states of to_basis in the states of
!! from_basis. Coefficients of a wave function transform as c_to = T^T c_from, and operators as
!! M_to = T^T M_from T. All transformations are real. Transformations exist only between bases of
!! the same space (can_transform).
!--------------------------------------------------------------------------------------------------
module evaluator_base_mod
    use global_defs
    use basis_mod, only : basis_known, basis_is_spin_free, basis_full, basis_like
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
    public :: overlap_warn_tol
    public :: evaluator_set_model_states
    public :: evaluator_save_state
    public :: evaluator_load_state


    real(dp), parameter :: soc_group_tol = 1.0e-10_dp !< Largest spin-orbit coupling within a group
                                                      !! of states without a warning.
    real(dp), parameter :: overlap_warn_tol = 0.9_dp !< Smallest singular value of the overlaps
    !! between the states of two points without a warning.


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
    !! link. The states of the previous point are reached through get_transformation and transform
    !! with the id of that point.
    type, abstract :: potential_evaluator
        type(state_layout) :: model_states !< States of the model of the potential (diabatic
                                           !! basis).
        type(state_layout) :: states !< States of the dynamics (adiabatic, mch and mch_sf bases):
                                     !! all states of the model, or a subset of them.
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
        integer :: id_diab = 0 !< Id of the diabatic states, the same at all points. 0 without
                               !! diabatic states.
        type(state_overlap) :: diab !< Spin-free MCH states of the dynamics in the spin-free
        !! diabatic states: bra 'diabatic_sf' (id_diab), ket 'mch_sf' (id), (n_sf of the model,
        !! n_sf). Its t is the transformation. Not set without diabatic states.
        real(dp), allocatable :: adiab_vec(:, :) !< Adiabatic states in the MCH states with their
        !! spin components, with the phases of the adiabatic states (n_state, n_state).
        type(state_overlap) :: link !< Overlaps between the states of the previous point and of
        !! this point, both in the basis link_basis. Not set without a previous point.
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
        procedure :: has_diabatic
        procedure :: can_transform
        procedure :: get_transformation
        procedure :: transform_vector
        procedure :: transform_operator
        generic :: transform => transform_vector, transform_operator
        procedure :: previous_point
        procedure :: link_basis
        procedure :: set_link
        procedure(get_energy_single), deferred :: get_energy_single
        procedure(get_energy_all), deferred :: get_energy_all
        procedure(get_gradient), deferred :: get_gradient
        procedure(get_nadv), deferred :: get_nadv
        procedure :: get_soc2
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
    !! Potentials with diabatic states get the id of these states here.
    !----------------------------------------------------------------------------------------------
    subroutine evaluator_set_model_states(self, nstate, multiplicity, spin_states)
        class(potential_evaluator), intent(inout) :: self
        integer, intent(in) :: nstate(:) !< Number of spin-free states of each multiplicity.
        integer, intent(in) :: multiplicity(:) !< Multiplicity of each block of states.
        character(len=*), intent(in) :: spin_states
        integer :: m, ci, c, k

        call self%model_states%set(nstate, multiplicity, spin_states)
        if (self%has_diabatic() .and. self%id_diab == 0) self%id_diab = request_id()

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
    !! The diabatic bases have the states of the model: all of them ('diabatic') or one for each
    !! spin-free state ('diabatic_sf'). The 'mch_sf' basis has one state for each spin-free state
    !! of the dynamics, and the other bases one state for each state of the dynamics.
    !----------------------------------------------------------------------------------------------
    pure function basis_size(self, basis) result(n)
        class(potential_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer :: n

        select case(basis)
        case('diabatic')
            n = self%model_states%n_state
        case('diabatic_sf')
            n = self%model_states%n_sf
        case('mch_sf')
            n = self%states%n_sf
        case default
            n = self%states%n_state
        end select
    end function basis_size


    !> Whether the potential has diabatic states (the 'diabatic' and 'diabatic_sf' bases). Not in
    !! the base class. Needs no evaluated point.
    pure function has_diabatic(self) result(has)
        class(potential_evaluator), intent(in) :: self
        logical :: has

        has = .false.
    end function has_diabatic


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: can_transform
    !> @brief Whether get_transformation and transform are available from from_basis to to_basis.
    !> @details
    !! Both must be bases of the same space (see basis_mod), and a diabatic basis needs a
    !! potential with diabatic states. Needs no evaluated point, so it can check the input before
    !! the first evaluation.
    !----------------------------------------------------------------------------------------------
    pure function can_transform(self, from_basis, to_basis) result(available)
        class(potential_evaluator), intent(in) :: self
        character(len=*), intent(in) :: from_basis
        character(len=*), intent(in) :: to_basis
        logical :: available

        available = basis_known(from_basis) .and. basis_known(to_basis)
        if (.not. available) return
        available = (basis_is_spin_free(from_basis) .eqv. basis_is_spin_free(to_basis))
        if (.not. available) return
        if (basis_full(from_basis) == 'diabatic' .or. basis_full(to_basis) == 'diabatic') then
            available = self%has_diabatic()
        end if
    end function can_transform


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: get_transformation
    !> @brief Matrix T(i, j) = <from_i|to_j>, whose columns are the states of to_basis in the
    !! states of from_basis (see the module description).
    !> @details
    !! Without from_point and to_point, both sets of states are those of this point, and the
    !! transformation is built from diab and adiab_vec:
    !! - 'diabatic_sf' -> 'mch_sf': diab%t.
    !! - 'diabatic' -> 'mch': diab%t for every spin component (P of evaluator_diabatic_mod).
    !! - 'mch' -> 'adiabatic': adiab_vec.
    !! - 'diabatic' -> 'adiabatic': the product of the two.
    !! The opposite directions are the transposes, and equal bases give the unit matrix. A
    !! transformation from or to a diabatic basis is not square when the dynamics runs on a subset
    !! of the states of the model.
    !!
    !! from_point (to_point) is the id of the point of the states of from_basis (to_basis). An id
    !! other than that of this point must be the id of the previous point, and the basis on that
    !! side the basis of link: the transformation is then the orthogonalised overlap matrix
    !! link%t, followed by the transformation from the basis of link to the other basis at this
    !! point. For from_basis = to_basis, this is the transformation between the locally diabatic
    !! basis of a step and the basis at its end (see sh_ld_mod).
    !!
    !! Stops if the transformation is not available (see can_transform), if the potential was not
    !! evaluated, or if link is not between the requested point and this point.
    !----------------------------------------------------------------------------------------------
    function get_transformation(self, from_basis, to_basis, from_point, to_point) result(trans)
        class(potential_evaluator), intent(in) :: self
        character(len=*), intent(in) :: from_basis
        character(len=*), intent(in) :: to_basis
        integer, intent(in), optional :: from_point !< Id of the point of the from_basis states.
        integer, intent(in), optional :: to_point !< Id of the point of the to_basis states.
        real(dp), allocatable :: trans(:, :)
        logical :: from_prev, to_prev

        from_prev = .false.
        to_prev = .false.
        if (present(from_point)) from_prev = (from_point /= self%id)
        if (present(to_point)) to_prev = (to_point /= self%id)
        if (from_prev .and. to_prev) then
            write(stderr, *) 'Error in evaluator_base_mod, get_transformation function.'
            write(stderr, *) '  Both sets of states are those of another point.'
            stop
        else if (from_prev) then
            trans = previous_transformation(self, from_basis, to_basis, from_point)
        else if (to_prev) then
            trans = transpose(previous_transformation(self, to_basis, from_basis, to_point))
        else
            trans = point_transformation(self, from_basis, to_basis)
        end if
    end function get_transformation


    !> Transformation between two bases at this point (see get_transformation).
    function point_transformation(self, from_basis, to_basis) result(trans)
        use matrix_mod, only : unit_mat
        class(potential_evaluator), intent(in) :: self
        character(len=*), intent(in) :: from_basis
        character(len=*), intent(in) :: to_basis
        real(dp), allocatable :: trans(:, :)

        if (.not. self%can_transform(from_basis, to_basis)) call not_available()
        if (from_basis == to_basis) then
            trans = unit_mat(self%basis_size(to_basis))
            return
        end if

        select case(from_basis//' '//to_basis)
        case('diabatic_sf mch_sf')
            call check_evaluated(.true., .false.)
            trans = self%diab%t
        case('mch_sf diabatic_sf')
            call check_evaluated(.true., .false.)
            trans = transpose(self%diab%t)
        case('diabatic mch')
            call check_evaluated(.true., .false.)
            trans = self%model_states%expand(self%diab%t, self%states)
        case('mch diabatic')
            call check_evaluated(.true., .false.)
            trans = transpose(self%model_states%expand(self%diab%t, self%states))
        case('mch adiabatic')
            call check_evaluated(.false., .true.)
            trans = self%adiab_vec
        case('adiabatic mch')
            call check_evaluated(.false., .true.)
            trans = transpose(self%adiab_vec)
        case('diabatic adiabatic')
            call check_evaluated(.true., .true.)
            trans = matmul(self%model_states%expand(self%diab%t, self%states), self%adiab_vec)
        case('adiabatic diabatic')
            call check_evaluated(.true., .true.)
            trans = transpose(matmul(self%model_states%expand(self%diab%t, self%states),          &
            &                        self%adiab_vec))
        end select

    contains

        subroutine not_available()
            write(stderr, *) 'Error in evaluator_base_mod, get_transformation function.'
            if (.not. basis_known(from_basis)) then
                write(stderr, *) '  Unknown basis: ', from_basis
            else if (.not. basis_known(to_basis)) then
                write(stderr, *) '  Unknown basis: ', to_basis
            else if (basis_is_spin_free(from_basis) .neqv. basis_is_spin_free(to_basis)) then
                write(stderr, *) '  No transformation from the ', from_basis, ' to the ',         &
                &                to_basis, ' basis.'
                write(stderr, *) '  One has the spin components of the states and the other ',   &
                &                'one state for each multiplet.'
            else
                write(stderr, *) '  No transformation from the ', from_basis, ' to the ',         &
                &                to_basis, ' basis.'
                write(stderr, *) '  The diabatic basis is not available for this potential.'
            end if
            stop
        end subroutine not_available

        subroutine check_evaluated(need_diab, need_vec)
            logical, intent(in) :: need_diab
            logical, intent(in) :: need_vec
            logical :: ok

            ok = .true.
            if (need_diab) ok = self%diab%valid()
            if (need_vec) ok = ok .and. allocated(self%adiab_vec)
            if (ok) return
            write(stderr, *) 'Error in evaluator_base_mod, get_transformation function.'
            write(stderr, *) '  The potential was not evaluated.'
            stop
        end subroutine check_evaluated

    end function point_transformation


    !> Transformation from the states of prev_basis at the previous point, with id point, to the
    !! states of basis at this point (see get_transformation).
    function previous_transformation(self, prev_basis, basis, point) result(trans)
        class(potential_evaluator), intent(in) :: self
        character(len=*), intent(in) :: prev_basis
        character(len=*), intent(in) :: basis
        integer, intent(in) :: point
        real(dp), allocatable :: trans(:, :)

        if (.not. self%link%matches(point, prev_basis, self%id, prev_basis)) then
            write(stderr, *) 'Error in evaluator_base_mod, get_transformation function.'
            write(stderr, *) '  The overlaps of the potential do not link the requested point ', &
            &                'to the current point.'
            write(stderr, '(3x,a,i0,3a,i0,a)') 'Requested: point ', point, ', ', prev_basis,     &
            &   ' basis -> point ', self%id, '.'
            if (self%link%valid()) then
                write(stderr, '(3x,a,i0,3a,i0,a)') 'Stored: point ', self%link%id_bra, ', ',     &
                &   self%link%basis_bra, ' basis -> point ', self%link%id_ket, '.'
            else
                write(stderr, *) '  No overlaps with a previous point are stored.'
            end if
            stop
        end if
        if (basis == prev_basis) then
            trans = self%link%t
        else
            trans = matmul(self%link%t, point_transformation(self, prev_basis, basis))
        end if
    end function previous_transformation


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: transform_vector
    !> @brief Coefficients of a wave function in the states of to_basis, c_to = T^T c, from its
    !! coefficients c in the states of from_basis (T from get_transformation).
    !> @details
    !! With the dynamics on a subset of the states of the model, the transformation from the
    !! diabatic basis is not square, and the result is the projection onto the states of the
    !! dynamics. from_point and to_point are those of get_transformation. Stops if the
    !! transformation is not available or c has the wrong size.
    !----------------------------------------------------------------------------------------------
    function transform_vector(self, from_basis, to_basis, c, from_point, to_point) result(c_to)
        class(potential_evaluator), intent(in) :: self
        character(len=*), intent(in) :: from_basis
        character(len=*), intent(in) :: to_basis
        complex(dp), intent(in) :: c(:) !< Coefficients in the states of from_basis.
        integer, intent(in), optional :: from_point !< Id of the point of the from_basis states.
        integer, intent(in), optional :: to_point !< Id of the point of the to_basis states.
        complex(dp), allocatable :: c_to(:)
        real(dp), allocatable :: t(:, :)

        t = self%get_transformation(from_basis, to_basis, from_point, to_point)
        call check_transform_size(from_basis, to_basis, shape(c), size(t, 1), 'transform_vector')
        c_to = matmul(transpose(t), c)
    end function transform_vector


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: transform_operator
    !> @brief Matrix of an operator in the states of to_basis, T^T mat T, from its matrix mat in
    !! the states of from_basis (T from get_transformation).
    !> @details
    !! A rank 2 argument is always an operator. A set of coefficient vectors is transformed one at
    !! a time with transform_vector. from_point and to_point are those of get_transformation.
    !! Stops if the transformation is not available or mat has the wrong shape.
    !----------------------------------------------------------------------------------------------
    function transform_operator(self, from_basis, to_basis, mat, from_point, to_point)           &
    &   result(mat_to)
        class(potential_evaluator), intent(in) :: self
        character(len=*), intent(in) :: from_basis
        character(len=*), intent(in) :: to_basis
        real(dp), intent(in) :: mat(:, :) !< Operator in the states of from_basis.
        integer, intent(in), optional :: from_point !< Id of the point of the from_basis states.
        integer, intent(in), optional :: to_point !< Id of the point of the to_basis states.
        real(dp), allocatable :: mat_to(:, :)
        real(dp), allocatable :: t(:, :)

        t = self%get_transformation(from_basis, to_basis, from_point, to_point)
        call check_transform_size(from_basis, to_basis, shape(mat), size(t, 1),                  &
        &                         'transform_operator')
        mat_to = matmul(matmul(transpose(t), mat), t)
    end function transform_operator


    !> Stop if an argument of transform with shape arg_shape is not over the n_from states of
    !! from_basis.
    subroutine check_transform_size(from_basis, to_basis, arg_shape, n_from, caller)
        character(len=*), intent(in) :: from_basis
        character(len=*), intent(in) :: to_basis
        integer, intent(in) :: arg_shape(:)
        integer, intent(in) :: n_from
        character(len=*), intent(in) :: caller

        if (all(arg_shape == n_from)) return
        write(stderr, *) 'Error in evaluator_base_mod, ', caller, ' function.'
        write(stderr, '(3x,a,1000(i0,:,1x))') 'Shape of the argument: ', arg_shape
        write(stderr, '(3x,5a,i0,a)') 'Transformation from the ', from_basis, ' to the ',        &
        &                             to_basis, ' basis needs ', n_from, ' states.'
        stop
    end subroutine check_transform_size


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: get_soc2
    !> @brief Squared magnitudes of the spin-orbit couplings between all states of a basis.
    !> @details
    !! The full space has the Hamiltonian H_sf + H_soc and the spin-free space H_sf (see
    !! basis_mod), so in a basis of the full space ('diabatic', 'mch') the spin-orbit couplings
    !! are get_Hamiltonian of that basis minus the Hamiltonian of the spin-free basis of the same
    !! kind, expanded to the spin components (expand of the layout). The spin-free part is
    !! subtracted rather than the blocks within the groups set to zero, which holds for complex
    !! couplings as well.
    !!
    !! In a basis of the spin-free space ('diabatic_sf', 'mch_sf'), element (a, b) is the sum of
    !! the squared couplings between all components of the spin-free states a and b
    !! (sum_components). For two manifolds of components, this is the root-sum-square coupling,
    !! squared. It depends neither on the signs of the states nor on the basis of the components.
    !!
    !! Stops for the 'adiabatic' basis, whose states include the spin-orbit couplings.
    !----------------------------------------------------------------------------------------------
    function get_soc2(self, basis) result(s2)
        class(potential_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        real(dp), allocatable :: s2(:, :)
        character(len=:), allocatable :: full

        if (.not. basis_known(basis)) then
            write(stderr, *) 'Error in evaluator_base_mod, get_soc2 function.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end if
        full = basis_full(basis)
        select case(full)
        case('diabatic')
            s2 = squared_couplings(self%model_states)
        case('mch')
            s2 = squared_couplings(self%states)
        case default
            write(stderr, *) 'Error in evaluator_base_mod, get_soc2 function.'
            write(stderr, *) '  No spin-orbit couplings in the ', basis, ' basis: its states ',  &
            &                'include them.'
            stop
        end select

    contains

        !> Squared couplings in the basis, for the states of layout.
        function squared_couplings(layout) result(s2)
            type(state_layout), intent(in) :: layout
            real(dp), allocatable :: s2(:, :)

            s2 = abs(self%get_Hamiltonian(full)                                                  &
            &        - layout%expand(self%get_Hamiltonian(basis_like(full, 'mch_sf'))))**2
            if (basis_is_spin_free(basis)) s2 = layout%sum_components(s2)
        end function squared_couplings

    end function get_soc2


    !> Id of the previous point, to which link connects this point. 0 without a link.
    pure function previous_point(self) result(id)
        class(potential_evaluator), intent(in) :: self
        integer :: id

        id = 0
        if (self%link%valid()) id = self%link%id_bra
    end function previous_point


    !> Basis in which set_link stores the overlaps with the states of the previous point: the
    !! basis of the request (the basis of the dynamics), or 'adiabatic' without a request.
    pure function link_basis(self) result(basis)
        class(potential_evaluator), intent(in) :: self
        character(len=:), allocatable :: basis

        if (allocated(self%request%basis)) then
            basis = self%request%basis
        else
            basis = 'adiabatic'
        end if
    end function link_basis


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: set_link
    !> @brief Set link, the overlaps between the states of link_basis of the previous point and of
    !! this point.
    !> @details
    !! The evaluator gives the overlaps between the two points in its anchor, a kind of states
    !! whose overlaps it knows: 'diabatic' for a model, whose diabatic states are the same at all
    !! points (s is not given), or 'mch' for states computed at each point, with the overlaps s
    !! of the spin-free states from their wave functions. The overlaps in the basis of the link
    !! are T_prev^T S T, where T is the transformation from the anchor to that basis at each
    !! point, and S is s (for every spin component in a basis of the full space) or the unit
    !! matrix. diab and adiab_vec of both points must be set, with their final phases.
    !!
    !! A warning is written if the smallest singular value of the overlaps is below
    !! overlap_warn_tol: the states of the two points then do not span the same space.
    !----------------------------------------------------------------------------------------------
    subroutine set_link(self, prev, anchor, s)
        use matrix_mod, only : unit_mat
        class(potential_evaluator), intent(inout) :: self
        class(potential_evaluator), intent(in) :: prev !< Evaluator of the previous point.
        character(len=*), intent(in) :: anchor !< 'diabatic' or 'mch'.
        real(dp), intent(in), optional :: s(:, :) !< Overlaps between the spin-free anchor states
                                                  !! of prev (rows) and of this point (columns).
        character(len=:), allocatable :: basis
        character(len=:), allocatable :: a !< Anchor basis in the space of the link basis.
        real(dp), allocatable :: s_link(:, :)
        real(dp) :: sigma_min

        basis = self%link_basis()
        a = basis_like(anchor, basis)
        if (a == basis) then
            if (present(s)) then
                s_link = anchor_overlap()
            else
                s_link = unit_mat(self%basis_size(basis))
            end if
        else
            s_link = self%get_transformation(a, basis)
            if (present(s)) s_link = matmul(anchor_overlap(), s_link)
            s_link = matmul(transpose(prev%get_transformation(a, basis)), s_link)
        end if
        call self%link%set(prev%id, basis, self%id, basis, s_link, sigma_min)
        if (sigma_min < overlap_warn_tol) then
            write(stderr, *) 'Warning in evaluator_base_mod, set_link subroutine.'
            write(stderr, '(3x,3a,f8.5,a)') 'Smallest singular value of the overlaps between '// &
            &   'the ', basis, ' states of two points: ', sigma_min, '.'
            write(stderr, *) '  A state of the dynamics changed character, for example at a ',   &
            &                'crossing with a state which is not in the dynamics.'
        end if

    contains

        !> s in the space of the link basis.
        function anchor_overlap() result(s_a)
            real(dp), allocatable :: s_a(:, :)

            if (basis_is_spin_free(basis)) then
                s_a = s
            else if (a == 'diabatic') then
                s_a = prev%model_states%expand(s, self%model_states)
            else
                s_a = prev%states%expand(s, self%states)
            end if
        end function anchor_overlap

    end subroutine set_link


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
    !! Adds the number of states of each of the five bases (n_state, see basis_mod), and
    !! spin_states and the spin manifolds of the states of the dynamics (see add_to_json of
    !! state_layout). The object 'model' holds spin_states and the spin manifolds of the states of
    !! the model, which are the states of the diabatic bases. Evaluators that are built from
    !! another potential override this to describe it as well.
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
        call json%add(nst, 'diabatic_sf', self%basis_size('diabatic_sf'))
        call json%add(nst, 'mch_sf', self%basis_size('mch_sf'))
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
    !! Only what depends on the point is saved: the quantities evaluated at it, its id, the
    !! relations between its bases (diab, adiab_vec) and its link to the previous point. The setup
    !! (model, spin manifolds, options) is built from the input when the backup is loaded. Each
    !! type adds its own object to obj (here 'base', with the numbers of states of the dynamics
    !! and of the model for checking), so a backup of another type of potential is detected on
    !! loading. Extensions override this and call the subroutine of their parent first.
    !----------------------------------------------------------------------------------------------
    subroutine evaluator_save_state(self, json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_add_array
        class(potential_evaluator), intent(in) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the state is added.
        type(json_value), pointer :: base, link, diab

        call json%create_object(base, 'base')
        call json%add(base, 'n_state', self%states%n_state)
        call json%add(base, 'n_sf', self%states%n_sf)
        call json%add(base, 'n_model_state', self%model_states%n_state)
        call json%add(base, 'n_model_sf', self%model_states%n_sf)
        call json%add(base, 'spin_states', self%states%spin_states)
        call json%add(base, 'id', self%id)
        call json%add(base, 'id_diab', self%id_diab)
        call json%create_object(diab, 'diab')
        call self%diab%save_state(json, diab)
        call json%add(base, diab)
        if (allocated(self%adiab_vec)) then
            call json_add_array(json, base, 'adiab_vec', self%adiab_vec)
        end if
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
        call json_get_required(json, base, 'id_diab', self%id_diab)
        call self%diab%load_state(json, json_get_child(json, base, 'diab'))
        if (allocated(self%adiab_vec)) deallocate(self%adiab_vec)
        if (json%valid_path(base, 'adiab_vec')) then
            call json_get_required(json, base, 'adiab_vec', self%adiab_vec,                      &
            &                      [self%states%n_state, self%states%n_state])
        end if
        call self%link%load_state(json, json_get_child(json, base, 'link'))
    end subroutine evaluator_load_state


end module evaluator_base_mod
