module evaluator_base_mod
    use global_defs
    implicit none

    private
    public :: potential_evaluator
    public :: spin_manifold
    public :: state_location
    public :: set_spin_manifolds
    public :: soc_group_tol


    real(dp), parameter :: soc_group_tol = 1.0e-10_dp !< Largest spin-orbit coupling within a group
                                                      !! of states without a warning.


    !> Spin-free states of one multiplicity.
    !> @details
    !! With spin_states = 'components', each spin component is a separate state, and the states of
    !! the manifold form mult consecutive groups of n states, one group for each component. With
    !! spin_states = 'effective', each spin-free state is a single state, in a single group.
    type :: spin_manifold
        integer :: mult = 1 !< Multiplicity, 2S+1.
        integer :: n = 0 !< Number of spin-free states.
        integer :: n_comp = 1 !< Number of groups: mult for components, 1 for effective states.
        integer :: sf0 = 0 !< First spin-free state of the manifold.
        integer :: i0 = 0 !< First state of the first group of the manifold.
    end type spin_manifold


    !> Base class of all potentials.
    !> @details
    !! Specific evaluators only implement set_geometry and the evaluation of the potential at the
    !! current geometry, and fill the variables of the base classes. The base classes do the rest.
    type, abstract :: potential_evaluator
        integer :: n_state = 0
        integer :: n_sf = 0 !< Number of spin-free states.
        type(spin_manifold), allocatable :: manifold(:)
        character(len=:), allocatable :: spin_states !< 'components' or 'effective'.
        logical :: mean_couplings = .false. !< Couplings between manifolds of effective states
                                            !! are given as the mean over the component pairs.
        character(len=:), allocatable :: classical_coordinates !< Coordinates expected by
                                                               !! update_geometry.
        integer :: n_source_states = 0 !< Number of states given by the source of the potential,
                                       !! checked by set_spin_manifolds. 0 if not fixed.
        character(len=12) :: component_order = 'by_component' !< Order of the spin components in
        !! a full Hamiltonian given by the source of the potential:
        !! - by_component - All states of the first component, then of the second, ... (ZagHop).
        !! - by_state - All components of the first state, then of the second, ... (Quantics).
        integer, allocatable :: state_map(:) !< Index in the source of each state (n_state).
    contains
        procedure(eval), deferred :: eval
        procedure(set_geometry), deferred :: set_geometry
        procedure, non_overridable :: update_geometry
        procedure :: save_previous
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
        procedure :: set_spin_manifolds
        procedure :: spin_states_keyword
        procedure :: basis_size
        procedure :: locate_state
        procedure :: expand
        procedure :: sum_components
        procedure :: to_rss_couplings
        procedure :: max_soc_in_groups
        procedure :: get_reference_geometry
        generic :: get_energy => get_energy_single, get_energy_all
    end type potential_evaluator


    !> Position of a state in the spin manifolds of a potential_evaluator.
    type :: state_location
        integer :: manifold = 0 !< Manifold containing the state.
        integer :: mult = 0 !< Multiplicity of the state.
        integer :: component = 0 !< Spin component (position of the group within the manifold).
        integer :: index = 0 !< Index of the state within its group.
        integer :: sf = 0 !< Spin-free state.
        integer :: group = 0 !< Group containing the state, numbered in state order.
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
        !> Set the geometry at which the potential is evaluated next. Coordinates gives the
        !! coordinate system of the geometry. Anything that get_derivative_target_coord needs to
        !! convert derivatives to these coordinates must be stored here.
        subroutine set_geometry(self, geometry, coordinates)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator), intent(inout) :: self
            real(dp), intent(in) :: geometry(:, :)
            character(len=*), intent(in) :: coordinates
        end subroutine set_geometry
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
        !! 'spin_free' basis, the squares are summed over the components of each pair of
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
    ! SUBROUTINE: update_geometry
    !> @brief Move the potential to a new geometry.
    !> @details
    !! Saves what the evaluator keeps from the previous geometry (save_previous) and passes the
    !! geometry to the evaluator (set_geometry). Gradients and nonadiabatic couplings are returned
    !! in the coordinates of this geometry.
    !----------------------------------------------------------------------------------------------
    subroutine update_geometry(self, geometry, coordinates)
        class(potential_evaluator), intent(inout) :: self
        real(dp), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates

        call self%save_previous()
        call self%set_geometry(geometry, coordinates)
    end subroutine update_geometry


    !> Save the quantities of the current geometry that are needed at the next geometry. Nothing
    !! is saved in the base class.
    subroutine save_previous(self)
        class(potential_evaluator), intent(inout) :: self
    end subroutine save_previous


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: get_derivative_target_coord
    !> @brief Convert a derivative with respect to the coordinates of the potential to the
    !! coordinates of the last geometry passed to update_geometry.
    !> @details
    !! The result is flat. In the base class the two coordinate systems are the same. An
    !! evaluator that accepts geometries in other coordinates overrides this function. Everything
    !! it needs must be stored in set_geometry, not in update_geometry.
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
    ! SUBROUTINE: set_spin_manifolds
    !> @brief Set the spin manifolds from the number of states of each multiplicity.
    !> @details
    !! Each block i of nstate(i) states with multiplicity multiplicity(i) gives one manifold.
    !! Blocks without states are skipped. Each multiplicity can be given only once, so couplings
    !! between manifolds are always spin-orbit couplings.
    !!
    !! spin_states is one of:
    !! - components - Each spin component is a separate state.
    !! - effective_rss - Each spin-free state is a single state. Couplings between manifolds are
    !!   the root-sum-square over the component pairs.
    !! - effective_mean - As effective_rss, but the squared couplings between manifolds are the
    !!   mean over the component pairs.
    !!
    !! If n_source_states is set, the number of states must match it. state_map is built from
    !! component_order.
    !----------------------------------------------------------------------------------------------
    subroutine set_spin_manifolds(self, nstate, multiplicity, spin_states)
        class(potential_evaluator), intent(inout) :: self
        integer, intent(in) :: nstate(:) !< Number of spin-free states of each multiplicity.
        integer, intent(in) :: multiplicity(:) !< Multiplicity of each block of states.
        character(len=*), intent(in) :: spin_states
        integer :: i, m, sf, ci, c, k

        select case(spin_states)
        case('components')
            self%spin_states = 'components'
            self%mean_couplings = .false.
        case('effective_rss')
            self%spin_states = 'effective'
            self%mean_couplings = .false.
        case('effective_mean')
            self%spin_states = 'effective'
            self%mean_couplings = .true.
        case default
            write(stderr, *) 'Error in evaluator_base_mod, set_spin_manifolds subroutine.'
            write(stderr, *) '  Unrecognized spin_states: ', spin_states
            write(stderr, *) '  Expected components, effective_rss or effective_mean.'
            stop
        end select

        do i = 1, size(nstate)
            if (nstate(i) == 0) cycle
            if (multiplicity(i) < 1) then
                write(stderr, *) 'Error in evaluator_base_mod, set_spin_manifolds subroutine.'
                write(stderr, '(3x,a,i0,a)') 'Invalid multiplicity: ', multiplicity(i), '.'
                stop
            end if
            if (any(multiplicity(:i-1) == multiplicity(i) .and. nstate(:i-1) > 0)) then
                write(stderr, *) 'Error in evaluator_base_mod, set_spin_manifolds subroutine.'
                write(stderr, '(3x,a,i0,a)') 'Multiplicity ', multiplicity(i), ' is given twice.'
                write(stderr, *) '  Give all states of the same multiplicity in one block.'
                stop
            end if
        end do

        if (allocated(self%manifold)) deallocate(self%manifold)
        allocate(self%manifold(count(nstate > 0)))
        m = 0
        sf = 1
        ci = 1
        do i = 1, size(nstate)
            if (nstate(i) == 0) cycle
            m = m + 1
            self%manifold(m)%mult = multiplicity(i)
            self%manifold(m)%n = nstate(i)
            self%manifold(m)%n_comp = 1
            if (self%spin_states == 'components') self%manifold(m)%n_comp = multiplicity(i)
            self%manifold(m)%sf0 = sf
            self%manifold(m)%i0 = ci
            sf = sf + nstate(i)
            ci = ci + nstate(i) * self%manifold(m)%n_comp
        end do
        self%n_sf = sf - 1
        self%n_state = ci - 1

        if (self%n_source_states > 0 .and. self%n_state /= self%n_source_states) then
            write(stderr, *) 'Error in evaluator_base_mod, set_spin_manifolds subroutine.'
            write(stderr, '(3x,a,i0,a)') 'Input nstate, multiplicity and spin_states give ',     &
            &                            self%n_state, ' states.'
            write(stderr, '(3x,a,i0,a)') 'The potential has ', self%n_source_states, ' states.'
            stop
        end if

        ! Index in the source of each state. With one component, both orders are the same.
        if (allocated(self%state_map)) deallocate(self%state_map)
        allocate(self%state_map(self%n_state))
        do m = 1, size(self%manifold)
            associate(mf => self%manifold(m))
                do c = 1, mf%n_comp
                    do k = 1, mf%n
                        ci = mf%i0 + (c - 1) * mf%n + k - 1
                        select case(self%component_order)
                        case('by_component')
                            self%state_map(ci) = ci
                        case('by_state')
                            self%state_map(ci) = mf%i0 + (k - 1) * mf%n_comp + c - 1
                        case default
                            write(stderr, *) 'Error in evaluator_base_mod, set_spin_manifolds ',  &
                            &                'subroutine.'
                            write(stderr, *) '  Unknown component_order: ', trim(self%component_order)
                            stop
                        end select
                    end do
                end do
            end associate
        end do
    end subroutine set_spin_manifolds


    !> Value of the spin_states keyword that gives the spin states of the potential.
    pure function spin_states_keyword(self) result(keyword)
        class(potential_evaluator), intent(in) :: self
        character(len=:), allocatable :: keyword

        if (self%spin_states == 'components') then
            keyword = 'components'
        else if (self%mean_couplings) then
            keyword = 'effective_mean'
        else
            keyword = 'effective_rss'
        end if
    end function spin_states_keyword


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: basis_size
    !> @brief Number of states of a basis.
    !> @details
    !! The 'spin_free' basis has one state for each spin-free state. The other bases have one
    !! state for each state of the potential.
    !----------------------------------------------------------------------------------------------
    pure function basis_size(self, basis) result(n)
        class(potential_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer :: n

        if (basis == 'spin_free') then
            n = self%n_sf
        else
            n = self%n_state
        end if
    end function basis_size


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: locate_state
    !> @brief Manifold, multiplicity, spin component, index within the group, spin-free state and
    !! group of a state.
    !----------------------------------------------------------------------------------------------
    pure function locate_state(self, state) result(loc)
        class(potential_evaluator), intent(in) :: self
        integer, intent(in) :: state
        type(state_location) :: loc
        integer :: m, k

        if (state < 1 .or. state > self%n_state) then
            error stop 'Error in evaluator_base_mod, locate_state function. State out of range.'
        end if
        loc%group = 0
        do m = 1, size(self%manifold)
            if (state < self%manifold(m)%i0 + self%manifold(m)%n * self%manifold(m)%n_comp) exit
            loc%group = loc%group + self%manifold(m)%n_comp
        end do
        k = state - self%manifold(m)%i0
        loc%manifold = m
        loc%mult = self%manifold(m)%mult
        loc%component = k / self%manifold(m)%n + 1
        loc%index = mod(k, self%manifold(m)%n) + 1
        loc%sf = self%manifold(m)%sf0 + loc%index - 1
        loc%group = loc%group + loc%component
    end function locate_state


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: expand
    !> @brief Copy each manifold block of a matrix over the spin-free states onto every group of the
    !! manifold.
    !> @details
    !! Elements between manifolds of the spin-free matrix are not copied, and the elements between
    !! groups of the result are zero.
    !----------------------------------------------------------------------------------------------
    pure function expand(self, sf_mat) result(mat)
        class(potential_evaluator), intent(in) :: self
        real(dp), intent(in) :: sf_mat(:, :) !< Matrix over the spin-free states.
        real(dp) :: mat(self%n_state, self%n_state)
        integer :: m, c, i0, s0, n

        mat = 0.0_dp
        do m = 1, size(self%manifold)
            n = self%manifold(m)%n
            s0 = self%manifold(m)%sf0
            do c = 1, self%manifold(m)%n_comp
                i0 = self%manifold(m)%i0 + (c - 1) * n
                mat(i0:i0+n-1, i0:i0+n-1) = sf_mat(s0:s0+n-1, s0:s0+n-1)
            end do
        end do
    end function expand


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: sum_components
    !> @brief Sum a matrix over the states to a matrix over the spin-free states.
    !> @details
    !! Element (a, b) of the result is the sum of the elements (i, j) for all components i of
    !! spin-free state a and all components j of spin-free state b.
    !----------------------------------------------------------------------------------------------
    pure function sum_components(self, mat) result(sf_mat)
        class(potential_evaluator), intent(in) :: self
        real(dp), intent(in) :: mat(:, :) !< Matrix over the states.
        real(dp) :: sf_mat(self%n_sf, self%n_sf)
        type(state_location) :: loc
        integer :: sf(self%n_state)
        integer :: i, j

        do i = 1, self%n_state
            loc = self%locate_state(i)
            sf(i) = loc%sf
        end do
        sf_mat = 0.0_dp
        do j = 1, self%n_state
            do i = 1, self%n_state
                sf_mat(sf(i), sf(j)) = sf_mat(sf(i), sf(j)) + mat(i, j)
            end do
        end do
    end function sum_components


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: to_rss_couplings
    !> @brief Convert couplings between manifolds from the mean to the root-sum-square convention.
    !> @details
    !! The squared couplings between states of manifolds a and b are multiplied by
    !! mult_a * mult_b, so the couplings by its square root.
    !----------------------------------------------------------------------------------------------
    pure subroutine to_rss_couplings(self, w)
        class(potential_evaluator), intent(in) :: self
        real(dp), intent(inout) :: w(:, :) !< Matrix over the states.
        integer :: a, b, ia, ib, na, nb

        do a = 1, size(self%manifold)
            ia = self%manifold(a)%i0
            na = self%manifold(a)%n * self%manifold(a)%n_comp
            do b = 1, size(self%manifold)
                if (a == b) cycle
                ib = self%manifold(b)%i0
                nb = self%manifold(b)%n * self%manifold(b)%n_comp
                w(ia:ia+na-1, ib:ib+nb-1) = w(ia:ia+na-1, ib:ib+nb-1)                            &
                &   * sqrt(real(self%manifold(a)%mult * self%manifold(b)%mult, dp))
            end do
        end do
    end subroutine to_rss_couplings


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: max_soc_in_groups
    !> @brief Largest magnitude of the spin-orbit couplings between states of the same group.
    !> @details
    !! For a real Hamiltonian, these couplings are zero, and a nonzero value points to an input
    !! error, such as components which are not ordered as assumed.
    !----------------------------------------------------------------------------------------------
    pure function max_soc_in_groups(self, soc) result(smax)
        class(potential_evaluator), intent(in) :: self
        real(dp), intent(in) :: soc(:, :) !< Spin-orbit couplings between the states.
        real(dp) :: smax
        integer :: m, c, i0, n

        smax = 0.0_dp
        do m = 1, size(self%manifold)
            n = self%manifold(m)%n
            do c = 1, self%manifold(m)%n_comp
                i0 = self%manifold(m)%i0 + (c - 1) * n
                smax = max(smax, maxval(abs(soc(i0:i0+n-1, i0:i0+n-1))))
            end do
        end do
    end function max_soc_in_groups


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
