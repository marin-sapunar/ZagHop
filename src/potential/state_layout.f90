!--------------------------------------------------------------------------------------------------
! MODULE: state_layout_mod
!> @brief Layout of a set of electronic states in spin manifolds.
!> @details
!! A state_layout gives the number of spin-free states of each multiplicity (the spin manifolds),
!! and the order of the states built from them. With spin_states = 'components', each spin
!! component is a separate state. The states of a manifold then form mult consecutive groups of n
!! states, one group for each component. With spin_states = 'effective_rss' or 'effective_mean',
!! each spin-free state is a single state, in a single group.
!!
!! The layout does not depend on the basis of the states. A potential has one layout for the
!! states of its model and one for the states of the dynamics, which can be a subset of the
!! first (subset, subset_index).
!--------------------------------------------------------------------------------------------------
module state_layout_mod
    use global_defs
    implicit none

    private
    public :: state_layout
    public :: spin_manifold
    public :: state_location


    !> Spin-free states of one multiplicity.
    type :: spin_manifold
        integer :: mult = 1 !< Multiplicity, 2S+1.
        integer :: n = 0 !< Number of spin-free states.
        integer :: n_comp = 1 !< Number of groups: mult for components, 1 for effective states.
        integer :: sf0 = 0 !< First spin-free state of the manifold.
        integer :: i0 = 0 !< First state of the first group of the manifold.
    end type spin_manifold


    !> Position of a state in the spin manifolds of a state_layout.
    type :: state_location
        integer :: manifold = 0 !< Manifold containing the state.
        integer :: mult = 0 !< Multiplicity of the state.
        integer :: component = 0 !< Spin component (position of the group within the manifold).
        integer :: index = 0 !< Index of the state within its group.
        integer :: sf = 0 !< Spin-free state.
        integer :: group = 0 !< Group containing the state, numbered in state order.
    end type state_location


    type :: state_layout
        character(len=:), allocatable :: spin_states !< 'components' or 'effective_*'.
        type(spin_manifold), allocatable :: manifold(:)
        integer :: n_state = 0 !< Number of states.
        integer :: n_sf = 0 !< Number of spin-free states.
    contains
        procedure :: set
        procedure :: is_set
        procedure :: subset
        procedure :: subset_index
        procedure :: locate_state
        procedure :: expand
        procedure :: sum_components
        procedure :: to_rss_couplings
        procedure :: max_soc_in_groups
        procedure :: add_to_json
    end type state_layout


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: set
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
    !----------------------------------------------------------------------------------------------
    subroutine set(self, nstate, multiplicity, spin_states)
        class(state_layout), intent(inout) :: self
        integer, intent(in) :: nstate(:) !< Number of spin-free states of each multiplicity.
        integer, intent(in) :: multiplicity(:) !< Multiplicity of each block of states.
        character(len=*), intent(in) :: spin_states
        integer :: i, m, sf, ci

        select case(spin_states)
        case('components', 'effective_rss', 'effective_mean')
            self%spin_states = spin_states
        case default
            write(stderr, *) 'Error in state_layout_mod, set subroutine.'
            write(stderr, *) '  Unrecognized spin_states: ', spin_states
            write(stderr, *) '  Expected components, effective_rss or effective_mean.'
            stop
        end select

        if (size(nstate) /= size(multiplicity)) then
            write(stderr, *) 'Error in state_layout_mod, set subroutine.'
            write(stderr, '(3x,a,i0,a,i0,a)') 'Number of states given for ', size(nstate),        &
            &   ' blocks, multiplicity for ', size(multiplicity), '.'
            stop
        end if
        do i = 1, size(nstate)
            if (nstate(i) == 0) cycle
            if (nstate(i) < 0) then
                write(stderr, *) 'Error in state_layout_mod, set subroutine.'
                write(stderr, '(3x,a,i0,a)') 'Invalid number of states: ', nstate(i), '.'
                stop
            end if
            if (multiplicity(i) < 1) then
                write(stderr, *) 'Error in state_layout_mod, set subroutine.'
                write(stderr, '(3x,a,i0,a)') 'Invalid multiplicity: ', multiplicity(i), '.'
                stop
            end if
            if (any(multiplicity(:i-1) == multiplicity(i) .and. nstate(:i-1) > 0)) then
                write(stderr, *) 'Error in state_layout_mod, set subroutine.'
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
    end subroutine set


    !> The spin manifolds are set.
    pure logical function is_set(self)
        class(state_layout), intent(in) :: self

        is_set = allocated(self%manifold)
    end function is_set


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: subset
    !> @brief Layout with only the first nstate(i) spin-free states of each multiplicity.
    !> @details
    !! The manifolds of the subset are in the order of this layout, also when the multiplicities
    !! are given in another order. A multiplicity which is not given has no states in the subset.
    !! Stops if a multiplicity with states is not in this layout, if more states are requested
    !! than this layout has, or if no state is left.
    !----------------------------------------------------------------------------------------------
    function subset(self, nstate, multiplicity) result(sub)
        class(state_layout), intent(in) :: self
        integer, intent(in) :: nstate(:) !< Number of spin-free states of each multiplicity.
        integer, intent(in) :: multiplicity(:) !< Multiplicity of each block of states.
        type(state_layout) :: sub
        integer :: n_sub(size(self%manifold))
        integer :: i, m

        if (size(nstate) /= size(multiplicity)) then
            call subset_error()
            write(stderr, '(3x,a,i0,a,i0,a)') 'Number of states given for ', size(nstate),        &
            &   ' blocks, multiplicity for ', size(multiplicity), '.'
            stop
        end if
        n_sub = 0
        do i = 1, size(nstate)
            if (nstate(i) == 0) cycle
            if (any(multiplicity(:i-1) == multiplicity(i) .and. nstate(:i-1) > 0)) then
                call subset_error()
                write(stderr, '(3x,a,i0,a)') 'Multiplicity ', multiplicity(i), ' is given twice.'
                stop
            end if
            m = findloc(self%manifold%mult, multiplicity(i), dim=1)
            if (m == 0) then
                call subset_error()
                write(stderr, '(3x,a,i0,a)') 'There are no states of multiplicity ',              &
                &                            multiplicity(i), '.'
                stop
            end if
            if (nstate(i) < 0 .or. nstate(i) > self%manifold(m)%n) then
                call subset_error()
                write(stderr, '(3x,i0,a,i0,a,i0,a)') nstate(i), ' states of multiplicity ',       &
                &   multiplicity(i), ' requested, ', self%manifold(m)%n, ' available.'
                stop
            end if
            n_sub(m) = nstate(i)
        end do
        if (all(n_sub == 0)) then
            call subset_error()
            write(stderr, *) '  No states requested.'
            stop
        end if
        call sub%set(n_sub, self%manifold%mult, self%spin_states)

    contains

        subroutine subset_error()
            write(stderr, *) 'Error in state_layout_mod, subset function.'
        end subroutine subset_error

    end function subset


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: subset_index
    !> @brief Position in this layout of each state of a subset of it.
    !> @details
    !! sub has, for each multiplicity, the first spin-free states of this layout (see subset).
    !! Component c of spin-free state k of a manifold is the same component of the same spin-free
    !! state in both layouts.
    !----------------------------------------------------------------------------------------------
    pure subroutine subset_index(self, sub, idx, idx_sf, rest)
        class(state_layout), intent(in) :: self
        type(state_layout), intent(in) :: sub !< Subset of this layout.
        integer, allocatable, intent(out) :: idx(:) !< State of this layout for each state of sub.
        integer, allocatable, intent(out) :: idx_sf(:) !< The same for the spin-free states.
        integer, allocatable, intent(out) :: rest(:) !< States of this layout which are not in sub.
        logical :: in_sub(self%n_state)
        integer :: ms, m, c, k, i

        allocate(idx(sub%n_state))
        allocate(idx_sf(sub%n_sf))
        in_sub = .false.
        do ms = 1, size(sub%manifold)
            m = findloc(self%manifold%mult, sub%manifold(ms)%mult, dim=1)
            if (m == 0) error stop 'Error in state_layout_mod, subset_index. Not a subset.'
            associate(mf => self%manifold(m), sf => sub%manifold(ms))
                if (sf%n > mf%n .or. sf%n_comp /= mf%n_comp) then
                    error stop 'Error in state_layout_mod, subset_index. Not a subset.'
                end if
                do k = 1, sf%n
                    idx_sf(sf%sf0 + k - 1) = mf%sf0 + k - 1
                    do c = 1, sf%n_comp
                        i = mf%i0 + (c - 1) * mf%n + k - 1
                        idx(sf%i0 + (c - 1) * sf%n + k - 1) = i
                        in_sub(i) = .true.
                    end do
                end do
            end associate
        end do
        rest = pack([(i, i = 1, self%n_state)], .not. in_sub)
    end subroutine subset_index


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: locate_state
    !> @brief Manifold, multiplicity, spin component, index within the group, spin-free state and
    !! group of a state.
    !----------------------------------------------------------------------------------------------
    pure function locate_state(self, state) result(loc)
        class(state_layout), intent(in) :: self
        integer, intent(in) :: state
        type(state_location) :: loc
        integer :: m, k

        if (state < 1 .or. state > self%n_state) then
            error stop 'Error in state_layout_mod, locate_state function. State out of range.'
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
        class(state_layout), intent(in) :: self
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
        class(state_layout), intent(in) :: self
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
        class(state_layout), intent(in) :: self
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
        class(state_layout), intent(in) :: self
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
    ! SUBROUTINE: add_to_json
    !> @brief Add spin_states and the spin manifolds to a JSON object.
    !> @details
    !! For each manifold: multiplicity, n (spin-free states), n_comp (groups: the spin components,
    !! or 1 for effective states), first_spin_free and first_state (1-based). State
    !!     first_state + (c - 1) * n + k - 1
    !! is component c of spin-free state first_spin_free + k - 1.
    !----------------------------------------------------------------------------------------------
    subroutine add_to_json(self, json, obj)
        use json_module, only : json_core, json_value
        class(state_layout), intent(in) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the layout is added.
        type(json_value), pointer :: arr, mf
        integer :: m

        call json%add(obj, 'spin_states', self%spin_states)
        call json%create_array(arr, 'manifolds')
        do m = 1, size(self%manifold)
            call json%create_object(mf, '')
            call json%add(mf, 'multiplicity', self%manifold(m)%mult)
            call json%add(mf, 'n', self%manifold(m)%n)
            call json%add(mf, 'n_comp', self%manifold(m)%n_comp)
            call json%add(mf, 'first_spin_free', self%manifold(m)%sf0)
            call json%add(mf, 'first_state', self%manifold(m)%i0)
            call json%add(arr, mf)
        end do
        call json%add(obj, arr)
    end subroutine add_to_json


end module state_layout_mod
