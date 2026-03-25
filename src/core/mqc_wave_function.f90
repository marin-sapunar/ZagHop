module mqc_wave_function_mod
    use global_defs
    use state_mod
    use rmat_mod
    implicit none
    private
    public :: mqc_wave_function

    type :: mqc_wave_function
        integer :: n_state_group = 0 !< Number of state groups in the wave function.
        integer :: n_state = 0 !< Total number of basis states in the wave function.
        integer, allocatable :: n_state_per_group(:) !< Number of states in each state group.
        type(state), allocatable :: qm_state(:) !< States as given by the quantum method. Dimension: n_state
        real(dp), allocatable :: en(:) !< Energies of the states. Dimension: n_state
        complex(dp), allocatable :: coeff(:) !< Coeffs of the el. states in the total wf.
                                             !! Dimension: n_state
        integer :: active_state !< Index of the active state.
        integer, allocatable :: phase(:) !< Phase factors for each state to keep wf continuity.
                                         !! Dimension: n_state
        real(dp), allocatable :: prob(:) !< Probability of hop to each state.
                                         !! Dimension: n_state
        real(dp), allocatable :: overlap(:, :) !< Overlap matrix between states.
                                               !! Dimension: (n_state, n_state)
        real(dp), allocatable :: adt(:, :) !< Adiabatic-to-diabatic transformation matrix.
                                           !! Dimension: (n_state, n_state)
        logical, allocatable :: need_gradient(:) !< Whether gradients are needed for each state.
                                        !! Dimension: n_state
        logical, allocatable :: need_nadv(:, :) !< Whether non-adiabatic vectors are needed.
                                        !! Dimension: (n_state, n_state)
        logical, allocatable :: need_soc(:, :) !< Whether spin-orbit couplings are needed.
                                        !! Dimension: (n_state, n_state)
    contains
        procedure :: initialize
        procedure :: index => index_group_state_to_full, index_full_to_group_state
    end type mqc_wave_function

contains

    subroutine initialize(self, n_state, multiplicity, active_state)
        class(mqc_wave_function), intent(inout) :: self
        integer, intent(in) :: n_state(:)
        integer, intent(in) :: multiplicity(:)
        integer, intent(in), optional :: active_state
        integer :: i, ms, igroup, istate, j
        integer :: cmult

        self%n_state = sum(n_state * multiplicity)
        allocate(self%qm_state(self%n_state))
        allocate(self%en(self%n_state), source=0.0_dp)
        allocate(self%coeff(self%n_state), source=(0.0_dp, 0.0_dp))
        allocate(self%prob(self%n_state), source=0.0_dp)
        allocate(self%phase(self%n_state), source=1)
        allocate(self%need_gradient(self%n_state), source=.false.)
        allocate(self%need_nadv(self%n_state, self%n_state), source=.false.)
        allocate(self%need_soc(self%n_state, self%n_state), source=.false.)
        self%n_state_group = sum(multiplicity, mask=(n_state > 0))
        allocate(self%n_state_per_group(self%n_state_group))
        igroup = 0
        istate = 0
        do i = 1, size(n_state)
            if (n_state(i) == 0) cycle
            if (mod(multiplicity(i), 2) == 1) then
                igroup = igroup + 1
                self%n_state_per_group(igroup) = n_state(i)
                do j = 1, n_state(i)
                    istate = istate + 1
                    call self%qm_state(istate)%initialize(igroup, j, istate, self%n_state, multiplicity(i)-1, 0)
                end do
                cmult = (multiplicity(i) - 1) / 2
            else
                cmult = multiplicity(i) / 2
            end if
            do ms = 1, cmult
                igroup = igroup + 1
                self%n_state_per_group(igroup) = n_state(i)
                do j = 1, n_state(i)
                    istate = istate + 1
                    call self%qm_state(istate)%initialize(igroup, j, istate, self%n_state, multiplicity(i)-1, ms)
                end do
                igroup = igroup + 1
                self%n_state_per_group(igroup) = n_state(i)
                do j = 1, n_state(i)
                    istate = istate + 1
                    call self%qm_state(istate)%initialize(igroup, j, istate, self%n_state, multiplicity(i)-1, -ms)
                end do
            end do
        end do

        if (present(active_state)) then
            if (active_state < 1 .or. active_state > self%n_state) then
                call errstop("mqc_wave_function_mod", "Active state index out of bounds.", 1)
            end if
            self%active_state = active_state
        end if
    end subroutine initialize


    function index_full_to_group_state(self, i_state) result(group_state)
        class(mqc_wave_function), intent(in) :: self
        integer, intent(in) :: i_state
        integer :: group_state(2)

        if (i_state < 1 .or. i_state > self%n_state) then
            call errstop("mqc_wave_function_mod", "State index out of bounds.", 1)
        end if

        group_state(1) = self%qm_state(i_state)%group
        group_state(2) = self%qm_state(i_state)%group_state
    end function index_full_to_group_state


    function index_group_state_to_full(self, group, group_state) result(i_state)
        class(mqc_wave_function), intent(in) :: self
        integer, intent(in) :: group
        integer, intent(in) :: group_state
        integer :: i_state

        if (group < 1 .or. group > self%n_state_group) then
            call errstop("mqc_wave_function_mod", "Group index out of bounds.", 1)
        end if
        if (group_state < 1 .or. group_state > self%n_state_per_group(group)) then
            call errstop("mqc_wave_function_mod", "Group state index out of bounds.", 1)
        end if

        do i_state = 1, self%n_state
            if (self%qm_state(i_state)%group /= group) cycle
            if (self%qm_state(i_state)%group_state == group_state) return
        end do

        call errstop("mqc_wave_function_mod", "State with given group and group_state not found.", 1)
    end function index_group_state_to_full



end module mqc_wave_function_mod