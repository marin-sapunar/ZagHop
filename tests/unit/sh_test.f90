module sh_test_mod
    use global_defs, only : dp
    use evaluator_diabatic_mod, only : diabatic_evaluator
    implicit none

    !> Evaluator with energies and spin-orbit couplings set directly by the test. The spin-free
    !! matrix is diagonal, so the energies of the spin-free states must be given in ascending order
    !! within each manifold.
    type, extends(diabatic_evaluator) :: scripted_evaluator
    contains
        procedure :: setup => scripted_setup
        procedure :: set_energies => scripted_set_energies
        procedure :: set_soc => scripted_set_soc
        procedure :: update_geometry => scripted_update_geometry
        procedure :: eval_diab => scripted_eval_diab
    end type scripted_evaluator

contains

    subroutine scripted_setup(self, nstate, multiplicity, spin_states)
        class(scripted_evaluator), intent(inout) :: self
        integer, intent(in) :: nstate(:)
        integer, intent(in) :: multiplicity(:)
        character(len=*), intent(in), optional :: spin_states

        if (present(spin_states)) then
            call self%set_spin_manifolds(nstate, multiplicity, spin_states)
        else
            call self%set_spin_manifolds(nstate, multiplicity, 'components')
        end if
        self%ndof = 1
        call self%allocate_hamiltonian()
    end subroutine scripted_setup

    subroutine scripted_set_energies(self, energy)
        class(scripted_evaluator), intent(inout) :: self
        real(dp), intent(in) :: energy(:)
        integer :: i

        self%sf_w = 0.0_dp
        do i = 1, self%n_sf
            self%sf_w(i, i) = energy(i)
        end do
        call self%eval()
    end subroutine scripted_set_energies

    !> Spin-orbit couplings between the states, used from the next call of set_energies.
    subroutine scripted_set_soc(self, soc)
        class(scripted_evaluator), intent(inout) :: self
        real(dp), intent(in) :: soc(:, :)

        self%soc_w = soc
    end subroutine scripted_set_soc

    subroutine scripted_update_geometry(self, geometry, coordinates)
        class(scripted_evaluator), intent(inout) :: self
        real(dp), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates
    end subroutine scripted_update_geometry

    subroutine scripted_eval_diab(self)
        class(scripted_evaluator), intent(inout) :: self
    end subroutine scripted_eval_diab

end module sh_test_mod


program sh_test
    implicit none
    integer :: ret


    ret = trajectory_buffer_test()
    if (ret /= 0) then
        write(*,*) "Error in trajectory_buffer_test."
        stop ret
    end if

    ret = lz_formula_test()
    if (ret /= 0) then
        write(*,*) "Error in lz_formula_test."
        stop ret
    end if

    ret = check_gap_test()
    if (ret /= 0) then
        write(*,*) "Error in check_gap_test."
        stop ret
    end if

    ret = lz_extra_minimum_test()
    if (ret /= 0) then
        write(*,*) "Error in lz_extra_minimum_test."
        stop ret
    end if

    ret = lz_guard_test()
    if (ret /= 0) then
        write(*,*) "Error in lz_guard_test."
        stop ret
    end if

    ret = lz_isc_test()
    if (ret /= 0) then
        write(*,*) "Error in lz_isc_test."
        stop ret
    end if

    contains


    !> Times of all stored points in time order.
    function stored_times() result(times)
        use global_defs, only : dp
        use system_type_mod, only : system_type, trajectory_n_point, trajectory_point
        real(dp), allocatable :: times(:)
        type(system_type), pointer :: pt
        integer :: i

        allocate(times(trajectory_n_point()))
        do i = 1, size(times)
            pt => trajectory_point(i)
            times(i) = pt%time
        end do
    end function stored_times


    integer function trajectory_buffer_test() result(ret)
        use global_defs, only : dp
        use system_type_mod, only : system_type, tr1, tr2, allocate_trajectory_data,             &
        &                           trajectory_next, trajectory_goto, trajectory_n_point,        &
        &                           trajectory_point
        real(dp), parameter :: dt_0 = 2.0_dp
        type(system_type), pointer :: pt
        integer :: i

        ret = 1
        call allocate_trajectory_data(6)
        tr1%step = 0
        tr1%time = 0.0_dp
        allocate(tr1%geom(1, 1))
        ! Mark each point by its geometry to check that the data moves with the point.
        do i = 1, 3
            call trajectory_next(dt_0)
            tr1%geom = real(i, dp)
        end do
        if (any(stored_times() /= [0.0_dp, 2.0_dp, 4.0_dp, 6.0_dp])) return
        if (tr1%step /= 3 .or. tr1%dt /= dt_0 .or. .not. tr1%on_grid) return
        if (tr2%time /= 4.0_dp) return

        ! Insert a point between the points at 2 and 4, then one between the points at 2 and 3.
        ret = 2
        call trajectory_goto(2.0_dp, keep_later=.true.)
        if (tr1%time /= 2.0_dp .or. tr1%geom(1, 1) /= 1.0_dp .or. tr2%time /= 0.0_dp) return
        call trajectory_next(dt_0, 1.0_dp)
        tr1%geom = 1.5_dp
        if (any(stored_times() /= [0.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 6.0_dp])) return
        if (tr1%on_grid .or. tr1%step /= 1 .or. tr1%dt /= 1.0_dp .or. tr1%written) return
        call trajectory_goto(2.0_dp, keep_later=.true.)
        call trajectory_next(dt_0, 0.5_dp)
        tr1%geom = 1.25_dp
        if (any(stored_times() /= [0.0_dp, 2.0_dp, 2.5_dp, 3.0_dp, 4.0_dp, 6.0_dp])) return
        if (tr2%time /= 2.0_dp) return

        ! Later points survive the insertions.
        ret = 3
        pt => trajectory_point(4)
        if (pt%geom(1, 1) /= 1.5_dp) return
        pt => trajectory_point(5)
        if (pt%geom(1, 1) /= 2.0_dp) return
        pt => trajectory_point(6)
        if (pt%geom(1, 1) /= 3.0_dp .or. .not. pt%on_grid .or. pt%step /= 3) return

        ! Look points up by time.
        ret = 4
        call trajectory_goto(4.0_dp, keep_later=.true.)
        if (tr1%geom(1, 1) /= 2.0_dp .or. tr2%time /= 3.0_dp) return
        call trajectory_goto(6.0_dp, keep_later=.false.)
        if (tr1%geom(1, 1) /= 3.0_dp .or. trajectory_n_point() /= 6) return

        ! Going back without keeping later points discards them. The next step from a point off the
        ! grid ends on the grid.
        ret = 5
        call trajectory_goto(2.5_dp, keep_later=.false.)
        if (trajectory_n_point() /= 3 .or. tr1%geom(1, 1) /= 1.25_dp) return
        call trajectory_next(dt_0)
        if (tr1%time /= 4.0_dp .or. tr1%dt /= 1.5_dp .or. tr1%step /= 2 .or. .not. tr1%on_grid) &
        &   return
        if (any(stored_times() /= [0.0_dp, 2.0_dp, 2.5_dp, 4.0_dp])) return

        ! When all slots are used, the oldest written point is dropped.
        ret = 6
        do i = 1, trajectory_n_point()
            pt => trajectory_point(i)
            pt%written = .true.
        end do
        do i = 1, 3
            call trajectory_next(dt_0)
            tr1%written = .true.
        end do
        if (any(stored_times() /= [2.0_dp, 2.5_dp, 4.0_dp, 6.0_dp, 8.0_dp, 10.0_dp])) return
        if (tr1%step /= 5) return

        ! More slots keep the stored points.
        ret = 7
        call allocate_trajectory_data(8)
        if (any(stored_times() /= [2.0_dp, 2.5_dp, 4.0_dp, 6.0_dp, 8.0_dp, 10.0_dp])) return
        if (tr1%time /= 10.0_dp .or. tr2%time /= 8.0_dp) return
        call trajectory_next(dt_0)
        call trajectory_next(dt_0)
        if (trajectory_n_point() /= 8 .or. tr1%time /= 14.0_dp) return
        ret = 0
    end function trajectory_buffer_test


    integer function lz_formula_test() result(ret)
        use global_defs, only : dp
        use sh_lz_mod, only : lz_prob, sec_deriv_3p, sec_deriv_3p_err
        real(dp) :: x(3), y(3), sd(3)
        real(dp), parameter :: tiny = 1.0e-12_dp

        ! exp(-pi/2 * sqrt(gap**3 / sd))
        ret = 1
        if (abs(lz_prob(0.01_dp, 1.0e-6_dp) - 0.20787957635076193_dp) > tiny) return
        if (abs(lz_prob(0.02_dp, 2.0e-6_dp) - 0.04321391826377226_dp) > tiny) return

        ! Exact for a quadratic on an uneven grid.
        ret = 2
        x = [0.3_dp, 1.0_dp, 2.7_dp]
        y = 1.7_dp * x**2 - 0.4_dp * x + 0.2_dp
        if (abs(sec_deriv_3p(x(1), x(2), x(3), y(1), y(2), y(3)) - 3.4_dp) > tiny) return

        ! Errors aligned for the smallest second derivative (smallest probability) first.
        ret = 3
        sd = sec_deriv_3p_err(0.0_dp, 1.0_dp, 2.5_dp, 0.05_dp, 0.01_dp, 0.04_dp, 1.0e-3_dp)
        if (abs(sd(2) - 0.048_dp) > tiny) return
        if (.not. (sd(1) <= sd(2) .and. sd(2) <= sd(3))) return
        ! Errors larger than the change of the gap.
        sd = sec_deriv_3p_err(0.0_dp, 1.0_dp, 2.0_dp, 0.0105_dp, 0.01_dp, 0.0105_dp, 1.0e-3_dp)
        if (sd(1) /= 0.0_dp .or. sd(2) > sd(3)) return
        ret = 0
    end function lz_formula_test


    integer function check_gap_test() result(ret)
        use global_defs, only : dp
        use sh_lz_mod, only : check_gap
        integer, parameter :: group(4) = [1, 1, 1, 2]
        real(dp) :: qe0(4), qe1(4), qe2(4)
        logical :: skip(4)

        ! Active state 3. Gap minima at the middle point with states 1 (not a neighbour), 2 and
        ! 4 (another group).
        qe0 = [-0.30_dp, -0.20_dp, 0.0_dp, 0.20_dp]
        qe1 = [-0.01_dp, -0.05_dp, 0.0_dp, 0.05_dp]
        qe2 = [-0.30_dp, -0.10_dp, 0.0_dp, 0.10_dp]
        skip = .false.
        ret = 1
        if (any(check_gap(3, group, qe0, qe1, qe2, skip) .neqv. [.false., .true., .false., .false.])) &
        &   return
        ! Skipped states are not checked.
        ret = 2
        skip(2) = .true.
        if (any(check_gap(3, group, qe0, qe1, qe2, skip))) return
        ! Minimum at the first point.
        ret = 3
        skip = .false.
        if (any(check_gap(3, group, qe1, qe0, qe2, skip))) return
        ! Minimum at the last point.
        ret = 4
        if (any(check_gap(3, group, qe0, qe2, qe1, skip))) return
        ! With state 2 active, state 1 is a neighbour.
        ret = 5
        qe0 = [-0.30_dp, 0.0_dp, 0.30_dp, 0.20_dp]
        qe1 = [-0.01_dp, 0.0_dp, 0.20_dp, 0.05_dp]
        qe2 = [-0.30_dp, 0.0_dp, 0.10_dp, 0.10_dp]
        if (any(check_gap(2, group, qe0, qe1, qe2, skip) .neqv. [.true., .false., .false., .false.])) &
        &   return
        ret = 0
    end function check_gap_test


    !> Gap minima with states 1 and 3 of a group of three, at neighbouring points of a bisected
    !! window. After the first minimum is evaluated without a hop, the second must also be
    !! evaluated. State 2 has zero energy, and the gaps to states 1 and 3 are set at each point.
    integer function lz_extra_minimum_test() result(ret)
        use global_defs, only : dp
        use sh_test_mod, only : scripted_evaluator
        use sh_lz_mod, only : tsh_lz
        use tsh_method_mod, only : request_none, request_goto, request_insert
        use random_mod, only : rng_type
        use random_pcg_mod, only : rng_pcg_xsh_rr
        type(scripted_evaluator) :: pot
        type(tsh_lz) :: lz
        class(rng_type), pointer :: rng

        allocate(rng_pcg_xsh_rr :: rng)
        call rng%init(1)
        call pot%setup([3], [1])
        call lz%initialize(rng, pot, 2)
        lz%prob_conv = 0.05_dp
        lz%min_dt = 60.0_dp
        lz%en_err = 1.0e-13_dp

        ret = 1
        call pot%set_energies([-0.30_dp, 0.0_dp, 0.40_dp])
        call lz%start(0.0_dp, pot)
        call pot%set_energies([-0.10_dp, 0.0_dp, 0.20_dp])
        call lz%propagate(100.0_dp, pot)
        if (lz%request%action /= request_none .or. lz%request%t_final /= 0.0_dp) return

        ! Minimum with state 1 at t = 100 needs a smaller step.
        ret = 2
        call pot%set_energies([-0.20_dp, 0.0_dp, 0.30_dp])
        call lz%propagate(200.0_dp, pot)
        if (lz%request%action /= request_insert) return
        if (lz%request%time /= 0.0_dp .or. lz%request%dt /= 50.0_dp) return
        if (lz%request%hop .or. lz%active_state /= 2) return

        ! Windows (0, 50, 100), minimum with state 1 at t = 50, and (50, 100, 200), minimum with
        ! state 3 at t = 100. No bisection, since the step is below min_dt, and no hop.
        ret = 3
        call pot%set_energies([-0.09_dp, 0.0_dp, 0.35_dp])
        call lz%propagate(50.0_dp, pot)
        if (lz%active_state /= 2) return
        if (lz%prob(1) /= 0.0_dp .or. .not. lz%prob(3) > 0.0_dp) return
        if (lz%request%action /= request_goto .or. lz%request%time /= 200.0_dp) return
        if (lz%request%hop .or. lz%request%t_final /= 100.0_dp) return
        ret = 0
    end function lz_extra_minimum_test


    !> A hop at a gap minimum, rejected by the driver. Continuing from the minimum does not
    !! evaluate the same event again.
    integer function lz_guard_test() result(ret)
        use global_defs, only : dp
        use sh_test_mod, only : scripted_evaluator
        use sh_lz_mod, only : tsh_lz
        use tsh_method_mod, only : request_none, request_goto
        use random_mod, only : rng_type
        use random_pcg_mod, only : rng_pcg_xsh_rr
        type(scripted_evaluator) :: pot
        type(tsh_lz) :: lz
        class(rng_type), pointer :: rng

        allocate(rng_pcg_xsh_rr :: rng)
        call rng%init(1)
        call pot%setup([2], [1])
        call lz%initialize(rng, pot, 2)
        lz%en_err = 1.0e-13_dp

        ! Narrow minimum at t = 100 with hopping probability close to 1.
        ret = 1
        call pot%set_energies([-0.05_dp, 0.0_dp])
        call lz%start(0.0_dp, pot)
        call pot%set_energies([-0.001_dp, 0.0_dp])
        call lz%propagate(100.0_dp, pot)
        call pot%set_energies([-0.05_dp, 0.0_dp])
        call lz%propagate(200.0_dp, pot)
        if (.not. lz%prob(1) > 0.98_dp) return
        if (lz%request%action /= request_goto .or. lz%request%time /= 100.0_dp) return
        if (.not. lz%request%hop .or. lz%request%previous_state /= 2) return
        if (lz%active_state /= 1 .or. lz%request%t_final /= 0.0_dp) return

        ret = 2
        call lz%hop_rejected()
        if (lz%active_state /= 2 .or. lz%request%hop) return

        ! The same point again from the minimum.
        ret = 3
        call lz%propagate(200.0_dp, pot)
        if (any(lz%prob /= 0.0_dp) .or. lz%active_state /= 2) return
        if (lz%request%action /= request_none .or. lz%request%t_final /= 100.0_dp) return
        ret = 0
    end function lz_guard_test


    !> Intersystem crossing between a singlet and a triplet, from either state, given as
    !! components and as an effective state. The squared coupling summed over the components is
    !! s2 = |h|**2 at t = 0 and 4 |h|**2 at t = 100, and the crossing is at t = 25, so the
    !! interpolated s2 is 1.75 |h|**2. The gap slope is 0.04 / 100.
    integer function lz_isc_test() result(ret)
        use global_defs, only : dp
        use sh_lz_mod, only : tsh_lz
        use tsh_method_mod, only : request_none
        type(tsh_lz) :: lz
        real(dp), parameter :: h(3) = [1.0e-4_dp, -2.0e-4_dp, 0.5e-4_dp]
        real(dp), parameter :: tiny = 1.0e-12_dp
        real(dp) :: s2

        s2 = 1.75_dp * sum(h**2)

        ! S -> T. The singlet has a single component, so initial_mean is the same as rss.
        ret = 1
        if (abs(isc_prob(lz, 'components', 'rss', 1, h) - p_isc(s2)) > tiny) return
        if (lz%n_state /= 2 .or. size(lz%prob) /= 2) return
        if (abs(isc_prob(lz, 'components', 'rms', 1, h) - p_isc(s2 / 3)) > tiny) return
        if (abs(isc_prob(lz, 'components', 'initial_mean', 1, h) - p_isc(s2)) > tiny) return

        ! T -> S.
        ret = 2
        if (abs(isc_prob(lz, 'components', 'rss', 2, h) - p_isc(s2)) > tiny) return
        if (abs(isc_prob(lz, 'components', 'rms', 2, h) - p_isc(s2 / 3)) > tiny) return
        if (abs(isc_prob(lz, 'components', 'initial_mean', 2, h) - p_isc(s2 / 3)) > tiny) return

        ! The effective state with the root-sum-square coupling gives the same probabilities.
        ret = 3
        if (abs(isc_prob(lz, 'effective_rss', 'rss', 1, h) - p_isc(s2)) > tiny) return
        if (abs(isc_prob(lz, 'effective_rss', 'rms', 2, h) - p_isc(s2 / 3)) > tiny) return
        if (abs(isc_prob(lz, 'effective_rss', 'initial_mean', 2, h) - p_isc(s2 / 3)) > tiny)    &
        &   return

        ! A strong coupling gives a hop at the latest point.
        ret = 4
        if (.not. isc_prob(lz, 'components', 'rss', 1, 100 * h) > 0.999_dp) return
        if (lz%active_state /= 2 .or. .not. lz%request%hop .or. .not. lz%request%isc) return
        if (lz%request%previous_state /= 1 .or. lz%request%action /= request_none) return
        ret = 0
    end function lz_isc_test


    !> Probability of an intersystem crossing with squared coupling s2 in lz_isc_test.
    real(dp) function p_isc(s2) result(prob)
        use global_defs, only : dp, pi
        real(dp), intent(in) :: s2

        prob = 1.0_dp - exp(-2.0_dp * pi * s2 / (0.04_dp / 100.0_dp))
    end function p_isc


    !> Run the crossing of lz_isc_test from spin-free state active (1 singlet, 2 triplet) and
    !! return the probability of the hop to the other state.
    real(dp) function isc_prob(lz, spin_states, coupling, active, h) result(prob)
        use global_defs, only : dp
        use sh_test_mod, only : scripted_evaluator
        use sh_lz_mod, only : tsh_lz
        use random_mod, only : rng_type
        use random_pcg_mod, only : rng_pcg_xsh_rr
        type(tsh_lz), intent(out) :: lz
        character(len=*), intent(in) :: spin_states
        character(len=*), intent(in) :: coupling
        integer, intent(in) :: active
        real(dp), intent(in) :: h(3) !< Couplings of the singlet and the triplet components at t = 0.
        type(scripted_evaluator) :: pot
        class(rng_type), pointer :: rng

        allocate(rng_pcg_xsh_rr :: rng)
        call rng%init(1)
        call pot%setup([1, 1], [1, 3], spin_states)
        call lz%initialize(rng, pot, active)
        lz%isc_coupling = coupling
        lz%en_err = 1.0e-13_dp

        call pot%set_soc(st_soc(pot%n_state, h))
        call pot%set_energies([0.0_dp, 0.01_dp])
        call lz%start(0.0_dp, pot)
        call pot%set_soc(st_soc(pot%n_state, 2 * h))
        call pot%set_energies([0.0_dp, -0.03_dp])
        call lz%propagate(100.0_dp, pot)
        prob = lz%prob(3 - active)
    end function isc_prob


    !> Spin-orbit couplings of a singlet and a triplet, with n_state = 4 for components and 2 for
    !! an effective triplet with the root-sum-square coupling.
    function st_soc(n_state, h) result(soc)
        use global_defs, only : dp
        integer, intent(in) :: n_state
        real(dp), intent(in) :: h(3) !< Couplings of the singlet and the triplet components.
        real(dp) :: soc(n_state, n_state)

        soc = 0.0_dp
        if (n_state == 4) then
            soc(1, 2:4) = h
        else
            soc(1, 2) = norm2(h)
        end if
        soc(:, 1) = soc(1, :)
    end function st_soc


end program sh_test
