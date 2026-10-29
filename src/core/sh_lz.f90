!--------------------------------------------------------------------------------------------------
! MODULE: sh_lz_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date May, 2020
!
! DESCRIPTION:
!> @brief Surface hopping algorithm based on Landau-Zener formula.
!> @details
!! A hop between two neighbouring states of the same group is considered at each minimum of the
!! gap between the active state and the other state. A minimum is found at the middle point of the
!! last three points. The hopping probability is given by the adiabatic Landau-Zener formula, with
!! the second derivative of the gap from the three points.
!!
!! If the estimated error of the probability is above prob_conv, the method requests a new point
!! between the stored points, halving the larger step next to the minimum. After a hop, the method
!! requests the trajectory to continue from the minimum.
!!
!! The dynamics runs in the group_adiabatic basis. For a single group of states this is the
!! adiabatic basis.
!--------------------------------------------------------------------------------------------------
module sh_lz_mod
    use global_defs
    use constants
    use tsh_method_mod, only : tsh_method, tsh_request, base_initialize => initialize,          &
    &                          request_goto, request_insert
    use evaluator_base_mod, only : potential_evaluator, state_location
    use random_mod, only : rng_type
    implicit none

    private
    public :: tsh_lz
    public :: check_gap
    public :: lz_prob
    public :: lz_prob_err_gap
    public :: lz_prob_err_both
    public :: sec_deriv_3p
    public :: sec_deriv_3p_err


    !> Point of the trajectory stored by the method.
    type :: lz_point
        real(dp) :: time = 0.0_dp
        real(dp), allocatable :: h(:, :) !< Hamiltonian in the basis of the dynamics.
    end type lz_point


    !> Gap minimum at which the hopping probability between two states was evaluated.
    type :: lz_event
        real(dp) :: time = 0.0_dp !< Time of the gap minimum.
        integer :: pair(2) = 0 !< The two states.
    end type lz_event


    type, extends(tsh_method) :: tsh_lz
        real(dp) :: prob_conv = 0.05_dp !< Convergence threshold for the estimated error of the
                                        !! probabilities.
        real(dp) :: min_dt = 0.0_dp !< Minimum time step when converging the probabilities.
        real(dp) :: en_err = 1.0e-6_dp !< Estimated error in the energies.
        integer, allocatable :: group(:) !< Group of each state.
        type(lz_point), allocatable :: point(:) !< Stored points in time order.
        integer :: i_cur = 0 !< Position of the current point in point.
        real(dp) :: t_window = 0.0_dp !< Time of the first point of the last bisected window.
        real(dp), allocatable :: gap_2deriv(:, :) !< Second derivatives of the gaps with each state
                                                  !! (see sec_deriv_3p_err) at the gap minimum
                                                  !! being evaluated.
        type(lz_event), allocatable :: guard(:) !< Events which are not evaluated again.
    contains
        procedure :: initialize => lz_initialize
        procedure :: start => lz_start
        procedure :: propagate => lz_propagate
        procedure :: max_points => lz_max_points
        procedure, private :: check_window
        procedure, private :: evaluate_window
    end type tsh_lz


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: lz_initialize
    !> @brief Set default options and the group of each state.
    !----------------------------------------------------------------------------------------------
    subroutine lz_initialize(self, rng, potential, active_state)
        class(tsh_lz), intent(inout) :: self
        class(rng_type), pointer, intent(in) :: rng
        class(potential_evaluator), intent(in) :: potential
        integer, intent(in) :: active_state
        type(state_location) :: loc
        integer :: i

        call base_initialize(self, rng, potential, active_state)
        self%basis = 'group_adiabatic'
        self%min_dt = 0.05_dp / aut_fs
        allocate(self%group(self%n_state))
        do i = 1, self%n_state
            loc = potential%locate_state(i)
            self%group(i) = loc%group
        end do
        allocate(self%gap_2deriv(3, self%n_state), source=0.0_dp)
        allocate(self%point(0))
        allocate(self%guard(0))
        if (potential%n_group > 1) then
            write(stderr, *) 'Warning in sh_lz_mod, lz_initialize subroutine.'
            write(stderr, *) '  Only hops between states of the same group are considered.'
        end if
    end subroutine lz_initialize


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: lz_max_points
    !> @brief Largest number of points stored at once, for output grid step dt_0.
    !> @details
    !! Each bisection adds a point to the three points of the window. The bisected step is at least
    !! 2/3 of the larger step in the window, so the larger step halves at least every second
    !! bisection, and steps are not bisected below min_dt.
    !----------------------------------------------------------------------------------------------
    integer function lz_max_points(self, dt_0) result(n)
        class(tsh_lz), intent(in) :: self
        real(dp), intent(in) :: dt_0 !< Time step of the output grid.
        integer, parameter :: max_halving = 50
        integer :: n_halving

        n_halving = max_halving
        if (self%min_dt > 0.0_dp) then
            n_halving = min(n_halving, max(0, ceiling(log(dt_0 / self%min_dt) / log(2.0_dp))))
        end if
        n = 2 * n_halving + 5
    end function lz_max_points


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: lz_start
    !> @brief Store the initial point.
    !----------------------------------------------------------------------------------------------
    subroutine lz_start(self, t, potential)
        class(tsh_lz), intent(inout) :: self
        real(dp), intent(in) :: t
        class(potential_evaluator), intent(in) :: potential

        self%t = t
        self%point = [lz_point(t, potential%get_Hamiltonian(self%basis))]
        self%i_cur = 1
        self%request = tsh_request(t_final=t)
    end subroutine lz_start


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: lz_propagate
    !> @brief Store the new point and evaluate the hopping probabilities at a gap minimum.
    !> @details
    !! Usually the new point is the latest one, and the window of the last three points is checked
    !! for a gap minimum. After a requested insertion, the new point lies between the four points
    !! q1 < q2 < q3 < q4 made of the bisected window and the new point. The windows (q1, q2, q3)
    !! and (q2, q3, q4) are then checked in turn. After a gap minimum in the first one is evaluated
    !! without a hop, the second one is also checked.
    !----------------------------------------------------------------------------------------------
    subroutine lz_propagate(self, t2, potential)
        class(tsh_lz), intent(inout) :: self
        real(dp), intent(in) :: t2 !< Time of the new point.
        class(potential_evaluator), intent(in) :: potential
        logical :: inserted !< The new point lies between stored points.
        logical :: min_first !< Gap minimum in the first window after an insertion.
        logical :: need_bisect
        logical :: check(self%n_state)
        integer :: iw !< Position of the first point of the window.
        integer :: i, cstate

        cstate = self%active_state
        self%request = tsh_request(previous_state=cstate)
        self%prob = 0.0_dp
        self%t = t2

        inserted = (self%i_cur < size(self%point))
        if (.not. inserted) then
            ! Only the last two points are needed for the next window.
            if (self%i_cur > 2) then
                self%point = self%point(self%i_cur - 1:)
                self%i_cur = 2
                self%guard = pack(self%guard, self%guard%time >= self%point(1)%time)
            end if
        end if
        self%point = [self%point(:self%i_cur),                                                    &
        &             lz_point(t2, potential%get_Hamiltonian(self%basis)),                        &
        &             self%point(self%i_cur + 1:)]
        self%i_cur = self%i_cur + 1
        if (size(self%point) < 3) then
            call finish()
            return
        end if

        min_first = .false.
        if (inserted) then
            do iw = 1, size(self%point)
                if (self%point(iw)%time == self%t_window) exit
            end do
            check = self%check_window(iw)
            min_first = any(check)
            if (.not. min_first) then
                iw = iw + 1
                check = self%check_window(iw)
                if (.not. any(check)) then
                    write(stderr, *) 'Warning. Gap minimum not found after adding an extra time step.'
                    do i = iw - 1, iw + 2
                        write(stderr, '(999(e24.16, 1x))') self%point(i)%time,                     &
                        &                                  diagonal(self%point(i)%h)
                    end do
                end if
            end if
        else
            iw = size(self%point) - 2
            check = self%check_window(iw)
        end if

        if (any(check)) then
            do
                call self%evaluate_window(iw, check, need_bisect)
                if (need_bisect) then
                    call bisect(iw)
                    call finish()
                    return
                end if

                ! Check if a hop should occur.
                if (sum(self%prob) > 1.0_dp) then
                    write(stderr, *) ' Warning. Sum of hopping probabilities for all states higher than 1.'
                    write(stderr, *) '   time:', self%point(iw + 1)%time
                    write(stderr, *) '   cstate:', cstate
                    write(stderr, *) '   fprob:', self%prob
                end if
                call self%select_state()
                if (self%active_state /= cstate) exit

                if (stdp2) write(stdout, '(3x,a)') 'Hop probability evaluated, no hop.'
                if (.not. min_first) exit
                min_first = .false.
                iw = iw + 1
                check = self%check_window(iw)
                if (.not. any(check)) exit
                if (stdp2) write(stdout, '(3x,a)') ' Extra gap minimum in same step.'
                if (stdp2) write(stdout, '(3x,a)') ' Running LZSH procedure for new steps.'
            end do
        end if

        if (self%active_state /= cstate) then
            ! Continue from the gap minimum in the new state.
            if (stdp1) then
                write(stdout, '(3x,a,f0.4,a)') 'Hop occurred, resuming trajectory from t=',        &
                &                               self%point(iw + 1)%time * aut_fs, '.'
            end if
            self%guard = [self%guard, lz_event(self%point(iw + 1)%time, [cstate, self%active_state])]
            self%point = self%point(:iw + 1)
            self%i_cur = iw + 1
            self%request%action = request_goto
            self%request%time = self%point(self%i_cur)%time
            self%request%hop = .true.
        else if (any(check) .or. inserted) then
            ! Continue from the latest point.
            if (stdp1) then
                write(stdout, '(3x,a,f0.4,a)') 'Resuming trajectory from t=',                      &
                &                               self%point(size(self%point))%time * aut_fs, '.'
            end if
            if (self%i_cur < size(self%point)) then
                self%i_cur = size(self%point)
                self%request%action = request_goto
                self%request%time = self%point(self%i_cur)%time
            end if
        end if
        self%gap_2deriv = 0.0_dp
        call finish()

    contains

        !> Points up to the one before the current point are final. A gap minimum at the current
        !! point is only found with the next point.
        subroutine finish()
            if (self%i_cur > 1) self%request%t_final = self%point(self%i_cur - 1)%time
        end subroutine finish

        !> Request a new point which halves the step before the middle of window iw, or the step
        !! after it if that one is more than 1.5 times longer.
        subroutine bisect(iw)
            integer, intent(in) :: iw
            real(dp) :: dt
            integer :: substep

            associate(t0 => self%point(iw)%time, t1 => self%point(iw + 1)%time,                    &
            &         t2 => self%point(iw + 2)%time)
                if (1.5_dp * (t1 - t0) > (t2 - t1)) then
                    dt = 0.5_dp * (t1 - t0)
                    substep = 1
                else
                    dt = 0.5_dp * (t2 - t1)
                    substep = 2
                end if
                self%t_window = t0
            end associate
            self%i_cur = iw + substep - 1
            self%request%action = request_insert
            self%request%time = self%point(self%i_cur)%time
            self%request%dt = dt
            if (stdp1) then
                write(stdout, '(3x,a, i0)') 'Adding new substep: ', substep
                write(stdout, '(5x,a,f10.4)') 't=', (self%request%time + dt) * aut_fs
            end if
        end subroutine bisect

    end subroutine lz_propagate


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: check_window
    !> @brief Check for gap minima with the active state at the middle of the window starting at
    !! position iw, skipping events in the guard list.
    !----------------------------------------------------------------------------------------------
    function check_window(self, iw) result(check)
        class(tsh_lz), intent(in) :: self
        integer, intent(in) :: iw
        logical :: check(self%n_state)
        logical :: skip(self%n_state)
        integer :: i, cstate

        cstate = self%active_state
        skip = .false.
        do i = 1, size(self%guard)
            if (self%guard(i)%time /= self%point(iw + 1)%time) cycle
            if (self%guard(i)%pair(1) == cstate) skip(self%guard(i)%pair(2)) = .true.
            if (self%guard(i)%pair(2) == cstate) skip(self%guard(i)%pair(1)) = .true.
        end do
        check = check_gap(cstate, self%group, diagonal(self%point(iw)%h),                         &
        &                 diagonal(self%point(iw + 1)%h), diagonal(self%point(iw + 2)%h), skip)
    end function check_window


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: evaluate_window
    !> @brief Calculate the hopping probabilities at the gap minima at the middle of the window
    !! starting at position iw, and decide whether the step should be reduced.
    !----------------------------------------------------------------------------------------------
    subroutine evaluate_window(self, iw, check, need_bisect)
        class(tsh_lz), intent(inout) :: self
        integer, intent(in) :: iw
        logical, intent(in) :: check(:) !< States with a gap minimum.
        logical, intent(out) :: need_bisect
        real(dp) :: qe0(self%n_state), qe1(self%n_state), qe2(self%n_state)
        real(dp) :: t0, t1, t2
        real(dp) :: prob(3)
        real(dp) :: g0, g1, g2, gap_err
        real(dp) :: dt !< Last time step.
        integer :: i, cstate

        cstate = self%active_state
        qe0 = diagonal(self%point(iw)%h)
        qe1 = diagonal(self%point(iw + 1)%h)
        qe2 = diagonal(self%point(iw + 2)%h)
        t0 = self%point(iw)%time
        t1 = self%point(iw + 1)%time
        t2 = self%point(iw + 2)%time
        dt = self%point(self%i_cur)%time - self%point(self%i_cur - 1)%time

        need_bisect = .false.
        self%prob = 0.0_dp
        ! Evaluate hopping probability and decide whether the time step should be reduced
        do i = 1, self%n_state
            if (.not. check(i)) cycle
            g0 = qe0(cstate) - qe0(i)
            g1 = qe1(cstate) - qe1(i)
            g2 = qe2(cstate) - qe2(i)
            gap_err = abs((g0 - 2*g1 + g2) / 2)
            if ((self%gap_2deriv(2, i) == 0.0_dp) .or. (gap_err > 20 * self%en_err)) then
                ! Calculate 2nd derivative of the gap if it hasn't already been calculated
                ! for this pair of states at this gap minimum.
                ! Also re-calculate the 2nd derivative of the gap if the energies of the states
                ! have changed significantly enough with respect to the convergence threshold
                ! for the energy.
                self%gap_2deriv(:, i) = sec_deriv_3p_err(t0, t1, t2, g0, g1, g2, self%en_err)
            else
                ! Otherwise, keep previously calculated 2nd derivative since the random errors
                ! in the energies due to the convergence threshold might cause larger errors
                ! in the calculated second derivative.
                continue
            end if
            call lz_prob_err_gap(g1, self%gap_2deriv(2, i), gap_err, self%en_err, prob)
            self%prob(i) = prob(2)
            if (stdp3) then
                write(stdout, '(3x,a)') 'LZSH probability estimates:'
                write(stdout, '(5x,a,e15.7)') 'Pmin = ', prob(1)
                write(stdout, '(5x,a,e15.7)') 'P = ', prob(2)
                write(stdout, '(5x,a,e15.7)') 'Pmax = ', prob(3)
            end if
            if (prob(3) - prob(1) > self%prob_conv) then
                if (self%min_dt >= dt) then
                    if (stdp3) write(stdout, '(3x,a)') 'LZSH not bisecting due to time step.'
                    cycle
                end if
                if (gap_err < 2 * self%en_err) then
                    if (stdp3) write(stdout, '(3x,a)') 'LZSH not bisecting due to qm_en_err.'
                    cycle
                end if
                need_bisect = .true.
                ! Not returning from the subroutine right away so we can also calculate the
                ! second derivative of the gap between the current state and another state if
                ! needed.
            end if
        end do
    end subroutine evaluate_window


    pure function diagonal(h) result(d)
        real(dp), intent(in) :: h(:, :)
        real(dp) :: d(size(h, 1))
        integer :: i

        do i = 1, size(h, 1)
            d(i) = h(i, i)
        end do
    end function diagonal


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: check_gap
    !> @brief Check for a gap minimum with the active state at the middle of three points.
    !> @details
    !! Only states in the same group as the active state, whose index differs from it by 1, are
    !! checked.
    !----------------------------------------------------------------------------------------------
    pure function check_gap(cstate, group, qe0, qe1, qe2, skip) result(check)
        integer, intent(in) :: cstate !< Current active state
        integer, intent(in) :: group(:) !< Group of each state.
        real(dp), intent(in) :: qe0(:) !< Energies at t0
        real(dp), intent(in) :: qe1(:) !< Energies at t1
        real(dp), intent(in) :: qe2(:) !< Energies at t2
        logical, intent(in) :: skip(:) !< States which are not checked.
        logical :: check(size(qe1))
        integer :: i
        real(dp) :: g0, g1, g2

        check = .false.
        do i = 1, size(qe1)
            if (i == cstate) cycle
            if (abs(i-cstate) /= 1) cycle ! Only allow hops to neighbouring states.
            if (group(i) /= group(cstate)) cycle
            if (skip(i)) cycle ! Prevent evaluating the same event again after rewinding.
            g0 = qe0(cstate) - qe0(i)
            g1 = qe1(cstate) - qe1(i)
            g2 = qe2(cstate) - qe2(i)
            if ((abs(g1) < abs(g2)) .and. (abs(g1) < abs(g0))) then
                check(i) = .true.
            end if
        end do
    end function check_gap


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: lz_prob_err_both
    !
    !> @brief Expected error interval for using the LZ formula.
    !> @details
    !! Due to discreet time steps, the true gap minimum is not the lowest calculated point so an
    !! error proportional to the change in the gap between the time steps is assumed.
    !! Additionally, the energies are usually calculated numerically only up to a specified
    !! precision which affects both the estimated gap minimum and the numerical evaluation of the
    !! second derivative of the gap.
    !> @note The sd variable passed to this subroutine is assumed to be calculated by the
    !!       function sec_deriv_3p_err.
    !----------------------------------------------------------------------------------------------
    pure subroutine lz_prob_err_both(g1, sd, gap_err, tol_err, prob)
        real(dp), intent(in) :: g1 !< Gap at time-step with lowest gap
        real(dp), intent(in) :: sd(3) !< Estimated range for the second deriv. at the gap minimum
        real(dp), intent(in) :: gap_err !< Assumed error in gap minimum due to finite step
        real(dp), intent(in) :: tol_err !< Assumed error in all energies due to precision
        real(dp), intent(out) :: prob(3) !< Probabilities calculated assuming:
                                         !< (1) errors aligning to give minimum probability,
                                         !< (2) no errors
                                         !< (3) errors aligning to give maximum probability.
        real(dp) :: err_gap
        real(dp), parameter :: tinydp = 1.0e-15_dp

        ! Minimum probability
        err_gap = abs(g1) + tol_err
        prob(1) = lz_prob(err_gap, abs(sd(1)))

        ! No error probability
        prob(2) = lz_prob(abs(g1), abs(sd(2)))

        ! Maximum probability
        err_gap = max(abs(g1) - gap_err - tol_err, tinydp)
        prob(3) = lz_prob(err_gap, abs(sd(3)))
    end subroutine lz_prob_err_both


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: lz_prob_err_gap
    !
    !> @brief Expected error interval for using the LZ formula due to accuracy of gap minimum.
    !> @details
    !! Due to discreet time steps, the true gap minimum is not the lowest calculated point so an
    !! error proportional to the change in the gap between the time steps is assumed.
    !----------------------------------------------------------------------------------------------
    pure subroutine lz_prob_err_gap(g1, sd, gap_err, tol_err, prob)
        real(dp), intent(in) :: g1 !< Gap at gap minimum
        real(dp), intent(in) :: sd !< Second derivative at gap minimum
        real(dp), intent(in) :: gap_err !< Assumed error in gap minimum due to finite step
        real(dp), intent(in) :: tol_err !< Assumed error in all energies due to precision
        real(dp), intent(out) :: prob(3) !< Probabilities calculated assuming (1) errors aligning to
            !! give min probability, (2) no errors and (3) errors aligning for max probability
        real(dp) :: err_gap
        real(dp), parameter :: tinydp = 1.0e-15_dp

        ! Minimum probability
        err_gap = abs(g1) + tol_err
        prob(1) = lz_prob(err_gap, abs(sd))

        ! No error probability
        prob(2) = lz_prob(abs(g1), abs(sd))

        ! Maximum probability
        err_gap = max(abs(g1) - gap_err - tol_err, tinydp)
        prob(3) = lz_prob(err_gap, abs(sd))
    end subroutine lz_prob_err_gap


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: lz_prob
    !> @brief Return the Landau-Zener hopping probability
    !----------------------------------------------------------------------------------------------
    pure function lz_prob(gap, sd) result(prob)
        real(dp), intent(in) :: gap !< Gap at gap minimum
        real(dp), intent(in) :: sd !< Second derivative at gap minimum
        real(dp) :: prob
        prob = exp(-pio2 * sqrt(gap**3 / sd))
    end function lz_prob


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: sec_deriv_3p_err
    !
    !> @brief Expected error interval for using the LZ formula.
    !> @details
    !! When using finite time steps and calculating the second derivative numerically, the energy
    !! at the gap minimum and the differences in the gap energies used for the finite difference
    !! method are not known exactly. For both, the error is determined by the precision of the
    !! numerical procedure used to calculate energies. Additionaly, the true gap minimum is not
    !! the lowest calculated point so an additional error proportional to the change in the gap
    !! between the time steps is assumed.
    !----------------------------------------------------------------------------------------------
    pure function sec_deriv_3p_err(t0, t1, t2, g0, g1, g2, tol_err) result(sd)
        real(dp), intent(in) :: t0 !< Time for which g0 is evaluated
        real(dp), intent(in) :: t1 !< Time at gap minimum
        real(dp), intent(in) :: t2 !< Time for which g2 is evaluated
        real(dp), intent(in) :: g0 !< Finite difference gap before minimum
        real(dp), intent(in) :: g1 !< Gap at gap minimum
        real(dp), intent(in) :: g2 !< Finite difference gap after minimum
        real(dp), intent(in) :: tol_err !< Assumed error in all energies due to precision
        real(dp) :: sd(3) !< Calculated numerical second derivatives assuming:
                          !< (1) errors aligning to give minimum 2nd deriv. (min probability),
                          !< (2) no errors in energies,
                          !< (3) errors aligning to give maximum 2nd deriv. (max probability).
        real(dp) :: err_g0, err_g1, err_g2
        real(dp), parameter :: tinydp = 1.0e-15_dp

        ! Minimum probability
        err_g0 = abs(g0) - tol_err
        err_g1 = abs(g1) + tol_err
        err_g2 = abs(g2) - tol_err
        if ((err_g0 < err_g1) .or. (err_g2 < err_g1)) then
            sd(1) = 0.0_dp
        else
            sd(1) = sec_deriv_3p(t0, t1, t2, err_g0, err_g1, err_g2)
        end if

        ! No error probability
        sd(2) = sec_deriv_3p(t0, t1, t2, g0, g1, g2)

        ! Maximum probability
        err_g0 = abs(g0) + tol_err
        err_g1 = max(abs(g1) - tol_err, tinydp)
        err_g2 = abs(g2) + tol_err
        sd(3) = sec_deriv_3p(t0, t1, t2, err_g0, err_g1, err_g2)
    end function sec_deriv_3p_err


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: sec_deriv_3p
    !> @brief Calculate numerical second derivative for three unevenly spaced points.
    !----------------------------------------------------------------------------------------------
    pure function sec_deriv_3p(x1, x2, x3, y1, y2, y3) result(sd)
        real(dp), intent(in) :: x1
        real(dp), intent(in) :: x2
        real(dp), intent(in) :: x3
        real(dp), intent(in) :: y1
        real(dp), intent(in) :: y2
        real(dp), intent(in) :: y3
        real(dp) :: sd
        real(dp) :: d21, d31, d32
        d31 = x3 - x1
        d32 = x3 - x2
        d21 = x2 - x1
        sd = 2 * y1 / d21 / d31 - 2 * y2 / d32 / d21 + 2 * y3 / d32 / d31
    end function sec_deriv_3p


end module sh_lz_mod
