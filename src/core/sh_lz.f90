!--------------------------------------------------------------------------------------------------
! MODULE: sh_lz_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date May, 2020
!
! DESCRIPTION:
!> @brief Surface hopping algorithm based on Landau-Zener formula.
!--------------------------------------------------------------------------------------------------
module sh_lz_mod
    use global_defs
    use system_var, only : trajtype
    use constants
    implicit none

    private
    public :: lzsh


contains

    subroutine lzsh(rng, dt, qm_en_err, prob_conv, min_dt, dt_0)
        use system_var
        use random_mod, only : rng_type
        class(rng_type), allocatable, intent(inout) :: rng
        real(dp), intent(inout) :: dt
        real(dp), intent(in) :: qm_en_err
        real(dp), intent(in) :: prob_conv
        real(dp), intent(in) :: min_dt
        real(dp), intent(in) :: dt_0
        type(trajtype), pointer :: t0 !< Trajectory variables before gap minimum.
        type(trajtype), pointer :: t1 !< Trajectory variables at gap minimum.
        type(trajtype), pointer :: t2 !< Trajectory variables after gap minimum.
        type(trajtype), pointer :: t_wrk !< Temporary pointer
        logical :: need_bisect
        logical, allocatable :: check(:)
        real(dp) :: prob(3)
        real(dp) :: gap_sd(3)
        real(dp) :: g0, g1, g2, gap_err
        real(dp) :: rnum
        integer :: i
        logical :: gap_min_before_bisect
        logical :: check_hop

        if (tr1%step < 2) return
        ! Prevent re-checking at the same gap minimum if the trajectory was rewound
        ! due to a hop which was rejected due to energy conservation.
        if (tr1%substep == -2) then
            tr1%substep = -1
            i = trajectory_data(index_offset(data_index_1, -2))%cstate
            if (tr2%cstate == i) return
        end if

        gap_min_before_bisect = .true.
        allocate(check(tr1%nstate))

        select case(tr1%substep)
        case(-1)
            t0 => trajectory_data(index_offset(data_index_1, -2))
            t1 => trajectory_data(data_index_2)
            t2 => trajectory_data(data_index_1)
        case(1)
            ! Now in step 0.5 (between steps 0 and 1 from previous call).
            ! First checking for gap between steps 0, 0.5 and 1.
            t0 => trajectory_data(data_index_2)
            t1 => trajectory_data(data_index_1)
            t2 => trajectory_data(index_offset(data_index_1, 1))
            check = check_gap(t0, t1, t2)
            if (.not. any(check)) then
                ! Now checking for gap between steps 0.5, 1 and 2.
                t_wrk => trajectory_data(data_index_2)
                t0 => trajectory_data(data_index_1)
                t1 => trajectory_data(index_offset(data_index_1, 1))
                t2 => trajectory_data(index_offset(data_index_1, 2))
                gap_min_before_bisect = .false.
            end if
        case(2)
            ! Now in step 1.5 (between steps 1 and 2 from previous call).
            ! First checking for gap between steps 0, 1 and 1.5.
            t0 => trajectory_data(index_offset(data_index_1, -2))
            t1 => trajectory_data(data_index_2)
            t2 => trajectory_data(data_index_1)
            check = check_gap(t0, t1, t2)
            if (.not. any(check)) then
                ! Now checking for gap between steps 1, 1.5 and 2.
                t_wrk => trajectory_data(index_offset(data_index_1, -2))
                t0 => trajectory_data(data_index_2)
                t1 => trajectory_data(data_index_1)
                t2 => trajectory_data(index_offset(data_index_1, 1))
                gap_min_before_bisect = .false.
            end if
        end select

        check = check_gap(t0, t1, t2)
        if ((tr1%substep > 0) .and. (.not. any(check))) then
            write(stderr, *) 'Warning. Gap minimum not found after adding an extra time step.'
            write(stderr, '(999(e24.16, 1x))') t_wrk%time, t_wrk%qe
            write(stderr, '(999(e24.16, 1x))') t0%time, t0%qe
            write(stderr, '(999(e24.16, 1x))') t1%time, t1%qe
            write(stderr, '(999(e24.16, 1x))') t2%time, t2%qe
        end if

        if (.not. any(check)) return

500     need_bisect = .false.

        ! Evaluate hopping probability and decide whether the time step should be reduced
        do i = 1, t1%nstate
            if (.not. check(i)) cycle
            g0 = t0%qe(t1%cstate) - t0%qe(i)
            g1 = t1%qe(t1%cstate) - t1%qe(i)
            g2 = t2%qe(t1%cstate) - t2%qe(i)
            gap_err = abs((g0 - 2*g1 + g2) / 2)
            if ((t1%gap_2deriv(2, i) == 0.0_dp) .or. (gap_err > 20 * qm_en_err)) then
                ! Calculate 2nd derivative of the gap if it hasn't already been calculated
                ! for this pair of states at this gap minimum.
                ! Also re-calculate the 2nd derivative of the gap if the energies of the states
                ! have changed significantly enough with respect to the convergence threshold
                ! for the energy.
                gap_sd = sec_deriv_3p_err(t0%time, t1%time, t2%time, g0, g1, g2, qm_en_err)
                t1%gap_2deriv(:, i) = gap_sd
            else
                ! Otherwise, keep previously calculated 2nd derivative since the random errors
                ! in the energies due to the convergence threshold might cause larger errors
                ! in the calculated second derivative.
                continue
            end if
            call lz_prob_err_gap(g1, t1%gap_2deriv(2, i), gap_err, qm_en_err, prob)
        !   call lz_prob_err_both(g1, t1%gap_2deriv(:, i), gap_err, qm_en_err, prob)
            t1%prob(i) = prob(2)
            if (stdp3) then
                write(stdout, '(3x,a)') 'LZSH probability estimates:'
                write(stdout, '(5x,a,e15.7)') 'Pmin = ', prob(1)
                write(stdout, '(5x,a,e15.7)') 'P = ', prob(2)
                write(stdout, '(5x,a,e15.7)') 'Pmax = ', prob(3)
            end if
            if (prob(3) - prob(1) > prob_conv) then
                if (min_dt >= dt) then
                    if (stdp3) write(stdout, '(3x,a)') 'LZSH not bisecting due to time step.'
                    cycle
                end if
                if (gap_err < 2 * qm_en_err) then
                    if (stdp3) write(stdout, '(3x,a)') 'LZSH not bisecting due to qm_en_err.'
                    cycle
                end if
                need_bisect = .true.
                ! Not returning from the subroutine right away so we can also calculate the
                ! second derivative of the gap between the current state and another state if
                ! needed.
            end if
        !   call lz_prob_err_both(g1, gap_sd, gap_err, qm_en_err, prob)
        end do

        if (need_bisect) then
            if (1.5*(t1%time - t0%time) > (t2%time - t1%time)) then
                ! Set t0 as the current step and slide t1 and t2 forward in the array
                ! so they don't get overwritten by the extra step.
                dt = 0.5_dp * (t1%time - t0%time)
                call trajectory_slide_forward(t0%step, 2, .true.)
                tr1%substep = 1
                tr1%step = maxval(trajectory_data(:)%step) + 1
            else
                ! Set t1 as the current step and slide t2 forward in the array
                ! so it doesn't get overwritten by the extra step.
                dt = 0.5_dp * (t2%time - t1%time)
                call trajectory_slide_forward(t1%step, 1, .true.)
                tr1%substep = 2
                tr1%step = maxval(trajectory_data(:)%step) + 1
            end if
            if (stdp1) then
                write(stdout, '(3x,a, i0)') 'Adding new substep: ', tr1%substep
                write(stdout, '(5x,a,f10.4)') 't=', (tr1%time + dt) * aut_fs
            end if
            return
        end if

        ! Check if a hop should occur.
        if (sum(t1%prob) > 1.0_dp) then
            write(stderr, *) ' Warning. Sum of hopping probabilities for all states higher than 1.'
            write(stderr, *) '   time:', t1%time
            write(stderr, *) '   cstate:', t1%cstate
            write(stderr, *) '   fprob:', t1%prob
        end if
        check_hop = .false.
        prob(2) = 0.0_dp
        call rng%uniform(rnum)
        do i = 1, size(t1%prob)
            prob(2) = prob(2) + t1%prob(i)
            if (rnum < prob(2)) then
                check_hop = .true.
                t1%cstate = i
                exit
            end if
        end do

        if (check_hop) then
            if (stdp1) then
                write(stdout, '(3x,a,f0.4,a)') 'Hop occurred, resuming trajectory from t=', &
                &                               t1%time * aut_fs, '.'
            end if
            call trajectory_set_current(t1%step)
            tr1%substep = -2
        else
            if (stdp2) write(stdout, '(3x,a)') 'Hop probability evaluated, no hop.'
            if (tr1%substep > 0) then
                if (gap_min_before_bisect) then
                    if (tr1%substep == 1) then
                        t0 => trajectory_data(data_index_1)
                        t1 => trajectory_data(index_offset(data_index_1, 1))
                        t2 => trajectory_data(index_offset(data_index_1, 2))
                    else if (tr1%substep == 2) then
                        t0 => trajectory_data(data_index_2)
                        t1 => trajectory_data(data_index_1)
                        t2 => trajectory_data(index_offset(data_index_1, 1))
                    end if
                    check = check_gap(t0, t1, t2)
                    if (any(check)) then
                        if (stdp2) write(stdout, '(3x,a)') ' Extra gap minimum in same step.'
                        if (stdp2) write(stdout, '(3x,a)') ' Running LZSH procedure for new steps.'
                        goto 500
                    end if
                end if
            end if
            if (stdp1) then
                write(stdout, '(3x,a,f0.4,a)') 'Resuming trajectory from t=', t2%time * aut_fs, '.'
            end if
            call trajectory_set_current(t2%step)
            tr1%substep = -1
        end if

        ! Reset parameters changed by the adaptive-step LZSH algorithm
        dt = dt_0
        tr1%gap_2deriv = 0.0_dp
    end subroutine lzsh


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: check_gap
    !> @brief Check for a gap minimum with the active state during the previous three time steps.
    !----------------------------------------------------------------------------------------------
    function check_gap(t0, t1, t2) result(check)
        type(trajtype), intent(in) :: t0
        type(trajtype), intent(inout) :: t1
        type(trajtype), intent(in) :: t2
        logical :: check(t1%nstate)
        integer :: i, cstate, pstate
        real(dp) :: g0, g1, g2

        pstate = t0%cstate
        cstate = t1%cstate
        check = .false.
        t1%prob = 0.0_dp
        do i = 1, t1%nstate
            if (i == cstate) cycle
            if (abs(i-cstate) /= 1) cycle ! Only allow hops to neighbouring states.
            if (i == pstate) cycle ! Prevent hop back to same state after rewinding.
            g0 = t0%qe(cstate) - t0%qe(i)
            g1 = t1%qe(cstate) - t1%qe(i)
            g2 = t2%qe(cstate) - t2%qe(i)
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
    subroutine lz_prob_err_both(g1, sd, gap_err, tol_err, prob)
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
    subroutine lz_prob_err_gap(g1, sd, gap_err, tol_err, prob)
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
    function lz_prob(gap, sd) result(prob)
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
    function sec_deriv_3p_err(t0, t1, t2, g0, g1, g2, tol_err) result(sd)
        real(dp), intent(in) :: t0 !< Time for which g0 is evaluated
        real(dp), intent(in) :: t1 !< Time at gap minimum
        real(dp), intent(in) :: t2 !< Time for which g2 is evaluated
        real(dp), intent(in) :: g0 !< Finite difference gap before minimum
        real(dp), intent(in) :: g1 !< Gap at gap minimum
        real(dp), intent(in) :: g2 !< Finite difference gap after minimum
        real(dp), intent(in) :: tol_err !< Assumed error in all energies due to precision
        real(dp) :: sd(3) !< Calculated numerical second derivatives assuming:
                          !< (1) errors aligning to give maximum 2nd deriv. (min probability),
                          !< (2) no errors in energies,
                          !< (3) errors aligning to give minimum 2nd deriv. (max probability).
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
    function sec_deriv_3p(x1, x2, x3, y1, y2, y3) result(sd)
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
