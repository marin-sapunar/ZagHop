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

    subroutine lzsh()
        use control_var, only : ctrl
        use system_var
        type(trajtype) :: t0, t1, t2
        logical :: need_bisect
        logical, allocatable :: check(:)
        logical, allocatable :: sd_converged(:)
        real(dp), allocatable :: gap_sd(:)
        real(dp) :: prob(3)
        real(dp) :: sd(3)
        real(dp) :: g0, g1, g2, gap_err
        integer :: i, istate
        logical :: err_check

        allocate(sd_converged(t1%nstate), source=.false.)
        allocate(check(t1%nstate), source=.false.)
        allocate(gap_sd(t1%nstate))
        if (t(1)%step < 2) return
        if (.not. any(check_gap(t(3), t(2), t(1)))) return

        istate = t(2)%cstate
        t0 = t(3)
        t1 = t(2)
500     t2 = t(1)

        do
            need_bisect = .false.
            check = check_gap(t(3), t(2), t(1))
            do i = 1, t1%nstate
                if (.not. check(i)) cycle
                g0 = t0%qe(istate) - t0%qe(i)
                g1 = t1%qe(istate) - t1%qe(i)
                g2 = t2%qe(istate) - t2%qe(i)
                gap_err = (g0 - 2*g1 + g2) / 2
                if (sd_converged(i)) then
                    call lz_prob_interval_gap_only(g1, gap_sd(i), gap_err, ctrl%qm_en_err, prob)
                else
                    call lz_prob_interval(t0%time, t1%time, t2%time, g0, g1, g2, gap_err, &
                    &                     ctrl%qm_en_err, sd, prob)
                    if (gap_err < 10 * ctrl%qm_en_err) then
                        sd_converged(i) = .true.
                        gap_sd(i) = sd(2)
                    end if
                end if
                t1%prob(i) = prob(2)
                if (prob(3) - prob(1) < ctrl%lz_prob_conv) then
                    if (sd_converged(i)) then
                        call lz_prob_interval(t0%time, t1%time, t2%time, g0, g1, g2, gap_err,      &
                        &                     0.0_dp, sd, prob)
                        if (prob(3) - prob(1) < ctrl%lz_prob_conv) then
                            need_bisect = .true.
                        end if
                    else
                        need_bisect = .true.
                    end if
                end if
            end do
            if (.not. need_bisect) exit
            call bisect_gap(t0, t1, t2, err_check)
        end do

        if (check_hop(t1)) then
            t(3) = t0
            t(2) = t1
            t(1) = t2
            write(stdout, '(3x,a,f8.4,a)') ' Hop occurred, resuming trajectory from t=', &
            &                               t(2)%time * aut_fs, '.'
            call trajectory_rewind(t, 1, .false.)
            t(2)%cstate = istate
        else
            write(stdout, '(3x,a)') 'Hop probability evaluated, no hop.'
            if (any(check_gap(t1, t2, t(1)))) then
                write(stdout, '(3x,a)') ' Extra gap minimum in same step.'
                t0 = t1
                t1 = t2
                goto 500
            else
                t(3) = t0
                t(2) = t1
            end if
        end if
    end subroutine lzsh

    function check_gap(t0, t1, t2) result(check)
      ! use control_var, only : ctrl
      ! use constants
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
          ! if ((i > cstate) .and. (pstate > cstate)) cycle ! Prevent hopping back up.
          ! if ((i < cstate) .and. (pstate < cstate)) cycle ! Prevent hopping back down.
            if (i == pstate) cycle ! Prevent hop back to same state after rewinding.
            g0 = t0%qe(cstate) - t0%qe(i)
            g1 = t1%qe(cstate) - t1%qe(i)
            g2 = t2%qe(cstate) - t2%qe(i)
            if ((abs(g1) < abs(g2)) .and. (abs(g1) < abs(g0))) then
                check(i) = .true.
            end if
        end do

1002 format (3(f11.5,1x),2x,i3,2x,i3,4f18.10)
    end function check_gap


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: check_hop
    !> @brief Evaluate if hop should occur based on hopping probability for t1.
    !----------------------------------------------------------------------------------------------
    function check_hop(t1) result(check)
        type(trajtype), intent(inout) :: t1
        logical :: check
        real(dp) :: prob, rnum
        integer :: i

        if (sum(t1%prob) > 1.0_dp ) then
            write(stderr, *) ' Warning. Sum of hopping probabilities for all states higher than 1.'
            write(stderr, *) '   time:', t1%time
            write(stderr, *) '   cstate:', t1%cstate
            write(stderr, *) '   fprob:', t1%prob
        end if
        check = .false.
        prob = 0.0_dp
        call random_number(rnum)
        do i = 1, size(t1%prob)
            prob = prob + t1%prob(i)
            if (rnum < prob) then
                check = .true.
                t1%cstate = i
                exit
            end if
        end do
    end function check_hop


    subroutine new_step(t0, dt, new_t)
        use control_var, only : ctrl
        use nuclear_dyn_mod
        use interface_mod
        type(trajtype), intent(in) :: t0
        real(dp), intent(in) :: dt
        type(trajtype), intent(out) :: new_t

        write(stdout, '(5x,a,f10.4,a)') ' Computing step at t=', (t0%time+dt) * aut_fs ,'.'
        new_t = t0
        new_t%time = t0%time + dt
        call dyn_updategeom(dt, t0%mass, t0%geom, t0%grad, t0%velo, new_t%geom, ctrl%cns, &
        &                   ctrl%orientlvl)
        call run_qm(new_t, .false.)
        call dyn_updatevelo(dt, t0%mass, new_t%geom, t0%grad, new_t%grad, t0%velo, new_t%velo, &
        &                   ctrl%cns, ctrl%orientlvl)
    end subroutine new_step


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: bisect_gap
    !> @brief Add extra time step between three steps to locate gap minimum
    !> @description
    !! Adds an extra time step halfway between t0 and t1 or t1 and t2 (depending on which time
    !! difference is larger). Checks between which three points the gap falls and updates t0, t1
    !! and t2 to be those steps.
    !----------------------------------------------------------------------------------------------
    subroutine bisect_gap(t0, t1, t2, check)
        type(trajtype), intent(inout) :: t0 !< Trajectory variables before gap minimum
        type(trajtype), intent(inout) :: t1 !< Trajectory variables at gap minimum
        type(trajtype), intent(inout) :: t2 !< Trajectory variables after gap minimum
        logical, intent(out) :: check
        type(trajtype) :: wrk_t !< Trajectory variables at bisection
        real(dp) :: dt1, dt2

        dt1 = t1%time - t0%time
        dt2 = t2%time - t1%time
        check = .false.
        if (1.5*dt1 > dt2) then
            call new_step(t0, dt1 / 2, wrk_t)
            if (any(check_gap(t0, wrk_t, t1))) then
                t2 = t1
                t1 = wrk_t
                check = .true.
            else
                if (any(check_gap(wrk_t, t1, t2))) then
                    t0 = wrk_t
                    check = .true.
                end if
            end if
        else
            call new_step(t1, dt2 / 2, wrk_t)
            if (any(check_gap(t0, t1, wrk_t))) then
                t2 = wrk_t
                check = .true.
            else
                if (any(check_gap(t1, wrk_t, t2))) then
                    t0 = t1
                    t1 = wrk_t
                    check = .true.
                end if
            end if
        end if

        if (.not. check) then
            write(stderr, *) 'Warning. Gap minimum not found after adding an extra time step.'
            write(stderr, '(999(e24.16, 1x))') t0%time, t0%qe
            write(stderr, '(999(e24.16, 1x))') t1%time, t1%qe
            write(stderr, '(999(e24.16, 1x))') t2%time, t2%qe
            write(stderr, '(999(e24.16, 1x))') wrk_t%time, wrk_t%qe
        end if
    end subroutine bisect_gap


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: lz_prob_interval
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
    subroutine lz_prob_interval(t0, t1, t2, g0, g1, g2, gap_err, tol_err, sd, prob)
        real(dp), intent(in) :: t0 !< Time for which g0 is evaluated
        real(dp), intent(in) :: t1 !< Time at gap minimum
        real(dp), intent(in) :: t2 !< Time for which g2 is evaluated
        real(dp), intent(in) :: g0 !< Finite difference gap before minimum
        real(dp), intent(in) :: g1 !< Gap at gap minimum
        real(dp), intent(in) :: g2 !< Finite difference gap after minimum
        real(dp), intent(in) :: gap_err !< Assumed error in gap minimum due to finite step
        real(dp), intent(in) :: tol_err !< Assumed error in all energies due to precision
        real(dp), intent(out) :: sd(3) !< Calculated numerical second derivs (same indexing as prob)
        real(dp), intent(out) :: prob(3) !< Probabilities calculated assuming (1) errors aligning to
            !! give min probability, (2) no errors and (3) errors aligning for max probability
        real(dp) :: err_gap, err_g0, err_g1, err_g2
        real(dp), parameter :: tinydp = 1.0e-15_dp

        ! Minimum probability
        err_gap = abs(g1) + tol_err
        err_g0 = abs(g0) - tol_err
        err_g1 = abs(g1) + tol_err
        err_g2 = abs(g2) - tol_err
        if ((err_g0 < err_g1) .or. (err_g2 < err_g1)) then
            prob(1) = 0.0_dp
            sd(1) = 0.0_dp
        else
            sd(1) = second_derivative(t0, t1, t2, err_g0, err_g1, err_g2)
            prob(1) = lz_prob(err_gap, abs(sd(1)))
        end if

        ! No error probability
        sd(2) = second_derivative(t0, t1, t2, g0, g1, g2)
        prob(2) = lz_prob(abs(g1), abs(sd(2)))

        ! Maximum probability
        err_gap = max(abs(g1) - gap_err - tol_err, tinydp)
        err_g0 = abs(g0) + tol_err
        err_g1 = max(abs(g1) - tol_err, tinydp)
        err_g2 = abs(g2) + tol_err
        sd(3) = second_derivative(t0, t1, t2, err_g0, err_g1, err_g2)
        prob(3) = lz_prob(err_gap, abs(sd(3)))
    end subroutine lz_prob_interval

    !----------------------------------------------------------------------------------------------
    ! FUNCTION: lz_prob_interval_gap_only
    !
    !> @brief Expected error interval for using the LZ formula due to accuracy of gap minimum.
    !> @details
    !! Additionaly, the true gap minimum is not the lowest calculated point so an error 
    !! proportional to the change in the gap between the time steps is assumed.
    !----------------------------------------------------------------------------------------------
    subroutine lz_prob_interval_gap_only(g1, sd, gap_err, tol_err, prob)
        real(dp), intent(in) :: g1 !< Gap at gap minimum
        real(dp), intent(in) :: sd !< Second derivative at gap minimum
        real(dp), intent(in) :: gap_err !< Assumed error in gap minimum due to finite step
        real(dp), intent(in) :: tol_err !< Assumed error in all energies due to precision
        real(dp), intent(out) :: prob(3) !< Probabilities calculated assuming (1) errors aligning to
            !! give min probability, (2) no errors and (3) errors aligning for max probability
        real(dp) :: err_gap
        real(dp), parameter :: tinydp = 1.0e-15_dp

        ! Minimum probability
        err_gap = abs(g1) + tol_err + gap_err
        prob(1) = lz_prob(err_gap, abs(sd))

        ! No error probability
        prob(2) = lz_prob(abs(g1), abs(sd))

        ! Maximum probability
        err_gap = max(abs(g1) - gap_err - tol_err, tinydp)
        prob(3) = lz_prob(err_gap, abs(sd))
    end subroutine lz_prob_interval_gap_only


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: lz_prob
    !> @brief Return the Landau-Zener hopping probability
    !----------------------------------------------------------------------------------------------
    function lz_prob(gap, sd) result(prob)
        real(dp), intent(in) :: gap !< Gap at gap minimum
        real(dp), intent(in) :: sd !< Second derivative at gap minimum
        real(dp) :: prob
        prob = exp(- pio2 * sqrt( gap**3 / sd ))
    end function lz_prob


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: second_derivative
    !> @brief Calculate numerical second derivative for unevenly spaced points.
    !----------------------------------------------------------------------------------------------
    function second_derivative(x1, x2, x3, y1, y2, y3) result(sd)
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
    end function second_derivative


end module sh_lz_mod
