!--------------------------------------------------------------------------------------------------
! MODULE: sh_fssh_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2016
!> @author Cristina Sanz, Autonoma University Madrid
!> @date May, 2024: sh_sosh and sh_interpolate_sovec added for the Fewest Switches Surface 
!! Hopping method using spin-orbit couplings
!
! DESCRIPTION:
!> @brief Surface hopping algorithm.
!--------------------------------------------------------------------------------------------------
module sh_fssh_mod
    use global_defs
    implicit none

    private
    public :: sh_adiabatic!, sh_sosh
    

contains

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: sh_fssh
    !
    ! DESCRIPTION:
    !> @brief Tully's Fewest Switches Surface Hopping Method.
    !> @details
    !! Propagates electronic wave function coefficients and determines hops for the FSSH method.
    !! The coefficients are propagated using an external subroutine for solving the system of
    !! ordinary differential equations.
    !! The nuclear time step (dt) is split into smaller time steps for which the coefficients and
    !! hopping probabilities are calculated.
    !----------------------------------------------------------------------------------------------
    subroutine sh_adiabatic(opt_clvl, interpolation_en, interpolation_tdc, t0, t1, t2, wf_t1, &
    &                       wf_t2, nstep, vel_t1, vel_t2, rng)
        use ode_call_mod ! Interface to Shampine/Gordon ODE solver.
        use mqc_wave_function_mod, only : mqc_wave_function
        use random_mod, only : rng_type
        use matrix_mod, only : diagonal_mat
        use tdc_mod, only : hst_tdc, nadvec2tdc, npi_tdc_integrated
        character(len=*), intent(in) :: opt_clvl !< Method for calculating time-derivative couplings.
        integer, intent(in) :: interpolation_en !< Method for interpolating energies during the time step.
        integer, intent(in) :: interpolation_tdc !< Method for interpolating overlaps during the time step.
        real(dp), intent(in) :: t0 !< Time at previous step.
        real(dp), intent(in) :: t1 !< Time at current step.
        real(dp), intent(in) :: t2 !< Time at next step.
        type(mqc_wave_function), intent(in) :: wf_t1 !< Wave function at t1.
        type(mqc_wave_function), intent(inout) :: wf_t2 !< Wave function at t2.
        integer, intent(in) :: nstep !< Number of substeps.
        real(dp), intent(in) :: vel_t1(:, :) !< Velocities at t1.
        real(dp), intent(in) :: vel_t2(:, :) !< Velocities at t2.
        class(rng_type), allocatable, intent(inout) :: rng !< Random number generator to use.
        
        real(dp) :: edt !< Time step for the propagation of the electronic WF.
        real(dp) :: tt !< Current time during propagation.
        real(dp) :: prob !< Probability of hopping into a state.
        real(dp) :: cprob !< Cumulative probability of hopping into any state.

        integer :: cstate !< Current state index.
        integer :: i
        integer :: st
        integer :: de_flag
        real(dp) :: rnum !< Random number for surface hopping.
        real(dp), allocatable :: nadv1(:, :, :) !< Nonadiabatic coupling vectors at t1.
        real(dp), allocatable :: nadv2(:, :, :) !< Nonadiabatic coupling vectors at t2.
        real(dp), allocatable :: en_t(:) !< Energies at current time.
        real(dp), allocatable :: tdc_t1(:, :) !< TDC matrix at t1.
        real(dp), allocatable :: tdc_t2(:, :) !< TDC matrix at t2.
        real(dp), allocatable :: tdc_t(:, :) !< TDC matrix at current time.
        real(dp), allocatable :: wrk_1(:, :) !< Work array for storing TDCs at half step.
        real(dp), allocatable :: wrk_2(:, :) !< Work array for storing TDCs at half step.
        complex(dp), allocatable :: soc_t1(:, :) !< Spin-orbit coupling matrix at t1.
        complex(dp), allocatable :: soc_t2(:, :) !< Spin-orbit coupling matrix at t2.
        complex(dp), allocatable :: soc_t(:, :) !< Spin-orbit coupling matrix at current time.
        integer, parameter :: im_i = cmplx(0.0_dp, 1.0_dp, kind=dp)

        odens = size(wf_t2%coeff)
        allocate(odecmat(odens, odens))

        ! Propagation time step.
        tt = t1
        edt = (t2 - t1) / nstep
        cstate = wf_t1%active_state

        select case(opt_clvl)
        case("hst")
            call hst_tdc(t1-t0, wf_t1%overlap(1)%c, wrk_1)
            call hst_tdc(t2-t1, wf_t2%overlap(1)%c, wrk_2)
            call sh_interpolate_tdc(interpolation_tdc, 0.5_dp*(t0 + t1), 0.5_dp*(t1 + t2), &
            &                       t1, wrk_1, wrk_2, tdc_t1)
            call sh_interpolate_tdc(interpolation_tdc, 0.5_dp*(t0 + t1), 0.5_dp*(t1 + t2), &
            &                       t2, wrk_1, wrk_2, tdc_t2)
        case("npi")
            call npi_tdc_integrated(t1-t0, wf_t1%overlap(1)%c, wrk_1)
            call npi_tdc_integrated(t2-t1, wf_t2%overlap(1)%c, wrk_2)
            call sh_interpolate_tdc(interpolation_tdc, 0.5_dp*(t0 + t1), 0.5_dp*(t1 + t2), &
            &                       t1, wrk_1, wrk_2, tdc_t1)
            call sh_interpolate_tdc(interpolation_tdc, 0.5_dp*(t0 + t1), 0.5_dp*(t1 + t2), &
            &                       t2, wrk_1, wrk_2, tdc_t2)
        case("nadvec")
            nadv1 = build_nadvec_matrix(wf_t1%need_nadv, wf_t1%qm_state)
            nadv2 = build_nadvec_matrix(wf_t2%need_nadv, wf_t2%qm_state)
            call nadvec2tdc(nadv1, vel_t1, tdc_t1)
            call nadvec2tdc(nadv2, vel_t2, tdc_t2)
        end select

        if (any(wf_t2%need_soc)) then
            soc_t1 = build_soc_matrix(wf_t1%need_soc, wf_t1%qm_state)
            soc_t2 = build_soc_matrix(wf_t2%need_soc, wf_t2%qm_state)
        end if

        wf_t2%prob = 0.0_dp

        do i = 1, nstep
            ! Get energies and TDCs for current substep.
            call sh_interpolate_energy(interpolation_en, t1, t2, tt, wf_t1%qm_state(:)%energy, wf_t2%qm_state(:)%energy, en_t)
            call sh_interpolate_tdc(interpolation_tdc, t1, t2, tt, tdc_t1, tdc_t2, tdc_t)
            odecmat = cmplx(0.0_dp, -diagonal_mat(en_t), kind=dp) - tdc_t
            if (any(wf_t2%need_soc)) then
                call sh_interpolate_soc(interpolation_tdc, t1, t2, tt, soc_t1, soc_t2, soc_t)
                odecmat = odecmat - im_i * soc_t
            end if
            
            ! Propagate wf coefficients.
            call callode(odens, wf_t2%coeff, tt, edt, de_flag)
            tt = tt + edt

            ! Determine hopping probabilities.
            call rng%uniform(rnum)
            cprob = 0.0_dp
            hop: do st = 1, odens
                if (st == cstate) cycle
                prob = - 2 * edt * odecmat(st, cstate) * real(conjg(wf_t2%coeff(st)) * wf_t2%coeff(cstate)) / &
                     & (abs(wf_t2%coeff(cstate))**2)
                if (prob > 0.0_dp) then ! Not actual probability, can be negative.
                    cprob = cprob + prob
                    wf_t2%prob(st) = wf_t2%prob(st) + prob
                    if (rnum < cprob) then
                        cstate = st
                        exit hop
                    end if
                end if
            end do hop
            
        end do

        wf_t2%active_state = cstate

        deallocate(odecmat)
    end subroutine sh_adiabatic

    function build_nadvec_matrix(need_nadv, states) result(nadvec)
        use state_mod, only : state
        logical, intent(in) :: need_nadv(:, :)
        type(state), intent(in) :: states(:)
        real(dp), allocatable :: nadvec(:, :, :)
        integer :: i, j, n_dof

        if (.not. any(need_nadv)) then
            call errstop("sh_fssh_mod", "No nonadiabatic coupling vectors requested.", 1)
        end if
        do i = 1, size(states)
            if (.not. allocated(states(i)%nadv)) then
                call errstop("sh_fssh_mod", "Nonadiabatic coupling vectors not allocated.", 1)
            end if
            do j = 1, size(states)
                if (need_nadv(i, j)) then
                    if (.not. allocated(states(i)%nadv(j)%c)) then
                        call errstop("sh_fssh_mod", "Nonadiabatic coupling vector requested but not allocated.", 1)
                    end if
                    n_dof = size(states(i)%nadv(j)%c)
                end if
            end do
        end do

        allocate(nadvec(n_dof, size(states), size(states)), source=0.0_dp)
        do i = 1, size(states)
            do j = 1, size(states)
                if (need_nadv(i, j)) nadvec(:, i, j) = states(i)%nadv(j)%c
            end do
        end do
    end function build_nadvec_matrix


    function build_soc_matrix(need_soc, states) result(soc_mat)
        use state_mod, only : state
        logical, intent(in) :: need_soc(:, :)
        type(state), intent(in) :: states(:)
        complex(dp), allocatable :: soc_mat(:, :)
        integer :: i, j, n_states

        if (.not. any(need_soc)) then
            call errstop("sh_fssh_mod", "No spin-orbit couplings requested.", 1)
        end if
        n_states = size(states)
        do i = 1, n_states
            if (any(need_soc(i, :))) then
                if (.not. allocated(states(i)%soc)) then
                    call errstop("sh_fssh_mod", "Spin-orbit couplings not allocated.", 1)
                end if
            end if
        end do

        allocate(soc_mat(n_states, n_states), source=(0.0_dp, 0.0_dp))
        do i = 1, n_states
            do j = 1, n_states
                if (need_soc(i, j)) then
                    soc_mat(i, j) = states(i)%soc(j)
                end if
            end do
        end do
    end function build_soc_matrix

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: SH_Interpolate_Energy
    !
    ! DESCRIPTION:
    !> @brief Energy interpolation during single nuclear time step in surface hopping.
    !> @details
    !! Returns electronic state energies during the propagation of wave function coefficients in a
    !! surface hopping calculation. The interpolation is done based on the calculated energies at
    !! the beginning and end of the time nuclear time step, E(t) and E(t + dt).
    !!
    !! The behaviour of the subroutine is determined by the value of opt:
    !! - 0 - At all substeps, the energy is constant and equal to the E(t + dt).
    !! - 1 - The energies are equal to E(t) until t + dt/2, and to E(t + dt) afterwards.
    !! - 2 - The energies are a linear interpolation between E(t) and E(t + dt).
    !----------------------------------------------------------------------------------------------
    subroutine sh_interpolate_energy(opt, t1, t2, tt, en1, en2, ien)
        integer, intent(in) :: opt !< Type of interpolation.
        real(dp), intent(in) :: t1 !< Time at which en1 is evaluated.
        real(dp), intent(in) :: t2 !< Time at which en2 is evaluated.
        real(dp), intent(in) :: tt !< Current time during propagation.
        real(dp), intent(in) :: en1(:) !< Energies in previous nuclear time step.
        real(dp), intent(in) :: en2(:) !< Energies in current nuclear time step.
        real(dp), allocatable, intent(out) :: ien(:) !< Interpolated energies.

        select case(opt)
        ! Constant value.
        case(0)
            ien = en2
        ! Step function in middle of time step.
        case(1)
            if (tt < (t1 + t2) / 2) then
                ien = en1
            else
                ien = en2
            end if
        ! Linear interpolation.
        case(2)
            ien = en1 + (en2 - en1) * (tt - t1) / (t2 - t1)
        end select
    end subroutine sh_interpolate_energy


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: sh_interpolate_tdc
    !
    ! DESCRIPTION:
    !> @brief TDC interpolation during single nuclear time step in surface hopping.
    !> @details
    !! Returns time-derivative couplings during the propagation of wave function coefficients in a
    !! surface hopping calculation.
    !----------------------------------------------------------------------------------------------
    subroutine sh_interpolate_tdc(opt, t1, t2, tt, tdc1, tdc2, itdc)
        integer, intent(in) :: opt !< Type of interpolation.
        real(dp), intent(in) :: t1 !< Time at which en1 is evaluated.
        real(dp), intent(in) :: t2 !< Time at which en2 is evaluated.
        real(dp), intent(in) :: tt !< Current time during propagation.
        real(dp), intent(in) :: tdc1(:, :) !< Overlaps between wfs at t0 - dt and t0.
        real(dp), intent(in) :: tdc2(:, :) !< Overlaps between wfs at t0 and t0 + dt.
        real(dp), allocatable, intent(inout) :: itdc(:, :) !< Time-derivative couplings.
        
        select case(opt)
        ! Constant value.
        case(0)
            itdc = tdc2
        ! Step function in middle of time step.
        case(1)
            if (tt < (t1 + t2) / 2) then
                itdc = tdc1
            else
                itdc = tdc2
            end if
        ! Linear interpolation.
        case(2)
            itdc = tdc1 + (tdc2 - tdc1) * (tt - t1) / (t2 - t1)
        end select
    end subroutine sh_interpolate_tdc


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: sh_interpolate_soc
    !
    ! DESCRIPTION:
    !> @brief SOC interpolation during single nuclear time step in surface hopping.
    !> @details
    !! Returns spin-orbit couplings during the propagation of wave function coefficients in a
    !! surface hopping calculation.
    !----------------------------------------------------------------------------------------------
    subroutine sh_interpolate_soc(opt, t1, t2, tt, soc1, soc2, isoc)
        integer, intent(in) :: opt !< Type of interpolation.
        real(dp), intent(in) :: t1 !< Time at which en1 is evaluated.
        real(dp), intent(in) :: t2 !< Time at which en2 is evaluated.
        real(dp), intent(in) :: tt !< Current time during propagation.
        complex(dp), intent(in) :: soc1(:, :) !< Overlaps between wfs at t0 - dt and t0.
        complex(dp), intent(in) :: soc2(:, :) !< Overlaps between wfs at t0 and t0 + dt.
        complex(dp), allocatable, intent(inout) :: isoc(:, :) !< Spin-orbit couplings.

        select case(opt)
        ! Constant value.
        case(0)
            isoc = soc2
        ! Step function in middle of time step.
        case(1)
            if (tt < (t1 + t2) / 2) then
                isoc = soc1
            else
                isoc = soc2
            end if
        ! Linear interpolation.
        case(2)
            isoc = soc1 + (soc2 - soc1) * (tt - t1) / (t2 - t1)
        end select
    end subroutine sh_interpolate_soc
 
 
end module sh_fssh_mod
