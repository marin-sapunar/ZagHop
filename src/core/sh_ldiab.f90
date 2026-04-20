!--------------------------------------------------------------------------------------------------
! MODULE: sh_ldiab_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2016
!
! DESCRIPTION:
!> @brief Surface hopping algorithm.
!--------------------------------------------------------------------------------------------------
module sh_ldiab_mod
    use global_defs
    implicit none

    private
    public :: sh_diabatic


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: SH_Diabatic
    !
    ! DESCRIPTION:
    !> @brief Fewest switches surface hopping in a locally diabatic basis.
    !> @details
    !! Propagates electronic wave function coefficients and determines hops for the SH method.
    !----------------------------------------------------------------------------------------------
    subroutine sh_diabatic(t1, t2, wf_t1, wf_t2, rng)
        use matrix_mod, only : diag, &
                               mat_sy_exp
        use orthog_mod, only : orthog_lowdin
        use linalg_wrapper_mod, only : gemm, gemv
        use random_mod, only : rng_type
        use mqc_wave_function_mod, only : mqc_wave_function
        real(dp), intent(in) :: t1 !< Initial time.
        real(dp), intent(in) :: t2 !< Final time.
        type(mqc_wave_function), intent(in) :: wf_t1 !< Wave function at time t1.
        type(mqc_wave_function), intent(inout) :: wf_t2 !< Wave function at time t2.
        class(rng_type), intent(inout) :: rng
        real(dp) :: t(wf_t1%n_state, wf_t1%n_state) !< Orthogonalized overlap matrix.
        real(dp) :: w1(wf_t1%n_state, wf_t1%n_state) !< Work array 1.
        real(dp) :: w2(wf_t1%n_state, wf_t1%n_state) !< Work array 2.
        complex(dp) :: u(wf_t1%n_state, wf_t1%n_state) !< Transformation matrix.
        complex(dp) :: w3(wf_t1%n_state, wf_t1%n_state) !< Work array 3.
        complex(dp) :: w4(wf_t1%n_state, wf_t1%n_state) !< Work array 4.
        complex(dp) :: pwf(wf_t1%n_state) !< Previous WF coefficients.
        complex(dp) :: cwf(wf_t1%n_state) !< Current WF coefficients.
        real(dp) :: b(wf_t1%n_state) !< Contributions of each state to change in population.
        real(dp) :: rnum !< Random number for surface hopping.
        real(dp) :: denom !< Denominator in calculation of b.
        real(dp), parameter :: thresh = 1.e-10_dp !< Threshold for denominator.
        real(dp) :: cprob !< Cumulative probability of hopping into any state.
        integer :: k, l

        t = wf_t2%overlap
        pwf = wf_t1%coeff

        ! Orthogonalize overlap matrix.
        call orthog_lowdin(t)

        ! Generate Z matrix. (Approx. Hamiltonian at half step.)
        call gemm(t, diag(wf_t2%en), w1)
        call gemm(w1, t, w2, transb='T') ! H(t+dt/2) = T.E(t+dt).Tt
        w2 = (diag(wf_t1%en) + w2) * 0.5_dp ! Z = (E(0) + H(t+dt))/2

        ! Generate U matrix. (Transformation matrix.)
        w3 = mat_sy_exp(w2, cmplx(0.0_dp, t1-t2, kind=dp)) ! exp(-i * Z * dt)
        w4 = cmplx(t, 0.0_dp, dp) !< Convert T to complex matrix.
        call gemm(w4, w3, u, transa='T') ! u = T^t exp(-i * Z * dt)

        ! Calculate new wf coefficients.
        cwf = matmul(u, pwf) ! A(t0 + dt) = U A(t0)

        ! Calculate hopping probabilities.
        k = wf_t2%active_state
        b = 0.0_dp
        denom = abs(cwf(k))**2 - real(u(k, k) * pwf(k) * conjg(cwf(k)))
        if (denom > thresh) then
            denom = (abs(cwf(k))**2 - abs(pwf(k))**2) / denom
            do l = 1, wf_t1%n_state
                if (l == k) cycle
                b(l) = - real(u(k, l) * pwf(l) * conjg(cwf(k))) * denom
            end do
        else
            ! If the denominator is close to zero, the property b(k, l) = - b(l, k) is used.
            do l = 1, wf_t1%n_state
                if (l == k) cycle
                denom = abs(cwf(l))**2 - real(u(l, l) * pwf(l) * conjg(cwf(l)))
                if (denom < thresh) cycle
                denom = (abs(cwf(l))**2 - abs(pwf(l))**2) / denom
                b(l) = real(u(l, k) * pwf(k) * conjg(cwf(l))) * denom
            end do
        end if
        b = b / abs(pwf(k))**2

        ! Determine if hop should occur.
        call rng%uniform(rnum)
        cprob = 0.0_dp
        wf_t2%prob = 0.0_dp
        hop: do l = 1, wf_t1%n_state
            if (l == k) cycle
            if (b(l) > 0.0_dp) then ! Not actual probability, can be negative.
                cprob = cprob + b(l)
                wf_t2%prob(l) = wf_t2%prob(l) + b(l)
                if (rnum < cprob) then
                    wf_t2%active_state = l
                    exit hop
                end if
            end if
        end do hop

        wf_t2%coeff = cwf
    end subroutine sh_diabatic


end module sh_ldiab_mod
