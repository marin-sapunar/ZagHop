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
    subroutine sh_diabatic(n, dt, qe1, qe2, cwf, tst, olp, fprob)
        use matrix_mod, only : diagonal_mat, &
                               mat_sy_exp
        use orthog_mod, only : orthog_lowdin
        use linalg_wrapper_mod, only : gemm, gemv
        integer, intent(in) :: n !< Number of states.
        real(dp), intent(in) :: dt !< Nuclear dynamics time step.
        real(dp), intent(in) :: olp(n, n) !< Overlap matrix for current step.
        real(dp), intent(in) :: qe1(n) !< Energies at t.
        real(dp), intent(in) :: qe2(n) !< Energies at t + dt.
        complex(dp), intent(inout) :: cwf(n) !< WF coefficients in the adiabatic basis.
        integer, intent(inout) :: tst !< Current state.
        real(dp), intent(out) :: fprob(n) !< Final probability for each state.
        real(dp) :: t(n, n) !< Orthogonalized overlap matrix.
        real(dp) :: w1(n, n) !< Work array 1.
        real(dp) :: w2(n, n) !< Work array 2.
        complex(dp) :: u(n, n) !< Transformation matrix.
        complex(dp) :: w3(n, n) !< Work array 3.
        complex(dp) :: w4(n, n) !< Work array 4.
        complex(dp) :: pwf(n) !< Previous WF coefficients.
        real(dp) :: b(n) !< Contributions of each state to change in population.
        real(dp) :: rnum !< Random number for surface hopping.
        real(dp) :: denom !< Denominator in calculation of b.
        real(dp), parameter :: thresh = 1.e-10_dp !< Threshold for denominator.
        real(dp) :: cprob !< Cumulative probability of hopping into any state.
        integer :: k, l

        t = olp
        pwf = cwf

        ! Orthogonalize overlap matrix.
        call orthog_lowdin(t)

        ! Generate Z matrix. (Approx. Hamiltonian at half step.)
        call gemm(t, diagonal_mat(qe2), w1)
        call gemm(w1, t, w2, transb='T') ! H(t+dt/2) = T.E(t+dt).Tt
        w2 = (diagonal_mat(qe1) + w2) * 0.5_dp ! Z = (E(0) + H(t+dt))/2

        ! Generate U matrix. (Transformation matrix.)
        w3 = mat_sy_exp(w2, cmplx(0.0_dp, -dt, kind=dp)) ! exp(-i * Z * dt)
        w4 = cmplx(t, 0.0_dp, dp) !< Convert T to complex matrix.
        call gemm(w4, w3, u, transa='T') ! u = T^t exp(-i * Z * dt)

        ! Calculate new wf coefficients.
        cwf = matmul(u, pwf) ! A(t0 + dt) = U A(t0)

        ! Calculate hopping probabilities.
        k = tst
        b = 0.0_dp
        denom = abs(cwf(k))**2 - real(u(k, k) * pwf(k) * conjg(cwf(k)))
        if (denom > thresh) then
            denom = (abs(cwf(k))**2 - abs(pwf(k))**2) / denom
            do l = 1, n
                if (l == k) cycle
                b(l) = - real(u(k, l) * pwf(l) * conjg(cwf(k))) * denom
            end do
        else
            ! If the denominator is close to zero, the property b(k, l) = - b(l, k) is used.
            do l = 1, n
                if (l == k) cycle
                denom = abs(cwf(l))**2 - real(u(l, l) * pwf(l) * conjg(cwf(l)))
                if (denom < thresh) cycle
                denom = (abs(cwf(l))**2 - abs(pwf(l))**2) / denom
                b(l) = real(u(l, k) * pwf(k) * conjg(cwf(l))) * denom
            end do
        end if
        b = b / abs(pwf(k))**2

        ! Determine if hop should occur.
        call random_number(rnum)
        cprob = 0.0_dp
        fprob = 0.0_dp
        hop: do l = 1, n
            if (l == k) cycle
            if (b(l) > 0.0_dp) then ! Not actual probability, can be negative.
                cprob = cprob + b(l)
                fprob(l) = fprob(l) + b(l)
                if (rnum < cprob) then
                    tst = l
                    exit hop
                end if
            end if
        end do hop
    end subroutine sh_diabatic


end module sh_ldiab_mod
