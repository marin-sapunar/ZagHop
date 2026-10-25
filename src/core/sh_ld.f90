!--------------------------------------------------------------------------------------------------
! MODULE: sh_ld_mod
!
! DESCRIPTION:
!> @brief Fewest switches surface hopping in the locally diabatic representation (LD-FSSH).
!> @details
!! The electronic wave function is propagated in the locally diabatic basis. The active state
!! and the coefficients refer to the adiabatic basis.
!--------------------------------------------------------------------------------------------------
module sh_ld_mod
    use global_defs
    use tsh_method_mod, only : tsh_method, base_initialize => initialize
    use evaluator_base_mod, only : potential_evaluator
    use orthog_mod, only : orthog_lowdin
    use matrix_mod, only : mat_sy_exp
    use random_mod, only : rng_type
    implicit none

    private
    public :: tsh_ld


    type, extends(tsh_method) :: tsh_ld
        character(len=:), allocatable :: representation
        character(len=:), allocatable :: prob_formula !< Formula for the hopping probabilities:
        !! 'sharc', 'ld01', 'ld01_l' or 'ld19'.
        integer :: n_substep = 1 !< Substeps for propagating the coefficients during a step.
        integer :: decohlvl = 1 !< Decoherence correction:
        !! - 0 - None.
        !! - 1 - Energy based decoherence scheme (Granucci, Persico; doi: 10.1063/1.2715585).
        complex(dp), allocatable :: coeff(:)
        real(dp), allocatable :: u_t1(:, :)
        real(dp), allocatable :: h_t1(:, :)
    contains
        procedure :: initialize => ld_initialize
        procedure :: start => ld_start
        procedure :: propagate => ld_propagate
        procedure :: decohere => ld_decohere
        procedure :: write_coefficients => ld_write_coefficients
    end type tsh_ld


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ld_initialize
    !> @brief Set default options and put the full wave function in the initial state.
    !----------------------------------------------------------------------------------------------
    subroutine ld_initialize(self, rng, potential, active_state)
        class(tsh_ld), intent(inout) :: self
        class(rng_type), pointer, intent(in) :: rng
        class(potential_evaluator), intent(in) :: potential
        integer, intent(in) :: active_state

        call base_initialize(self, rng, potential, active_state)
        self%representation = 'locally_diabatic'
        self%prob_formula = 'sharc'
        allocate(self%coeff(self%n_state), source=cmplx(0.0_dp, 0.0_dp, kind=dp))
        self%coeff(active_state) = cmplx(1.0_dp, 0.0_dp, kind=dp)
    end subroutine ld_initialize


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ld_start
    !> @brief Store the adiabatic Hamiltonian and states at the initial point.
    !----------------------------------------------------------------------------------------------
    subroutine ld_start(self, t, potential)
        class(tsh_ld), intent(inout) :: self
        real(dp), intent(in) :: t
        class(potential_evaluator), intent(in) :: potential

        self%t = t
        self%h_t1 = potential%get_Hamiltonian('adiabatic')
        self%u_t1 = potential%get_transformation('diabatic', 'adiabatic')
    end subroutine ld_start


    subroutine ld_propagate(self, t2, potential)
        use matrix_mod, only : unitary_transform, unit_mat
        class(tsh_ld), intent(inout) :: self
        class(potential_evaluator), intent(in) :: potential
        real(dp), intent(in) :: t2 !< Time after propagation
        real(dp) :: dt
        real(dp) :: t(self%n_state, self%n_state) !< Orthogonalized overlap matrix.
        real(dp) :: h_t(self%n_state, self%n_state) !< Hamiltonian during substep.
        real(dp) :: h_t2(self%n_state, self%n_state) !< Hamiltonian at time t2.
        real(dp) :: dh(self%n_state, self%n_state) !< Change of Hamiltonian over the step (LD basis).
        complex(dp) :: w3(self%n_state, self%n_state) !< Work array 3.
        complex(dp) :: w4(self%n_state, self%n_state) !< Work array 4.
        real(dp), allocatable :: u_t2(:, :)
        integer :: i, cstate
        complex(dp) :: c_propagation_t0(self%n_state) !< Coefficients at start of substep in propagation representation.
        complex(dp) :: c_propagation_dt(self%n_state) !< Coefficients at end of substep in propagation representation.
        complex(dp) :: c_hopping_t0(self%n_state) !< Coefficients at start of substep in hopping representation.
        complex(dp) :: c_hopping_dt(self%n_state) !< Coefficients at end of substep in hopping representation.

        cstate = self%active_state
        c_propagation_t0 = self%coeff
        c_hopping_t0 = self%coeff
        self%prob = 0.0_dp
        u_t2 = potential%get_transformation('diabatic', 'adiabatic')
        h_t2 = potential%get_Hamiltonian('adiabatic')

        ! Get overlap matrix.
        t = matmul(transpose(self%u_t1), u_t2) ! T = U(t)^t . U(t+dt)
        call orthog_lowdin(t)

        dh = matmul(matmul(t, h_t2), transpose(t)) ! H(t+dt) = T.H(t+dt).T^t
        dh = dh - self%h_t1

        dt = (t2 - self%t) / self%n_substep

        do i = 1, self%n_substep
            h_t = self%h_t1 + dh * real(i - 0.5_dp, dp) / self%n_substep
            w3 = mat_sy_exp(h_t, cmplx(0.0_dp, -dt, kind=dp))
            c_propagation_dt = matmul(w3, c_propagation_t0)
            c_hopping_dt = matmul(transpose(t), c_propagation_dt)
            select case(self%prob_formula)
            case('sharc')
                w4 = matmul(transpose(t), w3)
                self%prob = self%prob + prob_sharc(cstate, c_hopping_t0, c_hopping_dt, w4)
            case('ld01')
                w4 = matmul(transpose(t), w3)
                self%prob = self%prob + prob_ld01(cstate, c_hopping_t0, c_hopping_dt, w4)
            case('ld01_l')
                self%prob = self%prob + prob_ld01_l(cstate, c_hopping_t0, c_hopping_dt, t)
            case('ld19')
                self%prob = self%prob + prob_ld19(cstate, c_hopping_t0, c_hopping_dt, t)
            end select
            c_propagation_t0 = c_propagation_dt
            c_hopping_t0 = c_hopping_dt
        end do

        call self%select_state()

        self%coeff = c_hopping_dt
        self%u_t1 = u_t2
        self%h_t1 = h_t2
        self%t = t2
    end subroutine ld_propagate


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ld_decohere
    !> @brief Apply the decoherence correction selected by decohlvl to the coefficients.
    !----------------------------------------------------------------------------------------------
    subroutine ld_decohere(self, ekin, dt, potential)
        use decoherence_mod, only : edc
        class(tsh_ld), intent(inout) :: self
        real(dp), intent(in) :: ekin !< Kinetic energy.
        real(dp), intent(in) :: dt !< Time step.
        class(potential_evaluator), intent(in) :: potential

        select case(self%decohlvl)
        case(0) ! No decoherence correction.
        case(1) ! Energy based decoherence.
            if (stdp3) write(stdout, *) ' Dechoerence correction type EDC.'
            call edc(self%active_state, ekin, potential%get_energy('adiabatic'), dt, self%coeff)
        end select
    end subroutine ld_decohere


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ld_write_coefficients
    !> @brief Write real parts followed by imaginary parts of the coefficients on one line.
    !----------------------------------------------------------------------------------------------
    subroutine ld_write_coefficients(self, unit)
        class(tsh_ld), intent(in) :: self
        integer, intent(in) :: unit

        write(unit, '(1x,1000e22.12)') real(self%coeff), aimag(self%coeff)
    end subroutine ld_write_coefficients


    pure function prob_sharc(b, c_t, c_dt, p) result(prob)
        integer, intent(in) :: b
        complex(dp), intent(in) :: c_t(:)
        complex(dp), intent(in) :: c_dt(:)
        complex(dp), intent(in) :: p(:, :)
        real(dp) :: prob(size(c_t))
        real(dp) :: b_diff
        integer :: a
        real(dp), parameter :: tiny = 1.0e-10_dp

        prob = 0.0_dp
        b_diff = 1 - abs(c_dt(b))**2 / abs(c_t(b))**2
        if (abs(b_diff) <= tiny) return
        do a = 1, size(c_t)
            if (a == b) cycle
            prob(a) = max(b_diff * real(c_dt(a) * conjg(p(a, b)) * conjg(c_t(b))) / &
            &         (abs(c_t(b))**2 - real(c_dt(b) * conjg(p(b, b)) * conjg(c_t(b)))), 0.0_dp)
        end do
    end function prob_sharc


    pure function prob_ld01(m, c_t, c_dt, ubar) result(prob)
        integer, intent(in) :: m
        complex(dp), intent(in) :: c_t(:)
        complex(dp), intent(in) :: c_dt(:)
        complex(dp), intent(in) :: ubar(:, :)
        real(dp) :: prob(size(c_t))
        integer :: n

        prob = 0.0_dp
        do n = 1, size(c_t)
            if (m == n) cycle
            prob(n) = -real(ubar(m, n) * c_t(n) * conjg(c_dt(m))) / abs(c_t(m))**2 * &
                      (abs(c_dt(m))**2 - abs(c_t(m))**2) / &
                      (abs(c_dt(m))**2 - real(ubar(m, m) * c_t(m) * conjg(c_dt(m))))
        end do
    end function prob_ld01


    pure function prob_ld01_l(m, c_t, c_dt, u) result(prob)
        integer, intent(in) :: m
        complex(dp), intent(in) :: c_t(:)
        complex(dp), intent(in) :: c_dt(:)
        real(dp), intent(in) :: u(:, :)
        real(dp) :: prob(size(c_t))
        integer :: n

        prob = 0.0_dp
        do n = 1, size(c_t)
            if (m == n) cycle
            prob(n) = max(2.0_dp * real(u(m, n) * c_t(n) * conjg(c_dt(m))) / abs(c_t(m))**2, 0.0_dp)
        end do
    end function prob_ld01_l


    pure function prob_ld19(m, c_t, c_dt, u) result(prob)
        integer, intent(in) :: m
        complex(dp), intent(in) :: c_t(:)
        complex(dp), intent(in) :: c_dt(:)
        real(dp), intent(in) :: u(:, :)
        real(dp) :: prob(size(c_t))
        integer :: n
        real(dp) :: w_m

        prob = 0.0_dp
        w_m = (1 - abs(c_dt(m))**2 / abs(c_t(m))**2)
        if (w_m <= 0.0_dp) return

        do n = 1, size(c_t)
            if (m == n) cycle
            ! Numerator of x_mn
            prob(n) = sqrt(max(abs(c_dt(n))**2 - abs(c_t(n))**2, 0.0_dp)) * abs(u(m, n))
        end do

        prob = w_m * prob / sum(prob)
    end function prob_ld19


end module sh_ld_mod
