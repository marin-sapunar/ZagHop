module tdse_propagator_mod
    use global_defs, only : dp
    use evaluator_base_mod, only : potential_evaluator
    use orthog_mod, only : orthog_lowdin
    use matrix_mod, only : diag, mat_sy_exp
    use random_mod, only : rng_type
    use linalg_wrapper_mod, only : gemm
    implicit none

    type :: tdse_propagator
        class(rng_type), pointer :: rng
        character(len=:), allocatable :: representation
        integer :: active_state = 0
        integer :: n_state = 0
        integer :: n_substep = 1
        real(dp) :: t = 0.0_dp
        complex(dp), allocatable :: coeff(:)
        real(dp), allocatable :: prob(:)
        real(dp), allocatable :: u_t1(:, :)
        real(dp), allocatable :: h_t1(:, :)
    contains
        procedure :: initialize
        procedure :: propagate
    end type tdse_propagator


contains

    subroutine initialize(self, rng, n_state)
        use matrix_mod, only : unit_mat
        class(tdse_propagator), intent(inout) :: self
        class(rng_type), pointer, intent(in) :: rng
        integer, intent(in) :: n_state

        self%rng => rng
        self%n_state = n_state
        allocate(self%coeff(self%n_state), source=cmplx(0.0_dp, 0.0_dp, kind=dp))
        allocate(self%prob(self%n_state), source=0.0_dp)
        self%representation = 'locally_diabatic'
    end subroutine initialize


    subroutine propagate(self, t2, potential)
        use matrix_mod, only : unitary_transform
        class(tdse_propagator), intent(inout) :: self
        class(potential_evaluator), intent(in) :: potential
        real(dp), intent(in) :: t2 !< Time after propagation
        real(dp) :: dt
        real(dp) :: t(self%n_state, self%n_state) !< Orthogonalized overlap matrix.
        real(dp) :: h_t(self%n_state, self%n_state) !< Hamiltonian during substep.
        real(dp) :: h_t2(self%n_state, self%n_state) !< Hamiltonian at time t2.
        complex(dp) :: u(self%n_state, self%n_state) !< Transformation matrix.
        complex(dp) :: w3(self%n_state, self%n_state) !< Work array 3.
        complex(dp) :: w4(self%n_state, self%n_state) !< Work array 4.
        real(dp), allocatable :: u_t2(:, :)
        integer :: i, cstate
        real(dp) :: cprob
        real(dp) :: rnum
        complex(dp) :: coeff_t(self%n_state) !< Coefficients at time t.

        cstate = self%active_state
        coeff_t = self%coeff
        u_t2 = potential%get_transformation('diabatic', 'adiabatic')
        h_t2 = potential%get_Hamiltonian('adiabatic')

        ! if (t2 > 59.40_dp / 0.02418884326509_dp) then
        !     write(*,*) "here"
        ! end if

        ! Get overlap matrix.
        t = matmul(transpose(self%u_t1), u_t2) ! T = U(t)^t . U(t+dt)
        call orthog_lowdin(t)

        h_t2 = matmul(matmul(t, h_t2), transpose(t)) ! H(t+dt) = T.H(t+dt).T^t
        h_t2 = h_t2 - self%h_t1

        dt = (t2 - self%t) / self%n_substep

        do i = 1, self%n_substep
             h_t = self%h_t1 + h_t2 * real(i - 0.5_dp, dp) / self%n_substep
             w3 = mat_sy_exp(h_t, cmplx(0.0_dp, -dt, kind=dp))
             coeff_t = matmul(w3, coeff_t)
        end do
        
        w4 = cmplx(t, 0.0_dp, dp)
        coeff_t = matmul(transpose(w4), coeff_t)
        w3 = matmul(transpose(w4), w3)
        self%prob = prob_sharc(cstate, self%coeff, coeff_t, w3)

        cprob = 0.0_dp
        call self%rng%uniform(rnum)
        do i = 1, self%n_state
            if (i == cstate) cycle
            cprob = cprob + self%prob(i)
            if (rnum < cprob) then
                self%active_state = i
                exit
            end if
        end do

        self%coeff = coeff_t
        self%u_t1 = u_t2
        self%h_t1 = h_t2
        self%t = t2
    end subroutine propagate


    pure function prob_sharc(b, c_t, c_dt, p) result(prob)
        integer, intent(in) :: b
        complex(dp), intent(in) :: c_t(:)
        complex(dp), intent(in) :: c_dt(:)
        complex(dp), intent(in) :: p(:, :)
        real(dp) :: prob(size(c_t))
        integer :: a

        prob = 0.0_dp
        do a = 1, size(c_t)
            if (a == b) cycle
            prob(a) = max((1 - abs(c_dt(b))**2 / abs(c_t(b))**2) * &
            &         real(c_dt(a) * conjg(p(a, b)) * conjg(c_t(b))) / &
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
        complex(dp), intent(in) :: u(:, :)
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
        complex(dp), intent(in) :: u(:, :)
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


end module tdse_propagator_mod