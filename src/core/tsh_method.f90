!--------------------------------------------------------------------------------------------------
! MODULE: tsh_method_mod
!
! DESCRIPTION:
!> @brief Base type for surface hopping methods.
!> @details
!! A method holds the active state and the electronic data needed to choose it at each new point
!! of the trajectory. The driver calls start at the initial point and propagate at each new point,
!! after the energies and couplings at that point have been evaluated.
!!
!! After each call, the method tells the driver through its request what to do next. The driver
!! can make an earlier stored point of the trajectory current (request_goto), or make a stored
!! point current and add a new point after it with a shorter step (request_insert). The method
!! object itself is carried over to the new current point, so it keeps its own history.
!--------------------------------------------------------------------------------------------------
module tsh_method_mod
    use global_defs, only : dp, stdout, stderr, stdp1, stdp2
    use evaluator_base_mod, only : potential_evaluator, state_location
    use random_mod, only : rng_type
    implicit none

    private
    public :: tsh_method
    public :: tsh_none
    public :: tsh_request
    public :: initialize
    public :: set_initial_state
    public :: request_none
    public :: request_goto
    public :: request_insert
    public :: method_save_state
    public :: method_load_state


    integer, parameter :: request_none = 0 !< Continue from the current point.
    integer, parameter :: request_goto = 1 !< Make the stored point at time current and discard
                                           !! all later points.
    integer, parameter :: request_insert = 2 !< Make the stored point at time current, keep all
                                             !! later points and take the next step with dt.


    !> What the driver should do after the method processed a point.
    type :: tsh_request
        integer :: action = request_none
        real(dp) :: time = 0.0_dp !< Time of the stored point for request_goto/request_insert.
        real(dp) :: dt = 0.0_dp !< Step for request_insert.
        logical :: hop = .false. !< The active state changed at the current point after the action.
        logical :: isc = .false. !< The hop is an intersystem crossing.
        integer :: previous_state = 0 !< Active state before the hop.
        real(dp) :: t_final = -huge(1.0_dp) !< Points up to this time will not change anymore.
    end type tsh_request


    type, abstract :: tsh_method
        class(rng_type), pointer :: rng => null()
        character(len=:), allocatable :: basis !< Basis of active_state and of the nuclear dynamics.
        integer :: n_state = 0 !< Number of states of the basis.
        integer :: active_state = 0
        real(dp) :: t = 0.0_dp !< Time of the last processed point.
        real(dp), allocatable :: prob(:) !< Hopping probabilities from the last step (output).
        type(tsh_request) :: request !< Request from the last processed point.
    contains
        procedure :: initialize
        procedure :: set_initial_state
        procedure(start_iface), deferred :: start
        procedure(propagate_iface), deferred :: propagate
        procedure :: hop_rejected
        procedure :: decohere
        procedure :: write_coefficients
        procedure :: select_state
        procedure :: rescale_velocity
        procedure :: save_state => method_save_state
        procedure :: load_state => method_load_state
    end type tsh_method


    !> Surface hopping turned off. The trajectory stays in the initial state.
    type, extends(tsh_method) :: tsh_none
    contains
        procedure :: start => none_start
        procedure :: propagate => none_propagate
        procedure :: save_state => none_save_state
        procedure :: load_state => none_load_state
    end type tsh_none


    abstract interface
        !> Store the data needed from the initial point of the trajectory and set the request.
        subroutine start_iface(self, t, potential)
            import :: tsh_method, dp, potential_evaluator
            class(tsh_method), intent(inout) :: self
            real(dp), intent(in) :: t
            class(potential_evaluator), intent(in) :: potential
        end subroutine start_iface

        !> Process a new point of the trajectory at time t2, choose the active state and set the
        !! request.
        subroutine propagate_iface(self, t2, potential)
            import :: tsh_method, dp, potential_evaluator
            class(tsh_method), intent(inout) :: self
            real(dp), intent(in) :: t2
            class(potential_evaluator), intent(in) :: potential
        end subroutine propagate_iface
    end interface


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: initialize
    !> @brief Common part of the initialization. Extensions override it and call it first.
    !> @details
    !! The basis is 'adiabatic', unless the extension sets another one before the call. The
    !! active state, the number of states and the probabilities refer to this basis.
    !----------------------------------------------------------------------------------------------
    subroutine initialize(self, rng, potential, active_state)
        class(tsh_method), intent(inout) :: self
        class(rng_type), pointer, intent(in) :: rng
        class(potential_evaluator), intent(in) :: potential
        integer, intent(in) :: active_state

        self%rng => rng
        if (.not. allocated(self%basis)) self%basis = 'adiabatic'
        self%n_state = potential%basis_size(self%basis)
        if (active_state < 1 .or. active_state > self%n_state) then
            write(stderr, *) 'Error in tsh_method_mod, initialize subroutine.'
            write(stderr, '(3x,a,i0,a,i0,3a)') 'Initial state ', active_state,                   &
            &   ' is not between 1 and ', self%n_state, ', the number of states of the ',          &
            &   self%basis, ' basis.'
            stop
        end if
        self%active_state = active_state
        allocate(self%prob(self%n_state), source=0.0_dp)
    end subroutine initialize


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: set_initial_state
    !> @brief Set the active state at the initial point from a state of a given basis.
    !> @details
    !! state is a state of basis, 'diabatic' or 'adiabatic'. A diabatic state is given by its index
    !! in the source of the potential (see diabatic_from_source). If basis is the basis of the
    !! method, the state becomes the active state. Otherwise, the active state a is drawn with
    !! probability |T(i, a)|**2, where T is the transformation from basis to the basis of the
    !! method at the current point (the same weights as zagpy sample state). One random number is
    !! drawn only in this case. With the 'mch' basis of the method, a diabatic state is
    !! projected through its spin-free state.
    !! If the dynamics runs on a subset of the states of the model, a diabatic state is only
    !! partly in the states of the method basis. The weights are then divided by their sum (the
    !! retained weight), and a state without weight in the states of the dynamics is an error.
    !! The potential must be evaluated at the initial point. Extensions that hold more than the
    !! active state override this subroutine and call it first.
    !----------------------------------------------------------------------------------------------
    subroutine set_initial_state(self, potential, state, basis)
        class(tsh_method), intent(inout) :: self
        class(potential_evaluator), intent(in) :: potential
        integer, intent(in) :: state !< Initial state, a state of basis.
        character(len=*), intent(in) :: basis !< Basis of state.
        integer :: i !< Index of the state in basis.
        integer :: a
        real(dp), allocatable :: trans(:, :)
        real(dp), allocatable :: weight(:) !< Probability of each state of the method basis.
        real(dp) :: retained !< Weight of the initial state in the states of the method basis.
        real(dp), parameter :: weight_tol = 1.0e-12_dp !< Smallest retained weight.
        real(dp) :: cprob
        real(dp) :: rnum
        type(state_location) :: loc

        select case(basis)
        case('diabatic')
            i = potential%diabatic_from_source(state)
        case('adiabatic')
            i = state
            if (i < 1 .or. i > potential%basis_size(basis)) i = 0
        case default
            write(stderr, *) 'Error in tsh_method_mod, set_initial_state subroutine.'
            write(stderr, *) '  Unknown basis of the initial state: ', basis
            stop
        end select
        if (i == 0) then
            write(stderr, *) 'Error in tsh_method_mod, set_initial_state subroutine.'
            write(stderr, '(3x,a,i0,3a)') 'Initial state ', state, ' is not a state of the ',     &
            &                             basis, ' basis.'
            stop
        end if

        if (basis == self%basis) then
            self%active_state = i
            if (stdp1) then
                write(stdout, *)
                write(stdout, '(1x,a,i0,3a,i0,a)') 'Initial state ', state, ' of the ', basis,    &
                &                                  ' basis is state ', i, '.'
            end if
            return
        end if

        if (self%basis == 'mch') then
            if (basis /= 'diabatic') then
                write(stderr, *) 'Error in tsh_method_mod, set_initial_state subroutine.'
                write(stderr, *) '  The mch basis is only given in the diabatic basis.'
                write(stderr, *) '  Basis of the initial state: ', basis
                stop
            end if
            loc = potential%model_states%locate_state(i)
            i = loc%sf
            trans = potential%get_transformation('diabatic', 'mch')
        else
            trans = potential%get_transformation(basis, self%basis)
        end if
        weight = trans(i, :)**2
        retained = sum(weight)
        if (retained < weight_tol) then
            write(stderr, *) 'Error in tsh_method_mod, set_initial_state subroutine.'
            write(stderr, '(3x,a,i0,5a)') 'Initial state ', state, ' of the ', basis,             &
            &   ' basis has no weight in the states of the ', self%basis, ' basis of the dynamics.'
            stop
        end if
        weight = weight / retained

        call self%rng%uniform(rnum)
        a = findloc(weight > 0.0_dp, .true., dim=1, back=.true.) ! Guard against rounding.
        cprob = 0.0_dp
        do i = 1, self%n_state
            cprob = cprob + weight(i)
            if (rnum < cprob) then
                a = i
                exit
            end if
        end do
        self%active_state = a

        if (stdp1) then
            write(stdout, *)
            write(stdout, '(1x,a,i0,6a)') 'Initial state ', state, ' of the ', basis, ' basis ',  &
            &                             'projected onto the ', self%basis, ' basis.'
            if (abs(retained - 1.0_dp) > 1.0e-8_dp) then
                write(stdout, '(3x,a,f10.6)') 'Weight retained in the states of the dynamics:',   &
                &                             retained
            end if
            write(stdout, '(3x,a,1000f10.6)') 'Weights:', weight
            write(stdout, '(3x,a,i0,a)') 'Initial active state: ', a, '.'
        end if
    end subroutine set_initial_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: hop_rejected
    !> @brief Undo the hop of the last request, after the driver found it frustrated.
    !----------------------------------------------------------------------------------------------
    subroutine hop_rejected(self)
        class(tsh_method), intent(inout) :: self

        self%active_state = self%request%previous_state
        self%request%hop = .false.
    end subroutine hop_rejected


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: decohere
    !> @brief Apply a decoherence correction before propagating to the next point. No-op by
    !! default.
    !----------------------------------------------------------------------------------------------
    subroutine decohere(self, ekin, dt, potential)
        class(tsh_method), intent(inout) :: self
        real(dp), intent(in) :: ekin !< Kinetic energy.
        real(dp), intent(in) :: dt !< Time step.
        class(potential_evaluator), intent(in) :: potential
    end subroutine decohere


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: write_coefficients
    !> @brief Write the electronic wave function coefficients to unit. No-op by default.
    !----------------------------------------------------------------------------------------------
    subroutine write_coefficients(self, unit)
        class(tsh_method), intent(in) :: self
        integer, intent(in) :: unit
    end subroutine write_coefficients


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: select_state
    !> @brief Choose the new active state from the hopping probabilities in prob.
    !> @details
    !! Exactly one random number is drawn, also when all probabilities are zero.
    !----------------------------------------------------------------------------------------------
    subroutine select_state(self)
        class(tsh_method), intent(inout) :: self
        integer :: i, cstate
        real(dp) :: cprob
        real(dp) :: rnum

        cstate = self%active_state
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
    end subroutine select_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: rescale_velocity
    !> @brief Handle momentum correction after a hop from state pst to the active state.
    !> @details
    !! The behaviour of the subroutine is determined by two options, opt_mc and opt_fh.
    !! The first, opt_mc, determines the direction along which momentum should be rescaled:
    !! - 0 - No velocity rescaling (shouldn't be selected)
    !! - 1 - Rescale along the velocity vector (shouldn't be selected)
    !! - 2 - Rescale along the component along the gradient difference vector
    !! - 3 - Rescale along the component along the nonadiabatic coupling vector.
    !!
    !! The second option, opt_fh, determines what should be done for "frustrated hops" when there
    !! is not enough kinetic energy along the chosen direction to conserve the total energy after
    !! a hop:
    !! - 1 - Return to previous state without changing velocity.
    !! - 2 - Return to previous state, but also invert the velocity along the rescale direction.
    !!       (this option only makes sense with opt_mc=2 or 3)
    !! The return to the previous state is left to the caller, through hop_rejected.
    !----------------------------------------------------------------------------------------------
    subroutine rescale_velocity(self, opt_mc, opt_fh, amask, pst, potential, mass, velo, frustrated)
        class(tsh_method), intent(in) :: self
        integer, intent(in) :: opt_mc !< Type of momentum correction.
        integer, intent(in) :: opt_fh !< Behaviour at frustrated hop.
        integer, intent(in) :: amask(:) !< Atoms considered when rescaling.
        integer, intent(in) :: pst !< Previous state.
        class(potential_evaluator), intent(in) :: potential
        real(dp), intent(in) :: mass(:) !< Masses.
        real(dp), intent(inout) :: velo(:, :) !< Velocities.
        logical, intent(out) :: frustrated !< Not enough energy for the hop.
        real(dp) :: mvel(size(velo, 1), size(amask)) !< Mass weighted velocity.
        real(dp) :: pgrd(size(velo, 1), size(amask)) !< Gradient of the previous state.
        real(dp) :: cgrd(size(velo, 1), size(amask)) !< Gradient of the current state.
        real(dp) :: m(size(velo, 1), size(amask)) !< Temporary mass array.
        real(dp) :: rescale_dir(size(velo, 1), size(amask)) !< Direction along which to rescale.
        real(dp) :: nadv(size(velo, 1)*size(amask)) !< Nonadiabatic coupling vector between the two states.
        real(dp) :: mvel_dir !< Component of mass weighted velocity along rescale direction.
        real(dp) :: delta_e !< Required change in kinetic energy.

        frustrated = .false.
        if (stdp2) write(stdout, '(5x,a)') 'Ensuring energy conservation.'

        ! Work with temporary arrays and use mass-weighted coordinates.
        m = spread(mass(amask), 1, size(velo, 1))
        mvel = velo(:, amask) * sqrt(m)

        select case(opt_mc)
        case(0) ! No rescaling.
            return
        case(1) ! Rescale along velocity vector.
            rescale_dir = mvel
        case(2) ! Rescale along gradient difference vector.
            pgrd = reshape(potential%get_gradient(self%basis, pst), shape(pgrd))
            cgrd = reshape(potential%get_gradient(self%basis, self%active_state), shape(cgrd))
            rescale_dir = (pgrd - cgrd) / sqrt(m)
        case(3) ! Rescale along nonadiabatic coupling vector.
            nadv = potential%get_nadv(self%basis, pst, self%active_state)
            rescale_dir = reshape(nadv, shape(mvel)) / sqrt(m)
        end select

        ! Rescale velocity
        delta_e = potential%get_energy(self%basis, self%active_state)                             &
        &       - potential%get_energy(self%basis, pst)
        rescale_dir = rescale_dir / sqrt(sum(rescale_dir**2))
        mvel_dir = sum(rescale_dir * mvel)
        if (mvel_dir**2 > 2 * delta_e) then
            ! The sign here is based on the work of Herman (see chapter 6 of 10.1007/0-306-46949-9
            ! and references within). Choosing the other direction would also work, but this one
            ! ensures that the change in momentum is as small as possible (but also means that the
            ! direction of the change doesn't depend on the sign of the vector along which rescaling
            ! is performed).
            mvel = mvel + rescale_dir * (sign(sqrt(mvel_dir**2 - 2*delta_e), mvel_dir) - mvel_dir)
        else
            if (stdp1) then
                write(stdout, '(7x,a)') 'Insufficient energy for hop.'
                write(stdout, '(7x,a,e16.8)') 'Energy difference:', delta_e
                write(stdout, '(7x,a,e16.8)') 'Available energy:', mvel_dir**2 / 2
            end if

            frustrated = .true.

            select case(opt_fh)
            case(1) ! Just return to previous state.
                continue
            case(2) ! Invert velocity along rescale direction.
                mvel = mvel - 2 * mvel_dir * rescale_dir
            end select
        end if

        ! Save changes to velocity array.
        velo(:, amask) = mvel / sqrt(m)
    end subroutine rescale_velocity


    subroutine none_start(self, t, potential)
        class(tsh_none), intent(inout) :: self
        real(dp), intent(in) :: t
        class(potential_evaluator), intent(in) :: potential

        self%t = t
        self%request = tsh_request(t_final=t)
    end subroutine none_start


    subroutine none_propagate(self, t2, potential)
        class(tsh_none), intent(inout) :: self
        real(dp), intent(in) :: t2
        class(potential_evaluator), intent(in) :: potential

        self%t = t2
        self%request = tsh_request(previous_state=self%active_state, t_final=t2)
    end subroutine none_propagate


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: method_save_state
    !> @brief Add the state of the method to a JSON object, for a backup.
    !> @details
    !! Only what changes during the dynamics is saved, in the object 'tsh': the active state, the
    !! time and probabilities of the last processed point and the last request. The basis and the
    !! number of states are added for checking. Options come from the input when the backup is
    !! loaded. Extensions add their own object and call the subroutine of their parent first. The
    !! random number generator is saved separately.
    !----------------------------------------------------------------------------------------------
    subroutine method_save_state(self, json, obj)
        use json_module, only : json_core, json_value
        class(tsh_method), intent(in) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the state is added.
        type(json_value), pointer :: tsh, req

        call json%create_object(tsh, 'tsh')
        call json%add(tsh, 'basis', self%basis)
        call json%add(tsh, 'n_state', self%n_state)
        call json%add(tsh, 'active_state', self%active_state)
        call json%add(tsh, 't', self%t)
        call json%add(tsh, 'prob', self%prob)
        call json%create_object(req, 'request')
        call json%add(req, 'action', self%request%action)
        call json%add(req, 'time', self%request%time)
        call json%add(req, 'dt', self%request%dt)
        call json%add(req, 'hop', self%request%hop)
        call json%add(req, 'isc', self%request%isc)
        call json%add(req, 'previous_state', self%request%previous_state)
        call json%add(req, 't_final', self%request%t_final)
        call json%add(tsh, req)
        call json%add(obj, tsh)
    end subroutine method_save_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: method_load_state
    !> @brief Restore the state of the method saved by save_state.
    !> @details
    !! The method must be initialized from the same input as the saved one. Stops if the basis or
    !! the number of states differ.
    !----------------------------------------------------------------------------------------------
    subroutine method_load_state(self, json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_get_required, json_get_child, json_check_value
        class(tsh_method), intent(inout) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object holding the state.
        type(json_value), pointer :: tsh, req

        tsh => json_get_child(json, obj, 'tsh')
        call json_check_value(json, tsh, 'basis', self%basis)
        call json_check_value(json, tsh, 'n_state', self%n_state)
        call json_get_required(json, tsh, 'active_state', self%active_state)
        call json_get_required(json, tsh, 't', self%t)
        call json_get_required(json, tsh, 'prob', self%prob, self%n_state)
        req => json_get_child(json, tsh, 'request')
        call json_get_required(json, req, 'action', self%request%action)
        call json_get_required(json, req, 'time', self%request%time)
        call json_get_required(json, req, 'dt', self%request%dt)
        call json_get_required(json, req, 'hop', self%request%hop)
        call json_get_required(json, req, 'isc', self%request%isc)
        call json_get_required(json, req, 'previous_state', self%request%previous_state)
        call json_get_required(json, req, 't_final', self%request%t_final)
    end subroutine method_load_state


    !> State of the method without surface hopping, with an empty object 'none' marking the method.
    subroutine none_save_state(self, json, obj)
        use json_module, only : json_core, json_value
        class(tsh_none), intent(in) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the state is added.
        type(json_value), pointer :: none

        call method_save_state(self, json, obj)
        call json%create_object(none, 'none')
        call json%add(obj, none)
    end subroutine none_save_state


    !> Restore the state saved by none_save_state. Stops if the backup is of another method.
    subroutine none_load_state(self, json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_get_child
        class(tsh_none), intent(inout) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object holding the state.
        type(json_value), pointer :: none

        call method_load_state(self, json, obj)
        none => json_get_child(json, obj, 'none')
    end subroutine none_load_state


end module tsh_method_mod
