!--------------------------------------------------------------------------------------------------
! MODULE: external_json_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2026
!
! DESCRIPTION:
!> @brief Potential computed by an external program, through JSON files.
!> @details
!! At each point, the evaluator writes a request (request_file) with the geometry, the states and
!! the derivatives which are needed, runs the command given in the input, and reads the response
!! (response_file). The program runs in the working directory of the trajectory and keeps nothing
!! between calls: the evaluator of each point stores the wave functions of its states, and the
!! overlaps with the states of the previous point are computed here, from the wave functions of
!! both points.
!!
!! The signs of the states returned by the program are arbitrary. With a previous point, they are
!! matched to the signs of its states (match_phase), and the nonadiabatic coupling vectors get the
!! same signs. At the first point, the signs are taken as returned.
!!
!! Only singlet states are supported: the states of the program are then the adiabatic states,
!! and the 'mch' and 'adiabatic' bases are the same. The diabatic basis is not available.
!!
!! Protocol (version 1). All quantities are in atomic units, indices start at 1, and geometries,
!! gradients and nonadiabatic coupling vectors are nested arrays with one inner array per atom
!! (or per coordinate label): v[atom][dim].
!!
!! Request:
!! - protocol ("zaghop"), version (1), point_id.
!! - calculation: "new_point", or "same_point" for more derivatives at the geometry of the last
!!   new point (the program must not change the stored results of the new point).
!! - coordinates: coordinate system of the geometry (e.g. "cartesian"), atoms, geometry.
!! - states: [{"multiplicity", "nstate"}].
!! - request: gradients [{"multiplicity", "state"}], nacv [{"multiplicity", "states": [i, j]}],
!!   wave_function (logical), properties (names of extra quantities, not used by the evaluator).
!!
!! Response:
!! - protocol, version, point_id (the same as in the request).
!! - states: [{"multiplicity", "energies"}].
!! - gradients: [{"multiplicity", "state", "gradient"}].
!! - nacv: [{"multiplicity", "states": [i, j], "vector"}], vector = <i|d/dR j>.
!! - wave_function: [{"multiplicity", "format", ...}]. Format "vectors": "c", the coefficients of
!!   the states in an orthonormal basis which is the same at all points, c[basis][state]. The
!!   overlaps of two points are C(prev)^T C.
!! - properties: object with the requested extra quantities.
!--------------------------------------------------------------------------------------------------
module external_json_mod
    use global_defs
    use evaluator_base_mod
    implicit none

    private
    public :: external_json_evaluator
    public :: request_file
    public :: response_file


    character(len=*), parameter :: request_file = 'qm_in.json' !< Request written by the evaluator.
    character(len=*), parameter :: response_file = 'qm_out.json' !< Response of the program.
    character(len=*), parameter :: protocol_name = 'zaghop'
    integer, parameter :: protocol_version = 1
    real(dp), parameter :: overlap_warn_tol = 0.9_dp !< Smallest singular value of the overlaps
    !! between the states of two points without a warning.
    real(dp), parameter :: energy_warn_tol = 1.0e-6_dp !< Largest change of an energy in a
    !! same_point calculation without a warning (Eh).


    !----------------------------------------------------------------------------------------------
    !> TYPE: external_json_evaluator
    !> @brief Potential computed by an external program (see the module description).
    !----------------------------------------------------------------------------------------------
    type, extends(potential_evaluator) :: external_json_evaluator
        character(len=:), allocatable :: command !< Command which runs the program.
        character(len=:), allocatable :: atoms(:) !< Atom symbols (or coordinate labels).
        character(len=:), allocatable :: coordinates !< Coordinates of the geometry.
        real(dp), allocatable :: geometry(:, :) !< Current geometry (ndim, natom).
        integer :: ndof = 0 !< Number of coordinates, size(geometry).
        real(dp), allocatable :: energy(:) !< Energies of the states (n_state).
        real(dp), allocatable :: grad(:, :) !< Gradients (ndof, n_state).
        logical, allocatable :: have_grad(:) !< The gradient of a state was computed.
        real(dp), allocatable :: nadv(:, :, :) !< Nonadiabatic coupling vectors <i|d/dR j>,
                                                !! with the phases of wf (ndof, n_state, n_state).
        logical, allocatable :: have_nadv(:, :) !< The coupling of a pair was computed.
        character(len=:), allocatable :: wf_format !< Format of the wave functions.
        real(dp), allocatable :: wf(:, :) !< Wave functions of the states, with matched phases
                                          !! (format 'vectors': coefficients, (n_basis, n_state)).
    contains
        procedure :: init => external_init
        procedure :: set_atoms
        procedure :: set_model_states => external_set_model_states
        procedure :: update_geometry => external_update_geometry
        procedure :: eval => external_eval
        procedure :: complete => external_complete
        procedure :: get_Hamiltonian => external_get_Hamiltonian
        procedure :: get_transformation => external_get_transformation
        procedure :: get_energy_single => external_get_energy_single
        procedure :: get_energy_all => external_get_energy_all
        procedure :: get_gradient => external_get_gradient
        procedure :: get_nadv => external_get_nadv
        procedure :: get_soc2 => external_get_soc2
        procedure :: save_state => external_save_state
        procedure :: load_state => external_load_state
        procedure, private :: run_program
        procedure, private :: write_request
        procedure, private :: read_response
        procedure, private :: add_derivatives
        procedure, private :: check_basis
    end type external_json_evaluator


    !> Results read from a response.
    type :: qm_response
        real(dp), allocatable :: energy(:) !< Energies of the states.
        integer, allocatable :: grad_state(:) !< States of the gradients.
        real(dp), allocatable :: grad(:, :) !< Gradients (ndof, size(grad_state)).
        integer, allocatable :: nadv_pair(:, :) !< Pairs of the couplings (2, n_pair).
        real(dp), allocatable :: nadv(:, :) !< Coupling vectors <i|d/dR j> (ndof, n_pair).
        character(len=:), allocatable :: wf_format !< Format of the wave functions.
        real(dp), allocatable :: wf(:, :) !< Wave functions of the states.
    end type qm_response


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: external_init
    !> @brief Set the command which runs the program.
    !> @details
    !! The geometry is passed in the coordinates of the trajectory, by default Cartesian. The
    !! states are set with set_model_states, and the atoms with set_atoms.
    !----------------------------------------------------------------------------------------------
    subroutine external_init(self, command)
        class(external_json_evaluator), intent(inout) :: self
        character(len=*), intent(in) :: command

        if (len_trim(command) == 0) then
            write(stderr, *) 'Error in external_json_mod, init subroutine.'
            write(stderr, *) '  No command given for the external program.'
            stop
        end if
        self%command = trim(command)
        self%classical_coordinates = 'cartesian'
    end subroutine external_init


    !> Set the atom symbols (or coordinate labels) written in the request.
    subroutine set_atoms(self, sym)
        class(external_json_evaluator), intent(inout) :: self
        character(len=*), intent(in) :: sym(:)
        integer :: i, n

        n = 1
        do i = 1, size(sym)
            n = max(n, len_trim(sym(i)))
        end do
        allocate(character(len=n) :: self%atoms(size(sym)))
        do i = 1, size(sym)
            self%atoms(i) = trim(sym(i))
        end do
    end subroutine set_atoms


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: external_set_model_states
    !> @brief Set the states from the number of states of each multiplicity. Only singlets.
    !----------------------------------------------------------------------------------------------
    subroutine external_set_model_states(self, nstate, multiplicity, spin_states)
        class(external_json_evaluator), intent(inout) :: self
        integer, intent(in) :: nstate(:) !< Number of spin-free states of each multiplicity.
        integer, intent(in) :: multiplicity(:) !< Multiplicity of each block of states.
        character(len=*), intent(in) :: spin_states

        if (any(multiplicity /= 1 .and. nstate > 0)) then
            write(stderr, *) 'Error in external_json_mod, set_model_states subroutine.'
            write(stderr, *) '  Only singlet states are available for potentials computed by an ', &
            &                'external program.'
            stop
        end if
        call evaluator_set_model_states(self, nstate, multiplicity, spin_states)
    end subroutine external_set_model_states


    !> Store the geometry for the next calculation. Derivatives are returned in the same
    !! coordinates, so get_derivative_target_coord is the identity.
    subroutine external_update_geometry(self, geometry, coordinates)
        class(external_json_evaluator), intent(inout) :: self
        real(dp), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates

        self%geometry = geometry
        self%coordinates = coordinates
        self%ndof = size(geometry)
    end subroutine external_update_geometry


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: external_eval
    !> @brief Compute the energies, the requested derivatives and the wave functions at the
    !! current geometry, and the overlaps with the states of prev.
    !> @details
    !! With prev, the phases of the states are matched to those of prev, and link holds the
    !! overlaps C(prev)^T C between the states of prev and of this point, after matching.
    !----------------------------------------------------------------------------------------------
    subroutine external_eval(self, prev)
        use state_overlap_mod, only : match_phase
        class(external_json_evaluator), intent(inout) :: self
        class(potential_evaluator), intent(in), optional :: prev
        type(eval_request) :: req
        type(qm_response) :: res
        integer :: n
        integer, allocatable :: phase(:)
        real(dp), allocatable :: s(:, :)
        real(dp) :: sigma_min

        call self%new_point()
        if (.not. allocated(self%geometry)) then
            write(stderr, *) 'Error in external_json_mod, eval subroutine.'
            write(stderr, *) '  No geometry given.'
            stop
        end if
        n = self%states%n_state
        ! The evaluator of a new point is a copy of that of the previous point.
        if (allocated(self%energy)) deallocate(self%energy)
        if (allocated(self%grad)) deallocate(self%grad)
        if (allocated(self%have_grad)) deallocate(self%have_grad)
        if (allocated(self%nadv)) deallocate(self%nadv)
        if (allocated(self%have_nadv)) deallocate(self%have_nadv)
        if (allocated(self%wf)) deallocate(self%wf)
        allocate(self%grad(self%ndof, n), source=0.0_dp)
        allocate(self%have_grad(n), source=.false.)
        allocate(self%nadv(self%ndof, n, n), source=0.0_dp)
        allocate(self%have_nadv(n, n), source=.false.)

        if (allocated(self%request%basis)) then
            call self%check_basis(self%request%basis, 'eval')
            req = self%request
        else
            req = new_request('adiabatic')
        end if
        res = self%run_program('new_point', req, .true.)
        self%energy = res%energy
        self%wf_format = res%wf_format
        self%wf = res%wf
        allocate(phase(n), source=1)

        if (present(prev)) then
            select type(prev)
            class is (external_json_evaluator)
                if (.not. allocated(prev%wf)) then
                    write(stderr, *) 'Error in external_json_mod, eval subroutine.'
                    write(stderr, *) '  The previous point was not evaluated.'
                    stop
                end if
                if (any(shape(prev%wf) /= shape(self%wf))                                         &
                &   .or. prev%wf_format /= self%wf_format) then
                    write(stderr, *) 'Error in external_json_mod, eval subroutine.'
                    write(stderr, *) '  The wave functions of the previous point have a ',         &
                    &                'different format or size.'
                    stop
                end if
                call match_phase(prev%wf, self%wf, s, phase)
                call self%link%set(prev%id, self%id, 'adiabatic', s, sigma_min)
                if (sigma_min < overlap_warn_tol) then
                    write(stderr, *) 'Warning in external_json_mod, eval subroutine.'
                    write(stderr, '(3x,a,f8.5,a)') 'Smallest singular value of the overlaps '//   &
                    &   'between the states of two points: ', sigma_min, '.'
                    write(stderr, *) '  A state changed character, for example at a crossing ',   &
                    &                'with a state which is not computed.'
                end if
            class default
                write(stderr, *) 'Error in external_json_mod, eval subroutine.'
                write(stderr, *) '  The previous point is not computed by an external program.'
                stop
            end select
        end if
        call self%add_derivatives(res, phase)
    end subroutine external_eval


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: external_complete
    !> @brief Compute the requested derivatives which are missing at the evaluated point.
    !> @details
    !! The program is run with calculation "same_point". For nonadiabatic coupling vectors, the
    !! wave functions are requested as well, and the signs of the new couplings are matched to the
    !! stored states, whose phases don't change. The energies are kept; a change larger than
    !! energy_warn_tol gives a warning.
    !----------------------------------------------------------------------------------------------
    subroutine external_complete(self, req)
        use state_overlap_mod, only : match_phase
        class(external_json_evaluator), intent(inout) :: self
        type(eval_request), intent(in) :: req
        type(eval_request) :: missing
        type(qm_response) :: res
        integer, allocatable :: phase(:)
        integer :: i, k

        if (self%id == 0 .or. .not. allocated(self%energy)) then
            write(stderr, *) 'Error in external_json_mod, complete subroutine.'
            write(stderr, *) '  The potential was not evaluated.'
            stop
        end if
        call self%check_basis(req%basis, 'complete')
        missing = new_request(req%basis)
        do i = 1, size(req%gradient)
            if (.not. self%have_grad(req%gradient(i))) then
                missing%gradient = [missing%gradient, req%gradient(i)]
            end if
        end do
        do k = 1, size(req%nadv, 2)
            if (.not. self%have_nadv(req%nadv(1, k), req%nadv(2, k))) then
                missing%nadv = reshape([missing%nadv, req%nadv(:, k)],                            &
                &                      [2, size(missing%nadv, 2) + 1])
            end if
        end do
        if (size(missing%gradient) == 0 .and. size(missing%nadv, 2) == 0) return

        res = self%run_program('same_point', missing, size(missing%nadv, 2) > 0)
        if (maxval(abs(res%energy - self%energy)) > energy_warn_tol) then
            write(stderr, *) 'Warning in external_json_mod, complete subroutine.'
            write(stderr, '(3x,a,es10.3,a)') 'Energies changed by up to ',                        &
            &   maxval(abs(res%energy - self%energy)), ' Eh in a calculation at the same point.'
        end if
        allocate(phase(self%states%n_state), source=1)
        if (allocated(res%wf)) then
            if (any(shape(res%wf) /= shape(self%wf)) .or. res%wf_format /= self%wf_format) then
                write(stderr, *) 'Error in external_json_mod, complete subroutine.'
                write(stderr, *) '  The wave functions have a different format or size than ',    &
                &                'at the first calculation of the point.'
                stop
            end if
            call match_phase(self%wf, res%wf, phase=phase)
        end if
        call self%add_derivatives(res, phase)
    end subroutine external_complete


    !> Store the derivatives of a response, with the signs of the states changed by phase.
    subroutine add_derivatives(self, res, phase)
        class(external_json_evaluator), intent(inout) :: self
        type(qm_response), intent(in) :: res
        integer, intent(in) :: phase(:) !< Sign applied to each state of the response.
        integer :: k, i, j

        do k = 1, size(res%grad_state)
            self%grad(:, res%grad_state(k)) = res%grad(:, k)
            self%have_grad(res%grad_state(k)) = .true.
        end do
        do k = 1, size(res%nadv_pair, 2)
            i = res%nadv_pair(1, k)
            j = res%nadv_pair(2, k)
            self%nadv(:, i, j) = phase(i) * phase(j) * res%nadv(:, k)
            self%nadv(:, j, i) = -self%nadv(:, i, j)
            self%have_nadv(i, j) = .true.
            self%have_nadv(j, i) = .true.
        end do
    end subroutine add_derivatives


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: run_program
    !> @brief Write the request, run the program and read its response.
    !> @details
    !! The response file is removed before the program is run. Stops if the program fails (nonzero
    !! exit status) or the response is missing or incomplete. Both files are then left in the
    !! working directory.
    !----------------------------------------------------------------------------------------------
    function run_program(self, calculation, req, want_wf) result(res)
        class(external_json_evaluator), intent(in) :: self
        character(len=*), intent(in) :: calculation !< 'new_point' or 'same_point'.
        type(eval_request), intent(in) :: req !< Derivatives to compute.
        logical, intent(in) :: want_wf !< Request the wave functions.
        type(qm_response) :: res
        integer :: exitstat, cmdstat, u
        character(len=256) :: cmdmsg
        logical :: exists

        call self%write_request(calculation, req, want_wf)
        inquire(file=response_file, exist=exists)
        if (exists) then
            open(newunit=u, file=response_file)
            close(u, status='delete')
        end if

        if (stdp2) write(stdout, '(5x,a)') 'Running external program ('//calculation//').'
        exitstat = 0
        cmdstat = 0
        cmdmsg = ''
        call execute_command_line(self%command, wait=.true., exitstat=exitstat, cmdstat=cmdstat,  &
        &                         cmdmsg=cmdmsg)
        if (cmdstat /= 0 .or. exitstat /= 0) then
            write(stderr, *) 'Error in external_json_mod, run_program function.'
            write(stderr, *) '  External program failed: ', self%command
            if (cmdstat /= 0) write(stderr, '(3x,a,i0,2a)') 'cmdstat = ', cmdstat, ': ',          &
            &                                               trim(cmdmsg)
            if (cmdstat == 0) write(stderr, '(3x,a,i0)') 'Exit status: ', exitstat
            stop
        end if
        inquire(file=response_file, exist=exists)
        if (.not. exists) then
            write(stderr, *) 'Error in external_json_mod, run_program function.'
            write(stderr, *) '  File ', response_file, ' not found after running: ', self%command
            stop
        end if
        res = self%read_response(req, want_wf)
    end function run_program


    !> Write the request for a calculation at the current geometry to request_file.
    subroutine write_request(self, calculation, req, want_wf)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_add_array
        class(external_json_evaluator), intent(in) :: self
        character(len=*), intent(in) :: calculation
        type(eval_request), intent(in) :: req
        logical, intent(in) :: want_wf
        type(json_core) :: json
        type(json_value), pointer :: root, arr, obj, rq
        integer :: i, k

        if (.not. allocated(self%atoms)) then
            write(stderr, *) 'Error in external_json_mod, write_request subroutine.'
            write(stderr, *) '  Atoms not set.'
            stop
        end if
        if (size(self%atoms) /= size(self%geometry, 2)) then
            write(stderr, *) 'Error in external_json_mod, write_request subroutine.'
            write(stderr, '(3x,a,i0,a,i0,a)') 'Number of atoms (', size(self%atoms),            &
            &   ') does not match the geometry (', size(self%geometry, 2), ').'
            stop
        end if

        call json%initialize()
        call json%create_object(root, '')
        call json%add(root, 'protocol', protocol_name)
        call json%add(root, 'version', protocol_version)
        call json%add(root, 'point_id', self%id)
        call json%add(root, 'calculation', calculation)
        call json%add(root, 'coordinates', self%coordinates)
        call json%create_array(arr, 'atoms')
        do i = 1, size(self%atoms)
            call json%add(arr, '', trim(self%atoms(i)))
        end do
        call json%add(root, arr)
        call json_add_array(json, root, 'geometry', transpose(self%geometry))

        call json%create_array(arr, 'states')
        call json%create_object(obj, '')
        call json%add(obj, 'multiplicity', 1)
        call json%add(obj, 'nstate', self%states%n_state)
        call json%add(arr, obj)
        call json%add(root, arr)

        call json%create_object(rq, 'request')
        call json%create_array(arr, 'gradients')
        do k = 1, size(req%gradient)
            call json%create_object(obj, '')
            call json%add(obj, 'multiplicity', 1)
            call json%add(obj, 'state', req%gradient(k))
            call json%add(arr, obj)
        end do
        call json%add(rq, arr)
        call json%create_array(arr, 'nacv')
        do k = 1, size(req%nadv, 2)
            call json%create_object(obj, '')
            call json%add(obj, 'multiplicity', 1)
            call json%add(obj, 'states', req%nadv(:, k))
            call json%add(arr, obj)
        end do
        call json%add(rq, arr)
        call json%add(rq, 'wave_function', want_wf)
        call json%create_array(arr, 'properties')
        call json%add(rq, arr)
        call json%add(root, rq)

        call json%print(root, request_file)
        call json%destroy(root)
        if (json%failed()) then
            write(stderr, *) 'Error in external_json_mod, write_request subroutine.'
            call json%print_error_message(stderr)
            stop
        end if
    end subroutine write_request


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: read_response
    !> @brief Read and check the response to a request.
    !> @details
    !! Stops if the response is for another point, has a different number of states, or misses a
    !! requested quantity. Quantities which were not requested are ignored.
    !----------------------------------------------------------------------------------------------
    function read_response(self, req, want_wf) result(res)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_get_required, json_get_child, json_check_value
        class(external_json_evaluator), intent(in) :: self
        type(eval_request), intent(in) :: req
        logical, intent(in) :: want_wf
        type(qm_response) :: res
        type(json_core) :: json
        type(json_value), pointer :: root, arr, item
        integer :: i, k, n, mult, state
        integer, allocatable :: pair(:)
        real(dp), allocatable :: v(:, :)
        logical :: found

        n = self%states%n_state
        call json%initialize()
        call json%load(file=response_file, p=root)
        if (json%failed()) then
            write(stderr, *) 'Error in external_json_mod, read_response function.'
            write(stderr, *) '  Could not read ', response_file, ':'
            call json%print_error_message(stderr)
            stop
        end if
        call json_check_value(json, root, 'protocol', protocol_name)
        call json_check_value(json, root, 'version', protocol_version)
        call json_get_required(json, root, 'point_id', i)
        if (i /= self%id) then
            write(stderr, *) 'Error in external_json_mod, read_response function.'
            write(stderr, '(3x,2a,i0,a,i0,a)') response_file, ': response for point ', i,          &
            &                                  ', expected point ', self%id, '.'
            stop
        end if

        ! Energies.
        arr => json_get_child(json, root, 'states')
        found = .false.
        do i = 1, json%count(arr)
            call json%get_child(arr, i, item)
            call json_get_required(json, item, 'multiplicity', mult)
            if (mult /= 1) cycle
            call json_get_required(json, item, 'energies', res%energy, n)
            found = .true.
        end do
        if (.not. found) call response_error('No energies of the singlet states.')

        ! Gradients, in the order of the request.
        allocate(res%grad_state(size(req%gradient)), source=req%gradient)
        allocate(res%grad(self%ndof, size(req%gradient)))
        do k = 1, size(req%gradient)
            found = .false.
            if (json%valid_path(root, 'gradients')) then
                arr => json_get_child(json, root, 'gradients')
                do i = 1, json%count(arr)
                    call json%get_child(arr, i, item)
                    call json_get_required(json, item, 'multiplicity', mult)
                    call json_get_required(json, item, 'state', state)
                    if (mult /= 1 .or. state /= req%gradient(k)) cycle
                    call json_get_required(json, item, 'gradient', v,                             &
                    &                      shape(transpose(self%geometry)))
                    res%grad(:, k) = reshape(transpose(v), [self%ndof])
                    found = .true.
                    exit
                end do
            end if
            if (.not. found) call response_error('Missing gradient of a requested state.',        &
            &                                    req%gradient(k:k))
        end do

        ! Nonadiabatic coupling vectors, in the order of the request.
        allocate(res%nadv_pair(2, size(req%nadv, 2)), source=req%nadv)
        allocate(res%nadv(self%ndof, size(req%nadv, 2)))
        do k = 1, size(req%nadv, 2)
            found = .false.
            if (json%valid_path(root, 'nacv')) then
                arr => json_get_child(json, root, 'nacv')
                do i = 1, json%count(arr)
                    call json%get_child(arr, i, item)
                    call json_get_required(json, item, 'multiplicity', mult)
                    call json_get_required(json, item, 'states', pair, 2)
                    if (mult /= 1) cycle
                    if (all(pair == req%nadv(:, k))) then
                        call json_get_required(json, item, 'vector', v,                           &
                        &                      shape(transpose(self%geometry)))
                        res%nadv(:, k) = reshape(transpose(v), [self%ndof])
                    else if (all(pair == req%nadv([2, 1], k))) then
                        ! <j|d/dR i> = -<i|d/dR j>.
                        call json_get_required(json, item, 'vector', v,                           &
                        &                      shape(transpose(self%geometry)))
                        res%nadv(:, k) = -reshape(transpose(v), [self%ndof])
                    else
                        cycle
                    end if
                    found = .true.
                    exit
                end do
            end if
            if (.not. found) call response_error('Missing coupling vector of a requested pair.',  &
            &                                    req%nadv(:, k))
        end do

        ! Wave functions.
        if (want_wf) then
            found = .false.
            if (json%valid_path(root, 'wave_function')) then
                arr => json_get_child(json, root, 'wave_function')
                do i = 1, json%count(arr)
                    call json%get_child(arr, i, item)
                    call json_get_required(json, item, 'multiplicity', mult)
                    if (mult /= 1) cycle
                    call json_get_required(json, item, 'format', res%wf_format)
                    select case(res%wf_format)
                    case('vectors')
                        call json_get_required(json, item, 'c', res%wf)
                        if (size(res%wf, 2) /= n) then
                            call response_error('Wrong number of states in the wave functions.')
                        end if
                    case default
                        call response_error('Unknown wave function format: '//res%wf_format)
                    end select
                    found = .true.
                    exit
                end do
            end if
            if (.not. found) call response_error('No wave functions of the singlet states.')
        end if
        call json%destroy(root)

    contains

        subroutine response_error(msg, states)
            character(len=*), intent(in) :: msg
            integer, intent(in), optional :: states(:)

            write(stderr, *) 'Error in external_json_mod, read_response function.'
            write(stderr, *) '  ', response_file, ': ', msg
            if (present(states)) write(stderr, '(3x,a,100(1x,i0))') 'States:', states
            stop
        end subroutine response_error

    end function read_response


    !> Stop unless basis is one of the bases of the potential ('adiabatic' or 'mch').
    subroutine check_basis(self, basis, caller)
        class(external_json_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        character(len=*), intent(in) :: caller

        select case(basis)
        case('adiabatic', 'mch')
        case default
            write(stderr, *) 'Error in external_json_mod, ', caller, '.'
            write(stderr, *) '  Basis "', basis, '" is not available for potentials computed by ', &
            &                'an external program. Expected adiabatic or mch.'
            stop
        end select
    end subroutine check_basis


    !> Hamiltonian in the states of the program: the diagonal matrix of the energies.
    function external_get_Hamiltonian(self, basis) result(H)
        use matrix_mod, only : diag
        class(external_json_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        real(dp), allocatable :: H(:, :)

        call self%check_basis(basis, 'get_Hamiltonian')
        H = diag(self%energy)
    end function external_get_Hamiltonian


    !> Transformation between the 'adiabatic' and 'mch' bases: the unit matrix.
    function external_get_transformation(self, from_basis, to_basis) result(trans)
        use matrix_mod, only : unit_mat
        class(external_json_evaluator), intent(in) :: self
        character(len=*), intent(in) :: from_basis, to_basis
        real(dp), allocatable :: trans(:, :)

        call self%check_basis(from_basis, 'get_transformation')
        call self%check_basis(to_basis, 'get_transformation')
        trans = unit_mat(self%states%n_state)
    end function external_get_transformation


    pure function external_get_energy_single(self, basis, state) result(energy)
        class(external_json_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: state
        real(dp) :: energy

        select case(basis)
        case('adiabatic', 'mch')
            energy = self%energy(state)
        case default
            error stop 'Error in external_json_mod, get_energy_single function. Unknown basis.'
        end select
    end function external_get_energy_single


    pure function external_get_energy_all(self, basis) result(energy)
        class(external_json_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        real(dp), allocatable :: energy(:)

        select case(basis)
        case('adiabatic', 'mch')
            energy = self%energy
        case default
            error stop 'Error in external_json_mod, get_energy_all function. Unknown basis.'
        end select
    end function external_get_energy_all


    !> Gradient of a state. Stops if it was not computed at this point (see complete).
    function external_get_gradient(self, basis, istate) result(grad)
        class(external_json_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: istate
        real(dp), allocatable :: grad(:)

        call self%check_basis(basis, 'get_gradient')
        if (.not. self%have_grad(istate)) then
            write(stderr, *) 'Error in external_json_mod, get_gradient function.'
            write(stderr, '(3x,a,i0,a)') 'Gradient of state ', istate,                           &
            &                            ' was not computed at this point.'
            stop
        end if
        grad = self%grad(:, istate)
    end function external_get_gradient


    !> Nonadiabatic coupling vector <istate1|d/dR istate2>. Stops if it was not computed at this
    !! point (see complete).
    function external_get_nadv(self, basis, istate1, istate2) result(nadv)
        class(external_json_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: istate1, istate2
        real(dp), allocatable :: nadv(:)

        call self%check_basis(basis, 'get_nadv')
        if (.not. self%have_nadv(istate1, istate2)) then
            write(stderr, *) 'Error in external_json_mod, get_nadv function.'
            write(stderr, '(3x,a,i0,a,i0,a)') 'Coupling vector between states ', istate1, ' and ', &
            &                                 istate2, ' was not computed at this point.'
            stop
        end if
        nadv = self%nadv(:, istate1, istate2)
    end function external_get_nadv


    !> Squared spin-orbit couplings: zero, only singlet states are available.
    function external_get_soc2(self, basis) result(s2)
        class(external_json_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        real(dp), allocatable :: s2(:, :)

        call self%check_basis(basis, 'get_soc2')
        allocate(s2(self%states%n_state, self%states%n_state), source=0.0_dp)
    end function external_get_soc2


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: external_save_state
    !> @brief Add the state of the evaluated point to a JSON object, in the object
    !! 'external_json'.
    !> @details
    !! The computed quantities and the wave functions (with their matched phases) are saved, so a
    !! loaded point needs no calculation. Which gradients and couplings were computed is saved as
    !! 0/1 flags (have_grad, and have_nadv by columns).
    !----------------------------------------------------------------------------------------------
    subroutine external_save_state(self, json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_add_array
        class(external_json_evaluator), intent(in) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the state is added.
        type(json_value), pointer :: e

        call evaluator_save_state(self, json, obj)
        call json%create_object(e, 'external_json')
        call json%add(e, 'coordinates', self%coordinates)
        call json_add_array(json, e, 'geometry', self%geometry)
        call json%add(e, 'energy', self%energy)
        call json_add_array(json, e, 'grad', self%grad)
        call json%add(e, 'have_grad', merge(1, 0, self%have_grad))
        call json_add_array(json, e, 'nadv', self%nadv)
        call json%add(e, 'have_nadv', merge(1, 0, reshape(self%have_nadv, [size(self%have_nadv)])))
        call json%add(e, 'wf_format', self%wf_format)
        call json_add_array(json, e, 'wf', self%wf)
        call json%add(obj, e)
    end subroutine external_save_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: external_load_state
    !> @brief Restore the state of an evaluated point saved by save_state.
    !----------------------------------------------------------------------------------------------
    subroutine external_load_state(self, json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_get_required, json_get_child
        class(external_json_evaluator), intent(inout) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object holding the state.
        type(json_value), pointer :: e
        integer, allocatable :: flags(:)
        integer :: n

        call evaluator_load_state(self, json, obj)
        e => json_get_child(json, obj, 'external_json')
        n = self%states%n_state
        call json_get_required(json, e, 'coordinates', self%coordinates)
        call json_get_required(json, e, 'geometry', self%geometry)
        self%ndof = size(self%geometry)
        call json_get_required(json, e, 'energy', self%energy, n)
        call json_get_required(json, e, 'grad', self%grad, [self%ndof, n])
        call json_get_required(json, e, 'have_grad', flags, n)
        self%have_grad = (flags == 1)
        call json_get_required(json, e, 'nadv', self%nadv, [self%ndof, n, n])
        call json_get_required(json, e, 'have_nadv', flags, n * n)
        self%have_nadv = reshape(flags == 1, [n, n])
        call json_get_required(json, e, 'wf_format', self%wf_format)
        call json_get_required(json, e, 'wf', self%wf)
        if (size(self%wf, 2) /= n) then
            write(stderr, *) 'Error in external_json_mod, load_state subroutine.'
            write(stderr, *) '  Wrong number of states in the saved wave functions.'
            stop
        end if
    end subroutine external_load_state


end module external_json_mod
