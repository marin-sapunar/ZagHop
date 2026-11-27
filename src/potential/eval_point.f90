!--------------------------------------------------------------------------------------------------
! PROGRAM: eval_point
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2026
!
! DESCRIPTION:
!> @brief Evaluate a model potential at one point, for a request of the external_json evaluator.
!> @details
!! Usage: eval_point [--help] <potential> [args]
!! - vibronic_coupling [template] - Vibronic coupling model (default template LVC.template).
!! - model <name> [parameter value ...] - One-dimensional model (see model_1d_mod).
!!
!! Reads the request (qm_in.json) and writes the response (qm_out.json) in the working directory.
!!
!! Extra properties ("properties" of the request):
!! - diabatic_transformation - The adiabatic states in the diabatic states, t[diabatic][adiabatic].
!! - diabatic_hamiltonian - The Hamiltonian of the model in its diabatic states.
!--------------------------------------------------------------------------------------------------
program eval_point
    use global_defs
    use string_mod, only : string
    use evaluator_diabatic_mod, only : diabatic_evaluator
    use vc_evaluator_mod, only : vc_evaluator
    use vc_model_mod, only : vc_model
    use model_1d_mod, only : model_1d_evaluator
    use coord_transform_mod, only : normal_modes
    use external_json_mod, only : request_file, response_file
    use json_module, only : json_core, json_value
    use json_util_mod, only : json_add_array, json_get_required, json_get_child, json_check_value
    implicit none

    class(diabatic_evaluator), allocatable :: pot
    type(vc_model), target :: model
    type(normal_modes), target :: nm
    type(string), allocatable :: args(:)
    integer :: narg, i

    ! Request.
    type(json_core) :: json
    type(json_value), pointer :: req, res, arr, item, obj
    integer :: point_id, mult, nstate, state
    character(len=:), allocatable :: coordinates, str
    real(dp), allocatable :: geom(:, :) !< Geometry (ndim, natom).
    integer, allocatable :: grad_states(:)
    integer, allocatable :: nacv_pairs(:, :)
    integer, allocatable :: pair(:)
    logical :: want_wf
    type(string), allocatable :: properties(:)

    ! Command line.
    narg = command_argument_count()
    if (narg == 0 .or. get_arg(1) == '--help') call usage()

    allocate(args(narg))
    do i = 1, narg
        args(i)%s = get_arg(i)
    end do

    select case(args(1)%s)
    case('vibronic_coupling')
        if (size(args) > 1) then
            call model%init(args(2)%s, nm)
        else
            call model%init('LVC.template', nm)
        end if
        allocate(vc_evaluator :: pot)
        select type(pot)
        type is (vc_evaluator)
            call pot%initialize(model)
        end select
    case('model')
        if (size(args) < 2) call usage()
        allocate(model_1d_evaluator :: pot)
        select type(pot)
        type is (model_1d_evaluator)
            call pot%init(args(2:))
        end select
    case default
        call usage()
    end select

    ! Read the request.
    call json%initialize()
    call json%load(file=request_file, p=req)
    if (json%failed()) then
        write(stderr, '(a)') 'Error in eval_point: could not read '//request_file//'.'
        call json%print_error_message(stderr)
        stop 1
    end if
    call json_check_value(json, req, 'protocol', 'zaghop')
    call json_check_value(json, req, 'version', 1)
    call json_get_required(json, req, 'point_id', point_id)
    call json_get_required(json, req, 'coordinates', coordinates)
    call json_get_required(json, req, 'geometry', geom)
    geom = transpose(geom)

    arr => json_get_child(json, req, 'states')
    nstate = 0
    do i = 1, json%count(arr)
        call json%get_child(arr, i, item)
        call json_get_required(json, item, 'multiplicity', mult)
        if (mult /= 1) call fail('Only singlet states are available.')
        call json_get_required(json, item, 'nstate', nstate)
    end do
    if (nstate < 1) call fail('No states requested.')
    if (any(pot%model_states%manifold%mult /= 1)) then
        call pot%select_states([nstate], [1])
    else if (nstate > pot%model_states%n_state) then
        call fail('More states requested than the model has.')
    else if (nstate < pot%model_states%n_state) then
        call pot%select_states([nstate], [1])
    end if

    allocate(grad_states(0))
    arr => json_get_child(json, req, 'request.gradients')
    do i = 1, json%count(arr)
        call json%get_child(arr, i, item)
        call json_get_required(json, item, 'multiplicity', mult)
        call json_get_required(json, item, 'state', state)
        if (mult /= 1 .or. state < 1 .or. state > nstate) call fail('Invalid gradient request.')
        grad_states = [grad_states, state]
    end do
    allocate(nacv_pairs(2, 0))
    arr => json_get_child(json, req, 'request.nacv')
    do i = 1, json%count(arr)
        call json%get_child(arr, i, item)
        call json_get_required(json, item, 'multiplicity', mult)
        call json_get_required(json, item, 'states', pair, 2)
        if (mult /= 1 .or. any(pair < 1) .or. any(pair > nstate) .or. pair(1) == pair(2)) then
            call fail('Invalid nacv request.')
        end if
        nacv_pairs = reshape([nacv_pairs, pair], [2, size(nacv_pairs, 2) + 1])
    end do
    call json_get_required(json, req, 'request.wave_function', want_wf)
    arr => json_get_child(json, req, 'request.properties')
    allocate(properties(json%count(arr)))
    do i = 1, size(properties)
        call json%get_child(arr, i, item)
        call json%get(item, str)
        if (json%failed()) call fail('Properties must be strings.')
        properties(i)%s = str
    end do
    call json%destroy(req)

    ! Evaluate.
    call pot%update_geometry(geom, coordinates)
    call pot%eval()

    ! Write the response.
    call json%create_object(res, '')
    call json%add(res, 'protocol', 'zaghop')
    call json%add(res, 'version', 1)
    call json%add(res, 'point_id', point_id)

    call json%create_array(arr, 'states')
    call json%create_object(obj, '')
    call json%add(obj, 'multiplicity', 1)
    call json%add(obj, 'energies', pot%get_energy('adiabatic'))
    call json%add(arr, obj)
    call json%add(res, arr)

    call json%create_array(arr, 'gradients')
    do i = 1, size(grad_states)
        call json%create_object(obj, '')
        call json%add(obj, 'multiplicity', 1)
        call json%add(obj, 'state', grad_states(i))
        call json_add_array(json, obj, 'gradient',                                                &
        &   transpose(reshape(pot%get_gradient('adiabatic', grad_states(i)), shape(geom))))
        call json%add(arr, obj)
    end do
    call json%add(res, arr)

    call json%create_array(arr, 'nacv')
    do i = 1, size(nacv_pairs, 2)
        call json%create_object(obj, '')
        call json%add(obj, 'multiplicity', 1)
        call json%add(obj, 'states', nacv_pairs(:, i))
        call json_add_array(json, obj, 'vector', transpose(reshape(                               &
        &   pot%get_nadv('adiabatic', nacv_pairs(1, i), nacv_pairs(2, i)), shape(geom))))
        call json%add(arr, obj)
    end do
    call json%add(res, arr)

    if (want_wf) then
        call json%create_array(arr, 'wave_function')
        call json%create_object(obj, '')
        call json%add(obj, 'multiplicity', 1)
        call json%add(obj, 'format', 'vectors')
        call json_add_array(json, obj, 'c', pot%adiab_trans)
        call json%add(arr, obj)
        call json%add(res, arr)
    end if

    call json%create_object(obj, 'properties')
    do i = 1, size(properties)
        select case(properties(i)%s)
        case('oscillator_strengths')
            select type(pot)
            type is (vc_evaluator)
                if (pot%model%dm(1)%max_order < 0) then
                    call fail('The model has no dipole moments for oscillator strengths.')
                end if
                call json%add(obj, 'oscillator_strengths', pot%get_oscill('adiabatic', 1))
            class default
                call fail('Oscillator strengths are only available for vibronic coupling models.')
            end select
        case('diabatic_transformation')
            call json_add_array(json, obj, 'diabatic_transformation', pot%adiab_trans)
        case('diabatic_hamiltonian')
            call json_add_array(json, obj, 'diabatic_hamiltonian', pot%diab_w)
        case default
            call fail('Unknown property: '//properties(i)%s)
        end select
    end do
    call json%add(res, obj)

    call json%print(res, response_file)
    call json%destroy(res)
    if (json%failed()) then
        write(stderr, '(a)') 'Error in eval_point: could not write '//response_file//'.'
        call json%print_error_message(stderr)
        stop 1
    end if

contains

    function get_arg(i) result(arg)
        integer, intent(in) :: i
        character(len=:), allocatable :: arg
        integer :: n

        call get_command_argument(i, length=n)
        allocate(character(len=n) :: arg)
        call get_command_argument(i, arg)
    end function get_arg


    subroutine usage()
        write(stderr, '(a)') 'Usage: eval_point [--help] <potential> [args]'
        write(stderr, '(a)') '  vibronic_coupling [template]'
        write(stderr, '(a)') '  model <name> [parameter value ...]'
        write(stderr, '(a)') 'Reads '//request_file//' and writes '//response_file//'.'
        stop 1
    end subroutine usage


    subroutine fail(msg)
        character(len=*), intent(in) :: msg

        write(stderr, '(a)') 'Error in eval_point: '//msg
        stop 1
    end subroutine fail

end program eval_point
