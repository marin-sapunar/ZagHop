!--------------------------------------------------------------------------------------------------
! MODULE: state_overlap_mod
!> @brief Overlaps between two sets of states, stored with the ids of the sets.
!> @details
!! Each set of states (for example, the electronic states computed at one geometry) gets an id
!! from request_id, unique within the program. A state_overlap holds the overlaps between a set
!! of bra states and a set of ket states, s(i, j) = <bra_i|ket_j>, together with the ids of the
!! two sets and the basis of the states. Users check the ids and the basis (matches) before using
!! the overlaps, so overlaps computed for states other than the expected ones are never used.
!!
!! The overlaps refer to the states as they are when they are set, including their phases. A
!! change of the states afterwards, such as a change of sign, makes them invalid.
!!
!! t is the Lowdin (symmetric) orthogonalisation of the columns of s, the orthogonal matrix
!! closest to s. For orthonormal states spanning the same space, s is already orthogonal. For
!! states which span only a part of the space (a subset of the states of a model, or the states
!! computed on the fly), s is not orthogonal, and its smallest singular value measures how much
!! of the space of the bra states is left in the ket states.
!!
!! match_phase chooses the signs of a set of real states so that they follow those of another
!! set, from the overlaps between the two.
!--------------------------------------------------------------------------------------------------
module state_overlap_mod
    use global_defs
    implicit none

    private
    public :: state_overlap
    public :: request_id
    public :: last_requested_id
    public :: restore_last_id
    public :: match_phase
    public :: phasematch_assigned_rotation


    integer :: last_id = 0 !< Last id given by request_id.


    type :: state_overlap
        integer :: id_bra = 0 !< Id of the bra states. 0 if the overlaps are not set.
        integer :: id_ket = 0 !< Id of the ket states.
        character(len=:), allocatable :: basis !< Basis of both sets of states.
        real(dp), allocatable :: s(:, :) !< Overlaps <bra_i|ket_j>.
        real(dp), allocatable :: t(:, :) !< Lowdin orthogonalised s.
    contains
        procedure :: set
        procedure :: clear
        procedure :: valid
        procedure :: matches
        procedure :: save_state
        procedure :: load_state
    end type state_overlap


contains


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: request_id
    !> @brief New id for a set of states, unique within the program.
    !> @details
    !! All ids are given by this function, so changes to how they are given are made here only
    !! (e.g. restoring the counter on restart). The counter is incremented atomically with OpenMP,
    !! so ids stay unique when they are requested from several threads. Ids are only unique within
    !! one process.
    !----------------------------------------------------------------------------------------------
    function request_id() result(id)
        integer :: id

        !$omp atomic capture
        last_id = last_id + 1
        id = last_id
        !$omp end atomic
    end function request_id


    !> Last id given by request_id, saved in a backup.
    integer function last_requested_id() result(id)
        id = last_id
    end function last_requested_id


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: restore_last_id
    !> @brief Continue the ids after id, the last id given by request_id in a previous run.
    !> @details
    !! Called when a run continues from a backup, after all objects that request ids at setup
    !! were created, so that the restored ids stay unique and new ids follow the previous run.
    !----------------------------------------------------------------------------------------------
    subroutine restore_last_id(id)
        integer, intent(in) :: id

        last_id = id
    end subroutine restore_last_id


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: set
    !> @brief Set the overlaps between two sets of states and their Lowdin orthogonalisation.
    !> @details
    !! orthog_lowdin computes t from the singular value decomposition s = U Sigma V^t as U V^t,
    !! which equals s (s^t s)^(-1/2) for nonsingular s. The smallest singular value is returned in
    !! sigma_min.
    !----------------------------------------------------------------------------------------------
    subroutine set(self, id_bra, id_ket, basis, s, sigma_min)
        use orthog_mod, only : orthog_lowdin
        class(state_overlap), intent(inout) :: self
        integer, intent(in) :: id_bra !< Id of the bra states.
        integer, intent(in) :: id_ket !< Id of the ket states.
        character(len=*), intent(in) :: basis !< Basis of both sets of states.
        real(dp), intent(in) :: s(:, :) !< Overlaps <bra_i|ket_j>.
        real(dp), intent(out), optional :: sigma_min !< Smallest singular value of s.
        real(dp), allocatable :: sigma(:)

        self%id_bra = id_bra
        self%id_ket = id_ket
        self%basis = basis
        self%s = s
        self%t = s
        call orthog_lowdin(self%t, sigma)
        if (present(sigma_min)) sigma_min = minval(sigma)
    end subroutine set


    !> Remove the overlaps.
    subroutine clear(self)
        class(state_overlap), intent(inout) :: self

        self%id_bra = 0
        self%id_ket = 0
        if (allocated(self%basis)) deallocate(self%basis)
        if (allocated(self%s)) deallocate(self%s)
        if (allocated(self%t)) deallocate(self%t)
    end subroutine clear


    !> The overlaps are set.
    pure logical function valid(self)
        class(state_overlap), intent(in) :: self

        valid = self%id_bra > 0 .and. allocated(self%t)
    end function valid


    !> The overlaps are set, between the states with the given ids, in the given basis.
    pure logical function matches(self, id_bra, id_ket, basis)
        class(state_overlap), intent(in) :: self
        integer, intent(in) :: id_bra !< Id of the bra states.
        integer, intent(in) :: id_ket !< Id of the ket states.
        character(len=*), intent(in) :: basis !< Basis of both sets of states.

        matches = .false.
        if (.not. self%valid()) return
        matches = self%id_bra == id_bra .and. self%id_ket == id_ket .and. self%basis == basis
    end function matches


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: save_state
    !> @brief Add the overlaps to a JSON object. Without overlaps, only the ids (0) are added.
    !----------------------------------------------------------------------------------------------
    subroutine save_state(self, json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_add_array
        class(state_overlap), intent(in) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the overlaps are added.

        call json%add(obj, 'id_bra', self%id_bra)
        call json%add(obj, 'id_ket', self%id_ket)
        if (.not. self%valid()) return
        call json%add(obj, 'basis', self%basis)
        call json_add_array(json, obj, 's', self%s)
        call json_add_array(json, obj, 't', self%t)
    end subroutine save_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: load_state
    !> @brief Restore the overlaps saved by save_state.
    !> @details
    !! t is read as saved, not computed again from s, so it is the same to the last bit.
    !----------------------------------------------------------------------------------------------
    subroutine load_state(self, json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_get_required
        class(state_overlap), intent(inout) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object holding the overlaps.
        integer :: id_bra, id_ket

        call self%clear()
        call json_get_required(json, obj, 'id_bra', id_bra)
        call json_get_required(json, obj, 'id_ket', id_ket)
        if (id_bra <= 0) return
        self%id_bra = id_bra
        self%id_ket = id_ket
        call json_get_required(json, obj, 'basis', self%basis)
        call json_get_required(json, obj, 's', self%s)
        call json_get_required(json, obj, 't', self%t, shape(self%s))
    end subroutine load_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: match_phase
    !> @brief Match the phases of the vectors in ket to those in bra (phasematch_assigned_rotation).
    !> @details
    !! The overlaps bra^T ket between the vectors, after the phases of ket are matched, are
    !! returned in olap, and the sign applied to each vector of ket in phase.
    !----------------------------------------------------------------------------------------------
    subroutine match_phase(bra, ket, olap, phase)
        real(dp), intent(in) :: bra(:, :)
        real(dp), intent(inout) :: ket(:, :)
        real(dp), allocatable, intent(out), optional :: olap(:, :) !< bra^T ket.
        integer, allocatable, intent(out), optional :: phase(:) !< Sign applied to each vector.
        integer :: i
        real(dp), allocatable :: overlap(:, :)
        integer, allocatable :: ph(:)

        overlap = matmul(transpose(bra), ket)
        allocate(ph(size(ket, 2)))
        call phasematch_assigned_rotation(overlap, ph)
        do i = 1, size(ket, 2)
            if (ph(i) == -1) then
                ket(:, i) = -ket(:, i)
            end if
        end do
        if (present(olap)) call move_alloc(overlap, olap)
        if (present(phase)) call move_alloc(ph, phase)
    end subroutine match_phase


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: phasematch_assigned_rotation
    !
    !> @brief Match phase of assigned states while avoiding reflections.
    !> @details 
    !! Each adiabatic state in the current step is matched with a single state in the previous
    !! step and its' phase is matched with the previous step. To ensure that the wave functions
    !! can be smoothly varied from the previous to the current step, a check is performed to make
    !! sure that the overlap matrix has the (approximate) form of a rotation matrix by making sure
    !! that the determinant is positive (for a pure rotation, the determinant is exactly 1, for a
    !! combination of rotation and reflection, the determinant is -1).
    !
    !> @note Instead of only checking that the full matrix is a rotation, the subroutine ensures 
    !! that each set of states which exchanged positions forms a rotation. (For example, if two
    !! trivial crossings occur in a single step, the overall matrix is a rotation if all assigned 
    !! overlaps are positive. However, it makes more sense to make sure the submatrix for each pair
    !! is individually a rotation (each pair will have 1 positive and 1 negative phase).
    !----------------------------------------------------------------------------------------------
    subroutine phasematch_assigned_rotation(olap, phase)
        use assignment_problem_mod
        use matrix_mod, only: mat_ge_det
        use ivec_mod
        real(dp), intent(inout) :: olap(:, :)
        integer, intent(out) :: phase(:)
        integer, allocatable :: cmatch(:)
        integer, allocatable :: rmatch(:)
        type(ivec), allocatable :: loops(:)
        integer :: st, i

        phase = 1

        ! Match states of current and previous step.
        call assignment_problem(-olap**2, rmatch, cmatch)

        ! Match phase based on assignment.
        do st = 1, minval(shape(olap))
            if (olap(rmatch(st), cmatch(st)) < 0) then
                olap(:, cmatch(st)) = -olap(:, cmatch(st))
                phase(cmatch(st)) = -1
            end if
        end do

        ! Check which sets of states swapped positions.
        call getloops(cmatch, loops)
        do i = 1, size(loops)
            ! Phase of states that didn't change position is already matched.
            if (size(loops(i)%c) == 1) cycle
            ! For sets of states that interchanged positions, make sure that the subsection of the
            ! overlap matrix spaning these states is a rotation matrix.
            if (mat_ge_det(olap(loops(i)%c, loops(i)%c)) < 0.0_dp) then
                olap(:, loops(i)%c(1)) = -olap(:, loops(i)%c(1))
                phase(loops(i)%c(1)) = -phase(loops(i)%c(1))
            end if
        end do
    end subroutine phasematch_assigned_rotation


end module state_overlap_mod
