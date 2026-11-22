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
!! closest to s. For orthonormal states spanning the same space, s is already orthogonal.
!--------------------------------------------------------------------------------------------------
module state_overlap_mod
    use global_defs
    implicit none

    private
    public :: state_overlap
    public :: request_id


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


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: set
    !> @brief Set the overlaps between two sets of states and their Lowdin orthogonalisation.
    !> @details
    !! orthog_lowdin computes t from the singular value decomposition s = U Sigma V^t as U V^t,
    !! which equals s (s^t s)^(-1/2) for nonsingular s.
    !----------------------------------------------------------------------------------------------
    subroutine set(self, id_bra, id_ket, basis, s)
        use orthog_mod, only : orthog_lowdin
        class(state_overlap), intent(inout) :: self
        integer, intent(in) :: id_bra !< Id of the bra states.
        integer, intent(in) :: id_ket !< Id of the ket states.
        character(len=*), intent(in) :: basis !< Basis of both sets of states.
        real(dp), intent(in) :: s(:, :) !< Overlaps <bra_i|ket_j>.

        self%id_bra = id_bra
        self%id_ket = id_ket
        self%basis = basis
        self%s = s
        self%t = s
        call orthog_lowdin(self%t)
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


end module state_overlap_mod
