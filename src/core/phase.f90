!--------------------------------------------------------------------------------------------------
! MODULE: phase_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date December, 2017
!
!> @brief Modify elements of overlap matrix to keep phase unchanged between bra and ket states.
!--------------------------------------------------------------------------------------------------
module phase_mod
    use global_defs
    implicit none

    private
    public :: phasematch
    public :: phasematch_diagonal
    public :: phasematch_assigned
    public :: phaseMatch_assigned_rotation


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: PhaseMatch
    !
    ! DESCRIPTION:
    !> @brief Modify overlap matrix elements to keep phase of wfs in current and previous step.
    !> @details
    !! Controlled by ctrl%phaselvl, see input_mod manual and subroutine documentation for details.
    !----------------------------------------------------------------------------------------------
    subroutine phasematch()
        use system_var
        use control_var

        select case(ctrl%tdc_type)
        case(1)
            select case(ctrl%phaselvl)
            case(0)
                ! Do not change overlap matrix.
            case(1)
                call phasematch_diagonal(t(1)%olap)
            case(2)
                call phasematch_assigned_rotation(t(1)%olap)
            case default
                write(stderr, *) 'Error in PhaseMod module, OverlapPhase subroutine.'
                write(stderr, *) ' Unrecognized method.'
                stop
            end select
        case(2)
            !> @todo Do anything for analytic NAD vectors?
        end select
    end subroutine phasematch


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: phasematch_diagonal
    !
    !> @brief Match phase of adiabatic states.
    !> @details 
    !! Phase of adiabatic states is kept constant. Signs of the rows/columns of the overlap matrix
    !! are changed to ensure the diagonal elements (corresponding to the overlap between each
    !! adiabatic state with itself in the next step) are always positive.
    !! In case of trivial crossings some of the diagonal elements can be close to zero.
    !----------------------------------------------------------------------------------------------
    subroutine phasematch_diagonal(olap)
        real(dp), intent(inout) :: olap(:, :)
        integer :: st, n1, n2
        logical, allocatable, save :: swap(:)

        n1 = size(olap, 1)
        n2 = size(olap, 2)
        if (.not. allocated(swap)) then
            ! In first step, allocate swap array.
            allocate(swap(n1))
        else
            ! In subsequent steps, swap sign of row of overlap matrix corresponding to state whose
            ! sign was changed in the previous step.
            do st = 1, n1
                if (swap(st)) then
                    olap(st, :) = -olap(st, :)
                end if
            end do
        end if
        swap = .false. ! Reset swap array.

        do st = 1, min(n1, n2)
            if (olap(st, st) < 0) then ! Check if sign of diagonal element is negative.
                ! Change sign of overlap matrix column corresponding to wf with changed sign.
                olap(:, st) = -olap(:, st)
                ! Save index of changed column, in next step change sign of the corresponding row.
                swap(st) = .true.
            end if
        end do
    end subroutine phasematch_diagonal


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
    subroutine phasematch_assigned_rotation(olap)
        use assignment_problem_mod
        use matrix_mod, only: mat_ge_det
        use ivec_mod
        real(dp), intent(inout) :: olap(:, :)
        integer, allocatable :: cmatch(:)
        integer, allocatable :: rmatch(:)
        logical, allocatable, save :: swap(:)
        type(ivec), allocatable :: loops(:)
        integer :: st, i, n1

        n1 = size(olap, 1)
        if (.not. allocated(swap)) then
            ! In first step, allocate swap array.
            allocate(swap(n1))
        else
            ! In subsequent steps, swap sign of row of overlap matrix corresponding to state whose
            ! sign was changed in the previous step.
            do st = 1, n1
                if (swap(st)) then
                    olap(st, :) = -olap(st, :)
                end if
            end do
        end if
        swap = .false. ! Reset swap array.

        ! Match states of current and previous step.
        call assignment_problem(-olap**2, rmatch, cmatch)

        ! Match phase based on assignment.
        do st = 1, n1
            if (olap(st, cmatch(st)) < 0) then
                olap(:, cmatch(st)) = -olap(:, cmatch(st))
                swap(cmatch(st)) = .not. swap(cmatch(st))
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
                swap(loops(i)%c(1)) = .not. swap(loops(i)%c(1))
            end if
        end do
    end subroutine phasematch_assigned_rotation


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: phasematch_assigned
    !
    !> @brief Match phase of assigned states.
    !> @details 
    !! Phase of matched states is kept constant. Identities of the states are tracked by solving 
    !! the assignment problem for -abs(S). Signs of the rows/columns of the overlap matrix are 
    !! changed to ensure the elements corresponding to the overlaps of matched states are always
    !! positive.
    !
    !> @note
    !! The overlap matrix should be close to a rotation matrix. This version of the code doesn't
    !! ensure this. This is most important in the case of an unavoided crossing between two states
    !! where the subroutine gives both states a positive phase, which means one of the states is
    !! reflected instead of rotated.
    !----------------------------------------------------------------------------------------------
    subroutine phasematch_assigned(olap)
        use assignment_problem_mod
        real(dp), intent(inout) :: olap(:, :)
        integer, allocatable :: cmatch(:)
        integer, allocatable :: rmatch(:)
        integer :: st, n1
        logical, allocatable, save :: swap(:)

        n1 = size(olap, 1)
        if (.not. allocated(swap)) then
            ! In first step, allocate swap array.
            allocate(swap(n1))
        else
            ! In subsequent steps, swap sign of row of overlap matrix corresponding to state whose
            ! sign was changed in the previous step.
            do st = 1, n1
                if (swap(st)) then
                    olap(st, :) = -olap(st, :)
                end if
            end do
        end if
        swap = .false. ! Reset swap array.

        ! Match states of current and previous step.
        call assignment_problem(-olap**2, rmatch, cmatch)


        do st = 1, n1
            if (olap(st, cmatch(st)) < 0) then ! Check if sign of any wf changed between steps.
                ! Change sign of overlap matrix column corresponding to wf with changed sign.
                olap(:, cmatch(st)) = -olap(:, cmatch(st))
                ! Save index of changed column, in next step change sign of the corresponding row.
                swap(cmatch(st)) = .true.
            end if
        end do
    end subroutine phasematch_assigned


end module phase_mod
