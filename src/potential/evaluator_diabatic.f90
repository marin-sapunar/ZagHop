module evaluator_diabatic_mod
    use global_defs
    use evaluator_base_mod
    implicit none


    type, abstract, extends(potential_evaluator) :: diabatic_evaluator
        integer :: n_group = 1
        integer, allocatable :: group_nstate(:)
        integer, allocatable :: group_i0(:)
        real(dp), allocatable :: diab_w(:, :)
        real(dp), allocatable :: diab_dw(:, :, :)
        real(dp), allocatable :: adiab_w(:, :)
        real(dp), allocatable :: adiab_dw(:, :, :)
        real(dp), allocatable :: adiab_trans(:, :)
        real(dp), allocatable :: group_adiab_w(:, :)
        real(dp), allocatable :: group_adiab_dw(:, :, :)
        real(dp), allocatable :: group_adiab_trans(:, :)
        real(dp), allocatable :: ldiab_trans(:, :)
    contains
        procedure :: eval => eval_all
        procedure(eval_diab), deferred :: eval_diab
        procedure :: change_basis => change_basis_mat
        procedure :: eval_eigvec => eval_eigvec
        procedure :: get_Hamiltonian => get_Hamiltonian
        procedure :: get_transformation => get_transformation
        procedure :: get_energy_single => get_energy_single
        procedure :: get_energy_all => get_energy_all
        procedure :: get_gradient => get_gradient
    end type diabatic_evaluator


    abstract interface
        subroutine eval_diab(self)
            import diabatic_evaluator
            class(diabatic_evaluator), intent(inout) :: self
        end subroutine eval_diab
    end interface


contains

    subroutine eval_all(self)
        class(diabatic_evaluator), intent(inout) :: self
        integer :: imode

        call self%eval_diab()
        call self%eval_eigvec('adiabatic')
        call self%eval_eigvec('group_adiabatic')

        if (allocated(self%ldiab_trans)) then
            call match_phase(self%ldiab_trans, self%adiab_trans)
        end if

        self%adiab_w = self%diab_w
        self%adiab_dw = self%diab_dw
        call self%change_basis('adiabatic', self%adiab_w)
        do imode = 1, size(self%diab_dw, 1)
            call self%change_basis('adiabatic', self%adiab_dw(imode, :, :))
        end do

        self%group_adiab_w = self%diab_w
        self%group_adiab_dw = self%diab_dw
        call self%change_basis('group_adiabatic', self%group_adiab_w)
        do imode = 1, size(self%diab_dw, 1)
            call self%change_basis('group_adiabatic', self%group_adiab_dw(imode, :, :))
        end do
    end subroutine eval_all

    function get_Hamiltonian(self, basis) result(H)
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        real(dp), allocatable :: H(:, :)

        select case(basis)
        case('diabatic')
            H = self%diab_w
        case('adiabatic')
            H = self%adiab_w
        case('group_adiabatic')
            H = self%group_adiab_w
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_Hamiltonian function.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select
    end function get_Hamiltonian


    function get_transformation(self, from_basis, to_basis) result(trans)
        use matrix_mod, only : unit_mat
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: from_basis, to_basis
        real(dp), allocatable :: trans(:, :)

        if (from_basis == to_basis) then
            trans = unit_mat(self%n_state)
            return
        end if

        select case(from_basis)
        case('diabatic')
            trans = unit_mat(self%n_state)
        case('adiabatic')
            trans = transpose(self%adiab_trans)
        case('group_adiabatic')
            trans = transpose(self%group_adiab_trans)
        case('locally_diabatic')
            trans = transpose(self%ldiab_trans)
        end select

        select case(to_basis)
        case('diabatic')
            ! No change of basis needed for diabatic representation
        case('adiabatic')
            trans = matmul(trans, self%adiab_trans)
        case('group_adiabatic')
            trans = matmul(trans, self%group_adiab_trans)
        case('locally_diabatic')
            trans = matmul(trans, self%ldiab_trans)
        end select
    end function get_transformation


    function get_gradient(self, basis, istate) result(grad)
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: istate
        real(dp), allocatable :: grad(:, :)

        allocate(grad(1, size(self%diab_dw, 1)))

        select case(basis)
        case('diabatic')
            grad(1, :) = self%diab_dw(:, istate, istate)
        case('adiabatic')
            grad(1, :) = self%adiab_dw(:, istate, istate)
        case('group_adiabatic')
            grad(1, :) = self%group_adiab_dw(:, istate, istate)
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_gradient function.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select
    end function get_gradient


    pure function get_energy_single(self, state) result(energy)
        class(diabatic_evaluator), intent(in) :: self
        integer, intent(in) :: state
        real(dp) :: energy

        energy = self%adiab_w(state, state)
    end function get_energy_single

    pure function get_energy_all(self) result(energy)
        use matrix_mod, only : diag
        class(diabatic_evaluator), intent(in) :: self
        real(dp), allocatable :: energy(:)

        energy = diag(self%adiab_w)
    end function get_energy_all


    subroutine eval_eigvec(self, basis)
        class(diabatic_evaluator), intent(inout) :: self
        character(len=*), intent(in) :: basis
        integer :: i, i0, ns

        select case(basis)
        case('diabatic')
            ! No change of basis needed for diabatic representation
        case('adiabatic')
            self%adiab_trans = get_eigvec(self%diab_w)
        case('group_adiabatic')
            self%group_adiab_trans = 0.0_dp
            do i = 1, self%n_group
                i0 = self%group_i0(i)
                ns = self%group_nstate(i)
                self%group_adiab_trans(i0:(i0 + ns - 1), i0:(i0 + ns - 1)) = &
                &   get_eigvec(self%diab_w(i0:(i0 + ns - 1), i0:(i0 + ns - 1)))
            end do
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, eval_eigvec subroutine.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select
    end subroutine eval_eigvec


    subroutine match_phase(bra, ket)
        real(dp), intent(in) :: bra(:, :)
        real(dp), intent(inout) :: ket(:, :)
        integer :: i
        real(dp), allocatable :: overlap(:, :)
        integer, allocatable :: phase(:)

        overlap = matmul(transpose(bra), ket)
        allocate(phase(size(ket, 2)))
        call phasematch_assigned_rotation(overlap, phase)
        do i = 1, size(ket, 2)
            if (phase(i) == -1) then
                ket(:, i) = -ket(:, i)
            end if
        end do
    end subroutine match_phase


    function get_eigvec(mat) result(eigv)
        use linalg_wrapper_mod, only : syev
        real(dp), intent(in) :: mat(:, :)
        real(dp), allocatable :: eigv(:, :)
        real(dp), allocatable :: wrk_e(:)

        eigv = mat
        allocate(wrk_e(size(eigv, 1)))
        call syev(eigv, wrk_e, jobz='V', uplo='U')
    end function get_eigvec


    subroutine change_basis_mat(self, basis, mat)
        use matrix_mod, only : unitary_transform
        class(diabatic_evaluator) :: self
        character(len=*), intent(in) :: basis
        real(dp), intent(inout) :: mat(:, :)

        select case(basis)
        case('diabatic')
            ! No change of basis needed for diabatic representation
        case('adiabatic')
            call unitary_transform(self%adiab_trans, mat, trans=.true.)
        case('group_adiabatic')
            call unitary_transform(self%group_adiab_trans, mat, trans=.true.)
        case('locally_diabatic')
            call unitary_transform(self%ldiab_trans, mat, trans=.true.)
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, change_basis_mat subroutine.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select        
    end subroutine change_basis_mat



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




end module evaluator_diabatic_mod