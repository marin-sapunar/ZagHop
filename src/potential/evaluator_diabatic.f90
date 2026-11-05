module evaluator_diabatic_mod
    use global_defs
    use evaluator_base_mod
    implicit none


    !> Potential given by a diabatic Hamiltonian.
    !> @details
    !! eval_diab gives the spin-free Hamiltonian (sf_w, sf_dw), once for each spin-free state, and
    !! the spin-orbit couplings between the states (soc_w, soc_dw). The full diabatic Hamiltonian
    !! is assembled from the two. The 'spin_free' basis has one state for each spin-free state.
    !! It diagonalises the spin-free Hamiltonian within each manifold, and its energies, gradients
    !! and couplings come from the spin-free Hamiltonian alone. The spin-orbit couplings between
    !! spin-free states are only given as squared magnitudes, summed over the components.
    type, abstract, extends(potential_evaluator) :: diabatic_evaluator
        integer :: ndof = 0
        real(dp), allocatable :: sf_w(:, :) !< Spin-free Hamiltonian (n_sf, n_sf).
        real(dp), allocatable :: sf_dw(:, :, :) !< Derivatives of sf_w (ndof, n_sf, n_sf).
        real(dp), allocatable :: soc_w(:, :) !< Spin-orbit couplings (n_state, n_state).
        real(dp), allocatable :: soc_dw(:, :, :) !< Derivatives of soc_w (ndof, n_state, n_state).
        real(dp), allocatable :: diab_w(:, :)
        real(dp), allocatable :: diab_dw(:, :, :)
        real(dp), allocatable :: adiab_w(:, :)
        real(dp), allocatable :: adiab_dw(:, :, :)
        real(dp), allocatable :: adiab_trans(:, :)
        real(dp), allocatable :: sf_trans(:, :) !< Eigenvectors of each manifold block of sf_w.
        real(dp), allocatable :: lsf_trans(:, :) !< sf_trans at the previous point.
        real(dp), allocatable :: spin_free_w(:, :) !< Spin-free Hamiltonian in the 'spin_free'
                                                   !! basis (n_sf, n_sf).
        real(dp), allocatable :: spin_free_dw(:, :, :) !< Derivatives of spin_free_w.
        real(dp), allocatable :: ldiab_trans(:, :)
    contains
        procedure :: eval => eval_all
        procedure(eval_diab), deferred :: eval_diab
        procedure :: allocate_hamiltonian
        procedure :: split_hamiltonian
        procedure :: change_basis => change_basis_mat
        procedure :: eval_eigvec => eval_eigvec
        procedure :: get_Hamiltonian => get_Hamiltonian
        procedure :: get_transformation => get_transformation
        procedure :: get_energy_single => get_energy_single
        procedure :: get_energy_all => get_energy_all
        procedure :: get_gradient => get_gradient
        procedure :: get_nadv => get_nadv
        procedure :: get_soc2 => get_soc2
    end type diabatic_evaluator


    abstract interface
        subroutine eval_diab(self)
            import diabatic_evaluator
            class(diabatic_evaluator), intent(inout) :: self
        end subroutine eval_diab
    end interface


contains

    subroutine eval_all(self)
        use matrix_mod, only : unitary_transform
        class(diabatic_evaluator), intent(inout) :: self
        integer :: imode, m, s0, n
        real(dp), allocatable :: w(:, :)

        call self%allocate_hamiltonian()
        call self%eval_diab()

        ! Full diabatic Hamiltonian.
        self%diab_w = self%expand(self%sf_w) + self%soc_w
        do imode = 1, self%ndof
            self%diab_dw(imode, :, :) = self%expand(self%sf_dw(imode, :, :))                       &
            &                         + self%soc_dw(imode, :, :)
        end do

        call self%eval_eigvec('adiabatic')
        call self%eval_eigvec('spin_free')

        if (allocated(self%ldiab_trans)) then
            call match_phase(self%ldiab_trans, self%adiab_trans)
        else
            ! First point, without states to match the phase to.
            call fix_sign(self%adiab_trans)
        end if
        if (allocated(self%lsf_trans)) then
            do m = 1, size(self%manifold)
                s0 = self%manifold(m)%sf0
                n = self%manifold(m)%n
                call match_phase(self%lsf_trans(s0:s0+n-1, s0:s0+n-1),                            &
                &                self%sf_trans(s0:s0+n-1, s0:s0+n-1))
            end do
        else
            call fix_sign(self%sf_trans)
        end if

        self%adiab_w = self%diab_w
        self%adiab_dw = self%diab_dw
        call self%change_basis('adiabatic', self%adiab_w)
        do imode = 1, size(self%diab_dw, 1)
            call self%change_basis('adiabatic', self%adiab_dw(imode, :, :))
        end do

        ! The spin-free quantities come from the spin-free Hamiltonian alone.
        self%spin_free_w = self%sf_w
        call unitary_transform(self%sf_trans, self%spin_free_w, trans=.true.)
        if (.not. allocated(self%spin_free_dw)) then
            allocate(self%spin_free_dw(self%ndof, self%n_sf, self%n_sf))
        end if
        do imode = 1, self%ndof
            w = self%sf_dw(imode, :, :)
            call unitary_transform(self%sf_trans, w, trans=.true.)
            self%spin_free_dw(imode, :, :) = w
        end do
    end subroutine eval_all


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: allocate_hamiltonian
    !> @brief Allocate the spin-free Hamiltonian, the spin-orbit couplings and the full diabatic
    !! Hamiltonian, if they are not allocated yet.
    !> @details
    !! The spin manifolds must be set. Spin-orbit couplings are zero unless eval_diab sets them.
    !----------------------------------------------------------------------------------------------
    subroutine allocate_hamiltonian(self)
        class(diabatic_evaluator), intent(inout) :: self
        integer :: ns, nsf, nd

        if (.not. allocated(self%manifold)) then
            write(stderr, *) 'Error in evaluator_diabatic_mod, allocate_hamiltonian subroutine.'
            write(stderr, *) '  Spin manifolds are not set.'
            stop
        end if
        ns = self%n_state
        nsf = self%n_sf
        nd = self%ndof
        if (.not. allocated(self%sf_w)) allocate(self%sf_w(nsf, nsf), source=0.0_dp)
        if (.not. allocated(self%sf_dw)) allocate(self%sf_dw(nd, nsf, nsf), source=0.0_dp)
        if (.not. allocated(self%soc_w)) allocate(self%soc_w(ns, ns), source=0.0_dp)
        if (.not. allocated(self%soc_dw)) allocate(self%soc_dw(nd, ns, ns), source=0.0_dp)
        if (.not. allocated(self%diab_w)) allocate(self%diab_w(ns, ns))
        if (.not. allocated(self%diab_dw)) allocate(self%diab_dw(nd, ns, ns))
    end subroutine allocate_hamiltonian


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: split_hamiltonian
    !> @brief Set the spin-free Hamiltonian and the spin-orbit couplings from a full Hamiltonian.
    !> @details
    !! The states of the full Hamiltonian are in the order of the spin manifolds. Within a group,
    !! the spin-orbit couplings sum to zero over the components of a manifold (the trace of a spin
    !! matrix is zero), so the spin-free Hamiltonian is the average of the group blocks over the
    !! components. It is computed as block_1 + sum(block_c - block_1) / n_comp, which is exactly
    !! block_1 when all blocks are equal. The rest of the full Hamiltonian are the spin-orbit
    !! couplings. With effective states, the manifold blocks are the spin-free Hamiltonian, and
    !! the elements between manifolds are the spin-orbit couplings.
    !----------------------------------------------------------------------------------------------
    subroutine split_hamiltonian(self, w, dw)
        class(diabatic_evaluator), intent(inout) :: self
        real(dp), intent(in) :: w(:, :) !< Full Hamiltonian.
        real(dp), intent(in) :: dw(:, :, :) !< Derivatives of the full Hamiltonian.
        integer :: imode

        call self%allocate_hamiltonian()
        self%sf_w = component_average(w)
        self%soc_w = w - self%expand(self%sf_w)
        do imode = 1, size(dw, 1)
            self%sf_dw(imode, :, :) = component_average(dw(imode, :, :))
            self%soc_dw(imode, :, :) = dw(imode, :, :) - self%expand(self%sf_dw(imode, :, :))
        end do
        if (self%mean_couplings) then
            call self%to_rss_couplings(self%soc_w)
            do imode = 1, size(dw, 1)
                call self%to_rss_couplings(self%soc_dw(imode, :, :))
            end do
        end if

    contains

        function component_average(mat) result(sf_mat)
            real(dp), intent(in) :: mat(:, :)
            real(dp) :: sf_mat(self%n_sf, self%n_sf)
            real(dp), allocatable :: delta(:, :)
            integer :: m, c, i0, s0, n

            sf_mat = 0.0_dp
            do m = 1, size(self%manifold)
                n = self%manifold(m)%n
                s0 = self%manifold(m)%sf0
                i0 = self%manifold(m)%i0
                allocate(delta(n, n), source=0.0_dp)
                do c = 2, self%manifold(m)%n_comp
                    delta = delta + (mat(i0+(c-1)*n:i0+c*n-1, i0+(c-1)*n:i0+c*n-1)               &
                    &                - mat(i0:i0+n-1, i0:i0+n-1))
                end do
                sf_mat(s0:s0+n-1, s0:s0+n-1) = mat(i0:i0+n-1, i0:i0+n-1)                          &
                &                            + delta / self%manifold(m)%n_comp
                deallocate(delta)
            end do
        end function component_average

    end subroutine split_hamiltonian

    function get_Hamiltonian(self, basis) result(H)
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        real(dp), allocatable :: H(:, :)

        select case(basis)
        case('diabatic')
            H = self%diab_w
        case('adiabatic')
            H = self%adiab_w
        case('spin_free')
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_Hamiltonian function.'
            write(stderr, *) '  No Hamiltonian in the spin_free basis. The spin-orbit couplings ', &
            &                'between spin-free states are only given by get_soc2.'
            stop
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_Hamiltonian function.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select
    end function get_Hamiltonian


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: get_transformation
    !> @brief Matrix whose columns are the states of to_basis in the states of from_basis.
    !> @details
    !! The 'spin_free' basis has one state for each spin-free state. It is only given in the
    !! spin-free diabatic states (from 'diabatic', with n_sf rows), and back.
    !----------------------------------------------------------------------------------------------
    function get_transformation(self, from_basis, to_basis) result(trans)
        use matrix_mod, only : unit_mat
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: from_basis, to_basis
        real(dp), allocatable :: trans(:, :)

        if (from_basis == to_basis) then
            trans = unit_mat(self%basis_size(to_basis))
            return
        end if

        if (from_basis == 'spin_free' .or. to_basis == 'spin_free') then
            if (from_basis == 'diabatic') then
                trans = self%sf_trans
            else if (to_basis == 'diabatic') then
                trans = transpose(self%sf_trans)
            else
                write(stderr, *) 'Error in evaluator_diabatic_mod, get_transformation function.'
                write(stderr, *) '  The spin_free basis is only given in the diabatic basis.'
                write(stderr, *) '  Requested: ', from_basis, ' to ', to_basis
                stop
            end if
            return
        end if

        select case(from_basis)
        case('diabatic')
            trans = unit_mat(self%n_state)
        case('adiabatic')
            trans = transpose(self%adiab_trans)
        case('locally_diabatic')
            trans = transpose(self%ldiab_trans)
        end select

        select case(to_basis)
        case('diabatic')
            ! No change of basis needed for diabatic representation
        case('adiabatic')
            trans = matmul(trans, self%adiab_trans)
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
        case('spin_free')
            grad(1, :) = self%spin_free_dw(:, istate, istate)
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_gradient function.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select
    end function get_gradient


    pure function get_energy_single(self, basis, state) result(energy)
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: state
        real(dp) :: energy

        select case(basis)
        case('diabatic')
            energy = self%diab_w(state, state)
        case('adiabatic')
            energy = self%adiab_w(state, state)
        case('spin_free')
            energy = self%spin_free_w(state, state)
        case default
            error stop 'Error in evaluator_diabatic_mod, get_energy_single function. Unknown basis.'
        end select
    end function get_energy_single

    pure function get_energy_all(self, basis) result(energy)
        use matrix_mod, only : diag
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        real(dp), allocatable :: energy(:)

        select case(basis)
        case('diabatic')
            energy = diag(self%diab_w)
        case('adiabatic')
            energy = diag(self%adiab_w)
        case('spin_free')
            energy = diag(self%spin_free_w)
        case default
            error stop 'Error in evaluator_diabatic_mod, get_energy_all function. Unknown basis.'
        end select
    end function get_energy_all


    subroutine eval_eigvec(self, basis)
        class(diabatic_evaluator), intent(inout) :: self
        character(len=*), intent(in) :: basis
        integer :: m, s0, n

        select case(basis)
        case('diabatic')
            ! No change of basis needed for diabatic representation
        case('adiabatic')
            self%adiab_trans = get_eigvec(self%diab_w)
        case('spin_free')
            if (.not. allocated(self%sf_trans)) allocate(self%sf_trans(self%n_sf, self%n_sf))
            self%sf_trans = 0.0_dp
            do m = 1, size(self%manifold)
                s0 = self%manifold(m)%sf0
                n = self%manifold(m)%n
                self%sf_trans(s0:s0+n-1, s0:s0+n-1) = get_eigvec(self%sf_w(s0:s0+n-1, s0:s0+n-1))
            end do
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, eval_eigvec subroutine.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select
    end subroutine eval_eigvec


    function get_nadv(self, basis, istate1, istate2) result(nadv)
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: istate1, istate2
        real(dp), allocatable :: nadv(:)
        integer :: imode
        real(dp) :: edif
        real(dp), parameter :: tiny_hf = 1.0e-8_dp

        allocate(nadv(self%ndof), source=0.0_dp)
        select case(basis)
        case('diabatic')
            continue
        case('adiabatic')
            edif = self%adiab_w(istate2, istate2) - self%adiab_w(istate1, istate1)
            if (abs(edif) < tiny_hf) then
                if (stdp1) then
                    write(stderr, *) 'Warning in evaluator_diabatic_mod, get_nadv function.'
                    write(stderr, *) '  Near-degeneracy between states ', istate1, ' and ', istate2, '.'
                    write(stderr, *) '  Setting denominator to ', tiny_hf, ' Hartree.'
                end if
                edif = sign(tiny_hf, edif)
            end if
            do imode = 1, self%ndof
                nadv(imode) = self%adiab_dw(imode, istate1, istate2) / edif
            end do
        case('spin_free')
             edif = self%spin_free_w(istate2, istate2) - self%spin_free_w(istate1, istate1)
             if (abs(edif) < tiny_hf) then
                if (stdp1) then
                    write(stderr, *) 'Warning in evaluator_diabatic_mod, get_nadv function.'
                    write(stderr, *) '  Near-degeneracy between states ', istate1, ' and ', istate2, '.'
                    write(stderr, *) '  Setting denominator to ', tiny_hf, ' Hartree.'
                end if
                edif = sign(tiny_hf, edif)
            end if
            do imode = 1, self%ndof
                nadv(imode) = self%spin_free_dw(imode, istate1, istate2) / edif
            end do
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_nadv function.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select
    end function get_nadv


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: get_soc2
    !> @brief Squared magnitudes of the spin-orbit couplings between all states of a basis.
    !> @details
    !! In the 'spin_free' basis, element (a, b) is the sum of the squared couplings between all
    !! components of the spin-free states a and b. For two manifolds of components, this is the
    !! root-sum-square coupling, squared. It depends neither on the signs of sf_trans nor on the
    !! basis of the components.
    !----------------------------------------------------------------------------------------------
    function get_soc2(self, basis) result(s2)
        use matrix_mod, only : unitary_transform
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        real(dp), allocatable :: s2(:, :)
        real(dp), allocatable :: soc(:, :)

        select case(basis)
        case('diabatic')
            s2 = abs(self%soc_w)**2
        case('spin_free')
            ! The same vectors of sf_trans transform each component of a spin-free state.
            soc = self%soc_w
            call unitary_transform(self%expand(self%sf_trans), soc, trans=.true.)
            s2 = self%sum_components(abs(soc)**2)
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_soc2 function.'
            write(stderr, *) '  Unsupported basis: ', basis
            stop
        end select
    end function get_soc2


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


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: fix_sign
    !> @brief Choose the sign of each vector so that its component with the largest magnitude is
    !! positive.
    !> @details
    !! The signs of the eigenvectors returned by the eigensolver are arbitrary. They can change
    !! with rounding differences in the matrix, which depend on the compiler options and the
    !! linear algebra library.
    !----------------------------------------------------------------------------------------------
    subroutine fix_sign(vec)
        real(dp), intent(inout) :: vec(:, :)
        integer :: i

        do i = 1, size(vec, 2)
            if (vec(maxloc(abs(vec(:, i)), dim=1), i) < 0.0_dp) vec(:, i) = -vec(:, i)
        end do
    end subroutine fix_sign


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