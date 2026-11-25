!--------------------------------------------------------------------------------------------------
! MODULE: evaluator_diabatic_mod
!
! DESCRIPTION:
!> @brief Potentials defined by a Hamiltonian in a diabatic basis.
!> @details
!! The model gives the spin-free part of the Hamiltonian and the spin-orbit couplings in its
!! diabatic states. The names of the variables give their basis (see evaluator_base_mod):
!! - sf_, soc_ - The two parts of the Hamiltonian of the model in the diabatic basis.
!! - diab_ - The full Hamiltonian of the model in the diabatic basis.
!! - mch_ - MCH states of the model: eigenstates of the spin-free part within each multiplicity.
!! - adiab_ - Adiabatic states of the dynamics.
!!
!! The dynamics can run on a subset of the states of the model (select_states). With P the MCH
!! states of the model with all spin components, and sel the lowest MCH states of each
!! multiplicity, the adiabatic states are the eigenstates of the Hamiltonian in the selected MCH
!! states,
!!     H_sub = [P^T (H_sf + H_soc) P](sel, sel) = diag(e)(sel) + [P^T H_soc P](sel, sel),
!! with eigenvectors V and U = P(:, sel) V in the diabatic basis (adiab_trans). This is the
!! construction of on-the-fly dynamics, where only some MCH states are computed and the
!! spin-orbit couplings are diagonalised between them (Mai, Marquetand and Gonzalez, WIREs Comput.
!! Mol. Sci. 8, e1370 (2018), doi: 10.1002/wcms.1370, eqs. 7 and 8, where V is called U). With
!! all states of the model selected, the adiabatic states are the eigenstates of the full
!! Hamiltonian.
!!
!! U has a column for each state of the dynamics and a row for each state of the model. U^T U
!! is the unit matrix, and U U^T is the projector onto the selected MCH states. The overlaps of
!! the adiabatic states of two points are then not an orthogonal matrix.
!!
!! Derivatives of the adiabatic energies of a subset. The selected MCH states change with the
!! geometry, so the derivative of H_sub has a term from the derivative couplings A_x between
!! the unselected (uns) and the selected states,
!!     dH_sub/dx = [P^T (dH/dx) P](sel, sel) + C_x + C_x^T,
!!     C_x = [P^T H_soc P](sel, uns) A_x(uns, sel),
!!     A_x(u, s) = <u|d/dx s> = [P^T (dH_sf/dx) P](u, s) / (e(s) - e(u)).
!! The derivative couplings between selected states do not contribute. adiab_dw is V^T dH_sub/dx
!! V. Its diagonal gives the gradients, and its off-diagonal elements divided by the energy gaps
!! the nonadiabatic couplings within the states of the dynamics. Without exact_gradient, C_x is
!! left out, which is the Hellmann-Feynman expression U^T (dH/dx) U. In the MCH states this is
!! the gradient transformation of eq. 13 of the reference, with the derivative of the spin-orbit
!! couplings of the model added. (With K_mn = <m|d/dx n>, the second term of eq. 13 enters as
!! (E_n - E_m) K_mn, opposite to the sign printed in the reference.) C_x is zero when all states
!! are selected or when the selected and unselected states have no spin-orbit coupling.
!--------------------------------------------------------------------------------------------------
module evaluator_diabatic_mod
    use global_defs
    use evaluator_base_mod
    implicit none


    real(dp), parameter :: gap_tol = 1.0e-8_dp !< Smallest magnitude of an energy gap in the
                                               !! denominator of a derivative coupling (Eh).
    real(dp), parameter :: overlap_warn_tol = 0.9_dp !< Smallest singular value of the overlaps
    !! between the adiabatic states of two points without a warning.


    !----------------------------------------------------------------------------------------------
    !> TYPE: diabatic_evaluator
    !> @brief Evaluate potentials defined in the diabatic basis.
    !> @details
    !! The spin-free Hamiltonian (sf_w) and the spin-orbit couplings (soc_w) along with their
    !! derivatives (sf_dw and soc_dw) are given by eval_hamiltonian, for all states of the model.
    !! The full diabatic Hamiltonian is assembled from the two. The 'mch' basis has one state for
    !! each spin-free state of the dynamics. Its energies, gradients and couplings come from the
    !! spin-free Hamiltonian alone. The spin-orbit couplings between its states are only given as
    !! squared magnitudes, summed over the components.
    !----------------------------------------------------------------------------------------------
    type, abstract, extends(potential_evaluator) :: diabatic_evaluator
        integer :: ndof = 0
        real(dp), allocatable :: sf_w(:, :) !< Spin-free Hamiltonian of the model (n_sf, n_sf).
        real(dp), allocatable :: sf_dw(:, :, :) !< Derivatives of sf_w (ndof, n_sf, n_sf).
        real(dp), allocatable :: soc_w(:, :) !< Spin-orbit couplings of the model
                                             !! (n_state, n_state).
        real(dp), allocatable :: soc_dw(:, :, :) !< Derivatives of soc_w (ndof, n_state, n_state).
        real(dp), allocatable :: diab_w(:, :) !< Hamiltonian of the model.
        real(dp), allocatable :: diab_dw(:, :, :) !< Derivatives of diab_w.
        real(dp), allocatable :: mch_trans(:, :) !< MCH states of the model: eigenvectors of each
                                                 !! manifold block of sf_w (n_sf, n_sf).
        real(dp), allocatable :: mch_w(:, :) !< Spin-free Hamiltonian in the MCH states of the
                                             !! model (n_sf, n_sf).
        real(dp), allocatable :: mch_dw(:, :, :) !< Derivatives of sf_w in the MCH states of the
                                                 !! model.
        real(dp), allocatable :: adiab_trans(:, :) !< Adiabatic states of the dynamics in the
                                                   !! diabatic states (model states, states).
        real(dp), allocatable :: adiab_w(:, :) !< Hamiltonian in the adiabatic states.
        real(dp), allocatable :: adiab_dw(:, :, :) !< Derivatives of the Hamiltonian in the
        !! adiabatic states, with the term from the change of the selected MCH states.
        logical :: exact_gradient = .true. !< Include the term from the change of the selected MCH
        !! states in adiab_dw (see the module description).
        logical :: soc_checked = .false. !< Spin-orbit couplings within groups were checked.
    contains
        procedure :: eval => eval_all
        procedure(eval_diab), deferred :: eval_diab
        procedure :: eval_hamiltonian
        procedure :: eval_states
        procedure :: allocate_hamiltonian
        procedure :: allocate_bases
        procedure :: split_hamiltonian
        procedure :: set_full_hamiltonian
        procedure :: in_basis
        procedure :: mch_soc
        procedure :: get_Hamiltonian => get_Hamiltonian
        procedure :: get_transformation => get_transformation
        procedure :: get_energy_single => get_energy_single
        procedure :: get_energy_all => get_energy_all
        procedure :: get_gradient => get_gradient
        procedure :: get_nadv => get_nadv
        procedure :: get_soc2 => get_soc2
        procedure :: save_state => diabatic_save_state
        procedure :: load_state => diabatic_load_state
    end type diabatic_evaluator


    abstract interface
        !> Evaluate the potential at the current geometry. Fills sf_w and sf_dw, and soc_w and
        !! soc_dw if the potential has spin-orbit couplings, directly or with
        !! set_full_hamiltonian.
        subroutine eval_diab(self)
            import diabatic_evaluator
            class(diabatic_evaluator), intent(inout) :: self
        end subroutine eval_diab
    end interface


contains

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: eval_all
    !> @brief Evaluate the potential and build all bases at the current geometry.
    !> @details
    !! With prev, the evaluator of the previous point, the phases of the MCH states of the model
    !! and of the adiabatic states are matched to those of prev, and link holds the overlaps
    !! between the adiabatic states of prev and of this point, U(prev)^T U. Without prev, the
    !! signs of the states are fixed by fix_sign.
    !----------------------------------------------------------------------------------------------
    subroutine eval_all(self, prev)
        class(diabatic_evaluator), intent(inout) :: self
        class(potential_evaluator), intent(in), optional :: prev
        integer :: imode

        call self%new_point()
        call self%eval_hamiltonian()

        ! Full diabatic Hamiltonian.
        self%diab_w = self%model_states%expand(self%sf_w) + self%soc_w
        do imode = 1, self%ndof
            self%diab_dw(imode, :, :) = self%model_states%expand(self%sf_dw(imode, :, :))        &
            &                         + self%soc_dw(imode, :, :)
        end do

        call self%allocate_bases()
        if (present(prev)) then
            select type(prev)
            class is (diabatic_evaluator)
                call check_previous(prev)
                call self%eval_states(prev)
            class default
                write(stderr, *) 'Error in evaluator_diabatic_mod, eval_all subroutine.'
                write(stderr, *) '  The previous point is not given by a diabatic Hamiltonian.'
                stop
            end select
        else
            call self%eval_states()
        end if

    contains

        !> Check that prev was evaluated, with the same states.
        subroutine check_previous(prev)
            class(diabatic_evaluator), intent(in) :: prev

            if (.not. allocated(prev%adiab_trans) .or. .not. allocated(prev%mch_trans)) then
                write(stderr, *) 'Error in evaluator_diabatic_mod, eval_all subroutine.'
                write(stderr, *) '  The previous point was not evaluated.'
                stop
            end if
            if (any(shape(prev%adiab_trans) /= [self%model_states%n_state, self%states%n_state]) &
            &   .or. any(shape(prev%mch_trans) /= shape(self%mch_trans))) then
                write(stderr, *) 'Error in evaluator_diabatic_mod, eval_all subroutine.'
                write(stderr, *) '  The previous point has a different number of states.'
                stop
            end if
        end subroutine check_previous

    end subroutine eval_all


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: eval_states
    !> @brief Build the MCH states of the model and the adiabatic states of the dynamics from the
    !! diabatic Hamiltonian, with the Hamiltonian and its derivatives in both.
    !> @details
    !! See the module description for the adiabatic states of a subset of the states of the model
    !! and their derivatives. The spin-free part of H_sub is the diagonal matrix of the
    !! eigenvalues of the MCH states, so without spin-orbit couplings the adiabatic states are
    !! exactly the selected MCH states.
    !----------------------------------------------------------------------------------------------
    subroutine eval_states(self, prev)
        use matrix_mod, only : unitary_transform
        class(diabatic_evaluator), intent(inout) :: self
        class(diabatic_evaluator), intent(in), optional :: prev !< Evaluated previous point.
        integer :: imode, m, s0, n, i
        type(state_location) :: loc
        real(dp), allocatable :: w(:, :)
        real(dp), allocatable :: e(:) !< Energies of the MCH states of the model (n_sf).
        real(dp), allocatable :: e_c(:) !< The same for the MCH states with all components.
        real(dp), allocatable :: p(:, :) !< MCH states of the model with all components.
        real(dp), allocatable :: soc(:, :) !< Spin-orbit couplings in the states p.
        real(dp), allocatable :: h_sub(:, :) !< Hamiltonian in the selected MCH states.
        real(dp), allocatable :: v(:, :) !< Adiabatic states in the selected MCH states.
        real(dp), allocatable :: s(:, :) !< Overlaps between the adiabatic states of prev and self.
        real(dp) :: sigma_min !< Smallest singular value of s.

        ! MCH states of the model.
        allocate(e(self%model_states%n_sf))
        self%mch_trans = 0.0_dp
        do m = 1, size(self%model_states%manifold)
            s0 = self%model_states%manifold(m)%sf0
            n = self%model_states%manifold(m)%n
            self%mch_trans(s0:s0+n-1, s0:s0+n-1) = get_eigvec(self%sf_w(s0:s0+n-1, s0:s0+n-1),    &
            &                                                 e(s0:s0+n-1))
            if (present(prev)) then
                call match_phase(prev%mch_trans(s0:s0+n-1, s0:s0+n-1),                            &
                &                self%mch_trans(s0:s0+n-1, s0:s0+n-1))
            end if
        end do
        ! First point, without states to match the phase to.
        if (.not. present(prev)) call fix_sign(self%mch_trans)

        ! The spin-free quantities come from the spin-free Hamiltonian alone.
        self%mch_w = self%sf_w
        call unitary_transform(self%mch_trans, self%mch_w, trans=.true.)
        do imode = 1, self%ndof
            w = self%sf_dw(imode, :, :)
            call unitary_transform(self%mch_trans, w, trans=.true.)
            self%mch_dw(imode, :, :) = w
        end do

        ! Adiabatic states: eigenstates of the Hamiltonian in the selected MCH states.
        p = self%model_states%expand(self%mch_trans)
        soc = self%mch_soc()
        allocate(e_c(self%model_states%n_state))
        do i = 1, self%model_states%n_state
            loc = self%model_states%locate_state(i)
            e_c(i) = e(loc%sf)
        end do
        h_sub = soc(self%sel, self%sel)
        do i = 1, size(self%sel)
            h_sub(i, i) = h_sub(i, i) + e_c(self%sel(i))
        end do
        v = get_eigvec(h_sub)
        self%adiab_trans = matmul(p(:, self%sel), v)

        if (present(prev)) then
            call match_phase(prev%adiab_trans, self%adiab_trans, s)
            call self%link%set(prev%id, self%id, 'adiabatic', s, sigma_min)
            if (sigma_min < overlap_warn_tol) then
                write(stderr, *) 'Warning in evaluator_diabatic_mod, eval_states subroutine.'
                write(stderr, '(3x,a,f8.5,a)') 'Smallest singular value of the overlaps '//       &
                &   'between the adiabatic states of two points: ', sigma_min, '.'
                write(stderr, *) '  A state of the dynamics changed character, for example at '//  &
                &                'a crossing with a state which is not selected.'
            end if
        else
            call fix_sign(self%adiab_trans)
        end if

        self%adiab_w = self%in_basis('adiabatic', self%diab_w)
        do imode = 1, self%ndof
            self%adiab_dw(imode, :, :) = self%in_basis('adiabatic', self%diab_dw(imode, :, :))
        end do
        if (self%exact_gradient .and. size(self%uns) > 0) call add_subspace_term()

    contains

        !> Add the term from the change of the selected MCH states with the geometry,
        !! V^T (C_x + C_x^T) V, to adiab_dw.
        subroutine add_subspace_term()
            real(dp), allocatable :: soc_su(:, :) !< Couplings of selected and unselected states.
            real(dp), allocatable :: d(:, :) !< Derivative of sf_w in the states p.
            real(dp), allocatable :: a(:, :) !< Derivative couplings A_x(uns, sel).
            real(dp), allocatable :: c(:, :) !< C_x + C_x^T.
            integer :: is, iu
            logical :: clamped, any_clamped

            soc_su = soc(self%sel, self%uns)
            if (all(soc_su == 0.0_dp)) return
            ! V with the phases of adiab_trans.
            v = matmul(transpose(p(:, self%sel)), self%adiab_trans)
            allocate(a(size(self%uns), size(self%sel)))
            any_clamped = .false.
            do imode = 1, self%ndof
                d = self%model_states%expand(self%mch_dw(imode, :, :))
                a = 0.0_dp
                do is = 1, size(self%sel)
                    do iu = 1, size(self%uns)
                        if (d(self%uns(iu), self%sel(is)) == 0.0_dp) cycle
                        a(iu, is) = d(self%uns(iu), self%sel(is))                                 &
                        &   / clamped_gap(e_c(self%sel(is)) - e_c(self%uns(iu)), clamped)
                        any_clamped = any_clamped .or. clamped
                    end do
                end do
                c = matmul(soc_su, a)
                c = c + transpose(c)
                self%adiab_dw(imode, :, :) = self%adiab_dw(imode, :, :)                           &
                &                          + matmul(matmul(transpose(v), c), v)
            end do
            if (any_clamped .and. stdp1) then
                write(stderr, *) 'Warning in evaluator_diabatic_mod, eval_states subroutine.'
                write(stderr, *) '  Near-degeneracy between a selected and an unselected MCH state.'
                write(stderr, *) '  Setting denominator to ', gap_tol, ' Hartree.'
            end if
        end subroutine add_subspace_term

    end subroutine eval_states


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: allocate_hamiltonian
    !> @brief Allocate the spin-free Hamiltonian, the spin-orbit couplings and the full diabatic
    !! Hamiltonian of the model, if they are not allocated yet.
    !> @details
    !! The states of the model must be set. Spin-orbit couplings are zero unless eval_diab sets
    !! them.
    !----------------------------------------------------------------------------------------------
    subroutine allocate_hamiltonian(self)
        class(diabatic_evaluator), intent(inout) :: self
        integer :: ns, nsf, nd

        if (.not. self%model_states%is_set()) then
            write(stderr, *) 'Error in evaluator_diabatic_mod, allocate_hamiltonian subroutine.'
            write(stderr, *) '  Spin manifolds are not set.'
            stop
        end if
        ns = self%model_states%n_state
        nsf = self%model_states%n_sf
        nd = self%ndof
        if (.not. allocated(self%sf_w)) allocate(self%sf_w(nsf, nsf), source=0.0_dp)
        if (.not. allocated(self%sf_dw)) allocate(self%sf_dw(nd, nsf, nsf), source=0.0_dp)
        if (.not. allocated(self%soc_w)) allocate(self%soc_w(ns, ns), source=0.0_dp)
        if (.not. allocated(self%soc_dw)) allocate(self%soc_dw(nd, ns, ns), source=0.0_dp)
        if (.not. allocated(self%diab_w)) allocate(self%diab_w(ns, ns))
        if (.not. allocated(self%diab_dw)) allocate(self%diab_dw(nd, ns, ns))
    end subroutine allocate_hamiltonian


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: allocate_bases
    !> @brief Allocate the Hamiltonian and its derivatives in the adiabatic states of the dynamics
    !! and in the MCH states of the model, if they are not allocated yet.
    !----------------------------------------------------------------------------------------------
    subroutine allocate_bases(self)
        class(diabatic_evaluator), intent(inout) :: self
        integer :: ns, nsf, nd

        ns = self%states%n_state
        nsf = self%model_states%n_sf
        nd = self%ndof
        if (.not. allocated(self%adiab_w)) allocate(self%adiab_w(ns, ns))
        if (.not. allocated(self%adiab_dw)) allocate(self%adiab_dw(nd, ns, ns))
        if (.not. allocated(self%mch_trans)) allocate(self%mch_trans(nsf, nsf))
        if (.not. allocated(self%mch_w)) allocate(self%mch_w(nsf, nsf))
        if (.not. allocated(self%mch_dw)) allocate(self%mch_dw(nd, nsf, nsf))
    end subroutine allocate_bases


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: eval_hamiltonian
    !> @brief Evaluate the spin-free Hamiltonian and the spin-orbit couplings at the current
    !! geometry, without building the bases.
    !> @details
    !! At the first evaluation, warns about spin-orbit couplings within a group of states. These
    !! are zero for a real Hamiltonian, so they point to an input error, such as components which
    !! are not ordered as assumed.
    !----------------------------------------------------------------------------------------------
    subroutine eval_hamiltonian(self)
        class(diabatic_evaluator), intent(inout) :: self
        real(dp) :: smax

        call self%allocate_hamiltonian()
        call self%eval_diab()

        if (.not. self%soc_checked) then
            smax = self%model_states%max_soc_in_groups(self%soc_w)
            if (smax > soc_group_tol) then
                write(stderr, *) 'Warning in evaluator_diabatic_mod, eval_hamiltonian subroutine.'
                write(stderr, '(3x,a,e10.3,a)') 'Spin-orbit couplings within a group of states, ', &
                &                               'up to ', smax, ' Eh.'
                write(stderr, *) '  These are zero for a real Hamiltonian. Check the order of ',  &
                &                'the states and spin_states.'
            end if
            self%soc_checked = .true.
        end if
    end subroutine eval_hamiltonian


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: set_full_hamiltonian
    !> @brief Set the spin-free Hamiltonian and the spin-orbit couplings from a full Hamiltonian
    !! given in the order of the source of the potential (component_order).
    !----------------------------------------------------------------------------------------------
    subroutine set_full_hamiltonian(self, w, dw)
        class(diabatic_evaluator), intent(inout) :: self
        real(dp), intent(in) :: w(:, :) !< Full Hamiltonian, in the order of the source.
        real(dp), intent(in) :: dw(:, :, :) !< Derivatives of w, (ndof, n_state, n_state).

        call self%split_hamiltonian(w(self%state_map, self%state_map),                           &
        &                           dw(:, self%state_map, self%state_map))
    end subroutine set_full_hamiltonian


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
        self%soc_w = w - self%model_states%expand(self%sf_w)
        do imode = 1, size(dw, 1)
            self%sf_dw(imode, :, :) = component_average(dw(imode, :, :))
            self%soc_dw(imode, :, :) = dw(imode, :, :)                                            &
            &                        - self%model_states%expand(self%sf_dw(imode, :, :))
        end do
        if (self%model_states%spin_states == 'effective_mean') then
            call self%model_states%to_rss_couplings(self%soc_w)
            do imode = 1, size(dw, 1)
                call self%model_states%to_rss_couplings(self%soc_dw(imode, :, :))
            end do
        end if

    contains

        function component_average(mat) result(sf_mat)
            real(dp), intent(in) :: mat(:, :)
            real(dp) :: sf_mat(self%model_states%n_sf, self%model_states%n_sf)
            real(dp), allocatable :: delta(:, :)
            integer :: m, c, i0, s0, n

            sf_mat = 0.0_dp
            do m = 1, size(self%model_states%manifold)
                n = self%model_states%manifold(m)%n
                s0 = self%model_states%manifold(m)%sf0
                i0 = self%model_states%manifold(m)%i0
                allocate(delta(n, n), source=0.0_dp)
                do c = 2, self%model_states%manifold(m)%n_comp
                    delta = delta + (mat(i0+(c-1)*n:i0+c*n-1, i0+(c-1)*n:i0+c*n-1)               &
                    &                - mat(i0:i0+n-1, i0:i0+n-1))
                end do
                sf_mat(s0:s0+n-1, s0:s0+n-1) = mat(i0:i0+n-1, i0:i0+n-1)                          &
                &                            + delta / self%model_states%manifold(m)%n_comp
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
        case('mch')
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_Hamiltonian function.'
            write(stderr, *) '  No Hamiltonian in the mch basis. The spin-orbit couplings ', &
            &                'between its states are only given by get_soc2.'
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
    !! The diabatic basis has all states of the model, and the adiabatic and mch bases the states
    !! of the dynamics, so a transformation from or to the diabatic basis is not square when the
    !! dynamics runs on a subset of the states. The 'mch' basis has one state for each spin-free
    !! state. It is only given in the spin-free diabatic states (from 'diabatic', with a row for
    !! each spin-free state of the model), and back.
    !----------------------------------------------------------------------------------------------
    function get_transformation(self, from_basis, to_basis) result(trans)
        use matrix_mod, only : unit_mat
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: from_basis, to_basis
        real(dp), allocatable :: trans(:, :)

        select case(from_basis)
        case('diabatic', 'adiabatic', 'mch')
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_transformation function.'
            write(stderr, *) '  Unknown basis: ', from_basis
            stop
        end select

        if (from_basis == to_basis) then
            trans = unit_mat(self%basis_size(to_basis))
            return
        end if

        select case(from_basis//' '//to_basis)
        case('diabatic adiabatic')
            trans = self%adiab_trans
        case('adiabatic diabatic')
            trans = transpose(self%adiab_trans)
        case('diabatic mch')
            trans = self%mch_trans(:, self%sel_mch)
        case('mch diabatic')
            trans = transpose(self%mch_trans(:, self%sel_mch))
        case('adiabatic mch', 'mch adiabatic')
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_transformation function.'
            write(stderr, *) '  The mch basis is only given in the diabatic basis.'
            write(stderr, *) '  Requested: ', from_basis, ' to ', to_basis
            stop
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_transformation function.'
            write(stderr, *) '  Unknown basis: ', to_basis
            stop
        end select
    end function get_transformation


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: get_gradient
    !> @brief Gradient of a state, flat, in the coordinates of the last geometry passed to
    !! update_geometry.
    !----------------------------------------------------------------------------------------------
    function get_gradient(self, basis, istate) result(grad)
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: istate
        real(dp), allocatable :: grad(:)
        real(dp), allocatable :: d(:)

        select case(basis)
        case('diabatic')
            d = self%diab_dw(:, istate, istate)
        case('adiabatic')
            d = self%adiab_dw(:, istate, istate)
        case('mch')
            d = self%mch_dw(:, self%sel_mch(istate), self%sel_mch(istate))
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_gradient function.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select
        grad = self%get_derivative_target_coord(d)
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
        case('mch')
            energy = self%mch_w(self%sel_mch(state), self%sel_mch(state))
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
        case('mch')
            energy = diag(self%mch_w)
            energy = energy(self%sel_mch)
        case default
            error stop 'Error in evaluator_diabatic_mod, get_energy_all function. Unknown basis.'
        end select
    end function get_energy_all


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: get_nadv
    !> @brief Nonadiabatic coupling vector between two states, flat, in the coordinates of the last
    !! geometry passed to update_geometry.
    !> @details
    !! With a subset of the states of the model, the couplings of the adiabatic states are those
    !! within the states of the dynamics (see the module description).
    !----------------------------------------------------------------------------------------------
    function get_nadv(self, basis, istate1, istate2) result(nadv)
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        integer, intent(in) :: istate1, istate2
        real(dp), allocatable :: nadv(:)
        real(dp), allocatable :: d(:)
        integer :: imode, k1, k2
        real(dp) :: edif
        logical :: clamped

        allocate(d(self%ndof), source=0.0_dp)
        clamped = .false.
        select case(basis)
        case('diabatic')
            continue
        case('adiabatic')
            edif = clamped_gap(self%adiab_w(istate2, istate2) - self%adiab_w(istate1, istate1),   &
            &                  clamped)
            do imode = 1, self%ndof
                d(imode) = self%adiab_dw(imode, istate1, istate2) / edif
            end do
        case('mch')
            k1 = self%sel_mch(istate1)
            k2 = self%sel_mch(istate2)
            edif = clamped_gap(self%mch_w(k2, k2) - self%mch_w(k1, k1), clamped)
            do imode = 1, self%ndof
                d(imode) = self%mch_dw(imode, k1, k2) / edif
            end do
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_nadv function.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select
        if (clamped .and. stdp1) then
            write(stderr, *) 'Warning in evaluator_diabatic_mod, get_nadv function.'
            write(stderr, *) '  Near-degeneracy between states ', istate1, ' and ', istate2, '.'
            write(stderr, *) '  Setting denominator to ', gap_tol, ' Hartree.'
        end if
        nadv = self%get_derivative_target_coord(d)
    end function get_nadv


    !> Energy gap for the denominator of a derivative coupling: edif, or gap_tol with the sign of
    !! edif if its magnitude is smaller (clamped).
    function clamped_gap(edif, clamped) result(gap)
        real(dp), intent(in) :: edif
        logical, intent(out) :: clamped
        real(dp) :: gap

        clamped = (abs(edif) < gap_tol)
        gap = edif
        if (clamped) gap = sign(gap_tol, edif)
    end function clamped_gap


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: mch_soc
    !> @brief Spin-orbit couplings of the model in its MCH states, with all spin components.
    !> @details
    !! The same vectors of mch_trans transform each component of a spin-free state.
    !----------------------------------------------------------------------------------------------
    function mch_soc(self) result(soc)
        use matrix_mod, only : unitary_transform
        class(diabatic_evaluator), intent(in) :: self
        real(dp), allocatable :: soc(:, :)

        soc = self%soc_w
        call unitary_transform(self%model_states%expand(self%mch_trans), soc, trans=.true.)
    end function mch_soc


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: get_soc2
    !> @brief Squared magnitudes of the spin-orbit couplings between all states of a basis.
    !> @details
    !! In the 'mch' basis, element (a, b) is the sum of the squared couplings between all
    !! components of the spin-free states a and b. For two manifolds of components, this is the
    !! root-sum-square coupling, squared. It depends neither on the signs of mch_trans nor on the
    !! basis of the components.
    !----------------------------------------------------------------------------------------------
    function get_soc2(self, basis) result(s2)
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        real(dp), allocatable :: s2(:, :)

        select case(basis)
        case('diabatic')
            s2 = abs(self%soc_w)**2
        case('mch')
            s2 = self%model_states%sum_components(abs(self%mch_soc())**2)
            s2 = s2(self%sel_mch, self%sel_mch)
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, get_soc2 function.'
            write(stderr, *) '  Unsupported basis: ', basis
            stop
        end select
    end function get_soc2


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: match_phase
    !> @brief Match the phases of the vectors in ket to those in bra (phasematch_assigned_rotation).
    !> @details
    !! The overlaps bra^T ket between the vectors, after the phases of ket are matched, are
    !! returned in olap.
    !----------------------------------------------------------------------------------------------
    subroutine match_phase(bra, ket, olap)
        real(dp), intent(in) :: bra(:, :)
        real(dp), intent(inout) :: ket(:, :)
        real(dp), allocatable, intent(out), optional :: olap(:, :) !< bra^T ket.
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
        if (present(olap)) call move_alloc(overlap, olap)
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


    !> Eigenvectors of a symmetric matrix and, with e, its eigenvalues in ascending order.
    function get_eigvec(mat, e) result(eigv)
        use linalg_wrapper_mod, only : syev
        real(dp), intent(in) :: mat(:, :)
        real(dp), intent(out), optional :: e(:)
        real(dp), allocatable :: eigv(:, :)
        real(dp), allocatable :: wrk_e(:)

        eigv = mat
        allocate(wrk_e(size(eigv, 1)))
        call syev(eigv, wrk_e, jobz='V', uplo='U')
        if (present(e)) e = wrk_e
    end function get_eigvec


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: in_basis
    !> @brief Matrix over the diabatic states of the model, in the states of another basis.
    !> @details
    !! In the adiabatic basis, the result is U^T mat U for the states of the dynamics.
    !----------------------------------------------------------------------------------------------
    function in_basis(self, basis, mat) result(res)
        class(diabatic_evaluator), intent(in) :: self
        character(len=*), intent(in) :: basis
        real(dp), intent(in) :: mat(:, :) !< Matrix over the diabatic states.
        real(dp), allocatable :: res(:, :)

        select case(basis)
        case('diabatic')
            res = mat
        case('adiabatic')
            res = matmul(matmul(transpose(self%adiab_trans), mat), self%adiab_trans)
        case default
            write(stderr, *) 'Error in evaluator_diabatic_mod, in_basis function.'
            write(stderr, *) '  Unknown basis: ', basis
            stop
        end select
    end function in_basis



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

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: diabatic_save_state
    !> @brief Add the evaluated Hamiltonians, their derivatives and the transformations between
    !! the bases to a JSON object, in the object 'diabatic'.
    !> @details
    !! All evaluated arrays are saved, so a loaded point is the same as the evaluated one without
    !! evaluating the potential again. Arrays which are not allocated (e.g. the bases of the
    !! source of effective spin states, which are never built) are left out.
    !----------------------------------------------------------------------------------------------
    subroutine diabatic_save_state(self, json, obj)
        use json_module, only : json_core, json_value
        class(diabatic_evaluator), intent(in) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object to which the state is added.
        type(json_value), pointer :: d

        call evaluator_save_state(self, json, obj)
        call json%create_object(d, 'diabatic')
        call json%add(d, 'ndof', self%ndof)
        call json%add(d, 'exact_gradient', self%exact_gradient)
        call json%add(d, 'soc_checked', self%soc_checked)
        call save_2d(d, 'sf_w', self%sf_w)
        call save_3d(d, 'sf_dw', self%sf_dw)
        call save_2d(d, 'soc_w', self%soc_w)
        call save_3d(d, 'soc_dw', self%soc_dw)
        call save_2d(d, 'diab_w', self%diab_w)
        call save_3d(d, 'diab_dw', self%diab_dw)
        call save_2d(d, 'adiab_w', self%adiab_w)
        call save_3d(d, 'adiab_dw', self%adiab_dw)
        call save_2d(d, 'adiab_trans', self%adiab_trans)
        call save_2d(d, 'mch_trans', self%mch_trans)
        call save_2d(d, 'mch_w', self%mch_w)
        call save_3d(d, 'mch_dw', self%mch_dw)
        call json%add(obj, d)

    contains

        subroutine save_2d(d, name, a)
            use json_util_mod, only : json_add_array
            type(json_value), pointer, intent(in) :: d
            character(len=*), intent(in) :: name
            real(dp), allocatable, intent(in) :: a(:, :)

            if (allocated(a)) call json_add_array(json, d, name, a)
        end subroutine save_2d

        subroutine save_3d(d, name, a)
            use json_util_mod, only : json_add_array
            type(json_value), pointer, intent(in) :: d
            character(len=*), intent(in) :: name
            real(dp), allocatable, intent(in) :: a(:, :, :)

            if (allocated(a)) call json_add_array(json, d, name, a)
        end subroutine save_3d

    end subroutine diabatic_save_state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: diabatic_load_state
    !> @brief Restore the evaluated arrays saved by save_state.
    !> @details
    !! Arrays missing from the backup are deallocated. The shape of each array is checked
    !! against the numbers of states of the model and of the dynamics and the degrees of freedom
    !! of the current setup. Stops if the gradients were computed with another exact_gradient.
    !----------------------------------------------------------------------------------------------
    subroutine diabatic_load_state(self, json, obj)
        use json_module, only : json_core, json_value
        use json_util_mod, only : json_get_required, json_get_child, json_check_value
        class(diabatic_evaluator), intent(inout) :: self
        type(json_core), intent(inout) :: json
        type(json_value), pointer, intent(in) :: obj !< Object holding the state.
        type(json_value), pointer :: d
        integer :: ns, nm, nsf, nd

        call evaluator_load_state(self, json, obj)
        d => json_get_child(json, obj, 'diabatic')
        call json_check_value(json, d, 'ndof', self%ndof)
        call json_check_value(json, d, 'exact_gradient', self%exact_gradient)
        call json_get_required(json, d, 'soc_checked', self%soc_checked)
        ns = self%states%n_state
        nm = self%model_states%n_state
        nsf = self%model_states%n_sf
        nd = self%ndof
        call load_2d('sf_w', self%sf_w, [nsf, nsf])
        call load_3d('sf_dw', self%sf_dw, [nd, nsf, nsf])
        call load_2d('soc_w', self%soc_w, [nm, nm])
        call load_3d('soc_dw', self%soc_dw, [nd, nm, nm])
        call load_2d('diab_w', self%diab_w, [nm, nm])
        call load_3d('diab_dw', self%diab_dw, [nd, nm, nm])
        call load_2d('adiab_w', self%adiab_w, [ns, ns])
        call load_3d('adiab_dw', self%adiab_dw, [nd, ns, ns])
        call load_2d('adiab_trans', self%adiab_trans, [nm, ns])
        call load_2d('mch_trans', self%mch_trans, [nsf, nsf])
        call load_2d('mch_w', self%mch_w, [nsf, nsf])
        call load_3d('mch_dw', self%mch_dw, [nd, nsf, nsf])

    contains

        subroutine load_2d(name, a, expected)
            character(len=*), intent(in) :: name
            real(dp), allocatable, intent(inout) :: a(:, :)
            integer, intent(in) :: expected(2)

            if (allocated(a)) deallocate(a)
            if (json%valid_path(d, name)) call json_get_required(json, d, name, a, expected)
        end subroutine load_2d

        subroutine load_3d(name, a, expected)
            character(len=*), intent(in) :: name
            real(dp), allocatable, intent(inout) :: a(:, :, :)
            integer, intent(in) :: expected(3)

            if (allocated(a)) deallocate(a)
            if (json%valid_path(d, name)) call json_get_required(json, d, name, a, expected)
        end subroutine load_3d

    end subroutine diabatic_load_state


end module evaluator_diabatic_mod
